"""fit_cooling_data_dens.py

Create header files required for the Townsend+09 cooling algorithm implemented
in Athena. Will read in cooling function table, fit it (uniformly, i.e., non
optimally), and output header files as well as some plots.
"""
import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def fitloglog(x, y):
    x1 = np.log(x)
    y1 = np.log(y)
    p = np.polyfit(x1, y1, 1)
    A = p[0]
    B = np.exp(p[1])
    return A, B

def fit2point(x, y):
    x1 = np.log(x)
    y1 = np.log(y)
    A = (y1[1] - y1[0]) / (x1[1] - x1[0])
    B = np.exp(y1[1] - A * x1[1])
    return A, B

def plotfit(T0, T1, Lambda, index, **kwargs):
    x = np.array([T0, T1])
    y = Lambda * x**index
    plt.loglog(x, y, **kwargs)


########################################
# Configure here
out_prefix = "WSS09_z0_Z1/"
out_fn = "WSS09_z0_Z1.h"

# Reading in cooling data
f = h5py.File("WSS09_CoolingTables/z_0.000.hdf5")
d = f['Solar/Hydrogen_density_bins'][...]     # Densities
T = f['Solar/Temperature_bins'][...]          # Temperatures
Lambda_full = f['Solar/Net_cooling'][...].T   # Lambda(rho, T)
print("Loaded Lambda data with shape `%s`" %(str(Lambda_full.shape)))

########################################
s = out_prefix + "plots/"
if not os.path.exists(s):
    os.makedirs(s)


#Tmin = T[0]
Tmin = 5e3
Tmax = 1e8

binnum_T = 31
binnum_d = 81


# fitdata[:,0] -- temperature bins
# fitdata[:,1] -- Lambda
# fitdata[:,2] -- index
fitdata_full = np.zeros((binnum_d,binnum_T+1, 3))


T_bin = np.logspace(np.log10(Tmin), np.log10(Tmax), binnum_T+1,
	                endpoint=True, base=10.0)
#d_bin = np.logspace(np.log10(d[0]), np.log10(d[-1]), binnum_d+1,
#	                endpoint=True, base=10.0)
d_bin = d.copy()

assert len(d_bin) == binnum_d
assert len(T_bin) == binnum_T + 1

for i in tqdm(range(0,len(d))):
    plt.clf()

    fitdata = fitdata_full[i,:,:]
    Lambda = Lambda_full[i,:].copy()
    minL = 1e-30

    ## Try to make cut smoother --> failed!
    if Lambda[(T > Tmin) & (T < Tmax)].min() < minL:
        mclip = Lambda < minL
        Lambda = np.clip(Lambda, minL, np.inf)
        Ta = T_bin[(T_bin < T[mclip][-1])][-1]
        Tb = T_bin[(T_bin > T[mclip][-1])][0]
        La, Lb = [Lambda[np.min(np.abs(T - cT)) == np.abs(T - cT)] for cT in [Ta, Tb]]
        cm = (T >= Ta) & (T<=Tb)
        newv = (La - Lb) / (Ta - Tb) * (T[cm] - Ta) + La
        Lambda[cm] = newv


    fitdata[:,0] = T_bin

    plt.loglog(T, Lambda, 'k:')

    for n in range(binnum_T):
        if n % 2 == 0:
            mask = np.logical_and((T >= T_bin[n]), (T < T_bin[n+1]))
            index_fit, Lambda_fit = fitloglog(T[mask], Lambda[mask])
            fitdata[n,1] = Lambda_fit
            fitdata[n,2] = index_fit

            plotfit(T_bin[n], T_bin[n+1], Lambda_fit, index_fit)

    for n in range(binnum_T):
        if n % 2 == 1:
            T0, T1 = T_bin[n], T_bin[n+1]
            Lambda0 = fitdata[n-1,1] * T0**fitdata[n-1,2]
            Lambda1 = fitdata[n+1,1] * T1**fitdata[n+1,2]
            index_fit, Lambda_fit = fit2point(np.array([T0, T1]), 
                                              np.array([Lambda0, Lambda1]))
            #assert Lambda_fit
            fitdata[n,1] = Lambda_fit
            fitdata[n,2] = index_fit

            plotfit(T0, T1, Lambda_fit, index_fit)
    # Plot last bin again
    plotfit(T0, T1 * 5, Lambda_fit, index_fit, ls = '--', color = 'k')
    # Convert from Lambda * T**index to Lambda_k * (T / T_k)**index
    fitdata[:,1] *= fitdata[:,0]**fitdata[:,2]

    fitdata[np.isnan(fitdata)] = 9.949213290000138e-21

    # Note fitdata[binnum_T,1] = fitdata[binnum_T,2] = 0.
    # but the code will compute these based on powerlaw

    assert not np.any(np.isnan(fitdata)), "nan found!"

    # Change the units
    fitdata[:,0] *= 8.61733e-8 # K --> keV
    fitdata[:,1] /= 1e-23       # erg cm^3/s --> 1e-23 erg cm^3/s


    plt.xlabel(r"$T$ (code units)")
    plt.ylabel(r"$\Lambda_{\rm net}$ (code units)")
    plt.title(r"$n_{H} = %.3e\,{\rm cm}^{-3}$" %(d[i]))
    plt.savefig(out_prefix + "plots/%03d.png" %(i), bbox_inches = 'tight')


#  save
np.savetxt(out_prefix + 'cool_func.dat', fitdata)

# Now save for C
fp = open(out_prefix + out_fn, 'w')

fp.write("/* This file was automatically created. */\n\n")

fp.write("#define nfit_cool_T %d\n" %(binnum_T))
fp.write("#define nfit_cool_d %d\n\n" %(binnum_d))

fp.write("static Real sdT[nfit_cool_T] = {")
fp.write(", ".join(["%e" %(fitdata[i,0]) for i in range(binnum_T)]))
fp.write("};\n\n")

fp.write("static Real sdd[nfit_cool_d] = {")
fp.write(", ".join(["%e" %(i) for i in d_bin]))
fp.write("};\n\n")

fp.write("static const Real sdexpt[nfit_cool_d][nfit_cool_T] = {\n")
for idbin in range(binnum_d):
    fp.write("{")
    fp.write(", ".join(["%e" %(i) for i in fitdata_full[idbin,:binnum_T,2]]))
    fp.write("}")
    if idbin != binnum_d - 1:
        fp.write(",\n")
fp.write("};\n\n")

fp.write("static const Real sdL[nfit_cool_d][nfit_cool_T] = {\n")
for idbin in range(binnum_d):
    fp.write("{")
    fp.write(", ".join(["%e" %(i) for i in fitdata_full[idbin,:binnum_T,1]]))
    fp.write("}")
    if idbin != binnum_d - 1:
        fp.write(",\n")
fp.write("};\n\n")

fp.close()


print("Done. %s written." %(out_prefix + out_fn))


