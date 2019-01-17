cloud crushing version of Mike


# Configure

```
 # MHD
 ./configure --with-problem=cloud-3dtang --with-gas=mhd --with-order=3 --enable-mpi --with-flux=hlld --with-nscalars=1 

 # Hydrodynamical
 ./configure --with-problem=cloud-3dtang --with-order=2p --enable-mpi --with-flux=hllc --with-nscalars=1 --with-gas=hydro --with-integrator=vl
```

# Example input file

```
<job>
maxout          = 3         # Number of output files
num_domains     = 1
problem_id      = cloud

<output1>
dt              = 0.5
out_fmt         = hst

<output2>
dt              = 1.0
out             = prim
out_fmt         = vtk

<output3>
dt              = 10.0
out_fmt         = rst

<time>
cour_no         = 0.3
nlim            = 2000000.0
tlim            = 950.0

<domain1>
NGrid_x1        = 32.0
NGrid_x2        = 4.0
NGrid_x3        = 4.0
Nx1             = 480.0
Nx2             = 96.0
Nx3             = 96.0
bc_ix2          = 2
bc_ix3          = 2
bc_ox1          = 2.0
bc_ox2          = 2
bc_ox3          = 2
level           = 0.0
x1max           = 135.0
x1min           = -15.0
x2max           = 15.0
x2min           = -15.0
x3max           = 15.0
x3min           = -15.0

<problem>
Bz              = 0.365   # not needed in HD
betafloor       = 0.00   # not needed in HD1
betain          = 10.0   # not needed in HD
betaout         = 1e+20   # not needed in HD
d_MIN           = 0.001 
dp              = -0.6 # sets global pressure. P = Gamma - 1 + dp
dr              = 0.046875 # border of cloud in units of r_cl
drat            = 100.0 # density ratio
dtmin           = 1e-09 
gamma           = 1.66666666667
nterms          = 10.0    # not needed in HD
r_cloud         = 2.5 # radius of cloud
tnotcool        = 0.04 # switch off cooling above this T
vflow           = 0.5 # flow speed

<par_end>

```



