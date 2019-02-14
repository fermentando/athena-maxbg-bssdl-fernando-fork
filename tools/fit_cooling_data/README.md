
# fit_cooling_data_dens.py

Create header files required for the Townsend+09 cooling algorithm implemented
in Athena. Will read in cooling function table, fit it (uniformly, i.e., non
optimally), and output header files as well as some plots.

I used, for instance, the cooling tables provided by (Wiersma, R. P. C., Schaye,
J., and Smith, B. D. 2009; arxiv:0807.3748;
http://www.strw.leidenuniv.nl/WSS08/) which the script is currently configured
to run on. But any cooling table should work.

The generated template file should go in `src/prob/cooling_data/` and can be
included in the problem file.


