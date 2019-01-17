cloud crushing version of Mike


Configure

```
 # MHD
 ./configure --with-problem=cloud-3dtang --with-gas=mhd --with-order=3 --enable-mpi --with-flux=hlld --with-nscalars=1 

 # Hydrodynamical
 ./configure --with-problem=cloud-3dtang --with-order=2p --enable-mpi --with-flux=hllc --with-nscalars=1 --with-gas=hydro --with-integrator=vl
```


