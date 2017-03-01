# gasgiant

## Gas giant models

`data/grid3S` is a table of gas giant models from [Marleau & Cumming 2014](http://adsabs.harvard.edu/abs/2014MNRAS.437.1378M)

## Envelope Models

`envelope.py` calculates plane-parallel models of accreting envelopes (see section 3 of [Berardo et al. 2017](http://adsabs.harvard.edu/abs/2017ApJ...834..149B))

`python envelope.py L(LSun) Mdot(ME/yr) Ptop(bars) Ttop(K) plot_flag mass radius Starget`

If the first parameter `L` is negative, the code will search over luminosity to match the entropy `Starget` at the base of the envelope, otherwise `Starget` is ignored and the envelope is integrated with surface luminosity `L`

