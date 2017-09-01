# C-O-H
Code to calculate chemical speciation of an ideal C&ndash;O&ndash;H fluid in equilibrium with graphite at variable temperature, pressure, and O<sub>2</sub> fugacity.  Intended to simulate a carbon-saturated crustal fluid that might exist at depth in melt-derived gabbroic rocks beneath mid-ocean ridges.

Running the scripts `crustalfluidmodel_optim.R` and `plotresults2_FMQ.m` produces [Figure 3.7](https://github.com/dtexwang/thesis/blob/master/figures/Fig3.S1.pdf) in [my Ph.D. thesis](http://dx.doi.org/10.1575/1912/9052).

**System**
* R version 3.3.3 (x86_64) on Windows 7, with packages:
  * [CHNOSZ_1.0.8](http://www.chnosz.net/)
  * [rootSolve_1.7](https://cran.r-project.org/web/packages/rootSolve/index.html)
* MATLAB R2012b (8.0.0.783) 64-bit Windows 7
