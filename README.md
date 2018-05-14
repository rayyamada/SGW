# SGW
MATLAB code for surface-gravity wave topography interactions

Integrates surface gravity wave topography interaction equations (see Thomas and Yamada, 2018)
using the High-Order Spectral (HOS) Method.

To run the code, 
1) setup simulation parameters  (e.g. run script sample_HOS_init.m)
2) run HOS_spectral_ps.m to integrate forward in time (uses RK4 integration)

HOS_spectral_ps.m depends on HOS_getrhs_ps.m and HOS_getPhiz.m which depend on
k2ps.m and ps2k.m


