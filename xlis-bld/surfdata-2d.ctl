* use xdfopen
DSET /discover/nobackup/projects/nca/ytian/CESM_INPUT_DATA/lnd/clm2/surfdata_map/surfdata_0.9x1.25_simyr1850_c130415.nc
DTYPE netcdf
UNDEF  1.e+36
XDEF lsmlon     288 LINEAR    0   1.25
YDEF lsmlat     192 LINEAR    -90    0.942
TDEF time 1 LINEAR 00Z01Jan1901 1mo
VARS 2
PFTDATA_MASK=>pftdata_mask 0 y,x
TOPOT=>topo 0 y,x
ENDVARS

