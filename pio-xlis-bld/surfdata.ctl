* use xdfopen
DSET /discover/nobackup/projects/nca/ytian/CESM_INPUT_DATA/lnd/clm2/surfdata_map/surfdata_0.9x1.25_simyr1850_c130415.nc
DTYPE netcdf
UNDEF  1.e+36
XDEF lsmlon     288 LINEAR    0   1.25
YDEF lsmlat     192 LINEAR    -90    0.942
ZDEF lsmpft 17 LINEAR 1 1 
TDEF time 1 LINEAR 00Z01Jan1901 1mo
VARS 3
PCT_PFT=>pct_pft 17 z,y,x
MONTHLY_LAI=>monthly_lai 17 t,z,y,x
TOPOT=>topo 0 y,x
ENDVARS

