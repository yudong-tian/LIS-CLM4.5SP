#! /bin/bash 

source /usr/share/modules/init/bash
module purge
module load comp/intel-14.0.3.174 mpi/impi-5.0.3.048 other/cmake-2.8.11.2

mpif90  -c -I.  \
  -I/discover/nobackup/ytian/CESM_CODE/LIS-CLM4.5SP/bld/intel/impi/nodebug/nothreads/include \
  -I/discover/nobackup/ytian/CESM_CODE/LIS-CLM4.5SP/bld/intel/impi/nodebug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share \
  -I/usr/local/other/SLES11.1/netcdf/4.1.3/intel-12.1.0.233/include \
  -I/discover/nobackup/ytian/CESM_CODE/LIS-CLM4.5SP/bld/intel/impi/nodebug/nothreads/include \
  -I/home/ytian/CESM/LIS-CLM4.5SP/bld/lnd/obj \
  -I/home/ytian/CESM/LIS-CLM4.5SP/bld/atm/obj \
  -I/home/ytian/CESM/LIS-CLM4.5SP/bld/ocn/obj \
  -I/home/ytian/CESM/LIS-CLM4.5SP/bld/cpl/obj \
  -I/discover/nobackup/ytian/CESM_CODE/LIS-CLM4.5SP/bld/lib/include  \
  -O2 -fp-model precise -convert big_endian -assume byterecl -ftz -traceback -DMAXPATCH_PFT=17 -DLINUX  \
  -DNDEBUG -DMCT_INTERFACE -DHAVE_MPI -DFORTRANUNDERSCORE -DNO_R16 -DFORTRANUNDERSCORE -DNO_R16 \
  -DHAVE_F2008_CONTIGUOUS -DLINUX -DCPRINTEL  -DHAVE_SLASHPROC -free \
  xlis_main.F90  \
  my_code/clm_initializeMod.F90 \
  my_code/surfrdMod.F90 \
  my_code/ncdio_pio.F90 


mpif90  -o xlis_main xlis_main.o clm_initializeMod.o surfrdMod.o ncdio_pio.o \
  -L/discover/nobackup/ytian/CESM_CODE/LIS-CLM4.5SP/bld/lib/  \
  -latm -lice -llnd -locn -lrof -lglc -lwav \
  -L/discover/nobackup/ytian/CESM_CODE/LIS-CLM4.5SP/bld/intel/impi/nodebug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share \
  -lcsm_share \
  -L/discover/nobackup/ytian/CESM_CODE/LIS-CLM4.5SP/bld/intel/impi/nodebug/nothreads/lib -lpio \
  -lgptl -lmct -lmpeu \
  -L/usr/local/other/SLES11.1/netcdf/4.1.3/intel-12.1.0.233/lib \
  -lnetcdff  -lnetcdf -lm -lcurl  \
  -L/usr/local/other/SLES11/lapack/3.3.1/intel-12.1.0.233/lib  \
  -llapack -lblas

