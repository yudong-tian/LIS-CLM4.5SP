  program single_col
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use abortutils       , only : endrun
    use clm_time_manager , only : get_nstep, get_step_size, set_timemgr_init, &
                                  set_nextsw_cday
    use clmtype
    use clm_atmlnd       , only : clm_l2a, clm_a2l
    use clm_glclnd       , only : clm_s2x
    use clm_initializeMod, only : initialize1, initialize2
    use clm_varctl       , only : finidat,single_column, set_clmvarctl, iulog, noland, &
                                  inst_index, inst_suffix, inst_name, &
                                  create_glacier_mec_landunit 

! CLM driver related 
    use clm_driver      ,only : clm_drv
    use clm_varorb       , only : eccen, obliqr, lambm0, mvelpp
    use filterMod           , only : filter
    use dynlandMod          , only : dynland_hwcontent
  use inicPerpMod         , only : inicPerp
  use accFldsMod          , only : updateAccFlds
  use clm_driverInitMod   , only : clm_driverInit
  use BalanceCheckMod     , only : BeginWaterBalance, BalanceCheck
  use SurfaceRadiationMod , only : SurfaceRadiation
  use Hydrology1Mod       , only : Hydrology1
  use Hydrology2Mod       , only : Hydrology2
  use SLakeFluxesMod   , only : SLakeFluxes
  use SLakeTemperatureMod, only : SLakeTemperature
  use SLakeHydrologyMod, only : SLakeHydrology
  use Biogeophysics1Mod   , only : Biogeophysics1
  use BareGroundFluxesMod , only : BareGroundFluxes
  use CanopyFluxesMod     , only : CanopyFluxes
  use Biogeophysics2Mod   , only : Biogeophysics2
  use SurfaceAlbedoMod    , only : SurfaceAlbedo
  use pft2colMod          , only : pft2col
  use ActiveLayerMod      , only : alt_calc
  use UrbanMod            , only : UrbanAlbedo, UrbanRadiation, UrbanFluxes
  use clm_atmlnd          , only : clm_map2gcell
  use SNICARMod           , only : SnowAge_grain
  use STATICEcosysDynMod  , only : EcosystemDyn

    use clm_varcon          , only : zlnd, isturb
    use controlMod       , only : control_setNL
    use decompMod        , only : get_proc_bounds, get_clump_bounds
    use domainMod        , only : ldomain
    use shr_pio_mod      , only : shr_pio_init1, shr_pio_init2
    use shr_kind_mod     , only : r8 => shr_kind_r8
    use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                  shr_file_getLogUnit, shr_file_getLogLevel, &
                                  shr_file_getUnit, shr_file_setIO
    use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod  , only : seq_timemgr_EClockGetData
    use seq_infodata_mod , only : seq_infodata_type, seq_infodata_GetData, seq_infodata_PutData, &
                                  seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                  seq_infodata_start_type_brnch
    use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
    use spmdMod          , only : masterproc, spmd_init, npes, iam
    use clm_varctl       , only : nsrStartup, nsrContinue, nsrBranch
    use clm_cpl_indices  , only : clm_cpl_indices_set, nflds_l2x
    use seq_flds_mod
    use mct_mod
    use ESMF

    implicit none
  include 'mpif.h' 
!
!
! !LOCAL VARIABLES:
    integer                          :: ierr 	     ! YDT 
    integer                          :: LNDID	     ! Land identifyer
    integer                          :: mpicom_lnd   ! MPI communicator
    type(mct_gsMap),         pointer :: GSMap_lnd    ! Land model MCT GS map
    type(mct_gGrid),         pointer :: dom_l        ! Land model domain
    type(mct_gsMap),         pointer :: GSMap_sno
    type(mct_gGrid),         pointer :: dom_s
    type(seq_infodata_type), pointer :: infodata     ! CESM driver level info data
    integer  :: lsize                                ! size of attribute vector
    integer  :: g, l, c, p, i,j, nc                  ! indices
    integer  :: dtime_sync                           ! coupling time-step from the input synchronization clock
    integer  :: dtime_clm                            ! clm time-step
    logical  :: exists                               ! true if file exists
    logical  :: atm_aero                             ! Flag if aerosol data sent from atm model
    logical  :: samegrid_al                          ! true if atmosphere and land are on the same grid
    real(r8) :: scmlat                               ! single-column latitude
    real(r8) :: scmlon                               ! single-column longitude
    character(len=SHR_KIND_CL) :: caseid             ! case identifier name
    character(len=SHR_KIND_CL) :: ctitle             ! case description title
    character(len=SHR_KIND_CL) :: starttype          ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL) :: calendar           ! calendar type name
    character(len=SHR_KIND_CL) :: hostname           ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version            ! Model version
    character(len=SHR_KIND_CL) :: username           ! user running the model
    integer :: nsrest                                ! clm restart type
    integer :: perpetual_ymd                         ! perpetual date
    integer :: ref_ymd                               ! reference date (YYYYMMDD)
    integer :: ref_tod                               ! reference time of day (sec)
    integer :: start_ymd                             ! start date (YYYYMMDD)
    integer :: start_tod                             ! start time of day (sec)
    integer :: stop_ymd                              ! stop date (YYYYMMDD)
    integer :: stop_tod                              ! stop time of day (sec)
    logical :: brnch_retain_casename                 ! flag if should retain the case name on a branch start type
    logical :: perpetual_run                         ! flag if should cycle over a perpetual date or not
    integer :: lbnum                                 ! input to memory diagnostic
    integer :: begg, endg
    integer :: begl, endl
    integer :: begc, endc
    integer :: begp, endp
    character(len=32), parameter :: sub = 'lnd_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"


!YDT local variables for clm_drv() 

    logical :: doalb                      ! .true. ==> do albedo calculation on this time step
    logical :: rstwr_sync                 ! .true. ==> write restart file before returning
    logical :: rstwr                      ! .true. ==> write restart file before returning
    logical :: nlend_sync                 ! Flag signaling last time-step
    logical :: nlend                      ! .true. ==> last time-step
    logical :: dosend                     ! true => send data back to driver
    real(r8):: nextsw_cday                ! calday from clock of next radiation computation
    real(r8):: caldayp1                   ! clm calday plus dtime offset
    integer :: shrlogunit,shrloglev       ! old values for share log unit and log level
    logical,save :: first_call = .true.         ! first call work
    logical  :: glcrun_alarm          ! if true, sno data is averaged and sent to glc this step
    logical  :: update_glc2sno_fields ! if true, update glacier_mec fields
    real(r8) :: calday                ! calendar day for nstep
    real(r8) :: declin                ! solar declination angle in radians for nstep
    real(r8) :: declinp1              ! solar declination angle in radians for nstep+1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
    real(r8) :: recip                 ! reciprical
    character(len=32)            :: rdate       ! date char string for restart file names

    integer :: num_inst_total = 1  ! land only

  integer , pointer :: clandunit(:) ! landunit index associated with each column
  integer , pointer :: itypelun(:)  ! landunit type
!-----------------------------------------------------------------------

    LNDID = 1

    write(*, *) "MPI initialization starts ..."
    call MPI_INIT(ierr)
    if (ierr .ne. MPI_SUCCESS) then
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, ierr)
    end if

    mpicom_lnd = MPI_COMM_WORLD
    call spmd_init(mpicom_lnd, LNDID) 
    write(*, *) "MPI initialization done ..., npes=", npes

    call shr_pio_init1(num_inst_total, "drv_in", mpicom_lnd) 
    !call shr_pio_init2(comp_id,comp_name,comp_iamin,comp_comm,comp_comm_iam)
    ! YDT comp_name has to be upcase
     call shr_pio_init2((/1/),      (/"LND"/),    (/.true./), (/mpicom_lnd/), (/iam/) )

    call ESMF_Initialize() 

    nsrest=0  !0: initial run. 1: restart: 3: branch
    version = "clm 4.5" 
    hostname = "discover" 
    username = "ytian" 
    single_column = .false.
    inst_name ="lnd" 


    call control_setNL("lnd_in") 

    ! got settings from drv_in
    call set_timemgr_init( calendar_in='NO_LEAP', start_ymd_in=20000101, start_tod_in=0, &
                           ref_ymd_in=20000101, ref_tod_in=0, stop_ymd_in=99990101,         &
                           stop_tod_in=0,  perpetual_run_in=.false.,     &
                           perpetual_ymd_in=-999 )

    call set_clmvarctl(    caseid_in=caseid, ctitle_in=ctitle,                     &
                           brnch_retain_casename_in=brnch_retain_casename,         &
                           single_column_in=single_column, scmlat_in=scmlat,       &
                           scmlon_in=scmlon, nsrest_in=nsrest, version_in=version, &
                           hostname_in=hostname, username_in=username)

    doalb = .false.   !YDT  not do alb for now

    write(*, *) "set_clmvarctl() done ..."


    ! Read namelist, grid and surface data
    call initialize1( )

    ! the following are called within initialize1() 
    ! call initClmtype()
    ! call init_atm2lnd_type(begg, endg, clm_a2l)
    ! call init_lnd2atm_type(begg, endg, clm_l2a)

    write(*, *) "initialize1() done ..."

    ! Finish initializing clm
    call initialize2()

    write(*, *) "initialize2() done ..."

    call get_clump_bounds(1, begg, endg, begl, endl, begc, endc, begp, endp)



    !YDT pick one grid box at (35N, 95W) 
    !do g=begg, endg
      ! write(*, '(3I7, 2F10.2)') g, grc%luni(g), grc%lunf(g), grc%latdeg(g), grc%londeg(g)  
      ! at (35N, 95W) 
      ! 12762  19878  19880     35.34    265.00
    !end do 

     write(*, *) "==================  Single grid decomposition ================================" 
     g=12762
     nc = 288  ! number of columns in the original grid 
     write(*, *) " g       gindex   ic    ir   lat   lon   nlandunits  ncolumns  npfts" 
     write(*, *) "---------------------------------------------------------------------------" 
     write(*, '(4I7, 2F8.2, 3I7)') g, grc%gindex(g), mod(grc%gindex(g), nc), &
                  (grc%gindex(g) - mod(grc%gindex(g), nc) )/nc + 1, & 
                 grc%latdeg(g), grc%londeg(g), & 
            grc%nlandunits(g), grc%ncolumns(g), grc%npfts(g) 
 
     write(*, *) 
     write(*, *) " g    gindex    luni    lunf   coli  colf   pfti   pftf" 
     write(*, *) "---------------------------------------------------------------------------" 
     write(*, '(9I7)') g, grc%gindex(g), grc%luni(g), grc%lunf(g), grc%coli(g), grc%colf(g), &
                          grc%pfti(g), grc%pftf(g)

     write(*, *) 
     write(*, *) " l    gridcell  wtgcell%  coli   colf   pfti   pftf     npfts   itype_lun" 
     write(*, *) "---------------------------------------------------------------------------" 
     do l=grc%luni(g), grc%lunf(g)      
       write(*, '(2I7, F10.2, 6I7)') l, lun%gridcell(l), lun%wtgcell(l)*100, lun%coli(l), lun%colf(l), &
                                        lun%pfti(l), lun%pftf(l), lun%npfts(l), lun%itype(l) 
     end do 
                            
     write(*, *) 
     write(*, *) " c    landunit  wtlunit%  gridcell  wtgcell%  pfti   pftf  npfts   itype_col" 
     write(*, *) "---------------------------------------------------------------------------" 
     do c=grc%coli(g), grc%colf(g)      
       write(*, '(2I7, F10.2, I7, F10.2, 3X, 4I7)') c, col%landunit(c), col%wtlunit(c)*100, col%gridcell(c), & 
                 col%wtgcell(c)*100, col%pfti(c), col%pftf(c), col%npfts(c), col%itype(c) 
     end do 
    
     write(*, *) 
     write(*, *) " p    column     wtcol%  landunit  wtlunit%  gridcell  wtgcell%  itype_pft" 
     write(*, *) "---------------------------------------------------------------------------" 
     do p=grc%pfti(g), grc%pftf(g)      
       write(*, '(2I7, F10.2, I7, F10.2, I7, F10.2, I7)') p, pft%column(p), pft%wtcol(p)*100, pft%landunit(p), &
                                        pft%wtlunit(p)*100, &
                                        pft%gridcell(p), pft%wtgcell(p)*100, pft%itype(p) 
     end do 

     write(*, *) "===========================================================================" 

     write(*, *) 
     write(*, *) "================== Filter Settings ========================================" 
     write(*, *) " num_soilc num_soilp num_hydroc num_lakec num_lakep num_nolakec num_nolakep" 
     write(*, '(7I9)') filter(1)%num_soilc, filter(1)%num_soilp, filter(1)%num_hydrologyc, &
                       filter(1)%num_lakec, filter(1)%num_lakep, filter(1)%num_nolakec, filter(1)%num_nolakep
     write(*, *) "===========================================================================" 


    !===========  ranges for a single grid box  ------------------------------------
     nc = 1    ! clump index: always 1 for now 

     begg = g
     endg = g

     begl=grc%luni(g)
     endl=grc%lunf(g)      

     begc=grc%coli(g)
     endc=grc%colf(g)      

     begp=grc%pfti(g)
     endp=grc%pftf(g)      
     


     ! Assign local pointers to derived subtypes components (column-level)
     itypelun            => lun%itype
     clandunit           =>col%landunit  

    !===========  grid-level initialization  ------------------------------------
     
     ! initialize heat and water content and dynamic balance fields to zero, for the single grid box 
        gwf%qflx_liq_dynbal(g) = 0._r8
        gws%gc_liq2(g)         = 0._r8
        gws%gc_liq1(g)         = 0._r8
        gwf%qflx_ice_dynbal(g) = 0._r8
        gws%gc_ice2(g)         = 0._r8
        gws%gc_ice1(g)         = 0._r8
        gef%eflx_dynbal(g)     = 0._r8
        ges%gc_heat2(g)        = 0._r8
        ges%gc_heat1(g)        = 0._r8

    ! initialize input for the original clm_drv arguments
    ! subroutine clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
    ! logical ,        intent(in) :: doalb       ! true if time for surface albedo calc
    ! real(r8),        intent(in) :: nextsw_cday ! calendar day for nstep+1
    ! real(r8),        intent(in) :: declinp1    ! declination angle for next time step
    ! real(r8),        intent(in) :: declin      ! declination angle for current time step
    ! logical,         intent(in) :: rstwr       ! true => write restart file this step
    ! logical,         intent(in) :: nlend       ! true => end of run on this step
    ! character(len=*),intent(in) :: rdate       ! restart file time stamp for name
 
    doalb = .false. 
    nextsw_cday = -999 ! to be fixed 
    declinp1 = 0.0
    declin = 0.0
    rstwr = .false.  ! writing restart or not 
    nlend = .false.  ! not end of run 
    rdate = ""    ! not save restart for now

    ! Now ripping components out of clm_drv(). Skip 'CN' and 'CNDV" components for now 
    ! forget DUST and VOC too 

    !--- get initial heat,water content ---
       call dynland_hwcontent( begg, endg, gws%gc_liq1(begg:endg), &
                               gws%gc_ice1(begg:endg), ges%gc_heat1(begg:endg) )


    ! Determine decomp vertical profiles
        call alt_calc(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)

    ! Initialize the mass balance checks: water
     call BeginWaterBalance(begc, endc, begp, endp, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec, &
          filter(nc)%num_hydrologyc, filter(nc)%hydrologyc)

     ! Initialize variables from previous time step and
     ! Determine canopy interception and precipitation onto ground surface.
     ! Determine the fraction of foliage covered by water and the fraction
     ! of foliage that is dry and transpiring. Initialize snow layer if the
     ! snow accumulation exceeds 10 mm.

     pcf%cisun_z(begp:endp,:) = -999._r8
     pcf%cisha_z(begp:endp,:) = -999._r8

     ! initialize declination for current timestep
     do c = begc,endc
        cps%decl(c) = declin
     end do

     call clm_driverInit(begc, endc, begp, endp, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec)

     call Hydrology1(begc, endc, begp, endp, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_nolakep, filter(nc)%nolakep)


     ! Surface Radiation for non-urban columns

     call SurfaceRadiation(begp, endp, &
                           filter(nc)%num_nourbanp, filter(nc)%nourbanp)

     ! Surface Radiation for urban columns

     call UrbanRadiation(nc, begl, endl, begc, endc, begp, endp, &
                         filter(nc)%num_nourbanl, filter(nc)%nourbanl, &
                         filter(nc)%num_urbanl, filter(nc)%urbanl, &
                         filter(nc)%num_urbanc, filter(nc)%urbanc, &
                         filter(nc)%num_urbanp, filter(nc)%urbanp)

     ! Determine leaf temperature and surface fluxes based on ground
     ! temperature from previous time step.

     call Biogeophysics1(begg, endg, begc, endc, begp, endp, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)

     ! Determine bare soil or snow-covered vegetation surface temperature and fluxes
     ! Calculate Ground fluxes (frac_veg_nosno is either 1 or 0)

     ! BareGroundFluxes for all pfts except lakes and urban landunits

     call BareGroundFluxes(begp, endp, &
                           filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp)

     call UrbanFluxes(nc, begp, endp, begl, endl, begc, endc, &
                      filter(nc)%num_nourbanl, filter(nc)%nourbanl, &
                      filter(nc)%num_urbanl, filter(nc)%urbanl, &
                      filter(nc)%num_urbanc, filter(nc)%urbanc, &
                      filter(nc)%num_urbanp, filter(nc)%urbanp)

     ! Determine non snow-covered vegetation surface temperature and fluxes
     ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
     ! and leaf water change by evapotranspiration

     call CanopyFluxes(begg, endg, begc, endc, begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)

    ! Determine lake temperature and surface fluxes

     call SLakeFluxes(begc, endc, begp, endp, &
                         filter(nc)%num_lakec, filter(nc)%lakec, &
                         filter(nc)%num_lakep, filter(nc)%lakep)
     call SLakeTemperature(begc, endc, begp, endp, &
                              filter(nc)%num_lakec, filter(nc)%lakec, &
                              filter(nc)%num_lakep, filter(nc)%lakep)

     ! Determine soil/snow temperatures including ground temperature and
     ! update surface fluxes for new ground temperature.

     call Biogeophysics2(begl, endl, begc, endc, begp, endp, &
                         filter(nc)%num_urbanl,  filter(nc)%urbanl, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)

    ! Perform averaging from PFT level to column level
    ! YDT: potential end of pft loop? 

     call pft2col(begc, endc, filter(nc)%num_nolakec, filter(nc)%nolakec)

     ! Vertical (column) soil and surface hydrology
     call Hydrology2(begc, endc, begp, endp, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
                     filter(nc)%num_urbanc, filter(nc)%urbanc, &
                     filter(nc)%num_snowc, filter(nc)%snowc, &
                     filter(nc)%num_nosnowc, filter(nc)%nosnowc)

     ! Lake hydrology
     call SLakeHydrology(begc, endc, begp, endp, filter(nc)%num_lakec, filter(nc)%lakec, &
                            filter(nc)%num_lakep, filter(nc)%lakep)


     ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)

     do c = begc,endc
        l = clandunit(c)
        if (itypelun(l) == isturb) then
           ! Urban landunit use Bonan 1996 (LSM Technical Note)
           cps%frac_sno(c) = min( cps%snow_depth(c)/0.05_r8, 1._r8)
        end if
     end do


     ! Snow aging routine based on Flanner and Zender (2006), Linking snowpack
     ! microphysics and albedo evolution, JGR, and Brun (1989), Investigation of
     ! wet-snow metamorphism in respect of liquid-water content, Ann. Glaciol.
     ! ============================================================================
     ! Note the snow filters here do not include lakes; SnowAge_grain is called
     ! for lakes from SLakeHydrology.

     call SnowAge_grain(begc, endc, &
          filter(nc)%num_snowc, filter(nc)%snowc, &
          filter(nc)%num_nosnowc, filter(nc)%nosnowc)

     
     ! YDT now back to pft level? 
     ! Prescribed biogeography,
     ! prescribed canopy structure, some prognostic carbon fluxes
     call EcosystemDyn(begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep, &
                       doalb)

     ! Check the energy and water balance
     call BalanceCheck(begp, endp, begc, endc, begl, endl, begg, endg)

     ! Determine albedos for next time step
     if (doalb) then

        ! Albedos for non-urban columns
        call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, &
                           filter(nc)%num_nourbanc, filter(nc)%nourbanc, &
                           filter(nc)%num_nourbanp, filter(nc)%nourbanp, &
                           nextsw_cday, declinp1)

        ! Albedos for urban columns
        if (filter(nc)%num_urbanl > 0) then
           call UrbanAlbedo(nc, begl, endl, begc, endc, begp, endp,   &
                            filter(nc)%num_urbanl, filter(nc)%urbanl, &
                            filter(nc)%num_urbanc, filter(nc)%urbanc, &
                            filter(nc)%num_urbanp, filter(nc)%urbanp)
        end if

     end if

  ! Determine gridcell averaged properties to send to atm (l2as and l2af derived types)
  call clm_map2gcell( )

end program single_col
