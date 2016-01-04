
! Usage: 
! read_netcdf inputfile data_type varname dim1_name dim2_name ... 
! Three data_type: int, float, double 
! Max four dimensions 


  program read_netcdf 

  use netcdf

  implicit none
  include "netcdf.inc"

  character (len = 200) :: infile 
  character (len = 100) :: data_type, varname
  character (200) ::  dim_name(18) 
  !integer :: start(18), ocount(18), dimid(18), count(18)   
  integer :: start(4), ocount(4), dimid(4), count(4)   
  real*8, allocatable  :: ddata1d(:)
  real*8, allocatable  :: ddata4d(:, :, :, :)

  integer :: ncid, varid, iargc, ndims, nargs, i, j, k, l, len, ic 

  nargs = iargc() 
  if (nargs .gt. 7 ) then 
   write(*, *)"max dimension is 4! " 
   stop
  end if 

  ndims = nargs - 3 

  call getarg(1, infile) 
  call getarg(2, data_type) 
  call getarg(3, varname) 
  do i=1, ndims 
    call getarg(i+3, dim_name(i) ) 
  end do 
  
  call check( nf90_open(infile, NF90_NOWRITE, ncid) )

  ! get all the dimension sizes 
  start(:) = 1
  count(:) = 1
  len = 1 
  !do i=1, ndims 
  do i=ndims, 1, -1  ! swap the dim order for fortran for count() array
    ic = ndims - i + 1  ! reverse order
    call check( nf90_inq_dimid(ncid, trim(dim_name(i)), dimid(i) ) ) 
    call check( nf90_inquire_dimension(ncid, dimid(i), len=count(ic) ) ) 
    write(*, *) trim(dim_name(i)), " dimid= ", dimid(i), "  dim size=", count(ic) 
    len = len * count(ic) 
  end do 
  do i=1, ndims
    write(*, *) count(i) 
  end do 

  call check( nf90_inq_varid(ncid, trim(varname), varid) )

! forget data type matching. Convert everything to double when reading whatever types
! see https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/Type-Conversion.html#Type-Conversion
!  select case (trim(data_type)) 
!    case ('int') 
           allocate(ddata4d(count(1), count(2), count(3), count(4)))  
           !call check( nf90_get_var(ncid, varid, ddata1d, start=start, count=count ) )
           call check( nf90_get_var(ncid, varid, ddata4d, start=start, count=count ) )
           !call check( nf90_get_var(ncid, varid, ddata4d ) )
           !call check( nf90_get_var(ncid, varid, ddata1d ) )
           write(*, '(10F10.2)') real(ddata4d)
           !write(*, '(10F10.2)') real(ddata1d)

  ! other types to be implemented 
  !end select 

  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )
  deallocate(ddata1d)  


contains

  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

end program read_netcdf 

