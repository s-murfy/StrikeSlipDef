
module initiation
real,parameter :: pi = 3.14159265358979323846264338327950288419716939937510
real,parameter :: deg2rad=pi/180.
real,parameter :: rad2deg=180./pi



type obs_type
  integer :: Nobs
  double precision :: x0,y0,dip
  real :: lonref,latref     ! reference lat/lon used in coversion to/from UTM
  real  :: lon,lat
  double precision,allocatable,dimension(:):: x,y,z
  double precision,allocatable,dimension(:,:):: ux,uy,uz,p
  double precision,allocatable,dimension(:,:):: ux_ten,uy_ten,uz_ten,p_ten
end type

type source_type
  integer :: Nsources
  real :: lonref,latref    ! reference lat/lon used in coversion to/from UTM
  integer :: depth
  double precision,allocatable,dimension(:):: x,y,z
  real,allocatable,dimension(:):: latE,lonE,latW,lonW
  double precision,allocatable,dimension(:):: strike,rake,dip,dl,dw,slip


end type
private


public :: read_file,obs_type,source_type,output,output_piezo,output_geodetic

contains
!================================================================================
subroutine output_geodetic(obs,fname,fname_ten,N_sources)
  use georef

  implicit none
  type(obs_type) :: obs
  character(len=40) :: filename
  character(len=30) :: fname,fname_ten,file_id
 real :: x,y,lon,lat
  integer :: ion,n,i,j
  integer :: N_sources
!  type(geoid) :: gref
!  integer :: hemi,sens
  double precision :: geodetic_lon,geodetic_lat
  character(len=80)::fmt

! save strike slip deformation
  filename = trim(adjustl(fname)) // '.dat'
  open(14,file=trim(adjustl(filename)),form = 'formatted', action = 'write')
  geodetic_lon = obs%lon
  geodetic_lat = obs%lat
  ! convert from utm to lat/Long
  x = obs%x(1)+obs%x0
  y = obs%y(1)+obs%y0
  write(14,*) N_sources
  write(14,*) obs%Nobs
  write(14,*) obs%z(2)-obs%z(1)
  call convgeoutm(lat,lon,obs%latref,obs%lonref,x,y,-1)
  write(14,*) lon,lat
  write(fmt,'(a,i5,a)') '(',N_sources,'f20.6)'
  !write(*,'(a)') fmt
  write(14,fmt) (obs%ux(1,j),j = 1,N_sources)
  write(14,fmt) (obs%uy(1,j),j = 1,N_sources)
  write(14,fmt) (obs%uz(1,j),j = 1,N_sources)
  close(14)

! save tensile deformation
  filename = trim(adjustl(fname_ten)) // '.dat'
  open(15,file=trim(adjustl(filename)),form = 'formatted', action = 'write')
  write(15,*) N_sources
  write(15,*) 3  ! obs%Nobs
  write(15,*) obs%z(2)-obs%z(1)
  write(15,*) lon,lat
  write(fmt,'(a,i5,a)') '(',N_sources,'f20.6)'
  write(15,fmt) (obs%ux_ten(1,j),j = 1,N_sources)
  write(15,fmt) (obs%uy_ten(1,j),j = 1,N_sources)
  write(15,fmt) (obs%uz_ten(1,j),j = 1,N_sources)
  close(15)


return
end subroutine output_geodetic
!================================================================================
subroutine output(obs,fname,N_sources)
  use georef
  implicit none
  type(obs_type) :: obs
  character(len=40) :: filename
  character(len=30) :: fname,file_id
  real :: lon,lat,x,y
  integer :: ion,n,i,j
  integer :: N_sources

  do j = 1,N_sources
   write(file_id, '(i0)') j
   filename = trim(adjustl(fname)) // trim(adjustl(file_id)) // '.dat'
   open(14,file=trim(adjustl(filename)),form = 'formatted', action = 'write')
  ! convert from utm to lat/Long
   do i = 1,obs%Nobs
     x = real(obs%x(i)+obs%x0)
     y = real(obs%y(i)+obs%y0)

     call convgeoutm(lat,lon,obs%latref,obs%lonref,x,y,-1)
     ! lat = lat*rad2deg
     ! lon = lon*rad2deg

!save to file lat/Long/depth/pressure/surface
    write(14,*) lon,lat,obs%z(i),obs%p(i,j),obs%ux(i,j),obs%uy(i,j),obs%uz(i,j)
  enddo
 close(14)
enddo

return
end subroutine output
!================================================================================
subroutine read_file(source,surf,obs_tekr)
use georef

implicit none
double precision,allocatable,dimension(:,:) :: x,y,d
double precision :: x1dash,x2dash,x3dash,x4dash
double precision :: y1dash,y2dash,y3dash,y4dash
double precision :: d1dash,d2dash,d3dash,d4dash
double precision :: ang
double precision :: depth,xx,yy
real :: qlon,qlat,x0,y0
double precision :: dp3,xp3,dp4,xp4
double precision :: x1,x2,x3,x4
double precision :: y1,y2,y3,y4
double precision :: d1,d2,d3,d4
real :: lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4
double precision :: x_mean,y_mean,d_mean
double precision, parameter :: kilo = 1000.d0
double precision :: x1dd,x2dd,x3dd,x4dd, y1dd,y2dd,y3dd,y4dd
real :: obs_lat,obs_lon
double precision :: obs_len,obs_depth,obs_dx,obs_dz
integer :: surf_nx,surf_ny
double precision :: surf_dx,surf_dy,dz
integer :: ion,i,n,rc
integer :: icount,j,nl,nz
double precision :: obs_dl,obs_x0,obs_y0,obs_z0
double precision :: source_dy,source_dz
type(obs_type) :: surf,obs_tekr
real :: latref,lonref
real :: tekr_lon, tekr_lat
double precision :: ref_x0,ref_y0
double precision :: rake, slipping_width
integer :: N_sources
double precision :: slat2,slon2,slat1,slon1,stk,sdw,sdl,sdip,sstk
type(source_type) :: source
character(80) :: fault_file
logical :: file_exists
namelist /PARAMETERS/ depth, rake, slipping_width, fault_file
namelist /OBSERVATIONS/  tekr_lon, tekr_lat


! reference position is based on first coordinate trace
! latref = (40.72007535)
! lonref = (29.27885615)
! latref =40.6681
! lonref =28.1123
latref = 0.
lonref = 27.0

!================= set reference coordinate ==========
qlat = 40.8739
qlon = 27.911 !d0!27.750807d0!-360.d0

! convert from lat/lon to utm
call convgeoutm(qlat,qlon,latref,lonref,x0,y0,1)

ref_x0 = dble(x0)
ref_y0 = dble(y0)


!======================= set sources =============================
! Open and read Namelist file.
open (50,file='input_file.txt', action='read',  iostat=rc)
read (50,nml=PARAMETERS, iostat=rc)
if (rc /= 0) then
    write (*,*) 'Error: invalid PARAMETERS format'
    stop
endif
read (50,nml=OBSERVATIONS, iostat=rc)
if (rc /= 0) then
    write (*,*) 'Error: invalid OBSERVATIONS format'
    stop
endif
close (50)

inquire (file=fault_file, exist=file_exists)   ! check file exists
if (.not.file_exists) then
      write (*,*) 'Error: file "', trim(fault_file), '" does not exist.'
      stop
end if

!================== read in fault coordinates ===============
!open(12,file='fault_trace_coordinates.txt',action='read')
open(12,file=fault_file,action='read')

read(12,*) N_sources
source%Nsources = N_sources
allocate(source%latE(N_sources),source%lonE(N_sources),source%latW(N_sources),source%lonW(N_sources))
allocate(source%x(N_sources),source%y(N_sources),source%z(N_sources))
allocate(source%strike(N_sources),source%rake(N_sources),source%dip(N_sources),source%slip(N_sources))
allocate(source%dl(N_sources),source%dw(N_sources))


allocate(x(n,4),y(n,4),d(n,4))

source%latref = latref
source%lonref = lonref

open(55,file='../res/fault_latlon.txt',action='write')

do i = 1,N_sources
      read(12,*) slat2,slon2,slat1,slon1,stk,sdl,sdw
      source%latE(i) = slat1
      source%lonE(i) = slon1
      source%latW(i) = slat2
      source%lonW(i) = slon2
      source%strike(i) = stk
      source%rake(i) = rake !180.d0 !-140.d0!-147.d0 !180.d0
      source%dip(i) =  90.d0! 82.d0!69.d0!90.d0
      source%dl(i) = sdl
      source%dw(i) = -slipping_width! -2000.d0!10000.d0-sdw !-sdw ! don't change this!
      source%slip(i) = 1.d0

      call convgeoutm(source%latW(i),source%lonW(i),latref,lonref,x0,y0,1)

! define fault on x-y plane and then rotate along x-axis to place at correct dip
!  note: based on OKADA, fault needs to dip in negative y-axis direction
      source%depth = int(depth)
      source_dy = -depth*cos(source%dip(i)*deg2rad)
      source_dz = -depth*sin(source%dip(i)*deg2rad)
      source%x(i) = dble(x0)-ref_x0
      source%y(i) = dble(y0)-ref_y0+source_dy
      source%z(i) = abs(source_dz)!
      !          ^   Northing (y))
      !          |
      !          |
      !          .------>  Easting (x)
      !
      !                 length
      !     x1,y1   .---------.  x2,y2
      !             |         |
      !             |         |   width
      !             |         |
      !             |         |
      !     x4,y4   .---------.  x3,y3
      !
      x1 = 0.d0; y1 = 0.d0;       d1 = 0.d0
      x2 = 0.d0 ; y2 = source%dl(i);   d2 = 0.d0
      x3 = source%dw(i) ; y3 = source%dl(i);  d3 = 0.d0
      x4 = source%dw(i) ; y4 = 0.d0;     d4 = 0.d0


      ! rotate for dip
      sdip = source%dip(i)*deg2rad
      xp3 = x3*cos(sdip)+d3*sin(sdip)
      dp3 = -x3*sin(sdip)+d3*cos(sdip)

      xp4 = x4*cos(sdip)+d4*sin(sdip)
      dp4 = -x4*sin(sdip)+d4*cos(sdip)

      x1dash = x1; y1dash = y1; d1dash = d1
      x2dash = x2; y2dash = y2; d2dash = d2
      x3dash = xp3; y3dash = y3; d3dash = dp3
      x4dash = xp4; y4dash = y4; d4dash = dp4

      ! rotate for strike
      sstk = source%strike(i)*deg2rad;
      x1dd=y1dash*sin(sstk)+x1dash*cos(sstk)
      y1dd=y1dash*cos(sstk)-x1dash*sin(sstk)

      x2dd=y2dash*sin(sstk)+x2dash*cos(sstk)
      y2dd=y2dash*cos(sstk)-x2dash*sin(sstk)

      x3dd=y3dash*sin(sstk)+x3dash*cos(sstk)
      y3dd=y3dash*cos(sstk)-x3dash*sin(sstk)

      x4dd=y4dash*sin(sstk)+x4dash*cos(sstk)
      y4dd=y4dash*cos(sstk)-x4dash*sin(sstk)


      ! correction to place slip relative to reference coordinates
      x1dash = x1dd+source%x(i)+ref_x0; y1dash = y1dd+source%y(i)+ref_y0;
      x2dash = x2dd+source%x(i)+ref_x0; y2dash = y2dd+source%y(i)+ref_y0;
      x3dash = x3dd+source%x(i)+ref_x0; y3dash = y3dd+source%y(i)+ref_y0;
      x4dash = x4dd+source%x(i)+ref_x0; y4dash = y4dd+source%y(i)+ref_y0;


! convert back to lat/lon
      call convgeoutm(lat1,lon1,latref,lonref,real(x1dash),real(y1dash),-1)
      call convgeoutm(lat2,lon2,latref,lonref,real(x2dash),real(y2dash),-1)
      call convgeoutm(lat3,lon3,latref,lonref,real(x3dash),real(y3dash),-1)
      call convgeoutm(lat4,lon4,latref,lonref,real(x4dash),real(y4dash),-1)

      write(55,*) lon1,lat1,d1dash,         &
                  lon2,lat2,d2dash,         &
                  lon3,lat3,d3dash,         &
                  lon4,lat4,d4dash,source%strike(i) !source%slip(i)

enddo

close(55)

close(12)

! ============ coordinates of observation  =============
!===== produce surface deformation around  fault (used for debugging)
  obs_lat = 40.8739d0!   give a poistion in the area of interest
  obs_lon = 27.911d0 !
! convert from Lat/Long to UTM
  call convgeoutm(obs_lat,obs_lon,latref,lonref,x0,y0,1)
  obs_x0 = dble(x0)
  obs_y0 = dble(y0)
  obs_z0 =0.d0
  surf_dx = 500.d0
  surf_dy = 500.d0
  surf_nx = 200
  surf_ny = 200

!  =============  map of surface
  obs_z0 = 0.d0
  call set_surface(ref_x0,ref_y0,obs_x0,obs_y0,obs_z0,surf_nx,surf_ny,surf_dx,surf_dy,latref,lonref,surf,'xy')

  ! =============  pt observation @ tek
  ! obs_lon = 27.496
  ! obs_lat = 40.958
  dz = 100.d0
  nz = 150    ! go to a depth of 15km
  call convgeoutm(tekr_lat,tekr_lon,latref,lonref,x0,y0,1)
  obs_x0 = dble(x0)
  obs_y0 = dble(y0)
  call set_pt_obs(ref_x0,ref_y0,obs_x0,obs_y0,dz,nz,latref,lonref,obs_tekr)

return
end subroutine read_file
!========================================================================================
subroutine set_pt_obs(x0,y0,obs_x0,obs_y0,dz,nz,latref,lonref,obs)
  implicit none
  type(obs_type) :: obs
  real :: latref,lonref
  double precision :: obs_x0,obs_y0,x0,y0,dz,xx,yy,ang
  integer :: i,nz
  ang = -90.d0*deg2rad   ! check if correct!!

  obs%Nobs = nz
  obs%x0 = x0
  obs%y0 = y0
  obs%dip = 0.d0

  obs%latref = latref
  obs%lonref = lonref

  allocate(obs%x(obs%Nobs),obs%y(obs%Nobs),obs%z(obs%Nobs))
  do i = 1,obs%Nobs
        xx = obs_x0-x0
        yy = obs_y0-y0
        obs%z(i) = -dble(i-1)*dz
        obs%x(i) = xx
        obs%y(i) = yy
  enddo

return
end subroutine set_pt_obs
!========================================================================================
subroutine set_surface(x0,y0,obs_x0,obs_y0,obs_z0,nx,ny,dx,dy,latref,lonref,obs,plane)
implicit none
type(obs_type) :: obs
real :: latref,lonref

integer :: i,j,icount
double precision :: x0,y0,dx,dy,ang,xx,yy,obs_x0,obs_y0,obs_z0
double precision :: xx0,yy0
integer :: nx,ny
character(2) :: plane

!ang = 0.d0*deg2rad
ang = -90.d0*deg2rad   ! check if correct!!


obs%Nobs = nx*ny
obs%latref = latref
obs%lonref = lonref

allocate(obs%x(obs%Nobs),obs%y(obs%Nobs),obs%z(obs%Nobs))
icount = 0
obs%x0 = x0
obs%y0 = y0
obs%dip = 0.d0

if (plane == 'xy') then
  do i =1,nx
    do j = 1,ny
     icount = icount+1
     xx0 = obs_x0-x0
     yy0 = obs_y0-y0
     xx = xx0-float(nx)/2.d0*dx+float(i-1)*dx
     yy = yy0-float(ny)/2.d0*dy+float(j-1)*dy
     obs%z(icount) = 0.d0
     obs%x(icount) = xx*cos(ang)-yy*sin(ang)
     obs%y(icount) = xx*sin(ang)+yy*cos(ang)
    enddo
  enddo
elseif(plane == 'yz') then
  do i =1,nx
    do j = 1,ny
      icount = icount+1
      xx0 = obs_x0-x0
      yy0 = obs_y0-y0
      xx = xx0
      yy = yy0-float(nx)/2.d0*dx+float(i-1)*dx  ! make a N-S cross-section
      obs%z(icount) = float(j-1)*dy
      obs%x(icount) = xx !xx*cos(ang)-yy*sin(ang)
      obs%y(icount) = yy !xx*sin(ang)+yy*cos(ang)
    enddo
  enddo
else
  write(*,*) 'Error in definition of plane'
  stop
endif

return
end subroutine set_surface
!========================================================================================
end module initiation
