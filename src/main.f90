program  StrikeSlipDef
use initiation
use stiffness, only: dtau
implicit none
character(len=30):: filename,filename_ten
character(len=30) :: format_string,format_string2


type(obs_type) :: surf,obs_tekr
type(source_type) :: source
double precision,allocatable,dimension(:) :: rake,dip,strike,dl,dw,slip
double precision,allocatable,dimension(:,:) :: x,y,d

! read in file
call read_file(source,surf,obs_tekr)

! calculation
write(*,*) 'calc surface displacement'
call dtau(source,obs_tekr)

! save data
!filename = 'surf0_'
!call output(surf,filename,source%Nsources)

if (source%depth < 10000) then
  format_string = "(A16,I4.4)"
  format_string2 = "(A20,I4.4)"
else
  format_string = "(A16,I5.5)"
  format_string2 = "(A20,I5.5)"
endif

!filename = 'obs_tekr_14000'
write (filename,format_string) "../res/obs_tekr_", source%depth
write (filename_ten,format_string2) "../res/obs_tekr_ten_", source%depth

call output_geodetic(obs_tekr,filename,filename_ten,source%Nsources)


end program StrikeSlipDef
