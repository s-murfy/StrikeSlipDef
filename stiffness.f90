module stiffness
use initiation, only: obs_type,source_type
implicit none
real,parameter :: pi = 3.14159265358979323846264338327950288419716939937510
real,parameter :: deg2rad=pi/180.
real,parameter :: rad2deg=180./pi


private

public::  dtau

contains

!!===============================================
subroutine dtau(source,obs)
use okada
implicit none
type(obs_type) :: obs
type(source_type) :: source
double precision,allocatable,dimension(:) :: slip,rake,dip,strike,dl,dw
double precision,allocatable,dimension(:,:) :: x,y,d,x_l,y_l
double precision,dimension(3,3) :: lsig,f_sig
double precision :: AW1,AW2,AL1,AL2,alpha,z,c
integer :: iret
double precision :: UX1,UX2,UX3,UY1,UY2,UY3,UZ1,UZ2,UZ3
double precision :: UXX(3),UYX(3),UZX(3),UXY(3),UYY(3),UZY(3),UXZ(3),UYZ(3),UZZ(3)
double precision :: lam_mu
double precision :: dil_ss,sxy_ss,sxz_ss,syz_ss,sxx_ss,syy_ss,szz_ss
double precision :: dil_ds,sxy_ds,sxz_ds,syz_ds,sxx_ds,syy_ds,szz_ds
double precision :: dil_ten,sxy_ten,sxz_ten,syz_ten,sxx_ten,syy_ten,szz_ten
double precision :: dil,sxy,syz,sxz,sxx,syy,szz
double precision :: obs_x_l,obs_y_l
double precision :: x_obs,y_obs,z_obs,nu,mu,Fdepth
double precision :: slip_SS,slip_DS
double precision :: t_13,t_23,s_n,ang,ux_ss,uy_ss,ux_ds,uy_ds
double precision :: ux_tensile,uy_tensile,uz_tensile
integer :: i,j,n,no_slip_cells,no_aftershocks
character(3) :: mode
!real(dp),intent(out) :: k_tau,k_sig
!  Okada coordinate convention:
!   x is dirn along strke
!   y is dirn perpendicular to fault
!   z is height (i.e. neg in direction of depth)
!   c is depth (i.e. pos in dirn of depth)
!! Convention:
!       cell(i) = reciever cell of stress change
!       cell(j) = source cell of stress change
!
!  Okada Coordinate System
!
!      Z  y
!      | /
!      |/
!      . - -> x
!
!      _____________
!     |             |
!     |   fault     |
!     |             |
!      -------------

!      _______
!     |       |
!     |   o   |         x = ref pt for slipping cell
!     |       |         o = observation pt for okada
!     x ------
!
!
lam_mu = 0.d0
dil = 0.d0
sxy = 0.d0
sxz = 0.d0
syz = 0.d0
sxx = 0.d0
syy = 0.d0
szz = 0.d0

nu = 0.25d0

lam_mu = nu/(1.d0-2.d0*nu)
alpha = -0.5d0/(nu-1.d0)
allocate(x_l(no_slip_cells,4),y_l(no_slip_cells,4))

! deformation due to strike slip
allocate(obs%ux(obs%Nobs,source%Nsources),obs%uy(obs%Nobs,source%Nsources),obs%uz(obs%Nobs,source%Nsources))
allocate(obs%p(obs%Nobs,source%Nsources))


! deformation due to tensile motion
allocate(obs%ux_ten(obs%Nobs,source%Nsources),obs%uy_ten(obs%Nobs,source%Nsources),obs%uz_ten(obs%Nobs,source%Nsources))
allocate(obs%p_ten(obs%Nobs,source%Nsources))


obs%p = 0.d0; obs%ux = 0.d0;obs%uy = 0.d0;obs%uz = 0.d0
obs%p_ten = 0.d0; obs%ux_ten = 0.d0;obs%uy_ten = 0.d0;obs%uz_ten = 0.d0


do j = 1,source%Nsources
       slip_SS = source%slip(j)*cos(source%rake(j)*deg2rad) ! positive slip = left lateral
       slip_DS = source%slip(j)*sin(source%rake(j)*deg2rad) ! positive slip = reverse-fault
       if (abs(slip_SS) < epsilon(slip_SS)) slip_SS = 0.d0
       if (abs(slip_DS) < epsilon(slip_DS)) slip_DS = 0.d0
! positive tensile slip = open-type movement
       c = abs(source%z(j))
       ang = (source%strike(j)-90.d0)*deg2rad
       AL1 = 0.d0
       AL2 = source%dl(j)
       AW1 = source%dw(j)
       AW2 = 0.d0   ! reference point is top left hand corner of cell
!========================== observations ====================================
          do i = 1,obs%Nobs
 ! rotate observation point to correct position
            obs_x_l = obs%x(i) - source%x(j)
            obs_y_l = obs%y(i) - source%y(j)
            x_obs = obs_x_l*cos(ang)-obs_y_l*sin(ang)
            y_obs = obs_x_l*sin(ang)+obs_y_l*cos(ang)
! observation points (relative to the bottom left hand corner of the slipping patch)
             z_obs = obs%z(i)
              call  DCD3(alpha,x_obs,y_obs,z_obs,c,                            &
                     source%dip(j),AL1,AL2,AW1,AW2,                            &
                     UX1,UX2,UX3,UY1,UY2,UY3,UZ1,UZ2,UZ3,                      &
                     UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
              if (iret == 1) then
                     print*, 'error produced by DCD3'
                     stop
              endif
!    UXX(1) = strike slip, UXX(2) = dip slip , UXX(3) = tensile slip
! STRIKE SLIP (1 = STRIKE SLIP)
              dil_ss = lam_mu*(UXX(1)+UYY(1)+UZZ(1))
              sxy_ss = mu*(UYX(1)+UXY(1))
              sxz_ss = mu*(UXZ(1)+UZX(1))
              syz_ss = mu*(UYZ(1)+UZY(1))
              sxx_ss = 2.d0*mu*(dil+UXX(1))
              syy_ss = 2.d0*mu*(dil+UYY(1))
              szz_ss = 2.d0*mu*(dil+UZZ(1))

! DIP SLIP (2 = DIP SLIP)
              dil_ds = lam_mu*(UXX(2)+UYY(2)+UZZ(2))
              sxy_ds = mu*(UYX(2)+UXY(2))
              sxz_ds = mu*(UXZ(2)+UZX(2))
              syz_ds = mu*(UYZ(2)+UZY(2))
              sxx_ds = 2.d0*mu*(dil+UXX(2))
              syy_ds = 2.d0*mu*(dil+UYY(2))
              szz_ds = 2.d0*mu*(dil+UZZ(2))

! Tensile SLIP (3 = Tensile)
             dil_ten = lam_mu*(UXX(3)+UYY(3)+UZZ(3))
             sxy_ten = mu*(UYX(3)+UXY(3))
             sxz_ten = mu*(UXZ(3)+UZX(3))
             syz_ten = mu*(UYZ(3)+UZY(3))
             sxx_ten = 2.d0*mu*(dil+UXX(3))
             syy_ten = 2.d0*mu*(dil+UYY(3))
             szz_ten = 2.d0*mu*(dil+UZZ(3))

! combine components
              dil = dil_ss*slip_SS + dil_ds*slip_DS
              sxy = sxy_ss*slip_SS + sxy_ds*slip_DS
              syz = syz_ss*slip_SS + syz_ds*slip_DS
              sxz = sxz_ss*slip_SS + sxz_ds*slip_DS
              sxx = sxx_ss*slip_SS + sxx_ds*slip_DS
              syy = syy_ss*slip_SS + syy_ds*slip_DS
              szz = szz_ss*slip_SS + szz_ds*slip_DS

              lsig(1,:) = (/ sxx,sxy,sxz /)
              lsig(2,:) = (/ sxy,syy,syz /)
              lsig(3,:) = (/ sxz,syz,szz /)
!              f_sig = lsig
              f_sig = rotate(lsig,(source%strike(j)-90.d0)*deg2rad,obs%dip,-1)      ! str = zero okay

             obs%p(i,j) = (sxx+syy+szz)/3.d0
             obs%p_ten(i,j) = (sxx_ten+syy_ten+szz_ten)/3.d0

              ! calculation of displacement
              ang = (source%strike(j)-90.d0)*deg2rad
              ux_ss = UX1*dcos(ang)+UY1*dsin(ang)
              uy_ss = -UX1*dsin(ang)+UY1*dcos(ang)
              ux_ds = UX2*cos(ang)+UY2*sin(ang)
              uy_ds = -UX2*sin(ang)+UY2*cos(ang)
              ux_tensile = UX3*cos(ang)+UY3*sin(ang)
              uy_tensile = -UX3*sin(ang)+UY3*cos(ang)
              uz_tensile = UZ3
!
              obs%ux(i,j) = ux_ss*slip_SS + ux_ds*slip_DS
              obs%uy(i,j) = uy_ss*slip_SS + uy_ds*slip_DS
              obs%uz(i,j) = UZ1*slip_SS + UZ2*slip_DS

             obs%ux_ten(i,j) = ux_tensile
             obs%uy_ten(i,j) = uy_tensile
             obs%uz_ten(i,j) = UZ3

  enddo   ! loop over observation points
enddo   ! loop over active cells

return
end subroutine dtau
!!===============================================
!----------------------------------------------------------
pure function  rotate(mat,str,dip,i) result(sigma)
   implicit none
   double precision,dimension(3,3),intent(in) :: mat
   double precision,dimension(3,3) :: rot,trot,sigma,mat1
   double precision,intent(in) :: str,dip
   integer,intent(in) :: i

  rot(1,:) = (/ dsin(str)             , dcos(str)                ,0.d0       /)
  rot(2,:) = (/ -dcos(str)*dcos(dip)   , dsin(str)*dcos(dip)     ,dsin(dip)  /)
  rot(3,:) = (/ dcos(str)*dsin(dip)    , -dsin(str)*dsin(dip)    ,dcos(dip)  /)
   trot = transpose(rot)
   if (i >= 0) then
   ! this is a forward rotation (i.e. ) from global to local fault coordinates
    mat1 = matmult3(mat,trot)
    sigma = matmult3(rot,mat1)
   else
    ! this is a reverse rotation from local to global coordinates
    mat1 = matmult3(mat,rot)
    sigma = matmult3(trot,mat1)
   endif
end function rotate
!----------------------------------------------------------
! Function that multiplies two 3 x 3 matrices together
pure function matmult3(a,b)     result(c)
   implicit none
   double precision,dimension(3,3),intent(in) :: a,b
   double precision,dimension(3,3) :: c
   integer :: i,j,k

   c = 0.d0

   do i = 1,3
    do j = 1,3
     do k = 1,3
     c(i,j) = c(i,j) + a(i,k)*b(k,j)
     enddo
    enddo
   enddo

end function matmult3
!----------------------------------------------------------
end module stiffness
