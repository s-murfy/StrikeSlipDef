! Georef is a fortran90 module to convert coordinates between UTM and geographic
! This code is the application of the formulas in
! Map Projections: A Working Manual, p. 63 by J. P. Snyder (1987) published by the USGS.
! This has been take from the open source program georef at: https://github.com/andherit/georef
! procduced by André Herrero - andre.herrero@ingv.it
!
! No modifications have been made to the original code
!
!  Reference : doi : 10.5281/zenodo.840875
!
!
module georef
! Set of conversion tools between geographic coordinates
! and UTM ones :
! - function meriref : compute the reference meridian for a given UTM zone
! - function loadgeoid : load the parameter for a given geoid, e.g. wgs84
! - subroutine convgeoutm : principal frontend.
!   compute the relative position in meter of a geographic point. The function
!   is able to work in both way.
! - function guessfu : compute the UTM zone corresponding to a longitude
! - subroutine utmgeo : conversion from utm to geographic coordinates.
! - subroutine geoutm : conversion from geographic to utm coordinates.


   implicit none
   real*8, private, parameter :: pi=acos(-1.d0)
   real*8, private, parameter :: ak=0.9996d0


! parameter for wgs84 ellipsoide
!  real*8, private, parameter :: cc=6399593.6257585
!  real*8, private, parameter :: ep=0.08209443794984
   real*8, private, parameter :: cc=6378.137
   real*8, private, parameter :: ep=0.08181919132
   real*8, private, parameter :: e2=ep*ep

   type geoid
      real*8 :: ra,rb
      real*8 :: ef
   end type geoid

contains

!############################################################
function meriref(ifu)

      integer :: ifu
      real*8 :: meriref

      meriref=(dble(float(ifu))-.5)*(1.d0/30.d0*pi)-pi
   return
end function meriref
!############################################################
function loadgeoid(label)

      character*(*) :: label
      type(geoid) :: loadgeoid

      select case (label(1:len_trim(label)))
         case('wgs84')
           loadgeoid%ra=6378137.
           loadgeoid%rb=6356752.3142
         case('ed50')
           loadgeoid%ra=6378388.
           loadgeoid%rb=6356911.9
         case('clark')
            loadgeoid%ra=6378206.4
            loadgeoid%rb=6356583.8
         case default
           write(0,*) 'unlist geoid... using wgs84'
           loadgeoid%ra=6378137.
           loadgeoid%rb=6356752.3142
      end select
      loadgeoid%ef=sqrt(1.d0-loadgeoid%rb*loadgeoid%rb/loadgeoid%ra/loadgeoid%ra)
      return
end function loadgeoid

!############################################################
subroutine convgeoutm(lat,lon,latref,lonref,x,y,sens)

! front end
! the utm values are computed relatively to a reference point
! warning : the reference meridian is lonref and NOT the
! standard reference medirian of the zone where belong the
! reference point.
! sens is a switch to the direction of the conversion
!     sens=1  : absolute geo -> relative utm
!     sens=-1 : relative utm -> absolute geo
! units : angle in degrees
!         position in meters

      real*4 :: lat,lon,latref,lonref,x,y
      integer :: sens

      integer :: kh
      real*8 :: x0,y0,dlar,dlor,dla,dlo,xp,yp,lo
      type(geoid) :: gref

! position of the reference point in utm
      dlar=dble(latref)/180.*pi
      dlor=dble(lonref)/180.*pi
      gref=loadgeoid('wgs84')
!      gref=loadgeoid('clark')
      call geoutm(dlor,dlar,x0,y0,dlor,gref)
      if (sens > 0) then
         dla=dble(lat)/180.*pi
         dlo=dble(lon)/180.*pi
         call geoutm(dlo,dla,xp,yp,dlor,gref)
         x=sngl(xp-x0)
         y=sngl(yp-y0)
      else
         xp=dble(x)+x0
         yp=dble(y)+y0
         kh=sign(1,int(latref))
         call utmgeo(dlo,dla,xp,yp,dlor,kh,gref)
         lat=sngl(dla/pi)*180.
         lon=sngl(dlo/pi)*180.
      endif
    return
end subroutine convgeoutm

!############################################################
function guessfu(lon)

!  WARNING !!! lon is expressed in degree
! This function returns the reference meridian for a given UTM zone.
! Exceptions are not taking into account.

      integer :: guessfu
      real*4 :: lon

      if (lon > 180.) then
         guessfu=int(lon-360.)/6+30
      else
         guessfu=int(lon)/6+31
      endif
   return
end function guessfu

!############################################################
subroutine utmgeo(long,lat,xx,yy,lonref,hemi,gref)

! convert coordinates from utm (xx,yy) in meters (x -> East, y -> North) to geographic
! (long,lat) in RAD., using the geoid gref and the reference
! mederidian lonref (rad.). hemi is equal to 1 in northern hemisphere and
! -1 in the southern one

! Reference for formula :
! Map Projections: A Working Manual by J. P. Snyder (1987) published by
! the usgs (available online). page 63.

     real*8, intent(in) :: xx,yy,lonref
     integer, intent(in) :: hemi
     real*8, intent(out) :: long,lat
     type(geoid) :: gref


     real*8 :: ef2,ef4,ef6,ym
     real*8 :: m,mu,e1,p1,n1,r1,d,t1,c1,ep2

     ef2=gref%ef**2.d0
     ef4=gref%ef**4.d0
     ef6=gref%ef**6.d0
     ep2=ef2/(1.d0-ef2)
! false easting removal
     if (hemi < 0) then
        m=(yy-10000000.d0)/ak
     else
        m=yy/ak
     endif
     mu=m/(gref%ra*(1.d0-ef2/4.d0-(3.d0/64.d0)*ef4-(5.d0/256.d0)*ef6))
     e1=(1.d0-sqrt(1.d0-ef2))/(1.d0+sqrt(1.d0-ef2))
     p1=mu+(1.5d0*e1-(27.d0/32.d0)*e1**3.d0)*sin(2.d0*mu)
     p1=p1+((21.d0/16.d0)*e1*e1-(55.d0/32.d0)*e1**4.d0)*sin(4.d0*mu)
     p1=p1+(151.d0/96.d0)*e1**3.d0*sin(6.d0*mu)
     p1=p1+(1097.d0/512.d0)*e1**4.d0*sin(8.d0*mu)
     n1=gref%ra/sqrt(1.d0-ef2*sin(p1)*sin(p1))
     r1=gref%ra*(1.d0-ef2)/(1.d0-ef2*sin(p1)*sin(p1))**1.5d0
! false easting removal
     d=(xx-500000.d0)/n1/ak
     t1=tan(p1)**2.d0
     c1=ep2*cos(p1)**2.d0
     lat=d*d/2.d0
     lat=lat-(5.d0+3.d0*t1+10.d0*c1-4.d0*c1*c1-9.d0*ep2)*d**4.d0/24.d0
     lat=lat+(61.d0+90.d0*t1+298.d0*c1+45.d0*t1**2.d0-252*ep2-3.d0*c1**2.d0)*d**6.d0/720.d0
     lat=p1-lat*(n1*tan(p1)/r1)
     long=d-(1.d0+2.d0*t1+c1)*d**3.d0/6.d0
     long=long+(5.d0-2.d0*c1+28.d0*t1-3.d0*c1*c1+8.d0*ep2+24.d0*t1*t1)*d**5.d0/120.d0
     long=lonref+long/cos(p1)
     return
end subroutine utmgeo
!############################################################
subroutine geoutm(long,lat,xx,yy,lonref,gref)

! convert coordinates from geographic (long,lat) in rad. to
! utm (xx,yy) in meters (x -> East, y -> North), using the geoid gref and the reference
! mederidian lonref (rad.)

! Reference for formula :
! Map Projections: A Working Manual by J. P. Snyder (1987) published by
! the usgs (available online). page 61.

     real*8 :: long,lat,xx,yy,lonref
     type(geoid) :: gref

     real*8 :: ef2,ef4,ef6
     real*8 :: nu,sl,cl,a
     real*8 :: sp,s1,s2,s3,s4
     real*8 :: t,c

     ef2=gref%ef**2.d0
     ef4=gref%ef**4.d0
     ef6=gref%ef**6.d0
     sl=sin(lat)
     cl=cos(lat)
     nu=1./sqrt(1-ef2*sl*sl)
     a=(long-lonref)*cl
     s1=1.d0-ef2/4.d0-3.d0*ef4/64.d0-5.d0*ef6/256.d0
     s2=3.d0*ef2/8.d0+3.d0*ef4/32.d0+45.d0*ef6/1024.d0
     s3=15.d0*ef4/256.d0+45.d0*ef6/1024.d0
     s4=35.d0*ef6/3072.d0
     sp=s1*lat-s2*sin(2.d0*lat)+s3*sin(4.d0*lat)-s4*sin(6.d0*lat)
     t=tan(lat)*tan(lat)
     c=ef2/(1.d0-ef2)*cl*cl
     xx=500000.d0+ak*gref%ra*nu*(a+(1.d0-t+c)*a*a*a/6.d0+(5.d0-18.d0*t+t*t)*(a**5.d0)/120.d0)
     yy=ak*gref%ra*(sp+nu*tan(lat)*(a*a/2.d0+(5.d0-t+9.d0*c+4.d0*c*c)*a*a*a*a/24.d0+(61.d0-58.d0*t+t*t)*(a**6.d0)/720.d0))
! false northing for south hemisphere
     if (lat < 0.d0) yy=yy+10000000.d0
    return
end subroutine geoutm

end module georef
