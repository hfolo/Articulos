      implicit real*8 (a-h,o-z)
      real*8 r,rmin,rmax,dr,rho
      integer i,n
!
      open (3,file='density.in',status='replace')
!
      n = 100
!
      rmin = 0.0d0
      rmax = 1.0d0
      dr   = (rmax-rmin)/float(n)
!
      pi = 4.0d0*atan(1.0d0)
      rhoc = 1.0d0
      ck = pi/rmax
!
      do i = 1,n
       r = rmin+dr*float(i)
       rho = rhoc*sin(ck*r)/(ck*r)
!       rho = rhoc*(1.0d0-r*r)
!       rho = rhoc*(1.0d0-r)
       write (3,100) r,rho
      end do
!
 100  format (1p2d15.5)
!
      end
