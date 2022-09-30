      program graficador
      implicit real*8 (a-h,k-z)
      real*8 x(60000),y(60000)
      integer i,j,n
      open (1,file='simple10.dat',status='old')
      open (2,file='delta10.dat',status='replace')

      pi = 4.0d0*atan(1.0d0)
      i = 1
 3    read (1,*,err=4,end=4) xx,y1,zz,zz,zz,zz,zz
       x(i)=xx
       y(i)=y1
       i = i+1
      goto 3
 4    close (1)
      n = i-1
      do i=1,n
       exact = 6.0d0/(pi*x(i))**2+2.0d0/(pi*x(i)/dtan(pi*x(i))-1.0d0)
       write(2,100) x(i),dabs(exact-y(i))
      end do
      close(2)
 100  format (1p3e20.10)
      end
