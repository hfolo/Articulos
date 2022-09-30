      program graficador
      implicit real*8 (a-h,k-z)
      real*8 x(60000),y(60000)
      integer i,j,n
      open (1,file='datos-0.5.dat',status='old')
      open (2,file='densidad-0.5.in',status='replace')

       i = 1
 3    read (1,*,err=4,end=4) xx,y1,y2,y3
       x(i)=xx
       y(i)=y3
       i = i+1
      goto 3
 4    close (1)
      n = i-1
      r = x(i-1)
      do i=1,n
       write(2,100) x(i)/r,y(i)
      end do
      close(2)
 100  format (1p2e20.10)
      end
