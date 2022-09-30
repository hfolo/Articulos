      program graficador
      implicit real*8 (a-h,k-z)
      integer i,j
      open (2,file='3ene.dat',status='old')
      open (3,file='gf.dat',status='replace')
      open (4,file='app.dat',status='replace')


 13    read (2,*,err=14,end=14) ene,h,kf,c
       gf = 5.0d0*c/2.0d0
       exacto = (5.0d0-2.0d0*dsqrt(2.0d0/h-1.0d0))/3.0d0
       aproxi = (1.0d0+2.0d0*h)/3.0d0
       write(3,102)  h,gf
      goto 13
 14   close (2)

      hmin=0.4d0
      hmax=1.0d0
      dh  = (hmax-hmin)/1000.0d0
      do i=0,1000
       h=hmin+dh*dfloat(i)
       a = 0.4d0
       aa = dsqrt(a)
       b = 0.6d0
       bb = dsqrt(b)
       c = (1.0d0-aa*bb)/b
       d = 1.0d0/dsqrt(0.6d0*1.75d0)
       f1 = c*(h-a)+aa*dsqrt(h-a)/h
       f2 = dsqrt(0.25d0+h-0.25d0/h)
       f3 = d*dsqrt((h-a)*(h+0.75d0)/h)
       f4 = (5.0d0/3.0d0)-(2.0d0/3.0d0)/h
       write(4,104)  h,f1,f2,f3,f4
      end do
      close (4)

 102  format (1p2e20.10)
 103  format (1p3e20.10)
 104  format (1p5e20.10)
      end
