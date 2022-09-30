      program graficador
      implicit real*8 (a-h,k-z)
      integer i,j,itot
      open (1,file='two-layers.dat',status='replace')
c     constantes 
      cero   = 0.0d0
      uno    = 1.0d0
      dos    = 2.0d0
      tres   = 3.0d0
      cuatro = 4.0d0
      twopi  = dos*cuatro*datan(uno)
      pi     = twopi/dos
      rad    = pi/180.0d0
c     intervalo del grafico
      itot = 2000
      xmin = 0.0d0
      xmax = 1.0d0
      dx   = (xmax-xmin)/dfloat(itot)
c     parametros del grafico
      p = 0.2d0
      l = 0.2d0
      ca = (tres*l+dos)*(dos*l+5.0d0*(uno-l)*p**3)-9.0d0*l*(uno-l)*p**5
      k1 = 10.0d0*(l+(1.0d0-l)*p**3)**2/ca
      k2 = 2.0d0*(l+(1.0d0-l)*p**3)*(3.0d0*l+2.0d0+3.0d0*(1.0d0-l)*p**3)/ca

      do i=0,itot
       x  = xmin+dx*dfloat(i)
       if(x.lt.0.3) then
        y = 1.0d0
       else
        y = 0.6
       end if
       write(1,102) x,y
      end do

      close(1)
 102  format (1p2e20.10)
 103  format (1p3e20.10)
 104  format (1p4e20.10)
 105  format (1p5e20.10)
 106  format (1p6e20.10)
 107  format (1p7e20.10)
      end
