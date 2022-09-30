      program graficador
      implicit real*8 (a-h,k-z)
      real*8 y(4),y0(4),ymax(4),ymin(4)
      integer i,j,itot,imed
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /ext/   mext,msum
      common /body/  mtot,rtot,ctot,dtot,grav,ej,em
      common /orb/   a,e,n,per
!
      open (1,file='aprox-atractor-nu.dat')
      open (2,file='aprox-atractor-delta.dat')
      open (3,file='aprox-atractor-erho.dat')
      open (4,file='aprox-atractor-ez.dat')
      open (5,file='aprox-energia.dat')
      open (6,file='aprox-elem.dat')

!     constantes 
      cero   = 0.0d0
      uno    = 1.0d0
      dos    = 2.0d0
      tres   = 3.0d0
      cuatro = 4.0d0
      cinco  = 5.0d0
      twopi  = dos*cuatro*datan(uno)
      pi     = twopi/dos
      rad    = pi/180.0d0

!     constantes del problema
      k   = 1.720209895d-02                       ! Gaussian constant
      g   = k**2                                  ! Gravitational constant [AU^3/M_sun day^2]
      sty = uno/(60.0d0*60.0d0*24.0d0*365.2425d0) ! Second to year
      std = uno/(60.0d0*60.0d0*24.0d0)            ! Second to day
      uam = 1.495978707d11                        ! AU in m
      zmsolkg = 1.98911d30                        ! solar mass in kg
      gwatt = zmsolkg*uam**2*std**3*1.0d-9
      twatt = zmsolkg*uam**2*std**3*1.0d-12

!     basic data
      mext = 5.68326d26/zmsolkg             ! Mass [M_sun]
      mtot = 1.08d20/zmsolkg                ! Mass [M_sun]
      rtot = 252.1d3/uam                    ! Mean radius [AU]
      ctot = 0.4d0*mtot*rtot**2             ! Moment of inertial [M_sun AU^2]
      dtot = mtot/(cuatro*pi/tres*rtot**3)
      grav = g*mtot/rtot**2
      msum = mext+mtot                      ! m1+m2
      gm   = g*msum

!     orbital elements
      a   = 238.02d6/uam                   ! Semi-major axis [AU]
      e   = 0.0045d0                       ! Eccentricity
      n = dsqrt(gm/a**3)                   ! Mean motion [1/day]
      ej = 3.75d0*mext*rtot**3/mtot/a**3
      em = 1.25d0*n**2*rtot**3/(g*mtot)

!     intervalo del grafico
      itot = 2000
      letamin = 0.9d1
      letamax = 1.6d1
      dleta   = (letamax-letamin)/dfloat(itot)
      do i=0,itot
       eta  = 10.0d0**(letamin+dleta*dfloat(i))
       eta  = eta/zmsolkg/std*uam
       gama = dtot*grav*rtot/(dos*eta)
       p    = n/gama
       zet  = eta*zmsolkg*std/uam
       zga  = gama*std

       ymin = 1.0d10
       ymax =-1.0d10
       y0   = cero
       imed = 1000
       xmin = cero
       xmax = twopi
       dx   = (xmax-xmin)/dfloat(imed)
       do j=0,imed
        x = xmin+dx*dfloat(j)
        call atractor (x,gama,y,y0)
        do iy=1,4
         if(y(iy).gt.ymax(iy)) ymax(iy)=y(iy)
         if(y(iy).lt.ymin(iy)) ymin(iy)=y(iy)
        end do
       end do
       write(1,105) dlog10(zet),zga,y0(1)/rad,ymax(1)/rad,ymin(1)/rad
       write(2,105) dlog10(zet),zga,y0(2)/rad,ymax(2)/rad,ymin(2)/rad
       write(3,105) dlog10(zet),zga,y0(3),ymax(3),ymin(3)
       write(4,105) dlog10(zet),zga,y0(4),ymax(4),ymin(4)

       kapa = mext/msum
       al   = uno-tres*kapa*ej
       be   = 8.0d0*kapa*em
       div1 = uno+p**2
       div2 = uno+al**2*p**2

       E0 = g*mext*ctot*e**2/a**3
!      mechanical energy
       Etot0 = 1.5d0*E0*n*p*twatt
       zs1   = (7.0d0+(cuatro+tres*al**2)*p**2)*ej
       zs2   = dos*(uno-al)*(uno-al*p**2)*em
       Etot  = Etot0*(zs1+zs2)/div1/div2
       Etot2 = 10.5d0*E0*ej*n*p/div1*twatt

!      binding energy
       Eint0 = 9.0d0*E0*kapa**2*ej**2*em*n*p*twatt
       zr1   = tres*ej*(uno-p**2)*div2
       zr2   = 8.0d0*em*p**2*(tres-al*p**2)*(uno-al)
       Eint  = Eint0*(zr1-zr2)/(div1*div2)**2
       Eint2 = 1.8d1/7.0d0*kapa**2*ej*em*(uno-p**2)/div1**2*Etot2

!      orbital and rotational energies due dC/dt
       Eorb0 = E0/1.2d1*n*p*twatt
       zr1   = 1.8d1*ej**2*div1*div2**2
       zr2   = al*(1.9d1*al-1.6d1)*p**2+(4.3d1-al*(6.8d1-3.1d1*al))
       zr2   = zr2*p**2+(7.0d0-cuatro*al)*(9.0d0-8.0d0*al)
       zr2   = tres*ej*em*div2*zr2
       zr3   = tres*al**3*p**6-al*(uno-dos*al*(uno+al))*p**4
       zr3   = zr3+tres*(dos+al*(6.0d0+al*(dos+al)))*p**2
       zr3   = zr3-tres*(dos-al)
       zr3   =-8.0d0*em**2*(uno-al)*zr3
       Eorb  = Eorb0*(zr1+zr2+zr3)/(div1*div2)**2
       Eorb2 = 1.0d0/1.4d1*(dos*ej+em)*Etot2
       
       Erot0 = E0*ej*n*p*twatt
       zr1   = tres*ej*(uno+al*p**2)*div2
       zr2   =-8.0d0*em*(uno-al)*p**2*(al**2*p**2+(uno-al*(uno-al)))
       Erot  = Erot0*(zr1+zr2)/div1/div2**2
       Erot2 = 2.0d0/7.0d0*ej*Etot2
       write(5,110) zet,zga,Etot,Eint,Eorb,Erot,Etot2,Eint2,Eorb2,Erot2

       dadt = dos*a**2/(g*mext*mtot)*Etot/twatt*std/sty
       dedt = (uno-e**2)/(dos*a*e)*dadt
!        dadt = 2.1d1*ctot*ej*e**2/(mtot*a)*n*p/(uno+p**2)*std/sty
!        dedt = (uno-e**2)/(dos*a*e)*dadt
       write(6,105) zet,zga,dadt*uam/1.0d3,dedt
      end do

      close(1)
 102  format (1p2e20.10)
 103  format (1p3e20.10)
 104  format (1p4e20.10)
 105  format (1p5e20.10)
 106  format (1p6e20.10)
 107  format (1p7e20.10)
 108  format (1p8e20.10)
 110  format (1p10e20.10)
      end

 !****************************************************************************
!****************************************************************************
      subroutine atractor (ell,gama,y,y0)
      implicit real*8 (a-h,k-z)
      real*8 y(4),y0(4)
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /ext/   mext,msum
      common /body/  mtot,rtot,ctot,dtot,grav,ej,em
      common /orb/   a,e,n,per
!
      ell2 = dos*ell
!
      p    = n/gama
      kapa = mext/msum
      al   = uno-tres*kapa*ej
      be   = uno-(uno-al)/cuatro
      div  = uno+al**2*p**2
!
      b0 = 1.2d1*n*e**2/(uno+p**2)*(uno+al*p**2)/div
      b1 =-1.2d1*kapa*n*al*p**2*e*ej/div
      b2 = 1.2d1*kapa*n*p*e*ej/div
!
      d0 = tres*p*e**2/(uno+p**2)*(dos+(uno+al)*p**2)/div
      d1 =-dos*p*e/div
      d2 =-dos*p**2*al*e/div
!
      e0 = ej*(uno+tres*e**2/dos-d1**2-d2**2)
      e1 = tres*ej*e/(uno+p**2)
      e2 = tres*ej*p*e/(uno+p**2)
!
      z0 = ej/dos*(uno+tres*e**2/dos)
      z0 = z0+em*(uno+b0/n+(b1**2+b2**2)/(8.0d0*n**2))
      coe = gama/(n**2+gama**2)
      z1 = coe*(1.5d0*e*ej*gama+em*(b1*gama-b2*n)/n)
      z2 = coe*(1.5d0*e*ej*n+em*(b1*n+b2*gama)/n)
!
      f1 = (uno-al)*n/dos*(3*e*d2+(d1*e2+d2*e1)/ej)
      f2 =-(uno-al)*n/dos*(3*e*d1+(d1*e1-d2*e2)/ej)
      f3 = gama/dos*(3*e*d1-(d1*e1-d2*e2)/ej)
      f4 = gama/dos*(3*e*d2-(d1*e2+d2*e1)/ej)
      g1 = -2.5d0*n*e**2+f1/dos-f3
      g2 = f2/dos-f4
!
      d3 = (gama*g1-dos*be*g2)/(cuatro*be**2*n**2+gama**2)
      d4 = (gama*g2+dos*be*g1)/(cuatro*be**2*n**2+gama**2)
!
      b3 = (uno-al)*n*d4+f1
      b4 =-(uno-al)*n*d3+f2
!
      b3 = (uno-al)*n*d4+f1
      b4 =-(uno-al)*n*d3+f2
!
      coe = gama/(cuatro*n**2+gama**2)
      e3 = coe*ej*(gama*(4.5d0*e**2-d1**2+d2**2)+cuatro*n*d1*d2)
      e4 = coe*ej*(dos*n*(4.5d0*e**2-d1**2+d2**2)-dos*gama*d1*d2)
!
      coe1 = (gama*b3-dos*n*b4)/n
      coe1 = coe1 + (gama*(b1**2-b2**2)-cuatro*n*b1*b2)/(8.0d0*n**2)
      coe2 = (gama*b4+dos*n*b3)/n
      coe2 = coe2 + (n*(b1**2-b2**2)+cuatro*gama*b1*b2)/(8.0d0*n**2)
      z3 = coe*(2.25d0*e**2*ej*gama+coe1*em)
      z4 = coe*(4.5d0*e**2*ej*n+coe2*em)
      
!       y(1)=b0+b1*dcos(ell)+b2*dsin(ell)+b3*dcos(ell2)+b4*dsin(ell2)
!       y(2)=d0+d1*dcos(ell)+d2*dsin(ell)+d3*dcos(ell2)+d4*dsin(ell2)
!       y(3)=e0+e1*dcos(ell)+e2*dsin(ell)+e3*dcos(ell2)+e4*dsin(ell2)
!       y(4)=z0+z1*dcos(ell)+z2*dsin(ell)+z3*dcos(ell2)+z4*dsin(ell2)
      y(1)=b0+b1*dcos(ell)+b2*dsin(ell)
      y(2)=d0+d1*dcos(ell)+d2*dsin(ell)
      y(3)=e0+e1*dcos(ell)+e2*dsin(ell)
      y(4)=z0+z1*dcos(ell)+z2*dsin(ell)
!
      y0(1)=b0
      y0(2)=d0
      y0(3)=e0
      y0(4)=z0
!
      return     
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      Enceladus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      mtot = 1.08d20/zmsolkg          ! Mass [M_sun]
!      rtot = 252.1d3/uam              ! Mean radius [AU]
!      a   = 238.02d6/uam              ! Semi-major axis [AU]
!      e   = 0.0045d0                  ! Eccentricity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      Mimas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      mtot = 0.379d20/zmsolkg               ! Mass [M_sun]
!      rtot = 197.49d3/uam                   ! Mean radius [AU]
!      a   = 185.52d6/uam                    ! Semi-major axis [AU]
!      e   = 0.01986d0                       ! Eccentricity
