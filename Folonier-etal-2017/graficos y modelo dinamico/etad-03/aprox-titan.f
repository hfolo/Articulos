      implicit real*8 (a-h,k-z)
      real*8 zc(7),zs(7),worb(3),wrot(3),mzi(3)
      real*8 f1(-1:1,-1:1),f2(-1:1,-1:1),f3(-1:1,-1:1)
      real*8 f4(-1:1,-1:1),f5(-1:1,-1:1)
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /cte/   g,std,sty,gw,aukm
      common /ext/   mext,qext,next,perext,perdext,ldot0,ldot1
      common /body/  mtot,ejean,emacl
      common /elem/  a,e,n,per,perd,dadt0,dedt0
      common /core/  dc,rc,mc,gravc,cc,etac,gamac,ec,lamc
      common /shell/ ds,rs,ms,gravs,cs,etas,gamas,es,lams,ck,zhl
      common /torq/  tid,kag,mu,atm,tcc,tcs,tsc,tss

      open (1,file='nu-equil-12.dat')
      open (2,file='dwdt-12.dat')
      open (3,file='orbita-12.dat')
c      zhl = 1.5d1
c      zhl = 1.78d2
      zhl = 2.5d2
      call ctes
      gamac    = 1.0d-6
      p1       = n/gamac
      igama    = 1000
      lgamsmin = -1.0d1
      lgamsmax = -5.0d0
      dlgam    = (lgamsmax-lgamsmin)/dfloat(igama)
c
      do ilazo=0,igama
       gamas = 10.0d0**(lgamsmax-dfloat(ilazo)*dlgam)
       p2 = n/gamas
       zlp2 = dlog10(p2)
       zlg2 = dlog10(gamas)
c
       call aprox(gamac,gamas,zc,zs)
       zsequ    = zs(1)*e**2/sty/rad/dos
       zsamp1   = dsqrt(zs(2)**2+zs(3)**2)*e/sty/rad/dos
       zsamp2   = dsqrt(zs(4)**2+zs(5)**2)*e**2/sty/rad/dos
       zsampatm = dsqrt(zs(6)**2+zs(7)**2)/sty/rad/dos
       zssup    = zsequ+zsamp1
       zsinf    = zsequ-zsamp1
       zssupatm = zsequ+zsamp1+zsamp2+zsampatm
       zsinfatm = zsequ-zsamp1-zsamp2-zsampatm
       write(1,106) zlg2,zsequ,zssup,zsinf,zssupatm,zsinfatm
c
       call diss(gamac,gamas,zc,zs,worb,wrot,wgc,wfric,mzi)
       wred=worb(3)+wrot(3)
       wtot=worb(3)+wrot(3)+wgc+wfric
       write(2,109) zlp2,zlg2,worb(3),wrot(3),wgc,wfric,wred,wtot,-wtot
       da = dadt0*worb(3)/gw/std/aukm
       de = dedt0*(dsqrt(uno-e**2)*worb(3)/gw-n*mzi(3))/sty
       write(3,106) zlg2,da,de,-da,-de
      end do

      close(1)
      close(2)
      close(3)
 100  format (1p10d12.4)
 101  format (1p10d14.6)
 103  format (1p3d20.10)
 104  format (1p4d20.10)
 105  format (1p5d20.10)
 106  format (1p6d20.10)
 107  format (1p7d20.10)
 108  format (1p8d20.10)
 109  format (1p9d20.10)
      end program

C****************************************************************************
C****************************************************************************
      subroutine ctes
      implicit real*8 (a-h,k-z)
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /cte/   g,std,sty,gw,aukm
      common /ext/   mext,qext,next,perext,perdext,ldot0,ldot1
      common /body/  mtot,ejean,emacl
      common /elem/  a,e,n,per,perd,dadt0,dedt0
      common /core/  dc,rc,mc,gravc,cc,etac,gamac,ec,lamc
      common /shell/ ds,rs,ms,gravs,cs,etas,gamas,es,lams,ck,zhl
      common /torq/  tid,kag,mu,atm,tcc,tcs,tsc,tss

c     numeros
      cero   = 0.0d0
      uno    = 1.0d0
      dos    = 2.0d0
      tres   = 3.0d0
      cuatro = 4.0d0
      cinco  = 5.0d0
      twopi  = 8.0d0*datan(1.0d0)
      pi     = twopi/dos
      rad    = pi/180.0d0

c     constantes del problema
      g   = 6.67384d-20   ! Gravitational constant [km^3/kg s^2]
      sty = uno/(60.0d0*60.0d0*24.0d0*365.0d0)
      std = uno/(60.0d0*60.0d0*24.0d0)
      gw  = 1.0d-3                               ! dE/dt  [GW]
c      gw  = 1.0d6/(2.0*twopi*(2.575d6)**2)        ! dE/dt / Sup_sphere [W/m^2]
      aukm = 1.495978707d8    ! AU in kilometer

c     Particle basic data
      mext    = 5.6832592d26    ! mass [kg]
      qext    = 1.6d4           ! Dissipation
      next    = 6.713428d-9     ! Mean motion [1/s]
      perext  = twopi/next      ! Orbital period [s]
      perdext = perext/8.64d4   ! Orbital period [day]
      ldot0   = 1.65d19         ! [kg km^2/s]
      ldot1   = 1.35d19         ! [kg km^2/s]   ! Tokano
c      ldot1   = 0.82d18         ! [kg km^2/s]   ! Charnay

c     Extended body basic data
      mtot = 13.45d22       ! mass [kg]
      rtot = 2.575d3        ! mean radius [km]
      
c     Mean radius
c      zhl = 178.0d0
      ric = 2.084d3   ! mean radius [km] Rock core (inner core)
      roc = 2.286d3   ! mean radius [km] Ice core (outer core)
      rl  = roc+zhl   ! mean radius [km] Ocean (liquid)
      rk  = rtot      ! mean radius [km] Crust

c     Density
      doc = 1.30d12    ! density [kg/km^3] Ice core (outer core)
      dl  = 1.07d12    ! density [kg/km^3] Ocean (liquid)
      dk  = 0.951d12   ! density [kg/km^3] Crust
      dic = 13.45d22/(cuatro*pi/tres*ric**3)
      dic = dic-doc*(roc**3-ric**3)/ric**3-dl*(rl**3-roc**3)/ric**3
      dic = dic-dk*(rk**3-rl**3)/ric**3

c     Mass
      mic = cuatro*pi/tres*dic*ric**3            ! mass [kg] Rock core (inner core)
      moc = cuatro*pi/tres*doc*(roc**3-ric**3)   ! mass [kg] Ice core (outer core)
      ml  = cuatro*pi/tres*dl*(rl**3-roc**3)     ! mass [kg] Ocean (liquid)
      mk  = cuatro*pi/tres*dk*(rk**3-rl**3)      ! mass [kg] Crust

c     Moment of inertial
      cic = 0.4d0*mic*ric**2                            ! moment of inertial [kg km^2] Rock core (inner core)
      coc = 0.4d0*moc*(roc**5-ric**5)/(roc**3-ric**3)   ! moment of inertial [kg km^2] Ice core (outer core)
      cl  = 0.4d0*ml*(rl**5-roc**5)/(rl**3-roc**3)      ! moment of inertial [kg km^2] Ocean (liquid)
      ck  = 0.4d0*mk*(rk**5-rl**5)/(rk**3-rl**3)        ! moment of inertial [kg km^2] Crust

c     Core basic data
      rc    = roc                           ! mean radius [km]
      mc    = mic+moc                       ! mass [kg]
      dc    = mc/(cuatro*pi/tres*rc**3)     ! density [kg/km^3]
      gravc = g*mc/rc**2                    ! gravity [km/s^2]
      cc    = cic+coc                       ! moment of inertial [kg km^2]

c     Shell basic data
      rs    = rk                                  ! mean radius [km]
      ms    = ml+mk                               ! mass [kg]
      ds    = ms/(cuatro*pi/tres*(rs**3-rc**3))   ! density [kg/km^3]
      gravs = g*(mc+ms)/rs**2                     ! gravity [km/s^2]
      cs    = cl+ck                               ! moment of inertial [kg km^2]

c     Layers
      call layer4(ric/rs,roc/rs,rl/rs,doc/dic,dl/dic,dk/dic,hc,hs)

c     orbital elements
      a     = 1.22195d6                   ! Semi-major axis [km]
      e     = 2.7587d-2                   ! Eccentricity
      n     = dsqrt(g*(mtot+mext)/a**3)   ! Mean motion [1/s]
      per   = twopi/n                     ! Orbital period [s]
      perd  = per/8.64d4                  ! Orbital period [day]
      dadt0 = dos*a**2/(g*mtot*mext)
      dedt0 = a*dsqrt(uno-e**2)/(g*mtot*mext*e)

c     Flattenings
      ejean = 3.75d0*(mext/mtot)*(rs/a)**3    ! Jeans homogeneous flattening
      emacl = 1.25d0*n**2*rs**3/(g*mtot)      ! MacLaurin homogeneous flattening [s^2]
      ec    = hc*ejean                        ! equatorial flattening
      es    = hs*ejean                        ! equatorial flattening
      
c     parameters
      lams = cero
      lamc = cero
      eta  = 1.0d3*1.0d-3   ! viscosity [Pa s]  !!!!!  [Pa s]= 1.0d3*[kg/km s]
      etah = eta/zhl        ! (eta/d)_max
c
      tid = 5.625d0*g*mext**2*rs**3/(mtot*a**6)
      kag = 3.2d1*pi**2*g/7.5d1*(uno+lamc)*(uno+lams)*ec*es*dc*ds*rc**5
      mu  = 8.0d0*pi/tres*rc**4*etah
      atmk= dos*ldot1*next/ck
      atms= dos*ldot1*next/cs
      atm = atmk
c      kag  = cero
c      mu   = cero

c     coeficientes constantes en la marea
      tcc = tid*hc
      tcs = cero
      tsc = tid*hc*rc**5/(rs**5-rc**5)
      tss = tid*hs*rs**5/(rs**5-rc**5)
      
 201  format (1p11d20.10)
      return
      end

C****************************************************************************
C****************************************************************************
      subroutine aprox(gc,gs,zc,zs)
      implicit real*8 (a-h,k-z)
      real*8 zl0(2),zl1(4),zl2(4),zc(7),zs(7)
      real*8 zd(2,2),zt(2,2),zq(2,2)
      real*8 zd1(4,4),zd2(4,4)
      real*8 zp(2),zp1(4),zp2(4),zr(2),zr2(4)
      real*8 zdinv(2,2),zd1inv(4,4),zd2inv(4,4)
      real*8 zlatm(4),zdatm(4,4),zdatminv(4,4),zpatm(4)

      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /cte/   g,std,sty,gw,aukm
      common /ext/   mext,qext,next,perext,perdext,ldot0,ldot1
      common /elem/  a,e,n,per,perd,dadt0,dedt0
      common /core/  dc,rc,mc,gravc,cc,etac,gamac,ec,lamc
      common /shell/ ds,rs,ms,gravs,cs,etas,gamas,es,lams,ck,zhl
      common /torq/  tid,kag,mu,atm,tcc,tcs,tsc,tss
c
      p1 = n/gc
      p2 = n/gs
c
      t11  = dos*tcc/gc
      t12  = dos*tcs/gc
      t21  = dos*tsc/gs*cs/ck
      t22  = dos*tss/gs*cs/ck
      zk1  = dos*kag/gc/cc
      zk2  = dos*kag/gs/ck
      zmu1 = mu/gc/cc
      zmu2 = mu/gs/ck
c       
      zd(1,1)= (t11+zk1/(uno+lamc)+zmu1*gc)
      zd(1,2)=-(t12+zk1/(uno+lams)+zmu1*gs)
      zd(2,1)=-(t21+zk2/(uno+lamc)+zmu2*gc)
      zd(2,2)= (t22+zk2/(uno+lams)+zmu2*gs)
      call inv2(zd,zdinv)
c
      zt(1,1)= t11
      zt(1,2)=-t12
      zt(2,1)=-t21
      zt(2,2)= t22
c
      zq(1,1)= tres*(dos+p1**2+p1**4)/(dos*(uno+p1**2)**2)
      zq(1,2)=-tres*p1/((uno+p1**2)**2)
      zq(2,1)=-tres*(dos+p2**2+p2**4)/(dos*(uno+p2**2)**2)
      zq(2,2)= tres*p2/((uno+p2**2)**2)
c       
      zd1(1,1)= zd(1,1)
      zd1(1,2)= n
      zd1(1,3)= zd(1,2)
      zd1(1,4)= cero
      zd1(2,1)=-n
      zd1(2,2)= zd(1,1)
      zd1(2,3)= cero
      zd1(2,4)= zd(1,2)
      zd1(3,1)= zd(2,1)
      zd1(3,2)= cero
      zd1(3,3)= zd(2,2)
      zd1(3,4)= n
      zd1(4,1)= cero
      zd1(4,2)= zd(2,1)
      zd1(4,3)=-n
      zd1(4,4)= zd(2,2)
      call inv4(zd1,zd1inv)
c       
      zd2(1,1)= zd(1,1)
      zd2(1,2)= dos*n
      zd2(1,3)= zd(1,2)
      zd2(1,4)= cero
      zd2(2,1)=-dos*n
      zd2(2,2)= zd(1,1)
      zd2(2,3)= cero
      zd2(2,4)= zd(1,2)
      zd2(3,1)= zd(2,1)
      zd2(3,2)= cero
      zd2(3,3)= zd(2,2)
      zd2(3,4)= dos*n
      zd2(4,1)= cero
      zd2(4,2)= zd(2,1)
      zd2(4,3)=-dos*n
      zd2(4,4)= zd(2,2)
      call inv4(zd2,zd2inv)
c       
      zdatm(1,1)= zd(1,1)
      zdatm(1,2)= dos*next
      zdatm(1,3)= zd(1,2)
      zdatm(1,4)= cero
      zdatm(2,1)=-dos*next
      zdatm(2,2)= zd(1,1)
      zdatm(2,3)= cero
      zdatm(2,4)= zd(1,2)
      zdatm(3,1)= zd(2,1)
      zdatm(3,2)= cero
      zdatm(3,3)= zd(2,2)
      zdatm(3,4)= dos*next
      zdatm(4,1)= cero
      zdatm(4,2)= zd(2,1)
      zdatm(4,3)=-dos*next
      zdatm(4,4)= zd(2,2)
      call inv4(zdatm,zdatminv)
c       
      zp(1) = 1.2d1*p1/(uno+p1**2)
      zp(2) = 1.2d1*p2/(uno+p2**2)
c
      zp1(1) = t11*p1/(uno+p1**2)   -t12*p2/(uno+p2**2)
      zp1(2) = t11*p1**2/(uno+p1**2)-t12*p2**2/(uno+p2**2)
      zp1(3) =-t21*p1/(uno+p1**2)   +t22*p2/(uno+p2**2)
      zp1(4) =-t21*p1**2/(uno+p1**2)+t22*p2**2/(uno+p2**2)
      zp1 = cuatro*zp1
c
      zp2(1)= t11*p1/(uno+4.*p1**2)-t12*p2/(uno+4.*p2**2)
      zp2(2)=dos*(t11*p1**2/(uno+4.*p1**2)-t12*p2**2/(uno+4.*p2**2))
      zp2(3)=-t21*p1/(uno+4.*p1**2)+t22*p2/(uno+4.*p2**2)
      zp2(4)=dos*(-t21*p1**2/(uno+4.*p1**2)+t22*p2**2/(uno+4.*p2**2))
      zp2 = 1.7d1*zp2
c
      zpatm(1)=cero
      zpatm(2)=cero
      zpatm(3)=cero
      zpatm(4)=uno
      zpatm= dos*atm/gs*zpatm
c
      zl1 = matmul(zd1inv,zp1)
c
      zr(1) = zq(1,1)*zl1(1)+zq(1,2)*zl1(2)
      zr(2) = zq(2,1)*zl1(3)+zq(2,2)*zl1(4)
c
      zcoe1 = zq(1,1)*zl1(1)-zq(1,2)*zl1(2)
      zcoe2 = zq(2,1)*zl1(3)-zq(2,2)*zl1(4)
      zcoe3 = zq(1,2)*zl1(1)+zq(1,1)*zl1(2)
      zcoe4 = zq(2,2)*zl1(3)+zq(2,1)*zl1(4)
      zr2(1) = t11*zcoe1-t12*zcoe2
      zr2(2) = t11*zcoe3-t12*zcoe4
      zr2(3) =-t21*zcoe1+t22*zcoe2
      zr2(4) =-t21*zcoe3+t22*zcoe4
c
      zl2   = matmul(zd2inv,zp2)-matmul(zd2inv,zr2)
      zl0   = matmul(zdinv,matmul(zt,zp))-matmul(zdinv,matmul(zt,zr))
      zlatm = matmul(zdatminv,zpatm)
c
      zc(1) = gc*zl0(1)     ! nu_c = Bc0 e**2   + Cc1 e cos(l)   + Sc1 e sin(l)   + Cc2 e**2 cos(2l)   + Sc2 e**2 sin(2l)
      zc(2) = gc*zl1(1)     ! nu_c = zc(1) e**2 + zc(2) e cos(l) + zc(3) e sin(l) + zc(4) e**2 cos(2l) + zc(5) e**2 sin(2l)
      zc(3) = gc*zl1(2)
      zc(4) = gc*zl2(1)
      zc(5) = gc*zl2(2)
      zc(6) = gc*zlatm(1)
      zc(7) = gc*zlatm(2)
c
      zs(1) = gs*zl0(2)     ! nu_s = Bs0 e**2   + Cs1 e cos(l)   + Ss1 e sin(l)   + Cs2 e**2 cos(2l)   + Ss2 e**2 sin(2l)
      zs(2) = gs*zl1(3)     ! nu_s = zs(1) e**2 + zs(2) e cos(l) + zs(3) e sin(l) + zs(4) e**2 cos(2l) + zs(5) e**2 sin(2l)
      zs(3) = gs*zl1(4)
      zs(4) = gs*zl2(3)
      zs(5) = gs*zl2(4)
      zs(6) = gs*zlatm(3)
      zs(7) = gs*zlatm(4)

      return
      end

C****************************************************************************
C****************************************************************************

      subroutine func(gc,gs,zc,zs,f1,f2,f3,f4,f5)
      implicit real*8 (a-h,k-z)
      real*8 zc(7),zs(7)
      real*8 f1(-1:1,-1:1),f2(-1:1,-1:1),f3(-1:1,-1:1)
      real*8 f4(-1:1,-1:1),f5(-1:1,-1:1)
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /elem/  a,e,n,per,perd,dadt0,dedt0
c
      p1 = n/gc
      p2 = n/gs
      zp1 = uno+p1**2
      zp2 = uno+p2**2
c
      b0c = zc(1)
      c1c = zc(2)
      s1c = zc(3)
      b0s = zs(1)
      c1s = zs(2)
      s1s = zs(3)

c     f_{j,k}^{1}
      f1(0,0)  = dos*p1*b0c/n
      f1(0,1)  = p1/(dos*zp1)
      f1(0,-1) =-4.9d1*p1/(dos*zp1)
      f1(1,0)  =-p1*c1c/(dos*n)
      f1(1,1)  = cero
      f1(1,-1) = 3.5d0*p1*(c1c*(uno-p1**2)+dos*s1c*p1)/(n*zp1**2)
      f1(-1,0) = 7.0d0*p1*c1c/(dos*n)
      f1(-1,1) =-0.5d0*p1*(c1c*(uno-p1**2)+dos*s1c*p1)/(n*zp1**2)
      f1(-1,-1)= cero

c     f_{j,k}^{2}
      f2(0,0)  = dos*p2*b0s/n
      f2(0,1)  = p2/(dos*zp2)
      f2(0,-1) =-4.9d1*p2/(dos*zp2)
      f2(1,0)  =-p2*c1s/(dos*n)
      f2(1,1)  = cero
      f2(1,-1) = 3.5d0*p2*(c1s*(uno-p2**2)+dos*s1s*p2)/(n*zp2**2)
      f2(-1,0) = 7.0d0*p2*c1s/(dos*n)
      f2(-1,1) =-0.5d0*p2*(c1s*(uno-p2**2)+dos*s1s*p2)/(n*zp2**2)
      f2(-1,-1)= cero

c     f_{j,k}^{1,1}
      f3(0,0)  = p1*(c1c*c1c+s1c*s1c)/n
      f3(0,1)  = cero
      f3(0,-1) = cero
      f3(1,0)  =-0.5d0*s1c
      f3(1,1)  = cero
      f3(1,-1) = 3.5d0*(-c1c*p1+s1c)/zp1
      f3(-1,0) =-3.5d0*s1c
      f3(-1,1) = 0.5d0*(-c1c*p1+s1c)/zp1
      f3(-1,-1)= cero

c     f_{j,k}^{1,2}
      f4(0,0)  = p1*(c1c*c1s+s1c*s1s)/n
      f4(0,1)  = cero
      f4(0,-1) = cero
      f4(1,0)  =-0.5d0*s1s
      f4(1,1)  = cero
      f4(1,-1) = 3.5d0*(-c1s*p1+s1s)/zp1
      f4(-1,0) =-3.5d0*s1s
      f4(-1,1) = 0.5d0*(-c1s*p1+s1s)/zp1
      f4(-1,-1)= cero

c     f_{j,k}^{2,2}
      f5(0,0)  = p2*(c1s*c1s+s1s*s1s)/n
      f5(0,1)  = cero
      f5(0,-1) = cero
      f5(1,0)  =-0.5d0*s1s
      f5(1,1)  = cero
      f5(1,-1) = 3.5d0*(-c1s*p2+s1s)/zp2
      f5(-1,0) =-3.5d0*s1s
      f5(-1,1) = 0.5d0*(-c1s*p2+s1s)/zp2
      f5(-1,-1)= cero
c
      f1 = f1*e**2  ! f_{j,k}^{2}
      f2 = f2*e**2  ! f_{j,k}^{2}
      f3 = f3*e**2  ! f_{j,k}^{1,1}
      f4 = f4*e**2  ! f_{j,k}^{2,1}
      f5 = f5*e**2  ! f_{j,k}^{2,2}

!       write(*,201) f5(0,0),f5(0,1),f5(0,-1)
!       write(*,201) f5(1,0),f5(1,1),f5(1,-1)
!       write(*,201) f5(-1,0),f5(-1,1),f5(-1,-1)
!       stop
!  201  format (1p11d20.10)
      return
      end

C****************************************************************************
C****************************************************************************

      subroutine diss(gc,gs,zc,zs,worb,wrot,wgc,wfric,mzi)
      implicit real*8 (a-h,k-z)
      real*8 zc(7),zs(7),worb(3),wrot(3),mzi(3)
      real*8 f1(-1:1,-1:1),f2(-1:1,-1:1),f3(-1:1,-1:1)
      real*8 f4(-1:1,-1:1),f5(-1:1,-1:1)
      integer j,k
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /cte/   g,std,sty,gw,aukm
      common /elem/  a,e,n,per,perd,dadt0,dedt0
      common /core/  dc,rc,mc,gravc,cc,etac,gamac,ec,lamc
      common /shell/ ds,rs,ms,gravs,cs,etas,gamas,es,lams,ck,zhl
      common /torq/  tid,kag,mu,atm,tcc,tcs,tsc,tss
c
      p1 = n/gc
      p2 = n/gs
      zp1 = uno+p1**2
      zp2 = uno+p2**2
c
      g1 = cero
      g2 = cero
      g3 = cero
      g4 = cero
      g5 = cero
      g6 = cero
      g7 = cero
      call func(gc,gs,zc,zs,f1,f2,f3,f4,f5)
      do j=-1,1
       do k=-1,1
        g1=g1+dfloat(k-2)*n*f1(j,k)
        g2=g2+dfloat(k-2)*n*f2(j,k)
        g3=g3+dos*n*f1(j,k)+f3(j,k)
        g4=g4+dos*n*f1(j,k)+f4(j,k)
        g5=g5+dos*n*f2(j,k)+f5(j,k)
        g6=g6+f1(j,k)
        g7=g7+f2(j,k)
       end do
      end do
c
      g8 = tres*gc*n**2*e**2/(n**2+gc**2)
      g9 = tres*gs*n**2*e**2/(n**2+gs**2)
c
      worb(1) = -tcc*cc/cuatro*(g1+g8)*gw
      worb(2) = -tss*cs/cuatro*(g2+g9)*gw+tsc*cs/cuatro*(g1+g8)*gw
      worb(3) = (worb(1)+worb(2))
c
      wrot(1) = -tcc*cc/cuatro*g3*gw
      wrot(2) = -tss*cs/cuatro*g5*gw+tsc*cs/cuatro*g4*gw
      wrot(3) = (wrot(1)+wrot(2))
c
      coe1 = ((zc(2)-zs(2))*zc(2)+(zc(3)-zs(3))*zc(3))/(gc*(uno+lamc))
      coe2 = ((zc(2)-zs(2))*zs(2)+(zc(3)-zs(3))*zs(3))/(gs*(uno+lams))
      wgc  = -kag*e**2/cuatro*(coe1-coe2)*gw
c
      wfric = -mu*e**2/8.0d0*((zc(2)-zs(2))**2+(zc(3)-zs(3))**2)*gw
c
      mzi(1) = tcc*cc/dos*g6
      mzi(2) = tss*cs/dos*g7-tsc*cs/dos*g6
      mzi(3) = mzi(1)+mzi(2)

 201  format (1p11d20.10)

      return
      end

C****************************************************************************
C****************************************************************************
      subroutine angxi(nuc,nus,gc,gs,xi)
      implicit real*8 (a-h,k-z)
      integer iorder,i,j
      common /num/ cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /core/  dc,rc,mc,gravc,cc,etac,gamac,ec,lamc
      common /shell/ ds,rs,ms,gravs,cs,etas,gamas,es,lams,ck,zhl

      sin2sc = dos*gc*nuc/(gc**2+nuc**2)            ! sin(2sigma_c0)
      cos2sc = (gc**2-nuc**2)/(gc**2+nuc**2)        ! cos(2sigma_c0)
      eps0c  = datan2(sin2sc,uno+dos*lamc+cos2sc)   ! vartheta_c
      sin2ss = dos*gs*nus/(gs**2+nus**2)            ! sin(2sigma_s0)
      cos2ss = (gs**2-nus**2)/(gs**2+nus**2)        ! cos(2sigma_s0)
      eps0s  = datan2(sin2ss,uno+dos*lams+cos2ss)   ! vartheta_s
      xi     = (eps0s-eps0c)/dos                    ! xi

      return
      end

      subroutine inv2(E,F)
      implicit real*8(a-h,l-z)
      parameter (imax=2)
      integer n
      real*8 E(imax,imax),F(imax,imax)
      double precision aaa(imax,imax), ccc(imax,imax)
      double precision LLL(imax,imax), UUU(imax,imax)
      double precision coeff,bbb(imax),ddd(imax),xxx(imax)

!     calculation of the matrix F inversa de E
!============================================================
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
      n = 2
      do i=1,n
       do j=1,n
        aaa(i,j) = E(i,j)
       end do
      end do

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
      LLL=0.0
      UUU=0.0
      bbb=0.0

! step 1: forward elimination
      do k=1, n-1
       do i=k+1,n
        coeff=aaa(i,k)/aaa(k,k)
        LLL(i,k) = coeff
        do j=k+1,n
         aaa(i,j) = aaa(i,j)-coeff*aaa(k,j)
        end do
       end do
      end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
      do i=1,n
       LLL(i,i) = 1.0
      end do
! U matrix is the upper triangular part of A
      do j=1,n
       do i=1,j
        UUU(i,j) = aaa(i,j)
       end do
      end do

! Step 3: compute columns of the inverse matrix C
      do k=1,n
       bbb(k)=1.0
       ddd(1) = bbb(1)
! Step 3a: Solve Ld=b using the forward substitution
       do i=2,n
        ddd(i)=bbb(i)
        do j=1,i-1
         ddd(i) = ddd(i) - LLL(i,j)*ddd(j)
        end do
       end do
! Step 3b: Solve Ux=d using the back substitution
       xxx(n)=ddd(n)/UUU(n,n)
       do i = n-1,1,-1
        xxx(i) = ddd(i)
        do j=n,i+1,-1
         xxx(i)=xxx(i)-UUU(i,j)*xxx(j)
        end do
        xxx(i) = xxx(i)/uuu(i,i)
       end do
! Step 3c: fill the solutions x(n) into column k of C
       do i=1,n
        ccc(i,k) = xxx(i)
       end do
       bbb(k)=0.0
      end do 

      do i=1,n
       do j=1,n
        F(i,j) = ccc(i,j)
       end do
      end do

 301  format(1p5e15.5)
      return
      end

      subroutine inv4(E,F)
      implicit real*8(a-h,l-z)
      parameter (imax=4)
      integer n
      real*8 E(imax,imax),F(imax,imax)
      double precision aaa(imax,imax), ccc(imax,imax)
      double precision LLL(imax,imax), UUU(imax,imax)
      double precision coeff,bbb(imax),ddd(imax),xxx(imax)

!     calculation of the matrix F inversa de E
!============================================================
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
      n = 4
      do i=1,n
       do j=1,n
        aaa(i,j) = E(i,j)
       end do
      end do

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
      LLL=0.0
      UUU=0.0
      bbb=0.0

! step 1: forward elimination
      do k=1, n-1
       do i=k+1,n
        coeff=aaa(i,k)/aaa(k,k)
        LLL(i,k) = coeff
        do j=k+1,n
         aaa(i,j) = aaa(i,j)-coeff*aaa(k,j)
        end do
       end do
      end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
      do i=1,n
       LLL(i,i) = 1.0
      end do
! U matrix is the upper triangular part of A
      do j=1,n
       do i=1,j
        UUU(i,j) = aaa(i,j)
       end do
      end do

! Step 3: compute columns of the inverse matrix C
      do k=1,n
       bbb(k)=1.0
       ddd(1) = bbb(1)
! Step 3a: Solve Ld=b using the forward substitution
       do i=2,n
        ddd(i)=bbb(i)
        do j=1,i-1
         ddd(i) = ddd(i) - LLL(i,j)*ddd(j)
        end do
       end do
! Step 3b: Solve Ux=d using the back substitution
       xxx(n)=ddd(n)/UUU(n,n)
       do i = n-1,1,-1
        xxx(i) = ddd(i)
        do j=n,i+1,-1
         xxx(i)=xxx(i)-UUU(i,j)*xxx(j)
        end do
        xxx(i) = xxx(i)/uuu(i,i)
       end do
! Step 3c: fill the solutions x(n) into column k of C
       do i=1,n
        ccc(i,k) = xxx(i)
       end do
       bbb(k)=0.0
      end do 

      do i=1,n
       do j=1,n
        F(i,j) = ccc(i,j)
       end do
      end do

 401  format(1p5e15.5)
      return
      end

      subroutine layer4(x1,x2,x3,d2,d3,d4,h2,h4)
      implicit real*8(a-h,l-z)
      parameter (imax=4)
      real*8 x(0:imax),d(1:imax+1),h(1:imax)
      real*8 E(imax,imax),F(imax,imax)
      real*8 f1(imax),f2(imax),f3(imax),a(0:imax),b(0:imax),c(0:imax)
      real*8 kf,c20,c22,kl(imax)
      integer i,j,k,ii,n,ifull,iscreen
      double precision aaa(imax,imax), ccc(imax,imax)
      double precision LLL(imax,imax), UUU(imax,imax)
      double precision coeff,bbb(imax),ddd(imax),xxx(imax)
      common /num/ cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
     
c     radios normalizados
      x(0) = cero
      x(1) = x1
      x(2) = x2
      x(3) = x3
      x(4) = uno
c     densidades normalizadas
      d(1) = uno
      d(2) = d2
      d(3) = d3
      d(4) = d4
      d(5) = cero

c     calculo de fn
      n=4
      fn = cero
      do i=1,n
       fn = fn + (d(i)-d(i+1))*x(i)**3
      end do
      c32f = tres/(dos*fn)
      c52f = cinco/(dos*fn)
      c52  = cinco/dos
c     calculation of the matrix E
      E = 0.0d0
      do i=1,n
       do j=1,n        
        if(i.lt.j) then
         E(i,j) = -c32f*(d(j)-d(j+1))*x(i)**3
        end if

        if(i.eq.j) then
         E(i,j)=0.0d0
         do k=i,n
          E(i,j)=E(i,j)+(d(k)-d(k+1))*(x(k)**3-x(i)**3)
         end do
         E(i,j)=-c32f*(d(i)-d(i+1))*x(i)**3+c52-c52f*E(i,j)
        end if

        if(i.gt.j) then
         E(i,j) = -c32f*(d(j)-d(j+1))*x(j)**5/x(i)**2
        end if
       end do
      end do


!     calculation of the matrix F inversa de E
!============================================================
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
      do i=1,n
       do j=1,n
        aaa(i,j) = E(i,j)
       end do
      end do

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
      LLL=0.0
      UUU=0.0
      bbb=0.0

! step 1: forward elimination
      do k=1, n-1
       do i=k+1,n
        coeff=aaa(i,k)/aaa(k,k)
        LLL(i,k) = coeff
        do j=k+1,n
         aaa(i,j) = aaa(i,j)-coeff*aaa(k,j)
        end do
       end do
      end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
      do i=1,n
       LLL(i,i) = 1.0
      end do
! U matrix is the upper triangular part of A
      do j=1,n
       do i=1,j
        UUU(i,j) = aaa(i,j)
       end do
      end do

! Step 3: compute columns of the inverse matrix C
      do k=1,n
       bbb(k)=1.0
       ddd(1) = bbb(1)
! Step 3a: Solve Ld=b using the forward substitution
       do i=2,n
        ddd(i)=bbb(i)
        do j=1,i-1
         ddd(i) = ddd(i) - LLL(i,j)*ddd(j)
        end do
       end do
! Step 3b: Solve Ux=d using the back substitution
       xxx(n)=ddd(n)/UUU(n,n)
       do i = n-1,1,-1
        xxx(i) = ddd(i)
        do j=n,i+1,-1
         xxx(i)=xxx(i)-UUU(i,j)*xxx(j)
        end do
        xxx(i) = xxx(i)/uuu(i,i)
       end do
! Step 3c: fill the solutions x(n) into column k of C
       do i=1,n
        ccc(i,k) = xxx(i)
       end do
       bbb(k)=0.0
      end do 

      do i=1,n
       do j=1,n
        F(i,j) = ccc(i,j)
       end do
      end do

      do i=1,n
       h(i) = 0.0d0
       do j=1,n
        h(i)=h(i)+F(i,j)*x(j)**3
       end do
      end do

      h1 = h(1)
      h2 = h(2)
      h3 = h(3)
      h4 = h(4)
 301  format(1p5e15.5)
      return
      end

