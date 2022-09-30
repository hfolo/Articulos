      implicit real*8 (a-h,k-z)
      real*8 y(5),w(8),ymin(4),ymax(4),ysum(4),wsum(8)
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /cte/   g,k,gm,std,sty,uam,zmsolkg,gwatt,twatt
      common /ext/   mext,next,perext,msum,mred
      common /body/  mtot,rtot,ctot,dtot,grav,eta,gama,ej,em
      common /orb/   a,e,n,per

      open (1,file='atractor-nu.dat')
      open (2,file='atractor-delta.dat')
      open (3,file='atractor-erho.dat')
      open (4,file='atractor-ez.dat')
      open (5,file='energia.dat')

      ip      = 50
      letamin = 0.9d1
      letamax = 1.6d1
      dleta   = (letamax-letamin)/dfloat(ip)
      do i=0,50
       call ctes
       eta  = 10.0d0**(letamin+dleta*dfloat(i))
       eta  = eta/zmsolkg/std*uam
       gama = dtot*grav*rtot/(dos*eta)
       p    = n/gama

!      lazo principal de integracion y precision de la integracion.
       time0  = 0.0d0
       tfinal = 7.0d4*per      ! final integration time (initial is = 0)
       htime  = per/1.0d2  ! output interval
       h      = per/1.0d2  ! initial timestep
       ll     = 13

!      Condiciones iniciales
       tau = time0
       y(1) = 1.2d1*n*e**2/(uno+p**2)                                  ! nu inicial
       y(2) = 6.0d0*p*e**2/(uno+p**2)                                  ! delta inicial
       y(3) = ej*(uno+1.5d0*e**2-4.0d0*p**2*e**2/(uno+p**2))           ! varepsilon_rho inicial
       y(4) = ej/dos*(uno+1.5d0*e**2)+em*(uno+1.2d1*e**2/(uno+p**2))   ! varepsilon_z inicial
       y(5) = tau                                                      ! t0

       ymin = 1.0d10
       ymax =-1.0d10
       wsum = cero
       asum = cero
       icont = 0
       do 50 while (tau.lt.tfinal)
        call bs (y,htime,h,ll,5)
        tau = tau + htime
        y(2) = dmod(y(2),pi)
        if (y(2).lt.0.0d0)    y(2)=y(2)+pi
        if (y(2).gt.pi/dos)   y(2)=y(2)-pi
        if(tau.gt.tfinal-1.0d4*per) then
         call energia (y,w)
         do iy=1,4
          if(y(iy).gt.ymax(iy)) ymax(iy)=y(iy)
          if(y(iy).lt.ymin(iy)) ymin(iy)=y(iy)
          ysum(iy)=ysum(iy)+y(iy)
          wsum(iy)=wsum(iy)+w(iy)
          wsum(iy+4)=wsum(iy+4)+w(iy+4)
         end do
         icont=icont+1
        end if
 50    continue
       do iy=1,4
        ysum(iy)=ysum(iy)/dfloat(icont)
        wsum(iy)=wsum(iy)/dfloat(icont)
        wsum(iy+4)=wsum(iy+4)/dfloat(icont)
       end do
       zet = eta*zmsolkg*std/uam
       zga = gama*std
      write(1,110) zet,zga,ysum(1)/rad,ymax(1)/rad,ymin(1)/rad
      write(2,110) zet,zga,ysum(2)/rad,ymax(2)/rad,ymin(2)/rad
      write(3,110) zet,zga,ysum(3),ymax(3),ymin(3)
      write(4,110) zet,zga,ysum(4),ymax(4),ymin(4)
      write(5,110) zet,zga,-wsum
      write(*,101) zet,zga,ysum(1)/rad,ymax(1)/rad,ymin(1)/rad
      end do

      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
 100  format (1p10d12.4)
 101  format (1p10d14.6)
 109  format (1p9d20.10)
 110  format (1p10d20.10)
      end program

!****************************************************************************
!****************************************************************************
      subroutine force (yc,f)
      implicit real*8 (a-h,k-z)
      real*8 yc(5),f(5)
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /cte/   g,k,gm,std,sty,uam,zmsolkg,gwatt,twatt
      common /ext/   mext,next,perext,msum,mred
      common /body/  mtot,rtot,ctot,dtot,grav,eta,gama,ej,em
      common /orb/   a,e,n,per
!
      nu   = yc(1)
      del  = yc(2)
      ep   = yc(3)
      epz  = yc(4)
      time = yc(5)
!
      co1 = dcos(n*time)
      si1 = dsin(n*time)
      co2 = co1*co1-si1*si1
      si2 = si1*co1+co1*si1
      co3 = co1*co2-si1*si2
      si3 = si1*co2+co1*si2
      co4 = co1*co3-si1*si3
      si4 = si1*co3+co1*si3
!
      ar3 = uno+tres*e*co1
      ar3 = ar3+tres/dos*e**2*(uno+tres*co2)
      ar3 = ar3+e**3/8.0d0*(2.7d1*co1+5.3d1*co3)
      ar3 = ar3+e**4/8.0d0*(1.5d1+2.8d1*co2+7.7d1*co4)
      dtheta = -dos*n*e*co1-cinco*n/dos*e**2*co2
      dtheta = dtheta+n/cuatro*e**3*(co1-1.3d1*co3)
      dtheta = dtheta+n/1.2d1*e**4*(1.1d1*co2-5.15d1*co4)
!
      co = dcos(dos*del)
      si = dsin(dos*del)
!
      ejean = 3.75d0*(mext/mtot)*(rtot/a)**3*ar3
      emacl = 1.25d0*(nu/dos+n)**2*rtot**3/(g*mtot)
      epola = ejean/dos+emacl
!
      f(1) =-tres*g*mext/a**3*ar3*ep*si
      f(2) = nu/dos+dtheta-gama*ejean*si/ep/dos
      f(3) = gama*(ejean*co-ep)
      f(4) = gama*(epola-epz)
      f(5) = uno

 205  format (2i4,1p6d13.5)
 206  format (1p4d13.5)
 210  format (1p10d15.5)
 211  format (1p10d20.10)
      return     
      end

!****************************************************************************
!****************************************************************************
      subroutine bs (y,deltat,step,ll,neqt)
      implicit real*8 (a-h,k-z)
      parameter (imax=5)
      real*8 y(imax),dydx(imax),vv(imax)
      integer neqt
!
      eps  = 10.0**(-ll)
      htry = step
      t    = 0.0
!
      do while (t.lt.deltat)
       if (t+htry.gt.deltat) htry = deltat - t
       call force (y,dydx)
       call bstep (y,dydx,neqt,t,htry,eps,hdid,hnext)
       htry = hnext
      end do
      step = hnext
!
      return
      end

!****************************************************************************
!****************************************************************************
      subroutine bstep (y,dydx,nv,x,htry,eps,hdid,hnext)
      implicit real*8 (a-h,o-z)
      parameter (nmax=5,kmaxx=8,safe1=.25,safe2=.7,redmax=1.e-5,redmin=.7,tiny=1.e-30,scalmx=.1)
      integer nseq(kmaxx+1)
      real*8 y(nmax),dydx(nmax),yscal(nmax),a(kmaxx+1),alf(kmaxx,kmaxx)
      real*8 err(kmaxx),yerr(nmax),ysav(nmax),yseq(nmax)
      logical first,reduct
      save a,alf,epsold,first,kmax,kopt,nseq,xnew
      data first/.true./,epsold/-1./
      data nseq /2,4,6,8,10,12,14,16,18/
!
      if (eps.ne.epsold) then
       hnext = -1.0d29
       xnew  = -1.0d29
       eps1  = safe1*eps
       a(1) = nseq(1) + 1
       do 11 k = 1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
 11    continue
       do 13 iq = 2,kmaxx
        do 12 k = 1,iq-1
         alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+1)))
 12     continue
 13    continue
       epsold = eps
       do 14 kopt = 2,kmaxx-1
        if (a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt)) goto 1
 14    continue
 1     kmax = kopt
      end if
      h = htry
      do 15 i = 1,nv
       ysav(i) = y(i)
 15   continue
      if (h.ne.hnext.or.x.ne.xnew) then
       first = .true.
       kopt = kmax
      end if
      reduct = .false.
 2    do 17 k = 1,kmax
       xnew = x + h
       if (xnew.eq.x) return
       call mmid (ysav,dydx,nv,x,h,nseq(k),yseq)
       xest = (h/nseq(k))**2
       call pzextr (k,xest,yseq,y,yerr,nv)
       if (k.ne.1) then
        errmax = tiny
        do 16 i = 1,nv
         errmax = max(errmax,abs(yerr(i)))
 16     continue
        errmax = errmax/eps
        km = k - 1
        err(km) = (errmax/safe1)**(1./(2*km+1))
       end if
       if (k.ne.1.and.(k.ge.kopt-1.or.first)) then
        if (errmax.lt.1.) goto 4
        if (k.eq.kmax.or.k.eq.kopt+1) then
         red = safe2/err(km)
         goto 3
        else if (k.eq.kopt) then
         if (alf(kopt-1,kopt).lt.err(km)) then
          red = 1./err(km)
          goto 3
         end if
        else if (kopt.eq.kmax) then
         if (alf(km,kmax-1).lt.err(km)) then
          red = alf(km,kmax-1)*safe2/err(km)
          goto 3
         end if
        else if (alf(km,kopt).lt.err(km)) then
         red = alf(km,kopt-1)/err(km)
         goto 3
        end if
       end if
 17   continue
 3    red = min(red,redmin)
      red = max(red,redmax)
      h = h*red
      reduct = .true.
      goto 2
 4    x = xnew
      hdid = h
      first = .false.
      wrkmin = 1.d35
      do 18 kk = 1,km
       fact = max(err(kk),scalmx)
       work = fact*a(kk+1)
       if (work.lt.wrkmin) then
        scale = fact
        wrkmin = work
        kopt = kk + 1
       end if
 18   continue
      hnext = h/scale
      if (kopt.ge.k.and.kopt.ne.kmax.and..not.reduct) then
       fact = max(scale/alf(kopt-1,kopt),scalmx)
       if (a(kopt+1)*fact.le.wrkmin) then
        hnext = h/fact
        kopt = kopt + 1
       end if
      end if
!
      return
      end

!****************************************************************************
!****************************************************************************
      subroutine mmid (y,dydx,nvar,xs,htot,nstep,yout)
      implicit real*8 (a-h,o-z)
      parameter (nmax=5)
      real*8 dydx(nmax),y(nmax),ym(nmax),yn(nmax),yout(nmax)
      real*8 derivs(nmax)
!
      h  = htot/nstep
      do 11 i = 1,nvar
       ym(i) = y(i)
       yn(i) = y(i) + h*dydx(i)
 11   continue
      x  = xs + h
      call force (yn,derivs)
      h2 = 2.*h
      do 13 n = 2,nstep
       do 12 i = 1,nvar
        swap = ym(i) + h2*derivs(i)
        ym(i) = yn(i)
        yn(i) = swap
 12    continue
       x = x + h
       call force (yn,derivs)
 13   continue
      do 14 i = 1,nvar
       yout(i) = 0.5*(ym(i)+yn(i)+h*derivs(i))
 14   continue
!     
      return
      end

!****************************************************************************
!****************************************************************************
      subroutine pzextr (iest,xest,yest,yz,dy,nv)
      implicit real*8 (a-h,o-z)
      parameter (immax=23,nmax=5)
      real*8 dy(nmax),yest(nmax),yz(nmax),d(nmax),qcol(nmax,immax)
      real*8 x(immax)
      save qcol,x
!     
      x(iest) = xest
      do 11 j = 1,nv
       dy(j) = yest(j)
       yz(j) = yest(j)
 11   continue
      if (iest.eq.1) then
       do 12 j = 1,nv
        qcol(j,1) = yest(j)
 12    continue
      else
       do 13 j = 1,nv
        d(j) = yest(j)
 13    continue
       do 15 k1 = 1,iest-1
        delta = 1.0/(x(iest-k1)-xest)
        f1 = xest*delta
        f2 = x(iest-k1)*delta
        do 14 j = 1,nv
         q = qcol(j,k1)
         qcol(j,k1) = dy(j)
         delta = d(j) - q
         dy(j) = f1*delta
         d(j)  = f2*delta
         yz(j) = yz(j) + dy(j)
 14     continue
 15    continue
       do 16 j = 1,nv
        qcol(j,iest) = dy(j)
 16    continue
      end if
!     
      return
      end

!****************************************************************************
!****************************************************************************
      subroutine energia (yc,w)
      implicit real*8 (a-h,k-z)
      real*8 yc(5),w(8)
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /cte/   g,k,gm,std,sty,uam,zmsolkg,gwatt,twatt
      common /ext/   mext,next,perext,msum,mred
      common /body/  mtot,rtot,ctot,dtot,grav,eta,gama,ej,em
      common /orb/   a,e,n,per
!
      nu   = yc(1)
      del  = yc(2)
      ep   = yc(3)
      epz  = yc(4)
      time = yc(5)
!
      co1 = dcos(n*time)
      si1 = dsin(n*time)
      co2 = co1*co1-si1*si1
      si2 = si1*co1+co1*si1
      co3 = co1*co2-si1*si2
      si3 = si1*co2+co1*si2
      co4 = co1*co3-si1*si3
      si4 = si1*co3+co1*si3
!
      ar3 = uno+tres*e*co1
      ar3 = ar3+tres/dos*e**2*(uno+tres*co2)
      ar3 = ar3+e**3/8.0d0*(2.7d1*co1+5.3d1*co3)
      ar3 = ar3+e**4/8.0d0*(1.5d1+2.8d1*co2+7.7d1*co4)
      v = n*time+dos*e*si1+cinco/cuatro*e**2*si2
      v = v-e**3/cuatro*(si1-1.3d1/tres*si3)
      v = v-e**4/2.4d1*(1.1d1*si2-2.575d1*si4)
!
      omega = nu/dos+n
      ejean = 3.75d0*(mext/mtot)*(rtot/a)**3*ar3
      emacl = 1.25d0*omega**2*rtot**3/(g*mtot)
      epola = ejean/dos+emacl
!
      co = dcos(dos*del)
      si = dsin(dos*del)
!
      E0   = g*mext*ctot/(dos*a**3)*ar3
      E1   = 0.08d0*g*mtot**2*gama/rtot
      Etot = tres/dos*E0*gama*((ep*co-ejean)+dos/tres*(epz-epola))
      Erot =-tres*E0*omega*ep*si
      Eorb = Etot-Erot
      Eint =-E1*(ep*(ep-ejean*co)+cuatro/tres*epz*(epz-epola))
      dEorb= E0*(ep*co+dos/tres*epz)*gama*(epz-epola)
      dErot=-ctot/tres*omega**2*gama*(epz-epola)
      dadt = dos*a**2/(g*mtot*mext)*Eorb
      dedt = dadt/dos/a-n*a/(g*mtot*mext*dsqrt(uno-e**2))*Erot/omega
      dedt = (uno-e**2)/e*dedt
!
      w(1) = twatt*Eorb              ! E_orb  [TW]
      w(2) = twatt*Erot              ! E_rot  [TW]
      w(3) = twatt*Etot              ! E_tot  [TW]
      w(4) = twatt*Eint              ! E_int  [TW]
      w(5) = twatt*dEorb             ! dEorb  [TW]
      w(6) = twatt*dErot             ! dErot  [TW]
      w(7) = dadt*std/sty*uam/1.0d3  ! da/dt  [km/yr]
      w(8) = dedt*std/sty            ! de/dt  [1/yr]

 203  format (2i4,1p6d13.5)
 204  format (1p4d13.5)
 202  format (1p8d13.5)
 210  format (1p8d20.10)
      return     
      end

!****************************************************************************
!****************************************************************************
      subroutine atractor (time,nu,del,er,ez)
      implicit real*8 (a-h,k-z)
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /ext/   mext,next,perext,msum,mred
      common /body/  mtot,rtot,ctot,dtot,grav,eta,gama,ej,em
      common /orb/   a,e,n,per
!
      ell  = n*time
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
      d0 = tres*p*e**2/(uno+p**2)*(dos+(tres*al-uno)*p**2)/div
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
      
      nu  = b0+b1*dcos(ell)+b2*dsin(ell)+b3*dcos(ell2)+b4*dsin(ell2)
      del = d0+d1*dcos(ell)+d2*dsin(ell)+d3*dcos(ell2)+d4*dsin(ell2)
      er  = e0+e1*dcos(ell)+e2*dsin(ell)+e3*dcos(ell2)+e4*dsin(ell2)
      ez  = z0+z1*dcos(ell)+z2*dsin(ell)+z3*dcos(ell2)+z4*dsin(ell2)

      return     
      end

!****************************************************************************
!****************************************************************************
      subroutine ctes
      implicit real*8 (a-h,k-z)
      common /num/   cero,uno,dos,tres,cuatro,cinco,twopi,pi,rad
      common /cte/   g,k,gm,std,sty,uam,zmsolkg,gwatt,twatt
      common /ext/   mext,next,perext,msum,mred
      common /body/  mtot,rtot,ctot,dtot,grav,eta,gama,ej,em
      common /orb/   a,e,n,per

!     numeros
      cero   = 0.0d0
      uno    = 1.0d0
      dos    = 2.0d0
      tres   = 3.0d0
      cuatro = 4.0d0
      cinco  = 5.0d0
      twopi  = 8.0d0*datan(1.0d0)
      pi     = twopi/dos
      rad    = pi/180.0d0

!     constantes del problema
      k   = 1.720209895d-02                    ! Gaussian constant
      g   = k**2                               ! Gravitational constant [AU^3/M_sun day^2]
      sty = uno/(60.0d0*60.0d0*24.0d0*365.0d0) ! Second to year
      std = uno/(60.0d0*60.0d0*24.0d0)         ! Second to day
      uam = 1.495978707d11                     ! AU in m
      zmsolkg = 1.98911d30                     ! solar mass in kg
      gwatt = zmsolkg*uam**2*std**3*1.0d-9
      twatt = zmsolkg*uam**2*std**3*1.0d-12

!     Particle basic data
      mext    = 5.68326d26/zmsolkg    ! Mass [M_sun]
      next    = 6.713428d-9/std       ! Mean motion [1/day]
      perext  = twopi/next            ! Orbital period [day]

!     Extended body basic data
      mtot = 1.08d20/zmsolkg          ! Mass [M_sun]
      rtot = 252.1d3/uam              ! Mean radius [AU]
      ctot = 0.4d0*mtot*rtot**2       ! Moment of inertial [M_sun AU^2]
      dtot = mtot/(cuatro*pi/tres*rtot**3)  ! Density [M_sun/AU^3]
      grav = g*mtot/rtot**2                 ! Gravity [AU/day^2]
      msum = mext+mtot                ! m1+m2
      gm   = g*msum
      mred = mext*mtot/msum

!     orbital elements
      a   = 238.02d6/uam           ! Semi-major axis [AU]
      e   = 0.0045d0               ! Eccentricity
      n   = dsqrt(gm/a**3)         ! Mean motion [1/day]
      per = twopi/n                ! Orbital period [day]

!     flattenings
      ej = 3.75d0*(mext/mtot)*(rtot/a)**3
      em = 1.25d0*n**2*rtot**3/(g*mtot)

 201  format (1p11d20.10)
 203  format (1p6d13.5)
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      Enceladus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      mtot = 1.08d20/zmsolkg          ! Mass [M_sun]
!      rtot = 252.1d3/uam              ! Mean radius [AU]
!      a   = 238.02d6/uam           ! Semi-major axis [AU]
!      e   = 0.0045d0               ! Eccentricity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      Mimas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      mtot = 0.379d20/zmsolkg               ! Mass [M_sun]
!      rtot = 197.49d3/uam                   ! Mean radius [AU]
!      a   = 185.52d6/uam           ! Semi-major axis [AU]
!      e   = 0.01986d0              ! Eccentricity
