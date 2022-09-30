      implicit real*8(a-h,k-z)
      real*8 y(3),f(3)
      integer i,ih,iflag
      common /con/ iflag
!     2.0000000000E-01    5.8020377066E-01    4.5050942665E-01    2.4615384615E-01
!     1.0000000000E-01    5.7501982215E-01    4.3754955537E-01    2.4313725490E-01
!     1.0000000000E+00    6.1850156379E-01    7.6299687242E-01    1.2336215736E+00
!     2.0000000000E+00    6.5887182566E-01    6.8225634869E-01    1.0354917635E+00
!     5.0000000000E+00    7.4328638326E-01    5.1342723349E-01    6.9075291174E-01
!     1.0000000000E+01    8.1944118496E-01    3.6111763009E-01    4.4068767438E-01
!     1.0000000000E+02    9.7183483715E-01    5.6330325708E-02    5.7962859073E-02
      iflag = -1
      zhs  = 0.57501982215
      ih   = 10000
      hmin = 0.0d0
      hmax = 1.0d0
      dh   = (hmax-hmin)/dfloat(ih)
      y1   = 0.0d0
      y2   = 0.0d0
      do i=2,ih
       zh0 = hmin+dfloat(i)*dh
       call integra (zh0,zh1)
!        write(*,*) zh0,zh1,zhs
!        read(*,*)
       if(zh1.gt.zhs) then
        x1 = zh0-dh
        x2 = zh0
        y2 = zh1
c
        za = (y2-y1)/dh
        zb = y2-za*x2
c
        xc = (zhs-zb)/za
        if(zh1.gt.zhs) goto 49
       end if
       y1 = zh1
      end do
 49   continue
      iflag = 1
      call integra (xc,caca)
      end


      subroutine integra (zz0,zzs)
      implicit real*8 (a-h,k-z)
      real*8 y(3),f(3)
      integer iflag,poten
      common /con/ iflag

c     lazo principal de integracion y precision de la integracion.
      time   = 1.0d-6
!      time   = 1.0d-30
!      time   = 0.0d0
      tfinal = 1.0d0          ! final integration time (initial is = 0).
      htime  = 1.0d-3         ! output interval (in days)
      h      = 1.0d-3         ! initial timestep (in days)
      ll     = 13
c
      y(1) = zz0
      y(2) = 0.0d0
      y(3) = 0.0d0
      poten = 100
      densi = 1.0-time**poten
!      densi = 1.0-time**0.1
!      write(*,*) densi
!      stop
      if (dfloat(iflag).gt.0.0) then
       open (1,file='alpha100-k2.dat',status='replace')
       write(1,103) time,densi,y(1)
       write(*,103) time,densi,y(1)
      end if

      do 50 while (time.lt.tfinal)
       call bus (y,htime,h,ll,3)
!       write(*,*) 'y',time,y
!       read(*,*) 
       time = time + htime
         if (dfloat(iflag).gt.0.0) then
          densi = 1.0-time**poten
!          densi = 1.0-time**0.1
          write(1,103) time,densi,y(1)
          write(*,103) time,densi,y(1)
        end if
 50   continue
 51   continue
      zzs = y(1)

      if (dfloat(iflag).gt.0.0) then
      close(1)
      end if
 102  format (1p2e20.10)
 103  format (1p3e20.10)
 104  format (1p4e20.10)
      return     
      end

      subroutine force (y,f)
      implicit real*8 (a-h,k-z)
      real*8 y(3),f(3)
      integer a

      u = y(1)
      v = y(2)
      t = y(3)

c     equations of motion.
      a = 100
      b = dfloat(a)
      q = -6.0*b*t**(a-1)/(3.0+b-3.0*t**a)
!      q = -6.0*0.1*t**(0.1-1.0)/(3.0+0.1-3.0*t**0.1)
!      q=-6.0/(4.0-3.0*t**2)
!      q=-12.0d0*t/(5.0-3.0*t**2)

!        pi = 4.0d0*atan(1.0d0)
!       xt=pi*t
!       if(t.eq.0) then
!        q=0.
!       else
!        q=-6./t+2.*pi*xt*dsin(xt)/(dsin(xt)-xt*dcos(xt))
!       end if

!       pi = 4.0d0*atan(1.0d0)
!       xt=pi*t/2.0
!       xx=pi*t
!       if(t.eq.0) then
!        q=0.
!       else
!        q=-6./t+pi*xx**2/(4.*xx+dtan(xt)*(xx**2-8.))
!       end if

!      if(t.lt.0.0000001) then
      if(t.eq.0.0) then
       f(1) = 0.0
       f(2) = 0.0
      else
       f(1) = u*v/t**6
       f(2) = -v**2/t**6-q*v-q*t**5
      end if
      f(3) = 1.0d0
 !     write(*,*) 'f',f(1),f(2)
 !     read(*,*)

 101  format (1p4e20.10)
      return     
      end


      subroutine bus (y,deltat,step,ll,neqt)
      implicit real*8 (a-h,k-z)
      parameter (nmax=3)
      real*8 y(nmax),dydx(nmax)
c
      eps  = 10.0**(-ll)
      htry = step
      t    = 0.0
c
      do while (t.lt.deltat)
       if (t+htry.gt.deltat) htry = deltat - t
       call force (y,dydx)
       call bstep (y,dydx,neqt,t,htry,eps,hdid,hnext)
       htry = hnext
      end do
      step = hnext
c
      return
      end


      subroutine bstep (y,dydx,nv,x,htry,eps,hdid,hnext)
      implicit real*8 (a-h,o-z)
      parameter (nmax=3,kmaxx=8,safe1=.25,safe2=.7,
     *     redmax=1.e-5,redmin=.7,tiny=1.e-30,scalmx=.1)
      integer nseq(kmaxx+1)
      real*8 y(nmax),dydx(nmax),yscal(nmax),a(kmaxx+1),alf(kmaxx,kmaxx)
      real*8 err(kmaxx),yerr(nmax),ysav(nmax),yseq(nmax)
      logical first,reduct
      save a,alf,epsold,first,kmax,kopt,nseq,xnew
      data first/.true./,epsold/-1./
      data nseq /2,4,6,8,10,12,14,16,18/
c
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
c
      return
      end


      subroutine mmid (y,dydx,nvar,xs,htot,nstep,yout)
      implicit real*8 (a-h,o-z)
      parameter (nmax=3)
      real*8 dydx(nmax),y(nmax),ym(nmax),yn(nmax),yout(nmax)
      real*8 derivs(nmax)
c
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
c     
      return
      end

      subroutine pzextr (iest,xest,yest,yz,dy,nv)
      implicit real*8 (a-h,o-z)
      parameter (immax=23,nmax=3)
      real*8 dy(nmax),yest(nmax),yz(nmax),d(nmax),qcol(nmax,immax)
      real*8 x(immax)
      save qcol,x
c     
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
c     
      return
      end
