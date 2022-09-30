      implicit real*8 (a-h,k-z)
      real*8 y(3)
      common /ind/ n
      open (1,file='density0.5.dat',status='replace')
      n = 0.5d0
c     lazo principal de integracion y precision de la integracion.
      time   = 1.0d-6
      tfinal = 1.0d2         ! final integration time (initial is = 0).
      htime  = 5.0d-3        ! output interval (in days)
      h      = 5.0d-3        ! initial timestep (in days)
      ll     = 13
c
      y(1) = 0.0d0    ! phi
      y(2) = 1.0d0    ! theta
      y(3) = time     ! xi
      do 50 while (time.lt.tfinal)
c     integrate using Bulirsch-Stoer (precision given by ll).
       call bus (y,htime,h,ll,3)
       time = time + htime
       if(y(2).lt.0.0d0.or.y(2).ne.y(2)) goto 51
 50   continue
 51   continue
      rtime  = time
      time   = 1.0d-6
      tfinal = rtime                 ! final integration time (initial is = 0).
      htime  = rtime/2000.0d0        ! output interval (in days)
      h      = rtime/2000.0d0        ! initial timestep (in days)
      ll     = 13
c
      y(1) = 0.0d0    ! phi
      y(2) = 1.0d0    ! theta
      y(3) = time     ! xi
      write(1,102) time/rtime,y(2)**n
      do 60 while (time.lt.tfinal)
c     integrate using Bulirsch-Stoer (precision given by ll).
       call bus (y,htime,h,ll,3)
       time = time + htime
       if(y(2).gt.0.0d0) write(1,102) time/rtime,y(2)**n
 60   continue
 61   continue


 102  format (1p2e20.10)
 103  format (1p3e20.10)
 104  format (1p4e20.10)
 105  format (1p5e20.10)
c     
      end

      subroutine force (y,f)
      implicit real*8 (a-h,k-z)
      real*8 y(3),f(3)
      common /ind/ n

      u = y(1)
      v = y(2)
      t = y(3)


c     equations of motion.
      f(1) = v**n*t**2
      f(2) = -u/t**2
      f(3) = 1.0d0
!      f(4) = 1.0d0

 101  format (1p4e20.10)
      return     
      end


      subroutine bus (y,deltat,step,ll,neqt)
      implicit real*8 (a-h,k-z)
      parameter (imax=3)
      real*8 y(imax),dydx(imax)
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
      parameter (imax=3)
      parameter (nmax=4,kmaxx=8,safe1=.25,safe2=.7,
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
      parameter (imax=3)
      parameter (nmax=4)
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
      parameter (imax=3)
      parameter (immax=23,nmax=4)
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


