      implicit real*8(a-h,k-z)
      real*8 y(2),f(2)
!      integer a
      common /pot/ a
      open (1,file='2alpha.dat',status='replace')

!       amin = 0.1d0
!       amax = 1.0d3
!       da   = (amax-amin)/5000.0
!       do i=0,5000
      amin = -4.0d0
      amax = 4.0d0
      da   = (amax-amin)/100.0
      do i=0,100
       a=10.0d0**(amin+dfloat(i)*da)
c     lazo principal de integracion y precision de la integracion.
       time   = 0.0d0
!       time   = 1.0d-6
       tfinal = 1.0d0         ! final integration time (initial is = 0).
       htime  = 1.0d-4         ! output interval (in days)
       h      = 1.0d-4         ! initial timestep (in days)
       ll     = 13
c
       y(1) = 0.0d0    ! phi
       y(2) = 0.0d0    ! phi

       do 50 while (time.lt.tfinal)
        call bus (y,htime,h,ll,2)
        time = time + htime
 50    continue
 51    continue
       zhs  = 2.0/(2.0+y(1))
       zdhs = 2.0*y(1)/(2.0+y(1))
       zkf  = 5.0d0*zhs/2.0d0-1.0d0
       zcc  = 0.4d0*(3.0d0+a)/(5.0d0+a)

       write(1,104) dlog10(a),zhs,zkf,zcc
       write(*,104) dlog10(a),zhs,zkf,zcc
      end do

      close(1)
 102  format (1p2e20.10)
 103  format (1p3e20.10)
 104  format (1p4e20.10)
 105  format (1p5e20.10)
      end

      subroutine force (y,f)
      implicit real*8 (a-h,k-z)
      real*8 y(2),f(2)
!      integer a
      common /pot/ a

      u = y(1)
      t = y(2)
!      b=i(a)
!      b=a


c     equations of motion.
      q=-6.0*a*t**(a-1.0d0)/(3.0+a-3.0*t**a)
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

      if(t.lt.0.00001) then
       f(1) = 5.
      else
!       f(1) = -u**2/t**6-q*u-q*t**5
       f(1) = -u**2/t-(q+5.0d0/t)*u-q
      end if
      f(2) = 1.0d0

 101  format (1p4e20.10)
      return     
      end


      subroutine bus (y,deltat,step,ll,neqt)
      implicit real*8 (a-h,k-z)
      parameter (imax=2)
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
      parameter (imax=2)
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
      parameter (imax=2)
      parameter (nmax=2)
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
      parameter (imax=2)
      parameter (immax=23,nmax=2)
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
