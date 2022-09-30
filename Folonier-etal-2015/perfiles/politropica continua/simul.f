      implicit real*8(a-h,k-z)
      open (9,file='perf2.dat',status='replace')

      in   = 100
      nmin = 0.0d0
      nmax = 5.0d0
      dn   = (nmax-nmin)/dfloat(in)
      do i=0,in
       zn = nmin+dn*dfloat(i)
       call integracion (zn)
       call layer (zh0,zkf,zinerc)
       write(*,302) zn,zh0,0.4d0*(zkf+1.0d0),zkf,zinerc
       write(9,302) zn,zh0,0.4d0*(zkf+1.0d0),zkf,zinerc
      end do

      close(9)
 302  format (1p5e20.10)
      end

      subroutine integracion (zn)
      implicit real*8 (a-h,k-z)
      real*8 y(3),nindex
      common /ind/ nindex
      open (8,file='density.in',status='replace')
c     lazo principal de integracion y precision de la integracion.
      nindex = zn
      time   = 1.0d-6
      tfinal = 1.0d10         ! final integration time (initial is = 0).
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
      htime  = rtime/1000.0d0        ! output interval (in days)
      h      = rtime/1000.0d0        ! initial timestep (in days)
      ll     = 11
c
      y(1) = 0.0d0    ! phi
      y(2) = 1.0d0    ! theta
      y(3) = time     ! xi

      do 60 while (time.lt.tfinal)
c     integrate using Bulirsch-Stoer (precision given by ll).
       call bus (y,htime,h,ll,3)
       time = time + htime
       if(y(2).gt.0.0d0) write(8,102) time/rtime,y(2)**nindex
 60   continue
 61   continue
      close(8)
 102  format (1p2e20.10)
      return
      end

      subroutine force (y,f)
      implicit real*8 (a-h,k-z)
      real*8 y(3),f(3),nindex
      common /ind/ nindex

      u = y(1)
      v = y(2)
      t = y(3)


c     equations of motion.
      f(1) = v**nindex*t**2
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

      subroutine layer (zzh0,kf,zzctot)
      implicit real*8(a-h,l-z)
      parameter (imax=2010)
      real*8 r(0:imax+1),rho(0:imax+1),mlay(imax)
      real*8 sigma(imax),m(imax),mtot,h(0:imax)
      real*8 E(imax,imax),F(imax,imax),D(imax)
      real*8 kf,a(imax),b(imax),c(imax),f1(imax),f2(imax),f3(imax)
      integer i,j,k,ii,iin,ifull,iscreen
      double precision aaa(imax,imax), ccc(imax,imax)
      double precision LLL(imax,imax), UUU(imax,imax)
      double precision coeff,bbb(imax),ddd(imax),xxx(imax)
      common /cte/ g,gkg,mekg,mjkg,mskg,mkm,cmkm,aukm,ds
      common /num/ cero,uno,dos,tres,cinco,pi,twopi
      common /den/ r,rho,rs,iin
      common /out/ ifull,iscreen
      common /tyr/ mpar,rpar,w

      open (1,file='layers.in',status='old')
      call data
     
!     we change layers for ellipsoids
      m = 0.0d0
      mtot = 0.0d0
      do i=1,iin
       sigma(i) = rho(i)-rho(i+1)                                ! density of the ellipsoids
       m(i)     = 4.0d0*pi/3.0d0*sigma(i)*r(i)**3                ! mass of the ellipsoids
       mlay(i)  = 4.0d0*pi/3.0d0*rho(i)*(r(i)**3-r(i-1)**3)      ! mass of the layers
       mtot=mtot+m(i)                                            ! mass of the body
       D(i) = r(i)**3/r(iin)**3                                    ! auxiliary vector
      end do
      mrho = tres*mtot*1.0d-12/(4.0d0*pi*rs**3)                  ! mean density
!     calculation of the equivalent homogeneous solutions
      ej = 3.75d0*(mpar/mtot)*(rs/rpar)**3
      em = 1.25d0*rs**3*w**2/(mtot*g)
      er = ej+em

!     calculation of the matrix E
      E = 0.0d0
      do i=1,iin
       do j=1,iin
        if(i.eq.j) then
         E(i,j)=0.0d0
         do k=i,iin
          E(i,j)=E(i,j)+5.0d0*m(k)/(2.0d0*mtot)*(1.0d0-(r(i)/r(k))**3)
         end do
         E(i,j)=+1.0d0+3.0d0*(mtot-m(i))/(2.0d0*mtot)-E(i,j)
        end if

        if(i.gt.j) then
         E(i,j) = -3.0d0*m(j)/(2.0d0*mtot)*(r(j)/r(i))**2
        end if

        if(i.lt.j) then
         E(i,j) = -3.0d0*m(j)/(2.0d0*mtot)*(r(i)/r(j))**3
        end if
       end do
      end do


!     calculation of the matrix F

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
      do i=1,iin
       do j=1,iin
        aaa(i,j) = E(i,j)
       end do
      end do

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
      LLL=0.0
      UUU=0.0
      bbb=0.0

! step 1: forward elimination
      do k=1, iin-1
       do i=k+1,iin
        coeff=aaa(i,k)/aaa(k,k)
        LLL(i,k) = coeff
        do j=k+1,iin
         aaa(i,j) = aaa(i,j)-coeff*aaa(k,j)
        end do
       end do
      end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
      do i=1,iin
       LLL(i,i) = 1.0
      end do
! U matrix is the upper triangular part of A
      do j=1,iin
       do i=1,j
        UUU(i,j) = aaa(i,j)
       end do
      end do

! Step 3: compute columns of the inverse matrix C
      do k=1,iin
       bbb(k)=1.0
       ddd(1) = bbb(1)
! Step 3a: Solve Ld=b using the forward substitution
       do i=2,iin
        ddd(i)=bbb(i)
        do j=1,i-1
         ddd(i) = ddd(i) - LLL(i,j)*ddd(j)
        end do
       end do
! Step 3b: Solve Ux=d using the back substitution
       xxx(iin)=ddd(iin)/UUU(iin,iin)
       do i = iin-1,1,-1
        xxx(i) = ddd(i)
        do j=iin,i+1,-1
         xxx(i)=xxx(i)-UUU(i,j)*xxx(j)
        end do
        xxx(i) = xxx(i)/uuu(i,i)
       end do
! Step 3c: fill the solutions x(n) into column k of C
       do i=1,iin
        ccc(i,k) = xxx(i)
       end do
       bbb(k)=0.0
      end do 

      do i=1,iin
       do j=1,iin
        F(i,j) = ccc(i,j)
       end do
      end do
!============================================================
!============================================================

!     calculation of the flattening and semiaxes 
      zctot = cero           ! C ellipsoid
      do i=1,iin
       h(i) = 0.0d0
       do k=1,iin
        h(i)=h(i)+F(i,k)*D(k)
       end do
       ttq = 4.0d0*pi*(rho(i)-rho(i+1))*r(i)**5/15.0d0
       zctot = zctot + ttq
      end do

      zzctot = zctot/(mtot*rs**2)
!     calculation of the equivalent fluid Love number, C20 and C22
      h(0)=cero
      kf=0.0d0
      do i=1,iin
       kf=kf+m(i)*h(i)*r(i)**2
      end do
      kf = 1.5d0*kf/(mtot*r(iin)**2)            ! numero de Love fluido equivalente del cuerpo
      zzh0 = h(1)

      close (1)
 402  format (1p7e15.5)
      return
      end

      subroutine data
      implicit real*8(a-h,l-z)
      parameter (imax=2010)
      real*8 r(0:imax+1),rho(0:imax+1)
      integer i,iflag,ii,iin,ifull,iscreen
      character(len=5)  zout,askfull,askonscreen
      character(len=45) shortline
      character(len=50) line,archin,archout
      common /cte/ g,gkg,mekg,mjkg,mskg,mkm,cmkm,aukm,ds
      common /num/ cero,uno,dos,tres,cinco,pi,twopi
      common /den/ r,rho,rs,iin
      common /out/ ifull,iscreen
      common /tyr/ mpar,rpar,w

!     some constants
      g    = 6.693d-20        ! gravitational constant [km^3/kg*s^2]
      gkg  = 1.0d-3           ! gram in kilogram
      mekg = 5.952d24         ! mass of Earth in kilogram
      mjkg = 1.899d27         ! mass of Jupiter in kilogram
      mskg = 1.98911d30       ! mass of Sun in kilogram
      mkm  = 1.0d-3           ! meter in kilometer
      cmkm = 1.0d-5           ! centimeter in kilometer
      aukm = 1.495978707d8    ! AU in kilometer
      ds   = 8.64d4           ! day in second

!     numbers
      cero  = 0.0d0
      uno   = 1.0d0
      dos   = 2.0d0
      tres  = 3.0d0
      cinco = 5.0d0
      twopi = 8.0d0*datan(1.0d0)
      pi    = twopi/2.0d0

      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a50,d20.10)') line,mpar
      read (1,'(a50,d20.10)') line,rpar
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a50,d20.10)') line,w
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a)') dummy

      read (1,'(a50,a50)') line,archin
      read (1,'(a50,i6)') line,iin
      iflag=0
      r   = cero
      rho = cero
      if (iin.eq.0) then
       open (2,file=archin,status='old')
       ii = 1
 3     continue
       read (2,*,err=4,end=4) r(ii),rho(ii)
       rho(ii) = rho(ii)*1.0d12
       ii = ii+1
       goto 3
 4     iin = ii-1
       close (2)
       do while (iflag.eq.0)
        read (1,'(a5)') zout
        if(zout.eq.'-----') goto 5
       end do
      else
       do i=1,iin
        read (1,'(a50,d20.10)') line,r(i)
       end do
       do i=1,iin
        read (1,'(a50,d20.10)') line,rho(i)
        rho(i) = rho(i)*1.0d12
       end do
       do while (iflag.eq.0)
        read (1,'(a5)') zout
        if(zout.eq.'-----') goto 5
       end do
      end if
 5    continue
      rs = r(iin)
      r(0) = cero
      rho(iin+1) = cero

      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a50,a50)') line,archout
      read (1,'(a50,a5)') line,askfull
      read (1,'(a50,a5)') line,askonscreen

!      open (3,file=archout,status='replace')
      ifull=0
      if (askfull.eq.'yes') ifull=1
      iscreen=0
      if (askonscreen.eq.'yes') iscreen=1

      return
      end

