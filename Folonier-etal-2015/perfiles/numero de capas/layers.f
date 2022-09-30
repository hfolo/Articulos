      program layers
      implicit real*8(a-h,l-z)
      parameter (imax=3000)
      real*8 r(0:imax+1),rho(0:imax+1),mlay(imax)
      real*8 sigma(imax),m(imax),mtot,h(0:imax),xl(0:imax)
      real*8 za(imax),zb(imax),zc(imax)
      real*8 zla(imax),zlb(imax),zlc(imax)
      real*8 E(imax,imax),F(imax,imax),D(imax)
      real*8 f1(imax),f2(imax),f3(imax),a(0:imax),b(0:imax),c(0:imax)
      real*8 kf,c20,c22,kl(imax)
      character(len=7) ind,ind1
      character(len=15) l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,q7c,q7d
      character(len=50) q1,q2,q3,q4,q5,q6,q7,q8,q9,q10
      character(len=50) q11,q12,q13,q14,q15,q16,q7b
      character(len=15) coe1,coe2,coe3,coe4,sh1,sh2,sh3
      character(len=50) dline
      integer i,j,k,ii,n,ifull,iscreen
      double precision aaa(imax,imax), ccc(imax,imax)
      double precision LLL(imax,imax), UUU(imax,imax)
      double precision coeff,bbb(imax),ddd(imax),xxx(imax)
      common /cte/ g,gkg,mekg,mjkg,mskg,mkm,cmkm,aukm,ds
      common /num/ cero,uno,dos,tres,cinco,pi,twopi
      common /den/ r,rho,rs,n
      common /out/ ifull,iscreen
      common /tyr/ mpar,rpar,w
      common /car1/ coe1,coe2,coe3,coe4,sh1,sh2,sh3,dline,ind,ind1
      common /car2/ l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,q16,q7b,q7c,q7d
      common /car3/ q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15

      open (1,file='layers.in',status='old')
      open (8,file='simple3000.dat',status='replace')
      call data
      call caract
     
!     we change layers for ellipsoids
      m = 0.0d0
      mtot = 0.0d0
      do i=1,n
       sigma(i) = rho(i)-rho(i+1)                                ! density of the ellipsoids
       m(i)     = 4.0d0*pi/3.0d0*sigma(i)*r(i)**3                ! mass of the ellipsoids
       mlay(i)  = 4.0d0*pi/3.0d0*rho(i)*(r(i)**3-r(i-1)**3)      ! mass of the layers
       mtot=mtot+m(i)                                            ! mass of the body
       D(i) = r(i)**3/r(n)**3                                    ! auxiliary vector
      end do
      mrho = tres*mtot*1.0d-12/(4.0d0*pi*rs**3)                  ! mean density
!     calculation of the equivalent homogeneous solutions
      ej = 3.75d0*(mpar/mtot)*(rs/rpar)**3
      em = 1.25d0*rs**3*w**2/(mtot*g)
      er = ej+em

!     calculation of the matrix E
      E = 0.0d0
      do i=1,n
       do j=1,n
        if(i.eq.j) then
         E(i,j)=0.0d0
         do k=i,n
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
!============================================================
!============================================================

!     calculation of the flattening and semiaxes 
      a(0) = cero
      b(0) = cero
      c(0) = cero

      do i=1,n
       h(i) = 0.0d0
       do k=1,n
        h(i)=h(i)+F(i,k)*D(k)
       end do
       f1(i) = h(i)*(ej+em)
       f2(i) = h(i)*em
       f3(i) = h(i)*ej
       cc1 = (uno-f1(i))**(1.0/3.0)
       cc2 = (uno-f2(i))**(1.0/3.0)
       a(i) = r(i)*cc2/(cc1**2)
       b(i) = r(i)*cc1/(cc2**2)
       c(i) = r(i)*cc1*cc2
       xl(i) = (h(i)*r(i)**5-h(i-1)*r(i-1)**5)/(r(i)**5-r(i-1)**5)
      end do


!     calculation of the moments of inertia
c
      zla = cero             ! A_i layer
      zlb = cero             ! B_i layer
      zlc = cero             ! C_i layer
c
      zlatot = cero          ! A layer
      zlbtot = cero          ! B layer
      zlctot = cero          ! C layer
c
      za = cero              ! A_i ellipsoid
      zb = cero              ! B_i ellipsoid
      zc = cero              ! C_i ellipsoid
c
      zatot = cero           ! A ellipsoid
      zbtot = cero           ! B ellipsoid
      zctot = cero           ! C ellipsoid
c
      ztot = cero            ! I
c
      do i=1,n
       ttt = 4.0d0*pi*rho(i)/15.0d0
       zaa = r(i)**3*(b(i)**2+c(i)**2)-r(i-1)**3*(b(i-1)**2+c(i-1)**2)
       zbb = r(i)**3*(c(i)**2+a(i)**2)-r(i-1)**3*(c(i-1)**2+a(i-1)**2)
       zcc = r(i)**3*(a(i)**2+b(i)**2)-r(i-1)**3*(a(i-1)**2+b(i-1)**2)

       zla(i) = ttt*zaa
       zlb(i) = ttt*zbb
       zlc(i) = ttt*zcc

       zlatot = zlatot + zla(i)
       zlbtot = zlbtot + zlb(i)
       zlctot = zlctot + zlc(i)

       za(i) = uno*m(i)*(c(i)**2+b(i)**2)/cinco
       zb(i) = uno*m(i)*(c(i)**2+a(i)**2)/cinco
       zc(i) = uno*m(i)*(a(i)**2+b(i)**2)/cinco

       zatot = zatot + za(i)
       zbtot = zbtot + zb(i)
       zctot = zctot + zc(i)
      end do

      zzlatot = zlatot/(mtot*rs**2)
      zzlbtot = zlbtot/(mtot*rs**2)
      zzlctot = zlctot/(mtot*rs**2)

      zzatot = zatot/(mtot*rs**2)
      zzbtot = zbtot/(mtot*rs**2)
      zzctot = zctot/(mtot*rs**2)

      ztot = (zzlatot+zzlbtot+zzlctot)/2.0d0

!     calculation of the equivalent fluid Love number, C20 and C22
      h(0)=cero
      kf=0.0d0
      do i=1,n
       kf=kf+m(i)*h(i)*r(i)**2
c
       kl(i) = 3.75d0*zlc(i)*xl(i)*r(n)**3/(mtot*r(i)**5)                       ! numero de Love fluido de la capa i-esima
      end do
      kf = 1.5d0*kf/(mtot*r(n)**2)            ! numero de Love fluido equivalente del cuerpo
      c20 = 2.0d0*kf*(ej+2.0d0*em)/15.0d0
      c22 = kf*ej/15.0d0

!     output
      if (iscreen.eq.1) then
      write(*,*) 'GLOBAL PARAMETERS'
      write(*,*) 
      write(*,'(a50,1p1e11.3)') q1,mtot
      write(*,'(a50,1p1e11.3)') q2,rs
      write(*,'(a50,1p1e11.3)') q3,mrho
      write(*,'(a50,1p1e11.3)') q4,zzctot
      write(*,'(a50,1p1e11.3)') q5,c20
      write(*,'(a50,1p1e11.3)') q6,c22
      write(*,'(a50,1p1e11.3)') q7,kf
      write(*,'(a50,1p1e11.3)') q7b,h(n)
      write(*,*) 
      write(*,*) 'PARAMETERS ON SURFACE'
      write(*,*) 
      write(*,'(a50,1p1e11.3)') q8,a(n)-rs
      write(*,'(a50,1p1e11.3)') q9,b(n)-rs
      write(*,'(a50,1p1e11.3)') q10,c(n)-rs
      write(*,*) 
      write(*,'(a50,1p1e11.3)') q11,f1(n)
      write(*,'(a50,1p1e11.3)') q12,f2(n)
      write(*,'(a50,1p1e11.3)') q13,f3(n)
      write(*,*) 
      write(*,*) 'EQUIVALENT HOMOGENEOUS SOLUTIONS'
      write(*,*) 
      write(*,'(a50,1p1e11.3)') q14,ej
      write(*,'(a50,1p1e11.3)') q15,em
      write(*,'(a50,1p1e11.3)') q16,er
      write(*,*) 
      write(*,'(3a49)') dline,dline,dline
      write(*,*) 
      write(*,'(a8,7a13)') ind,l1,l2,l3,l4,q7c,q7d
       do i=1,n
        write(*,104) i,r(i),rho(i),mlay(i),zlc(i),h(i),kl(i)
       end do
      end if

      write(3,*) 'GLOBAL PARAMETERS'
      write(3,*) 
      write(3,'(a50,1p1e18.8)') q1,mtot
      write(3,'(a50,1p1e18.8)') q2,rs
      write(3,'(a50,1p1e18.8)') q3,mrho
      write(3,'(a50,1p1e18.8)') q4,zctot/(mtot*rs**2)
      write(3,'(a50,1p1e18.8)') q5,c20
      write(3,'(a50,1p1e18.8)') q6,c22
      write(3,'(a50,1p1e18.8)') q7,kf
      write(3,'(a50,1p1e18.8)') q7b,h(n)
      write(3,*) 
      write(3,*) 'PARAMETERS ON SURFACE'
      write(3,*) 
      write(3,'(a50,1p1e18.8)') q8,a(n)-rs
      write(3,'(a50,1p1e18.8)') q9,b(n)-rs
      write(3,'(a50,1p1e18.8)') q10,c(n)-rs
      write(3,*) 
      write(3,'(a50,1p1e18.8)') q11,f1(n)
      write(3,'(a50,1p1e18.8)') q12,f2(n)
      write(3,'(a50,1p1e18.8)') q13,f3(n)
      write(3,*) 
      write(3,*) 'EQUIVALENT HOMOGENEOUS SOLUTIONS'
      write(3,*) 
      write(3,'(a50,1p1e18.8)') q14,ej
      write(3,'(a50,1p1e18.8)') q15,em
      write(3,'(a50,1p1e18.8)') q16,er
      if (ifull.eq.1) then
       write(3,*) 
       write(3,'(3a49)') dline,dline,dline
       write(3,*) 
       write(3,'(a5,7a18)') ind1,l1,l2,l3,l4,q7c,q7d
       do i=1,n
        write(3,103) i,r(i),rho(i),mlay(i),zlc(i),h(i),kl(i)
        write(8,108) r(i)/rs,h(i),kl(i),rho(i)/1.0d12,mlay(i)/mtot,
     *  zlc(i)/(mlay(i)*r(i)**2),zlc(i)/(mtot*rs**2)
       end do
      end if

      close (1)
      close (3)
      close (8)
 100  format(9A)
 101  format(1p3e15.5)
 102  format(1p4e15.5)
 103  format(i4,1p7e18.8)
 104  format(i4,1p7e13.3)
 108  format(1p7e20.10)
 109  format(1p2e20.10)
      end program

      subroutine data
      implicit real*8(a-h,l-z)
      parameter (imax=3000)
      real*8 r(0:imax+1),rho(0:imax+1)
      integer i,iflag,ii,n,ifull,iscreen
      character(len=5)  zout,askfull,askonscreen
      character(len=45) shortline
      character(len=50) line,archin,archout
      common /cte/ g,gkg,mekg,mjkg,mskg,mkm,cmkm,aukm,ds
      common /num/ cero,uno,dos,tres,cinco,pi,twopi
      common /den/ r,rho,rs,n
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
      read (1,'(a50,i6)') line,n
      iflag=0
      r   = cero
      rho = cero
      if (n.eq.0) then
       open (2,file=archin,status='old')
       ii = 1
 3     continue
       read (2,*,err=4,end=4) r(ii),rho(ii)
       rho(ii) = rho(ii)*1.0d12
       ii = ii+1
       goto 3
 4     n = ii-1
       close (2)
       do while (iflag.eq.0)
        read (1,'(a5)') zout
        if(zout.eq.'-----') goto 5
       end do
      else
       do i=1,n
        read (1,'(a50,d20.10)') line,r(i)
       end do
       do i=1,n
        read (1,'(a50,d20.10)') line,rho(i)
        rho(i) = rho(i)*1.0d12
       end do
       do while (iflag.eq.0)
        read (1,'(a5)') zout
        if(zout.eq.'-----') goto 5
       end do
      end if
 5    continue
      rs = r(n)
      r(0) = cero
      rho(n+1) = cero

      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a50,a50)') line,archout
      read (1,'(a50,a5)') line,askfull
      read (1,'(a50,a5)') line,askonscreen

      open (3,file=archout,status='replace')
      ifull=0
      if (askfull.eq.'yes') ifull=1
      iscreen=0
      if (askonscreen.eq.'yes') iscreen=1

      return
      end

      subroutine caract
      character(len=7) ind,ind1
      character(len=15) l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,q7c,q7d
      character(len=50) q1,q2,q3,q4,q5,q6,q7,q8,q9,q10
      character(len=50) q11,q12,q13,q14,q15,q16,q7b
      character(len=15) coe1,coe2,coe3,coe4,sh1,sh2,sh3
      character(len=50) dline
      common /car1/ coe1,coe2,coe3,coe4,sh1,sh2,sh3,dline,ind,ind1
      common /car2/ l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,q16,q7b,q7c,q7d
      common /car3/ q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15

      coe1 = '    C20'
      coe2 = '    C22'
      coe3 = '    kf'
      coe4 = '    Ic'
      sh1 = '    f_Jeans'
      sh2 = '    f_MacLaurin'
      sh3 = '    f_Roche'
      ind  = '  i'
      ind1 = '   i'
      dline = '=================================================='
      l1 = 'R(i)'
      l2 = 'rho(i)'
      l3 = 'mlay(i)'
      l4 = 'Clay(i)'
      l5 = 'a(i)'
      l6 = 'b(i)'
      l7 = 'c(i)'
      l8 = 'f1=(a-c)/a'
      l9 = 'f2=(b-c)/b'
      l10= 'f3=(a-b)/a'
      q7c= 'H(i)                                            :'
      q7d= 'kl(i)                                           :'
      q1 = ' m [kg]                                          :'
      q2 = ' R [km]                                          :'
      q3 = ' Mean Density [g/cm^3]                           :'
      q4 = ' C/mR^2 (Total Moment of Inertia Respect to z)   :'
      q5 = ' C20                                             :'
      q6 = ' C22                                             :'
      q7 = ' ekf (Equivalent Fluid Love Number)              :'
      q7b= ' H(n)                                            :'
      q8 = ' a-R [km]                                        :'
      q9 = ' b-R [km]                                        :'
      q10= ' c-R [km]                                        :'
      q11= ' f1=(a-c)/a                                      :'
      q12= ' f2=(b-c)/b                                      :'
      q13= ' f3=(a-b)/a                                      :'
      q14= ' f_Jeans                                         :'
      q15= ' f_MacLaurin                                     :'
      q16= ' f_Roche                                         :'

      return
      end

      subroutine inverse(y,z,n,m)
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
      implicit none 
      integer n,m
      real*8 y(m,m),z(m,m)
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k

      do i=1,n
       do j=1,n
        a(i,j) = y(i,j)
       end do
      end do

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0

! step 1: forward elimination
      do k=1, n-1
       do i=k+1,n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
        end do
       end do
      end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
      do i=1,n
       L(i,i) = 1.0
      end do
! U matrix is the upper triangular part of A
      do j=1,n
       do i=1,j
        U(i,j) = a(i,j)
       end do
      end do

! Step 3: compute columns of the inverse matrix C
      do k=1,n
       b(k)=1.0
       d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
       do i=2,n
        d(i)=b(i)
        do j=1,i-1
         d(i) = d(i) - L(i,j)*d(j)
        end do
       end do
! Step 3b: Solve Ux=d using the back substitution
       x(n)=d(n)/U(n,n)
       do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
         x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
       end do
! Step 3c: fill the solutions x(n) into column k of C
       do i=1,n
        c(i,k) = x(i)
       end do
       b(k)=0.0
      end do 

      do i=1,n
       do j=1,n
        z(i,j) = c(i,j)
       end do
      end do

      end subroutine inverse

