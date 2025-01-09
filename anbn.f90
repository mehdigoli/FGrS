
    module fft_anbn

    use global_variables

    implicit none

    contains

    !anbn computes the lump coefficient for several gravity related quantities
    !input:
    ! nthread ..... integer, number of computational threads
    ! tip     ..... integer, flag for different gravity realted quantities
    ! cnm, snm .... 1d double array, Cosine/Sine SH coefficients
    ! cnmH/snmH ... 1d double array, SH coeff. only for quant = 'geoid'
    ! tf        ... 2d double array, eigenvalues of several gravity related quantities.
    ! t1, u1,   ... double scalar, sin(phi), cos(phi), phi: geocentric lat.
    ! ips       ... 1d integer array, the ix part of Xnumber of sectorial Pnm
    ! ps        ... 1d double array,  the  x part of Xnumber of sectorial Pnm
    ! xi, eta   ... 1d array, precomputed useful array for computation tesseral Pnm
    ! Lp        ... integer scalar, optimized SHS for polar region using Jekeli's rule.
    ! nmin, nmax .. integer scalar, min/max degree of SHS
    ! output
    ! an, bn    ... 2d double array, lumb coefficients
    ! example an(0:nmax,1:col) , col number of different gravity related quantities.
    !
    ! programmer: Mehdi Goli.

    subroutine anbn(cnm, snm, cnmH, snmH, tf, t1, u1, an, bn)

    use alf

    implicit none

    real(8), allocatable , intent(in):: cnm(:), snm(:), cnmH(:), snmH(:)

    real(8), intent(in) :: tf(0:,:)

    real(8), allocatable, intent(inout):: an(:,:), bn(:,:)

    real(8), intent (in):: u1, t1

    real(8) :: tdu   , tfi

    integer:: n, m0, mp , m

    integer(8):: k

    real(8), allocatable:: pnm(:), dpnm(:), d2pnm(:),vm(:),tempc(:), temps(:)

    allocate( pnm(0:nmax+1) )

    an=0d0;

    bn=0d0

    tdu = t1/u1

    pnm(0:nmax+1)=0d0;

    m0 = 8d-3*nmax+100

    k=1

    do  n = 0, nmin-1
        k=k+n
        call alf_FRM ( n, n, t1, u1, a_n, b_n, sqr, zsqr, mu(n), pnm)
    enddo

    select case (tip)

        !-------------------------  !W, T, sa_dg, sa_ga, dgdr, sa_dgadr, Z, height, 'none'

    case(1, 3)

        do  n = nmin, nmax

            k=k+n

            mp = u1*n + m0

            mp = min(n,mp)

            call alf_FRM ( n, mp, t1, u1, a_n, b_n, sqr, zsqr, mu(n), pnm)

            tfi=  tf(n,1)

            do m=0, mp

                an(m,1)=an(m,1)+   cnm(k+m)*pnm(m)  *tfi
                bn(m,1)=bn(m,1)+   snm(k+m)*pnm(m)  *tfi

            enddo

        enddo

    case(2) ! geoid

        allocate(tempc(0:nmax), temps(0:nmax))

        do n = nmin, nmax

            mp = u1*n + m0

            mp = min(n,mp)

            call alf_FRM( n, mp, t1, u1,  a_n, b_n, sqr, zsqr, mu(n), pnm)

            k=k+n

            tempc(0:mp)= cnm(k:k+mp)*pnm(0:mp)
            temps(0:mp)= snm(k:k+mp)*pnm(0:mp)

            an(0:mp,1)= an(0:mp,1)+tf(n,1)*tempc(0:mp) !C T(m)

            bn(0:mp,1)= bn(0:mp,1)+tf(n,1)*temps(0:mp)!S T(m)

            an(0:mp,2)=an(0:mp,2)+tf(n,2)*tempc(0:mp) !C T(m),dT/dr(m)

            bn(0:mp,2)=an(0:mp,2)+tf(n,2)*temps(0:mp) !S T(m),dT/dr(m)

            an(0:mp,3)=an(0:mp,1)+ cnmH(k:k+mp)*pnm(0:mp)  ! C H(m)

            bn(0:mp,3)=an(0:mp,1)+ snmH(k:k+mp)*pnm(0:mp)  ! S H(m)

        enddo

    case(4) !DOV

        allocate(dpnm(0:nmax), vm(0:nmax))

        vm = [ (n, n= 0 , nmax) ]

        dpnm=0d0

        do n = nmin, nmax

            mp = u1*n + m0

            mp = min(n,mp)

            call alf_FRM( n, n, t1, u1,  a_n, b_n, sqr, zsqr, mu(n), pnm)

            call dalfx_RW (n, sqr, pnm, dpnm)

            k=k+n

            an(0:n,1)= an(0:n,1)+tf(n,1)* cnm(k:k+n)*dpnm(0:n) !C dT/dphi

            bn(0:n,1)= bn(0:n,1)+tf(n,1)* snm(k:k+n)*dpnm(0:n) !S dT/dphi

            an(0:mp,2)= an(0:mp,2)+tf(n,2)* vm(0:mp) * snm(k:k+mp)*pnm(0:mp) !C dT/dl

            bn(0:mp,2)= bn(0:mp,2)+tf(n,2)* (-vm(0:mp) * cnm(k:k+mp)*pnm(0:mp)) !S dT/dl


        enddo

    case(5) !g, gg, dT, dW, dV

        allocate(dpnm(0:nmax), tempc(0:nmax), temps(0:nmax), vm(0:nmax))

        vm = [ (n, n= 0 , nmax) ]

        do n = nmin, nmax

            mp = u1*n + m0

            mp = min(n,mp)

            call alf_FRM( n, n, t1, u1,  a_n, b_n, sqr, zsqr, mu(n), pnm)
            call dalfx_RW (n, sqr, pnm, dpnm)

            k=k+n

            tempc(0:mp)= cnm(k:k+mp)*pnm(0:mp)

            temps(0:mp)= snm(k:k+mp)*pnm(0:mp)

            an(0:mp,1)= an(0:mp,1)+tf(n,1)* tempc(0:mp) !C d*/dr

            bn(0:mp,1)= bn(0:mp,1)+tf(n,1)* temps(0:mp) !S d*/dr

            an(0:mp,2)= an(0:mp,2)+tf(n,2)* ( vm(0:mp) * temps(0:mp)) !C   d*/dlam

            bn(0:mp,2)= bn(0:mp,2)+tf(n,2)* ( -vm(0:mp) * tempc(0:mp)) !S   d*/dlam

            an(0:n,3) = an(0:n,3)+tf(n,3) * ( cnm(k:k+n)*dpnm(0:n)) !C d*/dphi

            bn(0:n,3) = bn(0:n,3)+tf(n,3) * ( snm(k:k+n)*dpnm(0:n)) !S d*/dphi

        enddo

    case(7) ! dg

        allocate(dpnm(0:nmax), tempc(0:nmax), temps(0:nmax), vm(0:nmax))

        vm = [ (n, n= 0 , nmax) ]

        do n = nmin, nmax

            mp = u1*n + m0

            mp = min(n,mp)

            call alf_FRM( n, n, t1, u1,  a_n, b_n, sqr, zsqr, mu(n), pnm)
            call dalfx_RW (n, sqr, pnm, dpnm)
            k=k+n

            tempc(0:mp)= cnm(k:k+mp)*pnm(0:mp)

            temps(0:mp)= snm(k:k+mp)*pnm(0:mp)

            an(0:mp,1)= an(0:mp,1)+tf(n,1)* tempc(0:mp) !C d*/dr

            bn(0:mp,1)= bn(0:mp,1)+tf(n,1)* temps(0:mp) !S d*/dr

            an(0:mp,2)= an(0:mp,2)+tf(n,2)* tempc(0:mp) !C d*/dr

            bn(0:mp,2)= bn(0:mp,2)+tf(n,2)* temps(0:mp) !S d*/dr

            an(0:mp,3)= an(0:mp,3)+tf(n,3)* ( vm(0:mp) * temps(0:mp)) !C   d*/dlam

            bn(0:mp,3)= bn(0:mp,3)+tf(n,3)* ( -vm(0:mp) * tempc(0:mp)) !S   d*/dlam

            an(0:n,4) = an(0:n,4)+tf(n,4) * ( cnm(k:k+n)*dpnm(0:n)) !C d*/dphi

            bn(0:n,4) = bn(0:n,4)+tf(n,4) * ( snm(k:k+n)*dpnm(0:n)) !S d*/dphi


        enddo

    case(6) ! 2nd gradient d2W d2T

        allocate(d2pnm(0:nmax), dpnm(0:nmax), tempc(0:nmax), temps(0:nmax), vm(0:nmax))

        vm = [ (n, n= 0 , nmax) ]

        do n = 0, nmax

            mp = u1*n + m0

            mp = n!min(n,mp)

            call alf_FRM( n, n, t1, u1,  a_n, b_n, sqr, zsqr, mu(n), pnm)

            call dalfx_RW (n, sqr, pnm, dpnm)

            call dalfx_RW (n, sqr, dpnm, d2pnm)

            k=k+n
            tempc(0:mp)= cnm(k:k+mp)*pnm(0:mp)

            temps(0:mp)= snm(k:k+mp)*pnm(0:mp)

            an(0:mp,1)=an(0:mp,1)+ tf(n,1) * tempc(0:mp) ! C: dr
            bn(0:mp,1)=bn(0:mp,1)+ tf(n,1) * temps(0:mp) !S: dr

            an(0:mp,2)= an(0:mp,2)+tf(n,2) *  vm(0:mp) * temps(0:mp) !C   d*/dlam
            bn(0:mp,2)= bn(0:mp,2)+tf(n,2) *  (-vm(0:mp) * tempc(0:mp)) !S   d*/dlam

            an(0:mp,3)=an(0:mp,3)+ tf(n,3) * tempc(0:mp) ! C: drdr
            bn(0:mp,3)=bn(0:mp,3)+ tf(n,3) * temps(0:mp) !S: drdr

            an(0:mp,4)= an(0:mp,4)+tf(n,4) * ( -vm(0:mp)*vm(0:mp) * tempc(0:mp))      ! C:  dldl
            bn(0:mp,4)= bn(0:mp,4)+tf(n,4) * (- vm(0:mp)*vm(0:mp) * temps(0:mp)) !S:  dldl

            an(0:mp,5)= an(0:mp,5)+tf(n,5) *   vm(0:mp) * temps(0:mp) ! C:  dldr
            bn(0:mp,5)= bn(0:mp,5)+tf(n,5) * (-vm(0:mp) * tempc(0:mp)) !S:  dldr

            tempc(0:n)= cnm(k:k+n)*dpnm(0:n)
            temps(0:n)= snm(k:k+n)*dpnm(0:n)

            an(0:n,6)= an(0:n,6)+tf(n,6)* tempc(0:n) !C dphi
            bn(0:n,6)= bn(0:n,6)+tf(n,6)* temps(0:n) !S dphi

            an(0:n,7)= an(0:n,7)+tf(n,7)* tempc(0:n) !C  drdphi
            bn(0:n,7)= bn(0:n,7)+tf(n,7)* temps(0:n) !S  drdphi

            an(0:n,8)= an(0:n,8)+tf(n,8)*    vm(0:n) * temps(0:n)  !C  dldphi
            bn(0:n,8)= bn(0:n,8)+tf(n,8)* (-vm(0:n) * tempc(0:n))   !S  dldphi

            an(0:n,9)= an(0:n,9)+tf(n,9)* cnm(k:k+n)*d2pnm(0:n) !C dphidphi
            bn(0:n,9)= bn(0:n,9)+tf(n,9)* snm(k:k+n)*d2pnm(0:n) !S dphidphi

        enddo
        !stop
    end select

    end subroutine anbn

    !---------------------------------------------------------------------------------
    function eigengrav_ell (nx, func, r, aegm, GMegm) result(tf)

    ! eigengrav_ell returns the isotropic spectral transfer
    ! (or: eigenvalues) of several gravity related quantities.
    !
    ! input:
    !    nx .. integer argument, spherical harmonic degree
    !    func ..  integer argument , deing the functional under consideration:
    !    r .....  double scalar argument, geocentric radius [m].
    !    aegm ... double scalar argument, scale of SH coefficient (set aemg=1 for func='none' , 'height')
    !    GM   ... double scalar argument, GM of EGM file. (set gmemg=1 for func='none' , 'height')
    ! output:
    !    tf ..... double 2d array, eigenvalues of several gravity related quantities
    !
    ! programmer: M.Goli. translated to FORTRAN from SHbundle software

    !  use global_variables

    implicit none

    integer, intent(in):: nx

    integer, intent(in):: func

    real(8), intent(in) :: r, GMegm, aegm

    real(8), allocatable :: tf(:,:)

    integer:: n

    real(8) , allocatable:: tfW(:)

    !------------------------------------------------------

    if (.not.(allocated(tfw))) allocate(tfW(0:nx))
    tfw=0d0;

    do n=0, nx

        tfW(n) =  GMegm/aegm * (aegm/r)**(n+1)

    enddo

    select case (func)

    case( 1:3 )      ! V,W,T,Z

        if (.not.(allocated(tf))) allocate(tf(0:nx,1))
        tf=0d0

        tf(: , 1) = tfW

    case ( 4  )  !   geoid

        if (.not.(allocated(tf))) allocate(tf(0:nx,3))
        tf=0d0;

        tf(: , 1) = tfW

        do n=0, nx
            tf(n , 2) =  tfW(n) * (n+1)/r !  ! d*/dr
        enddo

        tf(: , 3) = 1.0d0 !

    case (5 )    !spherical approximation of dg

        if (.not.(allocated(tf))) allocate(tf(0:nx,1))
        tf=0d0
        do n=0, nx
            tf(n, 1) = tfW(n) * (n+1)/r *1d5 ! [mGal]
        enddo

    case (6)        !spherical approximation of dg

        if (.not.(allocated(tf)))allocate(tf(0:nx,1))
        tf=0d0
        do n=0, nx
            tf(n, 1) = tfW(n) * (n-1)/r *1d5 ! [mGal]
        enddo


    case(7,8)    !  = CS coef. defines output , height

        if (.not.(allocated(tf)))allocate(tf(0:nx,1))
        tf = 1.0d0

    case ( 9 ) !DOV

        if (.not.(allocated(tf)))allocate(tf(0:nx,2))
        tf=0d0
        tf(:, 1) =  - tfW * 180.d0 * 3600.d0 /pi !xi  [arcsec]

        tf(:, 2) =  - tfW  * 180.d0 * 3600.d0 /pi  !eta [arcsec]

    case (10,13) !gravity anomaly


        if (.not.(allocated(tf)))allocate(tf(0:nx,4))

        tf=0d0

        tf(:,1) = tfW

        do n=0, nx
            tf(n,2) = -tfW(n) * (n+1)/r*1d5        !dwdr  [mGal]
        enddo

        tf(:,3) = tfW *1d5                        !dwdl  [mGal]

        tf(:,4) = tfW*1d5                          !dwdphi   [mGal]

    case ( 11,12,14, 15) ! fisrt gradient of V

        if (.not.(allocated(tf)))allocate(tf(0:nx,3))

        tf=0d0

        do n=0, nx
            tf(n,1) = -tfW(n) * (n+1)/r*1d5        !dwdr  [mGal]
        enddo

        tf(:,2) = tfW *1d5                        !dwdl  [mGal]

        tf(:,3) = tfW*1d5                          !dwdphi   [mGal]

    case ( 16:17 )

        if (.not.(allocated(tf)))allocate(tf(0:nx, 9))

        tf=0d0

        do n=0, nx
            tf(n, 1) =  -tfW(n) * (n+1)/r   !dr
            tf(n, 3) =   tfW(n) * (n+1)*(n+2)/r/r   !drdr
        enddo

        tf(:, 2) = tfW (:)                                !dlam

        tf(:, 4)  = tfW (:)                              ! dlamdlam

        tf(:, 5) =  tf(:,1)   !  drdlam

        tf(:, 6) = tfW (:)     !dphi

        tf(:, 7) =  tf(:,1)    !drdphi

        tf(:, 8)  = tfW (:)    !, dphidlam

        tf(:, 9)  = tfW (:)    ! dphidphi

        tf=tf*1d9 ! Etovos

    end select

    end function eigengrav_ell

    !---------------------------------------------------------------------------------
    !
    !subroutine inv_spec_v(Am, Bm, nc, lam, f, c_shf, s_shf)
    !
    !! in_spec_v computes the function from Fourier coef an, bn
    !! before the apply the ifft, am,bm should be shifted using shift theorem
    !! for coarse grid ifft is replaced by direct summation.
    !! input:
    !! an, bn ....... 2d double array, lumb coefficients.
    !! nc     ....... integer parameter = 360/step of grid
    !! lam    ....... 1d double array, array of longitude (deg).
    !! c_chf,s_chf .. 1d double array (optional). if exist shift the an, bn, otherwise no
    !! output
    !! f (1:nlam, 1:col) ... output values for all point on longitude band
    !!
    !! programmer: M. Goli
    !
    !use Fast_Fourier
    !
    !implicit none
    !
    !real(8) :: am(0:,:), bm(0:,:)
    !
    !real(8), allocatable :: f(:,:)
    !
    !real(8), allocatable, intent(in) :: lam(:)
    !
    !integer, intent(in) :: nc
    !
    !real(8), allocatable , intent (in) , optional:: c_shf(:), s_shf(:)
    !
    !real(8), allocatable :: Am_(:,:), Bm_(:,:), temp(:,:)
    !
    !integer::  nlam, i, col, j, nmax
    !
    !real(8), allocatable :: Am1(:), Bm1(:)
    !
    !real(8), allocatable::  cosml(:),  sinml(:)
    !
    !integer, allocatable :: vm(:)
    !
    !real(8) :: rlami , s
    !
    !integer:: nc1, is , nm , m
    !
    !f=0d0
    !nlam=size(lam)
    !
    !nmax =size(Am,1)-1
    !col =size(Am, 2)
    !
    !nc1=nc
    !
    !do while (nmax+1>nc1/2)
    !    nc1=nc1*2
    !    nlam=nlam*2
    !enddo
    !
    !if (nc1/nc>30)then
    !
    !    allocate( cosml(0:nmax), sinml(0:nmax), vm(0:nmax))
    !
    !    do m=0, nmax
    !        vm(m)=m
    !    enddo
    !
    !    do i=1, col
    !
    !        do j=1, size(lam)
    !
    !            rlami = lam(j) * dr
    !
    !            cosml(vm) = cos(vm*rlami);
    !
    !            sinml(vm) = sin(vm*rlami)
    !
    !            s=0d0
    !
    !            do m=0,nmax
    !                s= s +    am(m,i)*cosml(m) +bm(m,i)*sinml(m)
    !            enddo
    !
    !            f (j, i) = s
    !
    !        enddo
    !
    !    enddo
    !
    !    return
    !
    !endif
    !!
    !!----- shift am bm if need
    !!
    !if (present(c_shf) .and. present(s_shf) ) then
    !
    !    if (.not.(allocated(temp))) allocate(temp(0:nmax,col))
    !
    !    do i=1,col
    !
    !        temp (0:nmax,i) = am(0:nmax,i) * c_shf(0:nmax) - &
    !            bm(0:nmax,i) * s_shf(0:nmax)
    !
    !        bm(0:nmax,i) =  bm(0:nmax,i) * c_shf(0:nmax) + &
    !            am(0:nmax,i) * s_shf(0:nmax)
    !
    !        am(0:nmax, i) = temp (0:nmax,i)
    !
    !    enddo
    !
    !endif
    !
    !nm=max(nmax+1,nc1)
    !
    !if (.not.(allocated(Am_))) allocate (Am_(nm,col) , Bm_(nm,col), Am1(nc1) , Bm1(nc1))
    !
    !Am_= 0.d0; Bm_ = 0.d0; Am1 = 0; Bm1 = 0; f = 0.d0
    !
    !do i=1, nmax+1
    !
    !    Am_ (i, 1:col)  = Am(i-1, 1:col)
    !    Bm_ (i, 1:col)  = Bm(i-1, 1:col)
    !
    !enddo
    !
    !Am_(1,:) = Am_(1,:) * 2.0d0;
    !
    !is= nlam/size(lam)
    !
    !do i = 1, col
    !
    !    Am1=0d0;
    !    Bm1=0d0;
    !
    !    AM1(1:nmax+1)=  Am_(1:nmax+1, i)  ;  Am1(1)=Am1(1)/2;
    !    BM1(1:nmax+1)=  Bm_(1:nmax+1, i)
    !    call fft1 (Am1, bm1, nc1, +1, j)
    !
    !    f(1:nlam,i)=AM1(1:nc1:is)
    !
    !enddo
    !
    !end subroutine inv_spec_v
    subroutine inv_spec_v(An, Bn, nc, lam, f, c_shf, s_shf )

    !in_spec_v computes the function from Fourier coef an, bn
    !before the apply the ifft, an,bn should be shifted using shift theorem
    !ifft performs using FFTW3 library.
    !for coarse grid ifft is replaced by direct summation.
    !input:
    !an, bn ....... 2d double array, lumb coefficients.
    !nc     ....... integer parameter = 360/step of grid
    !lam    ....... 1d double array, array of longitude (deg).
    !c_chf,s_chf .. 1d double array (optional). if exist shift the an, bn, otherwise no
    !output
    !f (1:nlam, 1:col) ... output values for all point on longitude band

    !programmer: M. Goli

    use global_variables

    use FFTW3

    implicit none

    real(8) ,intent(inout) :: an(0:,:), bn(0:,:) , f(:,:)

    real(8), allocatable, intent(in) :: lam(:)

    ! integer, intent(in) :: nc
    integer:: nc

    real(8), allocatable , intent (in) , optional:: c_shf(:), s_shf(:)

    real(8), allocatable :: An_(:,:), Bn_(:,:), d(:), temp(:,:)!, f1(:,:)

    integer::  nlam, i, nc1, nlam1, k, j

    complex(c_double_complex), allocatable:: cn(:)

    real(8) :: rlami

    type(c_ptr) :: plan

    !============================================================

    nmax =size(An,1)-1

    col =size(An, 2)

    nlam=size(lam)

    nlam1=nlam

    nc1=nc
    
    !f1=f;

    do

        if (nmax+1<=nc1/2) exit
        nc1=nc1*2
        nlam1=nlam1*2

    end do

    if ( nc1/nc >100 ) then

        do i=1, col

            do j=1, nlam

                f(j, i) = sum( an(0:nmax,i)*cosml(:,j)+bn(0:nmax,i)*sinml(:,j))

            enddo

        enddo

        return

    endif

    !----- shift an bn if need

    if (present(c_shf) .and. present(s_shf) ) then

        allocate(temp(0:nmax,col))

        do i=1,col

            temp (0:nmax,i) = an(0:nmax,i) * c_shf(0:nmax) - &
                bn(0:nmax,i) * s_shf(0:nmax)

            bn(0:nmax,i) =  bn(0:nmax,i) * c_shf(0:nmax) + &
                an(0:nmax,i) * s_shf(0:nmax)

            an(0:nmax, i) = temp (0:nmax,i)

        enddo

    endif

    An(0,:) = An(0,:) * 2.0d0;

    allocate( cn( nc1) , d(nc1) )

    allocate (An_(1:nc1, 1:col) , Bn_(1:nc1, 1:col))

    An_ =0.d0

    Bn_ =0.d0

     f = 0.d0

    do i=1, nmax+1

        An_ (i, 1:col)  = An(i-1 , 1:col)
        Bn_ (i, 1:col)  = Bn(i-1 , 1:col)

    enddo

    do i = 1, col

        cn  = cmplx(an_(:,i), -bn_(:,i), kind(1.d0))  /2.0d0

        d=0.0d0

        !$OMP CRITICAL

        plan = fftw_plan_dft_c2r_1d(nc1, cn(1:nc1), d(1:nc1), FFTW_ESTIMATE+FFTW_UNALIGNED)

        call dfftw_execute_dft_c2r(plan, cn(1:nc1), d(1:nc1))

        call dfftw_destroy_plan(plan)

        !$OMP end CRITICAL

        if (nc1==nc) then
          
            f(1:nlam, i) = d(1:nlam)
            
        else
            j=1
            do k = 1, nlam1, nc1/nc
                f(j, i) = d(k)
                j = j + 1
            end do 
            
        endif
        
        !print*, maxval(abs(f1-f))
        !!
        !!    j = 1
        !!    do k = 1, nc1, nc1/nc
        !!        f(j, i) = d(k)
        !!        j = j + 1
        !!    end do
        !!endif

    enddo

    return
    end subroutine inv_spec_v


    !--------------------------------------------------------------------------

    function nm2row(n, m)

    !convert (n,m)=> 1d array (row-wise)
    !00 ->1
    !10 ->2
    !11 ->3
    !20 ->4
    !....
    !programmer:M.Goli

    integer, intent(in):: m, n

    integer(8):: nm2row

    nm2row=n*(n+1)/2+m+1

    end function nm2row



    end module fft_anbn