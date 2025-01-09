    module alf

    use global_variables

    implicit none

    contains

    !--------------------------------------------------------------------------------

    subroutine Pn0 (t1, p)

    ! Pn0 computes the Legendre polynomials Pn(t1) up to degree nmax

    !input :
    !  - nmax : integer argument, maximum degree of expansion
    !  - t1   : double argument
    !output
    !  - p (0:nmax) , 1d double array, Legendre polynomials

    implicit none

    real(8), intent(in):: t1

    real(8),  allocatable, intent(inout)::  p(:)

    real(8):: a1, a2, a3

    integer:: i

    p(0)=1.0d0 ;

    p(1)=sqrt(3d0)*t1;

    do i=2, nmax

        a1=dble(2*i+1)*dble(2*i-1)

        a2=dble(2*i+1)

        a3=dble(2*i-3)

        p(i)=p(i-1)*sqrt(a1)/i*t1 - p(i-2)*(i-1)*sqrt(a2)/(i*sqrt(a3))

    END do

    return

    end subroutine pn0

    !--------------------------------------------------------------------------

    subroutine dalfx_RW (n, sqr, dp, ddp)

    integer, intent(in) :: n

    real(8), intent(in) :: sqr(-1:), dp(0:)

    real(8), intent(inout) :: ddp(0:)

    !real(8), allocatable, intent (in) :: sqr(:), dp(:)
    !
    !real(8), allocatable, intent (inout):: ddp(:)

    integer:: m


    if (n==0) then

        ddp(0) =0.d0

        return

    endif

    if (n==1) then

        ddp(0) = dp(1)

        ddp(1) = - dp(0)

        return

    endif

    ddp(0) = sqr(n)*sqr(n+1) * sqr(2) * dp(1)

    ddp(1) = sqr(n+2)*sqr(n-1) * dp(2) - &
        sqr(2) * sqr(n+1)*sqr(n)  * dp(0)

    do m=2, n-1
        ddp(m)  = sqr(n+m+1)*sqr(n-m)  *dp(m+1) - &
            sqr(n+m)  *sqr(n-m+1)*dp(m-1)
    enddo

    ddp(n) = - sqr(2)*sqr(n) * dp(m-1)

    ddp(0:n) = ddp(0:n) / 2d0

    end subroutine dalfx_RW

    !--------------------------------------------------------------------------


    subroutine alf_FRM( n, mp, t, u, a_n, b_n, sqr, zsqr, mun, pn)

    ! alf_RFM computes the P_n,0:n uisng 4 terms  RW recursive formula
    ! reference, Xing et al., 2019, Journal of Geodesy (JoG)
    ! this algorithm optimized for polar region. i.e., computaion is perform up to
    ! mp<=n mp is input argument, mp = n sin(theta) + mu.
    ! mu is the empirical constant, cf. Balmino et al., 2012 JoG for details.
    !
    ! input :
    ! n, mp<=n .......... integer argument , degree, optimized order
    ! t ................. double scalar argument, cos(theta)
    ! u ................. double scalar argument, sin(theta)
    ! an, bn, sqr, zsqr.. double 1d array
    ! mun ............... double scalar argument
    ! input/output
    ! pn ................ double 1d array,
    ! pn as input contains P_(n-1,0:n)
    ! pn as output contains P_(n,0:n)
    !
    ! programmer: M.Goli

    integer , intent (in) :: n

    integer, intent(in):: mp

    real(8), dimension(:), allocatable , intent(in) ::  sqr, zsqr, a_n, b_n

    real(8), dimension(0:), intent(inout) :: pn

    real(8), intent(in) :: t, u, mun

    real(8):: p1, p2, p3, pnm, pn0, pn1

    integer:: m

    real(8):: cnm, dnm, enm, ud2

    ud2 = u / 2.d0

    if (n==0) then

        pn(0) =1.d0

        return

    endif

    if (n==1) then

        pn(0) = sqr(3) * t

        pn(1) = sqr(3) * u

        return

    endif

    pn0 = a_n(n) * t * pn(0) - b_n(n) * ud2 * pn(1)

    m=1

    cnm = sqr(n+m)*sqr(n-m)

    dnm = sqr(n-m)*sqr(n-m-1)

    enm = sqr(n+m)*sqr(n+m-1)*sqr(2)

    pn1 = cnm * t * Pn(m)  - ud2 * ( dnm * pn(m+1) -enm *pn(m-1) )

    p1=pn(1)

    p2=pn(2)

    do m=2, mp

        p3=pn(m+1)

        cnm = sqr(n+m)*sqr(n-m)

        dnm =zsqr(n-m)! sqr(n-m)*sqr(n-m-1)

        enm =zsqr(n+m)! sqr(n+m)*sqr(n+m-1)

        pnm = cnm * t * P2  - ud2 * ( dnm * p3 -enm *p1 )

        p1=pn(m)

        pn(m)=pnm

        p2=p3

    enddo

    pn(1)=pn1

    pn(1:mp)=pn(1:mp)*mun

    pn(0)=pn0

    end subroutine alf_FRM

    !------------------------------------------------------------------------
    subroutine FRM_ini(nmax, mu, sqr, zsqr, an, bn)

    ! RFM_ini prepares some initial array for ALF computation by
    ! RFM method
    ! input
    !   nmax ..... integer argument, maximum degree of exapnsion
    ! output
    !   mu,sqr, zsqr, an, bn ..... double 1d array
    !programmer: M.Goli

    integer, intent(in) :: nmax

    real(kind = dp), allocatable, intent(out) ::  sqr(:), mu(:), an(:), bn(:), zsqr(:)

    integer:: n

    allocate( mu(0:nmax), sqr( - 1:2*nmax + 3), an(nmax), bn(nmax), zsqr( - 1:2*nmax + 3))

    mu(0) = 0

    do n = 1, nmax

        an(n) = sqrt(2.d0*n + 1.d0)/sqrt(2.d0*n - 1.d0)

        bn(n) = sqrt (2.d0*(n - 1.d0)/n) * an(n)

        mu(n) = 1.d0/dble(n) * an(n)

    enddo

    sqr( - 1) = 0.d0

    do n = 0, 2*nmax + 3

        sqr(n) = sqrt(dble(n))

    enddo

    do n = 1, 2*nmax + 3

        zsqr(n) = sqrt(dble(n)*dble(n - 1))

    enddo

    zsqr(0) = 0.d0

    end subroutine FRM_ini

    !--------------------------------------------------------------

    end module alf