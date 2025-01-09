
    module shs

    use global_variables
    use alf
    use normal_field_
    use lagInterp
    use fft_anbn
    use utility

    implicit none

    contains

    function shs_grid(show, cnm, snm, cnmH, snmH, lat, lam, height, GMegm, aegm) result(f)

    ! Spherical Harmonics Synthesis for grid
    ! INPUT:

    !     sub_normal ........ logical argument, subtract normal field true/false
    !     normal_field....... string argument, set normal filed ('WGS84' or 'GRS80')
    !                         false for func ='none' and 'height'
    !     nthread     ....... integer argument, Number of parallel threads used to run.
    !     cnm, snm    ....... 1d arrays argument, set of SH coefficients
    !     cnmH, snmH  ....... 1d arrays argument, set of SH coef of topography, only required for option%func = 'geoid'
    !     nmin, nmax   ...... integer argument, min/maximum degree/order
    !     option%func      .......  string argument, defining the field option%funcity.
    !     lat, lon   ... .... 1D array, vector of geodetic latitudes, longitudes
    !     height      ....... double scalar, ellipsoidal height of grid in meter
    !     GMegm, aegm ....... GM and scale of EGM. (set GM==1d0, aegm==1d0 for option%func== none and height)
    !
    ! OUTPUT:
    !    f .................. 3D double array, grid file with size(lat) X size (lon) X n (dim of output function)


    implicit none

    logical , intent(in) :: show

    real(8), allocatable , intent(in):: cnm(:), snm(:), cnmH(:), snmH(:)

    real(8), allocatable, intent(in):: lat(:), lam(:)

    real(8), intent(in) :: aegm , GMegm, height

    real(8), allocatable::  f(:,:,:)

    integer:: m, j

    integer(4):: wid

    integer:: done

    real(8), allocatable :: fi(:,:)  , fii(:,:,:)

    integer:: nlat, nlon

    real(8):: west, westm

    real(8), allocatable:: c_shf(:), s_shf(:)

    integer::  i, nc, nc1

    nlat = size(lat);

    nlon = size(lam)

    if (allocated(f)) deallocate(f)
    if (allocated(fi)) deallocate(fi)

    allocate( f(1:nlat, 1:nlon, iprint))
    allocate( fi(1:nlon, iprint))
    allocate( fii(1:2, 1:nlon, iprint))

    f=0d0; fi=0; fii=0;

    allocate( c_shf(0:nmax), s_shf(0:nmax))
    c_shf=0; s_shf=0

    West = minval(lam)

    do m=0, nmax

        westm = west * dble(m) * dr
        c_shf(m) = cos (westm)
        s_shf(m) = -sin(westm)

    enddo

    nc = nint(360.d0/option%step_lon)
    nc1=nc
    
    do

        if (nmax+1<=nc1/2) exit
        nc1=nc1*2

    end do

    if (nc1/nc>100) then

        allocate(cosml(1:nmax+1,nlon), sinml(1:nmax+1,nlon))

        do m = 0, nmax
            cosml(m+1, :) = cos(m *  lam*dr  )
            sinml(m+1, :) = sin(m *  lam*dr  )
        end do
        
    endif

    call  FRM_ini(option%nmax, mu, sqr, zsqr, a_n, b_n)

    if ( show) then
        done=0; wid=0;
    endif

    !$OMP PARALLEL PRIVATE(i,j, fi) &
    !$OMP SHARED(f,cnm,snm,cnmH,snmH,option,lat,lam,height,GMegm, aegm,c_shf, s_shf)
    !$OMP do    schedule (dynamic)

    do i=1, size(lat)

        fi = shs_parallel(cnm, snm, cnmH, snmH, lat(i), lam, height, nlon, GMegm, aegm, c_shf, s_shf)

        f(i,:,:)=fi

        if (show) then

            done=done+1

            call progress(wid, done, nlat)

        endif

    enddo

    !$OMP END do

    !$OMP END PARALLEL

    deallocate(c_shf, s_shf)

    end function shs_grid



    function shs_parallel(cnm, snm, cnmH, snmH, lat, lam, height, &
        nlon, GMegm, aegm, c_shf, s_shf ) result (f)

    implicit none

    real(8), allocatable , intent(in):: cnm(:), snm(:)

    real(8), allocatable , intent(in) :: cnmH(:), snmH(:)

    real(8), allocatable, intent(in):: lam(:), c_shf(:), s_shf(:)

    real(8), allocatable:: f(:,:)

    real(8), intent(in) :: aegm , GMegm, height, lat

    integer, intent(in)::  nlon

    real(8), allocatable :: bm(:, :), am(:, :) , pnm(:), tf(:, :)

    real(8):: u1, t1

    real(8):: ri, thetai

    integer:: nadd

    real(8), allocatable :: fi(:,:)

    real(8):: anu

    real(8) :: C2n(0:nmax)

    integer::  nc, n

    if ( nlon>1 ) then

        nc = nint(360.d0/option%step_lon)

    endif

    !- allocate output matrix

    if (.not.(allocated(fi))) allocate(fi(1:nlon,col))
    fi=0d0

    if (.not.(allocated(f))) allocate(f(1:nlon,iprint))
    f=0.d0

    nadd=60 + 8d-3*nmax

    !------------- geodetic to spherical coordinate transformation

    call geod2sphe(lat, height, ri, thetai)

    if (.not.(allocated(bm))) allocate(bm(0:nmax, col), am(0:nmax, col), pnm(0:nmax) )

    pnm=0

    tf = eigengrav_ell(nmax, option%func, ri, aegm, GMegm)

    am=0.0; bm=0.0; anu=0.d0

    u1=sin(thetai*dr)

    t1=cos(thetai*dr)
    fi=0.d0

    call anbn(cnm, snm, cnmH, snmH, tf, t1, u1 , am, bm)

    call inv_spec_v(Am, Bm, nc, lam, fi, c_shf= c_shf, s_shf= s_shf )

    aNu=0d0;

    if (option%func == 10)then

        do n=0,nmax
            C2n(n) = klm(n)*GMe/GMegm*(ae/aegm)**(n+1)
        enddo

        call pn0(t1, pnm)

        aNu= sum(C2n(nmin:nmax)*pnm(nmin:nmax)*tf(nmin:nmax,1)) ! normal potential

    endif

    f =  f2func (aNu, ri, omega, u1, t1, lat, option%func, fi, height, nlon)

    end function shs_parallel

    !------------------------------------------------------------------

    function shs_scatter(cnm, snm, cnmH, snmH, lat, lam, h, GMegm, aegm) result(f)

    ! Spherical Harmonics Synthesis for grid
    ! INPUT:

    !     sub_normal ........ logical argument, subtract normal field true/false
    !     normal_field....... string argument, set normal filed ('WGS84' or 'GRS80')
    !                         false for func ='none' and 'height'
    !     nthread     ....... integer argument, Number of parallel threads used to run.
    !     cnm, snm    ....... 1d arrays argument, set of SH coefficients
    !     cnmH, snmH  ....... 1d arrays argument, set of SH coef of topography, only required for func = 'geoid'
    !     lat, lon,h   ... .... 1D array, vector of geodetic latitudes, longitudes, height
    !
    !     GMegm, aegm ....... GM and scale of EGM. (set GM==1d0, aegm==1d0 for func== none and height)
    !
    ! OUTPUT:
    !    f .................. 3D double array, grid file with size(lat) X size (lon) X n (dim of output function)

    implicit none

    real(8), allocatable , intent(in):: cnm(:), snm(:), cnmH(:), snmH(:)

    real(8), allocatable, intent(in):: lat(:), lam(:) , h(:)

    real(8) , allocatable :: f(:,:)

    real(8), intent(in) :: aegm , GMegm

    integer(4):: wid

    integer:: done, npoint,   i

    real(8), allocatable :: fi(:)

    npoint = size(lat);

    allocate( f(npoint, iprint), fi(iprint) )

    call  FRM_ini(option%nmax, mu, sqr, zsqr, a_n, b_n)

    done=0; wid=0

    call progress(wid, 0, npoint)

    !$OMP PARALLEL PRIVATE(i, fi ) &
    !$OMP SHARED(npoint,f,nmin,nmax,cnm,snm,cnmH,snmH,option,lat,lam,h,GMegm, aegm )
    !$OMP do !schedule (dynamic)

    do i=1,npoint

        fi  = shs_point(cnm, snm, cnmH, snmH, lat(i), lam(i), h(i), GMegm, aegm )
        f(i,1:iprint) = fi(1:iprint)

        done=done+1

        call progress(wid, done, npoint)

    enddo

    !$OMP END do

    !$OMP END PARALLEL

    end function shs_scatter

    !-----------------------------------------------------------------

    function shs_point(cnm, snm, cnmH, snmH, lat, lam, height, GMegm, aegm) result(f)

    implicit none

    real(8), allocatable , intent(in):: cnm(:), snm(:), cnmH(:), snmH(:)

    real(8), intent(in) :: aegm , GMegm, height, lat,  lam

    real(8), allocatable :: bm(:, :), am(:, :), tf(:, :), pnm(:),  cosml(:), sinml(:), &

        C2n(:), f(:), vf(:,:),  vfi(:, :)

    real(8):: u1, t1, ri, thetai, anu

    integer:: nadd, n, m

    !- allocate output matrix

    allocate( cosml(0:nmax), sinml(0:nmax), C2n(0:nmax), f(iprint), vf(1,iprint),  vfi(1, col) )

    f=0.d0; vfi=0d0

    nadd=60 + 8d-3*nmax

    !------------- geodetic to spherical coordinate transformation

    call geod2sphe(lat, height, ri, thetai)

    if (.not.(allocated(bm))) then
        allocate(bm(0:nmax, col), am(0:nmax, col), pnm(0:nmax+1))
    endif

    tf = eigengrav_ell(nmax, option%func, ri, aegm, GMegm)

    am=0.0d0; bm=0.0d0;

    u1=sind(thetai); t1=cosd(thetai)

    call anbn(cnm, snm, cnmH, snmH, tf, t1, u1 , am, bm)

    do m=0, nmax

        cosml(m) = cosd(m* lam);
        sinml(m) = sind(m* lam)

    enddo

    do m=0, nmax

        vfi (1,:) = vfi (1,:) + am(m,:)*cosml(m) + bm(m,:)*sinml(m)

    enddo

    aNu=0d0;

    if (option%func == 10)then

        do n=0,nmax
            C2n(n) = klm(n)*GMe/GMegm*(ae/aegm)**(n+1)
        enddo

        call pn0(t1, pnm)

        aNu= sum(C2n(nmin:nmax)*pnm(nmin:nmax)*tf(nmin+1:nmax+1,1)) ! normal potential

    endif

    vf =  f2func (aNu, ri, omega, u1, t1, lat, option%func, vfi, height, 1)

    f = vf(1,:)

    end function shs_point

    !------------------------------------------------------------------

    function f2func (aNu, ri, omega, u, t, lat, func, fi, height, nlon) result(f)

    integer:: func, nlon, j

    real(8):: anU, ri, omega, u, t, lat, height, gamma

    real(8), allocatable, dimension(:):: gg_r, gg_l, gg_f, gg, gravity, z, wr, wf, wl, wll, wlf, wff, wrr, wfr, wlr

    real(8):: gc_r, gc_f, gc, tdu

    real(8), allocatable:: fi(:,:), f(:,:)

    if (.not.(allocated(fi))) allocate (fi(1:nlon,1:col))
    allocate (f(1:nlon,1:iprint))

    tdu=t/u;

    select case (func)

    case(1)   ! W,V

        f( 1:nlon, 1) = fi(1:nlon, 1)! V
        f( 1:nlon, 2) = fi(1:nlon, 1)+0.5d0 * (ri*omega*u)**2 !W

    case(2,5:8)! " T, -dT/dr, 2/r*T-dT/dr , H ,none

        f( 1:nlon, 1) = fi(1:nlon, 1)! T

    case(3)! "3 height anomaly: Z = T/gamma_q             "

        do j = 1 , nlon

            call normal_gravity_p( height-fi(j,1) /9.81d0 , lat, gamma)

            f(j,1) = fi(j,1) /gamma

        enddo

    case (4) !"4 geoid undulation: N = T0/gamma_0 +ADCB    "

        call normal_gravity_p( -fi(j,1) /9.81d0 , lat, gamma)

        f( 1:nlon,1) = fi(1:nlon,1)*( 1.d0 + fi(1:nlon,2) +  &

            fi(1:nlon,3)**2 * 2.d0 * pi * G * rho) / gamma

    case(9) ! DOV (xi , eta)

        call normal_gravity_p( height , lat, gamma)

        f(1:nlon,1) = fi(1:nlon,1) /gamma / ri

        f(1:nlon,2) = fi(1:nlon,2) /gamma / ri / u

    case(10,13)! norm (dg) , norm(ga)

        if (.not.(allocated(gg_r))) allocate(gg_r(nlon), gg_l(nlon), gg_f(nlon))

        gg_r=0d0;gg_l=0d0; gg_f=0d0

        gg_r(1:nlon) =   fi(1:nlon, 2)

        gg_l(1:nlon) = -  fi(1:nlon, 3) / ri / u

        gg_f(1:nlon) =  fi(1:nlon, 4) / ri

        gc_r =  ri*omega**2 * u *u  *1d5

        gc_f = - ri*omega**2 * u *t  *1d5

        gc   = sqrt(gc_r**2+gc_f **2)

        gg   = sqrt(gg_r(1:nlon)**2+gg_f(1:nlon)**2+gg_l(1:nlon)**2)

        if (.not.(allocated(gravity))) allocate(gravity(1:nlon), z(1:nlon))

        gravity=0d0; z=0d0;

        gravity(1:nlon) = sqrt( (gg_r+gc_r)**2 + (gg_f+gc_f)**2 + gg_l**2 )

        z(1:nlon)= (fi(1:nlon,1)-anu) / 9.81d0

        !normal gravity at height (h-z)

        do j=1, nlon

            if (func==13) z(j)=0

            call normal_gravity_p(height-z(j), lat , gamma)

            f(j, 1) = gravity(j) - gamma*1d5

        enddo

    case(11) !gravity: gx, gy, gz (LNOF)

        if (.not.(allocated(gg_r))) allocate(gg_r(nlon), gg_l(nlon), gg_f(nlon))

        gg_r=0d0; gg_l=0d0; gg_f=0d0

        gc_r =  ri*omega**2 * u *u  *1d5

        gc_f = - ri*omega**2 * u *t  *1d5

        gg_r(1:nlon) = fi(1:nlon, 1)

        gg_l(1:nlon) = - fi(1:nlon, 2) / ri / u

        gg_f(1:nlon) =  fi(1:nlon, 3) / ri

        f(1:nlon,1) =  gg_f  + gc_f       ! dx
        f(1:nlon,2) =  gg_l               ! dy
        f(1:nlon,3) =  gg_r  + gc_r       ! dz

    case(12, 14) !gravity disturbance: dgx,dgy,dgz , dV/dx, dV/dy, dV/dz (LNOF)

        f(1:nlon,1) =  fi(1:nlon, 3) / ri        ! dx
        f(1:nlon,2) =  - fi(1:nlon, 2) / ri / u  ! dy
        f(1:nlon,3) =  fi(1:nlon, 1)             ! dz

    case(15) !dV/dphi, dV/dlam, dV/dr

        f(1:nlon,1) =  fi(1:nlon, 3)   ! d_dphi
        f(1:nlon,2) =  fi(1:nlon, 2)   ! d_lam
        f(1:nlon,3) =  fi(1:nlon, 1)   ! d_dr

    case(16,17) !  2nd gradient  d2T

        if (.not.(allocated(wl))) allocate(Wl(nlon), Wf(nlon), Wr(nlon), Wll(nlon),Wlf(nlon), &
            Wlr(nlon), Wff(nlon), Wfr(nlon), Wrr(nlon))

        Wl=0d0; Wf=0d0; Wr=0d0; Wll=0d0; Wlf=0d0; Wlr=0d0; Wff=0d0; Wfr=0d0; Wrr=0d0

        Wl = fi(1:nlon, 2) ; Wf = fi(1:nlon,6)  ;  Wr= fi(1:nlon,1);

        Wll = fi(1:nlon,4) ; Wlf = fi(1:nlon,8) ;  Wlr = fi(1:nlon,5)

        Wff = fi(1:nlon,9) ;         Wfr = fi(1:nlon,7)

        Wrr = fi(1:nlon,3)

        if (func == 16) then

            f(1:nlon, 1) = Wr / ri  + Wff / ri / ri          !Wxx

            f(1:nlon,2) = -1.d0 / ri/ri / u * Wlf - tdu / ri /ri / u * Wl          !Wxy

            f(1:nlon,3) = -1.d0 / ri/ri * Wf + 1.0d0 / ri * Wfr          !Wxz

            f(1:nlon,4) =  Wr / ri - tdu /ri/ri * Wf + 1.d0 /ri/ri/u/u * Wll                !Wyy

            f(1:nlon,5) = +1.d0 / ri/ri/ u  * Wl - 1.0d0 / ri /u* Wlr                   !Wyz

            f(1:nlon,6) = Wrr     !Wzz

        else

            f(1:nlon, 1) = Wr
            f(1:nlon, 2) = Wrr
            f(1:nlon, 3) = Wf
            f(1:nlon, 4) = Wff
            f(1:nlon, 5) = Wl
            f(1:nlon, 6) = Wll
            f(1:nlon, 7) = Wfr
            f(1:nlon, 8) = Wlr
            f(1:nlon, 9) = Wlf

        endif

    end select

    end

    function shs_igrid (cnm, snm, cnmh, snmh, lat, lam, h, GMegm, aegm) result(f)

    ! Spherical Harmonics Synthesis for grid with irregular top
    ! INPUT:

    !     cnm, snm    ....... 1d arrays argument, set of SH coefficients
    !     cnmH, snmH  ....... 1d arrays argument, set of SH coef of topography, only required for quant = 'geoid'
    !     lat, lon   ... .... 1D array, vector of geodetic latitudes, longitudes
    !     height      ....... 2d array, conatins height of grid H(1:nlat, 1:nlon)
    !     GMegm, aegm ....... GM and scale of EGM. (set GM==1_dp, aegm==1_dp for quant== none and height)
    !
    ! OUTPUT:
    !    f .................. 3D double array, grid file with size(lat) X size (lon) X n (dim of output function)

    use lagInterp

    implicit none

    real(8), allocatable , intent(in):: cnm(:), snm(:), h(:,:)

    real(8), allocatable, intent(in):: lat(:), lam(:)

    real(8), allocatable:: f(:,:,:)

    real(8), intent(in) :: aegm , GMegm

    integer(4):: wid

    real(8), allocatable :: fh(:,:,:), cnmH(:), snmH(:), hi(:)

    integer:: nlat, nlon

    integer:: ii, i, j , done

    real(8) , allocatable:: h0(:)

    integer:: nh, m

    real(8), allocatable :: href(:)

    integer:: ih, i0, i1, i2

    real(8), allocatable :: c_shf(:), s_shf(:)

    real(8) :: west, westm
    
    integer:: nc, nc1

    allocate( c_shf(0:nmax), s_shf(0:nmax))
    c_shf=0; s_shf=0

    West = minval(lam)

    do m=0, nmax

        westm = west * dble(m) * dr
        c_shf(m) = cos (westm)
        s_shf(m) = -sin(westm)

    enddo

    call  FRM_ini(option%nmax, mu, sqr, zsqr, a_n, b_n)

    nlat = size(lat);   nlon = size(lam)

    !- allocate output matrix

    allocate(f(1:nlat, 1:nlon, iprint))

    f=0d0
    
    nc = nint(360.d0/option%step_lon)
    nc1=nc
    
    do

        if (nmax+1<=nc1/2) exit
        nc1=nc1*2

    end do

     if (nc1/nc>100) then

        allocate(cosml(1:nmax+1,nlon), sinml(1:nmax+1,nlon))

        do m = 0, nmax
            cosml(m+1, :) = cos(m *  lam*dr  )
            sinml(m+1, :) = sin(m *  lam*dr  )
        end do
        
    endif

    call h_ref(option%func, nmax, href)

    ih = 4

    allocate(fh(size(href),1:nlon, iprint), h0(size(href)))

    done=0; wid=0;

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(ii, i, j, hi, h0, fh, i0, i1, i2, nh) &
    !$OMP SHARED(h, f, cnm, snm, cnmH, snmH, lat, lam, GMegm, aegm, c_shf, s_shf) &
    !$OMP SHARED(ih, iprint, href, nlon, nlat, wid, done)

    !$OMP DO SCHEDULE(DYNAMIC)
    do ii = 1, nlat
        ! Select height for each parallel iteration
        hi = h(ii, 1:nlon)
        call select_h(hi, nlon, href, ih, h0, nh)

        ! allocate(fh(nh, nlon, iprint))
        fh = 0.0_dp

        do i = 1, nh
            fh(i, 1:nlon, 1:iprint) = shs_parallel(cnm, snm, cnmH, snmH, lat(ii), lam, h0(i), nlon, GMegm, aegm, c_shf, s_shf)
        end do

        if (nh == 1) then
            f(ii, :, :) = fh(1, :, :)
        else
            do i = 1, nlon
                i0 = minloc(abs(hi(i) - h0), dim=1)
                i1 = max(i0 - ih, 1)
                i2 = min(i0 + ih, nh)
                do j = 1, iprint
                    f(ii, i, j) = lagrange_interp_1d(i2 - i1 + 1, h0(i1:i2), fh(i1:i2, i, j), hi(i))
                end do
            end do
        end if

        ! deallocate(fh)

        !$OMP CRITICAL
        done = done + 1
        call progress(wid, done, nlat)
        !$OMP END CRITICAL

    end do
    !$OMP END DO
    !$OMP END PARALLEL

    end function shs_igrid

    !------------------------------------------------------------------

    subroutine h_ref (func , nmax, href)

    integer, intent(in):: nmax
    integer, intent(in) :: func
    real(8), allocatable , intent(out):: href(:)

    real(8):: h, v, height
    integer::     i, j, n

    if (func==1 .or. func == 2 .or. func == 3 .or. func == 8) then

        h=2d-2

    else

        h=1.5d-2

    endif

    h = h / (nmax/2160d0)

    if (h>0.05d0) h=0.05d0

    i=0; n=0

    do

        i=i+1

        v  =  1d0 - (dble(i-1)*h)**4 / 4d0

        height  = 6378137*(1d0/v-1d0)

        if (height > 544000) n=n+1

        if(n>4) exit

    end do

    allocate(href(i))

    do j=1,i

        v  =  1d0 - (dble(j-1)*h)**4 / 4d0

        href(j) = 6378137*(1d0/v-1d0)

    end do

    end subroutine h_ref

    subroutine select_h( h, nlon, href, ih, h0, nh)
    !
    real(8), allocatable :: h(:), href(:)

    real(8):: h0(:)

    integer:: nlon, nh, is, ie, ih

    real(8) :: maxh, minh, range, meanh

    !  if (allocated(h0)) deallocate(h0)

    meanh=sum(h)/nlon

    maxh=maxval(h )

    minh=minval(h )

    range= maxh - minh

    if (range <1d-2 )then

        nh = 1

        !  allocate( h0(1) )

        h0(1) = meanh

    else

        is = minloc( abs(minh - href ), 1  )

        if (minh < href(is) ) is=is-1

        is = is -ih

        is= max(1, is)

        ie = minloc( abs(maxh - href ), 1  )

        ie=ie+ ih

        ie= min(size(href),ie)

        nh= ie-is+1

        ! allocate( h0(nh) )

        h0(1:nh)= href(is:ie)

    endif

    end

    end module shs


