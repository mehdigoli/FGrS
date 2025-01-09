    module normal_field_

    use global_variables
    use read_files

    implicit none

    real(8):: ae, GMe, e2, Omega, flat, be, j2, omega2, EE, ae2, be2

    real(8), allocatable:: klm(:)

    contains

    !------------------------------------------------------------------

    subroutine geod2sphe(lat, h, r, theta)

    !   geod2sphe converts geodetic latitue and geodetic height to
    !   spherical geocentric radius and geocenteric co-latitude
    !
    !   input:
    !    lat ...... scalar double argument, geodetic latitude [deg],
    !    h   ...... scalar double argument, geodetic height   [m],
    !    ae  ...... scalar double argument, semi-major of normal ellipsoid [m],
    !    e2  ...... scalar double argument, squared eccentricity of normal ellipsoid [-],
    !
    ! output
    !  r      ....... scalar double argument, geocentric radius [m]
    !  theta  ....... scalar double argument, co geocentric latitude [deg]
    ! programmer: M.Goli

    use global_variables, only: dr, pi, dp

    implicit none

    real(8), intent(in):: lat, h

    real(8), intent(out):: r, theta

    real(8)::  Ne

    real(8) :: x, z

    Ne=ae/sqrt(1d0-e2*sin(lat*dr)**2)

    x=(Ne+h)*cos(lat*dr)

    z=(Ne*(1d0-e2)+h)*sin(lat*dr)

    r=sqrt(x*x+z*z)

    theta=acos(z/r)*180.0d0/pi !in degree

    end subroutine geod2sphe

    !------------------------------------------------------------------

    subroutine GRS(name)

    !      subroutine GRS returns geometrical and physical parameters of GRS80.
    ! output
    ! ae ........... scalar double argument, semi-major of GRS80 [m]
    ! e2 ........... scalar double argument, squared eccentricity of GRS80
    ! omega ........ scalar double argument, angular velocity of rotation of GRS80 [1/s]
    ! flat ......... scalar double argument, flattening form factor of GRS80
    ! j2 ........... scalar double argument, dynamical form factor of GRS80
    ! GMe .......... scalar double argument, Geocentric gravitational constant of GRS80 [m3/s2]
    ! programmer: M.Goli

    use global_variables, only : dp

    implicit none

    integer, intent(in):: name

    !  name_normal(1)=" 1   GRS80"
    !  name_normal(2)=" 2   WGS84"
    !  name_normal(3)=" 3   EGM2008"
    !  name_normal(4)=" 4   NIMA96"
    !  name_normal(5)=" 5   WGS72"
    !  name_normal(6)=" 6   GRS67"
    !  name_normal(7)=" 7   GRIM5"
    !  name_normal(8)=" 8   JGM3"
    !  name_normal(9)=" 9   MOON"
    !  name_normal(10)="10  MARS"
    !  name_normal(11)="11  VENUS"
    
    select case (name)
        
    case (2) !'WGS84' .or. 'wgs84'
        ae=           6378137.d0
        GMe=  3.986004418000000d+14
        J2=  1.082629821257275d-03
        Omega=  7.292115000000000d-05
        flat=  2.982572235630000d+02
    case (1) ! 'GRS80' .or. 'grs80'
        ae=           6378137.d0
        GMe=  3.986005000000000d+14
        J2=  1.082629999944205d-03
        Omega=  7.292115000000000d-05
        flat=  2.982572221010000d+02
    case (7) !'GRIM5' .or. 'grim5'
        ae=           6378136.46d0
        GMe=  3.986004415000000d+14
        J2=  1.082626920081630d-03
        Omega=  7.292115000000000d-05
        flat=  2.982576500000000d+02
    case (3) !'EGM2008' .or. 'egm2008'
        ae=           6378136.580d0
        GMe=  3.986004415000000d+14
        J2=  1.082626585714383d-03
        Omega=  7.292115000000000d-05
        flat=  2.982576860000000d+02
    case (6) !'GRS67' .or. 'grs67'
        ae=           6378160.d000
        GMe=  3.986030000000000d+14
        J2=  1.082699999922893d-03
        Omega=  7.292115146700000d-05
        flat=  2.982471674270000d+02
    case (8) !'JGM3' .or. 'jgm3'
        ae=           6378136.3d0
        GMe=  3.986004415000000d+14
        J2=  1.082631872225640d-03
        Omega=  7.292115000000000d-05
        flat=  2.982570000000000d+02
    case (5) !'WGS72' .or. 'wgs72'
        ae=           6378135d0
        GMe=  3.986005000000000d+14
        J2=  1.082632743299607d-03
        Omega=  7.292115000000000d-05
        flat=  2.982570000000000d+02
    case (4) !'NIMA96' .or. 'nima96'
        ae=           6378136.d0
        GMe=  3.986004415000000d+14
        J2=  1.082636412777191d-03
        Omega=  7.292115000000000d-05
        flat=  2.982564150990000d+02
    case (9) !'MOON' .or. 'moon'
        ae=           1738140.d000
        GMe=  4.902801076000000d+12
        J2=  2.032013788211070d-04
        Omega=  2.661621000000000d-06
        flat=  3.240000000000000d+03
    case (10) !'MARS' .or. 'mars'
        ae=           3397000d000
        GMe=  4.282836977393900d+13
        J2=  1.955402378489285d-03
        Omega=  7.088220000000000d-05
        flat=  1.911800000000000d+02
    case (11) !'VENUS' .or. 'venus'
        ae=           6051000d0
        GMe=  3.248585897260000d+14
        J2=  4.403426442298153d-06
        Omega=  2.992600000000000d-07
        flat=  1.507000000000000d+05
    end select

!    other parameters
       
    flat=1d0/flat
    be=ae-ae*flat
    e2=(ae*ae-be*be)/ae/ae
    omega2 = omega *omega
    ae2= ae*ae
    be2= be*be
    EE=   (ae2-be2 )

    end subroutine GRS

    !------------------------------------------------------------------

    subroutine normal_field_par(nmax, ellip)

    ! normal_filed_par returnes 4pi fully normalized
    ! spherical harmonics of the the Earth's normal potential
    !
    ! input
    !
    ! nmax ......... integer argument, maximum degree of expansion
    ! ellip ........ string argument, define the normal filed (grs80 or wgs84)
    !
    !output
    ! ae ........... scalar double argument, semi-major of normal ellipsoid [m]
    ! e2 ........... scalar double argument, squared eccentricity of normal ellipsoid
    ! GMe .......... scalar double argument, Geocentric gravitational constant
    !                of normal ellipsoid [m3/s2]
    ! kml .......... 1D double array, spherical harmonics of the the Earth's normal potential
    !                size 0:nmax
    !programmer: M. Goli

    use global_variables, only : dp

    implicit none

    integer, intent(in)  :: nmax

    integer, intent(in) :: ellip

    real(8):: DJ2N

    integer:: n, L

    klm=0.0d0

    call GRS(ellip)

    do n=0,nmax,2

        L=n/2

        DJ2N=(-1.0d0)**(L+1)* 3.0d0 * e2 ** L /(n+1.0d0)/ &

            (n+3.0d0)*(1.0d0 - L + 5.0d0* L * j2/ e2)

        klm(n)=-DJ2N * DSQRT(1.0d0/(2.0d0*dble(n)+1.0d0))

    end do

    end subroutine normal_field_par

    !------------------------------------------------------------------

    function dgamdh(gamma, slat) result( dgdh )

    real(8):: slat, dgdh, gamma

    dgdh= -2d0*gamma/ae*(1d0 + flat + omega2*ae2*be/GMe - 2d0* flat * slat*slat)

    end

    !------------------------------------------------------------------

    subroutine normal_gravity_p (h, lat, gamma, dgdh)

    ! normal_gravity_p computes the normal gravity at point (lat,h)

    ! input :
    ! GRS ..... character , the name of normal ellipsoid
    ! h   ..... double scalar, geodetic height [m]
    ! lat ..... double scalar, geodetic latitude [deg]
    ! output
    ! gamma ... double scalar, normal gravity at altitude h
    ! dgdh  ... double scalar, optional argument, vertical gradient of gamma at altitude h [mGal/m]
    !programmer : M. Goli

    use global_variables , only :  dp, pi, dr

    implicit none

    real(8), intent(in) :: h, lat

    real(8), intent(out):: gamma

    real(8), intent(out), optional:: dgdh

    real(8) :: slat, clat

    real(8) :: N

    real(8) :: x, z, r

    real(8) ::    E, u, u2pE2

    real(8) :: w, beta, cbeta, sbeta, q, q0, qp

    real(8) :: gamma_beta, gamma_u , uu

    ! reference : World Geodetic System 1984
    ! Its Definition and Relationships with Local Geodetic Systems
    ! NIMA report pp. 4-2 : 4-4

    slat= sin(lat*dr)

    clat= cos(lat*dr)

    N= ae / sqrt( 1.d0 - e2*slat**2 )

    x= ( N + h ) * clat

    z= ( N*(1-e2) +h) *slat

    E=sqrt(EE)

    r=  x*x+ z*z -E*E

    u =  sqrt( r * 0.5d0 * ( 1.d0 + sqrt ( 1.d0 + 4.d0 * EE * Z*Z / r**2 )))

    uu=u*u

    u2pE2 = uu + EE

    beta = atan ( z * sqrt (u2pE2 ) / u / x )

    sbeta= sin(beta)

    cbeta= cos(beta)

    q= ( ( 1.d0 +3.d0 *uu/EE) * atan (E/u) -3.d0 * u / E) /2.d0

    q0= ( ( 1.d0 +3.d0 *be2/EE) * atan (E/be) -3.d0 * be / E) /2.d0

    qp= 3.d0 * ( 1.d0 + uu/EE) * ( 1.d0- u/E * atan(E/u)) - 1.d0

    w = sqrt ( (uu +  EE *sbeta *sbeta ) / u2pE2)

    gamma_u = - GMe / u2pE2 - omega2 * ae2 * E / u2pE2 * qp /q0 * &
        ( sbeta *sbeta/ 2.d0 -1.d0/6.d0) + &
        omega2 * u * cbeta * cbeta
    gamma_u = gamma_u /w

    gamma_beta= omega2 * ae2 / sqrt(u2pE2) * q/q0 * sbeta * cbeta - &
        omega2 * sqrt (u2pE2) * sbeta * cbeta

    gamma_beta= gamma_beta / w

    gamma = sqrt ( gamma_beta **2 + gamma_u**2)

    if (present(dgdh)) dgdh= -2d0*gamma/ae*(1d0 + flat + omega2*ae2*be/GMe - 2d0* flat * slat*slat)

    end subroutine normal_gravity_p

    !------------------------------------------------------------------

    subroutine sub_normal( option, aegm, gmegm, cnm, klm)

    type(opt) option
    real(8), allocatable:: cnm(:), klm(:)
    real(8):: aegm, gmegm
    integer:: n, k

    if (allocated(klm)) deallocate (klm)

    allocate(klm(0:option%nmax))

    call normal_field_par(option%nmax, option%inormal)

    if ( option%func == 1 .or. option%func == 7  .or. option%func == 8  .or. option%func == 10 &
        .or. option%func == 11 .or. option%func == 13 ) return

    do n=0, 20, 2

        k=nm2row(n,0)

        cnm(k)=cnm(k) -klm(n)*GMe/GMegm*(ae/aegm)**n

    enddo

    end subroutine sub_normal

    !------------------------------------------------------------------

    end module normal_field_