    module lagInterp

    implicit none

    ! module for 1d/3d Lagrange interpolation

    contains

    !--------------------------------------------------------------------

    function lagrange_basis_function_1d ( mx, xd, i, xi) result( yi )

    ! basic lagrange function
    ! input:
    ! mx, ... integer, scalar
    ! xd  ... xd(mx) double array
    ! i   ... integer, scalar
    ! xi  ... double, scalar
    ! output
    ! yi .... double, scalar
    !
    ! programmer:  John Burkardt

    implicit none

    integer , intent (in)::  mx, i
    real (8), intent (in):: xd(mx)
    real (8), intent (in):: xi
    real (8)  :: yi
    integer:: j

    yi = 1.0D+00

    if ( xi /= xd(i) ) then

        do j = 1, i-1

            yi = yi * ( xi - xd(j) ) / ( xd(i) - xd(j) )

        end do
        
        do j = i+1, mx

            yi = yi * ( xi - xd(j) ) / ( xd(i) - xd(j) )

        end do

    endif

    end function lagrange_basis_function_1d

    !--------------------------------------------------------------------

    function lagrange_interp_3d ( mx, my, mz, x, y, z, f, xi) result(fi)

    ! 3D lagrange interpolation
    ! input:
    ! mx, my, mz         ..... integer, scalar, dims of function in x,y,z
    ! x(mx),y(my),z(mz)  ..... 1d array, double, nodes in x,y,z
    ! f(mx,my,mz)        ..... 3d array, double, data at x,y,z
    ! xi(3)              ..... 1d array, double, position of qurey point,
    !                           x(1)=x, x(2)=y, x(3)= z
    ! output
    ! fi                 ..... double, scalar, interpolated value at qurey point
    !
    ! programmer: John Burkardt
    !    modifed by Mehdi Goli 2d => 3d

    implicit none

    integer, intent(in):: mx, my, mz

    real(8) , intent (in) :: x(mx), y(my), z(mz), f(mx, my, mz)

    real(8), intent(in) :: xi(3)

    real(8) :: fi

    integer:: i, j, k

    real(8):: lx(mx), ly(my), lz(mz), lxi, lyj


    fi =0d0

    do i = 1, mx

        lx(i) = lagrange_basis_function_1d ( mx, x, i, xi(1))

    enddo

    do j = 1, my

        ly(j) = lagrange_basis_function_1d ( my, y, j, xi(2))

    enddo

    do k= 1, mz

        lz(k) = lagrange_basis_function_1d ( mz, z, k, xi(3))

    enddo

    do i = 1, mx

        lxi=lx(i)

        do j = 1, my

            lyj=ly(j)

            do k= 1, mz

                fi  = fi  + f(i, j, k) * lxi * lyj * lz(k)

            end do

        end do

    end do

    return

    end function lagrange_interp_3d

    !---------------------------------------------------------------------------

    function lagrange_interp_2d (nf, mx, my, x, y,  f, xi) result(fi)

    ! 3D lagrange interpolation
    ! input:
    ! mx, my, mz         ..... integer, scalar, dims of function in x,y,z
    ! x(mx),y(my)   ..... 1d array, double, nodes in x,y,z
    ! f(mx,my )        ..... 3d array, double, data at x,y,z
    ! xi(2)              ..... 1d array, double, position of qurey point,
    !                           x(1)=x, x(2)=y, x(3)= z
    ! output
    ! fi                 ..... double, scalar, interpolated value at qurey point
    !
    ! programmer: John Burkardt
    !    modifed by Mehdi Goli 2d => 3d

    implicit none

    integer, intent(in):: mx, my, nf

    real(8) , intent (in) :: x(mx), y(my) , f(mx, my, 0:nf)

    real(8), intent(in) :: xi(2)

    real(8) :: fi( 0:nf)

    integer::  i, j

    real(8):: lx(mx), ly(my) , lxi , lyj


    fi =0d0

    do i = 1, mx

        lx(i) = lagrange_basis_function_1d ( mx, x, i, xi(1))

    enddo

    do j = 1, my

        ly(j) = lagrange_basis_function_1d ( my, y, j, xi(2))

    enddo

    do i = 1, mx

        lxi=lx(i)

        do j = 1, my

            lyj=ly(j)


            fi(:)  = fi(:)  + f(i, j, : ) * lxi * lyj


        end do

    end do

    return

    end function lagrange_interp_2d

    !---------------------------------------------------------------------------

    function hermit_m1 (n, u,x,y,dy) result(H)  ! Hermite interpolation (Lagrange)

    ! u: discrete data points;
    ! vector x: [x_1,...,x_n]
    ! vector y: [y_1,...,y_n]
    ! vector dy: [y'_1,...,y'_n]

    integer, intent(in) :: n

    real(8) ,intent(in) :: u, x(n), y(n), dy(n)

    real(8) :: H

    real(8) :: li  ! Lagrange basis polynomials

    real(8) :: a ! basis polynomials alpha(x)

    real(8) :: b ! basis polynomials beta(x)

    integer:: i, j

    real(8):: dl, l2

    H = 0d0

    do i=1, n

        dl=0d0;       ! derivative of Lagrange basis

        li=1d0

        do j=1, n

            if (j == i) cycle

            dl=dl+1/(x(i)-x(j));

            li=li*(u-x(j))/(x(i)-x(j));

        end do

        l2 = li*Li;

        b = (u-x(i)) * l2;            ! basis polynomial alpha(x)

        a = (1d0-2d0*(u-x(i))*dl) *l2;  ! basis polynomial beta(x)

        H = H + a * y(i) + b * dy(i) ;  ! Hermite polynomial H(x)

    end do

    end function hermit_m1

    !---------------------------------------------------------------------------

    function lagrange_interp_1d ( mx, x, f, xi) result (fi)

    ! 3D lagrange interpolation
    ! input:
    ! mx          ..... integer, scalar, dims of function in x,y,z
    ! x(mx)       ..... 1d array, double, nodes in x,y,z
    ! f(mx)       ..... 3d array, double, data at x,y,z
    ! xi          ..... scalar, double, position of qurey point,
    !
    ! output
    ! fi                 ..... double, scalar, interpolated value at qurey point
    !
    ! programmer: John Burkardt
    !    modifed by Mehdi Goli 1d

    implicit none

    integer, intent(in):: mx

    real(8) , intent (in) :: x(mx) , f(mx )

    real(8), intent(in) :: xi

    real(8) :: fi

    integer:: l, i

    real(8):: lx(mx) , lxi

    l = 0

    fi = 0d0

    do i = 1, mx

        lxi=lagrange_basis_function_1d ( mx, x, i, xi)

        fi  = fi  + f(i) * lxi

    end do

    return

    end function lagrange_interp_1d

    end module lagInterp