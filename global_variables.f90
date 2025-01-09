    module global_variables

    ! module defines constants
    use, intrinsic :: iso_c_binding

    integer, parameter :: sp = selected_real_kind(6, 33)
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(8), parameter:: pi=3.14159265358979323846264338327950288d0
    real(8), parameter:: dr=pi/180.0d0 !degree2radius
    real(8), parameter:: rd=180.0d0/pi !radius2degree
    real(8), parameter:: Rmean=6371008.0d0
    real(8), parameter:: rho=2670.d0
    real(8), parameter:: G =6.67428d-11;         !gravitational constant [m^3 /(kg s^2)]

    Type opt

        character (128):: egmfile, normal, ofile, ifile, hnm
        real(8):: North, south, west, east, step_lat, step_lon, h0
        integer*4:: mod, nmin, nmax, nthread, func, inormal

    end type opt

    type (opt) :: option
    
    integer:: col, iprint, tip, nmin, nmax 
    !
    character(len=50), dimension(21):: name_func, name_normal
    !
    !character(len=1000) :: all_text
    !
    character (len=1) :: degree

    real(8), allocatable::  sqr(:), zsqr(:), a_n(:), b_n(:), mu(:) , cosml(:,:), sinml(:,:)

    end module   global_variables