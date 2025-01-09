    module FFTW3

    use, intrinsic :: iso_c_binding

    integer, parameter :: C_FFTW_R2R_KIND = C_INT32_T

    integer(C_INT), parameter :: FFTW_R2HC = 0
    integer(C_INT), parameter :: FFTW_HC2R = 1
    integer(C_INT), parameter :: FFTW_DHT = 2
    integer(C_INT), parameter :: FFTW_REDFT00 = 3
    integer(C_INT), parameter :: FFTW_REDFT01 = 4
    integer(C_INT), parameter :: FFTW_REDFT10 = 5
    integer(C_INT), parameter :: FFTW_REDFT11 = 6
    integer(C_INT), parameter :: FFTW_RODFT00 = 7
    integer(C_INT), parameter :: FFTW_RODFT01 = 8
    integer(C_INT), parameter :: FFTW_RODFT10 = 9
    integer(C_INT), parameter :: FFTW_RODFT11 = 10
    integer(C_INT), parameter :: FFTW_FORWARD = -1
    integer(C_INT), parameter :: FFTW_BACKWARD = +1
    integer(C_INT), parameter :: FFTW_MEASURE = 0
    integer(C_INT), parameter :: FFTW_DESTROY_INPUT = 1
    integer(C_INT), parameter :: FFTW_UNALIGNED = 2
    integer(C_INT), parameter :: FFTW_CONSERVE_MEMORY = 4
    integer(C_INT), parameter :: FFTW_EXHAUSTIVE = 8
    integer(C_INT), parameter :: FFTW_PRESERVE_INPUT = 16
    integer(C_INT), parameter :: FFTW_PATIENT = 32
    integer(C_INT), parameter :: FFTW_ESTIMATE = 64
    integer(C_INT), parameter :: FFTW_WISdoM_ONLY = 2097152
    integer(C_INT), parameter :: FFTW_ESTIMATE_PATIENT = 128
    integer(C_INT), parameter :: FFTW_BELIEVE_PCOST = 256
    integer(C_INT), parameter :: FFTW_NO_DFT_R2HC = 512
    integer(C_INT), parameter :: FFTW_NO_NONTHREADED = 1024
    integer(C_INT), parameter :: FFTW_NO_BUFFERING = 2048
    integer(C_INT), parameter :: FFTW_NO_INDIRECT_OP = 4096
    integer(C_INT), parameter :: FFTW_ALLOW_LARGE_GENERIC = 8192
    integer(C_INT), parameter :: FFTW_NO_RANK_SPLITS = 16384
    integer(C_INT), parameter :: FFTW_NO_VRANK_SPLITS = 32768
    integer(C_INT), parameter :: FFTW_NO_VRECURSE = 65536
    integer(C_INT), parameter :: FFTW_NO_SIMD = 131072
    integer(C_INT), parameter :: FFTW_NO_SLOW = 262144
    integer(C_INT), parameter :: FFTW_NO_FIXED_RADIX_LARGE_N = 524288
    integer(C_INT), parameter :: FFTW_ALLOW_PRUNING = 1048576

    type, bind(C) :: fftw_iodim
        integer(C_INT) n, is, os
    end type fftw_iodim
    type, bind(C) :: fftw_iodim64
        integer(C_INTPTR_T) n, is, os
    end type fftw_iodim64

    interface

    type(C_PTR) function fftw_plan_dft_c2r_1d(n,in,out,flags) bind(C, name='fftw_plan_dft_c2r_1d')
    import
    integer(C_INT), value :: n
    complex(C_doUBLE_COMPLEX), dimension(*), intent(out) :: in
    real(C_doUBLE), dimension(*), intent(out) :: out
    integer(C_INT), value :: flags
    end function fftw_plan_dft_c2r_1d

    subroutine fftw_execute_dft_c2r(p,in,out) bind(C, name='fftw_execute_dft_c2r')
    import
    type(C_PTR), value :: p
    complex(C_doUBLE_COMPLEX), dimension(*), intent(inout) :: in
    real(C_doUBLE), dimension(*), intent(out) :: out
    end subroutine fftw_execute_dft_c2r

    subroutine fftw_destroy_plan(p) bind(C, name='fftw_destroy_plan')
    import
    type(C_PTR), value :: p
    end subroutine fftw_destroy_plan

    end interface

    end module FFTW3
