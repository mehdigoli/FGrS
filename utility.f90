    module utility

    use global_variables

    implicit none

    contains

    !------------------------------------------------------------------

    function num2str (f, typ, fmt)

    ! num2str converts number to string with fmt format
    ! input
    ! f    ...... double argument
    ! typ  ...... character argument
    !           - 'int' integer intput
    !           - 'real' real input
    ! fmt    .... optional character argument, in the form of FORTRAN format style
    ! output
    ! num2str ... charater argument
    !programmer: M. Goli
    ! example:
    ! num2str( 1000.d0, '(f4.0)')

    real(8) , intent (in):: f

    character (len= *) , intent(in)  :: typ

    character (len= *) , intent(in) , optional :: fmt

    character(len=30) :: num2str

    if (present (fmt)) then

        if (typ == 'int') then
            write( num2str, fmt ) int(f)
        else
            write( num2str, fmt ) f
        endif

    else

        if (typ == 'int') then
            write( num2str, * ) int(f)
        else
            write( num2str, * ) f
        endif

    endif

    end function num2str


    subroutine progress(wid, done, all)
    ! Progress creates a progress bar

    ! Input:
    !     wid  - integer, percentage of completed iterations (in/out)
    !     done - integer, number of completed iterations
    !     all  - integer, total number of iterations

    use iso_fortran_env, only: output_unit
    use iso_c_binding, only: c_backspace, C_CARRIAGE_RETURN

    implicit none

    integer, intent(in) :: done, all
    integer(4), intent(inout) :: wid  ! in/out for tracking percentage
    integer(4) :: progress_percent    ! calculated percentage of progress
    integer :: k                      ! loop index for updating the bar
    character(len=13) :: bar          ! display bar for progress
    character(1) :: newline_char      ! newline character for output


    newline_char = C_CARRIAGE_RETURN
    ! Initialize
    bar = 'progress ???%'     ! Template for progress bar
    progress_percent = done * 100 / all  ! Calculate progress percentage
  

    ! Update percentage in the bar
    write(bar(10:12), fmt='(I0)') progress_percent

    ! If progress hasn't advanced, return early
    if (progress_percent <= wid) return

    ! Update the saved progress percentage
    wid = progress_percent

    ! Overwrite the previous progress bar
    write(output_unit, '(A)', advance='no') repeat(c_backspace, len_trim(bar))

    ! If not complete, display updated progress bar
    if (progress_percent /= 100) then
        write(output_unit, '(A1,A13)', advance='no') newline_char, bar
    else
        ! On completion, display final message
        write(output_unit, '(A1,A)', advance='no') newline_char, 'progress completed'
    end if

    flush(output_unit)  ! Ensure immediate output

    return
    end subroutine progress



    !------------------------------------------------------------------

    subroutine Tic(start, rate)

    ! subroutine to set srart of time
    ! input :(-)
    ! output :
    ! start ... integer, start of timer
    ! rate .... integer, rate of timer (miliseconds)
    ! output (-)

    implicit none

    integer(4), intent(out) :: start, rate

    call system_clock(count_rate=rate)

    call system_clock(start)

    end subroutine Tic

    !------------------------------------------------------------------

    function toc (start, rate) result(txt_time)

    ! subroutine to measure and display (on screen) the elapsed time
    ! since the stopwatch timer started.
    ! input:
    ! start ... integer, start of timer
    ! rate .... integer, rate of time (miliseconds)
    ! output (-)
    ! programmer: M.Goli

    implicit none

    integer(4), intent(in) :: start, rate

    integer(4):: finish

    real(4):: time

    character(80) :: txt_time

    call system_clock(finish)

    time=dble(finish-start)/rate

    call sec2hms(time, txt_time )

    end function toc

    !------------------------------------------------------------------

    SUBROUTINE sec2hms(T, text_out)

    !reference https://rosettacode.org/wiki/Convert_seconds_to_compound_duration#Fortran

    implicit none

    real(4), intent(in):: T      !The time, in seconds. Positive only, please.

    INTEGER NTYPES  !How many types of time?

    PARAMETER (NTYPES = 4)!This should do.

    INTEGER USIZE(NTYPES)!Size of the time unit.

    CHARACTER(len=3):: UNAME(NTYPES)!Name of the time unit.

    PARAMETER (USIZE = (/24*60*60, 60*60,   60,    1/)) !The compiler does some arithmetic.

    PARAMETER (UNAME = (/"day","hr ","min","sec"/)) !Approved names, with trailing spaces.

    CHARACTER(len=45) :: TEXT

    CHARACTER(len=80) :: text_out

    INTEGER I,L,N,S

    S = int(T) !A copy I can mess with.

    L = 0 !No text has been generated.

    do I = 1,NTYPES !Step through the types to do so.

        if (I/=NTYPES) then

            N = S/USIZE(I)

            IF (N.GT.0) THEN

                S = S - N*USIZE(I)

                IF (L.GT.0) THEN  !Is this the first text to be rolled?

                    L = L + 2

                    TEXT(L - 1:L) = ", "

                END IF   !Now ready for this count.

                WRITE (TEXT(L + 1:),'(I0,1X,A)') N,UNAME(I) !Place, with the unit name.

                L = LEN_TRIM(TEXT)  !Find the last non-blank resulting.

            END IF   !Since I'm not keeping track.

        else

            N = S/USIZE(I) !Largest first.

            S = S - N*USIZE(I)  !Yes! Remove its contribution.
            !No.

            if (L>0) then

                L = L + 2 ; TEXT(L - 1:L) = ","  !Cough forth some punctuation.

            endif

            WRITE (TEXT(L + 1:), '(f6.3,1x, A)' ) N+(T-int(T))  ,UNAME(I) !Place, with the unit name.

        endif

    END do   !On to the next unit.

    !Cast forth the result.
    WRITE (text_out,'(a)') 'computation time= ' // TEXT
    WRITE (111,'(a)') 'computation time= ' // TEXT

    END SUBROUTINE sec2hms

    !------------------------------------------------------------------

    subroutine print_error(msg, msg1, msg2)

    character(len = *) , intent(in) :: msg
    character(len = *) , intent(in), optional :: msg1
    character(len = *) , intent(in), optional :: msg2

    print*, msg
    write(111,*) msg

    if (present(msg1) ) then
        print*, msg1
        write(111,*) msg1
    endif

    if (present(msg2) ) then
        print*, msg2
        write(111,*) msg2
    endif

    stop

    end subroutine print_error

    !------------------------------------------------------------------

    subroutine tip_func (func)

    ! tip_function returns the number of functions required for different functions, and
    ! the number of output
    ! input:
    ! func ...... integer id of gravity functional
    ! output
    ! col ........ integer, number of reqiured function for computation of func
    ! iprintd .... integer, number of finally outputs
    ! tip     .... integer, flag shows that the category of func

    integer :: func

    iprint=1

    select case (func)

    case(1) !W,V

        col=1; tip = 1; iprint=2

    case (4) ! geoid

        col=3; tip = 2; iprint=1

    case(2:3, 5:8)!T, Z! dgsa, gasa, H, none

        col=1; tip = 3; iprint=1

    case(9)!    DoV

        col=2; tip = 4; iprint=2

    case(10,13) !norm(dg) norm(ga)

        col=4; tip = 7; iprint=1

    case( 11:12, 14,15 )   ! vector of gravity, gravitional,  dW, dT

        col=3; tip = 5; iprint=3

    case(16) !second derivative

        col = 9; tip = 6; iprint= 6

    case(17)

        col = 9; tip = 6; iprint= 9

    end select

    end subroutine tip_func

    !------------------------------------------------------------------

    subroutine stop_error(msg)

    character(*):: msg

    write(*, '(a)') ADJUSTL(trim(msg))
    write(111,  '(a)') ADJUSTL(trim(msg))

    stop

    end
    end module utility
