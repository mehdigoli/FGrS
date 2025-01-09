    module options

    use global_variables
    use utility

    contains

    subroutine print_logo(fid)
    integer:: fid
    write(fid,*) '                        _____ ____      ____   '
    write(fid,*) '                       ||  ___/ ___|_ __/ ___|  '
    write(fid,*) '                       || |_ | |  _| ''__\___ \  '
    write(fid,*) '                       ||  _|| |_| | |   ___) | '
    write(fid,*) '                       ||_|   \____|_|  |____/  '
    write(fid,*) '                       -------------------------'
    end subroutine print_logo

    !------------------------------------------------------------------------

    subroutine write_option(fid)
    integer:: fid
    character (len=8)  cdate
    character (len=10) ctime
    character (len=5) czone
    integer(4) ival(8)
    integer(4) yr,mon,day,hr,min,sec,ms

    call print_logo(fid)
    write(fid, '(20x,a)')  'FGrS: Fast Gravimetric Synthesis V1.0'
    call date_and_time(cdate,ctime,czone,ival)

    read (cdate(1:4),*) yr
    read (cdate(5:6),*) mon
    read (cdate(7:8),*) day
    read (ctime(1:2),*) hr
    read (ctime(3:4),*) min
    read (ctime(5:6),*) sec
    read (ctime(8:10),*) ms

    ! Get current date and time
    degree = char(176)
    ! Print date and time
    write(fid, '(a, 10x, i0,a,i0,a,i0)') "Date:                    ", yr, "-", mon, "-", day
    write(fid, '(a, 10x, i0,a,i0,a,i0)') "Time:                    ", hr, ":", min, ":", sec
    write(fid, '(a, 10x, a)') "Functional oputput:      ", name_func(option%func)(5:)
    write(fid, '(a, 10x, a)') "Reference model:         ", name_normal(option%inormal)(6:)
    write(fid, '(a, 10x, a)') "Geopotential model file: ", trim(option%egmfile)
    write(fid, '(a, 10x, i0)') "Minimum used degree:     ", option%nmin
    write(fid, '(a, 10x, i0)') "Maximum used degree:     ", option%nmax
    write(fid, '(a, 10x, i0)') "Number of threads:       ", option%nthread

    end subroutine write_option

    !=====================================================================

    subroutine check_option(fid)

    integer:: fid
    character(len=80)::txt
    LOGICAL :: file_exists
    txt='none'

    if (option%nmax < 0) txt="nmax<0"

    if (option%nmin < 0) txt= "nmin<0"

    if (option%nmin > option%nmax) txt="nmin>nmax"

    if (option%inormal < 1 .or. option%inormal >11) then
        txt= "select a valid normal field"
    endif

    if (option%mod ==1 ) then

        if (option%south <= -90 .OR. option%south >= 90) &
            txt= "south <= -90 .OR. south >= 90"

        if (option%north <= -90 .OR. option%north >= 90) &
            txt= "north <= -90 .OR. north >= 90"

        if (option%south > option%south) txt="south > north"

        if (option%west > 360 .OR. option%west < -360) &
            txt="west > 360 .OR. west < -360"

        if (option%east > 360 .OR. option%east < -360) &
            txt= "east > 360 .OR. east < -360"

        if (option%west > option%east) txt= "west > east"

        ! option%step_lat=option%step_lat/3600.d0
        if (option%step_lat < 0 ) txt="step lat < 0"

        ! option%step_lon =option%step_lon/3600.d0
        if (option%step_lon < 0 ) txt = "step lon < 0"

        if (option%h0 < 0 ) txt="h < 0"

    endif

    if (option%nthread < 0 ) txt="nthread < 0"


    ! Use INQUIRE to check if the egmfile exists

    inquire(FILE=option%egmfile, EXIST=file_exists)

    if (not(file_exists)) then
        txt= "EGM file does not exist."//trim(option%egmfile)
    end if

    ! Use INQUIRE to check if the inputfile exists

    if (option%mod /= 1 ) then
        inquire(FILE=option%ifile, EXIST=file_exists)

        if (not(file_exists)) then
            txt= "Input file does not exist: "//trim(option%ifile);
        end if
    endif

    
    if (txt/='none') then

        write(fid,'(a)') trim(txt)
        call stop_error(trim(txt))

    endif
    end

    end module options
