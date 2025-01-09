    module read_files

    use utility

    implicit none

    contains

    !------------------------------------------------------------------

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

    !------------------------------------------------------------------

    

    !------------------------------------------------------------------
    subroutine readHfile(option, lat, lon, h)
    
    implicit none
    type(opt) :: option
   
    real(8), allocatable, intent(out) :: lat(:), lon(:), h(:,:)
    
    integer :: i, j, io_stat, n_points, n_lat, n_lon,k
    real(8), allocatable :: temp_data(:,:)
    
    ! First pass: count the number of points and find unique lat/lon
    open(unit=10, file=option%ifile, status='old', action='read')
    n_points = 0
    do
      read(10, *, iostat=io_stat)
      if (io_stat /= 0) exit
      n_points = n_points + 1
    end do
    close(10)
    
    ! Allocate temporary array
    allocate(temp_data(n_points, 3))
    
    ! Second pass: read the data
    open(unit=10, file=option%ifile, status='old', action='read')
    do i = 1, n_points
      read(10, *) temp_data(i, :)
    end do
    close(10)
    
    option%north = maxval(temp_data(:,1))
    option%south = minval(temp_data(:,1))
    option%east = maxval(temp_data(:,2))
    option%west = minval(temp_data(:,2))
    
    n_lat = 1
    do i = 2, n_points
      if (temp_data(i, 1) /= temp_data(i-1, 1)) n_lat = n_lat + 1
    end do
    
    n_lon=n_points/n_lat
    
    option%step_lat = (option%north-option%south)/(n_lat-1)
    option%step_lon = (option%east-option%west)/(n_lon-1)
    
    option%step_lon = nint(option%step_lon * 3600) / 3600d0
    option%step_lat = nint(option%step_lat * 3600) / 3600d0
    
    ! Allocate output arrays
    allocate( h(n_lat, n_lon))
    
    option%north = nint((option%north + option%step_lat/2) * 3600 ) / 3600d0
    option%west = nint((option%west - option%step_lon/2) * 3600) / 3600d0
    option%east = nint((option%east  + option%step_lon/2) * 3600) / 3600d0
    option%south = nint((option%south  - option%step_lat/2) * 3600) / 3600d0
  
    call makegrid(option, lat, lon)
  
    ! Fill h array
    k=0
    do i = 1, n_lat
        do j=1, n_lon
            k=k+1
            h(i,j) = temp_data(k, 3)
        enddo
        
    end do
    
    ! Deallocate temporary array
    deallocate(temp_data)
  end subroutine readHfile
  
    !------------------------------------------------------------------

    subroutine ReadScatterFile (option, lat, lon, h)

    type(opt) option
    integer:: stat, npoint, i
    real(8):: x, y, z, minh
    real(8), allocatable:: lat(:), lon(:), h(:)
    
    open(10, file= option%ifile , action = 'read' , iostat=stat)

    if (stat /= 0) call stopbox1( trim(option%ifile) // " cannot be opened")

    npoint=0

    do

        read(10,*, iostat=stat) x, y, z
        if (stat<0) exit
        npoint=npoint+1

    enddo
    
    if (allocated(lat)) deallocate(lat)
    if (allocated(lon)) deallocate(lon)
    if (allocated(h)) deallocate(h)

    allocate (lat(npoint), lon(npoint), h(npoint))

    rewind(10)

    do i=1, npoint

        read(10,*, iostat=stat) lat(i), lon(i), h(i)

    enddo

    close(10)

    option%south= minval(lat)

    if (option%south <= -90d0) then

        i= minloc(lat,1)
        call stopbox1('lat of point #' // num2str(dble(i), 'int' , '(i7)') // '<= -90')

    endif

    option%north= maxval(lat)

    if (option%north >= 90d0) then

        i= maxloc(lat,1)
        call stopbox1('lat of point #' // num2str(dble(i), 'int', '(i7)') // '>= 90')

    endif

    option%west= minval(lon)

    if (option%west < -180d0) then

        i= minloc(lon,1)
        call stopbox1('lon of point #' // num2str(dble(i), 'int', '(i7)') // '< -180')

    endif
    option%east= maxval(lon)

    if (option%east > 360d0 ) then

        i= minloc(lon,1)
        call stopbox1('lon of point #' // num2str(dble(i), 'int', '(i7)') // '> -180')

    endif

    minh= minval(h)

    if (minh < 0d0 ) then

        i= minloc(h,1)
        call stopbox1('h of point #' // num2str(dble(i), 'int', '(i7)') // '< 0')

    endif

    end subroutine ReadScatterFile

    !------------------------------------------------------------------

    subroutine makegrid(option, lat, lon)

    real(8), allocatable:: lat(:), lon(:)
    integer:: nlat, nlon, i, centre
    type(opt) option
    
    centre = 1

    nlat = nint((option%north-option%south)/option%step_lat)+(1-centre)
    nlon = nint((option%east-option%west)/option%step_lon)+(1-centre)

    if (allocated(lat)) deallocate(lat)
    if (allocated(lon)) deallocate(lon)
    allocate ( lat(nlat) , lon(nlon) )

    lat = (/ (option%north - option%step_lat/2 - (i-1) * option%step_lat, i=1, nlat) /)

    lon = (/ (option%west +  option%step_lon/2 + (i-1) * option%step_lon, i=1, nlon) /)

    end subroutine makegrid

    !------------------------------------------------------------------

    subroutine stopbox1(msg)
    character(*):: msg
    stop msg
    end

    !------------------------------------------------------------------
    
    SUBROUTINE readEGM(fmt,  file, nmax, cnm, snm, scale, GM, tide)

    ! readEGM read Earth potential model cs file
    !
    !
    ! input
    !
    ! file   ........ string argument, file holding cs coefficient
    ! nmax  ......... integer argument, maximum degree of expansion
    !
    ! output
    ! cnm, snm ........... 2D double array, cs coeffcients, size (0:nmax, 0:lamx)
    ! scale (optional) ... scalar double argument, scale of EGM [m]
    ! GM(optional) ....... scalar double argument, Geocentric gravitational constant of EGM [m3/s2]
    !
    ! acceptable format
    !     1- GFZ ICGEM format
    !     2- free format:
    !
    !  free format:
    !  headerlines (optional)
    !  [n, m,  cnm, snm sigma_c(optional) sigma_s(optional)]
    !
    !  Programmer: M.Goli

    use global_variables, only: dp

    implicit none

    character(*), intent(in):: file, fmt

    integer, intent (in):: nmax

    real(8), allocatable , intent(out):: cnm(:), snm(:)

    real(8), intent(out) :: GM, scale

    character(len=4), intent(out):: tide

    character(:) , allocatable::  imsg

    character(len=100):: Q

    character(len=  6):: dummy

    integer:: nmaxr,  i, j, nall, k

    integer(8) :: ii

    real(8):: c, s

    integer :: istat

    open(file=file, unit=1010,  status='old', iostat=istat, iomsg=imsg)

    if (fmt=='gfc') then

        do

            read(1010, '(A)', err=10) Q

            i=index(Q, 'earth_gravity_constant' )

            if (i /=0 )then

                read(Q(i+22:80) , *, iostat=istat) GM

                if (istat/=0)  call stopbox1('error in readEGM.f90: invalid GM_egm')

            endif

            i=index(Q, 'radius')

            if (i /=0 )then

                read(Q(i+12:80) , *, iostat=istat) scale

                if (istat/=0)  call stopbox1( 'error in readEGM.f90: invalid a_egm')

            endif

            i=index(Q, 'max_degree')

            if (i /=0 )then

                read( Q(i+12:80) , *, iostat=istat) nmaxr

                if (istat/=0) call stopbox1( 'error in readEGM.f90: invalid nmax')


            endif

            i=  index(Q,'tide_system')

            if (i /=0 ) then

                imsg=trim(adjustl(trim( Q(i+11:80))))

                i=  index(Q,'free')
                if (i /= 0 ) tide='free'

                i=  index(Q,'none')
                if (i /= 0 ) tide='none'

                i=  index(Q,'unknown')
                if (i /= 0 ) tide='none'

                i=  index(Q,'zero')
                if (i /= 0 ) tide='zero'

                i=  index(Q,'mean')
                if (i /= 0 ) tide='mean'


            endif

            if ( Q(1:11) == 'end_of_head' ) exit

        end do

        if ( nmaxr < nmax) then

            call stopbox1('insufficient max_degree of model, reduce nmax')

            write (111, '(a25, i5, a3, i5)') adjustl('Error, max_degree of model = ') , nmaxr, ' < ', nmax

            write(111, *) '--- reduce nmax '

            stop

        endif

    else

        gm= 1d0
        scale=1d0
        tide = 'none'

    endif

    if (allocated(cnm)) deallocate(cnm)

    if (allocated(snm)) deallocate(snm)

    nall=(nmax+1)*(nmax+2)/2

    allocate( cnm(nall), snm(nall))

    cnm=0.d0; snm=0.d0

    k=0;

    if (fmt=='gfc') then

        do

            READ(1010,*, err=10, end=40) dummy, i, j, c, s

            if(i>nmax) cycle

            if (dummy/='gfc')cycle

            ii=nm2row(i, j)

            cnm(ii)=c

            snm(ii)=s

            if (i==nmax)then

                if(j==nmax) exit

            endif

        enddo

    else

        do

            READ(1010,*, err=10, end=40) i, j, c, s

            if(i>nmax) cycle

            ii=nm2row( i, j)

            cnm(ii)=c

            snm(ii)=s

            if (i==nmax .and. j==nmax) exit

        enddo

    endif

40  if (i<nmax) then

        print*,  'error in -readEGM:: expected nmax=' , nmax
        write(111, *)  'error in -readEGM:: expected nmax=' , nmax

        print*,  '....................: found  nmax=', i
        write(111,*)  '....................: found  nmax=', i

        stop

    endif

    close(1010)

    return

10  call stopbox1( 'error occures in reading EGMfile= ' // file )
    write(111, *) 'error occures in reading EGMfile=' , trim(file)

    stop

    END subroutine readEGM
    
    end module read_files
