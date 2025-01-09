    program FGrS

    use global_variables
    use shs
    use utility
    use read_files
    use normal_field_
    use options

    real(8):: aegm, gmegm

    character(len=4) :: egmtide

    real(8), allocatable:: cnm(:), snm(:)

    real(8), allocatable:: cnmH(:), snmH(:)

    real(8),allocatable:: lat(:), lon(:) , h(:), hg(:,:)

    real(8), allocatable:: f(:,:), f3(:,:,:)

    character(80) :: tempi

    integer  :: start, rate, i, j, wid

    character(len=20):: name_func1(17), mod
    
    logical:: file_exists

    !    character(len=3) :: fmt

    name_normal(1)=" 1   GRS80"
    name_normal(2)=" 2   WGS84"
    name_normal(3)=" 3   EGM2008"
    name_normal(4)=" 4   NIMA96"
    name_normal(5)=" 5   WGS72"
    name_normal(6)=" 6   GRS67"
    name_normal(7)=" 7   GRIM5"
    name_normal(8)=" 8   JGM3"
    name_normal(9)=" 9   MOON"
    name_normal(10)="10  MARS"
    name_normal(11)="11  VENUS"

    name_func(1)= " 1  gravity and gravitational potential: W, V "
    name_func(2)= " 2  anomolous potential: T = W -- U           "
    name_func(3)= " 3  height anomaly: Z = T/gamma_q             "
    name_func(4)= " 4  geoid undulation: N = T0/gamma_0 + ADCB   "
    name_func(5)= " 5  spherical app gravity disturbance: - dT/dr"
    name_func(6)= " 6  sph app gravity anomaly: 2/r T-dT/dr      "
    name_func(7)= " 7  height: H                                 "
    name_func(8)= " 8  CS coef defines output: none              "
    name_func(9)= " 9  DOV (xi , eta)                            "
    name_func(10)="10  gravity anomaly:: norm(g) - norm(gamma_q) "
    name_func(11)="11  gravity: gx, gy, gz (LNOF)                "
    name_func(12)="12  gravity disturbance: dgx,dgy,dgz (LNOF)   "
    name_func(13)="13  gravity disturbance: norm(g)- norm(gamma) "
    name_func(14)="14  dV/dx, dV/dy, dV/dz (LNOF)                "
    name_func(15)="15  dV/dphi, dV/dlam, dV/dr                   "
    name_func(16)="16  T_xx, T_yy, T_zz, T_xy, T_xz, T_yz        "
    name_func(17)="17  T_r, T_rr, T_p, T_pp, T_l, T_ll, T_rp, T_rl, T_pl"

    inquire(FILE='ioptions.txt', EXIST=file_exists)
    
    if (not(file_exists)) then
         call stop_error('file ioptions.txt does not exist.')
    endif
    
    open(1, file='ioptions.txt', action='read')
    
    read(1,'(a100)')option%egmfile
    read(1,*)option%nmin
    read(1,*)option%nmax
    read(1,*)option%inormal
    read(1,*)option%func
    read(1,*)mod
    read(1,*)option%south
    read(1,*)option%step_lat
    read(1,*)option%north
    read(1,*)option%west
    read(1,*)option%step_lon
    read(1,*)option%east
    read(1,*)option%h0

    if (mod=='grid') then
        option%mod=1;
    elseif(mod=='igrid')then
        option%mod=3;
    else
        option%mod=2;
    endif 
       
    read(1,'(a100)')option%ifile
    read(1,'(a100)')option%ofile
    read(1,*)option%nthread
    close(1)

    open(111, file='log_FGrS.txt', action='write')
   
    
    call write_option(111)
    call write_option(6)
    call check_option(111)

    !write(6, '(/a)')  'Reading SH_ceoff file ...'

    if (option%func==4) then

        print*, 'enter hnm: spherical harmonics of topography'
        read(*, *) option%hnm
        call readEGM( 'free' ,  option%hnm , option%nmax, cnmH, snmH, aegm, gmegm, egmtide)

    endif

    call readEGM( 'gfc' , option%egmfile, option%nmax, cnm, snm, aegm, gmegm, egmtide)

    write (6, '(a , 10x, *(a) )')  'Scale of EGM (m3/s2):    ', adjustl(num2str(aegm,'real', '(f12.3)'))
    write (6, '(a , 10x, *(a) )')  'GM of EGM (meter):       ' , adjustl(num2str(gmegm,'real', '(d20.12)'))
    write (6, '(a, 10x, *(a) )' )  'Tide system:             ', egmtide//'_tide'


    write (111, '(a , 10x, *(a) )')  'Scale of EGM (m3/s2):    ', adjustl(num2str(aegm,'real', '(f12.3)'))
    write (111, '(a , 10x, *(a) )')  'GM of EGM (meter):       ' , adjustl(num2str(gmegm,'real', '(d20.12)'))
    write (111, '(a, 10x, *(a) )' )  'Tide system:             ', egmtide//'_tide'

    if (option%mod == 1 .or. option%mod ==3 ) then

        if   (option%mod ==3)then
            call readHfile (option, lat, lon, hg)
        else
            option%step_lat=option%step_lat/3600d0
            option%step_lon=option%step_lon/3600d0
            call makegrid(option, lat, lon)
        endif

        write(tempi , '(f12.5)') option%south

        write( 6, '( a , 10x, a )' ) 'Latitude limit south:     ', trim(adjustl(tempi)) // ' deg'
        write( 111, '( a , 10x, a )' ) 'Latitude limit south:     ', trim(adjustl(tempi)) //' deg'

        write(tempi , '(f12.5)') option%north

        write (6, '( a , 10x, a )' ) 'Latitude limit north:     ', trim(adjustl(tempi)) //' deg'
        write (111, '( a , 10x, a )' ) 'Latitude limit north:     ', trim(adjustl(tempi)) //' deg'

        write(tempi , '(f12.5)') option%west
        write (6,  '( a , 10x, a )')'Longitude limit west:     ',trim(adjustl(tempi)) //' deg'
        write (111,  '( a , 10x, a )')'Longitude limit west:     ',trim(adjustl(tempi)) //' deg'

        write(tempi , '(f12.5)') option%east
        write (*, '( a , 10x, a )' )'Longitude limit east:     ',trim(adjustl(tempi)) //' deg'
        write (111, '( a , 10x, a )' )'Longitude limit east:     ',trim(adjustl(tempi)) //' deg'

        write(tempi, '(f7.2)')  option%step_lat*60

        write(*,  '(a , 10x, *(a) )' )  'Gridsize lat (arcmin):   ', trim(adjustl(tempi))
        write(111,  '(a , 10x, *(a) )' )  'Gridsize lat (arcmin):   ', trim(adjustl(tempi))

        write(tempi, '(f7.2)')  option%step_lon*60
        write( *,  '(a , 10x, *(a) )' )  'Gridsize lon (arcmin):   ', trim(adjustl(tempi))
        write( 111,  '(a , 10x, *(a) )' )  'Gridsize lon (arcmin):   ', trim(adjustl(tempi))

        write( tempi, *) size(lat,1)
        write(*, '(a , 10x, *(a) )' )  'Grid dimension in lat:   ', trim(adjustl(tempi))
        write(111, '(a , 10x, *(a) )' )  'Grid dimension in lat:   ', trim(adjustl(tempi))

        write(tempi, *) size(lon,1)
        write (*, '(a , 10x, *(a) )')  'Grid dimension in lon:   ',trim(adjustl(tempi))
        write (111, '(a , 10x, *(a) )')  'Grid dimension in lon:   ',trim(adjustl(tempi))

    endif

    if (option%mod == 3) then

        if ( option%func == 4 .or. option%func == 7  .or. option%func == 7) then

            write(111,*) 'error: output functional must be solid (3D) function'
            call stop_error ('error: output functional must be solid (3D) function')

        endif

        call readHfile (option, lat, lon, hg)

    elseif (option%mod == 2 .or. option%mod == 3) then

        call ReadScatterFile (option, lat, lon, h)

        if (option%mod == 2 )then
            write(*, '(/a)' )  'scattered points statistics:'
        else
            write(*, '(/a)' )  'irregular points statistics:'
        endif

        write(6, '(a, 10x, i0 )')  'Number computation point:    ', size(lat,1)
        write(111, '(a, 10x, i0 )')  'Number computation point:    ', size(lat,1)

        write(tempi , '(f12.5)') option%south

        write( *, '( a , 10x, a )' ) 'Latitude limit south:     ', trim(adjustl(tempi)) //' deg'
        write( 111, '( a , 10x, a )' ) 'Latitude limit south:     ', trim(adjustl(tempi)) //' deg'

        write(tempi , '(f12.5)') option%north

        write (*, '( a , 10x, a )' ) 'Latitude limit north:     ', trim(adjustl(tempi)) //' deg'
        write (111, '( a , 10x, a )' ) 'Latitude limit north:     ', trim(adjustl(tempi)) //' deg'

        write(tempi , '(f12.5)') option%west
        write (*,  '( a , 10x, a )')'Longitude limit west:     ',trim(adjustl(tempi)) //' deg'
        write (111,  '( a , 10x, a )')'Longitude limit west:     ',trim(adjustl(tempi)) //' deg'

        write(tempi , '(f12.5)') option%east
        write (*, '( a , 10x, a )' )'Longitude limit east:     ',trim(adjustl(tempi)) //' deg'
        write (111, '( a , 10x, a )' )'Longitude limit east:     ',trim(adjustl(tempi)) //' deg'

        write(tempi , '(f8.3)') minval(h)
        write (*, '( a , 10x, a )' )  'min height (meter):       ', trim(adjustl(tempi))
        write (111, '( a , 10x, a )' )  'min height (meter):       ', trim(adjustl(tempi))

        write(tempi , '(f8.3)') maxval(h)
        write (*, '(   a, 10x, a )' ) 'max height (meter):       ', trim(adjustl(tempi))
        write (111, '(   a, 10x, a )' ) 'max height (meter):       ', trim(adjustl(tempi))

    endif

    if (option%func == 4 .or. option%func == 7) then

        select case (option%mod)
        case(1)
            option%h0=0
        case(2)
            h=0
        case(3)
            hg=0
        end select

    endif

    call tip_func (option%func)

    nmin=option%nmin
    nmax=option%nmax

    call GRS( option%inormal)

    call sub_normal(option, aegm, gmegm, cnm, klm)

    write(*, '(/a)' ) 'computation start... '

    call tic(start, rate)

    !$ call OMP_SET_NUM_THREADS(option%nthread);

    select case (option%mod)
    case(1)
        f3 = shs_grid (.true., cnm, snm, cnmH, snmH, lat, lon, option%h0, GMegm, aegm )
    case(2)

        f = shs_scatter (cnm, snm, cnmH, snmH, lat, lon, h, GMegm, aegm )
    case(3)
        f3 = shs_igrid (cnm, snm, cnmH, snmH, lat, lon, hg, GMegm, aegm )
    end select

    !write(*, '(/a)') 'progress completed.'

    write(*, '(/a)')  toc(start, rate)
  !  write(1,*) option%nthread, toc(start, rate)

    write(*, '(/a)') 'write results into file= ' // trim(option%ofile)
     
    write(111, '(a, 10x, a)') 'Write results into file: ' , trim(option%ofile)

    open(10,file= option%ofile, action='write')
 
    wid=0 
 
    select case (option%mod)
     
    case(2)

        do i=1,size(f,1)

            write (10, '(f12.5 , 1x, f12.5 , 1x , f12.3, 1x , 10e25.15)') &
                lat(i) , lon(i), h(i), f(i,:)
            call progress(wid, i,size(f,1))
        enddo
        deallocate(f, lat, lon)

    case(1,3)

        do i=1,size(lat,1)
            
            do j=1, size(lon,1)
        
                if (option%mod == 3) option%h0 = hg(i,j)
        
                write (10, '(f25.15 , 1x, f25.15 , 1x , f12.3, 1x , 10e25.15)') &
                    lat(i), lon(j), option%h0, f3(i, j ,1:iprint)
                
            enddo
            call progress(wid, i,size(lat,1))
        enddo
        deallocate(f3, lat, lon)

    end select

    write(*, '(//a)') 'program terminate normally.'
    write(111, '(//a)') 'program terminate normally.'

    close(111)
    close(10)
    
    end
