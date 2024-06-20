module coredensity_lib
    implicit none

! List of elements to prepare
    integer :: ntype
    integer, allocatable :: zat_list(:)
    integer :: map_zat(118)

! Corresponding core densities and radial grids
    integer, allocatable :: ngrid(:)
    real*8,  allocatable :: radial_grid(:,:)
    real*8,  allocatable :: core_densities(:,:)

! For spline fitting
    real*8,  allocatable :: y2(:,:)

contains
    subroutine map_elements(z_of_atom)
        implicit none

        integer,intent(in) :: z_of_atom(:)
        integer :: index, i

        map_zat = 0
        index = 0
        do i=1,size(z_of_atom)
            if(map_zat(z_of_atom(i)) == 0 ) then
                index = index + 1
                map_zat(z_of_atom(i)) = index
            end if
        end do

        ntype = maxval(map_zat)

        allocate(zat_list(ntype))
        do i = 1,118
            if(map_zat(i) /= 0) zat_list(map_zat(i)) = i
        end do
    end subroutine

    subroutine prepare_coredensity
        implicit none

        integer :: ngrid_max, itype, io, nlines, iline
        real*8  :: r, rho, fourpi = 12.566370614359172954
        character(len=100) :: filename
        character(len=10)  :: file_id

        ngrid_max = 0
        do itype = 1,ntype
            write(file_id,'(i0)') zat_list(itype)
            filename = 'core_dens_' // trim(adjustl(file_id)) // '.dat'
            
            nlines = 0
            open(unit = 10, file = filename, status = 'old')
            do
                read(10,*,iostat = io)
                if(io/=0) exit
                nlines = nlines + 1
            enddo
            close(10)

            ngrid_max = max(ngrid_max,nlines)
        end do

        allocate(ngrid(ntype), radial_grid(ntype,ngrid_max), &
            core_densities(ntype,ngrid_max), y2(ntype,ngrid_max))

        radial_grid = 0.0
        core_densities = 0.0
        y2 = 0.0

        do itype = 1,ntype
            write(file_id,'(i0)') zat_list(itype)
            filename = 'core_dens_' // trim(adjustl(file_id)) // '.dat'

            open(unit = 10, file = filename, status = 'old')
            nlines = 0
            do
                read(10,*,iostat = io) r, rho
                if(io/=0) exit
                nlines = nlines + 1
                radial_grid(itype,nlines) = r
                core_densities(itype,nlines) = rho / r / r / fourpi
            enddo
            ngrid(itype) = nlines
            close(10)

            if(radial_grid(itype,ngrid(itype)) < 6) then
                write(*,*) 'radial grid too small for element: ',zat_list(itype)
                call exit(1)
            end if
            call SPLINE_dn(radial_grid(itype,:nlines),core_densities(itype,:nlines),&
                nlines,y2(itype,:nlines))
        enddo
    end subroutine

    subroutine get_core_dens(zat,r,val)
        implicit none
        
        integer :: zat, index
        real*8  :: r, val

        index = map_zat(zat)
        call SPLINT_dn(radial_grid(index,:ngrid(index)),&
            core_densities(index,:ngrid(index)),&
            y2(index,:ngrid(index)),ngrid(index),r,val)
    end subroutine
end module coredensity_lib