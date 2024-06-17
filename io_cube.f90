subroutine read_cube
    use coredens
    implicit none

    integer :: iatom, zval
    integer :: i1, i2, i3, iline, idir
    real*8  :: latvec_mat(3,3), recip_mat(3,3)
    logical :: stat
    real*8, allocatable :: tmp(:)
    integer :: index
    character(len = 200) :: line
    real*8  :: pos_d(3)

    nspin = 1

    inquire( file=trim('SPIN1_CHG.cube'), exist=stat )
    if(.not. stat) then
        write(*,*) 'Error: SPIN1_CHG.cube does not exist'
        call exit(1)
    end if

    inquire( file=trim('SPIN2_CHG.cube'), exist=stat )
    if(stat) then
        nspin = 2
    end if

    open(unit = 10, file = 'SPIN1_CHG.cube')

! read header
    read(10,*) ! Skip the first two lines
    read(10,*)
    read(10,*) natom
    read(10,*) n1, dx1
    read(10,*) n2, dx2
    read(10,*) n3, dx3

    dv = abs(dot_product(dx1, xprod(dx2,dx3)))

    latvec1 = n1 * dx1
    latvec2 = n2 * dx2
    latvec3 = n3 * dx3

    latvec_mat(:,1) = latvec1
    latvec_mat(:,2) = latvec2
    latvec_mat(:,3) = latvec3

    call m33inv(latvec_mat, recip_mat, stat)
    if(.not. stat) then
        write(*,*) 'Error: lattice vectors are not invertible'
        call exit(1)
    end if

    bvec1 = recip_mat(1,:)
    bvec2 = recip_mat(2,:)
    bvec3 = recip_mat(3,:)

! read atomic information
    allocate(z_of_atom(natom), ncore_of_atom(natom), pos_of_atom(natom,3))
    do iatom = 1, natom
        read(10,*) z_of_atom(iatom), zval, pos_of_atom(iatom,:)
        ncore_of_atom(iatom) = z_of_atom(iatom) - zval
    end do

! bring to ucell at origin
! seems not necessary as tau_d in abacus is already in the unit cell
!    do iatom = 1, natom
!        pos_d(1) = sum(pos_of_atom(iatom,:)*bvec1)
!        pos_d(2) = sum(pos_of_atom(iatom,:)*bvec2)
!        pos_d(3) = sum(pos_of_atom(iatom,:)*bvec3)
!        do idir = 1,3
!            do while(pos_d(idir) > 1.0d0)
!                pos_d(idir) = pos_d(idir) - 1.0d0
!            enddo
!            do while(pos_d(idir) < 0.0d0)
!                pos_d(idir) = pos_d(idir) + 1.0d0
!            enddo
!        end do
!
!        pos_of_atom(iatom,:) = pos_d(1) * latvec1 + pos_d(2) * latvec2 + pos_d(3) * latvec3
!    enddo

! read valence density
    allocate(rho(nspin,n1,n2,n3),tmp(n1*n2*n3))

    read(10,*) tmp

    index = 1
    do i1 = 1, n1
        do i2 = 1, n2
            do i3 = 1, n3
                rho(1,i1,i2,i3) = tmp(index)
                index = index + 1
            end do
        end do
    end do

    close(10)

    if(nspin .eq. 2) then

        open(unit = 10, file = 'SPIN2_CHG.cube')
        read(10,*) tmp

        do iline = 1, natom + 6
            read(10,'(A)') line
        enddo

        index = 1
        do i1 = 1, n1
            do i2 = 1, n2
                do i3 = 1, n3
                    rho(2,i1,i2,i3) = tmp(index)
                    index = index + 1
                end do
            end do
        end do

        close(10)
    endif

    deallocate(tmp)

contains
    function xprod(x,y)
    implicit none
    real*8, dimension(3) :: xprod, x, y
    xprod(1) = x(2)*y(3) - x(3)*y(2)
    xprod(2) =-x(3)*y(1) + x(1)*y(3)
    xprod(3) = x(1)*y(2) - x(2)*y(1)
    end function xprod
end subroutine read_cube

subroutine write_cube
    use coredens
    implicit none

    character(len = 200) :: line
    integer :: iline, index
    integer :: i1,i2,i3

    open(unit = 10, file = 'SPIN1_CHG.cube')
    open(unit = 11, file = 'SPIN1_CHG_core.cube')
    do iline = 1, natom + 6
        read(10,'(A)') line
        write(11,'(A)') trim(line)
    enddo

    index = 1
    do i1 = 1, n1
        do i2 = 1, n2
            do i3 = 1, n3
                write(11,"(e13.4)",advance="no") rho(1,i1,i2,i3)
                index = index + 1
                if(i3 == n3 .or. (mod(i3,6) == 0)) write(11,*)
            end do
        end do
    end do

    close(10)
    close(11)

    if(nspin .eq. 2) then
        open(unit = 10, file = 'SPIN2_CHG.cube')
        open(unit = 11, file = 'SPIN2_CHG_core.cube')
        do iline = 1, natom + 6
            read(10,'(A)') line
            write(11,'(A)') trim(line)
        enddo

        index = 1
        do i1 = 1, n1
            do i2 = 1, n2
                do i3 = 1, n3
                    write(11,"(e13.4)",advance="no") rho(2,i1,i2,i3)
                    index = index + 1
                    if(i3 == n3 .or. (mod(i3,6) == 0)) write(11,*)
                end do
            end do
        end do

        close(10)
        close(11)
    endif
end subroutine