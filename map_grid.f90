subroutine map_grid
    use coredens
    implicit none

    real*8  :: rcut = 6  ! cutoff radius for core densities
    ! at 6 Bohr, the Gaussians are already vanishingly small

    integer :: nboxext_1, nboxext_2, nboxext_3
    integer :: ibox1, ibox2, ibox3, iatom_ext, iatom
    integer :: i1,i2,i3,index
    real*8  :: pos_grid(3)
    integer :: natom_around_grid_max, ii

    ! Extend the supercell by rcut
    nboxext_1 = ceiling(rcut * norm2(bvec1))
    nboxext_2 = ceiling(rcut * norm2(bvec2))
    nboxext_3 = ceiling(rcut * norm2(bvec3))

    natom_ext = (2*nboxext_1+1) * (2*nboxext_2+1) * (2*nboxext_3+1) * natom
    allocate(pos_of_atom_ext(natom_ext,3))

    ! Calculate positions of atoms in the extended cell
    iatom_ext = 1
    do ibox1 = -nboxext_1,nboxext_1
        do ibox2 = -nboxext_2,nboxext_2
            do ibox3 = -nboxext_3,nboxext_3
                do iatom = 1,natom
                    pos_of_atom_ext(iatom_ext,:) = pos_of_atom(iatom,:) + ibox1*latvec1 + ibox2*latvec2 + ibox3*latvec3
                    iatom_ext = iatom_ext + 1
                enddo
            enddo
        enddo
    enddo

    ! Search for neighbouring grid points of atoms in the extended cell
    allocate(natom_around_grid(n1*n2*n3))
    index = 1
    do i1 = 1,n1
        do i2 = 1,n2
            do i3 = 1,n3
                pos_grid = (i1-1)*dx1 + (i2-1)*dx2 + (i3-1)*dx3
                natom_around_grid(index) = 0
                do iatom_ext = 1,natom_ext
                    if (norm2(pos_of_atom_ext(iatom_ext,:) - pos_grid) < rcut) then
                        natom_around_grid(index) = natom_around_grid(index) + 1
                    endif
                enddo
                index = index + 1
            enddo
        enddo
    enddo
    natom_around_grid_max = maxval(natom_around_grid)
    allocate(atom_around_grid(n1*n2*n3,natom_around_grid_max))

    index = 1
    do i1 = 1,n1
        do i2 = 1,n2
            do i3 = 1,n3
                pos_grid = (i1-1)*dx1 + (i2-1)*dx2 + (i3-1)*dx3
                ii = 0
                do iatom_ext = 1,natom_ext
                    if (norm2(pos_of_atom_ext(iatom_ext,:) - pos_grid) < rcut) then
                        ii = ii + 1
                        atom_around_grid(index,ii) = iatom_ext
                    endif
                enddo
                if(ii /= natom_around_grid(index)) then
                    print *, 'Error in map_grid: ii /= natom_around_grid(index)'
                    stop
                endif
                index = index + 1
            enddo
        enddo
    enddo
end subroutine map_grid