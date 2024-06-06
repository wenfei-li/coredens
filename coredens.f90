! This is a simple program to add core electron density defined in edflib.f90
! to the valence electron density given in .cube format
! and output the total electron density in .cube format

module coredens
    implicit none

! Information of the 3D FFT grid
    integer :: n1,n2,n3
    real*8  :: dx1(3), dx2(3), dx3(3), dv
    real*8  :: latvec1(3), latvec2(3), latvec3(3)
    real*8  :: bvec1(3), bvec2(3), bvec3(3)

! Information of the atoms
    integer :: natom
    integer, allocatable :: z_of_atom(:), ncore_of_atom(:)
    real*8,  allocatable :: pos_of_atom(:,:)

! Array containing the density
    real*8,  allocatable :: rho(:,:,:)

! Atoms in extended cell
    integer :: natom_ext
    real*8, allocatable :: pos_of_atom_ext(:,:)

! Keep record of neighbouring atoms around grid points
    integer, allocatable :: natom_around_grid(:), atom_around_grid(:,:)

end module coredens

program main
    use coredens
    implicit none

    call read_cube
    call map_grid
    call add_core
    call write_cube
end program
