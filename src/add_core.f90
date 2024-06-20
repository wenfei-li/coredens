subroutine add_core
    use coredens
    use coredensity_lib
    implicit none

    integer :: i1,i2,i3,index
    real*8  :: tmp,tmp1
    integer :: iat,iat_ucell
    real*8  :: pos_grid(3), pos_atom(3),r
    integer :: zat,ncore
    real*8  :: alf(100),coe(100)
    integer :: nfun,ifun
    real*8  :: tot_ncore
    
    real*8, allocatable :: rhocore(:,:,:)
    
    allocate(rhocore(n1,n2,n3))
    rhocore = 0.0
    
    index = 1
    do i1 = 1, n1
        do i2 = 1, n2
            do i3 = 1, n3
                pos_grid = (i1-1)*dx1 + (i2-1)*dx2 + (i3-1)*dx3
                tmp = 0.0
                do iat = 1, natom_around_grid(index)
                    pos_atom = pos_of_atom_ext(atom_around_grid(index,iat),:)
                    r = norm2(pos_grid - pos_atom)
                    iat_ucell = mod(atom_around_grid(index,iat)-1,natom)+1
                    zat = z_of_atom(iat_ucell)
                    ncore = ncore_of_atom(iat_ucell)

                    if(database_type == 0) then
                        call edflib(zat,ncore,nfun,alf,coe)
                        do ifun = 1, nfun
                            tmp = tmp + coe(ifun)*exp(-alf(ifun)*r**2)
                        enddo
                    else
                        call get_core_dens(zat,r,tmp1)
                        tmp = tmp + tmp1
                    endif
                enddo
                index = index + 1

                rhocore(i1,i2,i3) = rhocore(i1,i2,i3) + tmp
            end do
        end do
    end do

    tot_ncore = sum(ncore_of_atom)
    rhocore = rhocore * tot_ncore / (sum(rhocore) * dv) / dble(nspin)
    ! Need to renormalize; because the grid for valence density typically
    ! will not be sufficient for core density ...
    rho(1,:,:,:) = rho(1,:,:,:) + rhocore

    if(nspin == 2) then
        rho(2,:,:,:) = rho(2,:,:,:) + rhocore
    endif

    deallocate(rhocore)
end subroutine