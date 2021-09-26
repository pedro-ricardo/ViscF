! #############################################################################
! This module contains the set of functions used to initialize a variable
! in MFSim
!
! #############################################################################
module amr_startvariables
    use data_variable
    use data_boxes
    use amr_coords
    use amrex_amr_module
    use data_input

    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !All subroutines/variables are
    !accessible outside this module
    !Except the folowing ones
    public:: zeros, ones, exact
    !-------------------------------------
    
contains

! -----------------------------------------------------------------------------
! Fill variable with zeros
subroutine zeros(t,dt,patch,ngcl,Pvar)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    integer, intent(in):: ngcl
    type(Tpatch), intent(in):: patch
    type(Tpatchbox), target:: Pvar(:)
    !Local:
    integer:: c
    double precision, contiguous, pointer:: var(:,:,:)
    
    ! Loop in components
    do c = 1,Pvar(1)%ncomp
        ! Get variable in correct indexes
        var => Pvar(1)%patch3D(c)
        ! Set all values in one pass
        var = 0.d0
    end do
    
end subroutine zeros
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Fill variable with ones
subroutine ones(t,dt,patch,ngcl,Pvar)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    integer, intent(in):: ngcl
    type(Tpatch), intent(in):: patch
    type(Tpatchbox), target:: Pvar(:)
    !Local:
    integer:: c
    double precision, contiguous, pointer:: var(:,:,:)
    
    ! Loop in components
    do c = 1,Pvar(1)%ncomp
        ! Get variable in correct indexes
        var => Pvar(1)%patch3D(c)
        ! Set all values in one pass
        var = 1.d0
    end do

end subroutine ones
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Fill variable with cartesian function
subroutine exact(t,dt,patch,ngcl,Pvar)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    integer, intent(in):: ngcl
    type(Tpatch), intent(in):: patch
    type(Tpatchbox), target:: Pvar(:)
    !Local:
    integer:: c, i ,j ,k, sn(3)
    double precision:: x, y, z
    double precision, contiguous, pointer:: var(:,:,:)
    
    ! Loop in components
    do c = 1,Pvar(1)%ncomp
        ! Get variable in correct indexes
        var => Pvar(1)%patch3D(c)
        ! Face info
        sn = Pvar(1)%sn
        
        do k = patch%ilo(3)-ngcl, patch%ihi(3)+sn(3)+ngcl
            z = index2coord(amrex_problo(3),k,patch%dn(3),sn(3))
            do j = patch%ilo(2)-ngcl, patch%ihi(2)+sn(2)+ngcl
                y = index2coord(amrex_problo(2),j,patch%dn(2),sn(2))
                do i = patch%ilo(1)-ngcl, patch%ihi(1)+sn(1)+ngcl
                    x = index2coord(amrex_problo(1),i,patch%dn(1),sn(1))

                    var(i,j,k) = exact_fun(t,x,y,z)

                end do
            end do
        end do
        
    end do

end subroutine exact
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Temporary function
function exact_fun(t,x,y,z) result(exact)
    implicit none
    !Entrada:
    double precision, intent(in):: t, x, y, z
    !Saída:
    double precision:: exact
    
    exact = t + sqrt(x**2 + y**2 + z**2)
    
end function exact_fun
! -----------------------------------------------------------------------------

end module amr_startvariables