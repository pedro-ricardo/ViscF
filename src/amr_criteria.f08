! #############################################################################
! This module contains the sets of functions that represent the refinement
! criteria available.
! #############################################################################
module amr_criteria
    use data_variable
    use data_boxes
    use data_input
    use amr_coords
    use amrex_amr_module

    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: tag2AMReX, initial_tag
    !-------------------------------------
    
    double precision, parameter:: mark=10.d0, clear=0.d0
    double precision, parameter:: eps=1.d-2

contains

! -----------------------------------------------------------------------------
! Translate the tag variable marked with custom criteria to AMReX tagging
subroutine tag2AMReX(patch, settag, cleartag, tagp, taglo, taghi, atagp, ataglo, ataghi)
    implicit none
    !Entrada:
    type(Tpatch), intent(in):: patch
    character(kind=c_char), intent(in):: settag, cleartag
    integer, intent(in):: taglo(4), taghi(4), ataglo(4), ataghi(4)
    double precision, intent(in) :: tagp(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    character(kind=c_char), intent(inout) :: atagp(ataglo(1):ataghi(1),ataglo(2):ataghi(2),ataglo(3):ataghi(3))
    !Local:
    integer:: i, j, k
    
    do k = patch%ilo(3), patch%ihi(3)
        do j = patch%ilo(2), patch%ihi(2)
           do i = patch%ilo(1), patch%ihi(1)

                if (tagp(i,j,k)>(mark-eps)) then
                    ! Mark cell to refine
                    atagp(i,j,k) = settag
                elseif (tagp(i,j,k)<(clear+eps)) then
                    ! Mark cell to coarse
                    atagp(i,j,k) = cleartag
                endif

           enddo
        enddo
     enddo

end subroutine tag2AMReX
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Use setted values from input to fill 
subroutine initial_tag(t,dt,patch,ngcl,Pvar)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    integer, intent(in):: ngcl
    type(Tpatch), intent(in):: patch
    type(Tpatchbox), target:: Pvar(:)
    !Local:
    integer:: n, i, j ,k
    integer:: npatches, sn(3)
    type(Tpatch), allocatable:: init_patch(:)
    double precision:: x, y, z
    double precision, contiguous, pointer:: var(:,:,:)


    ! Extract from input to reduce name lenght
    npatches = size(input%msh%init_grid)
    allocate(init_patch(npatches))
    init_patch = input%msh%init_grid

    ! Get variable in correct indexes
    var => Pvar(1)%patch3D()
    ! Face info
    sn = Pvar(1)%sn
    
    ! Loop in input patches
    do n = 1,npatches
        ! Loop Z
        do k = patch%ilo(3)-ngcl, patch%ihi(3)+ngcl
            z = index2coord(amrex_problo(3),k,patch%dn(3),sn(3))
            ! Check for patch limits
            if ( (z<=init_patch(n)%hi(3)).and.(z>=init_patch(n)%lo(3)) ) then
                ! Loop in Y
                do j = patch%ilo(2)-ngcl, patch%ihi(2)+ngcl
                    y = index2coord(amrex_problo(2),j,patch%dn(2),sn(2))
                    ! Check for patch limits
                    if ( (y<=init_patch(n)%hi(2)).and.(y>=init_patch(n)%lo(2)) ) then
                        ! Loop in X
                        do i = patch%ilo(1)-ngcl, patch%ihi(1)+ngcl
                            x = index2coord(amrex_problo(1),i,patch%dn(1),sn(1))
                            ! Check for patch limits
                            if ( (x<=init_patch(n)%hi(1)).and.(x>=init_patch(n)%lo(1)) ) then

                                ! Mark it !!!
                                var(i,j,k) = mark
                                
                            end if
                        end do

                    end if
                end do

            end if
        end do

    end do
        
    
    deallocate(init_patch)

end subroutine initial_tag
! -----------------------------------------------------------------------------

end module amr_criteria