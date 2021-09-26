! #############################################################################
! This module holds the functions used to initialize the variables in this
! example.
! #############################################################################
module example_functions
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
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: exact_rhs
    public:: exact_d2pdx2, exact_d2pdy2, exact_d2pdz2, exact_p
    !-------------------------------------
    
    double precision, parameter:: w1 = 2.d0, w2 = 2.d0
    double precision, parameter:: w3 = 2.d0, w4 = 1.d0
    double precision, parameter:: pi = atan(1.d0)*4.d0  

contains

! -----------------------------------------------------------------------------
! Patch function to apply the selected exact function in the variable
subroutine exact_rhs(t,dt,patch,ngcl,Pvar)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    integer, intent(in):: ngcl
    type(Tpatch), intent(in):: patch
    type(Tpatchbox), target:: Pvar(:)
    !Local:
    integer:: c, i ,j ,k, sn(3)
    double precision:: x, y, z
    double precision, contiguous, pointer:: rhs(:,:,:)
    
    ! Loop in components
    do c = 1,Pvar(1)%ncomp
        ! Get variable in correct indexes
        rhs => Pvar(1)%patch3D(c)
        ! Face info
        sn = Pvar(1)%sn
        
        do k = patch%ilo(3)-ngcl, patch%ihi(3)+sn(3)+ngcl
            z = index2coord(amrex_problo(3),k,patch%dn(3),sn(3))
            do j = patch%ilo(2)-ngcl, patch%ihi(2)+sn(2)+ngcl
                y = index2coord(amrex_problo(2),j,patch%dn(2),sn(2))
                do i = patch%ilo(1)-ngcl, patch%ihi(1)+sn(1)+ngcl
                    x = index2coord(amrex_problo(1),i,patch%dn(1),sn(1))

                    rhs(i,j,k) = exact_d2pdx2(t,x,y,z)+ exact_d2pdy2(t,x,y,z)+ exact_d2pdz2(t,x,y,z)

                end do
            end do
        end do
        
    end do

end subroutine exact_rhs
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Standard manufactured solution in MFSim:
! Pressure
function exact_p(t,x,y,z)
    implicit none
    double precision, intent(in):: t, x, y, z
    double precision:: exact_p

    exact_p = sin(pi*w3*z+pi*w2*y+pi*w1*x+t*w4)**2.d+0
    
end function exact_p
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Standard manufactured solution in MFSim:
! Pressure second derivative in X
function exact_d2pdx2(t,x,y,z)
    implicit none
    double precision, intent(in):: t, x, y, z
    double precision:: exact_d2pdx2

    exact_d2pdx2 = 2.d+0*pi**2*w1**2*cos(pi*w3*z+pi*w2*y+pi*w1*x+t*w4)**2 &
        -2.d+0*pi**2*w1**2*sin(pi*w3*z+pi*w2*y+pi*w1*x+t*w4)**2

end function exact_d2pdx2
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Standard manufactured solution in MFSim: 
! Pressure second derivative in Y
function exact_d2pdy2(t,x,y,z)
    implicit none
    double precision, intent(in):: t, x, y, z
    double precision:: exact_d2pdy2

    exact_d2pdy2 = 2.d+0*pi**2*w2**2*cos(pi*w3*z+pi*w2*y+pi*w1*x+t*w4)**2 &
        -2.d+0*pi**2*w2**2*sin(pi*w3*z+pi*w2*y+pi*w1*x+t*w4)**2

end function exact_d2pdy2
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Standard manufactured solution in MFSim:
! Pressure second derivative in Z
function exact_d2pdz2(t,x,y,z)
    implicit none
    double precision, intent(in):: t, x, y, z 
    double precision:: exact_d2pdz2

    exact_d2pdz2 = 2.d+0*pi**2*w3**2*cos(pi*w3*z+pi*w2*y+pi*w1*x+t*w4)**2 &
        -2.d+0*pi**2*w3**2*sin(pi*w3*z+pi*w2*y+pi*w1*x+t*w4)**2

end function exact_d2pdz2
! -----------------------------------------------------------------------------

end module example_functions