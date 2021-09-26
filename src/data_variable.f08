! #############################################################################
! This module has the definition and functions of a AMR variable
!
! #############################################################################
module data_variable
    use amrex_amr_module
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: TamrVar, Tpatch
    public:: center, faceX, faceY, faceZ, nodal
    !-------------------------------------
    
    ! Face ID flags
    integer, parameter:: center(3)=[0,0,0]
    integer, parameter:: faceX(3) =[1,0,0]
    integer, parameter:: faceY(3) =[0,1,0]
    integer, parameter:: faceZ(3) =[0,0,1]
    integer, parameter:: nodal(3) =[1,1,1]

    ! Type to save variable information
    type, public:: TamrVar
        ! Number of components in a variable
        integer:: components
        ! Number of ghost cells in variable
        integer:: ghosts
        ! Face/Centered information (by component)
        integer:: face_info(3)
        ! Boundary condition information
        integer, allocatable:: bc(:,:,:)
        ! Variable data
        type(amrex_multifab), allocatable:: mfab(:)
        ! Functions
        contains
            ! Allocate variable with levels only
            procedure:: init => TamrVarInit
            ! Deallocate variable
            procedure:: destroy => TamrVarDestroy
            ! Allocate only a certain level
            procedure:: levelInit => TamrVarLevelInit
            ! Deallocate only a certain level
            procedure:: levelDestroy => TamrVarLevelDestroy
    end type

    ! Type to save patch info
    type Tpatch
        ! Patch index limit
        integer:: ilo(3), ihi(3)
        ! Level of this patch
        integer:: lev
        ! Patch physical corners
        double precision:: lo(3), hi(3)
        ! Cells size
        double precision:: dn(3)
    end type

contains

! -----------------------------------------------------------------------------
! Perform correct allocation of AMR variable. Still no data in levels.
subroutine TamrVarInit(var, comp, ngc, sn)
    implicit none
    !Entrada:
    integer, intent(in):: comp, ngc
    integer, optional, intent(in):: sn(3)
    class(TamrVar), intent(inout):: var
    !Local:
    integer:: s(3)

    ! Default is centered variable
    s = center
    if(present(sn))then
        s = sn
    end if
    
    ! Clear if previously init
    call var%destroy()
    ! Save components
    var%components = comp
    ! Save ghosts
    var%ghosts = ngc
    ! Save face/center index types
    var%face_info = s

    ! Allocate levels
    allocate(var%mfab(amrex_max_level+1))
    ! Allocate BCs for each component
    allocate(var%bc(3,comp,2))

    ! Set default BCs as Neumann
    var%bc = amrex_bc_reflect_even
    
end subroutine TamrVarInit
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Perform correct level allocation
subroutine TamrVarLevelInit(var, lev, ba, dm)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    type(amrex_boxarray), intent(in) :: ba
    type(amrex_distromap), intent(in) :: dm
    class(TamrVar), intent(inout):: var
    !Local:
    
    ! Clear the level to re-init if needed
    call var%levelDestroy(lev)
    ! Create the level
    call amrex_multifab_build(var%mfab(lev), ba, dm, var%components, var%ghosts, var%face_info==1)
    
end subroutine TamrVarLevelInit
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Perform correct deallocation of AMR variable
subroutine TamrVarDestroy(var)
    implicit none
    !Entrada:
    class(TamrVar), intent(inout):: var
    !Local:
    integer:: lev

    ! Check init
    if (allocated(var%mfab)) then
        ! Clear each level
        do lev=1,size(var%mfab)
            call var%levelDestroy(lev)
        end do
        ! Free memory
        deallocate(var%mfab)
    end if

    if (allocated(var%bc)) deallocate(var%bc)
    
end subroutine TamrVarDestroy
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Perform correct level deallocation
subroutine TamrVarLevelDestroy(var,lev)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    class(TamrVar), intent(inout):: var
    !Local:
    
    ! Check init
    if (allocated(var%mfab)) then
        ! If level exists in variable clear it
        if (lev<=size(var%mfab)) then
            call amrex_multifab_destroy(var%mfab(lev))
        end if
    end if
    
end subroutine TamrVarLevelDestroy
! -----------------------------------------------------------------------------



end module data_variable