! #############################################################################
! This module has the types and functions to compress a list of variables and
! constants into a box. This box is needed for any mesh operation done with
! `AMRmesh%forCells` calls.
! #############################################################################
module data_boxes
    use data_variable
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: Tbox, Tpatchbox
    !-------------------------------------
    
    ! Type to represent a box full of variables and constants
    ! that can be unpack inside a patch function
    type, public:: Tbox
        ! Number of constans in box
        integer:: ncte = 0
        ! Number of variables in box
        integer:: nvar = 0
        ! Pile of variable boxes
        type(TboxVar), allocatable:: Vars(:)
        ! Pile of variable boxes
        type(Tpatchbox), allocatable:: Patch(:)
        ! Pile of contant boxes
        type(TboxCte), allocatable:: Cbox(:)
        ! Function
        contains
            ! Insert variable in box
            procedure:: addVar => TboxVarAdd
            ! Insert constant in box
            procedure:: addCte => TboxCteAdd
            ! Deallocate when finished
            procedure:: cleanVars => TboxCleanVars
            ! Deallocate when finished
            procedure:: cleanCtes => TboxCleanCtes
            ! Deallocate when finished
            procedure:: destroy => TboxDestroy
            final:: TboxEnd
    end type

    ! Type to compress a variables
    type, public:: TboxVar
        ! Pointer to the desired variable
        type(TamrVar), pointer:: var => Null()
    end Type

    ! Type to compress a variable in patch index
    type, public:: Tpatchbox
        ! Pointer to the desired variable with patch indexes
        double precision, contiguous, pointer:: var(:,:,:,:) => Null()
        ! Number of components
        integer:: ncomp
        ! Variable limits
        integer:: ilim(3,2)
        ! Variable type (centered/faced)
        integer:: sn(3)
        ! Function
        contains
            ! Returns a 3D pointer with correct patch indexes
            procedure:: patch3D => TpatchboxPatch3D
    end Type

    ! Type to compress a constant
    type, public:: TboxCte
        ! Value of desired constant
        double precision:: cte = 0.d0
    end type

contains

! -----------------------------------------------------------------------------
! Insert a variable inside the box
subroutine TboxVarAdd(box,var)
    implicit none
    !Entrada:
    class(Tbox), intent(inout):: box
    type(TamrVar), target, intent(in):: var
    !Local:
    type(Tbox):: tmp
    integer:: i, n
    
    if(box%nvar==0) then ! First allocation
        n = 1
        allocate(box%Vars(n))
        ! Point to variable
        box%Vars(n)%var => var

    else ! Change in size
        tmp = box
        ! Clean up box
        call box%cleanVars()
        n = tmp%nvar+1
        allocate(box%Vars(n))
        ! Transfer old pointers
        do i=1,n-1
            box%Vars(i)%var => tmp%Vars(i)%var
        end do
        ! Point last one to variable
        box%Vars(n)%var => var
    end if

    ! Reallocate patch pointers
    allocate(box%Patch(n))

    ! Update counters
    if (allocated(box%Vars)) box%nvar = size(box%Vars)
    
end subroutine TboxVarAdd
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Insert a constant inside the box
subroutine TboxCteAdd(box,cte)
    implicit none
    !Entrada:
    class(Tbox), intent(inout):: box
    double precision, intent(in):: cte
    !Local:
    type(Tbox):: tmp
    integer:: i, n
    
    if(box%ncte==0) then ! First allocation
        allocate(box%Cbox(1))
        ! Set to constant
        box%Cbox(1)%cte = cte

    else ! Change in size
        tmp = box
        ! Clean up box
        call box%cleanCtes()
        n = tmp%ncte+1
        allocate(box%Cbox(n))
        ! Transfer old values
        do i=1,n-1
            box%Cbox(i)%cte = tmp%Cbox(i)%cte
        end do
        ! Set last one to constant
        box%Cbox(n)%cte = cte
    end if

    ! Update counters
    if (allocated(box%Cbox)) box%ncte = size(box%Cbox)
    
end subroutine TboxCteAdd
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Remove everything from box
subroutine TboxDestroy(box)
    implicit none
    !Entrada:
    class(Tbox), intent(inout):: box
    !Local:
    
    call box%cleanVars()
    call box%cleanCtes()
    
end subroutine TboxDestroy
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Same as 'TboxDestroy' but for final procedure
subroutine TboxEnd(box)
    implicit none
    !Entrada:
    type(Tbox), intent(inout):: box
    !Local:
    
    call box%destroy()
    
end subroutine TboxEnd
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Remove variables from box
subroutine TboxCleanVars(box)
    implicit none
    !Entrada:
    class(Tbox), intent(inout):: box
    !Local:
    integer:: i

    ! Check Variable box
    if (allocated(box%Vars)) then
        do i=1, size(box%Vars)
            ! Reset pointers
            nullify(box%Vars(i)%var)
        end do
        ! Deallocate it
        deallocate(box%Vars)
    end if

    ! Check Patch box
    if (allocated(box%Patch)) then
        do i=1, size(box%Patch)
            ! Reset pointers
            nullify(box%Patch(i)%var)
        end do
        ! Deallocate it
        deallocate(box%Patch)
    end if
    ! Zero counter
    box%nvar = 0

end subroutine TboxCleanVars
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Remove constants from box
subroutine TboxCleanCtes(box)
    implicit none
    !Entrada:
    class(Tbox), intent(inout):: box
    !Local:
    
    ! Check constant box
    if (allocated(box%Cbox)) then
        ! Deallocate it
        deallocate(box%Cbox)
    end if
    ! Zero counter
    box%ncte = 0
    
end subroutine TboxCleanCtes
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Routine executed inside Patch variable.
! Returns a pointer with correct access indexes
function TpatchboxPatch3D(Pvar,component) result(point3D)
    implicit none
    double precision, contiguous, pointer:: point3D(:,:,:)
    !Entrada:
    class(Tpatchbox), intent(inout):: Pvar
    integer, optional, intent(in):: component
    !Local:
    integer:: comp

    ! Default is first component
    comp = 1
    ! Get correct one if passed
    if (present(component)) then
        comp=component
    end if

    ! Transform in 3D with correct indexes
    point3D( Pvar%ilim(1,1):Pvar%ilim(1,2),&
             Pvar%ilim(2,1):Pvar%ilim(2,2),&
             Pvar%ilim(3,1):Pvar%ilim(3,2) ) => Pvar%var(:,:,:,comp)
    
end function TpatchboxPatch3D
! -----------------------------------------------------------------------------

end module data_boxes