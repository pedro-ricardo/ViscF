! #############################################################################
! This module has the object to store the pair name/variable in order to
! print them correctly
!
! #############################################################################
module output_names
    use amrex_string_module
    use data_variable

    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: Toutput
    !-------------------------------------

    ! Output pair type
    type Toutput
        ! Variable name
        type(amrex_string):: name
        ! Variable data pointer
        type(TamrVar), pointer:: var => null()
        ! Flag to print it
        logical:: print = .false.
        ! Functions
        contains
            ! Set a pair name/var
            procedure:: set => TpairAdd
            ! Choose to print or not
            procedure:: setPrint =>TpairSetPrint
            ! Clear the type
            procedure:: clear => TpairDestroy
            final:: destroyTpair
    end type
    
contains

! -----------------------------------------------------------------------------
! Insert a name/var pair, overwritting the previous one if any
subroutine TpairAdd(out, name, var)
    implicit none
    !Entrada:
    class(Toutput), intent(inout):: out
    character(len=*), intent(in):: name
    type(TamrVar), target, intent(in):: var
    !Local:
    
    ! Clean type
    call out%clear()
    
    ! Transform the Fortran String in AMReX type
    call amrex_string_build(out%name,name)
    ! Point to the variable
    out%var => var
    ! Start with print disabled
    out%print = .false.
    
    
end subroutine TpairAdd
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Clean up the name/var pair
subroutine TpairDestroy(out)
    implicit none
    !Entrada:
    class(Toutput), intent(inout):: out
    !Local:

    if (allocated(out%name%data)) then
        deallocate(out%name%data)
    end if
    ! Clear pointer
    nullify(out%var)
    
end subroutine TpairDestroy
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Clean up the name/var pair.
! But with a different declaration needed by 'final' procedure
subroutine destroyTpair(out)
    implicit none
    !Entrada:
    type(Toutput), intent(inout):: out
    !Local:

    ! Call the 'TpairDestroy' function.
    call out%clear()

end subroutine destroyTpair
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Set the pair to be printable
subroutine TpairSetPrint(out,bool)
    implicit none
    !Entrada:
    class(Toutput), intent(inout):: out
    logical, intent(in):: bool
    !Local:
    
    out%print = bool
    
end subroutine TpairSetPrint
! -----------------------------------------------------------------------------

end module output_names