! #############################################################################
! This module stores the Remesh parameters to MFSim software.
!
! #############################################################################
module data_remesh
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: Tremesh
    !-------------------------------------
    
    ! Type for regrid parameters
    type Tremesh
        ! Remesh frequency in iterations
        integer, allocatable:: freq(:)
        ! Remesh on restart
        logical:: on_restart
        ! Berger-Rigoutsis Remesh algorithm CutOff
        double precision:: cut_off
        ! Buffer Zone from flagged point
        integer:: buff_fp
        ! Functions
        contains
            ! Insert default values in type
            procedure:: init => TremeshDefault
            ! Deallocate memory of type
            procedure:: destroy => TremeshDeallocate
    end type

contains


! -----------------------------------------------------------------------------
! Deallocate components in Tremesh type
pure subroutine TremeshDeallocate( remsh )
    implicit none
    !Entrada:
    class(Tremesh), intent(inout) :: remsh
    
    if (allocated(remsh%freq)) then
        deallocate(remsh%freq)
    end if
    
end subroutine TremeshDeallocate
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Initialize Remesh default parameters
pure subroutine TremeshDefault( remsh, lev )
    implicit none
    !Entrada:
    class(Tremesh), intent(inout) :: remsh
    integer, intent(in):: lev
    
    ! Default sets
    allocate(remsh%freq(max(lev,1)))
    remsh%freq(:) = -1
    remsh%on_restart = .false.
    remsh%cut_off = 0.7d0
    remsh%buff_fp = 2
    
end subroutine TremeshDefault
! -----------------------------------------------------------------------------


end module data_remesh