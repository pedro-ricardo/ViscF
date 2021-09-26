! #############################################################################
! This module stores the Input/Output parameters to MFSim software.
!
! #############################################################################
module data_io
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: Tio
    !-------------------------------------
    
    ! Type for output control parameters
    type Tio
        ! Maximum number restart files maintained
        integer:: rst_max
        ! Name prefix in restart files
        character(:), allocatable:: rst_prefix
        ! Save restart frequency in iterations
        integer:: rst_it
        ! Save restart frequency in seconds
        double precision:: rst_time
        ! Name prefix in results files
        character(:), allocatable:: out_prefix
        ! Save results frequency in iterations
        integer:: out_it
        ! Save results frequency in seconds
        double precision:: out_time
        ! Output variable name list
        character(len=32), allocatable:: out_names(:)
        ! Functions
        contains
            ! Insert default values in type
            procedure:: init => TioDefault
            ! Deallocate memory of type
            procedure:: destroy => TioDeallocate
    end type

contains

! -----------------------------------------------------------------------------
! Set default values for Tio type
pure subroutine TioDefault( io )
    implicit none
    !Entrada:
    class(Tio), intent(inout) :: io
    
    ! Default sets
    io%rst_max = 2
    io%rst_prefix = 'chk'
    io%rst_it = 500
    io%rst_time = -1.d0

    io%out_prefix = 'plt'
    io%out_it = 50
    io%out_time = -1.d0

    allocate(io%out_names(4))
    io%out_names(1:3) = ['U','V','W']
    io%out_names(4) = 'Pressure'
    
end subroutine TioDefault
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Deallocate components in Tio type
pure subroutine TioDeallocate( io )
    implicit none
    !Entrada:
    class(Tio), intent(inout) :: io
    
    if (allocated(io%rst_prefix)) then
        deallocate(io%rst_prefix)
    end if
    
    if (allocated(io%out_prefix)) then
        deallocate(io%out_prefix)
    end if

    if (allocated(io%out_names)) then
        deallocate(io%out_names)
    end if

end subroutine TioDeallocate
! -----------------------------------------------------------------------------



end module data_io