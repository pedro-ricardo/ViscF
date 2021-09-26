! #############################################################################
! This module stores the Simulation parameters to MFSim software.
!
! #############################################################################
module data_simulation
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: Tsimulation
    !-------------------------------------
    
    ! Type for simulation control parameters
    type Tsimulation
        ! Final simulation Time
        double precision:: end_time
        ! Maximum number of iterations
        integer:: max_it
        ! Functions
        contains
            ! Insert default values in type
            procedure:: init => TsimulationDefault
            ! Deallocate memory of type
            procedure:: destroy => TsimulationDeallocate
    end type

contains

! -----------------------------------------------------------------------------
! Set default values for Tsimulation type
pure subroutine TsimulationDefault( sim )
    implicit none
    !Entrada:
    class(Tsimulation), intent(inout) :: sim
    
    ! Default sets
    sim%end_time = 1.d0
    sim%max_it = -1
    
    
end subroutine TsimulationDefault
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Deallocate components in Tsimulation type
pure subroutine TsimulationDeallocate( sim )
    implicit none
    !Entrada:
    class(Tsimulation), intent(inout) :: sim
    
    ! No allocatable data in Tsimulation.
    ! Dummy function 
    sim%max_it = -1   
    
end subroutine TsimulationDeallocate
! -----------------------------------------------------------------------------

end module data_simulation