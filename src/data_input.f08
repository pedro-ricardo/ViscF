! #############################################################################
! This module stores the input parameters to MFSim software.
!
!   To add a new input parameters create a type like `Tmyclass`
!   and declare your parameters in it.
!   Every new input type must have:
!       * A `init()` pure subroutine that initialize a default
!         value to the parameters.
!       * A `destroy()` pure subroutine that deallocate the memory
!         of allocatable parameters. This might be a dummy subroutine.
!   Your type needs to be added than to `Tinput` as a parameter. And
!   the `init()` and `destroy()` methods of it must be called inside
!   `Tinput_default()` and `Tinput_deallocate` subroutines. 
!       
! #############################################################################
module data_input
    use data_io
    use data_mesh
    use data_remesh
    use data_simulation
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: input
    !-------------------------------------

    ! Type with all input parameters in it
    type Tinput
        ! Mesh parameters
        type(Tmeshparams):: msh
        ! Remesh parameters
        type(Tremesh):: remsh
        ! Simulation parameters
        type(Tsimulation):: sim
        ! Input/Output parameters
        type(Tio):: io
        ! Functions
        contains
            ! Insert default values in type
            procedure:: init => TinputDefault
            ! Deallocate memory of type
            procedure:: destroy => TinputDeallocate
    end type


    ! Module variable
    type(Tinput):: input

contains


! -----------------------------------------------------------------------------
! Deallocate components in Tinput type
pure subroutine TinputDeallocate( input )
    implicit none
    !Entrada:
    class(Tinput), intent(inout) :: input

    ! Mesh data
    call input%msh%destroy()
    ! Remesh data
    call input%remsh%destroy()
    ! Simulation data
    call input%sim%destroy()
    ! Input/Output data
    call input%io%destroy()

end subroutine TinputDeallocate
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Initialize Tinput type with default values for it
pure subroutine TinputDefault( input )
    implicit none
    !Entrada:
    class(Tinput), intent(inout) :: input
    
    ! Mesh data
    call input%msh%init()
    ! Remesh data
    call input%remsh%init(input%msh%levels)
    ! Simulation data
    call input%sim%init()
    ! Input/Output data
    call input%io%init()
    
end subroutine TinputDefault
! -----------------------------------------------------------------------------

end module data_input