! #############################################################################
! This module defines all variables that use the in AMR structure. As well as
! lists all these variables in the functions needed by remesh
!   1. To add a variable first declare it with: `type(TamrVar), public:: myvar`
!   2. In `initAmrVariables` subroutine, add it's declaration like:
!      `call myvar%init(comp=1,ngc=2,sn=center)` where `comp` is the number of
!      components, `ngc` the number of ghosts in it and `sn` is the variable face
!      type that can be: `center`, `faceX`, `faceY`, `faceZ` and `nodal`.
!   3. Insert a line in `startAmrLevel` for `myvar` with:
!      `call AMRmesh%initCells(lev, t, dt, myvar, func)` where `lev` is the
!       level to apply the operation, `t` is the time, `dt` the time delta and
!       `func` is a function defined in **amr_startvariables** file following
!       the template there.
!   4. Insert a line for `myvar` in `destroyAmrVariables`, `buildAmrLevel`,
!      `remakeAmrLevel`, `interpAmrLevel` and `clearAmrLevel` by just
!       following the templates in each subroutine.
!   5. Add the variable to the print list if necessary. Follow the function 
!      template by calling `call out%set("My Variable Print Name",myvar)` then
!      add a `call var_list%append(out)` to insert this info to the list.
!
! #############################################################################

module amr_variables
    use data_variable
    use data_mesh
    use output_list
    use amrex_amr_module
    use amr_startvariables
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines are the only
    !ones accessible outside this module
    public:: initAmrVariables, destroyAmrVariables
    public:: startAmrLevel, buildAmrLevel, clearAmrLevel
    public:: remakeAmrLevel, interpAmrLevel
    public:: initAmrOutputs
    !---------------------------

    !The following variables are accessible outside this module
    type(TamrVar), target, public:: tag
    type(TamrVar), target, public:: phi, velu, velv, velw
    
contains

! -----------------------------------------------------------------------------
! This function holds the list of all AMR variables to allocate
subroutine initAmrVariables()
    implicit none
    !Local:
    
    ! Mesh topology variable
    call tag%init(comp=1,ngc=0,sn=center)
    ! True Variables
    call phi%init(comp=1,ngc=1,sn=center)
    call velu%init(comp=1,ngc=2,sn=faceX)
    call velv%init(comp=1,ngc=2,sn=faceY)
    call velw%init(comp=1,ngc=2,sn=faceZ)
    
end subroutine initAmrVariables
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! This function holds the list of all AMR variables do deallocate
subroutine destroyAmrVariables()
    implicit none
    !Local:
    
    ! Mesh topology variable
    call tag%Destroy()
    ! True Variables
    call phi%Destroy()
    call velu%Destroy()
    call velv%Destroy()
    call velw%Destroy()
    
end subroutine destroyAmrVariables
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Initialize the list of AMR variables in a mesh level during remesh
subroutine buildAmrLevel(lev, ba, dm)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    type(amrex_boxarray), intent(in) :: ba
    type(amrex_distromap), intent(in) :: dm
    !Local:

    ! Mesh topology variable
    call tag%levelInit(lev,ba,dm)
    ! True Variables
    call phi%levelInit(lev,ba,dm)
    call velu%levelInit(lev,ba,dm)
    call velv%levelInit(lev,ba,dm)
    call velw%levelInit(lev,ba,dm)
    
end subroutine buildAmrLevel
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Initialize the list of AMR variables in a mesh level during remesh
subroutine startAmrLevel(lev)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    !Local:
    double precision:: t, dt
    
    ! Set dummy variables
    t = 0.d0
    dt = 0.d0

    ! Mesh topology variable
    call AMRmesh%initCells(lev,t,dt,tag,zeros)
    ! True Variables
    call AMRmesh%initCells(lev,t,dt,phi,zeros)
    call AMRmesh%initCells(lev,t,dt,velu,zeros)
    call AMRmesh%initCells(lev,t,dt,velv,zeros)
    call AMRmesh%initCells(lev,t,dt,velw,zeros)
    
end subroutine startAmrLevel
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Interpolate the list of AMR variables between coarse to fine levels
subroutine interpAmrLevel(lev, ba, dm)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    type(amrex_boxarray), intent(in):: ba
    type(amrex_distromap), intent(in):: dm
    !Local:
    
    ! Mesh topology variable
    call AMRmesh%coarse2Fine(lev, tag, ba, dm)
    ! True Variables
    call AMRmesh%coarse2Fine(lev, phi, ba, dm)
    call AMRmesh%coarse2Fine(lev, velu, ba, dm)
    call AMRmesh%coarse2Fine(lev, velv, ba, dm)
    call AMRmesh%coarse2Fine(lev, velw, ba, dm)
    
end subroutine interpAmrLevel
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Finish interpolating procedure on the list of AMR variables.
! Joins the older fine cells with the ones interpolated from coarser level.
subroutine remakeAmrLevel(lev, ba, dm)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    type(amrex_boxarray), intent(in):: ba
    type(amrex_distromap), intent(in):: dm
    !Local:
    
    ! Mesh topology variable
    call AMRmesh%levelRemake(lev, tag, ba, dm)
    ! True Variables
    call AMRmesh%levelRemake(lev, phi, ba, dm)
    call AMRmesh%levelRemake(lev, velu, ba, dm)
    call AMRmesh%levelRemake(lev, velv, ba, dm)
    call AMRmesh%levelRemake(lev, velw, ba, dm)
    
end subroutine remakeAmrLevel
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Destroy the list of AMR variables in a mesh level during remesh
subroutine clearAmrLevel(lev)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    !Local:
    
    ! Mesh topology variable
    call tag%levelDestroy(lev)
    ! True Variables
    call phi%levelDestroy(lev)
    call velu%levelDestroy(lev)
    call velv%levelDestroy(lev)
    call velw%levelDestroy(lev)
    
end subroutine clearAmrLevel
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Initialize the list of variables to be available at output
subroutine initAmrOutputs()
    implicit none
    !Local:
    type(Toutput):: out

    ! Initialize list - Don't touch this
    call var_list%init()

    ! Mesh topology variable
    call out%set("Mesh Tagging",tag); call var_list%append(out)
    ! True Variables
    call out%set("Scalar",phi); call var_list%append(out)
    call out%set("U",velu); call var_list%append(out)
    call out%set("V",velv); call var_list%append(out)
    call out%set("W",velw); call var_list%append(out)
    
end subroutine initAmrOutputs
! -----------------------------------------------------------------------------

end module amr_variables