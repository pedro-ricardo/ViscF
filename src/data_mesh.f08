! #############################################################################
! This module stores the Mesh parameters to MFSim software.
!
! #############################################################################
module data_mesh
    use mpi_f08
    use data_boxes
    use data_variable
    use amr_ghosts
    use amrex_amr_module
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: Tmesh, Tmeshparams, op_in_patch
    public:: AMRmesh
    !Variables from other modules made accesible from here
    public:: Tbox, Tpatchbox, TamrVar
    !-------------------------------------
    
    ! Type for mesh parameters
    type Tmeshparams
        ! Number of cells in each direction
        integer:: bot_cells(3)
        ! Maximun number of refinement levels
        integer:: levels
        ! Refinement Ratio
        integer, allocatable:: ratio(:)
        ! Maximun Lbot patch size - default is 32
        integer:: max_grid_size
        ! Domain Geometry
        double precision:: dom_lo(3), dom_hi(3)
        ! Periodic Directions
        logical:: periodic(3)
        ! Initial patches in fine level
        type(Tpatch), allocatable:: init_grid(:)
        ! Functions
        contains
            ! Insert default values in type
            procedure:: init => TmeshparamsDefault
            ! Deallocate memory of type
            procedure:: destroy => TmeshparamsDeallocate
    end type

    ! Type for mesh operations
    type Tmesh
        ! Mesh input Parameters
        type(Tmeshparams), pointer:: params
        ! Number of levels in the mesh
        integer:: nlev
        ! Number of visible cells in the mesh
        integer(8):: ncell
        ! Mesh attributes
        type(amrex_geometry), allocatable:: geo(:)
        type(amrex_boxarray), allocatable:: ba(:)
        type(amrex_distromap), allocatable:: dm(:)
        ! Functions
        contains
            ! Insert default values in type
            procedure:: init => TmeshSetParams
            ! Clear mesh characteristics
            procedure:: updateAttributes => TmeshUpdateAttributes
            ! Deallocate memory of type
            procedure:: destroy => TmeshDeallocate
            ! Iterate over mesh cells
            procedure:: forCells => TmeshForCells
            ! Iterate Initialize value in mesh cells
            procedure:: initCells => TmeshInitCells
            ! Interpolate the coarse cells into fine ones
            procedure:: coarse2Fine => TmeshCoarse2Fine
            ! Interpolate fine cells into coarse ones
            procedure:: fine2Coarse => TmeshFine2Coarse
            ! Interpolate center variable to 3 face variables
            procedure:: center2Face => TmeshCenter2Face
            ! Level remake
            procedure:: levelRemake => TmeshLevelRemake
    end type

    ! User subroutine template in forCells
    interface
        subroutine op_in_patch(t,dt,patch,ngcl,box)
            import Tpatch, Tpatchbox
            implicit none
            double precision, intent(in):: t, dt
            integer, intent(in):: ngcl
            type(Tpatch), intent(in):: patch
            type(Tpatchbox), target:: box(:)
        end subroutine
    end interface

    ! Global mesh variable
    type(Tmesh):: AMRmesh

contains

! -----------------------------------------------------------------------------
! Deallocate components in Tmeshparams type
pure subroutine TmeshparamsDeallocate( msh )
    implicit none
    !Entrada:
    class(Tmeshparams), intent(inout) :: msh
    
    if (allocated(msh%ratio)) then
        deallocate(msh%ratio)
    end if
    if (allocated(msh%init_grid)) then
        deallocate(msh%init_grid)
    end if
    
end subroutine TmeshparamsDeallocate
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Initialize Mesh parameters default values
pure subroutine TmeshparamsDefault( msh )
    implicit none
    !Entrada:
    class(Tmeshparams), intent(inout) :: msh

    ! Default sets
    msh%bot_cells = [32, 32, 32]
    msh%levels = 0
    allocate(msh%ratio(max(msh%levels,1)))
    msh%ratio(:) = 2
    msh%max_grid_size = 32
    msh%dom_hi = [1.d0, 1.d0, 1.d0]
    msh%dom_lo = [0.d0, 0.d0, 0.d0]
    allocate(msh%init_grid(1))
    msh%init_grid(1)%hi=[1.d0, 1.d0, 1.d0]
    msh%init_grid(1)%lo=[0.d0, 0.d0, 0.d0]
    msh%periodic = [.true.,.true.,.true.]
    
end subroutine TmeshparamsDefault
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Initialize Tmesh type with values
subroutine TmeshSetParams(msh,mparam)
    implicit none
    !Entrada:
    class(Tmesh), intent(inout):: msh
    type(Tmeshparams), target, intent(in):: mparam
    !Local:
    
    ! Set pointer to input mesh parameters
    msh%params => mparam
    ! Set number of levels
    msh%nlev = msh%params%levels+1
    ! Set number of cell initially
    msh%ncell = 0

end subroutine TmeshSetParams
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Deallocate components in Tmesh type
subroutine TmeshDeallocate(msh)
    implicit none
    !Entrada:
    class(Tmesh), intent(inout):: msh
    !Local:
    integer:: i

    ! Get rid of pointer to parameters
    nullify(msh%params)
    ! Deallocate BoxArray
    if (allocated(msh%ba)) then
        do i=1,size(msh%ba)
            call amrex_boxarray_destroy(msh%ba(i))
        end do
        deallocate(msh%ba)
    end if
    ! Deallocate distribution mapping
    if (allocated(msh%dm)) then
        do i=1,size(msh%dm)
            call amrex_distromap_destroy(msh%dm(i))
        end do
        deallocate(msh%dm)
    end if
    ! Deallocate Geometry
    if (allocated(msh%geo)) then
        do i=1,size(msh%geo)
            call amrex_geometry_destroy(msh%geo(i))
        end do
        deallocate(msh%geo)
    end if
    ! Set levels to zero
    msh%nlev = 0
    
end subroutine TmeshDeallocate
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Routine executed inside AMR mesh.
! Loop though all patches in level and calls the user operation
subroutine TmeshForCells(msh, lev, t, dt, ngc, user_function, box)
    implicit none
    !Entrada:
    class(Tmesh), intent(inout) :: msh
    integer, intent(in):: lev, ngc
    double precision, intent(in):: t, dt
    procedure(op_in_patch):: user_function
    type(Tbox), target, intent(inout):: box
    !Local:
    integer:: v
    integer:: vlo(4), vhi(4)
    type(Tpatch):: patch
    type(amrex_box):: bx
    type(Tpatchbox), pointer:: VarPatch
    type(amrex_mfiter) :: mfi
    
    ! Perform integrity checks
    if (box%nvar==0) then
        write(*,*)'Warning: Empty variable box passed to "forCells".'
        ! Do Nothing
        return
    end if
    
    ! Build the iterator for patches
    call amrex_mfiter_build(mfi, box%Vars(1)%var%mfab(lev), tiling=.true.)

    ! Iterate over patches
    do while(mfi%next())

        ! Get first box
        bx = mfi%tilebox()
        ! Ensure a centered
        call bx%cellize(1)  ! weird name but it
        call bx%cellize(2)  ! sets a center index
        call bx%cellize(3)  ! on each direction

        ! Ensure only phisical cells in patch box
        call bx%grow( -box%Vars(1)%var%mfab(lev)%nghost() )

        ! Important Information:
        ! Center index will always be in `ilo` and `ihi` variables
        !
        !     i      i+1     i+2      -> Center indexes (default loop)
        ! |-------|-------|-------|
        ! |   o   |   o   |   o   |
        ! |-------|-------|-------|
        ! i      i+1     i+2     i+3  -> Equivalent Face index
        !
        
        ! Set patch info type
        patch%ilo = bx%lo
        patch%ihi = bx%hi
        patch%lev = lev
        ! Watchout for C like index in amrex_geom
        patch%dn = amrex_geom(lev-1)%dx

        ! Loop through variables in box
        do v=1,box%nvar
            ! Just an alias
            VarPatch => box%Patch(v)
            ! Set pointer to it's data
            ! The variable patch has 4 indexes (i,i,k,component)
            VarPatch%var => box%Vars(v)%var%mfab(lev)%dataPtr(mfi)
            ! Get variable index limits
            vlo = lbound(VarPatch%var)
            vhi = ubound(VarPatch%var)
            ! Separate 3D index from components
            VarPatch%ilim(:,1) = vlo(1:3)
            VarPatch%ilim(:,2) = vhi(1:3)
            VarPatch%ncomp = vhi(4)
            ! Copy variable type (centered/faced)
            VarPatch%sn = box%Vars(v)%var%face_info

        end do
        ! Unset alias
        nullify(VarPatch)

        ! Call the user defined function in current patch
        call user_function(t,dt,patch,ngc,box%Patch)

    end do

    ! Deallocate iterator
    call amrex_mfiter_destroy(mfi)

end subroutine TmeshForCells
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Routine executed inside AMR mesh. For initialization purposes
! Loop though all patches in level and calls the user operation for a
! single variable
subroutine TmeshInitCells(msh, lev, t, dt, var, user_function)
    implicit none
    !Entrada:
    class(Tmesh), intent(inout) :: msh
    integer:: lev
    double precision:: t, dt
    type(TamrVar), intent(inout):: var 
    procedure(op_in_patch):: user_function
    !Local:
    type(Tbox):: box

    ! Encapsulate variable in box
    call box%addVar(var)
    ! Call forCells to use user function through all cells
    call msh%forCells(lev, t, dt, var%ghosts, user_function, box)
    
end subroutine TmeshInitCells
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Interpolate data from the level below. Can be used in
! the same mesh or a new one.
subroutine TmeshCoarse2Fine(msh, lev, var, ba, dm, interp)
    implicit none
    !Entrada:
    class(Tmesh), intent(inout):: msh
    type(TamrVar), intent(inout):: var
    integer, intent(in):: lev
    integer, optional, intent(in):: interp
    type(amrex_boxarray), optional, intent(in):: ba
    type(amrex_distromap), optional, intent(in):: dm
    !Local:
    integer:: Clev, int, c
    double precision:: t, t0
    integer:: def_lo_bc(3,var%components), def_hi_bc(3,var%components)
    
    ! Interpolation scheme
    int = amrex_interp_cell_cons
    if (present(interp)) int = interp

    ! Dummy variables to bypass am AMReX time interpolation
    t0= 0.d0
    t = t0+1.d0
    ! Level to C index
    Clev = lev-1
    ! Default BCs (Neumann)
    def_lo_bc = amrex_bc_reflect_even
    def_hi_bc = amrex_bc_reflect_even
    ! Set module variable to be used inside `domainBC`
    call VarBCs%set(var)

    ! Remake level before interpolation
    if (present(ba).and.present(dm)) then
        call var%levelInit(lev,ba,dm)
    end if

    ! Interpolate from coarse
    do c = 1,var%components
        call amrex_fillcoarsepatch(var%mfab(lev), t0, var%mfab(lev-1), t, var%mfab(lev-1), &
            amrex_geom(Clev-1), domainBC, amrex_geom(Clev), domainBC, t, c, c, var%components,&
            amrex_ref_ratio(Clev-1), int, def_lo_bc, def_hi_bc)
    end do
    
end subroutine TmeshCoarse2Fine
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Average data in all components from level `lev+1` to lev.
! Can perform an average and store it in a different variable if needed.
subroutine TmeshFine2Coarse(msh, lev, var, var2)
    implicit none
    !Entrada:
    class(Tmesh), intent(inout):: msh
    type(TamrVar), target, intent(inout):: var
    integer, intent(in):: lev
    !Optional:
    type(TamrVar), optional, target, intent(inout):: var2
    !Local:
    integer:: Clev, c
    type(TamrVar), pointer:: average_var => Null()

    ! Level to C index
    Clev = lev-1

    ! Set variable to store the average
    ! Default is the same variable
    average_var => var
    if (present(var2)) then
        average_var => var2
    end if

    ! Check for level above
    if (lev<amrex_get_numlevels()) then
        ! Average all components
        do c = 1,var%components
            call amrex_average_down(var%mfab(lev+1), var2%mfab(lev), amrex_geom(Clev+1), &
                amrex_geom(Clev), c, var%components, amrex_ref_ratio(Clev))
        end do
    end if
    
end subroutine TmeshFine2Coarse
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Average a CellCenter variable to 3 face variables
subroutine TmeshCenter2Face(msh, lev, var, varX, varY, varZ)
    implicit none
    !Entrada:
    class(Tmesh), intent(inout):: msh
    type(TamrVar), target, intent(in):: var
    type(TamrVar), target, optional, intent(inout):: varX
    type(TamrVar), target, optional, intent(inout):: varY
    type(TamrVar), target, optional, intent(inout):: varZ
    integer, intent(in):: lev
    !Local:
    integer:: Clev, c
    type(amrex_multifab), allocatable:: faces(:,:)

    ! Level to C index
    Clev = lev-1

    ! Allocate faces variable
    allocate(faces(3,size(var%mfab)))

    ! Copy variable to faces
    if (present(varX)) faces(1,:) = varX%mfab
    if (present(varY)) faces(2,:) = varY%mfab
    if (present(varZ)) faces(3,:) = varZ%mfab

    ! Check for level above
    if (lev<amrex_get_numlevels()) then
        ! Average all components
        do c = 1,var%components
            call amrex_average_cellcenter_to_face(faces(:,lev), &
                var%mfab(lev), amrex_geom(Clev))
        end do
    end if

    ! Set values back
    if (present(varX)) varX%mfab = faces(1,:)
    if (present(varY)) varY%mfab = faces(2,:)
    if (present(varZ)) varZ%mfab = faces(3,:)
    ! Free faces variable
    deallocate(faces)
    
end subroutine TmeshCenter2Face
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Rebuild a fine level with given level and mesh info
subroutine TmeshLevelRemake(msh, lev, var, ba, dm, interp)
    implicit none
    !Entrada:
    class(Tmesh), intent(inout):: msh
    type(TamrVar), intent(inout):: var
    integer, intent(in):: lev
    type(amrex_boxarray), intent(in):: ba
    type(amrex_distromap), intent(in):: dm
    integer, optional, intent(in):: interp
    !Local:
    type(TamrVar):: var_tmp
    integer:: Clev, int
    
    ! Interpolation scheme
    int = amrex_interp_cell_cons
    if (present(interp)) int = interp

    ! Create a temp variable with the same kind
    call var_tmp%init(comp=var%components,ngc=var%ghosts,sn=var%face_info)
    ! Create level in temp variable
    call var_tmp%levelInit(lev,ba,dm)

    ! Translate C-index
    Clev = lev-1

    ! Fill data to temporary variables
    call fillAllBoundaries(lev, var_tmp, var, int)
    ! Clear the level
    call var%levelDestroy(lev)
    ! Move temporary to variable
    call var%mfab(lev)%move(var_tmp%mfab(lev))
    ! Delete temporary
    call var_tmp%destroy()

end subroutine TmeshLevelRemake
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Update mesh characteristics
subroutine TmeshUpdateAttributes(msh)
    implicit none
    !Entrada:
    class(Tmesh), intent(inout):: msh
    !Local:
    integer:: l, nlev, Clev, i, ierr
    integer(8):: cells
    type(TamrVar):: dummy
    
    ! Get number of levels
    nlev = amrex_get_numlevels()

    ! Remove old BoxArray
    if (allocated(msh%ba)) then
        do i=1,size(msh%ba)
            call amrex_boxarray_destroy(msh%ba(i))
        end do
        deallocate(msh%ba)
    end if
    ! Remove old Distribution Map
    if (allocated(msh%dm)) then
        do i=1,size(msh%dm)
            call amrex_distromap_destroy(msh%dm(i))
        end do
        deallocate(msh%dm)
    end if
    ! Remove old Geometry
    if (allocated(msh%geo)) then
        do i=1,size(msh%geo)
            call amrex_geometry_destroy(msh%geo(i))
        end do
        deallocate(msh%geo)
    end if

    ! Allocate new
    allocate(msh%dm(nlev),msh%ba(nlev),msh%geo(nlev))

    ! Get new data from AMReX
    do l=1,nlev
        Clev = l-1
        msh%dm(l) = amrex_get_distromap(Clev)
        msh%ba(l) = amrex_get_boxarray(Clev)
        msh%geo(l) = amrex_get_geometry(Clev)
    end do

    ! Update number of cells
    msh%ncell = 0
    ! Create a dummy variable
    call dummy%init(1,0,center)
    ! Build it
    do l = 1,msh%nlev
        call dummy%levelInit(l,msh%ba(l),msh%dm(l))
        ! Calculate cells on each level
        call msh%initCells(l,0.d0,0.d0,dummy,countCells)
    end do
    ! Deallocate dummy
    call dummy%destroy()
    ! Sum cells in all processes
    cells = msh%ncell
    call MPI_Allreduce(cells, msh%ncell, 1, MPI_Integer8, MPI_Sum, MPI_COMM_WORLD, ierr)

    ! Update number of levels
    msh%nlev = nlev

end subroutine TmeshUpdateAttributes
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Go through patches counting the number of visible cells
subroutine countCells(t,dt,patch,ngcl,Pvar)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    integer, intent(in):: ngcl
    type(Tpatch), intent(in):: patch
    type(Tpatchbox), target:: Pvar(:)
    !Local:
    integer:: cells(3)
    
    ! Get number of cells in each direction
    cells = patch%ihi-patch%ilo+1
    ! Sums in the ncell global type
    AMRmesh%ncell = AMRmesh%ncell + (cells(1)*cells(2)*cells(3))
    
end subroutine countCells
! -----------------------------------------------------------------------------

end module data_mesh