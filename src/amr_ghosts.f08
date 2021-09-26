! #############################################################################
! This module holds subroutine/functions to calculate the patch ghosts cells.
! Can fill Sisters/Internal, Coarse Fine and Domain Boundary.
! #############################################################################
module amr_ghosts
    use data_variable
    use amrex_amr_module

    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: domainBC, fillSisterBoundary, fillAllBoundaries
    public:: VarBCs
    !-------------------------------------
    
    ! Type to store boundary information
    type TboundaryType
        ! Pointer data to BC array in the variable
        integer, pointer:: val(:,:,:) => Null()
        ! Number of components in variable
        integer:: ncomp = 1
        ! Functions
        contains
            ! Set the BC for a variable
            procedure:: set => setTboundaryType
            ! Clear the variable
            procedure:: clear => clearTboundaryType
            ! Get the min faces BC types
            procedure:: min => minTboundaryType
            ! Get the max faces BC types
            procedure:: max => maxTboundaryType
    end type TboundaryType

    type(TboundaryType):: VarBCs

contains

! -----------------------------------------------------------------------------
! Set boundary pointer to variable BCs
subroutine setTboundaryType(bt,var)
    implicit none
    !Entrada:
    class(TboundaryType), intent(inout):: bt
    type(TamrVar), target, intent(in):: var
    !Local:

    ! Save components
    bt%ncomp = var%components
    ! Set pointer
    bt%val => var%bc
    
end subroutine setTboundaryType
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Clear boundary pointer
subroutine clearTboundaryType(bt)
    implicit none
    !Entrada:
    class(TboundaryType), intent(inout):: bt
    !Local:
    
    ! Clear components
    bt%ncomp = 1
    ! Clear pointer
    nullify(bt%val)
    
end subroutine clearTboundaryType
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Return the minimun of boundary types (west,south,bottom)
function minTboundaryType(bt) result(minbc)
    implicit none
    !Entrada:
    class(TboundaryType), intent(inout):: bt
    !Saida:
    integer:: minbc(3,bt%ncomp)
    !Local:
    
    ! Default type is Neumann
    minbc = amrex_bc_reflect_even
    ! Check pointer to get real BC info
    if (associated(bt%val)) minbc = bt%val(:,:,1)
    
end function minTboundaryType
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Return the maximun of boundary types (east,north,top)
function maxTboundaryType(bt) result(maxbc)
    implicit none
    !Entrada:
    class(TboundaryType), intent(inout):: bt
    !Saida:
    integer:: maxbc(3,bt%ncomp)
    !Local:
    
    ! Default type is Neumann
    maxbc = amrex_bc_reflect_even
    ! Check pointer to get real BC info
    if (associated(bt%val)) maxbc = bt%val(:,:,2)
    
    
end function maxTboundaryType
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! This function applies the domain boundary conditions to a variable
! It is a `amrex_physbc_proc` procedure needed by interpolation functions.
subroutine domainBC(Cvar, src_comp, ncomp, time, Cgeom) bind(c)
    ! Entrada:
    type(c_ptr), value:: Cvar, Cgeom
    integer(c_int), value:: src_comp, ncomp
    real(c_double), value:: time
    ! Local:
    type(amrex_geometry):: geom
    type(amrex_multifab):: var
    type(amrex_mfiter):: mfi
    double precision, contiguous, pointer:: v(:,:,:,:)
    integer:: lo(4), hi(4)
    
    ! Check if problem isn't fully periodic
    if (.not. amrex_is_all_periodic()) then

        ! Translate C types to Fortran
        geom = Cgeom
        var = Cvar

        ! Start iterator (without tiling for performance boost in ghost sweep)
        call amrex_mfiter_build(mfi, var, tiling=.false.)
        do while(mfi%next())
            ! Get patch
            v => var%dataptr(mfi)
            ! Check if it touches the domain
            if (.not. geom%domain%contains(v)) then
                ! Get patch limits
                lo = lbound(v)
                hi = ubound(v)

                ! Get boundary function for standard AMReX types
                ! VarBCs is a Module variable that need to be set
                ! othewise it will be considered Neumann
                call amrex_filcc(v, lo, hi, &
                    geom%domain%lo, geom%domain%hi, geom%dx, &
                    geom%get_physical_location(lo), VarBCs%min(), VarBCs%max() )

                ! Get boundary function for custom dirichlet type
                ! (...)

            end if
        end do

    end if

end subroutine domainBC
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Inject sister patches ghosts, including the periodic at same level.
subroutine fillSisterBoundary(lev,var)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    type(TamrVar), intent(inout):: var
    !Local:
    
    ! Call sister injection
    call var%mfab(lev)%fill_boundary(amrex_geom(lev-1))
    
end subroutine fillSisterBoundary
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Apply ghosts calculations to sister, coarse/fine and domain boundary cells.
! Equal `var` and `srcvar` to fill `var` own ghosts.
subroutine fillAllBoundaries(lev, var, srcvar, interp)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    type(TamrVar), intent(inout):: var
    type(TamrVar), intent(in):: srcvar
    integer, intent(in), optional:: interp
    !Local:
    integer:: Clev, int, c
    double precision:: t, t0
    integer:: def_lo_bc(3,srcvar%components), def_hi_bc(3,srcvar%components)
    
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
    call VarBCs%set(srcvar)

    ! Check level in variable
    if (lev <= size(var%mfab) .and. lev <= size(srcvar%mfab)) then

        if (lev == 1) then
            ! Single level ghosts
            do c = 1,srcvar%components
                call amrex_fillpatch(var%mfab(lev), t0, srcvar%mfab(lev), t, srcvar%mfab(lev), &
                    amrex_geom(Clev), domainBC, t, c, c, srcvar%components)
            end do
        else
            ! Multiple level ghosts
            do c = 1,srcvar%components
                call amrex_fillpatch(var%mfab(lev), t0, srcvar%mfab(lev-1), t, srcvar%mfab(lev-1), &
                    amrex_geom(Clev-1), domainBC, t0, srcvar%mfab(lev), t, &
                    srcvar%mfab(lev), amrex_geom(Clev), domainBC, t, c, c, &
                    srcvar%components, amrex_ref_ratio(Clev-1), int, def_lo_bc, def_hi_bc)
            end do
        end if

    end if
    
end subroutine fillAllBoundaries
! -----------------------------------------------------------------------------

end module amr_ghosts