! #############################################################################
! This module holds the functions needed to set and perform mesh changes
! using AMReX library
! #############################################################################
module amr_remesh
    use data_variable
    use data_mesh
    use amr_variables
    use amr_criteria
    use amrex_amr_module
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: initAmrMesh
    !-------------------------------------
    
contains

! -----------------------------------------------------------------------------
! Set the functions needed by AMReX mesh operations
subroutine initAmrMesh()
    implicit none
    !Local:

    ! Set the function pointers in AMReX to our functions
    call amrex_init_virtual_functions( &
        newLevelFromScratch, &
        newLevelFromCoarse, &
        remakeLevel, &
        clearLevel, &
        buildCriteria &
    )

    ! Set the function pointer in AMReX to be executed after regrid
    call amrex_init_post_regrid_function( &
        afterRemesh &
    )
    
    
end subroutine initAmrMesh
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Starts a level from nothing (initial condition)
! Only used during initialization.
subroutine newLevelFromScratch(Clev, time, pba, pdm) bind(c)
    implicit none
    !Entrada:
    integer(c_int), intent(in), value :: Clev
    real(c_double), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm
    !Local:
    integer:: lev
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm

    ! Translate C-index to Fortran
    lev = Clev+1
    ! Set C pointers to have a defined type
    ba = pba
    dm = pdm
    
    ! Reallocate all AMR variables in current level
    call buildAmrLevel(lev, ba, dm)
    ! Start all AMR variables with initial condition
    call startAmrLevel(lev)

end subroutine newLevelFromScratch
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Make a new level and fill with interpolated data from coarser level.
subroutine newLevelFromCoarse(Clev, time, pba, pdm) bind(c)
    implicit none
    !Entrada:
    integer(c_int), intent(in), value :: Clev
    real(c_double), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm
    !Local:
    integer:: lev
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm

    ! Translate C-index to Fortran
    lev = Clev+1
    ! Set C pointers to have a defined type
    ba = pba
    dm = pdm
    
    ! Interpolate from coarser level
    call interpAmrLevel(lev,ba,dm)

end subroutine newLevelFromCoarse
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Remake an existing level and fill with data from existing fine and 
! coarse levels.
subroutine remakeLevel(Clev, time, pba, pdm) bind(c)
    implicit none
    !Entrada:
    integer(c_int), intent(in), value :: Clev
    real(c_double), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm
    !Local:
    integer:: lev
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm

    ! Translate C-index to Fortran
    lev = Clev+1
    ! Set C pointers to have a defined type
    ba = pba
    dm = pdm
    
    ! Rebuild the whole level using old mesh and coarse info
    call remakeAmrLevel(lev,ba,dm)

end subroutine remakeLevel
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Mark cell to refine/coarse
subroutine buildCriteria(Clev, cp, time, settag, cleartag) bind(c)
    implicit none
    !Entrada:
    integer(c_int), intent(in), value :: Clev
    real(c_double), intent(in), value :: time
    type(c_ptr), intent(in), value :: cp
    character(kind=c_char), intent(in), value :: settag, cleartag
    !Local:
    integer:: lev
    type(amrex_tagboxarray):: amrex_tag

    ! Translate C-index to Fortran
    lev = Clev+1
    ! Set C pointers to have a defined type
    amrex_tag = cp
    
    ! Choose a criterion
    if (time<=0.d0) then
        ! Initial mesh
        call AMRmesh%initCells(lev, time, 0.d0, tag, initial_tag)
    else

        !TODO: Insert other mesh criteria here
        ! using the box and forCells function
    
    end if

    ! Set tagged values to AMReX system
    call markTagOnAMReX(lev,tag%mfab(lev),amrex_tag,settag,cleartag)

end subroutine buildCriteria
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Translate tagged values to AMReX structure
subroutine markTagOnAMReX(lev, tagvar, amrex_tag, settag, cleartag)
    implicit none
    !Entrada:
    integer, intent(in):: lev
    type(amrex_multifab), intent(in):: tagvar
    type(amrex_tagboxarray), intent(inout):: amrex_tag
    character(kind=c_char), intent(in):: settag, cleartag
    !Local:
    integer:: Clev
    type(Tpatch):: patch
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    double precision, contiguous, pointer :: tagarr(:,:,:,:)
    character(kind=c_char), contiguous, pointer :: amrex_tagarr(:,:,:,:)
    
    ! Translate to C-index
    Clev = lev-1

    ! Custom patch loop oppen to set tag variable into amrex structure
    ! Get iterator
    call amrex_mfiter_build(mfi, tagvar, tiling=.true.)
    ! Loop in patches
    do while(mfi%next())

        ! Get limits
        bx = mfi%tilebox()

        ! Set patch info type
        patch%ilo = bx%lo
        patch%ihi = bx%hi
        patch%lev = lev
        ! Watchout for C like index in amrex_geom
        patch%dn = amrex_geom(Clev)%dx

        ! Open both tag variables
        tagarr => tagvar%dataptr(mfi)
        amrex_tagarr => amrex_tag%dataptr(mfi)

        ! Call translation funtion
        call tag2AMReX(patch, settag, cleartag, &
            tagarr, lbound(tagarr), ubound(tagarr), &
            amrex_tagarr, lbound(amrex_tagarr), ubound(amrex_tagarr) )
    end do

    ! Deallocate Iterator
    call amrex_mfiter_destroy(mfi)
    
end subroutine markTagOnAMReX
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Deallocate AMR variables in the level
subroutine clearLevel(Clev) bind(c)
    implicit none
    !Entrada:
    integer(c_int), intent(in), value :: Clev
    !Local:
    integer:: lev

    ! Translate C-index to Fortran
    lev = Clev+1

    ! Deallocate all AMR variables in current level
    call clearAmrLevel(lev)
    
end subroutine clearLevel
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! This is executed automatically after regrid
subroutine afterRemesh()
    implicit none
    !Local:
    
    ! Mesh type Update
    call AMRmesh%updateAttributes()
    
end subroutine afterRemesh
! -----------------------------------------------------------------------------

end module amr_remesh