! #############################################################################
! This module contains functions to operate/convert indexes in several ways
!
! #############################################################################
module amr_coords
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: coord2index, index2coord
    !-------------------------------------
    
contains

! -----------------------------------------------------------------------------
! Transform from coodinate (x,y,z) to amr-index
! Old amr3d: index = idnint((x-gax)/dx+0.5d0*dble(1+sx))
elemental function coord2index(problo,x,dx,faced) result(ix)
    implicit none
    !Entrada:
    integer, intent(in):: faced
    double precision, intent(in):: problo, x, dx
    !Saída:
    integer:: ix
    !Local:

    ix = nint( (x - problo - (1-faced)*0.5d0*dx)/dx )

end function coord2index
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Transform from amr-index to coodinate (x,y,z)
! Old amr3d: coord = gax + (dble(ix)-0.5d0*dble(1+sx))*dx
elemental function index2coord(problo,ix,dx,faced) result(x)
    implicit none
    !Entrada:
    integer, intent(in):: ix, faced
    double precision, intent(in):: problo, dx
    !Saída:
    double precision:: x
    !Local:
    
    x = problo + dble(ix)*dx + (1-faced)*0.5d0*dx

end function index2coord
! -----------------------------------------------------------------------------


end module amr_coords