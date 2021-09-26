! #############################################################################
! This module has a list for the pairs name/var needed to print. This list 
! template follows the "Modern Fortran Heterogeneous Linked List" implementation
! at GitHub 
!
! Check the original at: https://github.com/pedro-ricardo/LinkedList
!
! #############################################################################
module output_list
    use output_names

implicit none

!#####################################
!            Module Usage:
!#####################################
private
  
!The following subroutines/variables are the only
!ones accessible outside this module
public:: Tnode, Tlist
!Variables from other modules made accesible from here
public:: Toutput, var_list
!-------------------------------------

!-------------------------------------
type Tnode
    type(Tnode), pointer :: next => null()
    type(Tnode), pointer :: prev => null()
    type(Toutput), allocatable :: item
    contains
    procedure, private :: destroy => node_finalizer
    procedure, private :: destroy_all => node_finalizer_snowball
end type Tnode
!-------------------------------------

!-------------------------------------
type Tlist
    integer, private :: num_items = 0
    type(Tnode), pointer :: head => null()
    type(Tnode), pointer :: tail => null()
    contains
    procedure:: init => start_list
    procedure:: append => list_append_item
    procedure:: destroy => list_finalizer
    procedure:: foreach => list_foreach
    procedure:: pop_node => list_pop_node_n
    procedure:: pop_this => list_pop_node
    procedure:: len => list_lenght
    procedure:: len_print => list_lenght_printable
    final:: destroy_list
end type Tlist
!-------------------------------------

!-------------------------------------
! Definition of a subroutine that operates on node item
interface
subroutine process(item)
    import Toutput
    implicit none
    !Entrada:
    type(Toutput), intent(inout) :: item
    
end subroutine process
end interface
!-------------------------------------


! Module variable: List of variables to print
type(Tlist):: var_list

contains

! -----------------------------------------------------------------------------
! Retuns a initialized node. Can have an object in it or be empty.
pure function start_node( item ) result( val )
    implicit none
    type(Tnode) :: val
    !Entrada:
    type(Toutput), intent(in), optional :: item
    !Local:
    
    !Insert a item to the node if present
    if (present(item)) allocate(val%item, source=item)
    
end function start_node
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Returns a initialized list. Seems useless at this point due to standard values
! present on declarations.
pure subroutine start_list( this_list )
    implicit none
    !Entrada:
    class(Tlist), intent(inout) :: this_list
    !Local:

    this_list%num_items = 0

end subroutine start_list
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Delete all nodes in sequence from the list and frees the memory in the items.
pure recursive subroutine node_finalizer_snowball( this_node )
    implicit none
    !Entrada:
    class(Tnode), intent(inout) :: this_node
    !Local:
    
    !Deallocate it's item
    if (allocated(this_node%item)) deallocate(this_node%item)
    !Nullify it's pointers
    if (associated(this_node%next)) then
        call this_node%next%destroy_all()
        deallocate(this_node%next)
        nullify(this_node%next)
    end if
    nullify(this_node%prev)
    
end subroutine node_finalizer_snowball
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Pop out a node from the list, by a given number.
pure subroutine list_pop_node_n( this_list, node_num )
    implicit none
    !Entrada:
    class(Tlist), intent(inout) :: this_list
    integer, intent(in):: node_num
    !Local:
    type(Tnode), pointer:: curr
    integer:: cont

    !Foward sweep
    curr => this_list%head
    cont = 1
    do while ( associated(curr) )
        if (cont==node_num) then
            
            call this_list%pop_this(curr)
            
            !Exit when found
            return
        end if
        curr => curr%next
        cont = cont+1
    end do
    
    
end subroutine list_pop_node_n
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Pop out a node from the list, by the given node.
pure subroutine list_pop_node(this_list, this_node)
    implicit none
    !Entrada:
    class(Tlist), intent(inout) :: this_list
    type(Tnode), pointer :: this_node
    !Local:
    
    if (associated(this_node%prev).and.associated(this_node%next)) then
        !In List middle
        this_node%next%prev => this_node%prev
        this_node%prev%next => this_node%next

    else if (associated(this_node%prev)) then
        !In List tail
        nullify(this_node%prev%next)
        this_list%tail => this_node%prev

    else if (associated(this_node%next)) then
        !In List head
        nullify(this_node%next%prev)
        this_list%head => this_node%next
    end if

    !Remove node from count
    this_list%num_items = this_list%num_items - 1

    !Destroy node content
    call this_node%destroy()            
    !Free it's memmory
    deallocate(this_node)
    
    
end subroutine list_pop_node
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Delete a node from the list and frees the memory in the item.
pure subroutine node_finalizer( this_node )
    implicit none
    !Entrada:
    class(Tnode), intent(inout) :: this_node
    !Local:
    
    !Deallocate it's item
    if (allocated(this_node%item)) deallocate(this_node%item)
    !Nullify it's pointers
    nullify(this_node%next)
    nullify(this_node%prev)
        
end subroutine node_finalizer
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Delete the entire list. Nullifing the head and triggering one node_finalizer
! after the other
pure subroutine list_finalizer( this_list )
    implicit none
    !Entrada:
    class(Tlist), intent(inout) :: this_list
    !Local:

    this_list%num_items = 0
    if (associated(this_list%head)) then
        call this_list%head%destroy_all()
        deallocate(this_list%head)
        nullify(this_list%head)
        nullify(this_list%tail)
    end if
    
end subroutine list_finalizer
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Delete the entire list. Nullifing the head and triggering one node_finalizer
! after the other
pure subroutine destroy_list( this_list )
    implicit none
    !Entrada:
    type(Tlist), intent(inout) :: this_list
    !Local:

    this_list%num_items = 0
    if (associated(this_list%head)) then
        call this_list%head%destroy_all()
        deallocate(this_list%head)
        nullify(this_list%head)
        nullify(this_list%tail)
    end if
    
end subroutine destroy_list
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Insert an item by creating a new node on the list.
pure subroutine list_append_item( this_list, item )
    implicit none
    !Entrada:
    class(Tlist), intent(inout) :: this_list
    type(Toutput), intent(in) :: item
    !Local:

    if (associated(this_list%tail)) then
        allocate(this_list%tail%next, source=start_node(item))
        this_list%tail%next%prev => this_list%tail
        this_list%tail => this_list%tail%next
    else
        allocate(this_list%head, source=start_node(item))
        this_list%tail => this_list%head
    end if

    this_list%num_items = this_list%num_items + 1

end subroutine list_append_item
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Loop through nodes executing the subroutine on items
subroutine list_foreach(this_list, subr)
    implicit none
    !Entrada:
    class(Tlist), intent(inout) :: this_list
    procedure(process):: subr
    !Local:
    type(Tnode), pointer:: curr
    
    !Foward sweep
    curr => this_list%head
    do while ( associated(curr) )
        if (allocated(curr%item)) call subr(curr%item)
        curr => curr%next
    end do
    
end subroutine list_foreach
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Get the number of items in list
pure function list_lenght(this_list) result(len)
    implicit none
    !Saída:
    integer:: len
    !Entrada:
    class(Tlist), intent(in) :: this_list
    !Local:
    
    len = this_list%num_items
    
end function list_lenght
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Get the number of printable variables in the output
function list_lenght_printable(this_list) result(len)
    implicit none
    !Saída:
    integer:: len
    !Entrada:
    class(Tlist), intent(in) :: this_list
    !Local:
    type(Tnode), pointer:: curr
    
    !Initialize
    len=0
    !Foward sweep
    curr => this_list%head
    do while ( associated(curr) )
        if (allocated(curr%item)) then
            
            if (curr%item%print) len = len+1

        end if
        curr => curr%next
    end do
    
    
end function list_lenght_printable
! -----------------------------------------------------------------------------

end module output_list