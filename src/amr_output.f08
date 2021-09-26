! #############################################################################
! This module contains the functions to print amr variables into files for
! post-processing tools like visit, paraview, tecplot, etc...
!
! #############################################################################
module amr_output
    use data_mesh
    use data_input
    use data_boxes
    use output_list
    use data_variable
    use amr_variables
    use amr_startvariables
    use amrex_amr_module

    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: selectAmrOutputs, writeAmrOutput
    !Variables from other modules made accesible from here
    public:: var_list
    !-------------------------------------
    
contains

! -----------------------------------------------------------------------------
! Select which variable to print in output
subroutine selectAmrOutputs()
    implicit none
    !Local:
    type(Tnode), pointer:: curr
    character(len=32):: outname
    
    ! Loop output list forward
    curr => var_list%head
    do while ( associated(curr) )
        ! Get variable name
        outname = trim(amrex_string_c_to_f(curr%item%name%data))
        ! Check if name is in input list
        if (any(input%io%out_names == outname)) then
            ! Print it
            call curr%item%setPrint(.true.)
        else
            ! Don't print it
            call curr%item%setPrint(.false.)
        end if
        ! Go to next item
        curr => curr%next
    end do
    
end subroutine selectAmrOutputs
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Write AMR variables file to disk
subroutine writeAmrOutput(t,dt,ct)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    integer, intent(in):: ct
    !Local:
    integer:: nlev, nprint, ngc
    integer:: i, lev
    integer:: sn(3), nfile
    character(len=80):: fname='', tname=''
    logical:: exist_timef
    type(Tnode), pointer:: curr
    type(Tbox):: box
    type(TamrVar):: mfabAll
    type(amrex_string), allocatable:: varnames(:)
    type(amrex_multifab), pointer:: var1(:) => null()

    ! Get levels
    nlev = amrex_get_numlevels()
    ! Only one ghost in print variable
    ngc = 1
    ! Print variable is centered
    sn = center
    
    ! Get list size
    nprint = var_list%len_print()

    if (nprint>0) then
        ! Initialize print variable
        call mfabAll%init(nprint,ngc,sn)

        ! Allocate all names
        allocate(varnames(nprint))

        ! Loop output list forward
        curr => var_list%head
        do while ( associated(curr) )
            ! Check print
            if (curr%item%print) then
                ! First printable var
                var1 => curr%item%var%mfab
                ! Exit the list loop
                exit
            end if
            ! Go to next item
            curr => curr%next
        end do
        
        ! Print variable declaration
        do lev=1, nlev
            call mfabAll%levelInit(lev,var1(lev)%ba,var1(lev)%dm)
            call AMRmesh%initCells(lev,t,dt,mfabAll,zeros)
        end do

        ! Put print variable first in the box
        call box%addVar(mfabAll)

        ! TODO: Use mfab(lev)%copy() to fill mfabAll instead of box+forcells
        ! ---------------------------------
        ! Loop output list forward
        curr => var_list%head
        i = 0
        do while ( associated(curr) )
            ! Check print
            if (curr%item%print) then
                i=i+1
                ! Add variable to box
                call box%addVar(curr%item%var)
                ! Get name information
                varnames(i) = curr%item%name
            end if
            ! Go to next item
            curr => curr%next
        end do
        ! ---------------------------------

        ! Copy infomation 
        do lev=1, nlev
            ! Copy values from all variables to the print one
            call AMRmesh%forcells(lev,t,dt,ngc,fillPrintVar,box)
        end do

        ! File name
        write(fname,'(a,i9.9)')input%io%out_prefix,ct
        ! Write to plotfile
        call amrex_write_plotfile (trim(fname), nlev, mfabAll%mfab, &
            varnames, amrex_geom, t, [ct], amrex_ref_ratio)

        ! Write file with time plot data
        if(amrex_parallel_ioprocessor()) then
            ! Set definitions
            nfile = 4654
            tname = 'time_'//input%io%out_prefix//'.visit'
            ! Check existance
            inquire(file=tname, exist=exist_timef)
            ! Open File (new when ct=1)
            if (ct==1) then
                open(unit=nfile, file=tname, status='replace', access="sequential",&
                    position='append', action='write')
            else
                open(unit=nfile, file=tname, status='unknown', access="sequential",&
                    position='append', action='write')
            end if
            ! Write filename in it
            write(nfile,'(a)')trim(fname)//'/Header'
            ! Close file
            close(nfile)
        end if

        ! Clear name info
        do i=1,nprint
            deallocate(varnames(i)%data)
        end do
        deallocate(varnames)

    else
        write(*,'(a)')'##########'
        write(*,'(a)')' WARNING: No variable set for output.'
        write(*,'(a)')'##########'
    end if
    
end subroutine writeAmrOutput
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! This function is executed inside the mesh to copy data
! from variables 2:n to the first box variable. This copy
! is needed to print data from multiple variables.
subroutine fillPrintVar(t,dt,patch,ngcl,Pvar)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    integer, intent(in):: ngcl
    type(Tpatch), intent(in):: patch
    type(Tpatchbox), target:: Pvar(:)
    !Local:
    integer:: n
    double precision, contiguous, pointer:: var(:,:,:)
    double precision, contiguous, pointer:: print_var(:,:,:)
    integer sn(3), ix(2), iy(2), iz(2)
    integer:: comp

    ! Default component to print is the first
    comp = 1
    
    ! Loop in variables to print
    ! The first is the print variable
    do n=2,size(Pvar)
        ! Get print varirable
        print_var => Pvar(1)%patch3D(n-1)
        ! Get variable
        var => Pvar(n)%patch3D(comp)
        ! Get face info
        sn = Pvar(n)%sn
        ! Get patch limits
        ix = [patch%ilo(1)-ngcl, patch%ihi(1)+ngcl]
        iy = [patch%ilo(2)-ngcl, patch%ihi(2)+ngcl]
        iz = [patch%ilo(3)-ngcl, patch%ihi(3)+ngcl]
        
        ! Copy information, making average to center
        print_var(ix(1):ix(2), iy(1):iy(2), iz(1):iz(2)) = &
            ( &
            var(ix(1):ix(2), iy(1):iy(2), iz(1):iz(2)) + &
            var(ix(1)+sn(1):ix(2)+sn(1), iy(1)+sn(2):iy(2)+sn(2), iz(1)+sn(3):iz(2)+sn(3)) &
            ) * 0.5d0

    end do

    
end subroutine fillPrintVar
! -----------------------------------------------------------------------------

end module amr_output