! #############################################################################
! This module contain functions to read and save information from input files
!
! #############################################################################
module input_reader
    use data_input
    use data_mesh
    use amrex_amr_module
    use json_module
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: readInputFile
    !-------------------------------------
    
contains


! -----------------------------------------------------------------------------
! Read input file and fill input structure.
! This routine is executed by AMReX during initialization in `amrex_init` function.
subroutine readInputFile() bind(c)
    implicit none
    !Local:
    type(json_file):: inputfile
    character(kind=c_char,len=:),allocatable:: svec
    character(32):: sint
    double precision, allocatable:: dvec(:)
    integer, allocatable:: ivec(:)
    logical, allocatable:: lvec(:)
    double precision:: dval
    integer:: ival, ival2, n_items
    integer:: i
    logical:: found

    ! Initialize Json object
    call inputfile%initialize(comment_char=json_CK_'//', path_separator=json_CK_'.')
    ! Set file
    call inputfile%load_file(filename = 'main.jsonc')
    ! File check
    if (inputfile%failed() .and. amrex_parallel_ioprocessor()) then
        write(*,*)'No input file found. Using default parameters.'
    end if

    ! Initialize MFSim input
    call input%init()

    ! Mesh
    call inputfile%get('Mesh.Domain Mininum',dvec,found)
    if (found) input%msh%dom_lo=dvec
    if(allocated(dvec))deallocate(dvec)
    call inputfile%get('Mesh.Domain Maximum',dvec,found)
    if (found) input%msh%dom_hi=dvec
    if(allocated(dvec))deallocate(dvec)
    call inputfile%get('Mesh.Lbot Cells',ivec,found)
    if (found) input%msh%bot_cells=ivec
    if(allocated(ivec))deallocate(ivec)
    call inputfile%get('Mesh.Domain Periodicity',lvec,found)
    if (found) input%msh%periodic=lvec
    if(allocated(lvec))deallocate(lvec)
    call inputfile%get('Mesh.Levels above Lbot',ival,found)
    if (found) input%msh%levels=ival
    call inputfile%get('Mesh.Refinement Ratio',ival,found)
    ival2 = input%msh%ratio(1)
    if(allocated(input%msh%ratio)) deallocate(input%msh%ratio)
    allocate(input%msh%ratio(max(input%msh%levels,1)))
    input%msh%ratio = ival2
    if (found) input%msh%ratio=ival
    call inputfile%info('Mesh.Fine Mesh Initial Topology',found,n_children=n_items)
    if (found) then
        if(allocated(input%msh%init_grid))deallocate(input%msh%init_grid)
        allocate(input%msh%init_grid(n_items))
        do i=1, n_items
            write(sint,*)i
            call inputfile%get('Mesh.Fine Mesh Initial Topology('//trim(sint)//').Box Min',dvec,found)
            if (found) input%msh%init_grid(i)%lo = dvec
            if(allocated(dvec))deallocate(dvec)
            call inputfile%get('Mesh.Fine Mesh Initial Topology('//trim(sint)//').Box Max',dvec,found)
            if (found) input%msh%init_grid(i)%hi = dvec
            if(allocated(dvec))deallocate(dvec)
        end do
    end if

    ! Remesh
    call inputfile%get('Regrid.Frequency',ival,found)
    ival2 = input%remsh%freq(1)
    if(allocated(input%remsh%freq)) deallocate(input%remsh%freq)
    allocate(input%remsh%freq(max(input%msh%levels,1)))
    input%remsh%freq = ival2
    if (found) input%remsh%freq=ival
    call inputfile%get('Regrid.Cut Off',dval,found)
    if (found) input%remsh%cut_off=dval
    call inputfile%get('Regrid.Buffer Zone',ival,found)
    if (found) input%remsh%buff_fp=ival

    ! Simulation
    call inputfile%get('Simulation Controls.Maximum Iterations',ival,found)
    if (found) input%sim%max_it=ival
    call inputfile%get('Simulation Controls.End Time',dval,found)
    if (found) input%sim%end_time=dval

    ! Input/Output
    call inputfile%get('I/O.Restart Files.Number of Files',ival,found)
    if (found) input%io%rst_max=ival
    call inputfile%get('I/O.Restart Files.File Prefix',svec,found)
    if (found) then
        if(allocated(input%io%rst_prefix))deallocate(input%io%rst_prefix)
        input%io%rst_prefix=svec
    end if
    if(allocated(svec))deallocate(svec)
    call inputfile%get('I/O.Restart Files.Save by Iteration',ival,found)
    if (found) input%io%rst_it=ival
    call inputfile%get('I/O.Restart Files.Save by Time',dval,found)
    if (found) input%io%rst_time=dval
    call inputfile%get('I/O.Output Files.File Prefix',svec,found)
    if (found) then
        if(allocated(input%io%out_prefix))deallocate(input%io%out_prefix)
        input%io%out_prefix=svec
    end if
    if(allocated(svec))deallocate(svec)
    call inputfile%get('I/O.Output Files.Save by Iteration',ival,found)
    if (found) input%io%out_it=ival
    call inputfile%get('I/O.Output Files.Save by Time',dval,found)
    if (found) input%io%out_time=dval
    call inputfile%info('I/O.Output Files.Print Variables',found,n_children=n_items)
    if (found) then
        if(allocated(input%io%out_names))deallocate(input%io%out_names)
        allocate(input%io%out_names(n_items))
        do i=1, n_items
            write(sint,*)i
            call inputfile%get('I/O.Output Files.Print Variables('//trim(sint)//')',svec,found)
            if (found) input%io%out_names(i) = svec
            if(allocated(svec))deallocate(svec)
        end do
    end if
    

    ! End Json object
    call inputfile%destroy()

    ! Insert the parameters in AMReX
    call setAmrexParameters()

    ! Copy parameters to AMRmesh variable
    call AMRmesh%init(input%msh)

end subroutine readInputFile
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Set the input parameters into AMReX
subroutine setAmrexParameters()
    implicit none
    !Local:
    type(amrex_parmparse) :: pp
    
    ! Global Paramaters
    call amrex_parmparse_build(pp)
    call pp%add("max_step", input%sim%max_it)
    call pp%add("stop_time", input%sim%end_time)
    call amrex_parmparse_destroy(pp)

    ! Geometry Parameters
    call amrex_parmparse_build(pp,"geometry")
    call pp%add("coord_sys", 0)
    call pp%addarr("prob_lo", input%msh%dom_lo)
    call pp%addarr("prob_hi", input%msh%dom_hi)
    call pp%addarr("is_periodic", transfer(input%msh%periodic,1,3))
    call amrex_parmparse_destroy(pp)

    ! AMR Parameters
    call amrex_parmparse_build(pp,"amr")
    call pp%addarr("n_cell", input%msh%bot_cells)
    call pp%add("max_level", input%msh%levels)
    call pp%addarr("ref_ratio", input%msh%ratio)
    call pp%add("max_grid_size", input%msh%max_grid_size)

    call pp%addarr("regrid_int", input%remsh%freq)
    call pp%add("grid_eff", input%remsh%cut_off)
    call pp%add("regrid_on_restart", input%remsh%on_restart)
    call pp%add("n_error_buf", input%remsh%buff_fp)

    call pp%add("checkpoint_nfiles", input%io%rst_max)
    call pp%add("check_file", input%io%rst_prefix)
    call pp%add("check_int", input%io%rst_it)
    call pp%add("check_per", input%io%rst_time)
    call pp%add("plot_file", input%io%out_prefix)
    call pp%add("plot_int", input%io%out_it)
    call pp%add("plot_per", input%io%out_time)

    call amrex_parmparse_destroy(pp)

    ! AMReX System Parameters
    call amrex_parmparse_build(pp,"amrex")
    call pp%add("verbose", 0) ! No debug prints
    call amrex_parmparse_destroy(pp)
    
end subroutine setAmrexParameters
! -----------------------------------------------------------------------------


end module input_reader