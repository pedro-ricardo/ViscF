program main
    use mpi_f08
    use data_input
    use input_reader
    use data_mesh
    use amr_remesh
    use amr_variables
    use amr_output
    use amrex_amr_module
    
    implicit none
    integer :: ierr, i
    double precision:: t1, t2, t = 0.d0
    
    ! Initialize MPI
    call mpi_init(ierr)

    t1 = mpi_Wtime()

    ! Initialize AMReX
    call amrex_init(comm=MPI_COMM_WORLD%mpi_val,proc_parmparse=readInputFile)
    call amrex_amrcore_init()
    ! Initialize global variables
    call initAmrMesh()
    call initAmrVariables()
    call initAmrOutputs()
    ! Initialize mesh
    call amrex_init_from_scratch(t)

    ! Mesh info
    if(amrex_parallel_ioprocessor()) then
        write(*,'(a)')'Domain Information'
        write(*,'(a,i3)')'Total Levels: ',amrex_get_numlevels()
        write(*,'(a,8(i3,1x))')'Refinement Ratios: ',amrex_ref_ratio
        do i=1,amrex_get_numlevels()
            write(*,'(a,i3)')'Level : ',i
            write(*,'(a,3(e12.5,1x))')'Deltas: ',amrex_geom(i-1)%dx
            write(*,'(a,3(i12,1x))')'Cells : ',amrex_geom(i-1)%domain%hi-amrex_geom(i-1)%domain%lo+1
        end do
    end if
    
    ! PLot it
    call selectAmrOutputs()
    call writeAmrOutput(1.2d0,1.d-4,15)

    ! Deallocate input
    call input%destroy()
    ! Deallocate output vars list
    call var_list%destroy()
    ! Deallocate global variables
    call destroyAmrVariables()
    ! Deallocate Mesh structure
    call AMRmesh%destroy()
    ! Finalize AMReX
    call amrex_amrcore_finalize() ! TODO: REPORT ERROR IN DEALLOCATION HERE
    call amrex_finalize()

    t2 = mpi_Wtime()

    if(amrex_parallel_ioprocessor()) write(*,'(a,f8.4,a)')'Time enlapsed: ',t2-t1,' seconds'
    ! Finalize MPI
    call mpi_finalize(ierr)

end program main

