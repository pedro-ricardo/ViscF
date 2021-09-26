program main
    use mpi_f08
    use data_input
    use data_mesh
    use data_boxes
    use data_variable
    use input_reader
    use amr_remesh
    use amr_variables
    use amr_output
    use amr_ghosts
    use amrex_amr_module
    use example_functions
    use example_solver
    
    implicit none
    integer:: ierr, proc_id
    double precision:: t1, t2, t11=0.d0, t21=0.d0
    double precision:: t, dt
    double precision:: residue
    integer:: ct, l
    type(TamrVar):: rhs
    
    ! Initialize MPI
    call mpi_init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, proc_id, ierr)
    t1 = mpi_Wtime()

    ! Initialize AMReX
    call amrex_init(comm=MPI_COMM_WORLD%mpi_val,proc_parmparse=readInputFile)
    call amrex_amrcore_init()
    ! Initialize global variables
    call initAmrMesh()
    call initAmrVariables()
    call initAmrOutputs()
    ! Set printable variables
    call selectAmrOutputs()
    ! Initialize mesh
    call amrex_init_from_scratch(0.d0)

    ! Mesh info
    if(proc_id==0) then
        write(*,'(a)')'Domain Information'
        write(*,'(a,i12)')'Total Cells: ',AMRmesh%ncell
        write(*,'(a,i3)')'Total Levels: ',AMRmesh%nlev
        write(*,'(a,8(i3,1x))')'Refinement Ratios: ',amrex_ref_ratio
        do l=1,AMRmesh%nlev
            write(*,'(a,i3)')'Level : ',l
            write(*,'(a,3(e12.5,1x))')'Deltas: ',AMRmesh%geo(l)%dx
            write(*,'(a,3(i12,1x))')'Cells : ',AMRmesh%geo(l)%domain%hi-AMRmesh%geo(l)%domain%lo+1
        end do
    end if
    
    ! Loop definitions
    dt = 1.d-2 !s
    t = 0.d0
    ct = 0
    residue = 0.d0

    ! ------------------------------------------
    ! Time loop
    do while (t < input%sim%end_time)
        ! Advance counters
        ct = ct +1
        t = t + dt
        ! Print timestep header
        if(proc_id==0) then
            if (ct==1 .or. mod(ct,10)==0) write(*,'(a10,4(a1,a10))')'Iter','|','Time','|','Delta','|','Residue','|','It. Time'
            write(*,'(i10,4(a1,es10.3))')ct,'|',t,'|',dt,'|',residue,'|',t21-t11
        end if

        t11 = mpi_Wtime()
        ! RHS calculation
        call rhs%init(comp=1,ngc=1,sn=center)
        ! Build it
        do l = 1,AMRmesh%nlev
            call rhs%levelInit(l,AMRmesh%ba(l),AMRmesh%dm(l))
            ! Calculate on level
            call AMRmesh%initCells(l,t,dt,rhs,exact_rhs)
        end do
        ! Ghost fill
        do l = 1,AMRmesh%nlev
            call fillAllBoundaries(l,rhs,rhs)
        end do

        ! Calculate
        call solve_elliptic(t,dt,phi,rhs,residue)
        
        ! Deallocate rhs
        call rhs%destroy()
        
        ! Check plot step
        if ( ct==1 .or. (mod(ct,input%io%out_it)==0 .and. input%io%out_it>0) ) then
            ! Plot it
            call writeAmrOutput(t,dt,ct)
        end if

        t21 = mpi_Wtime() 
        ! Max iteration condition
        if (ct >= input%sim%max_it .and. input%sim%max_it > 0) exit
    end do
    ! ------------------------------------------
    
    if(proc_id==0) write(*,'(a)')'Simulation Ended'

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

    if(proc_id==0) write(*,'(a,f8.4,a)')'Time enlapsed: ',t2-t1,' seconds'
    ! Finalize MPI
    call mpi_finalize(ierr)

end program main

