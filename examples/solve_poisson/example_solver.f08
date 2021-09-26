! #############################################################################
! Header
!
! #############################################################################
module example_solver
    use data_variable
    use data_boxes
    use data_mesh
    use amrex_base_module
    use amrex_linear_solver_module
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: press_rhs, solve_elliptic
    !-------------------------------------
    
contains

! -----------------------------------------------------------------------------
! Calcuate pressure correction equation RHS
subroutine press_rhs(t,dt,patch,ngcl,Pvar)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    integer, intent(in):: ngcl
    type(Tpatch), intent(in):: patch
    type(Tpatchbox), target:: Pvar(:)
    !Local:
    integer:: i, j, k, sn(3)
    double precision :: ue, uw, vn, vs, wt, wb, dudx, dvdy, dwdz
    double precision, contiguous, pointer:: rhs(:,:,:),&
        rho(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:)
    
    ! Face info
    sn = Pvar(1)%sn
    ! Get variable in correct indexes
    rhs => Pvar(1)%patch3D()
    u => Pvar(2)%patch3D()
    v => Pvar(3)%patch3D()
    w => Pvar(4)%patch3D()
    rho => Pvar(5)%patch3D()

    do k= patch%ilo(3), patch%ihi(3)+sn(3)
        do j= patch%ilo(2), patch%ihi(2)+sn(2)
            do i= patch%ilo(1), patch%ihi(1)+sn(1)

                ue = u(i+1,j  ,k  ); uw = u(i  ,j  ,k  )
                vn = v(i  ,j+1,k  ); vs = v(i  ,j  ,k  )
                wt = w(i  ,j  ,k+1); wb = w(i  ,j  ,k  )

                dudx = (ue-uw)/patch%dn(1)
                dvdy = (vn-vs)/patch%dn(2)
                dwdz = (wt-wb)/patch%dn(3)

                rhs(i,j,k) = rho(i,j,k)*(dudx + dvdy + dwdz)/dt
            end do
        end do
    end do

end subroutine press_rhs
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Uses MLMG to solve the ellictic equation:
! ∇²(phi) = rhs
subroutine solve_elliptic(t,dt,phi,rhs,res)
    implicit none
    !Entrada:
    double precision, intent(in):: t, dt
    type(TamrVar), intent(in):: rhs
    type(TamrVar), intent(inout):: phi
    double precision, intent(out):: res
    !Local:
    integer:: l, Clev
    type(amrex_poisson) :: poisson
    type(amrex_multigrid) :: multigrid
    
    ! Build equation type as Poisson
    call amrex_poisson_build(poisson, AMRmesh%geo, AMRmesh%ba, AMRmesh%dm, &
        metric_term=.false., agglomeration=.true., consolidation=.true., &
        max_coarsening_level=30)
    ! Define discretization Order on laplacian
    call poisson%set_maxorder(2)
    ! Set Boundary condition types solver
    call poisson%set_domain_bc([amrex_lo_periodic,amrex_lo_periodic,amrex_lo_periodic], &
        [amrex_lo_periodic,amrex_lo_periodic,amrex_lo_periodic])
    
    ! ???? Don't know yet ????
    do l = 1, AMRmesh%nlev
        Clev = l-1
        call poisson%set_level_bc(Clev, phi%mfab(l))
    end do

    ! Build solver for poisson solver type
    call amrex_multigrid_build(multigrid, poisson)
    ! Information on MLMG iterations
    call multigrid%set_verbose(0)
    ! Information on bottom solver
    call multigrid%set_bottom_verbose(0)
    ! Maximun number of iterations
    call multigrid%set_max_iter(100)
    ! F-cycle MG iteration before V-Cycle
    call multigrid%set_max_fmg_iter(0)
    ! Set Bottom solver type (0-smooth 1-bcgstab 2-cg 3-hypre 4-petsc)
    call multigrid%set_bottom_solver(1)

    ! WARNING - Fixed iterations for SpeedUp evaluation
    call multigrid%set_fixed_iter(5)

    ! Solve equation      solution  rhs    abs_tol, rel_tol 
    res = multigrid%solve(phi%mfab, rhs%mfab, 0.d0, 0.d0)

    ! Free equation and solver
    call amrex_poisson_destroy(poisson)
    call amrex_multigrid_destroy(multigrid)

end subroutine solve_elliptic
! -----------------------------------------------------------------------------

end module example_solver