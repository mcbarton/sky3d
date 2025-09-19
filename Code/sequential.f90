!------------------------------------------------------------------------------
! MODULE: Parallel
!------------------------------------------------------------------------------
! DESCRIPTION:
!> @brief
!!contains the same routines as in file \c parallel.f90, but with mostly empty
!!functions to enable sequential running.
!------------------------------------------------------------------------------
MODULE Parallel
  USE Params, ONLY: wflag,db
  USE Levels, ONLY: nstmax,npsi,nstloc
  IMPLICIT NONE
  SAVE
  LOGICAL,PARAMETER :: tmpi=.FALSE.
  INTEGER, ALLOCATABLE :: node(:),localindex(:),globalindex(:)
  INTEGER :: mpi_nprocs,mpi_ierror,mpi_myproc, &
       processor_name,proc_namelen
  INTEGER :: mpi_comm_world,mpi_sum,mpi_double_precision
CONTAINS     !  all dummy subroutines to run on a sequential machine
  !************************************************************************
  SUBROUTINE alloc_nodes
    ALLOCATE(node(nstmax),localindex(nstmax),globalindex(nstmax))
  END SUBROUTINE alloc_nodes
  !************************************************************************
  SUBROUTINE init_all_mpi
    WRITE(*,*) '***** Running sequential version *****'
    mpi_nprocs=1
    mpi_myproc=0
    wflag=.TRUE.
  END SUBROUTINE init_all_mpi
  !************************************************************************
  !> dummy function for the MPI routine
  SUBROUTINE mpi_init(ierror)
    INTEGER, intent(in) :: ierror
    STOP ' MPI_INIT: parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_init
  !************************************************************************
  !> dummy function for the MPI routine
  SUBROUTINE mpi_comm_size(comm_world,nprocs,ierror)
    INTEGER, intent(in) :: ierror, nprocs, comm_world
    STOP ' MPI_COMM_SIZE: parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_comm_size
  !************************************************************************
  !> dummy function for the MPI routine
  SUBROUTINE mpi_comm_rank(comm_world,myproc,ierror)
    INTEGER, intent(in) :: ierror, myproc, comm_world
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_comm_rank
  !************************************************************************
  !> dummy function for the MPI routine
  SUBROUTINE mpi_get_processor_name(processor_name,proc_namelen,ierror)
    INTEGER, intent(in) :: ierror, processor_name, proc_namelen
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_get_processor_name
  !************************************************************************
  !> dummy function for the MPI routine
  SUBROUTINE mpi_barrier (comm_world, ierror)
    INTEGER , intent(in) :: ierror, comm_world
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_barrier
  !************************************************************************
  SUBROUTINE associate_nodes
    INTEGER :: i
    node=0
    nstloc=nstmax
    FORALL(i=1:nstmax)
       globalindex(i)=i
       localindex(i)=i
    END FORALL
  END SUBROUTINE associate_nodes
  !************************************************************************
  !> dummy function for the MPI routine
  SUBROUTINE mpi_allreduce(rho,tmp_rho,length,        &
       i_double_precision,sum,  &
       comm_world,ierror)
    INTEGER, intent(in) :: ierror, comm_world, i_double_precision, length, sum
    REAL(db), DIMENSION(:), INTENT(IN) :: rho,tmp_rho
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE mpi_allreduce
  !************************************************************************
  SUBROUTINE collect_densities
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_densities
  !************************************************************************
  SUBROUTINE collect_sp_properties
    STOP ' parallel calls inhibited '
    RETURN
  END SUBROUTINE collect_sp_properties
  !************************************************************************
  SUBROUTINE finish_mpi
  END SUBROUTINE finish_mpi
END MODULE Parallel
