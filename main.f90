Program main

    Use global
    Use mc_sampling
    Use hamiltonian
    Use ensemble_average

    implicit none
    include 'mpif.h'
    
    character*500:: infile, outfile
    integer   :: rseed
    integer*4 :: nMCIter, recordFreq
    integer*4 :: inn, inc, imc, ine, iter

    integer*4, allocatable :: IsingVar(:, :), IsingVarTmp(:)
    double precision, allocatable :: ensembleAvg(:), ensembleAvgTot(:), &
                                     expConstr(:), params(:), delta(:), &
                                     etot(:), dx_save(:)


    double precision :: rnum, diff, learning_rate
    
    double precision :: temp_min, temp_max, dT, deltaU, deltaBeta, etmp, ratio
    integer          :: rpl_targetT, accept, naccept, ntrial, irpl, irpl_start

    ! MPI stuff
    integer          :: ierr, proc_id, no_procs
    integer          :: buf_int(10), buflen
    double precision :: buf_dble(10)

    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, no_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, proc_id, ierr)
    buflen = 10

    if (proc_id .eq. 0) then 
        ! read input file
        infile = 'input.txt'
        open(unit=10, file=trim(infile))
        read(10, *) nStates
        read(10, *) nSites
        read(10, *) nreplica

        read(10, *) temp_min
        read(10, *) temp_max

        read(10, *) nMCIter
        read(10, *) recordFreq
        close(10)

        ! ensemble average
        nExpConstr = nSites*nStates + nSites*(nSites-1)/2*nStates*nStates

        allocate(expConstr(nExpConstr))
        open(unit=10,file='experimental_constraints.txt')
        do ine = 1, nExpConstr
            read(10, *) expConstr(ine)
        enddo
        close(10)

        buf_dble(1) = temp_min
        buf_dble(2) = temp_max

        buf_int(1)  = nStates
        buf_int(2)  = nSites
        buf_int(3)  = nReplica
        buf_int(4)  = nMCIter
        buf_int(5)  = recordFreq
        buf_int(6)  = nExpConstr
        open(unit=10, file='log.txt', status='replace')
        close(10)
    endif
    !
    !     Broadcast the input data
    !
    Call MPI_BCAST(buf_dble,buflen,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    Call MPI_BCAST(buf_int, buflen,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)

    temp_min     = buf_dble(1)
    temp_max     = buf_dble(2)

    nStates  = buf_int(1)
    nSites = buf_int(2)
    nReplica     = buf_int(3)
    nMCIter      = buf_int(4)
    recordFreq   = buf_int(5)
    nExpConstr   = buf_int(6)

    ! output the system setup for debugging purposes
    if (proc_id .eq. 0) then
        open(unit=10, file='log.txt', position='append')
        write(10,*) "++++++++++ Setup ++++++++++"
        write(10,*) "nStates       : ", nStates 
        write(10,*) "nSites      : ", nSites 
        write(10,*) "temp              : ", temp_min, temp_max
        write(10,*) "nMCIter           : ", nMCIter
        write(10,*) "recordFreq        : ", recordFreq
        write(10,*) "no_procs          : ", no_procs

        close(10)
    endif

    ! allocate variables
    allocate(IsingVar(nSites, nreplica))
    allocate(IsingVarTmp(nSites))

    allocate(etot(nreplica))
    allocate(temp(nReplica))
    allocate(beta(nReplica))
    dT = (temp_max-temp_min)/(nreplica-1)
    do irpl = 1, nReplica
        temp(irpl) = temp_min + (irpl-1)*dT
        beta(irpl) = 1.0 / temp(irpl) / kB
        if(abs(temp(irpl)-1.0) < 1e-16) then
            rpl_targetT = irpl
        endif
    enddo

    if (proc_id .eq. 0)  then
        open(unit=10, file='log.txt', position='append')
        write(10, *) "target: ", rpl_targetT 
        write(10, *) temp
        close(10)
    endif

    ! allocate variables
    allocate(params(nExpConstr))
    allocate(field(nSites, nStates))
    allocate(coupling(nSites, nSites, nStates, nStates))
    allocate(pos_field(nSites, nStates))
    allocate(pos_coupling(nSites, nSites, nStates, nStates))

    if (proc_id .eq. 0) then
        call init_ham()
        call update_params_backward(params, nExpConstr)
    endif
    Call MPI_BCAST(params,nExpConstr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    Call update_params_forward(params, nExpConstr)

    call map_field()
    call map_coupling()


    ! initialize the random number
    rseed = 12345 + proc_id
    call init_random_seed(rseed) 

    ! initialize the Ising variable
    do irpl = 1, nReplica
        do inn = 1, nSites
            call random_number(rnum)
            inc = floor(rnum*nStates) + 1
            IsingVar(inn, irpl) = inc
        enddo
    enddo

    ! initialize the average
    allocate(ensembleAvg(nExpConstr))
    allocate(ensembleAvgTot(nExpConstr))
    allocate(delta(nExpConstr))
    allocate(dx_save(nExpConstr))
    dx_save=0.0

    diff = 99999  ! set to a large initial number
    learning_rate = 0.01  ! SGD optimization
    iter = 0
    do while (diff > 0.05)
        iter = iter + 1
        ensembleAvg = 0.0
        ensembleAvgTot = 0.0

        do irpl = 1, nReplica
            call eval_tot_energy(IsingVar(:,irpl), etot(irpl))
        enddo

        naccept = 0
        ntrial  = 0
        do imc = 1, nMCIter

            do irpl = 1, nreplica
                call mc_move(IsingVar(:,irpl), etot(irpl), beta(irpl))
            enddo

            if (mod(imc, recordFreq) .eq. 0) then

                ! save the configurations first
                call calc_ensemble_avg(IsingVar(:,rpl_targetT), ensembleAvg)

                call random_number(rnum)
                if (rnum <= 0.50) then
                    irpl_start = 1
                else
                    irpl_start = 2
                endif

                do irpl = irpl_start, nreplica, 2
                    if ((irpl+1) <= nreplica) then
                        deltaBeta = beta(irpl+1) - beta(irpl)
                        deltaU = etot(irpl+1) - etot(irpl)

                        accept = 0
                        if (deltaBeta*deltaU > 0) then
                            accept = 1
                        else
                            call random_number(rnum)
                            if (rnum < exp(deltaBeta*deltaU)) then
                                accept = 1
                            endif
                        endif
                        if (accept .eq. 1) then
                            IsingVarTmp = IsingVar(:,irpl)
                            IsingVar(:,irpl) = IsingVar(:,irpl+1)
                            IsingVar(:,irpl+1) = IsingVarTmp

                            etmp = etot(irpl)
                            etot(irpl) = etot(irpl+1)
                            etot(irpl+1) = etmp
                        endif

                        naccept = naccept + accept
                        ntrial  = ntrial + 1

                    endif
                enddo

            endif
        enddo

        Call MPI_Reduce(ensembleAvg, ensembleAvgTot, nExpConstr, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (proc_id .eq. 0) then
            ensembleAvg = ensembleAvgTot / no_procs / (nMCIter/recordFreq)

            diff = sum(abs(ensembleAvg-expConstr)) / sum(expConstr)

            !add L2-regularizer with alpha=0.01
            delta = (expConstr-ensembleAvg) + 0.01*2*params

            !add momentum step with beta=0.9
            params = params - learning_rate * delta + 0.9*dx_save
            dx_save = - learning_rate * delta + 0.9*dx_save

        endif

        Call MPI_BCAST(params,nExpConstr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        Call update_params_forward(params, nExpConstr)

        Call MPI_BCAST(diff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        if (proc_id .eq. 0) then

            open(unit=10, file='log.txt', position='append')
            write(10, '(i8, 1x, f12.6)') iter, diff
            close(10)


            if (mod(iter, 10) .eq. 0) then
                outfile = './output/ensemble_average.txt'
                call output_ensemble_avg(ensembleAvg, outfile)
                call output_params(params, nExpConstr)
            endif
        endif
    enddo

    ! conclude the code and free memory
    call free_ham()

end program

SUBROUTINE init_random_seed(rseed)
    implicit none

    INTEGER :: i, n, clock, rseed
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /) + rseed
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
END SUBROUTINE
