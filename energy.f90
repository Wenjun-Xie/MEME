Program energy

    Use global
    Use mc_sampling
    Use hamiltonian

    implicit none

    character*500:: infile
    integer*4 :: inn, inc

    integer*4, allocatable :: IsingVar(:)

    double precision :: totE

    integer :: i, mut_num

    ! read input file
    infile = 'input.txt'
    open(unit=10, file=trim(infile))
    read(10, *) nStates
    read(10, *) nSites
    close(10)

    ! allocate variables
    allocate(IsingVar(nSites))
    allocate(field(nSites, nStates))
    allocate(coupling(nSites, nSites, nStates, nStates))

    call init_ham()

    mut_num=20000
    open(unit=1, file='msa_mut_non_gap.txt')
    open(unit=11, file='msa_mut_MaxEnt_energy.txt', status='replace')

    do i=1,mut_num
        do inn=1,nSites
            if(inn==nSites)then
                read(1,'(i2)') IsingVar (inn)
            else
                read(1,'(i3)',advance='no') IsingVar (inn)
            endif
        enddo

        call eval_tot_energy(IsingVar, totE)
        write(11,*) totE

    enddo

    call free_ham()

end program

