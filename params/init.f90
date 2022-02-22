program init
    implicit none
    integer*4 :: inn, inc, jnn, jnc, nSites, nStates
    real:: p

    nSites = 287
    nStates = 21

    open(1, file='../experimental_constraints.txt')

    open(10, file='IsingHamiltonian_field.txt')
    do inn = 1, nSites
        do inc = 1, nStates
            read(1, *) p
            write(10, '(f8.5)') -log(p+0.0001)
        enddo
    enddo
    close(10)

    open(10, file='IsingHamiltonian_coupling.txt')
    do inn = 1, nSites
        do jnn = inn+1, nSites
            do inc = 1, nStates
                do jnc = 1, nStates
                    write(10, '(f8.6)') 0.0
                enddo
            enddo
        enddo
    enddo
    close(10)
end
