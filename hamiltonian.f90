module hamiltonian

    Use global

    implicit none

    double precision, allocatable :: field(:,:), coupling(:,:,:,:)

contains

! Initialize the hamiltonian and read in the parameters from files.
    subroutine init_ham()

        implicit none
        integer*4 :: inn, inc, jnn, jnc

        field = 0.0
        coupling = 0.0

        open(10, file='./params/IsingHamiltonian_field.txt')
        do inn = 1, nSites
            do inc = 1, nStates
                read(10, *) field(inn,inc)
            enddo
        enddo
        close(10)

        open(10, file='./params/IsingHamiltonian_coupling.txt')
        do inn = 1, nSites
            do jnn = inn+1, nSites
                do inc = 1, nStates
                    do jnc = 1, nStates
                        read(10, *) coupling(inn,jnn,inc,jnc)
                        coupling(jnn,inn,jnc,inc) = coupling(inn,jnn,inc,jnc)
                    enddo
                enddo
            enddo
        enddo
        close(10)

    end subroutine


! Free the hamiltonian 
    subroutine free_ham()

        implicit none

        deallocate(field)
        deallocate(coupling)

    end subroutine


! Update params to field+coupling
    subroutine update_params_forward(params, np)

        implicit none
        integer          :: np
        double precision :: params(np)
        integer          :: inn, inc, jnn, jnc, ip

        ip = 0
        do inn = 1, nSites
            do inc = 1, nStates
                ip = ip + 1
                field(inn,inc) = params(ip)
            enddo
        enddo

        do inn = 1, nSites
            do jnn = inn+1, nSites
                do inc = 1, nStates
                    do jnc = 1, nStates
                        ip = ip + 1
                        coupling(inn,jnn,inc,jnc) = params(ip)
                        coupling(jnn,inn,jnc,inc) = coupling(inn,jnn,inc,jnc)
                    enddo
                enddo
            enddo
        enddo
        
    end subroutine


! Output parameters
    subroutine output_params(params, np)

        implicit none
        integer          :: inn, inc, jnn, jnc, ine, np
        double precision :: params(np)

        open(unit=10, file='./params/params.txt')
        do ine = 1, np
            write(10, '(e15.8)') params(ine)
        enddo
        close(10)

        open(unit=10, file='./params/IsingHamiltonian_field.txt')
        do inn = 1, nSites
            do inc = 1, nStates
                write(10, '(e15.8)') field(inn,inc)  
            enddo
        enddo
        close(10)

        open(unit=10, file='./params/IsingHamiltonian_coupling.txt')
        do inn = 1, nSites
            do jnn = inn+1, nSites
                do inc = 1, nStates
                    do jnc = 1, nStates
                        write(10, '(e15.8)') coupling(inn,jnn,inc,jnc)
                    enddo
                enddo
            enddo
        enddo
        close(10)

    end subroutine


! Update field+coupling to params
    subroutine update_params_backward(params, np)

        implicit none
        integer          :: np
        double precision :: params(np)
        integer          :: inn, inc, jnn, jnc, ip

        ip = 0
        do inn = 1, nSites
            do inc = 1, nStates
                ip = ip +1
                params(ip) = field(inn,inc)  
            enddo
        enddo

        do inn = 1, nSites
            do jnn = inn+1, nSites
                do inc = 1, nStates
                    do jnc = 1, nStates
                        ip = ip + 1
                        params(ip) = coupling(inn,jnn,inc,jnc)
                    enddo
                enddo
            enddo
        enddo

    end subroutine


! Evaluate the energy change after a spin flip
    subroutine eval_energy(IsingVar, flipN, flipC, IsingEner)

        implicit none
        integer*4, intent(in) :: IsingVar(nSites)
        double precision      :: IsingEner
        integer*4             :: flipN, flipC

        integer*4             :: inn

        IsingEner = 0.0

        IsingEner = IsingEner + field(flipN, flipC)

        do inn = 1, nSites
            if (inn<flipN) then
                IsingEner = IsingEner + coupling(inn,flipN,IsingVar(inn),flipC)
            elseif (inn>flipN) then
                IsingEner = IsingEner + coupling(flipN,inn,flipC,IsingVar(inn))
            endif
        enddo
    end subroutine


! Evaluate the total energy given a configuration
    subroutine eval_tot_energy(IsingVar, IsingEner)

        implicit none
        integer*4, intent(in) :: IsingVar(nSites)
        double precision      :: IsingEner

        integer*4             :: inn, jnn

        IsingEner = 0.0

        do inn = 1, nSites
            IsingEner = IsingEner + field(inn,IsingVar(inn))
        enddo

        do inn = 1, nSites
            do jnn = inn+1, nSites
                IsingEner = IsingEner + coupling(inn,jnn,IsingVar(inn),IsingVar(jnn))
            enddo
        enddo

    end subroutine

end module
