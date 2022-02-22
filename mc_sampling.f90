module mc_sampling

    Use global
    Use hamiltonian

    implicit none

contains

! main monte carlo routine
!
! MC move accepted based on the Metropolis acceptance criterion

    subroutine mc_move(IsingVar, etot, beta_rpl)

        implicit none

        integer*4 :: IsingVar(nSites)

        double precision :: rnum, oldE, newE, deltaE, beta_rpl, etot
        integer*4        :: inn_iter, inn, inc, accept

        do inn_iter = 1, nSites

            call random_number(rnum)
            inn = floor(rnum*nSites) + 1

            call eval_energy(IsingVar, inn, IsingVar(inn), oldE)

            call random_number(rnum)
            inc = floor(rnum*nStates) + 1

            call eval_energy(IsingVar, inn, inc, newE)

            deltaE = newE - oldE

            accept = 0
            if (deltaE < 0) then
                accept = 1
            else
                call random_number(rnum)
                if (rnum < exp(-beta_rpl*deltaE)) then
                    accept = 1
                endif
            endif

            if (accept .eq. 1) then
                IsingVar(inn) = inc
                etot = etot + deltaE
            endif

        enddo

    end subroutine

end module
