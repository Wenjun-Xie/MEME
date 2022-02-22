module ensemble_average

    Use global

    implicit none

    double precision, save :: nDump = 0

    integer*4, allocatable :: pos_field(:,:), pos_coupling(:,:,:,:)

contains

! Used to get index in params
    subroutine map_field()
        implicit none
        integer*4        :: inn, inc, ip
        ip = 0
        do inn = 1, nSites
            do inc = 1, nStates
                ip = ip +1
                pos_field(inn,inc) = ip
            enddo
        enddo
    end


! Used to get index in params
    subroutine map_coupling()
        implicit none
        integer*4        :: inn, inc, jnn, jnc, ip
        ip = nSites*nStates
        do inn = 1, nSites
            do jnn = inn+1, nSites
                do inc = 1, nStates
                    do jnc = 1, nStates
                        ip = ip +1
                        pos_coupling(inn,jnn,inc,jnc) = ip
                    enddo
                enddo
            enddo
        enddo
    end


! Calculate the ensemble average
    subroutine calc_ensemble_avg(IsingVar, ensembleAvg)

        implicit none

        integer*4        :: IsingVar(nSites)
        double precision :: ensembleAvg(nExpConstr)

        integer*4        :: inn, jnn, ip

        nDump = nDump + 1    ! record how many averages that we performed

        do inn = 1, nSites
            ip = pos_field(inn, IsingVar(inn))
            ensembleAvg(ip) = ensembleAvg(ip) + 1
        enddo

        do inn = 1, nSites
            do jnn = inn+1, nSites
                ip = pos_coupling(inn, jnn, IsingVar(inn), IsingVar(jnn))
                ensembleAvg(ip) = ensembleAvg(ip) + 1
            enddo
        enddo

    end subroutine


! Output ensemble average
    subroutine output_ensemble_avg(ensembleAvg, outfile)
        implicit none

        character*500    :: outfile
        double precision :: ensembleAvg(nExpConstr)

        integer*4        :: ine

        open(unit=10, file=trim(outfile))
        do ine = 1, nExpConstr
            write(10, '(e15.8)') ensembleAvg(ine)
        enddo
        close(10)

    end subroutine

end module
