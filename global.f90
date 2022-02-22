module global
    
    implicit none

    integer          :: nStates, nSites, nReplica, nExpConstr
    double precision, parameter :: kB = 1.0
    double precision, allocatable :: beta(:), temp(:)

end module
