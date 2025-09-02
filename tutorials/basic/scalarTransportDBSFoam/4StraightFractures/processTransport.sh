#!/bin/bash

###### USERS INPUT ############################################################


#### END OF USER INPUT #######################################################


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    echo -e "processSpeciesTransfer in parallel on $NP processors"
    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun processSpeciesTransfer -parallel > processSpeciesTransferT.out
    else
	mpirun -np $NP processSpeciesTransfer -parallel > processSpeciesTransferT.out
    fi
else
    echo -e "processSpeciesTransfer"
    processSpeciesTransfer > processSpeciesTransferT.out
fi
