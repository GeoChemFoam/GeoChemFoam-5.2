#!/bin/bash

###### USERS INPUT ############################################################


#### END OF USER INPUT #######################################################


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    echo -e "processHeatTransfer in parallel on $NP processors"
    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun processHeatTransfer -parallel > processHeatTransferH.out
    else
	mpirun -np $NP processHeatTransfer -parallel > processHeatTransferH.out
    fi
else
    echo -e "processHeatTransfer"
    processHeatTransfer > processHeatTransferH.out
fi
