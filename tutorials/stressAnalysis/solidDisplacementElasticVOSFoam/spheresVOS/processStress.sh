#!/bin/bash

###### USERS INPUT ############################################################


#### END OF USER INPUT #######################################################


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    echo -e "processStressStrain in parallel on $NP processors"
    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun processStressStrain -parallel > processStress.out
    else
	mpirun -np $NP processStressStrain -parallel > processStress.out
    fi
else
    echo -e "processStressStrain"
    processStressStrain > processStress.out
fi
