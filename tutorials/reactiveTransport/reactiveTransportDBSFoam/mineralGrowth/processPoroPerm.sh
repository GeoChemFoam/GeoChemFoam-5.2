#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        # Run simpleFoam in parallel
        echo -e "Run processPoroPerm in parallel on $NP processors"
        srun processPoroPerm -parallel  > processPoroPermRT.out

    else
        # Run processPoroPerm in parallel
        echo -e "Run processPoroPerm in parallel on $NP processors"
        mpirun -np $NP processPoroPerm -parallel  > processPoroPermRT.out
    fi
else
    echo -e "Run processPoroPerm"
    processPoroPerm > processPoroPermRT.out
fi




