#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then

        # Run processSolidArea in parallel
        echo -e "Run processSoldiArea in parallel on $NP processors"
        srun processSolidArea -parallel  > processSolidAreaRT.out

    else
        # Run processPoroPerm in parallel
        echo -e "Run processSoldiArea in parallel on $NP processors"
        mpirun -np $NP processSolidArea -parallel  > processSolidAreaRT.out
    fi
else
    echo -e "Run processSolidArea"
    processSolidArea > processSolidArea.out
fi




