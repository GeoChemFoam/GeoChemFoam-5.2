#!/bin/bash

###### USERS INPUT ############################################################


#### END OF USER INPUT #######################################################


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    echo -e "processDisp in parallel on $NP processors"
    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun processDispersion -parallel  > processDispD.out
    else 
        mpirun -np $NP processDispersion -parallel  > processDispD.out
    fi
else
    echo -e "processDisp"
    processDispersion > processDispD.out
fi
