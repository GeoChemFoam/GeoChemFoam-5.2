#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun python $GCFOAM_DIR/applications/utilities/pyTools/moveResultTo0.py "U" "p" "phi" 

        # Run simpleFoam in parallel
        echo -e "Run processPoroPerm in parallel on $NP processors"
        srun processPoroPerm -parallel  > processPoroPermFlow.out

    else
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/moveResultTo0.py "U" "p" "phi"

        # Run processPoroPerm in parallel
        echo -e "Run processPoroPerm in parallel on $NP processors"
        mpirun -np $NP processPoroPerm -parallel  > processPoroPermFlow.out
    fi
else
    for j in [1-9]*; do rm -rf $j/uniform; mv $j/U 0/.; mv $j/p 0/.; mv $j/phi 0/.; done
    rm -rf [1-9]*
    echo -e "Run processPoroPerm"
    processPoroPerm > processPoroPermFlow.out
fi




