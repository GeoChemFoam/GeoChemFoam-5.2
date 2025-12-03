#!/bin/bash

###### USERS INPUT ############################################################

TotalTime=800
WriteTimestep=800
initTimestep=0.1
maxTimestep=50

#### END OF USER INPUT #######################################################

cp system/controlDict1 system/controlDict
cp system/fvSolution1 system/fvSolution

sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/initTimestep/$initTimestep/g" system/controlDict
sed -i "s/maxTimestep/$maxTimestep/g" system/controlDict
sed -i "s/NCORR/100/g" system/fvSolution



if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "run reactiveTransportDBSFoam in parallel on $NP processors"
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun reactiveTransportDBSFoam -parallel > reactiveTransportDBSFoamRT.out
    else
        mpirun -np $NP reactiveTransportDBSFoam -parallel > reactiveTransportDBSFoamRT.out
    fi
else
    echo -e "run reactiveTransportDBSFoam"
    reactiveTransportDBSFoam > reactiveTransportDBSFoamRT.out 
fi

