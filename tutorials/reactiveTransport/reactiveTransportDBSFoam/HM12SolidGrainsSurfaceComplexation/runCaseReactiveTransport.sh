#!/bin/bash

###### USERS INPUT ############################################################

## Define the total time of the simulation and how often to output concentration fields
TotalTime=0.1
WriteTimestep=0.005
runTimestep=0.0002

#### END OF USER INPUT #######################################################

cp system/controlDict1 system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/runTimestep/$runTimestep/g" system/controlDict

cp system/fvSolution1 system/fvSolution
sed -i "s/cSmooth/0/g" system/fvSolution
sed -i "s/nSmooth/0/g" system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # Run reactiveTransportFoam in parallel
    echo -e "Run reactiveTransportFoam in parallel on $NP processors"
    mpirun -np $NP reactiveTransportDBSFoam -parallel  > reactiveTransportDBSFoamRT.out
else
    echo -e "Run reactiveTransportFoam"
    reactiveTransportDBSFoam > reactiveTransportDBSFoamRT.out
fi

