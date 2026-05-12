#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

cp system/topoSetDict1 system/topoSetDict


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "Run topoSet in parallel on $NP processors"
    mpirun -np $NP topoSet -parallel  > topoSet.out

    mpirun -np $NP subsetMesh epsCells -overwrite -parallel > subsetMesh.out

    mpirun -np $NP createPatch -overwrite -parallel > createPatch.out

    rm processor*/0/eps

else
    echo "Run topoSet"
    topoSet > topoSet.out

    subsetMesh  epsCells -overwrite > subsetMesh.out

    createPatch -overwrite > createPatch.out

    rm 0/eps

fi

