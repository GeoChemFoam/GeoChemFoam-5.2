#!/bin/bash

###### USERS INPUT ############################################################

#Initial Tracer and boundary condition (dimensionless)
T0=0.0
Tin=1.0
#wall_boundary_type="fixedValue" # This value will be T0
wall_boundary_type="zeroGradient"

#Diffusion coefficient (m^2/s)
Diff=1e-9

#model='molDiff' ##Deff=Diff

model='tortuosity'
## D=Diff/tau

#model='PecletDependent' 
## Dipsersivity constant
## Pe<1
## Dx=Diff/tau*(1+betax*Pe^alpha1x)
## Dy=Diff/tau*(1+betay*Pe^alpha1y)
## Dz=Diff/tau*(1+betaz*Pe^alpha1z)
## Pe>1
## Dx=Diff/tau*(1+betax**Pe^alpha2x)
## Dy=Diff/tau*(1+betay*Pe^alpha2y)
## Dz=Diff/tau*(1+betaz*Pe^alpha2z)


#define the labels of the phases
phases=(0,1,2)

#define the porosity of each phase, note that the porosity of the solid phase CANNOT be 0, default to 0.0001
micro_por=(0.0001,0.4,1)

##tortuosity factor of micropores
tau=(2.5,2.5,1)

##micro pore size (m) - Pe = UL/Diff/eps
#Lpore=(5e-6,5e-6,5e-6)

#betax=(0.5,0.5,0)
#alpha1x=(1,1,1)
#alpha2x=(1,1,1)
#betay=(0.5,0.5,0)
#alpha1y=(1,1,1)
#alpha2y=(1,1,1)
#betaz=(0.5,0.5,0)
#alpha1z=(1,1,1)
#alpha2z=(1,1,1)

#### END OF USER INPUT

echo -e "set flow and transport properties"
cp constant/transportProperties2 constant/transportProperties
sed -i "s/Diff/$Diff/g" constant/transportProperties

NPX="$(tail -n 1 system/NPX)"
NPY="$(tail -n 1 system/NPY)"
NPZ="$(tail -n 1 system/NPZ)"

if { [ -f 0/eps ] && grep -q "frontAndBack" 0/eps; } || \
   { [ -f processor0/0/eps ] && grep -q "frontAndBack" processor0/0/eps; }; then
   dimension="2D"
else dimension="3D"

fi

echo -e "create T"
if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createT.py $dimension $NPX $NPY $NPZ $T0 $Tin $wall_boundary_type "zeroGradient"

	if [[ "$model" == "PecletDependent" ]];then
          echo "Calculate micro dispersivity"
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependent.py $dimension $Diff $micro_por $tau $Lpore $betax $alpha1x $alpha2x  $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ
	elif [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro dispersivity"
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension $Diff $micro_por $tau $NPX $NPY $NPZ
	fi
    else
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createT.py $dimension $NPX $NPY $NPZ $T0 $Tin $wall_boundary_type "zeroGradient"

        if [[ "$model" == "PecletDependent" ]];then
          echo "Calculate micro dispersivity"
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependent.py $dimension $Diff $micro_por $tau $Lpore $betax $alpha1x $alpha2x  $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ 
        elif [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro dispersivity"
          mpirun -np $NP  python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension $Diff $micro_por $tau $NPX $NPY $NPZ
	fi
    fi
else
    python $GCFOAM_DIR/applications/utilities/pyTools/createT.py $dimension $NPX $NPY $NPZ $T0 $Tin $wall_boundary_type "zeroGradient"

    if [[ "$model" == "PecletDependent" ]];then
      echo "Calculate micro dispersivity"
      python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependent.py $dimension $Diff $micro_por $tau $Lpore $betax $alpha1x $alpha2x  $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ
    elif [[ "$model" == "tortuosity" ]];then
      echo "Calculate micro dispersivity"
      python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension $Diff $micro_por $tau $NPX $NPY $NPZ
    fi
fi

echo -e "Case initialised"

