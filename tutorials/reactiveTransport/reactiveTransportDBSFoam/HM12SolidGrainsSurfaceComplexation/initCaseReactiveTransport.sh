#!/bin/bash

###### USERS INPUT ############################################################

#Species in the domain
species=("Ni+2","NiOH+","NiO2H2","Cl-","H+","OH-")

#Diffusion coefficients (m^2/s), one per species
Diff=(1e-9,1e-9,1e-9,1e-9,1e-9,1e-9)

#Initial concentration
C0=(9.838e-06,1.223e-07,9.722e-09,1.994e-05,1.003e-08,1.010e-06)
#inlet concentration
C_in=(0,0,0,0,9.974e-08,1.0095e-07)

#Surface master species and site density (kmol/m2)= 2e-3 mol/kg / 1 m^2/g
Surface_masters=("Surf_sO")
master_density=(2e-9)

#Surface species and initial concentrations (kmol/m2)
Surf_species=("Surf_sO-","Surf_sONi+","Surf_sOH","Surf_sONiOH")
Surf_C0=(4.631e-10,4.715e-10,8.997e-13,1.0644e-09)

# Minimum porosity included in the fixed PHREEQC chemistry mapping
eps_chemistry_min=1e-3

model='molDiff' ##Deff=Diff

## model='tortuosity'
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
#phases=(0,1,2)

#define the porosity of each phase, note that the porosity of the solid phase CANNOT be 0, default to 0.0001
#micro_por=(0.0001,0.4,1)

##tortuosity factor of micropores
#tau=(2.5,2.5,1)

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


IFS=',' read -ra species_array <<< "${species[0]}"
IFS=',' read -ra Diff_array <<< "${Diff[0]}"
IFS=',' read -ra C0_array <<< "${C0[0]}"
IFS=',' read -ra C_in_array <<< "${C_in[0]}"
IFS=',' read -ra Surface_masters_array <<< "${Surface_masters[0]}"
IFS=',' read -ra master_density_array <<< "${master_density[0]}"
IFS=',' read -ra Surf_species_array <<< "${Surf_species[0]}"
IFS=',' read -ra Surf_C0_array <<< "${Surf_C0[0]}"

python $GCFOAM_DIR/applications/utilities/pyTools/createThermoPhysicalProperties.py $species $Diff $Surf_species $Surface_masters $eps_chemistry_min

NPX="$(tail -n 1 system/NPX)"
NPY="$(tail -n 1 system/NPY)"
NPZ="$(tail -n 1 system/NPZ)"

if { [ -f constant/polyMesh/boundary ] && grep -q "frontAndBack" constant/polyMesh/boundary; } || \
   { [ -f processor0/constant/polyMesh/boundary ] && grep -q "frontAndBack" processor0/constant/polyMesh/boundary; }; then
   dimension="2D"
else dimension="3D"

fi

echo -e "create concentration"
if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        for i in "${!species_array[@]}"; do
            srun python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ "${species_array[$i]}" "${C0_array[$i]}" "${C_in_array[$i]}"
        done
        for i in "${!Surface_masters_array[@]}"; do
            srun python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surface_masters_array[$i]}" "${master_density_array[$i]}"
        done
        for i in "${!Surf_species_array[@]}"; do
            srun python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surf_species_array[$i]}" "${Surf_C0_array[$i]}"
        done

	if [[ "$model" == "PecletDependent" ]];then
          echo "Calculate micro dispersivity"
          for i in "${!species_array[@]}"; do
              srun python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependent.py $dimension "${Diff_array[$i]}" $micro_por $tau $Lpore $betax $alpha1x $alpha2x $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ "DT_${species_array[$i]}"
          done
	elif [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro dispersivity"
          for i in "${!species_array[@]}"; do
              srun python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension "${Diff_array[$i]}" $micro_por $tau $NPX $NPY $NPZ "D_${species_array[$i]}"
          done
	fi
    else
        for i in "${!species_array[@]}"; do
            mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ "${species_array[$i]}" "${C0_array[$i]}" "${C_in_array[$i]}"
        done

        for i in "${!Surface_masters_array[@]}"; do
            mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surface_masters_array[$i]}" "${master_density_array[$i]}"
        done
        for i in "${!Surf_species_array[@]}"; do
            mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surf_species_array[$i]}" "${Surf_C0_array[$i]}"
        done
        if [[ "$model" == "PecletDependent" ]];then
          echo "Calculate micro dispersivity"
            for i in "${!species_array[@]}"; do
                mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependent.py $dimension "${Diff_array[$i]}" $micro_por $tau $Lpore $betax $alpha1x $alpha2x $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ "DT_${species_array[$i]}"
            done
        elif  [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro dispersivity"
            for i in "${!species_array[@]}"; do
                mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension "${Diff_array[$i]}" $micro_por $tau $NPX $NPY $NPZ "D_${species_array[$i]}"
            done
	fi
    fi
else
        for i in "${!species_array[@]}"; do
            python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ "${species_array[$i]}" "${C0_array[$i]}" "${C_in_array[$i]}"
        done
        for i in "${!Surface_masters_array[@]}"; do
            python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surface_masters_array[$i]}" "${master_density_array[$i]}"
        done
        for i in "${!Surf_species_array[@]}"; do
            python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surf_species_array[$i]}" "${Surf_C0_array[$i]}"
        done


    if [[ "$model" == "PecletDependent" ]];then
      echo "Calculate micro dispersivity"
      for i in "${!species_array[@]}"; do
          python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependent.py $dimension "${Diff_array[$i]}" $micro_por $tau $Lpore $betax $alpha1x $alpha2x $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ "DT_${species_array[$i]}"
      done
    elif [[ "$model" == "tortuosity" ]];then
      echo "Calculate micro dispersivity"
      for i in "${!species_array[@]}"; do
          python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension "${Diff_array[$i]}" $micro_por $tau $NPX $NPY $NPZ "D_${species_array[$i]}"
      done
    fi
fi



echo -e "Case initialised"
