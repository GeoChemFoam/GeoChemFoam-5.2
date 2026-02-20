#!/bin/bash

###### USERS INPUT ############################################################

#Define image name
Image_name="calcitePost"
#Image_name="HM12_2_4"

#Image directory location ($PWD if current directory)
dir="$GCFOAM_IMG/raw"

#Choose image format
format='raw'

#Choose if the image is compressed or not
compressed='yes'

# Define image dimensions
x_dim=536
y_dim=300
z_dim=40

segmentation='grayscale'
#Values of solid, pore, and minimum porosity value for the solid phase (note: if the image contains solid voxels, this CANNOT be 0)
pores_value=255
solid_value=0
eps_min=0.0001

# Define cropping parameters
x_min=0
x_max=536
y_min=0
y_max=300
z_min=0
z_max=40

#padding for inlet/outlet
#This adds voxels in the image directly
padWidth=0

# define resolution (m)
res=0.000005

# number of cells of initial mesh
## They need to be a factor of the number of voxels in the image
n_x=134
n_y=75
n_z=10

dynamicMeshRef='true'
#Mesh refinement level
nlevel=1
nRef=200
refineStokes=0

dimension="3D"
#dimension="2D"
# flow direction 0, 1 or 2 
direction=0

# Number of processors in each direction
NPX=2
NPY=2
NPZ=2

#### END OF USER INPUT #######################################################


if [ $format != 'raw' ]
then
        echo "ERROR: only raw format is implemented for this solver"
        exit
fi

if [[ "${dimension}" == "2D" ]]
then
  if [ $direction == 2 ]
  then
     echo "ERROR: direction cannot be 2 in 2D simulations"
     exit
  elif [ $n_z -gt 1 ]
  then
     echo "ERROR: n_z= $n_z must be 1 in 2D simulations"
     exit
  elif [ $NPZ -gt 1 ]
  then
     echo "EROR: NPZ must be 1 in 2D simulations"
     exit
  elif [ $dynamicMeshRef == 'true' ]
  then
     echo "ERROR: no dynamicmeshRef for 2D simulations"
     exit
  fi


fi

#Insert dimensions in postProcessDict
x_1=0
y_1=0
z_1=0

xSize=$(echo "$x_max - $x_min" | bc)
ySize=$(echo "$y_max - $y_min" | bc)
zSize=$(echo "$z_max - $z_min" | bc)

x_2=$(expr $xSize*$res | bc)
y_2=$(expr $ySize*$res | bc)
z_2=$(expr $zSize*$res | bc)

cp system/postProcessDict1 system/postProcessDict
sed -i "s/x_1/$x_1/g" system/postProcessDict
sed -i "s/y_1/$y_1/g" system/postProcessDict
sed -i "s/z_1/$z_1/g" system/postProcessDict

sed -i "s/x_2/$x_2/g" system/postProcessDict
sed -i "s/y_2/$y_2/g" system/postProcessDict
sed -i "s/z_2/$z_2/g" system/postProcessDict

sed -i "s/flowdir/$direction/g" system/postProcessDict

mkdir -p constant/polyMesh

mkdir -p 0

filename=$Image_name\.$format

if [ $compressed == 'yes' ]
then
        filename=$Image_name\.$format\.tar.gz
fi

mkdir -p constant/triSurface
cp $dir\/$filename constant/triSurface/.
cd constant/triSurface
if [ $compressed == 'yes' ]
then
        tar -xf $filename
fi
cd ../..

NP=$(echo "$NPX*$NPY*$NPZ" | bc)

rm -f system/NP*

echo $NPX >> system/NPX
echo $NPY >> system/NPY
echo $NPZ >> system/NPZ

if [ $segmentation == 'grayscale' ]
then
	micro_por=0
	phases=0
elif [ $segmentation == 'phases' ]
then
     echo "ERROR: only grayscale segmentation implemented for this solver"
     exit
fi

cyclic='no'

if [[ $NP > 1 ]]; then
    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then

        # error check that number of process files matches number of MPI tasks in parallel job
        if [[ ${NP} != ${SLURM_NTASKS} ]]
        then
           echo "ERROR: Number of MPI tasks does not equal number of processor files"
           exit
        fi

        # prepend spindle to srun command to pre-load python modules, otherwise comment-out sprindle line to simply use srun
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createMesh.py $x_dim  $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $NPX $NPY $NPZ

    else
       mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createMesh.py $x_dim  $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $NPX $NPY $NPZ

    fi
else
    python $GCFOAM_DIR/applications/utilities/pyTools/createMesh.py $x_dim  $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $NPX $NPY $NPZ

fi

cp constant/dynamicMeshDict1 constant/dynamicMeshDict

cp system/controlDict1 system/controlDict
sed -i "s/TotalTime/1e-06/g" system/controlDict
sed -i "s/WriteTimestep/1e-06/g" system/controlDict
sed -i "s/initTimestep/1e-06/g" system/controlDict
sed -i "s/maxTimestep/1e-06/g" system/controlDict

if [ $dynamicMeshRef == 'true' ]
then

  #Dummy fluid properties
  flowRate=1e-30
  Visc=1e-6
  Diff=1e-30
  kreac=0
  scoeff=1
  rhos=2710
  Mws=100
  cinlet=0
  kf=0 

  cp constant/transportProperties1 constant/transportProperties
  sed -i "s/Visc/$Visc/g" constant/transportProperties
  sed -i "s/rho_s/$rhos/g" constant/transportProperties
  sed -i "s/Mw_s/$Mws/g" constant/transportProperties
  sed -i "s/k_f/$kf/g" constant/transportProperties

  cp constant/thermoPhysicalProperties1 constant/thermoPhysicalProperties
  sed -i "s/Diff/$Diff/g" constant/thermoPhysicalProperties
  sed -i "s/s_coeff/$scoeff/g" constant/thermoPhysicalProperties
  sed -i "s/k_reac/$kreac/g" constant/thermoPhysicalProperties

  cp system/fvSolution1 system/fvSolution
  sed -i "s/nSmooth/1/g" system/fvSolution
  sed -i "s/cSmooth/1/g" system/fvSolution
  sed -i "s/NCORR/1/g" system/fvSolution

  cp constant/dynamicMeshDict2 constant/dynamicMeshDict

  upRef=$(echo "scale=4; 1 - (1-$eps_min)*0.005 + 0.01 * $refineStokes" | bc)
  lowRef=$(echo "scale=4; $eps_min+(1-$eps_min) * 0.005" | bc)

  sed -i "s/nRef/1/g" constant/dynamicMeshDict
  sed -i "s/lowRef/$lowRef/g" constant/dynamicMeshDict
  sed -i "s/upRef/$upRef/g" constant/dynamicMeshDict
  sed -i "s/refLevel/$nlevel/g" constant/dynamicMeshDict

  dt=$(echo "scale=10; 0.000001 / $nlevel" | bc -l)

  cp system/controlDict1 system/controlDict
  sed -i "s/TotalTime/1e-06/g" system/controlDict
  sed -i "s/WriteTimestep/1e-06/g" system/controlDict
  sed -i "s/initTimestep/$dt/g" system/controlDict
  sed -i "s/maxTimestep/1e-06/g" system/controlDict

  if [[ $NP > 1 ]]; then
      # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
      if [[ "${PLATFORM}" == "ARCHER2" ]]; then

          # error check that number of process files matches number of MPI tasks in parallel job
          if [[ ${NP} != ${SLURM_NTASKS} ]]
          then
             echo "ERROR: Number of MPI tasks does not equal number of processor files"
             exit
          fi
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createU.py $dimension $NPX $NPY $NPZ 'flow_rate' $flowRate
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createP.py $dimension $NPX $NPY $NPZ 'flow_rate' 0
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ 'C' $cinlet

	  echo "refine mesh at interface by running reactiveTransportDBSFoam for small time"

          echo "run reactiveTransportDBSFoam in parallel on $NP processor "
          srun reactiveTransportDBSFoam -parallel> reactiveTransportDBSFoam.out

          echo "move new mesh to 0"
          srun python $GCFOAM_DIR/applications/utilities/pyTools/moveMeshTo0.py

          srun python $GCFOAM_DIR/applications/utilities/pyTools/moveResultTo0.py "eps"

          echo "process mesh centers"
          srun processMeshCellCenters -parallel > processMeshCellCenters.out

          echo "create refined eps"
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createEpsDynamic.py $x_dim $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $nlevel $NPX $NPY $NPZ


      else
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createU.py $dimension $NPX $NPY $NPZ 'flow_rate' $flowRate
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createP.py $dimension $NPX $NPY $NPZ 'flow_rate' 0
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ 'C' $cinlet

          echo "refine mesh at interface by running reactiveTransportDBSFoam for small time"

          echo "run reactiveTransportDBSFoam in parallel on $NP processor "
          mpirun -np $NP reactiveTransportDBSFoam -parallel> reactiveTransportDBSFoam.out

	  echo "move new mesh to 0"
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/moveMeshTo0.py

          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/moveResultTo0.py "eps"

          echo "process mesh centers"
          mpirun -np $NP processMeshCellCenters -parallel > processMeshCellCenters.out

          echo "create refined eps"
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createEpsRefinement.py $x_dim $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $nlevel $refineStokes $NPX $NPY $NPZ


      fi
  else
      python $GCFOAM_DIR/applications/utilities/pyTools/createU.py $dimension $NPX $NPY $NPZ 'flow_rate' $flowRate
      python $GCFOAM_DIR/applications/utilities/pyTools/createP.py $dimension $NPX $NPY $NPZ 'flow_rate' 0
      python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ 'C' $cinlet

      echo "refine mesh at interface by running reactiveTransportDBSFoam for small time"

      echo "run reactiveTransportDBSFoam"
      reactiveTransportDBSFoam > reactiveTransportDBSFoam.out

      echo "move new mesh to 0"
      mv 1e-06/polyMesh/* constant/polyMesh/.
      mv 1e-06/eps 0/.
      rm -rf 1e-06

  
      echo "process mesh centers"
      processMeshCellCenters > processMeshCellCenters.out

      python $GCFOAM_DIR/applications/utilities/pyTools/createEpsDynamic.py $x_dim $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $nlevel $NPX $NPY $NPZ
  fi

  cp constant/dynamicMeshDict2 constant/dynamicMeshDict

  sed -i "s/nRef/$nRef/g" constant/dynamicMeshDict
  sed -i "s/lowRef/$lowRef/g" constant/dynamicMeshDict
  sed -i "s/upRef/$upRef/g" constant/dynamicMeshDict
  sed -i "s/refLevel/$nlevel/g" constant/dynamicMeshDict


fi

rm -rf constant/triSurface
  
echo -e "Mesh created. It is advised to check in paraview to confirm mesh and 0/eps are reasonable before running flow"



