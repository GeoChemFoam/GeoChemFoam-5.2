import sys
from mpi4py import MPI
import numpy as np
import array
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dimension=sys.argv[1]

NPX=int(sys.argv[2])
NPY=int(sys.argv[3])
NPZ=int(sys.argv[4])

direction=int(sys.argv[5])
pconf=sys.argv[6]
epsilon1=sys.argv[7]

isSolidWalls=int(sys.argv[8])

fluidPressure=sys.argv[9]

NP=NPX*NPY*NPZ

ipz = rank // (NPX * NPY)  # Get ipz by integer division
remainder = rank % (NPX * NPY)
ipy = remainder // NPX     # Get ipy
ipx = remainder % NPX      # Get ipx

                
###################################################################
###### p ##########################################################
###################################################################
 
if NP>1:
    output_path = "processor"+str(rank)+"/"
else:
    output_path = ""

data = [
'/*--------------------------------*- C++ -*----------------------------------*\\\\n',
'| =========                 |                                                 |\n',
'| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n',
'|  \\\\    /   O peration     | Version:  2212                                  |\n',
'|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n',
'|    \\\\/     M anipulation  |                                                 |\n',
'\\*---------------------------------------------------------------------------*/\n',
'FoamFile\n',
"{\n",
"    version     2.0;\n",
"    format      ascii;\n",
"    arch        \"LSB;label=32;scalar=64\";\n",
"    class       volVectorField;\n",
"    location    \"0\";\n",
"    object      D;\n",
"}\n",
"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
"\n",
"dimensions      [0 1 0 0 0 0 0];\n",
"\n",
"internalField   uniform (0 0 0);\n",
'\n'
"boundaryField"'\n',
"{"'\n',
]
if direction==0:
  data.extend([
"    wall_left\n"
"    {\n",
"        type            fixedValue;\n",
"        value           uniform (0 0 0);\n",
"    }\n"
"    wall_right\n"
"    {\n",
"        type            fixedValue;\n",
f"        value           uniform ({epsilon1} 0 0);\n",
"    }\n"
])
else:
  data.extend([
"    wall_left\n"
"    {\n",
"        type            tractionDisplacement;\n",
"        traction        uniform (0 0 0);\n",
f"        pressure        uniform {pconf};\n",
"        value           uniform (0 0 0);\n",
"    }\n"
"    wall_right\n"
"    {\n",
"        type            tractionDisplacement;\n",
"        traction        uniform (0 0 0);\n",
f"        pressure        uniform {pconf};\n",
"        value           uniform (0 0 0);\n",
"    }\n"
])


if direction==1:
  data.extend([
"    wall_bottom\n"
"    {\n",
"        type            fixedValue;\n",
"        value           uniform (0 0 0);\n",
"    }\n"
"    wall_top\n"
"    {\n",
"        type            fixedValue;\n",
f"        value           uniform (0 {epsilon1} 0);\n",
"    }\n"
])
else:
  data.extend([
"    wall_bottom\n"
"    {\n",
"        type            tractionDisplacement;\n",
"        traction        uniform (0 0 0);\n",
f"        pressure        uniform {pconf};\n",
"        value           uniform (0 0 0);\n",
"    }\n"
"    wall_top\n"
"    {\n",
"        type            tractionDisplacement;\n",
"        traction        uniform (0 0 0);\n",
f"        pressure        uniform {pconf};\n",
"        value           uniform (0 0 0);\n",
"    }\n"
])

if (dimension=="2D"):
   data.extend([
   "    frontAndBack\n"
   "    {\n",
   "        type            empty;\n",
   "    }\n"
   ])

else:
  if direction==2:
    data.extend([
  "    wall_front\n"
  "    {\n",
  "        type            fixedValue;\n",
  "        value           uniform (0 0 0);\n",
  "    }\n"
  "    wall_back\n"
  "    {\n",
  "        type            fixedValue;\n",
  f"        value           uniform (0 0 {epsilon1});\n",
  "    }\n"
  ])
  else:
    data.extend([
  "    wall_front\n"
  "    {\n",
  "        type            tractionDisplacement;\n",
  "        traction        uniform (0 0 0);\n",
  f"        pressure        uniform {pconf};\n",
  "        value           uniform (0 0 0);\n",
  "    }\n"
  "    wall_back\n"
  "    {\n",
  "        type            tractionDisplacement;\n",
  "        traction        uniform (0 0 0);\n",
  f"        pressure        uniform {pconf};\n",
  "        value           uniform (0 0 0);\n",
  "    }\n"
  ])


##proc k to k-1
if (ipz>0):
   data.extend([
   f"    procBoundary{rank}to{rank - NPX * NPY}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform (0 0 0);\n",
   "    }\n"
   ])
            
##proc j to j-1                
if (ipy>0):
   data.extend([
   f"    procBoundary{rank}to{rank - NPX}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform (0 0 0);\n",
   "    }\n"
   ])

##proc i to i-1                
if (ipx>0):
   data.extend([
   f"    procBoundary{rank}to{rank - 1}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform (0 0 0);\n",
   "    }\n"
   ])

##proc i to i+1                
if (ipx<NPX-1):
   data.extend([
   f"    procBoundary{rank}to{rank + 1}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform (0 0 0);\n",
   "    }\n"
   ])

##proc j to j+1                
if (ipy<NPY-1):
   data.extend([
   f"    procBoundary{rank}to{rank + NPX}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform (0 0 0);\n",
   "    }\n"
   ])

##proc k to k+1                
if (ipz<NPZ-1):
   data.extend([
   f"    procBoundary{rank}to{rank + NPX * NPY}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform (0 0 0);\n",
   "    }\n"
   ])

if isSolidWalls==1:
   data.extend([
   f"    solidwalls\n",
   "    {\n",
   "        type            tractionDisplacement;\n",
   "        traction        uniform (0 0 0);\n",
  f"        pressure        uniform {fluidPressure};\n",
   "        value           uniform (0 0 0);\n",
   "    }\n"
   ])

data.extend([
"}\n",
"\n",
"\n",
"// ************************************************************************* //"
])


#io_start_time = time.time()

# Write data to file
with open(output_path+'0/D', 'w') as f:
  f.writelines(data)
  f.close()
