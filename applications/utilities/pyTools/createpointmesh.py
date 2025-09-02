import numpy as np
import time

def main(xDim, yDim, zDim, xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, res, padWidth, direction, NPX, NPY, NPZ, rank, output_path):
 
 #start_time = time.time()
 
 if direction==0:
 
     #bounding box
     xmin=-padWidth
     xmax=xMax-xMin+padWidth
     ymin=0
     ymax=yMax-yMin
     zmin=0
     zmax=zMax-zMin
     
 elif direction==1:
    
     #bounding box
     xmin=0
     xmax=xMax-xMin
     ymin=-padWidth
     ymax=yMax-yMin+padWidth
     zmin=0
     zmax=zMax-zMin
     
 else:
    
     #bounding box
     xmin=0
     xmax=xMax-xMin
     ymin=0
     ymax=yMax-yMin
     zmin=-padWidth
     zmax=zMax-zMin+padWidth
     
 #number of cells
 p=int((xMax-xMin)/nX)
 q=int((yMax-yMin)/nY)
 r=int((zMax-zMin)/nZ)
 
 nx1=int((xmax-xmin)/p)
 ny1=int((ymax-ymin)/q)
 nz1=int((zmax-zmin)/r)
 
 resX=res*(xmax-xmin)/nx1
 resY=res*(ymax-ymin)/ny1
 resZ=res*(zmax-zmin)/nz1
 
 ipz = rank // (NPX * NPY)  # Get ipz by integer division
 remainder = rank % (NPX * NPY)
 ipy = remainder // NPX     # Get ipy
 ipx = remainder % NPX      # Get ipx
 
 baseZ = nz1 // NPZ
 remZ = nz1 % NPZ

 if ipz < remZ:
    nzp = baseZ + 1
    startZ = ipz * nzp
 else:
    nzp = baseZ
    startZ = remZ * (baseZ + 1) + (ipz - remZ) * baseZ

 baseY = ny1 // NPY
 remY = ny1 % NPY

 if ipy < remY:
    nyp = baseY + 1
    startY = ipy * nyp
 else:
    nyp = baseY
    startY = remY * (baseY + 1) + (ipy - remY) * baseY

 baseX = nx1 // NPX
 remX = nx1 % NPX

 if ipx < remX:
    nxp = baseX + 1
    startX = ipx * nxp 
 else:
    nxp = baseX
    startX = remX * (baseX + 1) + (ipx - remX) * baseX

 
 #nzp=int(np.rint(nz1/NPZ))
 #nyp=int(np.rint(ny1/NPY))
 #nxp=int(np.rint(nx1/NPX))
                                                      
 ###################################################################
 ###### points #####################################################
 ###################################################################
 npoints=(nxp+1)*(nyp+1)*(nzp+1)
                 
 ##So that the files are exactly the same than with blockMesh - remove later
 data = [
 '/*--------------------------------*- C++ -*----------------------------------*\\\\n',
 '| =========                 |                                                 |\n',
 '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n',
 '|  \\\\    /   O peration     | Version:  2212                                  |\n',
 '|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n',
 '|    \\\\/     M anipulation  |                                                 |\n',
 '\\*---------------------------------------------------------------------------*/\n',
 "FoamFile\n",
 "{\n",
 "    version     2.0;\n",
 "    format      ascii;\n",
 "    arch        \"LSB;label=32;scalar=64\";\n",
 "    class       vectorField;\n",
 "    location    \"constant/polyMesh\";\n",
 "    object      points;\n",
 "}\n",
 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
 "\n",
 "\n",
 f"{npoints}\n",
 "(\n"
 ]
 # Collect data in a list
 for k in range(nzp + 1):
  z = round(startZ * resZ + zmin * res + k * resZ, 9)
  num_z = '0' if abs(z) < 1e-9 else str(z)
  for j in range(nyp + 1):
    y = round(startY * resY + ymin * res + j * resY, 9)
    num_y = '0' if abs(y) < 1e-9 else str(y)
    for i in range(nxp + 1):
      x = round(startX * resX + xmin * res + i * resX, 9)
      num_x = '0' if abs(x) < 1e-9 else str(x)
      data.append(f"({num_x} {num_y} {num_z})\n")
 data.append(")\n")
 ##So that the files are exactly the same than with blockMesh - remove later
 data.extend(['\n', '\n', "// ************************************************************************* //"])
 
 #io_start_time = time.time()

 # Write data to file
 with open(output_path+"constant/polyMesh/points", 'w') as f:
     f.writelines(data)

 #io_end_time = time.time()
 #io_elapsed_time = io_end_time - io_start_time
                 
 if NPX*NPY*NPZ>1:
    ###################################################################
   ###### pointProcAddressing ########################################
   ###################################################################
                 
   data = [
   ##So that the files are exactly the same than with blockMesh - remove later
   '/*--------------------------------*- C++ -*----------------------------------*\\\\n',
   '| =========                 |                                                 |\n',
   '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n',
   '|  \\\\    /   O peration     | Version:  2212                                  |\n',
   '|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n',
   '|    \\\\/     M anipulation  |                                                 |\n',
   '\\*---------------------------------------------------------------------------*/\n',
   "FoamFile\n",
   "{\n",
   "    version     2.0;\n",
   "    format      ascii;\n",
   "    arch        \"LSB;label=32;scalar=64\";\n",
   "    class       labelList;\n",
   "    location    \"constant/polyMesh\";\n",
   "    object      pointProcAddressing;\n",
   "}\n",
   "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
   "\n",
   "\n",
   f"{npoints}\n",
   "(\n"
   ]
               
   # Precompute invariants
   ny_factor = (ny1 + 1) * (nx1 + 1)
   nx_factor = (nx1 + 1)
   base_index = (
   startZ * ny_factor +
   startY * nx_factor +
   startX
   )
                
   # Generate data
   for k in range(nzp + 1):
     k_offset = k * ny_factor
     for j in range(nyp + 1):
       j_offset = j * nx_factor
       for i in range(nxp + 1):
         index = base_index + k_offset + j_offset + i
         data.append(f"{index}\n")
                 
   # Add footer
   data.append(")\n")
   data.extend(['\n', '\n', "// ************************************************************************* //"])
                 
   #io_start_time = time.time()
                 
   # Write data to file
   with open(f'processor{rank}/constant/polyMesh/pointProcAddressing', 'w') as f:
     f.writelines(data)
     f.close()
 
   #io_end_time = time.time()
   #io_elapsed_time = io_elapsed_time + (io_end_time - io_start_time)

 #end_time = time.time()
 #elapsed_time = end_time - start_time
 #if(rank==0):
   #print("Elapsed time and io time pointmesh:", rank, elapsed_time, io_elapsed_time, "seconds")
 
