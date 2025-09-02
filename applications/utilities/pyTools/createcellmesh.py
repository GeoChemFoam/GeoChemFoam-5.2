import numpy as np
import time

def main(xDim, yDim, zDim, xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, padWidth, direction, NPX, NPY, NPZ, rank):
 
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

 if rank==0:
  with open('system/Nx1', 'w') as f: f.write(f"{nx1}\n")
  with open('system/Ny1', 'w') as f: f.write(f"{ny1}\n")
  with open('system/Nz1', 'w') as f: f.write(f"{nz1}\n")
         
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
 ncells = nxp*nyp*nzp  
                                                      
 ###################################################################
 ###### cellProcAddressing #########################################
 ###################################################################              
 
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
 "    class       labelList;\n",
 "    location    \"constant/polyMesh\";\n",
 "    object      cellProcAddressing;\n",
 "}\n",
 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
 "\n",
 "\n",
 f"{ncells}\n",
 "(\n"
 ]
                 
 # Precompute invariants
 ny_factor = ny1 * nx1
 nx_factor = nx1
 base_index = (
 startZ * ny_factor +
 startY * nx_factor +
 startX
 )
                 
 # Generate data
 for k in range(0,nzp):
   k_offset = k * ny_factor
   for j in range(0,nyp):
     j_offset = j * nx_factor
     for i in range(0,nxp):
       index = base_index + k_offset + j_offset + i
       data.append(f"{index}\n")
                
 # Add footer
 data.append(")\n")
 data.extend(['\n', '\n', "// ************************************************************************* //"])
              
 #io_start_time = time.time()
 
 # Write data to file
 with open(f'processor{rank}/constant/polyMesh/cellProcAddressing', 'w') as f:
   f.writelines(data)
   f.close()
 
 #io_end_time = time.time()
 #io_elapsed_time = io_end_time - io_start_time
 #end_time = io_end_time
 #elapsed_time = end_time - start_time
 #if(rank==0):
   #print("Elapsed time and io time cellmesh:", rank, elapsed_time, io_elapsed_time, "seconds")
 
 
 
