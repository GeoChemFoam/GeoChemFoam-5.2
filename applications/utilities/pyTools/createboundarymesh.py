import numpy as np
import time

def main(xDim, yDim, zDim, xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, padWidth, dimension,direction, NPX, NPY, NPZ, rank, output_path):
   
 #start_time = time.time()
 
 if direction == 0:
 
     # bounding box
     xmin = -padWidth
     xmax = xMax - xMin + padWidth
     ymin = 0
     ymax = yMax - yMin
     zmin = 0
     zmax = zMax - zMin
 
 elif direction == 1:
 
     # bounding box
     xmin = 0
     xmax = xMax - xMin
     ymin = -padWidth
     ymax = yMax - yMin + padWidth
     zmin = 0
     zmax = zMax - zMin
 
 else:
 
     # bounding box
     xmin = 0
     xmax = xMax - xMin
     ymin = 0
     ymax = yMax - yMin
     zmin = -padWidth
     zmax = zMax - zMin + padWidth
 
 # number of cells
 p = int((xMax - xMin) / nX)
 q = int((yMax - yMin) / nY)
 r = int((zMax - zMin) / nZ)
 
 nx1 = int((xmax - xmin) / p)
 ny1 = int((ymax - ymin) / q)
 nz1 = int((zmax - zmin) / r)
 
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

 #nzp = int(np.rint(nz1 / NPZ))
 #nyp = int(np.rint(ny1 / NPY))
 #nxp = int(np.rint(nx1 / NPX))
 
 nIfaces = nxp * nyp * (nzp - 1) + nxp * (nyp - 1) * nzp + (nxp - 1) * nyp * nzp
 ###################################################################
 ###### boundary   #################################################
 ###################################################################

 # Prepare the data to be written
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
 "    class       polyBoundaryMesh;\n",
 "    location    \"constant/polyMesh\";\n",
 "    object      boundary;\n",
 "}\n",
 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
 '\n'
 ]
 
 nbound = 5
 if (dimension=="3D"):
   nbound +=1
 if ipx > 0:
   nbound += 1
 if ipx < NPX - 1:
   nbound += 1
 if ipy > 0:
   nbound += 1
 if ipy < NPY - 1:
   nbound += 1
 if ipz > 0:
   nbound += 1
 if ipz < NPZ - 1:
   nbound += 1
 
 data.append(f"{nbound}\n")
 data.append("(\n")
 
 startFace = nIfaces

 # Add boundary data directly
 if ipx == 0:
   nbfaces = nyp * nzp
 else:
   nbfaces = 0
 if direction == 0:
   data.append("    inlet\n")
 else:
   data.append("    wall_left\n")
 data.append("    {\n")
 data.append("        type            patch;\n")
 data.append(f"        nFaces          {nbfaces};\n")
 data.append(f"        startFace       {startFace};\n")
 data.append("    }\n")
 startFace += nbfaces
 
 if ipx == NPX - 1:
   nbfaces = nyp * nzp
 else:
   nbfaces = 0
 if direction == 0:
   data.append("    outlet\n")
 else:
   data.append("    wall_right\n")
 data.append("    {\n")
 data.append("        type            patch;\n")
 data.append(f"        nFaces          {nbfaces};\n")
 data.append(f"        startFace       {startFace};\n")
 data.append("    }\n")
 startFace += nbfaces
 
 if ipy == 0:
   nbfaces = nxp * nzp
 else:
   nbfaces = 0
 if direction == 1:
   data.append("    inlet\n")
 else:
   data.append("    wall_bottom\n")
 data.append("    {\n")
 data.append("        type            patch;\n")
 data.append(f"        nFaces          {nbfaces};\n")
 data.append(f"        startFace       {startFace};\n")
 data.append("    }\n")
 startFace += nbfaces
 
 if ipy == NPY - 1:
   nbfaces = nxp * nzp
 else:
   nbfaces = 0
 if direction == 1:
   data.append("    outlet\n")
 else:
   data.append("    wall_top\n")
 data.append("    {\n")
 data.append("        type            patch;\n")
 data.append(f"        nFaces          {nbfaces};\n")
 data.append(f"        startFace       {startFace};\n")
 data.append("    }\n")
 startFace += nbfaces
 
 if (dimension=="3D"):
   if ipz == 0:
     nbfaces = nxp * nyp
   else:
     nbfaces = 0
   if direction == 2:
     data.append("    inlet\n")
   else:
     data.append("    wall_front\n")
   data.append("    {\n")
   data.append("        type            patch;\n")
   data.append(f"        nFaces          {nbfaces};\n")
   data.append(f"        startFace       {startFace};\n")
   data.append("    }\n")
   startFace += nbfaces

   if ipz == NPZ - 1:
     nbfaces = nxp * nyp
   else:
     nbfaces = 0
   if direction == 2:
     data.append("    outlet\n")
   else:
     data.append("    wall_back\n")
   data.append("    {\n")
   data.append("        type            patch;\n")
   data.append(f"        nFaces          {nbfaces};\n")
   data.append(f"        startFace       {startFace};\n")
   data.append("    }\n")
   startFace += nbfaces
 else:
   nbfaces = 2*nxp*nyp
   data.append("    frontAndBack\n")
   data.append("    {\n")
   data.append("        type            empty;\n")
   data.append("        inGroups        1(empty);\n")
   data.append(f"        nFaces          {nbfaces};\n")
   data.append(f"        startFace       {startFace};\n")
   data.append("    }\n")
   startFace += nbfaces


 # Add processor boundaries
 if ipz > 0:
   data.append(f"    procBoundary{rank}to{rank - NPX * NPY}\n")
   data.append("    {\n")
   data.append("        type            processor;\n")
   data.append("        inGroups        1(processor);\n")
   data.append(f"        nFaces          {nxp * nyp};\n")
   data.append(f"        startFace       {startFace};\n")
   data.append("        matchTolerance  0.0001;\n")
   data.append("        transform       unknown;\n")
   data.append(f"        myProcNo        {rank};\n")
   data.append(f"        neighbProcNo    {rank - NPX * NPY};\n")
   data.append("    }\n")
   startFace = startFace + nxp * nyp
                 
 if ipy > 0:
   data.append(f"    procBoundary{rank}to{rank - NPX}\n")
   data.append("    {\n")
   data.append("        type            processor;\n")
   data.append("        inGroups        1(processor);\n")
   data.append(f"        nFaces          {nxp * nzp};\n")
   data.append(f"        startFace       {startFace};\n")
   data.append("        matchTolerance  0.0001;\n")
   data.append("        transform       unknown;\n")
   data.append(f"        myProcNo        {rank};\n")
   data.append(f"        neighbProcNo    {rank - NPX};\n")
   data.append("    }\n")
   startFace = startFace + nxp * nzp
                   
 if ipx > 0:
   data.append(f"    procBoundary{rank}to{rank - 1}\n")
   data.append("    {\n")
   data.append("        type            processor;\n")
   data.append("        inGroups        1(processor);\n")
   data.append(f"        nFaces          {nyp * nzp};\n")
   data.append(f"        startFace       {startFace};\n")
   data.append("        matchTolerance  0.0001;\n")
   data.append("        transform       unknown;\n")
   data.append(f"        myProcNo        {rank};\n")
   data.append(f"        neighbProcNo    {rank - 1};\n")
   data.append("    }\n")
   startFace = startFace + nyp * nzp
                   
 if ipx < NPX - 1:
   data.append(f"    procBoundary{rank}to{rank + 1}\n")
   data.append("    {\n")
   data.append("        type            processor;\n")
   data.append("        inGroups        1(processor);\n")
   data.append(f"        nFaces          {nyp * nzp};\n")
   data.append(f"        startFace       {startFace};\n")
   data.append("        matchTolerance  0.0001;\n")
   data.append("        transform       unknown;\n")
   data.append(f"        myProcNo        {rank};\n")
   data.append(f"        neighbProcNo    {rank + 1};\n")
   data.append("    }\n")
   startFace = startFace + nyp * nzp
                   
 if ipy < NPY - 1:
   data.append(f"    procBoundary{rank}to{rank + NPX}\n")
   data.append("    {\n")
   data.append("        type            processor;\n")
   data.append("        inGroups        1(processor);\n")
   data.append(f"        nFaces          {nxp * nzp};\n")
   data.append(f"        startFace       {startFace};\n")
   data.append("        matchTolerance  0.0001;\n")
   data.append("        transform       unknown;\n")
   data.append(f"        myProcNo        {rank};\n")
   data.append(f"        neighbProcNo    {rank + NPX};\n")
   data.append("    }\n")
   startFace = startFace + nxp * nzp
                   
 if ipz < NPZ - 1:
   data.append(f"    procBoundary{rank}to{rank + NPX * NPY}\n")
   data.append("    {\n")
   data.append("        type            processor;\n")
   data.append("        inGroups        1(processor);\n")
   data.append(f"        nFaces          {nxp * nyp};\n")
   data.append(f"        startFace       {startFace};\n")
   data.append("        matchTolerance  0.0001;\n")
   data.append("        transform       unknown;\n")
   data.append(f"        myProcNo        {rank};\n")
   data.append(f"        neighbProcNo    {rank + NPX * NPY};\n")
   data.append("    }\n")
                     
 data.append(")\n")
 data.append('\n')
 data.append("// ************************************************************************* //\n")
                     
 #io_start_time = time.time()
 
 # Write all data at once
 with open(output_path+'constant/polyMesh/boundary', 'a') as f:
   f.seek(0)  # get to the first position
   f.writelines(data)
   f.close()
 
 #io_end_time = time.time()
 #io_elapsed_time = io_end_time - io_start_time
 
 if NPX*NPY*NPZ>1:
   ###################################################################
   ###### boundaryProcAddressing #####################################
   ###################################################################
   # Prepare the data to be written
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
   "    object      boundaryProcAddressing;\n",
   "}\n",
   "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
   '\n'
    ]
    
   if nbound < 11:
     data.append(f"{nbound}(0 1 2 3 4")
     if (dimension=="3D"):
         data.append(" 5")
     if ipz > 0:
       data.append(" -1")
     if ipy > 0:
       data.append(" -1")
     if ipx > 0:
       data.append(" -1")
     if ipx < NPX - 1:
       data.append(" -1")
     if ipy < NPY - 1:
       data.append(" -1")
     if ipz < NPZ - 1:
       data.append(" -1")
     data.append(")\n")
   else:
     data.append('\n')
     data.append(f"{nbound}\n")
     data.append("(\n")
     data.append("0\n")
     data.append("1\n")
     data.append("2\n")
     data.append("3\n")
     data.append("4\n")
     if (dimension=="3D"):
        data.append("5\n")
     if ipz > 0:
       data.append("-1\n")
     if ipy > 0:
       data.append("-1\n")
     if ipx > 0:
       data.append("-1\n")
     if ipx < NPX - 1:
       data.append("-1\n")
     if ipy < NPY - 1:
       data.append("-1\n")
     if ipz < NPZ - 1:
       data.append("-1\n")
     data.append(")\n")

   data.append('\n')
   data.append("// ************************************************************************* //\n")
    
   #io_start_time = time.time()
 
   # Write data to file at once
   with open(f'processor{rank}/constant/polyMesh/boundaryProcAddressing', 'a') as f:
     f.seek(0)  # get to the first position
     f.writelines(data)
     f.close()
 
   #io_end_time = time.time()
   #io_elapsed_time = io_elapsed_time + (io_end_time - io_start_time)
 #end_time = time.time()
 #elapsed_time = end_time - start_time
 #if(rank==0):
   #print("Elapsed time and io time boundarymesh:", rank, elapsed_time, io_elapsed_time, "seconds")
 
