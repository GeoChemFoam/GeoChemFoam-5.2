import numpy as np
import time

def main(xDim, yDim, zDim, xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, padWidth, direction, NPX, NPY, NPZ, rank, output_path):
 
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
 npoints=(nxp+1)*(nyp+1)*(nzp+1)                                                     
 ###################################################################
 ###### faces  #####################################################
 ###################################################################

 nfaces=(nxp+1)*nyp*nzp+nxp*(nyp+1)*nzp+nxp*nyp*(nzp+1)
 nIfaces=nxp*nyp*(nzp-1)+nxp*(nyp-1)*nzp+(nxp-1)*nyp*nzp
 
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
 "    class       faceList;\n",
 "    location    \"constant/polyMesh\";\n",
 "    object      faces;\n",
 "}\n",
 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
 "\n",
 "\n",
 f"{nfaces}\n",
 "(\n"
 ]
 ##internal faces
 for k in range(0,nzp):
   for j in range(0,nyp):
     for i in range (0,nxp):
       if (i<nxp-1):
         ## i->i+1 face
         x1=k*(nxp+1)*(nyp+1)+j*(nxp+1)+i+1
         x2=k*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i+1
         x3=(k+1)*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i+1
         x4=(k+1)*(nxp+1)*(nyp+1)+j*(nxp+1)+i+1
         data.append(f"4({x1} {x2} {x3} {x4})\n")
       if (j<nyp-1):
         ## j->j+1 face
         x1=k*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i
         x2=(k+1)*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i
         x3=(k+1)*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i+1
         x4=k*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i+1
         data.append(f"4({x1} {x2} {x3} {x4})\n")
       if (k<nzp-1):
         x1=(k+1)*(nxp+1)*(nyp+1)+j*(nxp+1)+i
         x2=(k+1)*(nxp+1)*(nyp+1)+j*(nxp+1)+i+1
         x3=(k+1)*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i+1
         x4=(k+1)*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i
         data.append(f"4({x1} {x2} {x3} {x4})\n")
 if (ipx==0):
   ##Left boundary
   for k in range(0,nzp):
     for j in range(0,nyp):
     ##left face of (jth,kth) column
       x1=k*(nxp+1)*(nyp+1)+j*(nxp+1)
       x2=(k+1)*(nxp+1)*(nyp+1)+j*(nxp+1)
       x3=(k+1)*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)
       x4=k*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)
       data.append(f"4({x1} {x2} {x3} {x4})\n") 
 if (ipx==NPX-1):
   ##right boundary
   for k in range(0,nzp):
     for j in range(0,nyp):
       ##left face of (jth,kth) column
       x1=k*(nxp+1)*(nyp+1)+j*(nxp+1)+nxp
       x2=k*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+nxp
       x3=(k+1)*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+nxp
       x4=(k+1)*(nxp+1)*(nyp+1)+j*(nxp+1)+nxp
       data.append(f"4({x1} {x2} {x3} {x4})\n")
 if (ipy==0):
   ##bottom boundary
   for i in range(0,nxp):
     for k in range(0,nzp):
     ##bottom face of (ith,kth) column
       x1=k*(nxp+1)*(nyp+1)+i
       x2=k*(nxp+1)*(nyp+1)+i+1
       x3=(k+1)*(nxp+1)*(nyp+1)+i+1
       x4=(k+1)*(nxp+1)*(nyp+1)+i
       data.append(f"4({x1} {x2} {x3} {x4})\n")
 if (ipy==NPY-1):
   ##top boundary
   for i in range(0,nxp):
     for k in range(0,nzp):
       ##bottom face of (ith,kth) column
       x1=k*(nxp+1)*(nyp+1)+nyp*(nxp+1)+i
       x2=(k+1)*(nxp+1)*(nyp+1)+nyp*(nxp+1)+i
       x3=(k+1)*(nxp+1)*(nyp+1)+nyp*(nxp+1)+i+1
       x4=k*(nxp+1)*(nyp+1)+nyp*(nxp+1)+i+1
       data.append(f"4({x1} {x2} {x3} {x4})\n")
                           
 if (ipz==0):
   ##front boundary
   for i in range(0,nxp):
     for j in range(0,nyp):
       ##bottom face of (ith,kth) column
       x1=j*(nxp+1)+i
       x2=(j+1)*(nxp+1)+i
       x3=(j+1)*(nxp+1)+i+1
       x4=j*(nxp+1)+i+1
       data.append(f"4({x1} {x2} {x3} {x4})\n")

             
 if (ipz==NPZ-1):        
   ##back boundary
   for i in range(0,nxp):
     for j in range(0,nyp):
       ##bottom face of (ith,kth) column
       x1=nzp*(nxp+1)*(nyp+1)+j*(nxp+1)+i
       x2=nzp*(nxp+1)*(nyp+1)+j*(nxp+1)+i+1
       x3=nzp*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i+1
       x4=nzp*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i
       data.append(f"4({x1} {x2} {x3} {x4})\n")
 
 if (ipz>0):
   ##proc in front
   for j in range(0,nyp):
     for i in range(0,nxp):
       ##bottom face of (ith,kth) column
       x1=j*(nxp+1)+i
       x2=(j+1)*(nxp+1)+i
       x3=(j+1)*(nxp+1)+i+1
       x4=j*(nxp+1)+i+1
       data.append(f"4({x1} {x2} {x3} {x4})\n")
                             
 if (ipy>0):
   ##proc in bottom
   for k in range(0,nzp):
     for i in range(0,nxp):
       ##bottom face of (ith,kth) column
       x1=k*(nxp+1)*(nyp+1)+i
       x2=k*(nxp+1)*(nyp+1)+i+1
       x3=(k+1)*(nxp+1)*(nyp+1)+i+1
       x4=(k+1)*(nxp+1)*(nyp+1)+i
       data.append(f"4({x1} {x2} {x3} {x4})\n")
 
 if (ipx>0):
   ##proc on the left
   for k in range(0,nzp):
     for j in range(0,nyp):
       ##left face of (jth,kth) column
       x1=k*(nxp+1)*(nyp+1)+j*(nxp+1)
       x2=(k+1)*(nxp+1)*(nyp+1)+j*(nxp+1)
       x3=(k+1)*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)
       x4=k*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)
       data.append(f"4({x1} {x2} {x3} {x4})\n")
 if (ipx<NPX-1):
   ##proc on the right
   for k in range(0,nzp):
     for j in range(0,nyp):
       ##left face of (jth,kth) column
       x1=k*(nxp+1)*(nyp+1)+j*(nxp+1)+nxp
       x2=k*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+nxp
       x3=(k+1)*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+nxp
       x4=(k+1)*(nxp+1)*(nyp+1)+j*(nxp+1)+nxp
       data.append(f"4({x1} {x2} {x3} {x4})\n")
 
 if (ipy<NPY-1):
   ##proc on top
   for k in range(0,nzp):
     for i in range(0,nxp):
       ##bottom face of (ith,kth) column
       x1=k*(nxp+1)*(nyp+1)+nyp*(nxp+1)+i
       x2=(k+1)*(nxp+1)*(nyp+1)+nyp*(nxp+1)+i
       x3=(k+1)*(nxp+1)*(nyp+1)+nyp*(nxp+1)+i+1
       x4=k*(nxp+1)*(nyp+1)+nyp*(nxp+1)+i+1
       data.append(f"4({x1} {x2} {x3} {x4})\n")
                         
 if (ipz<NPZ-1):        
   ##proc on back
   for j in range(0,nyp):
     for i in range(0,nxp):
       ##bottom face of (ith,kth) column
       x1=nzp*(nxp+1)*(nyp+1)+j*(nxp+1)+i
       x2=nzp*(nxp+1)*(nyp+1)+j*(nxp+1)+i+1
       x3=nzp*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i+1
       x4=nzp*(nxp+1)*(nyp+1)+(j+1)*(nxp+1)+i
       data.append(f"4({x1} {x2} {x3} {x4})\n")

 data.append(")\n")
 data.extend(['\n', '\n', "// ************************************************************************* //"])

 #io_start_time = time.time()

 # Write all data to file at once
 with open(output_path+'constant/polyMesh/faces', 'w') as f:
   f.writelines(data)
   f.close()
 
 #io_end_time = time.time()
 #io_elapsed_time = io_end_time - io_start_time
 
 if NPX*NPY*NPZ>1:
                 
   ###################################################################
   ###### faceProcAddressing  ######################################
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
   "    object      faceProcAddressing;\n",
   "}\n",
   "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
   "\n",
   "\n",
   f"{nfaces}\n",
   "(\n"
   ]

   ##internal faces
   for k in range(0,nzp-1):
     for j in range(0,nyp-1):
       for i in range (0,nxp-1):
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(i+startX)*3+1}\n")
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(i+startX)*3+2}\n")
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(i+startX)*3+3}\n")
       ##i=nxp
       if (ipx==NPX-1):
         ### right boundary is external in global mesh
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+1}\n")
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+2}\n")
       else:
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+2}\n")
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+3}\n")
     ##j=nyp
     if (ipy==NPY-1):
       ### top boundary is external in global mesh
       for i in range (0,nxp-1):
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*2+1}\n")
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*2+2}\n")
       ##i=nxp
       if (ipx==NPX-1):
         ### right boundary is external in global mesh
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+((nyp-1)+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*2+1}\n")
       else:
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+((nyp-1)+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*2+2}\n")
     else:
       for i in range (0,nxp-1):
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*3+1}\n")
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*3+3}\n")
       ##i=nxp
       if (ipx==NPX-1):
         ### right boundary is external in global mesh
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+((nyp-1)+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+2}\n")
       else:
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+((nyp-1)+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+3}\n")
   ###k=nzp-1
   if (ipz==NPZ-1):
     ###k=back boundary is external in global mesh
     for j in range (0,nyp-1):
       for i in range (0,nxp-1):
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+nx1)+(i+startX)*2+1}\n")
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+nx1)+(i+startX)*2+2}\n")
       ### i=nxp-1
       if (ipx==NPX-1):
         ### right boundary is external in global mesh
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+nx1)+((nxp-1)+startX)*2+1}\n")
       else:
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+nx1)+((nxp-1)+startX)*2+2}\n")
     ### j=nyp-1
     if (ipy==NPY-1):
       for i in range (0,nxp-1):
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+((nyp-1)+startY)*((nx1-1)+nx1)+(i+startX)+1}\n")
     else:
       for i in range (0,nxp-1):
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+((nyp-1)+startY)*((nx1-1)+nx1)+(i+startX)*2+1}\n")
   else:
     for j in range (0,nyp-1):
       for i in range (0,nxp-1):
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(i+startX)*3+1}\n")
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(i+startX)*3+2}\n")
       ### i=nxp-1
       if (ipx==NPX-1):
         ### right boundary is external in global mesh
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+1}\n")
       else:
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+2}\n")
     ### j=nyp-1
     if (ipy==NPY-1):
       for i in range (0,nxp-1):
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+((nyp-1)+startY)*((nx1-1)+2*nx1)+(i+startX)*2+1}\n")
     else:
       for i in range (0,nxp-1):
         data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+((nyp-1)+startY)*((nx1-1)+2*nx1)+(i+startX)*3+1}\n")
                           
   if (ipx==0):
     ##Left boundary
     for k in range(0,nzp):
       for j in range(0,nyp):
         data.append(f"{nx1*ny1*(nz1-1)+nx1*(ny1-1)*nz1+(nx1-1)*ny1*nz1+(k+startZ)*ny1+j+startY+1}\n")
             
   if (ipx==NPX-1):
     ##right boundary
     for k in range(0,nzp):
       for j in range(0,nyp):
         data.append(f"{nx1*ny1*(nz1-1)+nx1*(ny1-1)*nz1+(nx1-1)*ny1*nz1+ny1*nz1+(k+startZ)*ny1+j+startY+1}\n")
   if (ipy==0):
     ##bottom boundary
     for i in range(0,nxp):
       for k in range(0,nzp):
         data.append(f"{nx1*ny1*(nz1-1)+nx1*(ny1-1)*nz1+(nx1-1)*ny1*nz1+2*ny1*nz1+(i+startX)*nz1+k+startZ+1}\n")
                           
   if (ipy==NPY-1):
     ##top boundary
     for i in range(0,nxp):
       for k in range(0,nzp):
         data.append(f"{nx1*ny1*(nz1-1)+nx1*(ny1-1)*nz1+(nx1-1)*ny1*nz1+2*ny1*nz1+nx1*nz1+(i+startX)*nz1+(k+startZ)+1}\n")
                          
   if (ipz==0):
     ##front boundary
     for i in range(0,nxp):
       for j in range(0,nyp):
         data.append(f"{nx1*ny1*(nz1-1)+nx1*(ny1-1)*nz1+(nx1-1)*ny1*nz1+2*ny1*nz1+2*nx1*nz1+(i+startX)*ny1+j+startY+1}\n")
                           
   if (ipz==NPZ-1):
     ##back boundary
     for i in range(0,nxp):
       for j in range(0,nyp):
         data.append(f"{nx1*ny1*(nz1-1)+nx1*(ny1-1)*nz1+(nx1-1)*ny1*nz1+2*ny1*nz1+2*nx1*nz1+nx1*ny1+(i+startX)*ny1+j+startY+1}\n")

   if (ipz>0):
     ##proc k to proc k-1
     for j in range(0,nyp-1):
       for i in range(0,nxp-1):
         data.append(f"{-((startZ-1)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(i+startX)*3+3)}\n")
       ##i=nxp-1
       if (ipx==NPX-1):
         ##right boundary is external in global mesh
         data.append(f"{-((startZ-1)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+2)}\n")
       else:
         data.append(f"{-((startZ-1)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+3)}\n")
     #j=nyp-1
     if (ipy==NPY-1):
       ##top boundary is external in global mesh
       for i in range(0,nxp-1):
         data.append(f"{-((startZ-1)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*2+2)}\n")
       ##i=nxp=1
       if (ipx==NPX-1):
         ##right boundary is external in global mesh
         data.append(f"{-((startZ-1)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*2+1)}\n")
       else:
         data.append(f"{-((startZ-1)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*2+2)}\n")
     else:
       for i in range(0,nxp-1):
         data.append(f"{-((startZ-1)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*3+3)}\n")
       ##i=nxp=1
       if (ipx==NPX-1):
        ##right boundary is external in global mesh
        data.append(f"{-((startZ-1)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+2)}\n")
       else:
        data.append(f"{-((startZ-1)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+3)}\n")
   if (ipy>0):
     ##proc j to proc j-1
     for k in range(0,nzp-1):
       for i in range(0,nxp-1):
         data.append(f"{-((k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(startY-1)*((nx1-1)+2*nx1)+(i+startX)*3+2)}\n")
       ##i=nxp=1
       if (ipx==NPX-1):
         ##right boundary is external in global mesh
         data.append(f"{-((k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(startY-1)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+1)}\n")
       else:
         data.append(f"{-((k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(startY-1)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+2)}\n")
     #j=nyp-1                    
     if (ipz==NPZ-1):
       ##back boundary is external in global mesh                    
       for i in range(0,nxp-1):
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(startY-1)*((nx1-1)+nx1)+(i+startX)*2+2)}\n")
       ##i=nxp=1                            
       if (ipx==NPX-1):
         ##right boundary is external in global mesh
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(startY-1)*((nx1-1)+nx1)+(nxp-1+startX)*2+1)}\n")
       else:
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(startY-1)*((nx1-1)+nx1)+(nxp-1+startX)*2+2)}\n")
     else:
       for i in range(0,nxp-1):
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(startY-1)*((nx1-1)+2*nx1)+(i+startX)*3+2)}\n")
       ##i=nxp=1 
       if (ipx==NPX-1):
         ##right boundary is external in global mesh
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(startY-1)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+1)}\n")
       else:
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(startY-1)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+2)}\n")
   if (ipx>0):
     ##proc i to i-1
     for k in range(0,nzp-1):
       for j in range(0,nyp-1):
         data.append(f"{-((k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(startX-1)*3+1)}\n")
       #j=nyp-1
       if (ipy==NPY-1):
         ##top boundary is external in global mesh
         data.append(f"{-((k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(startX-1)*2+1)}\n")
       else:
         data.append(f"{-((k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(startX-1)*3+1)}\n")
     #k=nzp-1
     if (ipz==NPZ-1):
       #back boundary is external in global mesh
       for j in range(0,nyp-1):
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+nx1)+(startX-1)*2+1)}\n")
       #j=nyp-1
       if (ipy==NPY-1):
         ##top boundary is external in global mesh
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+nx1)+(startX-1)+1)}\n")
       else:
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+nx1)+(startX-1)*2+1)}\n")
     else:
       for j in range(0,nyp-1):
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(startX-1)*3+1)}\n")
       #j=nyp-1
       if (ipy==NPY-1):
         ##top boundary is external in global mesh
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(startX-1)*2+1)}\n")
       else:
         data.append(f"{-((nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(startX-1)*3+1)}\n")
                
   if (ipx<NPX-1):
     ##proc i to proc i+1
     for k in range(0,nzp-1):
       for j in range(0,nyp-1):
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+1}\n")
       #j=nyp-1
       if (ipy==NPY-1):
         ##top boundary is external in global mesh
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*2+1}\n")
       else:
         data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+1}\n")
     #k=nzp-1
     if (ipz==NPZ-1):
       #back boundary is external in global mesh
       for j in range(0,nyp-1):
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+nx1)+(nxp-1+startX)*2+1}\n")
       #j=nyp-1
       if (ipy==NPY-1):
         ##top boundary is external in global mesh
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+nx1)+(nxp-1+startX)+1}\n")
       else:
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+nx1)+(nxp-1+startX)*2+1}\n")
     else:
       for j in range(0,nyp-1):
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+1}\n")
       #j=nyp-1
       if (ipy==NPY-1):
         ##top boundary is external in global mesh
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*2+1}\n")
       else:
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(nxp-1+startX)*3+1}\n")
                   
   if (ipy<NPY-1):
     ##proc j to proc j+1
     for k in range(0,nzp-1):
        for i in range (0,nxp-1):
          data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*3+2}\n")
        ##i=nxp=1                            
        if (ipx==NPX-1):
           ##right boundary is external in global mesh
           data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+1}\n")
        else:
           data.append(f"{(k+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+2}\n")
     #k=nzp-1
     if (ipz==NPZ-1):
       #back boundary is external in global mesh
       for i in range (0,nxp-1):
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+nx1)+(i+startX)*2+2}\n")
       ##i=nxp=1                            
       if (ipx==NPX-1):
         ##right boundary is external in global mesh
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+nx1)+((nxp-1)+startX)*2+1}\n")
       else:
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+nx1)+((nxp-1)+startX)*2+2}\n")
     else:
       for i in range (0,nxp-1):
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*3+2}\n")
       ##i=nxp=1                            
       if (ipx==NPX-1):
         ##right boundary is external in global mesh
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+1}\n")
       else:
         data.append(f"{(nzp-1+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+2}\n")
   if (ipz<NPZ-1):
      ##proc k to proc k+1
      for j in range (0,nyp-1):
        for i in range (0,nxp-1):
          data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+(i+startX)*3+3}\n")
          ##i=nxp=1                            
        if (ipx==NPX-1):
          ##right boundary is external in global mesh
          data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+2}\n")
        else:
          data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(j+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+3}\n")
      #j=nyp-1
      if (ipy==NPY-1):
        ##top boundary is external in global mesh
        for i in range (0,nxp-1):
          data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*2+2}\n")
          ##i=nxp=1                            
        if (ipx==NPX-1):
          ##right boundary is external in global mesh
          data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*2+1}\n")
        else:
          data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*2+2}\n")
      else:
        for i in range (0,nxp-1):
          data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+(i+startX)*3+3}\n")
        ##i=nxp=1                            
        if (ipx==NPX-1):
          ##right boundary is external in global mesh
          data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+2}\n")
        else:
          data.append(f"{((nzp-1)+startZ)*(ny1*(nx1-1)+(ny1-1)*nx1+ny1*nx1)+(nyp-1+startY)*((nx1-1)+2*nx1)+((nxp-1)+startX)*3+3}\n")
 
   data.extend([
   ")\n",
   "\n",
   "\n",
   "// ************************************************************************* //"
   ])

   #io_start_time = time.time()

   with open(f'processor{rank}/constant/polyMesh/faceProcAddressing', 'w') as f:
     f.writelines(data)
     f.close()
 
   #io_end_time = time.time()
   #io_elapsed_time = io_elapsed_time + (io_end_time - io_start_time)
 
 ###################################################################
 ###### owner  #####################################################
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
 f"    note        \"nPoints:{npoints}  nCells:{ncells}  nFaces:{nfaces}  nInternalFaces:{nIfaces}\";\n",
 "    class       labelList;\n",
 "    location    \"constant/polyMesh\";\n",
 "    object      owner;\n",
 "}\n",
 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
 "\n",
 "\n",
 f"{nfaces}\n",
 "(\n"
 ]
             
 # Internal faces
 for k in range(nzp):
   for j in range(nyp):
     for i in range(nxp):
       if (i < nxp - 1):
         ## i->i+1 face
         data.append(f"{k * nxp * nyp + j * nxp + i}\n")
       if (j < nyp - 1):
         ## j->j+1 face
         data.append(f"{k * nxp * nyp + j * nxp + i}\n")
       if (k < nzp - 1):
         ## k->k+1 face
         data.append(f"{k * nxp * nyp + j * nxp + i}\n")

 # Boundary faces
 if (ipx == 0):
   ##Left boundary
   for k in range(nzp):
     for j in range(nyp):
       data.append(f"{k * nxp * nyp + j * nxp}\n")
 if (ipx == NPX - 1):
   ##right boundary
   for k in range(nzp):
     for j in range(nyp):
       data.append(f"{k * nxp * nyp + j * nxp + nxp - 1}\n")
 if (ipy == 0):
   ##bottom boundary
   for i in range(nxp):
     for k in range(nzp):
       data.append(f"{k * nxp * nyp + i}\n")
 if (ipy == NPY - 1):
   ##top boundary
   for i in range(nxp):
     for k in range(nzp):
       data.append(f"{k * nxp * nyp + (nyp - 1) * nxp + i}\n")
 if (ipz == 0):
   ##front boundary
   for i in range(nxp):
     for j in range(nyp):
       data.append(f"{j * nxp + i}\n")
 if (ipz == NPZ - 1):
   ##back boundary
   for i in range(nxp):
     for j in range(nyp):
       data.append(f"{(nzp - 1) * nxp * nyp + j * nxp + i}\n")
             
 # Processor interfaces
 if (ipz > 0):
   #proc in front
   for j in range(nyp):
     for i in range(nxp):
       data.append(f"{j * nxp + i}\n")
 if (ipy > 0):
   ##proc in bottom
   for k in range(nzp):
     for i in range(nxp):
       data.append(f"{k * nxp * nyp + i}\n")
 if (ipx > 0):
   ##proc on the left
   for k in range(nzp):
     for j in range(nyp):
       data.append(f"{k * nxp * nyp + j * nxp}\n")
 if (ipx < NPX - 1):
   ##proc on the right
   for k in range(nzp):
     for j in range(nyp):
       data.append(f"{k * nxp * nyp + j * nxp + nxp - 1}\n")
 if (ipy < NPY - 1):
   ##proc on the top
   for k in range(nzp):
     for i in range(nxp):
       data.append(f"{k * nxp * nyp + (nyp - 1) * nxp + i}\n")
 if (ipz < NPZ - 1):
   ##proc on back
   for j in range(nyp):
     for i in range(nxp):
       data.append(f"{(nzp - 1) * nxp * nyp + j * nxp + i}\n")
              
 data.append(")\n")
 data.extend(['\n', '\n', "// ************************************************************************* //"])
                 
 #io_start_time = time.time()

 # Write data to file all at once
 with open(output_path+'constant/polyMesh/owner', 'w') as f:
   f.writelines(data)
   f.close()
 
 #io_end_time = time.time()
 #io_elapsed_time = io_elapsed_time + (io_end_time - io_start_time)
 
 ###################################################################
 ###### neighbour  #################################################
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
 f"    note        \"nPoints:{npoints}  nCells:{ncells}  nFaces:{nfaces}  nInternalFaces:{nIfaces}\";\n",
 "    class       labelList;\n",
 "    location    \"constant/polyMesh\";\n",
 "    object      neighbour;\n",
 "}\n",
 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
 "\n",
 "\n",
 f"{nIfaces}\n",
 "(\n"
 ]
              
 ##internal faces
 for k in range(nzp):
   for j in range(nyp):
     for i in range(nxp):
       if (i < nxp - 1):
         data.append(f"{k * nxp * nyp + j * nxp + i + 1}\n")
       if (j < nyp - 1):
         data.append(f"{k * nxp * nyp + (j + 1) * nxp + i}\n")
       if (k < nzp - 1):
         data.append(f"{(k + 1) * nxp * nyp + j * nxp + i}\n")
             
 data.append(")\n")
 data.extend(['\n', '\n', "// ************************************************************************* //"])
               
 #io_start_time = time.time()
 
 # Write data to file all at once
 with open(output_path+'constant/polyMesh/neighbour', 'w') as f:
   f.writelines(data)
   f.close()
 
 #io_end_time = time.time()
 #io_elapsed_time = io_elapsed_time + (io_end_time - io_start_time)
 #end_time = io_end_time
 #elapsed_time = end_time - start_time
 #if(rank==0):
   #print("Elapsed time and io time facemesh:", rank, elapsed_time, io_elapsed_time, "seconds")
 
     
 
