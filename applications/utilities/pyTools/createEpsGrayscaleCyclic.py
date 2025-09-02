import numpy as np
import array
import os
import time

def main(xDim, yDim, zDim, xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, Image_name, padWidth, pores_value,solid_value, eps_min,dimension,direction, NPX, NPY, NPZ, rank, output_path):

 #start_time = time.time()
 
 
 # NPX, NPY and NPZ are the number of processors in the x, y, and z directions.
 # xDim, yDim and zDim are the extents of the image file
 # xMin, xMax, etc., are the lower/upper bounds of the cropped image xMin <= x < xMax
 # nX, nY and nZ are are the number of cells of the initial mesh
 # GJP I guess nX=xMax-xMin but then why define them separately in the calling routine
 
 if direction==0:
    #bounding box
    xmin=-padWidth
    xmax=xMax-xMin+padWidth
    ymin=0
    ymax=yMax-yMin
    zmin=0
    zmax=zMax-zMin
    pad_in_z=0
 elif direction==1:
    #bounding box
    xmin=0
    xmax=xMax-xMin
    ymin=-padWidth
    ymax=yMax-yMin+padWidth
    zmin=0
    zmax=zMax-zMin
    pad_in_z=0
 else:
    direction=2
    #bounding box
    xmin=0
    xmax=xMax-xMin
    ymin=0
    ymax=yMax-yMin
    zmin=-padWidth
    zmax=zMax-zMin+padWidth
    direction=2
    pad_in_z=padWidth
 
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
 ncells=nxp*nyp*nzp
 
 # if padding in z prepare a slice of pore_values
 if direction ==2:
    my_local_pad_array = np.full((yDim*xDim), pores_value, dtype=np.uint8)
 
 io_start_time = time.time()
 
 # open image
 f = open('constant/triSurface/'+Image_name+'.raw', 'rb')
 file_path = f.name
 num_bytes = int(os.path.getsize(file_path) / float(zDim))
 
 io_end_time = time.time()
 io_elapsed_time = io_end_time - io_start_time
 
 # run through every layer, included potential padding in z
 for global_layer in range(0,zmax-zmin):

   ##get subarray for rank

   # ipz is the z-coord of the proc
   # ipy is the y-coord of the proc
   # ipz is the z-coord of the proc
   # nzp is the number of cells in the z-dir for each proc
   # nyp is the number of cells in the z-dir for each proc
   # nxp is the nubmer of cells in the x-dir for each proc

   if (ipz==NPZ-1) and (NPZ>1):
       if (0 <= global_layer < r):

         if global_layer < pad_in_z or global_layer > zMax-zMin+pad_in_z-1:
           my_array = my_local_pad_array
         else:
           offset=(zMin+global_layer-pad_in_z)*xDim*yDim

           f.seek(offset)
           #io_start_time = time.time()
           # read layer from file into local my_array slice
           my_array = np.fromfile(f, dtype=np.uint8, count=num_bytes)

           #io_end_time = time.time()

           #io_elapsed_time = io_elapsed_time + (io_end_time - io_start_time)

           # convert 1D my_array to a 2D array
           my_array = np.reshape(my_array,(yDim,xDim))
           # crop that 2D array
           my_array = my_array[yMin:yMax,xMin:xMax]

           if direction==0:
             #add pad in direction 0
             my_array = np.pad(my_array,pad_width=((0,0),(padWidth,padWidth)),mode='constant',constant_values=pores_value)
           elif direction==1:
             #add pad in direction 1
             my_array = np.pad(my_array,pad_width=((padWidth,padWidth),(0,0)),mode='constant',constant_values=pores_value)

           my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q, startX*p:(startX+nxp)*p], (1,nyp*q,nxp*p))
           if (global_layer==0):
             my_array_p_z01 = my_array_3d
           else:
             my_array_p_z01 = np.concatenate((my_array_p_z01, my_array_3d), axis=0)

   if (ipz==0) and (NPZ>1):
       if (zmax-zmin-r <= global_layer < zmax-zmin):

         if global_layer < pad_in_z or global_layer > zMax-zMin+pad_in_z-1:
           my_array = my_local_pad_array
         else:
           offset=(zMin+global_layer-pad_in_z)*xDim*yDim

           f.seek(offset)
           #io_start_time = time.time()
           # read layer from file into local my_array slice
           my_array = np.fromfile(f, dtype=np.uint8, count=num_bytes)

           #io_end_time = time.time()

           #io_elapsed_time = io_elapsed_time + (io_end_time - io_start_time)

           # convert 1D my_array to a 2D array
           my_array = np.reshape(my_array,(yDim,xDim))
           # crop that 2D array
           my_array = my_array[yMin:yMax,xMin:xMax]

           if direction==0:
             #add pad in direction 0
             my_array = np.pad(my_array,pad_width=((0,0),(padWidth,padWidth)),mode='constant',constant_values=pores_value)
           elif direction==1:
             #add pad in direction 1
             my_array = np.pad(my_array,pad_width=((padWidth,padWidth),(0,0)),mode='constant',constant_values=pores_value)

           my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q, startX*p:(startX+nxp)*p], (1,nyp*q,nxp*p))
           if (global_layer==zmax-zmin-r):
             my_array_p_z10 = my_array_3d
           else:
             my_array_p_z10 = np.concatenate((my_array_p_z10, my_array_3d), axis=0)



   # np.concatenate cannot add layer onto pile if pile is empty, thus the first layer is a simply copy.
   if (startZ-1)*r <= global_layer < (startZ+nzp+1)*r:

     if global_layer < pad_in_z or global_layer > zMax-zMin+pad_in_z-1:
       my_array = my_local_pad_array
     else:
       offset=(zMin+global_layer-pad_in_z)*xDim*yDim

       f.seek(offset)
       #io_start_time = time.time()
       # read layer from file into local my_array slice
       my_array = np.fromfile(f, dtype=np.uint8, count=num_bytes)

       #io_end_time = time.time()

       #io_elapsed_time = io_elapsed_time + (io_end_time - io_start_time)

     # convert 1D my_array to a 2D array
     my_array = np.reshape(my_array,(yDim,xDim))
     # crop that 2D array
     my_array = my_array[yMin:yMax,xMin:xMax]

     if direction==0:
       #add pad in direction 0
       my_array = np.pad(my_array,pad_width=((0,0),(padWidth,padWidth)),mode='constant',constant_values=pores_value)
     elif direction==1:
       #add pad in direction 1 
       my_array = np.pad(my_array,pad_width=((padWidth,padWidth),(0,0)),mode='constant',constant_values=pores_value)


   # np.concatenate cannot add layer onto pile if pile is empty, thus the first layer is a simply copy.
   if startZ*r <= global_layer < (startZ+nzp)*r:
     my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q, startX*p:(startX+nxp)*p], (1,nyp*q,nxp*p))
     if global_layer == startZ*r:
       my_array_p = my_array_3d
     else:
       my_array_p = np.concatenate((my_array_p, my_array_3d), axis=0)

     #proc j to j-1
     if (ipy>0):
       my_array_3d = np.reshape(my_array[(startY-1)*q:startY*q,startX*p:(startX+nxp)*p], (1,q,nxp*p))
       if global_layer == startZ*r:
         my_array_p_y0 = my_array_3d
       else:
         my_array_p_y0 = np.concatenate((my_array_p_y0, my_array_3d), axis=0)

     if (ipy==NPY-1) and (NPY>1):
       my_array_3d = np.reshape(my_array[0:q,startX*p:(startX+nxp)*p], (1,q,nxp*p))
       if global_layer == startZ*r:
         my_array_p_y01 = my_array_3d
       else:
         my_array_p_y01 = np.concatenate((my_array_p_y01, my_array_3d), axis=0)

     #proc i to i-1
     if (ipx>0):
       my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q,(startX-1)*p:startX*p], (1,nyp*q,p))
       if global_layer == startZ*r:
         my_array_p_x0 = my_array_3d
       else:
         my_array_p_x0 = np.concatenate((my_array_p_x0, my_array_3d), axis=0)

     if (ipx==NPX-1) and (NPX>1):
       my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q,0:p], (1,nyp*q,p))
       if global_layer == startZ*r:
         my_array_p_x01 = my_array_3d
       else:
         my_array_p_x01 = np.concatenate((my_array_p_x01, my_array_3d), axis=0)

     #proc i to i+1
     if (ipx<NPX-1):
       my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q,(startX+nxp)*p:(startX+nxp+1)*p], (1,nyp*q,p))
       if global_layer == startZ*r:
         my_array_p_x1 = my_array_3d
       else:
         my_array_p_x1 = np.concatenate((my_array_p_x1, my_array_3d), axis=0)

     if (ipx==0) and (NPX>1):
       my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q,xmax-xmin-p:xmax-xmin], (1,nyp*q,p))
       if global_layer == startZ*r:
         my_array_p_x10 = my_array_3d
       else:
         my_array_p_x10 = np.concatenate((my_array_p_x10, my_array_3d), axis=0)

     #proc j to j+1
     if (ipy<NPY-1):
       my_array_3d = np.reshape(my_array[(startY+nyp)*q:(startY+nyp+1)*q,startX*p:(startX+nxp)*p], (1,q,nxp*p))
       if global_layer == startZ*r:
         my_array_p_y1 = my_array_3d
       else:
         my_array_p_y1 = np.concatenate((my_array_p_y1, my_array_3d), axis=0)

     if (ipy==0) and (NPY>1):
       my_array_3d = np.reshape(my_array[ymax-ymin-q:ymax-ymin,startX*p:(startX+nxp)*p], (1,q,nxp*p))
       if global_layer == startZ*r:
         my_array_p_y10 = my_array_3d
       else:
         my_array_p_y10 = np.concatenate((my_array_p_y10, my_array_3d), axis=0)

   if (startZ-1)*r <= global_layer < startZ*r:
     #proc k to k-1
     if (ipz>0):
       my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q,startX*p:(startX+nxp)*p], (1,nyp*q,nxp*p))
       if global_layer == (startZ-1)*r:
         my_array_p_z0 = my_array_3d
       else:
         my_array_p_z0 = np.concatenate((my_array_p_z0, my_array_3d), axis=0)

   if (startZ+nzp)*r <= global_layer < (startZ+nzp+1)*r:
     #proc k to k+1
     if (ipz<NPZ-1):
       my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q,startX*p:(startX+nxp)*p], (1,nyp*q,nxp*p))
       if global_layer == (startZ+nzp)*r:
         my_array_p_z1 = my_array_3d
       else:
         my_array_p_z1 = np.concatenate((my_array_p_z1, my_array_3d), axis=0)

 ###################################################################
 ###### eps ########################################################
 ###################################################################

 ##get subarray for rank

 eps=np.zeros((nxp,nyp,nzp),dtype="float64")
 for i in range (0,nxp):
   for ii in range (0,p):
     for j in range (0,nyp):
       for jj in range (0,q):
         for k in range (0,nzp):
           for kk in range (0,r):
             eps[i,j,k]+=((my_array_p[k*r+kk,j*q+jj,i*p+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p[k*r+kk,j*q+jj,i*p+ii])/(pores_value-solid_value))/p/q/r
 ##get subarray for proc i-1
 eps_bx0 =np.zeros((nyp,nzp),dtype="float64")
 ##get subarray for proc right 
 eps_bx01 =np.zeros((nyp,nzp),dtype="float64")
 ##get subarray for proc i+1
 eps_bx1 =np.zeros((nyp,nzp),dtype="float64")
 ##get subarray for proc left 
 eps_bx10 =np.zeros((nyp,nzp),dtype="float64")
 ##get subarray for proc j-1
 eps_by0 =np.zeros((nxp,nzp),dtype="float64")
 ##get subarray for proc trhough top 
 eps_by01=np.zeros((nxp,nzp),dtype="float64")
 ##get subarray for proc j+1
 eps_by1 =np.zeros((nxp,nzp),dtype="float64")
 ##get subarray for proc trhough bottom 
 eps_by10=np.zeros((nxp,nzp),dtype="float64")
 ##get subarray for proc k-1
 eps_bz0 =np.zeros((nxp,nyp),dtype="float64")
 ##get subarray for proc trhough front 
 eps_bz01=np.zeros((nxp,nyp),dtype="float64")
 ##get subarray for proc k-1
 eps_bz1 =np.zeros((nxp,nyp),dtype="float64")
 ##get subarray for proc trhough back
 eps_bz10=np.zeros((nxp,nyp),dtype="float64")


 #proc k to k-1
 if (ipz>0):
   for i in range (0,nxp):
     for ii in range (0,p):
       for j in range (0,nyp):
         for jj in range (0,q):
           for kk in range (0,r):
             eps_bz0[i,j]+=((my_array_p_z0[kk,j*q+jj,i*p+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_z0[kk,j*q+jj,i*p+ii])/(pores_value-solid_value))/p/q/r

 #proc through front
 if (ipz==NPZ-1) and (NPZ>1):
   for i in range (0,nxp):
     for ii in range (0,p):
       for j in range (0,nyp):
         for jj in range (0,q):
           for kk in range (0,r):
             eps_bz01[i,j]+=((my_array_p_z01[kk,j*q+jj,i*p+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_z01[kk,j*q+jj,i*p+ii])/(pores_value-solid_value))/p/q/r

 #proc j to j-1
 if (ipy>0):
   for i in range (0,nxp):
     for ii in range (0,p):
       for jj in range (0,q):
         for k in range (0,nzp):
           for kk in range (0,r):
             eps_by0[i,k]+=((my_array_p_y0[k*r+kk,jj,i*p+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_y0[k*r+kk,jj,i*p+ii])/(pores_value-solid_value))/p/q/r

 #proc through top 
 if (ipy==NPY-1) and (NPY>1):
   for i in range (0,nxp):
     for ii in range (0,p):
       for jj in range (0,q):
         for k in range (0,nzp):
           for kk in range (0,r):
             eps_by01[i,k]+=((my_array_p_y01[k*r+kk,jj,i*p+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_y01[k*r+kk,jj,i*p+ii])/(pores_value-solid_value))/p/q/r

 #proc i to i-1
 if (ipx>0):
   for ii in range (0,p):
     for j in range (0,nyp):
       for jj in range (0,q):
         for k in range (0,nzp):
           for kk in range (0,r):
             eps_bx0[j,k]+=((my_array_p_x0[k*r+kk,j*q+jj,ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_x0[k*r+kk,j*q+jj,ii])/(pores_value-solid_value))/p/q/r

 #proc trough right 
 if (ipx==NPX-1) and (NPX>1):
   for ii in range (0,p):
     for j in range (0,nyp):
       for jj in range (0,q):
         for k in range (0,nzp):
           for kk in range (0,r):
             eps_bx01[j,k]+=((my_array_p_x01[k*r+kk,j*q+jj,ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_x01[k*r+kk,j*q+jj,ii])/(pores_value-solid_value))/p/q/r

 #proc i to i+1
 if (ipx<NPX-1):
   for ii in range (0,p):
     for j in range (0,nyp):
       for jj in range (0,q):
         for k in range (0,nzp):
           for kk in range (0,r):
             eps_bx1[j,k]+=((my_array_p_x1[k*r+kk,j*q+jj,ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_x1[k*r+kk,j*q+jj,ii])/(pores_value-solid_value))/p/q/r

 #proc trough left 
 if (ipx==0) and (NPX>1):
   for ii in range (0,p):
     for j in range (0,nyp):
       for jj in range (0,q):
         for k in range (0,nzp):
           for kk in range (0,r):
             eps_bx10[j,k]+=((my_array_p_x10[k*r+kk,j*q+jj,ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_x10[k*r+kk,j*q+jj,ii])/(pores_value-solid_value))/p/q/r


 #proc j to j+1
 if (ipy<NPY-1):
   for i in range (0,nxp):
     for ii in range (0,p):
       for jj in range (0,q):
         for k in range (0,nzp):
           for kk in range (0,r):
             eps_by1[i,k]+=((my_array_p_y1[k*r+kk,jj,i*p+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_y1[k*r+kk,jj,i*p+ii])/(pores_value-solid_value))/p/q/r

 #proc through bottom 
 if (ipy==0) and (NPY>1):
   for i in range (0,nxp):
     for ii in range (0,p):
       for jj in range (0,q):
         for k in range (0,nzp):
           for kk in range (0,r):
             eps_by10[i,k]+=((my_array_p_y10[k*r+kk,jj,i*p+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_y10[k*r+kk,jj,i*p+ii])/(pores_value-solid_value))/p/q/r


 #proc k to k+1
 if (ipz<NPZ-1):
   for i in range (0,nxp):
     for ii in range (0,p):
       for j in range (0,nyp):
         for jj in range (0,q):
           for kk in range (0,r):
             eps_bz1[i,j]+=((my_array_p_z1[kk,j*q+jj,i*p+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_z1[kk,j*q+jj,i*p+ii])/(pores_value-solid_value))/p/q/r

 #proc through front
 if (ipz==0) and (NPZ>1):
   for i in range (0,nxp):
     for ii in range (0,p):
       for j in range (0,nyp):
         for jj in range (0,q):
           for kk in range (0,r):
             eps_bz10[i,j]+=((my_array_p_z10[kk,j*q+jj,i*p+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_z10[kk,j*q+jj,i*p+ii])/(pores_value-solid_value))/p/q/r


 # Initialize list for file content
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
 "    class       volScalarField;\n",
 "    location    \"0\";\n",
 "    object      eps;\n",
 "}\n",
 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
 "\n",
 "dimensions      [0 0 0 0 0 0 0];\n",
 "\n",
 "internalField   nonuniform List<scalar> \n",
 f"{ncells}\n",
 "(\n"
 ]
 for k in range(nzp):
   for j in range(nyp):
     for i in range(nxp):
       num_string = format(eps[i, j, k], ".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
 data.extend([
 ")\n",
 ";\n",
 "\n",
 "boundaryField\n",
 "{\n"
 "    left\n"
 "    {\n",
 "        type            cyclic;\n",
 "    }\n"
 "    right\n"
 "    {\n",
 "        type            cyclic;\n",
 "    }\n"
 "    bottom\n"
 "    {\n",
 "        type            cyclic;\n",
 "    }\n"
 "    top\n"
 "    {\n",
 "        type            cyclic;\n",
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
     data.extend([
     "    front\n"
     "    {\n",
     "        type            cyclic;\n",
     "    }\n"
     "    back\n"
     "    {\n",
     "        type            cyclic;\n",
     "    }\n"
     ])

 ##proc rank to rank-NPX*NPY*(NPZ-1) throughback
 if (ipz==NPZ-1) and (NPZ>2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX*NPY*(NPZ-1)}throughback\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nyp}\n",
   "(\n"
   ])
   for i in range(nxp):
     for j in range(nyp):
       num_string = format(eps_bz01[i,j],".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])

 ##proc k to k-1
 if (ipz>0):
   data.extend([
   f"    procBoundary{rank}to{rank - NPX * NPY}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nyp}\n",
   "(\n"
   ])
   for j in range(nyp):
     for i in range(nxp):
       num_string = format(eps_bz0[i, j], ".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])

  ##proc rank to rank-NPX*NPY*(NPZ-1) throughback
 if (ipz==NPZ-1) and (NPZ==2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX*NPY*(NPZ-1)}throughback\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nyp}\n",
   "(\n"
   ])
   for i in range(nxp):
     for j in range(nyp):
       num_string = format(eps_bz01[i,j],".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])


  ##proc rank to rank-NPX*(NPY-1) throughtop
 if (ipy==NPY-1) and (NPY>2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX*(NPY-1)}throughtop\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nzp}\n",
   "(\n"
   ])
   for i in range(nxp):
     for k in range(nzp):
       num_string = format(eps_by01[i,k],".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])

 ##proc j to j-1
 if (ipy>0):
   data.extend([
   f"    procBoundary{rank}to{rank - NPX}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nzp}\n",
   "(\n"
   ])
   for k in range(nzp):
     for i in range(nxp):
       num_string = format(eps_by0[i, k], ".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])

  ##proc rank to rank-NPX*(NPY-1) throughtop
 if (ipy==NPY-1) and (NPY==2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX*(NPY-1)}throughtop\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nzp}\n",
   "(\n"
   ])
   for i in range(nxp):
     for k in range(nzp):
       num_string = format(eps_by01[i,k],".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])


 ##proc rank to rank-NPX+1 through right
 if (ipx==NPX-1) and (NPX>2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX+1}throughright\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nyp * nzp}\n",
   "(\n"
   ])
   for k in range(nzp):
     for j in range(nyp):
       num_string = format(eps_bx01[j,k],".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])


 ##proc i to i-1
 if (ipx>0):
   data.extend([
   f"    procBoundary{rank}to{rank - 1}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nyp * nzp}\n",
   "(\n"
   ])
   for k in range(nzp):
     for j in range(nyp):
       num_string = format(eps_bx0[j, k], ".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])

 #proc rank to rank-NPX+1 through right
 if (ipx==NPX-1) and (NPX==2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX+1}throughright\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nyp * nzp}\n",
   "(\n"
   ])
   for k in range(nzp):
     for j in range(nyp):
       num_string = format(eps_bx01[j,k],".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])

 ##proc i to i+1
 if (ipx<NPX-1):
   data.extend([
   f"    procBoundary{rank}to{rank + 1}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nyp * nzp}\n",
   "(\n"
   ])
   for k in range(nzp):
     for j in range(nyp):
       num_string = format(eps_bx1[j, k], ".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])
 
 ##proc rank to rank+NPX-1 through left
 if (ipx==0) and (NPX>1):
   data.extend([
   f"    procBoundary{rank}to{rank+NPX-1}throughleft\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nyp * nzp}\n",
   "(\n"
   ])
   for k in range(nzp):
     for j in range(nyp):
       num_string = format(eps_bx10[j,k],".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])


 ##proc j to j+1
 if (ipy<NPY-1):
   data.extend([
   f"    procBoundary{rank}to{rank + NPX}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nzp}\n",
   "(\n"
   ])
   for k in range(nzp):
     for i in range(nxp):
       num_string = format(eps_by1[i, k], ".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])

 ##proc rank to rank+NPX*(NPY-1) through bottom 
 if (ipy==0) and (NPY>1):
   data.extend([
   f"    procBoundary{rank}to{rank+NPX*(NPY-1)}throughbottom\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nzp}\n",
   "(\n"
   ])
   for i in range(nxp):
     for k in range(nzp):
       num_string = format(eps_by10[i,k],".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])

 ##proc k to k+1
 if (ipz<NPZ-1):
   data.extend([
   f"    procBoundary{rank}to{rank + NPX * NPY}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nyp}\n",
   "(\n"
   ])
   for j in range(nyp):
     for i in range(nxp):
       num_string = format(eps_bz1[i, j], ".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
   "    }\n"
   ])

 ##proc rank to rank+NPX*NPY*(NPZ-1) through front 
 if (ipz==0) and (NPZ>1):
   data.extend([
   f"    procBoundary{rank}to{rank+NPX*NPY*(NPZ-1)}throughfront\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nxp * nyp}\n",
   "(\n"
   ])
   for i in range(nxp):
     for j in range(nyp):
       num_string = format(eps_bz10[i,j],".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
   data.extend([
   ")\n",
   ";\n",
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
 with open(output_path+'0/eps', 'w') as f:
   f.writelines(data)
   f.close()

 #io_end_time = time.time()

 #io_elapsed_time = io_elapsed_time + (io_end_time - io_start_time)

 #end_time = io_end_time
 #elapsed_time = end_time - start_time

 #if(rank==0):
 #  print("Elapsed time and io time createEps:", rank, elapsed_time, io_elapsed_time, "seconds")

