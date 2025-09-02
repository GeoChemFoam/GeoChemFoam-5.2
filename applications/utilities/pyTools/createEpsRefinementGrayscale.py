import numpy as np
import array
import os
import time
import sys
from mpi4py import MPI
import h5py

def main(xDim, yDim, zDim, xMin, xMax, yMin, yMax, zMin, zMax, nX, nY, nZ, res,Image_name, padWidth, pores_value,solid_value, eps_min,dimension,direction, nlevel,refineStokes,NPX, NPY, NPZ, rank, output_path):

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
 
 p1=int(p/(2**nlevel))
 q1=int(q/(2**nlevel))
 r1=int(r/(2**nlevel))

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
       my_array_3d = np.reshape(my_array[startY*q-q1:startY*q,startX*p:(startX+nxp)*p], (1,q1,nxp*p))
       if global_layer == startZ*r:
         my_array_p_y0 = my_array_3d
       else:
         my_array_p_y0 = np.concatenate((my_array_p_y0, my_array_3d), axis=0)
     #proc i to i-1
     if (ipx>0):
       my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q,startX*p-p1:startX*p], (1,nyp*q,p1))
       if global_layer == startZ*r:
         my_array_p_x0 = my_array_3d
       else:
         my_array_p_x0 = np.concatenate((my_array_p_x0, my_array_3d), axis=0)
     #proc i to i+1
     if (ipx<NPX-1):
       my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q,(startX+nxp)*p:(startX+nxp)*p+p1], (1,nyp*q,p1))
       if global_layer == startZ*r:
         my_array_p_x1 = my_array_3d
       else:
         my_array_p_x1 = np.concatenate((my_array_p_x1, my_array_3d), axis=0)
     #proc j to j+1
     if (ipy<NPY-1):
       my_array_3d = np.reshape(my_array[(startY+nyp)*q:(startY+nyp)*q+q1,startX*p:(startX+nxp)*p], (1,q1,nxp*p))
       if global_layer == startZ*r:
         my_array_p_y1 = my_array_3d
       else:
         my_array_p_y1 = np.concatenate((my_array_p_y1, my_array_3d), axis=0)


   if startZ*r-r1 <= global_layer < startZ*r:
     #proc k to k-1
     if (ipz>0):
       my_array_3d = np.reshape(my_array[startY*q:(startY+nyp)*q,startX*p:(startX+nxp)*p], (1,nyp*q,nxp*p))
       if global_layer == startZ*r-r1:
         my_array_p_z0 = my_array_3d
       else:
         my_array_p_z0 = np.concatenate((my_array_p_z0, my_array_3d), axis=0)

   if (startZ+nzp)*r <= global_layer < (startZ+nzp)*r+r1:
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

 dx = res*(xmax-xmin)/nx1/(2**nlevel)
 dy = res*(ymax-ymin)/ny1/(2**nlevel)
 dz = res*(zmax-zmin)/nz1/(2**nlevel)

 NP=NPX*NPY*NPZ

 if NP>1:
   output_path = "processor"+str(rank)+"/"
 else:
   output_path = ""
 with open(output_path+"0/eps", 'r') as f:
     lines = f.readlines()

 # Find line with 'internalField' and get the line with the number of entries
 for i, line in enumerate(lines):
     if 'internalField' in line and 'nonuniform' in line:
         uniform=0
         n = int(lines[i + 1].strip())  # number of scalar values
         eps= np.empty(n, dtype=np.float64)
         start_idx = i + 3  # Skip the '(', go to the data
         for j in range(n):
             eps[j] = float(lines[start_idx + j].strip())
         i=start_idx+n
         break
     elif 'internalField' in line and 'uniform' in line:
         uniform=1
         n=1
         eps= np.empty(n, dtype=np.float64)
         parts = lines[i].replace(';', '').split()
         eps[0] = float(parts[-1].rstrip(';'))
         i +=1
         break

 if (ipz>0):
     found_ipz0=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank-NPX*NPY}" in line:
             found_ipz0=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i] and 'nonuniform' in lines[i]:
                     uniform_bz0=0
                     nbz0 = int(lines[i + 1].strip())
                     start_idx = i + 3
                     eps_bz0 = np.empty(nbz0, dtype=np.float64)
                     for j in range(nbz0):
                         eps_bz0[j] = float(lines[start_idx + j].strip())
                     i = start_idx + nbz0  # skip ahead
                     break
                 elif 'value' in lines[i] and 'uniform' in lines[i]:
                     uniform_bz0=1
                     nbz0=1
                     eps_bz0 = np.empty(nbz0, dtype=np.float64)
                     parts = lines[i].replace(';', '').split()
                     eps_bz0[0] = float(parts[-1].rstrip(';'))
                     break
                 i += 1
         i +=1
         if (found_ipz0==1):
             break

 if (ipy>0):
     found_ipy0=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank-NPX}" in line:
             found_ipy0=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i] and 'nonuniform' in lines[i]:
                     uniform_by0=0
                     nby0 = int(lines[i + 1].strip())
                     start_idx = i + 3
                     eps_by0 = np.empty(nby0, dtype=np.float64)
                     for j in range(nby0):
                         eps_by0[j] = float(lines[start_idx + j].strip())
                     i = start_idx + nby0  # skip ahead
                     break
                 elif 'value' in lines[i] and 'uniform' in lines[i]:
                     uniform_by0=1
                     nby0=1
                     eps_by0 = np.empty(nby0, dtype=np.float64)
                     parts = lines[i].replace(';', '').split()
                     eps_by0[0] = float(parts[-1].rstrip(';'))
                     break
                 i += 1
         i +=1
         if (found_ipy0==1):
             break

 if (ipx>0):
     found_ipx0=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank-1}" in line:
             found_ipx0=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i] and 'nonuniform' in lines[i]:
                     uniform_bx0=0
                     nbx0 = int(lines[i + 1].strip())
                     start_idx = i + 3
                     eps_bx0 = np.empty(nbx0, dtype=np.float64)
                     for j in range(nbx0):
                         eps_bx0[j] = float(lines[start_idx + j].strip())
                     i = start_idx + nbx0  # skip ahead
                     break
                 elif 'value' in lines[i] and 'uniform' in lines[i]:
                     uniform_bx0=1
                     nbx0=1
                     eps_bx0 = np.empty(nbx0, dtype=np.float64)
                     parts = lines[i].replace(';', '').split()
                     eps_bx0[0] = float(parts[-1].rstrip(';'))
                     break
                 i += 1
         i +=1
         if (found_ipx0==1):
             break

 if (ipx<NPX-1):
     found_ipx1=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank+1}" in line:
             found_ipx1=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i] and 'nonuniform' in lines[i]:
                     uniform_bx1=0
                     nbx1 = int(lines[i + 1].strip())
                     start_idx = i + 3
                     eps_bx1 = np.empty(nbx1, dtype=np.float64)
                     for j in range(nbx1):
                         eps_bx1[j] = float(lines[start_idx + j].strip())
                     i = start_idx + nbx1  # skip ahead
                     break
                 elif 'value' in lines[i] and 'uniform' in lines[i]:
                     uniform_bx1=1
                     nbx1=1
                     eps_bx1 = np.empty(nbx1, dtype=np.float64)
                     parts = lines[i].replace(';', '').split()
                     eps_bx1[0] = float(parts[-1].rstrip(';'))
                     break
                 i += 1
         i +=1
         if (found_ipx1==1):
             break

 if (ipy<NPY-1):
     found_ipy1=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank+NPX}" in line:
             found_ipy1=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i] and 'nonuniform' in lines[i]:
                     uniform_by1=0
                     nby1 = int(lines[i + 1].strip())
                     start_idx = i + 3
                     eps_by1 = np.empty(nby1, dtype=np.float64)
                     for j in range(nby1):
                         eps_by1[j] = float(lines[start_idx + j].strip())
                     i = start_idx + nby1  # skip ahead
                     break
                 elif 'value' in lines[i] and 'uniform' in lines[i]:
                     uniform_by1=1
                     nby1=1
                     eps_by1 = np.empty(nby1, dtype=np.float64)
                     parts = lines[i].replace(';', '').split()
                     eps_by1[0] = float(parts[-1].rstrip(';'))
                     break
                 i += 1
         i +=1
         if (found_ipy1==1):
             break

 if (ipz<NPZ-1):
     found_ipz1=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank+NPX*NPY}" in line:
             found_ipz1=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i] and 'nonuniform' in lines[i]:
                     uniform_bz1=0
                     nbz1 = int(lines[i + 1].strip())
                     start_idx = i + 3
                     eps_bz1 = np.empty(nbz1, dtype=np.float64)
                     for j in range(nbz1):
                         eps_bz1[j] = float(lines[start_idx + j].strip())
                     i = start_idx + nbz1  # skip ahead
                     break
                 elif 'value' in lines[i] and 'uniform' in lines[i]:
                     uniform_bz1=1
                     nbz1=1
                     eps_bz1 = np.empty(nbz1, dtype=np.float64)
                     parts = lines[i].replace(';', '').split()
                     eps_bz1[0] = float(parts[-1].rstrip(';'))
                     break
                 i += 1
         i +=1
         if (found_ipz1==1):
             break
     i += 1

 with open(output_path+"0/cellCenters", 'r') as f:
     lines = f.readlines()

 # Find line with 'internalField' and get the line with the number of entries
 for i, line in enumerate(lines):
     if 'internalField' in line:
         n = int(lines[i + 1].strip())  # number of scalar values
         if (uniform==1):
             uniform=0
             epsVal=eps[0]
             eps=np.full(n,epsVal,dtype=np.float64)
         ix= np.empty(n, dtype=np.int64)
         iy= np.empty(n, dtype=np.int64)
         iz= np.empty(n, dtype=np.int64)      
         start_idx = i + 3  # Skip the '(', go to the data
         for j in range(n):
             vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
             parts = vec_line.split()
             x, y, z = map(float, parts)
             ix[j]=np.floor((x-res*(xmin+startX*p))/dx);
             iy[j]=np.floor((y-res*(ymin+startY*q))/dy);
             iz[j]=np.floor((z-res*(zmin+startZ*r))/dz);
         i=start_idx+n
         break

 if (ipz>0):
     found_ipz0=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank-NPX*NPY}" in line:
             found_ipz0=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i]:
                     nbz0 = int(lines[i + 1].strip())
                     if (uniform_bz0==1):
                         uniform_bz0=0
                         epsVal=eps_bz0[0]
                         eps_bz0=np.full(nbz0,epsVal,dtype=np.float64)
                     start_idx = i + 3
                     ix_bz0 = np.empty(nbz0, dtype=np.int64)
                     iy_bz0 = np.empty(nbz0, dtype=np.int64)
                     #iz_bz0 = np.empty(nbz0, dtype=np.int64)
                     for j in range(nbz0):
                         vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                         parts = vec_line.split()
                         x, y, z = map(float, parts)
                         ix_bz0[j]=np.floor((x-res*(xmin+startX*p))/dx);
                         iy_bz0[j]=np.floor((y-res*(ymin+startY*q))/dy);
                         #iz_bz0[j]=np.floor((z-res*(zmin+(startZ-1)*r))/dz);
                     i = start_idx + nbz0  # skip ahead
                     break
                 i +=1
         i +=1
         if (found_ipz0==1):
             break

 if (ipy>0):
     found_ipy0=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank-NPX}" in line:
             found_ipy0=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i]:
                     nby0 = int(lines[i + 1].strip())
                     if (uniform_by0==1):
                         uniform_by0=0
                         epsVal=eps_by0[0]
                         eps_by0=np.full(nby0,epsVal,dtype=np.float64)
                     start_idx = i + 3
                     ix_by0 = np.empty(nby0, dtype=np.int64)
                     #iy_by0 = np.empty(nby0, dtype=np.int64)
                     iz_by0 = np.empty(nby0, dtype=np.int64)
                     for j in range(nby0):
                         vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                         parts = vec_line.split()
                         x, y, z = map(float, parts)
                         ix_by0[j]=np.floor((x-res*(xmin+startX*p))/dx);
                         #iy_by0[j]=np.floor((y-res*(ymin+(startY-1)*q))/dy);
                         iz_by0[j]=np.floor((z-res*(zmin+startZ*r))/dz);
                     i = start_idx + nby0  # skip ahead
                     break
                 i += 1
         i +=1
         if (found_ipy0==1):
             break

 if (ipx>0):
     found_ipx0=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank-1}" in line:
             found_ipx0=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i]:
                     nbx0 = int(lines[i + 1].strip())
                     if (uniform_bx0==1):
                         uniform_bx0=0
                         epsVal=eps_bx0[0]
                         eps_bx0=np.full(nbx0,epsVal,dtype=np.float64)
                     start_idx = i + 3
                     #ix_bx0 = np.empty(nbx0, dtype=np.int64)
                     iy_bx0 = np.empty(nbx0, dtype=np.int64)
                     iz_bx0 = np.empty(nbx0, dtype=np.int64)
                     for j in range(nbx0):
                         vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                         parts = vec_line.split()
                         x, y, z = map(float, parts)
                         #ix_bx0[j]=np.floor((x-res*(xmin+(startX-1)*p))/dx);
                         iy_bx0[j]=np.floor((y-res*(ymin+startY*q))/dy);
                         iz_bx0[j]=np.floor((z-res*(zmin+startZ*r))/dz);
                     i = start_idx + nbx0  # skip ahead
                     break
                 i += 1
         i +=1
         if (found_ipx0==1):
             break

 if (ipx<NPX-1):
     found_ipx1=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank+1}" in line:
             found_ipx1=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i]:
                     nbx1 = int(lines[i + 1].strip())
                     if (uniform_bx1==1):
                         uniform_bx1=0
                         epsVal=eps_bx1[0]
                         eps_bx1=np.full(nbx1,epsVal,dtype=np.float64)
                     start_idx = i + 3
                     #ix_bx1 = np.empty(nbx1, dtype=np.int64)
                     iy_bx1 = np.empty(nbx1, dtype=np.int64)
                     iz_bx1 = np.empty(nbx1, dtype=np.int64)
                     for j in range(nbx1):
                         vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                         parts = vec_line.split()
                         x, y, z = map(float, parts)
                         #ix_bx1[j]=np.floor((x-res*(xmin+(startX+nxp)*p))/dx);
                         iy_bx1[j]=np.floor((y-res*(ymin+startY*q))/dy);
                         iz_bx1[j]=np.floor((z-res*(zmin+startZ*r))/dz);
                     i = start_idx + nbx1  # skip ahead
                     break
                 i += 1
         i +=1
         if (found_ipx1==1):
             break

 if (ipy<NPY-1):
     found_ipy1=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank+NPX}" in line:
             found_ipy1=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i]:
                     nby1 = int(lines[i + 1].strip())
                     if (uniform_by1==1):
                         uniform_by1=0
                         epsVal=eps_by1[0]
                         eps_by1=np.full(nby1,epsVal,dtype=np.float64)
                     start_idx = i + 3
                     ix_by1 = np.empty(nby1, dtype=np.int64)
                     #iy_by1 = np.empty(nby1, dtype=np.int64)
                     iz_by1 = np.empty(nby1, dtype=np.int64)
                     for j in range(nby1):
                         vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                         parts = vec_line.split()
                         x, y, z = map(float, parts)
                         ix_by1[j]=np.floor((x-res*(xmin+startX*p))/dx);
                         #iy_by1[j]=np.floor((y-res*(ymin+(startY+nyp)*q))/dy);
                         iz_by1[j]=np.floor((z-res*(zmin+startZ*r))/dz);
                     i = start_idx + nby1  # skip ahead
                     break
                 i += 1
         i +=1
         if (found_ipy1==1):
             break

 if (ipz<NPZ-1):
     found_ipz1=0
     while (i<len(lines)):
         line=lines[i].strip()
         if f"procBoundary{rank}to{rank+NPX*NPY}" in line:
             found_ipz1=1
             # Search for internalField inside this patch
             while (i<len(lines)):
                 if 'value' in lines[i]:
                     nbz1 = int(lines[i + 1].strip())
                     if (uniform_bz1==1):
                         uniform_bz1=0
                         epsVal=eps_bz1[0]
                         eps_bz1=np.full(nbz1,epsVal,dtype=np.float64)
                     start_idx = i + 3
                     ix_bz1 = np.empty(nbz1, dtype=np.int64)
                     iy_bz1 = np.empty(nbz1, dtype=np.int64)
                     #iz_bz1 = np.empty(nbz1, dtype=np.int64)
                     for j in range(nbz1):
                         vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                         parts = vec_line.split()
                         x, y, z = map(float, parts)
                         ix_bz1[j]=np.floor((x-res*(xmin+startX*p))/dx);
                         iy_bz1[j]=np.floor((y-res*(ymin+startY*q))/dy);
                         #iz_bz1[j]=np.floor((z-res*(zmin+(startZ+nzp)*r))/dz);
                     i = start_idx + nbz1  # skip ahead
                     break
                 i += 1
         i +=1
         if (found_ipz1==1):
             break
     i += 1

 ##get subarray for rank

 #eps=np.zeros((nxp,nyp,nzp),dtype="float64")
 for j in range(0,n):
    epsVal=eps[j]
    is_refined=0
    if (epsVal>1.0+1e-3-refineStokes*3e-3):
      is_refined=1
    else:
       if (epsVal<1-(1-eps_min)*5e-3) and (epsVal> eps_min+(1-eps_min)*5e-3):
         is_refined=1
    if (is_refined):
      eps[j]=0
      for ii in range (0,p1):
        for jj in range (0,q1):
          for kk in range (0,r1):
             eps[j]+=((my_array_p[iz[j]*r1+kk,iy[j]*q1+jj,ix[j]*p1+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p[iz[j]*r1+kk,iy[j]*q1+jj,ix[j]*p1+ii])/(pores_value-solid_value))/p1/q1/r1

 if (ipz>0):
   for j in range(0,nbz0):
     epsVal=eps_bz0[j]
     is_refined=0
     if (epsVal>1.0+1e-3-refineStokes*3e-3):
       is_refined=1
     else:
       if (epsVal<1-(1-eps_min)*5e-3) and (epsVal> eps_min+(1-eps_min)*5e-3):
         is_refined=1
     if (is_refined):
       eps_bz0[j]=0
       for ii in range (0,p1):
         for jj in range (0,q1):
           for kk in range (0,r1):
             eps_bz0[j]+=((my_array_p_z0[kk,iy_bz0[j]*q1+jj,ix_bz0[j]*p1+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_z0[kk,iy_bz0[j]*q1+jj,ix_bz0[j]*p1+ii])/(pores_value-solid_value))/p1/q1/r1


 if (ipy>0):
   for j in range(0,nby0):
     epsVal=eps_by0[j]
     is_refined=0
     if (epsVal>1.0+1e-3-refineStokes*3e-3):
       is_refined=1
     else:
       if (epsVal<1-(1-eps_min)*5e-3) and (epsVal> eps_min+(1-eps_min)*5e-3):
         is_refined=1
     if (is_refined):
       eps_by0[j]=0
       for ii in range (0,p1):
         for jj in range (0,q1):
           for kk in range (0,r1):
               eps_by0[j]+=((my_array_p_y0[iz_by0[j]*r1+kk,jj,ix_by0[j]*p1+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_y0[iz_by0[j]*r1+kk,jj,ix_by0[j]*p1+ii])/(pores_value-solid_value))/p1/q1/r1

 if (ipx>0):
   for j in range(0,nbx0):
     epsVal=eps_bx0[j]
     is_refined=0
     if (epsVal>1.0+1e-3-refineStokes*3e-3):
       is_refined=1
     else:
       if (epsVal<1-(1-eps_min)*5e-3) and (epsVal> eps_min+(1-eps_min)*5e-3):
         is_refined=1
     if (is_refined):
       eps_bx0[j]=0
       for ii in range (0,p1):
         for jj in range (0,q1):
           for kk in range (0,r1):
               eps_bx0[j]+=((my_array_p_x0[iz_bx0[j]*r1+kk,iy_bx0[j]*q1+jj,ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_x0[iz_bx0[j]*r1+kk,iy_bx0[j]*q1+jj,ii])/(pores_value-solid_value))/p1/q1/r1

 if (ipx<NPX-1):
   for j in range(0,nbx1):
     epsVal=eps_bx1[j]
     is_refined=0
     if (epsVal>1.0+1e-3-refineStokes*3e-3):
       is_refined=1
     else:
       if (epsVal<1-(1-eps_min)*5e-3) and (epsVal> eps_min+(1-eps_min)*5e-3):
         is_refined=1
     if (is_refined):
       eps_bx1[j]=0
       for ii in range (0,p1):
         for jj in range (0,q1):
           for kk in range (0,r1):
               eps_bx1[j]+=((my_array_p_x1[iz_bx1[j]*r1+kk,iy_bx1[j]*q1+jj,ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_x1[iz_bx1[j]*r1+kk,iy_bx1[j]*q1+jj,ii])/(pores_value-solid_value))/p1/q1/r1

 if (ipy<NPY-1):
   for j in range(0,nby1):
     epsVal=eps_by1[j]
     is_refined=0
     if (epsVal>1.0+1e-3-refineStokes*3e-3):
       is_refined=1
     else:
       if (epsVal<1-(1-eps_min)*5e-3) and (epsVal> eps_min+(1-eps_min)*5e-3):
         is_refined=1
     if (is_refined):
       eps_by1[j]=0
       for ii in range (0,p1):
         for jj in range (0,q1):
           for kk in range (0,r1):
               eps_by1[j]+=((my_array_p_y1[iz_by1[j]*r1+kk,jj,ix_by1[j]*p1+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_y1[iz_by1[j]*r1+kk,jj,ix_by1[j]*p1+ii])/(pores_value-solid_value))/p1/q1/r1
  
 if (ipz<NPZ-1):
   for j in range(0,nbz1):
     epsVal=eps_bz1[j]
     is_refined=0
     if (epsVal>1.0+1e-3-refineStokes*3e-3):
       is_refined=1
     else:
       if (epsVal<1-(1-eps_min)*5e-3) and (epsVal> eps_min+(1-eps_min)*5e-3):
         is_refined=1
     if (is_refined):
       eps_bz1[j]=0
       for ii in range (0,p1):
         for jj in range (0,q1):
           for kk in range (0,r1):
               eps_bz1[j]+=((my_array_p_z1[kk,iy_bz1[j]*q1+jj,ix_bz1[j]*p1+ii]-solid_value)/(pores_value-solid_value)+ eps_min*(pores_value-my_array_p_z1[kk,iy_bz1[j]*q1+jj,ix_bz1[j]*p1+ii])/(pores_value-solid_value))/p1/q1/r1

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
 f"{n}\n",
 "(\n"
 ]
 for j in range(n):
       num_string = format(eps[j], ".8f").rstrip('0').rstrip('.')
       data.append(f"{num_string}\n")
 data.extend([
 ")\n",
 ";\n",
 "\n",
 "boundaryField\n",
 "{\n"
 "    inlet\n"
 "    {\n",
 "        type            zeroGradient;\n",
 "    }\n"
 "    outlet\n"
 "    {\n",
 "        type            zeroGradient;\n",
 "    }\n"
 "    \"wall_.*\"""\n"
 "    {\n",
 "        type            zeroGradient;\n",
 "    }\n"
 ])

 if (dimension=="2D"):
     data.extend([
     "    frontAndBack\n"
     "    {\n",
     "        type            empty;\n",
     "    }\n"
     ]) 


 ##proc k to k-1
 if (ipz>0):
   data.extend([
   f"    procBoundary{rank}to{rank - NPX * NPY}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           nonuniform List<scalar> \n",
   f"{nbz0}\n",
   "(\n"
   ])
   for j in range(nbz0):
       num_string = format(eps_bz0[j], ".8f").rstrip('0').rstrip('.')
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
   f"{nby0}\n",
   "(\n"
   ])
   for j in range(nby0):
       num_string = format(eps_by0[j], ".8f").rstrip('0').rstrip('.')
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
   f"{nbx0}\n",
   "(\n"
   ])
   for j in range(nbx0):
       num_string = format(eps_bx0[j], ".8f").rstrip('0').rstrip('.')
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
   f"{nbx1}\n",
   "(\n"
   ])
   for j in range(nbx1):
       num_string = format(eps_bx1[j], ".8f").rstrip('0').rstrip('.')
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
   f"{nby1}\n",
   "(\n"
   ])
   for j in range(nby1):
       num_string = format(eps_by1[j], ".8f").rstrip('0').rstrip('.')
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
   f"{nbz1}\n",
   "(\n"
   ])
   for j in range(nbz1):
       num_string = format(eps_bz1[j], ".8f").rstrip('0').rstrip('.')
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
