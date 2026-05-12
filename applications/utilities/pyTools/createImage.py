import sys
import os
from mpi4py import MPI
import createImgGrayscale as createImgGrayscale
import createImgGrayscaleCyclic as createImgGrayscaleCyclic
import createImgPhases as createImgPhases
import createImgPhasesCyclic as createImgPhasesCyclic
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

x_dim = int(sys.argv[1])
y_dim = int(sys.argv[2])
z_dim = int(sys.argv[3])
x_min = int(sys.argv[4])
x_max = int(sys.argv[5])
y_min = int(sys.argv[6])
y_max = int(sys.argv[7])
z_min = int(sys.argv[8])
z_max = int(sys.argv[9])
n_x = int(sys.argv[10])
n_y = int(sys.argv[11])
n_z = int(sys.argv[12])
res = float(sys.argv[13])
Image_name = sys.argv[14]
padWidth = int(sys.argv[15])
pad_value = float(sys.argv[16])
value1 = float(sys.argv[17])
value0 = float(sys.argv[18])
field_min = float(sys.argv[19])
dimension = sys.argv[20]
direction = int(sys.argv[21])
cyclic = sys.argv[22]
segmentation = sys.argv[23]
phases_value=sys.argv[24]
phases=sys.argv[25]
NPX = int(sys.argv[26])
NPY = int(sys.argv[27])
NPZ = int(sys.argv[28])
fieldName=sys.argv[29]
#start_time = time.time()
NP=NPX*NPY*NPZ


if NP>1:

    output_path = "processor"+str(rank)+"/"
else:

    output_path = ""


if rank == 0:
    print(f"create {fieldName}")

if (segmentation=='grayscale'):
  if (cyclic=='yes'):
    createImgGrayscaleCyclic.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, Image_name, padWidth, pad_value,value1, value0, field_min, dimension,direction, NPX, NPY, NPZ, rank, output_path,fieldName)
  else:
    createImgGrayscale.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, Image_name, padWidth, pad_value,value1, value0, field_min, dimension,direction, NPX, NPY, NPZ, rank, output_path,fieldName)
elif (segmentation=='phases'):
  if (cyclic=='yes'):
    createImgPhasesCyclic.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, Image_name, padWidth, pad_value,value1, value0, phases_value, phases, dimension,direction, NPX, NPY, NPZ, rank, output_path,fieldName)
  else:
    createImgPhases.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, Image_name, padWidth, pad_value,value1, value0, phases_value, phases, dimension,direction, NPX, NPY, NPZ, rank, output_path, fieldName)
else:
  raise TypeError("only grayscale and phases segmentation accepted")


#end_time = time.time()
#elapsed_time = end_time - start_time

#if rank ==0:
    #print("Elapsed time createMesh.py:", elapsed_time, "seconds")

MPI.Finalize()

