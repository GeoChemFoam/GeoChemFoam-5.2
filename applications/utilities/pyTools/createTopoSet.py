import numpy as np
import array
import os
import time
import sys

segmentation = sys.argv[1]
eps_min=float(sys.argv[2])
micro_por=sys.argv[3]
refineStokes=int(sys.argv[4])

f=open("system/topoSetDict",'w')
f.seek(0) #get to the first position
f.write("FoamFile"+'\n')
f.write("{"+'\n')
f.write("    version     2.0;"+'\n')
f.write("    format      ascii;"+'\n')
f.write("    class       dictionary;"+'\n')
f.write("    object      topoSetDict;"+'\n')
f.write("}"+'\n')
f.write(""+'\n')
f.write("actions"+'\n')
f.write("("+'\n')


if (segmentation=='phases'):
  micro_por_array = micro_por.split(',')
  nlist=len(micro_por_array);
  poro=np.float64(sorted(micro_por_array,key=float))

  f.write("{"+'\n')
  f.write("name refinementRegion;"+'\n')
  f.write("type cellSet;"+'\n')
  f.write("action new;"+'\n')
  f.write("source fieldToCell;"+'\n')
  f.write("field eps;"+'\n')
  f.write("min "+str(poro[nlist-1]+1e-3-3e-3*refineStokes)+";"+'\n')
  f.write("max "+str(poro[nlist-1]+3e-3)+";"+'\n')
  f.write("}"'\n')

  for j in range (0,nlist-1):
    f.write("{"+'\n')
    f.write("name refinementRegion;"+'\n')
    f.write("type cellSet;"+'\n')
    f.write("action add;"+'\n')
    f.write("source fieldToCell;"+'\n')
    f.write("field eps;"+'\n')
    f.write("min "+str(poro[j]+(poro[j+1]-poro[j])*5e-3)+";"+'\n')
    f.write("max "+str(poro[j+1]-(poro[j+1]-poro[j])*5e-3)+";"+'\n')
    f.write("}"+'\n')
  f.write(");")
  f.close()

else:
  f.write("{"+'\n')
  f.write("name refinementRegion;"+'\n')
  f.write("type cellSet;"+'\n')
  f.write("action new;"+'\n')
  f.write("source fieldToCell;"+'\n')
  f.write("field eps;"+'\n')
  f.write("min "+str(eps_min+(1-eps_min)*5e-3)+";"+'\n')
  f.write("max "+str(1.0-(1-eps_min)*5e-3 +refineStokes*1e-2)+";"+'\n')
  f.write("}"+'\n')
  f.write(");")
  f.close()

