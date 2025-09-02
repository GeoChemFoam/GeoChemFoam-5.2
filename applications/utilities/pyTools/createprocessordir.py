import os
import time

def main(rank):
    #start_time = time.time()
    
    
    base_path = f'processor{rank}'

    dirs = [f"{base_path}/0", f"{base_path}/constant/polyMesh"]

    for d in dirs:
      os.makedirs(d, exist_ok=True)

    #elapsed_time = time.time() - start_time
    #if rank == 0:
        #print(f"Elapsed time and IO time create processor: {rank} {elapsed_time:.6f} 0.0 seconds")
 
 
 
