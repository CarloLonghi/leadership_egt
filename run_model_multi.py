import numpy as np
import mechanics.evoEGT as evo
from mechanics.heterogeneous4 import calcH, calcWCD
import os
import concurrent
import multiprocessing as mp

from ctypes import c_double

class SharedArrayManager:
    def __init__(self, n):
        self.n = n
        self.shared_array = mp.Array(c_double, len(rv)*len(deltaLv)*len(pSv)*4)
        self.shape = (MAT.shape)
    
    def get_array(self):
        return np.frombuffer(self.shared_array.get_obj()).reshape(self.shape)

def worker_function(iddl):
    """
    Worker function that performs calculations
    """
    result_array = shared_array_manager.get_array()

    deltaL = deltaLv[iddl]
    deltaF = deltaL
    for idps in range(0, pSv.shape[0]):
        pS = pSv[idps]

        pF=np.zeros((2,2))
        pF[0,0] = 1/(1+np.exp(-betaF*(f)))
        pF[1,1] = 1/(1+np.exp(-betaF*(f)))
        pF[0,1] = 1/(1+np.exp(-betaF*(f+deltaF)))
        pF[1,0] = 1/(1+np.exp(-betaF*(f-deltaF)))   

        WCD=calcWCD(N,eps,pF,deltaL,pS,M)
        print(deltaL, pS)
        
        for idr, r, in enumerate(rv):
            SD,fixM = evo.Wgroup2SD(WCD,H,[r,-1.],beta,infocheck=False)
            
            result_array[idr, iddl, idps] = SD[:,0]            

def coop_pF_r():
    
    # Use ProcessPoolExecutor for parallel execution
    with concurrent.futures.ProcessPoolExecutor(max_workers=mp.cpu_count()) as executor:
        # Submit all tasks
        futures = [executor.submit(worker_function, i) for i in range(0, deltaLv.shape[0])]
        
        # Wait for all tasks to complete
        for future in futures:
            future.result()
    
    # Get final results
    final_result = shared_array_manager.get_array()
    return final_result



if __name__ == "__main__":

    import time

    t0=time.time()

    eps=0.01 #0.01
    Z=100
    k=1
    N=9*k
    beta=1.
    M=0
    f=0
    betaF=1

    H=calcH(N-1,Z-1)

    deltaLv=np.array([0, 1, 2, 4, 8])
    pSv=np.linspace(0,1.,num=50)
    rv=[3*k, 6*k, 8*k]
    
    MAT = np.zeros((len(rv), len(deltaLv), len(pSv), 4))
    # Create shared array manager
    n = 5
    shared_array_manager = SharedArrayManager(n)
    # Initialize array with zeros
    result_array = shared_array_manager.get_array()
    result_array[:] = 0    

    folder = './newtests/2bits/strengthstrat/multileader/'
    labfilenpy= folder + 'res_k4'
    backup = folder + 'backup_k4'
    MAT = coop_pF_r()
    np.save(labfilenpy,MAT)       
    print('data saved to file!')