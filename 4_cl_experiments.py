import numpy as np
import evoEGT as evo
from heterogeneous4 import calcH, calcWCD
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


def plotCOOPheat(MAT,f1v,f2v,sv,label):
# Input: MAT (matrix from "coop_pF_r" function), pFv, rv ,Mv (vectors with values of pF, r, and M), label (name for the output file)
# Output: heatmap plot of the fraction of cooperators as a function of pF and r, for different M
    import matplotlib.pyplot as plt
    fntsize=12
    nr=5
    nc=2
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)
    k=-1
    for idx in range(len(sv)):
        i = idx // nc
        j = idx % nc

        ax=axs[i,j]
        k=k+1
        h=ax.imshow(MAT[:,:,k],origin='lower', interpolation='none',aspect='auto',vmin=0,vmax=4)
        nticksY=5
        nticksX=3
        ax.set_xticks(np.linspace(0, MAT.shape[1]-1, nticksX))
        ax.set_yticks(np.linspace(0, MAT.shape[0]-1, nticksY))
        ax.set_xticklabels(np.linspace(f2v[0],f2v[-1],nticksX))
        ax.set_yticklabels(np.linspace(f2v[0],f2v[-1],nticksY))
        ax.text(25,50,"$S=%.2f$" % sv[k], size=10 )
        if i==nr-1: ax.set_xlabel(r'$p_{F_{ww,ss}}$', fontsize=fntsize)
        if j==0: ax.set_ylabel(r'$p_{F_{ws}}$', fontsize=fntsize)
    
    f.subplots_adjust(right=0.8)
    cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
    cb = f.colorbar(h, cax=cbar_ax)
    cb.set_ticks([1,2,3,4])
    cb.set_ticklabels(['ALLD','WCSD','WDSC','ALLC'])

#cb=f.colorbar(h, fraction=0.1,format='%.2f')
    #cb.set_label(label=r'$f_C$')
    f.savefig(label+'.png',bbox_inches='tight',dpi=300)
    f.clf()     
    return

def plotsingleheat(MAT,fv,rv,label):
# Input: MAT (matrix from "coop_pF_r" function), pFv, rv ,Mv (vectors with values of pF, r, and M), label (name for the output file)
# Output: heatmap plot of the fraction of cooperators as a function of pF and r, for different M
    import matplotlib.pyplot as plt
    fntsize=12
    f,ax=plt.subplots()
    h=ax.imshow(MAT,origin='lower', interpolation='none',aspect='auto')
    nticksY=5
    nticksX=3
    ax.set_xticks(np.linspace(0, MAT.shape[1]-1, nticksX))
    ax.set_yticks(np.linspace(0, MAT.shape[0]-1, nticksY))
    ax.set_xticklabels(np.linspace(fv[0],fv[-1],nticksX))
    ax.set_yticklabels(np.linspace(rv[0],rv[-1],nticksY))
    ax.set_xlabel(r'$f$', fontsize=fntsize)
    ax.set_ylabel(r'$r$', fontsize=fntsize)
#cb=f.colorbar(h, fraction=0.1,format='%.2f')
    #cb.set_label(label=r'$f_C$')
    f.savefig(label+'.png',bbox_inches='tight',dpi=300)
    f.clf()     
    return


if __name__ == "__main__":

    import time

    t0=time.time()

####### Plot heatmap #########################################
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
    np.save(labfilenpy,MAT)             # save matrix for heatmap
    print('data saved to file!')
    
    # MAT=np.load(labfilenpy+'.npy')      # load matrix for heatmap 
    # plotCOOPheat(MAT,f1v,f2v,sv,labfilenpy)      # plot heatmap
    #plotsingleheat(MAT,fv,rv,labfilenpy)
#####################################################