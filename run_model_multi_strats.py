import numpy as np
import mechanics.evoEGT as evo
from mechanics.heterogeneous16 import calcH, calcWCD

def coop_pF_r(rv,M,N,HZ,beta,eps,pSv,f,betaF,deltav):
# Input: pFv, rv, Mv (vectors with values of pF, r, and M), N, HZ (H or Z), beta, eps
# Output: matrix with the fraction of cooperators as a function of pF and r
    if np.isscalar(HZ):
        H=calcH(N-1,HZ-1)

    MAT = np.zeros((len(deltav), len(pSv), len(rv), 16, 16))

    for iddl, deltaL in enumerate(deltav):
        deltaF = deltaL
        pF=np.zeros((2,2))
        pF[0,0] = 1/(1+np.exp(-betaF*(f)))
        pF[1,1] = 1/(1+np.exp(-betaF*(f)))
        pF[0,1] = 1/(1+np.exp(-betaF*(f+deltaF)))
        pF[1,0] = 1/(1+np.exp(-betaF*(f-deltaF)))
        for idps, pS in enumerate(pSv):
            WCD=calcWCD(N,eps,pF,deltaL,pS,M)
            print(deltaL, pS)
            for idr, r in enumerate(rv):
                SD,fixM = evo.Wgroup2SD(WCD,H,[r,-1.],beta,infocheck=False)
                MAT[iddl, idps, idr] = fixM
    return MAT



if __name__ == "__main__":

    import time

    t0=time.time()

    eps=0.01 #0.01
    Z=100
    N=9
    beta=1.
    M=0
    betaF=1

    f=0

    deltaLv=np.array([1.5, 4])
    pSv=np.array([0.1, 0.5, 0.7])
    rv=np.array([6,])
    
    labfilenpy='./newtests/4bits/multileader/fm_ip'
    MAT=coop_pF_r(rv,M,N,Z,beta,eps,pSv,f,betaF,deltaLv)
    np.save(labfilenpy,MAT)             # save matrix for heatmap
    print('data saved to file!')
    