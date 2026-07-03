import numpy as np
import mechanics_singleleader.evoEGT as evo
from mechanics_singleleader.heterogeneous4 import calcH, calcWCD

def coop_pF_r(rv,M,N,HZ,beta,eps,pSv,deltaLv,f,betaF):
# Input: pFv, rv, Mv (vectors with values of pF, r, and M), N, HZ (H or Z), beta, eps
# Output: matrix with the fraction of cooperators as a function of pF and r
    if np.isscalar(HZ):
        H=calcH(N-1,HZ-1)

    MAT = np.zeros((len(rv), len(deltaLv), len(pSv), 4))

    for iddl, deltaL in enumerate(deltaLv):
        deltaF=deltaL
        for idps, pS in enumerate(pSv):
            pF=np.zeros((2,2))
            pF[0,0] = 1/(1+np.exp(-betaF*(f)))
            pF[1,1] = 1/(1+np.exp(-betaF*(f)))
            pF[0,1] = 1/(1+np.exp(-betaF*(f+deltaF)))
            pF[1,0] = 1/(1+np.exp(-betaF*(f-deltaF)))

            WCD=calcWCD(N,eps,pF,deltaL,pS,M)
            print(deltaL, pS)
            for idr, r, in enumerate(rv):
                SD,fixM = evo.Wgroup2SD(WCD,H,[r,-1.],beta,infocheck=False)
                MAT[idr, iddl, idps] = SD[:,0]
    return MAT


if __name__ == "__main__":

    import time

    t0=time.time()

    eps=0.01 #0.01
    Z=100
    N=9
    beta=1.
    M=0
    f=0
    betaF=1

    deltaLv=[0, 1, 2, 4, 8]
    pSv=np.linspace(0,1.,num=50)
    rv=np.linspace(1,10,num=10)
    
    labfilenpy='newtests/res_2bits_singleleader'
    MAT=coop_pF_r(rv,M,N,Z,beta,eps,pSv,deltaLv,f,betaF)
    np.save(labfilenpy,MAT)             # save matrix for heatmap
    print('data saved to file!')
    