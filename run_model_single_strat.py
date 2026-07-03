import numpy as np
import mechanics_singleleader.evoEGT as evo
from mechanics_singleleader.heterogeneous4_leaderstrat import calcH, calcWCD

def coop_pF_r(rv,M,N,HZ,beta,eps,pSv,f,betaF,deltav):
# Input: pFv, rv, Mv (vectors with values of pF, r, and M), N, HZ (H or Z), beta, eps
# Output: matrix with the fraction of cooperators as a function of pF and r
    if np.isscalar(HZ):
        H=calcH(N-1,HZ-1)

    MAT = np.zeros((len(deltav), len(pSv), len(rv), 4))

    for iddf, deltaF in enumerate(deltav):
        deltaL = deltaF
        pF=np.zeros((2,2))
        pF[0,0] = 1/(1+np.exp(-betaF*(f)))
        pF[1,1] = 1/(1+np.exp(-betaF*(f)))
        pF[0,1] = 1/(1+np.exp(-betaF*(f+deltaF)))
        pF[1,0] = 1/(1+np.exp(-betaF*(f-deltaF)))
        for idps, pS in enumerate(pSv):
            WCD=calcWCD(N,eps,pF,deltaL,pS,M)
            print(deltaF,pS)
            for idr, r in enumerate(rv):
                SD,fixM = evo.Wgroup2SD(WCD,H,[r,-1.],beta,infocheck=False)
                best_s = np.argmax(SD)
                if SD[best_s] >= 0.5:
                    MAT[iddf, idps, idr, best_s] = SD[best_s]
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

    deltaLv=np.linspace(0,8,num=50)
    pSv=np.linspace(0.,1.,num=50)
    rv=np.linspace(1, 10, 10)
    
    labfilenpy='newtests/2bits/leadstrat/singleleader/res_2bits_singleleader_deltaps'
    MAT=coop_pF_r(rv,M,N,Z,beta,eps,pSv,f,betaF,deltaLv)
    np.save(labfilenpy,MAT)             # save matrix for heatmap
    print('data saved to file!')
    