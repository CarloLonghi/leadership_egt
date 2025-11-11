import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import math

def gaussian(x, mu, sig):
    return np.exp(-((x - mu) ** 2) / (2 * (sig ** 2))) / (np.sqrt(2 * np.pi) * sig)

nr = 3
nc = 3
fntsize=15

k=2
pSv=np.linspace(0.,1.,num=50)
deltaLv=[0, 1, 2, 4, 8]
f=0
betaF=1.
N = 9
eps = 0.01
eps1 = 1 - eps
rv=np.array([3*k,6*k,8*k])


fig,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all', figsize=(10,9))
fig.subplots_adjust(hspace=0.4, wspace=0.2)
nticksY=6
nticksX=3

cmap = plt.get_cmap('viridis')

file = './newtests/1bit/leadership/multileader/res_k2'
data = np.load(file + '.npy')

for idr, r in enumerate(rv):
    ax=axs[idr,0]

    for iddl, deltaL in enumerate(deltaLv):
        deltaF = deltaL
        pleadS = 1/(1+np.exp(-deltaL))
        pleadW=1-pleadS
        res = np.zeros((pSv.shape[0]))
        for strat in range(2):
            s1 = strat
            for idps, pS in enumerate(pSv):
                pW = 1 - pS
                pF = np.zeros((2,2))

                pF[0,0] = 1/(1+np.exp(-betaF*(f)))
                pF[1,1] = 1/(1+np.exp(-betaF*(f)))
                pF[0,1] = 1/(1+np.exp(-betaF*(f+deltaF)))
                pF[1,0] = 1/(1+np.exp(-betaF*(f-deltaF)))
                
                benefit = 0
                cost = 0
                for Ns in range(N+1):
                    ns1 = Ns
                    Nw = N - Ns
                    nw1 = Nw
                    pNs = math.factorial(N)/(math.factorial(Ns)*math.factorial(Nw)) * pS**Ns * pW**Nw
                    for ns1l in range(ns1+1):
                        probs1l = math.factorial(ns1)/(math.factorial(ns1l)*math.factorial(ns1-ns1l)) * pleadS**ns1l * (1-pleadS)**(ns1-ns1l)
                        for nw1l in range(nw1+1):
                            probw1l = math.factorial(nw1)/(math.factorial(nw1l)*math.factorial(nw1-nw1l)) * pleadW**nw1l * (1-pleadW)**(nw1-nw1l)
                            
                            probl = probs1l*probw1l

                            Nwc = (nw1-nw1l)*s1
                            Nsc = (ns1-ns1l)*s1
                            Nwd = (nw1-nw1l)*(1-s1)
                            Nsd = (ns1-ns1l)*(1-s1)

                            Nwcl = nw1l*s1
                            Nscl = ns1l*s1
                            Nwdl = nw1l*(1-s1)
                            Nsdl = ns1l*(1-s1)
                            Nwl = Nwcl + Nwdl
                            Nsl = Nscl + Nsdl

                            Nwnl = Nw - Nwl
                            Nsnl = Ns - Nsl

                            if nw1l + ns1l == 0:
                                benefit += pNs*probl*(
                                    (Nwc+Nsc)*eps1 + (Nwd+Nsd)*eps
                                )
                            else:
                                follow_s = (pleadS * Nsl) / (pleadS * Nsl + pleadW * Nwl)
                                follow_w = (pleadW * Nwl) / (pleadS * Nsl + pleadW * Nwl)
                                benefit += pNs*probl*(
                                    (Nwcl + Nscl) * eps1 + (Nwdl + Nsdl) * eps # leaders
                                )                         
                                if Nwl > 0:
                                    benefit += pNs*probl*(
                                        follow_w*(
                                            (1-pF[0,0])*(Nwc*eps1 + Nwd*eps)+
                                            (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                            pF[0,0]*Nwnl*((Nwcl/Nwl)*(eps1**2+eps**2) + (Nwdl/Nwl)*(2*eps1*eps))+
                                            pF[1,0]*Nsnl*((Nwcl/Nwl)*(eps1**2+eps**2) + (Nwdl/Nwl)*(2*eps1*eps))
                                        )
                                    )

                                if Nsl > 0:
                                    benefit += pNs*probl*(
                                        follow_s*(
                                            (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                            (1-pF[1,1])*(Nsc*eps1 + Nsd*eps)+
                                            pF[0,1]*Nwnl*((Nscl/Nsl)*(eps1**2+eps**2) + (Nsdl/Nsl)*(2*eps1*eps))+
                                            pF[1,1]*Nsnl*((Nscl/Nsl)*(eps1**2+eps**2) + (Nsdl/Nsl)*(2*eps1*eps))                                                
                                        )
                                    )

                benefit /= N

                res[idps] += benefit * data[r-1, iddl, idps, strat]
            

        ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
        ax.set_xticklabels(np.linspace(pSv[0],pSv[-1],nticksX), fontsize=12)
        ax.set_yticks(np.linspace(0, 1, 3))
        ax.set_yticklabels(np.linspace(0,1,3), fontsize=12)
        ax.set_ylim(0.0, 1.0)
        ax.plot(res, label='$\Delta_f=\Delta_f=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

        ax.set_ylabel(r'cooperation level', fontsize=fntsize)
        ax.text(-40,0.46,"$r$=%d" % rv[idr], size=13)
        if idr==nr-1: ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idr==0: ax.set_title("1 bit" % rv[idr], size=13)

file = './newtests/2bits/strengthstrat/multileader/res_2bits_multileader'
data = np.load(file + '.npy')

for idr, r in enumerate(rv):
    ax=axs[idr,1]

    for iddl, deltaL in enumerate(deltaLv):
        deltaF = deltaL
        pleadS = 1/(1+np.exp(-deltaL))
        pleadW=1-pleadS
        res = np.zeros((pSv.shape[0]))
        for strat in range(4):
            s1 = [strat%2, strat//2]
            for idps, pS in enumerate(pSv):
                pW = 1 - pS
                pF = np.zeros((2,2))

                pF[0,0] = 1/(1+np.exp(-betaF*(f)))
                pF[1,1] = 1/(1+np.exp(-betaF*(f)))
                pF[0,1] = 1/(1+np.exp(-betaF*(f+deltaF)))
                pF[1,0] = 1/(1+np.exp(-betaF*(f-deltaF)))
                
                benefit = 0
                cost = 0
                for Ns in range(N+1):
                    ns1 = Ns
                    Nw = N - Ns
                    nw1 = Nw
                    pNs = math.factorial(N)/(math.factorial(Ns)*math.factorial(Nw)) * pS**Ns * pW**Nw
                    for ns1l in range(ns1+1):
                        probs1l = math.factorial(ns1)/(math.factorial(ns1l)*math.factorial(ns1-ns1l)) * pleadS**ns1l * (1-pleadS)**(ns1-ns1l)
                        for nw1l in range(nw1+1):
                            probw1l = math.factorial(nw1)/(math.factorial(nw1l)*math.factorial(nw1-nw1l)) * pleadW**nw1l * (1-pleadW)**(nw1-nw1l)
                            
                            probl = probs1l*probw1l

                            Nwc = (nw1-nw1l)*s1[0]
                            Nsc = (ns1-ns1l)*s1[1]
                            Nwd = (nw1-nw1l)*(1-s1[0])
                            Nsd = (ns1-ns1l)*(1-s1[1])

                            Nwcl = nw1l*s1[0]
                            Nscl = ns1l*s1[1]
                            Nwdl = nw1l*(1-s1[0])
                            Nsdl = ns1l*(1-s1[1])
                            Nwl = Nwcl + Nwdl
                            Nsl = Nscl + Nsdl

                            Nwnl = Nw - Nwl
                            Nsnl = Ns - Nsl

                            if nw1l + ns1l == 0:
                                benefit += pNs*probl*(
                                    (Nwc+Nsc)*eps1 + (Nwd+Nsd)*eps
                                )
                            else:
                                follow_s = (pleadS * Nsl) / (pleadS * Nsl + pleadW * Nwl)
                                follow_w = (pleadW * Nwl) / (pleadS * Nsl + pleadW * Nwl)
                                benefit += pNs*probl*(
                                    (Nwcl + Nscl) * eps1 + (Nwdl + Nsdl) * eps # leaders
                                )                         
                                if Nwl > 0:
                                    benefit += pNs*probl*(
                                        follow_w*(
                                            (1-pF[0,0])*(Nwc*eps1 + Nwd*eps)+
                                            (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                            pF[0,0]*Nwnl*((Nwcl/Nwl)*(eps1**2+eps**2) + (Nwdl/Nwl)*(2*eps1*eps))+
                                            pF[1,0]*Nsnl*((Nwcl/Nwl)*(eps1**2+eps**2) + (Nwdl/Nwl)*(2*eps1*eps))
                                        )
                                    )

                                if Nsl > 0:
                                    benefit += pNs*probl*(
                                        follow_s*(
                                            (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                            (1-pF[1,1])*(Nsc*eps1 + Nsd*eps)+
                                            pF[0,1]*Nwnl*((Nscl/Nsl)*(eps1**2+eps**2) + (Nsdl/Nsl)*(2*eps1*eps))+
                                            pF[1,1]*Nsnl*((Nscl/Nsl)*(eps1**2+eps**2) + (Nsdl/Nsl)*(2*eps1*eps))                                                
                                        )
                                    )

                benefit /= N

                res[idps] += benefit * data[r-1, iddl, idps, strat]
            

        ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
        ax.set_xticklabels(np.linspace(pSv[0],pSv[-1],nticksX), fontsize=12)
        ax.set_yticks(np.linspace(0, 1, 3))
        ax.set_yticklabels(np.linspace(0,1,3), fontsize=12)
        ax.set_ylim(0.0, 1.0)
        ax.plot(res, label='$\Delta_f=\Delta_f=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

        if idr==nr-1: ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idr==0: ax.set_title("2 bits strength" % rv[idr], size=13)

file = './newtests/2bits/leadstrat/multileader/res_new'
data = np.load(file + '.npy')

for idr, r in enumerate(rv):
    ax=axs[idr,2]

    for iddl, deltaL in enumerate(deltaLv):
        deltaF = deltaL
        pleadS = 1/(1+np.exp(-deltaL))
        pleadW=1-pleadS
        res = np.zeros((pSv.shape[0]))
        for strat in range(4):
            s1 = [strat%2, strat//2]
            for idps, pS in enumerate(pSv):
                pW = 1 - pS
                pF = np.zeros((2,2))

                pF[0,0] = 1/(1+np.exp(-betaF*(f)))
                pF[1,1] = 1/(1+np.exp(-betaF*(f)))
                pF[0,1] = 1/(1+np.exp(-betaF*(f+deltaF)))
                pF[1,0] = 1/(1+np.exp(-betaF*(f-deltaF)))
                
                benefit = 0
                cost = 0
                for Ns in range(N+1):
                    ns1 = Ns
                    Nw = N - Ns
                    nw1 = Nw
                    pNs = math.factorial(N)/(math.factorial(Ns)*math.factorial(Nw)) * pS**Ns * pW**Nw
                    for ns1l in range(ns1+1):
                        probs1l = math.factorial(ns1)/(math.factorial(ns1l)*math.factorial(ns1-ns1l)) * pleadS**ns1l * (1-pleadS)**(ns1-ns1l)
                        for nw1l in range(nw1+1):
                            probw1l = math.factorial(nw1)/(math.factorial(nw1l)*math.factorial(nw1-nw1l)) * pleadW**nw1l * (1-pleadW)**(nw1-nw1l)
                            
                            probl = probs1l*probw1l

                            Nwc = (nw1-nw1l)*s1[0]
                            Nsc = (ns1-ns1l)*s1[0]
                            Nwd = (nw1-nw1l)*(1-s1[0])
                            Nsd = (ns1-ns1l)*(1-s1[0])

                            Nwcl = nw1l*s1[1]
                            Nscl = ns1l*s1[1]
                            Nwdl = nw1l*(1-s1[1])
                            Nsdl = ns1l*(1-s1[1])
                            Nwl = Nwcl + Nwdl
                            Nsl = Nscl + Nsdl

                            Nwnl = Nw - Nwl
                            Nsnl = Ns - Nsl

                            if nw1l + ns1l == 0:
                                benefit += pNs*probl*(
                                    (Nwc+Nsc)*eps1 + (Nwd+Nsd)*eps
                                )
                            else:
                                follow_s = (pleadS * Nsl) / (pleadS * Nsl + pleadW * Nwl)
                                follow_w = (pleadW * Nwl) / (pleadS * Nsl + pleadW * Nwl)
                                benefit += pNs*probl*(
                                    (Nwcl + Nscl) * eps1 + (Nwdl + Nsdl) * eps # leaders
                                )                         
                                if Nwl > 0:
                                    benefit += pNs*probl*(
                                        follow_w*(
                                            (1-pF[0,0])*(Nwc*eps1 + Nwd*eps)+
                                            (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                            pF[0,0]*Nwnl*((Nwcl/Nwl)*(eps1**2+eps**2) + (Nwdl/Nwl)*(2*eps1*eps))+
                                            pF[1,0]*Nsnl*((Nwcl/Nwl)*(eps1**2+eps**2) + (Nwdl/Nwl)*(2*eps1*eps))
                                        )
                                    )

                                if Nsl > 0:
                                    benefit += pNs*probl*(
                                        follow_s*(
                                            (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                            (1-pF[1,1])*(Nsc*eps1 + Nsd*eps)+
                                            pF[0,1]*Nwnl*((Nscl/Nsl)*(eps1**2+eps**2) + (Nsdl/Nsl)*(2*eps1*eps))+
                                            pF[1,1]*Nsnl*((Nscl/Nsl)*(eps1**2+eps**2) + (Nsdl/Nsl)*(2*eps1*eps))                                                
                                        )
                                    )

                benefit /= N

                res[idps] += benefit * data[r-1, iddl, idps, strat]
            

        ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
        ax.set_xticklabels(np.linspace(pSv[0],pSv[-1],nticksX), fontsize=12)
        ax.set_yticks(np.linspace(0, 1, 3))
        ax.set_yticklabels(np.linspace(0,1,3), fontsize=12)
        ax.set_ylim(0.0, 1.0)
        ax.plot(res, label='$\Delta_f=\Delta_f=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))
        if idr==nr-1: ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idr==0: ax.set_title("2 bits leader" % rv[idr], size=13)


file = './newtests/4bits/multileader/res_new'
data = np.load(file + '.npy')

for idr, r in enumerate(rv):
    ax=axs[idr,3]

    for iddl, deltaL in enumerate(deltaLv):
        deltaF = deltaL
        pleadS = 1/(1+np.exp(-deltaL))
        pleadW=1-pleadS
        res = np.zeros((pSv.shape[0]))
        for strat in range(16):
            s1=[strat%8%4%2//1, strat%8%4//2, strat%8//4, strat//8]
            for idps, pS in enumerate(pSv):
                pW = 1 - pS
                pF = np.zeros((2,2))

                pF[0,0] = 1/(1+np.exp(-betaF*(f)))
                pF[1,1] = 1/(1+np.exp(-betaF*(f)))
                pF[0,1] = 1/(1+np.exp(-betaF*(f+deltaF)))
                pF[1,0] = 1/(1+np.exp(-betaF*(f-deltaF)))
                
                benefit = 0
                cost = 0
                for Ns in range(N+1):
                    ns1 = Ns
                    Nw = N - Ns
                    nw1 = Nw
                    pNs = math.factorial(N)/(math.factorial(Ns)*math.factorial(Nw)) * pS**Ns * pW**Nw
                    for ns1l in range(ns1+1):
                        probs1l = math.factorial(ns1)/(math.factorial(ns1l)*math.factorial(ns1-ns1l)) * pleadS**ns1l * (1-pleadS)**(ns1-ns1l)
                        for nw1l in range(nw1+1):
                            probw1l = math.factorial(nw1)/(math.factorial(nw1l)*math.factorial(nw1-nw1l)) * pleadW**nw1l * (1-pleadW)**(nw1-nw1l)
                            
                            probl = probs1l*probw1l

                            Nwc = (nw1-nw1l)*s1[0]
                            Nsc = (ns1-ns1l)*s1[1]
                            Nwd = (nw1-nw1l)*(1-s1[0])
                            Nsd = (ns1-ns1l)*(1-s1[1])

                            Nwcl = nw1l*s1[2]
                            Nscl = ns1l*s1[3]
                            Nwdl = nw1l*(1-s1[2])
                            Nsdl = ns1l*(1-s1[3])
                            Nwl = Nwcl + Nwdl
                            Nsl = Nscl + Nsdl

                            Nwnl = Nw - Nwl
                            Nsnl = Ns - Nsl

                            if nw1l + ns1l == 0:
                                benefit += pNs*probl*(
                                    (Nwc+Nsc)*eps1 + (Nwd+Nsd)*eps
                                )
                            else:
                                follow_s = (pleadS * Nsl) / (pleadS * Nsl + pleadW * Nwl)
                                follow_w = (pleadW * Nwl) / (pleadS * Nsl + pleadW * Nwl)      
                                benefit += pNs*probl*(
                                    (Nwcl + Nscl) * eps1 + (Nwdl + Nsdl) * eps # leaders
                                )                      
                                if Nwl > 0:
                                    benefit += pNs*probl*(
                                        follow_w*(
                                            (1-pF[0,0])*(Nwc*eps1 + Nwd*eps)+
                                            (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                            pF[0,0]*Nwnl*((Nwcl/Nwl)*(eps1**2+eps**2) + (Nwdl/Nwl)*(2*eps1*eps))+
                                            pF[1,0]*Nsnl*((Nwcl/Nwl)*(eps1**2+eps**2) + (Nwdl/Nwl)*(2*eps1*eps))
                                        )
                                    )

                                if Nsl > 0:
                                    benefit += pNs*probl*(
                                        follow_s*(
                                            (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                            (1-pF[1,1])*(Nsc*eps1 + Nsd*eps)+
                                            pF[0,1]*Nwnl*((Nscl/Nsl)*(eps1**2+eps**2) + (Nsdl/Nsl)*(2*eps1*eps))+
                                            pF[1,1]*Nsnl*((Nscl/Nsl)*(eps1**2+eps**2) + (Nsdl/Nsl)*(2*eps1*eps))                                                
                                        )
                                    )
                benefit /= N

                res[idps] += benefit * data[r-1, iddl, idps, strat]
            

        ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
        ax.set_xticklabels(np.linspace(pSv[0],pSv[-1],nticksX), fontsize=12)
        ax.set_yticks(np.linspace(0, 1, 3))
        ax.set_yticklabels(np.linspace(0,1,3), fontsize=12)
        ax.set_ylim(0.0, 1.0)
        ax.plot(res, label='$\Delta_f=\Delta_f=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))
        if idr==nr-1: ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idr==0: ax.set_title("4 bits" % rv[idr], size=13)


legend_elements = [Line2D([], [], marker='None', label='$\Delta_l=\Delta_f$', linestyle='None')]
legend_elements += [Line2D([], [], marker='s', color=cmap((idx)/(len(deltaLv))), label='%d'%deltaLv[idx],
                          markerfacecolor=cmap((idx)/(len(deltaLv))), markersize=10, linestyle='None') for idx in range(len(deltaLv))]
plt.legend( loc='upper center', bbox_to_anchor=(-1.3, -0.4),
          fancybox=True, shadow=False, ncol=6, columnspacing=0.0, handles=legend_elements,handletextpad=-0.3,fontsize=13)
plt.savefig('./newtests/summary_cl_plots.png', bbox_inches='tight', dpi=300)

# plt.show()