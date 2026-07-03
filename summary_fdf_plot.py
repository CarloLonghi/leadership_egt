import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import math

def gaussian(x, mu, sig):
    return np.exp(-((x - mu) ** 2) / (2 * (sig ** 2))) / (np.sqrt(2 * np.pi) * sig)

def aeps(pact,eps):
# Input: pact probability to perform the desired action (without error), eps probability of comitting an error 
# Output: action actually performed 
    return pact*(1.-2.*eps)+eps    

r = 6
idr = 5

nr = 2
nc = 4
fntsize=18

pSv=np.linspace(0.,1.,num=50)
deltaLv=[0, 1, 2, 4]
betaF=1.
N = 9
eps = 0.01
eps1 = 1 - eps
rv=np.array([3,6,8])


fig,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all', figsize=(10,6))
fig.subplots_adjust(hspace=0.4, wspace=0.2)
nticksY=6
nticksX=3

cmap = plt.get_cmap('viridis')

file = './newtests/1bit/leadership/multileader/res_fm4_df0'
data = np.load(file + '.npy')

ax = axs[0,0]
f=-4
for iddl, deltaL in enumerate(deltaLv):
    deltaF = 0
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

            res[idps] += benefit * data[idr, iddl, idps, strat]

    ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
    ax.set_xticklabels(["0", "0.5", "1"], fontsize=fntsize-2)
    ax.set_yticks(np.linspace(0, 1, nticksY))
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
    ax.set_ylim(0.0, 1.0)
    ax.plot(res, label='$\Delta_l=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

    ax.set_ylabel(r'cooperation level', fontsize=fntsize)
    ax.text(-60,0.53,r"$f=-4$", size=fntsize)
    ax.text(-60,0.39,r"$\Delta_f=0$", size=fntsize)
    # ax.set_xlabel(r'$p_s$', fontsize=fntsize)
    ax.set_title("B", size=fntsize, fontweight="bold")

file = './newtests/1bit/leadership/multileader/res_fm4_df8'
data = np.load(file + '.npy')

ax = axs[1,0]
f=-4
for iddl, deltaL in enumerate(deltaLv):
    deltaF = 8
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

            res[idps] += benefit * data[idr, iddl, idps, strat]

    ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
    ax.set_xticklabels(["0", "0.5", "1"], fontsize=fntsize-2)
    ax.set_yticks(np.linspace(0, 1, nticksY))
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
    ax.set_ylim(0.0, 1.0)
    ax.plot(res, label='$\Delta_l=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

    ax.set_ylabel(r'cooperation level', fontsize=fntsize)
    ax.text(-60,0.53,r"$f=-4$", size=fntsize)
    ax.text(-60,0.39,r"$\Delta_f=8$", size=fntsize)
    ax.set_xlabel(r'$p_s$', fontsize=fntsize)
    # ax.set_title("1 bit" % rv[idr], size=13)



file = './newtests/2bits/strengthstrat/multileader/res_fm4_df0'
data = np.load(file + '.npy')

ax = axs[0,1]
f=-4
for iddl, deltaL in enumerate(deltaLv):
    deltaF = 0
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

            res[idps] += benefit * data[idr, iddl, idps, strat]

    ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
    ax.set_xticklabels(["0", "0.5", "1"], fontsize=fntsize-2)
    ax.set_yticks(np.linspace(0, 1, nticksY))
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
    ax.set_ylim(0.0, 1.0)
    ax.plot(res, label='$\Delta_l=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

    # ax.set_xlabel(r'$p_s$', fontsize=fntsize)
    ax.set_title("S", size=fntsize, fontweight="bold")

file = './newtests/2bits/strengthstrat/multileader/res_fm4_df8'
data = np.load(file + '.npy')

ax = axs[1,1]
f=-4
for iddl, deltaL in enumerate(deltaLv):
    deltaF = 8
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

            res[idps] += benefit * data[idr, iddl, idps, strat]

    ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
    ax.set_xticklabels(["0", "0.5", "1"], fontsize=fntsize-2)
    ax.set_yticks(np.linspace(0, 1, nticksY))
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
    ax.set_ylim(0.0, 1.0)
    ax.plot(res, label='$\Delta_l=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

    ax.set_xlabel(r'$p_s$', fontsize=fntsize)
    # ax.set_title("1 bit" % rv[idr], size=13)



file = './newtests/2bits/leadstrat/multileader/res_fm4_df0'
data = np.load(file + '.npy')

ax = axs[0,2]
f=-4
for iddl, deltaL in enumerate(deltaLv):
    deltaF = 0
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

            res[idps] += benefit * data[idr, iddl, idps, strat]

    ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
    ax.set_xticklabels(["0", "0.5", "1"], fontsize=fntsize-2)
    ax.set_yticks(np.linspace(0, 1, nticksY))
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
    ax.set_ylim(0.0, 1.0)
    ax.plot(res, label='$\Delta_l=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

    # ax.set_xlabel(r'$p_s$', fontsize=fntsize)
    ax.set_title("L", size=fntsize, fontweight="bold")

file = './newtests/2bits/leadstrat/multileader/res_fm4_df8'
data = np.load(file + '.npy')

ax = axs[1,2]
f=-4
for iddl, deltaL in enumerate(deltaLv):
    deltaF = 8
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

            res[idps] += benefit * data[idr, iddl, idps, strat]
            
    ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
    ax.set_xticklabels(["0", "0.5", "1"], fontsize=fntsize-2)
    ax.set_yticks(np.linspace(0, 1, nticksY))
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
    ax.set_ylim(0.0, 1.0)
    ax.plot(res, label='$\Delta_l=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

    ax.set_xlabel(r'$p_s$', fontsize=fntsize)
    # ax.set_title("1 bit" % rv[idr], size=13)


file = './newtests/4bits/multileader/res_fm4_df0'
data = np.load(file + '.npy')

ax = axs[0,3]
f=-4
for iddl, deltaL in enumerate(deltaLv):
    deltaF = 0
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

            res[idps] += benefit * data[idr, iddl, idps, strat]

    ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
    ax.set_xticklabels(["0", "0.5", "1"], fontsize=fntsize-2)
    ax.set_yticks(np.linspace(0, 1, nticksY))
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
    ax.set_ylim(0.0, 1.0)
    ax.plot(res, label='$\Delta_l=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

    # ax.set_xlabel(r'$p_s$', fontsize=fntsize)
    ax.set_title("S + L", size=fntsize, fontweight="bold")

file = './newtests/4bits/multileader/res_fm4_df8'
data = np.load(file + '.npy')

ax = axs[1,3]
f=-4
for iddl, deltaL in enumerate(deltaLv):
    deltaF = 8
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

            res[idps] += benefit * data[idr, iddl, idps, strat]

    ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
    ax.set_xticklabels(["0", "0.5", "1"], fontsize=fntsize-2)
    ax.set_yticks(np.linspace(0, 1, nticksY))
    ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
    ax.set_ylim(0.0, 1.0)
    ax.plot(res, label='$\Delta_l=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

    ax.set_xlabel(r'$p_s$', fontsize=fntsize)



legend_elements = [Line2D([], [], marker='None', label='$\Delta_l:$', linestyle='None')]
legend_elements += [Line2D([], [], marker='s', color=cmap((idx)/(len(deltaLv))), label='%d'%deltaLv[idx],
                          markerfacecolor=cmap((idx)/(len(deltaLv))), markersize=10, linestyle='None') for idx in range(len(deltaLv))]
plt.legend( loc='upper center', bbox_to_anchor=(-1.3, -0.4),
          fancybox=True, shadow=False, ncol=6, columnspacing=0.0, handles=legend_elements,handletextpad=-0.3,fontsize=fntsize)
plt.savefig('./newtests/summary_fdf_plots_d4.png', bbox_inches='tight', dpi=300)

# plt.show()