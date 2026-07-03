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

nr = 4
nc = 3
fntsize=18

pSv=np.linspace(0.,1.,num=50)
deltaLv=[0, 1, 2, 4]
f=0
betaF=1.
N = 9
eps = 0.01
eps1 = 1 - eps
rv=np.array([3,6,8])


fig,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all', figsize=(12,13))
fig.subplots_adjust(hspace=0.2, wspace=0.2)
nticksY=6
nticksX=5

cmap = plt.get_cmap('viridis')

file = './newtests/1bit/leadership/singleleader/res_2bits_singleleader'
data = np.load(file + '.npy')


x_range = pSv.shape[0] - 1  # e.g. 100
y_range = 1.0  # ylim is 0 to 1

desired_physical_ratio = 12/13  # height/width of the box you want

for idr, r in enumerate(rv):
    ax=axs[0,idr]
    ax.set_aspect(desired_physical_ratio * x_range / y_range)

    for iddl, deltaL in enumerate(deltaLv):
        deltaF = deltaL
        fs = 1/(1+np.exp(-deltaL))
        fw=1-fs
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
                    nw1 = N - ns1
                    Nw = N - Ns
                    pNs = math.factorial(N)/(math.factorial(Ns)*math.factorial(Nw)) * pS**Ns * pW**Nw     
                    Nwc = (N-ns1)*s1
                    Nsc = ns1*s1
                    Nwd = (N-ns1)*(1-s1)
                    Nsd = ns1*(1-s1)

                    benefit_s = 0
                    benefit_w = 0
                    cost_s = 0
                    cost_w = 0

                    if Nw > 0:
                        benefit_w = (
                            (Nwc/Nw)*( # leader is a cooperator
                                eps1 + 
                                (1-pF[0,0])*((Nwc-1)*eps1 + Nwd*eps)+
                                (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                pF[0,0]*(Nw-1)*(eps1**2+eps**2)+pF[1,0]*Ns*(eps1**2+eps**2)
                            ) + (Nwd/Nw)*( # leader is a defector
                                eps + 
                                (1-pF[0,0])*(Nwc*eps1 + (Nwd-1)*eps)+
                                (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                pF[0,0]*(Nw-1)*(2*eps*eps1) + pF[1,0]*Ns*(2*eps*eps1)
                            )
                        )

                    if Ns > 0:

                        benefit_s = (
                            (Nsc/Ns)*( # leader is a cooperator
                                eps1 + 
                                (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                (1-pF[1,1])*((Nsc-1)*eps1 + Nsd*eps)+
                                pF[0,1]*Nw*(eps1**2+eps**2)+pF[1,1]*(Ns-1)*(eps1**2+eps**2)
                            ) + (Nsd/Ns)*( # leader is a defector
                                eps + 
                                (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                (1-pF[1,1])*(Nsc*eps1 + (Nsd-1)*eps)+
                                pF[0,1]*Nw*(2*eps*eps1) + pF[1,1]*(Ns-1)*(2*eps*eps1)
                            )
                        )

                    benefit += pNs*(((Nw*fw)/(Nw*fw+Ns*fs))*benefit_w + ((Ns*fs)/(Nw*fw+Ns*fs))*benefit_s)

                res[idps] += (benefit/N) * data[r-1, iddl, idps, strat]            

        ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fntsize-2)
        ax.set_yticks(np.linspace(0, 1, nticksY))
        ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
        ax.set_ylim(0.0, 1.0)
        ax.plot(res, label='$\Delta_f=\Delta_f=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

        if idr == 0: ax.set_ylabel(r'cooperation level', fontsize=fntsize)
        if idr==nr-1: ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idr==0: ax.text(-30, 0.465, "B", fontsize=fntsize, fontweight="bold")
        ax.text(18, 1.1, 'r=%d' % rv[idr], fontsize=fntsize)

file = './newtests/2bits/strengthstrat/singleleader/res_2bits_singleleader'
data = np.load(file + '.npy')

for idr, r in enumerate(rv):
    ax=axs[1,idr]
    ax.set_aspect(desired_physical_ratio * x_range / y_range)

    for iddl, deltaL in enumerate(deltaLv):
        deltaF = deltaL
        fs = 1/(1+np.exp(-deltaL))
        fw=1-fs
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
                for Ns in range(N+1):
                    ns1 = Ns
                    nw1 = N - ns1
                    Nw = N - Ns
                    pNs = math.factorial(N)/(math.factorial(Ns)*math.factorial(Nw)) * pS**Ns * pW**Nw     
                    Nwc = (N-ns1)*s1[0]
                    Nsc = ns1*s1[1]
                    Nwd = (N-ns1)*(1-s1[0])
                    Nsd = ns1*(1-s1[1])

                    benefit_s = 0
                    benefit_w = 0
                    cost_s = 0
                    cost_w = 0

                    if Nw > 0:
                        benefit_w = (
                            (Nwc/Nw)*( # leader is a cooperator
                                eps1 + 
                                (1-pF[0,0])*((Nwc-1)*eps1 + Nwd*eps)+
                                (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                pF[0,0]*(Nw-1)*(eps1**2+eps**2)+pF[1,0]*Ns*(eps1**2+eps**2)
                            ) + (Nwd/Nw)*( # leader is a defector
                                eps + 
                                (1-pF[0,0])*(Nwc*eps1 + (Nwd-1)*eps)+
                                (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                pF[0,0]*(Nw-1)*(2*eps*eps1) + pF[1,0]*Ns*(2*eps*eps1)
                            )
                        )

                    if Ns > 0:

                        benefit_s = (
                            (Nsc/Ns)*( # leader is a cooperator
                                eps1 + 
                                (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                (1-pF[1,1])*((Nsc-1)*eps1 + Nsd*eps)+
                                pF[0,1]*Nw*(eps1**2+eps**2)+pF[1,1]*(Ns-1)*(eps1**2+eps**2)
                            ) + (Nsd/Ns)*( # leader is a defector
                                eps + 
                                (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                (1-pF[1,1])*(Nsc*eps1 + (Nsd-1)*eps)+
                                pF[0,1]*Nw*(2*eps*eps1) + pF[1,1]*(Ns-1)*(2*eps*eps1)
                            )
                        )

                    benefit += pNs*(((Nw*fw)/(Nw*fw+Ns*fs))*benefit_w + ((Ns*fs)/(Nw*fw+Ns*fs))*benefit_s)

                res[idps] += (benefit/N) * data[r-1, iddl, idps, strat]
            

        ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fntsize-2)
        ax.set_yticks(np.linspace(0, 1, nticksY))
        ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
        ax.set_ylim(0.0, 1.0)
        ax.plot(res, label='$\Delta_f=\Delta_f=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))

        if idr==nr-1: ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idr == 0: ax.set_ylabel(r'cooperation level', fontsize=fntsize)
        if idr==0: ax.text(-30, 0.465, "S", fontsize=fntsize, fontweight="bold")

file = './newtests/2bits/leadstrat/singleleader/res_2bits_singleleader'
data = np.load(file + '.npy')

for idr, r in enumerate(rv):
    ax=axs[2,idr]
    ax.set_aspect(desired_physical_ratio * x_range / y_range)

    for iddl, deltaL in enumerate(deltaLv):
        deltaF = deltaL
        fs = 1/(1+np.exp(-deltaL))
        fw=1-fs
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
                    nw1 = N - ns1
                    Nw = N - Ns
                    pNs = math.factorial(N)/(math.factorial(Ns)*math.factorial(Nw)) * pS**Ns * pW**Nw     

                    Nwcnl = (N-ns1)*s1[0]
                    Nscnl = ns1*s1[0]
                    Nwdnl = (N-ns1)*(1-s1[0])
                    Nsdnl = ns1*(1-s1[0])

                    benefit_s = 0
                    benefit_w = 0
                    cost_s = 0
                    cost_w = 0

                    if Nw > 0:
                        benefit_w += (
                            (nw1/Nw)*( # leader is of strategy 1
                                aeps(s1[1], eps) + 
                                (1-pF[0,0])*((Nwcnl-s1[0])*eps1 + (Nwdnl-(1-s1[0]))*eps)+
                                (1-pF[1,0])*(Nscnl*eps1 + Nsdnl*eps)+
                                pF[0,0]*(Nw-1)*(s1[1]*(eps1**2+eps**2) + (1-s1[1])*(2*eps1*eps))
                                +pF[1,0]*Ns*(s1[1]*(eps1**2+eps**2) + (1-s1[1])*(2*eps1*eps))
                            )
                        )

                    if Ns > 0:
                        benefit_s += (
                            (ns1/Ns)*( # leader is of strategy 1
                                aeps(s1[1], eps) + 
                                (1-pF[0,1])*(Nwcnl*eps1 + Nwdnl*eps)+
                                (1-pF[1,1])*((Nscnl-s1[0])*eps1 + (Nsdnl-(1-s1[0]))*eps)+
                                pF[0,1]*Nw*(s1[1]*(eps1**2+eps**2) + (1-s1[1])*(2*eps1*eps))
                                +pF[1,1]*(Ns-1)*(s1[1]*(eps1**2+eps**2) + (1-s1[1])*(2*eps1*eps))
                            )
                        )

                    benefit += pNs*(((Nw*fw)/(Nw*fw+Ns*fs))*benefit_w + ((Ns*fs)/(Nw*fw+Ns*fs))*benefit_s)

                res[idps] += (benefit/N) * data[r-1, iddl, idps, strat]
            

        ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fntsize-2)
        ax.set_yticks(np.linspace(0, 1, nticksY))
        ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
        ax.set_ylim(0.0, 1.0)
        ax.plot(res, label='$\Delta_f=\Delta_f=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))
        if idr==nr-1: ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idr == 0: ax.set_ylabel(r'cooperation level', fontsize=fntsize)
        if idr==0: ax.text(-30, 0.465, "L", fontsize=fntsize, fontweight="bold")

file = './newtests/4bits/singleleader/res_singleleader'
data = np.load(file + '.npy')

for idr, r in enumerate(rv):
    ax=axs[3,idr]
    ax.set_aspect(desired_physical_ratio * x_range / y_range)

    for iddl, deltaL in enumerate(deltaLv):
        deltaF = deltaL
        fs = 1/(1+np.exp(-deltaL))
        fw=1-fs
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
                    Nw = N - Ns
                    ns1 = Ns
                    nw1 = N - ns1
                    pNs = math.factorial(N)/(math.factorial(Ns)*math.factorial(Nw)) * pS**Ns * pW**Nw

                    Nwc = nw1*s1[2]
                    Nsc = ns1*s1[3]
                    Nwd = nw1*(1-s1[2])
                    Nsd = ns1*(1-s1[3])

                    N1wcl = nw1*s1[0]
                    N1wdl = nw1*(1-s1[0])
                    N1scl = ns1*s1[1]
                    N1sdl = ns1*(1-s1[1])

                    benefit_s = 0
                    benefit_w = 0
                    cost_s = 0
                    cost_w = 0

                    if Nw > 0:

                        if nw1 > 0:
                            benefit_w += (nw1/Nw)*( # leader is of strategy 1
                                    (N1wcl/nw1)*( # leader is a cooperator
                                        eps1+
                                        (1-pF[0,0])*((Nwc-s1[2])*eps1 + (Nwd-(1-s1[2]))*eps)+
                                        (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                        pF[0,0]*(Nw-1)*(eps1**2+eps**2)+pF[1,0]*Ns*(eps1**2+eps**2)
                                    )+
                                    (N1wdl/nw1)*( # leader is a defector
                                        eps+
                                        (1-pF[0,0])*((Nwc-s1[2])*eps1 + (Nwd-(1-s1[2]))*eps)+
                                        (1-pF[1,0])*(Nsc*eps1 + Nsd*eps)+
                                        pF[0,0]*(Nw-1)*(2*eps1*eps)+pF[1,0]*Ns*(2*eps1*eps)
                                    )
                                )
                            
                    if Ns > 0:

                        if ns1 > 0:
                            benefit_s += (ns1/Ns)*( # leader is of strategy 1
                                (N1scl/ns1)*( # leader is a cooperator
                                    eps1+
                                    (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                    (1-pF[1,1])*((Nsc-s1[3])*eps1 + (Nsd-(1-s1[3]))*eps)+
                                    pF[0,1]*Nw*(eps1**2+eps**2)+pF[1,1]*(Ns-1)*(eps1**2+eps**2)
                                )+
                                (N1sdl/ns1)*( # leader is a defector
                                    eps+
                                    (1-pF[0,1])*(Nwc*eps1 + Nwd*eps)+
                                    (1-pF[1,1])*((Nsc-s1[3])*eps1 + (Nsd-(1-s1[3]))*eps)+
                                    pF[0,1]*Nw*(2*eps1*eps)+pF[1,1]*(Ns-1)*(2*eps1*eps)
                                )
                            )

                    benefit += pNs*(((Nw*fw)/(Nw*fw+Ns*fs))*benefit_w + ((Ns*fs)/(Nw*fw+Ns*fs))*benefit_s)
                    cost += pNs*(((Nw*fw)/(Nw*fw+Ns*fs))*cost_w + ((Ns*fs)/(Nw*fw+Ns*fs))*cost_s)


                res[idps] += (benefit/N) * data[r-1, iddl, idps, strat]
            

        ax.set_xticks(np.linspace(0, pSv.shape[0]-1, nticksX))
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fntsize-2)
        ax.set_yticks(np.linspace(0, 1, nticksY))
        ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=fntsize-2)
        ax.set_ylim(0.0, 1.0)
        ax.plot(res, label='$\Delta_f=\Delta_f=%d$'%deltaF, color=cmap((iddl)/(len(deltaLv))))
        if idr==nr-1: ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idr == 0: ax.set_ylabel(r'cooperation level', fontsize=fntsize)
        if idr==0: ax.text(-35, 0.465, "S+L", fontsize=fntsize, fontweight="bold")

legend_elements = [Line2D([], [], marker='None', label='$\Delta:$', linestyle='None')]
legend_elements += [Line2D([], [], marker='s', color=cmap((idx)/(len(deltaLv))), label='%d'%deltaLv[idx],
                          markerfacecolor=cmap((idx)/(len(deltaLv))), markersize=10, linestyle='None') for idx in range(len(deltaLv))]
plt.legend( loc='upper center', bbox_to_anchor=(-0.85, -0.28),
          fancybox=True, shadow=False, ncol=6, columnspacing=0.0, handles=legend_elements,handletextpad=-0.3,fontsize=fntsize)
plt.savefig('./newtests/summary_cl_plots_singleleader_d4.png', bbox_inches='tight', dpi=300)

# plt.show()