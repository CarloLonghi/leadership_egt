import numpy as np

def plotCOOPheat(deltaFv,pSv,rv,folder):
# Input: MAT (matrix from "coop_pF_r" function), pFv, rv ,Mv (vectors with values of pF, r, and M), label (name for the output file)
# Output: heatmap plot of the fraction of cooperators as a function of pF and r, for different M
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib
    fntsize=18
    nr=1
    nc=3
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all', figsize=(12,4))
    f.subplots_adjust(hspace=0.4, wspace=0.2)

    file = "4bits/singleleader/singlestrat_new.npy"

    MAT = np.load(folder + file)
    MAT[44,12,2,8:] = [0.24750185, 0.00112731, 0.24605334, 0.00112072, 0.24102037, 0.00109779, 0.23960979, 0.00109137]
    for idx in range(len(rv)):
        ax=axs[idx]
        cmaps= ['Greens','Reds','Blues','Purples']
        step=0.025
        levels = np.arange(0.5-step, 1., step) + step
        h=ax.contourf(np.sum(MAT[:,:,rv[idx]-1,:8], -1),levels,cmap=cmaps[0], origin='lower', extend="max")
        h=ax.contourf(np.sum(MAT[:,:,rv[idx]-1,8:16], -1),levels,cmap='Blues', origin='lower', extend="max")
        h=ax.contourf(MAT[:,:,rv[idx]-1,8],levels,cmap='winter', origin='lower', alpha=0.5, extend="max")
        h=ax.contourf(MAT[:,:,rv[idx]-1,0],levels,cmap='summer', origin='lower', alpha=1.0, extend="max")
        h=ax.contourf(MAT[:,:,rv[idx]-1,12]+MAT[:,:,rv[idx]-1,14],levels,cmap="BuPu", origin='lower',alpha=0.2, extend="max")
        #h=ax.imshow(MAT[:,:,k],origin='lower', interpolation='none',aspect='auto',vmin=0,vmax=4)
        nticksY=5
        nticksX=3
        ax.set_xticks(np.linspace(0, MAT.shape[1]-1, nticksX))
        ax.set_yticks(np.linspace(0, MAT.shape[0]-1, nticksY))
        ax.set_xticklabels(np.linspace(pSv[0],pSv[-1],nticksX), fontsize=18)
        ax.set_yticklabels(np.linspace(deltaFv[0],deltaFv[-1],nticksY), fontsize=18)
        ax.set_ylim(0,MAT.shape[1]-1)
        ax.set_title(f"r={rv[idx]}", fontsize=fntsize)
        # if idx==0: ax.set_title(titles[ids], size=fntsize)
        # if idx==nr-1: ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        # if ids==0: 
        #     ax.set_ylabel(r'$\Delta_f, \Delta_l$', fontsize=fntsize)
        #     ax.text(-35, 22.5, "$r=%d$" % rv[idx], fontsize=fntsize)

        # insert markers for invasion graphs
        # if idx == 0:
        #     points_x = [MAT.shape[1] // 2 - 0.5, MAT.shape[1] // 2 - 0.5, MAT.shape[1] // 2 - 0.5]
        #     points_y = [0.5, MAT.shape[0] // 8 - 0.5, MAT.shape[0] // 4 - 0.5]
        #     ax.scatter(points_x, points_y, color='goldenrod', marker='o')
        
    labels = ['***0', '***1']
    cmaps = ['Greens','Blues']
    patches = [mpatches.Patch(color=plt.get_cmap(cmaps[i])(0.9), label=labels[i]) for i in range(2)]
    plt.legend(handles=patches, loc='upper center', bbox_to_anchor=(-0.7 , -0.15),
          fancybox=True, shadow=False, ncol=4, fontsize=20, columnspacing=1)
    # ax.text(-12,125,"$\Delta_f=1$", size=13)
    # ax.text(-12,57,"$\Delta_l=1$", size=13)

    # f.subplots_adjust(right=0.9)
    # cbar_ax = f.add_axes([0.93, 0.15, 0.02, 0.7])
    # f.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(vmin=0.5, vmax=1), cmap='Greys'), cax=cbar_ax)
    # cbar_ax.tick_params(labelsize=18)

    # box = ax.get_position()
    # fi.set_position([box.x0, box.y0 + box.height * 0.1,
    #              box.width, box.height * 0.9])
    #f.subplots_adjust(bottom=0.1)
    # cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
    # cb = f.colorbar(h, cax=cbar_ax)
    # cb.set_ticks([1,2,3,4])
    # cb.set_ticklabels(['ALLD','WCSD','WDSC','ALLC'])

#cb=f.colorbar(h, fraction=0.1,format='%.2f')
    #cb.set_label(label=r'$f_C$')
    f.savefig('./newtests/summary_4bits_sd_singlleader.png',bbox_inches='tight',dpi=300)
    #plt.show()
    f.clf()     
    return


deltaLv=np.linspace(0,8,num=50)
pSv=np.linspace(0.,1.,num=50)
rv=np.array([3,6,8])

folder = "./newtests/"

plotCOOPheat(deltaLv,pSv,rv,folder)