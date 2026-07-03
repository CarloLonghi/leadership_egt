import numpy as np

def plotCOOPheat(deltaFv,pSv,rv,folder):
# Input: MAT (matrix from "coop_pF_r" function), pFv, rv ,Mv (vectors with values of pF, r, and M), label (name for the output file)
# Output: heatmap plot of the fraction of cooperators as a function of pF and r, for different M
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib
    from matplotlib.colors import LinearSegmentedColormap

    fntsize=18
    nr=4
    nc=3
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all', figsize=(12,13))
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    files = ["1bit/leadership/singleleader/res_2bits_singleleader_deltaps.npy", 
             "2bits/strengthstrat/singleleader/res_2bits_singleleader_deltaps.npy",
             "2bits/leadstrat/singleleader/res_2bits_singleleader_deltaps.npy"]
    
    titles = ["B", "S", "L"]

    for ids in range(3):
        MAT = np.load(folder + files[ids])
        for idx in range(len(rv)):
            ax=axs[ids, idx]
            nstrats = 2 if ids == 0 else 4
            cmaps=['Greens','Purples'] if ids==0 else ['Greens','Reds','Blues','Purples']
            for strat in range(nstrats):
                step=0.025
                levels = np.arange(0.5-step, 1., step) + step
                h=ax.contourf(MAT[:,:,rv[idx]-1,strat],levels,cmap=cmaps[strat], origin='lower',)
            #h=ax.imshow(MAT[:,:,k],origin='lower', interpolation='none',aspect='auto',vmin=0,vmax=4)
            nticksY=6
            nticksX=5
            ax.set_xticks(np.linspace(0, MAT.shape[1]-1, nticksX))
            ax.set_yticks(np.linspace(0, MAT.shape[0]-1, nticksY))
            ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fntsize-2)
            ax.set_yticklabels(['0', '1', '2', '3', '4', '5'], fontsize=fntsize-2)
 
            ax.set_ylim(0,MAT.shape[1]-1)
            if idx==0: 
                ax.text(-23, 14, titles[ids], fontsize=fntsize, fontweight="bold")
                ax.set_ylabel(r'$\Delta$', fontsize=fntsize) 
            if ids==0: 
                ax.text(18, 35, 'r=%d' % rv[idx], fontsize=fntsize)

            # insert markers for invasion graphs
            # if idx == 0:
            #     points_x = [MAT.shape[1] // 2 - 0.5, MAT.shape[1] // 2 - 0.5, MAT.shape[1] // 2 - 0.5]
            #     points_y = [0.5, MAT.shape[0] // 8 - 0.5, MAT.shape[0] // 4 - 0.5]
            #     ax.scatter(points_x, points_y, color='goldenrod', marker='o')

            # add strategies labels
            if ids == 1:
                if idx == 0:
                    ax.text(37, 27, "[*0]", fontsize=16)
                if idx == 1:
                    ax.text(20, 27, "[*1]", fontsize=16)
                if idx == 2:
                    ax.text(26, 27, "[*1]", fontsize=16)

            if ids == 2:
                if idx == 2:
                    ax.text(9, 27, "[*1]", fontsize=16) 

    file = "4bits/singleleader/singlestrat_new.npy"

    MAT = np.load(folder + file)
    MAT[44,12,2,8:] = [0.24750185, 0.00112731, 0.24605334, 0.00112072, 0.24102037, 0.00109779, 0.23960979, 0.00109137]

    data = MAT[:,:,:,[0,8,12,14]]

    colors = ['#f7f7f7', '#d4d4f7', "#c890ee", "#AD3896", "#9c288d", "#6B0861"]
    n_bins = 256
    cmap_name = 'Purples_custom'
    LC_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

    D = data[...,0]
    SLC = data[...,1]
    LC = data[...,2] + data[...,3]
    LC_overlap = np.where((D < 0.5) & (LC < 0.5) & (SLC < 0.5), LC, 0)
    LC_no_overlap = np.where(LC > 0.5, LC, 0)
    SLC_overlap = np.where((D < 0.5) & (LC < 0.5) & (SLC < 0.5), SLC, 0)
    SLC_no_overlap = np.where(SLC > 0.5, SLC, 0) 
    D_overlap = np.where((D < 0.5) & (LC < 0.5) & (SLC < 0.5), D, 0)
    D_no_overlap = np.where(D > 0.5, D, 0)      

    LC_overlap = LC
    SLC_overlap = SLC
    D_overlap = D     
    
    for idx in range(len(rv)):
        ax=axs[3, idx]
        cmaps= ['Greens','Reds','Blues','Purples']
        step=0.025
        levels = np.arange(0.5-step, 0.66666, step) + step
        levels_overlap = np.arange(0.05, 0.5, step) + step
        # h=ax.contourf(np.sum(MAT[:32,:,rv[idx]-1,:8], -1),levels,cmap=cmaps[0], origin='lower', extend="max")
        # h=ax.contourf(np.sum(MAT[:32,:,rv[idx]-1,8:16], -1),levels,cmap='Blues', origin='lower', extend="max")
        # h=ax.contourf(MAT[:32,:,rv[idx]-1,8],levels,cmap='winter', origin='lower', alpha=0.4, extend="max")
        # h=ax.contourf(MAT[:32,:,rv[idx]-1,0],levels,cmap='summer', origin='lower', alpha=0.5, extend="max")
        # h=ax.contourf(MAT[:32,:,rv[idx]-1,12]+MAT[:32,:,rv[idx]-1,14],levels,cmap="BuPu", origin='lower',alpha=0.2, extend="max")
        #h=ax.imshow(MAT[:,:,k],origin='lower', interpolation='none',aspect='auto',vmin=0,vmax=4)

        h=ax.contourf(D_no_overlap[...,rv[idx]-1],levels,cmap="Greens", origin='lower', extend="max", alpha=1.0)
        h=ax.contourf(LC_no_overlap[...,rv[idx]-1],levels,cmap=LC_cmap, origin='lower', extend="max", alpha=1.0)        
        h=ax.contourf(SLC_no_overlap[...,rv[idx]-1],levels,cmap="Blues", origin='lower', extend="max", alpha=1.0)

        h=ax.contourf(D_overlap[...,rv[idx]-1],levels_overlap,cmap="Greens", origin='lower', extend="max", alpha=0.15)         
        h=ax.contourf(LC_overlap[...,rv[idx]-1],levels_overlap,cmap=LC_cmap, origin='lower', extend="max", alpha=0.15)
        h=ax.contourf(SLC_overlap[...,rv[idx]-1],levels_overlap,cmap="Blues", origin='lower', extend="max", alpha=0.2)        
                
        nticksY=6
        nticksX=5
        ax.set_xticks(np.linspace(0, MAT.shape[1]-1, nticksX))
        ax.set_yticks(np.linspace(0, 32, nticksY))
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fntsize-2)
        ax.set_yticklabels(['0', '1', '2', '3', '4', '5'], fontsize=fntsize-2)
        ax.set_ylim(0,32)
        ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idx == 0: 
            ax.text(-28, 14, "S+L", fontsize=fntsize, fontweight="bold")
            ax.set_ylabel(r'$\Delta$', fontsize=fntsize)

        if idx == 0:
            ax.text(30, 26, "[***1]", fontsize=16)
        if idx == 1:
            ax.text(30, 26, "[***1]", fontsize=16)
        if idx == 2:
            ax.text(30, 26, "[***1]", fontsize=16)

    legend_prop = {"size": fntsize}

    labels = ['$\\bf{AllD}$ [0]', '$\\bf{AllC}$ [1]']
    cmaps=["Greens", "Purples"]
    patches = [mpatches.Patch(color=plt.get_cmap(cmaps[i])(0.9), label=labels[i]) for i in range(2)]    
    first_legend = ax.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.92 , 4.36),
          fancybox=True, shadow=False, ncol=1, fontsize=20, columnspacing=1.0, title_fontsize=20, frameon=True, prop=legend_prop)

    # patches += [mpatches.Patch(color="none")]
    # vertical_line = matplotlib.lines.Line2D([], [], color="black", marker='|', linestyle='None', markersize=20, markeredgewidth=1.5)
    # patches += [vertical_line,]
    labels = ['$\\bf{AllD}$ [00]', '$\\bf{WDSC}$ [01]', '$\\bf{AllC}$ [11]']
    cmaps=["Greens", "Blues", "Purples"]
    patches = [mpatches.Patch(color=plt.get_cmap(cmaps[i])(0.9), label=labels[i]) for i in range(3)]
    second_legend = plt.legend(handles=patches, loc='upper right', bbox_to_anchor=(2.1 , 3.20),
          fancybox=True, shadow=False, ncol=1, fontsize=20, columnspacing=1.0, title_fontsize=20, frameon=True, prop=legend_prop)
    
    labels = ['$\\bf{AllD}$ [00]', '$\\bf{NDLC}$ [01]']
    cmaps=["Greens", "Blues"]
    patches = [mpatches.Patch(color=plt.get_cmap(cmaps[i])(0.9), label=labels[i]) for i in range(2)]
    third_legend = plt.legend(handles=patches, loc='upper right', bbox_to_anchor=(2.06 , 1.95),
          fancybox=True, shadow=False, ncol=1, fontsize=20, columnspacing=1.0, title_fontsize=20, frameon=True, prop=legend_prop)    

    labels = ['$\\bf{AllD}$ [0000]', '$\\bf{SLC}$ [0001]', '$\\bf{LC}$ [0*11]']
    cmaps=["Greens", "Blues", LC_cmap]
    patches = [mpatches.Patch(color=plt.get_cmap(cmaps[i])(0.9), label=labels[i]) for i in range(3)]
    fourth_legend = plt.legend(handles=patches, loc='upper right', bbox_to_anchor=(2.1 , 0.8),
          fancybox=True, shadow=False, ncol=1, fontsize=20, columnspacing=1.0, title_fontsize=20, frameon=True, prop=legend_prop)        

    # cmaps=['Greens', 'summer','Blues', 'winter', 'BuPu']
    # colors = [plt.get_cmap("Greens")(0.9), (26/255, 199/255, 57/255, 0.85), plt.get_cmap("Blues")(0.9), (34/255, 51/255, 215/255, 1), (15/255, 12/255, 120/255, 0.85)]

    ax.add_artist(third_legend)
    ax.add_artist(second_legend)
    ax.add_artist(first_legend)

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
    f.savefig('./newtests/summary_single_sd_4bits_d5.png',bbox_inches='tight',dpi=300)
    #plt.show()
    f.clf()     
    return


deltaLv=np.linspace(0,8,num=50)
pSv=np.linspace(0.,1.,num=50)
rv=np.array([3,6,8])

folder = "./newtests/"

plotCOOPheat(deltaLv,pSv,rv,folder)