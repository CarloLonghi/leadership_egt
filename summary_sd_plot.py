import numpy as np

def plotCOOPheat(deltaFv,pSv,rv,folder):
# Input: MAT (matrix from "coop_pF_r" function), pFv, rv ,Mv (vectors with values of pF, r, and M), label (name for the output file)
# Output: heatmap plot of the fraction of cooperators as a function of pF and r, for different M
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap

    fntsize=18
    nr=4
    nc=3
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all', figsize=(10,13))
    f.subplots_adjust(hspace=0.4, wspace=0.2)

    files = ["1bit/leadership/multileader/multisd_2bits_multileader.npy", 
             "2bits/strengthstrat/multileader/multisd_2bits_multileader.npy",
             "2bits/leadstrat/multileader/multisd_2bits_multileader.npy"]
    
    titles = ["B", "S", "L"]

    n_bins = 256

    # Define color transitions from white to purple
    # colors = ['#f7f7f7', '#d4d4f7', "#99b0ee", "#6180d6", "#304dad", "#08298B"]
    # cmap_name = 'c1'
    # blue_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

    colors = ['#f7f7f7', "#f7d4f4", "#ee99e3", "#d661bd", "#a62faa", "#790B83"]
    cmap_name = 'c2'
    violet_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

    # colors = ['#f7f7f7', "#f7f2d4", "#eee199", "#d6c061", "#ad9732", "#926F07"]
    # cmap_name = 'c1'
    # yellow_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

    # colors = ["#f7f7f7", "#d4f7da", "#99eea4", "#61d671", "#2b9b34", "#0A8629"]
    # cmap_name = 'c1'
    # green_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)    

    colors = ['#f7f7f7', '#d4d4f7', "#99b0ee", "#6180d6", "#3353bb", "#0F3197"]
    cmap_name = 'c1'
    blue_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

    colors = ['#f7f7f7', "#f7d4d9", "#ee99a4", "#d6616b", "#bb3345", "#920B1D"]
    cmap_name = 'c2'
    red_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

    colors = ['#f7f7f7', "#dff7d4", "#bbee99", "#77c559", "#409626", "#257E08"]
    cmap_name = 'c1'
    green_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)    

    for ids in range(3):
        MAT = np.load(folder + files[ids])
        for idx in range(len(rv)):
            ax=axs[ids, idx]
            nstrats = 2 if ids == 0 else 4
            cmaps=["Greens", "Purples"] if ids==0 else ["Greens", "Reds", "Blues", "Purples"]
            for strat in range(nstrats):
                step=0.025
                levels = np.arange(0.5-step, 1., step) + step
                h=ax.contourf(MAT[:32,:,rv[idx]-1,strat],levels,cmap=cmaps[strat], origin='lower',)                
                # alpha_array = np.clip(MAT[:32,:,rv[idx]-1,strat], 0, 1)
                # ax.pcolormesh(MAT[:32,:,rv[idx]-1,strat], cmap=cmaps[strat], alpha=alpha_array, shading='auto')

            #h=ax.imshow(MAT[:,:,k],origin='lower', interpolation='none',aspect='auto',vmin=0,vmax=4)
            nticksY=6
            nticksX=5
            ax.set_xticks(np.linspace(0, MAT.shape[1]-1, nticksX))
            ax.set_yticks(np.linspace(0, 32, nticksY))
            ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fntsize-2)
            ax.set_yticklabels(['0', '1', '2', '3', '4', '5'], fontsize=fntsize-2)
            
            ax.set_ylim(0,32)
            if idx==0: 
                ax.text(-30, 15, titles[ids], fontsize=fntsize, fontweight="bold")
                ax.set_ylabel(r'$\Delta$', fontsize=fntsize) 
            if ids==0: 
                ax.text(18, 40, 'r=%d' % rv[idx], fontsize=fntsize)

            # insert markers for invasion graphs
            # if idx == 0:
            #     points_x = [MAT.shape[1] // 2 - 0.5, MAT.shape[1] // 2 - 0.5, MAT.shape[1] // 2 - 0.5]
            #     points_y = [0.5, MAT.shape[0] // 8 - 0.5, MAT.shape[0] // 4 - 0.5]
            #     ax.scatter(points_x, points_y, color='goldenrod', marker='o')

            # add strategies labels
            if ids == 1:
                if idx == 0:
                    ax.text(40, 26, "[*0]", fontsize=16)
                if idx == 1:
                    ax.text(40, 26, "[*0]", fontsize=16)
                    ax.text(20, 26, "[*1]", fontsize=16)
                if idx == 2:
                    ax.text(40, 26, "[*0]", fontsize=16)
                    ax.text(26, 26, "[*1]", fontsize=16)

            if ids == 2:
                if idx == 0:
                    ax.text(34, 27, "[*0]", fontsize=16)
                if idx == 1:
                    ax.text(38, 26, "[*0]", fontsize=16)
                    ax.text(20, 26, "[*1]", fontsize=16)
                if idx == 2:
                    ax.text(40, 26, "[*0]", fontsize=16)
                    ax.text(26, 26, "[*1]", fontsize=16) 

    file = "4bits/multileader/singlestrat_new.npy"

    MAT = np.load(folder + file)

    MAT[44,12,2,8:] = [0.24750185, 0.00112731, 0.24605334, 0.00112072, 0.24102037, 0.00109779, 0.23960979, 0.00109137]

    data = MAT[:32,:,:,[0,8,12,14]]
    # ms = np.argmax(data, axis=-1)

    # for i in range(data.shape[0]):
    #     for j in range(data.shape[1]):
    #         for k in range(data.shape[2]):
    #             mask = np.arange(data.shape[-1]) != ms[i,j,k]
    #             data[i,j,k,mask] = 0

    # Define color transitions from white to purple
    # colors = ['#f7f7f7', '#d4d4f7', "#909aee", "#4575db", "#2955cf", "#0737D3"]
    # colors = ['#f7f7f7', '#d4d4f7', "#90a6ee", "#4575db", "#3661d8", "#0026A3"]
    # colors = ['#f7f7f7', '#d4d4f7', "#c890ee", "#b33abe", "#8231c5", "#8800A3"]
    # colors = ['#f7f7f7', '#d4d4f7', "#c890ee", "#AD3896", "#9c288d", "#6B0861"]
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
        step=0.025
        levels = np.arange(0.5-step, 1.0, step) + step
        levels_overlap = np.arange(0.0, 0.5, step) + step
        # h=ax.contourf(np.sum(MAT[:32,:,rv[idx]-1,:8], -1),levels,cmap=cmaps[0], origin='lower', extend="max")
        # h=ax.contourf(np.sum(MAT[:32,:,rv[idx]-1,8:16], -1),levels,cmap='Blues', origin='lower', extend="max")
        # h=ax.contourf(data[...,rv[idx]-1,1],levels,cmap="Blues", origin='lower', extend="max")
        # h=ax.contourf(data[...,rv[idx]-1,0],levels,cmap="Greens", origin='lower', extend="max")
        # h=ax.contourf(data[...,rv[idx]-1,2]+data[...,rv[idx]-1,3],levels,cmap=LC_cmap, origin='lower', extend="max")

        h=ax.contourf(D_no_overlap[...,rv[idx]-1],levels,cmap="Greens", origin='lower', extend="max", alpha=1.0)
        h=ax.contourf(LC_no_overlap[...,rv[idx]-1],levels,cmap=LC_cmap, origin='lower', extend="max", alpha=1.0)        
        h=ax.contourf(SLC_no_overlap[...,rv[idx]-1],levels,cmap="Blues", origin='lower', extend="max", alpha=1.0)

        h=ax.contourf(D_overlap[...,rv[idx]-1],levels_overlap,cmap="Greens", origin='lower', extend="max", alpha=0.15)         
        h=ax.contourf(LC_overlap[...,rv[idx]-1],levels_overlap,cmap=LC_cmap, origin='lower', extend="max", alpha=0.15)
        h=ax.contourf(SLC_overlap[...,rv[idx]-1],levels_overlap,cmap="Blues", origin='lower', extend="max", alpha=0.15)        
        
        # alpha_array = np.clip(D[...,rv[idx]-1], 0, 1)
        # ax.pcolormesh(D[...,rv[idx]-1], cmap='Greens', alpha=alpha_array, shading='auto')
        # alpha_array = np.clip(LC[...,rv[idx]-1], 0, 1)
        # ax.pcolormesh(LC[...,rv[idx]-1], cmap=LC_cmap, alpha=alpha_array, shading='auto')
        # alpha_array = np.clip(SLC[...,rv[idx]-1], 0, 1)
        # ax.pcolormesh(SLC[...,rv[idx]-1], cmap='Blues', alpha=alpha_array, shading='auto')                

        nticksY=6
        nticksX=5
        ax.set_xticks(np.linspace(0, MAT.shape[1]-1, nticksX))
        ax.set_yticks(np.linspace(0, 32, nticksY))
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fntsize-2)
        ax.set_yticklabels(['0', '1', '2', '3', '4', '5'], fontsize=fntsize-2)
        ax.set_ylim(0,32)
        ax.set_xlabel(r'$p_s$', fontsize=fntsize)
        if idx == 0: 
            ax.text(-35, 15, "S+L", fontsize=fntsize, fontweight="bold")
            ax.set_ylabel(r'$\Delta$', fontsize=fntsize)

        if idx == 0:
            ax.text(30, 27, "[***0]", fontsize=16)
            ax.text(1, 27, "[***1]", fontsize=16)
        if idx == 1:
            ax.text(34, 26, "[***0]", fontsize=16)
            ax.text(9, 26, "[***1]", fontsize=16)
        if idx == 2:
            ax.text(35.5, 26, "[***0]", fontsize=16)
            ax.text(14, 26, "[***1]", fontsize=16)             

    legend_prop = {"size": fntsize}

    labels = ['$\\bf{AllD}$ [0]', '$\\bf{AllC}$ [1]']
    cmaps=["Greens", "Purples"]
    patches = [mpatches.Patch(color=plt.get_cmap(cmaps[i])(0.9), label=labels[i]) for i in range(2)]    
    first_legend = ax.legend(handles=patches, loc='upper right', bbox_to_anchor=(2.2 , 4.96),
          fancybox=True, shadow=False, ncol=1, fontsize=20, columnspacing=1.0, title_fontsize=20, frameon=True, prop=legend_prop)

    # patches += [mpatches.Patch(color="none")]
    # vertical_line = matplotlib.lines.Line2D([], [], color="black", marker='|', linestyle='None', markersize=20, markeredgewidth=1.5)
    # patches += [vertical_line,]
    labels = ['$\\bf{AllD}$ [00]', '$\\bf{WDSC}$ [01]', '$\\bf{AllC}$ [11]']
    cmaps=["Greens", "Blues", "Purples"]
    patches = [mpatches.Patch(color=plt.get_cmap(cmaps[i])(0.9), label=labels[i]) for i in range(3)]
    second_legend = plt.legend(handles=patches, loc='upper right', bbox_to_anchor=(2.4 , 3.68),
          fancybox=True, shadow=False, ncol=1, fontsize=20, columnspacing=1.0, title_fontsize=20, frameon=True, prop=legend_prop)
    
    labels = ['$\\bf{AllD}$ [00]', '$\\bf{NDLC}$ [01]']
    cmaps=["Greens", "Blues"]
    patches = [mpatches.Patch(color=plt.get_cmap(cmaps[i])(0.9), label=labels[i]) for i in range(2)]
    third_legend = plt.legend(handles=patches, loc='upper right', bbox_to_anchor=(2.36 , 2.18),
          fancybox=True, shadow=False, ncol=1, fontsize=20, columnspacing=1.0, title_fontsize=20, frameon=True, prop=legend_prop)    

    labels = ['$\\bf{AllD}$ [0000]', '$\\bf{SLC}$ [0001]', '$\\bf{LC}$ [0*11]']
    cmaps=["Greens", "Blues", LC_cmap]
    patches = [mpatches.Patch(color=plt.get_cmap(cmaps[i])(0.9), label=labels[i]) for i in range(3)]
    fourth_legend = plt.legend(handles=patches, loc='upper right', bbox_to_anchor=(2.43 , 0.87),
          fancybox=True, shadow=False, ncol=1, fontsize=20, columnspacing=1.0, title_fontsize=20, frameon=True, prop=legend_prop)        

    # cmaps=['Greens', 'summer','Blues', 'winter', 'BuPu']
    # colors = [plt.get_cmap("Greens")(0.9), (26/255, 199/255, 57/255, 0.85), plt.get_cmap("Blues")(0.9), (34/255, 51/255, 215/255, 1), (15/255, 12/255, 120/255, 0.85)]

    ax.add_artist(third_legend)
    ax.add_artist(second_legend)
    ax.add_artist(first_legend)

    # add colorbar
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
    f.savefig('./newtests/summary_multi_sd_4bits_d5_alt.png',bbox_inches='tight',dpi=300)
    #plt.show()
    f.clf()     
    return


deltaLv=np.linspace(0,8,num=50)
l = int(5 / (8 / 50))
deltaLv = deltaLv[:l]
pSv=np.linspace(0.,1.,num=50)
rv=np.array([3,6,8])

folder = "./newtests/"

plotCOOPheat(deltaLv,pSv,rv,folder)