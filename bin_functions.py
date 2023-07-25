# -*- coding: utf-8 -*-
"""
Created on Wed May 10 12:59:17 2023

@author: morri
"""
import numpy as np
import matplotlib.pyplot as plt
from h5_functions import *
from scipy.optimize import curve_fit
import time as t

#%% READING DATA

def ReadDist(filename, dir = r'C:\Users\morri\Desktop\Thesis Data\bin data', ns = 4, numberel = 101): 
    file = dir + filename
             
    fxy = np.zeros((numberel, numberel, ns))
    fxz = np.zeros((numberel, numberel, ns))
    fyz = np.zeros((numberel, numberel, ns))
    fxyz = np.zeros((numberel,numberel,numberel,ns))
    axes = np.zeros((numberel,ns))

    with open(file, 'rb') as file:
        
        header = np.fromfile(file, 'i8', count = 1)
        last_pos = file.tell() - 4 
        for i in range(1, ns+1):
            
            file.seek(last_pos)
            axes[:, i-1] = np.fromfile(file, dtype='float32', count=numberel)
            last_pos = file.tell()
        
        last_pos = file.tell()
        file.seek(last_pos)
        x1 = np.fromfile(file, dtype = 'f4', count = 1)
        
        last_pos = file.tell()
        file.seek(last_pos)
        x2 = np.fromfile(file, dtype = 'f4', count = 1)
        
        last_pos = file.tell()
        file.seek(last_pos)
        z1 = np.fromfile(file, dtype = 'f4', count = 1)
        
        last_pos = file.tell()
        file.seek(last_pos)
        z2 = np.fromfile(file, dtype = 'f4', count = 1)
        
        last_pos = file.tell()
        file.seek(last_pos)
        ic = np.fromfile(file, dtype = 'i4', count = ns)
        
        last_pos = file.tell()
        for i in range(1, ns+1):
            file.seek(last_pos)
            fxyz[:,:,:,i-1] = np.fromfile(file, dtype = '<f4', count = numberel*numberel*numberel).reshape((numberel,numberel,numberel))
            last_pos = file.tell()
       
        for i in range(1, ns+1):
           
            file.seek(last_pos)
            fxy[:, :, i-1] = np.fromfile(file, dtype='<f4', count = numberel**2).reshape((numberel,numberel))
            last_pos = file.tell()
            
        for i in range(1, ns+1):
            
            file.seek(last_pos)
            fxz[:, :, i-1] = np.fromfile(file, dtype='<f4', count = numberel**2).reshape((numberel,numberel))
            last_pos = file.tell()
        
        for i in range(1, ns+1):
            
            file.seek(last_pos)
            fyz[:, :, i-1] = np.fromfile(file, dtype='<f4', count = numberel**2).reshape((numberel,numberel))
            last_pos = file.tell()
        
        file.seek(last_pos)   
        vxa = np.fromfile(file, dtype = 'f4', count = ns)
        
        last_pos = file.tell()
        file.seek(last_pos)
        vya = np.fromfile(file, dtype = 'f4', count = ns)
        
        last_pos = file.tell()
        file.seek(last_pos)
        vza = np.fromfile(file, dtype = 'f4', count = ns)

    return axes,x1,x2,z1,z2,fxyz,fxy,fxz,fyz

#%% ANALYSIS FUNCTIONS

def Vinterval(data,maskarray,lower,upper):
    data = np.ma.masked_where(maskarray < lower, data)
    data = np.ma.masked_where(maskarray > upper, data)
    return data

def CalcMaskArray(dim, species = 'i', vmaxi = 5, vmaxe = 15, bins = 100, spacing = 'lin'):
    if dim == 1:
        vx = np.zeros(bins+1)
    elif dim == 2:
        vxy = np.zeros((bins+1,bins+1))
    elif dim == 3:
        vxyz = np.zeros((bins+1,bins+1,bins+1))

    if spacing == 'lin':
        if species == 'i':
            vs = np.linspace(-vmaxi,vmaxi, num = bins+1)
        elif species == 'e':
            vs = np.linspace(-vmaxe,vmaxe, num = bins+1)
        for i in range(len(vs)):
            if dim == 1:
                vx[i] = vs[i]**2
            else:
                for j in range(len(vs)):
                    if dim == 2:
                        vxy[i,j] = vs[i]**2 + vs[j]**2  
                    elif dim == 3:
                        for k in range(len(vs)):
                            vxyz[i,j,k] = vs[i]**2 + vs[j]**2 + vs[k]**2
    elif spacing == 'log':
        if species == 'i':
            vs2 = np.logspace(np.log10(0.00001), np.log10(vmaxi**2), num = bins+1)
        elif species == 'e':
            vs2 = np.logspace(np.log10(0.00001), np.log10(vmaxe**2), num = bins+1)
        for i in range(len(vs2)):
            if dim == 1:
                vx[i] = vs2[i]
            else:
                for j in range(len(vs2)):
                    if dim == 2:
                        vxy[i,j] = vs2[i] + vs2[j] 
                    elif dim == 3:
                        for k in range(len(vs2)):
                            
                            vxyz[i,j,k] = vs2[i] + vs2[j] + vs2[k]
    if dim == 1:
        return vx
    elif dim == 2:
        return vxy
    else:
        return vxyz

def CalcBinedges(minv = 0, vmaxe = 20, vmaxi = 7.5, species = 'i', spacing = 'lin', num = 101):
    if spacing == 'lin' and species == 'i':
        return np.linspace(0,(vmaxi**2)*3, num = num)
        #return np.linspace(0,(vmaxi**2)*3, num = num+1)
    elif spacing == 'lin' and species == 'e':
        return np.linspace(0,(vmaxe**2)*3, num = num)
        #return np.linspace(0,(vmaxe**2)*3, num = num+1)
    elif spacing == 'log' and species == 'i':
        return np.logspace(np.log10(0.00001),np.log10(3*vmaxi**2), num = num) 
    elif spacing == 'log' and species == 'e':
        return np.logspace(np.log10(0.00001),np.log10(3*vmaxe**2), num = num)

def CalcdVSumCount(data,maskarray,binedges):
    data_dv_list = []
    data_dv_count_list = []
    data_dv_sum_list = []
    data_dv_cumsum_list = []
    v_bin_mid_list = []
    
    for i in range(len(binedges)-1):    
        huh = Vinterval(data, maskarray, binedges[i], binedges[i+1])
        data_dv_list.append(huh)
        data_dv = data_dv_list[i]
        data_dv_count_list.append(data_dv.count())
        data_dv_sum_list.append(np.sum(data_dv, dtype = 'float64'))
        v_bin_mid = (binedges[i] + binedges[i+1])/2
        v_bin_mid_list.append(v_bin_mid)
        
        if i != 0:
            
            data_dv_cumsum1 = data_dv_cumsum_list[i-1]
            cumsumel = data_dv_cumsum1 + np.sum(data_dv, dtype = 'float64')
            data_dv_cumsum_list.append(cumsumel)
        else:
            data_dv_cumsum_list.append(np.sum(data_dv, dtype = 'float64'))
        
    return np.array(data_dv_sum_list), np.array(data_dv_count_list)#, np.array(data_dv_list), np.array(data_dv_cumsum_list)

def CalcStairs(filenamelist, dir = dir, vmaxi = 7.5, vmaxe = 20, bins = 100, spacing = 'log'):
    bincount_i_list = []
    bincount_e_list = []
    
    vxyz_i = CalcMaskArray(3, species = 'i', vmaxi = vmaxi, bins = bins, spacing = spacing)
    binedges_i = CalcBinedges(species = 'i', vmaxi = vmaxi, num = bins, spacing = spacing)
    vxyz_e = CalcMaskArray(3, species = 'e', vmaxe = vmaxe, bins = bins, spacing = spacing)
    binedges_e = CalcBinedges(species = 'e', vmaxe = vmaxe, num = bins, spacing = spacing)
    
    for file in filenamelist:
        #SUS
        dox = [0,5]
        axes = ReadDist(file, dir, numberel = bins+1)[0]
        fxyz = ReadDist(file, dir, numberel = bins+1)[5]
        vi = axes[:,0]
        ve = axes[:,1]
        fxyz_i = (fxyz[:,:,:,0] + fxyz[:,:,:,2])/((vi[1]-vi[0])**3)
        fxyz_e = (fxyz[:,:,:,1] + fxyz[:,:,:,3])/((ve[1]-ve[0])**3)
        
        sumlist3D_i, countlist3D_i = CalcdVSumCount(fxyz_i, vxyz_i, binedges = binedges_i)
        v_freq_norm_i = sumlist3D_i/countlist3D_i
        bincount_i_list.append(v_freq_norm_i)
        
        sumlist3D_e, countlist3D_e = CalcdVSumCount(fxyz_e, vxyz_e, binedges = binedges_e)
        v_freq_norm_e = sumlist3D_e/countlist3D_e
        bincount_e_list.append(v_freq_norm_e)
        
        
        
    return bincount_i_list, bincount_e_list, binedges_i, binedges_e, vxyz_i, vxyz_e
    
#%% PLOTTING ROUTINE FOR V**2
filelist = [r'\\varyingrun\whole_sim\twpe03000_vmax_i7.dat'
            ,r'\\varyingrun\whole_sim\twpe04000_vmax_i7.dat'
            ,r'\\varyingrun\whole_sim\twpe05000_vmax_i7.dat'
            ,r'\\varyingrun\whole_sim\twpe06000_vmax_i7.dat'
            ,r'\\varyingrun\whole_sim\twpe07000_vmax_i7.dat']

def PlotWholeSimV2(filelist):
    fig,ax = plt.subplots(2,1)
    for file in filelist:
        axes,x1,x2,z1,z2,fxyz,fxy,fxz,fyz = ReadDist(file, ns = 4)
        fxyz_i = fxyz[:,:,:,0]+fxyz[:,:,:,2]
        fxyz_e = fxyz[:,:,:,1]+fxyz[:,:,:,3]
        vi = axes[:,0]
        ve = axes[:,1]
        fxyz_i = (fxyz_i / (vi[1]-vi[0])**3)
        fxyz_e = (fxyz_e / (ve[1]-ve[0])**3)
        
        vxyz_i = CalcMaskArray(3, species = 'i', vmaxi = 7.5)
        binedges_i = CalcBinedges(species = 'i', vmaxi = 7.5)
        sumlist3D_i, countlist3D_i, masked_list, cumsum_i = CalcdVSumCount(fxyz_i, vxyz_i, binedges = binedges_i)
        v_freq_norm_i = sumlist3D_i/countlist3D_i
        
        vxyz_e = CalcMaskArray(3, species = 'e', vmaxe = 20)
        binedges_e = CalcBinedges(species = 'e', vmaxe = 20)
        sumlist3D_e, countlist3D_e, masked_list, cumsum_e = CalcdVSumCount(fxyz_e, vxyz_e, binedges = binedges_e)
        v_freq_norm_e = sumlist3D_e/countlist3D_e
        
        maxi = 62.9177
        maxe = 3.3043
        
        ax[0].stairs(100*(cumsum_i/maxi),binedges_i)
        ax[1].stairs(100*(cumsum_e/maxe),binedges_e)
        
    #ax[0].set_yscale('log')
    #ax[1].set_yscale('log')
    
    ax[0].grid()
    ax[1].grid()
    
    ax[0].set_ylabel('$f(v^2) ions$')
    ax[1].set_ylabel('$f(v^2) electrons$')
    
    ax[1].set_xlabel('$v^2$')
    
    ax[0].legend(['t=3000','t=4000','t=5000','t=6000','t=7000'])
    
    return cumsum_i, cumsum_e

#%% PLOTTING FUNCTIONS

def MakeStairsPlots(bincounts_i,bincounts_e,binedges_i,binedges_e, yscale = 'log', xscale = 'linear'):  
    fig, ax = plt.subplots(2,1)
    
    for j in range(len(bincounts_i)):
        ax[0].stairs(bincounts_i[j], binedges_i)
        ax[1].stairs(bincounts_e[j], binedges_e)
    
    ax[0].set_yscale(yscale)
    ax[1].set_yscale(yscale)
    
    ax[0].set_xscale(xscale)
    ax[1].set_xscale(xscale)
    
    ax[0].legend(['varying', 'baseline'])
    
#%% MAXWELLIAN FITTING ROUTINE FOR VELOCITY DISTRIBUTION

def Maxwellian(x,a,b):
    return (x**2)*a*np.exp(-b*((x-0.5)**2))

def FitMaxwellian(filenamelist, bins = 100, spacing = 'log'):
    y_i, y_e, edge_i, edge_e, vxyzi, vxyze = CalcStairs(filenamelist, dir = r'C:\Users\morri\Desktop\Thesis Data\bin data', spacing = spacing)
    
    edge_i2 = np.insert(edge_i[1:], -1, 0)
    mids_i = (edge_i + edge_i2)/2 #[sum(x) for x in zip(edge_i, edge_i2)]/2
    mids_i = mids_i[:len(mids_i)-1]
    
    y_i = np.nan_to_num(y_i)
    
    first = 1 
    for ys in y_i:
        mids_i = mids_i
        ys = ys*1000
        params_i = curve_fit(Maxwellian, mids_i[:50], ys[:50])[0]
        if first == 1:
            params_i_list = np.array(params_i)
            first = 0 
        else:
            params_i_list = np.vstack((params_i_list, params_i))
    return params_i_list, mids_i, y_i, edge_i

def PlotMaxwellian(a,b,c,d,x,y,edges):
    Maxy = Maxwellian(x,a,b,c,d)
    plt.figure()
    plt.stairs(y,edges)
    plt.plot(x,Maxy)
    plt.yscale('log')
    
def MaxPlotRoutine(filenamelist, spacing = 'log'):
    params_i_list, mids_i, y_i, edges_i = FitMaxwellian(filenamelist, spacing = spacing)
    for i in range(len(y_i)):
        print(params_i_list)
        params_i_list = [params_i_list]
        PlotMaxwellian(params_i_list[i][0],params_i_list[i][1],params_i_list[i][2],params_i_list[i][3],mids_i,y_i[0]*1000, edges_i)

f3000 = [r'\\varyingrun\whole_sim\twpe03000_vmax_i7.dat', r'\\baselinerun\twpe03000_baseline.dat']
f4000 = [r'\\varyingrun\whole_sim\twpe04000_vmax_i7.dat', r'\\baselinerun\twpe04000_baseline.dat']
f5000 = [r'\\varyingrun\whole_sim\twpe05000_vmax_i7.dat', r'\\baselinerun\twpe05000_baseline.dat']
f6000 = [r'\\varyingrun\whole_sim\twpe06000_vmax_i7.dat', r'\\baselinerun\twpe06000_baseline.dat']
f7000 = [r'\\varyingrun\whole_sim\twpe07000_vmax_i7.dat', r'\\baselinerun\twpe07000_baseline.dat']
flists = [f3000,f4000,f5000,f6000,f7000]

def PlotFullSim1(figurename, flists = flists):
    #Whole Sim all times per run
    fig,ax = plt.subplots(2,2)
    for i,species in zip([0,1],['e','i']):
        for j in [1,0]:
            for flist in flists:
                fname = flist[j]
                binc,edges,mids  = Bincountnorm(species, 'lin', fname, dir = r"C:\Users\morri\Desktop\Thesis Data\bin data", numberel = 101)
                ax[j,i].stairs(binc,edges)
            ax[j,i].set_yscale('log')
            ax[j,i].grid()
            if species == 'e':
                ax[j,i].set_xlim(0,1500)
                ax[j,i].set_ylim(0,1)
            elif species == 'i':
                ax[j,i].set_xlim(0,(7.5**2)*3)
                ax[j,i].set_ylim(0,1)
    
    ax[0,1].legend(labels = ['t=30','t=40','t=50','t=60','t=70'], fontsize = 13)
            
    ax[0,0].text(30,0.01, 'varying run', fontsize = 15)
    ax[0,1].text(10,0.05, 'varying run', fontsize = 15)
    ax[1,0].text(30,0.02, 'baseline run', fontsize = 15)
    ax[1,1].text(10,0.05, 'baseline run', fontsize = 15)    
    
    ax[0,0].set_title('Electron distributions', fontsize = 20)
    ax[0,1].set_title('Ion distributions', fontsize = 20)
    
    ax[1,0].set_xlabel('Velocity squared ($v_A^2$)', fontsize = 15)
    ax[1,1].set_xlabel('Velocity squared ($v_A^2$)', fontsize = 15)
    ax[0,0].set_ylabel('$\mathit{f}(v^2)$', fontsize = 15)
    ax[1,0].set_ylabel('$\mathit{f}(v^2)$', fontsize = 15)
    #plt.savefig(r"C:\Users\morri\Desktop\Thesis Figures" + figurename)
    ax[0,0].set_xlim(0,1200)
    ax[1,0].set_xlim(0,1200)
    ax[0,1].set_xlim(0,140)
    ax[1,1].set_xlim(0,140)
    
    ax[0,0].set_ylim(0,1)
    ax[1,0].set_ylim(0,1)
    ax[0,1].set_ylim(0,1)
    ax[1,1].set_ylim(0,1)
    
def PlotFullSim2(timef = [f3000,f5000,f7000]):
    fig,ax = plt.subplots(2,1)
    for it,color in zip(timef, ['#1f77b4','#ff7f0e','#2ca02c']):
        for i,species in zip([0,1], ['e','i']):
            for fname,linest in zip(it,['--','-']):
                binc,edges,mids  = Bincountnorm(species, 'lin', fname, dir = r"C:\Users\morri\Desktop\Thesis Data\bin data", numberel = 101)
                ax[i].stairs(binc,edges,color = color, linestyle = linest)
                ax[i].grid('major')
                ax[i].set_yscale('log')
                ax[i].set_ylabel('$\mathit{f}(v^2)$', fontsize = 15)
    
    ax[1].set_xlabel('Velocity squared ($v_A^2$)', fontsize = 15)
    ax[0].legend(labels = ['varying,t=30','baseline,t=30','varying,t=50','baseline,t=50','varying,t=70','baseline,t=70'], fontsize = 17)
    #ax[0].text(210,0.003,'electrons', fontsize = 15 )
    #ax[1].text(26,0.018,'ions', fontsize = 15)
    ax[0].set_title('Electron distributions', fontsize = 20)
    ax[1].set_title('Ion distributions', fontsize = 20)
    ax[0].set_xlim(0,1200)
    ax[0].set_ylim(0,1)
    ax[1].set_xlim(0,140)
    ax[1].set_ylim(0,1)
    
    
def MakeByvsnoBy(flow):
    datfilelistby, xlist, zlist = CoordstoBoxes('gridcoords90_45_t5000', '10000',
                                                ns = 'north2', mask = 'both', flow = flow)
    datfilelistby = list(map(str, datfilelistby))
    
    #PlotDistboxes(xlist,zlist,dX = 2.386, iteration = '10000', dir = 'D:\\', filename = 'field varying.h5',
    #              title = 'anything else brudda')
    
    datfilelistnoby, xlist, zlist = CoordstoBoxes('gridcoords90_45_t5000', '10000',
                                                  ns = 'north2', mask = 'nonboth', flow = flow)
    datfilelistnoby = list(map(str, datfilelistnoby))
    
    #PlotDistboxes(xlist,zlist,dX = 2.386, iteration = '10000', dir = 'D:\\', filename = 'field varying.h5',
    #              title = 'anything else brudda')
    
    bince = np.array([])
    binci = np.array([])
    t0 = t.time()
    for datfilelist in [datfilelistby,datfilelistnoby]:
        bincesum = np.zeros(120)
        bincisum = np.zeros(120)
        for f in datfilelist:
            t1 = t.time()
            find = datfilelist.index(f)
            f = f.zfill(4)
            f = r'\\' + f + '.dat' 
            
            bince1,edgese,mids  = Bincountnorm('e', 'lin', f, dir = r"D:\Data_Sim_MR\DistsBoxes\t5000\\", numberel = 121, binnum = 121)
            bincesum += bince1
            binci1,edgesi,mids  = Bincountnorm('i', 'lin', f, dir = r"D:\Data_Sim_MR\DistsBoxes\t5000\\", numberel = 121, binnum = 121)
            bincisum += binci1
            t2 = t.time()
            #print(t2-t1) 
            lend = len(datfilelist)
            print(f"{find}/{lend}")
            
        bincesumnorm = bincesum / len(datfilelist)
        bincisumnorm = bincisum / len(datfilelist)
        
        if datfilelist == datfilelistby:
            bince = np.append(bince, bincesumnorm)
            binci = np.append(binci, bincisumnorm)
        elif datfilelist == datfilelistnoby:
            bince = np.vstack((bince,bincesumnorm))
            binci = np.vstack((binci,bincisumnorm))
            
            
    tf = t.time()  
    print(tf-t0)
    return bince, binci, edgese, edgesi

def RefMaxwellian(flow,fitbin,lastbin):
    datfilelistnoby, xlist, zlist = CoordstoBoxes('gridcoords90_45_t5000', '10000',
                                                  ns = 'north2', mask = 'nonboth', flow = flow)
    datfilelistnoby = list(map(str, datfilelistnoby))
    
    #PlotDistboxes(xlist,zlist,dX = 2.386, iteration = '10000', dir = 'D:\\', filename = 'field varying.h5',
    #              title = 'anything else brudda')
    
    bince = np.array([])
    binci = np.array([])
    t0 = t.time()
    for datfilelist in [datfilelistnoby]:
        bincesum = np.zeros(99)
        bincisum = np.zeros(120)
        for f in datfilelist:
            find = datfilelist.index(f)
            t1 = t.time()
            f = f.zfill(4)
            f = r'\\' + f + '.dat' 
            
            bince1,edgese,midse = Bincountnorm('e', 'lin', f, dir = r"D:\Data_Sim_MR\DistsBoxes\t5000\\", numberel = 121, binnum = 100)
            paramse = curve_fit(Maxwellian, 1000*bince1[fitbin:lastbin], midse[fitbin:lastbin])[0]
            fite = Maxwellian(midse[fitbin:lastbin], paramse[0], paramse[1])/1000
            print(paramse)
            
            print(midse[fitbin:][0])
            print(midse[fitbin:lastbin][-1])
            bincesum += bince1
            binci1,edgesi,midsi = Bincountnorm('i', 'lin', f, dir = r"D:\Data_Sim_MR\DistsBoxes\t5000\\", numberel = 121, binnum = 50)
            #bincisum += binci1
            t2 = t.time()
            print(t2-t1)
            lend = len(datfilelist)
            print(f"{find} / {lend}")
        
            #plt.figure()
            #plt.title(f"{f}")
            
            #plt.stairs(bince1,edgese)
            #plt.plot(midse[fitbin:lastbin], fite)
            #plt.stairs(binci1,edgesi)
            
            #plt.yscale('log')
            #plt.xscale('log')
            #plt.show()
        bincesumnorm = bincesum / len(datfilelist)
        #bincisumnorm = bincisum / len(datfilelist)
        
        #bince = np.append(bince, bincesumnorm)
        #binci = np.append(binci, bincisumnorm)
    plt.figure()
    plt.plot(bincesumnorm, edgese)
    tf = t.time()  
    print((tf-t0)/60)
    return bince, binci, edgese, edgesi
    
#%%
def MakeV2Array(vmax,axis,numberel):
    
    dv = (axis[1]-axis[0])

    array_size = (numberel, numberel, numberel)
    
    i, j, k = np.indices(array_size)
    i = (i-((numberel-1)/2))*(dv)
    j = (j-((numberel-1)/2))*(dv)
    k = (k-((numberel-1)/2))*(dv)
    
    result_array = i**2 + j**2 + k**2
    return result_array

def Bincountnorm(species, spacing, filename, dir, numberel = 101, binnum = 101):
    axes,_,_,_,_,fxyz,_,_,_ = ReadDist(filename, dir = dir, numberel = numberel)
    if species == 'e':
        v2_3D = MakeV2Array(np.amax(axes), axes[:,1], numberel)
        fxyz = fxyz[:,:,:,1] + fxyz[:,:,:,3]
    elif species == 'i':
        v2_3D =  MakeV2Array(np.amax(axes), axes[:,0], numberel)
        fxyz = fxyz[:,:,:,0] + fxyz[:,:,:,2]
    
    v2_1D = np.reshape(v2_3D, numberel**3)
    fxyz_1D = np.reshape(fxyz, numberel**3)
    binedgesv2 = CalcBinedges(species = species, spacing = spacing, num = binnum)
    bincount = np.histogram(v2_1D, binedgesv2, weights = fxyz_1D)[0]
    
    binedgesv = np.sqrt(binedgesv2)
    binvmid = binedgesv[1:] - 0.5 * (binedgesv[1:] - binedgesv[:-1])
    bindv = np.ediff1d(binedgesv)
    binvol = binvmid**2 * bindv
    

    bincountnorm = bincount/binvol
    return bincountnorm, binedgesv2, binvmid**2

def PullingThrough(filename, dir, numberel):
    axes,_,_,_,_,fxyz,_,_,_ = ReadDist(filename, dir = dir, numberel = numberel)

    v2e_3D = MakeV2Array(np.amax(axes), axes[:,1], numberel)
    v2i_3D = MakeV2Array(np.amax(axes[0][:]), axes[:,0], numberel)
    fxyz_e = fxyz[:,:,:,1] + fxyz[:,:,:,3]
    fxyz_i = fxyz[:,:,:,0] + fxyz[:,:,:,2]
    
    fxyze_1D = np.reshape(fxyz_i, numberel**3)
    fxyzi_1D = np.reshape(fxyz_e, numberel**3)
    v2e_1D = np.reshape(v2e_3D, numberel**3)
    v2i_1D = np.reshape(v2i_3D, numberel**3)
        
    logedgesv2_e = CalcBinedges(species = 'e', spacing = 'log')
    logedgesv2_i = CalcBinedges(species = 'i', spacing = 'log')
    linedgesv2_e = CalcBinedges(species = 'e', spacing = 'lin')
    linedgesv2_i = CalcBinedges(species = 'i', spacing = 'lin')
    
    logedgesv_e = np.sqrt(logedgesv2_e)
    logedgesv_i = np.sqrt(logedgesv2_i)
    linedgesv_e = np.sqrt(linedgesv2_e)
    linedgesv_i = np.sqrt(linedgesv2_i)
    
    vmidlog_e = logedgesv_e[1:] - 0.5 * (logedgesv_e[1:] - logedgesv_e[:-1])
    vmidlog_i = logedgesv_i[1:] - 0.5 * (logedgesv_i[1:] - logedgesv_i[:-1])
    vmidlin_e = linedgesv_e[1:] - 0.5 * (linedgesv_e[1:] - linedgesv_e[:-1])
    vmidlin_i = linedgesv_i[1:] - 0.5 * (linedgesv_i[1:] - linedgesv_i[:-1])
    
    dvlog_e = np.ediff1d(logedgesv_e)
    dvlog_i = np.ediff1d(logedgesv_i)
    dvlin_e = np.ediff1d(linedgesv_e)
    dvlin_i = np.ediff1d(linedgesv_i)
    
    vollog_e = vmidlog_e**2 * dvlog_e
    vollog_i = vmidlog_i**2 * dvlog_i
    vollin_e = vmidlin_e**2 * dvlin_e
    vollin_i = vmidlin_i**2 * dvlin_i
    
    binc_lin_e = np.histogram(v2e_1D, linedgesv2_e, weights = fxyze_1D)[0]
    bincnorm_lin_e = binc_lin_e/vollin_e
    
    binc_log_e = np.histogram(v2e_1D, logedgesv2_e, weights = fxyze_1D)[0]
    bincnorm_log_e = binc_log_e/vollog_e
    print(bincnorm_log_e)
    params_log_e = curve_fit(Maxwellian, vmidlog_e[59:], bincnorm_log_e[59:])[0]
    params_lin_e = curve_fit(Maxwellian, vmidlin_e[5:], bincnorm_lin_e[5:])[0]
    
    Maxlog = Maxwellian(vmidlog_e, params_log_e[0], params_log_e[1], params_log_e[2], params_log_e[3])
    Maxlin = Maxwellian(vmidlin_e, params_log_e[0], params_lin_e[1], params_lin_e[2], params_lin_e[3])
    
    plt.figure()
    plt.stairs(bincnorm_log_e,logedgesv2_e)
    #plt.stairs(bincnorm_lin_e,linedgesv2_e)
    
    plt.plot(vmidlog_e[:], Maxlog[:])
    #plt.plot(vmidlin_e, Maxlin)
    
    plt.yscale('log')
    plt.xscale('log')
    print(np.sum(binc_lin_e))
    print(np.sum(fxyz_e))
    

    