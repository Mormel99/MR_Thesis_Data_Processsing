#%% IMPORTING 
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 11:28:27 2023

@author: morri
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time as t
#%% DICTIONARIES

units = {'space': '($c/\u03C9_{pi}$) ', 'density': '[$n_{0}$)]', 'v': '[$v_{A}$]', 
         'B': '($B_{0}$)', 'E': '($E_{0}$)', 'rho': '[$en_{0}$]', 'J': '[$v_{A}en_{0}$]'}

#%% READING DATA + NORMALIZING
def GetData(datakey, iteration,
            dir = r'C:\Users\morri\Desktop\Thesis Data\h5 data\\', file = 'fields.h5'):
    iteration = iteration.zfill(10)
    
    with h5py.File(dir + file, 'r') as file:
        if datakey == 'nnx' or datakey == 'nnx' or datakey == 'nnz' or datakey == 'teti' or datakey == 'wpewce' or datakey == 'xe' or datakey == 'xmax' or datakey == 'ze' or datakey == 'zmax':
            sim_info = file['simulaton_information']
            return sim_info[datakey][0][0]
        else:
            data = file['data'][iteration]
            wpewce = file['simulation_information']['wpewce'][0][0]
            mime = 25
            #fix mime plz
            
            if datakey == 'ex' or datakey == 'ey' or datakey == 'ez':
                return np.array(data[datakey]) * wpewce**2 * np.sqrt(mime)
            
            elif datakey == 'bx' or datakey == 'by' or datakey == 'bz':
                return np.array(data[datakey]) * wpewce
            
            elif datakey == 'A':
                return np.array(data[datakey])
            elif datakey == 'vxsi' or datakey == 'vysi' or datakey == 'vzsi' or datakey == 'vxxi' or datakey == 'vyyi' or datakey == 'vzzi' or datakey == 'vxse' or datakey == 'vyse' or datakey == 'vzse' or datakey == 'vxxe' or datakey == 'vyye' or datakey == 'vzze' or datakey == 'dnsi' or datakey == 'dnse':
                dfacn = data['dns']['2'].attrs['dfac'][0]
                dfacs = data['dns']['3'].attrs['dfac'][0]
                
                if datakey == 'dnsi':
                    return np.array(data['dns']['1'])*dfacn + np.array(data['dns']['3'])*dfacs
                elif datakey == 'dnse':
                    return np.array(data['dns']['2'])*dfacn + np.array(data['dns']['4'])*dfacs
                elif datakey == 'vxsi' or datakey == 'vysi' or datakey == 'vzsi':
                    datakey = datakey.replace('i','')
                    return (np.array(data[datakey]['1'])*dfacn + np.array(data[datakey]['3'])*dfacs) * wpewce * np.sqrt(mime)
                elif datakey == 'vxse' or datakey == 'vyse' or datakey == 'vzse':
                    datakey = datakey.replace('e','')
                    return (np.array(data[datakey]['2'])*dfacn + np.array(data[datakey]['4'])*dfacs) * wpewce * np.sqrt(mime)
                elif datakey == 'vxxi' or datakey == 'vyyi' or datakey == 'vzzi':
                    datakey = datakey.replace('i','')
                    return (np.array(data[datakey]['1'])*dfacn + np.array(data[datakey]['3'])*dfacs) * (wpewce**2) * mime
                elif datakey == 'vxxe' or datakey == 'vyye' or datakey == 'vzze':
                    datakey = datakey.replace('e','')
                    return (np.array(data[datakey]['2'])*dfacn + np.array(data[datakey]['4'])*dfacs) * (wpewce**2) * mime
            elif datakey == 'Axline' or datakey == 'RE' or datakey == 'UB' or datakey == 'UK_electrons' or datakey == 'UK_ions' or datakey == 'UT_electrons' or datakey == 'UT_ions' or datakey == 'dt' or datakey == 'time' or datakey == 'xline_position':
                    return data.attrs[datakey]
                
def SubsetBy(data, iteration, threshold = 0.005):
             by = GetData('by', iteration)
             data = np.ma.masked_where(by < threshold, data)
             #data = np.ma.masked_outside(data, 0.2, 0.7)
             return data
         
'''
LIJST VAN EENHEDEN EN DAN DYNAMISCH MAKEN MET DE FUCKING LABELS VAN DE ASSEN MAAR HUH
'''
#%% PLOTTING FUNCTIONS

#Makes multiple subplots, given the datakeys and iteration of the quantities 
#you want to plot
def MakeSubplots(datakeys, iteration): 
    A = GetData('A', iteration)
    
    dist = 265/5 
    x = np.linspace(0, 4 * dist, num = 3200)
    z = np.linspace(-dist, dist, num = 3200)
    n_figs = len(datakeys)
    fig,ax = plt.subplots(n_figs, 1)
    if n_figs == 1: 
        
        cmp = ax.pcolormesh(x, z, GetData(datakeys[0], iteration))
        #ax.contour(x, z, A, levels = 25, colors = 'black', linewidths = 0.5)
        ax.set_xlabel('X ' + units['space'])
        ax.set_ylabel('Z ' + units['space'])
        fig.colorbar(cmp, ax = ax)
    else:
        ctopbound = 0
        cbotbound = 0
        
        for i in range(n_figs):
            #PLZ FIX THIS EVENTUALLY APPEND IS DUMBASS SLOW I THINK ATLEAST
            data = np.array([])
            data = np.append(data, GetData(datakeys[i], iteration))
            datmax = np.amax(data[i])
            datmin = np.amin(data[i])
            if datmax > ctopbound:
                ctopbound = datmax
            if datmin < cbotbound:
                cbotbound = datmin
        
        for i in range(n_figs):   
            cmp = ax[i].pcolormesh(x, z, data[i], vmin = cbotbound, vmax = ctopbound)
            ax[i].contour(x, z, A, levels = 25, colors = 'black', linewidths = 0.5)
            if i == n_figs -1:
                ax[i].set_xlabel('X ' + units['space'])
            ax[i].set_ylabel('Z ' + units['space'])
            fig.colorbar(cmp, ax = ax[i])
            
#Without GetData, quicker for small adjustments of plots
def MakeSubplots2(data, iteration, vmax = 10, xlim = 'none', ylim = 'none', dir = 'D:\\', file = 'testfields2.h5', unit = 'insert unit', title ='', cbarl = None, cbarh = None): 
    A = GetData('A', iteration, dir = dir, file = file)
    
    dist = 265/5 
    x = np.linspace(0, 4 * dist, num = 3200)
    z = np.linspace(-dist, dist, num = 3200)
    n_figs = len(data)
    fig,ax = plt.subplots(n_figs, 1)
    if n_figs == 1: 
        
        cmp = ax.pcolormesh(x, z, data[0], vmin = -0.5, vmax = 0.5, cmap = 'bwr')
        ax.contour(x, z, A, levels = 10, colors = 'black', linewidths = 0.5)
        ax.set_xlabel('X ' + units['space'])
        ax.set_ylabel('Z ' + units['space'])
        #fig.colorbar(cmp, ax = ax)
    else:
        ctopbound = 0
        cbotbound = 0
        
        for i in range(n_figs):
            datmax = np.amax(data[i])
            datmin = np.amin(data[i])
            
            if datmax > ctopbound:
                ctopbound = datmax
            if datmin < cbotbound:
                cbotbound = datmin
            
        for i in range(n_figs):   
            cmp = ax[i].pcolormesh(x, z, data[i], vmin = cbotbound, vmax = ctopbound)
            ax[i].contour(x, z, A, levels = 25, colors = 'black', linewidths = 0.5)
            
            if i == n_figs -1:
                ax[i].set_xlabel('X ' + units['space'])
            ax[i].set_ylabel('Z ' + units['space'])
            
            
            
            #fig.colorbar(cmp, ax = ax[i])
        fig.colorbar(cmp, ax = ax.ravel().tolist()).set_label(units['J'])
    if xlim != 'none':
        plt.xlim(xlim[0], xlim[1])
    if ylim != 'none':
        plt.ylim(ylim[0], ylim[1])

    fig.suptitle(title)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes('right', size = '5%', pad = 1)
    fig.colorbar(cmp, ax = ax).set_label(unit)

    #plt.gca().set_aspect(aspect = 0.75)
#Makes cut through the z-axis at a specified x coordinate. 
#So coordinate wise it shows (z,q), q being the plotted quantity.
def MakeZcut(datakeys, x, iteration):
    
    dist = 265/5 
    #x = np.linspace(0, 4 * dist, num = 3200)
    z = np.linspace(-dist, dist, num = 3200)
    n_figs = len(datakeys)
    fig,ax = plt.subplots(n_figs, 1)
    if n_figs == 1: 
        
        y = GetData(datakeys[0], iteration = iteration)[:,x]
        ax.plot(z, y)
        ax.set_xlabel('Z ' + units['space'])
        ax.set_ylabel('insert correct unit here')
        
    else:
        for i in range(n_figs):   
            y = GetData(datakeys[i], iteration = iteration)[:,x]
            ax[i].plot(z, y)
            #DOES SOMETHING WEIRD WITH DASHED AND NON DASHED DISPLAY
            #ALSO FOR DIFFERENT PLOTS DIFFERENT DEFINITIONS OF MAJOR/MINOR
            #CAUSE OF DIFF DATA RANGES
            ax[i].set_yticks([-1, -0.5, 0, 0.5, 1], minor = 'True')
            ax[i].grid(which = 'major', c = 'grey')
            ax[i].grid(which = 'minor', c = 'grey', linestyle = '--')
            
            if i == n_figs -1:
                ax[i].set_xlabel('Z ' + units['space'])
            ax[i].set_ylabel('insert correct uni here')
            
def MakeZcut2(data, x):
    
    dist = 265/5 
    #x = np.linspace(0, 4 * dist, num = 3200)
    z = np.linspace(-dist, dist, num = 3200)
    n_figs = len(data)
    fig,ax = plt.subplots(n_figs, 1)
    if n_figs == 1: 
        
        y = data[0][:,x]
        ax.plot(z, y)
        ax.set_xlabel('Z ' + units['space'])
        ax.set_ylabel('insert correct unit here')
        
    else:
        for i in range(n_figs):   
            y = data[i][:,x]
            ax[i].plot(z, y)
    
            if i == n_figs -1:
                ax[i].set_xlabel('Z ' + units['space'])
            ax[i].set_ylabel('insert correct unit here')
            
def MakeXcut(datakeys, z, iteration):
    
    dist = 265/5 
    x = np.linspace(0, 4 * dist, num = 3200)
    #z = np.linspace(-dist, dist, num = 3200)
    n_figs = len(datakeys)
    fig,ax = plt.subplots(n_figs, 1)
    if n_figs == 1: 
        
        y = GetData(datakeys[0], iteration = iteration)[z,:]
        ax.plot(x, y)
        ax.set_xlabel('X ' + units['space'])
        ax.set_ylabel('insert correct unit here')
        
    else:
        for i in range(n_figs):   
            y = GetData(datakeys[i], iteration = iteration)[z,:]
            ax[i].plot(x, y)
    
            if i == n_figs -1:
                ax[i].set_xlabel('X ' + units['space'])
            ax[i].set_ylabel('insert correct unit here')
            
def MakeXcut2(data, z, iteration):
    
    dist = 265/5 
    x = np.linspace(0, 4 * dist, num = 3200)
    #z = np.linspace(-dist, dist, num = 3200)
    n_figs = len(data)
    fig,ax = plt.subplots(n_figs, 1)
    if n_figs == 1: 
        
        y = GetData(data[0])[z,:]
        ax.plot(x, y)
        ax.set_xlabel('hoi')
        ax.set_ylabel('hai')
        
    else:
        for i in range(n_figs):   
            y = GetData(data[i])[z,:]
            ax[i].plot(x, y)
    
            if i == n_figs -1:
                ax[i].set_xlabel('X ' + units['space'])
            ax[i].set_ylabel('insert correct unit here')
            
def PlotDistboxes(xcoor, zcoor, dX, iteration = '8000', dir = 'smth', filename = 'smthn', title = 'placeholder'):
    dim = 265/5
    x = np.linspace(0,4*dim, num = 3200)
    z = np.linspace(-dim,dim, num = 3200)
    
    fig, ax = plt.subplots()
    ax.set_xlim(0, dim*4)
    ax.set_ylim(-dim,dim)
    ax.set_aspect('equal')
    
    for i in range(len(xcoor)):
        ax.add_patch(plt.Rectangle((xcoor[i]-dX/2,zcoor[i]-dX/2), dX, dX, fill=False, color='black'))
    by = GetData('by', iteration = iteration, file = filename, dir = dir)
    A = GetData('A', iteration = iteration, file = filename, dir = dir)
    cmp = ax.pcolormesh(x,z,by)
    ax.contour(x, z, A, levels = 20, colors = 'black', linewidths = 0.5) #levels = [37.85,40.29]
    plt.xlabel('X'+units['space'], fontsize = 13)
    plt.ylabel('Z'+units['space'], fontsize = 13)
    ax.set_ylim(-30,30)
    ax.title.set_text(title) #fontsize = 13)
    cb = plt.colorbar(cmp, fraction = 0.05, aspect = 10)
    cb.set_label( label = '$B_{y}$($B_{0}$)', fontsize = 14)
    
def PlotMacroUDivision(sim = 'baseline', figurename = 'anytin'):
    
    if sim == 'baseline':
        rangemax = 16000
        dirlist = [r'C:\Users\morri\Desktop\Thesis Data\h5 data\\']
        filelist = ['fields_baseline.h5']
    elif sim == 'varying':
        rangemax = 16000
        dir = [r'D:\\']
        filelist = ['field varying.h5']
    elif sim == 'both':
        rangemax = 16000
        dirlist = [r'C:\Users\morri\Desktop\Thesis Data\h5 data\\',r'D:\\']
        filelist = ['fields_baseline.h5','field varying.h5']
        
    
    for dir,file,linestyle in zip(dirlist,filelist,['-','--']):
        U_matrix = np.array([])
        UBlist = []
        UK_elist = []
        UK_ilist = []
        UT_elist = []
        UT_ilist = []
        
        itlist = range(0,rangemax+400, 400)
        for iteration in itlist:
            if iteration == 0:
                iteration = 1
            iteration = str(iteration)
            UB =   GetData('UB', iteration, dir = dir, file = file)[0]
            UK_e = GetData('UK_electrons', iteration, dir = dir, file = file)[0]
            UK_i = GetData('UK_ions', iteration, dir = dir, file = file)[0]
            UT_e = GetData('UT_electrons', iteration, dir = dir, file = file)[0]
            UT_i = GetData('UT_ions', iteration, dir = dir, file = file)[0]
            #U_mat_el = np.array([[UB, UK_e, UK_i, UT_e, UT_i]])
            #U_matrix = np.concatenate((U_matrix, U_mat_el), axis = 0)
            UBlist.append(UB)
            UK_elist.append(UK_e)
            UK_ilist.append(UK_i)
            UT_elist.append(UT_e)
            UT_ilist.append(UT_i)
        U_matrix = np.array([UBlist,UK_elist,UK_ilist,UT_elist,UT_ilist])
        U_sum = np.array([sum(x) for x in zip(*U_matrix)])
        tlist = np.array(itlist) / (4*25) # wpewce * mime
        #for i in range(len(U_matrix)):
        #    U_matrix[i] = (U_matrix[i]/U_sum)*100
    
    #'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
        plt.plot(tlist,U_matrix[0]-U_matrix[0][0], color = '#1f77b4',  linestyle = linestyle)
        plt.plot(tlist,U_matrix[1]-U_matrix[1][0], color = '#ff7f0e',  linestyle = linestyle)
        plt.plot(tlist,U_matrix[2]-U_matrix[2][0], color = '#2ca02c',  linestyle = linestyle)
        plt.plot(tlist,U_matrix[3]-U_matrix[3][0], color = '#d62728',        linestyle = linestyle)
        plt.plot(tlist,U_matrix[4]-U_matrix[4][0], color = '#9467bd',    linestyle = linestyle)
        plt.plot(tlist,U_sum - U_sum[0],           color = '#8c564b',    linestyle = linestyle)
    
    plt.title(f'Temporal change in energy division \n of {sim} runs')
    plt.xlabel('Time ($\u03A9_{ci}^{-1}$)')
    plt.ylabel('Change in energy density ($AED$)')
    
    legend1 = plt.legend(labels = ['$\u0394 u_{b}$','$\u0394 u_{k,e}$','$\u0394 u_{k,i}$','$\u0394 u_{t,e}$','$\u0394 u_{t,i}$','$\u0394 u_{sum}$'])

    plt.grid('major')
    plt.xlim(0,160)
    
    line1 = plt.Line2D([], [], color='black', linestyle='-', label='baseline run')
    line2 = plt.Line2D([], [], color='black', linestyle='--', label='varying run')
    
    # Add the second legend
    legend2 = plt.legend(handles=[line1, line2], loc='lower left', fontsize='medium')
    
    # Add the second legend to the current axes
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)
    
    plt.show()  
    
    return U_matrix


def PlotDeltaUb(sim = 'both', figurename = 'anytin'):
    
    if sim == 'both':
        rangemax = 16000
        dirlist = [r'C:\Users\morri\Desktop\Thesis Data\h5 data\\',r'D:\\']
        filelist = ['fields_baseline.h5','field varying.h5']
        
        itlist = range(0,rangemax+400, 400)
        plt.figure()
        for dir,file in zip(dirlist,filelist):
            UBlist = []

            for iteration in itlist:
                if iteration == 0:
                    iteration = 1
                iteration = str(iteration)
                UB =   GetData('UB', iteration, dir = dir, file = file)[0]
                UBlist.append(UB)
            tlist = np.array(itlist) / (4*25) # wpewce * mime
        
            plt.plot(tlist,UBlist-UBlist[0])
        plt.grid('major')
        plt.yticks(ticks = [-1200000,-1100000,-1000000,-900000,-800000,-700000,-600000,-500000,-400000,-300000,-200000,-100000,0])
        plt.grid('both', axis = 'y')
        
        plt.title('Temporal change in magnetic energy density \n of both runs')
        plt.legend(labels= ['baseline run', 'varying run'])
        plt.ylabel('$\u0394 u_{b}$')
        plt.xlabel('Time ($\u03A9_{ci}^{-1}$)')
        plt.ylabel('$\u0394 u_{b}$ ($AED$)')
        plt.savefig(r"C:\Users\morri\Desktop\Thesis Figures" + figurename)
    else:
        print('Why tho?')
        
def PlotUbUtUk(sim = 'both', figurename = 'anytin'):
    check = 1
    if sim == 'both':
        rangemax = 16000
        dirlist = [r'C:\Users\morri\Desktop\Thesis Data\h5 data\\',r'D:\\']
        filelist = ['fields_baseline.h5','field varying.h5']
        
        itlist = range(0,rangemax+400, 400)
        plt.figure()
        for dir,file in zip(dirlist,filelist):
            U_matrix = np.array([])
            UBlist = []
            UK_list = []
            UT_list = []
            
            itlist = range(0,rangemax+400, 400)
            for iteration in itlist:
                if iteration == 0:
                    iteration = 1
                iteration = str(iteration)
                UB =   GetData('UB', iteration, dir = dir, file = file)[0]
                UK_e = GetData('UK_electrons', iteration, dir = dir, file = file)[0]
                UK_i = GetData('UK_ions', iteration, dir = dir, file = file)[0]
                UT_e = GetData('UT_electrons', iteration, dir = dir, file = file)[0]
                UT_i = GetData('UT_ions', iteration, dir = dir, file = file)[0]
                #U_mat_el = np.array([[UB, UK_e, UK_i, UT_e, UT_i]])
                #U_matrix = np.concatenate((U_matrix, U_mat_el), axis = 0)
                UBlist.append(UB)
                UK_list.append(UK_e + UK_i)
                UT_list.append(UT_e + UT_i)
            U_matrix = np.array([UBlist,UK_list,UT_list])
            U_sum = np.array([sum(x) for x in zip(*U_matrix)])
            tlist = np.array(itlist) / (4*25) # wpewce * mime
            #for i in range(len(U_matrix)):
            #    U_matrix[i] = (U_matrix[i]/U_sum)*100
            
            if check == True:
                plt.plot(tlist,U_matrix[0]-U_matrix[0][0], label = 'baseline $\u0394 u_{b}$')
                plt.plot(tlist,U_matrix[1]-U_matrix[1][0], label = 'baseline $\u0394 u_{k}$')
                plt.plot(tlist,U_matrix[2]-U_matrix[2][0], label = 'baseline $\u0394 u_{t}$')
                plt.plot(tlist,U_sum - U_sum[0], 'black', label = 'baseline $\u0394 u_{sum}$')
                check = 0
            else: 
                plt.plot(tlist,U_matrix[0]-U_matrix[0][0], color = '#1f77b4', linestyle = '--', label = 'var $\u0394 u_{b}$')
                plt.plot(tlist,U_matrix[1]-U_matrix[1][0], color = '#ff7f0e', linestyle = '--', label = 'var $\u0394 u_{k}$')
                plt.plot(tlist,U_matrix[2]-U_matrix[2][0], color = '#2ca02c',linestyle = '--', label = 'var $\u0394 u_{t}$')
                plt.plot(tlist,U_sum - U_sum[0], 'black', linestyle = '--', label = 'var $\u0394 u_{sum}$')
    
    plt.title('Temportal changes in energy division \n of both runs')
    plt.xlabel('Time ($\u03A9_{ci}^{-1}$)')
    plt.ylabel('Change in energy density ($AED$)')
    plt.xlim(0,160)
    
    legend1 = plt.legend(labels = ['$\u0394 u_{b}$', '$\u0394 u_{k}$', '$\u0394 u_{t}$','$\u0394 u_{sum}$'], loc = 'upper left')
    
    line1 = plt.Line2D([], [], color='black', linestyle='-', label='baseline run')
    line2 = plt.Line2D([], [], color='black', linestyle='--', label='varying run')
    
    # Add the second legend
    legend2 = plt.legend(handles=[line1, line2], loc='lower left', fontsize='medium')
    
    # Add the second legend to the current axes
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)
    
    plt.grid('major')
    plt.show()  
    plt.savefig(r"C:\Users\morri\Desktop\Thesis Figures" + figurename)
    return U_matrix

def PlotERec(sim = 'baseline', figurename = 'anytin'):
    if sim == 'baseline':
        rangemax = 20000
        dir = r'C:\Users\morri\Desktop\Thesis Data\h5 data\\'
        file = 'fields_baseline.h5'
    elif sim == 'varying':
        rangemax = 16000
        dir = r'D:\\'
        file = 'field varying.h5'
    elif sim == 'both':
        rangemax = 16000
        fig,ax = plt.subplots(2,1)
        dirlist = [r'C:\Users\morri\Desktop\Thesis Data\h5 data\\', r'D:\\']
        filelist = ['fields_baseline.h5', 'field varying.h5']
        
        
        itlist = range(0,rangemax+400, 400)
        for dir,file in zip(dirlist,filelist):
            UBlist = []

            for iteration in itlist:
                if iteration == 0:
                    iteration = 1
                iteration = str(iteration)
                UB =   GetData('UB', iteration, dir = dir, file = file)[0]
                UBlist.append(UB)
            tlist = np.array(itlist) / (4*25) # wpewce * mime
        
            ax[0].plot(tlist,UBlist-UBlist[0])
        ax[0].set_yticks(ticks = [-1200000,-1100000,-1000000,-900000,-800000,-700000,-600000,-500000,-400000,-300000,-200000,-100000,0])
        ax[0].grid('both', axis = 'y')
        
        ax[0].set_title('Temporal change in magnetic energy density \n of both runs')
        ax[0].legend(labels= ['baseline run', 'varying run'])
        ax[0].set_ylabel('$\u0394 u_{b}$ ($AED$)')
        
        
    for dir, file in zip(dirlist, filelist):
        ERlist = np.array([])
        
        itlist = range(0,rangemax+400, 400)
        for iteration in itlist:
            if iteration == 0:
                iteration = 1
            iteration = str(iteration)
            ER = GetData('RE', iteration, dir = dir, file = file)[0]
            ERlist = np.append(ERlist, ER)
        
        
        ax[1].plot(np.array(itlist)/(4*25), ERlist)
    
    ax[1].set_title('Reconnection rate')
    ax[0].grid('both')
    ax[1].grid('both')
    ax[1].legend(labels = ['baseline run', 'varying run'])
    
    ax[1].set_xlabel('Time ($\u03A9_{ci}^{-1}$)')
    ax[1].set_ylabel(r'$ \frac{\partial}{\partial t} \phi_{rec}$')
    
    ax[0].set_xlim(4.5,160)
    ax[1].set_xlim(4.5,160)
    ax[1].set_ylim(0,0.25)
    plt.show()
    plt.savefig(r"C:\Users\morri\Desktop\Thesis Figures" + figurename)
    
#%% CALCULATING QUANTITIES FROM READ DATA

def CalcE(iteration, dir = 'something', file = 'something'):
    E = np.array([GetData('ex', iteration, dir = dir, file = file),GetData('ey', iteration, dir = dir, file = file), GetData('ez', iteration, dir = dir, file = file)])
    return E,np.linalg.norm(E, axis = 0)

def CalcB(iteration):
    B = np.array([GetData('bx', iteration),GetData('by', iteration), GetData('bz', iteration)])
    return B,np.linalg.norm(B, axis = 0)

def Calc_u(iteration):
    return 1/2 * CalcB[1]**2 #+ CalcE[1]**2

def CalcS(iteration):
    Sx = np.multiply(CalcE(iteration)[0][1],CalcB(iteration)[0][2]) - np.multiply(CalcE(iteration)[0][2],CalcB(iteration)[0][1])
    Sy = np.multiply(CalcE(iteration)[0][2],CalcB(iteration)[0][0]) - np.multiply(CalcE(iteration)[0][0],CalcB(iteration)[0][2])
    Sz = np.multiply(CalcE(iteration)[0][0],CalcB(iteration)[0][1]) - np.multiply(CalcE(iteration)[0][1],CalcB(iteration)[0][0])
    S = np.array([Sx, Sy, Sz])
    return S, np.linalg.norm(S, axis = 0)

def CalcJ(iteration, species = 'both', dir = 'something', file = 'something2'):
    if species == 'both' or species == 'ion':
        vxsi = GetData('vxsi', iteration, dir = dir, file = file)
        vysi = GetData('vysi', iteration, dir = dir, file = file)
        vzsi = GetData('vzsi', iteration, dir = dir, file = file)
        
        Vsi = np.array([vxsi, vysi, vzsi])
        Ji = Vsi 
    if species == 'both' or species == 'electron':
        vxse = GetData('vxse', iteration, dir = dir, file = file)
        vyse = GetData('vyse', iteration, dir = dir, file = file)
        vzse = GetData('vzse', iteration, dir = dir, file = file)
     
        Vse = np.array([vxse, vyse, vzse])
        Je = -1*Vse
    if species == 'ion':
        return Ji
    elif species == 'electron':
        return Je
    else:
        return Ji + Je
    
def JdotE(iteration, dir = 'something', file = 'something', species = 'both'):
    return np.sum(np.multiply(CalcJ(iteration, dir = dir, file = file, species = species),CalcE(iteration, dir = dir, file = file)[0]), axis = 0)

def CalcPT(iteration, species = 'both'):
    if species == 'both' or species == 'ion':
        vxsi = GetData('vxsi', iteration)
        vysi = GetData('vysi', iteration)
        vzsi = GetData('vzsi', iteration) 
        dnsi = GetData('dnsi', iteration)
        Vbi = np.array([vxsi, vysi, vzsi])/dnsi
        
        vxxi = GetData('vxxi', iteration)
        vyyi = GetData('vyyi', iteration)
        vzzi = GetData('vzzi', iteration)
        pdiai = np.array([vxxi, vyyi, vzzi])
        
        for i in range(3):
            Vbi[i][dnsi < 0.005] = 0
            pdiai[i] = pdiai[i] - Vbi[i]**2 * dnsi 
        
        pi = (1/3)*np.sum(pdiai, axis=0)
        Ti = (pdiai / dnsi)
        for i in range(3):
            Ti[i][dnsi < 0.005] = 0 #hmhmhmmh 2 for loops
        
    if species == 'both' or species == 'electron':
        vxse = GetData('vxse', iteration)
        vyse = GetData('vyse', iteration)
        vzse = GetData('vzse', iteration)
        dnse = GetData('dnse', iteration)
        Vbe = np.array([vxse, vyse, vzse])/dnse
        
        vxxe = GetData('vxxe', iteration)
        vyye = GetData('vyye', iteration)
        vzze = GetData('vzze', iteration)
        pdiae = np.array([vxxe, vyye, vzze])
        
        for i in range(3):
            Vbi[i][dnsi < 0.005] = 0
            pdiae[i] = (1/25)*(pdiae[i] - Vbe[i]**2 * dnse)
        
        pe = (1/3)*np.sum(pdiae, axis=0)
        Te = (pdiae / dnse)
        for i in range(3):
            Te[i][dnse < 0.00005] = 0 #hmhmhmmh 2 for loops
        
    #the pdia's are the diagonal components of the p tensor shifted to the bulk frame
    if species == 'both':
        pdia = pdiai + pdiae
        p = pi + pe
        T = Ti + Te
        return pdia, p, T
    elif species == 'ion': 
        return pdiai, pi, Ti
    else: 
        return pdiae, pe, Te
    
def CalcTs(iteration, species = 'both'):
    T = CalcPT(iteration, species = species)[2]
    B,b = CalcB(iteration)
    t = np.linalg.norm(T, axis = 0)
    tpar = np.sum(np.multiply(T,B), axis = 0) / b
    
    tper = t - tpar
    return T, t, tpar, tper

#%% DISTRIBUTION BOXES FUNCTIONS

def GridDists(xdens = 20, zdens = 10):
    dim = 265/5
    X = np.linspace(0,4*dim, num = 3200)
    Z = np.linspace(-dim, dim, num = 3200)
    
    dDx = 3200/xdens
    dDz = 3200/zdens
    Xlist = []
    Zlist = []
    with open('gridcoords' + str(xdens) + '_' + str(zdens), 'w') as file:
        for i in range(xdens-1):
            for j in range(zdens-1):
                #xbotleft = X[int(i*dDx)]
                #xbotright = X[int((i+1)*dDx)]
                #xtopleft = xbotleft
                #xtopright = xbotright
                
                #zbotleft = Z[int(j*dDz)]
                #zbotright = zbotleft
                #ztopleft = Z[int((j+1)*dDz)]
                #ztopright = ztopleft
            
                #file.write(f"({xtopleft},{ztopleft}) ({xtopright},{ztopright})\n({xbotleft},{zbotleft}) ({xbotright},{zbotright})\n\n")
                
                xlow = X[int(i*(dDx))]
                xlow = np.round(xlow, decimals = 3)
                xhigh = X[int((i+1)*dDx)]
                xhigh = np.round(xhigh, decimals = 3)
                zlow = Z[int(j*dDz)]
                zlow = np.round(zlow, decimals = 3)
                zhigh = Z[int((j+1)*dDz)]
                zhigh = np.round(zhigh, decimals = 3)
                
                Xc = (xlow + xhigh)/2
                Zc = (zlow + zhigh)/2
                if Zc < 19 and Zc > -19 and Xc > 30 and Xc < 107:
                    if Xc > 63 and (Zc < -16 or Zc > 16):
                        pass     
                    else:    
                        file.write(f"{xlow} {xhigh} {zlow} {zhigh}\n")
                        Xlist.append(Xc)
                        Zlist.append(Zc)
        print(len(Xlist))
        
        ds = xhigh - xlow
        return Xlist, Zlist, ds
                
def SelectFieldline(iteration = '8000', levels = 10):
    A = GetData('A', iteration)
    by = GetData('by', iteration)
    dim = 265/5
    X = np.linspace(0,4*dim, num = 3200)
    Z = np.linspace(-dim, dim, num = 3200)
    fig,ax = plt.subplots()
    cunt = ax.contour(X, Z, A, levels = levels, colors = 'black', linewidths = 0.5)
    
    coords = []
    for i in range(1,len(cunt.collections)-1):    
       for j in [0,1]:
        p = cunt.collections[i].get_paths()
        v = p[j].vertices
        x = v[:,0]
        y = v[:,1]
        coords.append([x,y])
    rightline = 'no'
    while rightline == 'no':
        line = input('What contour line do you want distribution boxes for?')
        plt.figure()
        for i in range(len(coords)):    
            if i == int(line):
                plt.plot(coords[i][0],coords[i][1], c = 'red')
            else:
                plt.plot(coords[i][0],coords[i][1], c = 'black')
        plt.pcolormesh(X,Z,by)
        plt.show()
        rightline = input('Is that the line you want to make dist. boxes for?')

    return coords, line

def FieldlineDists(coords, line, num):
    xc = coords[int(line)][0]
    yc = coords[int(line)][1]
    dx = xc[+1:]-xc[:-1]
    dy = yc[+1:]-yc[:-1]    
    ds = np.array((0,*np.sqrt(dx*dx + dy*dy)))
    s = np.cumsum(ds)
    print()
    xinter =  np.interp(np.linspace(0, s[-1], num), s, xc)
    yinter =  np.interp(np.linspace(0, s[-1], num), s, yc)
    
    '''
    plt.figure()
    plt.scatter(xinter, yinter, c = 'green', s = 2)
    for i in range(len(coords)):    
        if i == int(line):
            plt.plot(coords[i][0],coords[i][1], c = 'red')
        else:
            plt.plot(coords[i][0],coords[i][1], c = 'black')
    plt.show()
    '''
    
    dsinter = np.square((xinter[1] - xinter[0])**2 + (yinter[1]-yinter[2])**2) 
    dX = dsinter/4
    dZ = dX
    
    with open('fieldcoords' + str(line) + '_' + str(num), 'w') as file:
        for i in range(len(xinter)):
            for j in range(len(yinter)):
                #xbotleft = xinter[i] - dX
                #xbotright = xinter[i] + dX
                #xtopleft = xbotleft
                #xtopright = xbotright
                
                #zbotleft = yinter[j] - dZ
                #zbotright = zbotleft
                #ztopleft = yinter[j] + dZ
                #ztopright = ztopleft
    
                #file.write(f"({xtopleft},{ztopleft}) ({xtopright},{ztopright})\n({xbotleft},{zbotleft}) ({xbotright},{zbotright})\n\n")
                
                xlow = xinter[i] - dX
                xhigh = xinter[i] + dX
                zlow = yinter[i] - dZ
                zhigh = yinter[i] + dZ
                file.write(f"({xlow},{xhigh},{zlow},{zhigh})\n")
    PlotDistboxes(xinter, yinter, dX)

def ByDists(iteration, dS = 1, threshold = 0.10, dir = r'C:\Users\morri\Desktop\Thesis Code\dist boxes\\'):
    dim = 265/5
    X = np.linspace(0,4*dim, num = 3200)
    Z = np.linspace(-dim, dim, num = 3200)
    dD = int(3200/ ((265/5)*4) * (dS*2))
    XX,ZZ = np.meshgrid(X,Z)
    XXby = SubsetBy(XX, iteration, threshold = threshold)
    ZZby = SubsetBy(ZZ, iteration, threshold = threshold)
    
    Xlist = np.array([])
    Zlist = np.array([])
    thresholdstr = str(threshold).replace('.','_')
    counter = 0
    with open(f'bycoords_{iteration}_{thresholdstr}', 'w') as file:
        for i in range(int(3200/dD)):
            i = i*dD
            for j in range(int(3200/(np.floor(dD/2)))):
                j = int(j*(np.floor(dD/2)))
                X = XXby[i,j]
                Xlist = np.append(Xlist, X)
                Z = ZZby[i,j]
                Zlist = np.append(Zlist, Z)
                
                xlow = X - dS/2
                xhigh = X + dS/2
                zlow = Z - dS/2
                zhigh = Z + dS/2
                if xlow < 0 or xlow > 0:
                    file.write(f"({xlow},{xhigh},{zlow},{zhigh})\n")
                    counter += 1
    file.close()
    PlotDistboxes(Xlist, Zlist, dS)
    print(counter)

def SubsetA(data, iteration, lowthreshold = 37.8501, highthreshold = 40.2875, threshold = 0):
    A = GetData('A', iteration)
    if threshold != 0:    
        data = np.ma.masked_where(A < threshold, data)
    else:
        data = np.ma.masked_where(A<lowthreshold, data)
        data = np.ma.masked_where(A>highthreshold, data)
    return data
    
def ADists(iteration, dS = 1, lowthreshold = 37.8501, highthreshold = 40.2875, threshold = 0, dir = r'C:\Users\morri\Desktop\Thesis Code\dist boxes\\'):
    dim = 265/5
    X = np.linspace(0,4*dim, num = 3200)
    Z = np.linspace(-dim, dim, num = 3200)
    dD = int(3200/ ((265/5)*4) * (dS*2))
    XX,ZZ = np.meshgrid(X,Z)
    XXA = SubsetA(XX, iteration, lowthreshold = lowthreshold, highthreshold = highthreshold, threshold = threshold)
    ZZA = SubsetA(ZZ, iteration, lowthreshold = lowthreshold, highthreshold = highthreshold, threshold = threshold)
    
    
    Xlist = np.array([])
    Zlist = np.array([])
    counter = 0
    #thresholdstr = str(threshold).replace('.','_')
    with open(f'Acoords_{iteration}', 'w') as file:
        for i in range(int(3200/dD)):
            i = i*dD
            for j in range(int(3200/(np.floor(dD/2)))):
                j = int(j*(np.floor(dD/2)))
                X = XXA[i,j]
                Xlist = np.append(Xlist, X)
                Z = ZZA[i,j]
                Zlist = np.append(Zlist, Z)
                
                xlow = X - dS/2
                xhigh = X + dS/2
                zlow = Z - dS/2
                zhigh = Z + dS/2
                
                if xlow < 0 or xlow > 0:
                    counter += 1
                    file.write(f"({xlow},{xhigh},{zlow},{zhigh})\n")
    file.close()
    PlotDistboxes(Xlist, Zlist, dS)
    print(counter)
    
#%% From coordsfile to correct distboxes
def CoordstoBoxes(filename, iteration, Aliml = 37.8501, Alimh = 40.2875, dir = r'C:\Users\morri\Desktop\Thesis Code\\',
                  dir2 = 'D:\\', filename2 = 'field varying.h5', ns = 'both', mask = 'by2', flow = 'both'):
    Axline = GetData('Axline', iteration = iteration, dir = r'D:\\', file = 'field varying.h5')
    with open(filename, 'r') as file:
        Coordsarray = np.loadtxt(file, delimiter = '/n', dtype = 'str')
        A = GetData('A', iteration, dir = dir2, file = filename2)
        dim = 265/5
        X = np.around(np.linspace(0,4*dim, num = 3200), decimals = 3)
        Z = np.around(np.linspace(-dim,dim, num = 3200), decimals = 3)
        datfilelist = []
        centerlistx = []
        centerlistz = []
        counter = 1
        for coords in Coordsarray:
            fltcoords = list(map(float, coords.split()))
            xli = np.squeeze(np.where(X == fltcoords[0])[0])
            xhi = np.squeeze(np.where(X == fltcoords[1])[0])
            zli = np.squeeze(np.where(Z == fltcoords[2])[0])
            zhi = np.squeeze(np.where(Z == fltcoords[3])[0])       
            
            xc = round((fltcoords[0]+fltcoords[1])/2, 3)
            zc = round((fltcoords[2]+fltcoords[3])/2, 3)
            
            xcrel = (dim*4)/xc
            xci = round(3200/xcrel)
            zcrel = (dim*2)/(zc+dim)
            zci = round(3200/zcrel)
            
            Atl = A[xli,zhi]
            Atr = A[xhi,zhi]
            Abl = A[xli,zli]
            Abr = A[xhi,zli]
            Ac = A[zci,xci]
            
            #and (Ac < 43.9454 - 0.12 or Ac > 46.3935 + 0.12)
            #if (Atl > Aliml and Atl < Alimh) and (Atr > Aliml and Atr < Alimh) and (Abl > Aliml and Abl < Alimh) and (Abr > Aliml and Abr < Alimh):
            Aliml2 = 37.8501
            Alimh2 = 40.2875
            Aliml1 = 43.9454
            Alimh1 = 46.3935
            if mask == 'by2':
                Aliml = Aliml2
                Alimh = Alimh2
            elif mask == 'by1':
                Aliml = Aliml1
                Alimh = Alimh1

            if flow == 'inflow' and Ac < Axline:
                
                if mask == 'both' and ((Ac > Aliml2-0.09 and Ac < Alimh2+0.09) or (Ac > Aliml1 - 0.09 and Ac < Alimh1 + 0.09)):
                    if ns == 'north2' and (Ac > Aliml2-0.09 and Ac < Alimh2+0.09) and zc > 0:
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                    elif ns == 'north2' and (Ac > Aliml1-0.09 and Ac < Alimh1+0.09):
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                    elif ns == 'both':
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                elif mask == 'nonboth' and (Ac < Aliml2-0.2 or Ac > Alimh2+0.2) and (Ac < Aliml1-0.2 or Ac > Alimh1+0.2):
                    datfilelist.append(counter)
                    centerlistx.append(xc)
                    centerlistz.append(zc)
                        
                elif (mask == 'by1' or mask == 'by2') and (Ac > Aliml-0.09 and Ac < Alimh+0.09):
                    if ns == 'north' and zc > 0:
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                    elif ns == 'both':
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                        
            elif flow == 'outflow' and Ac > Axline:
                
                if mask == 'both' and ((Ac > Aliml2-0.09 and Ac < Alimh2+0.09) or (Ac > Aliml1 - 0.09 and Ac < Alimh1 + 0.09)):
                    if ns == 'north2' and (Ac > Aliml2-0.09 and Ac < Alimh2+0.09) and zc > 0:
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                    elif ns == 'north2' and (Ac > Aliml1-0.09 and Ac < Alimh1+0.09):
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                    elif ns == 'both':
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                elif mask == 'nonboth' and (Ac < Aliml2-0.2 or Ac > Alimh2+0.2) and (Ac < Aliml1-0.2 or Ac > Alimh1+0.2):
                    datfilelist.append(counter)
                    centerlistx.append(xc)
                    centerlistz.append(zc)
                        
                elif (mask == 'by1' or mask == 'by2') and (Ac > Aliml-0.09 and Ac < Alimh+0.09):
                    if ns == 'north' and zc > 0:
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                    elif ns == 'both':
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
            elif flow == 'both':
                    
                if mask == 'both' and ((Ac > Aliml2-0.09 and Ac < Alimh2+0.09) or (Ac > Aliml1 - 0.09 and Ac < Alimh1 + 0.09)):
                    if ns == 'north2' and (Ac > Aliml2-0.09 and Ac < Alimh2+0.09) and zc > 0:
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                    elif ns == 'north2' and (Ac > Aliml1-0.09 and Ac < Alimh1+0.09):
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                    elif ns == 'both':
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                elif mask == 'nonboth' and (Ac < Aliml2-0.2 or Ac > Alimh2+0.2) and (Ac < Aliml1-0.2 or Ac > Alimh1+0.2):
                    datfilelist.append(counter)
                    centerlistx.append(xc)
                    centerlistz.append(zc)
                       
                elif (mask == 'by1' or mask == 'by2') and (Ac > Aliml-0.09 and Ac < Alimh+0.09):
                    if ns == 'north' and zc > 0:
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)
                    elif ns == 'both':
                        datfilelist.append(counter)
                        centerlistx.append(xc)
                        centerlistz.append(zc)     
            counter += 1
            
        file.close()
    return datfilelist, centerlistx, centerlistz
       
        