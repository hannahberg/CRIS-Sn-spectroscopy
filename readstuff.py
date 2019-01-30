import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import glob
import re
import seaborn as sns
from scipy.constants import pi, speed_of_light, hbar, Boltzmann
from scipy import signal

# Read in the necessary files
def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]"""
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    #Sort the given list in the way that humans expect
    l.sort(key=alphanum_key)
    return l

def name_dat_file(folder):
    filenames = glob.glob('/home/hannahcb/cern/ode-solver/' + folder +'/' + '*.*')
    filenames = sort_nicely(filenames)
    return filenames

# Put all the data into dataframes
def read_to_df(filename):
    colnames = ['observed wavelength','unc','ritz wavelength','r_unc','rel. int.','Aki','Acc','Ei','l_configuration','l_term','l_J','u_configuration','u_term','u_J','type','tp ref.','line ref.']
    df = pd.read_csv(filename,skiprows=6,names=colnames,delimiter='\s*\|',engine='python',comment='#',index_col=False,)
    df[['Ei','Ek']] = df['Ei'].str.split('\s*-',expand=True)
    df = df.dropna(thresh=2).reset_index(drop=True) # drops empty rows
    return df

def read_to_df2(filename):
    colnames = ['observed wavelength','unc','ritz wavelength','r_unc','t. wavenumber','rel. int.','Aki','Acc','Ei','l_configuration','l_term','l_J','u_configuration','u_term','u_J','type','tp ref.','line ref.']
    df = pd.read_csv(filename,skiprows=6,names=colnames,delimiter='\s*\|',engine='python',comment='#',index_col=False,)
    df[['Ei','Ek']] = df['Ei'].str.split('\s*-',expand=True)
    df = df.dropna(thresh=4).reset_index(drop=True) # drops empty rows
    df[['Aki','Ei','Ek']] = df[['Aki','Ei','Ek']].apply(pd.to_numeric,errors='coerce')
    return df

def read_to_df3(filename):
    colnames = ['observed wavelength','ritz wavelength','rel. int.','Aki','Acc','Ei','l_configuration','l_term','l_J','u_configuration','u_term','u_J','type','tp ref.','line ref.']
    df = pd.read_csv(filename,skiprows=6,names=colnames,delimiter='\s*\|',engine='python',comment='#',index_col=False,)
    df[['Ei','Ek']] = df['Ei'].str.split('\s*-',expand=True)
    df = df.dropna(thresh=2).reset_index(drop=True) # drops empty rows
    return df

# read in initial or final cross section and normalize to get population

def read_population(filename,colnames=['level','xs']):
    pop = pd.read_csv(filename,names=colnames,delimiter=';',float_precision='high').apply(pd.to_numeric,errors='coerce')
    pop = renorm(pop).sort_values(by='level').reset_index(drop=True)
    return pop
# Functions to normalize arrays and dataframes
def normd(array):
    suma = np.sum(array)
    if suma == 0:
        return array
    else:
        newarray = array/suma
        return newarray
# test normd()
def test_norm():
    ta = np.random.rand(1,50)
    tanew = normd(ta)
    if (np.sum(tanew) - 1) != 0:
        print(np.sum(tanew))
        raise ValueError

def renorm(df):
    sum_xs = sum(df['xs'])
    df['pop'] = df['xs']/sum_xs
    return df

# Make df out of Einstein coefficients and add energy indexes
def get_ki(df_A, p):
    indexlist = []
#     for i, (index, row) in enumerate(df.iterrows())
    for i,(ai,arow) in enumerate(df_A.iterrows()):
        ei = round(arow.Ei,2)
        ek = round(arow.Ek,2)
        iindex = p.index[p['level'].between(ei-0.1,ei+0.1)].tolist()[0]
        kindex =  p.index[p['level'].between(ek-0.1,ek+0.1)].tolist()[0]
        indexlist.append([kindex,iindex])
    return indexlist#pd.DataFrame(indexes, columns=list('ki'), dtype=int)

def einsteindf(df,populationdf):
    emc2 = df[['Ek','Ei','Aki']].sort_values(by='Ek').reset_index(drop=True).apply(pd.to_numeric)
    indexes = get_ki(emc2,populationdf)
    aindexes = pd.DataFrame(indexes, columns=list('ki'), dtype=int)
    adf = pd.concat([emc2,aindexes],axis=1)
    return adf

# Short program for transmission without pumping

def shortA(df_A,initialpop,stop,steps):
    time = np.linspace(0,stop,steps)
    dt = -(time[1]-time[0])
    population = initialpop['pop'].copy()
    for t in range(len(time)):
        for i, (ai, arow) in enumerate(df_A.iterrows()):
            pk = population[int(arow.k)]
            if pk > 0:
                diff = pk*(1-np.exp(dt*arow.Aki))
                population[int(arow['k'])] -= diff
                population[int(arow['i'])] += diff
    return population

# Make an on/off laser array
def laser_bool(time,steps,pumpstart,pumpduration,pumpoff):
    dt = time/steps
    '''if pumpduration.any() < dt:
        raise ValueError
    elif pumpoff < dt:
        raise ValueError'''
    print(steps)
#     time = np.linspace(0,t_stop,steps)
#     dt = time[1]-time[0]
    tstrt = np.tile([False],int(pumpstart/dt))
    print(int((pumpstart)/dt))
    ton = np.tile([True],int(pumpduration/dt))
    print(int(pumpduration/dt))
    toff = np.tile([False],int(round(pumpoff/dt,0)))
#     print(toff)
    onoff = np.append(ton,toff)
#     print(ton/toff)
    print(steps/len(onoff))
    onoffl = np.tile(onoff,int(steps/len(onoff))+1)
    laserb = np.append(tstrt,onoffl)
    return laserb[0:steps]
