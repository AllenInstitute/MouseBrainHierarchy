from __future__ import division
from functools import partial
import logging

import numpy as np
import pandas as pd

"""
Created on Mon Aug 20 20:30:31 2018

@author: Hannah Choi
"""
"""
This code contains functions needed for mapping TC clusters to FF/FB by 
maximizing the global hierarchy score of the given thalamo-cortical connectivity data, 
and is called in "run_TC.py" and "run_TC_shuffled.py"

"""

def confidence(df):
     """Returns multiplier which prevents all thalamic regions placed below cortical regions, by biasing towards equal # of FF and FB connections"""   
     count_ff = len(df[df.ffb==1])
     count_fb = len(df[df.ffb==-1])
     confnew = min(count_ff, count_fb)/(count_ff+count_fb)
 
     return confnew


def hierarchy_values(df, val='ffb'):
    """Returns hierarchy values for talamic regions as source"""
    source = df.groupby('source')[val].mean()
    levels = -source

    return df.source.map(levels) 


def assign_hierarchy_levels(df, j):
    """Returns global hierarchy scores with (hgc)/without (hg) TC confidence 
    (= equal # of FF/FB bias)"""
    def ff_or_fb(cluster):
        b = bin(c0 + j)[-(cluster)]
        return 2*int(b)-1

    n = len(df.clu.unique())
    #n = len(df.clu.unique()-1)
    c0 = 2**n

    df.loc[:, "ffb"] = df.clu.apply(ff_or_fb)
    confg = confidence(df)

    ##################################################################
    '''Change the file name accordingly'''
    df_cortex=pd.read_excel(r'./Output/CC_conf_iter.xls') 
    # df_cortex=pd.read_excel(r'./Output/CC_noconf_iter.xls') 
    ##################################################################
    
    df = df.join(df_cortex.set_index('areas'), on='target')    
    ##################################################################
    '''Take the last hierarchy scores of cortical regions at the end of iterations (20 steps)'''
    levels_t = df.groupby('target')[20].mean()
    ##################################################################
    hr_t = df.target.map(levels_t)
    hrc_t = hr_t
    hr_s = hierarchy_values(df)
    hrc_s= hr_s 


    hg = np.mean(df.ffb*(hr_t - hr_s))
    hgc = np.mean(df.ffb*(hrc_t - hrc_s))*confg
    
    logging.debug("{:0{n}b}  (hg, hgc) = ({:.3f}, {:.3f})".format(j, hg, hgc, n=n))

    return hg, hgc


def fit(df, parallel=True, n_procs=-1):
    """iterates through assign_hierarchy_levels with parallel support"""
    def generator():
        n_possible = 2**(len(df.clu.unique()) - 1)
        return range(n_possible)

    hierarchy_vals = []
    func = partial(assign_hierarchy_levels, df)
    if parallel:
        from multiprocessing import Pool, cpu_count
        if n_procs < 1:
            n_procs = cpu_count()
        pool = Pool(processes=n_procs)

        for result in pool.imap(func, generator()):
            hierarchy_vals.append(result)
    else:
        for result in map(func, generator()):
            hierarchy_vals.append(result)

    return hierarchy_vals
