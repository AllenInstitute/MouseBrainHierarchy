from __future__ import division
from functools import partial
import logging
import numpy as np
import pandas as pd


"""
Created on Fri Nov 16 17:45:36 2018

@author: Hannah Choi
"""
"""
This code contains functions needed for mapping TC clusters to FF/FB by 
maximizing the global hierarchy score of the shuffled thalamo-cortical connectivity data, 
and is called in "run_TC_shuffled.py"

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


def assign_hierarchy_levels(df, i_shuffle,  j):
    """Returns global hierarchy scores with (hgc)/without (hg) TC confidence 
    (= equal # of FF/FB bias)"""
    def ff_or_fb(cluster):
        b = bin(c0 + j)[-(cluster)]
        return 2*int(b)-1

    n = len(df.clu.unique())
    c0 = 2**n

    df.loc[:, "ffb"] = df.clu.apply(ff_or_fb)

    confg = confidence(df)

    ##################################################################
    '''Change the file name accordingly'''
    df_cortex=pd.read_excel(r'./Output/shuffled/CCshuffled_conf_iter'+str(i_shuffle)+'.xls')
    # df_cortex=pd.read_excel(r'./Output/shuffled/CCshuffled_noconf_iter'+str(i_shuffle)+'.xls')
    ##################################################################

    df = df.join(df_cortex.set_index('areas'), on='target')    
    ##################################################################
    '''Take the last hierarchy scores of cortical regions at the end of iterations (10 steps)'''
    levels_t = df.groupby('target')[10].mean()
    ##################################################################
    hr_t = df.target.map(levels_t)
    hrc_t = hr_t
    hr_s = hierarchy_values(df)
    hrc_s= hr_s 


    hg = np.mean(df.ffb*(hr_t - hr_s))
    hgc = np.mean(df.ffb*(hrc_t - hrc_s))*confg

    # Signed
    # hg = np.mean(df.ffb*np.sign(hr_t - hr_s))
    
    logging.debug("{:0{n}b}  (hg, hgc) = ({:.3f}, {:.3f})".format(j, hg, hgc, n=n))

    return hg, hgc 


def fit(df, i_shuffle, parallel=True, n_procs=-1):   
    """iterates through assign_hierarchy_levels with parallel support"""
    def generator():
        n_possible = 2**(len(df.clu.unique()) - 1)
        return range(n_possible)

    hierarchy_vals = []
    func = partial(assign_hierarchy_levels, df, i_shuffle)
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
