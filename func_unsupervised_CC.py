from __future__ import division
from functools import partial
import logging
import numpy as np

"""
Created on Mon Aug 20 20:30:31 2018

@author: Hannah Choi
"""
"""
This code contains functions needed for mapping CC clusters to FF/FB by 
maximizing the global hierarchy score of the given cortico-cortical connectivity data, 
and is called in "run_CC.py" and "run_CC_shuffled.py"

"""

def CC_confidence(df):
  """Returns confidence of cre lines"""
  func = lambda x: 1 - np.abs(x.mean())

  return df.groupby('creline')['ffb'].transform(func)

def hierarchy_values(df, val='ffb'):
    """Returns hierarchy values for sources and targets"""
    source = df.groupby('source')[val].mean()
    target = df.groupby('target')[val].mean()
    levels = 0.5*(target - source)

    return df.source.map(levels), df.target.map(levels)


def hierarchy_values_with_conf(df):
    """Returns hierarchy values for sources and targets when weighted by Cre confidence"""
    df.loc[:, "conf_ffb"] = df.ffb * df.conf

    return hierarchy_values(df, val='conf_ffb')


def assign_hierarchy_levels(df, j):
    """Assigns a consistency val"""
    def ff_or_fb(cluster):
        b = bin(c0 + j)[-(cluster)]
        return 2*int(b)-1

    n = len(df.clu.unique())
    c0 = 2**n

    df.loc[:, "ffb"] = df.clu.apply(ff_or_fb)
    df.loc[:, "conf"] = CC_confidence(df)
    confg = df.conf.mean() + 1e-12 

    hr_s, hr_t = hierarchy_values(df)
    hrc_s, hrc_t = hierarchy_values_with_conf(df)


    hg = np.mean(df.ffb*(hr_t - hr_s))
    hgc = np.mean(df.ffb*df.conf*(hrc_t - hrc_s))  / (confg**2)  
    
    # Signed: if computing the global hierarchy score by taking "sign", use this.
    #hg = np.mean(df.ffb*np.sign(hr_t - hr_s))
    #hgc = np.mean(df.ffb*df.conf*np.sign(hrc_t - hrc_s))  / (confg)  
    
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
