# -*- coding: utf-8 -*-
"""
Created on Tue Apr 02 18:48:28 2019

@author: HannahChoi
"""

# In[]:
# This code takes the "final hierarchy" found from self-consistency and iteration from C-C and T-C connections.
# Then, based on this hierarchy scores of all cortical and thalamic regions, assign +/-1 (FF/FB) for each C->T pairs.

# In[]:
from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# In[]:

CreConf = 1    # 1 if using CC hierarhcy without Cre-confidence; 0 if not

# In[]: cortico-thalamic source-line-target pairs
if CreConf == 1:
    input_dir1 =r'./Input_Conf/'
    output_dir = r'./results_Conf/results/'     
    dft2=pd.read_excel(input_dir1+"TC_CCconf_iter.xls")
elif CreConf == 0:
    input_dir1 = r'./Input_NoConf/'        
    output_dir = r'./results_NoConf/results/'
    dft2=pd.read_excel(input_dir1+"TC_CCnoconf_iter.xls")
    
dft1=pd.read_excel(input_dir1+"L5 v L6 NPV values.xlsx")
dfV=dft1[(dft1.target != "VISC")&(dft1.source != "SSp-un")&(dft1.target != "SSp-un")&(dft1.source != "VISC")]

# In[]: cortical & thalamic hierarchy from thalamo-cortical & cortico-cortical connections

dfh=dft2[["CortexThalamus","areas",20]]
dfh=dfh.rename(columns={20: "h"})

# In[]:
dfV['FFB']=""

for i_pair in dfV.index:
    h_source = np.array(dfh[dfh.areas == dfV['source'][i_pair]].h)
    h_target = np.array(dfh[dfh.areas == dfV['target'][i_pair]].h)
    if len(h_source)>0 and len(h_target)>0:
        if h_source>h_target:
            dfV['FFB'][i_pair]=-1
        elif h_source<h_target:
            dfV['FFB'][i_pair]=1
        else:
            dfV['FFB'][i_pair]=0

# In[]:

dfV_FB=dfV[dfV['FFB']==-1]
dfV_FF=dfV[dfV['FFB']==1]
dfV_LAT=dfV[dfV['FFB']==0]            


# In[]:

f = plt.figure()
plt.plot(dfV_FB['LOG_NPV_L5'],dfV_FB['LOG_NPV_L6'],'ro',label='FB')
plt.plot(dfV_FF['LOG_NPV_L5'],dfV_FF['LOG_NPV_L6'],'bo',label='FF')
#plt.plot(dfV_LAT['LOG_NPV_L5'],dfV_LAT['LOG_NPV_L6'],'ko',label='Lat')
plt.xlabel('Log NPV L5')
plt.ylabel('Log NPV L6')
plt.legend()
#plt.title('')
plt.show()
f.savefig(output_dir+"LOG_CT_sourcelayer_FFB.pdf", bbox_inches='tight')



dfV=dfV.rename(columns={"FFB": "FFB_from CC+TC"})
dfV.to_excel(output_dir+'CT_sourcelayer_FFB.xls')


# In[]:

f = plt.figure()
plt.plot(10**np.array(dfV_FB['LOG_NPV_L5']),10**np.array(dfV_FB['LOG_NPV_L6']),'ro',label='FB')
plt.plot(10**np.array(dfV_FF['LOG_NPV_L5']),10**np.array(dfV_FF['LOG_NPV_L6']),'bo',label='FF')
plt.xlabel('NPV L5')
plt.ylabel('NPV L6')
plt.legend()
#plt.title('')
plt.show()
f.savefig(output_dir+"NoLOG_CT_sourcelayer_FFB.pdf", bbox_inches='tight')

