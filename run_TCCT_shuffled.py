from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IterativeMethod import iterativeTCCT

"""
Created on Sat Nov 17 19:48:14 2018

@author: Hannah Choi
"""
"""
This code genereates shuffled versions of the CT connectitivity data and 
computes the global hierarchy scores of the shuffled data.
"""


# In[]: Set input and output directories 

input_dir = r'./Input/'        # Directory with the file "CC_TC_CT_clusters.xlsx" & "CT_sourcelayer_FFB.xls"
input_dir2 = r'./Output/'               # Directory with the file  "TC_CCconf_iter.xls" or "TC_CCnoconf_iter.xls", "ghc_TC.xls"
output_dir = r'./Output/shuffled/'   # Directory to save the ouputs from the shuffled experimental data

''' ATTENTION! Change the "df_cortex" accordingly in func_unsupervised_TC as well! '''
CreConf = 1                 # 1 if using CC hierarhcy with Cre-confidence; 0 if not

# In[]: Read the excel file with source-target-creline pairs and their cluster numbers. 
    
xls=pd.ExcelFile(input_dir+"CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'CC+TC clusters with MD avged')
df=df[(df.hemi == "ipsi")&(df.creline != "C57BL/6J / Emx1")&(df.target != "VISC")&(df.source != "SSp-un")&(df.target != "SSp-un")&(df.source != "VISC")]

# In[]: Cortical and thalamic hierarchy from CC+TC iteration

if CreConf == 0:
    h_CC_TC=pd.read_excel(input_dir2+"TC_CCnoconf_iter.xls")
elif CreConf == 1:
    h_CC_TC=pd.read_excel(input_dir2+"TC_CCconf_iter.xls")
    
C_areas = h_CC_TC[h_CC_TC["CortexThalamus"]=="C"]["areas"].unique()
T_areas = h_CC_TC[h_CC_TC["CortexThalamus"]=="T"]["areas"].unique()

# In[]: 2 clusters (FF vs FB) based on LDA of thalamo-cortical source-line-target pairs (dfCT)

dfCT = pd.read_excel(input_dir+"CT_sourcelayer_FFB.xls") 
dfCT = dfCT[['source','target','FFB_LDA']]
dfCT = dfCT.rename(columns={"FFB_LDA":"ffb"})


# In[]: 9 clusters of thalamo-cortical source-line-target pairs (dfTC)

dfTC=pd.read_excel(input_dir2+"inputexpanded_TC9.xls")
dfTC=dfTC[['source','target','ffb_c']]
dfTC = dfTC.rename(columns={"ffb_c":"ffb"})

# In[]: Merge dataframes of CT and TC connections 

dfVT0 = pd.concat([dfTC, dfCT], ignore_index=True)

source_areas = dfVT0["source"].unique()
target_areas = dfVT0["target"].unique()

num_TC = np.shape(dfTC)[0]
num_CT = np.shape(dfCT)[0]


# In[ ]: Find global hierarchy scores of shuffled TC+CT connectivity data

n_iter = 10
by_creline = 0
n_shuffle = 100

hr_ct_shuffled = np.zeros(n_shuffle)
hr_iter_c_shuffled = np.zeros(n_shuffle)
hr_iter_ct_shuffled = np.zeros(n_shuffle)


for i_shuffle in range(0,n_shuffle):
    
    print('i_shuffle='+str(i_shuffle))
    
    dfVT_shuffled = dfVT0
    
    if by_creline == 0:
        source_list= dfVT_shuffled.source
        target_list= dfVT_shuffled.target
        source_shuffled = source_list.sample(frac=1).reset_index(drop=True)
        target_shuffled = target_list.sample(frac=1).reset_index(drop=True)
        source_shuffled.index = source_list.index
        target_shuffled.index = target_list.index 
        dfVT_shuffled.loc[source_list.index,"source"]=np.array(source_shuffled)
        dfVT_shuffled.loc[target_list.index,"target"]=np.array(target_shuffled)
    
    
    ###########################################################################
    '''Find initial hierarchy score of each thalamic area'''   
    
    n_T_areas=len(T_areas) 
    
    hrs=range(0,n_T_areas)
    hrt=range(0,n_T_areas)
    hrc=range(0,n_T_areas)
    
    for i in range(0,n_T_areas):
        hrs[i]=-np.mean(dfVT0[dfVT0.source == T_areas[i]].ffb)
        if len(dfVT0[dfVT0.target == T_areas[i]]) != 0:
            hrt[i]=np.mean(dfVT0[dfVT0.target == T_areas[i]].ffb)
        else:
            hrt[i]=0
            
        hrc[i]=(hrs[i]+hrt[i])/2
    
    data=[T_areas,hrc]
    data=np.transpose(data)
    columns = ['areas','h']
    dfiT = pd.DataFrame(data,columns=columns)
    ###########################################################################

    ###########################################################################
    '''Iterate thalamic + cortical hierarhcy scores'''
    
    if CreConf == 1:
        dfiC = pd.read_excel(output_dir+'CCshuffled_conf_iter'+str(i_shuffle)+'.xls') 
    elif CreConf == 0:
        dfiC = pd.read_excel(output_dir+'CCshuffled_noconf_iter'+str(i_shuffle)+'.xls') 
    
    dfiC['h'] = dfiC[5]
    dfVC = pd.read_excel(output_dir+'inputexpanded_CC_shuffled'+str(i_shuffle)+'.xls')
    dfVC = dfVC[["source","target","ffb_nc","ffb_c"]]
    if CreConf == 0:
        dfVC["ffb"] = dfVC["ffb_nc"]
    elif CreConf == 1:
        dfVC["ffb"] = dfVC["ffb_c"]
    
    
    dfiT = dfiT[["areas","h"]]
    dfiC = dfiC[["areas","h"]]
    dfVT = dfVT0[["source","target","ffb"]]
    dfVC = dfVC[["source","target","ffb"]]
    
    hr_iter = iterativeTCCT(dfiC, dfVC, dfiT, dfVT, n_iter)
    
    
    iteration=np.arange(0,n_iter+1,1)
    n_area=np.shape(hr_iter)[0]
    allareas = hr_iter["areas"].unique()
    
    for i_area in range(0,n_area):
        if hr_iter['areas'][i_area] in list(dfiC['areas']):
            hr_iter.loc[i_area,'CortexThalamus'] = 'C'
        else:
            hr_iter.loc[i_area,'CortexThalamus'] = 'T'
    
    hr_iter = hr_iter[['areas','CortexThalamus', 0,n_iter] ]  
   
    #    if CreConf == 1:
    #        hr_iter.to_excel(output_dir+'TCCT_CCconf_iter'+str(i_shuffle)+'.xls')
    #    elif CreConf == 0:
    #        hr_iter.to_excel(output_dir+'TCCT_CCnoconf_iter'+str(i_shuffle)+'.xls')
    ###########################################################################


    ###########################################################################
    '''global hierarchy score of the shuffled data before & after iteration'''
    
    dfi_TCCT = hr_iter[["CortexThalamus","areas",0,n_iter]]
    dfi_TCCT = dfi_TCCT.rename(columns={0: "h0", n_iter:"h_iter"})
    
    dfV_CC = dfVC[['source','target','ffb']]
    dfV_TCCT = dfVT[["source","target","ffb"]]
    
    
    dfi_cortex1 = dfi_TCCT[(dfi_TCCT.CortexThalamus == 'C')]
    dfi_cortex1 = dfi_cortex1[['areas','h_iter']]
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='source')
    dfV_CC=dfV_CC.rename(columns={"h_iter": "hs"})
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_CC=dfV_CC.rename(columns={"h_iter": "ht"})
    dfV_CC = dfV_CC.dropna() 
    hg_CC_1 = dfV_CC.ffb*(dfV_CC.ht- dfV_CC.hs)
    
    dfi_thalamus1=dfi_TCCT[(dfi_TCCT.CortexThalamus == 'T')]
    dfi_thalamus1 = dfi_thalamus1[['areas','h_iter']]
    dfV_TCCT = dfV_TCCT.join(dfi_thalamus1.set_index('areas'), on ='source')
    dfV_TCCT=dfV_TCCT.rename(columns={"h_iter": "hs"})
    dfV_TCCT = dfV_TCCT.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_TCCT=dfV_TCCT.rename(columns={"h_iter": "ht"})
    dfV_TCCT = dfV_TCCT.dropna() 
    hg_TCCT_1 = dfV_TCCT.ffb*(dfV_TCCT.ht- dfV_TCCT.hs)
    
    hg_cortex_TCCT_iter = np.mean(hg_CC_1)
    hg_TCCT_iter = np.mean(hg_CC_1.append(hg_TCCT_1))
    
    dfV_CC = dfVC[['source','target','ffb']]
    dfV_TCCT = dfVT[["source","target","ffb"]]
    
    dfi_cortex1 = dfi_TCCT[(dfi_TCCT.CortexThalamus == 'C')]
    dfi_cortex1 = dfi_cortex1[['areas','h0']]
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='source')
    dfV_CC=dfV_CC.rename(columns={"h0": "hs"})
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_CC=dfV_CC.rename(columns={"h0": "ht"})
    dfV_CC = dfV_CC.dropna() 
    hg_CC_1 = dfV_CC.ffb*(dfV_CC.ht- dfV_CC.hs)
    
    dfi_thalamus1=dfi_TCCT[(dfi_TCCT.CortexThalamus == 'T')]
    dfi_thalamus1 = dfi_thalamus1[['areas','h0']]
    dfV_TCCT = dfV_TCCT.join(dfi_thalamus1.set_index('areas'), on ='source')
    dfV_TCCT=dfV_TCCT.rename(columns={"h0": "hs"})
    dfV_TCCT = dfV_TCCT.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_TCCT=dfV_TCCT.rename(columns={"h0": "ht"})
    dfV_TCCT = dfV_TCCT.dropna() 
    hg_TCCT_1 = dfV_TCCT.ffb*(dfV_TCCT.ht- dfV_TCCT.hs)
    
    hg_cortex_TCCT_init = np.mean(hg_CC_1)
    hg_TCCT_init  = np.mean(hg_CC_1.append(hg_TCCT_1))
    
  
    hr_ct_shuffled[i_shuffle] = hg_TCCT_init
    hr_iter_c_shuffled[i_shuffle] = hg_cortex_TCCT_iter
    hr_iter_ct_shuffled[i_shuffle] = hg_TCCT_iter


pd.DataFrame(hr_ct_shuffled).to_excel(output_dir+'shuffled_hg_TCCT_init_cortexthalamus.xls')
pd.DataFrame(hr_iter_c_shuffled).to_excel(output_dir+'shuffled_hg_TCCT_iter_cortex.xls')
pd.DataFrame(hr_iter_ct_shuffled).to_excel(output_dir+'shuffled_hg_TCCT_iter_cortexthalamus.xls')


# In[]: Plot global hierarchy scores of 100 shuffled data with the global hierarchy score of the original data

"""Global hierarchy scores of the original cortico-thalamic + thalamo-cortical connectivity"""

df_hg_TC = pd.read_excel(input_dir+'ghs_TCCT.xls')

hg_all_init = df_hg_TC["hg_TCCT_init"][0]
hg_cortex_iter = df_hg_TC["hg_cortex_TCCT_iter"][0]
hg_all_iter = df_hg_TC["hg_TCCT_iter"][0] 


### No Cre-conf
#hg_all_init = 0.140955843202832
#hg_cortex_iter = 0.10400187059214
#hg_all_iter = 0.12977039470228

#### Cre Conf
#hg_all_init = 0.103134133157259
#hg_cortex_iter = 0.0991790409540279
#hg_all_iter = 0.127804097236645

### WT
#hg_all_init = 0.224612114930841
#hg_cortex_iter = 0.116412327238314
#hg_all_iter = 0.190188780227


hm1 = (hg_all_init-np.mean(hr_ct_shuffled))/np.std(hr_ct_shuffled)              # Z-score for thalamus+cortex before iteration
hm2 = (hg_cortex_iter-np.mean(hr_iter_c_shuffled))/np.std(hr_iter_c_shuffled)   # Z-score for cortex after iteration
hm3 = (hg_all_iter-np.mean(hr_iter_ct_shuffled))/np.std(hr_iter_ct_shuffled)    # Z-score for thalamus+cortex after iteration

""" Figure showing global hierarchy scores of shuffled data & original CT+TC data before iteration, for cortex+thalamus """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_ct_shuffled, bins=bins, label='confidence adjusted')
ax.axvline(x=hg_all_init,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm1))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffledgh_TCCT_all_init.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original CT+TC data after iteration, for cortex only """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_iter_c_shuffled, bins=bins, label='confidence adjusted')
ax.axvline(x=hg_cortex_iter,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm2))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffledgh_TCCT_cortex_iter.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original CT+TC data after iteration, for cortex+thalamus """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_iter_ct_shuffled, bins=bins, label='confidence adjusted')
ax.axvline(x=hg_all_iter,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm3))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffledgh_TCCT_all_iter.pdf", bbox_inches='tight')



