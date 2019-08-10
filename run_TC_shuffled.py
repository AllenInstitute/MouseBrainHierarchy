from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_TC_shuffled import fit
from IterativeMethod import iterativeTC

"""
Created on Mon Oct 29 13:48:02 2018

@author: Hannah Choi
"""
"""
This code genereates shuffled versions of the TC connectivity data and 
computes the global hierarchy scores of the shuffled data.
"""

# In[]: Set input and output directories 

input_dir = r'./Input/'        # Directory with the file "CC_TC_CT_clusters.xlsx"
input_dir2 = r'./Output/'               # Directory with the file "ghc_TC.xls"
output_dir = r'./Output/shuffled/'      # Directory to save the ouputs from the shuffled experimental data

''' ATTENTION! Change the "df_cortex" accordingly in func_unsupervised_TC as well! '''
CreConf = 1                    # 1 if using CC hierarhcy with Cre-confidence; 0 if not

# In[]: Read in the excel file with source-target-creline pairs and their cluster numbers. Construct a dataframe using only the thalamo-cortical connections. 

xls=pd.ExcelFile(input_dir+"CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'CC+TC clusters with MD avged')
df=df[(df.hemi == "ipsi")&(df.target != "VISC")&(df.source != "SSp-un")&(df.target != "SSp-un")&(df.source != "VISC")]
df = df[(df["Target Major Division"] == "isocortex")&(df["Source Major Division"] == "thalamus")] # T-C connections

dfV1 = df[['source','target','creline','Cluster ID']]
dfV1=dfV1.rename(columns={"Cluster ID": "clu"})
dfV1 = dfV1.reset_index(drop=True)
dfV1['clu'] = dfV1['clu'].apply(np.int64)
    
source_areas = dfV1["source"].unique()
target_areas1 = dfV1["target"].unique()   
print target_areas1[~np.in1d(target_areas1, source_areas)]

dfVT0 = dfV1 
source_areas = dfVT0["source"].unique()
target_areas = dfVT0["target"].unique()

num_clu = len(dfVT0["clu"].unique())

dfVT0 = dfVT0[["source","target","creline","clu"]]#.copy()



# In[ ]: Find global hierarchy scores of shuffled TC connectivity data

n_iter = 10
by_creline = 0
line_list = dfVT0["creline"].unique()
n_shuffle = 100

hr_ct_shuffled = np.zeros(n_shuffle)
hr_iter_c_shuffled = np.zeros(n_shuffle)
hr_iter_ct_shuffled = np.zeros(n_shuffle)

conf_shuffled = np.zeros(n_shuffle)

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
        
    elif by_creline == 1:
        for i in range(0,len(line_list)):
            source_list = dfVT_shuffled[(dfVT_shuffled.creline == str(line_list[i]))].source
            target_list = dfVT_shuffled[(dfVT_shuffled.creline == str(line_list[i]))].target
            source_shuffled = source_list.sample(frac=1).reset_index(drop=True)
            target_shuffled = target_list.sample(frac=1).reset_index(drop=True)
            source_shuffled.index = source_list.index
            target_shuffled.index = target_list.index        
            dfVT_shuffled["source"][source_list.index]=np.array(source_shuffled)
            dfVT_shuffled["target"][target_list.index]=np.array(target_shuffled)



    hierarchy_vals = fit(dfVT_shuffled,i_shuffle)

    jmax_raw, jmax = np.argmax(hierarchy_vals, axis=0)
    jmax_raw_val = hierarchy_vals[jmax_raw][0]
    jmax_val = hierarchy_vals[jmax][1]
    logging.debug("RESULTS")
    n = len(dfVT0.clu.unique())
    logging.debug("(jmax_raw, val) = ({:0{n}b}, {:.3f})".format(jmax_raw, jmax_raw_val, n=n))
    logging.debug("(jmax,     val) = ({:0{n}b}, {:.3f})".format(jmax, jmax_val, n=n))
    
    results = dict(jmax=bin(2**n+jmax),
                   jmax_val=jmax_val,
                   jmax_raw=bin(2**n+jmax_raw),
                   jmax_raw_val=jmax_raw_val)

    hrc_original = jmax_val
    hr_original = jmax_raw_val

    ###########################################################################
    """Define functions needed"""
    
    c0=2**num_clu;
    def ffb_c (cls):
        """Direction of each cluster with confidence"""
        b=(bin(c0+jmax)[-(cls)])
        return -(2*int(b)-1)
    
    def ffb_nc (cls):
        """Direction of each cluster without confidence"""
        b=(bin(c0+jmax_raw)[-(cls)])
        return -(2*int(b)-1) #-(2*int(b)-1)
    
    def confidence(df):
         """Returns multiplier which biases towards roughly equal # of FF and FB connections"""   
         count_ff = len(df[df.ffb_c==1])
         count_fb = len(df[df.ffb_c==-1])
         confnew = min(count_ff, count_fb)/(count_ff+count_fb)
         return confnew
        
    def hrf (area):
        '''Hierarchy score of each area without confidence'''
        return -np.mean(dfVT_shuffled[dfVT_shuffled.source == area].ffb_nc)
    
    def hrcf (area):
        '''Hierarchy score of each area with confidence'''
        return -np.mean(dfVT_shuffled[dfVT_shuffled.source == area].ffb_c)
    ###########################################################################
    
    
    ###########################################################################
    """Produce an expanded data frame with  FF/FB, hierarchy values as source & target 
    for each pair of TC connections"""
    
    dfVT_shuffled["ffb_c"]=dfVT_shuffled["clu"].apply(ffb_c)
    dfVT_shuffled["ffb_nc"]=dfVT_shuffled["clu"].apply(ffb_nc)
    dfVT_shuffled["hrc_s"]=dfVT_shuffled["source"].apply(hrcf)
    dfVT_shuffled["hr_s"]=dfVT_shuffled["source"].apply(hrf)
    
    conf = confidence(dfVT_shuffled)
    
    ###########################################################################

    ###########################################################################
    '''Finding initial hierarchy score of each thalamic area (21)'''    
    areas = source_areas
    n_areas=len(areas) # 21 thalamic regions
    hr=range(0,n_areas)
    hrc=range(0,n_areas)
    
    for i in range(0,n_areas):
        hr[i]=-np.mean(dfVT_shuffled[dfVT_shuffled.source == areas[i]].ffb_nc) 
        hrc[i]=conf*(-np.mean(dfVT_shuffled[dfVT_shuffled.source == areas[i]].ffb_c))
    
    data=[areas,hrc]
    data=np.transpose(data)
    columns = ['areas','hrc']
    dfiT = pd.DataFrame(data,columns=columns)  
    
    ###########################################################################
    
    ###########################################################################
    '''Iterate thalamic + cortical hierarhcy scores'''
    if CreConf == 1:
        dfiC = pd.read_excel(output_dir+'CCshuffled_conf_iter'+str(i_shuffle)+'.xls')  
    elif CreConf == 0:
        dfiC = pd.read_excel(output_dir+'CCshuffled_noconf_iter'+str(i_shuffle)+'.xls') 

    ''' Note that n_iter=10 in run_CC_shuffled.py'''
    dfiC['h'] = dfiC[10]  
    dfVC = pd.read_excel(output_dir+'inputexpanded_CC_shuffled'+str(i_shuffle)+'.xls')
    dfVC = dfVC[["source","target","ffb_nc","ffb_c"]]
    if CreConf == 0:
        dfVC["ffb"] = dfVC["ffb_nc"]
    elif CreConf == 1:
        dfVC["ffb"] = dfVC["ffb_c"]
        
    dfiT['h'] = dfiT['hrc']
    dfVT0["ffb"] = dfVT0["ffb_c"]
    
    
    dfiT = dfiT[["areas","h"]]
    dfiC = dfiC[["areas","h"]]
    dfVT = dfVT0[["source","target","ffb"]]
    dfVC = dfVC[["source","target","ffb"]]
    
    hr_iter = iterativeTC(dfiC, dfVC, dfiT, dfVT, n_iter)
    
    
    iteration=np.arange(0,n_iter+1,1)
    n_area=np.shape(hr_iter)[0]
    allareas = hr_iter["areas"].unique()
    
       
    for i_area in range(0,n_area):
        if hr_iter['areas'][i_area] in list(dfiC['areas']):
            hr_iter.loc[i_area,'CortexThalamus'] = 'C'
        else:
            hr_iter.loc[i_area,'CortexThalamus'] = 'T'
    
    hr_iter = hr_iter[['areas','CortexThalamus', 0, n_iter] ]     
    
    #    if CreConf == 1:
    #        hr_iter.to_excel(output_dir+'TC_CCconf_iter'+str(i_shuffle)+'.xls')
    #    elif CreConf == 0:
    #        hr_iter.to_excel(output_dir+'TC_CCnoconf_iter'+str(i_shuffle)+'.xls')
    ###########################################################################


    ###########################################################################
    '''global hierarchy score of the shuffled data before & after iteration'''
     
    dfi_TC = hr_iter[["CortexThalamus","areas",0, n_iter]]
    dfV_TC = dfVT[["source","target","ffb"]]
    dfV_CC = dfVC[['source','target','ffb']]
    dfi_TC = dfi_TC.rename(columns={0: "h0", n_iter:"h_iter"})
    
    dfi_cortex1 = dfi_TC[(dfi_TC.CortexThalamus == 'C')]
    dfi_cortex1 = dfi_cortex1[['areas','h_iter']]
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='source')
    dfV_CC=dfV_CC.rename(columns={"h_iter": "hs"})
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_CC=dfV_CC.rename(columns={"h_iter": "ht"})
    dfV_CC = dfV_CC.dropna() 
    hg_CC_1 = dfV_CC.ffb*(dfV_CC.ht- dfV_CC.hs)
    
    dfi_thalamus1=dfi_TC[(dfi_TC.CortexThalamus == 'T')]
    dfi_thalamus1 = dfi_thalamus1[['areas','h_iter']]
    dfV_TC = dfV_TC.join(dfi_thalamus1.set_index('areas'), on ='source')
    dfV_TC=dfV_TC.rename(columns={"h_iter": "hs"})
    dfV_TC = dfV_TC.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_TC=dfV_TC.rename(columns={"h_iter": "ht"})
    dfV_TC = dfV_TC.dropna() 
    hg_TC_1 = dfV_TC.ffb*(dfV_TC.ht- dfV_TC.hs)
    
    hg_cortex_TC_iter = np.mean(hg_CC_1)
    hg_TC_iter = np.mean(hg_CC_1.append(hg_TC_1))
    
    
    dfV_TC = dfVT[["source","target","ffb"]]
    dfV_CC = dfVC[['source','target','ffb']]
    
    dfi_cortex1 = dfi_TC[(dfi_TC.CortexThalamus == 'C')]
    dfi_cortex1 = dfi_cortex1[['areas','h0']]
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='source')
    dfV_CC=dfV_CC.rename(columns={"h0": "hs"})
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_CC=dfV_CC.rename(columns={"h0": "ht"})
    dfV_CC = dfV_CC.dropna() 
    hg_CC_1 = dfV_CC.ffb*(dfV_CC.ht- dfV_CC.hs)
    
    dfi_thalamus1=dfi_TC[(dfi_TC.CortexThalamus == 'T')]
    dfi_thalamus1 = dfi_thalamus1[['areas','h0']]
    dfV_TC = dfV_TC.join(dfi_thalamus1.set_index('areas'), on ='source')
    dfV_TC=dfV_TC.rename(columns={"h0": "hs"})
    dfV_TC = dfV_TC.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_TC=dfV_TC.rename(columns={"h0": "ht"})
    dfV_TC = dfV_TC.dropna() 
    hg_TC_1 = dfV_TC.ffb*(dfV_TC.ht- dfV_TC.hs)
    
    hg_thalamus_TC_init = np.mean(hg_CC_1)
    hg_cortex_TC_init = np.mean(hg_CC_1)
    hg_TC_init = np.mean(hg_CC_1.append(hg_TC_1))
    

    
    hr_ct_shuffled[i_shuffle] = hg_TC_init
    hr_iter_c_shuffled[i_shuffle] = hg_cortex_TC_iter
    hr_iter_ct_shuffled[i_shuffle] = hg_TC_iter



pd.DataFrame(hr_ct_shuffled).to_excel(output_dir+'TC_hg_shuffled_all_init.xls')
pd.DataFrame(hr_iter_c_shuffled).to_excel(output_dir+'TC_hg_shuffled_cortex_iter.xls')
pd.DataFrame(hr_iter_ct_shuffled).to_excel(output_dir+'TC_hg_shuffled_all_iter.xls')

# In[]: Plot global hierarchy scores of 100 shuffled data with the global hierarchy score of the original data

"""Global hierarchy scores of the original thalamo-cortical connectivity"""

df_hg_TC = pd.read_excel(input_dir2+'ghs_TC.xls')

hg_all_init = df_hg_TC["hg_TC_init"][0]
hg_cortex_iter = df_hg_TC["hg_cortex_TC_iter"][0]
hg_all_iter = df_hg_TC["hg_TC_iter"][0] 

### No conf
#hg_all_init = 0.0905605156930546
#hg_cortex_iter = 0.102959136598683
#hg_all_iter = 0.119694181094546

### conf
#hg_all_init = 0.0905605156930546
#hg_cortex_iter = 0.102959136598683
#hg_all_iter = 0.119694181094546

#### WT
#hg_all_init = 0.186626390380453
#hg_cortex_iter = 0.147195762950937
#hg_all_iter = 0.171537027881485

hm1 = (hg_all_init-np.mean(hr_ct_shuffled))/np.std(hr_ct_shuffled)              # Z-score for thalamus+cortex before iteration
hm2 = (hg_cortex_iter-np.mean(hr_iter_c_shuffled))/np.std(hr_iter_c_shuffled)   # Z-score for cortex after iteration
hm3 = (hg_all_iter-np.mean(hr_iter_ct_shuffled))/np.std(hr_iter_ct_shuffled)    # Z-score for thalamus+cortex after iteration


""" Figure showing global hierarchy scores of shuffled data & original TC data before iteration, for cortex+thalamus """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_ct_shuffled, bins=bins, label='confidence adjusted')
ax.axvline(x=hg_all_init,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm1))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffledgh_TC_all_init.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original TC data after iteration, for cortex only """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_iter_c_shuffled, bins=bins, label='confidence adjusted')
ax.axvline(x=hg_cortex_iter,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm2))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffledgh_TC_cortex_iter.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original TC data after iteration, for cortex+thalamus """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_iter_ct_shuffled, bins=bins, label='confidence adjusted')
ax.axvline(x=hg_all_iter,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm3))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffledgh_TC_all_iter.pdf", bbox_inches='tight')

