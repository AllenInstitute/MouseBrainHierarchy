from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_CC import fit
from IterativeMethod import iterativeCC

"""
Created on Sun Mar 10 17:16:12 2019

@author: Hannah Choi
"""
"""
This code genereates shuffled versions of the CC connectitivity data and 
computes the global hierarchy scores of the shuffled data.
"""


# In[]: Set input and output directories 

input_dir = r'./Input/'                     # Directory with the file "CC_TC_CT_clusters.xlsx"
input_dir2 = r'./Output/'                   # Directory with the file "ghc_CC.xls"
output_dir = r'./Output/shuffled/'          # Directory to save the ouputs from the shuffled experimental data

# In[]: Read in the excel file with source-target-creline pairs and their cluster numbers. Construct a dataframe using only the cortico-cortical connections. 

xls=pd.ExcelFile(input_dir+"CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'CC+TC clusters with MD avged')

df=df[(df.hemi == "ipsi")&(df.creline != "C57BL/6J / Emx1")&(df.target != "VISC")&(df.source != "SSp-un")&(df.target != "SSp-un")&(df.source != "VISC")] # Consider ipsilateral, Cre connectivity data only

df = df[(df["Target Major Division"] == "isocortex")&(df["Source Major Division"] == "isocortex")]  # Consider C-C connections only

dfV1 = df[['source','target','creline','Cluster ID']]
dfV1=dfV1.rename(columns={"Cluster ID": "clu"})

dfV1 = dfV1.reset_index(drop=True)
dfV1['clu'] = dfV1['clu'].apply(np.int64)

source_areas = dfV1["source"].unique()
target_areas1 = dfV1["target"].unique()   
print target_areas1[~np.in1d(target_areas1, source_areas)]
 
dfV = dfV1 
source_areas = dfV["source"].unique()
target_areas = dfV["target"].unique()

areas = source_areas
n_areas=len(areas)

num_clu = len(dfV["clu"].unique())

dfV_concise = dfV[["source","target","creline","clu"]].copy()


# In[]: Global hierarchy scores of shuffled CC connectivity data

n_iter = 10 # or 5, for fast process
n_shuffle = 100
by_creline = 1   # 1 if shuffle within each Cre-line; 0 if shuffle across Cre-lines 
line_list = dfV_concise["creline"].unique()

hrc_shuffled = np.zeros(n_shuffle)
hr_shuffled = np.zeros(n_shuffle)

hr_init_shuffled = np.zeros(n_shuffle)
hr_iter_shuffled = np.zeros(n_shuffle)
hrc_iter_shuffled = np.zeros(n_shuffle)
hrc_init_shuffled = np.zeros(n_shuffle)

for i_shuffle in range(0,n_shuffle):
    dfV_shuffled = dfV_concise#.copy()
    
    if by_creline == 0:
        source_list= dfV_shuffled.source
        target_list= dfV_shuffled.target
        source_shuffled = source_list.sample(frac=1).reset_index(drop=True)
        target_shuffled = target_list.sample(frac=1).reset_index(drop=True)
        source_shuffled.index = source_list.index
        target_shuffled.index = target_list.index 
        dfV_shuffled["source"][source_list.index]=np.array(source_shuffled)
        dfV_shuffled["target"][target_list.index]=np.array(target_shuffled)
        
    elif by_creline == 1:
        for i in range(0,len(line_list)):
            source_list = dfV_shuffled[(dfV_shuffled.creline == str(line_list[i]))].source
            target_list = dfV_shuffled[(dfV_shuffled.creline == str(line_list[i]))].target
            source_shuffled = source_list.sample(frac=1).reset_index(drop=True)
            target_shuffled = target_list.sample(frac=1).reset_index(drop=True)
            source_shuffled.index = source_list.index
            target_shuffled.index = target_list.index        
            dfV_shuffled["source"][source_list.index]=np.array(source_shuffled)
            dfV_shuffled["target"][target_list.index]=np.array(target_shuffled)
    

    hierarchy_vals = fit(dfV_shuffled)
        
    jmax_raw, jmax = np.argmax(hierarchy_vals, axis=0)
    jmax_raw_val = hierarchy_vals[jmax_raw][0]
    jmax_val = hierarchy_vals[jmax][1]
    logging.debug("RESULTS")
    n = num_clu #len(df.clu.unique())
    logging.debug("(jmax_raw, val) = ({:0{n}b}, {:.3f})".format(jmax_raw, jmax_raw_val, n=n))
    logging.debug("(jmax,     val) = ({:0{n}b}, {:.3f})".format(jmax, jmax_val, n=n))
    
    results = dict(jmax=bin(2**n+jmax),
                   jmax_val=jmax_val,
                   jmax_raw=bin(2**n+jmax_raw),
                   jmax_raw_val=jmax_raw_val)
    
    print('i_shuffle='+str(i_shuffle))
    if jmax_val<0:
        print('jmax='+str(bin(jmax)))
    
    hrc_shuffled[i_shuffle]=jmax_val
    hr_shuffled[i_shuffle]=jmax_raw_val
    
    ###########################################################################    
    """Define functions needed"""
    
    c0=2**num_clu
    
    def ffb_c (cls):
        """Direction of each cluster with Cre-confidence"""
        b=(bin(c0+jmax)[-(cls)])
        return -(2*int(b)-1)    # or (2*int(b)-1)
    
    
    def ffb_nc (cls):
        """Direction of each cluster without Cre-confidence"""
        b=(bin(c0+jmax_raw)[-(cls)])
        return -(2*int(b)-1)
    
    def cre_confidence1(df):
        """Returns confidence of cre lines"""
        func = lambda x: 1 - np.abs(x.mean())
        return df.groupby('creline')['ffb_c'].transform(func)
    
    def hrf (area):
        '''Hierarchy score of each area without Cre-confidence'''
        return ((-np.mean(dfV_shuffled[dfV_shuffled.source == area].ffb_nc)
                 +np.mean(dfV_shuffled[dfV_shuffled.target == area].ffb_nc))/2)
    
    def hrcf (area):
        '''Hierarchy score of each area with Cre-confidence'''
        return ((-np.mean(dfV_shuffled[dfV_shuffled.source == area].ffb_c*dfV_shuffled[dfV_shuffled.source == area].conf)
                 +np.mean(dfV_shuffled[dfV_shuffled.target == area].ffb_c*dfV_shuffled[dfV_shuffled.target == area].conf))/2)
    ###########################################################################
    
    ###########################################################################
    """Produce expanded data frame with  FF/FB, Cre-confidence, hierarchy values as source & target 
    for each pair of CC connections"""
    
    dfV_shuffled.loc[:,"ffb_c"]=dfV_shuffled["clu"].apply(ffb_c)
    dfV_shuffled.loc[:,"ffb_nc"]=dfV_shuffled["clu"].apply(ffb_nc)
    dfV_shuffled.loc[:, "conf"] = cre_confidence1(dfV_shuffled)
    dfV_shuffled.loc[:,"hrc_s"]=dfV_shuffled["source"].apply(hrcf)
    dfV_shuffled.loc[:,"hrc_t"]=dfV_shuffled["target"].apply(hrcf)
    dfV_shuffled.loc[:,"hr_s"]=dfV_shuffled["source"].apply(hrf)
    dfV_shuffled.loc[:,"hr_t"]=dfV_shuffled["target"].apply(hrf)
    
    dfV_shuffled.to_excel(output_dir+'inputexpanded_CC_shuffled'+str(i_shuffle)+'.xls')
    ###########################################################################
    
    hrs=range(0,n_areas)
    hrt=range(0,n_areas)
    hr=range(0,n_areas)
    hrc=range(0,n_areas)
    for i in range(0,n_areas):
        hrs[i]=-np.mean(dfV_shuffled[dfV_shuffled.source == areas[i]].ffb_nc)
        hrt[i]=np.mean(dfV_shuffled[dfV_shuffled.target == areas[i]].ffb_nc)
        hr[i]=(hrs[i]+hrt[i])/2
        hrc[i]=0.5*(-np.mean(dfV_shuffled[dfV_shuffled.source == areas[i]].ffb_c*dfV_shuffled[dfV_shuffled.source == areas[i]].conf)
        +np.mean(dfV_shuffled[dfV_shuffled.target == areas[i]].ffb_c*dfV_shuffled[dfV_shuffled.target == areas[i]].conf))
     
    data=[areas,hrc,hr]
    data=np.transpose(data)
    columns = ['areas','hrc','hr']
    dfi_shuffled = pd.DataFrame(data,columns=columns) 
    
    hr_iter, hrc_iter = iterativeCC(dfi_shuffled,dfV_shuffled,n_iter)
    
   
    hrc_iter = hrc_iter[['areas',0,n_iter]]
    hr_iter = hr_iter[['areas',0,n_iter]]  
    
    ###########################################################################
    '''Save hierarchy scores of cortical areas in the shuffled data'''
    hr_iter.to_excel(output_dir+'CCshuffled_noconf_iter'+str(i_shuffle)+'.xls')
    hrc_iter.to_excel(output_dir+'CCshuffled_conf_iter'+str(i_shuffle)+'.xls')
    ###########################################################################
    
    dfV_temp = dfV_shuffled[['source','target','clu','ffb_c','ffb_nc','conf']]
    dfi_temp = hr_iter[['areas',0,n_iter]]
    dfi_temp = dfi_temp.rename(columns={0: "h0", n_iter:"h_iter"})
    dfi_temp_conf = hrc_iter[['areas',0,n_iter]]
    dfi_temp_conf = dfi_temp_conf.rename(columns={0: "h0", n_iter:"h_iter"})

    dfi_t = dfi_temp[['areas','h0','h_iter']]
    dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
    dfV_temp =dfV_temp.rename(columns={"h_iter": "hs_iter"})
    dfV_temp =dfV_temp.rename(columns={"h0": "hs0"})
    dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
    dfV_temp =dfV_temp.rename(columns={"h_iter": "ht_iter"})
    dfV_temp =dfV_temp.rename(columns={"h0": "ht0"})
    dfV_temp = dfV_temp.dropna()
    
    ###########################################################################
    '''global hierarchy score of the shuffled data before & after iteration, without Cre-confidence'''
#    hr_init_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_nc*np.sign(dfV_temp.ht0 - dfV_temp.hs0))
#    hr_iter_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_nc*np.sign(dfV_temp.ht_iter - dfV_temp.hs_iter))
    hr_init_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_nc*(dfV_temp.ht0 - dfV_temp.hs0))
    hr_iter_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_nc*(dfV_temp.ht_iter - dfV_temp.hs_iter))
    
    dfV_temp = dfV_shuffled[['source','target','clu','ffb_c','ffb_nc','conf']]
    dfi_t = dfi_temp_conf[['areas','h0','h_iter']]
    dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
    dfV_temp =dfV_temp.rename(columns={"h_iter": "hs_iter"})
    dfV_temp =dfV_temp.rename(columns={"h0": "hs0"})
    dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
    dfV_temp =dfV_temp.rename(columns={"h_iter": "ht_iter"})
    dfV_temp =dfV_temp.rename(columns={"h0": "ht0"})
    dfV_temp = dfV_temp.dropna() 
    
     ###########################################################################
    '''global hierarchy score of the shuffled data before & after iteration, with Cre-confidence'''   
#    hrc_init_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_c*np.sign(dfV_temp.ht0 - dfV_temp.hs0))
#    hrc_iter_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_c*np.sign(dfV_temp.ht_iter - dfV_temp.hs_iter))
    hrc_init_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_c*(dfV_temp.ht0 - dfV_temp.hs0))
    hrc_iter_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_c*(dfV_temp.ht_iter - dfV_temp.hs_iter))
    ###########################################################################    
    
    
pd.DataFrame(hrc_init_shuffled).to_excel(output_dir +'CC_hg_init_shuffled.xls')
pd.DataFrame(hrc_iter_shuffled).to_excel(output_dir +'CC_hg_iter_shuffled.xls')

# In[]: Plot global hierarchy scores of 100 shuffled data with the global hierarchy score of the original data

"""Global hierarchy scores of the original cortico-cortical connectivity"""
df_hg_CC = pd.read_excel(input_dir2+'ghc_CC.xls')

hg_CC_conf_init = df_hg_CC["hg_CC_conf_init"][0]
hg_CC_conf_iter = df_hg_CC["hg_CC_conf_iter"][0]

hmc_init = (hg_CC_conf_init-np.mean(hrc_init_shuffled))/np.std(hrc_init_shuffled) # Z-score before iteration
hmc_iter = (hg_CC_conf_iter-np.mean(hrc_iter_shuffled))/np.std(hrc_iter_shuffled) # Z-score after iteration

""" Figure showing global hierarchy scores of shuffled data & original CC data before & after iteration """
fig,ax=plt.subplots()
ax.hist(hrc_init_shuffled, bins=25, label='before iterate')
ax.axvline(x=hg_CC_conf_init,linestyle='--')
ax.hist(hrc_iter_shuffled, bins=25, label='after iterate',color='red')
ax.axvline(x=hg_CC_conf_iter,linestyle='--',color='red')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hmc(init)='+str(hmc_init)+'; hmc(iter)='+str(hmc_iter))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffled_globalhierarchy_conf_CC.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original CC data before iteration """
fig,ax=plt.subplots()
ax.hist(hrc_init_shuffled, bins=10, label='before iterate')
ax.axvline(x=hg_CC_conf_init,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm(init)='+str(hmc_init))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffled_globalhierarchy_CC_init.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original CC data after iteration """
fig,ax=plt.subplots()
ax.hist(hrc_iter_shuffled, bins=25, label='after iterate',color='red')
ax.axvline(x=hg_CC_conf_iter,linestyle='--',color='red')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm(iter)='+str(hmc_iter))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffled_globalhierarchy_CC_iter.pdf", bbox_inches='tight')
