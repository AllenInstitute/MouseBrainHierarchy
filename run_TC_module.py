from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_TC import fit
from IterativeMethod import iterativeTC

"""
Created on Mon Jun 24 16:04:07 2019

@author: Hannah Choi
"""
"""
This code finds intra- or inter-module hierarchies of cortical network 
by maximizing the global hierarchy score followed by interations, 
based on the cortico-cortical (CC) & thalamo-cortical (TC) connectivity data.
"""

# In[]: Set input and output directories 

CreConf = 1                    # 1 if using CC hierarhcy with Cre-confidence; 0 if not

input_dir = r'./Input/'        # Directory with the file "CC_TC_CT_clusters.xlsx" & "clustermapping.xlsx"
output_dir = r'./Output/module/'      # Directory to save the ouputs from the experimental data


# In[]: Define the module. 

module = 'inter_predefined' #'VisualMedial' #'inter_predefined'

# In the paper, we used the following: 'VisualMedial' & 'inter_predefined'
# Possible modules: 'VisualMedial', 'Visual', 'Medial', 'Auditory', 'Somatomotor', 'PFC', 'Lateral', 'inter_predefined', 'inter'

# In[]: 9 clusters of thalamo-courtical source-line-target pairs

xls=pd.ExcelFile(input_dir+"CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'CC+TC clusters with MD avged')
df=df[(df.hemi == "ipsi")&(df.target != "VISC")&(df.source != "SSp-un")&(df.target != "SSp-un")&(df.source != "VISC")]

df = df[(df["Target Major Division"] == "isocortex")&(df["Source Major Division"] == "thalamus")] # T-C connections

list_module = df["Cortical Target Module"].unique()
clu_ffb=pd.read_excel(input_dir+"clustermapping.xlsx")
clu_ffb = clu_ffb.rename(columns={'TC': 'ffb'})

# In[]: For intra medial, lump togeter Visual and Medial & select only the areas within the chosen module

if (module != 'inter') and (module != 'inter_predefined'):   

    if module == 'VisualMedial':
        df.loc[((df["Cortical Target Module"] == "Visual")|(df["Cortical Target Module"] == "Medial")),"Cortical Target Module"]='VisualMedial'
        df = df[(df["Cortical Target Module"] == module)]
        df = df[(df["source"]!='RSPagl')&(df["target"]!='RSPagl')&(df["source"]!='RSPd')&(df["target"]!='RSPd')
        &(df["source"]!='RSPv')&(df["target"]!='RSPv')]
    else:
        df = df[(df["Cortical Target Module"] == module)]

# In[]: If inter-module, change all the target & source area names to the module name

if (module == 'inter') or (module == 'inter_predefined'): 
    for i_module in range(0,len(list_module)):
        df.loc[df["Cortical Target Module"] == list_module[i_module],'target'] = list_module[i_module]
        df.loc[df["Cortical Source Module"] == list_module[i_module],'source'] = list_module[i_module] 


 
# In[]: Trim the dataframe 

dfV1 = df[['source','target','creline','Cluster ID','Cortical Source Module', 'Cortical Target Module']]
dfV1=dfV1.rename(columns={"Cluster ID": "clu"})
dfV1 = dfV1.reset_index(drop=True)
dfV1['clu'] = dfV1['clu'].apply(np.int64)
    
source_areas = dfV1["source"].unique()
target_areas1 = dfV1["target"].unique()   
print target_areas1[~np.in1d(target_areas1, source_areas)]

dfVT = dfV1 
source_areas = dfVT["source"].unique()
target_areas = dfVT["target"].unique()

num_clu = len(dfVT["clu"].unique())

dfVT_concise = dfVT[["source","target","creline","clu","Cortical Source Module", "Cortical Target Module"]].copy()

# In[ ]: If inter-module, we may want to find the mapping rule, that is not pre-defined.

if module == 'inter':
    print dfVT.shape
    logging.debug("performing initial hierarchy assignment")
    hierarchy_vals = fit(dfVT)#, parallel=True, n_procs=-1)

    jmax_raw, jmax = np.argmax(hierarchy_vals, axis=0)
    jmax_raw_val = hierarchy_vals[jmax_raw][0]
    jmax_val = hierarchy_vals[jmax][1]
    logging.debug("RESULTS")
    n = len(dfVT.clu.unique())
    logging.debug("(jmax_raw, val) = ({:0{n}b}, {:.3f})".format(jmax_raw, jmax_raw_val, n=n))
    logging.debug("(jmax,     val) = ({:0{n}b}, {:.3f})".format(jmax, jmax_val, n=n))
    
    results = dict(jmax=bin(2**n+jmax),
                   jmax_val=jmax_val,
                   jmax_raw=bin(2**n+jmax_raw),
                   jmax_raw_val=jmax_raw_val)
    hrc_original = jmax_val
    hr_original = jmax_raw_val


# In[ ]: Define functions needed.

c0=2**num_clu;
def ffb_c (cls):
    """Direction of each cluster with confidence"""
    if module == 'inter':
        b=(bin(c0+jmax)[-(cls)])
        return -(2*int(b)-1) # or -(2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb)   

def ffb_nc (cls):
    """Direction of each cluster without confidence"""
    if module == 'inter':
        b=(bin(c0+jmax_raw)[-(cls)])
        return -(2*int(b)-1)  # or (2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb)

def confidence(df):
     """Returns multiplier which biases towards roughly equal # of FF and FB connections"""    
     count_ff = len(df[df.ffb_c==1])
     count_fb = len(df[df.ffb_c==-1])
     confnew = min(count_ff, count_fb)/(count_ff+count_fb)
     return confnew


def hrf (area):
    '''Hierarchy score of each area without confidence'''
    return -np.mean(dfVT[dfVT.source == area].ffb_nc)

def hrcf (area):
    '''Hierarchy score of each area with confidence'''
    return -np.mean(dfVT[dfVT.source == area].ffb_c)

# In[ ]: Produce expanded data frame with  FF/FB, hierarchy values as source & target for each pair of TC connections
'''ATTENTION! Use the confidence-weighted (biased) mapping of the TC clusters ("ffb_c" & "hrc_s");
Otherwise, the algorithm places all thalamic regions below cortical regions'''

dfVT["ffb_c"]=dfVT["clu"].apply(ffb_c)
dfVT["ffb_nc"]=dfVT["clu"].apply(ffb_nc)
dfVT["hrc_s"]=dfVT["source"].apply(hrcf)
dfVT["hr_s"]=dfVT["source"].apply(hrf)

conf = confidence(dfVT)

dfVT.to_excel(output_dir+'inputexpanded_TC9_'+module+'.xls')

# In[ ]: Find hierarchy scores of thalamic areas within a module or of modules

areas = source_areas
n_areas=len(areas) # 21 thalamic regions
hr=range(0,n_areas)
hrc=range(0,n_areas)

for i in range(0,n_areas):
    hr[i]=-np.mean(dfVT[dfVT.source == areas[i]].ffb_nc) 
    hrc[i]=conf*(-np.mean(dfVT[dfVT.source == areas[i]].ffb_c))

data=[areas,hrc,hr]
data=np.transpose(data)
columns = ['areas','hrc','hr']
dfiT = pd.DataFrame(data,columns=columns)
dfiT.head()

#dfiT.to_excel(output_dir+'initialhierarchy_TC9_'+module+'.xls')

# In[ ]: Iterate thalamic + cortical hierarhcy scores

n_iter = 20

if CreConf == 1:
    dfiC = pd.read_excel(output_dir+"CC_conf_iter_"+module+".xls")
elif CreConf == 0:
    dfiC = pd.read_excel(output_dir+"CC_noconf_iter_"+module+".xls")

dfiC['h'] = dfiC[n_iter]
dfVC = pd.read_excel(output_dir+"inputexpanded_CC9_"+module+".xls")
dfVC = dfVC[["source","target","ffb_nc","ffb_c"]]

if CreConf == 1:
    dfVC["ffb"] = dfVC["ffb_c"]
elif CreConf == 0:
    dfVC["ffb"] = dfVC["ffb_nc"]

dfiT['h'] = dfiT['hrc']
dfVT["ffb"] = dfVT["ffb_c"]


dfiT = dfiT[["areas","h"]]
dfiC = dfiC[["areas","h"]]
dfVT = dfVT[["source","target","ffb"]]
dfVC = dfVC[["source","target","ffb"]]

hr_iter = iterativeTC(dfiC, dfVC, dfiT, dfVT, n_iter)


iteration=np.arange(0,n_iter+1,1)
n_area=np.shape(hr_iter)[0]
allareas = hr_iter["areas"].unique()

##########################################################
""" Figure of hierarchy score iterations """
fig,ax=plt.subplots()
for i in range(0,n_area):
    y=np.squeeze(np.asarray(hr_iter[hr_iter.areas==allareas[i]].ix[:,1::1]))
    ax.plot(iteration,y)
ax.set_xlim([0, n_iter])
ax.set_xticks(np.arange(0, n_iter, step=5))
ax.set_xlabel('iter')
ax.set_ylabel('hierarchy value')
ax.set_title('confidence adjusted')
plt.show()
#fig.savefig(output_dir+"TC_CCnoconf_iter_"+module+".pdf", bbox_inches='tight')
##########################################################

##########################################################
""" Figure showing correlation between hierarchy scores before & after iterations"""
hr_final = hr_iter[:][n_iter]
hr_initial = hr_iter[:][0]   
f = plt.figure()
plt.plot(hr_initial,hr_final,'ro')
plt.xlabel('initial hierarchy (conf)')
plt.ylabel('final hierarchy (conf)')
plt.title('r='+str(np.corrcoef(hr_initial, hr_final)[0, 1]))
plt.show()
#f.savefig(output_dir+"TC_init_vs_fin_CCnoconf_"+module+".pdf", bbox_inches='tight')
##########################################################


##########################################################
'''Save hierarchy scores before and after iteration'''

for i_area in range(0,n_area):
    if hr_iter['areas'][i_area] in list(dfiC['areas']):
        hr_iter.loc[i_area,'CortexThalamus'] = 'C'
    else:
        hr_iter.loc[i_area,'CortexThalamus'] = 'T'

hr_iter = hr_iter[['areas','CortexThalamus', 0,n_iter] ]  
hr_iter_save = hr_iter #hr_iter[(hr_iter.CortexThalamus=='C')]

if CreConf == 1:
    hr_iter_save.to_excel(output_dir+'TC_CCconf_iter_'+module+'.xls')
elif CreConf == 0: 
    hr_iter_save.to_excel(output_dir+'TC_CCnoconf_iter_'+module+'.xls')

##########################################################

 # In[]: Print out global hierarchy scores of CC + TC connectivity data before and after iteration
 
dfi_TC = hr_iter[["CortexThalamus","areas",0,n_iter]]
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

hg_thalmus_TC_init = np.mean(hg_TC_1)
hg_cortex_TC_init = np.mean(hg_CC_1)
hg_TC_init = np.mean(hg_CC_1.append(hg_TC_1))


print('hg of TC thalamus='+str(hg_thalmus_TC_init))
print('hg of CC+TC before iterate cortex & thalamus='+str(hg_TC_init)) 
print('hg of CC+TC iterate cortex='+str(hg_cortex_TC_iter))
print('hg of CC+TC iterate cortex & thalamus='+str(hg_TC_iter))


'''Save global hierarchy scores'''
newDF= pd.DataFrame([])
newDF=newDF.append(pd.DataFrame({'hg_thalmus_TC_init':hg_thalmus_TC_init, 'hg_TC_init':hg_TC_init,
                                 'hg_cortex_TC_iter':hg_cortex_TC_iter, 'hg_TC_iter':hg_TC_iter},index=[0]))
newDF.to_excel(output_dir+'ghs_TC_'+module+'.xls')
