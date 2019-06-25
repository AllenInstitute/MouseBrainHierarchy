from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_TC import fit
from IterativeMethod import iterativeTC

"""
Created on Mon Oct 29 13:48:02 2018

@author: Hannah Choi
"""
"""
This code finds hierarchy scores of cortical & thalamic areas by maximizing the global hierarchy score 
followed by interations, based on the cortico-cortical (CC) AND the thalamo-cortical (TC) connectivity data.
"""

# In[]: Set input and output directories 

input_dir = r'./Input/'     # Directory with the file "CC_TC_CT_clusters.xlsx" 
output_dir = r'./Output/'   # Directory to save the ouputs from the experimental data

''' ATTENTION! Change the "df_cortex" accordingly in func_unsupervised_TC as well! '''
CreConf = 1                 # 1 if using CC hierarhcy with Cre-confidence; 0 if not


# In[]: Read the excel file with source-target-creline pairs and their cluster numbers. Use only the thalamo-cortical connections. 

xls=pd.ExcelFile(input_dir+"CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'CC+TC clusters with MD avged')

df=df[(df.hemi == "ipsi")&(df.target != "VISC")&(df.source != "SSp-un")&(df.target != "SSp-un")&(df.source != "VISC")]
df = df[(df["Target Major Division"] == "isocortex")&(df["Source Major Division"] == "thalamus")] # Consider T-C connections only

dfV1 = df[['source','target','creline','Cluster ID']]
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


# In[]: Map 9 clusters of thalamo-cortical connections to FF/FB directions, by maximizing the global hierarchy score

hierarchy_vals = fit(dfVT)

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


# In[ ]: Define functions needed

c0=2**num_clu

def ffb_c (cls):
    """Direction of each cluster with confidence"""
    b=(bin(c0+jmax)[-(cls)])
    return -(2*int(b)-1)

def ffb_nc (cls):
    """Direction of each cluster without confidence"""
    b=(bin(c0+jmax_raw)[-(cls)])
    return -(2*int(b)-1) 

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

dfVT.to_excel(output_dir+'inputexpanded_TC9.xls')

# In[ ]: Find hierarchy score for each of thalamic areas (21)

areas = source_areas
n_areas=len(areas) 
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

#dfiT.to_excel(output_dir+'initialhierarchy_TC9.xls')

# In[ ]: Iterate thalamic + cortical hierarchy scores

n_iter = 20

##########################################################
'''Load results from CC hierarchy'''

if CreConf == 0:
    dfiC = pd.read_excel(output_dir+"CC_noconf_iter.xls")
elif CreConf == 1:
    dfiC = pd.read_excel(output_dir+"CC_conf_iter.xls")
    
dfiC['h'] = dfiC[n_iter]
dfVC = pd.read_excel(output_dir+"inputexpanded_CC9.xls")
dfVC = dfVC[["source","target","ffb_nc","ffb_c"]]

if CreConf == 0:
    dfVC["ffb"] = dfVC["ffb_nc"]
elif CreConf == 1:
    dfVC["ffb"] = dfVC["ffb_c"]
##########################################################


##########################################################
'''Use initial hierarchy scores and cluster directions (FF/FB) 
for TC connections with the "equal #s of FF/FB"-confidence weight'''  

dfiT['h'] = dfiT['hrc']
dfVT["ffb"] = dfVT["ffb_c"]
##########################################################

dfiT = dfiT[["areas","h"]]
dfiC = dfiC[["areas","h"]]
dfVT = dfVT[["source","target","ffb"]]
dfVC = dfVC[["source","target","ffb"]]


##########################################################
'''Iterations ''' 

hr_iter = iterativeTC(dfiC, dfVC, dfiT, dfVT, n_iter)
##########################################################

iteration=np.arange(0,n_iter+1,1)
n_area=np.shape(hr_iter)[0]
allareas = hr_iter["areas"].unique()


##########################################################
'''Figure of hierarchy score iterations'''

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
#if CreConf == 0:
#    fig.savefig(output_dir+"TC_CCnoconf_iter.pdf", bbox_inches='tight')
#elif CreConf == 1:
#    fig.savefig(output_dir+"TC_CCconf_iter.pdf", bbox_inches='tight')
##########################################################


##########################################################
'''Figure showing correlation between hierarchy scores before & after iterations'''

hr_final = hr_iter[:][n_iter]
hr_initial = hr_iter[:][0]   
f = plt.figure()
plt.plot(hr_initial,hr_final,'ro')
plt.xlabel('initial hierarchy (conf)')
plt.ylabel('final hierarchy (conf)')
plt.title('r='+str(np.corrcoef(hr_initial, hr_final)[0, 1]))
plt.show()
#if CreConf == 0:
#    f.savefig(output_dir+"TC_init_vs_fin_CCnoconf.pdf", bbox_inches='tight')
#elif CreConf == 1:
#    f.savefig(output_dir+"TC_init_vs_fin_CCconf.pdf", bbox_inches='tight')
##########################################################

##########################################################
'''Save hierarchy scores before and after iteration'''

for i_area in range(0,n_area):
    if hr_iter['areas'][i_area] in list(dfiC['areas']):
        hr_iter.loc[i_area,'CortexThalamus'] = 'C'
    else:
        hr_iter.loc[i_area,'CortexThalamus'] = 'T'

hr_iter = hr_iter[['areas','CortexThalamus', 0,n_iter] ]  

if CreConf == 0:
    hr_iter.to_excel(output_dir+'TC_CCnoconf_iter.xls')
elif CreConf == 1:
    hr_iter.to_excel(output_dir+'TC_CCconf_iter.xls')
##########################################################

 # In[]: Print out global hierarchy scores for the TC connectivity data before and after iteration
 
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
newDF.to_excel(output_dir+'ghs_TC.xls')
