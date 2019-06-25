from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IterativeMethod import iterativeTCCT

"""
Created on Mon Oct 29 13:48:02 2018

@author: Hannah Choi
"""
"""
This code finds hierarchy scores of cortical & thalamic areas via iterations, based on the cortico-thalamic (CT) connections 
in addition to the cortico-cortical (CC) AND the thalamo-cortical (TC) connectivity data.
"""


# In[]: Set input and output directories 
    
input_dir = r'./Input/'     # Directory with the file "CC_TC_CT_clusters.xlsx" & "CT_sourcelayer_FFB.xls"
output_dir = r'./Output/'   # Directory to save the ouputs from the experimental data

''' ATTENTION! Change the "df_cortex" accordingly in func_unsupervised_TC as well! '''
CreConf = 1                 # 1 if using CC hierarhcy with Cre-confidence; 0 if not

# In[]: Read the excel file with source-target-creline pairs and their cluster numbers. 
    
xls=pd.ExcelFile(input_dir+"CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'CC+TC clusters with MD avged')
df=df[(df.hemi == "ipsi")&(df.creline != "C57BL/6J / Emx1")&(df.target != "VISC")&(df.source != "SSp-un")&(df.target != "SSp-un")&(df.source != "VISC")]

# In[]: Cortical and thalamic hierarchy from CC+TC iteration

if CreConf == 0:
    h_CC_TC=pd.read_excel(output_dir+"TC_CCnoconf_iter.xls")
elif CreConf == 1:
    h_CC_TC=pd.read_excel(output_dir+"TC_CCconf_iter.xls")
    
C_areas = h_CC_TC[h_CC_TC["CortexThalamus"]=="C"]["areas"].unique()
T_areas = h_CC_TC[h_CC_TC["CortexThalamus"]=="T"]["areas"].unique()

# In[]: 2 clusters (FF vs FB) based on LDA of thalamo-cortical source-target pairs (dfCT)

dfCT = pd.read_excel(input_dir+"CT_sourcelayer_FFB.xls") 
dfCT = dfCT[['source','target','FFB_LDA']]
dfCT = dfCT.rename(columns={"FFB_LDA":"ffb"})

# In[]: 9 clusters of thalamo-cortical source-line-target pairs (dfTC)

dfTC=pd.read_excel(output_dir+"inputexpanded_TC9.xls")
dfTC=dfTC[['source','target','clu','ffb_c']]
dfTC = dfTC.rename(columns={"ffb_c":"ffb"})

# In[]: Merge dataframes of C-T and T-C connections 

dfVT = pd.concat([dfTC, dfCT], ignore_index=True)

source_areas = dfVT["source"].unique()
target_areas = dfVT["target"].unique()

num_TC = np.shape(dfTC)[0]
num_CT = np.shape(dfCT)[0]

dfVT.to_excel(output_dir+'inputexpanded_TC9CT2.xls')

# In[ ]: Find initial hierarchy scores of thalamic areas (21)

n_T_areas=len(T_areas) 

hrs=range(0,n_T_areas)
hrt=range(0,n_T_areas)
hrc=range(0,n_T_areas)

for i in range(0,n_T_areas):
    hrs[i]=-np.mean(dfVT[dfVT.source == T_areas[i]].ffb)
    if len(dfVT[dfVT.target == T_areas[i]]) != 0:
        hrt[i]=np.mean(dfVT[dfVT.target == T_areas[i]].ffb)
    else:
        hrt[i]=0
        
    hrc[i]=(hrs[i]+hrt[i])/2

data=[T_areas,hrc]
data=np.transpose(data)
columns = ['areas','h']
dfiT = pd.DataFrame(data,columns=columns)
dfiT.head()


# In[ ]: Iterate thalamic + cortical hierarhcy scores

n_iter = 20

##########################################################
'''Load results from CC hierarchy'''

if CreConf == 0:
    dfiC = pd.read_excel(output_dir+"CC_noconf_iter.xls")
elif CreConf == 1:
    dfiC = pd.read_excel(output_dir+"CC_conf_iter.xls")
##########################################################
   
dfiC['h'] = dfiC[n_iter]
dfVC = pd.read_excel(output_dir+"inputexpanded_CC9.xls")
dfVC = dfVC[["source","target","ffb_nc","ffb_c"]]

if CreConf == 0:
    dfVC["ffb"] = dfVC["ffb_nc"]
elif CreConf == 1:
    dfVC["ffb"] = dfVC["ffb_c"]
    
dfiT = dfiT[["areas","h"]]
dfiC = dfiC[["areas","h"]]
dfVT = dfVT[["source","target","ffb"]]
dfVC = dfVC[["source","target","ffb"]]


##########################################################
'''Iterations ''' 
hr_iter = iterativeTCCT(dfiC, dfVC, dfiT, dfVT, n_iter)
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
#    fig.savefig(output_dir+"TCCT_CCnoconf_iter.pdf", bbox_inches='tight')
#elif CreConf == 1:
#    fig.savefig(output_dir+"TCCT_CCconf_iter.pdf", bbox_inches='tight')
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
#    f.savefig(output_dir+"TCCT_init_vs_fin_CCnoconf.pdf", bbox_inches='tight')
#elif CreConf == 1:
#    f.savefig(output_dir+"TCCT_init_vs_fin_CCconf.pdf", bbox_inches='tight')
##########################################################


'''Save hierarchy scores before and after iteration'''
for i_area in range(0,n_area):
    if hr_iter['areas'][i_area] in list(dfiC['areas']):
        hr_iter.loc[i_area,'CortexThalamus'] = 'C'
    else:
        hr_iter.loc[i_area,'CortexThalamus'] = 'T'

hr_iter = hr_iter[['areas','CortexThalamus', 0,n_iter] ]  

if CreConf == 0:
    hr_iter.to_excel(output_dir+'TCCT_CCnoconf_iter.xls')
elif CreConf == 1:
    hr_iter.to_excel(output_dir+'TCCT_CCconf_iter.xls')

 # In[]: Print out global hierarchy scores for the TC + CT connectivity data before and after iteration

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


print('hg of CC+TC+CT before iterate cortex & thalamus='+str(hg_TCCT_init))
print('hg of CC+TC+CT iterate cortex='+str(hg_cortex_TCCT_iter)) 
print('hg of CC+TC+CT iterate cortex & thalamus='+str(hg_TCCT_iter))

'''Save global hierarchy scores'''
newDF= pd.DataFrame([])
newDF=newDF.append(pd.DataFrame({'hg_TCCT_init':hg_TCCT_init, 'hg_cortex_TCCT_iter':hg_cortex_TCCT_iter,
                                 'hg_TCCT_iter':hg_TCCT_iter},index=[0]))
newDF.to_excel(output_dir+'ghs_TCCT.xls')
