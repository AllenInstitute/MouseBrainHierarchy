from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IterativeMethod import iterativeTCCT

"""
Created on Mon Jun 24 16:04:07 2019

@author: Hannah Choi
"""
"""
This code finds intra- or inter-module hierarchies of cortical network, 
based on the cortico-cortical (CC), thalamo-cortical (TC), and cortico-thalamic (CT) connectivity data.
"""

# In[]:

CreConf = 1                    # 1 if using CC hierarhcy with Cre-confidence; 0 if not

input_dir = r'./Input/'        # Directory with the file "CC_TC_CT_clusters.xlsx" & "clustermapping.xlsx"
output_dir = r'./Output/module/'      # Directory to save the ouputs from the experimental data


# In[]: Define the module. 

module = 'VisualMedial' #'VisualMedial' #'inter_predefined'

# In the paper, we used the following: 'VisualMedial' & 'inter_predefined'
# Possible modules: 'VisualMedial', 'Visual', 'Medial', 'Auditory', 'Somatomotor', 'PFC', 'Lateral', 'inter_predefined', 'inter'

# In[]: 9 clusters of cortico-cortico & thalamo-courtical source-line-target pairs
    
xls=pd.ExcelFile(input_dir+"CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'CC+TC clusters with MD avged')
df=df[(df.hemi == "ipsi")&(df.creline != "C57BL/6J / Emx1")&(df.target != "VISC")&(df.source != "SSp-un")&(df.target != "SSp-un")&(df.source != "VISC")]

list_module = df["Cortical Target Module"].unique()
clu_ffb=pd.read_excel(input_dir+"clustermapping.xlsx")

# In[]: If inter-module, change all the target & source area names to the module name

if (module == 'inter') or (module == 'inter_predefined'): 
    for i_module in range(0,len(list_module)):
        df.loc[df["Cortical Target Module"] == list_module[i_module],'target']= list_module[i_module]
        df.loc[df["Cortical Source Module"] == list_module[i_module],'source'] = list_module[i_module] 


# In[]: Cortical and thalamic hierarchy from CC+TC iteration

if CreConf == 1:
    h_CC_TC=pd.read_excel(output_dir+"TC_CCconf_iter_"+module+".xls")
elif CreConf == 0:
    h_CC_TC=pd.read_excel(output_dir+"TC_CCnoconf_iter_"+module+".xls")
C_areas = h_CC_TC[h_CC_TC["CortexThalamus"]=="C"]["areas"].unique()
T_areas = h_CC_TC[h_CC_TC["CortexThalamus"]=="T"]["areas"].unique()



# In[]: 2 clusters (FF/FB) based on LDA of cortico-thalamic source-line-target pairs (dfCT)

dfCT = pd.read_excel(input_dir+"CT_sourcelayer_FFB.xls") 
dfCT=dfCT[['source','target','FFB_LDA','Cortical Source Module']]
dfCT = dfCT.rename(columns={"FFB_LDA":"ffb"})

# In[]: For intra medial, select only the areas within the chosen module for cortico-thalamic connections (dfCT)

if (module != 'inter') and (module != 'inter_predefined'):     

    if module == 'VisualMedial':
        dfCT.loc[((dfCT["Cortical Source Module"] == "Visual")|(dfCT["Cortical Source Module"] == "Medial")),"Cortical Source Module"]='VisualMedial'
        dfCT = dfCT[(dfCT["Cortical Source Module"] == module)]
        dfCT = dfCT[(dfCT["source"]!='RSPagl')&(dfCT["target"]!='RSPagl')&(dfCT["source"]!='RSPd')&(dfCT["target"]!='RSPd')
        &(dfCT["source"]!='RSPv')&(dfCT["target"]!='RSPv')]
    else:
        dfCT = dfCT[(dfCT["Cortical Source Module"] == module)]

# In[]: For intra medial, select only the areas within the chosen module for thalamo-cortical source-line-target pairs (dfTC)
 
dfTC=pd.read_excel(output_dir+"inputexpanded_TC9_"+module+".xls")
dfTC = dfTC.rename(columns={"ffb_c":"ffb"})

if (module != 'inter') and (module != 'inter_predefined'):    
   
    if module == 'VisualMedial':
        dfTC.loc[((dfTC["Cortical Target Module"] == "Visual")|(dfTC["Cortical Target Module"] == "Medial")),"Cortical Target Module"]='VisualMedial'
        dfTC = dfTC[(dfTC["Cortical Target Module"] == module)]
        dfTC = dfTC[(dfTC["source"]!='RSPagl')&(dfTC["target"]!='RSPagl')&(dfTC["source"]!='RSPd')&(dfTC["target"]!='RSPd')
        &(dfTC["source"]!='RSPv')&(dfTC["target"]!='RSPv')]
    else:
        dfTC = dfTC[(dfTC["Cortical Target Module"] == module)]
        
# In[]: Merge dataframes of CT and TC connections 

if (module == 'inter') or (module == 'inter_predefined'):
    dfCT = dfCT[['Cortical Source Module','target','ffb']]
    dfCT = dfCT.rename(columns={"Cortical Source Module": "source"})
else:
    dfCT = dfCT[['source','target','ffb']]    
    
dfTC = dfTC[['source','target','ffb']]

dfVT = pd.concat([dfTC, dfCT], ignore_index=True)

# In[]: Produce expanded data frame with  FF/FB, hierarchy values as source & target for each pair of TC+CT connections

source_areas = dfVT["source"].unique()
target_areas = dfVT["target"].unique()

num_TC = np.shape(dfTC)[0]
num_CT = np.shape(dfCT)[0]

dfVT.to_excel(output_dir+'inputexpanded_TC9CT2_'+module+'.xls')


# In[ ]: Find initial hierarchy scores of thalamic areas (21)

n_T_areas=len(T_areas) # 21 thalamic regions

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
#dfiT.to_excel(output_dir+'initialhierarchy_TC9CT2_'+module+'.xls')
#dfiT.head()


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



dfiT = dfiT[["areas","h"]]
dfiC = dfiC[["areas","h"]]
dfVT = dfVT[["source","target","ffb"]]
dfVC = dfVC[["source","target","ffb"]]

hr_iter = iterativeTCCT(dfiC, dfVC, dfiT, dfVT, n_iter)


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
#fig.savefig(output_dir+"TCCT_CCnoconf_iter_"+module+".pdf", bbox_inches='tight')
##########################################################

##########################################################
""" Figure showing correlation between hierarchy scores before & after iterations"""
hr_final = hr_iter[:][n_iter]
hr_initial = hr_iter[:][0]   
f = plt.figure()
plt.plot(hr_initial,hr_final,'ro')
plt.xlabel('initial hierarchy')
plt.ylabel('final hierarchy')
plt.title('r='+str(np.corrcoef(hr_initial, hr_final)[0, 1]))
plt.show()
#f.savefig(output_dir+"TCCT_init_vs_fin_CCnoconf_"+module+".pdf", bbox_inches='tight')
##########################################################

##########################################################
'''Save hierarchy scores before and after iteration'''

for i_area in range(0,n_area):
    if hr_iter['areas'][i_area] in list(dfiC['areas']):
        hr_iter.loc[i_area,'CortexThalamus'] = 'C'
    else:
        hr_iter.loc[i_area,'CortexThalamus'] = 'T'

hr_iter = hr_iter[['areas','CortexThalamus', 0,n_iter] ]  
hr_iter_save = hr_iter[(hr_iter.CortexThalamus=='C')]

if CreConf == 1:
    hr_iter_save.to_excel(output_dir+'TCCT_CCconf_iter_'+module+'.xls')
elif CreConf == 0:
    hr_iter_save.to_excel(output_dir+'TCCT_CCnoconf_iter_'+module+'.xls')

##########################################################

 # In[]: Print out global hierarchy scores of CC+TC+CT connectivity data

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

#########################################
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
newDF.to_excel(output_dir+'ghs_TCCT_'+module+'.xls')

