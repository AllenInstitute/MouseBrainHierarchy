from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_CC import fit
from IterativeMethod import iterativeCC

"""
Created on Sat Oct 27 17:45:37 2018

@author: Hannah Choi
"""
"""
This code finds hierarchy scores of cortical areas by maximizing the global hierarchy score 
followed by interations, based on the cortico-cortical (CC) connectivity data.
"""

# In[]: Set input and output directories 

input_dir = r'./Input/'        # Directory with the file "CC_TC_CT_clusters.xlsx"
output_dir = r'./Output/'      # Directory to save the ouputs from the experimental data

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

num_clu = len(dfV["clu"].unique())


# In[]: Map 9 clusters of cortico-cortical connections to either of FF/FB directions, by maximizing the global hierarchy score

hierarchy_vals = fit(dfV)
    
jmax_raw, jmax = np.argmax(hierarchy_vals, axis=0)
jmax_raw_val = hierarchy_vals[jmax_raw][0]
jmax_val = hierarchy_vals[jmax][1]
logging.debug("RESULTS")
n = len(dfV.clu.unique())
logging.debug("(jmax_raw, val) = ({:0{n}b}, {:.3f})".format(jmax_raw, jmax_raw_val, n=n))
logging.debug("(jmax,     val) = ({:0{n}b}, {:.3f})".format(jmax, jmax_val, n=n))

results = dict(jmax=bin(2**n+jmax),
               jmax_val=jmax_val,
               jmax_raw=bin(2**n+jmax_raw),
               jmax_raw_val=jmax_raw_val)
hrc_original = jmax_val         # with Cre-confidence
hr_original = jmax_raw_val      # without Cre-confidence

# In[ ]: Define functions needed

c0=2**num_clu

def ffb_c (cls):
    """Direction of each cluster with Cre-confidence"""
    b=(bin(c0+jmax)[-(cls)])
    return (2*int(b)-1)    # or -(2*int(b)-1), depending on the final direction of the hierarchy


def ffb_nc (cls):
    """Direction of each cluster without Cre-confidence"""
    b=(bin(c0+jmax_raw)[-(cls)])
    return -(2*int(b)-1)

def cre_confidence1(df):
    """Returns confidence of Cre lines"""
    func = lambda x: 1 - np.abs(x.mean())
    return df.groupby('creline')['ffb_c'].transform(func)


def cre_confidence2(df):
     """Returns an alternative confidence of Cre lines"""   
     cre_list = df.creline.unique()
     areas = df.source.unique()
     n=len(areas)
     confnew =pd.DataFrame(np.zeros(len(df.index)),index=df.index, columns=['conf'])
 
     for kk in range(0,len(cre_list)):
         df_cre = df[df.creline == cre_list[kk]]
         print cre_list[kk]
         count_sym = 0
         count_sameffb = 0
         for ii in range(0,n):
             for jj in range(ii+1,n):
                 ij_ffb = np.array(df_cre[(df.source == areas[ii])&(df.target == areas[jj])].ffb_c)
                 ji_ffb = np.array(df_cre[(df.source == areas[jj])&(df.target == areas[ii])].ffb_c)
                 if len(ij_ffb)==1 and len(ji_ffb)==1:
                     count_sym = count_sym+1
                     if ij_ffb == ji_ffb:
                         count_sameffb = count_sameffb+1 
         confnew[df.creline == cre_list[kk]] = 1-count_sameffb/count_sym
     return confnew


def hrf (area):
    '''Hierarchy score of each area without Cre-confidence'''
    return (-np.mean(dfV[dfV.source == area].ffb_nc)+np.mean(dfV[dfV.target == area].ffb_nc))/2

def hrcf (area):
    '''Hierarchy score of each area with Cre-confidence'''
    return (-np.mean(dfV[dfV.source == area].ffb_c*dfV[dfV.source == area].conf)+np.mean(dfV[dfV.target == area].ffb_c*dfV[dfV.target == area].conf))/2



# In[ ]: Produce an expanded data frame (dfV) with Cre confidence values, confidence-weighted & unweighted FF/FB, hierarchy values as source & target for each pair of CC connections

dfV.loc[:,"ffb_c"]=dfV["clu"].apply(ffb_c)
dfV.loc[:,"ffb_nc"]=dfV["clu"].apply(ffb_nc)
dfV.loc[:, "conf"] = cre_confidence1(dfV)
dfV.loc[:,"hrc_s"]=dfV["source"].apply(hrcf)
dfV.loc[:,"hrc_t"]=dfV["target"].apply(hrcf)
dfV.loc[:,"hr_s"]=dfV["source"].apply(hrf)
dfV.loc[:,"hr_t"]=dfV["target"].apply(hrf)

dfV.to_excel(output_dir+'inputexpanded_CC9.xls')

# In[ ]: Find the initial hierarchy score of each of cortical areas (37 cortical areas)

areas = source_areas
n_areas=len(areas)
hrs=range(0,n_areas)
hrt=range(0,n_areas)
hr=range(0,n_areas)
hrc=range(0,n_areas)

for i in range(0,n_areas):
    hrs[i]=-np.mean(dfV[dfV.source == areas[i]].ffb_nc)
    hrt[i]=np.mean(dfV[dfV.target == areas[i]].ffb_nc)
    hr[i]=(hrs[i]+hrt[i])/2
    hrc[i]=0.5*(-np.mean(dfV[dfV.source == areas[i]].ffb_c*dfV[dfV.source == areas[i]].conf)+np.mean(dfV[dfV.target == areas[i]].ffb_c*dfV[dfV.target == areas[i]].conf))


data=[areas,hrc,hr]
data=np.transpose(data)
columns = ['areas','hrc','hr']
dfi = pd.DataFrame(data,columns=columns)
dfi.head()

# dfi.to_excel(output_dir+'initialhierarchy_CC9.xls')

# In[] : Iterate cortical hierarhcy scores, to refine the hierarhcy levels

##########################################################
'''Iterations '''

n_iter = 20
hr_iter, hrc_iter = iterativeCC(dfi,dfV,n_iter)
##########################################################

iteration=np.arange(0,n_iter+1,1)
n_area=len(source_areas)


##########################################################
""" Figure of hierarchy score iterations with Cre-confidence """

fig,ax=plt.subplots()
for i in range(0,n_area):
    y=np.squeeze(np.asarray(hrc_iter[hrc_iter.areas==areas[i]].ix[:,1::1]))
    ax.plot(iteration,y)

ax.set_xlim([0, n_iter])
ax.set_xticks(np.arange(0, n_iter, step=5))
ax.set_xlabel('iter')
ax.set_ylabel('hierarchy value')
ax.set_title('confidence adjusted')
plt.show()
#fig.savefig(output_dir+"CC_conf_iter.pdf", bbox_inches='tight')
##########################################################


##########################################################
""" Figure of hierarchy score iterations without Cre-confidence """

fig,ax=plt.subplots()
for i in range(0,n_area):
    y=np.squeeze(np.asarray(hr_iter[hr_iter.areas==areas[i]].ix[:,1::1]))
    ax.plot(iteration,y)

ax.set_xlim([0, n_iter])
ax.set_xticks(np.arange(0, n_iter, step=5))
ax.set_xlabel('iter')
ax.set_ylabel('hierarchy value')
ax.set_title('confidence not adjusted')
plt.show()
#fig.savefig(output_dir+"CC_noconf_iter.pdf", bbox_inches='tight')
##########################################################

##########################################################
""" Figure showing correlation between hierarchy scores before & after iterations with Cre-confidence """

hrc_final = hrc_iter[:][n_iter]
hrc_initial = hrc_iter[:][0]    
f = plt.figure()
plt.plot(hrc_initial,hrc_final,'ro')
plt.xlabel('initial hierarchy (conf)')
plt.ylabel('final hierarchy (conf)')
plt.title('r='+str(np.corrcoef(hrc_initial, hrc_final)[0, 1]))
plt.show()    
#f.savefig(output_dir+"CC_init_vs_fin_conf.pdf", bbox_inches='tight')
##########################################################

##########################################################
""" Figure showing correlation between hierarchy scores before & after iterations without Cre-confidence """

hr_final = hr_iter[:][n_iter]
hr_initial = hr_iter[:][0]   
f = plt.figure()
plt.plot(hr_initial,hr_final,'ro')
plt.xlabel('initial hierarchy (noconf)')
plt.ylabel('final hierarchy (noconf)')
plt.title('r='+str(np.corrcoef(hr_initial, hr_final)[0, 1]))
plt.show()
#f.savefig(output_dir+"CC_init_vs_fin_noconf.pdf", bbox_inches='tight')
##########################################################


'''Save hierarchy scores before and after iteration'''
hrc_iter = hrc_iter[['areas',0,n_iter]]
hr_iter = hr_iter[['areas',0,n_iter] ]  
hrc_iter.to_excel(output_dir+'CC_conf_iter.xls')  # Save before & after iteration hierarchy scores with Cre-confidence
hr_iter.to_excel(output_dir+'CC_noconf_iter.xls') # Save before & after iteration hierarchy scores without Cre-confidence
    
 # In[]: Print out global hierarchy scores for the CC connectivity data before and after iteration
 
dfV_temp = dfV[['source','target','clu','ffb_c','ffb_nc','conf']]
dfi_temp = hr_iter[['areas',0,n_iter]]
dfi_temp_conf = hrc_iter[['areas',0,n_iter]]
dfi_temp = dfi_temp.rename(columns={0: "h0", n_iter:"h_iter"})
dfi_temp_conf = dfi_temp_conf.rename(columns={0: "h0", n_iter:"h_iter"})

dfi_t = dfi_temp[['areas','h0']]
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
dfV_temp =dfV_temp.rename(columns={"h0": "hs_init"})
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
dfV_temp =dfV_temp.rename(columns={"h0": "ht_init"})
dfV_temp = dfV_temp.dropna()

dfi_t = dfi_temp[['areas','h_iter']]
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
dfV_temp =dfV_temp.rename(columns={"h_iter": "hs_iter"})
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
dfV_temp =dfV_temp.rename(columns={"h_iter": "ht_iter"})
dfV_temp = dfV_temp.dropna()

#hg_CC_init = np.mean(dfV_temp.ffb_nc*np.sign(dfV_temp.ht_init - dfV_temp.hs_init))
#hg_CC_iter = np.mean(dfV_temp.ffb_nc*np.sign(dfV_temp.ht_iter - dfV_temp.hs_iter))
hg_CC_init = np.mean(dfV_temp.ffb_nc*(dfV_temp.ht_init - dfV_temp.hs_init))
hg_CC_iter = np.mean(dfV_temp.ffb_nc*(dfV_temp.ht_iter - dfV_temp.hs_iter))


dfV_temp = dfV[['source','target','clu','ffb_c','ffb_nc','conf']]
dfi_t= dfi_temp_conf[['areas','h0']]
dfV_temp= dfV_temp.join(dfi_t.set_index('areas'), on ='source')
dfV_temp =dfV_temp.rename(columns={"h0": "hs_init"})
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
dfV_temp =dfV_temp.rename(columns={"h0": "ht_init"})
dfV_temp = dfV_temp.dropna()

dfi_t = dfi_temp_conf[['areas','h_iter']]
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
dfV_temp =dfV_temp.rename(columns={"h_iter": "hs_iter"})
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
dfV_temp =dfV_temp.rename(columns={"h_iter": "ht_iter"})
dfV_temp = dfV_temp.dropna()

#hg_CC_conf_init = np.mean(dfV_temp.ffb_c*np.sign(dfV_temp.ht_init - dfV_temp.hs_init))
#hg_CC_conf_iter = np.mean(dfV_temp.ffb_c*np.sign(dfV_temp.ht_iter - dfV_temp.hs_iter))
hg_CC_conf_init = np.mean(dfV_temp.ffb_c*(dfV_temp.ht_init - dfV_temp.hs_init))
hg_CC_conf_iter = np.mean(dfV_temp.ffb_c*(dfV_temp.ht_iter - dfV_temp.hs_iter))


print('hg of CC cortex='+str(hg_CC_init))
print('hg of CC iterate cortex='+str(hg_CC_iter))

print('hgc of CC cortex='+str(hg_CC_conf_init))
print('hgc of CC iterate cortex='+str(hg_CC_conf_iter))

'''Save global hierarchy scores'''
newDF= pd.DataFrame([])
newDF=newDF.append(pd.DataFrame({'hg_CC_init':hg_CC_init, 'hg_CC_iter':hg_CC_iter,
                                 'hg_CC_conf_init':hg_CC_conf_init, 'hg_CC_conf_iter':hg_CC_conf_iter},index=[0]))
newDF.to_excel(output_dir+'ghs_CC.xls')



