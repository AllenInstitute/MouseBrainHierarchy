from __future__ import division
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

"""
Created on Fri Aug 09 16:19:11 2019

@author: Hannah Choi
"""

"""
This code takes previously generated hierachy scores & global hierarchy score data from the "Results" folder, and make figures.
"""
# In[]:

input_dir = r'./Results/' 
xls1=pd.ExcelFile(input_dir+"hierarchy_summary_CreConf.xlsx")
xls2=pd.ExcelFile(input_dir+"gh_comparison_shuffled_CreConf.xlsx")

# In[]:
''' Hierarchy of all cortical & thalamic regions'''

df=pd.read_excel(xls1,'hierarchy_all_regions')
df = df.sort_values('CC+TC+CT iterated')

areas = df['areas']
areas_vec = np.arange(0,len(areas))
hs = df['CC+TC+CT iterated']

f = plt.figure()
plt.plot(hs,areas_vec,'ro')
plt.xlabel('hierarchy score')
plt.ylabel('areas')
plt.xlim(-1.1,1)
plt.yticks(areas_vec, areas)
plt.title('Hierarchy based on CC+TC+CT data')
f.set_size_inches(6, 10.5)
plt.show()    

# In[]:
''' Inter-module hierarchy of VisualMedial module'''

df=pd.read_excel(xls1,'hierarchy_visualmedial(VIS)')
df = df.sort_values('CC+TC+CT iterated')

areas = df['areas']
areas_vec = np.arange(0,len(areas))
hs = df['CC+TC+CT iterated']

f = plt.figure()
plt.plot(hs,areas_vec,'ro')
plt.xlabel('hierarchy score')
plt.ylabel('areas')
plt.xlim(-1.1,1)
plt.yticks(areas_vec, areas)
plt.title('Viual+Medial module hierarchy')
f.set_size_inches(4, 6)
plt.show()    

# In[]:
''' Intra-module cortical hierarchy'''

df=pd.read_excel(xls1,'hierarchy_intermodule')
df = df.sort_values('CC+TC+CT iterated')

areas = df['areas']
areas_vec = np.arange(0,len(areas))
hs = df['CC+TC+CT iterated']

f = plt.figure()
plt.plot(hs,areas_vec,'ro')
plt.xlabel('hierarchy score')
plt.ylabel('areas')
plt.xlim(-1.1,1)
plt.yticks(areas_vec, areas)
plt.title('Inter-module hierarchy')
f.set_size_inches(4, 6)
plt.show()  

# In[]:
''' Global hierarchy scores: comparison to shuffled data'''

df_CC = pd.read_excel(xls2,'CC')
CC_shuffled = df_CC['iter shuffled hg (cortex)']
CC_data = df_CC['iter data hg (cortex)'][0]
zscore_CC = (CC_data-np.mean(CC_shuffled))/np.std(CC_shuffled) 

df_TC = pd.read_excel(xls2,'CC+TC')
TC_shuffled = df_TC['iter shuffled hg (cortex+thal)']
TC_data = df_TC['iter data hg (cortex+thal)'][0]
zscore_TC = (TC_data-np.mean(TC_shuffled))/np.std(TC_shuffled) 

df_CT = pd.read_excel(xls2,'CC+TC+CT')
CT_shuffled = df_CT['iter shuffled hg (cortex+thal)']
CT_data = df_CT['iter data hg (cortex+thal)'][0]
zscore_CT = (CT_data-np.mean(CT_shuffled))/np.std(CT_shuffled) 

fig,ax=plt.subplots()
ax.hist(CC_shuffled, bins=25, label='CC',color='red')
ax.axvline(x=CC_data,linestyle='--',color='red')
ax.hist(TC_shuffled, bins=25, label='CC+TC',color='orange')
ax.axvline(x=TC_data,linestyle='--',color='orange')
ax.hist(CT_shuffled, bins=25, label='CC+TC+CT')
ax.axvline(x=CT_data,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
#ax.set_title('z-score(CC)='+str(zscore_CC),'z-score(CC+TC)='+str(zscore_TC),'z-score(CC+TC+CT)='+str(zscore_CT))
ax.legend(loc='upper right')