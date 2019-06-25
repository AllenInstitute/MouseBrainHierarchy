import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


"""
Created on Mon Jun 24 16:04:07 2019

@author: Hannah Choi
"""
"""
This code is called in main run files to refine hierarchy scores via iterations.
"""

# In[]: Cortico-cortical iterations

def iterativeCC(dfi,dfV,n_iter):

    df=dfi.sort_values('areas')
    
    hrc0=df["hrc"]
    hr0=df["hr"]
    
    areas=np.asarray(df.areas)
    n_area=np.shape(dfi)[0]
    
    hrc_iter_vals=np.zeros([n_area,n_iter])
    hr_iter_vals=np.zeros([n_area,n_iter])
    
    hrc_iter_vals[:,0]=hrc0
    hr_iter_vals[:,0]=hr0
    
    hrc_iter=pd.DataFrame(hrc_iter_vals)
    hrc_iter['areas'] = areas
    cols = hrc_iter.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    hrc_iter = hrc_iter[cols]  
    
    hr_iter=pd.DataFrame(hr_iter_vals)
    hr_iter['areas'] = areas
    cols = hr_iter.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    hr_iter = hr_iter[cols]

    
    for i_iter in range(1,n_iter+1):
        print i_iter
        for i_area in range(0,n_area):    
            current_target=np.array(dfV[dfV.source == areas[i_area]].target)
            current_source=np.array(dfV[dfV.target == areas[i_area]].source)
            
            hrs_ffbnc=np.array(dfV[dfV.source == areas[i_area]].ffb_nc)
            hrt_ffbnc=np.array(dfV[dfV.target == areas[i_area]].ffb_nc) 
            hrs_ffbc=np.array(dfV[dfV.source == areas[i_area]].ffb_c)
            hrt_ffbc=np.array(dfV[dfV.target == areas[i_area]].ffb_c)  
            cfs=np.array(dfV[dfV.source == areas[i_area]].conf)
            cft=np.array(dfV[dfV.target == areas[i_area]].conf)
            
            # i_area as source; compute influences from its targets
            hrc_asSource=np.zeros([np.size(current_target),1])
            hr_asSource=np.zeros([np.size(current_target),1])
            for j_area in range(0,np.size(current_target)):
                if any(hrc_iter.areas==current_target[j_area])==True:
                    hc_current_target = hrc_iter[hrc_iter.areas==current_target[j_area]][i_iter-1]
                    h_current_target = hr_iter[hr_iter.areas==current_target[j_area]][i_iter-1]
                    hrc_asSource[j_area]=cfs[j_area]*(-hc_current_target+hrs_ffbc[j_area])
                    hr_asSource[j_area]=(-h_current_target+hrs_ffbnc[j_area])
                else:
                    hrc_asSource[j_area]=0
                    hr_asSource[j_area]=0
                 
            
            # i_area as target; compute influences from its sources
            hrc_asTarget=np.zeros([np.size(current_source),1])
            hr_asTarget=np.zeros([np.size(current_source),1])
            for k_area in range(0,np.size(current_source)):
                if any(hrc_iter.areas==current_source[k_area])==True:
                    hc_current_source = hrc_iter[hrc_iter.areas==current_source[k_area]][i_iter-1]
                    h_current_source = hr_iter[hr_iter.areas==current_source[k_area]][i_iter-1]
                    hrc_asTarget[k_area]=cft[k_area]*(hc_current_source+hrt_ffbc[k_area])
                    hr_asTarget[k_area]=(h_current_source+hrt_ffbnc[k_area])
                else:
                    hrc_asTarget[k_area]=0
                    hr_asTarget[k_area]=0
                    
        
            hrc_iter.loc[hrc_iter.areas==areas[i_area],i_iter]=0.5*(-np.mean(hrc_asSource)+np.mean(hrc_asTarget))
            hr_iter.loc[hr_iter.areas==areas[i_area],i_iter]=0.5*(-np.mean(hr_asSource)+np.mean(hr_asTarget))

        hrc_iter[:][i_iter]=hrc_iter[:][i_iter]-np.mean(hrc_iter[:][i_iter])
        hr_iter[:][i_iter]=hr_iter[:][i_iter]-np.mean(hr_iter[:][i_iter])
    
    return hr_iter, hrc_iter



# In[]: Thalamo-cortical iterations

def iterativeTC(dfiC, dfVC, dfiT, dfVT, n_iter):
    
    dfi = dfiT.append(dfiC)
    dfV = dfVT.append(dfVC)

    df=dfi.sort_values('areas')
    
    hrc0=df["h"]
    
    areas=np.asarray(df.areas)
    n_area=len(areas) #np.shape(dfi)[0]
    
    hr_iter_vals=np.zeros([n_area,n_iter])
    hr_iter_vals[:,0]=hrc0
    
    hr_iter=pd.DataFrame(hr_iter_vals)
    hr_iter['areas'] = areas
    cols = hr_iter.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    hr_iter = hr_iter[cols]
    hr_iter.head()

    
    for i_iter in range(1,n_iter+1):
        print i_iter
        for i_area in range(0,n_area):    
            current_target=np.array(dfV[dfV.source == areas[i_area]].target)
            current_source=np.array(dfV[dfV.target == areas[i_area]].source)
   
            hrs_ffbc=np.array(dfV[dfV.source == areas[i_area]].ffb)
            hrt_ffbc=np.array(dfV[dfV.target == areas[i_area]].ffb)  
            
            # i_area as source; compute influences from its targets       
            if len(current_target)>0:
                hr_asSource=np.zeros([np.size(current_target),1])
                for j_area in range(0,np.size(current_target)):
                    if any(hr_iter.areas==current_target[j_area])==True:
                        h_current_target = hr_iter[hr_iter.areas==current_target[j_area]][i_iter-1]
                        hr_asSource[j_area]=(-h_current_target+hrs_ffbc[j_area])
                    else:
                        hr_asSource[j_area]=0  
            else:
                hr_asSource = 0
                
            # i_area as target; compute influences from its sources
            if len(current_source)>0:
                hr_asTarget=np.zeros([np.size(current_source),1])
                for k_area in range(0,np.size(current_source)):
                    if any(hr_iter.areas==current_source[k_area])==True:
                        h_current_source = hr_iter[hr_iter.areas==current_source[k_area]][i_iter-1]
                        hr_asTarget[k_area]=(h_current_source+hrt_ffbc[k_area])
                    else:
                        hr_asTarget[k_area]=0
            else:
                hr_asTarget = 0
                
            hr_iter.loc[hr_iter.areas==areas[i_area],i_iter]=0.5*(-np.mean(hr_asSource)+np.mean(hr_asTarget))


        hr_iter[:][i_iter]=hr_iter[:][i_iter]-np.mean(hr_iter[:][i_iter])
    
    return hr_iter


# In[]: Thalamo-cortical & cortico-thalamic iterations

def iterativeTCCT(dfiC, dfVC, dfiT, dfVT, n_iter):
    
    dfi = dfiT.append(dfiC)
    dfV = dfVT.append(dfVC)

    df=dfi.sort_values('areas')
    
    hrc0=df["h"]
    
    areas=np.asarray(df.areas)
    n_area=len(areas) #np.shape(dfi)[0]
    
    hr_iter_vals=np.zeros([n_area,n_iter])
    hr_iter_vals[:,0]=hrc0
    
    hr_iter=pd.DataFrame(hr_iter_vals)
    hr_iter['areas'] = areas
    cols = hr_iter.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    hr_iter = hr_iter[cols]
    hr_iter.head()

    
    for i_iter in range(1,n_iter+1):
        print i_iter
        for i_area in range(0,n_area):    
            current_target=np.array(dfV[dfV.source == areas[i_area]].target)
            current_source=np.array(dfV[dfV.target == areas[i_area]].source)
   
            hrs_ffbc=np.array(dfV[dfV.source == areas[i_area]].ffb)
            hrt_ffbc=np.array(dfV[dfV.target == areas[i_area]].ffb)  
            
            # i_area as source; compute influences from its targets       
            if len(current_target)>0:
                hr_asSource=np.zeros([np.size(current_target),1])
                for j_area in range(0,np.size(current_target)):
                    if any(hr_iter.areas==current_target[j_area])==True:
                        h_current_target = hr_iter[hr_iter.areas==current_target[j_area]][i_iter-1]
                        hr_asSource[j_area]=(-h_current_target+hrs_ffbc[j_area])
                    else:
                        hr_asSource[j_area]=0  
            else:
                hr_asSource = 0
                
            # i_area as target; compute influences from its sources
            if len(current_source)>0:
                hr_asTarget=np.zeros([np.size(current_source),1])
                for k_area in range(0,np.size(current_source)):
                    if any(hr_iter.areas==current_source[k_area])==True:
                        h_current_source = hr_iter[hr_iter.areas==current_source[k_area]][i_iter-1]
                        hr_asTarget[k_area]=(h_current_source+hrt_ffbc[k_area])
                    else:
                        hr_asTarget[k_area]=0
            else:
                hr_asTarget = 0
                
            hr_iter.loc[hr_iter.areas==areas[i_area],i_iter]=0.5*(-np.mean(hr_asSource)+np.mean(hr_asTarget))


        hr_iter[:][i_iter]=hr_iter[:][i_iter]-np.mean(hr_iter[:][i_iter])
    
    return hr_iter