# CorticalThalamic_Hierarchy

* Hierarchy searching algorithm based on layer-projection patterns of cortico-cortical, thalamo-cortical, and cortico-thalamic connections in Cre-transgenic mice
* Codes used in in Harris et al (2019) “The organization of corticocortical, thalamocortical, and corticothalamic connections by layer and cell class in the mouse brain”
* Written by Hannah Choi (hannahc@alleninstitute.org), 3/11/2019.
************************************************************************************
These Python codes build the hierarchy of the mouse cortical and thalamic regions based on clustered layer-projection patterns by maximizing the self-consistency represented by global hierarchy score, followed by iterative corrections. The hierarchy is first constructed based on the cortico-cortical connections only, which can be updated by taking thalamo-cortical and cortico-thalamic connections in to account. 

## Level of Support
We are releasing this code to the public as a tool we expect others to use. Questions concerning bugs and related issues are welcomed. We expect to address them promptly, pull requests will vetted by our staff before inclusion.

## Cortico-cortical hierarchy

### ```run_CC.py```
  * The main hierarchy generating code based on cortico-cortical connection info
  * Calls:
    - ```func_unsupervised_CC.py```
    - ```IterativeMethod.py```
  *	Input files:
    -	```CC_TC_CT_clusters.xlsx```
  *	Output files:
    -	```CC_conf_iter.xls``` & ```CC_noconf_iter.xls```which include hierarchy scores of all cortical regions before & after iterations with        or without Cre-confidence. This filed is read in ```run_TC.py``` and ```func_unsupervised_TC.py```.
    -	```inputexpanded_CC9.xls```, which include information on CC connections such as sources, targets, and FF/FB mapping. This file is          read in ```run_TC.py``` for iteration with TC connections. 
    -	```ghs_CC.xls``` which contains global hierarchy scores of TC connectivity data

### ```func_unsupervised_CC.py```
  * Contains functions needed for mapping CC clusters to FF/FB by maximizing the global hierarchy score of the given cortico-cortical       connectivity data.

### ```IterativeMethod.py```
  * Iterate hierarchy scores of all cortical regions to refine the hierarchy
  *	Import ```iterativeCC```

## Adding thalamo-cortical hierarchy
The cortico-cortical hierarchy can be updated by including the thalamo-cortical connectivity. Since the initial thalamic hierarchical positions are determined based on CC hierarchy scores of their target cortical regions, ```run_CC.py``` should be called before running ```run_TC.py```.

### ```run_TC.py```
  * The main hierarchy generating code based on thalamo-cortical connection info.
  * Calls:
    -	```func_unsupervised_TC.py```
    -	```IterativeMethod.py```
  *	Input files:
    -	```CC_TC_CT_clusters.xlsx```
    -	```CC_conf_iter.xls``` or ```CC_noconf_iter.xls```
    - ```inputexpanded_CC9.xls```
  *	Output files:
   - ```TC_CCconf_iter.xls``` or ```TC_CCnoconf_iter.xls```, which include hierarchy scores of all cortical and thalamic regions before & after iterations. 
   - ```inputexpanded_TC9.xls```, which include information on TC connections such as sources, targets, and FF/FB mapping. This file is          read in ```run_TCCT.py``` for iteration with TC connections. 
   - ```ghs_TC.xls``` which contains global hierarchy scores of TC connectivity data

### ```func_unsupervised_TC.py``` 
  *	Contains functions needed for mapping TC clusters to FF/FB by maximizing the global hierarchy score of the given thalamo-cortical       connectivity data. As the initial thalamic hierarchical positions and the FF/FB mapping are based on the cortical hierarchy             constructed from “run_CC.py”, the cortical hierarchy scores from CC connectivity should be read & changed accordingly (in               ```df_cortex=pd.read_excel(r'./Output/CC_conf_iter.xls')``` )  

###	```IterativeMethod.py```
  *	Iterate hierarchy scores of all cortical & thalamic regions to refine the hierarchy
  *	Import ```iterativeTC```

##	Adding cortico-thalamic hierarchy
The cortico-cortical +thalamo-cortical hierarchy can be further updated by including cortico-thalamic connectivity.

### ```run_CT.py```
  *	Calls:
    -	```IterativeMethod.py```
  * Input files:
    -	```CC_TC_CT_clusters.xlsx```
    - ```CT_sourcelayer_FFB.xls``` has FF/FB assignment for each source-target pair of cortico-thalamic connections based on 1) Linear           Discriminant Analysis (LDA) results on L5/L6 source NPV, and 2) FF/FB-ness for each cortico-thalamic connection ```run_CC.py``` with Cre-confidence.  If not using Cre-confidence for CC hierarchy, use the newly obtained FF/FB list for all CT connections to perform       LDA, and generate a similar spreadsheet to use here.
    -	```inputexpanded_TC9.xls```
    -	```TC_CCnoconf_iter.xls``` or ```TC_CCconf_iter.xls```
    -	```inputexpanded_CC9.xls```
    -	```CC_noconf_iter.xls``` or ```CC_conf_iter.xls```
  * Output files:
    -```inputexpanded_TC9CT2.xls```, which include information on TC+CT connections such as sources, targets, and FF/FB mapping.
    -```TCCT_CCnoconf_iter.xls``` or ```TCCT_CCconf_iter.xls```, which include hierarchy scores of all cortical and thalamic regions before & after iterations
    -```ghs_TCCT.xls``` which contains global hierarchy scores of TC + CT connectivity data

### ```IterativeMethod.py```
  * Iterate hierarchy scores of all cortical & thalamic regions to refine the hierarchy
  * Import ```iterativeTCCT```

To evaluate how hierarchical the mouse brain is, compare the global hierarchy score of the original connectivity data to the global hierarchy scores of shuffled connectivity data. To do this, use ```run_CC_shuffled.py```, ```run_TC_shuffled.py```, and ```run_TCCT_shuffled.py```. 


