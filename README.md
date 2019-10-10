# Mouse Brain Hierarchy

* Hierarchy searching algorithm based on layer-projection patterns of cortico-cortical, thalamo-cortical, and cortico-thalamic connections constructed from viral tracing experiments using Cre-driver lines
* Code used in Harris et al (2019) “Hierarchical organization of cortical and thalamic connectivity”
* If you have any questions or suggestions on the code, please contact Hannah Choi (hannahc@alleninstitute.org). 
* Written 8/9/2019. Latest update 10/10/2019 by Hannah Choi
************************************************************************************
This Python code builds the hierarchy of the mouse cortical and thalamic regions based on clustered layer-projection patterns by maximizing the self-consistency of the obtained hierarchy measured by global hierarchy score, followed by iterative corrections. The hierarchy is first constructed based on the cortico-cortical connections only, which can be updated by taking thalamo-cortical and cortico-thalamic connections in to account. 

## Level of Support
We are releasing this code to the public as a tool we expect others to use. Questions concerning bugs and related issues are welcomed. We expect to address them promptly, pull requests will vetted by our staff before inclusion.

## Cortico-cortical hierarchy

### ```run_CC.py```
  * The main hierarchy generating code based on cortico-cortical connectivity data.
  * Calls:
    - ```func_unsupervised_CC.py```
    - ```IterativeMethod.py```
  *	Input files:
    -	```CC_TC_CT_clusters.xlsx```
  *	Output files:
    -	```CC_conf_iter.xls``` & ```CC_noconf_iter.xls```which include hierarchy scores of all cortical regions before & after iterations with or without Cre-confidence. Either of these files is read in ```run_TC.py``` and ```func_unsupervised_TC.py```.
    -	```inputexpanded_CC9.xls```, which include information on CC connections such as sources, targets, and FF/FB mapping. This file is          read in ```run_TC.py``` for iteration with TC connections. 
    -	```ghs_CC.xls``` which contains global hierarchy scores of CC connectivity data

### ```func_unsupervised_CC.py```
  * Contains functions needed to map CC clusters to FF/FB by maximizing the global hierarchy score of the given cortico-cortical       connectivity data.

### ```IterativeMethod.py```
  * Iterate hierarchy scores of all cortical regions to refine the hierarchy
  *	Import ```iterativeCC```

## Adding thalamo-cortical hierarchy
The cortico-cortical hierarchy can be updated by including the thalamo-cortical connectivity data. Since the initial thalamic hierarchical positions are determined based on CC hierarchy scores of their target cortical regions, ```run_CC.py``` should be called before running ```run_TC.py```.

### ```run_TC.py```
  * The main hierarchy generating code based on thalamo-cortical connectivity data.
  * Calls:
    -	```func_unsupervised_TC.py```
    -	```IterativeMethod.py```
  *	Input files:
    -	```CC_TC_CT_clusters.xlsx```
    -	```CC_conf_iter.xls``` or ```CC_noconf_iter.xls```
    - ```inputexpanded_CC9.xls```
  * Output files:
    - ```TC_CCconf_iter.xls``` (used Cre-confidence for CC hierarchy) or ```TC_CCnoconf_iter.xls``` (no Cre-confidence for CC hierarchy), which include hierarchy scores of all cortical and thalamic regions before & after iterations. 
    - ```inputexpanded_TC9.xls```, which include information on TC connections such as sources, targets, and FF/FB mapping. This file is          read in ```run_TCCT.py``` for iteration with TC connections. 
    - ```ghs_TC.xls``` which contains global hierarchy scores of TC connectivity data

### ```func_unsupervised_TC.py``` 
  *	Contains functions needed for mapping TC clusters to FF/FB by maximizing the global hierarchy score of the given thalamo-cortical       connectivity data. As the initial thalamic hierarchical positions and the FF/FB mapping are based on the cortical hierarchy             constructed from ```run_CC.py```, the cortical hierarchy scores from CC connectivity should be read accordingly. (In ```run_CC.py```, use either ```df_cortex=pd.read_excel(r'./Output/CC_conf_iter.xls')``` or ```df_cortex=pd.read_excel(r'./Output/CC_noconf_iter.xls')```, depending on whether you use CC hierarchy computed with Cre-confidence or not. )  

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
    - ```CT_sourcelayer_FFB.xls``` has the FF/FB assignments for cortico-thalamic connections based on Linear Discriminant Analysis (LDA) results based on 1) L5/L6 source NPV, and 2) FF/FB-direction of each cortico-thalamic connection determined from CC & TC connectivity data. The FF/FB assignments of CT connections in this excel file was obtained based on the hierarchy constructed CC & TC connectivity data with Cre-confidence for CC hierarhcy. If not using Cre-confidence for CC hierarchy, use the accordingly updated FF/FB assignments for all CT connections to perform LDA, and generate a similar spreadsheet to use here.
    -	```inputexpanded_TC9.xls```
    -	```TC_CCnoconf_iter.xls``` or ```TC_CCconf_iter.xls```
    -	```inputexpanded_CC9.xls```
    -	```CC_noconf_iter.xls``` or ```CC_conf_iter.xls``` 
  * Output files:
    - ```inputexpanded_TC9CT2.xls```, which include information on TC+CT connections such as sources, targets, and FF/FB assignments.
    - ```TCCT_CCnoconf_iter.xls``` or ```TCCT_CCconf_iter.xls```, which include hierarchy scores of all cortical and thalamic regions before & after iterations
    - ```ghs_TCCT.xls``` which contains global hierarchy scores of CC+TC+CT connectivity data

### ```IterativeMethod.py```
  * Iterate hierarchy scores of all cortical & thalamic regions
  * Import ```iterativeTCCT```

## Comparison to shuffled connectivity data
To evaluate how hierarchical the mouse brain is, compare the global hierarchy scores of the original connectivity data to the global hierarchy scores of shuffled connectivity data. To do this, use ```run_CC_shuffled.py```, ```run_TC_shuffled.py``` (+```func_unsupervised_TC_shuffled.py```), and ```run_TCCT_shuffled.py```(+```func_unsupervised_TC_shuffled.py```).

## Inter-module and intra-module hierarchy
To find hierarchies of different cortical modules or within each module, run ```run_CC_module.py```, ```run_TC_module.py```, and ```run_TCCT_module.py``` to obtain CC, CC+TC, and CC+TC+CT-based hierarchies. For inter/intra-module hierarchy, instead of searching for the optimal mapping based on the limited numbers of intra-module or inter-module connections every time, use the pre-generated mapping of each of 9 clusters to either FF or FB direction constructed from all CC and TC connections, available in ```clustermapping.xlsx```. In the paper, we show intra-module hierarchy for visual-medial module (module = 'VisualMedial') as well as inter-module hierarchy (module ='inter_predefined').

## Figures and results
The summerized results of pre-generated hierarchy scores and global hierarchy scores are available in the ```\Results``` folder. You can use ```FigureGenerator.py``` to produce summary figures based on the pre-generated results.

## Terms of Use
https://alleninstitute.org/legal/terms-use/
