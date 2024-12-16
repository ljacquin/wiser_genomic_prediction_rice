[<img src="img/rice.jpg" width="600"/>]()

# genomic prediction for WISER, LS-means and BLUP phenotypes associated to rice traits

### 🎯 Objective

This repository contains R scripts designed for reproducible data analysis and results, aligned with the FAIR principles. The scripts perform data reformatting and phenotypic estimation using WISER, LS-means, and BLUP. The BLUP specifically integrate principal coordinates of genotypes, derived from genomic data, as fixed effects to account for population structure.

### 💻 Instructions

Download the ```wiser_genomic_prediction_rice``` repository in the current user's directory on a computing cluster or personal computer using one of the following commands :

  *  ```git clone git@github.com:ljacquin/wiser_genomic_prediction_rice.git``` <p> </p>
    or
  * ```git clone https://github.com/ljacquin/wiser_genomic_prediction_rice.git``` 
  <p> </p>
  
  ⚠️ Make sure``` git``` is installed beforehand; if not, install it with ```sudo apt install git```.
  <p> </p>

* Given that ```R ≥ 4.1.2``` is already installed, within the ```wiser_genomic_prediction_rice``` folder use the following command to install and test ```wiser_genomic_prediction_rice``` required ```R``` libraries : 

  * ```R -q --vanilla < src/requirements.R```
  * ```R -q --vanilla < src/test_requirements.R```
  <p> </p>
  
* The ```R``` script ```src/rice_data_treatment_and_analysis/rice_data_reformatting_and_blups_lsmeans_computation.R``` performs data reformatting and phenotypic estimation using LS-means and BLUP. The BLUP incorporate principal coordinates of genotypes, derived from genomic data, as fixed effects to account for population structure.

* The ```R``` script ```src/rice_genomic_prediction_and_analysis/rice_wiser_genomic_prediction_trait.R``` performs, for each trait, the genomic prediction tasks and analyses for the phenotypes estimated using WISER, LS-means, and BLUP. Note that this script also computes WISER's phenotypic estimates prior to the genomic prediction tasks.

* For genomic prediction tasks and analyses, execute the following commands to make scripts and programs executable:

  *  ```chmod u+rwx src/rice_genomic_prediction_and_analysis/*.sh```
  <p> </p>

* Finally, execute one of the following commands for executing the genomic prediction tasks and analyses :

  *  ```sbatch src/execute_rice_wiser_genomic_prediction_all_traits.sh```<p> </p>
    or
  * ```./src/execute_rice_wiser_genomic_prediction_all_traits.sh``` (i.e., interactive execution)
  <p> </p>

⚠️ The tasks and analyses performed by the ```R``` scripts in the ```wiser_genomic_prediction_rice``` repository can be run in either ```Unix/Linux``` or ```Windows``` environments, as long as ```R``` and the necessary libraries are installed. For local computations in ```RStudio```, ensure that the ```computation_mode``` variable is set to "local" in the ```R``` scripts located in ```src/```.

