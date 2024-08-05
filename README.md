# Identifying molecular actionable exposure that are potentially causative in neurodegenerative diseases

A repository to reproduce the analyses explained in Constantino de la Vega's dissertation titled "Identifying molecular actionable exposure that are potentially causative in neurodegenerative diseases", which was submitted for a masters in philosophy to the University of Cambridge in August 2024. 

This repository has the following folders:
- 1_data_transformation: This folder contains guidance on how to perform the data transformation steps mentioned in Methods section 3.3.1.
- 2_quality_control&formatting: This folder contains guidance on how to perform the quality control and formatting steps mentioned in Methods section 3.3.2 and 3.3.3, respectively.
- 3_LD_matrix_&_SuSiE: This folder contains guidance on how to perform LD matrix computation and sum of single effects (SuSiE) fine-mapping for each gene we analysed, explained in 3.4. 
- 4_HyPrColoc: This folder contains the data and code necessary to perform our HyPrColoc analyses on CR1, BST1, DDR1, HIP1R, and PILRA, expalined in 3.5. We excluded the ADAMTS4 and STX6 HyPrColoc analyses because they require the PSP data which is not publicly available. 
- 5_MR: This folder contains the code and data necessary to perform all of our two-sample MR analyses explained in 3.6. Importantly, this folder also contains the code to create Figures 15, 17, 19, 20D, and 21E. 

Additionally, each folder contains its own README file going into further detail about how to carry out each part. It is important to note that since the PSP data is not publicly avialble, I will not include it in the code. 
