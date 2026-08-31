# PB-NGS-in-cytopenia-investigation
This project is a retrospective audit of the use of PB NGS to assess capability as a screening tool in the investigation of cytopenias. The analysis obtained from these scripts characterises the cohort, including statistically significant changes between groups, assesses the clinically significant variants detected in the cohort, and determines the predictive power of PB NGS results. 

# Repository structure
R scripts are provided for: univariate analysis of cohort characterisation, oncoplot creation for NGS data assessment, geom_smooth plot of overall variant detection trend and predictive value assessments as per the labelled scripts. 

# Requirements
Each script installs/loads its own required packages at the top. R version and full package details are documented in the accompanying dissertation.

# How to run
1. Place the cleaned data file an appropriate folder and set as working directory— see Data availability note below.
2. Run scripts as needed for each analysis (they are independent rather than sequential):
3. Outputs (tables, oncoplot) are saved locally when scripts are run and are presented as results in the accompanying dissertation; they are not stored in this repository.

# Data availability
Patient-level data is **not included** in this repository for confidentiality reasons. The data comprises identifiable genomic information for 163 patients and cannot be shared publicly. Scripts are provided for methodological transparency and to allow reproduction of the analysis pipeline with appropriately governed access to the underlying data.

# Author
Jennifer Back
