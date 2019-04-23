## stunning-fiesta ##
# TCGA Assemble RJB #

Using adopted scripts from Jack's 'TCGA_Assemble' repo.

# Functions Script
In the Functions folder is a script with all key functions written by Jack for processing.

# LUAD_Read_V1.0
The 'luadread_tcga_pancancer.R' has all the pre-processing scripts for TCGA PanCancer Atlas for each of the following data:
1. Clinical
2. CNA (this data file is too big to be opened using Excel)
3. RNA-Seq - both median RPKM values and z-scores collated into one file (also too big to be opened in Excel)
4. Mutations
5. Fusions

# Analyse_Data_V1.0/Mutations
The 'luad_pancancer_mut_analysis.R' script has the following analyses which were done in April 2019:
1. Rough types and number of total mutations/patient
2. Detailed type of mutation
3. Protein change positions
4. Amino acid changes

The 'luad_pancancer_mut_subsetting.R' script was used in April 2019:
Used to subset predominatly KRAS single, STK11 single and KRAS/STK11 double mutation patients.

