# Group_Lasso_and_Multitask

Cox model with group lasso and multitask coupling term penalty.

To generate data: first run gettingTCGAsamplesWithExpressionAndClinicalData.R (checks all available TCGA data sets with RNA Seq V2) and then generate_data.r (downloads csv files for all available cancer types and creates 4 files for each: cancer_type_x.csv (gene expression data), cancer_type_os_months.csv (survival dates), cancer_type_os_status.csv (censoring/deceased status data), and cancer_type_groups.csv (pathway index for each of the genes)). The pathways used are in the file pathways_no_dupl.csv.

To run the survival model: move all files created in the previous step to the surv_data folder. Run Cox_GD_group_multitask.py for a desired combination of cancers. Output: beta and c-values for both single and multi-tasking runs (4 files).

To analyze: Run pathways_single.py (produces bar plot over most common pathways), ttest_multi.py (produces heatmaps over all cancer type combinations), ttest_single.py (compares grouped lasso versus naive lasso)

To test: Generate synthetic data mimicking COADREAD and STAD by running create_X.py and test them by running ttest_toy.py (need to rename toy_coadread and toy_stad in the file names to T1 and T2 first). The "real" betas used in the process are in dataforG.csv.
