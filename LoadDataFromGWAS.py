
import DataLoadingFunctions as dl
import logging
logging.basicConfig(filename='LoadDataFromGWAS.log', 
                    encoding='utf-8', 
                    level=logging.DEBUG,
                    format='%(asctime)s %(levelname)-8s %(message)s \n',
                    datefmt='%Y-%m-%d %H:%M:%S')

'''Python script to load data on pathway identification and GWAS 
data on the relationship to BMI and cancer outcome for the same SNPs.'''

# SECTION 1 - INPUT DATA
path='../Python/ICEP/'

# Inputs for results from: 
# "Harnessing tissue-specific genetic variation to dissect putative causal pathways between 
# body mass index and cardiometabolic phenotypes"
# Table saved from the supplementary material of the paper.
# This data is taken from: https://www.nature.com/articles/s41416-022-02060-6

col_labs={'rsid':'SNP','a1':'ea','a2':'nea','bx':'beta','bxse':'se','by':'beta','byse':'se',
          'Pathway1_suffix':'adipose','Pathway2_suffix':'brain','outcome':'trait','outcome_id':'id.outcome',
          'expo_lab':'BMI','coloc_lab':'PPA4'}
alltables               ='HarnessingExcelTables'
insheettyp              ='.xlsx'
# Sheets within the spreadsheet containing pathway colocization scores
pathway1_tab_list=[dl.sheet_table('Table 4',RS=4),dl.sheet_table('Table 7',RS=5,UC='A:G')]
pathway2_tab_list=[dl.sheet_table('Table 5',RS=4),dl.sheet_table('Table 7',RS=5,UC='K:Q')]
# Sheets within the spreadsheet containing the outcome_ids
outcome_id_tab_list=[dl.sheet_table('Table 13',4,24),dl.sheet_table('Table 14', 5,37),
                     dl.sheet_table('Table 18',3,41),dl.sheet_table('Table 18',43,69)]
expo_tab_list      =[dl.sheet_table('Table 16',3)]

# Set the input parameters based on the setup of the spreadsheet in the supplementary data
pathway_output_tab_list=[dl.sheet_table('Table 22',5,18),dl.sheet_table('Table 22',20,29)]
method_str='2SMR_Steiger_ivw'
col_labels={'outcome_lab':'outcome','bx_lab':'bx','bxse_lab':'se','nsnp_lab':'nsnp','path_lab':'exposure'}
logging.info('Manual inputs made for CAD paper')
logging.info('-----------------------')
# GET THE SNP LEVEL DATA FOR CLUSTERING
fullres_cad_df=dl.Load_dfs_for_clustering(path+alltables,col_labs,pathway1_tab_list,pathway2_tab_list,
                           outcome_id_tab_list,expo_tab_list)
# GET THE PATHWAY MR RESULTS FOR HYPOTHESIS TESTING
hypo_cad_df=dl.load_dfs_for_hypothesis_testing(path+alltables,pathway_output_tab_list,col_labels,method_str)

#--------------------------------------------------------------------------------------------
# Altnerative Inputs for results from: 
# "Disentangling the aetiological pathways between body mass index and site-specific cancer 
# risk using tissue-partitioned Mendelian randomisation"
# Supplementary data available from: https://www.nature.com/articles/s41416-022-02060-6 

col_labs={'rsid':'SNP','a1':'effect allele','a2':'other allele','bx':'beta','bxse':'se','by':'beta','byse':'se',
          'Pathway1_suffix':'adipose','Pathway2_suffix':'brain','outcome':'trait','outcome_id':'id.outcome',
          'coloc_lab':'PPA4','expo_lab':'BMI'}
alltables               ='41416_2022_2060_MOESM3_ESM'
insheettyp              ='.xlsx'
# Sheets within the spreadsheet containing pathway colocization scores
pathway1_tab_list=[dl.sheet_table('Table 1',5),dl.sheet_table('Table 2',5)]
pathway2_tab_list=[dl.sheet_table('Table 1',5),dl.sheet_table('Table 2',5)]
# Sheets within the spreadsheet containing the outcome_ids
outcome_id_tab_list=[dl.sheet_table('Table 4',4,30),dl.sheet_table('Table 5',19),
                     dl.sheet_table('Table 5',36,86),dl.sheet_table('Table 7',3,12),
                     dl.sheet_table('Table 7',14,44),dl.sheet_table('Table 8',3,19),
                     dl.sheet_table('Table 8',24,82)]
# Sheets within the spreadshet containing the exposure snp analysis
expo_tab_list=[dl.sheet_table('Table 1',5),dl.sheet_table('Table 2',5)]

# Set the input parameters based on the setup of the spreadsheet in the supplementary data
pathway_output_tab_list=[dl.sheet_table('Table 4',4,29),dl.sheet_table('Table 5',35,85),
                         dl.sheet_table('Table 8',24,82)]
col_labels={'outcome_lab':'outcome','bx_lab':'bx','bxse_lab':'se','nsnp_lab':'nsnp','path_lab':'exposure'}

logging.info('Manual inputs made for Cancer paper')
logging.info('-----------------------')
fullres_cancer_df=dl.Load_dfs_for_clustering(path+alltables,col_labs,pathway1_tab_list,pathway2_tab_list,
                           outcome_id_tab_list,expo_tab_list)
# GET THE PATHWAY MR RESULTS FOR HYPOTHESIS TESTING
hypo_cancer_df=dl.load_dfs_for_hypothesis_testing(path+alltables,pathway_output_tab_list,col_labels,method_str)


