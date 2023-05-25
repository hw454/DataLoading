import ieugwaspy as igd
import pandas as pd
from icecream import ic
import logging
logging.basicConfig(filename='LoadDataFromGWAS.log', 
                    encoding='utf-8', 
                    level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s \n',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

'''Python script to load data on pathway identification and GWAS 
data on the relationship to BMI and cancer outcome for the same SNPs.'''

# SECTION 1 - GET THE SNP LEVEL DATA. 
# REQUIRED FORM: SNP, CHR.POSITION, A1, A2, BX, BXSE, BY, BYSE, ADIPOSE_PPA4, BRAIN_PPA4

# The colocalization dataset for the adipose and brain pathway scores
# # is saved outside of the GWAS. 
# This data for pathway identification is taken from: https://www.nature.com/articles/s41416-022-02060-6
# Load this with: rsid, chr.pos, exposure, outcome, a1, a2, bx, bxse, Adipose_P, Brain_P
path='../Python/ICEP/'
labels=['rsid','chr.pos','exposure','outcome','a1','a2','bx','bxse','by','byse',
        'Adipose_P','Brain_P']


# Inputs for results from: 
# "Harnessing tissue-specific genetic variation to dissect putative causal pathways between 
# body mass index and cardiometabolic phenotypes"
# Table saved from the supplementary material of the paper.
# outcome_ids=['ebi-a-GCST006464','ieu-a-1120','ieu-a-1126','ieu-b-85','QyUAvp','kL0KZ3',
# 'KBEFfw','ebi-a-GCST006464','ieu-a-1120','ieu-a-1126']
col_labs={'rsid':'SNP','a1':'ea','a2':'nea','bx':'beta','bxse':'se','by':'beta','byse':'se',
          'Adipose_P':'PPA4_adipose','Brain_P':'PPA4_brain','outcome':'trait','outcome_id':'id.outcome',
          'expo_lab':'BMI','coloc_lab':'PPA4'}
alltables               ='HarnessingExcelTables'
insheettyp='.xlsx'
# Sheets within the spreadsheet containing pathway colocization scores
pathway1_tab_list=[sheet_table('Table 4',RS=4),sheet_table('Table 7',RS=5,UC='A:G')]
pathway2_tab_list=[sheet_table('Table 5',RS=4),sheet_table('Table 7',RS=5,UC='K:Q')]
# Sheets within the spreadsheet containing the outcome_ids
outcome_id_tab_list=[sheet_table('Table 13',4,24),sheet_table('Table 14', 5,37),sheet_table('Table 18',3,41),sheet_table('Table 18',43,69)]
expo_tab_list      =[sheet_table('Table 16',3)]
class sheet_table:
    def __init__(s,SN,RS=None,NT=None,UC=None):
        '''
        :param SN: string indicating sheetname
        :param RS: Integer indicating number of rows to skip
        :param NT: Total number of rows being searched for
        :param UC: Collumns to be used'''
        s.sheet_name=SN
        s.row_skip  =RS
        s.nrows     =NT
        s.cols      =UC
    def load(s,fname):
        df=pd.read_excel(fname,sheet_name=s.SN,skiprows=s.row_skip,nrows=s.nrows,usecols=s.cols)
        return df
def load_tab_list(tab_list,fname):
    '''' Iterate through the items in tab_list. Load the dataframe for each then merge. Duplicate rows will be deleted.'''
    df_tot=pd.DataFrame() 
    for tl in tab_list:
        df_single=tl.load(fname)
        df_tot=pd.concat([df_single,df_tot],ignore_index=True)
    df_tot.drop_duplicates(keep='first',ignore_index=True,inplace=True)
    return df_tot

def Load_dfs(col_labs,pathway1_tab_list,pathway2_tab_list,outcome_id_tab_list,expo_tab_list):
    '''
    :param col_loabs: dictionary of strings matching the column labels to desired information.
    :param pathway1_tab_list: list of the information for loading each table related to pathway1 in :py:class:`sheet_table` type.
    :param pathway2_tab_list: list of the information for loading each table related to pathway2 in :py:class:`sheet_table` type
    :param outcome_id_tab_list: list of the information for loading each table related to outcome ids in :py:class:`sheet_table` type.
    :param expo_tab_list: list of the information for loading each table related to exposures in :py:class:`sheet_table` type'''
# Set the labels based on the inputs
a1_lab      =col_labs['a1']
a2_lab      =col_labs['a2']
by_lab      =col_labs['by']
byse_l      =col_labs['byse']
rsid_lab    =col_labs['rsid']
adi_lab     =col_labs['Adipose_P']
bri_lab     =col_labs['Brain_P']
outcome_lab =col_labs['outcome_id']
bxlab       =col_labs['bx']
bxselab     =col_labs['bxse']


logging.info('Manual inputs made')
logging.info('-----------------------')
#Load Adipose and Brain colocalization scores.
adi_pathway_df  =pd.read_excel(path+alltables+insheettyp,sheet_name=Adipose_Sheet,skiprows=nskippath)[[rsid_lab,coloc_lab]]
brain_pathway_df=pd.read_excel(path+alltables+insheettyp,sheet_name=Brain_Sheet,skiprows=nskippath)[[rsid_lab,coloc_lab]]
pathway_df      =pd.merge(adi_pathway_df,brain_pathway_df,how='inner',on=rsid_lab,suffixes=('_brain','_adipose'))
pathway_df      =pathway_df[['SNP',adi_lab,bri_lab]]
pathway_df.rename(columns={rsid_lab:'rsid'},inplace=True)

logging.info('Data loaded from excel spreadshet '+path+alltables+insheettyp)
# Load the adipose and brain borderline scores. 
adi_bord_path_df=pd.read_excel(path+alltables+insheettyp,sheet_name=pathway_borderline_sheet,
                              skiprows=nskip_bor,usecols='A:G')[[rsid_lab,coloc_lab]]
bri_bord_path_df=pd.read_excel(path+alltables+insheettyp,sheet_name=pathway_borderline_sheet,
                              skiprows=nskip_bor,usecols='K:Q')[[rsid_lab+'.1',coloc_lab+'.1']]
bri_bord_path_df.rename(columns={rsid_lab+'.1':rsid_lab,coloc_lab+'.1':coloc_lab},inplace=True)
# Combine the borderline scores
border_path_df  =pd.merge(adi_bord_path_df,bri_bord_path_df,how='inner',on=rsid_lab,suffixes=('_brain','_adipose'))
border_path_df=border_path_df[[rsid_lab,adi_lab,bri_lab]]
border_path_df.rename(columns={rsid_lab:'rsid'},inplace=True)
# Make one full dataframe for colocization scores on both pathways
pathway_df=pd.concat([pathway_df,border_path_df])

logging.info('Pathway snp level data loaded')

# Retrieve the outcome_ids for getting the outcome results from open GWAS.
outcome_id_df=pd.read_excel(path+alltables+insheettyp,sheet_name=Outcome_id_Sheet,skiprows=nskip_outcome,nrows=noutcomes-nskip_outcome)[outcome_lab]
outcome_id_mo=pd.read_excel(path+alltables+insheettyp,sheet_name=Outcome_id_Sheet2,skiprows=nskip_outcome2)[outcome_lab]
outcome_id_th=pd.read_excel(path+alltables+insheettyp,sheet_name=Outcome_id_Sheet3,skiprows=nskip_outcome3,nrows=noutcomes3-nskip_outcome3)[outcome_lab]
outcome_id_fo=pd.read_excel(path+alltables+insheettyp,sheet_name=Outcome_id_Sheet3,skiprows=nskip_outcome4)[outcome_lab]
# Combine all the outcome ids found in the spreadsheet and remove duplicates.
outcome_id_df=pd.concat([outcome_id_df,outcome_id_mo,outcome_id_th,outcome_id_fo])
outcome_id_df.drop_duplicates(keep='first',inplace=True)
outcome_ids=outcome_id_df.values.tolist()
logging.info('Outcome ids loaded')

# Load the SNP level data for the exposure
expo_df=pd.read_excel(path+alltables+insheettyp,sheet_name=expo_sheet,skiprows=nskip_expo)[[rsid_lab,bxlab,bxselab]]
expo_df.rename(columns={rsid_lab:'rsid',bxlab:'bx',bxselab:'bxse'},inplace=True)
logging.info('Exposure snp level data loaded')

ic(expo_df.head())
ic(outcome_ids)
# Load the snp level data for the outcomes from open gwas
logging.info('-----------------------')
logging.info('Accessing Open GWAS')
assocs=pd.DataFrame.from_dict(igd.associations(expo_df['rsid'],outcome_ids,align_alleles=1, 
                                                 palindromes=1))
logging.info('Outcome snp level data loaded from open GWAS')
logging.info('-----------------------')
ic(assocs.head())
# Create the output dataframe for the snp level data.
# Required form: 
# rsid, chr.position, exposure, outcome, a1, a2, bx, xse, by, byse, adipose_ppa4, brain_ppa4 
logging.info('Reformat the data into readible csv')
logging.info('CSV column labels: '+str(labels))
out_df=pd.DataFrame(columns=labels)
out_df['rsid']=assocs['rsid']
out_df['exposure']=expo_lab
out_df['chr.pos']=['chr'+str(c)+':'+str(p) for c,p in zip(assocs['chr'],assocs['position'])]
out_df['a1']    =assocs[a1_lab]
out_df['a2']    =assocs[a2_lab]
out_df['by']    =assocs[by_lab]
out_df['byse']  =assocs[byse_l]
out_df['outcome']=[w.replace(" ","_") for w in assocs['trait']]

# Matche the outcome and exposure pairs to fill in the bx, bxse and PPA4 scores.
for rid in expo_df['rsid']:
    out_df.loc[out_df.rsid==rid,'bx']       =expo_df.loc[expo_df.rsid==rid,    'bx'].item()
    out_df.loc[out_df.rsid==rid,'bxse']     =expo_df.loc[expo_df.rsid==rid,    'bxse'].item()
for rid in pathway_df['rsid']:
    try: 
        out_df.loc[out_df.rsid==rid,'Adipose_P']=pathway_df.loc[pathway_df.rsid==rid,'PPA4_adipose'].item()
    except: pass
    try: 
        out_df.loc[out_df.rsid==rid,'Brain_P']  =pathway_df.loc[pathway_df.rsid==rid,'PPA4_brain'].item()
    except: pass
print(out_df.head())
out_df.to_csv(path+alltables+'FullResults.csv',index=False)
print('file at',path+alltables+'FullResults.csv')

# GET THE PATHWAY MR RESULTS FOR HYPOTHESIS TESTING
# Set the input parameters based on the setup of the spreadsheet in the supplementary data
path_steig_sheet='Table 22'
nskip_path_steig_1=5
nskip_path_steig_2=20
ntot_path_steig_1=18
ntot_path_steig_2=29
method_steig='2SMR_Steiger_ivw'
# Load the results computed with steiger filtering
path_steig_df1=pd.read_excel(path+alltables+insheettyp,sheet_name=path_steig_sheet,
                         skiprows=nskip_path_steig_1,nrows=ntot_path_steig_1-nskip_path_steig_1)
path_steig_df2=pd.read_excel(path+alltables+insheettyp,sheet_name=path_steig_sheet,
                         skiprows=nskip_path_steig_2,nrows=ntot_path_steig_2-nskip_path_steig_2)
path_steig_df            =pd.concat([path_steig_df1,path_steig_df2],ignore_index=True)
path_steig_df['Method']  =method_steig
# Get the results from the results without Steiger filtering.
pathway_df=path_steig_df#pd.concat([path_steig_df,path_nofilt_df],ignore_index=True)
pathway_df.drop_duplicates(keep='first',inplace=True)
print(pathway_df.keys())
print(pathway_df)
# Load the results into the desired format for the output dataframe
hypo_col_labs=['Label','Outcome','bx','OR','se','nSNPs','Method']
hypo_df=pd.DataFrame(columns=hypo_col_labs)
hypo_df.Outcome =[w.replace(" ","_") for w in pathway_df.outcome]
hypo_df.Label   =[c+'_'+str(p) for c,p in zip(pathway_df.exposure,pathway_df.Method)]
hypo_df.OR      =pathway_df.OR
hypo_df.Method  =pathway_df.Method
hypo_df.bx      =pathway_df.b 
hypo_df.se      =pathway_df.se
hypo_df.nSNPs   =pathway_df.nsnp
print(hypo_df)
hypo_df.to_csv(path+alltables+'ForHypoTest.csv',index=False)
print('file at',path+alltables+'ForHypoTest.csv')

# Altnerative Inputs for results from: 
# "Disentangling the aetiological pathways between body mass index and site-specific cancer 
# risk using tissue-partitioned Mendelian randomisation"
# col_labs={a1:'ea',a2:'nea',bx:'bx',bxse:'bxse',by:'beta',byse:'se',Adipose_P:'PPA4_adipose',
# Brain_P:'PPA4_brain',outcome:'trait'}
col_labs={'rsid':'SNP','a1':'effect allele','a2':'other allele','bx':'beta','bxse':'se','by':'beta','byse':'se',
          'Adipose_P':'PPA4_adipose','Brain_P':'PPA4_brain','outcome':'trait','outcome_id':'id.outcome'}
alltables               ='41416_2022_2060_MOESM3_ESM'
insheettyp='.xlsx'
# Sheets within the spreadsheet containing pathway colocization scores
Pathway_Sheet1           ='Table 1'
Pathway_Sheet2           ='Table 2'
# Sheets within the spreadsheet containing the outcome_ids
Outcome_id_Sheet        ='Table 4'
Outcome_id_Sheet2       ='Table 5'
Outcome_id_Sheet3       ='Table 7'
Outcome_id_Sheet4       ='Table 8'
# Sheets within the spreadshet containing the exposure snp analysis
expo_sheet1              ='Table 1'
expo_sheet2              ='Table 2'
coloc_lab_adi           =coloc_lab+'_adipose'
coloc_lab_bri           =coloc_lab+'_brain'
nskippath=5
nskip_outcome=4
nskip_outcome2=5
nskip_outcome3=3
nskip_outcome4=24
nskip_expo=5
nskip_bor=5
noutcomes=29
noutcomes2=18
noutcomes3=14
expo_lab='BMI'

#Load Adipose and Brain colocalization scores.
pathway_df1  =pd.read_excel(path+alltables+insheettyp,sheet_name=Pathway_Sheet1,skiprows=nskippath)[[rsid_lab,coloc_lab_adi,coloc_lab_bri]]
pathway_df2  =pd.read_excel(path+alltables+insheettyp,sheet_name=Pathway_Sheet2,skiprows=nskippath)[[rsid_lab,coloc_lab_adi,coloc_lab_bri]]
pathway_df      =pd.concat([pathway_df1,pathway_df2],ignore_index=True)
pathway_df.drop_duplicates(keep='first',inplace=True)
pathway_df.rename(columns={rsid_lab:'rsid'},inplace=True)



# Retrieve the outcome_ids for getting the outcome results from open GWAS.
outcome_id_df=pd.read_excel(path+alltables+insheettyp,sheet_name=Outcome_id_Sheet,skiprows=nskip_outcome,nrows=noutcomes-nskip_outcome)[outcome_lab]
ic(outcome_id_df.head())
outcome_id_mo=pd.read_excel(path+alltables+insheettyp,sheet_name=Outcome_id_Sheet2,skiprows=nskip_outcome2,nrows=noutcomes-nskip_outcome2)[outcome_lab]
outcome_id_th=pd.read_excel(path+alltables+insheettyp,sheet_name=Outcome_id_Sheet3,skiprows=nskip_outcome3,nrows=noutcomes3-nskip_outcome3)[outcome_lab]
outcome_id_fo=pd.read_excel(path+alltables+insheettyp,sheet_name=Outcome_id_Sheet4,skiprows=nskip_outcome4)

# Combine all the outcome ids found in the spreadsheet and remove duplicates.
outcome_id_df=pd.concat([outcome_id_df,outcome_id_mo,outcome_id_th,outcome_id_fo])
outcome_id_df.drop_duplicates(keep='first',inplace=True)
outcome_ids=outcome_id_df.values.tolist()

# Load the SNP level data for the exposure
expo_df_1=pd.read_excel(path+alltables+insheettyp,sheet_name=expo_sheet1,skiprows=nskip_expo)[[rsid_lab,bxlab,bxselab]]
expo_df_2=pd.read_excel(path+alltables+insheettyp,sheet_name=expo_sheet2,skiprows=nskip_expo)[[rsid_lab,bxlab,bxselab]]
expo_df=pd.concat([expo_df_1,expo_df_2],ignore_index=True)
expo_df.drop_duplicates(keep='first',inplace=True)
expo_df.rename(columns={rsid_lab:'rsid',bxlab:'bx',bxselab:'bxse'},inplace=True)

# Load the snp level data for the outcomes from open gwas
assocs=pd.DataFrame.from_dict(igd.associations(expo_df['rsid'],outcome_ids,align_alleles=1, 
                                                 palindromes=1))
# Create the output dataframe for the snp level data.
# Required form: 
# rsid, chr.position, exposure, outcome, a1, a2, bx, xse, by, byse, adipose_ppa4, brain_ppa4 
out_df=pd.DataFrame(columns=labels)
out_df['rsid']=assocs['rsid']
out_df['exposure']=expo_lab
out_df['chr.pos']=['chr'+str(c)+':'+str(p) for c,p in zip(assocs['chr'],assocs['position'])]
out_df['a1']    =assocs[a1_lab]
out_df['a2']    =assocs[a2_lab]
out_df['by']    =assocs[by_lab]
out_df['byse']  =assocs[byse_l]
out_df['outcome']=[w.replace(" ","_") for w in assocs['trait']]

# Match the outcome and exposure pairs to fill in the bx, bxse and PPA4 scores.
for rid in expo_df['rsid']:
    out_df.loc[out_df.rsid==rid,'bx']       =expo_df.loc[expo_df.rsid==rid,    'bx'].item()
    out_df.loc[out_df.rsid==rid,'bxse']     =expo_df.loc[expo_df.rsid==rid,    'bxse'].item()
for rid in pathway_df['rsid']:
    try: 
        out_df.loc[out_df.rsid==rid,'Adipose_P']=pathway_df.loc[pathway_df.rsid==rid,'PPA4_adipose'].item()
    except: pass
    try: 
        out_df.loc[out_df.rsid==rid,'Brain_P']  =pathway_df.loc[pathway_df.rsid==rid,'PPA4_brain'].item()
    except: pass
print(out_df.head())
out_df.to_csv(path+alltables+'FullResults.csv',index=False)
print('file at',path+alltables+'FullResults.csv')

# GET THE PATHWAY MR RESULTS FOR HYPOTHESIS TESTING
# Set the input parameters based on the setup of the spreadsheet in the supplementary data
path_sheet1='Table 4'
path_sheet2='Table 5'
path_sheet3='Table 8'
nskip_path_1=4
nskip_path_2=35
nskip_path_3=24
ntot_path_1=29
ntot_path_2=85
ntot_path_3=82
# Load the results computed with steiger filtering
path_df1=pd.read_excel(path+alltables+insheettyp,sheet_name=path_sheet1,
                         skiprows=nskip_path_1,nrows=ntot_path_1-nskip_path_1)
path_df2=pd.read_excel(path+alltables+insheettyp,sheet_name=path_sheet2,
                         skiprows=nskip_path_2,nrows=ntot_path_2-nskip_path_2)
path_df3=pd.read_excel(path+alltables+insheettyp,sheet_name=path_sheet3,
                         skiprows=nskip_path_3,nrows=ntot_path_3-nskip_path_3)
pathway_df            =pd.concat([path_df1,path_df2,path_df3],ignore_index=True)
pathway_df.method=[w.replace(" ","_") for w in pathway_df.method]
# Get the results from the results without Steiger filtering.
pathway_df.drop_duplicates(keep='first',inplace=True)
print(pathway_df.keys())
print(pathway_df)
# Load the results into the desired format for the output dataframe
hypo_col_labs=['Label','Outcome','bx','OR','se','nSNPs','Method']
hypo_df=pd.DataFrame(columns=hypo_col_labs)
hypo_df.Outcome =[w.replace(" ","_") for w in pathway_df.outcome]
hypo_df.Label   =[c+'_'+str(p) for c,p in zip(pathway_df.exposure,pathway_df.method)]
hypo_df.OR      =pathway_df.OR
hypo_df.Method  =pathway_df.Method
hypo_df.bx      =pathway_df.b 
hypo_df.se      =pathway_df.se
hypo_df.nSNPs   =pathway_df.nsnp
print(hypo_df)
hypo_df.to_csv(path+alltables+'ForHypoTest.csv',index=False)
print('file at',path+alltables+'ForHypoTest.csv')

