import ieugwaspy as igd
import pandas as pd
from icecream import ic
import logging
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
        ''' Loads a dataframe from the excel file at filename, sheet s.sheet_name.
        Loaded using :py:func:`Pandas.read_excel()`. 
        The dataframe will ignore the first s.row_skip rows and go until row s.nrows.

        The only columns loaded will be s.cols.

        :returns: df
        '''
        # To allow for row skip defined and not nrows check the case. NONE for rows and cols do not cause issues
        if s.nrows:
            df=pd.read_excel(fname,sheet_name=s.sheet_name,skiprows=s.row_skip,nrows=s.nrows-s.row_skip,usecols=s.cols)
        else:
            df=pd.read_excel(fname,sheet_name=s.sheet_name,skiprows=s.row_skip,usecols=s.cols)
        return df
    def load_extra_tab(s,fname,i):
        ''' Loads a dataframe from the excel file at filename, sheet s.sheet_name.
        Loaded using :py:func:`Pandas.read_excel()`. 
        The dataframe will ignore the first s.row_skip rows and go until row s.nrows.

        The only columns loaded will be s.cols.

        :returns: df
        '''
        df=s.load(fname)
        for col_lab in df.keys():
            if '.%d'%i in col_lab:
                df.rename(columns={col_lab:col_lab.replace('.%d'%i,'')},inplace=True)
        return df
def load_tab_list(tab_list,fname):
    '''' 
    :param tab_list: list of items in :py:class:`sheet_table` format.
    :param fname:    The datafile name including path and type.

    -Iterate through the items in tab_list. 
    -Load the dataframe for each then merge. 
    -Duplicate rows will be deleted.
    
    :rtype: :py:mod:`Pandas.DataFrame`.
    :returns: Dataframe containing all the information of the tables indicates by the :py:class:`sheet_table` items in the list. '''
    df_tot=pd.DataFrame() 
    sn0=''
    i=0
    for tl in tab_list:
        sn=tl.sheet_name
        if sn==sn0:
            # Two table with the same structure on the same sheet will have a suffix in the column names. 
            i+=1
            df_single=tl.load_extra_tab(fname,i)
        else:
            df_single=tl.load(fname)
            i=0
        df_tot=pd.concat([df_single,df_tot],ignore_index=True)
    df_tot.drop_duplicates(keep='first',ignore_index=True,inplace=True)
    return df_tot

def Load_dfs_for_clustering(datafilename,ftype,col_labs,pathway1_tab_list,pathway2_tab_list,outcome_id_tab_list,expo_tab_list):
    '''
    :param datafilename: path and filename of excel spreedsheet containing data to load.
    :param col_loabs: dictionary of strings matching the column labels to desired information.
    :param pathway1_tab_list: list of the information for loading each table related to pathway1 in :py:class:`sheet_table` type.
    :param pathway2_tab_list: list of the information for loading each table related to pathway2 in :py:class:`sheet_table` type
    :param outcome_id_tab_list: list of the information for loading each table related to outcome ids in :py:class:`sheet_table` type.
    :param expo_tab_list: list of the information for loading each table related to exposures in :py:class:`sheet_table` type
    
    columns in output dataframe:=['rsid','chr.pos','exposure','outcome','a1','a2','bx','bxse','by','byse',
        'pathway1_PPA4','pathway2_PPA4']
        '''

    # Set the labels based on the inputs
    a1_lab      =col_labs['a1']
    a2_lab      =col_labs['a2']
    by_lab      =col_labs['by']
    byse_l      =col_labs['byse']
    rsid_lab    =col_labs['rsid']
    p1_lab      =col_labs['Pathway1_suffix']
    p2_lab      =col_labs['Pathway2_suffix']
    outcome_lab =col_labs['outcome_id']
    bxlab       =col_labs['bx']
    bxselab     =col_labs['bxse']
    coloc_lab   =col_labs['coloc_lab']
    expo_lab    =col_labs['expo_lab']

    # Load Pathway 1 and Pathway 2 colocalization scores.
    pathway1_df     =load_tab_list(pathway1_tab_list,datafilename+ftype)[[rsid_lab,coloc_lab]]
    pathway2_df     =load_tab_list(pathway2_tab_list,datafilename+ftype)[[rsid_lab,coloc_lab]]
    if coloc_lab+'_'+p1_lab in pathway1_df.keys() and coloc_lab+'_'+p2_lab in pathway2_df.keys():
        pathway_df      =pd.merge(pathway1_df,pathway2_df,how='inner',on=rsid_lab)
    else:
        pathway_df      =pd.merge(pathway1_df,pathway2_df,how='inner',on=rsid_lab,suffixes=('_'+p1_lab,'_'+p2_lab))
    pathway_df      =pathway_df[[rsid_lab,coloc_lab+'_'+p1_lab,coloc_lab+'_'+p2_lab]]
    pathway_df.dropna(subset=[rsid_lab,coloc_lab+'_'+p1_lab,coloc_lab+'_'+p2_lab],inplace=True)
    pathway_df.rename(columns={rsid_lab:'rsid'},inplace=True)

    logging.info('Pathway snp level data loaded')
    logging.info(ic(pathway_df.head()))
    # Retrieve the outcome_ids for getting the outcome results from open GWAS.
    outcome_id_df=load_tab_list(outcome_id_tab_list,datafilename+ftype)[outcome_lab]
    outcome_id_df.dropna(inplace=True)
    outcome_ids=outcome_id_df.values.tolist()
    logging.info('Outcome ids loaded')
    logging.info(ic(outcome_ids))

    # Load the SNP level data for the exposure
    expo_df=load_tab_list(expo_tab_list,datafilename+ftype)[[rsid_lab,bxlab,bxselab]]
    expo_df.rename(columns={rsid_lab:'rsid',bxlab:'bx',bxselab:'bxse'},inplace=True)
    expo_df.dropna(subset=['rsid','bx','bxse'],inplace=True)
    logging.info('Exposure snp level data loaded')
    logging.info(ic(expo_df.head()))

    # Load the snp level data for the outcomes from open gwas
    logging.info('-----------------------')
    logging.info('Accessing Open GWAS')
    assocs=pd.DataFrame.from_dict(igd.associations(expo_df['rsid'],outcome_ids,align_alleles=1, 
                                                 palindromes=1))
    logging.info('Outcome snp level data loaded from open GWAS')
    logging.info('-----------------------')
    logging.info(ic(assocs.head()))
    
    # Create the output dataframe for the snp level data.
    # Required form: 
    # rsid, chr.position, exposure, outcome, a1, a2, bx, xse, by, byse, pathway1_ppa4, pathway2_ppa4 
    labels=['rsid','chr.pos','exposure','outcome','a1','a2','bx','bxse','by','byse',
        p1_lab+'_PPA4',p2_lab+'_PPA4']
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

    logging.info('Data loaded for a1,a2,rsid,by,byse,chr.position')
    logging.info('Load the association between SNPs and exposure')
    # Match the outcome and exposure pairs to fill in the bx, bxse and PPA4 scores.
    for rid in expo_df['rsid']:
        out_df.loc[out_df.rsid==rid,'bx']       =expo_df.loc[expo_df.rsid==rid,    'bx'].item()
        out_df.loc[out_df.rsid==rid,'bxse']     =expo_df.loc[expo_df.rsid==rid,    'bxse'].item()
    for rid in pathway_df['rsid']:
        try: 
            out_df.loc[out_df.rsid==rid,p1_lab+'_PPA4']=pathway_df.loc[pathway_df.rsid==rid,coloc_lab+'_'+p1_lab].item()
        except: pass
        try: 
            out_df.loc[out_df.rsid==rid,p2_lab+'_PPA4']  =pathway_df.loc[pathway_df.rsid==rid,coloc_lab+'_'+p2_lab].item()
        except: pass
    logging.info('Full csv ready for clustering')
    logging.info(ic(out_df.head()))
    logging.info('---------------------------------------------------------------')
    out_df.to_csv(datafilename+'FullResults.csv',index=False)
    logging.info('File saved at '+datafilename+'FullResults.csv')
    return out_df

def load_dfs_for_hypothesis_testing(datafilename,ftype, pathway_tab_list,col_labels,method=None):
    '''Function to generalise creating a csv for hypothesis testing from an excel spreadsheet with tables. 
    The inputs need to be created by investigating the spreadsheet and matching the sheet names and row numbers. 

    If the method varies then retrieve the methods from the spreadsheet. Otherwise set the method to a input 
    string for the method used.

    :param datafilename: path and filename of excel spreedsheet containing data to load.
    :param pathway1_tab_list: list of the information for loading each table related to pathway1 in :py:class:`sheet_table` type.
    :param method: string for method if all are the same.

    The pathway dataframe is a exposure association table with pathway labelled SNPs.

    columns in output dataframe:=['Label','Outcome','bx','OR','se','nSNPs','Method']
    
    The label is a string indicating the dominent pathway and the method of computation. 
    
    :rtype: :py:mod:`Pandas.Dataframe`.'''
    outcome_lab =col_labels['outcome_lab']
    bx_lab      =col_labels['bx_lab']
    bxse_lab    =col_labels['bxse_lab']
    nsnp_lab    =col_labels['nsnp_lab']
    path_lab    =col_labels['path_lab']

    # Load the results for the pathway outcome level
    pathway_df=load_tab_list(pathway_tab_list,datafilename+ftype)
    # Initialise the desired format for the output dataframe
    hypo_col_labs=['Label','Outcome','bx','OR','se','nSNPs','Method']
    hypo_df=pd.DataFrame(columns=hypo_col_labs)
    if not method:
        hypo_df.method=[w.replace(" ","_") for w in pathway_df.method_lab]
    else:
        hypo_df.Method  =method
    hypo_df.Outcome =[w.replace(" ","_") for w in pathway_df[outcome_lab]]
    hypo_df.Label   =[c+'_'+str(p) for c,p in zip(pathway_df[path_lab],hypo_df.Method)]
    hypo_df.bx      =pathway_df[bx_lab]
    hypo_df.se      =pathway_df[bxse_lab]
    hypo_df.nSNPs   =pathway_df[nsnp_lab]
    logging.info('Created hypothesis test csv.')
    logging.info(ic(hypo_df))
    hypo_df.to_csv(datafilename+'ForHypoTest.csv',index=False)
    logging.info('file at '+datafilename+'ForHypoTest.csv')
    return hypo_df


