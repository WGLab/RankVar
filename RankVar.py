import argparse
import pandas as pd
import numpy as np
import os
import joblib
import torch

#define exonic function
def exonic_fun(fun):
    if (fun=="nonsynonymous SNV"):
        return 0
    else:
        return 1

def impute_missing(df):
    # Set conditions for splicing and loss-of-function (lof) mutations
    condition = (
        (df['Func.refGene'] == 'splicing') |
        (((df['Func.refGene'] == 'exonic') | (df['Func.refGene'] == 'exonic;splicing')) &
         ((df['ExonicFunc.refGene'] == 'stopgain') |
          (df['ExonicFunc.refGene'] == 'stoploss') |
          (df['ExonicFunc.refGene'] == 'startloss') |
          (df['ExonicFunc.refGene'] == 'frameshift deletion') |
          (df['ExonicFunc.refGene'] == 'frameshift insertion')))
    )
    
    # Apply conditions to fill columns from 13 to the last column
    columns_to_fill =[col for col in df.columns if col.endswith('score')]   # Select columns starting from index 13 to the end
    for col in columns_to_fill:
        df.loc[condition & (df[col] == '.'), col] = 0.8  
        df.loc[~condition & (df[col] == '.'), col] = 0.0
    
    return df

def add_feature(mydata):
    dirname = os.path.dirname(__file__)
    #add LOEUF, PLI, and Haploinsufficiency
    #add LOEUF & PLI score
    LOEUF_PLI=pd.read_csv(os.path.join(dirname, 'database/PLI_LOEUF.bed'),sep='\t')
    LOEUF_PLI=LOEUF_PLI[['geneName','_loeuf','_pli']]
    LOEUF_PLI.columns=['Gene.refGene','LOEUF','PLI']
    # Calculate the rank of each score (1 being the highest)
    LOEUF_PLI['LOEUF_rank'] = LOEUF_PLI['LOEUF'].rank(method='min', ascending=False)
    # Calculate the rankscore as rank divided by the total number of scores
    LOEUF_PLI['LOEUF_rankscore'] = LOEUF_PLI['LOEUF_rank'] / len(LOEUF_PLI) 
    LOEUF_rankscore=LOEUF_PLI[['Gene.refGene','PLI','LOEUF','LOEUF_rankscore']]
    LOEUF_rankscore = LOEUF_rankscore.drop_duplicates(subset=['Gene.refGene'])

    # add evolutionary constrainr (fracCdsCons, fracConsPr)
    constraint=pd.read_csv(os.path.join(dirname, 'database/constraint.csv'), sep=',')
    constraint.columns=['Gene.refGene','fracCdsCons','fracConsPr']
    constraint=constraint.drop_duplicates(subset=['Gene.refGene'])

    #add omim gene
    omim=pd.read_csv(os.path.join(dirname, 'database/genemap2_mim.txt'), sep='\t')
    omim.columns=['omim','Gene.refGene','Phenotypes']
    #filter with both gene and phenotype
    omim_gene=omim[(omim['Gene.refGene'].notna()) & (omim['Phenotypes'].notna())].copy()
    omim_gene['omim_score']=1
    omim_gene=omim_gene.drop_duplicates(subset=['Gene.refGene'])

    #merge data and LOEUF_rankscore
    data_LOEUF=pd.merge(mydata, LOEUF_rankscore, on='Gene.refGene',how='left')
    data_LOEUF['LOEUF_rankscore'] = data_LOEUF['LOEUF_rankscore'].fillna(0)
    data_LOEUF['PLI_score'] = data_LOEUF['PLI'].fillna(0)

    # add Haploinsufficiency score
    hap_score=pd.read_csv(os.path.join(dirname, 'database/haploinsufficiency.csv'), sep='\t')
    hap_score=hap_score[['name','prediction']]
    hap_score.columns=['Gene.refGene','Haploinsufficiency_score']
    #merge data and haploinsufficiency score
    data_Haplo=pd.merge(data_LOEUF, hap_score, on='Gene.refGene',how='left')
    #fill nan with '.' (used for the next step to replace)
    data_Haplo['Haploinsufficiency_score']=data_Haplo['Haploinsufficiency_score'].fillna(0)

    #merge data and constraint score
    data_constraint=pd.merge(data_Haplo, constraint, on='Gene.refGene',how='left')
    #fill na with 0
    data_constraint['fracCdsCons_score'] = data_constraint['fracCdsCons'].fillna(0)
    data_constraint['fracConsPr_score'] = data_constraint['fracConsPr'].fillna(0)

    # add disease constraint score (pHaplo, pTriplo)
    pHaplo=pd.read_csv(os.path.join(dirname,'database/pHaplo.csv'), sep=',')
    pHaplo=pHaplo[['Gene','pHaplo','pTriplo']]
    pHaplo.columns=['Gene.refGene','pHaplo','pTriplo']
    #merge data and disease constraint score
    data_pHaplo=pd.merge(data_constraint, pHaplo, on='Gene.refGene',how='left')
    #fill na with 0
    data_pHaplo['pHaplo_score'] = data_pHaplo['pHaplo'].fillna(0)
    data_pHaplo['pTriplo_score'] = data_pHaplo['pTriplo'].fillna(0)

    #add RVIS score
    RVIS=pd.read_csv(os.path.join(dirname, 'database/RVIS.csv'), sep=',')
    RVIS.columns=['Gene.refGene','RVIS','RVIS_rank']
    #merge data and RVIS score
    data_RVIS=pd.merge(data_pHaplo, RVIS, on='Gene.refGene',how='left')
    # Calculate the rankscore as rank divided by the total number of scores
    data_RVIS['RVIS_rankscore'] = 1 - ((data_RVIS['RVIS_rank'] - 1) / (data_RVIS['RVIS_rank'].max() - 1))
    #fill na with 0
    data_RVIS['RVIS_rankscore']=data_RVIS['RVIS_rankscore'].fillna(0)
    
    #add clinvarNum
    clinvarNum=pd.read_csv(os.path.join(dirname, 'database/clinvarNum.csv'), sep='\t')
    #merge data and clinvarNum #########
    data_clinvarNum=pd.merge(data_RVIS, clinvarNum, on='Gene.refGene',how='left')
    #file na with 0
    data_clinvarNum['clinvarNumB_LB_score']=data_clinvarNum['clinvarNumB_LB_score'].fillna(0.3)
    data_clinvarNum['clinvarNumP_LP_score']=data_clinvarNum['clinvarNumP_LP_score'].fillna(0)
    
    #add spliceai score
    spliceai=pd.read_csv(os.path.join(dirname, 'database/SpliceAI_score.txt'),sep='\t')
    # Merge spliceai with data
    data_spliceai = pd.merge(data_clinvarNum, spliceai, on=['Chr','Start'], how='left')
    data_spliceai['spliceAImax_score'] = data_spliceai['spliceAImax_score'].fillna(0)
    
    #add mis_ratio
    mis_ratio=pd.read_csv(os.path.join(dirname, 'database/mis_ratio_clinvar.txt'),sep='\t')
    # Merge mis_ratio with data
    data_mis=pd.merge(data_spliceai, mis_ratio, on='Gene.refGene',how='left')
    data_mis['mis_ratio_score'] = data_mis['mis_ratio_score'].fillna(0)
    
    #add inheritance pattern
    inheritance=pd.read_csv(os.path.join(dirname, 'database/omim_inheriate.csv'),sep=',')
    inheritance=inheritance[['Approved Gene Symbol','recessive','dominant']]
    inheritance.columns=['Gene.refGene','recessive','dominant']
    
    # Merge inheritance pattern with data
    data_inheritance=pd.merge(data_mis,inheritance,on='Gene.refGene',how='left')
    data_inheritance['Recessive_score']=data_inheritance['recessive'].fillna(0)
    data_inheritance['Dominant_score']=data_inheritance['dominant'].fillna(0)
    
    #merge data and omim gene
    data_omim=pd.merge(data_inheritance, omim_gene, on='Gene.refGene',how='left')
    #file na with 0
    data_omim['omim_score']=data_omim['omim_score'].fillna(0)

    #map def
    data_omim['exonic_fun']=data_omim['ExonicFunc.refGene'].map(exonic_fun)

    return data_omim

# Function to map multiple genes to their corresponding HPO terms
def map_genes_to_hpo(gene_str, gene_hpo_mapping):
    if pd.isna(gene_str):  # Check if the value is NaN
        return []  # Return an empty list for NaN values
    genes = gene_str.split('|')  # Split the gene names by '|'
    hpo_terms = [gene_hpo_mapping.get(gene, []) for gene in genes]  # Get HPO terms for each gene
    # Flatten the list and remove duplicates
    hpo_terms = {hpo for sublist in hpo_terms for hpo in sublist}
    return hpo_terms

def predict(df):
    dirname = os.path.dirname(__file__)
    rf_model_loaded = joblib.load(os.path.join(dirname, 'models/random_forest_model_with_weights.pkl'))
    X=torch.tensor(df.values.astype(np.float32), dtype=torch.float32)
    rf_predictions = rf_model_loaded.predict_proba(X)[:,1]
        
    return rf_predictions

def extract_info_columns(df):
    # Split Otherinfo12 and Otherinfo13 into lists
    df['info_keys'] = df['Otherinfo12'].str.split(':')
    df['info_values'] = df['Otherinfo13'].str.split(':')
    
    # Extract values for AD and GQ
    for col in ['AD', 'GQ']:
        df[col] = df.apply(
            lambda row: row['info_values'][row['info_keys'].index(col)]
            if col in row['info_keys'] else None,
            axis=1
        )
    
    # Extract AD_alt from AD
    df['AD_alt'] = df['AD'].apply(lambda x: x.split(',')[1] if x and ',' in x else None)

    # Convert columns to integers
    df['GQ'] = pd.to_numeric(df['GQ'], errors='coerce').fillna(0).astype(int)
    df['AD_alt'] = pd.to_numeric(df['AD_alt'], errors='coerce').fillna(0).astype(int)

    # Drop helper columns
    df = df.drop(columns=['info_keys', 'info_values'])
    
    return df

def process_file(annotate_df, gq, ad, gnomad):
    dirname = os.path.dirname(__file__)
    white_list_path=os.path.join(dirname, 'database/white_list.txt')
    
    white_list = pd.read_csv(white_list_path, sep=",")  # Adjust `sep` based on your file format
    white_list['variant_id'] = "chr" + white_list['Chr'] + "_" + white_list['Start'].astype(str) + "_" + white_list['Ref'] + "_" + white_list['Alt']
    # Load annotation file
    
    # Extract AD, GQ, and AD_alt from Otherinfo12 and Otherinfo13
    annotate_df = extract_info_columns(annotate_df)
    
    # Replace '.' in gnomad41_exome_AF_grpmax with 0 and convert to float
    annotate_df['gnomad41_exome_AF_grpmax'] = annotate_df['gnomad41_exome_AF_grpmax'].replace('.', 0).astype(float)
    # Replace '.' in gnomad41_genome_AF_grpmax with 0 and convert to float
    annotate_df['gnomad41_genome_AF_grpmax'] = annotate_df['gnomad41_genome_AF_grpmax'].replace('.', 0).astype(float)
    
    # Create unique variant identifier
    annotate_df['variant_id'] = annotate_df['Chr'].astype(str) + "_" + annotate_df['Start'].astype(str) + "_" + annotate_df['Ref'] + "_" + annotate_df['Alt']
    
    # Define the functional filter
    functional_terms = ["stopgain", "stoploss", "startloss","frameshift deletion", "frameshift substitution", "frameshift insertion", "nonsynonymous SNV"]
    functional_filter = annotate_df['ExonicFunc.refGene'].isin(functional_terms)
    splicing_filter = annotate_df['Func.refGene'] == "splicing"

    # Separate whitelist and other variants
    white_list_variants = annotate_df[
        annotate_df['variant_id'].isin(white_list['variant_id']) &
        (annotate_df['GQ'] > gq) & (annotate_df['AD_alt'] > ad)
    ]

    functional_variants = annotate_df[
        (~annotate_df['variant_id'].isin(white_list['variant_id'])) &
        (annotate_df['gnomad41_exome_AF_grpmax'] <= gnomad) & #exome
        functional_filter &
        (annotate_df['GQ'] > gq) & (annotate_df['AD_alt'] > ad) 
    ]
    
    splicing_variants = annotate_df[
        (~annotate_df['variant_id'].isin(white_list['variant_id'])) &
        (annotate_df['gnomad41_genome_AF_grpmax'] <= gnomad) & # genome
        splicing_filter &
        (annotate_df['GQ'] > gq) & (annotate_df['AD_alt'] > ad)
    ]

    # Combine results
    filtered_variants = pd.concat([white_list_variants, functional_variants, splicing_variants])
    
    # Remove duplicates based on 'Chr' and 'Start'
    filtered_variants = filtered_variants.drop_duplicates(subset=['Chr', 'Start'])
    
    return filtered_variants

def main(annovar, phen2gene, hpo_ids, gq, ad, gnomad, output):
    dirname = os.path.dirname(__file__)
    print(f"Annotated VCF: {annovar}")
    print(f"Phen2Gene Score: {phen2gene}")
    print(f"HPO Ids: {hpo_ids}")
    print(f"GQ: {gq}")
    print(f"AD: {ad}")
    print(f"GnomAD: {gnomad}")
    
    gt=pd.read_csv(annovar, sep='\t')
    gt=process_file(gt, gq, ad, gnomad)
    
    phen2gene_score=pd.read_csv(phen2gene, sep='\t')
    phen2gene_score=phen2gene_score[['Gene','Score','Rank']]

    #merge phen2gene score to annotated txt
    variants_gt=pd.merge(gt, phen2gene_score, left_on='Gene.refGene', right_on='Gene', how='left')

    #fill nan with 0
    variants_gt['phen2gene_score'] = variants_gt['Score'].fillna(0)
    variants_gt['phen2gene_score_rank'] = variants_gt['Rank'].fillna(20000)
    # Step 2: Normalize by the total number of genes (20,000)
    total_genes = 20000
    variants_gt['phen2gene_rankscore'] = -np.log10(variants_gt['phen2gene_score_rank']/total_genes)
    
    #select columns from annotated data
    info_columns=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 
                  'gnomad41_exome_AF_grpmax', 'gnomad41_genome_AF_grpmax', 'phen2gene_score']
    rankscore_columns = [col for col in variants_gt.columns if col.endswith('rankscore') or col.endswith('QTL_gene')]
    select_columns=info_columns+rankscore_columns
    variants_data=variants_gt[select_columns].copy()

    #add feature
    variants_final_omim=add_feature(variants_data)
    variants_final_omim=variants_final_omim.drop(['PLI', 'LOEUF', 'fracCdsCons', 'fracConsPr', 'pHaplo', 'pTriplo', 'RVIS', 'RVIS_rank', 'omim', 'Phenotypes', 'recessive', 'dominant'], axis=1)
    column_info_from_anno=variants_final_omim.columns.tolist()
    variants_final_omim=impute_missing(variants_final_omim)
    
    #gene mapped hpo terms
    gene_to_pheno=pd.read_csv(os.path.join(dirname,'database/genes_to_phenotype.txt'), sep='\t')
    # Drop duplicates based on 'gene_symbol' and 'hpo_id', keeping the first occurrence
    gene_to_pheno = gene_to_pheno.drop_duplicates(subset=['gene_symbol', 'hpo_id'], keep='first')
    gene_to_pheno_df=pd.DataFrame(gene_to_pheno.groupby('gene_symbol')['hpo_id'].apply(list)).reset_index()
    gene_to_pheno_df.columns=['Gene.refGene','hpo']

    # Aggregate HPO terms for each gene as lists
    #gene_to_pheno with information for gene and mapped hpo
    gene_hpo_mapping = gene_to_pheno.groupby('gene_symbol')['hpo_id'].apply(list).to_dict()

    
    # Iterate through each gene column in the second DataFrame and create corresponding HPO columns
    qtl_df = pd.DataFrame(index=variants_final_omim.index)
    for col in variants_final_omim.columns:
        if col.endswith('QTL_gene'):  # Check if it's a gene-related column
            hpo_col = col.replace('QTL_gene', 'QTL_hpo')  # New column name for HPO
            qtl_df[hpo_col]=pd.Series(variants_final_omim[col].apply(map_genes_to_hpo, gene_hpo_mapping=gene_hpo_mapping)) # Correctly apply the function
    variants_final_omim=pd.concat([variants_final_omim, qtl_df],axis=1)
    variants_final_omim=variants_final_omim.drop(columns=[col for col in variants_final_omim.columns if col.endswith('QTL_gene')])
    
    with open(hpo_ids, 'r') as f:
        hpo_ids=set(f.read().splitlines())

   # Create a temporary DataFrame to store new columns
    new_columns = pd.DataFrame(index=variants_final_omim.index)

    # Loop through the '_hpo' columns and calculate the corresponding '_score' column
    for col in variants_final_omim.columns:
        if col.endswith('_hpo'):  # Check if it's an '_hpo' column
            count_col = col.replace('_hpo', '_score')  # New column name for count
            # Calculate the number of coinciding HPO terms
            new_columns[count_col] = variants_final_omim.apply(
                lambda row: len(hpo_ids.intersection(row[col])), axis=1
            )
    # Concatenate the new columns DataFrame to the original DataFrame
    variants_final_omim = pd.concat([variants_final_omim, new_columns], axis=1)

    # Step 1: Identify columns ending with 'eQTL_score' and 'sQTL_score'
    eQTL_columns = [col for col in variants_final_omim.columns if col.endswith('eQTL_score')]
    sQTL_columns = [col for col in variants_final_omim.columns if col.endswith('sQTL_score')]
    eQTL_hpo_columns = [col for col in variants_final_omim.columns if col.endswith('eQTL_hpo')]
    sQTL_hpo_columns = [col for col in variants_final_omim.columns if col.endswith('sQTL_hpo')]

    # Step 2: Combine the values (sum across rows) for the identified columns
    variants_final_omim['eQTL_tissue_score'] = variants_final_omim[eQTL_columns].sum(axis=1)
    variants_final_omim['sQTL_tissue_score'] = variants_final_omim[sQTL_columns].sum(axis=1)
    variants_final_omim=variants_final_omim.drop(columns=eQTL_columns + sQTL_columns+eQTL_hpo_columns+sQTL_hpo_columns)
    
    variants_final_omim = variants_final_omim[~variants_final_omim['Gene.refGene'].str.startswith(('ZNF', 'MUC', 'HLA', 'FAM', 'LOC'), na=False)]
    
    columns=['MetaRNN_rankscore', 'PLI_score', 'phen2gene_rankscore', 'omim_score', 'LOEUF_rankscore', 'exonic_fun', 'SIFT_converted_rankscore', 'DANN_rankscore', 'Haploinsufficiency_score', 'fracCdsCons_score', 'AlphaMissense_rankscore', 'clinvarNumB_LB_score', 'spliceAImax_score', 'eQTL_tissue_score', 'sQTL_tissue_score', 'RVIS_rankscore', 'clinvarNumP_LP_score', 'FATHMM_converted_rankscore', 'Recessive_score', 'Dominant_score']
    
    variants_final_omim_select=variants_final_omim[columns]
    variants_final_omim_info=variants_final_omim[info_columns].copy()
    predicitions=predict(variants_final_omim_select)
    variants_final_omim_info['pathogenecity_score']=predicitions
    variants_final_omim_info=variants_final_omim_info.sort_values('pathogenecity_score', ascending=False)
    variants_final_omim_info['rank']=variants_final_omim_info['pathogenecity_score'].rank(ascending=False, method='dense')
    
    variants_final_omim_info.to_csv(os.path.join(output, 'rank_var.tsv'), sep='\t', index=False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Parse arguments for genetic analysis software.")

    parser.add_argument("--annovar", type=str, help="Path to the annotated VCF file", required=True)
    parser.add_argument("--output", type=str, help="Path to the output folder", required=True)
    parser.add_argument("--hpo_ids", type=str, help="Path to the HPO id file", required=True)
    parser.add_argument("--phen2gene", type=str, help="Path to the Phen2Gene score file", required=True)
    parser.add_argument("--gq", type=float, help="Genotype Quality threshold (float)", default=20)
    parser.add_argument("--ad", type=float, help="Allelic Depth threshold (float)", default=15)
    parser.add_argument("--gnomad", type=float, help="GnomAD frequency threshold (float)", default=0.0001)

    args = parser.parse_args()
    
    if not args.output:
        args.output=os.getcwd()
    
    os.makedirs(args.output, exist_ok=True)

    main(args.annovar, args.phen2gene, args.hpo_ids, args.gq, args.ad, args.gnomad, args.output)
