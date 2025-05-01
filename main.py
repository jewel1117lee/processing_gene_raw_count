import pandas as pd
import scipy
from openpyxl import Workbook
from pathlib import Path
from src.clean import(filter_noncoding, filter_zeros, dropna_from_lists,)
from src.normalization import(compute_rpkm, compute_cpm,)
from src.deg_analysis import(run_deg_analysis)
from src.analysis import(match_database,ttest,make_box_figure,zip_data)

# I added this because i am testing on my own server,
# You don't need this if it is online 
BASE = Path(__file__).parent

# load needed data, raw counts are necessary
# you don't need transcript_file nor protein_coding if you are not analyzing them
raw_counts_path = BASE / "data" / "table_raw_count_for_all_samples_with_txID_gene_name.xlsx"
protein_coding_path = BASE / "data" / "protein coding gene.tsv"
transcript_file_path = BASE / "data" / "mouse_table_DF.csv"

# raw_df is general .xlsx, protein_coding is .tsv, transcript_file is .csv
raw_df = pd.read_excel(raw_counts_path)
protein_coding = pd.read_csv(protein_coding_path, sep='\t')
transcript_file = pd.read_csv(transcript_file_path)

# Set_output path to export everything you run
output_path = BASE / "result"

# calling filter_noncoding
# This step is necessary if you want to calculate for rpkm and etc
df = filter_noncoding(raw_df, protein_coding, "gene_name", 	"external_gene_name")

# calculating rpkm 
"""
This function will calculate for all of columns in the rpkm if you do not specify it
RPKM : One of the normalization rate where  RPKM = raw_count / ( geneLength/1000 * total_reads/1,000,000 )
Input: (df= dataframe, 
        gene_id_cols = the columns that you want to keep for new rpkm file, eg. gene_id, gene_name
        gene_length_col= name of the gene length col,
        sample_cols= the column or col lists you want to run, if non it will run all)
"""
# 1. Unfortunately, I do not have transcript length in my processing file, so I am going to merge the two excels and then run it
df_w_transcript = ensembl_df1= pd.merge(df, transcript_file[["transcript_id","Transcript_Length"]], on='transcript_id')
# 2. Run rpkm_df
rpkm_df = compute_rpkm(df_w_transcript , gene_id_cols= ["gene_name","gene_id"], gene_length_col="Transcript_Length")
# 3. Export
#rpkm_df.to_excel(output_path + "/rpkm.xlsx")

# calculating DEG Analysis
"""
This function will calculate the Differential express of the gene using package from pydeseq2
Input (
    counts_df = pd.DataFrame, the data you want to run 
    gene_id_col = str, name of the gene_id_col
    group1_cols = List[str], the list of columns you want to run for DEG for the first group
    group2_cols = List[str], the list of columns you want to run for DEG for the first group
    g1_label: str = str, the variable name for the first group
    g2_label: str = str, the variable name for the second group
    gene_name_col: str = 'gene_name',
    normalize_factor = False, this calls for a built in normalize function from pydeseq2)
        specifically it  mplements the median-of-ratios normalization
"""
# 1. Run the Test
test_DEG = run_deg_analysis(
    df,
    "transcript_id",
    ["FMR1_Granule_rep1","FMR1_Granule_rep2","FMR1_Granule_rep3"],
    ["WT_Granule_rep1","WT_Granule_rep2","WT_Granule_rep3"],
    g1_label = "FMR1",
    g2_label = "WT",
)

# 2. Export the DEG Results 
# test_DEG.to_excel(output_path + "/DEG.xlsx")

# Matching to online database 
"""
Today you are trying to see if your data matches to specific genes that you found online, you can use this function to if your corresponding gene is above or below
Input = (
    gene_names (list of str): List of gene names.
    values (list of numeric): Corresponding values for each gene (e.g., log-RPKM).
    comp_df (pd.DataFrame): Published data containing a column comp_col with gene names.
    comp_col (str): Column name in comp_df to match against gene_names.
    protein_coding_df (pd.DataFrame): DataFrame with 'gene_name' column listing protein-coding genes.    
    )
Output = (
    matched_genes: list(str), your gene list that matched to the comparison gene list
    matched_vals:  list(int, float), your correspondent gene value list that matched to the comparison gene list
    unmatched_genes: list(str),  your gene list that did not matched to the comparison gene list
    unmatched_vals: list(int, float), your correspondent gene value list that did not matched to the comparison gene list
)
"""
# 1. Importing comparison files
path_comparison = "C:/Users/lijew/Desktop/Sossin Lab/FMR1 project/data from Teodora/"
comparison_filename = "Teodora encompass.xlsx"
comparison_file = pd.read_excel(path_comparison + comparison_filename, sheet_name = "former gene list")

# 2. Run comparison
comparison_test = match_database(test_DEG, "gene_name", "log2FoldChange", comparison_file, "Darnell")

"""
The following runs t-test with unequal variance. It returns the t_stat, p_value, and how many stars for graphing
Input = (
    array1 = list or arrays,
    array2 = list or arrays
)
Output = (
    t_stat
    p_value
    star: statistical number of stars used in publication
)
"""
# 3. Run t-test
t_stat, p_value, star = ttest(comparison_test [1], comparison_test [3])

# 4. graphs
# 4.1 you should drop all of the NaN values before you make a graph, and since my graph code is in list of list, I wrote a function for that
data_for_graph = dropna_from_lists([comparison_test[1], comparison_test[3]])

# 4.2 graphing
"""
This creates boxplot for all the data that you had compared to the online genes
Inputs = (
    title (str): Plot title.
    data (list of lists): Two lists of numeric values.
    labels (list of str): Labels for the two boxes.
    xlabel, ylabel (str): Axis labels.
    colors (list of str): Face colors for the boxes.
    median_color (str): Color for median lines.
    ttest (float or str): If provided add the p-value or stars on top of the box plot
    output_path (str): If provided, saves plot to this path.
    )

"""
make_box_figure(
    title = "test", data = data_for_graph, labels = ["Darnell", "Darnell non"] , xlabel = "dataset", ylabel = "enrichment",
    colors=None, median_color='black', ttest = [star, ""])

#5. Export the matched and non_matched data
#5.1 because the matched and non_matched data are often unequal in length so you need to make it equal with NaN to export
"""
def zip_data
Input = (
    data = data you want to zip together
    col_names= if you have it it adds column names
)
"""
test_zip_df = (data_for_graph, ["matched", "unmatched"])
# test_zip_df.to_excel(output_path + "/comparison.xlsx")

