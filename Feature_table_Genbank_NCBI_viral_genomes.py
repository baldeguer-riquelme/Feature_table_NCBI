from Bio import SeqIO
from functools import reduce
import pandas as pd
import argparse
import logging
import os

pd.options.mode.chained_assignment = None

# command to create conda environment with all needed packages: 
# conda create --name feature_table -c bioconda -c conda-forge prodigal trnascan-se pandas biopython openjdk=11.0.8

# To install interproscan follow the instructions on the web site https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html

# 0. Calculate genome length
def calculate_genome_length(genome_file):
    for record in SeqIO.parse(genome_file, "fasta"):
        genome_length = len(record.seq)
        return(genome_length)

# 1. Predict tRNAs with tRNAScan-SE if not provided by user
def predict_tRNAs(genome, genome_file):
    logging.info("Running tRNAScan-SE...")
    trnascan = os.system('tRNAscan-SE -G -j {0}_tRNAScan-SE.gff -Q {1}  2>&1 | tee -a Feature_table_{0}.log > /dev/null'.format(genome, genome_file))
    logging.info("tRNAScan-SE finished successfully")

# 2. Analyze tRNAs
def calculate_tRNAs(trna_file, genome):
    logging.info("Reading tRNAScan-SE output file...")
    if os.stat(trna_file).st_size > 0:
        logging.info("tRNAs found in genome")
        trnascan_res=pd.read_csv(f"{genome}_tRNAScan-SE.gff",sep="\t",header=None, skiprows=[0])
        trnascan_res.columns=["contig", "tool", "type", "start", "stop", "score", "strand", "frame", "attribute"]
        only_trnas=trnascan_res[trnascan_res['type'] == "tRNA"]
        only_trnas['product'] = only_trnas['attribute'].str.split(";", expand=True)[2].str.replace("isotype=","")
        only_trnas['product'] = "tRNA-" + only_trnas['product'].astype(str)
        summary_trnas = only_trnas[["type","start","stop","strand","product"]]
        summary_trnas.insert(0,"Gene","tRNA")
        for idx in 0,(len(summary_trnas)-1):
            summary_trnas.loc[idx, ["Gene"]] = "tRNA-" + str(idx+1)
        summary_trnas.insert(5,"partial","00")
        summary_trnas['strand'] = summary_trnas['strand'].str.replace("+","1", regex=True)
        summary_trnas['strand'] = summary_trnas['strand'].str.replace("-","-1", regex=True)
        return(summary_trnas)
    else:
        logging.info("No tRNAs detected, skipping step...")
        summary_trnas = pd.DataFrame(columns=['Gene', 'type', 'start', 'stop', 'strand', 'partial', 'product'])
        return(summary_trnas)

# 3. Predict CDS with prodigal if not provided by user
def predict_CDS(genome, genome_file):
    logging.info("Running Prodigal...")
    prodigal = os.system('prodigal -i {0} -a {1}.faa -p meta  2>&1 | tee -a Feature_table_{1}.log > /dev/null'.format(genome_file, genome))
    rm_asterisk = os.system("sed 's/*//g' {0}.faa > tmp ; mv tmp {0}.faa".format(genome))
    logging.info("Prodigal finished successfully")

# 4. Annotate CDS with interproscan
def annotate(fasta_file, genome):
    logging.info("Running interproscan...")
    ipro = os.system('interproscan.sh -i {0} -b ipro_{1} -appl CDD,Pfam,SMART,TIGRFAM -t p -T ./temp_{1}/. >> Feature_table_{1}.log'.format(fasta_file, genome))
    raw_tab = pd.read_table("ipro_{0}.tsv".format(genome), sep="\t", header=None)
    raw_tab_clean = raw_tab.drop([1,2,4,6,7,8,9,10,11,13,14],axis=1)
    pfam = raw_tab_clean[raw_tab_clean[3] == "Pfam"].drop(12, axis=1)
    pfam.columns=["Tag","DB","Product"]
    CDD = raw_tab_clean[raw_tab_clean[3] == "CDD"].drop(12, axis=1)
    CDD.columns=["Tag","DB","Product"]
    TIGRFAM = raw_tab_clean[raw_tab_clean[3] == "TIGRFAM"].drop(12, axis=1)
    TIGRFAM.columns=["Tag","DB","Product"]
    SMART = raw_tab_clean[raw_tab_clean[3] == "SMART"].drop(5, axis=1)
    SMART.columns=["Tag","DB","Product"]
    to_merge = [pfam, CDD, TIGRFAM, SMART]
    merged_tab = pd.concat(to_merge)
    detected_dbs = merged_tab["DB"].unique().tolist()
    detected_dbs.insert(0, 'Tag')
    final_annot = merged_tab.pivot(columns='DB', index='Tag', values='Product').reset_index()
    final_annot = final_annot[detected_dbs]
    final_annot.to_csv("ipro_{0}_final.tsv".format(genome),sep='\t', index=False)
    logging.info("Interproscan finished successfully")

# 5. Analyze IDs from prodigal output to get coordinates and merge with annotation table
def calculate_CDS(fasta_file, annotation_file):
    logging.info("Reading gene/protein fasta file...")
    list_ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        list_ids.append(record.description)

    logging.info("Processing gene/protein fasta file to get coordinates...")
    gene_df = pd.DataFrame(list_ids)
    gene_df_expand = gene_df[0].str.split(" |;", expand=True)

    summary_genes = gene_df_expand[[0,2,4,6,9]]
    summary_genes.columns = ['Gene','start', 'stop','strand','partial']
    summary_genes.insert(1, "type", "CDS")
    summary_genes['partial'] = summary_genes['partial'].str.replace("partial=","")
    logging.info("Processing gene/protein fasta file done")

    # 5.1 Merge prodigal genes with annotation table
    logging.info("Reading annotation file...")
    annot = pd.read_table(annotation_file, header = 0, sep = '\t')
    logging.info("Processing anotation file...")
    annot['product'] = annot.bfill(axis=1).iloc[:, 1]
    final_annot = annot[['Tag','product']]
    summary_genes_product = summary_genes.merge(final_annot, left_on="Gene", right_on="Tag", how="left")
    summary_genes_product['product'] = summary_genes_product['product'].fillna("Hypothetical protein")
    summary_genes_product = summary_genes_product.drop('Tag', axis=1)
    logging.info("Processing annotation file done")
    return(summary_genes_product)
    

# 6. Merge both tRNAs and prodigal dataframes
def merged_dfs(summary_genes, summary_trnas):
    logging.info("Merging tRNAs and CDS dataframes...")
    final_df = pd.concat([summary_genes,summary_trnas], ignore_index = True)
    final_df['start'] = final_df['start'].astype(int)
    final_df_sorted = final_df.sort_values(by = ['start'])
    logging.info("Merging done")
    return(final_df_sorted)


def complete_genes(row):
    gene = row.Gene
    type=row.type
    start = row.start
    stop = row.stop
    product=row.product
    direction = row.strand
    if direction == "1" and type == "CDS":
        print('\t'.join([str(start), str(stop), type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'transl_table\t', '11']))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "-1" and type == "CDS":
        print('\t'.join([str(stop), str(start), type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'transl_table\t', '11']))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "1" and type == "tRNA":
        print('\t'.join([str(start), str(stop), type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "-1" and type == "tRNA":
        print('\t'.join([str(stop), str(start), type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))


def partial_genes_left(row):
    gene = row.Gene
    type=row.type
    start = row.start
    stop = row.stop
    product=row.product
    direction = row.strand
    if direction == "1" and type == "CDS":
        print(''.join(["<1\t", str(stop), "\t", type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'transl_table\t', '11']))
        print(''.join(['\t', '\t', '\t', 'codon_start\t', str(start)]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "-1" and type == "CDS":
        print(''.join([str(stop),"\t", ">1\t", type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'transl_table\t', '11']))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "1" and type == "tRNA":
        print(''.join(["<1\t", str(stop), "\t", type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'codon_start\t', str(start)]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "-1" and type == "tRNA":
        print(''.join([str(stop),"\t", ">1\t", type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))


def partial_genes_right(row, genome_length, fasta_file):
    gene = row.Gene
    type=row.type
    start = row.start
    stop = row.stop
    product=row.product
    direction = row.strand
    if direction == "1" and type == "CDS":
        print(''.join([str(start),"\t>", str(genome_length), "\t",type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'transl_table\t', '11']))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "-1" and type == "CDS":
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.name == gene:
                nt_start = (start -1) + len(record.seq)*3 + 3 # start position include the stop codon, this is why I sum + 3
                codon_start = genome_length - nt_start + 1
        print(''.join(["<", str(genome_length), "\t", str(start), "\t",type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'transl_table\t', '11']))
        print(''.join(['\t', '\t', '\t', 'codon_start\t', str(codon_start)]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "1" and type == "tRNA":
        print(''.join([str(start),"\t>", str(genome_length), "\t",type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "-1" and type == "tRNA":
        codon_start = 1
        print(''.join(["<", str(genome_length), "\t", str(start), "\t",type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'codon_start\t', str(codon_start)]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))


def partial_genes_both_ends(row, genome_length, fasta_file):
    gene = row.Gene
    type=row.type
    start = row.start
    stop = row.stop
    product=row.product
    direction = row.strand
    if direction == "1" and type == "CDS":
        print(''.join(["<1\t>", str(genome_length), "\t", type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'transl_table\t', '11']))
        print(''.join(['\t', '\t', '\t', 'codon_start\t', str(start)]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "-1" and type == "CDS":
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.name == gene:
                codon_start = start + len(record.seq)*3
        print(''.join(["<1\t>", str(genome_length), "\t", type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'transl_table\t', '11']))
        print(''.join(['\t', '\t', '\t', 'codon_start\t', str(codon_start)]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "1" and type == "tRNA":
        print(''.join(["<1\t>", str(genome_length), "\t", type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'codon_start\t', str(start)]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))
    elif direction == "-1" and type == "tRNA":
        print(''.join(["<1\t>", str(genome_length), "\t", type]))
        print(''.join(['\t', '\t', '\t', 'gene\t', gene]))
        print(''.join(['\t', '\t', '\t', 'codon_start\t', str(stop)]))
        print(''.join(['\t', '\t', '\t', 'product\t', product]))


def main():
    parser = argparse.ArgumentParser("Generates 'Feature tables' for submission of viral sequences to NCBI. It predicts proteins using Prodigal and annotates them using interproscan. tRNAs are predicted with tRNAScan-SE.")
    parser.add_argument("--genome", 
                        type=str, help="Genome fasta file", required=True)
    parser.add_argument("--prodigal",
                        type=str,
                        help="Prodigal output file (amino acid fasta format)", 
                        required=False)
    parser.add_argument("--trnas",
                        type=str,help="tRNAs annotation file", required=False)
    parser.add_argument("--annotation_cds",
                        type=str,help="CDS annotation file", required=False)
    args = parser.parse_args()

    genome_file = args.genome
    genome = genome_file.split(".")[0]
    logging.basicConfig(filename='Feature_table_{0}.log'.format(genome), encoding='utf-8', 
                        level=logging.DEBUG, format='\n[%(asctime)s] %(message)s\n')

    # Check user-provided arguments and get missing files
    if args.trnas is None:
        predict_tRNAs(genome, genome_file)
        trna_file = f"{genome}_tRNAScan-SE.gff"
    else:
        trna_file = args.trnas
    if args.prodigal is None:
        predict_CDS(genome, genome_file)
        fasta_file = f"{genome}.faa"
    else:
        fasta_file = args.prodigal
    if args.annotation_cds is None:
        annotate(fasta_file, genome)
        annotation_file = "ipro_{0}_final.tsv".format(genome)
    else:
        annotation_file = args.annotation_cds

    # Analyze files to get coordinates and annotations for all genes
    length = calculate_genome_length(genome_file)
    tRNAs = calculate_tRNAs(trna_file, genome)
    CDS = calculate_CDS(fasta_file, annotation_file)
    Merged_dataframe = merged_dfs(tRNAs, CDS)

    # Print feature table to terminal
    logging.info("Printing results to stdout...")
    print(' '.join(['>Feature', genome]))
    for row in Merged_dataframe.itertuples(index=False):
        if row.partial == "00":
            complete_genes(row)
        elif row.partial == "10":
            partial_genes_left(row)
        elif row.partial == "01":
            partial_genes_right(row, length, fasta_file)
        elif row.partial == "11":
            partial_genes_both_ends(row, length, fasta_file)
    logging.info("Feature table built successfully!")


if __name__ == '__main__':
    main()