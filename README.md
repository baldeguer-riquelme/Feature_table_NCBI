# Build "Feature table" for NCBI BankIt submission of viral genomes
Contains scripts and information to build "Feature Tables" using viral genome sequences in fasta format as the only input. 

NCBI asks users who submit viral genome sequences to provide a "Feature Table" (https://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html) but I haven't found an easy and fast way to build them. If you are facing the same problem, in this repository you will find a script and instructions to get your "Feature Tables" in a matter of seconds using your genome fasta file as the unique input.


### Dependencies

[tRNAScan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE) (tested with v2.0.12) 

[Prodigal](https://github.com/hyattpd/Prodigal) (tested with v2.6.3) 

[Interproscan](https://interproscan-docs.readthedocs.io/en/latest/index.html) (tested with v5.65-97.0) 

[Pandas](https://pandas.pydata.org/) (tested with v2.1.4) 

[Biopython](https://biopython.org/) (tested with v1.81) 

[Java](https://www.oracle.com/java/technologies/javase/jdk11-archive-downloads.html) (required v11.0.8) 

 
### Installation of dependencies

The easiest way to install dependencies is creating a new conda environment:

```conda create --name feature_table -c bioconda -c conda-forge prodigal trnascan-se pandas biopython openjdk=11.0.8```

We will also need to install Interproscan following the instructions available at https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html .
For your convenience, I have copied them below so you can just copy and paste them.

```
mkdir my_interproscan
cd my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.65-97.0/interproscan-5.65-97.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.65-97.0/interproscan-5.65-97.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.65-97.0-64-bit.tar.gz.md5
# Must return *interproscan-5.65-97.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.

# Uncompress
tar -pxvzf interproscan-5.65-97.0-64-bit.tar.gz
cd interproscan-5.65-97.0
python3 setup.py -f interproscan.properties
```

Then, add interproscan to the path by simply copying all the files to the "bin" folder within the "feature_table" conda environment. For example:

cp -r * /soft/anaconda3/envs/feature_table/bin/

```cp -r * /${path_to_env}/bin/```

### Usage
The easiest way to create your feature table is:

```python Feature_table_NCBI.py --genome virus.fasta > Output.txt```

A ".log" file will be generated which contains information about the process. If you face any error, check it.

### Usage with user-supplied files
Alternatively you might want to provide all or some of the files needed to build the feature table (see the "examples" folder), such as:
1. Proteins predicted from Prodigal in fasta format (amino acids)
2. GFF file with tRNA prediction info from tRNAScan-SE
3. Annotation table (At least two columns, first should be named as "Tag". If you provide multiple columns with the annotation from different databases, the script will give priority to columns from the left to the right)

The script will check if you provide any of these additional files and predict the missing ones.

```python Feature_table_NCBI.py --genome virus.fasta --prodigal virus.faa --trnas virus_tRNAScan-SE.gff --annotation_cds ipro_virus_final.tsv > Output.txt```
