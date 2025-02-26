# RankVar
RankVar is an AI-driven pipeline that integrates phenotype data and sequencing profiles to prioritize disease-causing genes and variants.

# Installation

We recommend using Conda to set up the environment. If Conda is not installed, run the following commands in Linux to install it.

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
After conda is installed successfully, RankVar sources can be downloaded:

```bash
git clone https://github.com/WGLab/RankVar.git
cd RankVar
```

# Inference
### Step 1: Install and run ANNOVAR
ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes.

#### a) Install ANNOVAR

Typically you will go to the [ANNOVAR website](https://annovar.openbioinformatics.org/en/latest/), fill in a registration form, and download the package there. When you have requested the ANNOVAR from the website and downloaded it, you will receive a compressed file ```annovar.latest.tar.gz```, you will need to unzip it. Then follow the user guide to install ANNOVAR. 

#### b) Run ANNOVAR

Input files to ANNOVAR refer to VCF file (example.vcf)

```bash
perl table_annovar.pl example.vcf humandb/ -buildver hg38 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp47a,gnomad41_exome,gnomad41_genome,clinvar_20240917,eQTL,sQTL -operation gx,r,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish
```
After that, you will find the result files ```myanno.hg38_multianno.txt``` and ```myanno.hg38_multianno.vcf```.

### Step 2: Install and run Phen2Gene
Phen2Gene is a phenotype-driven gene prioritization tool, that takes HPO (Human Phenotype Ontology) IDs as inputs, searches and prioritizes candidate causal disease genes.

#### a) Install Phen2gene

```bash
git clone https://github.com/WGLab/Phen2Gene.git
cd Phen2Gene
conda env create -f environment.yml
conda activate phen2gene
bash setup.sh
```
#### b) Run Phen2Gene

Input files to Phen2Gene should contain HPO IDs, separated by UNIX-recognized new line characters (i.e., \n). Alternatively you can use a space separated list of HPO IDs on the command line.

Here is an example file called ```hpo_list.txt```
```bash
HP:0000358
HP:0000039
HP:0008438
HP:0000891
HP:0000252
```
simply run:
```bash
python3 phen2gene.py -f example/hpo_list.txt -out phen2gene_out
```
After that, you will find the result files ```phen2gene_out/output_file.associated_gene_list```

### Step 3: Run RanVar in linux

Input files to RankVar are annotated VCF file (```myanno.hg38_multianno.txt```) and HPO related Phen2gene score file (```phen2gene_out/output_file.associated_gene_list```)

Type python RankVar/RankVar.py -help to see all options.
```bash
usage: RankVar.py [-h] --annovar ANNOVAR --output OUTPUT --hpo_ids HPO_IDS --phen2gene PHEN2GENE [--gq GQ] [--ad AD]
                  [--gnomad GNOMAD]

Parse arguments for genetic analysis software.

options:
  -h, --help            show this help message and exit
  --annovar ANNOVAR     Path to the annotated VCF file (default: None)
  --output OUTPUT       Path to the output folder (default: None)
  --hpo_ids HPO_IDS     Path to the HPO id file (default: None)
  --phen2gene PHEN2GENE
                        Path to the Phen2Gene score file (default: None)
  --gq GQ               Genotype Quality threshold (float) (default: 20)
  --ad AD               Allelic Depth threshold (float) (default: 15)
  --gnomad GNOMAD       GnomAD frequency threshold (float) (default: 0.0001)
```

#### Example

Download the example VCF files:
```bash
wget https://pmc.ncbi.nlm.nih.gov/articles/instance/5111005/bin/supp_mcs.a001131_Supp_File_2_KBG_family_Utah_VCF_files.zip
unzip supp_mcs.a001131_Supp_File_2_KBG_family_Utah_VCF_files.zip
```
After that, you will find the ```proband.vcf```

I am using GATK to convert ```proband.vcf``` from hg19 to hg38 and generate the ```proband.hg38.vcf```:

```bash
gatk --java-options "-Xmx16g" LiftoverVcf -I proband.vcf -O proband.hg38.vcf -CHAIN hg19ToHg38.over.chain.gz -REJECT unmapped_variants.vcf -R Homo_sapiens_assembly38.fasta
```

then, run annovar on ```proband.hg38.vcf``` and Phen2gene on ```hpo_list.txt``` to generate the files ```myanno.proband.hg38_multianno.txt``` and ```phen2gene_out/output_file.associated_gene_list```

run RankVar:
```bash
python RankVar.py --annovar myanno.proband.hg38_multianno.txt --phen2gene phen2gene_out/output_file.associated_gene_list  --hpo_ids hpo_list.txt --output $PWD
```
