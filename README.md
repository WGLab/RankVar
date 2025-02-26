# RankVar: Machine Learning-Based Variant Ranking and Reinterpretation for Rare Genetic Diseases
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
conda create -n rankvar python=3.10
conda activate rankvar
pip install numpy pandas joblib scikit-learn==1.3 torch
python RankVar.py --help
```

# Inference
### Step 1: Install and run ANNOVAR
ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes.

#### a) Install ANNOVAR

Typically you will go to the [ANNOVAR website](https://annovar.openbioinformatics.org/en/latest/), fill in a registration form, and download the package there. When you have requested the ANNOVAR from the website and downloaded it, you will receive a compressed file ```annovar.latest.tar.gz```, you will need to unzip it. Then follow the user guide to install ANNOVAR. 

#### b) Run ANNOVAR

Input files to ANNOVAR refer to VCF file (example.vcf)

```bash
perl table_annovar.pl example.vcf humandb/ -buildver hg38 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp47a,gnomad41_exome,gnomad41_genome,clinvar_20240917,GTEx_v8_eQTL,GTEx_v8_sQTL -operation gx,r,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish
```
After that, you will find the result files ```myanno.hg38_multianno.txt``` and ```myanno.hg38_multianno.vcf```.

### Step 2: Install and run Phen2Gene
Phen2Gene is a phenotype-driven gene prioritization tool, that takes HPO (Human Phenotype Ontology) IDs as inputs, searches and prioritizes candidate causal disease genes.

#### a) Install Phen2gene
Please follow [Phen2gene](https://github.com/WGLab/Phen2Gene) repository for instructions on how to install Phen2gene.

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
python3 phen2gene.py -f hpo_list.txt -out phen2gene_out
```
After that, you will find the result files ```phen2gene_out/output_file.associated_gene_list```

### Step 3: Run RanVar

Input files to RankVar are annotated VCF file (```myanno.hg38_multianno.txt```) and HPO terms (```hpo_list.txt```) and related Phen2gene score file (```phen2gene_out/output_file.associated_gene_list```)

Type ```python RankVar/RankVar.py -help``` to see all options.

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

Download the example VCF file:

```bash
wget https://pmc.ncbi.nlm.nih.gov/articles/instance/5111005/bin/supp_mcs.a001131_Supp_File_2_KBG_family_Utah_VCF_files.zip
unzip supp_mcs.a001131_Supp_File_2_KBG_family_Utah_VCF_files.zip
```

After that, you will find the ```proband.vcf```

If the input VCF is in hg19, you need to convert it to hg38 using GATK:

```bash
gatk --java-options "-Xmx16g" LiftoverVcf -I proband.vcf -O proband.hg38.vcf -CHAIN hg19ToHg38.over.chain.gz -REJECT unmapped_variants.vcf -R Homo_sapiens_assembly38.fasta
```
After that, you will find the result file ```proband.hg38.vcf```

Then, run annovar on ```proband.hg38.vcf``` and Phen2gene on ```hpo_list.txt``` to generate the files ```myanno.proband.hg38_multianno.txt``` and ```phen2gene_out/output_file.associated_gene_list```

run RankVar:
```bash
python RankVar.py --annovar myanno.proband.hg38_multianno.txt --phen2gene phen2gene_out/output_file.associated_gene_list  --hpo_ids hpo_list.txt --output output/
```
RankVar will write output in `output/rank_var.tsv` that will look like:
```
Chr    Start      End        Ref  Alt  Func.refGene  Gene.refGene  ExonicFunc.refGene       gnomad41_exome_AF_grpmax  phen2gene_score  pathogenicity_score  rank
chr16  89280526   89280526   -    T    exonic       ANKRD11       frameshift insertion     0.0                      1.0               1.0                  1.0
chr16  27701625   27701625   C    -    exonic       KATNIP        frameshift deletion      0.0                      0.0               1.0                  1.0
chr19  13298600   13298600   G    -    exonic       CACNA1A       frameshift deletion      0.0                      0.085824          0.77                 2.0
chr2   202555679  202555679  -    A    exonic       BMPR2         frameshift insertion     0.0                      0.108457          0.65                 3.0
chr19  53889950   53889950   -    G    exonic       PRKCG         frameshift insertion     0.0                      0.158115          0.59                 4.0
chr12  48966504   48966504   C    -    exonic       WNT10B        frameshift deletion      0.0                      0.1344            0.49                 5.0
...
```
