# AtomSeqTools
Run AtomSeq Target Pipeline Toolkit

## Requirements
- At least 16 cores CPU and 40G memory, and only analyse for GRCh38/hg38
- Singularity has been installed, [singularity user guide](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps) 
- Downloaded AtomSeqToolsDatabase
- Downloadeded AtomSeqTools's image, command:
```
singularity build AtomSeqTools_image_v2.8.sif docker://leemzheng/atomseqtools:v2.8
```
- Need a config file, fill it out like this:
```
Hg38_Fasta_Path=/PATH/to/GRCh38/hg38.fasta
Fusion_data_library=/PATH/to/AtomSeqToolsDatabase/Fusion_library
Variant_library=/PATH/to/AtomSeqToolsDatabase/Variant_library
```
- Before using AtomSeqTools for the first time, you need to establish hg38 alignment index using the command line:
```
python3 bin/make_index.py \
-image atomSeqTools_images_v2.8.sif \
-fasta /PATH/to/GRCh38/hg38.fasta
```

## Mode one(call variant or cnv or msi):

### 1) only variant calling:
```
python3 AtomSeqTools.py \
-image atomSeqTools_images_v2.8.sif \
-fq-dir /PATH/to/Rawdata \
-fq-prefix 0320-ZCZ-T27-A1 0327-ZLY-T27-A1 \
-outdir /PATH/to/Output \
-config configure_file \
-bed T27_DNA_25geneAndMSI_Primer_V1.1.bed \
-variant -sample-type Tissue
```
### 2) only cnv analysis:
```
python3 AtomSeqTools.py \
-image atomSeqTools_images_v2.8.sif \
-fq-dir /PATH/to/Rawdata \
-fq-prefix 0320-ZCZ-T27-A1 0327-ZLY-T27-A1 \
-outdir /PATH/to/Output \
-config configure_file \
-bed T27_DNA_25geneAndMSI_Primer_V1.1.bed \
-cnv -cnv-baseline CNV_baseline_file
```
### 3) only msi analysis:
```
python3 AtomSeqTools.py \
-image atomSeqTools_images_v2.8.sif \
-fq-dir /PATH/to/Rawdata \
-fq-prefix 0320-ZCZ-T27-A1 0327-ZLY-T27-A1 \
-outdir /PATH/to/Output \
-config configure_file \
-bed T27_DNA_25geneAndMSI_Primer_V1.1.bed \
-msi -mis-baseline MSI_baseline_file
```
### 4) variant, cnv and msi can be used in combination:
```
python3 AtomSeqTools.py \
-image atomSeqTools_images_v2.8.sif \
-fq-dir /PATH/to/Rawdata \
-fq-prefix 0320-ZCZ-T27-A1 0327-ZLY-T27-A1 \
-outdir /PATH/to/Output \
-config configure_file \
-bed T27_DNA_25geneAndMSI_Primer_V1.1.bed \
-variant -sample-type Tissue -cnv -cnv-baseline CNV_baseline_file -msi -mis-baseline MSI_baseline_file
```
## Mode two(methylation analysis):
```
python3 AtomSeqTools.py \
-image atomSeqTools_images_v2.8.sif \
-fq-dir /PATH/to/Rawdata  \
-fq-prefix 0408-V17-2-M5-A1 0408-V17-5-M5-A1 \
-outdir /PATH/to/Output \
-config configure_file \
-bed M5_LC_CRC_Methylation_Primer_V2.3.bed \
-meth
```
## Mode three(fusion analysis):
```
python3 AtomSeqTools.py \
-image atomSeqTools_images_v2.8.sif \
-fq-dir /PATH/to/Rawdata \
-fq-prefix 0320-ZJQ-T12-A1 \
-outdir /PATH/to/Output \
-config configure_file \
-bed T12_RNA_15gene_Primer_V2.0.bed \
-fusion
```

## Help Message:
```
python3 AtomSeqTools.py -h
usage: AtomSeqTools.py [-h] -image IMAGE -fq-dir PATH -fq-prefix STR [STR ...] -outdir PATH -config STR -bed STR
                       [-variant] [-sample-type STR] [-cnv] [-cnv-baseline STR] [-msi] [-msi-baseline STR]
                       [-meth] [-fusion] [-threads INT] [-data-to-analyse FLOAT] [-supporting-reads INT]
                       [-auto-remove]

Run AtomSeq Target Pipeline Toolkit

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -image IMAGE              Singularity image [*.sif]
  -fq-dir PATH              The path to directory of raw fastq
  -fq-prefix STR [STR ...]
                            Specify fastq prefix to be analyzed(can be many, delimited by space)
                            [eg: '0101-XXX-M3-A1_1.fq.gz' prefix is '0101-XXX-M3-A1']
  -outdir PATH              The path to directory of all output files
  -config STR               Config file. Some parameters and paths that need to be set
  -bed STR                  Input bed file of primers

Variant options:
  -variant                  Enable pipeline to call variant
  -sample-type STR          If call variant, it`s necessary to select sample type. Possible values: {cfDNA, Tissue}
  -cnv                      Enable pipeline to call cnv (cfDNA can`t call cnv)
  -cnv-baseline STR         If call cnv, it`s best to give a cnv baseline, but can be None. [default: None]
  -msi                      Enable pipeline to call msi (cfDNA can`t call msi)
  -msi-baseline STR         If call msi, it must give a msi baseline

Methylation options:
  -meth                     Enable pipeline to analyse methylation data

Fusion options:
  -fusion                   Enable pipeline to analyse RNA fusion data

Other options:
  -threads INT              Number of threads to use [default:12]
  -data-to-analyse FLOAT
                            Specify how many data size(G) to be analysed 
                            [default:0] [0 means analyse all fastq data]
  -supporting-reads INT
                            Only output consensus reads/pairs that merged by >= <supporting-reads> reads/pairs.
                            The valud should be 1~10, and the default value is 1. [default: 1]
  -auto-remove              Remove all temporary files (will remove Mid-files)

The bed file of primers should contain at least first six columns(delimited by table):
1. Chormosome ID  (consistent with reference genome);
2. Start of primer (0-based);
3. End of primer (1-based);
4. Primer ID;
5. Primer length;
6. Targeted strand (+/-).
```
## Contacts
If you have any questions or feedback, please contact us at: \
**Email:** mengzheng-li@ebiotron.com
