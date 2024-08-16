FROM condaforge/mambaforge as mamba_setup

#ENV HTTP_PROXY="socks5://172.14.60.101:9090"
#ENV HTTPS_PROXY="socks5://172.14.60.101:9090"

RUN mamba create -y -c conda-forge -c bioconda \
    -n vep python=3.8.10 ensembl-vep=108.2 star-fusion=1.11.1
RUN mamba create -y -c bioconda -c conda-forge \
    -n samtools samtools=1.16.1 bismark=0.24.1
RUN mamba install -y -c bioconda -c conda-forge \
	bedtools=2.30.0 \
	fastp=0.23.2 \
	snakemake=7.18.1 \
	gatk4=4.3.0.0 \
	cnvkit=0.9.10 \
	snpeff=5.1 \
	varscan=2.4.4 \
	bwa=0.7.17 \
	msisensor-pro=1.2.0
RUN mamba clean -a

FROM ubuntu:20.04
COPY --from=mamba_setup /opt/conda /opt/conda
LABEL maintainer="mengzheng-li@ebiotron.com"

ENV PATH=/opt/conda/envs/vep/bin:$PATH

RUN apt-get -y update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
                       aria2 \
					   python3 \
					   python3-pip 
RUN apt-get clean autoclean
RUN apt-get autoremove -y --purge

RUN pip3 install pysam==0.19.0 \
                 biopython==1.79 \
                 gffutils==0.11.0 \
                 natsort==8.1.0 \
	             pyyaml==5.3.1 \
                 pandas==1.4.2 \
                 igv-reports==1.6.1 \
				 matplotlib==3.5.2 \
	             -i https://pypi.tuna.tsinghua.edu.cn/simple/

# https://pypi.tuna.tsinghua.edu.cn/simple/
WORKDIR /opt/conda/bin
RUN aria2c http://opengene.org/gencore/gencore
RUN chmod a+x ./gencore

ARG CTAT_SPLICING_VERSION=v0.0.2
RUN aria2c https://github.com/NCIP/CTAT-SPLICING/releases/download/CTAT-SPLICING-${CTAT_SPLICING_VERSION}/CTAT-SPLICING.${CTAT_SPLICING_VERSION}.FULL.tar.gz && tar zxf CTAT-SPLICING.${CTAT_SPLICING_VERSION}.FULL.tar.gz && rm -f CTAT-SPLICING.${CTAT_SPLICING_VERSION}.FULL.tar.gz

RUN aria2c https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/CANCER_SPLICING_LIB_SUPPLEMENT/cancer_introns.GRCh38.Jun232020.tsv.gz

RUN aria2c https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz

RUN aria2c https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_assembly_report.txt

ENV PATH="${PATH}":/opt/conda/bin:/opt/conda/envs/samtools/bin:/opt/conda/bin/CTAT-SPLICING.${CTAT_SPLICING_VERSION}:/opt/conda/bin/CTAT-SPLICING.${CTAT_SPLICING_VERSION}/prep_genome_lib:/snakemake
#ENV PATH="${PATH}":/opt/conda/envs/vep/bin:/opt/conda/bin:/opt/conda/envs/samtools/bin:/opt/conda/bin/CTAT-SPLICING.${CTAT_SPLICING_VERSION}:/opt/conda/bin/CTAT-SPLICING.${CTAT_SPLICING_VERSION}/prep_genome_lib:/snakemake
ENV PATH=/opt/conda/envs/vep/bin:$PATH

COPY . /snakemake

RUN chmod 755 -R /snakemake

# CMD ["AtomSeqTools.py", "-h"]
