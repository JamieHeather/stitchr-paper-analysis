# stitchr-paper-analysis
## v 0.1.0
### 2021, Jamie Heather, MGH

This repo contains the raw data and scripts required to generate the plots for the Heather *et al.* manuscript describing [stitchr](https://github.com/JamieHeather/stitchr), the Python tool for generating full-length, coding T cell receptor (TCR) sequences out of the minimally-reported TCR V/J/CDR3 that is most often available.

Note that this version of the repo details the analysis for the `Stitchr` manuscript **before peer review**.

# Requirements

This analysis is designed to be run in a bash environment, with the following installed:

* python 3 (>= 3.6)
  * seaborn, numpy, pandas
* git 
* fastq-dump

# Raw data

There are four sources of bulk TCRseq data used to demonstrate the high-throughput capabilities of the `stitchr` suite of tools:

### VDJdb
* Manually annotated TCRs of known specificity, annotated from across the literature and various other databases.
* Represents a variety of high- and low-throughput TCR sequencing techniques, performed by a variety of labs, and thus should be as unbiased as one could hope to get with respect to potential methodological confounders.
* TRA and TRB chains.
* V gene, J gene, and CDR3 junction amino acid sequence only available.
* Links:
  * [Database: vdjdb.cdr3.net](https://vdjdb.cdr3.net/)
  * [Release used: 2021-02-02](https://github.com/antigenomics/vdjdb-db/releases/tag/2021-02-02) (accessed on 2021-08-19)
  * [Paper DOI: 10.1093/nar/gkz874](https://doi.org/10.1093/nar/gkz874)
  
### Emerson et al
* Healthy donor peripheral blood bulk TCRseq.
* Adaptive Biotechnologies' proprietary multiplex PCR from gDNA.
* TRB only.
* V gene, J gene, CDR3 junction amino acid and nucleotide sequence (87 nt total) available.
* Randomly selected 5x samples (HIP02873, HIP02805, HIP02811, HIP02820, and HIP02855).
* Links:
    * [Paper DOI: 10.1038/ng.3822](https://doi.org/10.1038/ng.3822)
    * [Data DOI: 10.21417/B7001Z](https://doi.org/10.21417/B7001Z)

### immuneSIM
* Artificially simulated TCR sequences.
* Uses IMGT-GENE/DB data and *in silico* simulated V(D)J recombination, producing TCRs of known properties lacking any sequencing errors.
* Links:
    * [Paper DOI: 10.1093/bioinformatics/btaa158](https://doi.org/10.1093/bioinformatics/btaa158)
    * [Github repo: github.com/GreiffLab/immuneSIM](https://github.com/GreiffLab/immuneSIM)
    * [CRAN package: cran.r-project.org/package=immuneSIM](https://cran.r-project.org/package=immuneSIM)
    
### Heather et al
* Heather donor peripheral blood bulk TCRseq.
* 5'RACE, UMI-barcoded amplification from cDNA, paired-end sequenced 2x250 on a MiSeq.
* A fraction of the amplicons where R1 and R2 overlapped (presumably those where the RT reaction ran far enough) were merged using [FLASH](https://doi.org/10.1093/bioinformatics/btr507), producing long-read variable-domain spanning reads.
* Most from a previous publication of the same lead author (JMH), with some additional samples from a later paper.
* Links:
  * [Original paper DOI: 10.3389/fimmu.2015.00644](https://doi.org/10.3389/fimmu.2015.00644)
  * [Follow up paper DOI: 10.3389/fimmu.2017.01267](https://doi.org/10.3389/fimmu.2017.01267)
* Raw merged FASTQ data available from SRA
    * See table below.

| Donor | Locus | SRA container accession | SRA run accession | 
|:-----:|:-----:|:-----------------------:|:-----------------:|
| HV01  |  TRA  |       SRX2454454        |    SRR5137043     |
| HV01  |  TRB  |       SRX2454419        |    SRR5137008     |
| HV02  |  TRA  |       SRX2454424        |    SRR5137013     |
| HV02  |  TRB  |       SRX2454417        |    SRR5137006     |
| HV03  |  TRA  |       SRX2454426        |    SRR5137015     |
| HV03  |  TRB  |       SRX2454449        |    SRR5137038     |
| HV04  |  TRA  |       SRX2454433        |    SRR5137022     |
| HV04  |  TRB  |       SRX2454448        |    SRR5137037     |
| HV07  |  TRA  |       SRX2454440        |    SRR5137029     |
| HV07  |  TRB  |       SRX2454406        |    SRR5136995     |
| HV08  |  TRA  |       SRX2454414        |    SRR5137003     |
| HV08  |  TRB  |       SRX2454434        |    SRR5137023     |
| HV09  |  TRA  |       SRX2454431        |    SRR5137020     |
| HV09  |  TRB  |       SRX2454415        |    SRR5137004     |
| HV10  |  TRA  |       SRX2454412        |    SRR5137001     |
| HV10  |  TRB  |       SRX2454459        |    SRR5137048     |
| HV11  |  TRA  |       SRX2454402        |    SRR5136991     |
| HV11  |  TRB  |       SRX2454444        |    SRR5137033     |
| HV13  |  TRA  |       SRX2454410        |    SRR5136999     |
| HV13  |  TRB  |       SRX2454442        |    SRR5137031     |
| HV15  |  TRA  |       SRX11859184       |    SRR15561076    |
| HV15  |  TRB  |       SRX11859192       |    SRR15561068    |
| HV16  |  TRA  |       SRX11859185       |    SRR15561075    |
| HV16  |  TRB  |       SRX11859193       |    SRR15561067    |
| HV17  |  TRA  |       SRX11859188       |    SRR15561072    |
| HV17  |  TRB  |       SRX11859194       |    SRR15561066    |
| HV18  |  TRA  |       SRX11859189       |    SRR15561071    |
| HV18  |  TRB  |       SRX11859195       |    SRR15561065    |
| HV19  |  TRA  |       SRX11859190       |    SRR15561070    |
| HV19  |  TRB  |       SRX11859186       |    SRR15561074    |
| HV20  |  TRA  |       SRX11859191       |    SRR15561069    |
| HV20  |  TRB  |       SRX11859187       |    SRR15561073    |

# Other scripts used

Note that the most up-to-date versions of these scripts are downloaded using git.  

### immunoseq2airr
* v1.2.0
* Converts Adaptive's immunoSEQ data into the AIRR-Community standard data format, including correcting its novel TCR gene nomenclature.
* [Release DOI: 10.5281/zenodo.3770611](https://doi.org/10.5281/zenodo.3770611)
* [Repo: github.com/JamieHeather/immunoseq2airr](https://github.com/JamieHeather/immunoseq2airr)

### autoDCR
* v0.1.0
* Annotates recombined TCR sequences.
  * Inspired by the core functionality of Decombinator (the TCR analysis software used in the above linked Chain group TCRseq studies).
  * Works off a similar conceptual framework (using Aho-Corasick tries to search for 'tag' sequences linked to V/J genes). 
  * Sacrifices the speed of Decombinator to allow for maximum-resolution allele level distinction across the length of a read. Given a FASTA of all V and J genes for a set of TRA/TRB loci, overlapping tags are generated across the entirety of all genes, and V/J genes used determined based on the maximum number of uniquely identifying tags, rather than the specific discovery of single manually-defined tags which uniquely identify a given gene proximal to the CDR3. This approach not only allows for greater inferences about CDR3-distal features, but automating tag (hence '*auto*DCR') production makes the process more robust to the discovery of novel alleles. 
  * [Original paper DOI: 10.1093/bioinformatics/btt004](https://doi.org/10.1093/bioinformatics/btt004)
  * [Recent update DOI: 10.1093/bioinformatics/btaa758](https://doi.org/10.1093/bioinformatics/btaa758)
  * [Repo: github.com/JamieHeather/autoDCR](https://github.com/JamieHeather/autoDCR)
  
# Running the analysis

Upon cloning this repo, there should be the following directories:
* Data/
  * Contains the raw VDJdb, Emerson, and immuneSIM data files.
* Outputs/
  * Empty; results produced by running the analysis will go here.
* Scripts/
  * Contains the various Python scripts to gather/generate/analyse/plot the data.

Scripts then need to be executed in the following order:

* `preparation.py`
  * This downloads the Heather *et al.* raw data from SRA, and converts the raw data from all data sources into a suitable format for processing with `thimble`, the high-throughput wrapper for `stitchr`.
  * Note that this must be run in a bash (or similar) environment, with [`fastq-dump` from `sra-tools`](https://ncbi.github.io/sra-tools/fastq-dump.html) installed
  * This script also uses `git` to download both the `stitchr` files, and the other required supplementary scripts described above.
  
* `stitching.py`
  * This submits the `thimble` format files generated above for processing.
  * It also repeatedly subsamples the third party datasets (VDJdb/Emerson) to generate a dynamic range of filesizes for benchmarking `thimble`'s speed.
  
* `plotting.py`
  * This reads in the stitched data and run metrics to produce the bulk of the plots that relate to the performance of `stitchr` and `thimble`

* `alleles.py`
  * This automates a rudimentary novel allele detection process on the Heather *et al*. TCRseq data, using additional information available from the output of `autoDCR`. 

* `replotting.py`
  * This uses the various functions established in `stitching.py` and `plotting.py` to repeat the analysis and plotting of the Heather *et al*. data using the IMGT germline TCR gene reference supplemented with potential novel alleles detected by `alleles.py`.

Scripts need to be executed from inside the Scripts directory, as relative paths to there are used throughout. Everything needs to be run in Python3, and has been tested on >= 3.6.9. (Note that the `functions.py` script should not be directly executed, as it simply contains objects used by the other three scripts.)
