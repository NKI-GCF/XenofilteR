# XenofilteR
XenofilteR: computational deconvolution of mouse and human reads in tumor xenograft sequence data

Human tumour samples or cancer cell lines, transplanted into mice are widely used as a 
model to study cancer. However, genomic analysis of tumour material derived from these 
xenografts is challenging. The sequence data not only contains reads that originate from 
the graft (human tumour or cell line) but also reads from host (mouse). We developed 
XenofilteR, an R-package for filtering host from graft reads in next generation sequence 
data. XenofilteR is a novel method that utilizes the edit distance of each read for 
classification. XenofilteR outperforms existing methods as validated on artificially 
mixed mouse/human samples and sets of patient derived xenograft samples. 

The paper that accompanies XenofilteR has been published in BMC Bioinformatics:
[XenoFilteR paper] (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2353-5)

Kluin, R. J. C. et al. XenofilteR: computational deconvolution of mouse and human reads in tumor xenograft sequence data. BMC Bioinformatics 19, 366 (2018).


## Installation

Installation of XenofilteR should be performed as
follows:

    > source("http://bioconductor.org/biocLite.R")
    > biocLite(c("Rsamtools", "GenomicAlignments", "BiocParallel", "futile.logger"))

As the last step in the installation process, the latest XenofilteR package can
be downloaded from the
[XenoFilteR releases webpage](https://github.com/PeeperLab/XenoFilteR/releases)
and installed using the following command:

    $ R CMD INSTALL XenofilteR*.tar.gz

Now you are all set to start your analysis.

## XenofilteR usage:


Load the XenofilteR package in R using:

    > library("XenofilteR")

XenofilteR fully supports paralel computing and is implemented in such a way
that every sample is processed on a single core. XenofilteR uses the
BiocParallel for parallel computing. This package requires the user to
specify which parallel environment is used by creating an instance of
BiocParallelParam. Here, we use the 'snow' environment and specify a SnowParam. 
The number of workers represents the number of CPUs used for the analysis and thereby 
the number of samples analyses simultaneously. Hence, it is wise to keep the number of 
CPUs low when analysing large samples as memory usage might be high. 

	> bp.param <- SnowParam(workers = 1, type = "SOCK")

XenofilteR requires a dataframe or matrix, named 'sample.list', with in the first 
column the bam file names as mapped to the graft reference. The second column contains the 
file names and paths to the bam files as mapped to the host reference. Each row in 
'sample.list' represents a single sequence run or sample. An optional list may be provided with 
alternative names for output files, 'output.names	'. Especially for RNAseq samples aligned with for example 
Tophat this may be convenient since all .bam files are named identical. 
The XenofilteR package and data are run in the following way: 

	> XenofilteR(sample.list, destination.folder = "./", bp.param = bp.param, output.names)


## Important note for STAR users

XenofilteR runs without issues on RNAseq mapped with STAR. However, the 'NM'-tag in the BAM file
is not generated by STAR with default settings. This tag is mandatory for XenofilteR to run. 
Please add the following to your STAR command when mapping to the human as well as the mouse reference genome:

--outSAMattributes NM

Additional explanation about the 'NM'-tag and optional output from STAR can be found here: 
STAR manual (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
(section 9.11)

## MM-threshold and Unmapped_penalty

For each sequence read an edit-distance is calculated which is used to classify a read as 
human or of mouse origin. Sequence reads with an edit distance higher than the MM-threshold 
(default == 4) will not be retained in the filtered output. This to remove sequence reads 
of mouse origin that do not map onto the mouse reference genome but do map onto the human 
reference genome. This problem is likely caused by gaps in the mouse reference genome and/or 
differences between mouse strains. The default of 3 mismatches has been tested for PE and 
SE sequence data of 75 basepairs. When applying XenofilteR to samples with longer 
sequence reads raising the default value should be considered. For example, setting the MM-threshold 
to 8 for samples with PE150. 

When analysing paired-end sequence data it happens that one of the reads of a pair is mapped 
to the reference genome while the other is not. The sequence read that is not mapped to the 
reference genome will be given an edit-distance as defined in 'Unmapped_penalty' (default == 8). 

## Independent performance assessment of XenofilteR 

The performance of XenofilteR was recently assessed by an independent group:

- Jo SY, Kim E, Kim S. Impact of mouse contamination in genomic profiling of 
	patient-derived models and best practice for robust analysis. Genome Biol. 2019 Nov 
	11;20(1):231. doi: 10.1186/s13059-019-1849-2. 
[[Pubmed]](http://www.ncbi.nlm.nih.gov/pubmed/31707992)    

## Contact
## 
We have tried to make the XenofilteR code readable and its use as easy
as possible. If any questions arise regarding the package, or if you
want to report any bugs, please do not hesitate and contact:

- [Roel Kluin](mailto:r.kluin@nki.nl)

Oscar and Roel both work at the NKI, Oscar as a postdoc in the lab of Prof. Daniel
Peeper, Roel as bioinformatician in the Genome Core Facility.


## Changes and additions we are currently working on
## 
- [ ] Catch STAR mapped samples missing the 'NM'-tag
- [ ] Optional output bam file with mouse reads
- [ ] Improve the vignette	
- [x] Add read-count statistics to the log file
- [x] Change output name and structure for RNAseq data
- [x] Making XenofilteR into an R-package 
- [x] Initial loop running the basic function of XenofilteR

