\name{XenoFilteR}
\alias{XenoFilteR}
\alias{list.files}
\title{
XenoFilteR: filter reads from xenograft sequencing
}
\description{
Xenograft filters DNA and RNA sequence data for reads that originate from the host in xenograft experiments. 
}
\usage{
XenoFilteR(sample.list, destination.folder, bp.param)
}
\arguments{
\item{sample.list}{
a data.frame or matrix that contains the path and name of the bam files of the reads alligned to the human reference (column 1) and the mouse reference (column 2).
}
\item{destination.folder}{
The path to the folder to which output should be written. The path can be
either absolute or relative.
}
\item{bp.param}{
a BiocParallelParam instance (see \code{BiocParallel} Bioconductor package)
that determines the settings used for parallel computing. Please refer to the
vignette for more information.
}
}
\details{
XenofilteR reads five parameters from the human and mouse bam file: 1) read name, 2) mapping quality 3) flag, 4) CIGAR and 5) mismatches (NM tag). Based on these parameters XenofilteR calculates two scores for each read, one for teh mapping against the human (graft) reference and a second for the mapping against teh mouse (host) reference. The read, or read-pair when performing paired-end experiments, is assigned to either derived from teh host or graft based on the score. 

}
\value{
\item{BamFiles}{
a folder with the .bam files filtered for mouse derived reads that are created during
the \code{XenofilteR}() run.
}
\item{XenofilteR.log}{
log file of XenofilteR. This file also contains the statistics for each sample processed. 
}
}
\references{
XenofilteR: DNA copy number detection from off-target sequence data. Roel Kluin, Oscar Krijgsman, et al. To be submitted for
publication. 
}
\author{
Oscar Krijgsman (o.krijgsman@nki.nl)
}
\examples{
\dontrun{
# 
setwd("/PATH/TO/BAMFILES/")
human_bams <- list.files(pattern = "hg19.bam$", full.names = TRUE)
mouse_bams <- list.files(pattern = "mm10.bam$", full.names = TRUE)

sample.list <- cbind(human_bams,mouse_bams)

CopywriteR(sample.list = sample.list, destination.folder =
        "/PATH/TO/DESTINATIONFOLDER/", ncpu = 1)
}
}
