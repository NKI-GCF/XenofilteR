\name{XenofilteR}
\alias{XenofilteR}
\alias{list.files}
\title{
XenofilteR: excluding host reads from xenograft sequencing data
}
\description{
Excludes DNA or RNA sequence reads that originate from contaminating species
(like the host in xenograft experiments. 
}
\usage{
XenofilteR(sample.list, destination.folder, bp.param, output.names, MM_threshold)
}
\arguments{
\item{sample.list}{
A data.frame or matrix that contains the full name of the bam files (including
path) of the reads aligned to the human reference (column 1) and the mouse
reference (column 2).
}
\item{destination.folder}{
The path to the folder to which output should be written. The path can be
either absolute or relative.
}
\item{bp.param}{
A BiocParallelParam instance (see \code{BiocParallel} Bioconductor package)
that determines the settings used for parallel computing. Please refer to the
vignette for more information.
}
\item{output.names}{
Optional list with alternative samples names for the generated filtered bam
file. This might become useful when analysing RNAseq data. Bam files are often
name identical and sample names are defined in the folder name, not in he bam
file.
}
\item{MM_threshold}{
The edit distance of reads that map to the human reference only should be below 
this threshold. Default = 5.
}
\item{Unmapped_penalty}{
Only applicable for paired-end sequence data. If one of the reads of a pair is
not mapped the edit distance for that read is set to this value.  Default = 4.
}
}
\details{
XenofilteR reads four parameters from the human and mouse bam file: 1) read
name, 2) mapping quality, 3) CIGAR, and 4) mismatches (NM tag). Based
on these parameters XenofilteR calculates two edit distances for each read, one
for the mapping to the human (graft) reference and a second for the mapping to
the mouse (host) reference. The read or read-pair is assigned to the host or
graft .bam files based on these edit distances.
}
\value{
\item{BamFiles}{
A folder with the .bam files from which reads mapping to xenogemones have been
removed.
}
\item{XenofilteR.log}{
Log file of XenofilteR. This file also contains the statistics for each sample
processed. 
}
}
\references{
XenofilteR: seperate graft from host reads in xenograft sequencing. Roel Kluin,
Oscar Krijgsman, et al. To be submitted for publication. 
}
\author{
Oscar Krijgsman (o.krijgsman@nki.nl)
Roel Kluin (r.kluin@nki.nl)
}
\examples{
\dontrun{
# 
library(XenofilteR)
path <- dir(system.file(package = "XenofilteR"), "extdata", full.names = TRUE)
samples <- list.files(path = path, pattern = ".bam$", full.names = TRUE)

sample.list <- data.frame(samples[1], samples[2])
bp.param <- SnowParam(workers = 1, type = "SOCK")

XenofilteR(sample.list = sample.list, destination.folder =
        "/PATH/TO/DESTINATIONFOLDER/", bp.param = bp.param, output.names = NULL, MM_threshold=5)
}
}
