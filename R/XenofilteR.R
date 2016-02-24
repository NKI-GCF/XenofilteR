XenofilteR <- function(sample.list, destination.folder, bp.param, output.names = NULL) {

    ##########################
    ## Check and initialise ##
    ##########################

    start.time <- Sys.time()

    ## Restore work directory upon exit
    wd.orig <- getwd()
    on.exit(setwd(wd.orig))
  
    ## Make folder paths absolute
    sample.list <- apply(sample.list, c(1, 2),
                         tools::file_path_as_absolute)
    sample.list <- data.frame(sample.list, stringsAsFactors = FALSE)
    colnames(sample.list) <- c("Graft", "Host")
    destination.folder <- tools::file_path_as_absolute(destination.folder)

    ## Create lists with graft bam files, paths and names
    sample.paths.graft <- unlist(sample.list[, 1])
    sample.paths.graft <- unique(sample.paths.graft[!is.na(sample.paths.graft)])
    sample.files.graft <- basename(sample.paths.graft)

    ## Create lists with graft bam files, paths and names
    sample.paths.host <- unlist(sample.list[, 2])
    sample.paths.host <- unique(sample.paths.host[!is.na(sample.paths.host)])
    sample.files.host <- basename(sample.paths.host)
    
    ## Check length output.names and if unique
    if (length(output.names) != 0){

				## All alternative output names should end with '.bam'
				output.names <- gsub(".bam", "", output.names)
				output.names <- paste0(output.names, ".bam")
		
				if (length(output.names) != nrow(sample.list["Graft"])){
						stop(.wrap("The number of provided names does not match the number of samples.", 
								       "Please correct the file names in:", sQuote(output.names)))
				}
			
				if (length(output.names) != length(unique(output.names))){
						stop(.wrap("Identical samples names are used for multiple samples", 
					          	 "Please correct the file names in:", sQuote(output.names)))
				}
    }

    ## Check whether destination folder exists
    if (!file.exists(destination.folder)) {
        stop(.wrap("The destination folder could not be found. Please change",
                   "the path specified in", sQuote(destination.folder)))
    }

    ## Check for write permissions in the destination folder
    if (file.access(destination.folder, 2) == -1) {
        stop(.wrap("You do not have write permission in the destination",
                   "folder."))
    }    
    
    ## Create folder
    destination.folder <- file.path(destination.folder, "Filtered_bams")
    tryCatch({
        if (!file.exists(file.path(destination.folder))) {
            dir.create(file.path(destination.folder),
                       recursive = TRUE)
        } else {
            stop(.wrap("The folder",
                       sQuote(file.path(destination.folder, "Filtered_bams")),
                       "already exists. Please remove it, or (in case you",
                       "still need it), rename it to prevent files from being",
                       "overwritten."))
        }
    }, warning = function(e) {
        stop(.wrap("You do not have write permissions in the destination",
                   "folder. Stopping execution of the remaining part of the",
                   "script..."))
    })

    ## Calculate the maximal number of CPUs to be used
    ncpu <- bpworkers(bp.param)

    ## Provide output to log
    flog.appender(appender.file(file.path(destination.folder,
                                          "XenofilteR.log")))
    flog.info(paste("Running XenofilteR version",
                    as(packageVersion("XenofilteR"), "character"), "..."))
    flog.info(paste0("XenofilteR was run using the following commands:", "\n\n",
                     "XenofilteR(sample.list = sample.list, ",
                     "destination.folder = \"", dirname(destination.folder),
                     "\", BPPARAM = bp.param", ")"))
    flog.info("The value of bp.param was:", getClass(bp.param), capture = TRUE)
    flog.info(paste("This analysis will be run on", ncpu, "cpus"))
    flog.info("The value of sample.list was:", sample.list, capture = TRUE)
    if (length(output.names)!=0){
				for (i in seq_along(output.names)) {
						flog.info(paste0("Alternative sample name for ",
						sample.files.graft[i],":", output.names[i]))
				}
    }

    cat(.wrap("The following samples will be analyzed:"), "\n")
    cat(paste("graft:", sample.list[,1], ";", "\t", "matching",
               "host:", sample.list[,2]), sep = "\n")
    cat(.wrap("This analysis will be run on", ncpu, "cpus"), "\n")


    ## Check if graft bam files are indexed (.bai needed for filter step)
    chr.sort.mode <- NULL

    tryCatch({
        for (samp in sample.paths.graft) {
            header <- scanBamHeader(samp)
            chr.sort.mode <- c(chr.sort.mode, list(header[[1]]$text$'@HD'))
        }
    }, error = function(e) {
        stop(.wrap("The BAM file header of file", sQuote(samp), "is corrupted",
                   "or truncated. Please rebuild this BAM file or exclude it",
                   "from analysis. Stopping execution of the remaining part of",
                   "the script..."))
    })
    chr.sort.mode <- unlist(lapply(chr.sort.mode, function(x) {
				length(grep("SO:coordinate", x))
    }))

    if (any(chr.sort.mode == 0)) {
        stop(.wrap("The following .bam files are unsorted:"), "\n",
             paste(sample.paths.graft[which(chr.sort.mode == 0)],
             collapse = "\n"), "\n",
             "Please sort these .bam files based on coordinates")
    }

    ## Index graft .bam files
    if (!all(file.exists(gsub("$", ".bai", sample.paths.graft)))) {
        if (file.access(".", 2) == -1) {
            stop(.wrap("The .bam files are not indexed and you do not have",
                       "write permission in (one of) the folder(s) where the",
                       ".bam files are located."))
        }
        IndexBam <- function(sample.paths.graft) {
            indexBam(sample.paths.graft)
            paste0("indexBam(\"", sample.paths.graft, "\")")
        }
        to.log <- bplapply(sample.paths.graft, IndexBam, BPPARAM = bp.param)
        lapply(to.log, flog.info)
    }

    ## Check whether BAMs are paired-end
    NumberPairedEndReads <- function(sample.paths.graft) {
    
        bam <- open(BamFile(sample.paths.graft, yieldSize = 1))
        close(bam)
        what <- c("flag")
        param <- ScanBamParam(what = what)
        bam <- readGAlignments(bam, param = param)
        intToBits(mcols(bam)$flag)[1] == 01
    }
    is.paired.end <- bplapply(sample.paths.graft, NumberPairedEndReads,
                              BPPARAM = bp.param)
    is.paired.end <- unlist(is.paired.end)
    for (i in seq_along(sample.paths.graft)) {
        flog.info(paste0("Paired-end sequencing for sample ", sample.paths.graft[i],
                         ": ", is.paired.end[i]))
    }


    ##############################################
    ## Assigning reads to either mouse or human ##
    ##############################################

    i <- c(seq_along(sample.paths.graft))
    ActualFilter <- function(i, destination.folder, sample.list, is.paired.end, 
                             sample.paths.graft, sample.paths.host, bp.param){

				## Settings for scanBam
				p4 <- ScanBamParam(tag=c("NM"), what=c("qname", "flag", "cigar"), 
				flag=scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE))

				## Read human data (mapped and primary alignment only)
				Human <- scanBam(paste(sample.paths.graft[i]), param=p4)
		
				## Read mouse data (mapped and primary alignment only)
				Mouse <- scanBam(paste(sample.paths.host[i]), param=p4)


				# Get human reads that also map to mouse (TRUE if reads also maps to mouse)
				set <- Human[[1]]$qname %in% Mouse[[1]]$qname

				## Check if overlap exists
				if (sum(set)==0){
						stop(.wrap("No reads names overlap between graft and host BAM.",
						           "Either nothing maps to the host reference or the BAM files do not match.",
						           " Please double check the bam files."))
				}
		

				## Get read names mapped to human reference only
				HumanSet <- unique(Human[[1]]$qname[set==FALSE])
				

				#######################
				## The actual filter ##
				#######################

				filt <- list(setStart=function(x) x$qname %in% HumanSet)
		
				if (length(output.names)==0){
						filterBam(paste(sample.paths.graft[i]), file.path(destination.folder,
							gsub(".bam","_Strict_Filtered.bam",sample.files.graft[i])), filter=FilterRules(filt))
				}
				if (length(output.names)!=0){
						filterBam(paste(sample.paths.graft[i]), file.path(destination.folder,
							gsub(".bam", "_Strict_Filtered.bam",output.names[i])), filter=FilterRules(filt))
				}

    }
   
     #############
     ## Wrap-up ##
     #############
  
    to.log <- bplapply(i, ActualFilter, destination.folder, sample.list, BPPARAM = bp.param, 
                       is.paired.end, sample.paths.graft, sample.paths.host)
    flog.appender(appender.file(file.path(destination.folder,"XenofilteR.log")))
    #lapply(to.log, flog.info)

    ## Report calculation time to log file
   flog.info(paste("Total calculation time of XenofilteR was",
                    round(difftime(Sys.time(), start.time, units = "hours"), 2),
                    "hours"))
    cat("Total calculation time of XenofilteR was: ",
        round(difftime(Sys.time(), start.time, units = "hours"), 2), "\n\n")
        
    flog.info(paste(sessionInfo()))

}

