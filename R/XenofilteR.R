XenofilteR<-function(Sample_list, destination.folder, bp.param){

	##########################
    ## Check and initialise ##
    ##########################

    start.time <- Sys.time()

    ## Restore work directory upon exit
    wd.orig <- getwd()
    on.exit(setwd(wd.orig))

	## Make folder paths absolute
    Sample_list <- apply(Sample_list, c(1, 2),
                            tools::file_path_as_absolute)
    Sample_list <- data.frame(Sample_list, stringsAsFactors = FALSE)
    colnames(Sample_list) <- c("Graft", "Host")
    destination.folder <- tools::file_path_as_absolute(destination.folder)

    ## Create lists with graft bam files, paths and names
    sample.paths.graft <- unlist(Sample_list[,1])
    sample.paths.graft <- unique(sample.paths.graft[!is.na(sample.paths.graft)])
    sample.files.graft <- basename(sample.paths.graft)

    ## Create lists with all bam files and paths
    sample.paths <- unlist(Sample_list)
    sample.paths <- unique(sample.paths[!is.na(sample.paths)])


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
                     "XenofilteR(Sample_list = Sample_list, ",
                     "destination.folder = \"", dirname(destination.folder),
                     "\", BPPARAM = bp.param", ")"))
    flog.info("The value of bp.param was:", getClass(bp.param), capture = TRUE)
    flog.info("The value of Sample_list was:", Sample_list,
              capture = TRUE)
    flog.info(paste("This analysis will be run on", ncpu, "cpus"))

	## Check if graft bam files are sorted and index if not (.bai needed for filter step)
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
    NumberPairedEndReads <- function(sample.paths) {
    
        bam <- open(BamFile(sample.paths, yieldSize = 1))
        close(bam)
        what <- c("flag")
        param <- ScanBamParam(what = what)
        bam <- readGAlignments(bam, param = param)
        intToBits(mcols(bam)$flag)[1] == 01
    }
    is.paired.end <- bplapply(sample.paths, NumberPairedEndReads,
                              BPPARAM = bp.param)
    is.paired.end <- unlist(is.paired.end)
    for (i in seq_along(sample.paths)) {
        flog.info(paste0("Paired-end sequencing for sample ", sample.paths[i],
                         ": ", is.paired.end[i]))
    }


	###################
    ## Actual filter ##
    ###################

    i <- c(seq_along(sample.paths))
	ActualFilter<-function(i, destination.folder, Sample_list, is.paired.end){

		## Read human data (all reads)
		p4 <- ScanBamParam(tag=c("NM"), what=c("qname", "mapq", "flag", "cigar"), flag=scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE)))
		Human <- scanBam(paste(sample.paths[i]), param=p4)
		cat("Finished reading human sample", Sample_list[i,1], "\n")
		
		## Read Mouse data (mapped only)
		p5 <- ScanBamParam(tag=c("NM"), what=c("qname", "mapq", "flag", "cigar"), flag=scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE))
		Mouse <- scanBam(paste(Sample_list[i,2]), param=p5)
		cat("Finished reading mouse sample", Sample_list[i,2], "\n")

		# Get human reads that also map to mouse (TRUE if reads also maps to mouse)
		set<-Human[[1]]$qname%in%Mouse[[1]]$qname

		## Check if overlap exists
		if (sum(set)==0){
			stop(.wrap("No reads names overlap between graft and host BAM. Either nothing maps to the host reference or the BAM files do not match. Execution stopped for this sample"))
		}

		ToHumanOnly<-unique(Human[[1]]$qname[set==FALSE])

		## Get the Clips + inserts + MisMatches (Mouse)
		Cigar.matrix<-cigarOpTable(Mouse[[1]]$cigar)
		Inserts<-Cigar.matrix[,colnames(Cigar.matrix)=="I"]
		Clips<-Cigar.matrix[,colnames(Cigar.matrix)=="S"]
		MM_I_mouse<-Clips+Inserts+Mouse[[1]]$tag$NM

		## Get the Clips + inserts + MisMatches (Human)
		Cigar.matrix<-cigarOpTable(Human[[1]]$cigar)
		Inserts<-Cigar.matrix[,colnames(Cigar.matrix)=="I"]
		Clips<-Cigar.matrix[,colnames(Cigar.matrix)=="S"]
		MM_I_human<-Clips+Inserts+Human[[1]]$tag$NM


		## Check per read-pair where read-pair belongs (mouse or human)

		## Filter for human reads that also map to mouse
		Human_qname_set<-Human[[1]]$qname[set==TRUE]
		Human_mapq_set<-Human[[1]]$mapq[set==TRUE]
		MM_I_human_set<-MM_I_human[set==TRUE]


		uni.name<-unique(Human_qname_set)
		Map_info<-matrix(data=0, ncol=8, nrow=length(uni.name))
		row.names(Map_info)<-uni.name
		colnames(Map_info)<-c("MM_mouse_F","MM_mouse_R","MM_human_F","MM_R_human", "Mq_mouse_F", "Mq_mouse_R", "Mq_human_F", "Mq_human_R")
	
	
		### Extensive data table with the number of mismatches (+ clips) and mapping quality.
		### This table is usefull for de-buging and checking data, But should be remove in 
		### final version of XenofilteR
 
		## Match for forward and reverse reads and get MM+I (mouse)
		FR_mouse<-lapply(Mouse[[1]]$flag, .FirstInPair)
		RR_mouse<-lapply(Mouse[[1]]$flag, .SecondInPair)

		## Fill dataframe with mismatches and mappin quality for mouse
		Map_info[,1]<-MM_I_mouse[unlist(FR_mouse)][match(uni.name, Mouse[[1]]$qname[unlist(FR_mouse)])]
		Map_info[,2]<-MM_I_mouse[unlist(RR_mouse)][match(uni.name, Mouse[[1]]$qname[unlist(RR_mouse)])]

		Map_info[,5]<-Mouse[[1]]$mapq[unlist(FR_mouse)][match(uni.name, Mouse[[1]]$qname[unlist(FR_mouse)])]
		Map_info[,6]<-Mouse[[1]]$mapq[unlist(RR_mouse)][match(uni.name, Mouse[[1]]$qname[unlist(RR_mouse)])]

		## Match for forward and reverse reads and get MM+I (human)
		FR_human<-lapply(Human[[1]]$flag[set==TRUE], .FirstInPair)
		RR_human<-lapply(Human[[1]]$flag[set==TRUE], .SecondInPair)
	
		## Fill dataframe with mismatches and mappin quality for human
		Map_info[,3]<-MM_I_human_set[unlist(FR_human)][match(uni.name, Human_qname_set[unlist(FR_human)])]
		Map_info[,4]<-MM_I_human_set[unlist(RR_mouse)][match(uni.name, Human_qname_set[unlist(RR_mouse)])]

		Map_info[,7]<-Human_mapq_set[unlist(FR_human)][match(uni.name, Human_qname_set[unlist(FR_human)])]
		Map_info[,8]<-Human_mapq_set[unlist(RR_mouse)][match(uni.name, Human_qname_set[unlist(RR_mouse)])]

		## Calculate average score for human and mouse reads	
		Score_mouse<-rowMeans(cbind(Map_info[,"MM_mouse_F"]/Map_info[,"Mq_mouse_F"], Map_info[,"MM_mouse_R"]/Map_info[,"Mq_mouse_R"]), na.rm=T)
		Score_human<-rowMeans(cbind(Map_info[,"MM_human_F"]/Map_info[,"Mq_human_F"], Map_info[,"MM_human_F"]/Map_info[,"Mq_human_F"]), na.rm=T)
	
		# Determine where reads fit better
		# Score human lower than mouse or no score for mouse at all (mapq==0)
		BetterToHuman<-row.names(Map_info)[which(Score_human<Score_mouse | (is.na(Score_mouse)==TRUE & is.na(Score_human)==FALSE))]

		Filt<-c(ToHumanOnly, BetterToHuman)

		cat("Finished calculating which reads can be assigned to human - Start writing filtered Bam files", "\n")

		###########################################################################
		############################ The actual filter ############################
		###########################################################################

		filt <- list(setStart=function(x) x$qname %in% Filt)
		filterBam(paste(Sample_list[i,1]), paste0(destination.folder,"/", gsub(".bam","_Filtered.bam",sample.files.graft[i])), 
			filter=FilterRules(filt))
		cat("Finished writing",gsub(".bam","_Filtered.bam",Sample_list[i,1]), " ---  sample", i, "out of", nrow(Sample_list), "\n")

	}
   
    to.log <- bplapply(i, destination.folder, ActualFilter, Sample_list, is.paired.end, BPPARAM = bp.param)
    lapply(to.log, flog.info)


	## Calculation time etc. 
 	flog.info(paste("Total calculation time of XenofilteR was",
                    round(difftime(Sys.time(), start.time, units = "hours"), 2),
                    "hours"))
    cat("Total calculation time of XenofilteR was: ",
        Sys.time() - start.time, "\n\n")

}

