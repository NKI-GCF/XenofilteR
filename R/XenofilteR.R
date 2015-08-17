XenofilteR<-function(Sample_list, destination.folder, bp.param){

	##########################
    ## Check and initialise ##
    ##########################

    start.time <- Sys.time()

    ## Restore work directory upon exit
    wd.orig <- getwd()
    on.exit(setwd(wd.orig))

	## Make folder paths absolute
    sample.control <- apply(Sample_list.control, c(1, 2),
                            tools::file_path_as_absolute)
    Sample_list <- data.frame(Sample_list, stringsAsFactors = FALSE)
    colnames(Sample_list) <- c("samples", "controls")
    destination.folder <- tools::file_path_as_absolute(destination.folder)

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

	## Provide output to log
    flog.appender(appender.file(file.path(destination.folder,
                                          "XenofilteR.log")))
    flog.info(paste("Running XenofilteR version",
                    as(packageVersion("XenofilteR"), "character"), "..."))
    flog.info(paste0("XenofilteR was run using the following commands:", "\n\n",
                     "XenofilteR(Sample_list = Sample_list, ",
                     "destination.folder = \"", dirname(destination.folder),
                     "\", BPPARAM = bp.param"))
    flog.info("The value of bp.param was:", getClass(bp.param), capture = TRUE)
    flog.info("The value of Sample_list was:", Sample_list,
              capture = TRUE)
    flog.info(paste("This analysis will be run on", ncpu, "cpus"))




	###################
    ## Actual filter ##
    ###################

	for (i in 1:nrow(Sample_list)){

		## Create list of .bam files
		flog.info("XenofilteR will analyze the following (unique) samples:",
				  Sample_list, capture = TRUE)


		## Read human data (all reads)
		p4 <- ScanBamParam(tag=c("NM"), what=c("qname", "mapq", "flag", "cigar"), flag=scanBamFlag(isUnmappedQuery=FALSE))
		Human <- scanBam(paste(Sample_list[i,1]), param=p4)
		cat("Finished reading human sample", Sample_list[i,1], "\n")
	
		# Filter Human data for 'Multi mappers'
	
		##### Still to do
	
		## Read Mouse data (mapped only)
		p5 <- ScanBamParam(tag=c("NM"), what=c("qname", "mapq", "flag", "cigar"), flag=scanBamFlag(isUnmappedQuery=FALSE))
		Mouse <- scanBam(paste(Sample_list[i,2]), param=p5)
		cat("Finished reading mouse sample", Sample_list[i,2], "\n")

		# Filter Mouse data for 'Multi mappers'

		##### Still to do

		# Get human reads that also map to mouse (TRUE if reads also maps to mouse)
		set<-Human[[1]]$qname%in%Mouse[[1]]$qname

		# Table with the classification of each read (based on human bam)
		Filter_table<-rep(0,length(Human[[1]]$qname))
	
		## If not mapped to mouse set as 1
		Filter_table[set==FALSE]<-1
	
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
		Filter_table_set<-Filter_table[set==TRUE]
		MM_I_human_set<-MM_I_human[set==TRUE]


		uni.name<-unique(Human_qname_set)
		Map_info<-matrix(data=0, ncol=8, nrow=length(uni.name))
		row.names(Map_info)<-uni.name
		colnames(Map_info)<-c("MM_mouse_F","MM_mouse_R","MM_human_F","MM_R_human", "Mq_mouse_F", "Mq_mouse_R", "Mq_human_F", "Mq_human_R")


		######################################################
		## Move functions to .private folder in R-package
		## Function to check which reads are first in pair 
		FirstInPair<-function(x){
			intToBits(x)[7]=="01"
		}
	
		## Function to check which reads are second in pair
		SecondInPair<-function(x){
			intToBits(x)[8]=="01"
		}
		######################################################
	
	
		### Extensive data table with the number of mismatches (+ clips) and mapping quality.
		### This table is usefull for de-buging and checking data, But should be remove in 
		### final version of XenofilteR
 
		## Match for forward and reverse reads and get MM+I (mouse)
		FR_mouse<-lapply(Mouse[[1]]$flag, FirstInPair)
		RR_mouse<-lapply(Mouse[[1]]$flag, SecondInPair)

		## Fill dataframe with mismatches and mappin quality for mouse
		Map_info[,1]<-MM_I_mouse[unlist(FR_mouse)][match(uni.name, Mouse[[1]]$qname[unlist(FR_mouse)])]
		Map_info[,2]<-MM_I_mouse[unlist(RR_mouse)][match(uni.name, Mouse[[1]]$qname[unlist(RR_mouse)])]

		Map_info[,5]<-Mouse[[1]]$mapq[unlist(FR_mouse)][match(uni.name, Mouse[[1]]$qname[unlist(FR_mouse)])]
		Map_info[,6]<-Mouse[[1]]$mapq[unlist(RR_mouse)][match(uni.name, Mouse[[1]]$qname[unlist(RR_mouse)])]

		## Match for forward and reverse reads and get MM+I (human)
		FR_human<-lapply(Human[[1]]$flag[set==TRUE], FirstInPair)
		RR_human<-lapply(Human[[1]]$flag[set==TRUE], SecondInPair)
	
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

		Filt<-c(unique(Human[[1]]$qname[set==FALSE]), BetterToHuman)

		cat("Finished calculating which reads can be assigned to human - Start writing filtered Bam files", "\n")

		###########################################################################
		############################ The actual filter ############################
		###########################################################################

		filt <- list(setStart=function(x) x$qname %in% Filt)
		filterBam(paste(Sample_list[i,1]), gsub(".bam","_Filtered.bam",Sample_list[i,1]), 
			filter=FilterRules(filt))
		cat("Finished writing",gsub(".bam","_Filtered.bam",Sample_list[i,1]), " ---  sample", i, "out of", nrow(Sample_list), "\n")

	}



	## Calculation time etc. 
 	flog.info(paste("Total calculation time of XenofilteR was",
                    round(difftime(Sys.time(), start.time, units = "hours"), 2),
                    "hours"))
    cat("Total calculation time of XenofilteR was: ",
        Sys.time() - start.time, "\n\n")

}

