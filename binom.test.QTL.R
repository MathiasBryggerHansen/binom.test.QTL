getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

gffRead <- function(gffFile) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = -1,
                   colClasses=c("character", "character", "character", "integer","integer","character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end","score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  return(gff)
}

#' Binomial test of SNP QTL
#'
#' This function runs the binomial test for each phenotype in global QTL dataframe against
#' SNP QTL. It can be used to support QTLs found in a study, and find novel interactions between different traits and the trait used in the study.
#'
#' @param gff Path to the global QTL file including multiple phenotypes
#' @param s is boolean for subsetting of breeds
#' @param s_v is a vector with the breeds to be subsetted on
#' @param genome_length  is the length in base pairs of the species analysed
#' @param chr_col chromosome name column, must be in gff format (Chr.nr)
#' @param pos_col SNP position column
#' @param snp_col SNP id column
#' @param alternative binomial alternative test hypothesis
#' @param qtl_table Path to the SNP QTL file for a given phenotype, csv format
#' @return dataframe of frequency overlap and binomial probability
#' @export
binom.test.QTL <- function(gff,qtl_table,s,s_v,genome_length,chr_col,pos_col,snp_col,alternative = c("two.sided", "less", "greater")) {
  #loading data
  qtls <- read.csv(file = qtl_table,header = T)
  qtls <- qtls[order(qtls[,chr_col]),]
  #Creating genomic ranges for both gff and qtl table files
  gff2 <- gffRead(gff)
  gff2$breed <- getAttributeField(x = gff2$attributes,field = "breed")
  if(s){
    gff2 <- gff2[grepl(pattern = paste(s_v,collapse = "|"),gff2$breed,ignore.case = T),]
  }
  gff2$trait <- getAttributeField(x = gff2$attributes,field = "trait")
  gff2$qtl_id <- getAttributeField(x = gff2$attributes,field = "QTL_ID")
  gff2 <- gff2[order(gff2$seqname,gff2$start),]
  gff2 <- gff2[!is.na(gff2$start)&!is.na(gff2$end),]

  gr <- GenomicRanges::GRanges(Rle(unique(gff2$seqname), table(gff2$seqname)),IRanges::IRanges(gff2$start, width=gff2$end - gff2$start + 1, names=gff2$trait))
  singles <- GenomicRanges::split(gr, names(gr))
  pheno_id <- names(singles)
  g_qtl <- GenomicRanges::GRanges(Rle(unique(qtls[,chr_col]), table(qtls[,chr_col])),IRanges::IRanges(qtls[,pos_col], width=1, names=as.vector(qtls[,snp_col])))
  eqtl_nr <- length(g_qtl)
  db_qtls <- length(singles)

  #result dataframe
  res <- matrix(data=list(),nrow = db_qtls,ncol = 5)

  for(i in 1:db_qtls){ #reduce(pintersect(..))?
    print(paste("testing nr ",i," out of ",db_qtls," db QTLs"))
    phen <- unlist(singles[i])
    phen <- GenomicRanges::reduce(phen) #reduce merges the QTLs that overlap, see ?reduce
    temp_test <- GenomicRanges::subsetByOverlaps(g_qtl, phen)
    intersect_count <- length(temp_test)
    qtl_width <- sum(as.numeric(GenomicRanges::width(GenomicRanges::union(phen,g_qtl))))
    if(intersect_count>eqtl_nr){
      intersect_count <- eqtl_nr
    }
    exp <- qtl_width/genome_length
    temp <- binom.test(x = intersect_count,n = eqtl_nr,p = exp,alternative = alternative)

    snps <- list(names(temp_test))
    res[[i,1]] <- snps
    res[[i,2]] <- exp
    res[[i,3]] <- intersect_count/eqtl_nr
    res[[i,4]] <- temp$p.value
  }
  res[,5] <- pheno_id
  res <- data.frame(res)
  colnames(res) <- c("SNP_ID","expected_frequency","actual_frequency","probability","phenotype")
  return(res)
}
