#!/usr/bin/Rscript

# Take a VCF file
# Extract REF-ALT counts --> frequencies
library(reshape2)
library(ggplot2)
setwd('~/mount2/scratch/jackson/ia327/first_exp/')
all.dirs <- dir(pattern = 'Del*')

for (d in all.dirs) {
  setwd(d)
  system('mkdir -p ploidy_estimate')
  all.files <- dir('calling/', pattern = 'vcf.gz')
  
  for (infile in all.files){
    ersno <- strsplit(infile, '\\.')[[1]][1]
    ncfile <- read.table('../name conversion.tsv', skip = 1 ,sep = '\t')
    sdno <- as.character(ncfile[ncfile[,6] == ersno, 5])
    
    ##### DATA INPUT and FORMATTING
    
    if (grepl('.gz', infile)){
      system(paste0('gunzip -c calling/', infile, ' > tmp.vcf'))
      datain <- read.table('tmp.vcf')
      system('rm tmp.vcf')
    } else {
      datain <- read.table(infile)
    }
    
    names(datain) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                       'INFO', 'FORMAT', 'VALUE')
    
    # Only SNPs and not from mitochondria
    datain <- datain[grepl('TYPE=snp;', datain$INFO),]
    datain <- datain[datain$CHROM != 'Mito',]
    
    values <- lapply(datain$VALUE, function (x) strsplit(as.character(x), ':')[[1]])
    values <- do.call("rbind", values)
    values <- data.frame(values)
    colnames(values) <- strsplit(as.character(datain$FORMAT[1]), ':')[[1]]
    
    alleles <- lapply(values$AD, function (x) {
      res <- strsplit(as.character(x), ',')[[1]]
      while (length(res) < 4){
        res <- c(res, NA)
      }
      return(res)
    })
    alleles <- do.call("rbind", alleles)
    alleles <- matrix(as.numeric(alleles), ncol = 4)
    colnames(alleles) <- c(paste0('A', 0:3))
    
    values <- cbind(values, alleles)
    values[,c('AD','RO','AO','QR','QA')] <- NULL
    values[,2:7] <- sapply(values[,2:7], function(x) as.numeric(as.character(x)))
    
    datain <- cbind(datain, values)
    datain$FORMAT <- NULL
    datain$VALUE <- NULL
    
    ##### FILTERING
    # Minimum coverage
    datain <- datain[datain$DP >= 10,]
    # Quality
    datain <- datain[datain$QUAL >= 10,]
    
    
    
    #### Fractions
    total.count <- rowSums(datain[,12:15], na.rm = TRUE)
    datain$A0 <- datain$A0 / total.count
    datain$A1 <- datain$A1 / total.count
    datain$A2 <- datain$A2 / total.count
    datain$A3 <- datain$A3 / total.count
    
    toplot <- list(A0 = datain$A0, A1 = datain$A1, A2 = datain$A2, A3 = datain$A3)
    
    pdf(paste0('ploidy_estimate/', ersno, '.pdf'))
    histo <- ggplot(melt(toplot), aes(value, fill = L1)) + 
      ggtitle(paste0(ersno, ' - ', sdno)) +
      scale_x_continuous(limits = c(0.05, 0.95)) +
      geom_histogram(aes(y=..density..),position = "stack", binwidth=0.015) +
      geom_density(stat = 'density', position = 'stack', alpha = 0.3) + 
      theme_minimal() + 
      labs(x = 'Allele frequency', y = 'Density') +
      geom_vline(xintercept = c(1/2, 2/3, 1/3, 3/4, 1/4, 4/5, 1/5, 5/6, 1/6),
                 linetype = 3)
    print(histo)
    dev.off()
  }
  setwd('..')
}
