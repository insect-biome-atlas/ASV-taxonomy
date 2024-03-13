# Assign the most specific rank with a confidence (bootstrap) of at least 0.5.
# Adapted from https://github.com/nf-core/ampliseq/blob/master/modules/local/dada2_taxonomy.nf
library(dada2)
seq<-getSequences(snakemake@input$qry)
taxLevels <- snakemake@params$taxLevels
ref<- snakemake@input$ref
threads <- snakemake@threads
set.seed(snakemake@params$seed)
minBoot <- snakemake@params$minBoot

taxa <- assignTaxonomy(seq,ref, minBoot = minBoot, taxLevels = taxLevels, multithread = threads, verbose=TRUE, outputBootstraps = TRUE)
tx <- data.frame(ASV_ID = names(seq), taxa, row.names = names(seq))
colnames(tx) <- gsub("tax.","", colnames(tx))
# (2) Set confidence to the bootstrap for the most specific taxon
# extract columns with taxonomic values
#tax <- tx[ , grepl( "tax." , names( tx ) ) ]
colnames(tx) <- gsub("tax.","", colnames(tx))

# find first occurrence of NA
#res <- max.col(is.na(tax), ties = "first")
# correct if no NA is present in column to NA
#if(any(res == 1)) is.na(res) <- (res == 1) & !is.na(tax[[1]])
# find index of last entry before NA, which is the bootstrap value
#res <- res-1
# if NA choose last entry
#res[is.na(res)] <- ncol(tax)
# extract bootstrap values
#boot <- tx[ , grepl( "boot." , names( tx ) ) ]
#boot$last_tax <- res
#valid_boot <- apply(boot,1,function(x) x[x[length(x)]][1]/100 )
# replace missing bootstrap values (NA) with 0
#valid_boot[is.na(valid_boot)] <- 0
# add bootstrap values to column confidence
#tx$confidence <- valid_boot

# (3) Reorder columns before writing to file
expected_order <- c("ASV_ID",taxLevels)
expected_order <- intersect(expected_order,colnames(tx))
taxa_export <- subset(tx, select = expected_order)
#colnames(taxa_export) <- sub("tax.", "", colnames(taxa_export))
#rownames(taxa_export) <- names(seq)

write.table(taxa_export, file = snakemake@output[[1]], sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')