sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
library(dada2)
seq<-getSequences(snakemake@input$qry)
taxLevels <- snakemake@params$taxLevels
ref<- snakemake@input$ref
threads <- snakemake@threads
minBoot <- as.integer(snakemake@wildcards$boot)

set.seed(snakemake@params$seed)
taxa <- assignTaxonomy(seq, ref, minBoot = minBoot, taxLevels = taxLevels, multithread = threads, tryRC = TRUE, verbose=TRUE, outputBootstraps = FALSE)
tx <- data.frame(ASV_ID = names(seq), taxa, row.names = names(seq))
#expected_order <- c("ASV_ID",taxLevels)
#expected_order <- intersect(expected_order,colnames(tx))
#taxa_export <- subset(tx, select = expected_order)
write.table(tx, file = snakemake@output[[1]], sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')
sink()