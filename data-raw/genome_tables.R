## code to prepare `DATASET` dataset goes here

cattle_genome_table <- readr::read_tsv("inst/extdata/cattle_genome_table.txt")

cattle_genome_table <- cattle_genome_table[, c("chr", "length", "genetic_length")]

usethis::use_data(cattle_genome_table, overwrite = TRUE)
