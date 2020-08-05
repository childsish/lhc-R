read_featureCounts = function(filename, n_samples) {
  read_tsv(filename, comment = '#', progress = FALSE, col_types = paste(c('cccccn', array('n', n_samples)), collapse = '')) %>%
    rename_at(vars(7:ncol(.)), paths_to_id) %>%
    rename(gene_id = Geneid,
             chr = Chr,
             start = Start,
             end = End,
             strand = Strand,
             length = Length) %>%
    mutate(chr = sapply(strsplit(chr, ';'), unique),
           start = sapply(lapply(strsplit(start, ';'), as.numeric), min),
           end = sapply(lapply(strsplit(end, ';'), as.numeric), max),
           strand = sapply(strsplit(strand, ';'), unique)) %>%
    return
}


