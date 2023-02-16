resume = function(inner_function, path, verbose = FALSE) {
  library(digest)

  checksum_path = paste0(path, '.checksum')
  checksum = digest(inner_function, algo = 'crc32')
  if (file.exists(path) && file.exists(checksum_path) && readLines(checksum_path) == checksum) {
    if (verbose) {
      message(paste('resuming from', path))
    }
    return(readRDS(path))
  }
  if (verbose) {
    message(paste('persisting to', path))
  }
  writeLines(checksum, checksum_path)
  data = inner_function()
  saveRDS(data, path)
  return(data)
}
