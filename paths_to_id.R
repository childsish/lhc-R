get_common_prefix = function(paths) {
  max_length = max(vapply(paths, nchar, 0))
  for (i in 1:max_length) {
    if (length(table(substr(paths, i, i))) > 1)
      return(substr(paths[1], 1, i - 1));
  }
  return(paths[1])
}

get_common_suffix = function(paths) {
  max_length = max(vapply(paths, nchar, 0))
  for (i in 0:max_length) {
    if (length(table(substr(paths, nchar(paths) - i, nchar(paths)))) > 1)
      return(substr(paths[1], nchar(paths[1]) - i + 1, nchar(paths[1])));
  }
  return(paths[1])
}

paths_to_id = function(paths) {
  prefix = get_common_prefix(paths)
  suffix = get_common_suffix(paths)
  return(substr(paths, nchar(prefix) + 1, nchar(paths) - nchar(suffix)))
}
