closest_join = function(left, right, by = 'chr') {
  left_join(left, right, by = by) %>%
    mutate(d = abs(start.y - start.x)) %>%
    group_by(chr, start.x, end.x) %>%
    filter(is.na(d) | d == min(d, na.rm = TRUE)) %>%
    ungroup %>%
    rename(start = start.x, end = end.x) %>%
    select(-start.y, -end.y) %>%
    return
}

