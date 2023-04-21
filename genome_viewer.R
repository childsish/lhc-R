add_feature_heights = function(features) {
  stopifnot(is_grouped_df(features))
  
  parents = features %>%
    arrange(start, end) %>%
    summarise(
      start = min(start, end),
      end = max(start, end),
      .groups = 'drop'
    )
  
  ends = parents$end[1]
  parents$height = NA
  parents$height[1] = 1
  if (nrow(parents) > 1) {
    for (i in 2:nrow(parents)) {
      for (j in 1:length(ends)) {
        if (ends[j] <= parents$start[i]) {
          parents$height[i] = j
          break
        }
      }
      if (is.na(parents$height[i])) {
        parents$height[i] = length(ends) + 1
      }
      ends[parents$height[i]] = parents$end[i]
    }
  }
  
  features %>%
    left_join(parents %>% select(-start, -end), by = group_vars(features)) %>%
    return
}

scale_coordinates = function(coordinates, lower_bound, upper_bound) {
  stopifnot(lower_bound < upper_bound)
  return((coordinates - lower_bound) / (upper_bound - lower_bound))
}

scale_and_transform_coordinates = function(coordinates, lower_bound, upper_bound) {
  stopifnot(lower_bound < upper_bound)
  return(coordinates * (upper_bound - lower_bound) + lower_bound)
}

transform_coordinates = function(coordinates, lower_bound, upper_bound) {
  stopifnot(lower_bound < upper_bound)
  return(coordinates / (upper_bound - lower_bound) + lower_bound)
}

plot_features = function(features, title, bounds = c(0, 0, 1, 1), labels = NULL) {
  legend(min(bounds[1], bounds[3]), max(bounds[2], bounds[4]), title, box.col = '#ffffff', adj = c(0, 1))
  discard = features %>%
    group_map(~ plot_feature(.x, .y, bounds, labels))
}

plot_feature = function(features, groups, bounds = c(0, 0, 1, 1), labels = NULL) {
  if (nrow(features) == 0) {
    return()
  }
  
  x_lower = bounds[1]
  y_lower = bounds[2]
  x_upper = bounds[3]
  y_upper = bounds[4]
  
  scaled_features = features %>%
    mutate(
      x1 = scale_and_transform_coordinates(x1, x_lower, x_upper),
      y1 = scale_and_transform_coordinates(y1, y_lower, y_upper),
      x2 = scale_and_transform_coordinates(x2, x_lower, x_upper),
      y2 = scale_and_transform_coordinates(y2, y_lower, y_upper),
    )
  
  plot_box(
    scaled_features$x1[1],
    scaled_features$y1[1],
    scaled_features$x2[1],
    scaled_features$y2[1],
    scaled_features$feature[1]
  )
  if (nrow(features) > 1) {
    for (i in 2:nrow(scaled_features)) {
      plot_box(
        scaled_features$x1[i],
        scaled_features$y1[i],
        scaled_features$x2[i],
        scaled_features$y2[i],
        scaled_features$feature[i]
      )
      ymid = (scaled_features$y1[i] + scaled_features$y2[i]) / 2
      segments(
        scaled_features$x2[i - 1],
        ymid,
        scaled_features$x1[i],
        ymid
      )
    }
  }
  label = ifelse(is.null(labels), paste0(groups[1,], collapse = ', '), scaled_features[[labels]][1])
  text(
    max(scaled_features$x1[1], 0),
    scaled_features$y1[1],
    label,
    adj = c(0, 0),
    cex = 0.8
  )
}

plot_box = function(x1, y1, x2, y2, feature=NULL) {
  yadj = ifelse(!is.null(feature) & feature %in% c('CDS'), 0.9, 0.8) * abs(y2 - y1)
  polygon(
    c(x1, x1, x2, x2),
    c(y1 + yadj, y2 - yadj, y2 - yadj, y1 + yadj)
  )
}

add_heights = function(reads) {
  ends = reads$x2[1]
  reads$height = NA
  reads$height[1] = 1
  for (i in 2:nrow(reads)) {
    for (j in 1:length(ends)) {
      if (ends[j] <= reads$x1[i]) {
        reads$height[i] = j
        break
      }
    }
    if (is.na(reads$height[i])) {
      reads$height[i] = length(ends) + 1
    }
    ends[reads$height[i]] = reads$x2[i]
  }
  return(reads)
}

plot_reads = function(reads, bounds = c(0, 0, 1, 1)) {
  if (nrow(reads) == 0) {
    return()
  }
  
  x_lower = bounds[1]
  y_lower = bounds[2]
  x_upper = bounds[3]
  y_upper = bounds[4]
  
  scaled_reads = reads %>%
    mutate(
      x1_read1 = scale_and_transform_coordinates(x1_read1, x_lower, x_upper),
      x2_read1 = scale_and_transform_coordinates(x2_read1, x_lower, x_upper),
      x1_read2 = scale_and_transform_coordinates(x1_read2, x_lower, x_upper),
      x2_read2 = scale_and_transform_coordinates(x2_read2, x_lower, x_upper),
      y1 = scale_and_transform_coordinates(y1, y_lower, y_upper),
      y2 = scale_and_transform_coordinates(y2, y_lower, y_upper),
    )
  
  xdiff = bounds[3] - bounds[1]
  ydiff = bounds[4] - bounds[2]
  
  for (i in 1:nrow(reads)) {
    if (!is.na(scaled_reads$pos_read1[i])) {
      draw_read(
        scaled_reads$x1_read1[i],
        scaled_reads$y1[i],
        scaled_reads$x2_read1[i],
        scaled_reads$y2[i],
        reads$strand_read1[i],
        col = reads$colour[i]
      )
    }
    if (!is.na(scaled_reads$pos_read2[i])) {
      draw_read(
        scaled_reads$x1_read2[i],
        scaled_reads$y1[i],
        scaled_reads$x2_read2[i],
        scaled_reads$y2[i],
        reads$strand_read2[i],
        col = reads$colour[i]
      )
    }
    if (!is.na(scaled_reads$pos_read1[i]) && !is.na(scaled_reads$pos_read2[i])) {
      ymid = (scaled_reads$y1[i] + scaled_reads$y2[i]) / 2
      if (reads$x2_read1[i] < reads$x1_read2[i]) {
        segments(reads$x2_read1[i], ymid, reads$x1_read2[i], ymid)
      } else if (reads$x2_read2[i] < reads$x1_read1[i]) {
        segments(reads$x2_read2[i], ymid, reads$x1_read1[i], ymid)
      }
    }
  }
}

draw_read = function(x1, y1, x2, y2, direction, col) {
  ymid = (y1 + y2) / 2
  ydiff = abs(ymid - y1)
  if (direction == '+') {
    polygon(
      c(x1, x1, x2 - ydiff, x2, x2 - ydiff),
      c(y1, y2, y2, ymid, y1),
      col = col,
      border = NA
    )
  } else {
    polygon(
      c(x1, x1 + ydiff, x2, x2, x1 + ydiff),
      c(ymid, y1, y1, y2, y2),
      col = col,
      border = NA
    )
  }
}

get_read_coordinates = function(reads, lower, upper) {
  reads %>%
    `[[`(1) %>%
    `[`(c('qname', 'flag', 'rname', 'strand', 'pos', 'qwidth', 'mapq')) %>%
    as_tibble %>%
    filter(
      mapq > 20,
      pos > lower,
      pos + qwidth < upper
    ) %>%
    mutate(
      read = ifelse(bitwAnd(flag, 64) > 0, 'read1', 'read2'),
    ) %>%
    select(-flag, -mapq) %>%
    pivot_wider(names_from = read, values_from = c(strand, pos, qwidth)) %>%
    mutate(
      colour = ifelse(strand_read1 == '+', '#d6604d7f', '#4393c37f'),
      x1_read1 = scale_coordinates(pos_read1, lower, upper),
      x2_read1 = scale_coordinates(pos_read1 + qwidth_read1, lower, upper),
      x1_read2 = scale_coordinates(pos_read2, lower, upper),
      x2_read2 = scale_coordinates(pos_read2 + qwidth_read2, lower, upper),
      x1 = pmin(x1_read1, x1_read2, na.rm = TRUE),
      x2 = pmax(x2_read1, x2_read2, na.rm = TRUE),
    ) %>%
    add_heights %>%
    mutate(
      y1 = 1 - scale_coordinates(height - 1, 0, max(height)),
      y2 = 1 - scale_coordinates(height, 0, max(height)),
    ) %>%
    return
}

get_features_in_range = function(features, chromosome, lower, upper) {
  feature_coordinates = features %>%
    filter(
      chr == chromosome,
      lower < end,
      start < upper,
    )
  
  if (nrow(feature_coordinates) == 0) {
    return(feature_coordinates)
  }
  
  feature_coordinates %>%
    group_by(transcript_id) %>%
    add_feature_heights %>%
    ungroup %>%
    mutate(
      x1 = scale_coordinates(start, lower, upper),
      y1 = scale_coordinates(height - 1, 0, max(height)),
      x2 = scale_coordinates(end, lower, upper),
      y2 = scale_coordinates(height, 0, max(height)),
    ) %>%
    group_by(transcript_id) %>%
    return
}