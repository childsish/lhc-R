rotate_array = function(array, n=-1) {
  if (n < 0) {
    return(c(array[(length(array) + n + 1):length(array)], array[1:(length(array) + n)]))
  }
  else {
    return(c(array[(n + 1):length(array)], array[1:n]))
  }
}
