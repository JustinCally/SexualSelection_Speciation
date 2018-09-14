gheatmap2<-function (p, data, offset = 0, width = 0.1, low = "green", high = "red", 
                    border_color = NULL, colnames = TRUE, colnames_position = "bottom", 
                    colnames_angle = 0, colnames_level = NULL, colnames_offset_x = 0, 
                    colnames_offset_y = 0, font.size = 4, hjust = 0.5) 
{
#the function is modified from ggtree::gheatmap
  require(magrittr)
  require(tidyr)
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  width <- width * (p$data$x %>% range %>% diff)/ncol(data)
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  df <- p$data
  df <- df[df$isTip, ]
  start <- max(df$x) + offset
  dd <- as.data.frame(data)
  lab <- df$label[order(df$y)]
  dd <- dd[lab, , drop = FALSE]
  dd$y <- sort(df$y)
  dd$lab <- lab
  dd <- gather(dd, variable, value, -c(lab, y))
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels = colnames(data))
  }
  else {
    dd$variable <- factor(dd$variable, levels = colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from = dd$variable, to = V2)
  mapping <- unique(mapping)
  dd$x <- V2
  dd$width <- width
  if (is.null(border_color)) {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), 
                        width = width, inherit.aes = FALSE, size=0)
  }
  else {
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), color = border_color,
                        width = width, inherit.aes = FALSE , size=0)
  }
  if (is(dd$value, "numeric")) {
    p2 <- p2 + scale_fill_gradient(low = low, high = high, 
                                   na.value = NA)
  }
  else {
    p2 <- p2 + scale_fill_discrete(na.value = NA)
  }
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    }
    else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    p2 <- p2 + geom_text(data = mapping, aes(x = to, y = y, 
                                             label = from), size = font.size, inherit.aes = FALSE, 
                         angle = colnames_angle, nudge_x = colnames_offset_x, 
                         nudge_y = colnames_offset_y, hjust = hjust)
  }
  p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
  attr(p2, "mapping") <- mapping
  return(p2)
}