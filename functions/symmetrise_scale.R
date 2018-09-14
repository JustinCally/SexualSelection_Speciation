#Symmetrise Scale when using facet_wrap function

symmetrise_scale <- function(p, axis = "x") {
  gb <- ggplot2::ggplot_build(p)
  type <- switch(axis, "x" = "x.range", "y" = "y.range")
  
  fname <- setdiff(names(gb$layout$panel_layout), c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y"))
  lims <- do.call(cbind, lapply(gb$layout$panel_ranges, "[[", type))
  facets <- gb$layout$panel_layout[, fname, drop = FALSE]
  
  
  # dev version ggplot2_2.2.1.9000 breaks everything yet again
  if (utils::packageVersion("ggplot2") >= "2.2.1.9") {
    fname <- setdiff(names(gb$layout$layout), c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y"))
    lims <- do.call(cbind, lapply(gb$layout$panel_params, "[[", type))
    facets <- gb$layout$layout[, fname, drop = FALSE]
  }
  
  
  lims2 <- as.vector(t(tcrossprod(apply(abs(lims), 2, max), c(-1, 1))))
  
  dummy <- stats::setNames(data.frame(facets[rep(seq_len(nrow(facets)), each = 2), ], lims2), c(fname, axis))
  
  switch(axis,
         "x" = p + geom_blank(data = dummy, aes_string(x = "x", y = "Inf"), inherit.aes = FALSE),
         "y" = p + geom_blank(data = dummy, aes_string(x = "Inf", y = "y"), inherit.aes = FALSE)
  )
}