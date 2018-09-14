#This function writes in the correlation value
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#Use this syntax to evaluate the correlation in your data. In this case we look at the 2nd, 4th, 5th, 6th, 7th and 8th columns of the mtcars data set
#pairs(mtcars[c(2,4,5,6,7,8)], lower.panel=panel.smooth, upper.panel=panel.cor)