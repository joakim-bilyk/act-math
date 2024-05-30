plotDiff <- function(path) {
  plotdf <- path %>% as.data.frame()
  ggplot(plotdf) + geom_line(aes(x= t,y=X)) + theme_custom()
}