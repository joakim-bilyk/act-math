gridplot <- function(g,n) {
  df <- data.frame(
    x = g(0:n/n)
  )
  ggplot(df) + geom_point(aes(x=x,y=0)) + theme_custom() + 
    theme(axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_line(),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA)) + labs(y="Grid")
}