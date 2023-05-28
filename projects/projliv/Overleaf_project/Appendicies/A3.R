#### Appendix A.3 ####
plot_function1 <- function(paths,num_states,s= 0, j = 1, debug = TRUE) {
  #Markov estimate
  total <- Estimate(paths,num_states,s= s, j = j, as_if = FALSE, debug = debug)
  #As-If-Markov estimate
  total2 <- Estimate(paths,num_states,s= s, j = j, as_if = TRUE, debug = debug)
  plotdf <- total$p_con
  colnames(plotdf)[2:4] <- paste0("j=",1:3,", Markov")
  plotdf <- plotdf%>% reshape2::melt(., id = "Time")
  plotdf2 <- total2$p_con
  colnames(plotdf2)[2:4] <- paste0("j=",1:3,", As-If")
  plotdf2 <- plotdf2%>% reshape2::melt(., id = "Time")
  times <- s+0:1000*(10-s)/1000
  plotdf3 <- data.frame(Time = times,
                        matrix(unlist(lapply(times, function(t) ((1:num_states == j)*1)%*%expm::expm(2*M*(log(1+0.5*t)-log(1+0.5*s))))),ncol=3,byrow=TRUE))
  colnames(plotdf3)[2:4] <- paste0("j=",1:3,", True")
  plotdf3 <- plotdf3 %>% reshape2::melt(., id = "Time")
  asif_n <- sum(total2$I[1,])-total2$I[1,1]
  markov_n <- sum(total$I[1,])-total$I[1,1]
  ggplot() +
    geom_step(data = plotdf ,mapping = aes(x=Time, y = value,col = variable)) + 
    geom_step(data = plotdf2,mapping = aes(x=Time, y = value,col = variable),linetype = "dashed") +
    geom_line(data = plotdf3,mapping = aes(x=Time, y = value,col = variable),linetype = "dotted",size=1) +
    theme_bw() +
    labs(title = TeX(paste0("Occupation probabilities under assumption $Z_",s,"=",j,"$")),
         y =  TeX(paste0("$P(Z_t=j|Z_",s,"=",j,")$")),
         subtitle = paste0("As-if estimate based on ",asif_n," observations,\nMarkov based on ",markov_n," observations")) +
    theme(legend.title = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic")) +
    scale_color_manual(values=c("#FA8072", "#FF0000","#8B0000",
                                "#7B68EE", "#1E90FF","#00008B",
                                "#3CB371", "#32CD32","#006400"))
}
plot1 <- plot_function1(paths,3,s= 0, j = 1)
plot2 <- plot_function1(paths,3,s= 1, j = 1)
plot3 <- plot_function1(paths,3,s= 3, j = 1)
plot4 <- plot_function1(paths,3,s= 6, j = 1)
ggarrange(plotlist = list(plot1,plot2,plot3,plot4),ncol = 2,nrow=2)