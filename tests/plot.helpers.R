## Helpers
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
myggsave = function(filename,plot,width=10,height=10){
  ggsave(filename=filename,plot=plot,width=width,height=height)
}
paper.theme = theme_set(theme_grey(base_size = 36)) +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1))
plot_preds = function(obs, true = NULL, fit = NULL, title = NULL, filename = NULL, d = 30){
  cols = gg_color_hue(nrow(true))
  pp = ggplot(data.frame(x=c(0, 1)), aes(x)) + paper.theme + #ylim(-6,6) +
    xlim(0,1) + labs(x = "time", y = "value")
  for (i in 1:nrow(true)){
    df = data.frame(x=0:(d-1)/(d-1),y=true[i,],yobs=obs[i,], yfit = fit[i,])
    if (!is.null(true))
      pp = pp + geom_segment(data=df, aes(x=x,y=y,xend=dplyr::lead(x),yend=dplyr::lead(y)),
                             color=cols[i], size=1.25, linetype="dashed")
    if (!is.null(fit))
      pp = pp + geom_segment(data=df, aes(x=x,y=yfit,xend=dplyr::lead(x),yend=dplyr::lead(yfit)),
                             color=cols[i], size=3)
    df = na.omit(df)
    if (!is.null(true))
      pp = pp + geom_segment(data=df, aes(x=x,y=y,xend=x,yend=yobs),
                             color=cols[i], size=1.75, linetype="dotted")
    pp = pp + geom_point(data=df, aes(x=x,y=yobs),
                         color=cols[i], size=5)
  }
  if (!is.null(title))
    pp = pp + labs(title = title)
  print(pp)
  myggsave(filename=paste0("docs/plots/",filename,".pdf"), plot=pp)
}
