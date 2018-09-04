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

plot_preds = function(obs, true = NULL, fit = NULL, title = NULL, filename = NULL, ylim = NULL, width = 10, height = 10){
  ylim = c(min(obs),max(obs))
  d = ncol(obs)
  cols = gg_color_hue(nrow(obs))
  pp = ggplot(data.frame(x=c(0, 1)), aes(x)) + paper.theme + #ylim(-6,6) +
    xlim(0,1) + labs(x = "time", y = "value")
  if (!is.null(ylim)){
    pp = pp + ylim(ylim)
  }
  for (i in 1:nrow(obs)){
    yt = NA
    yo = NA
    yf = NA
    if (!is.null(true))
      yt = true[i,]
    if (!is.null(obs))
      yo = obs[i,]
    if (!is.null(fit))
      yf = fit[i,]

    df = list(x=0:(d-1)/(d-1),y=yt,yobs=yo, yfit = yf)
    df = data.frame(df)
    if (!is.null(true))
      pp = pp + geom_segment(data=df, aes(x=x,y=y,xend=dplyr::lead(x),yend=dplyr::lead(y)),
                             color=cols[i], size=1.25, linetype="dashed")
    if (!is.null(fit))
      pp = pp + geom_segment(data=df, aes(x=x,y=yfit,xend=dplyr::lead(x),yend=dplyr::lead(yfit)),
                             color=cols[i], size=3)
    df_yyobs = na.omit(df[,c("x","y","yobs")])
    if (!is.null(true))
      pp = pp + geom_segment(data=df_yyobs, aes(x=x,y=y,xend=x,yend=yobs),
                             color=cols[i], size=1.75, linetype="dotted")
    df_yobs = na.omit(df[,c("x","yobs")])
    pp = pp + geom_point(data=df_yobs, aes(x=x,y=yobs),
                         color=cols[i], size=5)
  }
  if (!is.null(title))
    pp = pp + labs(title = title)
  print(pp)
  myggsave(filename=paste0("docs/plots/",filename,".pdf"), plot=pp, width = width, height = height)
  pp
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

weartals_theme = theme_classic() + theme(text = element_text(size=18), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
