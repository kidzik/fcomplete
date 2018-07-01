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
