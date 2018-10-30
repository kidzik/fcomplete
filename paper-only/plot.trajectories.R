source("tests/ggplot-theme.R")
labels = c("pelvic tilt", "pelvic obliquity", "pelvic rotation",
"hip flexion", "hip adduction", "hip rotation",
"knee flexion", "knee adduction", "knee rotation",
"ankle dorsiflexion", "foot progression")

gait.cycles = t(read.csv("/home/kidzik/Dropbox/DATA/CP/G_avg_CP.csv"))
head(gait.cycles)

mean.cycles = colMeans(gait.cycles)

library(ggplot2)
library(ggpubr)
library(gtable)
library(grid)


plots = list()
for (i in 1:11){
  y = mean.cycles[1:51 + (i-1)*51]
  x = 0:50/50
  df = data.frame(x=x, y=y)

  plots[[i]] = ggplot(df, aes(x, y)) + weartals_theme +
    geom_line(size=1.5, color="darkgreen") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0.1,0.1)) +
    theme(plot.margin=unit(rep(0.6,4), "cm")) +
    xlab("") +
    ylab(labels[[i]])
}

full.plot = ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
          plots[[5]], plots[[6]], plots[[7]], grid::rectGrob(gp=gpar(col=NA)),
          plots[[8]], plots[[9]], plots[[10]], plots[[11]],
          ncol = 4, nrow = 3)
ggsave(full.plot, filename = "docs/paper/images/kinematics.pdf",width=12,height=8)
