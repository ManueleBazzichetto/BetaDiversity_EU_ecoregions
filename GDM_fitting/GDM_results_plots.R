

#This script includes the R code to process and plot GDMs output.

library(ggplot2)
library(ggpubr)


# ------------------------------ grasslands


#deviance explained

expl_dev_grass_long <- data.frame(Expl_dev = as.vector(expl_dev_grass),
                                  ECO_NM = rownames(expl_dev_grass), #this gets recycled
                                  Period = rep(colnames(expl_dev_grass), each = nrow(expl_dev_grass))) 



ggplot(data = expl_dev_grass_long, aes(x = ECO_NM, y = Expl_dev, fill = Period)) +
  geom_col(position = 'identity') 







