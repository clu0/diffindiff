library(tidyverse)
library(RColorBrewer)

args=(commandArgs(TRUE))


colors = brewer.pal(8, "Dark2")

breaks = c("Samples means","OLS","direct_learner",
           "TR")
colors = c(colors[6],colors[3],colors[1],colors[3])
shapes = c(12,11,15,19)
sizes = c(2,1.5,2.2,2)

plotsize = function(x,y) options(repr.plot.width=x, repr.plot.height=y)

out = read_csv(paste0("output_const_full", ".csv"))

label_wrap <- function(variable, value) {
  paste0("Setup ", value)
}
out %>%
  select(-n, -p, -X1) %>%
  gather(alg,mse, -setup, -oracle_learner, -constant_eff, -naive_learner, -direct_learner, -constant_eff) %>%
  ggplot(aes(x=log(oracle_learner), y=log(mse), color=alg, shape=alg, size=alg,breaks=alg)) +
  scale_size_manual(breaks=breaks, values=sizes)+
  scale_shape_manual(breaks=breaks, values=shapes)+
  scale_colour_manual(breaks=breaks, values = colors) +
  theme(legend.title = element_blank()) +
  geom_point() +
  geom_abline(slope=1) +
  facet_wrap(~setup, ncol=2, scales="free", labeller = label_wrap) +
  labs(x = "log oracle mean-squared error", y = "log mean-squared error")
ggsave(paste0("const_test", ".png"), height = 10, width = 10)


