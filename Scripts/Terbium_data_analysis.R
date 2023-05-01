# Analyzing Terbium data

require(ggplot2)
require(tidyr)
require(dplyr)
require(lazyeval)


ans.l <- read.csv("C:/path/Terbium_NAME.l.csv",header = TRUE)

theme_mds_leg <-  theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line.x = element_line(colour = "black"), 
                        axis.line.y = element_line(colour = "black"), 
                        axis.ticks = element_line(colour = "black"), 
                        axis.text = element_text(color="black", size = 8),
                        strip.background = element_blank(),
                        legend.background = element_blank(),
                        legend.key = element_blank())

ntPalette <- c("#E69F00", "#56B4E9", "#F0E442","#009E73") # colorblind palette, but roughly same as old one
#ntPalette <- c("#C40233", "blue", "#FFCC00", "#009900") #This is the older palette

# want to get rid of mostly unused categories, Pin, Pdel, PtoX; comment out if needed
temp.l <- ans.l %>%
  filter(stat %in% c("Pstop"))
# weighted average not needed since these are not replicates

# now make column with values normalized to the top 10% of values; doing with nstat

temp.l <- temp.l %>% mutate(norm_scale = nstat/mean(temp.l$nstat[temp.l$nstat>=quantile(temp.l$nstat, 0.9, na.rm=TRUE)], na.rm=TRUE))
write.csv(temp.l, file = "Terbium_norm_NAME.l.csv")

# Make line plot of RT from different targets and reagents
v1 = 0
v2 = X
temp.l %>%
  filter(stat %in% c("Pstop")) %>%
  ggplot(aes(x=nt, y = RT, fill = ref)) +
  facet_grid(reagent + name ~ stat) +
  geom_line(stat = "identity") + scale_y_continuous(trans = 'log10')
  #coord_cartesian(ylim = c(0, 100000), xlim = c(v1,v2)) + 
  theme_mds_leg
#ggsave("NAME.pdf", width = 15, height = 10, units = "cm") # Saving the plot to about the right size is the key to get the fonts the way we want them.

# Bar plot of Pmod at nts
v1 = 0
v2 = X
temp.l %>%
  ggplot(aes(x=nt, y = norm_scale, fill = ref)) +
  facet_grid(reagent + name ~ .) + # , scales = "free"
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ntPalette,  # NT palette defined in functions file, conventional colors for acgt.
                    breaks = c("A", "C", "G", "T"),
                    labels = c("A", "C", "G", "U")) + # redefine the names to make the Ts Us.
  scale_color_manual(values = ntPalette,# Have to define this twice because the labels and the fill are different.
                     breaks = c("A", "C", "G", "T"),
                     labels = c("A", "C", "G", "U")) +
  geom_text(aes(label=ifelse(norm_scale > 0.5,as.character(nt),'')),size = 2.5, hjust=0,vjust=-1) +
  coord_cartesian(xlim = c(v1, v2), ylim = c(-1, 10))+
  theme_mds_leg
  ggsave("NAME.eps", width = 25, height = 15, units = "cm", dpi = 400) # Saving the plot to about the right size is the key to get the fonts the way we want them.


