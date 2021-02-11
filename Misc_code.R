#Phyloseq PCOAs

ps.mtl.sans.phasco <- subset_samples(ps, Population != "Mainland SA"
                                     & Population != "Kangaroo Island"
                                     & Koala != "Phasco")

ord.unw.uni <- ordinate(ps.mtl.sans.phasco, "PCoA", "unifrac", weighted=F)

#Weighted UniFrac
ord.w.uni <- ordinate(ps.mtl.sans.phasco, "PCoA", "unifrac", weighted=T)

#Axes 1/2
unwt.unifrac.1.2 <- plot_ordination(ps.mtl.sans.phasco, ord.unw.uni, color="Koala", axes = c(1, 2))
w.unifrac.1.2 <- plot_ordination(ps.mtl.sans.phasco, ord.w.uni, color="Koala", axes = c(1, 2))
#Axes 1/3
unwt.unifrac.1.3 <- plot_ordination(ps.mtl.sans.phasco, ord.unw.uni, color="Koala", axes = c(1, 3))
w.unifrac.1.3 <- plot_ordination(ps.mtl.sans.phasco, ord.w.uni, color="Koala", axes = c(1, 3))

#Plot it
FigX1a <- unwt.unifrac.1.2 +
  #stat_ellipse() +
  scale_x_reverse() +
  geom_point(size=9) +
  geom_text(aes(label=Collection_number), color="black", hjust=2.5, vjust=-1) +
  geom_path(size = 1) +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

FigX1b <- unwt.unifrac.1.3 +
  #stat_ellipse() +
  scale_x_reverse() +
  geom_point(size=9) +
  geom_text(aes(label=Collection_number), color="black", hjust=2.5, vjust=-1) +
  geom_path(size = 1) +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

#Create figure X:

ggarrange(FigX1a, FigX1b, nrow = 2, common.legend = TRUE, legend="right",
          labels = c("A)", "B)"), font.label = list(size=20, face="bold", color="black"))








