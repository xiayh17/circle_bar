library(clusterProfiler)
library(org.Hs.eg.db)
library(forcats)
library(ggplot2)
library(shadowtext)
library(ggstatsplot)
library(ggfun)
library(ggnewscale)
library(grid)
library(dplyr)

# https://community.rstudio.com/t/how-to-start-text-label-from-extreme-left-in-bar-plot-using-geom-col-and-geom-text-in-ggplot2/77256

data(geneList, package = "DOSE")
de <- names(geneList)[1:100]

go_CC <- enrichGO(de, 'org.Hs.eg.db', ont="CC", pvalueCutoff=0.01)

go_CC <- clusterProfiler::simplify(go_CC)

cc_dat <- go_CC@result
cc_dat$ONTOLOGY <- rep("CC",nrow(cc_dat))

cc_dat %>% 
  mutate(Description = fct_reorder(Description,Count)) %>% 
  ggplot(aes(Count, Description, fill = Description)) +
  geom_col(alpha = .6, width = .7) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(aes(label = Description),
            x = 0,
            hjust = 0, 
            size = 7) + theme_minimal()+
  theme(legend.position = "none",
        # plot.margin=margin(b=200,l=-2,unit="pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.y = element_blank())

cc_dat$group <- -1
cc_dat2 <- cc_dat
cc_dat2$Description <- paste0(cc_dat2$Description,"_fake")
cc_dat2$group <- 1

dat = rbind(cc_dat,cc_dat2)

dat$pvalue <- (-log10(dat$pvalue))*dat$group
dat=dat[order(dat$pvalue,decreasing = F),]

ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)),
                y=pvalue, 
                fill=group)) + 
  geom_bar(stat="identity",
           width=0.8) + 
  scale_fill_gradient(low="#3685af",high="#af4236",guide = "none") + 
  scale_x_discrete(name ="Pathway names") +
  scale_y_continuous(name ="-log10Pvalue",
                     expand = c(0,0)) +
  coord_flip() +
  theme_ggstatsplot()+
  theme(axis.text=element_text(face = "bold",size = 15),
        axis.title = element_text(face = 'bold',size = 15),
        plot.title = element_text(size = 20,hjust = 0.3),
        legend.position = "top",
        panel.grid = element_line(colour = 'white'))+
  ggtitle("The most enriched KEGG") 

dat <- dat %>% 
  mutate(Description = fct_reorder(Description,pvalue))

ggplot(dat, aes(pvalue, Description, fill = group)) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low="#3685af",high="#af4236",guide = "none") +
  geom_text(aes(label = stringr::str_wrap(Description,width = 30)),data = dat[which(dat["pvalue"] >= 0),],
            x = 0,
            hjust = 0, 
            size = 3) + 
  geom_text(aes(label = Description),data = dat[which(dat["pvalue"] < 0),],
            x = 0,
            hjust = 1, 
            size = 3) + 
  theme_minimal()+
  theme(legend.position = "none",
        # plot.margin=margin(b=200,l=-2,unit="pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.y = element_blank())

go_ALL <- enrichGO(de, 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.01)

split_color = RColorBrewer::brewer.pal(3,"Dark2")
bar_color = rev(RColorBrewer::brewer.pal(5,"GnBu")[3:5])

dat <- go_ALL@result %>%
  group_by(ONTOLOGY) %>%
  mutate(Description = fct_reorder(Description,Count)) %>%
  top_n(10) %>%
  mutate(myY = as.numeric(Description))

p <- ggplot(dat,aes(Count, myY, fill = p.adjust)) +
  geom_shadowtext(aes(color = ONTOLOGY, label = Description),
                  x = 0,
                  nudge_y = -0.5,
                  hjust = -0.01, show.legend = F,
                  # position=position_dodge2(width=0.9),
                  size = 4,
                  bg.colour='#ffffff',bg.r = NA) +
  scale_color_manual(values = split_color)+
  new_scale_color() +
  geom_segment(aes(x=0,y=myY-1,xend = Count,yend = myY-1,color = p.adjust),size = 1,lineend = 'round') +
  scale_y_continuous(name = "Description", breaks = dat$myY, labels = dat$Description, expand = c(0, 0.5)) +
  # geom_col(width = .15,alpha = 1,aes(color = p.adjust),
  #          position=position_dodge2(padding = 0.9)) +
  scale_color_gradientn(colours = bar_color) +
  # geom_vline(aes(xintercept = 0,colour = ONTOLOGY),size = 0.5,show.legend = F)+
  # geom_text(aes(label = Description,color = ONTOLOGY),
  #           x = 0,
  #           hjust = 0,
  #           size = 4) +
  scale_x_continuous(expand = c(0, 0))+
  theme_minimal()+
  theme(legend.position = "right",
        # plot.margin=margin(t= 100,b=2,l=-2,r= 2,unit="pt"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted"),
        panel.grid.minor = element_blank(),
        # axis.line = element_blank(),
        axis.text.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.y = element_blank())

p <- p +
  facet_grid(ONTOLOGY~., scale="free", switch="both") +
  theme(strip.background=element_roundrect(fill=NA, color=NA, r=0.31415,size = 0.5,linetype = "dotted"))

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))

k <- 1

for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  m <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- split_color[k]
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$lty <- "solid"
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- split_color[k]
  g$grobs[[i]]$grobs[[1]]$children[[m]]$children[[1]]$gp$col <- "white"
  k <- k+1
}

grid.draw(g)

ggsave(filename="ab.pdf", plot=g)
