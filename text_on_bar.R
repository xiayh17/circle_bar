library(clusterProfiler)
library(org.Hs.eg.db)
library(forcats)
library(ggplot2)
library(shadowtext)
library(ggstatsplot)

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
