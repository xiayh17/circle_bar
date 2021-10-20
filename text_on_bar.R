library(clusterProfiler)
library(org.Hs.eg.db)
library(forcats)
library(ggplot2)
library(shadowtext)

# https://community.rstudio.com/t/how-to-start-text-label-from-extreme-left-in-bar-plot-using-geom-col-and-geom-text-in-ggplot2/77256

data(geneList, package = "DOSE")
de <- names(geneList)[1:100]

go_CC <- enrichGO(de, 'org.Hs.eg.db', ont="CC", pvalueCutoff=0.01)

go_CC <- clusterProfiler::simplify(go_CC)

cc_dat <- go_CC@result
cc_dat$ONTOLOGY <- rep("CC",nrow(cc_dat))

cc_dat %>% 
  mutate(task = fct_reorder(Description,Count)) %>% 
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
