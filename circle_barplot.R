library(clusterProfiler)
library(org.Hs.eg.db)
library(forcats)
library(ggplot2)
library(shadowtext)

data(geneList, package = "DOSE")
de <- names(geneList)[1:100]

go_CC <- enrichGO(de, 'org.Hs.eg.db', ont="CC", pvalueCutoff=0.01)

go_CC <- clusterProfiler::simplify(go_CC)

cc_dat <- go_CC@result
cc_dat$ONTOLOGY <- rep("CC",nrow(cc_dat))

kong <- function(data,len = 2) {
  df2 <- data[1:len,]
  df2$Description <- rep(" ",len)
  df2$Count <- rep(0, len)
  res <- rbind(data,df2)
  return(res)
}

data2 <- kong(data = cc_dat)
data3 <- data2 %>%
  mutate(Description = fct_reorder(Description,Count))

color = "#f27e12"
label = "CC"
ggplot(data3, aes(x = Description, y = Count)) +
  geom_bar(width = 0.5, stat="identity",
           fill = color) +
  # geom_text(data = df2, aes(label = chars, x = x, y = y, angle = -angle, cex = 1), fontface = 2, family = "sans")+
  coord_polar(theta = "y",start = 0,clip = 'off' ) +
  xlab("") + ylab("") +
  ylim(c(0,max(data3$Count)*4/3)) +
  # ggtitle("Top Product Categories Influenced by Internet") +
  geom_text(data = cc_dat, size = 5, angle = 0,colour = color,hjust = 1.01,
            aes(x = Description, y = 0, label = Description)) +
  geom_shadowtext(data = cc_dat,aes(label=Count), position=position_dodge2(width=0.9),family = 'Times',
                  vjust=0.5,colour = color, bg.colour='#19b8fa',bg.r = 0.03,size = 4) + 
  geom_text(label=label, x=.5, y=.5, colour = color,size=8) +
  theme_minimal()+
  theme(legend.position = "none",
        # plot.margin=margin(b=200,l=-2,unit="pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())
