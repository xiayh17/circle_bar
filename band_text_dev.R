mtcars
df <- data.frame(des = rep(NA,3))
df$des <- c(" ",rownames(mtcars)[1:2])
df$value <- c(0:2)
library(ggplot2)
ggplot(df,aes(x = value,y = value)) + geom_bar(stat = "identity") +
  geom_text(aes(label = des),cex=5) +
  coord_polar('y')

labels <- df$des
chars = strsplit(labels, "")
strw = lapply(chars, strwidth,units =  "inches",cex=1,font = 2 ,family= "sans")

strwidth(labels[[2]],units =  "inches")

w_2 <- strw[[2]]

ls <- list()
for (i in seq_along(w_2)) {
  
  if (i == 1) {
    ls[[i]] <- w_2[[i]]
  }else {
    ls[[i]] <- sum(w_2[1:i])
  }
}

y_theta <- (unlist(ls) - min(unlist(ls)))/pi

df2 <- data.frame(
  chars = chars[[2]],
  x = rep(1,9),
  y = y_theta,
  angle = -y_theta * 180
)

ggplot(df,aes(x = value,y = value)) + geom_bar(stat = "identity") +
  geom_text(data = df2, aes(label = chars, x = x, y = y, cex = 1,angle = angle), fontface = 2, family = "sans", size= 3) +
  coord_polar('y')


