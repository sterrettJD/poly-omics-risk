

set.seed(1)
df <- data.frame(y=rnorm(10),x=c(1:5,1:3,1:2),group=c(rep("a",5),rep("b",3),rep("c",2)),name=c(paste("a",1:5,sep=""),paste("b",1:3,sep=""),paste("c",1:2,sep="")))
df$ymin <- df$y-runif(10,0.5,0.7)
df$ymax <- df$y+runif(10,0.5,0.7)

p1 <- ggplot(df,aes(y = name, x = y))+
  geom_point()+
  facet_wrap(~group,ncol=3,scales="free")+
  geom_segment(aes(x = ymin, xend = ymax, yend = name))+
  geom_vline(lty=2, aes(xintercept=0), colour = 'red')

p1