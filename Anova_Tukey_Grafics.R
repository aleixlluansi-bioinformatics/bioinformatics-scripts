library("ggplot2")
library("multcompView");library("dplyr")

View(my_data)
anova<-aov(IBS.SSS~Group,data=my_data)
summary(anova)
tukey<-TukeyHSD(anova)
print(tukey)
cld<-multcompLetters4(anova,tukey)
print(cld)

TK<-group_by(my_data,Group)%>%
  summarise(mean=mean(IBS.SSS),quant=quantile(IBS.SSS,probs=0.75))%>%
  arrange(desc(mean))

ggplot(my_data,aes(Group,IBS.SSS))+
  geom_boxplot()+
  labs(x="Group", y= "IBS-SSS")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_text(data=)

df1 <- data.frame(a = c(1, 1:2,2), b = c(350, 360, 360, 350))
df2 <- data.frame(a = c(1, 1,3, 3), b = c(370, 380, 380, 380, 370))


pp + geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 365, label = "***", size = 5) +
     geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 385, label = "***", size = 5)
     