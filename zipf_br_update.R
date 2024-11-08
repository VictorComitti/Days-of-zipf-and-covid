library(coronabr)
library(dplyr)
library(stats4)
library(ggplot2)
library(data.table)
library(officer)
library(flextable)
library(latex2exp)
library(gravity)
library(broom)
library(lmtest)
library(tidyr)

rm(list=ls(all=TRUE))
dados_br<-get_corona_br()

dados_br$date <- as.Date(dados_br$date)
dados_br<-subset( dados_br,  dados_br$city!="Importados/Indefinidos")
dados_br_covid <- dados_br %>% 
  filter(place_type == "city") %>% select("date", "city", "last_available_confirmed") %>%
  filter(date >= "2020-03-30") %>%
  filter(last_available_confirmed>0) %>% 
  setDT()

datas <- unique(dados_br_covid$date)
###zipf brasil

dados_br_cities <- dados_br %>% 
  filter(place_type=="city") %>% 
  select("estimated_population_2019") %>% 
  unique() %>% 
  arrange(., desc(.)) %>%
  mutate(rank = 1:nrow(.)) %>% 
  na.omit() %>% 
  setDT()


pop <- dados_br_cities$estimated_population_2019
rank <- dados_br_cities$rank

dados <- as.data.frame(cbind(pop, rank))
#setwd("C:/Users/Uriel/Desktop/Pesquisa/Shikida")

write.csv(dados,"covidbr.csv", row.names = TRUE)

# OLS Gabaix

zipf_pop_OLS <- lm(log(rank-0.5)~log(pop))

se <- sqrt(diag(vcov(zipf_pop_OLS)))[2]
pVal <- anova(zipf_pop_OLS)$'Pr(>F)'[1]

bptest(zipf_pop_OLS)$p.value


#PPML
pop_ppml <- ppml("estimated_population_2019", "rank",  
                 additional_regressors = c(), 
                 data = dados_br_cities)
summary(pop_ppml)
sqrt(diag(vcov(pop_ppml)))[2]
tidy(pop_ppml)$p.value[2]


extract_data <- function(dia){
  extracted_data = dados_br_covid %>% 
    filter(date == dia) %>% 
    na.omit() %>% 
    arrange(., desc(last_available_confirmed)) %>%
    mutate(rank = 1:nrow(.))
  return(extracted_data)
}


##estimacao via OLS with gabaix correction
zipf_OLS<-function(dia){
  covid <- extract_data(dia)
  casos <- covid$last_available_confirmed
  rank <- covid$rank-0.5
  fit <- lm(log(rank)~log(casos))
  coef <- coef(fit)[2]
  return(coef)
}

OLS<-lapply(datas, zipf_OLS) %>% 
  unlist() %>% as.numeric() %>% 
  abs()



##PPML estimation



zipf_PPML <- function(dia){
  covid <- extract_data(dia)
  fit <- ppml("last_available_confirmed", "rank",  additional_regressors = c(), data=covid)
  coef <-coef(fit)[2]
  return(coef)
}

PPML <- lapply(datas, zipf_PPML) %>% 
  unlist() %>% 
  as.numeric() %>% 
  abs()


days <- seq(1,length(PPML),1)
df <- as.data.frame(cbind(days, PPML, OLS))
df$days <- datas

meltdfPPML <-  reshape2::melt(df,id="days",
                              variable.name = "Estimator",  
                              value.name = "alpha")
ggplot(meltdfPPML,
       aes(x=days,y=alpha,colour=Estimator,group=Estimator)) + 
  geom_line(size=1.4)+xlab("Date")+
  ylab(TeX(" $\\hat{\\alpha}$"))+
  geom_hline(yintercept=0.873, 
             linetype="dashed", 
             color = "red")

####plot apenas do PPML

df_PPML<-as.data.frame(cbind(days, PPML))
df_PPML$days<-datas
meltdf_PPML <-  reshape2::melt(df_PPML,id="days", 
                               variable.name = "Estimator",  
                               value.name = "alpha")

ggplot(meltdf_PPML,aes(x=days,
                       y=alpha,
                       colour=Estimator,
                       group=Estimator)) + 
  geom_line(size=1.4)+
  xlab("Date")+
  ylab(TeX(" $\\hat{\\alpha}$")) +
  geom_hline(yintercept=0.873,
             linetype="dashed", 
             color = "blue")+
  ylim(0.7,1.4)+
  theme_grey(base_size = 20)+ 
  theme(legend.title=element_blank())



df_PPML_melt <- df_PPML %>% 
  gather(key = "PPML Estimate", value = "alpha", -days)

ggplot(df_PPML_melt, 
       aes(x = days, y = alpha)) + 
  geom_line(aes(color = "PPML Estimate"), size = 1) + 
  scale_color_manual(values = c("red"))+
  xlab("Date")+ylab(TeX(" $\\hat{\\alpha}$"))+ 
  theme(legend.title=element_blank())+
  theme_grey(base_size = 20)+
  theme(legend.title=element_blank())+
  ylim(0.78,1.4)+
  geom_hline(yintercept=0.873, linetype="dashed", color = "blue")+
  annotate("text", x=df_PPML_melt$days[110], y=0.87, 
           label = "Population alpha", vjust = -0.5)



#alpha_pop <- data.frame(yintercept=0.873, alpha_pop=factor("Population alpha"))

#p + geom_hline(aes(yintercept=yintercept, linetype=alpha_pop), data=alpha_pop, show.legend=TRUE) 

###tabela dias, casos, alpha

tabela<-function(dia){
  covid<-extract_data(dia)
  casos<-covid$last_available_confirmed
  n<-length(casos)
  rank<-covid$rank
  fit_G<-lm(log(rank-0.5)~log(casos))
  alpha_G<-round(abs(coef(fit_G)[2]), digits=4)
  se_G <- round(sqrt(diag(vcov(fit_G)))[2], digits=4)
  bp<-round(bptest(fit_G)$p.value, digits=4)
  pvalue_G<-tidy(fit_G)$p.value[2]
  fit_ppml<-ppml("last_available_confirmed", "rank",  additional_regressors = c(), data=covid)
  alpha_ppml<-round(abs(coef(fit_ppml)[2]),  digits=4)
  se_ppml <- round(sqrt(diag(vcov(fit_ppml)))[2], digits=4)
  pvalue_ppml<-tidy(fit_ppml)$p.value[2]
  day<-dia
  saida<-list(day,n,alpha_G, se_G, pvalue_G, bp, alpha_ppml, se_ppml, pvalue_ppml)
  return(saida)
}

table<-lapply(datas, tabela)


table<- matrix(unlist(table), ncol = 9, byrow = TRUE)
table #todas as estimativas de alpha são significativas a 1%

colnames<-c("Date", "Number of municipalities", "alpha_Gabaix", "SE", "pvalue_G", "bpteste", "alpha_PPML",  "SE_PPML", "pvalue_ppml")
colnames(table)<-colnames
table<-as.data.frame(table)
table$Date<-datas
table<-subset(table, select = -c(pvalue_G, pvalue_ppml))
table$alpha_Gabaix<-paste0(table$alpha_Gabaix, "***")
table$alpha_PPML<-paste0(table$alpha_PPML, "***")
table
ft <- flextable(data = table) %>% 
  theme_booktabs() %>% 
  autofit
ft



tmp <- tempfile(fileext = ".docx")

read_docx() %>% 
  body_add_flextable(ft) %>% 
  print(target = tmp)

browseURL(tmp)



