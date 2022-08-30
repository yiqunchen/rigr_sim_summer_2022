library(rigr)
library(tidyverse)
library(latex2exp)

set.seed(2022)
data(mri)
n_sim <- 1000
beta_0 <- 1
beta_1 <- 0
beta_2_seq <- seq(0,0.8,by=0.1)
beta_3 <- 0

p_val_collection_735_wald <- rep(NA, n_sim)
p_val_collection_735_naive <- rep(NA, n_sim)
p_val_collection_735_LRT <- rep(NA, n_sim)
p_val_collection_100_wald <- rep(NA, n_sim)
p_val_collection_100_naive <- rep(NA, n_sim)
p_val_collection_100_LRT <- rep(NA, n_sim)
p_val_collection_30_wald <- rep(NA, n_sim)
p_val_collection_30_naive <- rep(NA, n_sim)
p_val_collection_30_LRT <- rep(NA, n_sim)
p_val_collection_10_wald <- rep(NA, n_sim)
p_val_collection_10_naive <- rep(NA, n_sim)
p_val_collection_10_LRT <- rep(NA, n_sim)

df_p0_summary_wald <- data.frame()
df_p0_summary_naive<- data.frame()
df_p0_summary_LRT <- data.frame()

for (beta_2 in beta_2_seq){
  # full 735 samples
  for (i in c(1:n_sim)){
    cat(i,"\n")
    y_sim <- rpois(n = nrow(mri), beta_0 + mri$weight*beta_1 +
                     mri$height*beta_2+mri$alcoh*beta_3)
    curr_reg <- regress("rate", y_sim~weight+height+alcoh,data=mri)
    null_reg <- regress("rate", y_sim~weight,data=mri)
    curr_p_wald <- anova(null_reg, curr_reg,test = "Wald",robustSE = TRUE)$printMat[,4]
    curr_p_naive <- anova(null_reg, curr_reg,test = "Wald",robustSE = FALSE,
                          useFdstn = FALSE)$printMat[,3]
    curr_p_LRT <- anova(null_reg, curr_reg,test = "LRT")$printMat[,3]
    p_val_collection_735_wald[i] <- curr_p_wald
    p_val_collection_735_naive[i] <- curr_p_naive
    p_val_collection_735_LRT[i] <- curr_p_LRT
  }
  
  # sub 100
  for (i in c(1:n_sim)){
    mri_sub <- mri %>% slice_sample(n=100)
    cat(i,"\n")
    y_sim <- rpois(n = nrow(mri_sub), beta_0 + mri_sub$weight*beta_1 +
                     mri_sub$height*beta_2+mri_sub$alcoh*beta_3)
    curr_reg <- regress("rate", y_sim~weight+height+alcoh,data=mri_sub)
    null_reg <- regress("rate", y_sim~weight,data=mri_sub)
    curr_p_wald <- anova(null_reg, curr_reg,test = "Wald",robustSE = TRUE)$printMat[,4]
    curr_p_naive <- anova(null_reg, curr_reg,test = "Wald",robustSE = FALSE,
                          useFdstn = FALSE)$printMat[,3]
    curr_p_LRT <- anova(null_reg, curr_reg,test = "LRT")$printMat[,3]
    p_val_collection_100_wald[i] <- curr_p_wald
    p_val_collection_100_naive[i] <- curr_p_naive
    p_val_collection_100_LRT[i] <- curr_p_LRT
  }
  
  # sub 30
  for (i in c(1:n_sim)){
    mri_sub <- mri %>% slice_sample(n=30)
    cat(i,"\n")
    y_sim <- rpois(n = nrow(mri_sub), beta_0 + mri_sub$weight*beta_1 +
                     mri_sub$height*beta_2+mri_sub$alcoh*beta_3)
    curr_reg <- regress("rate", y_sim~weight+height+alcoh,data=mri_sub)
    null_reg <- regress("rate", y_sim~weight,data=mri_sub)
    curr_p_wald <- anova(null_reg, curr_reg,test = "Wald",robustSE = TRUE)$printMat[,4]
    curr_p_naive <- anova(null_reg, curr_reg,test = "Wald",robustSE = FALSE,
                          useFdstn = FALSE)$printMat[,3]
    curr_p_LRT <- anova(null_reg, curr_reg,test = "LRT")$printMat[,3]
    p_val_collection_30_wald[i] <- curr_p_wald
    p_val_collection_30_naive[i] <- curr_p_naive
    p_val_collection_30_LRT[i] <- curr_p_LRT
  }
  
  
  df_p0_summary_wald <- rbind(df_p0_summary_wald, 
                         data.frame(
                           #p_val_collection_10 = p_val_collection_10_wald,
                                                   p_val_collection_30 = p_val_collection_30_wald,
                                                   p_val_collection_100 = p_val_collection_100_wald,
                                                   p_val_collection_735 = p_val_collection_735_wald,
                                                   beta_2=beta_2,
                                                   type="Wald robust"))
  
  df_p0_summary_naive <- rbind(df_p0_summary_naive, 
                              data.frame(
                                #p_val_collection_10 = p_val_collection_10_naive,
                                         p_val_collection_30 = p_val_collection_30_naive,
                                         p_val_collection_100 = p_val_collection_100_naive,
                                         p_val_collection_735 = p_val_collection_735_naive,
                                         beta_2=beta_2,
                                         type="Wald naive"))
  
  df_p0_summary_LRT <- rbind(df_p0_summary_LRT, 
                              data.frame(
                                #p_val_collection_10 = p_val_collection_10_LRT,
                                         p_val_collection_30 = p_val_collection_30_LRT,
                                         p_val_collection_100 = p_val_collection_100_LRT,
                                         p_val_collection_735 = p_val_collection_735_LRT,
                                         beta_2=beta_2,
                                         type="LRT"))
}

df_p0_plot_data <- rbind(df_p0_summary_wald, 
      df_p0_summary_naive,
      df_p0_summary_LRT)

save(df_p0_plot_data,
     file="~/Desktop/summer_2022/power_anova_glm_poisson.RData")

df_p0_plot <- df_p0_plot_data %>% 
  pivot_longer(-c(beta_2, type), names_to = "p_type", values_to = "p_values") %>%
  mutate(p_type = as.factor(p_type)) %>%
  group_by(beta_2, p_type, type) %>%
  summarise(power=mean(p_values<=0.05),
            sd_power = sqrt(power*(1-power)/n()))


p_power <- df_p0_plot %>% 
  ggplot(aes(x=beta_2,y=power,
             colour=p_type,
             linetype=type))+
  geom_point()+
  geom_pointrange(aes(ymin = power+sd_power,
                      ymax = power-sd_power))+
  geom_line()+
  geom_abline(aes(slope=0,intercept=0.05))+
  ylab(TeX('Power at $\\alpha$=$0.05$'))+
  xlab(TeX('$\\beta_1$'))+
  theme_bw(base_size=18) +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        legend.position="bottom",
        legend.title = element_text(size=12),
        legend.key.size = unit(2.5, "lines"),
        axis.text = element_text(size=18),
        legend.text = element_text(size=13,hjust = 0),
        axis.title=element_text(size=18))+
  guides(colour = guide_legend(override.aes = list(size=1),nrow=2,byrow=TRUE),
         linetype = guide_legend(nrow=2,byrow=TRUE,override.aes = list(length=2)))+
  scale_colour_brewer(name="",palette="Dark2",
                      labels = unname(TeX(c('$n=100$','$n=30$','$n=735$'))))+
  scale_linetype(name = "",labels=c("LRT","Wald w/ Naive SE","Wald w/ robust SE & F dist."))


png("~/Desktop/summer_2022/glm_poisson_anova_power.png",
    width=10,height=9, units= "in", res=400)
p_power
dev.off()



