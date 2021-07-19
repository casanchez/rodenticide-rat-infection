# Figure 2
# run calculateAgeSMI.R and analysis.R scripts first

plotDat <- clean %>% 
  pivot_longer(cols = lepto:ecoli, names_to = "pathogen", 
               values_to = "infection_status") %>% 
  mutate_at(vars(infection_status), ~as.numeric(as.character(.x))) %>% 
  mutate_at(vars(poisoned, pathogen), as.factor) %>% 
  select(poisoned, pathogen, infection_status, ageClass)

# calculate standard errors
AR_SEs <- plotDat %>% 
  group_by(pathogen, poisoned) %>% 
  summarise(meanInf = mean(infection_status),
            se = std.error(infection_status))

age_SEs <- plotDat %>% 
  group_by(pathogen, ageClass) %>% 
  summarise(meanInf = mean(infection_status),
            se = std.error(infection_status))


myColors <- c("#9e9ac8", "#3f007d")

mytheme <-
  theme_bw() +
  theme(legend.position = c(0.75, 0.7), 
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        legend.box.background = element_rect(colour = "black")) 

# infection by AR status
plotDat %>%
  ggplot(aes(x = pathogen, y = infection_status, color = poisoned, 
             group = poisoned)) +
  scale_color_manual(values = myColors, 
                     labels = c("Unexposed", "AR-exposed"),
                     name = "") +
  geom_point(position = position_jitterdodge(jitter.height = 0.03, 
                                             jitter.width = 0.6),
             size = 2, alpha = 0.4) +
  geom_point(data = AR_SEs, aes(x = pathogen, y = meanInf, group = poisoned), 
             position = position_dodge(width = 0.75),
             size = 2) +
  geom_errorbar(data = AR_SEs, 
                mapping = aes(x = pathogen, y = meanInf, ymin = meanInf - se, 
                              ymax = meanInf + se, group = poisoned), 
                position = "dodge", width = 0.75, size = 1.25) +
  mytheme +
  scale_x_discrete(labels = c("E. coli", "Leptospira spp."), name = "") +
  scale_y_continuous(name = "Infection status", 
                     breaks = seq(0, 1, 0.2), 
                     labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")) -> p1

# infection by age class
plotDat %>%
  ggplot(aes(x = pathogen, y = infection_status, color = ageClass, 
             group = ageClass)) +
  scale_color_manual(values = myColors, 
                     labels = c("30-65 days", "> 65 days"),
                     name = "") +
  geom_point(position = position_jitterdodge(jitter.height = 0.03, 
                                             jitter.width = 0.4),
             size = 2, alpha = 0.4) +
  geom_point(data = age_SEs, aes(x = pathogen, y = meanInf, group = ageClass), 
             position = position_dodge(width = 0.75),
             size = 2) +
  geom_errorbar(data = age_SEs, 
                mapping = aes(x = pathogen, y = meanInf, ymin = meanInf - se, 
                              ymax = meanInf + se, group = ageClass), 
                position = "dodge", width = 0.75, size = 1.25) +
  mytheme +
  scale_x_discrete(labels = c("E. coli", "Leptospira spp."), name = "") +
  scale_y_continuous(name = "Infection status", 
                     breaks = seq(0, 1, 0.2), 
                     labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")) -> p2

# forest plots
mytheme2 <-
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        title = element_text(size = 10)) 

# lepto
plot_model(m.lepto, vline.color = "gray90", show.values = TRUE, 
           value.offset = 0.3,
           colors = "Paired", line.size = 1.25, dot.size = 3,
           title = "Predictors of Leptospira infection status",
           axis.labels = rev(c("AR-exposed", "Sex: male", "Age: > 65 d", 
                               "SMI"))) +
  mytheme2 +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous(trans = "log10", limits = c(0.01, 100),
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c(0.01, 0.1, 1, 10, 100)) -> f1

# E coli
plot_model(m.ecoli, vline.color = "gray90", show.values = TRUE, 
           value.offset = 0.3,
           colors = "Paired", line.size = 1.25, dot.size = 3,
           title = "Predictors of E. coli shedding status",
           axis.labels = rev(c("AR-exposed", "Sex: male", "Age: > 65 d", 
                               "SMI"))) +
  mytheme2 +
  scale_y_continuous(trans = "log10", limits = c(0.1, 10),
                     breaks = c(0.1, 1, 10),
                     labels = c(0.1, 1, 10)) -> f2

# save the plot
tiff("./Figures/Fig2.tiff", height = 7, width = 8, units = "in", res = 300)
grid.arrange(p1, p2, f1, f2, nrow = 2)
dev.off()