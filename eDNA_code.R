library(Matrix)
library(dplyr)
library(tidyverse)
library(tidybayes)
library(modelr)
library(ggplot2)
library(readxl)
library(janitor)
library(brms)
library(rstan)
library(StanHeaders)


eDNA_compiled_results <- read_excel("~/Desktop/eDNA/qPCR results/eDNA_compiled_results.xlsx") %>% 
  clean_names()

clean_results <- eDNA_compiled_results %>% replace_na(list(ct_mean=0, ct_sd=0,quantity=0,quantity_mean=0,quantity_sd=0)) %>% 
  ungroup %>% 
  mutate(quant_cat=case_when(quantity==0~0, TRUE~1),
         above_below=case_when(location=="above"~0,
                               location=="below"~1),
         river_cat=case_when(river=="big sioux"~0,
                             river=="vermillion"~1),
         quant_1000=quantity/1000,
         quant_mean_s=scale(quantity_mean))


##### extra unnecessary stuff #####

# mod_prob_detect <- brm(quant_cat~quant_mean_s,
#     data=clean_results,
#     family=bernoulli(),
#     prior=c(prior(normal(-10,1), class= "Intercept"),
#             prior(normal(0,1), class= "b")),
#             chains=4, iter = 2000)
# 
# plot(conditional_effects(mod_prob_detect, re_formula=NULL), points = T)
# 
# pp_check(mod_prob_detect)
# 
# cond_effect_mod_prob_detect <- conditional_effects(mod_prob_detect)
# cond_effect_mod_prob_detect$quant_mean_s
# 
# cond_effect_mod_prob_detect$quant_mean_s %>% 
#   ggplot(aes(x=quant_mean_s)) +
#   geom_pointrange(aes(y=estimate__, ymin=lower__, ymax=upper__))+
#   geom_point(data = mod_prob_detect$data, aes(x=quant_mean_s, y=quant_cat))+
#   theme_default()
# 
# cond_data_mod_prob_detect <- mod_prob_detect$data %>% distinct(quant_mean_s)
# 
# posts_mod_prob_detect <- add_epred_draws(mod_prob_detect, newdata= mod_prob_detect$data %>% 
#                                  distinct(quant_mean_s) , re_formula = NA)
# 
# 
# posts_mod_prob_detect <- add_predicted_draws(mod_prob_detect, newdata=mod_prob_detect$clean_results %>% 
#                                        distinct(mod_prob_detect,quant_cat) , re_formula = NA)

# 
# d_quant_mean_s <- d %>% distinct(length_mm, length_s)
# 
# PosteriorGSIlength <- posts_gsi_all %>%
#   group_by(length_s, lab_sex) %>%
#   left_join(d_lengthgsi) %>%
#   median_qi(.prediction) %>%
#   mutate(length_mm = (length_s*sd(d$length_mm)) + mean(d$length_mm)) %>%
#   ggplot(aes(x =length_mm, y = .prediction, fill = lab_sex)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = .lower, ymax = .upper),
#               alpha = 0.2) +
#   geom_point(data = d,
#              aes(y = gsi)) +
#   labs(title= "Blue Sucker GSI Prediction",
#        subtitle="Blue and pink bars incorporate the variation in individuals",
#        x="Length (mm)",
#        y="Predicted GSI")
# 
# ggsave(PosteriorGSIlength, file = "plots/PosteriorGSIlength.png", dpi = 750, width = 7, height = 5, units = "in")


##### real model #####

# next step: same model plus below/above dam after adjusting for quantity, 
# river as random effect and below/above in separate column

get_prior(quant_cat~quant_1000*location*river + (1|sample_name) + (1|target_name),
          data=clean_results,
          family=bernoulli())

mod_prob_ab <- brm(quant_cat~quant_1000*location*river + (1|sample_name) + (1|target_name),
                   data=clean_results,
                   family=bernoulli(),
                   prior=c(prior(normal(-10,1), class= "Intercept"),
                           prior(normal(0,1), class= "b")),
                   # sample_prior = "only",
                   chains=4, iter = 2000)


preds = mod_prob_ab$data %>% 
  distinct(location, river, target_name) %>% 
  expand_grid(quant_1000 = mean(clean_results$quant_1000), 
              sample_name = "new") %>% 
  add_epred_draws(mod_prob_ab, allow_new_levels = T) 


preds_graph <- preds %>% 
  ggplot(aes(x = location, y = .epred, fill = river)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~target_name)+
  labs(y="Probability of a positive sample") +  
  geom_point(data = mod_prob_ab$data, 
             aes(y = quant_cat, color = river, 
                 group = interaction(river, location)), 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0),
             shape = 21, alpha = 0.2)

ggsave(preds_graph, file = "plots/positive_predictions.png", dpi = 750, 
       width = 7, height = 5, units = "in")

conditional_effects(mod_prob_ab, conditions = tibble(target_name = "BHC"))

pp_check(mod_prob_ab)

ggplot(data=clean_results,aes(x=quant_1000, y=quant_cat))+
  geom_point()+
  facet_wrap(location~river)

clean_results %>% 
  group_by(location, river, target_name) %>% 
  summarize(positive = sum(quant_cat))

clean_results %>% filter(location == "above") %>% 
  filter(river == "vermillion")

clean_results %>% filter(is.na(river)) %>% view



pred_quant = mod_prob_ab$data %>% 
  distinct(location, river, target_name) %>% 
  expand_grid(quant_1000 = clean_results$quant_1000, 
              sample_name = "new") %>% 
  add_epred_draws(mod_prob_ab, allow_new_levels = T) 

pred_quant_summarized <- pred_quant %>% 
  group_by(location,target_name,river,quant_1000) %>% 
  median_qi(.epred,.width=c(0.75))

pred_quant_plot <- pred_quant_summarized %>% 
  mutate(river_location=paste(river,target_name))%>% 
  ggplot(aes(x = quant_1000, y = .epred, fill = location)) +
  geom_line(aes(color=location))+
  geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.2)+
  facet_wrap(river_location~location,ncol=2)+
  labs(y="Probability of a positive sample")

# how do i backcalculate the quantities so they are not at a magnitude of 1000?

ggsave(pred_quant_plot, file = "plots/loc_spec_probs.png", dpi = 750, 
       width = 7, height = 5, units = "in")


# Sensitivity analysis below, doubled the standard deviation 

mod_prob_ab_sens <- brm(quant_cat~quant_1000*location*river + (1|sample_name) + (1|target_name),
                        data=clean_results,
                        family=bernoulli(),
                        prior=c(prior(normal(-10,2), class= "Intercept"),
                                prior(normal(0,2), class= "b"),
                                prior(student_t(3,0,5), class="sd")),
                        # sample_prior = "only",
                        chains=4, iter = 2000)

preds_sens = mod_prob_ab_sens$data %>% 
  distinct(location, river, target_name) %>% 
  expand_grid(quant_1000 = mean(clean_results$quant_1000), 
              sample_name = "new") %>% 
  add_epred_draws(mod_prob_ab_sens, allow_new_levels = T) 


preds_sens_graph <- preds_sens %>% 
  ggplot(aes(x = location, y = .epred, fill = river)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~target_name)+
  labs(y="Probability of a positive sample",
       title="Sensitivity Analysis Result") +  
  geom_point(data = mod_prob_ab_sens$data, 
             aes(y = quant_cat, color = river, 
                 group = interaction(river, location)), 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0),
             shape = 21, alpha = 0.2)

ggsave(preds_sens_graph, file = "plots/sensitivity_analysis.png", dpi = 750, 
       width = 7, height = 5, units = "in")

conditional_effects(mod_prob_ab_sens, conditions = tibble(target_name = "BHC"))

pp_check(mod_prob_ab_sens)


# what's next? 
# 1) Model for differences between field filtered and lab filtered
# 2) Model for pilot study samples: how many samples do you need to get a positive detection?
