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

######### data ########

eDNA_compiled_results <- read_csv("data/average_sample_results.csv") %>% 
  clean_names()

clean_results <- eDNA_compiled_results %>% replace_na(list(ct_mean=0, ct_sd=0,quantity_mean=0,quantity_sd=0)) %>% 
  ungroup %>% 
  mutate(quant_cat=case_when(quantity_mean==0~0, TRUE~1),
         above_below=case_when(location=="above"~0,
                               location=="below"~1),
         river_cat=case_when(river=="big sioux"~0,
                             river=="vermillion"~1),
         quant_1000=quantity_mean/1000,
         quant_mean_s=scale(quantity_mean),
         target_color=case_when(target_name=="SVC"~"deeppink1",
                                target_name=="BHC"~"darkslategray3"),
         filter_color=case_when(filter_method=="field"~"lightgrey",
                                 filter_method=="lab"~"gray40"))

pilot_averages <- read_csv("data/pilot_averages.csv") %>% 
  replace_na(list(ct_mean=0, ct_sd=0,quantity_mean=0,quantity_sd=0))

##### modeling probability of a positive sample above and below the barriers #####

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
bayes_R2(mod_prob_ab)

ggplot(data=clean_results,aes(x=quant_1000, y=quant_cat))+
  geom_point()+
  facet_wrap(location~river)

clean_results %>% 
  group_by(location, river, target_name) %>% 
  summarize(positive = sum(quant_cat))

clean_results %>% filter(location == "above") %>% 
  filter(river == "vermillion")

# clean_results %>% filter(is.na(river)) %>% view

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


# Sensitivity analysis; doubled the standard deviation 

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
bayes_R2(mod_prob_ab_sens)


####### Field Filtered vs. Lab Filtered ######

get_prior(quant_cat~filter_method*target_name + (1|sample_name),
          data=clean_results,
          family=gaussian())

mod_prob_filter <- brm(quant_cat~filter_method*target_name + (1|sample_name),
                   data=clean_results,
                   family=gaussian(),
                   # prior=c(prior(normal(-10,1), class= "Intercept"),
                   #         prior(normal(0,1), class= "b")),
                   # sample_prior = "only",
                   chains=1, iter = 2000)

preds_filter = mod_prob_filter$data %>% 
  distinct(target_name, filter_method) %>% 
  expand_grid(quant_cat = mean(clean_results$quant_cat), 
              sample_name = "new") %>% 
  add_epred_draws(mod_prob_filter, allow_new_levels = T) 


preds_filter_graph <- preds_filter %>% 
  ggplot(aes(x = filter_method, y = .epred, fill = target_name)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~target_name)+
  labs(y="Probability of detection") +  
  geom_point(data = mod_prob_filter$data, 
             aes(y = quant_cat, color = filter_method),
                 # group = interaction(filter_method, target_name)), 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0),
             shape = 21, alpha = 0.2)

ggsave(preds_filter_graph, file = "plots/positive_predictions.png", dpi = 750, 
       width = 7, height = 5, units = "in")

conditional_effects(mod_prob_filter, conditions = tibble(target_name = "BHC"))

pp_check(mod_prob_filter)
bayes_R2(mod_prob_filter)


######## pilot study model #########

# 2 separate models: 
# first question: how many samples needed to get positive detection
# second question: is there a difference in quantity between the two protocols
# (kristie's and qiagen)


get_prior(quantity_mean~,
          data=pilot_averages,
          family=bernoulli())

mod_prob_pilot <- brm(quantity_mean~,
                   data=pilot_averages,
                   family=bernoulli(),
                   # prior=c(prior(normal(-10,1), class= "Intercept"),
                   #         prior(normal(0,1), class= "b")),
                   # sample_prior = "only",
                   chains=1, iter = 2000)


preds_pilot = mod_prob_pilot$data %>% 
  distinct() %>% 
  expand_grid(quant_mean = mean(pilot_averages$quant_1000), 
              sample_name = "new") %>% 
  add_epred_draws(mod_prob_ab, allow_new_levels = T) 


preds_pilot_graph <- preds_pilot %>% 
  ggplot(aes(x = , y = .epred, fill = river)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~target_name)+
  labs(y="Probability of a positive sample") +  
  geom_point(data = mod_prob_pilot$data, 
             aes(y = quant_cat, color = river, 
                 group = interaction(river, location)), 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0),
             shape = 21, alpha = 0.2)

ggsave(preds_pilot_graph, file = "plots/positive_predictions.png", dpi = 750, 
       width = 7, height = 5, units = "in")

conditional_effects(mod_prob_pilot, conditions = tibble(target_name = "BHC"))

pp_check(mod_prob_pilot)
bayes_R2(mod_prob_pilot)

