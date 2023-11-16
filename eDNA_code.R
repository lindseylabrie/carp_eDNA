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
library(scales)
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
  replace_na(list(ct_mean=0, ct_sd=0,quantity_mean=0,quantity_sd=0)) %>% 
  mutate(quant_cat=case_when(quantity_mean==0~0, TRUE~1))

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

# probabilities of differences between above and below samples in each river
preds %>% 
  group_by(river, location) %>% 
  median_qi(.epred)

# probability of below spillway samples always being higher than above spillway samples
preds %>% 
  ungroup() %>% 
  select(-.row, -.chain, -.iteration) %>% 
  pivot_wider(names_from = location, values_from = .epred) %>% 
  mutate(diff = below - above) %>% 
  summarize(total_greater = sum(diff>0)/nrow(.))
  median_qi(diff) 

preds_graph <- preds %>% 
  ggplot(aes(x = location, y = .epred, fill = river)) +
  facet_wrap(~target_name)+
  labs(y="Probability of a positive sample") +  
  geom_point(data = mod_prob_ab$data, 
             aes(y = quant_cat, color = river, 
                 group = interaction(river, location)), 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0),
             shape = 21, alpha = 0.2)+ 
  geom_boxplot(outlier.shape = NA)

ggsave(preds_graph, file = "plots/positive_predictions.png", dpi = 750, 
       width = 7, height = 5, units = "in")

ggsave(preds_graph, file = "plots/positive_predictions_compare.png", dpi = 750, 
       width = 3.5, height = 3, units = "in")

conditional_effects(mod_prob_ab, conditions = tibble(target_name = "BHC"))

pp_check(mod_prob_ab)
bayes_R2(mod_prob_ab)
summary(mod_prob_ab)

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
  ggplot(aes(x = location, y = .epred, fill = river))  +
  facet_wrap(~target_name)+
  labs(y="Probability of a positive sample") +  
  geom_point(data = mod_prob_ab_sens$data, 
             aes(y = quant_cat, color = river, 
                 group = interaction(river, location)), 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0),
             shape = 21, alpha = 0.2)+
  geom_boxplot(outlier.shape = NA)

ggsave(preds_sens_graph, file = "plots/sensitivity_analysis.png", dpi = 750, 
       width = 3.5, height = 3, units = "in")

conditional_effects(mod_prob_ab_sens, conditions = tibble(target_name = "BHC"))

pp_check(mod_prob_ab_sens)
bayes_R2(mod_prob_ab_sens)
summary(mod_prob_ab_sens)


####### Field Filtered vs. Lab Filtered ######

get_prior(quant_cat~filter_method*target_name + (1|sample_name),
          data=clean_results,
          family=bernoulli())

mod_prob_filter <- brm(quant_cat~filter_method*target_name + (1|sample_name),
                   data=clean_results,
                   family=bernoulli(),
                   # prior=c(prior(normal(-10,1), class= "Intercept"),
                   #         prior(normal(0,1), class= "b")),
                   # sample_prior = "only",
                   chains=4, iter = 2000)

preds_filter = mod_prob_filter$data %>% 
  distinct(target_name, filter_method) %>% 
  expand_grid(quant_cat = mean(clean_results$quant_cat), 
              sample_name = "new") %>% 
  add_epred_draws(mod_prob_filter, allow_new_levels = T) 

# probabilities of differences between lab and field filtered
preds_filter %>% 
  group_by(filter_method) %>% 
  median_qi(.epred)

# probability of lab vs field filtered
preds_filter %>% 
  ungroup() %>% 
  select(-.row, -.chain, -.iteration) %>% 
  pivot_wider(names_from = filter_method, values_from = .epred) %>% 
  mutate(diff = field - lab) %>% 
  # summarize(total_greater = sum(diff>0)/nrow(.))
  median_qi(diff)

# graph on a non log scale
preds_filter_graph_reg <- preds_filter %>% 
  ggplot(aes(x = filter_method, y = .epred, fill = target_name)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~target_name)+
  labs(y="Probability of detection",
       x="Filter method") +
  geom_point(data = mod_prob_filter$data, 
             aes(y = quant_cat),
             position = position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0),
             shape = 21, alpha = 0.5)+
  theme(legend.position="none")
  # scale_y_log10(breaks=c(0,0.0001,1), labels = function(x) format(x, scientific = FALSE))

ggsave(preds_filter_graph_reg, file = "plots/filter_location.png", dpi = 750, 
       width = 7, height = 5, units = "in")

# graph on a log 10 scale
preds_filter_graph_log10 <- preds_filter %>% 
  ggplot(aes(x = filter_method, y = .epred, fill = target_name)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~target_name)+
  labs(y="Probability of detection",
       x="Filter method") +
  geom_point(data = mod_prob_filter$data, 
             aes(y = quant_cat),
             position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0),
             shape = 21, alpha = 0.5)+
  theme(legend.position="none")+
  scale_y_log10(breaks=c(0,0.0001,0.01, 1), labels = function(x) format(x, scientific = FALSE))

ggsave(preds_filter_graph_log10, file = "plots/filter_location_log10.png", dpi = 750, 
       width = 6, height = 4.5, units = "in")


conditional_effects(mod_prob_filter, conditions = tibble(target_name = "BHC"))

pp_check(mod_prob_filter)
bayes_R2(mod_prob_filter)
summary(mod_prob_filter)

# no need for sensitivity analysis here -Jeff

######## pilot study #########

# first question: how many samples needed to get positive detection
# second question: is there a difference in quantity between the two protocols
# (kristie's and qiagen)

### How many samples? ###

pilot_averages %>%
  ungroup %>% 
  group_by(river,target_name) %>%
  sample_n(1) %>% 
  distinct(quant_cat) 

# at 1 sample for each location with known carp, 4/40 samples are false negatives
# i.e. a 10% chance that any sample is a false positive. This number halves as the
# number of samples increases.

### Difference in protocols? ###

Protocol_Plot <- ggplot(data=pilot_averages,
       aes(y=quantity_mean, x=protocol))+
  geom_point()+
  geom_boxplot(aes(group=protocol))+
  facet_wrap(~target_name)+
  ylab("Mean eDNA quantity (copies/μL)")+
  xlab("Protocol type")+
  scale_y_log10(labels = function(x) format(x, scientific = FALSE))

ggsave(Protocol_Plot, file = "plots/Protocol_Plot.png", dpi = 750, 
       width = 7, height = 5, units = "in")


# probably don't need to run a model based on this. Can use it as a justification for 
# using Kristie's protocol during extraction. Qiagen and kristie's get large amounts of DNA 
# as evidenced by the BHC model


###### other sanity checks to consider ######

# volume of water filtered and quantity of DNA

Volume_Plot_pilot <- ggplot(data=pilot_averages,
       aes(y=quantity_mean, x=volume_filtered_ml))+
  geom_point()+
  geom_boxplot(aes(group=volume_filtered_ml))+
  facet_wrap(~protocol)+
  ylab("Mean eDNA quantity (in copies/μL)")+
  xlab("Volume filtered (in mL)")

max(clean_results$volume_filtered_ml)

ggsave(Volume_Plot_pilot, file = "plots/Volume_Plot_pilot.png", dpi = 750, 
       width = 7, height = 5, units = "in")

Volume_Plot_data <- ggplot(data=clean_results,
                      aes(y=quantity_mean, x=volume_filtered_ml))+
  geom_point(aes(color=location))+
  # geom_boxplot(aes(group=location, fill=location))+
  facet_wrap(~target_name)+
  ylab("Mean eDNA quantity (copies/μL)")+
  xlab("Volume filtered (mL)")+
  scale_y_log10(labels = label_comma())

ggsave(Volume_Plot_data, file = "plots/Volume_Plot_data.png", dpi = 750, 
       width = 7, height = 5, units = "in")

# time between sampling and filtering affecting the amount of DNA found in the sample

Time_Plot <- ggplot(data=clean_results,
       aes(y=quantity_mean, x=time_from_sample_to_filter_hm))+
  geom_point(aes(color=location))+
  # geom_boxplot(aes(group=location,fill=location))+
  facet_wrap(~target_name)+
  ylab("Mean eDNA quantity (copies/μL)")+
  xlab("Time between sample and filter (hours)")+
  scale_y_log10(labels = label_comma())+ 
  theme(panel.spacing.x = unit(2, "lines"))
    

ggsave(Time_Plot, file = "plots/Time_Plot.png", dpi = 750, 
       width = 7, height = 5, units = "in")


# days between filtering and extraction affecting the amount of DNA found in the sample

Filter_to_Extraction <- ggplot(data=clean_results,
                    aes(y=quantity_mean, x=days_from_filter_to_extract))+
  geom_point(aes(color=location))+
  # geom_boxplot(aes(group=location,fill=location))+
  facet_wrap(~target_name)+
  ylab("Mean eDNA quantity (copies/μL)")+
  xlab("Days between filtering and extraction")+
  scale_y_log10(labels = label_comma())+ 
  theme(panel.spacing.x = unit(2, "lines"))

ggsave(Filter_to_Extraction, file = "plots/Filter_to_Extraction.png", dpi = 750, 
       width = 5, height = 3, units = "in")

# days between extraction and qPCR affecting the amount of DNA found in the sample

Extraction_to_qPCR <- ggplot(data=clean_results,
       aes(y=quantity_mean, x=days_from_extract_to_qpcr))+
  geom_point(aes(color=location))+
  # geom_boxplot(aes(group=location,fill=location))+
  facet_wrap(~target_name)+
  ylab("Mean eDNA quantity (copies/μL)")+
  xlab("Days between extraction and qPCR")+
  scale_y_log10(labels = label_comma())

ggsave(Extraction_to_qPCR, file = "plots/Extraction_to_qPCR.png", dpi = 750, 
       width = 5, height = 3, units = "in")

# does Month play a role?

level_order <- c('may', 'july', 'august','september','october') 

Month_Plot <- ggplot(data=clean_results,
       aes(y=quantity_mean, x=month))+
  geom_boxplot(aes(group=month),outlier.shape = NA)+
  geom_point(aes(color=location),
             shape=16,
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0.1))+
  facet_wrap(~target_name)+
  ylab("Mean eDNA quantity (copies/μL)")+
  xlab("Sample month")+
  scale_y_log10(labels = label_comma())+
  scale_x_discrete(limits=level_order)


ggsave(Month_Plot, file = "plots/Month_Plot.png", dpi = 750, 
       width = 7, height = 5, units = "in")



##### general statistics ######

clean_results %>% 
  group_by(location, river, target_name,year) %>% 
  summarize(positive = sum(quant_cat),
            average = mean(quantity_mean),
            sdmean=mean(quantity_sd))



