#---------------------------------------------------------------------------------#
# list of data files created here:

#random_m: a tibble with:
  #id = unique identifier
  #mat = list-column of matrix objects (matrix$matrix_A, matrix_U, etc.)
  #life-history traits (lam, GenT, Fec, rupr, rlwr, xt, etc.)

#model_random: a model summarizing how life-history traits predict resilience metrics

#data_random: cleaned posterior samples, but not linked to individual IDs.
#----------------------------------------------------------------------------------#

library(Rage)
library(popdemo)
library(tidyverse)
library(tidybayes)
library(tidyr)
library(brms)
library(bayestestR)
library(cowplot)
library(patchwork)
library(Rcompadre)
library(doSNOW)
library(parallel)
library(popbio)
library(Matrix)
library(Rage)
library(furrr)
library(here)

#library(mpmtools)

#Working directories

path <- gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath <- paste0(path,"/Data")
ResultPath <- paste0(path, "/Results") 
CodePath <- paste0(path, "/Code") 

# Load functions to simulate the data 

source(paste0(CodePath, "/SimulationFunctions.R"))
source(paste0(CodePath, "/Functions.R"))

# Store the range of matrix sizes

dimensions <- range(c(2,59))

# Random matrices --------------------------------------------------------------
# Here we simulate matrix population models with random non-negative elements, 
# keeping the Umat column sums to <= 1 (=survival cannot be higher than 1). 

# We will use a for each loop 

clus  <- makeCluster(detectCores() - 1)

# register the cluster with doParallel

registerDoSNOW(clus)

# Build random matrices with stasis, progression and retrogression

random <- foreach(i = c(3:max(dimensions)),
                  .packages = c("tidyverse", "popbio", "popdemo", "Rage"),
                  .combine = "rbind") %dopar% {
                    matrices <- replicate(n=10, 
                                          random_matrices(dimension=i),
                                          simplify = F)
                   tibble(id = i, mat = matrices)
                    return(matrices)
                  }

# Reactivity analysis of matrix A to estimate the life history traits----------------
random_m <- foreach(i = seq_along(random),
                    .combine = "rbind",
                    .packages = c("tidyverse", "popdemo", "popbio", "Rage"))  %dopar% {
                      matrix <- random[[i]]
                      id <- i
                      
                      if(!is.null(matrix$matrix_A) &&
                         is.matrix(matrix$matrix_A) &&
                         nrow(matrix$matrix_A) == ncol(matrix$matrix_A) &&
                         isErgodic(matrix$matrix_A)) {
                        
          #use tryCatch() to catch errors
                        lam <- tryCatch(lambda(matrix$matrix_A), error = function(e) NA)
                        GenT <- tryCatch(generation.time(matrix$matrix_A), error = function(e) NA)
                        Fec <- tryCatch(vital_rates(matrix$matrix_U, matrix$matrix_F)$fec, error = function(e) NA)
                        
                        #This is the maximum amplification after a small pulse shock:
                        rupr <- tryCatch(reac(A = matrix$matrix_A, bound = "upper"), error = function(e) NA)
                        
                        #This is the minimum possible growth immediately after a small disturbance:
                        rlwr <- tryCatch(reac(A = matrix$matrix_A, bound = "lower"), error = function(e) NA)
                        
                        #measures how quickly transients decay. Even without specifying a  disturbance vector, return.time() uses norms of projection matrices to estimate how long a typical deviation takes to fade.
                        xt <- tryCatch(return.time(matrix$matrix_A), error = function(e) NA)
                        
                        Dimension <- tryCatch(nrow(matrix$matrix_A), error = function(e) NA)
                        
                      } else {
                        lam <- NA
                        GenT <- NA
                        Fec <- NA
                        rupr <- NA
                        rlwr <- NA
                        xt <- NA
                        Dimension <- NA
                      }
                      
                      tibble(id = id, mat = list(matrix), lam, GenT, Fec, rupr, rlwr, xt, Dimension)
                    }

#Stop the cluster 

stopCluster(clus)

random_m<- random_m %>% 
  filter(!is.infinite(GenT)) %>% 
  mutate(xt=log(xt+1),
         rupr=log(rupr+1),
         rlwr=abs(log(rlwr+1)),
         Fec=log(Fec+1),
         GenT=log(GenT+1))

setwd(DataPath)

#save(random_m, file="SimData1.RData")


# Population time series --------------------------------------------------------

# Use furrr for parallel efficiency
plan(multisession)

source(here("C:/Users/qt37rapo/Desktop/PostDoc/Demography/DemographicResilience/functions/PopulationDynamics_disturbance.R"))


# Simulate for each population (control + disturbance cases)
sim_all <- random_m %>%
  drop_na(rlwr, rupr, mat) %>%
  mutate(
    control = future_map(mat, ~ simulate_disturbed_pop(.x, case = "control")),
    disturbed = future_pmap(list(mat, rlwr, rupr), 
                            ~ simulate_disturbed_pop(..1, rlwr = ..2, rupr = ..3, 
                            case = "resilience"))
  )


# Extract time series and compute RR
rr_data <- sim_all %>%
  mutate(
    rr_series = map2(disturbed, control, ~ {
      left_join(.x, .y, by = c("time","stage" ), suffix = c("_treat", "_ctrl")) %>%
        mutate(response_ratio = (total_treat - total_ctrl) / (total_treat + total_ctrl))
    })
  ) %>%
  select(id,rupr,rlwr,xt,lam,Dimension, rr_series) %>%
  unnest(rr_series)%>%
  filter(time >40 )
  

dum = rr_data %>%filter(id==1) 
dum%>%
  ggplot(.)+
  geom_point(aes(x = time, y = response_ratio))+
  geom_line(aes(x = time, y = response_ratio))+
  facet_grid(~stage)#+
  geom_line(aes(x = time, y = total_treat), color = 'darkred')+
  geom_line(aes(x = time, y = total_ctrl), color = 'black')



# Building communities --------------------------------------------------------------

# Number of communities to simulate
dplyr::n_distinct(random_m$id) #check how many spp

n_communities <-500

# Range of species per community
min_species <- 10

# Set seed for reproducibility
set.seed(42)

# Simulate random communities
sim_communities <- tibble(
  comm_id = 1:n_communities,
  n_species = min_species,
  members = map(n_species, ~ sample_n(sim_all, size = .x)))


# Add species-level time series to each community
sim_communities_ts <- sim_communities %>%
  mutate(
    species_ts = map(members, ~ pmap_dfr(.x, function(id, control, disturbed, rlwr, rupr, xt, lam, Dimension, ...) {
      # Join pre-computed time series
      treat <- disturbed %>%
        mutate(treatment = "treat", species_id = id)
      ctrl <- control %>%
        mutate(treatment = "control", species_id = id)
      
      bind_rows(treat, ctrl) %>%
        mutate(
          rlwr = rlwr,
          rupr = rupr,
          xt = xt,
          lam = lam,
          Dimension = Dimension
        )
    }))
  )



# Combine and aggregate to community level
community_ts <- sim_communities_ts %>%
  select(comm_id, species_ts) %>%
  unnest(species_ts) %>%
  group_by(comm_id, treatment, time) %>%
  reframe(
    mean_lambda = mean(lam, na.rm = TRUE),
    mean_xt = mean(xt, na.rm = TRUE),
    mean_rupr = mean(rupr, na.rm = TRUE),
    mean_rlwr = mean(rlwr, na.rm = TRUE),
    sd_rlwr = sd(rlwr),
    community_abundance = sum(abundance, na.rm = TRUE)  ) %>%
  pivot_wider(names_from = treatment, values_from = community_abundance) %>%
  mutate(totalRR = (treat - control) / (treat + control))

#plots
community_ts%>%
  filter(mean_lambda <3)%>%
ggplot(., aes(x = time, y = totalRR, color = as.factor(mean_lambda)))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = 0)

community_ts%>%
  filter(mean_lambda >6.5)%>%
  ggplot(., aes(x = time, y = totalRR, color = mean_lambda))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = 0)


#Stability calculation ------------------------------------------------------------------------
comm_resistance <- community_ts%>% 
  filter(time == 50) %>%
  select(comm_id, totalRR)%>%
  rename(resist = totalRR)

comm_recov <- community_ts%>% 
  filter(time == 200) %>%
  select(comm_id, totalRR)%>%
  rename(recov = totalRR)

stability<- community_ts%>%
  group_by(comm_id, mean_xt,sd_rlwr, mean_lambda, mean_rupr, mean_rlwr)%>%
  reframe(OEV = MESS::auc(y=totalRR, x=time, from=min(time), to=max(time), absolutearea = F),
          absOEV = MESS::auc(y=totalRR, x=time, from=min(time), to=max(time), absolutearea = T),
          CV = mean(totalRR)/sd(totalRR)) %>%
  left_join(., comm_resistance)%>%
  left_join(., comm_recov)



stability %>%
  pivot_longer(c(mean_xt, mean_lambda, mean_rupr, mean_rlwr,sd_rlwr),names_to = 'demography', values_to='values')%>%
ggplot(., aes(x = values, y = recov))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_grid(~demography, scales = 'free_x')+  
  geom_hline(yintercept = 0)

ggsave(plot = last_plot(), 
  file = here("C:/Users/qt37rapo/Desktop/PostDoc/Demography/DemographicResilience/output/DemographicResilience_recov.png"),
  width = 10, height = 4)
