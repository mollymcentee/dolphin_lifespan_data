#############################################################################
#                                                                             
# Sex bias in mortality risk changes over the lifespan of bottlenose dolphins
#                                                
#                                                                           
# Molly H.F. McEntee, Ewa Krzyszczyk, Vivienne Foroughirad, Janet Mann    
#                                  
#                                                                           
# R-script                                                    
#                                                                           
#                                                              June 2023  
#############################################################################

############################################################
### Data and packages     
############################################################

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(viridis)
library(patchwork)
library(mgcv)
library(pammtools)
library(BaSTA)
library(lubridate)
library(gridExtra)

lifespan.data <-
  read.csv(("lifespan_data.csv"), row.names = 1)%>%
  mutate_at(vars(Birth.Date, Death.Date, last.sight.date, first.sight.date), as.Date)

sexed.sample <- lifespan.data %>%
  filter(Sex != "Unknown") %>%
  filter(time > 3) %>%
  mutate(trunc.time = ifelse(trunc.time > 3, trunc.time, 3)) %>%
  mutate(Sex = case_when(Sex == "Female" ~ 0,
                         Sex == "Male" ~ 1))

############################################################
### Create KM plots              
############################################################

theme_set(theme_classic())

##KM plot for whole sample
KMplotall <- ggsurvplot(
  survfit(Surv(trunc.time, time, event) ~ 1,
          data = lifespan.data),
  risk.table = FALSE,
  conf.int = FALSE,
  break.time.by = 5,
  xlim = c(0, 52),
  font.x = 18,
  font.y = 18,
  xlab = NULL,
  font.tickslab = 14,
  font.title = 18,
  palette = c(adjustcolor("grey5", alpha.f = .6)),
  size = .9,
  censor.size = 6,
  axes.offset = TRUE,
  legend = "none",
  ggtheme = theme(
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 10))
  )
)


##KM plot by sex
KMplotsexed <- ggsurvplot(
  survfit(Surv(trunc.time, time, event) ~ Sex,
          data = lifespan.data),
  palette = adjustcolor(c("goldenrod2", "#039f88", "#440154"),
                        alpha.f = .6),
  risk.table = FALSE,
  conf.int = FALSE,
  break.time.by = 5,
  xlim = c(0, 52),
  xlab = "Age (years)",
  font.x = 18,
  font.y = 18,
  font.legend = 14,
  font.tickslab = 14,
  font.title = 18,
  legend = "bottom",
  size = 1,
  censor.size = 6,
  axes.offset = TRUE,
  legend.labs = c("Female", "Male", "Unknown sex"),
  legend.title	= "",
  ggtheme = theme(
    legend.key.width = unit(2, "line"),
    legend.background = element_blank(),
    legend.box.background = NULL,
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
    
  )
)

#extract plots
a <- KMplotall$plot
b <- KMplotsexed$plot


kmplot <- a / b +
  plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size = 16))

kmplot

############################################################
### Cox models                                                
############################################################

### Cox proportional hazards model with sex as a fixed predictor
cox.ph.mod <-
  coxph(Surv(trunc.time, time, event) ~ Sex, data = sexed.sample)

cox.ph.mod

### Test ph assumption
zph <- cox.zph(cox.ph.mod)

zph

### Plot scaled Schoenfeld residuals
plot(zph[1], lwd = 2) +
  abline(0, 0, col = 1, lty = 3, lwd = 2) +
  abline(
    h = cox.ph.mod$coef[1],
    col = 3,
    lwd = 2,
    lty = 2
  )

### Cox model with linear interaction between time and sex

cox.sex.time.mod <-
  coxph(
    Surv(trunc.time, time, event) ~ Sex + tt(Sex),
    data = sexed.sample,
    tt = function(x, t, ...)
    {
      x * t
    }
  )

cox.sex.time.mod

############################################################
### PAM models and plots              
############################################################

#format full data set for PAM
ped.data.full.sample <-
  as_ped(
    formula = Surv(trunc.time, time, event) ~ 1,
    id = "Dolphin.ID",
    data = lifespan.data
  ) %>%
  mutate(Dolphin.ID = as.factor(Dolphin.ID))

###PAM baseline hazards model with  full sample
PAM.model.full.sample <- gam(
  formula = ped_status ~ s(tend),
  data = ped.data.full.sample,
  family = poisson(),
  offset = offset,
  method = "REML"
)

summary(PAM.model.full.sample)

#plot baseline hazard
new.data.full.sample <-   ped.data.full.sample %>% 
  make_newdata(tend=unique(tend)) %>% 
  add_hazard(PAM.model.full.sample, type = "link") 

###log hazard plot
full.sample.plot <- ggplot(new.data.full.sample, aes(x = tend)) +
  geom_stephazard(aes(y = hazard)) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2)  +
  xlab("Age (years)") +
  ylab("Log-hazard") +
  theme_classic() +
  theme(
    text = element_text(size = 18),
    legend.position = "bottom",
    legend.key.width = unit(2, "line"),
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )

full.sample.plot

###format sexed data for PAM models
ped.data.sexed.sample <-
  as_ped(
    formula = Surv(trunc.time, time, event) ~ Sex,
    id = "Dolphin.ID",
    data = sexed.sample
  ) %>%
  mutate(Sex = as.factor(Sex)) %>%
  mutate(Dolphin.ID = as.factor(Dolphin.ID))


###PAM model stratified by sex
pam.model.strata.by.sex <- gam(
  formula = ped_status ~ Sex + s(tend, by = Sex),
  data = ped.data.sexed.sample,
  family = poisson(),
  offset = offset,
  method = "REML"
)

summary(pam.model.strata.by.sex)

###Plot stratified baseline hazards by sex

new.data.sexed.sample <-   ped.data.sexed.sample %>% 
  make_newdata(tend=unique(tend), Sex=unique(Sex)) %>% 
  group_by(Sex) %>%
  add_hazard(pam.model.strata.by.sex, type = "link")  %>%
  mutate(Sex = ifelse(Sex == "0", "Female", "Male"))


###log hazard plot 
log.haz.plot.strat.by.sex <- ggplot(new.data.sexed.sample,
                                    aes(x = tend, col = Sex)) +
  geom_step(aes(y = hazard, col = Sex)) +
  scale_color_manual(values = c(
    adjustcolor("goldenrod2", alpha.f = .6),
    adjustcolor("#039f88", alpha.f = .6)
  )) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Sex),
                  alpha = 0.3,
                  linetype = 0) +
  scale_fill_manual(values = c(
    adjustcolor("goldenrod2", alpha.f = .6),
    adjustcolor("#039f88", alpha.f = .6)
  )) +
  xlab(NULL) +
  ylab("Log-hazard") +
  theme_classic() +
  scale_x_continuous(breaks = seq(5, 50, 5),  limits = c(3, 52)) +
  theme(
    text = element_text(size = 18),
    legend.position = "bottom",
    legend.key.width = unit(2, "line"),
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )

###Plot effect of sex with age
new.data.diff.over.time <-  ped.data.sexed.sample %>% 
  make_newdata(tend=unique(tend),  Sex = c(1)) %>%
  add_term(pam.model.strata.by.sex, term = "Sex", reference = list(Sex = 0))

sex.haz.diff.plot <-
  ggplot(new.data.diff.over.time, aes(x = tend, y = fit)) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper),
                  alpha = 0.3,
                  linetype = 0) +
  geom_step() +
  xlab("Age (years)") +
  ylab("Effect of sex") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 51, 5), limits = c(3, 52)) +
  theme(
    text = element_text(size = 18),
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  geom_hline(
    yintercept = 0,
    size = .8,
    linetype = "dashed",
    color = "red"
  )

pammplot <- log.haz.plot.strat.by.sex / sex.haz.diff.plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

pammplot

############################################################
### PAM sensitivity analysis         
############################################################


mylist <- list()

for (i in 1:100) {
  
  randomized.death.date <- lifespan.data %>%
    ###calculate difference between  assigned death date and last sighting date
    mutate(death.diff = (as.numeric(Death.Date - last.sight.date)) / 365.25)    %>%
    ###remove individuals who are still alive 
    filter(!is.na(death.diff)) %>%
    ###calculate a new death offset from a gamma distribution, setting the 
    #mean to the difference between the individual's age at death and their 
    #age at last sighting. This means that the distribution will have the 
    #estimated age at death as the mean value, the left side will be bounded 
    #by the date of last sight and the right side will have a long tail 
    #proportional to the delay between the last sighting and the new death
    #date, which is based on maximum sighting gap.
    rowwise() %>%
    mutate(random.death.offset = rgamma(1, shape = death.diff, scale = 1)) %>%
    #apply new offset to generate randomized death date
    mutate(random.death.date = as.Date(as.Date(last.sight.date) + dyears(random.death.offset))) %>%
    select(Dolphin.ID, random.death.date)
  
  randomized.birth.date <- lifespan.data %>%
    ###calculate estimated age at first sighting
    mutate(birth.diff = (as.numeric(first.sight.date - Birth.Date)) / 365.25)    %>%
    ###remove individuals who were never seen in surveys and focals
    filter(total.num.sightings != 0) %>%
    ###calculate new age at first sight from a gamma distribution with the
    #estimated age at first sight as the mean. The new birth date will be
    #bounded on the right by first sighting date and have a long tail to left
    #proportional to the estimated age at first sighting.
    rowwise() %>%
    mutate(random.birth.offset = rgamma(1, shape = birth.diff, scale = 1)) %>%
    mutate(random.birth.date = as.Date(as.Date(first.sight.date) - dyears(random.birth.offset))) %>%
    select(Dolphin.ID, random.birth.date)
  
  ###bring new birth and death dates into one df
  randomized.life.table <- lifespan.data %>%
    left_join(randomized.birth.date) %>%
    left_join(randomized.death.date) %>%
    
    ###remove individuals with no sightings and therefore no new birth.dates
    filter(!is.na(random.birth.date)) %>%
    
    ###calculate  truncated time - if an individual was born
    ###before 1985, then trunc.time is age in 1985
    mutate(random.trunc.time = if_else(
      random.birth.date > as.Date('1985-01-01'),
      0,
      as.numeric((
        as.Date('1985-01-01') - random.birth.date
      ) / 365.25)
    )) %>%
    
    ###calculate age at death and age at censor; setting age at censor to 1 Jan, 2020
    mutate(random.death.age = as.numeric(random.death.date - random.birth.date) / 365.25)  %>%
    mutate(random.censor.age = as.numeric(as.Date('2020-01-01') - random.birth.date) / 365.25) %>%
    ###Merge the death age and censor age into one category, time
    ###if an individual has a randomized death date after Jan 1, 2020, use their censor age
    mutate(
      random.time = case_when(
        random.censor.age < random.death.age ~ random.censor.age,
        random.censor.age >= random.death.age ~ random.death.age,
        is.na(random.death.age) ~ random.censor.age
      )
    ) %>%
    
    ###assign event status based on the same rules as time
    mutate(
      random.event =
        case_when(
          random.censor.age < random.death.age ~ 0,
          random.censor.age >= random.death.age ~ 1,
          is.na(random.death.age) ~ 0
        )
    ) %>%
    ###remove death and age at censor columns
    select(-random.death.age,-random.censor.age) %>%
    mutate(Dolphin.ID = as.factor(Dolphin.ID)) %>%
    
    ###select sample for juvenile/adult sexed analysis
    
    filter(random.time > 3) %>%
    filter(Sex != "Unknown") %>%
    ###reformat sex
    mutate(Sex = case_when(Sex == "Female" ~ 0,
                           Sex == "Male" ~ 1)) %>%
    
    ###update everyone's left truncation to minimum 3 years old
    mutate(random.trunc.time = ifelse(random.trunc.time < 3, 3, random.trunc.time))
  
  
  ###format sexed data for PAM models
  ped.data.sexed.sample.rand <-
    as_ped(
      formula = Surv(random.trunc.time, random.time, random.event) ~ Sex,
      id = "Dolphin.ID",
      data = randomized.life.table
    ) %>%
    mutate(Sex = as.factor(Sex)) %>%
    mutate(Dolphin.ID = as.factor(Dolphin.ID))
  
  ###PAM model stratified by sex
  pam.model.strata.by.sex.rand <- gam(
    formula = ped_status ~ Sex + s(tend, by = Sex),
    data = ped.data.sexed.sample.rand,
    family = poisson(),
    offset = offset,
    method = "REML"
  )
  
  new.data.diff.over.time.rand <-  ped.data.sexed.sample.rand %>% 
    make_newdata(tend=unique(tend), Sex = c(1)) %>%
    add_term(pam.model.strata.by.sex.rand, term = "Sex", reference = list(Sex = 0))
  
  max.tend <- max(new.data.diff.over.time.rand$tend)
  
  sex.haz.diff.plot.rand <-
    ggplot(new.data.diff.over.time.rand, aes(x = tend, y = fit)) +
    geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper),
                    alpha = 0.3,
                    linetype = 0) +
    geom_step() +
    xlab("Age (years)") +
    ylab("Effect of sex") +
    theme_classic() +
    scale_x_continuous(breaks = seq(0, max.tend, 5), limits = c(3, max.tend)) +
    theme(
      text = element_text(size = 10),
      axis.text = element_text(color = "black"),
      axis.title.y = element_text(margin = margin(r = 0)),
      axis.title.x = element_text(margin = margin(t = 0))
    ) +
    geom_hline(
      yintercept = 0,
      size = .8,
      linetype = "dashed",
      color = "red")
  
  mylist[[i]] <- sex.haz.diff.plot.rand
}


###set layout for plot
layout <- rbind(c(1,2,3),
                c(4,5,6),
                c(7,8,9),
                c(10,11,12))


p <- marrangeGrob(grobs=mylist, layout_matrix=layout)



############################################################
### BaSTA models           
############################################################

BaSTA.data <- read.csv("BaSTA_data.csv", row.names = 1)


###separate the sexes into two dataframes
BaSTA.data.females <- BaSTA.data %>%
  filter(SexFemale == 1) %>%
  select(-SexFemale, -SexMale)

BaSTA.data.males <- BaSTA.data %>%
  filter(SexMale == 1) %>%
  select(-SexFemale, -SexMale)

###run multibasta models
basta.model <- multibasta(BaSTA.data, studyStart = 1985, studyEnd = 2019,
                          minAge = 3,
                          recaptTrans = c(1985, 1996),
                          niter = 500000,
                          burnin = 50000,
                          thinning = 2000,
                          nsim = 12,
                          parallel = TRUE,
                          ncpus = 4)


basta.model.female <- multibasta(BaSTA.data.females, studyStart = 1985, studyEnd = 2019,
                                 minAge = 3,
                                 recaptTrans = c(1985, 1996),
                                 niter = 500000,
                                 burnin = 50000,
                                 thinning = 2000,
                                 nsim = 12,
                                 parallel = TRUE,
                                 ncpus = 4)


basta.model.male <- multibasta(BaSTA.data.males, studyStart = 1985, studyEnd = 2019,
                               minAge = 3,
                               recaptTrans = c(1985, 1996),
                               niter = 500000,
                               burnin = 50000,
                               thinning = 2000,
                               nsim = 12,
                               parallel = TRUE,
                               ncpus = 4)

###Extract simple Gompertz model
go.si <- basta.model[["runs"]][["Go.Si"]]

summary(go.si)

plot(go.si)

plot(go.si, plot.trace = FALSE, fancy = TRUE,
     col = c(
       adjustcolor("goldenrod2", alpha.f = 1),
       adjustcolor("#039f88", alpha.f = 1)),
     names.legend = c("Females", "Males"))
