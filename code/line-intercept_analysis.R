#####Understory vegetation data: Line-intercept analysis#####

##Goal: Calculate and graph group-wise comparisons of veg conditions before and after treatment 
##DATE: 03 Aug 2023
##Author: Dr. Marc Mayes, SIG-NAL/USU/UCSB

#R packages
library(tidyverse)
library(sf)
library(lubridate)
library(nlme)
library(car)
library(rcompanion)
library(lsmeans)
library(emmeans)
library(multcompView)
library(multcomp)

#####1. Load line-intercept data#####

#data input directory:
indir <- "./output/"

#output directory
outdir <- "./output/line_int_analysis/"

#line-intercept data file:
line_intercept_infile <- "transect_height_weightedavg_allrecords_labels.csv"

tsect_lin <- read.csv(paste0(indir, line_intercept_infile))
head(tsect_lin)

#factor level: control timepoints
tsect_lin$timepoint <- factor(tsect_lin$timepoint,
                                                    levels=c("pre-grazing", "post-1mo", "post-2mo"))

#####2. Basic line-incercept dataframe orientation#####

head(tsect_lin)
#Transects each have grouping by time-point, plot type in observational experiment (grazed vs. control), overstory type, understory type and weighted-average height of all vegetation.

#####3. Summarize, graph, analysis:  Height differences of grazed vs. control, pre- and post-grazing#####

#Summarize by group
tsects_graze_v_control_alltimes <- tsect_lin %>%
  group_by(timepoint, plot_type) %>%
  summarize(N = length(ht_wtavg_m),
            ht_m = mean(ht_wtavg_m),
            ht_stdv = sd(ht_wtavg_m),
            ht_se = ht_stdv/sqrt(N))

#factor level: control timepoints
tsects_graze_v_control_alltimes$timepoint <- factor(tsects_graze_v_control_alltimes$timepoint,
                                                    levels=c("pre-grazing", "post-1mo", "post-2mo"))

#graph: boxplots of heights:
tsects_htcomp_graze_v_control <- ggplot(tsect_lin,
                                        aes(x=timepoint, y=ht_wtavg_m))+
  geom_boxplot(aes(fill=plot_type))+
  scale_fill_manual(values=c("wheat1", "wheat4"), name="Plot Type")+
  theme_bw()+
  labs(x="Time", y="Understory Height - m")+
  theme(axis.text=element_text(color="black"),
        axis.text.x=element_text(vjust=0.5, size=24),
        axis.text.y=element_text(size=20),
        axis.title=element_text(color="black", size=28),
        axis.title.x=element_text(margin=margin(t=20,r=0,b=0,l=0)),
        axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0)),
        legend.title=element_text(color="black", size=28),
        legend.text=element_text(color="black", size=24),
        #legend.position=c(0.9, 0.85),
        axis.line=element_line(color="black"),
        panel.border=element_rect(color="black")
        
  )
tsects_htcomp_graze_v_control
#export this graph at minimum:
png(filename=(paste0(outdir, "grazing_understory_impacts_v1.png")),
    width=10, height=6.75, units="in", res=600)
print(tsects_htcomp_graze_v_control)
dev.off()

#trial ANOVA: effect of timepoint and plot_type on understory height
names(tsect_lin)
graze_mod1 <- lm(ht_wtavg_m ~ plot_type * timepoint + tsect,
                 data=tsect_lin)
summary(graze_mod1)
anova(graze_mod1)

#####4. Repeated measures ANOVA: rcompanion version: https://rcompanion.org/handbook/I_09.html #####

#recode timepoint as integer.
tsect_lin$time <- as.integer(dplyr::recode(tsect_lin$timepoint,
                                "pre-grazing" = "1",
                                "post-1mo" = "2",
                                "post-2mo" = "3"))

#assess autocorrelation for lags in the time variable:
graze_mod2.a <- gls(ht_wtavg_m ~ plot_type + time + plot_type*time,
                    data=tsect_lin)
ACF(graze_mod2.a,
    form = ~ time | tsect)
#time-transect first-order correlation is 0.307, not accounting for random effects of transect.

graze_mod2.b <- lme(ht_wtavg_m ~ plot_type + time + plot_type*time,
                    random = ~ 1|tsect,
                    data=tsect_lin)
ACF(graze_mod2.b)
#time-transect first-order correlation is -0.555.  Use this value in the below..


#Mixed effects model assessing grazing effects on veg height, taking into account repeated measures and random effects.

graze_mod2 <- lme(ht_wtavg_m ~ plot_type + time + plot_type*time,
                  random = ~1|tsect,
                  correlation = corAR1(form = ~time | tsect,
                                       value=-0.555),
                  data=tsect_lin,
                  method="REML")
Anova(graze_mod2)

graze_mod2_v2 <- lme(ht_wtavg_m ~ plot_type + timepoint + plot_type*timepoint,
                     random = ~1|tsect,
                     correlation = corAR1(form = ~time | tsect,
                                          value=-0.555),
                     data=tsect_lin,
                     method="REML")
Anova(graze_mod2_v2)
##ANOVA result: time vs. grazing effects: grass height after grazing significantly different than pre-grazing, but effects differ by plot.  Between 

graze_mod2.fixed <- gls(ht_wtavg_m ~ plot_type + time + plot_type*time,
                        data=tsect_lin,
                        method="REML")
anova(graze_mod2,
      graze_mod2.fixed)
#No significant difference in models with and without random effects; random effects not significant
# Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# graze_mod2           1  7 17.98768 25.93614 -1.993843                        
# graze_mod2.fixed     2  5 19.57199 25.24946 -4.785995 1 vs 2 5.584305  0.0613


#calculate p-value and pseudo-r2 for model, using model with random effects
graze_mod2.null = lme(ht_wtavg_m ~ 1,
                      random = ~1|tsect,
                      data=tsect_lin)
nagelkerke(graze_mod2,
           graze_mod2.null)

#calculate p-value and pseudo-r2 for model, using null model without random effects
graze_mod2.null.2 <- gls(ht_wtavg_m~1,
                         data=tsect_lin)
nagelkerke(graze_mod2,
           graze_mod2.null.2)
##take the Cox and Snell (ML) pseudo-r2 for model (0.527) and p-value for model p = 0.001.


###post-hoc analysis by groups:
names(tsect_lin) #use emmeans
graze_mod2
#https://benwhalley.github.io/just-enough-r/contrasts-lmer.html
?emmeans
m.emm <- emmeans(graze_mod2_v2, ~plot_type*timepoint)
pairs(m.emm)
#this list shows pair-wise contrasts: Tukey tests

#targeted means comparisons: time-point comparisons across plot types.
m.emm2 <- emmeans(graze_mod2_v2, ~timepoint|plot_type)
pairs(m.emm2)
avg_pregraze_ht_all_plots <- sum(0.548, 0.719)/2
avg_pregraze_ht_all_plots
