library(rsq)
library(lme4)
library(vegan)
library(fitdistrplus)
library(dplyr)
library(car)
library(MASS)
library(interactions)
library(DHARMa)
library(MuMIn)
library(effectsize)
library(ggplot2)
library(jtools)
library(performance)
library(see)
library(qqplotr)
library(blmeco)
library(sjPlot)
library(bbmle) 
library(optimx)
library(nloptr)
library(dfoptim)
library(tidyr)
library(emmeans)
library(magrittr)
library(ggeffects)
library(viridis)
library(ggforce)
library(ggdist)
library(gghalves)
library(cowplot)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(broom.mixed)
library(jtools)
library(MASS)
library(stargazer)

X <- read.table("DATE.txt", h = T, sep="\t")

#Testing the premise whether the number of inflorescence in the 
#blocks could influence the observed visits

#Normality test
shapiro.test(X$N.VISITS)
shapiro.test(X$N.OF.INFLORESCENCE)

#Model
inflores <- glm.nb(N.VISITS~N.OF.INFLORESCENCE, data = X)

#Investigating the behavior of model residues
residuos.inflor <- simulateResiduals(fittedModel = inflores., n=1000)
plot(residuos.inflor) #OK

#Results
summary(inflores)
Anova(inflores)


##ANALYSIS OF THE WORK HYPOTHESES:
#Planned contrasts
contrasts(X$treatment)

#First line = Inflorescences without ants vs. any insect on the inflorescences (Any organism)
#Second line = inflorescences with pincushioned ants vs. with pinned beetle (Ant presence)
#Third line = inflorescence with live ants vs. with pinned ant (Ant activity)
mat.con.2<-as.matrix(cbind(c(-1,-1,-1,3),
                           c(0,1,-1,0),
                           c(-1,0,1,0)))

###LOOPS
response_vars <- c("Dipterans", "Wasps", "Bees", "Coleoptera", "Cockroaches", "Others.organism", "All.visitors")
# Create an empty list to store the models
model_list <- list()
model_list.TMB <- list()
Anova_list<-list()
effectsize.list<-list()
effectsize.list.vec<-list()
effectsize.list.ci<-list()

# Perform GLMs for each response variable and store in separate objects
for (response_var in response_vars) {
  formula <- as.formula(paste(response_var, " ~ treatment*+(1|blocks)"))
  glm_model <- glmer(formula, data = X, family = poisson,
                     contrasts=list(tratamentos=mat.con.2), 
                     glmerControl(optCtrl = list(maxfun = 1000000)))
  # Assign the model to a separate object
  assign(paste("Model_", response_var, sep = ""), glm_model)
  # Store the model in the list
  model_list[[response_var]] <- glm_model
  Anova_table<-Anova(model_list[[response_var]])
  # Assign the model to a separate object
  assign(paste("Anova_", response_var, sep = ""), Anova_table)
  # Store the Anova in the list
  Anova_list[[response_var]] <- Anova_table
  efectsize_table<-effectsize::effectsize(model_list[[response_var]],method = "refit", two_sd=F, 
                                          exponentiate=T)
  # Assign the model to a separate object
  assign(paste("Effect_", response_var, sep = ""), efectsize_table)
  # Store the Anova in the list
  effectsize.list[[response_var]] <- efectsize_table
  effectsize.list.vec[[response_var]]<-efectsize_table[,2]
  names(effectsize.list.vec[[response_var]])<-efectsize_table$Parameter
  effectsize.list.ci[[response_var]]<-matrix(c(efectsize_table[,4], efectsize_table[,5]), nrow = 4, ncol = 2)
  rownames(effectsize.list.ci[[response_var]])<-efectsize_table$Parameter
}

residuosdip <- simulateResiduals(fittedModel = glm_model, n=1000)
plot(residuosdip) #OK

# Access and inspect the GLM models
for (response_var in response_vars) {
  glm_model <- get(paste("Model_", response_var, sep = ""))
  print(paste("GLM model for", response_var))
  print(summary(glm_model))
  cat("\n")
}

model_list
Anova_list

#To generate tables with model summary
stargazer(Model_Dipterans,  Model_Wasps, Model_Bees, Model_Coleoptera, Model_Cockroaches, Model_Others.organism, Model_All.visitors,
          type = "text",style= "all",ci=T, 
          ci.level=0.95,out = "Results.models.html", single.row=TRUE, p.auto=F, digits=3,
          model.numbers=F, coef=effectsize.list.vec, ci.custom=effectsize.list.ci,report= "vcsp*" )

#Prediction graph:
Prediction_lists.mos<- ggpredict(model_list[[1]], type = "fixed", terms = c("treatment"))
Prediction_lists.mos

model_list
response_vars <- c("Dipterans", "Wasps", "Bees", "Coleoptera", "Cockroaches", "Others.organism", "All.visitors")
response_vars.list<-list("Dipterans"="Dipterans", "Wasps"= "Wasps", "Bees"="Bees", 
                         "Coleoptera"="Coleoptera", "Cockroaches"="Cockroaches", "Others.organism"="Others.organism", 
                         "All.visitors"="All.visitors")

Theme.coef.2<- theme(axis.text.y = element_text(size=8, angle=0,face="bold"),
                     legend.text = element_text(size=10),plot.tag = element_text(size=12, face="bold"),
                     axis.text.x =element_text(size=10),plot.title = element_text(size=12),
                     legend.position = "none", axis.line.y = element_line(linewidth = 1, colour = 1),
                     axis.line.x = element_line(linewidth = 1, colour = 1), 
                     axis.title.x = element_blank(),axis.title.y = element_blank(), 
                     #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     plot.margin = margin(0,-0.2,0,0,unit="cm"))

Dip.plots.loop<-list()
Dip.plots.loop.TMB<-list()

for (response_var in response_vars) {
  Dip.plots.loop[[response_var]]<-  plot_summs(model_list[[response_var]], scale = T, exp = F, 
                                               colors="Qual3", omit.coefs= c("(Intercept)", "Intercept"))+
    xlab("Standarized coefficient")+
    xlim(-0.6, 0.6)+
    scale_y_discrete(labels=c("treatment1" = "Manipulation", "treatment2" = "Ants/Other organism", "treatment3" = "Mobility"))+
    #scale_color_manual(values=c("NaFor"="#F0E442", "Crop"="#0072B2", "TRATMO"="#D55E00"))+
    ggtitle(response_vars.list[[response_var]])+theme_minimal()+
    Theme.coef.2+
    theme(axis.title.x = element_blank())
}

plot.grid2<- cowplot::plot_grid(plotlist=Dip.plots.loop, ncol = 2, nrow = 4,labels="auto",
                                hjust=c(0,0), vjust=1, scale = 0.9, rel_widths=c(1,1), axis="r")
plot.grid2
