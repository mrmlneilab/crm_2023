

library(survival)
library(survminer)

fit <- survfit(Surv(time, os) ~ Cluster, data = a4)
# fit <- survfit(Surv(time, os) ~ Cluster, data = subset(a4, Grade != "benign"))

survdiff(Surv(time, os) ~ Cluster, data = a4)
# survdiff(Surv(time, os) ~ Cluster, data = subset(a4, Grade != "benign"))

ggsurvplot(fit,
           pval = TRUE, 
           pval.size = 5,
           conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw()) # Change ggplot2 theme






