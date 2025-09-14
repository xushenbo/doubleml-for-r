library(ggplot2)
library(ggpubr)
library(pammtools)
library(DescTools)
library(data.table)
library(dplyr)
library(SuperLearner)
library(survival)

followup <- 6 # years
tau <- 6

setwd("")

load("cancer_incidence_20231228_newcode_newdata.RData")
cancer.incidence.cif <- cancer.incidence

# load("cancer_incidence_20231226_newcode_newdata_sep.RData")
# cancer.incidence.cif <- cancer.incidence
# 
# load("cancer_incidence_20231226_newcode_newdata_all.RData")
# cancer.incidence.cif$controlled.direct.cifj.driptw <- cancer.incidence$controlled.direct.cifj.driptw
# cancer.incidence.cif$composite.cifj.driptw <- cancer.incidence$composite.cifj.driptw
# cancer.incidence.cif$controlled.direct.cifj.drow <- cancer.incidence$controlled.direct.cifj.drow
# cancer.incidence.cif$composite.cifj.drow <- cancer.incidence$composite.cifj.drow

# add bound
bound.var <- c(1, grep("se", names(cancer.incidence.cif$sep.direct.astar1.cifj.driptw)), grep("diff", names(cancer.incidence.cif$sep.direct.astar1.cifj.driptw)))
cancer.incidence.cif$controlled.direct.cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$controlled.direct.cifj.driptw[, -bound.var])
cancer.incidence.cif$composite.cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$composite.cifj.driptw[, -bound.var])
cancer.incidence.cif$total.cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$total.cifj.driptw[, -bound.var])
cancer.incidence.cif$sep.direct.astar1.cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$sep.direct.astar1.cifj.driptw[, -bound.var])
cancer.incidence.cif$sep.direct.astar0.cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$sep.direct.astar0.cifj.driptw[, -bound.var])
cancer.incidence.cif$sep.indirect.astar1.cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$sep.indirect.astar1.cifj.driptw[, -bound.var])
cancer.incidence.cif$sep.indirect.astar0.cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$sep.indirect.astar0.cifj.driptw[, -bound.var])

cancer.incidence.cif$controlled.direct.cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$controlled.direct.cifj.drow[, -bound.var])
cancer.incidence.cif$composite.cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$composite.cifj.drow[, -bound.var])
cancer.incidence.cif$total.cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$total.cifj.drow[, -bound.var])
cancer.incidence.cif$sep.direct.astar1.cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$sep.direct.astar1.cifj.drow[, -bound.var])
cancer.incidence.cif$sep.direct.astar0.cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$sep.direct.astar0.cifj.drow[, -bound.var])
cancer.incidence.cif$sep.indirect.astar1.cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$sep.indirect.astar1.cifj.drow[, -bound.var])
cancer.incidence.cif$sep.indirect.astar0.cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.cif$sep.indirect.astar0.cifj.drow[, -bound.var])

cancer.incidence.cif$controlled.direct.cifj.driptw$t <- cancer.incidence.cif$controlled.direct.cifj.driptw$t/365.25
cancer.incidence.cif$composite.cifj.driptw$t <- cancer.incidence.cif$composite.cifj.driptw$t/365.25
cancer.incidence.cif$total.cifj.driptw$t <- cancer.incidence.cif$total.cifj.driptw$t/365.25
cancer.incidence.cif$sep.direct.astar1.cifj.driptw$t <- cancer.incidence.cif$sep.direct.astar1.cifj.driptw$t/365.25
cancer.incidence.cif$sep.direct.astar0.cifj.driptw$t <- cancer.incidence.cif$sep.direct.astar0.cifj.driptw$t/365.25
cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$t <- cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$t/365.25
cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$t <- cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$t/365.25

cancer.incidence.cif$controlled.direct.cifj.drow$t <- cancer.incidence.cif$controlled.direct.cifj.drow$t/365.25
cancer.incidence.cif$composite.cifj.drow$t <- cancer.incidence.cif$composite.cifj.drow$t/365.25
cancer.incidence.cif$total.cifj.drow$t <- cancer.incidence.cif$total.cifj.drow$t/365.25
cancer.incidence.cif$sep.direct.astar1.cifj.drow$t <- cancer.incidence.cif$sep.direct.astar1.cifj.drow$t/365.25
cancer.incidence.cif$sep.direct.astar0.cifj.drow$t <- cancer.incidence.cif$sep.direct.astar0.cifj.drow$t/365.25
cancer.incidence.cif$sep.indirect.astar1.cifj.drow$t <- cancer.incidence.cif$sep.indirect.astar1.cifj.drow$t/365.25
cancer.incidence.cif$sep.indirect.astar0.cifj.drow$t <- cancer.incidence.cif$sep.indirect.astar0.cifj.drow$t/365.25

### multiple causes of failure
## total
# ATE
total.ate.diff <- ggplot(cancer.incidence.cif$total.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, total effect") +
  theme_classic()

controlled.direct.ate.diff <- ggplot(cancer.incidence.cif$controlled.direct.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, controlled direct effect") +
  theme_classic()

sep.direct.astar1.ate.diff <- ggplot(cancer.incidence.cif$sep.direct.astar1.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(A1) ATE, separable direct effect, a"[d]*"=1")) +
  theme_classic()

sep.indirect.astar1.ate.diff <- ggplot(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(B1) ATE, separable indirect effect, a"[c]*"=0")) +
  theme_classic()

sep.direct.astar0.ate.diff <- ggplot(cancer.incidence.cif$sep.direct.astar0.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(C1) ATE, separable direct effect, a"[d]*"=0")) +
  theme_classic()

sep.indirect.astar0.ate.diff <- ggplot(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(D1) ATE, separable indirect effect, a"[c]*"=1")) +
  theme_classic()

# ATO
total.ato.diff <- ggplot(cancer.incidence.cif$total.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) OWATE, total effect") +
  theme_classic()

controlled.direct.ato.diff <- ggplot(cancer.incidence.cif$controlled.direct.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) OWATE, controlled direct effect") +
  theme_classic()

sep.direct.astar1.ato.diff <- ggplot(cancer.incidence.cif$sep.direct.astar1.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(A2) OWATE, separable direct effect, a"[d]*"=1")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.7))

sep.indirect.astar1.ato.diff <- ggplot(cancer.incidence.cif$sep.indirect.astar1.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(B2) OWATE, separable indirect effect, a"[c]*"=0")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.7))

sep.direct.astar0.ato.diff <- ggplot(cancer.incidence.cif$sep.direct.astar0.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(C2) OWATE, separable direct effect, a"[d]*"=0")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.7))

sep.indirect.astar0.ato.diff <- ggplot(cancer.incidence.cif$sep.indirect.astar0.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(D2) OWATE, separable indirect effect, a"[c]*"=1")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.7))

#
ggarrange(sep.direct.astar1.ate.diff, sep.direct.astar1.ato.diff,
          sep.indirect.astar1.ate.diff, sep.indirect.astar1.ato.diff,
          sep.direct.astar0.ate.diff, sep.direct.astar0.ato.diff,
          sep.indirect.astar0.ate.diff, sep.indirect.astar0.ato.diff,
          nrow = 4, ncol = 2)
ggsave("cancer_incidence_sep_cif_effect.pdf", width=8, height=16)

ggarrange(total.ate.diff, total.ato.diff, nrow = 1, ncol = 2)
ggsave("cancer_incidence_cif_total_effect.pdf", width=8, height=4)

ggarrange(controlled.direct.ate.diff, controlled.direct.ato.diff, nrow = 1, ncol = 2)
ggsave("cancer_incidence_cif_controlled_direct_effect.pdf", width=8, height=4)

# aje, each panel in a separate file
sep.direct.astar1.ate.diff
ggsave("cancer_incidence_sep_cif_effect_A1.pdf", width=4, height=4)

sep.direct.astar1.ato.diff
ggsave("cancer_incidence_sep_cif_effect_A2.pdf", width=4, height=4)

sep.indirect.astar1.ate.diff
ggsave("cancer_incidence_sep_cif_effect_B1.pdf", width=4, height=4)

sep.indirect.astar1.ato.diff
ggsave("cancer_incidence_sep_cif_effect_B2.pdf", width=4, height=4)

sep.direct.astar0.ate.diff
ggsave("cancer_incidence_sep_cif_effect_C1.pdf", width=4, height=4)

sep.direct.astar0.ato.diff
ggsave("cancer_incidence_sep_cif_effect_C2.pdf", width=4, height=4)

sep.indirect.astar0.ate.diff
ggsave("cancer_incidence_sep_cif_effect_D1.pdf", width=4, height=4)

sep.indirect.astar0.ato.diff
ggsave("cancer_incidence_sep_cif_effect_D2.pdf", width=4, height=4)

#
total.ate.diff
ggsave("cancer_incidence_cif_total_effect_A1.pdf", width=4, height=4)

total.ato.diff
ggsave("cancer_incidence_cif_total_effect_A2.pdf", width=4, height=4)

#
controlled.direct.ate.diff
ggsave("cancer_incidence_cif_controlled_direct_effect_A1.pdf", width=4, height=4)

controlled.direct.ato.diff
ggsave("cancer_incidence_cif_controlled_direct_effect_A2.pdf", width=4, height=4)

# defense
ggarrange(sep.direct.astar1.ate.diff, sep.direct.astar1.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_incidence_sep_direct_astar1_cif_effect.pdf", width=8, height=4)

ggarrange(sep.indirect.astar1.ate.diff, sep.indirect.astar1.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_incidence_sep_indirect_astar1_cif_effect.pdf", width=8, height=4)

ggarrange(sep.direct.astar0.ate.diff, sep.direct.astar0.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_incidence_sep_direct_astar0_cif_effect.pdf", width=8, height=4)

ggarrange(sep.indirect.astar0.ate.diff, sep.indirect.astar0.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_incidence_sep_indirect_astar0_cif_effect.pdf", width=8, height=4)



### treatment specific
## ATE
# a0
total.ate.a0 <- ggplot(cancer.incidence.cif$total.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.1, 0.2, 0.3)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, total effect, \n        Sulphonylurea (A=0)") +
  theme_classic()

controlled.direct.ate.a0 <- ggplot(cancer.incidence.cif$controlled.direct.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, controlled direct effect, \n        Sulphonylurea (A=0)") +
  theme_classic()

sep.direct.astar1.ate.a0 <- ggplot(cancer.incidence.cif$sep.direct.astar1.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(A1) ATE, separable direct effect, a"[d]*"=1,",
  "Sulphonylurea (A=0)              "))) +
  theme_classic()

sep.indirect.astar1.ate.a0 <- ggplot(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(B1) ATE, separable indirect effect, a"[c]*"=0,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

sep.direct.astar0.ate.a0 <- ggplot(cancer.incidence.cif$sep.direct.astar0.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(C1) ATE, separable direct effect, a"[d]*"=0,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

sep.indirect.astar0.ate.a0 <- ggplot(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(D1) ATE, separable indirect effect, a"[c]*"=1,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

# a1
total.ate.a1 <- ggplot(cancer.incidence.cif$total.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) ATE, total effect, \n        Metformin (A=1)") +
  theme_classic()

controlled.direct.ate.a1 <- ggplot(cancer.incidence.cif$controlled.direct.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) ATE, controlled direct effect, \n        Metformin (A=1)") +
  theme_classic()

sep.direct.astar1.ate.a1 <- ggplot(cancer.incidence.cif$sep.direct.astar1.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(A2) ATE, separable direct effect, a"[d]*"=1,",
                             "Metformin (A=1)                     "))) +
  theme_classic()

sep.indirect.astar1.ate.a1 <- ggplot(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(B2) ATE, separable indirect effect, a"[c]*"=0,",
                             "Metformin (A=1)                     "))) +
  theme_classic()

sep.direct.astar0.ate.a1 <- ggplot(cancer.incidence.cif$sep.direct.astar0.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(C2) ATE, separable direct effect, a"[d]*"=0,",
                             "Metformin (A=1)                     "))) +
  theme_classic()

sep.indirect.astar0.ate.a1 <- ggplot(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(D2) ATE, separable indirect effect, a"[c]*"=1,",
                             "Metformin (A=1)                     "))) +
  theme_classic()

## ATO
# a0
total.ato.a0 <- ggplot(cancer.incidence.cif$total.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) OWATE, total effect, \n        Sulphonylurea (A=0)") +
  theme_classic()

controlled.direct.ato.a0 <- ggplot(cancer.incidence.cif$controlled.direct.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) OWATE, controlled direct effect, \n        Sulphonylurea (A=0)") +
  theme_classic()

sep.direct.astar1.ato.a0 <- ggplot(cancer.incidence.cif$sep.direct.astar1.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(A1) OWATE, separable direct effect, a"[d]*"=1,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

sep.indirect.astar1.ato.a0 <- ggplot(cancer.incidence.cif$sep.indirect.astar1.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(B1) OWATE, separable indirect effect, a"[c]*"=0,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

sep.direct.astar0.ato.a0 <- ggplot(cancer.incidence.cif$sep.direct.astar0.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(C1) OWATE, separable direct effect, a"[d]*"=0,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

sep.indirect.astar0.ato.a0 <- ggplot(cancer.incidence.cif$sep.indirect.astar0.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(D1) OWATE, separable indirect effect, a"[c]*"=1,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

# a1
total.ato.a1 <- ggplot(cancer.incidence.cif$total.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B2) OWATE, total effect, \n        Metformin (A=1)") +
  theme_classic()

controlled.direct.ato.a1 <- ggplot(cancer.incidence.cif$controlled.direct.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B2) OWATE, controlled direct effect, \n        Metformin (A=1)") +
  theme_classic()

sep.direct.astar1.ato.a1 <- ggplot(cancer.incidence.cif$sep.direct.astar1.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(A2) OWATE, separable direct effect, a"[d]*"=1,",
                             "Metformin (A=1)                     "))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.7))

sep.indirect.astar1.ato.a1 <- ggplot(cancer.incidence.cif$sep.indirect.astar1.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(B2) OWATE, separable indirect effect, a"[c]*"=0,",
                             "Metformin (A=1)                     "))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.7))

sep.direct.astar0.ato.a1 <- ggplot(cancer.incidence.cif$sep.direct.astar0.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(C2) OWATE, separable direct effect, a"[d]*"=0,",
                             "Metformin (A=1)                     "))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.7))

sep.indirect.astar0.ato.a1 <- ggplot(cancer.incidence.cif$sep.indirect.astar0.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(D2) OWATE, separable indirect effect, a"[c]*"=1,",
                             "Metformin (A=1)                     "))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.7))

# ATE
ggarrange(sep.direct.astar1.ate.a0, sep.direct.astar1.ate.a1,
          sep.indirect.astar1.ate.a0, sep.indirect.astar1.ate.a1,
          sep.direct.astar0.ate.a0, sep.direct.astar0.ate.a1,
          sep.indirect.astar0.ate.a0, sep.indirect.astar0.ate.a1,
          nrow = 4, ncol = 2)
ggsave("cancer_incidence_sep_cif_ate_treatment_specific.pdf", width=8, height=16)

ggarrange(sep.direct.astar1.ato.a0, sep.direct.astar1.ato.a1,
          sep.indirect.astar1.ato.a0, sep.indirect.astar1.ato.a1,
          sep.direct.astar0.ato.a0, sep.direct.astar0.ato.a1,
          sep.indirect.astar0.ato.a0, sep.indirect.astar0.ato.a1,
          nrow = 4, ncol = 2)
ggsave("cancer_incidence_sep_cif_ato_treatment_specific.pdf", width=8, height=16)

ggarrange(total.ate.a0, total.ate.a1, total.ato.a0, total.ato.a1,
          nrow = 2, ncol = 2)
ggsave("cancer_incidence_cif_total_treatment_specific.pdf", width=8, height=8)

ggarrange(controlled.direct.ate.a0, controlled.direct.ate.a1, controlled.direct.ato.a0, controlled.direct.ato.a1,
          nrow = 2, ncol = 2)
ggsave("cancer_incidence_cif_controlled_treatment_specific.pdf", width=8, height=8)

## aje, each panel in a separate file
# ATE
sep.direct.astar1.ate.a0
ggsave("cancer_incidence_sep_cif_ate_treatment_specific_A1.pdf", width=4, height=4)

sep.direct.astar1.ate.a1
ggsave("cancer_incidence_sep_cif_ate_treatment_specific_A2.pdf", width=4, height=4)

sep.indirect.astar1.ate.a0
ggsave("cancer_incidence_sep_cif_ate_treatment_specific_B1.pdf", width=4, height=4)

sep.indirect.astar1.ate.a1
ggsave("cancer_incidence_sep_cif_ate_treatment_specific_B2.pdf", width=4, height=4)

sep.direct.astar0.ate.a0
ggsave("cancer_incidence_sep_cif_ate_treatment_specific_C1.pdf", width=4, height=4)

sep.direct.astar0.ate.a1
ggsave("cancer_incidence_sep_cif_ate_treatment_specific_C2.pdf", width=4, height=4)

sep.indirect.astar0.ate.a0
ggsave("cancer_incidence_sep_cif_ate_treatment_specific_D1.pdf", width=4, height=4)

sep.indirect.astar0.ate.a1
ggsave("cancer_incidence_sep_cif_ate_treatment_specific_D2.pdf", width=4, height=4)

# ATO
sep.direct.astar1.ato.a0
ggsave("cancer_incidence_sep_cif_ato_treatment_specific_A1.pdf", width=4, height=4)

sep.direct.astar1.ato.a1
ggsave("cancer_incidence_sep_cif_ato_treatment_specific_A2.pdf", width=4, height=4)

sep.indirect.astar1.ato.a0
ggsave("cancer_incidence_sep_cif_ato_treatment_specific_B1.pdf", width=4, height=4)

sep.indirect.astar1.ato.a1
ggsave("cancer_incidence_sep_cif_ato_treatment_specific_B2.pdf", width=4, height=4)

sep.direct.astar0.ato.a0
ggsave("cancer_incidence_sep_cif_ato_treatment_specific_C1.pdf", width=4, height=4)

sep.direct.astar0.ato.a1
ggsave("cancer_incidence_sep_cif_ato_treatment_specific_C2.pdf", width=4, height=4)

sep.indirect.astar0.ato.a0
ggsave("cancer_incidence_sep_cif_ato_treatment_specific_D1.pdf", width=4, height=4)

sep.indirect.astar0.ato.a1
ggsave("cancer_incidence_sep_cif_ato_treatment_specific_D2.pdf", width=4, height=4)

#
total.ate.a0
ggsave("cancer_incidence_cif_total_treatment_specific_A1.pdf", width=4, height=4)

total.ate.a1
ggsave("cancer_incidence_cif_total_treatment_specific_A2.pdf", width=4, height=4)

total.ato.a0
ggsave("cancer_incidence_cif_total_treatment_specific_B1.pdf", width=4, height=4)

total.ato.a1
ggsave("cancer_incidence_cif_total_treatment_specific_B2.pdf", width=4, height=4)

#
controlled.direct.ate.a0
ggsave("cancer_incidence_cif_controlled_treatment_specific_A1.pdf", width=4, height=4)

controlled.direct.ate.a1
ggsave("cancer_incidence_cif_controlled_treatment_specific_A2.pdf", width=4, height=4)

controlled.direct.ato.a0
ggsave("cancer_incidence_cif_controlled_treatment_specific_B1.pdf", width=4, height=4)

controlled.direct.ato.a1
ggsave("cancer_incidence_cif_controlled_treatment_specific_B2.pdf", width=4, height=4)

## cancer mortality
load("../mortality/natural_mortality_20231228.RData")
load("../mortality/cancer_mortality_20231228.RData")
load("../mortality/breast_cancer_mortality_20231228.RData")

natural.mortality.cif <- natural.mortality
cancer.mortality.cif <- cancer.mortality
breast.cancer.mortality.cif <- breast.cancer.mortality

# add bound
bound.var <- c(1, grep("se", names(cancer.mortality.cif$cifj.driptw)), grep("diff", names(cancer.mortality.cif$cifj.driptw)))
natural.mortality.cif$cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, natural.mortality.cif$cifj.driptw[, -bound.var])
natural.mortality.cif$cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, natural.mortality.cif$cifj.drow[, -bound.var])

cancer.mortality.cif$cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.mortality.cif$cifj.driptw[, -bound.var])
cancer.mortality.cif$cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.mortality.cif$cifj.drow[, -bound.var])

breast.cancer.mortality.cif$cifj.driptw[, -bound.var] <- Winsorize(minval = 0, maxval = 1, breast.cancer.mortality.cif$cifj.driptw[, -bound.var])
breast.cancer.mortality.cif$cifj.drow[, -bound.var] <- Winsorize(minval = 0, maxval = 1, breast.cancer.mortality.cif$cifj.drow[, -bound.var])

natural.mortality.cif$cifj.driptw$t <- natural.mortality.cif$cifj.driptw$t/365.25
natural.mortality.cif$cifj.drow$t <- natural.mortality.cif$cifj.drow$t/365.25

cancer.mortality.cif$cifj.driptw$t <- cancer.mortality.cif$cifj.driptw$t/365.25
cancer.mortality.cif$cifj.drow$t <- cancer.mortality.cif$cifj.drow$t/365.25

breast.cancer.mortality.cif$cifj.driptw$t <- breast.cancer.mortality.cif$cifj.driptw$t/365.25
breast.cancer.mortality.cif$cifj.drow$t <- breast.cancer.mortality.cif$cifj.drow$t/365.25

# ATE
natural.mortality.ate.diff <- ggplot(natural.mortality.cif$cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.01), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, all-cause mortality") +
  theme_classic()

composite.ate.diff <- ggplot(cancer.incidence.cif$composite.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.01), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) ATE, composite outcome") +
  theme_classic()

cancer.mortality.ate.diff <- ggplot(cancer.mortality.cif$cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0.01), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, mortality after cancer") +
  theme_classic()

breast.cancer.mortality.ate.diff <- ggplot(breast.cancer.mortality.cif$cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.8, 0.25), breaks = c(-0.75, -0.5, -0.25, 0, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) ATE, mortality after breast cancer") +
  theme_classic()

# ATO
natural.mortality.ato.diff <- ggplot(natural.mortality.cif$cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.01), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) OWATE, all-cause mortality") +
  theme_classic()

composite.ato.diff <- ggplot(cancer.incidence.cif$composite.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.01), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B2) OWATE, composite outcome") +
  theme_classic()

cancer.mortality.ato.diff <- ggplot(cancer.mortality.cif$cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0.01), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) OWATE, mortality after cancer") +
  theme_classic()

breast.cancer.mortality.ato.diff <- ggplot(breast.cancer.mortality.cif$cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.8, 0.25), breaks = c(-0.75, -0.5, -0.25, 0, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B2) OWATE, mortality after breast cancer") +
  theme_classic()

ggarrange(natural.mortality.ate.diff, natural.mortality.ato.diff,
          composite.ate.diff, composite.ato.diff,
          nrow = 2, ncol = 2)
ggsave("../figure/mortality_cif.pdf", width=8, height=8)

ggarrange(cancer.mortality.ate.diff, cancer.mortality.ato.diff,
          breast.cancer.mortality.ate.diff, breast.cancer.mortality.ato.diff,
          nrow = 2, ncol = 2)
ggsave("../figure/cancer_mortality_cif.pdf", width=8, height=8)

# aje
natural.mortality.ate.diff
ggsave("../figure/mortality_cif_A1.pdf", width=4, height=4)

natural.mortality.ato.diff
ggsave("../figure/mortality_cif_A2.pdf", width=4, height=4)

composite.ate.diff
ggsave("../figure/mortality_cif_B1.pdf", width=4, height=4)

composite.ato.diff
ggsave("../figure/mortality_cif_B2.pdf", width=4, height=4)

#
cancer.mortality.ate.diff
ggsave("../figure/cancer_mortality_cif_A1.pdf", width=4, height=4)

cancer.mortality.ato.diff
ggsave("../figure/cancer_mortality_cif_A2.pdf", width=4, height=4)

breast.cancer.mortality.ate.diff
ggsave("../figure/cancer_mortality_cif_B1.pdf", width=4, height=4)

breast.cancer.mortality.ato.diff
ggsave("../figure/cancer_mortality_cif_B2.pdf", width=4, height=4)

# defense
ggarrange(natural.mortality.ate.diff, natural.mortality.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/natural_mortality_cif.pdf", width=8, height=4)

ggarrange(composite.ate.diff, composite.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/composite_mortality_cif.pdf", width=8, height=4)

ggarrange(cancer.mortality.ate.diff, cancer.mortality.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_mortality_cif_only.pdf", width=8, height=4)

ggarrange(breast.cancer.mortality.ate.diff, breast.cancer.mortality.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/breast_cancer_mortality_cif.pdf", width=8, height=4)

### treatment specific
## ATE
# a0
natural.mortality.ate.a0 <- ggplot(natural.mortality.cif$cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, all-cause mortality, \n        Sulphonylurea (A=0)") +
  theme_classic()

composite.ate.a0 <- ggplot(cancer.incidence.cif$composite.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) ATE, composite outcome, \n        Sulphonylurea (A=0)") +
  theme_classic()

cancer.mortality.ate.a0 <- ggplot(cancer.mortality.cif$cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, mortality after cancer, \n        Sulphonylurea (A=0)") +
  theme_classic()

breast.cancer.mortality.ate.a0 <- ggplot(breast.cancer.mortality.cif$cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) ATE, mortality after breast cancer, \n        Sulphonylurea (A=0)") +
  theme_classic()

# a1
natural.mortality.ate.a1 <- ggplot(natural.mortality.cif$cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) ATE, all-cause mortality, \n        Metformin (A=1)") +
  theme_classic()

composite.ate.a1 <- ggplot(cancer.incidence.cif$composite.cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B2) ATE, composite outcome, \n        Metformin (A=1)") +
  theme_classic()

cancer.mortality.ate.a1 <- ggplot(cancer.mortality.cif$cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(C2) ATE, mortality after cancer, \n        Metformin (A=1)") +
  theme_classic()

breast.cancer.mortality.ate.a1 <- ggplot(breast.cancer.mortality.cif$cifj.driptw, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) ATE, mortality after breast cancer, \n        Metformin (A=1)") +
  theme_classic()

## ATO
# a0
natural.mortality.ato.a0 <- ggplot(natural.mortality.cif$cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) OWATE, all-cause mortality, \n        Sulphonylurea (A=0)") +
  theme_classic()

composite.ato.a0 <- ggplot(cancer.incidence.cif$composite.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) OWATE, composite outcome, \n        Sulphonylurea (A=0)") +
  theme_classic()

cancer.mortality.ato.a0 <- ggplot(cancer.mortality.cif$cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) OWATE, mortality after cancer, \n        Sulphonylurea (A=0)") +
  theme_classic()

breast.cancer.mortality.ato.a0 <- ggplot(breast.cancer.mortality.cif$cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) OWATE, mortality after breast cancer, \n        Sulphonylurea (A=0)") +
  theme_classic()

# a1
natural.mortality.ato.a1 <- ggplot(natural.mortality.cif$cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) OWATE, all-cause mortality, \n        Metformin (A=1)") +
  theme_classic()

composite.ato.a1 <- ggplot(cancer.incidence.cif$composite.cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B2) OWATE, composite outcome, \n        Metformin (A=1)") +
  theme_classic()

cancer.mortality.ato.a1 <- ggplot(cancer.mortality.cif$cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) OWATE, mortality after cancer, \n        Metformin (A=1)") +
  theme_classic()

breast.cancer.mortality.ato.a1 <- ggplot(breast.cancer.mortality.cif$cifj.drow, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B2) OWATE, mortality after breast cancer, \n        Metformin (A=1)") +
  theme_classic()

#
ggarrange(natural.mortality.ate.a0, natural.mortality.ate.a1,
          composite.ate.a0, composite.ate.a1, nrow = 2, ncol = 2)
ggsave("../figure/mortality_cif_ate_treatment_specific.pdf", width=8, height=8)

ggarrange(cancer.mortality.ate.a0, cancer.mortality.ate.a1,
          breast.cancer.mortality.ate.a0, breast.cancer.mortality.ate.a1, nrow = 2, ncol = 2)
ggsave("../figure/cancer_mortality_cif_ate_treatment_specific.pdf", width=8, height=8)

ggarrange(natural.mortality.ato.a0, natural.mortality.ato.a1,
          composite.ato.a0, composite.ato.a1, nrow = 2, ncol = 2)
ggsave("../figure/mortality_cif_ato_treatment_specific.pdf", width=8, height=8)

ggarrange(cancer.mortality.ato.a0, cancer.mortality.ato.a1,
          breast.cancer.mortality.ato.a0, breast.cancer.mortality.ato.a1, nrow = 2, ncol = 2)
ggsave("../figure/cancer_mortality_cif_ato_treatment_specific.pdf", width=8, height=8)

# aje
natural.mortality.ate.a0
ggsave("../figure/mortality_cif_ate_treatment_specific_A1.pdf", width=4, height=4)

natural.mortality.ate.a1
ggsave("../figure/mortality_cif_ate_treatment_specific_A2.pdf", width=4, height=4)

composite.ate.a0
ggsave("../figure/mortality_cif_ate_treatment_specific_B1.pdf", width=4, height=4)

composite.ate.a1
ggsave("../figure/mortality_cif_ate_treatment_specific_B2.pdf", width=4, height=4)

#
cancer.mortality.ate.a0
ggsave("../figure/cancer_mortality_cif_ate_treatment_specific_A1.pdf", width=4, height=4)

cancer.mortality.ate.a1
ggsave("../figure/cancer_mortality_cif_ate_treatment_specific_A2.pdf", width=4, height=4)

breast.cancer.mortality.ate.a0
ggsave("../figure/cancer_mortality_cif_ate_treatment_specific_B1.pdf", width=4, height=4)

breast.cancer.mortality.ate.a1
ggsave("../figure/cancer_mortality_cif_ate_treatment_specific_B2.pdf", width=4, height=4)

#
natural.mortality.ato.a0
ggsave("../figure/mortality_cif_ato_treatment_specific_A1.pdf", width=4, height=4)

natural.mortality.ato.a1
ggsave("../figure/mortality_cif_ato_treatment_specific_A2.pdf", width=4, height=4)

composite.ato.a0
ggsave("../figure/mortality_cif_ato_treatment_specific_B1.pdf", width=4, height=4)

composite.ato.a1
ggsave("../figure/mortality_cif_ato_treatment_specific_B2.pdf", width=4, height=4)

#
cancer.mortality.ato.a0
ggsave("../figure/cancer_mortality_cif_ato_treatment_specific_A1.pdf", width=4, height=4)

cancer.mortality.ato.a1
ggsave("../figure/cancer_mortality_cif_ato_treatment_specific_A2.pdf", width=4, height=4)

breast.cancer.mortality.ato.a0
ggsave("../figure/cancer_mortality_cif_ato_treatment_specific_B1.pdf", width=4, height=4)

breast.cancer.mortality.ato.a1
ggsave("../figure/cancer_mortality_cif_ato_treatment_specific_B2.pdf", width=4, height=4)


## crude incidence
load("crude/cancer_incidence_crude.RData")

cancer.incidence.crude <- cancer.incidence

# add bound
bound.var <- c(1, grep("se", names(cancer.incidence.crude$sep.direct.astar1.cifj.crude)), grep("diff", names(cancer.incidence.crude$sep.direct.astar1.cifj.crude)))
cancer.incidence.crude$controlled.direct.cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.crude$controlled.direct.cifj.crude[, -bound.var])
cancer.incidence.crude$composite.cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.crude$composite.cifj.crude[, -bound.var])
cancer.incidence.crude$total.cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.crude$total.cifj.crude[, -bound.var])
cancer.incidence.crude$sep.direct.astar1.cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.crude$sep.direct.astar1.cifj.crude[, -bound.var])
cancer.incidence.crude$sep.direct.astar0.cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.crude$sep.direct.astar0.cifj.crude[, -bound.var])
cancer.incidence.crude$sep.indirect.astar1.cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.crude$sep.indirect.astar1.cifj.crude[, -bound.var])
cancer.incidence.crude$sep.indirect.astar0.cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.incidence.crude$sep.indirect.astar0.cifj.crude[, -bound.var])

cancer.incidence.crude$controlled.direct.cifj.crude$t <- cancer.incidence.crude$controlled.direct.cifj.crude$t/365.25
cancer.incidence.crude$composite.cifj.crude$t <- cancer.incidence.crude$composite.cifj.crude$t/365.25
cancer.incidence.crude$total.cifj.crude$t <- cancer.incidence.crude$total.cifj.crude$t/365.25
cancer.incidence.crude$sep.direct.astar1.cifj.crude$t <- cancer.incidence.crude$sep.direct.astar1.cifj.crude$t/365.25
cancer.incidence.crude$sep.direct.astar0.cifj.crude$t <- cancer.incidence.crude$sep.direct.astar0.cifj.crude$t/365.25
cancer.incidence.crude$sep.indirect.astar1.cifj.crude$t <- cancer.incidence.crude$sep.indirect.astar1.cifj.crude$t/365.25
cancer.incidence.crude$sep.indirect.astar0.cifj.crude$t <- cancer.incidence.crude$sep.indirect.astar0.cifj.crude$t/365.25

### multiple causes of failure
## total
# ATE
total.ate.diff <- ggplot(cancer.incidence.crude$total.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, total effect") +
  theme_classic()

controlled.direct.ate.diff <- ggplot(cancer.incidence.crude$controlled.direct.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, controlled direct effect") +
  theme_classic()

sep.direct.astar1.ate.diff <- ggplot(cancer.incidence.crude$sep.direct.astar1.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(A1) ATE, separable direct effect, a"[d]*"=1")) +
  theme_classic()

sep.indirect.astar1.ate.diff <- ggplot(cancer.incidence.crude$sep.indirect.astar1.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(B1) ATE, separable indirect effect, a"[c]*"=0")) +
  theme_classic()

sep.direct.astar0.ate.diff <- ggplot(cancer.incidence.crude$sep.direct.astar0.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(C1) ATE, separable direct effect, a"[d]*"=0")) +
  theme_classic()

sep.indirect.astar0.ate.diff <- ggplot(cancer.incidence.crude$sep.indirect.astar0.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.04, 0.02), breaks = c(-0.04, -0.02, 0.00, 0.02)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression("(D1) ATE, separable indirect effect, a"[c]*"=1")) +
  theme_classic()

#
ggarrange(sep.direct.astar1.ate.diff, sep.direct.astar1.ato.diff,
          sep.indirect.astar1.ate.diff, sep.indirect.astar1.ato.diff,
          sep.direct.astar0.ate.diff, sep.direct.astar0.ato.diff,
          sep.indirect.astar0.ate.diff, sep.indirect.astar0.ato.diff,
          nrow = 4, ncol = 2)
ggsave("cancer_incidence_sep_cif_effect_newcode_newdata.pdf", width=8, height=16)

ggarrange(total.ate.diff, total.ato.diff, nrow = 1, ncol = 2)
ggsave("cancer_incidence_cif_total_effect_newcode_newdata.pdf", width=8, height=4)

ggarrange(controlled.direct.ate.diff, controlled.direct.ato.diff, nrow = 1, ncol = 2)
ggsave("cancer_incidence_cif_controlled_direct_effect_newcode_newdata.pdf", width=8, height=4)



# defense
ggarrange(sep.direct.astar1.ate.diff, sep.direct.astar1.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_incidence_sep_direct_astar1_cif_effect.pdf", width=8, height=4)

ggarrange(sep.indirect.astar1.ate.diff, sep.indirect.astar1.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_incidence_sep_indirect_astar1_cif_effect.pdf", width=8, height=4)

ggarrange(sep.direct.astar0.ate.diff, sep.direct.astar0.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_incidence_sep_direct_astar0_cif_effect.pdf", width=8, height=4)

ggarrange(sep.indirect.astar0.ate.diff, sep.indirect.astar0.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_incidence_sep_indirect_astar0_cif_effect.pdf", width=8, height=4)



### treatment specific
## ATE
# a0
total.ate.a0 <- ggplot(cancer.incidence.crude$total.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, total effect, \n        Sulphonylurea (A=0)") +
  theme_classic()

controlled.direct.ate.a0 <- ggplot(cancer.incidence.crude$controlled.direct.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, controlled direct effect, \n        Sulphonylurea (A=0)") +
  theme_classic()

sep.direct.astar1.ate.a0 <- ggplot(cancer.incidence.crude$sep.direct.astar1.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(A1) ATE, separable direct effect, a"[d]*"=1,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

sep.indirect.astar1.ate.a0 <- ggplot(cancer.incidence.crude$sep.indirect.astar1.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(B1) ATE, separable indirect effect, a"[c]*"=0,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

sep.direct.astar0.ate.a0 <- ggplot(cancer.incidence.crude$sep.direct.astar0.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(C1) ATE, separable direct effect, a"[d]*"=0,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

sep.indirect.astar0.ate.a0 <- ggplot(cancer.incidence.crude$sep.indirect.astar0.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(D1) ATE, separable indirect effect, a"[c]*"=1,",
                             "Sulphonylurea (A=0)              "))) +
  theme_classic()

# a1
total.ate.a1 <- ggplot(cancer.incidence.crude$total.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) ATE, total effect, \n        Metformin (A=1)") +
  theme_classic()

controlled.direct.ate.a1 <- ggplot(cancer.incidence.crude$controlled.direct.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) ATE, controlled direct effect, \n        Metformin (A=1)") +
  theme_classic()

sep.direct.astar1.ate.a1 <- ggplot(cancer.incidence.crude$sep.direct.astar1.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(A2) ATE, separable direct effect, a"[d]*"=1,",
                             "Metformin (A=1)                     "))) +
  theme_classic()

sep.indirect.astar1.ate.a1 <- ggplot(cancer.incidence.crude$sep.indirect.astar1.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(B2) ATE, separable indirect effect, a"[c]*"=0,",
                             "Metformin (A=1)                     "))) +
  theme_classic()

sep.direct.astar0.ate.a1 <- ggplot(cancer.incidence.crude$sep.direct.astar0.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(C2) ATE, separable direct effect, a"[d]*"=0,",
                             "Metformin (A=1)                     "))) +
  theme_classic()

sep.indirect.astar0.ate.a1 <- ggplot(cancer.incidence.crude$sep.indirect.astar0.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title=expression(atop("(D2) ATE, separable indirect effect, a"[c]*"=1,",
                             "Metformin (A=1)                     "))) +
  theme_classic()

# ATE
ggarrange(sep.direct.astar1.ate.a0, sep.direct.astar1.ate.a1,
          sep.indirect.astar1.ate.a0, sep.indirect.astar1.ate.a1,
          sep.direct.astar0.ate.a0, sep.direct.astar0.ate.a1,
          sep.indirect.astar0.ate.a0, sep.indirect.astar0.ate.a1,
          nrow = 4, ncol = 2)
ggsave("cancer_incidence_sep_cif_ate_treatment_specific_newcode_newdata.pdf", width=8, height=16)

ggarrange(sep.direct.astar1.ato.a0, sep.direct.astar1.ato.a1,
          sep.indirect.astar1.ato.a0, sep.indirect.astar1.ato.a1,
          sep.direct.astar0.ato.a0, sep.direct.astar0.ato.a1,
          sep.indirect.astar0.ato.a0, sep.indirect.astar0.ato.a1,
          nrow = 4, ncol = 2)
ggsave("cancer_incidence_sep_cif_ato_treatment_specific_newcode_newdata.pdf", width=8, height=16)

# ATO
ggarrange(total.ate.a0, total.ate.a1, total.ato.a0, total.ato.a1,
          nrow = 2, ncol = 2)
ggsave("cancer_incidence_cif_total_treatment_specific_newcode_newdata.pdf", width=8, height=8)

ggarrange(controlled.direct.ate.a0, controlled.direct.ate.a1, controlled.direct.ato.a0, controlled.direct.ato.a1,
          nrow = 2, ncol = 2)
ggsave("cancer_incidence_cif_controlled_treatment_specific_newcode_newdata.pdf", width=8, height=8)


## crude mortality
load("../mortality/crude/natural_mortality_crude.RData")
load("../mortality/crude/cancer_mortality_crude.RData")
load("../mortality/crude/breast_cancer_mortality_crude.RData")

natural.mortality.crude <- natural.mortality
cancer.mortality.crude <- cancer.mortality
breast.cancer.mortality.crude <- breast.cancer.mortality

# add bound
bound.var <- c(1, grep("se", names(cancer.mortality.crude$cifj.crude)), grep("diff", names(cancer.mortality.crude$cifj.crude)))
natural.mortality.crude$cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, natural.mortality.crude$cifj.crude[, -bound.var])
cancer.mortality.crude$cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, cancer.mortality.crude$cifj.crude[, -bound.var])
breast.cancer.mortality.crude$cifj.crude[, -bound.var] <- Winsorize(minval = 0, maxval = 1, breast.cancer.mortality.crude$cifj.crude[, -bound.var])

natural.mortality.crude$cifj.crude$t <- natural.mortality.crude$cifj.crude$t/365.25
cancer.mortality.crude$cifj.crude$t <- cancer.mortality.crude$cifj.crude$t/365.25
breast.cancer.mortality.crude$cifj.crude$t <- breast.cancer.mortality.crude$cifj.crude$t/365.25

# ATE
natural.mortality.ate.diff <- ggplot(natural.mortality.crude$cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0.05), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, all-cause mortality") +
  theme_classic()

composite.ate.diff <- ggplot(cancer.incidence.crude$composite.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0.05), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) ATE, composite outcome") +
  theme_classic()

cancer.mortality.ate.diff <- ggplot(cancer.mortality.crude$cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0.05), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(C1) ATE, mortality after cancer") +
  theme_classic()

breast.cancer.mortality.ate.diff <- ggplot(breast.cancer.mortality.crude$cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-1, 0.5), breaks = c(-1, -0.5, 0, 0.5)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(D1) ATE, mortality after breast cancer") +
  theme_classic()

ggarrange(natural.mortality.ate.diff, natural.mortality.ato.diff,
          composite.ate.diff, composite.ato.diff,
          cancer.mortality.ate.diff, cancer.mortality.ato.diff,
          breast.cancer.mortality.ate.diff, breast.cancer.mortality.ato.diff,
          nrow = 4, ncol = 2)
ggsave("../figure/cancer_mortality_cif_20231217.pdf", width=8, height=16)

# defense
ggarrange(natural.mortality.ate.diff, natural.mortality.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/natural_mortality_cif.pdf", width=8, height=4)

ggarrange(composite.ate.diff, composite.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/composite_mortality_cif.pdf", width=8, height=4)

ggarrange(cancer.mortality.ate.diff, cancer.mortality.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/cancer_mortality_cif_only.pdf", width=8, height=4)

ggarrange(breast.cancer.mortality.ate.diff, breast.cancer.mortality.ato.diff, nrow = 1, ncol = 2)
ggsave("figure/breast_cancer_mortality_cif.pdf", width=8, height=4)

### treatment specific
## ATE
# a0
natural.mortality.ate.a0 <- ggplot(natural.mortality.crude$cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A1) ATE, all-cause mortality, \n        Sulphonylurea (A=0)") +
  theme_classic()

composite.ate.a0 <- ggplot(cancer.incidence.crude$composite.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B1) ATE, composite outcome, \n        Sulphonylurea (A=0)") +
  theme_classic()

cancer.mortality.ate.a0 <- ggplot(cancer.mortality.crude$cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(C1) ATE, mortality after cancer, \n        Sulphonylurea (A=0)") +
  theme_classic()

breast.cancer.mortality.ate.a0 <- ggplot(breast.cancer.mortality.crude$cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(D1) ATE, mortality after cancer, \n        Sulphonylurea (A=0)") +
  theme_classic()

# a1
natural.mortality.ate.a1 <- ggplot(natural.mortality.crude$cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(A2) ATE, all-cause mortality, \n        Metformin (A=1)") +
  theme_classic()

composite.ate.a1 <- ggplot(cancer.incidence.crude$composite.cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(B2) ATE, composite outcome, \n        Metformin (A=1)") +
  theme_classic()

cancer.mortality.ate.a1 <- ggplot(cancer.mortality.crude$cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(C2) ATE, mortality after cancer, \n        Metformin (A=1)") +
  theme_classic()

breast.cancer.mortality.ate.a1 <- ggplot(breast.cancer.mortality.crude$cifj.crude, aes(x=t)) +
  geom_stepribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(D2) ATE, mortality after cancer, \n        Metformin (A=1)") +
  theme_classic()


#
ggarrange(natural.mortality.ate.a0, natural.mortality.ate.a1,
          composite.ate.a0, composite.ate.a1,
          cancer.mortality.ate.a0, cancer.mortality.ate.a1,
          breast.cancer.mortality.ate.a0, breast.cancer.mortality.ate.a1,
          nrow = 4, ncol = 2)
ggsave("figure/cancer_mortality_cif_ate_treatment_specific_20231217.pdf", width=8, height=16)

ggarrange(cancer.mortality.ate.diff, cancer.mortality.ate.a0, cancer.mortality.ate.a1, nrow = 1, ncol = 3)
ggsave("compare.pdf", width=12, height=4)

# direct kaplan-meier CIF plot
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$extraction_date <- max(data.df$lcd_n, na.rm=TRUE) # extraction_date <- as.Date(c("2016-12-31"))
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))

data.df %<>%
  # rowwise() %>%
  filter(prescription_date>cancer_eventdate) %>%
  mutate(event1 = death) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(deathdate_n,extraction_date,lcd_n,na.rm=TRUE),prescription_date,units="days"))) %>% # /365.242
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]
# data.df <- data.df[1:2000,]
time.grid <- 1
admin.cens <- max(data.df$end.time1)
n <- nrow(data.df)
data.df$a <- data.df$ref_sulf
data.df$end.time <- data.df$end.time1
data.df$event <- data.df$event1
data.df$id <- data.df$patid
data.df$naive <- 1
data.df <- data.df[order(data.df$patid),]


cancer.mortality.km <- data.frame(t=sort(unique(data.df$end.time1)))

cancer.mortality.km$estimate.a0 <- 1-summary(survfit(Surv(end.time1[data.df$a==0], 1*(event1[data.df$a==0]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$surv
cancer.mortality.km$estimate.a0.lower.ci <- 1-summary(survfit(Surv(end.time1[data.df$a==0], 1*(event1[data.df$a==0]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$upper
cancer.mortality.km$estimate.a0.upper.ci <- 1-summary(survfit(Surv(end.time1[data.df$a==0], 1*(event1[data.df$a==0]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$lower

cancer.mortality.km$estimate.a1 <- 1-summary(survfit(Surv(end.time1[data.df$a==1], 1*(event1[data.df$a==1]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$surv
cancer.mortality.km$estimate.a1.lower.ci <- 1-summary(survfit(Surv(end.time1[data.df$a==1], 1*(event1[data.df$a==1]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$upper
cancer.mortality.km$estimate.a1.upper.ci <- 1-summary(survfit(Surv(end.time1[data.df$a==1], 1*(event1[data.df$a==1]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$lower

cancer.mortality.km$estimate.diff <- cancer.mortality.km$estimate.a1-cancer.mortality.km$estimate.a0
cancer.mortality.km$estimate.diff.se <- sqrt(summary(survfit(Surv(end.time1[data.df$a==0], 1*(event1[data.df$a==0]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$std.err^2+
                                               summary(survfit(Surv(end.time1[data.df$a==1], 1*(event1[data.df$a==1]!=0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$std.err)
cancer.mortality.km$estimate.diff.lower.ci <- with(cancer.mortality.km, estimate.diff+qnorm(0.975)*estimate.diff.se) 
cancer.mortality.km$estimate.diff.upper.ci <- with(cancer.mortality.km, estimate.diff+qnorm(0.025)*estimate.diff.se)

cancer.mortality.km$t <- cancer.mortality.km$t/365.25

cancer.mortality.ate.diff <- ggplot(cancer.mortality.km, aes(x=t)) +
  geom_step(aes(y=estimate.diff), color="black") +
  geom_step(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0.05), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(C1) ATE, mortality after cancer") +
  theme_classic()

cancer.mortality.ate.a0 <- ggplot(cancer.mortality.km, aes(x=t)) +
  geom_step(aes(y=estimate.a0), color="black") +
  geom_step(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(C1) ATE, mortality after cancer, \n        Sulphonylurea (A=0)") +
  theme_classic()

cancer.mortality.ate.a1 <- ggplot(cancer.mortality.km, aes(x=t)) +
  geom_step(aes(y=estimate.a1), color="black") +
  geom_step(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_step(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = c(0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(C1) ATE, mortality after cancer, \n        Sulphonylurea (A=0)") +
  theme_classic()

ggarrange(cancer.mortality.ate.diff, cancer.mortality.ate.a0, cancer.mortality.ate.a1, nrow = 1, ncol = 3)
ggsave("compare_km.pdf", width=12, height=4)





### curves plot
## all-cause mortality
load("mortality/natural_mortality.RData")
natural.mortality$cifj.driptw <- natural.mortality$cifj.driptw/365.25
natural.mortality$cifj.drow <- natural.mortality$cifj.drow/365.25

# ATE
natural.mortality.ate.diff <- ggplot(natural.mortality$cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-4, 0), breaks = c(0, -1, -2, -3, -4)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)            ATE, all-cause mortality") +
  theme_classic()

# ATO
natural.mortality.ato.diff <- ggplot(natural.mortality$cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-4, 0), breaks = c(0, -1, -2, -3, -4)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)            OWATE, all-cause mortality") +
  theme_classic()

ggarrange(natural.mortality.ate.diff, natural.mortality.ato.diff,
          nrow = 1, ncol = 2)
ggsave("figure/natural_mortality_effect.pdf", width=8, height=4)

## cancer incidence
load("cancer incidence/cancer_incidence.RData")
cancer.incidence$controlled.direct.cifj.driptw <- cancer.incidence$controlled.direct.cifj.driptw/365.25
cancer.incidence$composite.cifj.driptw <- cancer.incidence$composite.cifj.driptw/365.25
cancer.incidence$total.cifj.driptw <- cancer.incidence$total.cifj.driptw/365.25
cancer.incidence$sep.direct.cifj.driptw <- cancer.incidence$sep.direct.cifj.driptw/365.25
cancer.incidence$sep.indirect.cifj.driptw <- cancer.incidence$sep.indirect.cifj.driptw/365.25

cancer.incidence$controlled.direct.cifj.drow <- cancer.incidence$controlled.direct.cifj.drow/365.25
cancer.incidence$composite.cifj.drow <- cancer.incidence$composite.cifj.drow/365.25
cancer.incidence$total.cifj.drow <- cancer.incidence$total.cifj.drow/365.25
cancer.incidence$sep.direct.cifj.drow <- cancer.incidence$sep.direct.cifj.drow/365.25
cancer.incidence$sep.indirect.cifj.drow <- cancer.incidence$sep.indirect.cifj.drow/365.25


####################################################################################
## breast cancer incidence
load("breast cancer incidence/breast_cancer_incidence.RData")
breast.cancer.incidence$controlled.direct.cifj.driptw <- breast.cancer.incidence$controlled.direct.cifj.driptw/365.25
breast.cancer.incidence$composite.cifj.driptw <- breast.cancer.incidence$composite.cifj.driptw/365.25
breast.cancer.incidence$total.cifj.driptw <- breast.cancer.incidence$total.cifj.driptw/365.25
breast.cancer.incidence$sep.direct.cifj.driptw <- breast.cancer.incidence$sep.direct.cifj.driptw/365.25
breast.cancer.incidence$sep.indirect.cifj.driptw <- breast.cancer.incidence$sep.indirect.cifj.driptw/365.25

breast.cancer.incidence$controlled.direct.cifj.drow <- breast.cancer.incidence$controlled.direct.cifj.drow/365.25
breast.cancer.incidence$composite.cifj.drow <- breast.cancer.incidence$composite.cifj.drow/365.25
breast.cancer.incidence$total.cifj.drow <- breast.cancer.incidence$total.cifj.drow/365.25
breast.cancer.incidence$sep.direct.cifj.drow <- breast.cancer.incidence$sep.direct.cifj.drow/365.25
breast.cancer.incidence$sep.indirect.cifj.drow <- breast.cancer.incidence$sep.indirect.cifj.drow/365.25

### multiple causes of failure
## total
# ATE
total.ate.diff <- ggplot(breast.cancer.incidence$total.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)               ATE, total effect") +
  theme_classic()

# ATO
total.ato.diff <- ggplot(breast.cancer.incidence$total.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)               OWATE, total effect") +
  theme_classic()

## sep direct
# ATE
sep.direct.ate.diff <- ggplot(breast.cancer.incidence$sep.direct.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)         ATE, separable direct effect") +
  theme_classic()

# ATO
sep.direct.ato.diff <- ggplot(breast.cancer.incidence$sep.direct.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)         OWATE, separable direct effect") +
  theme_classic()

## sep indirect
# ATE
sep.indirect.ate.diff <- ggplot(breast.cancer.incidence$sep.indirect.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)       ATE, separable indirect effect") +
  theme_classic()

# ATO
sep.indirect.ato.diff <- ggplot(breast.cancer.incidence$sep.indirect.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       OWATE, separable indirect effect") +
  theme_classic()

### one cause of failure
## controlled.direct
# ATE
controlled.direct.ate.diff <- ggplot(breast.cancer.incidence$controlled.direct.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)       ATE, controlled direct effect") +
  theme_classic()

# ATO
controlled.direct.ato.diff <- ggplot(breast.cancer.incidence$controlled.direct.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       OWATE, controlled direct effect") +
  theme_classic()

### mortality included
## composite
# ATE
composite.ate.diff <- ggplot(breast.cancer.incidence$composite.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)                 ATE, composite") +
  theme_classic()

# ATO
composite.ato.diff <- ggplot(breast.cancer.incidence$composite.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)                 OWATE, composite") +
  theme_classic()

## breast.cancer mortality
load("mortality/breast_cancer_mortality.RData")
breast.cancer.mortality$cifj.driptw <- breast.cancer.mortality$cifj.driptw/365.25
breast.cancer.mortality$cifj.drow <- breast.cancer.mortality$cifj.drow/365.25

# ATE
breast.cancer.mortality.ate.diff <- ggplot(breast.cancer.mortality$cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-15, 3), breaks = c(-15, -12, -9, -6, -3, 0, 3)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       ATE, mortality after breast cancer") +
  theme_classic()

# ATO
breast.cancer.mortality.ato.diff <- ggplot(breast.cancer.mortality$cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-15, 3), breaks = c(-15, -12, -9, -6, -3, 0, 3)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       OWATE, mortality after breast cancer") +
  theme_classic()

## composite
# ATE
composite.ate.diff <- ggplot(breast.cancer.incidence$composite.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-14, 3), breaks = c(-14, -12, -9, -6, -3, 0, 3)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)         ATE, composite outcome") +
  theme_classic()

# ATO
composite.ato.diff <- ggplot(breast.cancer.incidence$composite.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-14, 3), breaks = c(-14, -12, -9, -6, -3, 0, 3)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)         OWATE, composite outcome") +
  theme_classic()

ggarrange(total.ate.diff, total.ato.diff,
          sep.direct.ate.diff, sep.direct.ato.diff,
          sep.indirect.ate.diff, sep.indirect.ato.diff,
          controlled.direct.ate.diff, controlled.direct.ato.diff,
          nrow = 4, ncol = 2)
ggsave("figure/breast_cancer_incidence_effect.pdf", width=8, height=16)

ggarrange(breast.cancer.mortality.ate.diff, breast.cancer.mortality.ato.diff,
          composite.ate.diff, composite.ato.diff,
          nrow = 2, ncol = 2)
ggsave("figure/breast_cancer_mortality_effect.pdf", width=8, height=8)

#
ggarrange(breast.cancer.mortality.ate.diff, breast.cancer.mortality.ato.diff,
          nrow = 1, ncol = 2)
ggsave("figure/breast_cancer_mortality_only_effect.pdf", width=8, height=4)

####################################################################################
### table
## cancer incidence
## all-cause mortality
# ATE
paste0(round(natural.mortality.cif$cifj.driptw$estimate.diff[max(which(natural.mortality.cif$cifj.driptw$t<tau))], digits=3), " (", round(natural.mortality.cif$cifj.driptw$estimate.diff.se[max(which(natural.mortality.cif$cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(natural.mortality.cif$cifj.driptw$estimate.diff.lower.ci[max(which(natural.mortality.cif$cifj.driptw$t<tau))], digits=3), ", ", round(natural.mortality.cif$cifj.driptw$estimate.diff.upper.ci[max(which(natural.mortality.cif$cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(natural.mortality.cif$cifj.driptw$estimate.diff.lower.cb[max(which(natural.mortality.cif$cifj.driptw$t<tau))], digits=3), ", ", round(natural.mortality.cif$cifj.driptw$estimate.diff.upper.cb[max(which(natural.mortality.cif$cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(natural.mortality.cif$cifj.drow$t, natural.mortality.cif$cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(natural.mortality.cif$cifj.drow$t, natural.mortality.cif$cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(natural.mortality.cif$cifj.drow$t, natural.mortality.cif$cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(natural.mortality.cif$cifj.drow$t, natural.mortality.cif$cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(natural.mortality.cif$cifj.drow$t, natural.mortality.cif$cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(natural.mortality.cif$cifj.drow$t, natural.mortality.cif$cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(natural.mortality.crude$cifj.crude$estimate.diff[max(which(natural.mortality.crude$cifj.crude$t<tau))], digits=3), " (", round(natural.mortality.crude$cifj.crude$estimate.diff.se[max(which(natural.mortality.crude$cifj.crude$t<tau))], digits=3), ")",
       " & (", round(natural.mortality.crude$cifj.crude$estimate.diff.lower.ci[max(which(natural.mortality.crude$cifj.crude$t<tau))], digits=3), ", ", round(natural.mortality.crude$cifj.crude$estimate.diff.upper.ci[max(which(natural.mortality.crude$cifj.crude$t<tau))], digits=3), ")",
       " & (", round(natural.mortality.crude$cifj.crude$estimate.diff.lower.cb[max(which(natural.mortality.crude$cifj.crude$t<tau))], digits=3), ", ", round(natural.mortality.crude$cifj.crude$estimate.diff.upper.cb[max(which(natural.mortality.crude$cifj.crude$t<tau))], digits=3), ")")

## composite outcome
# ATE
paste0(round(cancer.incidence.cif$composite.cifj.driptw$estimate.diff[max(which(cancer.incidence.cif$composite.cifj.driptw$t<tau))], digits=3), " (", round(cancer.incidence.cif$composite.cifj.driptw$estimate.diff.se[max(which(cancer.incidence.cif$composite.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$composite.cifj.driptw$estimate.diff.lower.ci[max(which(cancer.incidence.cif$composite.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$composite.cifj.driptw$estimate.diff.upper.ci[max(which(cancer.incidence.cif$composite.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$composite.cifj.driptw$estimate.diff.lower.cb[max(which(cancer.incidence.cif$composite.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$composite.cifj.driptw$estimate.diff.upper.cb[max(which(cancer.incidence.cif$composite.cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(cancer.incidence.cif$composite.cifj.drow$t, cancer.incidence.cif$composite.cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.cif$composite.cifj.drow$t, cancer.incidence.cif$composite.cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$composite.cifj.drow$t, cancer.incidence.cif$composite.cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$composite.cifj.drow$t, cancer.incidence.cif$composite.cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$composite.cifj.drow$t, cancer.incidence.cif$composite.cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$composite.cifj.drow$t, cancer.incidence.cif$composite.cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(cancer.incidence.crude$composite.cifj.crude$estimate.diff[max(which(cancer.incidence.crude$composite.cifj.crude$t<tau))], digits=3), " (", round(cancer.incidence.crude$composite.cifj.crude$estimate.diff.se[max(which(cancer.incidence.crude$composite.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$composite.cifj.crude$estimate.diff.lower.ci[max(which(cancer.incidence.crude$composite.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$composite.cifj.crude$estimate.diff.upper.ci[max(which(cancer.incidence.crude$composite.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$composite.cifj.crude$estimate.diff.lower.cb[max(which(cancer.incidence.crude$composite.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$composite.cifj.crude$estimate.diff.upper.cb[max(which(cancer.incidence.crude$composite.cifj.crude$t<tau))], digits=3), ")")

## total effect
# ATE
paste0(round(cancer.incidence.cif$total.cifj.driptw$estimate.diff[max(which(cancer.incidence.cif$total.cifj.driptw$t<tau))], digits=3), " (", round(cancer.incidence.cif$total.cifj.driptw$estimate.diff.se[max(which(cancer.incidence.cif$total.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$total.cifj.driptw$estimate.diff.lower.ci[max(which(cancer.incidence.cif$total.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$total.cifj.driptw$estimate.diff.upper.ci[max(which(cancer.incidence.cif$total.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$total.cifj.driptw$estimate.diff.lower.cb[max(which(cancer.incidence.cif$total.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$total.cifj.driptw$estimate.diff.upper.cb[max(which(cancer.incidence.cif$total.cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(cancer.incidence.cif$total.cifj.drow$t, cancer.incidence.cif$total.cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.cif$total.cifj.drow$t, cancer.incidence.cif$total.cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$total.cifj.drow$t, cancer.incidence.cif$total.cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$total.cifj.drow$t, cancer.incidence.cif$total.cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$total.cifj.drow$t, cancer.incidence.cif$total.cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$total.cifj.drow$t, cancer.incidence.cif$total.cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(cancer.incidence.crude$total.cifj.crude$estimate.diff[max(which(cancer.incidence.crude$total.cifj.crude$t<tau))], digits=3), " (", round(cancer.incidence.crude$total.cifj.crude$estimate.diff.se[max(which(cancer.incidence.crude$total.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$total.cifj.crude$estimate.diff.lower.ci[max(which(cancer.incidence.crude$total.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$total.cifj.crude$estimate.diff.upper.ci[max(which(cancer.incidence.crude$total.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$total.cifj.crude$estimate.diff.lower.cb[max(which(cancer.incidence.crude$total.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$total.cifj.crude$estimate.diff.upper.cb[max(which(cancer.incidence.crude$total.cifj.crude$t<tau))], digits=3), ")")

## controlled direct
# ATE
paste0(round(cancer.incidence.cif$controlled.direct.cifj.driptw$estimate.diff[max(which(cancer.incidence.cif$controlled.direct.cifj.driptw$t<tau))], digits=3), " (", round(cancer.incidence.cif$controlled.direct.cifj.driptw$estimate.diff.se[max(which(cancer.incidence.cif$controlled.direct.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$controlled.direct.cifj.driptw$estimate.diff.lower.ci[max(which(cancer.incidence.cif$controlled.direct.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$controlled.direct.cifj.driptw$estimate.diff.upper.ci[max(which(cancer.incidence.cif$controlled.direct.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$controlled.direct.cifj.driptw$estimate.diff.lower.cb[max(which(cancer.incidence.cif$controlled.direct.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$controlled.direct.cifj.driptw$estimate.diff.upper.cb[max(which(cancer.incidence.cif$controlled.direct.cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(cancer.incidence.cif$controlled.direct.cifj.drow$t, cancer.incidence.cif$controlled.direct.cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.cif$controlled.direct.cifj.drow$t, cancer.incidence.cif$controlled.direct.cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$controlled.direct.cifj.drow$t, cancer.incidence.cif$controlled.direct.cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$controlled.direct.cifj.drow$t, cancer.incidence.cif$controlled.direct.cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$controlled.direct.cifj.drow$t, cancer.incidence.cif$controlled.direct.cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$controlled.direct.cifj.drow$t, cancer.incidence.cif$controlled.direct.cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(cancer.incidence.crude$controlled.direct.cifj.crude$estimate.diff[max(which(cancer.incidence.crude$controlled.direct.cifj.crude$t<tau))], digits=3), " (", round(cancer.incidence.crude$controlled.direct.cifj.crude$estimate.diff.se[max(which(cancer.incidence.crude$controlled.direct.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$controlled.direct.cifj.crude$estimate.diff.lower.ci[max(which(cancer.incidence.crude$controlled.direct.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$controlled.direct.cifj.crude$estimate.diff.upper.ci[max(which(cancer.incidence.crude$controlled.direct.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$controlled.direct.cifj.crude$estimate.diff.lower.cb[max(which(cancer.incidence.crude$controlled.direct.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$controlled.direct.cifj.crude$estimate.diff.upper.cb[max(which(cancer.incidence.crude$controlled.direct.cifj.crude$t<tau))], digits=3), ")")

## separable
# ATE
paste0(round(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$estimate.diff[max(which(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$t<tau))], digits=3), " (", round(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$estimate.diff.se[max(which(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$estimate.diff.lower.ci[max(which(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$estimate.diff.upper.ci[max(which(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$estimate.diff.lower.cb[max(which(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$estimate.diff.upper.cb[max(which(cancer.incidence.cif$sep.direct.astar1.cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(cancer.incidence.cif$sep.direct.astar1.cifj.drow$t, cancer.incidence.cif$sep.direct.astar1.cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.cif$sep.direct.astar1.cifj.drow$t, cancer.incidence.cif$sep.direct.astar1.cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$sep.direct.astar1.cifj.drow$t, cancer.incidence.cif$sep.direct.astar1.cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$sep.direct.astar1.cifj.drow$t, cancer.incidence.cif$sep.direct.astar1.cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$sep.direct.astar1.cifj.drow$t, cancer.incidence.cif$sep.direct.astar1.cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$sep.direct.astar1.cifj.drow$t, cancer.incidence.cif$sep.direct.astar1.cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(cancer.incidence.crude$sep.direct.astar1.cifj.crude$estimate.diff[max(which(cancer.incidence.crude$sep.direct.astar1.cifj.crude$t<tau))], digits=3), " (", round(cancer.incidence.crude$sep.direct.astar1.cifj.crude$estimate.diff.se[max(which(cancer.incidence.crude$sep.direct.astar1.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$sep.direct.astar1.cifj.crude$estimate.diff.lower.ci[max(which(cancer.incidence.crude$sep.direct.astar1.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$sep.direct.astar1.cifj.crude$estimate.diff.upper.ci[max(which(cancer.incidence.crude$sep.direct.astar1.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$sep.direct.astar1.cifj.crude$estimate.diff.lower.cb[max(which(cancer.incidence.crude$sep.direct.astar1.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$sep.direct.astar1.cifj.crude$estimate.diff.upper.cb[max(which(cancer.incidence.crude$sep.direct.astar1.cifj.crude$t<tau))], digits=3), ")")

# ATE
paste0(round(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$estimate.diff[max(which(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$t<tau))], digits=3), " (", round(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$estimate.diff.se[max(which(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$estimate.diff.lower.ci[max(which(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$estimate.diff.upper.ci[max(which(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$estimate.diff.lower.cb[max(which(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$estimate.diff.upper.cb[max(which(cancer.incidence.cif$sep.indirect.astar1.cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(cancer.incidence.cif$sep.indirect.astar1.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar1.cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.cif$sep.indirect.astar1.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar1.cifj.drow$estimate.diff.se, xout=tau)$y, digits=4), ")",
       " & (", round(approx(cancer.incidence.cif$sep.indirect.astar1.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar1.cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$sep.indirect.astar1.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar1.cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$sep.indirect.astar1.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar1.cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$sep.indirect.astar1.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar1.cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$estimate.diff[max(which(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$t<tau))], digits=3), " (", round(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$estimate.diff.se[max(which(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$t<tau))], digits=4), ")",
       " & (", round(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$estimate.diff.lower.ci[max(which(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$estimate.diff.upper.ci[max(which(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$estimate.diff.lower.cb[max(which(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$estimate.diff.upper.cb[max(which(cancer.incidence.crude$sep.indirect.astar1.cifj.crude$t<tau))], digits=3), ")")

# ATE
paste0(round(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$estimate.diff[max(which(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$t<tau))], digits=3), " (", round(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$estimate.diff.se[max(which(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$estimate.diff.lower.ci[max(which(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$estimate.diff.upper.ci[max(which(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$estimate.diff.lower.cb[max(which(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$estimate.diff.upper.cb[max(which(cancer.incidence.cif$sep.direct.astar0.cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(cancer.incidence.cif$sep.direct.astar0.cifj.drow$t, cancer.incidence.cif$sep.direct.astar0.cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.cif$sep.direct.astar0.cifj.drow$t, cancer.incidence.cif$sep.direct.astar0.cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$sep.direct.astar0.cifj.drow$t, cancer.incidence.cif$sep.direct.astar0.cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$sep.direct.astar0.cifj.drow$t, cancer.incidence.cif$sep.direct.astar0.cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$sep.direct.astar0.cifj.drow$t, cancer.incidence.cif$sep.direct.astar0.cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$sep.direct.astar0.cifj.drow$t, cancer.incidence.cif$sep.direct.astar0.cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(cancer.incidence.crude$sep.direct.astar0.cifj.crude$estimate.diff[max(which(cancer.incidence.crude$sep.direct.astar0.cifj.crude$t<tau))], digits=3), " (", round(cancer.incidence.crude$sep.direct.astar0.cifj.crude$estimate.diff.se[max(which(cancer.incidence.crude$sep.direct.astar0.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$sep.direct.astar0.cifj.crude$estimate.diff.lower.ci[max(which(cancer.incidence.crude$sep.direct.astar0.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$sep.direct.astar0.cifj.crude$estimate.diff.upper.ci[max(which(cancer.incidence.crude$sep.direct.astar0.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$sep.direct.astar0.cifj.crude$estimate.diff.lower.cb[max(which(cancer.incidence.crude$sep.direct.astar0.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$sep.direct.astar0.cifj.crude$estimate.diff.upper.cb[max(which(cancer.incidence.crude$sep.direct.astar0.cifj.crude$t<tau))], digits=3), ")")

# ATE
paste0(round(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$estimate.diff[max(which(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$t<tau))], digits=3), " (", round(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$estimate.diff.se[max(which(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$estimate.diff.lower.ci[max(which(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$estimate.diff.upper.ci[max(which(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$estimate.diff.lower.cb[max(which(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$t<tau))], digits=3), ", ", round(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$estimate.diff.upper.cb[max(which(cancer.incidence.cif$sep.indirect.astar0.cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(cancer.incidence.cif$sep.indirect.astar0.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar0.cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.cif$sep.indirect.astar0.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar0.cifj.drow$estimate.diff.se, xout=tau)$y, digits=4), ")",
       " & (", round(approx(cancer.incidence.cif$sep.indirect.astar0.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar0.cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$sep.indirect.astar0.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar0.cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.incidence.cif$sep.indirect.astar0.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar0.cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.cif$sep.indirect.astar0.cifj.drow$t, cancer.incidence.cif$sep.indirect.astar0.cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$estimate.diff[max(which(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$t<tau))], digits=3), " (", round(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$estimate.diff.se[max(which(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$t<tau))], digits=4), ")",
       " & (", round(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$estimate.diff.lower.ci[max(which(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$estimate.diff.upper.ci[max(which(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$estimate.diff.lower.cb[max(which(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$t<tau))], digits=3), ", ", round(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$estimate.diff.upper.cb[max(which(cancer.incidence.crude$sep.indirect.astar0.cifj.crude$t<tau))], digits=3), ")")

## cancer mortality
# ATE
paste0(round(cancer.mortality.cif$cifj.driptw$estimate.diff[max(which(cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), " (", round(cancer.mortality.cif$cifj.driptw$estimate.diff.se[max(which(cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.mortality.cif$cifj.driptw$estimate.diff.lower.ci[max(which(cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ", ", round(cancer.mortality.cif$cifj.driptw$estimate.diff.upper.ci[max(which(cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(cancer.mortality.cif$cifj.driptw$estimate.diff.lower.cb[max(which(cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ", ", round(cancer.mortality.cif$cifj.driptw$estimate.diff.upper.cb[max(which(cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(cancer.mortality.cif$cifj.drow$t, cancer.mortality.cif$cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.mortality.cif$cifj.drow$t, cancer.mortality.cif$cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.mortality.cif$cifj.drow$t, cancer.mortality.cif$cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality.cif$cifj.drow$t, cancer.mortality.cif$cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(cancer.mortality.cif$cifj.drow$t, cancer.mortality.cif$cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality.cif$cifj.drow$t, cancer.mortality.cif$cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(cancer.mortality.crude$cifj.crude$estimate.diff[max(which(cancer.mortality.crude$cifj.crude$t<tau))], digits=3), " (", round(cancer.mortality.crude$cifj.crude$estimate.diff.se[max(which(cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.mortality.crude$cifj.crude$estimate.diff.lower.ci[max(which(cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ", ", round(cancer.mortality.crude$cifj.crude$estimate.diff.upper.ci[max(which(cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ")",
       " & (", round(cancer.mortality.crude$cifj.crude$estimate.diff.lower.cb[max(which(cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ", ", round(cancer.mortality.crude$cifj.crude$estimate.diff.upper.cb[max(which(cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ")")

## breast cancer mortality
# ATE
paste0(round(breast.cancer.mortality.cif$cifj.driptw$estimate.diff[max(which(breast.cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), " (", round(breast.cancer.mortality.cif$cifj.driptw$estimate.diff.se[max(which(breast.cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(breast.cancer.mortality.cif$cifj.driptw$estimate.diff.lower.ci[max(which(breast.cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ", ", round(breast.cancer.mortality.cif$cifj.driptw$estimate.diff.upper.ci[max(which(breast.cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ")",
       " & (", round(breast.cancer.mortality.cif$cifj.driptw$estimate.diff.lower.cb[max(which(breast.cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ", ", round(breast.cancer.mortality.cif$cifj.driptw$estimate.diff.upper.cb[max(which(breast.cancer.mortality.cif$cifj.driptw$t<tau))], digits=3), ")")
# ATO
paste0(round(approx(breast.cancer.mortality.cif$cifj.drow$t, breast.cancer.mortality.cif$cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(breast.cancer.mortality.cif$cifj.drow$t, breast.cancer.mortality.cif$cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(breast.cancer.mortality.cif$cifj.drow$t, breast.cancer.mortality.cif$cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(breast.cancer.mortality.cif$cifj.drow$t, breast.cancer.mortality.cif$cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(breast.cancer.mortality.cif$cifj.drow$t, breast.cancer.mortality.cif$cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(breast.cancer.mortality.cif$cifj.drow$t, breast.cancer.mortality.cif$cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# crude
paste0(round(breast.cancer.mortality.crude$cifj.crude$estimate.diff[max(which(breast.cancer.mortality.crude$cifj.crude$t<tau))], digits=3), " (", round(breast.cancer.mortality.crude$cifj.crude$estimate.diff.se[max(which(breast.cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ")",
       " & (", round(breast.cancer.mortality.crude$cifj.crude$estimate.diff.lower.ci[max(which(breast.cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ", ", round(breast.cancer.mortality.crude$cifj.crude$estimate.diff.upper.ci[max(which(breast.cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ")",
       " & (", round(breast.cancer.mortality.crude$cifj.crude$estimate.diff.lower.cb[max(which(breast.cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ", ", round(breast.cancer.mortality.crude$cifj.crude$estimate.diff.upper.cb[max(which(breast.cancer.mortality.crude$cifj.crude$t<tau))], digits=3), ")")

# ## cancer mortality
# # ATE
# paste0(round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a0.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")
# 
# paste0(round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a1.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")
# 
# paste0(round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.diff.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.driptw$t, cancer.mortality$cifj.driptw$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# 
# # ATO
# paste0(round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a0.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")
# 
# paste0(round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a1.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")
# 
# paste0(round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.mortality$cifj.drow$t, cancer.mortality$cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# 
# ## composite
# # ATE
# paste0(round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a0.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")
# 
# paste0(round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a1.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")
# 
# paste0(round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.diff.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.driptw$t, cancer.incidence.crude$composite.cifj.driptw$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")
# 
# # ATO
# paste0(round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a0.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")
# 
# paste0(round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a1.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")
# 
# paste0(round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
#        " & (", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(cancer.incidence.crude$composite.cifj.drow$t, cancer.incidence.crude$composite.cifj.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")

####################################################################################
### CIF curves plot
## all-cause mortality
load("mortality/natural_mortality.RData")

natural.mortality$cifj.driptw <- as.data.frame(apply(natural.mortality$cifj.driptw, 2, diff))
natural.mortality$cifj.driptw[, !names(natural.mortality$cifj.driptw) %in% c("t")] <- natural.mortality$cifj.driptw[, !names(natural.mortality$cifj.driptw) %in% c("t")]/natural.mortality$cifj.driptw$t
natural.mortality$cifj.driptw$t <- natural.mortality$cifj.driptw$t[1:(length(natural.mortality$cifj.driptw$t)-1)]/365.25

natural.mortality$cifj.drow <- as.data.frame(apply(natural.mortality$cifj.drow, 2, diff))
natural.mortality$cifj.drow[, !names(natural.mortality$cifj.drow) %in% c("t")] <- natural.mortality$cifj.drow[, !names(natural.mortality$cifj.drow) %in% c("t")]/natural.mortality$cifj.drow$t
natural.mortality$cifj.drow$t <- natural.mortality$cifj.drow$t[1:(length(natural.mortality$cifj.drow$t)-1)]/365.25

# ATE
natural.mortality.ate.diff <- ggplot(natural.mortality$cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)            ATE, all-cause mortality") +
  theme_classic()

# ATO
natural.mortality.ato.diff <- ggplot(natural.mortality$cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)            OWATE, all-cause mortality") +
  theme_classic()

ggarrange(natural.mortality.ate.diff, natural.mortality.ato.diff,
          nrow = 1, ncol = 2)
ggsave("figure/mortality/natural_mortality_cif_effect.pdf", width=8, height=4)

## cancer incidence
load("cancer incidence/cancer_incidence.RData")

cancer.incidence$controlled.direct.cifj.driptw <- as.data.frame(apply(cancer.incidence$controlled.direct.cifj.driptw, 2, diff))
cancer.incidence$controlled.direct.cifj.driptw[, !names(cancer.incidence$controlled.direct.cifj.driptw) %in% c("t")] <- cancer.incidence$controlled.direct.cifj.driptw[, !names(cancer.incidence$controlled.direct.cifj.driptw) %in% c("t")]/cancer.incidence$controlled.direct.cifj.driptw$t
cancer.incidence$controlled.direct.cifj.driptw$t <- cancer.incidence$controlled.direct.cifj.driptw$t[1:(length(cancer.incidence$controlled.direct.cifj.driptw$t)-1)]/365.25

cancer.incidence$composite.cifj.driptw <- as.data.frame(apply(cancer.incidence$composite.cifj.driptw, 2, diff))
cancer.incidence$composite.cifj.driptw[, !names(cancer.incidence$composite.cifj.driptw) %in% c("t")] <- cancer.incidence$composite.cifj.driptw[, !names(cancer.incidence$composite.cifj.driptw) %in% c("t")]/cancer.incidence$composite.cifj.driptw$t
cancer.incidence$composite.cifj.driptw$t <- cancer.incidence$composite.cifj.driptw$t[1:(length(cancer.incidence$composite.cifj.driptw$t)-1)]/365.25

cancer.incidence$total.cifj.driptw <- as.data.frame(apply(cancer.incidence$total.cifj.driptw, 2, diff))
cancer.incidence$total.cifj.driptw[, !names(cancer.incidence$total.cifj.driptw) %in% c("t")] <- cancer.incidence$total.cifj.driptw[, !names(cancer.incidence$total.cifj.driptw) %in% c("t")]/cancer.incidence$total.cifj.driptw$t
cancer.incidence$total.cifj.driptw$t <- cancer.incidence$total.cifj.driptw$t[1:(length(cancer.incidence$total.cifj.driptw$t)-1)]/365.25

cancer.incidence$sep.direct.cifj.driptw <- as.data.frame(apply(cancer.incidence$sep.direct.cifj.driptw, 2, diff))
cancer.incidence$sep.direct.cifj.driptw[, !names(cancer.incidence$sep.direct.cifj.driptw) %in% c("t")] <- cancer.incidence$sep.direct.cifj.driptw[, !names(cancer.incidence$sep.direct.cifj.driptw) %in% c("t")]/cancer.incidence$sep.direct.cifj.driptw$t
cancer.incidence$sep.direct.cifj.driptw$t <- cancer.incidence$sep.direct.cifj.driptw$t[1:(length(cancer.incidence$sep.direct.cifj.driptw$t)-1)]/365.25

cancer.incidence$sep.indirect.cifj.driptw <- as.data.frame(apply(cancer.incidence$sep.indirect.cifj.driptw, 2, diff))
cancer.incidence$sep.indirect.cifj.driptw[, !names(cancer.incidence$sep.indirect.cifj.driptw) %in% c("t")] <- cancer.incidence$sep.indirect.cifj.driptw[, !names(cancer.incidence$sep.indirect.cifj.driptw) %in% c("t")]/cancer.incidence$sep.indirect.cifj.driptw$t
cancer.incidence$sep.indirect.cifj.driptw$t <- cancer.incidence$sep.indirect.cifj.driptw$t[1:(length(cancer.incidence$sep.indirect.cifj.driptw$t)-1)]/365.25

cancer.incidence$controlled.direct.cifj.drow <- as.data.frame(apply(cancer.incidence$controlled.direct.cifj.drow, 2, diff))
cancer.incidence$controlled.direct.cifj.drow[, !names(cancer.incidence$controlled.direct.cifj.drow) %in% c("t")] <- cancer.incidence$controlled.direct.cifj.drow[, !names(cancer.incidence$controlled.direct.cifj.drow) %in% c("t")]/cancer.incidence$controlled.direct.cifj.drow$t
cancer.incidence$controlled.direct.cifj.drow$t <- cancer.incidence$controlled.direct.cifj.drow$t[1:(length(cancer.incidence$controlled.direct.cifj.drow$t)-1)]/365.25

cancer.incidence$composite.cifj.drow <- as.data.frame(apply(cancer.incidence$composite.cifj.drow, 2, diff))
cancer.incidence$composite.cifj.drow[, !names(cancer.incidence$composite.cifj.drow) %in% c("t")] <- cancer.incidence$composite.cifj.drow[, !names(cancer.incidence$composite.cifj.drow) %in% c("t")]/cancer.incidence$composite.cifj.drow$t
cancer.incidence$composite.cifj.drow$t <- cancer.incidence$composite.cifj.drow$t[1:(length(cancer.incidence$composite.cifj.drow$t)-1)]/365.25

cancer.incidence$total.cifj.drow <- as.data.frame(apply(cancer.incidence$total.cifj.drow, 2, diff))
cancer.incidence$total.cifj.drow[, !names(cancer.incidence$total.cifj.drow) %in% c("t")] <- cancer.incidence$total.cifj.drow[, !names(cancer.incidence$total.cifj.drow) %in% c("t")]/cancer.incidence$total.cifj.drow$t
cancer.incidence$total.cifj.drow$t <- cancer.incidence$total.cifj.drow$t[1:(length(cancer.incidence$total.cifj.drow$t)-1)]/365.25

cancer.incidence$sep.direct.cifj.drow <- as.data.frame(apply(cancer.incidence$sep.direct.cifj.drow, 2, diff))
cancer.incidence$sep.direct.cifj.drow[, !names(cancer.incidence$sep.direct.cifj.drow) %in% c("t")] <- cancer.incidence$sep.direct.cifj.drow[, !names(cancer.incidence$sep.direct.cifj.drow) %in% c("t")]/cancer.incidence$sep.direct.cifj.drow$t
cancer.incidence$sep.direct.cifj.drow$t <- cancer.incidence$sep.direct.cifj.drow$t[1:(length(cancer.incidence$sep.direct.cifj.drow$t)-1)]/365.25

cancer.incidence$sep.indirect.cifj.drow <- as.data.frame(apply(cancer.incidence$sep.indirect.cifj.drow, 2, diff))
cancer.incidence$sep.indirect.cifj.drow[, !names(cancer.incidence$sep.indirect.cifj.drow) %in% c("t")] <- cancer.incidence$sep.indirect.cifj.drow[, !names(cancer.incidence$sep.indirect.cifj.drow) %in% c("t")]/cancer.incidence$sep.indirect.cifj.drow$t
cancer.incidence$sep.indirect.cifj.drow$t <- cancer.incidence$sep.indirect.cifj.drow$t[1:(length(cancer.incidence$sep.indirect.cifj.drow$t)-1)]/365.25

### multiple causes of failure
## total
# ATE
total.ate.diff <- ggplot(cancer.incidence$total.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)               ATE, total effect") +
  theme_classic()

# ATO
total.ato.diff <- ggplot(cancer.incidence$total.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)               OWATE, total effect") +
  theme_classic()

## sep direct
# ATE
sep.direct.ate.diff <- ggplot(cancer.incidence$sep.direct.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)         ATE, separable direct effect") +
  theme_classic()

# ATO
sep.direct.ato.diff <- ggplot(cancer.incidence$sep.direct.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)         OWATE, separable direct effect") +
  theme_classic()

## sep indirect
# ATE
sep.indirect.ate.diff <- ggplot(cancer.incidence$sep.indirect.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)       ATE, separable indirect effect") +
  theme_classic()

# ATO
sep.indirect.ato.diff <- ggplot(cancer.incidence$sep.indirect.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       OWATE, separable indirect effect") +
  theme_classic()

### one cause of failure
## controlled.direct
# ATE
controlled.direct.ate.diff <- ggplot(cancer.incidence$controlled.direct.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)       ATE, controlled direct effect") +
  theme_classic()

# ATO
controlled.direct.ato.diff <- ggplot(cancer.incidence$controlled.direct.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       OWATE, controlled direct effect") +
  theme_classic()

### mortality included
## composite
# ATE
composite.ate.diff <- ggplot(cancer.incidence$composite.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)                 ATE, composite") +
  theme_classic()

# ATO
composite.ato.diff <- ggplot(cancer.incidence$composite.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.2, 0.05), breaks = c(-0.2, -0.15, -0.1, -0.05, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)                 OWATE, composite") +
  theme_classic()

## cancer mortality
load("mortality/cancer_mortality.RData")

cancer.mortality$cifj.driptw <- as.data.frame(apply(cancer.mortality$cifj.driptw, 2, diff))
cancer.mortality$cifj.driptw[, !names(cancer.mortality$cifj.driptw) %in% c("t")] <- cancer.mortality$cifj.driptw[, !names(cancer.mortality$cifj.driptw) %in% c("t")]/cancer.mortality$cifj.driptw$t
cancer.mortality$cifj.driptw$t <- cancer.mortality$cifj.driptw$t[1:(length(cancer.mortality$cifj.driptw$t)-1)]/365.25

cancer.mortality$cifj.drow <- as.data.frame(apply(cancer.mortality$cifj.drow, 2, diff))
cancer.mortality$cifj.drow[, !names(cancer.mortality$cifj.drow) %in% c("t")] <- cancer.mortality$cifj.drow[, !names(cancer.mortality$cifj.drow) %in% c("t")]/cancer.mortality$cifj.drow$t
cancer.mortality$cifj.drow$t <- cancer.mortality$cifj.drow$t[1:(length(cancer.mortality$cifj.drow$t)-1)]/365.25

# ATE
cancer.mortality.ate.diff <- ggplot(cancer.mortality$cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       ATE, mortality after cancer") +
  theme_classic()

# ATO
cancer.mortality.ato.diff <- ggplot(cancer.mortality$cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       OWATE, mortality after cancer") +
  theme_classic()

## composite
# ATE
composite.ate.diff <- ggplot(cancer.incidence$composite.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)         ATE, composite outcome") +
  theme_classic()

# ATO
composite.ato.diff <- ggplot(cancer.incidence$composite.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)         OWATE, composite outcome") +
  theme_classic()

ggarrange(total.ate.diff, total.ato.diff,
          sep.direct.ate.diff, sep.direct.ato.diff,
          sep.indirect.ate.diff, sep.indirect.ato.diff,
          controlled.direct.ate.diff, controlled.direct.ato.diff,
          nrow = 4, ncol = 2)
ggsave("figure/incidence/cancer_incidence_cif_effect.pdf", width=8, height=16)

ggarrange(cancer.mortality.ate.diff, cancer.mortality.ato.diff,
          composite.ate.diff, composite.ato.diff,
          nrow = 2, ncol = 2)
ggsave("figure/mortality/cancer_mortality_cif_effect.pdf", width=8, height=8)

## breast cancer incidence
load("breast cancer incidence/breast_cancer_incidence.RData")

breast.cancer.incidence$controlled.direct.cifj.driptw <- as.data.frame(apply(breast.cancer.incidence$controlled.direct.cifj.driptw, 2, diff))
breast.cancer.incidence$controlled.direct.cifj.driptw[, !names(breast.cancer.incidence$controlled.direct.cifj.driptw) %in% c("t")] <- breast.cancer.incidence$controlled.direct.cifj.driptw[, !names(breast.cancer.incidence$controlled.direct.cifj.driptw) %in% c("t")]/breast.cancer.incidence$controlled.direct.cifj.driptw$t
breast.cancer.incidence$controlled.direct.cifj.driptw$t <- breast.cancer.incidence$controlled.direct.cifj.driptw$t[1:(length(breast.cancer.incidence$controlled.direct.cifj.driptw$t)-1)]/365.25

breast.cancer.incidence$composite.cifj.driptw <- as.data.frame(apply(breast.cancer.incidence$composite.cifj.driptw, 2, diff))
breast.cancer.incidence$composite.cifj.driptw[, !names(breast.cancer.incidence$composite.cifj.driptw) %in% c("t")] <- breast.cancer.incidence$composite.cifj.driptw[, !names(breast.cancer.incidence$composite.cifj.driptw) %in% c("t")]/breast.cancer.incidence$composite.cifj.driptw$t
breast.cancer.incidence$composite.cifj.driptw$t <- breast.cancer.incidence$composite.cifj.driptw$t[1:(length(breast.cancer.incidence$composite.cifj.driptw$t)-1)]/365.25

breast.cancer.incidence$total.cifj.driptw <- as.data.frame(apply(breast.cancer.incidence$total.cifj.driptw, 2, diff))
breast.cancer.incidence$total.cifj.driptw[, !names(breast.cancer.incidence$total.cifj.driptw) %in% c("t")] <- breast.cancer.incidence$total.cifj.driptw[, !names(breast.cancer.incidence$total.cifj.driptw) %in% c("t")]/breast.cancer.incidence$total.cifj.driptw$t
breast.cancer.incidence$total.cifj.driptw$t <- breast.cancer.incidence$total.cifj.driptw$t[1:(length(breast.cancer.incidence$total.cifj.driptw$t)-1)]/365.25

breast.cancer.incidence$sep.direct.cifj.driptw <- as.data.frame(apply(breast.cancer.incidence$sep.direct.cifj.driptw, 2, diff))
breast.cancer.incidence$sep.direct.cifj.driptw[, !names(breast.cancer.incidence$sep.direct.cifj.driptw) %in% c("t")] <- breast.cancer.incidence$sep.direct.cifj.driptw[, !names(breast.cancer.incidence$sep.direct.cifj.driptw) %in% c("t")]/breast.cancer.incidence$sep.direct.cifj.driptw$t
breast.cancer.incidence$sep.direct.cifj.driptw$t <- breast.cancer.incidence$sep.direct.cifj.driptw$t[1:(length(breast.cancer.incidence$sep.direct.cifj.driptw$t)-1)]/365.25

breast.cancer.incidence$sep.indirect.cifj.driptw <- as.data.frame(apply(breast.cancer.incidence$sep.indirect.cifj.driptw, 2, diff))
breast.cancer.incidence$sep.indirect.cifj.driptw[, !names(breast.cancer.incidence$sep.indirect.cifj.driptw) %in% c("t")] <- breast.cancer.incidence$sep.indirect.cifj.driptw[, !names(breast.cancer.incidence$sep.indirect.cifj.driptw) %in% c("t")]/breast.cancer.incidence$sep.indirect.cifj.driptw$t
breast.cancer.incidence$sep.indirect.cifj.driptw$t <- breast.cancer.incidence$sep.indirect.cifj.driptw$t[1:(length(breast.cancer.incidence$sep.indirect.cifj.driptw$t)-1)]/365.25

breast.cancer.incidence$controlled.direct.cifj.drow <- as.data.frame(apply(breast.cancer.incidence$controlled.direct.cifj.drow, 2, diff))
breast.cancer.incidence$controlled.direct.cifj.drow[, !names(breast.cancer.incidence$controlled.direct.cifj.drow) %in% c("t")] <- breast.cancer.incidence$controlled.direct.cifj.drow[, !names(breast.cancer.incidence$controlled.direct.cifj.drow) %in% c("t")]/breast.cancer.incidence$controlled.direct.cifj.drow$t
breast.cancer.incidence$controlled.direct.cifj.drow$t <- breast.cancer.incidence$controlled.direct.cifj.drow$t[1:(length(breast.cancer.incidence$controlled.direct.cifj.drow$t)-1)]/365.25

breast.cancer.incidence$composite.cifj.drow <- as.data.frame(apply(breast.cancer.incidence$composite.cifj.drow, 2, diff))
breast.cancer.incidence$composite.cifj.drow[, !names(breast.cancer.incidence$composite.cifj.drow) %in% c("t")] <- breast.cancer.incidence$composite.cifj.drow[, !names(breast.cancer.incidence$composite.cifj.drow) %in% c("t")]/breast.cancer.incidence$composite.cifj.drow$t
breast.cancer.incidence$composite.cifj.drow$t <- breast.cancer.incidence$composite.cifj.drow$t[1:(length(breast.cancer.incidence$composite.cifj.drow$t)-1)]/365.25

breast.cancer.incidence$total.cifj.drow <- as.data.frame(apply(breast.cancer.incidence$total.cifj.drow, 2, diff))
breast.cancer.incidence$total.cifj.drow[, !names(breast.cancer.incidence$total.cifj.drow) %in% c("t")] <- breast.cancer.incidence$total.cifj.drow[, !names(breast.cancer.incidence$total.cifj.drow) %in% c("t")]/breast.cancer.incidence$total.cifj.drow$t
breast.cancer.incidence$total.cifj.drow$t <- breast.cancer.incidence$total.cifj.drow$t[1:(length(breast.cancer.incidence$total.cifj.drow$t)-1)]/365.25

breast.cancer.incidence$sep.direct.cifj.drow <- as.data.frame(apply(breast.cancer.incidence$sep.direct.cifj.drow, 2, diff))
breast.cancer.incidence$sep.direct.cifj.drow[, !names(breast.cancer.incidence$sep.direct.cifj.drow) %in% c("t")] <- breast.cancer.incidence$sep.direct.cifj.drow[, !names(breast.cancer.incidence$sep.direct.cifj.drow) %in% c("t")]/breast.cancer.incidence$sep.direct.cifj.drow$t
breast.cancer.incidence$sep.direct.cifj.drow$t <- breast.cancer.incidence$sep.direct.cifj.drow$t[1:(length(breast.cancer.incidence$sep.direct.cifj.drow$t)-1)]/365.25

breast.cancer.incidence$sep.indirect.cifj.drow <- as.data.frame(apply(breast.cancer.incidence$sep.indirect.cifj.drow, 2, diff))
breast.cancer.incidence$sep.indirect.cifj.drow[, !names(breast.cancer.incidence$sep.indirect.cifj.drow) %in% c("t")] <- breast.cancer.incidence$sep.indirect.cifj.drow[, !names(breast.cancer.incidence$sep.indirect.cifj.drow) %in% c("t")]/breast.cancer.incidence$sep.indirect.cifj.drow$t
breast.cancer.incidence$sep.indirect.cifj.drow$t <- breast.cancer.incidence$sep.indirect.cifj.drow$t[1:(length(breast.cancer.incidence$sep.indirect.cifj.drow$t)-1)]/365.25

### multiple causes of failure
## total
# ATE
total.ate.diff <- ggplot(breast.cancer.incidence$total.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)               ATE, total effect") +
  theme_classic()

# ATO
total.ato.diff <- ggplot(breast.cancer.incidence$total.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)               OWATE, total effect") +
  theme_classic()

## sep direct
# ATE
sep.direct.ate.diff <- ggplot(breast.cancer.incidence$sep.direct.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)         ATE, separable direct effect") +
  theme_classic()

# ATO
sep.direct.ato.diff <- ggplot(breast.cancer.incidence$sep.direct.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)         OWATE, separable direct effect") +
  theme_classic()

## sep indirect
# ATE
sep.indirect.ate.diff <- ggplot(breast.cancer.incidence$sep.indirect.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)       ATE, separable indirect effect") +
  theme_classic()

# ATO
sep.indirect.ato.diff <- ggplot(breast.cancer.incidence$sep.indirect.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       OWATE, separable indirect effect") +
  theme_classic()

### one cause of failure
## controlled.direct
# ATE
controlled.direct.ate.diff <- ggplot(breast.cancer.incidence$controlled.direct.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)       ATE, controlled direct effect") +
  theme_classic()

# ATO
controlled.direct.ato.diff <- ggplot(breast.cancer.incidence$controlled.direct.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       OWATE, controlled direct effect") +
  theme_classic()

### mortality included
## composite
# ATE
composite.ate.diff <- ggplot(breast.cancer.incidence$composite.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a)                 ATE, composite") +
  theme_classic()

# ATO
composite.ato.diff <- ggplot(breast.cancer.incidence$composite.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.1, 0.05), breaks = c(-0.1, -0.075, -0.05, -0.025, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)                 OWATE, composite") +
  theme_classic()

## cancer mortality
load("mortality/breast_cancer_mortality.RData")

breast.cancer.mortality$cifj.driptw <- as.data.frame(apply(breast.cancer.mortality$cifj.driptw, 2, diff))
breast.cancer.mortality$cifj.driptw[, !names(breast.cancer.mortality$cifj.driptw) %in% c("t")] <- breast.cancer.mortality$cifj.driptw[, !names(breast.cancer.mortality$cifj.driptw) %in% c("t")]/breast.cancer.mortality$cifj.driptw$t
breast.cancer.mortality$cifj.driptw$t <- breast.cancer.mortality$cifj.driptw$t[1:(length(breast.cancer.mortality$cifj.driptw$t)-1)]/365.25

breast.cancer.mortality$cifj.drow <- as.data.frame(apply(breast.cancer.mortality$cifj.drow, 2, diff))
breast.cancer.mortality$cifj.drow[, !names(breast.cancer.mortality$cifj.drow) %in% c("t")] <- breast.cancer.mortality$cifj.drow[, !names(breast.cancer.mortality$cifj.drow) %in% c("t")]/breast.cancer.mortality$cifj.drow$t
breast.cancer.mortality$cifj.drow$t <- breast.cancer.mortality$cifj.drow$t[1:(length(breast.cancer.mortality$cifj.drow$t)-1)]/365.25

# ATE
breast.cancer.mortality.ate.diff <- ggplot(breast.cancer.mortality$cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       ATE, mortality after cancer") +
  theme_classic()

# ATO
breast.cancer.mortality.ato.diff <- ggplot(breast.cancer.mortality$cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)       OWATE, mortality after cancer") +
  theme_classic()

## composite
# ATE
composite.ate.diff <- ggplot(breast.cancer.incidence$composite.cifj.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)         ATE, composite outcome") +
  theme_classic()

# ATO
composite.ato.diff <- ggplot(breast.cancer.incidence$composite.cifj.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(-0.3, 0), breaks = c(-0.3, -0.2, -0.1, 0)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in CIF") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b)         OWATE, composite outcome") +
  theme_classic()

ggarrange(total.ate.diff, total.ato.diff,
          sep.direct.ate.diff, sep.direct.ato.diff,
          sep.indirect.ate.diff, sep.indirect.ato.diff,
          controlled.direct.ate.diff, controlled.direct.ato.diff,
          nrow = 4, ncol = 2)
ggsave("figure/incidence/breast_cancer_incidence_cif_effect.pdf", width=8, height=16)

ggarrange(breast.cancer.mortality.ate.diff, breast.cancer.mortality.ato.diff,
          composite.ate.diff, composite.ato.diff,
          nrow = 2, ncol = 2)
ggsave("figure/mortality/breast_cancer_mortality_cif_effect.pdf", width=8, height=8)





sum(cancer.incidence$total.cifj.drow$estimate.a0.upper.cb<cancer.incidence$total.cifj.drow$estimate.a0.upper.ci)
sum(cancer.incidence$total.cifj.drow$estimate.a1.upper.cb<cancer.incidence$total.cifj.drow$estimate.a1.upper.ci)

sum(cancer.incidence$total.cifj.drow$estimate.a0.lower.cb>cancer.incidence$total.cifj.drow$estimate.a0.lower.ci)
sum(cancer.incidence$total.cifj.drow$estimate.a1.lower.cb>cancer.incidence$total.cifj.drow$estimate.a1.lower.ci)

sum(cancer.incidence$total.cifj.driptw$estimate.a0.upper.cb<cancer.incidence$total.cifj.driptw$estimate.a0.upper.ci)
sum(cancer.incidence$total.cifj.driptw$estimate.a1.upper.cb<cancer.incidence$total.cifj.driptw$estimate.a1.upper.ci)

sum(cancer.incidence$total.cifj.driptw$estimate.a0.lower.cb>cancer.incidence$total.cifj.driptw$estimate.a0.lower.ci)
sum(cancer.incidence$total.cifj.driptw$estimate.a1.lower.cb>cancer.incidence$total.cifj.driptw$estimate.a1.lower.ci)

## cancer mortality
load("cancer_mortality_rv.RData")
f.sim$dr.direct.mc.ha.driptw <- f.sim$dr.direct.mc.ha.driptw/365.25
f.sim$dr.direct.mc.ha.drow <- f.sim$dr.direct.mc.ha.drow/365.25
# ATE
cancer.mortality.ate.diff <- ggplot(f.sim$dr.direct.mc.ha.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(c)                    ATE, survival") +
  theme_classic()

# ATO
cancer.mortality.ato.diff <- ggplot(f.sim$dr.direct.mc.ha.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(d)                    OWATE, survival") +
  theme_classic()

## treatment-specific
# IPTW
cancer.mortality.ate.a0 <- ggplot(f.sim$dr.direct.mc.ha.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a0), color="black") +
  geom_line(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  xlab("Follow-up time (years)") +
  ylab("RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(e)") +
  theme_classic()

cancer.mortality.ate.a1 <- ggplot(f.sim$dr.direct.mc.ha.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a1), color="black") +
  geom_line(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  xlab("Follow-up time (years)") +
  ylab("RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(f)") +
  theme_classic()

# OW
cancer.mortality.ato.a0 <- ggplot(f.sim$dr.direct.mc.ha.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a0), color="black") +
  geom_line(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  xlab("Follow-up time (years)") +
  ylab("RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(g)") +
  theme_classic()

cancer.mortality.ato.a1 <- ggplot(f.sim$dr.direct.mc.ha.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a1), color="black") +
  geom_line(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  xlab("Follow-up time (years)") +
  ylab("RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(h)") +
  theme_classic()

## arrange

ggarrange(cancer.incidence.ate.a0, cancer.incidence.ate.a1,
          cancer.incidence.ato.a0, cancer.incidence.ato.a1,
          cancer.mortality.ate.a0, cancer.mortality.ate.a1,
          cancer.mortality.ato.a0, cancer.mortality.ato.a1, nrow = 4, ncol = 2)
ggsave("app_treatment_specific.pdf", width=8, height=16)



ggarrange(cancer.incidence.ate.a0, cancer.incidence.ate.a1, cancer.incidence.ate.diff,
          cancer.incidence.ato.a0, cancer.incidence.ato.a1, cancer.incidence.ato.diff,
          cancer.mortality.ate.a0, cancer.mortality.ate.a1, cancer.mortality.ate.diff,
          cancer.mortality.ato.a0, cancer.mortality.ato.a1, cancer.mortality.ato.diff,
          nrow = 4, ncol = 3)
ggsave("application.pdf", width=12, height=16)

### table
## cancer incidence
load("cancer_incidence_rv.RData")
f.sim$dr.direct.mc.ha.driptw <- f.sim$dr.direct.mc.ha.driptw/365.25
f.sim$dr.direct.mc.ha.drow <- f.sim$dr.direct.mc.ha.drow/365.25
# ATE
paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")

# ATO
paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")

## cancer mortality
load("cancer_mortality_rv.RData")
f.sim$dr.direct.mc.ha.driptw <- f.sim$dr.direct.mc.ha.driptw/365.25
f.sim$dr.direct.mc.ha.drow <- f.sim$dr.direct.mc.ha.drow/365.25

# ATE
paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")

# ATO
paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")

### density plot
# cancer incidence
a.sl.lib <- c("SL.mean","SL.glm","SL.nnet","SL.rpartPrune","SL.xgboost","SL.ranger","SL.step","SL.glmnet","SL.earth") # a.sl.lib <- c("SL.mean","SL.glm") # ,"SL.kernelKnn","SL.gam"
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

####################################################################################
## cancer incidence curve
data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(cancer_eventdate)|(prescription_date<cancer_eventdate)) %>%
  filter(!(cancer==1 & is.na(cancer_eventdate))) %>%
  mutate(event1 = cancer) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(cancer_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(cancer, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & cancer==0), 2, cancer)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = cancer) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(cancer), cancer_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n, na.rm=TRUE))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(cancer), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(cancer), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n, na.rm=TRUE))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|cancer)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]
# data.df <- data.df[1:2000,]
time.grid <- 1
admin.cens <- max(data.df$end.time3)
n <- nrow(data.df)
data.df$a <- data.df$ref_sulf
# data.df$end.time <- data.df$end.time3
# data.df$event <- data.df$event3
data.df$id <- data.df$patid

## cancer incidence cumulative risk of censoring
cancer.incidence.G <- data.frame(t=sort(unique(data.df$end.time1))/365.25, cif.a0=1-summary(survfit(Surv(end.time1[data.df$a==0], 1*(event1[data.df$a==0]==0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$surv,
                    cif.a1=1-summary(survfit(Surv(end.time1[data.df$a==1], 1*(event1[data.df$a==1]==0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$surv)

cancer.incidence.G.plot <- ggplot(cancer.incidence.G, aes(x=t)) +
  geom_step(aes(y=cif.a1, linetype="solid"), color="black") +
  geom_step(aes(y=cif.a0, linetype="dotted"), color="black") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1))+
  xlab("Follow-up time (years)") +
  ylab("Cumulative risk of censoring") +
  # geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a) Cancer incidence") +
  theme_classic() +
  theme(legend.position = "none")

###
train.ind <- sample(1:n, ceiling(n/2), replace=FALSE) # sample training set
train.df <- data.df[train.ind,]
train.df <- train.df[order(train.df$patid),]
estimation.df <- data.df[setdiff(1:n, train.ind),]
estimation.df <- estimation.df[order(estimation.df$patid),] # estimation.df <- estimation.df[order(estimation.df$end.time,-estimation.df$event),]

train.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = estimation.df$a, X = estimation.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
estimation.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = train.df$a, X = train.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred) # estimation.df$ps <- predict(glm(a~x1+x2+x3+x4,data=estimation.df,family=binomial()),type="response")
data.df <- rbind(train.df, estimation.df)
# cancer.incidence.ps <- ggplot(data.df, aes(x=ps, group=factor(a)))+
#   geom_density(aes(linetype=ifelse(a, "dotted", "solid")))+
#   xlab("Propensity score")+
#   ylab("Density")+
#   scale_x_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1))+
#   labs(fill = "A:")+
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   guides(fill = guide_legend(title.position = "left"))+
#   scale_y_continuous(limits=c(0, 25), breaks=c(0, 5, 10, 15, 20, 25))+
#   ggtitle("(a) Cancer incidence")

# data.df <- cancer.incidence.hist$data

cancer.incidence.hist <- ggplot(data.df, aes(ps, fill=factor(a))) +
  geom_histogram(binwidth=0.05, position="identity", alpha = 0.6, col=I("black")) +
  xlab("Propensity score")+
  ylab("Frequency")+
  ylim(0, 35000)+
  scale_fill_manual(values=c("grey60", "grey20"))+
    theme(legend.position = "none", plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("(a) Cancer incidence")

# cancer mortality
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>%
  filter(prescription_date > cancer_eventdate) %>%
  mutate(event1 = death) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]
# data.df <- data.df[1:2000,]
time.grid <- 1
admin.cens <- max(data.df$end.time1)
n <- nrow(data.df)
data.df$a <- data.df$ref_sulf
data.df$end.time <- data.df$end.time1
data.df$event <- data.df$event1
data.df$id <- data.df$patid

## cancer mortality cumulative risk of censoring
cancer.mortality.G <- data.frame(t=sort(unique(data.df$end.time1))/365.25, cif.a0=1-summary(survfit(Surv(end.time1[data.df$a==0], 1*(event1[data.df$a==0]==0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$surv,
                                 cif.a1=1-summary(survfit(Surv(end.time1[data.df$a==1], 1*(event1[data.df$a==1]==0)) ~ 1, data=data.df), times=sort(unique(data.df$end.time1)), extend=TRUE)$surv)

cancer.mortality.G.plot <- ggplot(cancer.mortality.G, aes(x=t)) +
  geom_step(aes(y=cif.a1, linetype="solid"), color="black") +
  geom_step(aes(y=cif.a0, linetype="dotted"), color="black") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/3, 2*followup/3, followup)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1))+
  xlab("Follow-up time (years)") +
  ylab("Cumulative risk of censoring") +
  # geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b) Mortality after cancer diagnosis") +
  theme_classic() +
  theme(legend.position = "none")

ggarrange(cancer.incidence.G.plot, cancer.mortality.G.plot, nrow = 1, ncol = 2)
ggsave("censor_cif.pdf", width=8, height=4)

# aje
ggarrange(cancer.incidence.G.plot, nrow = 1, ncol = 1)
ggsave("censor_cif_a.pdf", width=4, height=4)

ggarrange(cancer.incidence.G.plot, nrow = 1, ncol = 1)
ggsave("censor_cif_b.pdf", width=4, height=4)

###
train.ind <- sample(1:n, ceiling(n/2), replace=FALSE) # sample training set
train.df <- data.df[train.ind,]
train.df <- train.df[order(train.df$patid),]
estimation.df <- data.df[setdiff(1:n, train.ind),]
estimation.df <- estimation.df[order(estimation.df$patid),] # estimation.df <- estimation.df[order(estimation.df$end.time,-estimation.df$event),]

train.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = estimation.df$a, X = estimation.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
estimation.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = train.df$a, X = train.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred) # estimation.df$ps <- predict(glm(a~x1+x2+x3+x4,data=estimation.df,family=binomial()),type="response")

data.df <- rbind(train.df, estimation.df)
# cancer.mortality.ps <- ggplot(data.df, aes(x=ps, group=factor(a)))+
#   geom_density(aes(linetype=ifelse(a, "dotted", "solid")))+
#   xlab("Propensity score")+
#   ylab("Density")+
#   scale_x_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1))+
#   labs(fill = "A:")+
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   guides(fill = guide_legend(title.position = "left"))+
#   scale_y_continuous(limits=c(0, 25), breaks=c(0, 5, 10, 15, 20, 25))+
#   ggtitle("(b) Mortality after cancer diagnosis")

data.df <- cancer.mortality.hist$data

cancer.mortality.hist <- ggplot(data.df, aes(ps, fill=factor(a))) +
  geom_histogram(binwidth=0.05, position="identity", alpha = 0.6, col=I("black")) +
  xlab("Propensity score")+
  ylab("Frequency")+
  ylim(0, 3500)+
  scale_fill_manual(values=c("grey60", "grey20"))+
  theme(legend.position = "none", plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("(b) Mortality after cancer diagnosis")

ggarrange(cancer.incidence.hist, cancer.mortality.hist, nrow = 1, ncol = 2)
ggsave("application_hist.pdf", width=8, height=4)

# aje
cancer.incidence.hist
ggsave("application_hist_a.pdf", width=4, height=4)

cancer.mortality.hist
ggsave("application_hist_b.pdf", width=4, height=4)

