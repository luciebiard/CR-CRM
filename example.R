load('list_scenarios.Rdata')
source('functions.R')

# Nsim=1
# n=70
# correct=2
# K=5
# n1=35
# scenar="plateau2_bis"
# tstar=8
# acc_typ='fixslow'
# delta=0.10
# graine=1234
# rando="adapt"
# typ="exponential"
# half=0.06

xx <- get_nsim_surv(Nsim=1, n=70, correct=2, K=5,n1=35, scenar="plateau2_bis", tstar=8, acc_typ='fixslow', delta=0.10, graine=1234, rando="adapt", typ="exponential", half=0.06)
design <- xx[1:2] # trial simulations: trial datasets and summarized performance metrics
bck <- xx[3:6] # benchmark: latent CI DLT, observed CI DLT, latent CI prog, observed CI prog

