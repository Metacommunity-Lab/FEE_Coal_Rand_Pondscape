library(dplyr); library(tibble); library(ggplot2); library(gridExtra)

# Charging functoins 
source("FEE_Random_landscape_Ponderful.R")
source("FEE_IterationModel_Ponderful.R")

load("Espacio_Retromed.RData")

#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
# Albera network iterations _______________________________________________________________####

Albera<-cbind(ALBERA,100)
Albera<-Albera[,-4]
colnames(Albera)[8]<-"J"
Albera<-Albera[,c(1:3, 5,7,6,4,8)]
M.distance.Albera<-as.matrix(dist(Albera[,2:3]))
Meta.comm.Albera<-Meta.comm.Albera
J.Albera<-Albera[,8]
D50s<-10^seq(log10(10),log10(3000),,50)
Gr<-Albera$grado
Cc<-Albera$cc
Bt<-Albera$bc
N<-nrow(M.distance.Albera)

# Model iterations
out.d50.obs.landscape.Albera<-NULL
for (d in D50s){
  out.t<-iteration.model(replicas=112, Meta.pool=Meta.comm.Albera, m.pool=0.01, Js=J.Albera, id.module=Albera[,7],
                         M.dist=M.distance.Albera, D50=d, m.max=1,
                         id.fixed=NULL, D50.fixed=1000, m.max.fixed=1, comm.fixed=Meta.comm.Albera,
                         Lottery=F, it=2000, prop.dead.by.it=0.05, id.obs=1:ncol(M.distance.Albera))
  out.t<-resume.out(out.t)
  S<-out.t[[1]][11:(10+N)]
  Be.all<-out.t[[1]][(10+N+1):(10+N+N)]
  gamma<-out.t[[1]][10]
  print(c(length(d),length(S), length(Be.all), length(gamma),length(Bt), length(Gr), length(Cc)))
  out.t2<-cbind(d,S, Be.all, gamma,Bt, Gr, Cc)
  
  print(tail(out.t2))
  out.d50.obs.landscape.Albera<-rbind(out.d50.obs.landscape.Albera, out.t2)
}

out.d50.obs.landscape.Albera
save.image("~/Dropbox/H2020/RETROMED/Frontiers_Diciembre_2022/Frontiers_Respuestas/Retromed_H2020_febrero_2023.RData")

# Albera NULL MODEL _______________________________________________________________####
# sequence of D50
random_vs_real_D50.Albera_Febrero_2023<-H2020_Random_Landscape_local_ponds.D50(n.random = 140,Meta.pool=Meta.comm.Albera, m.pool=0.01,  Js=Albera[,8],
                                                                               id.module=Albera[,7],
                                                                               D50.min=1, D50.max=2000, m.max=1,M.X.Y=Albera[,2:3],
                                                                               id.fixed=0, D50.fixed=100, m.max.fixed=1, comm.fixed=Meta.comm.Albera,
                                                                               Lottery=F, it=0, prop.dead.by.it=0, id.obs=1:nrow(Albera))

#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
# Giara network iterations _______________________________________________________________####

Giara<-cbind(GIARA,100)
Giara<-Giara[,-4]
colnames(Giara)[8]<-"J"
Giara<-Giara[,c(1:3, 5,7,6,4,8)]
M.distance.Giara<-as.matrix(dist(Giara[,2:3]))
Meta.comm.Giara<-Meta.comm
J.Giara<-Giara[,8]
#D50s<-seq(50,3000,,50)
D50s<-10^seq(log10(10),log10(3000),,50)
Gr<-Giara$grado
Cc<-Giara$cc
Bt<-Giara$bc
N<-nrow(M.distance.Giara)

# Model iterations
out.d50.obs.landscape.Giara<-NULL
for (d in D50s){
  out.t<-iteration.model(replicas=112, Meta.pool=Meta.comm.Giara, m.pool=0.01, Js=J.Giara, id.module=Giara[,7],
                         M.dist=M.distance.Giara, D50=d, m.max=1,
                         id.fixed=NULL, D50.fixed=1000, m.max.fixed=1, comm.fixed=Meta.comm.Giara,
                         Lottery=F, it=2000, prop.dead.by.it=0.05, id.obs=1:ncol(M.distance.Giara))
  out.t<-resume.out(out.t)
  S<-out.t[[1]][11:(10+N)]
  Be.all<-out.t[[1]][(10+N+1):(10+N+N)]
  gamma<-out.t[[1]][10]
  print(c(length(d),length(S), length(Be.all), length(gamma),length(Bt), length(Gr), length(Cc)))
  out.t2<-cbind(d,S, Be.all, gamma,Bt, Gr, Cc)
  
  print(tail(out.t2))
  out.d50.obs.landscape.Giara<-rbind(out.d50.obs.landscape.Giara, out.t2)
}

out.d50.obs.landscape.Giara

# Giara NULL MODEL_______________________________________________________________####
# sequence of D50
random_vs_real_D50.Giara_febrero_2023<-H2020_Random_Landscape_local_ponds.D50(n.random = 140,Meta.pool=Meta.comm.Giara, m.pool=0.01,  Js=Giara[,8],
                                                                              id.module=Giara[,7],
                                                                              D50.min=1, D50.max=2000, m.max=1,M.X.Y=Giara[,2:3],
                                                                              id.fixed=0, D50.fixed=100, m.max.fixed=1, comm.fixed=Meta.comm.Giara,
                                                                              Lottery=F, it=0, prop.dead.by.it=0, id.obs=1:nrow(Giara))

#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
# Guils network iterations _______________________________________________________________####
Guils<-cbind(GUILS,100)
Guils<-Guils[,-4]
colnames(Guils)[8]<-"J"
Guils<-Guils[,c(1:3, 5,7,6,4,8)]
M.distance.Guils<-as.matrix(dist(Guils[,2:3]))
Meta.comm.Guils<-Meta.comm.Albera
J.Guils<-Guils[,8]
#D50s<-seq(50,3000,,50)
D50s<-10^seq(log10(10),log10(3000),,50)
Gr<-Guils$grado
Cc<-Guils$cc
Bt<-Guils$bc
N<-nrow(M.distance.Guils)

# Model iterations
out.d50.obs.landscape.Guils<-NULL
for (d in D50s){
  out.t<-iteration.model(replicas=112, Meta.pool=Meta.comm.Guils, m.pool=0.01, Js=J.Guils, id.module=Guils[,7],
                         M.dist=M.distance.Guils, D50=d, m.max=1,
                         id.fixed=NULL, D50.fixed=1000, m.max.fixed=1, comm.fixed=Meta.comm.Guils,
                         Lottery=F, it=2000, prop.dead.by.it=0.05, id.obs=1:ncol(M.distance.Guils))
  out.t<-resume.out(out.t)
  S<-out.t[[1]][11:(10+N)]
  Be.all<-out.t[[1]][(10+N+1):(10+N+N)]
  gamma<-out.t[[1]][10]
  print(c(length(d),length(S), length(Be.all), length(gamma),length(Bt), length(Gr), length(Cc)))
  out.t2<-cbind(d,S, Be.all, gamma,Bt, Gr, Cc)
  
  print(tail(out.t2))
  out.d50.obs.landscape.Guils<-rbind(out.d50.obs.landscape.Guils, out.t2)
}

out.d50.obs.landscape.Guils

# Guils NULL MODEL _______________________________________________________________####
# sequence of D50
random_vs_real_D50.Guils_febrero_2023<-H2020_Random_Landscape_local_ponds.D50(n.random = 1000,Meta.pool=Meta.comm.Guils, m.pool=0.01,  Js=Guils[,8],
                                                                                     id.module=Guils[,7],
                                                                                     D50.min=1, D50.max=2000, m.max=1,M.X.Y=Guils[,2:3],
                                                                                     id.fixed=0, D50.fixed=100, m.max.fixed=1, comm.fixed=Meta.comm.Guils,
                                                                                     Lottery=F, it=0, prop.dead.by.it=0, id.obs=1:nrow(Guils))

#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
# PNAE network iterations _______________________________________________________________####
Pnae<-cbind(PNAE,100)
Pnae<-Pnae[,-4]
colnames(Pnae)[8]<-"J"
Pnae<-Pnae[,c(1:3, 5,7,6,4,8)]
M.distance.Pnae<-as.matrix(dist(Pnae[,2:3]))
Meta.comm.Pnae<-Meta.comm
J.Pnae<-Pnae[,8]
#D50s<-seq(50,3000,,50)
D50s<-10^seq(log10(10),log10(3000),,50)
Gr<-Pnae$grado
Cc<-Pnae$cc
Bt<-Pnae$bc
N<-nrow(M.distance.Pnae)

# Model iterations
out.d50.obs.landscape.Pnae<-NULL
for (d in D50s){
  out.t<-iteration.model(replicas=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01, Js=J.Pnae, id.module=Pnae[,7],
                         M.dist=M.distance.Pnae, D50=d, m.max=1,
                         id.fixed=NULL, D50.fixed=1000, m.max.fixed=1, comm.fixed=Meta.comm.Pnae,
                         Lottery=F, it=2000, prop.dead.by.it=0.05, id.obs=1:ncol(M.distance.Pnae))
  out.t<-resume.out(out.t)
  S<-out.t[[1]][11:(10+N)]
  Be.all<-out.t[[1]][(10+N+1):(10+N+N)]
  gamma<-out.t[[1]][10]
  print(c(length(d),length(S), length(Be.all), length(gamma),length(Bt), length(Gr), length(Cc)))
  out.t2<-cbind(d,S, Be.all, gamma,Bt, Gr, Cc)
  
  print(tail(out.t2))
  out.d50.obs.landscape.Pnae<-rbind(out.d50.obs.landscape.Pnae, out.t2)
}

out.d50.obs.landscape.Pnae

# PNAE NULL MODEL _______________________________________________________________####
# sequence of D50
random_vs_real_D50.Pnae_Febrero_2023<-H2020_Random_Landscape_local_ponds.D50(n.random = 140,Meta.pool=Meta.comm.Pnae, m.pool=0.01,  Js=Pnae[,8],
                                                                id.module=Pnae[,7],
                                                                D50.min=1, D50.max=2000, m.max=1,M.X.Y=Pnae[,2:3],
                                                                id.fixed=0, D50.fixed=100, m.max.fixed=1, comm.fixed=Meta.comm.Pnae,
                                                                Lottery=F, it=0, prop.dead.by.it=0, id.obs=1:nrow(Pnae))

#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
# Vila nova network iterations ____________________________________________________________####
Vilavona<-cbind(VILAVONA,100)
Vilavona<-Vilavona[,-4]
colnames(Vilavona)[8]<-"J"
Vilavona<-Vilavona[,c(1:3, 5,7,6,4,8)]
M.distance.Vilavona<-as.matrix(dist(Vilavona[,2:3]))
Meta.comm.Vilavona<-Meta.comm.Albera
J.Vilavona<-Vilavona[,8]
#D50s<-seq(50,3000,,50)
D50s<-10^seq(log10(10),log10(3000),,50)
Gr<-Vilavona$grado
Cc<-Vilavona$cc
Bt<-Vilavona$bc
N<-nrow(M.distance.Vilavona)

# Model iterations
out.d50.obs.landscape.Vilavona<-NULL
for (d in D50s){
  out.t<-iteration.model(replicas=112, Meta.pool=Meta.comm.Vilavona, m.pool=0.01, Js=J.Vilavona, id.module=Vilavona[,7],
                         M.dist=M.distance.Vilavona, D50=d, m.max=1,
                         id.fixed=NULL, D50.fixed=1000, m.max.fixed=1, comm.fixed=Meta.comm.Vilavona,
                         Lottery=F, it=2000, prop.dead.by.it=0.05, id.obs=1:ncol(M.distance.Vilavona))
  out.t<-resume.out(out.t)
  S<-out.t[[1]][11:(10+N)]
  Be.all<-out.t[[1]][(10+N+1):(10+N+N)]
  gamma<-out.t[[1]][10]
  print(c(length(d),length(S), length(Be.all), length(gamma),length(Bt), length(Gr), length(Cc)))
  out.t2<-cbind(d,S, Be.all, gamma,Bt, Gr, Cc)
  
  print(tail(out.t2))
  out.d50.obs.landscape.Vilavona<-rbind(out.d50.obs.landscape.Vilavona, out.t2)
}

out.d50.obs.landscape.Vilavona

# Vila Nova NULL MODEL ____________________________________________________________####
# sequence of D50
random_vs_real_D50.Vilavona_Febrero_2023<-H2020_Random_Landscape_local_ponds.D50(n.random = 140,Meta.pool=Meta.comm.Vilavona, m.pool=0.01,  Js=Vilavona[,8],
                                                                    id.module=Vilavona[,7],
                                                                    D50.min=1, D50.max=2000, m.max=1,M.X.Y=Vilavona[,2:3],
                                                                    id.fixed=0, D50.fixed=100, m.max.fixed=1, comm.fixed=Meta.comm.Vilavona,
                                                                    Lottery=F, it=0, prop.dead.by.it=0, id.obs=1:nrow(Vilavona))

#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
#_##############################################################################################
# Rocha network iterations ____________________________________________________________####
CHARCOS_ROCHA <- read.delim("~/Dropbox/H2020/RETROMED/CHARCOS_ROCHA.txt")
Rocha<-CHARCOS_ROCHA
Rocha <- cbind(Rocha,0,1,100)

#Rocha<-cbind(Rocha,100)
#Rocha<-Rocha[,-4]
colnames(Rocha)[6:8]<-c("bc","modules", "J")
Rocha[1:10,7]<-2
#Rocha<-Rocha[,c(1:3, 5,7,6,4,8)]
M.distance.Rocha<-as.matrix(dist(Rocha[,2:3]))
Meta.comm.Rocha<-Meta.comm.Albera
J.Rocha<-Rocha[,8]
#D50s<-seq(50,3000,,50)
D50s<-10^seq(log10(10),log10(3000),,50)
Gr<-Rocha$grado
Cc<-Rocha$cc
Bt<-Rocha$bc
N<-nrow(M.distance.Rocha)

# Model iterations
out.d50.obs.landscape.Rocha<-NULL
for (d in D50s){
  out.t<-iteration.model(replicas=112, Meta.pool=Meta.comm.Rocha, m.pool=0.01, Js=J.Rocha, id.module=Rocha[,7],
                         M.dist=M.distance.Rocha, D50=d, m.max=1,
                         id.fixed=NULL, D50.fixed=1000, m.max.fixed=1, comm.fixed=Meta.comm.Rocha,
                         Lottery=F, it=2000, prop.dead.by.it=0.05, id.obs=1:ncol(M.distance.Rocha))
  out.t<-resume.out(out.t)
  S<-out.t[[1]][11:(10+N)]
  Be.all<-out.t[[1]][(10+N+1):(10+N+N)]
  gamma<-out.t[[1]][10]
  print(c(length(d),length(S), length(Be.all), length(gamma),length(Bt), length(Gr), length(Cc)))
  out.t2<-cbind(d,S, Be.all, gamma,Bt, Gr, Cc)
  
  print(tail(out.t2))
  out.d50.obs.landscape.Rocha<-rbind(out.d50.obs.landscape.Rocha, out.t2)
}

out.d50.obs.landscape.Rocha

# Rocha NULL MODEL ____________________________________________________________####
# sequence of D50
random_vs_real_D50.Rocha_febrero2023<-H2020_Random_Landscape_local_ponds.D50(n.random = 140,Meta.pool=Meta.comm.Rocha, m.pool=0.01,  Js=Rocha[,8],
                                                                 id.module=Rocha[,7],
                                                                 D50.min=1, D50.max=2000, m.max=1,M.X.Y=Rocha[,2:3],
                                                                 id.fixed=0, D50.fixed=100, m.max.fixed=1, comm.fixed=Meta.comm.Rocha,
                                                                 Lottery=F, it=0, prop.dead.by.it=0, id.obs=1:nrow(Rocha))





