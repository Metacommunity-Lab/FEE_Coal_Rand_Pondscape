library(cowplot)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

# Figure 1
# Geographic distances
library(geosphere) # Charging package to caluclate distances
library(sp) # Other packages to deal with geospatial data 
library(stars) # Other packages to deal with geospatial data 
library(ggspatial) # Other packages to deal with geospatial data 

# Plots to obtain a gg network plot
library(GGally)
library(network)
library(ggnetwork)
library(dplyr)
library(ggplot2)
library(ggforce)
library(gtable)    
library(grid)
library(gridExtra) 
library(png)
library(viridis)
library(cowplot)
library(magick)

# Calculation of shortest path
library(ape)
xy.ALB <- read.table(file = "C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/00. URUGUAY/Ponderful/RETRO_H2020Pond/Figure 1/CHARCOS_ALBERA.txt",
                     sep = "\t", dec =",",header = T)
xy.GIA <- read.table(file = "C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/00. URUGUAY/Ponderful/RETRO_H2020Pond/Figure 1/CHARCOS_GIARA.txt",
                     sep = "\t", dec =",",header = T)
xy.GUI <- read.table(file = "C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/00. URUGUAY/Ponderful/RETRO_H2020Pond/Figure 1/CHARCOS_GUILS.txt",
                     sep = "\t", dec =",",header = T)
xy.PNAE <- read.table(file = "C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/00. URUGUAY/Ponderful/RETRO_H2020Pond/Figure 1/CHARCOS_PNAE.txt",
                      sep = "\t", dec =",",header = T)
xy.VMIL <- read.table(file = "C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/00. URUGUAY/Ponderful/RETRO_H2020Pond/Figure 1/CHARCOS_VILANOVAMILFONTES.txt",
                      sep = "\t", dec =",",header = T)
xy.ROC <- read.table(file = "C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/00. URUGUAY/Ponderful/RETRO_H2020Pond/Figure 1/CHARCOS_ROCHA.txt",
                     sep = "\t", dec =".",header = T)
# Check structure and correct PNAE
str(xy.ALB)
str(xy.GIA)
str(xy.GUI)# Remember that they are repeated
xy.GUI<- xy.GUI[-which(duplicated(xy.GUI[,2:3])==T),]

str(xy.PNAE)# Remember that they are repeated
xy.PNAE<- xy.PNAE[-which(duplicated(xy.PNAE[,2:3])==T),]

str(xy.VMIL)
str(xy.ROC)
c(nrow(xy.ALB),nrow(xy.GIA),nrow(xy.GUI),nrow(xy.PNAE),nrow(xy.VMIL),nrow(xy.ROC))
c(542,160,86,255,146)

# Creation of minimum spanning tree (one link per pond)
GG_ALB <- mst(as.matrix(dist(xy.ALB)))
GG_GIA <- mst(as.matrix(dist(xy.GIA)))
GG_GUI <- mst(as.matrix(dist(xy.GUI)))
GG_PNAE <- mst(as.matrix(dist(xy.PNAE)))
GG_VMIL <- mst(as.matrix(dist(xy.VMIL)))
GG_ROC <- mst(as.matrix(dist(cbind(xy.ROC[,2],xy.ROC[,3]))))


# List creation
xy.values <- list(xy.ALB[,2:3],xy.PNAE[,2:3],xy.GIA[,2:3],xy.VMIL[,1:2], xy.GUI[,2:3], xy.ROC[,2:3])
GG.values <- list(GG_ALB,GG_PNAE, GG_GIA, GG_VMIL, GG_GUI, GG_ROC)
Pondscapes_names <- c("Albera","PNAE","Giara", "Vila Nova","Guils", "Rocha")

# Nice colours? CUNILLERA_palette is what you need
source("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")

library(sna)
gg_plot_network <- list()
for (u in 1:length(xy.values)) {
  detach("package:sna", unload = TRUE)
  library(igraph)
  # In order to plot the network in a ggplot way do the follwing
  weight_degree <- Pondscapes_outputs[[u]][which(Pondscapes_outputs[[u]][,1]==max_diff[u]),7]
  
  detach("package:igraph", unload = TRUE)
  library("sna")
  n<- network(GG.values[[u]], directed=F, diag=F) 
  n %v% "family" <- weight_degree # Family is an standard name... 
  
  # ALL THESE LINES MUST BE RUNNED!!!
  #from here_____________________________________________________________________
  
  gg_plot_network[[u]] <- ggplot(n, layout=as.matrix(xy.values[[u]]),
                                 aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(color = "grey30", size=0.5) +
    geom_nodes(aes(fill=family, size=family) ,color ="black" ,shape=21, alpha=.55)+
    scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
    labs(x="",y="",title = Pondscapes_names[[u]])+
    theme(axis.line = element_blank(),
          plot.title = element_text(size=35, colour = "black",hjust = 0.5),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background = element_blank())
  
  legend_plots<- get_legend(ggplot(n, layout=as.matrix(xy.values[[u]]),
                                   aes(x = x, y = y, xend = xend, yend = yend))+
                              geom_edges(color = "black", size=0.5) +
                              geom_nodes(aes(fill=family) ,color ="grey50" ,shape=21, alpha=.55)+
                              scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F, name="Weighted degree",
                                                   guide = guide_colorbar(ticks = T,
                                                                          barheight = 5,
                                                                          barwidth = 1))+
                              theme_classic()+
                              theme(legend.direction = "vertical",
                                    legend.box="vertical"))
  
}

png(filename ="C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/00. URUGUAY/Ponderful/RETRO_H2020Pond/Figure 1/Figure_1_Pondscapes_MST.png",
    width =940*5 ,height =600*5 ,units ="px",res = 300)
grid.arrange(gg_plot_network[[1]],gg_plot_network[[2]],gg_plot_network[[3]],legend_plots,
             gg_plot_network[[4]],gg_plot_network[[5]],gg_plot_network[[6]],
             ncol=4,nrow=2)
dev.off()


# Figure 3
load("Retromed_H2020_febrero_2023.RData")

# Nice colours? CUNILLERA_palette is what you need
source("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")

Pondscapes_outputs <- list(out.d50.obs.landscape.Albera,
                           out.d50.obs.landscape.Pnae,
                           out.d50.obs.landscape.Giara,
                           out.d50.obs.landscape.Vilavona,
                           out.d50.obs.landscape.Guils,
                           out.d50.obs.landscape.Rocha) 

Pondscapes__rand_outputs <- list(random_vs_real_D50.Albera_Febrero_2023,
                                 random_vs_real_D50.Pnae_Febrero_2023,
                                 random_vs_real_D50.Giara_febrero_2023,
                                 random_vs_real_D50.Vilavona_Febrero_2023,
                                 random_vs_real_D50.Guils_Febrero_2023,
                                 random_vs_real_D50.Rocha_febrero2023) 

Pondscapes_names <- c("Albera", "PNAE", "Giara", "Vila Nova", "Guils","Rocha")

final_plot_v1 <- list()
final_plot_v2 <- list()
supp_CV_plot <- list()

#Albera________________________________________________________________________________________________________####
p=1

a<- Pondscapes_outputs[[p]]
############# Figures Albera
#########################################3
### Plot results: alpha and beta diversity in species` dispersal gradient
### also considering the gradient in local ponds isolation

# Calculation of the minumum and maximum quantiles 
ii.min.g<-which(a[,6]<quantile(a[,6],.05))
ii.max.g<-which(a[,6]>quantile(a[,6],.95))


# Minimum quantile line
d<-a[ii.min.g,1]
S<-a[ii.min.g,2]
gg<-a[ii.min.g,6]
cc<-a[ii.min.g,7]
B<-a[ii.min.g,3]

# Maximum quantile line
d2<-a[ii.max.g,1]
S2<-a[ii.max.g,2]
gg2<-a[ii.max.g,6]
cc2<-a[ii.max.g,7]
B2<-a[ii.max.g,3]

### ALPHA - S
#Function to descrive the minimum curve
M.hill<-nls(S~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=4, Smax=30, q=3, d50=150))
summary(M.hill)
ph<-coefficients(M.hill)
f.hill<-function(x) ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3])

#Function to descrive the maximum curve
M.hill.2<-nls(S2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=4, Smax=32, q=3, d50=100))
summary(M.hill.2)
ph2<-coefficients(M.hill.2)
f.hill.2<-function(x) ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3])

# Differential delta between the two lines and small plot of it
f.hill.delta<-function(x)(ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3]))/(ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3]))
S.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Alpha ratio",subtitle = "Alpha ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))


# Plot for alpha 
S.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>% ggplot(aes(d,S))+
            #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
            geom_jitter(aes(fill=Gr,colour=Gr),height =0.2 ,shape=21, alpha=0.2,size=2)+
            scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
            scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
            labs(title="A)",subtitle = Pondscapes_names[p], y="Alpha diversity", x=NULL)+
  
            geom_function(fun=f.hill, col="white",  size=1.8)+ 
            geom_function(fun=f.hill, col="black", linetype="dashed", size=1.5)+
            
            geom_function(fun=f.hill.2, col="white",  size=1.8)+ 
            geom_function(fun=f.hill.2, col="black", linetype="solid", size=1.5)+
  
            xlab("Dispersal (d50)")+
  
            scale_x_log10()+theme_bw()+theme(legend.position = "none")

S.o.Albera<-S.o.Albera+annotation_custom(ggplotGrob(S.ratio), xmin = log10(40), xmax = log10(250),  ymin = 20, ymax = 37)

#### Beta
#Function to descrive the minimum curve
B.M.hill<-nls(B~B0+(Bmax-B0)*(d^q)/(d50^q+d^q), start = list(B0=4.396e-01, Bmax=9.179e-01, q=-2.56, d50=500))
summary(B.M.hill)
B.ph<-coefficients(B.M.hill)
B.f.hill<-function(x) B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3])

#Function to descrive the maximum curve
B.M.hill.2<-nls(B2~B0+(Bmax-B0)*(d2^q)/(d50^q+d2^q), start = list(B0=0, Bmax=10, q=-3, d50=150))
summary(B.M.hill.2)
B.ph2<-coefficients(B.M.hill.2)
B.f.hill.2<-function(x) B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3])

# Differential delta between the two lines and small plot of it
B.f.hill.delta<-function(x)(B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3]))/(B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3]))
B.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=B.f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Beta ratio",subtitle = "Beta ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))


# Plot for beta 
B.o.Albera<-data.frame(a) %>% arrange(desc(d)) %>% ggplot()+aes(d,Be.all)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
   geom_jitter(aes(fill=Gr,colour=Gr),height =0.02 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(title="B)",subtitle = " ", y="Beta diversity", x=NULL)+
  
  geom_function(fun=B.f.hill, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=B.f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

B.o.Albera<-B.o.Albera+annotation_custom(ggplotGrob(B.ratio), xmin = log10(40), xmax = log10(250),  ymin = .5, ymax = .8)

#Legend extraction of the gradient of Closeness for both alpha and beta (it is the same for both)
legend <- get_legend(ggplot(data.frame(a))+aes(d,Be.all, fill=Gr)+
                       geom_point(shape=21, alpha=0.1,size=3, colour="black")+
                       scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F, name="Degree")+
                       scale_x_log10()+theme_bw())

# Obtention of Z-plots responding to the differential between alpha diversity and beta diversity
#Check function "plotea.null.landscape.D50" for more information
b<-data.frame(Pondscapes__rand_outputs[[p]])
b<-mutate(b, Z.S=(mean.S.real.-mean.S.null.)/sd.S.null.)
b$Z.S[which(b$Z.S==Inf)] <- 0
b<-mutate(b, Z.Bet=(mean.Be.all.real.-mean.Be.all.null.)/sd.Be.all.null.)
b<-mutate(b, CV.oe.S=(sd.S.real./mean.S.real.-sd.S.null./mean.S.null.))
b<-mutate(b, CV.oe.Bet=(sd.Be.all.real./mean.Be.all.real.-sd.Be.all.null./mean.Be.all.null.))
b<-mutate(b, Z.Gama=(Gamma.real-Gamma.null)/sd(Gamma.null))

S.oe<-ggplot(b, aes(D50,Z.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  labs(y="Z-Alpha", x="Dispersal (d50)", size=0.25, title = "C)", subtitle = " ")+
  geom_hline(yintercept=0, linetype="dashed")+scale_x_log10()+
  theme_bw()
B.oe<-ggplot(b, aes(D50,Z.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  labs(y="Z-Beta", x="Dispersal (d50)", size=0.25, title="D)", subtitle = " ")+
  theme_bw()+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")


S.cv.oe<-ggplot(b, aes(D50,CV.oe.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  ylab("CV alpha real- CV alpha null")+
  xlab("Dispersal (d50)")+
  labs(title=Pondscapes_names[p])+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()
B.cv.oe<-ggplot(b, aes(D50,CV.oe.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  ylab("CV beta real- CV beta null")+
  xlab("Dispersal (d50)")+
  labs(title=" ")+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()

supp_CV_plot[[p]]<- grid.arrange(S.cv.oe,B.cv.oe, ncol=2)

output_Figure3_data <- data.frame()
output_Figure3_data <- output_Figure3_data%>%bind_rows(
data.frame("Site" = "Albera",
"Max_RatioAlpha"=max(f.hill.delta(d)),
"Dist_Max_RatioAlpha"=unique(d[which(f.hill.delta(d)==max(f.hill.delta(d)))]),
"Min_RatioBeta"=min(B.f.hill.delta(d)),
"Dist_Min_RatioBeta"=unique(d[which(B.f.hill.delta(d)==min(B.f.hill.delta(d)))]),
"Max_Z_Alpha"=max(b$Z.S),
"Dist_Max_Z_Alpha"=b$D50[which(b$Z.S==max(b$Z.S))],
"Min_Z_Alpha"=min(b$Z.S),
"Dist_Min_Z_Alpha"=b$D50[which(b$Z.S==min(b$Z.S))],
"Max_Z_Beta"=max(b$Z.Bet),
"Dist_Max_Z_Beta"=b$D50[which(b$Z.Bet==max(b$Z.Bet))],
"Min_Z_Beta"=min(b$Z.Bet),
"Dist_Min_Z_Beta"=b$D50[which(b$Z.Bet==min(b$Z.Bet))]
))

# Final plot assembly
final_plot_v1[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.oe,      x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.oe,      x = 0.85,y = 0,   width = 0.15, height = 1)

final_plot_v2[[p]] <- ggdraw() +
              draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
              draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
              draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
              draw_plot(S.cv.oe,   x = 0.7, y = 0,   width = 0.15, height = 1)+
              draw_plot(B.cv.oe,   x = 0.85,y = 0,   width = 0.15, height = 1)



#PNAE________________________________________________________________________________________________________####
p=2

a<- Pondscapes_outputs[[p]]
############# Figures Albera
#########################################3
### Plot results: alpha and beta diversity in species` dispersal gradient
### also considering the gradient in local ponds isolation

# Calculation of the minumum and maximum quantiles 
ii.min.g<-which(a[,6]<quantile(a[,6],.01))
ii.max.g<-which(a[,6]>quantile(a[,6],.99))


# Minimum quantile line
d<-a[ii.min.g,1]
S<-a[ii.min.g,2]
gg<-a[ii.min.g,6]
cc<-a[ii.min.g,7]
B<-a[ii.min.g,3]

# Maximum quantile line
d2<-a[ii.max.g,1]
S2<-a[ii.max.g,2]
gg2<-a[ii.max.g,6]
cc2<-a[ii.max.g,7]
B2<-a[ii.max.g,3]

### ALPHA - S
#Function to descrive the minimum curve
M.hill<-nls(S~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=0, Smax=30, q=3, d50=500))
summary(M.hill)
ph<-coefficients(M.hill)
f.hill<-function(x) ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3])

#Function to descrive the maximum curve
M.hill.2<-nls(S2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=0, Smax=20, q=3, d50=100))
summary(M.hill.2)
ph2<-coefficients(M.hill.2)
f.hill.2<-function(x) ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3])

# Differential delta between the two lines and small plot of it
f.hill.delta<-function(x)(ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3]))/(ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3]))
S.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Alpha ratio",subtitle = "Alpha ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))

# Plot for alpha 
S.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,S)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.2 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(subtitle=Pondscapes_names[p], y="Alpha diversity", x=NULL)+
  
  geom_function(fun=f.hill, col="white",  size=1.8)+ 
  geom_function(fun=f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

S.o.Albera<-S.o.Albera+annotation_custom(ggplotGrob(S.ratio), xmin = log10(520), xmax = log10(2900),  ymin = 1, ymax = 22)

#### Beta
#Function to descrive the minimum curve
B.M.hill<-nls(B~B0+(Bmax-B0)*(d^q)/(d50^q+d^q), start = list(B0=4.396e-01, Bmax=9.179e-01, q=-2.56, d50=250) )
summary(B.M.hill)
B.ph<-coefficients(B.M.hill)
B.f.hill<-function(x) B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3])

#Function to descrive the maximum curve
B.M.hill.2<-nls(B2~B0+(Bmax-B0)*(d2^q)/(d50^q+d2^q), start = list(B0=0, Bmax=10, q=-1, d50=200)   )
summary(B.M.hill.2)
B.ph2<-coefficients(B.M.hill.2)
B.f.hill.2<-function(x) B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3])

# Differential delta between the two lines and small plot of it
B.f.hill.delta<-function(x)(B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3]))/(B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3]))
B.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=B.f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Beta ratio",subtitle = "Beta ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))


# Plot for beta 
B.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,Be.all)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.02 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(subtitle="", y="Beta diversity", x=NULL)+
  
  geom_function(fun=B.f.hill, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=B.f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

B.o.Albera<-B.o.Albera+annotation_custom(ggplotGrob(B.ratio), xmin = log10(510), xmax = log10(3100),  ymin = .7, ymax = .99)

#Legend extraction of the gradient of Closeness for both alpha and beta (it is the same for both)
legend <- get_legend(ggplot(data.frame(a))+aes(d,Be.all, fill=Gr)+
                       geom_point(shape=21, alpha=0.1,size=3, colour="black")+
                       scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F, name="Degree")+
                       scale_x_log10()+theme_bw())

# Obtention of Z-plots responding to the differential between alpha diversity and beta diversity
#Check function "plotea.null.landscape.D50" for more information
b<-data.frame(Pondscapes__rand_outputs[[p]])
b<-mutate(b, Z.S=(mean.S.real.-mean.S.null.)/sd.S.null.)
b$Z.S[which(b$Z.S==Inf)] <- 0
b<-mutate(b, Z.Bet=(mean.Be.all.real.-mean.Be.all.null.)/sd.Be.all.null.)
b<-mutate(b, CV.oe.S=(sd.S.real./mean.S.real.-sd.S.null./mean.S.null.))
b<-mutate(b, CV.oe.Bet=(sd.Be.all.real./mean.Be.all.real.-sd.Be.all.null./mean.Be.all.null.))
b<-mutate(b, Z.Gama=(Gamma.real-Gamma.null)/sd(Gamma.null))

S.oe<-ggplot(b, aes(D50,Z.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  labs(y="Z-Alpha", x="Dispersal (d50)", size=0.25, title = " ")+
  geom_hline(yintercept=0, linetype="dashed")+#scale_x_log10()+
  theme_bw()
B.oe<-ggplot(b, aes(D50,Z.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  labs(y="Z-Beta", x="Dispersal (d50)", size=0.25, title=" ")+
  theme_bw()+#scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")


S.cv.oe<-ggplot(b, aes(D50,CV.oe.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  ylab("CV alpha real- CV alpha null")+
  xlab("Dispersal (d50)")+
  labs(title=Pondscapes_names[p])+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()
B.cv.oe<-ggplot(b, aes(D50,CV.oe.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  ylab("CV beta real- CV beta null")+
  xlab("Dispersal (d50)")+
  labs(title=" ")+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()

supp_CV_plot[[p]]<- grid.arrange(S.cv.oe,B.cv.oe, ncol=2)

output_Figure3_data <- output_Figure3_data%>%bind_rows(
  data.frame("Site" = "PNAE",
             "Max_RatioAlpha"=max(f.hill.delta(d)),
             "Dist_Max_RatioAlpha"=unique(d[which(f.hill.delta(d)==max(f.hill.delta(d)))]),
             "Min_RatioBeta"=min(B.f.hill.delta(d)),
             "Dist_Min_RatioBeta"=unique(d[which(B.f.hill.delta(d)==min(B.f.hill.delta(d)))]),
             "Max_Z_Alpha"=max(b$Z.S),
             "Dist_Max_Z_Alpha"=b$D50[which(b$Z.S==max(b$Z.S))],
             "Min_Z_Alpha"=min(b$Z.S),
             "Dist_Min_Z_Alpha"=b$D50[which(b$Z.S==min(b$Z.S))],
             "Max_Z_Beta"=max(b$Z.Bet),
             "Dist_Max_Z_Beta"=b$D50[which(b$Z.Bet==max(b$Z.Bet))],
             "Min_Z_Beta"=min(b$Z.Bet),
             "Dist_Min_Z_Beta"=b$D50[which(b$Z.Bet==min(b$Z.Bet))]
  ))

# Final plot assembly
final_plot_v1[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.oe,      x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.oe,      x = 0.85,y = 0,   width = 0.15, height = 1)

final_plot_v2[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.cv.oe,   x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.cv.oe,   x = 0.85,y = 0,   width = 0.15, height = 1)



#Giara________________________________________________________________________________________________________####
p=3

a<- Pondscapes_outputs[[p]]
############# Figures Albera
#########################################3
### Plot results: alpha and beta diversity in species` dispersal gradient
### also considering the gradient in local ponds isolation

# Calculation of the minumum and maximum quantiles 
ii.min.g<-which(a[,6]<quantile(a[,6],.05))
ii.max.g<-which(a[,6]>quantile(a[,6],.95))


# Minimum quantile line
d<-a[ii.min.g,1]
S<-a[ii.min.g,2]
gg<-a[ii.min.g,6]
cc<-a[ii.min.g,7]
B<-a[ii.min.g,3]

# Maximum quantile line
d2<-a[ii.max.g,1]
S2<-a[ii.max.g,2]
gg2<-a[ii.max.g,6]
cc2<-a[ii.max.g,7]
B2<-a[ii.max.g,3]

### ALPHA - S
#Function to descrive the minimum curve
M.hill<-nls(S~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=0, Smax=25, q=3, d50=500)   )
summary(M.hill)
ph<-coefficients(M.hill)
f.hill<-function(x) ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3])

#Function to descrive the maximum curve
M.hill.2<-nls(S2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=4, Smax=35, q=1, d50=100)   )
summary(M.hill.2)
ph2<-coefficients(M.hill.2)
f.hill.2<-function(x) ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3])

# Differential delta between the two lines and small plot of it
f.hill.delta<-function(x)(ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3]))/(ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3]))
S.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Alpha ratio",subtitle = "Alpha ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))

# Plot for alpha 
S.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,S)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.2 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(subtitle = Pondscapes_names[p], y="Alpha diversity", x=NULL)+
  
  geom_function(fun=f.hill, col="white",  size=1.8)+ 
  geom_function(fun=f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

S.o.Albera<-S.o.Albera+annotation_custom(ggplotGrob(S.ratio), xmin = log10(900), xmax = log10(3200),  ymin = 0, ymax = 22)

#### Beta
#Function to descrive the minimum curve
B.M.hill<-nls(B~B0+(Bmax-B0)*(d^q)/(d50^q+d^q), start = list(B0=4.396e-01, Bmax=9.179e-01, q=-2.56, d50=500) )
summary(B.M.hill)
B.ph<-coefficients(B.M.hill)
B.f.hill<-function(x) B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3])

#Function to descrive the maximum curve
B.M.hill.2<-nls(B2~B0+(Bmax-B0)*(d2^q)/(d50^q+d2^q), start = list(B0=10, Bmax=1, q=-3, d50=50)   )
summary(B.M.hill.2)
B.ph2<-coefficients(B.M.hill.2)
B.f.hill.2<-function(x) B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3])

# Differential delta between the two lines and small plot of it
B.f.hill.delta<-function(x)(B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3]))/(B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3]))
B.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=B.f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Beta ratio",subtitle = "Beta ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))


# Plot for beta 
B.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,Be.all)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.02 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(subtitle="", y="Beta diversity", x=NULL)+
  
  geom_function(fun=B.f.hill, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=B.f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

B.o.Albera<-B.o.Albera+annotation_custom(ggplotGrob(B.ratio), xmin = log10(40), xmax = log10(250),  ymin = .42, ymax = .74)

#Legend extraction of the gradient of Closeness for both alpha and beta (it is the same for both)
legend <- get_legend(ggplot(data.frame(a))+aes(d,Be.all, fill=Gr)+
                       geom_point(shape=21, alpha=0.1,size=3, colour="black")+
                       scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F, name="Degree")+
                       scale_x_log10()+theme_bw())

# Obtention of Z-plots responding to the differential between alpha diversity and beta diversity
#Check function "plotea.null.landscape.D50" for more information
b<-data.frame(Pondscapes__rand_outputs[[p]])
b<-mutate(b, Z.S=(mean.S.real.-mean.S.null.)/sd.S.null.)
b$Z.S[which(b$Z.S==Inf)] <- 0
b<-mutate(b, Z.Bet=(mean.Be.all.real.-mean.Be.all.null.)/sd.Be.all.null.)
b<-mutate(b, CV.oe.S=(sd.S.real./mean.S.real.-sd.S.null./mean.S.null.))
b<-mutate(b, CV.oe.Bet=(sd.Be.all.real./mean.Be.all.real.-sd.Be.all.null./mean.Be.all.null.))
b<-mutate(b, Z.Gama=(Gamma.real-Gamma.null)/sd(Gamma.null))

S.oe<-ggplot(b, aes(D50,Z.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  labs(y="Z-Alpha", x="Dispersal (d50)", size=0.25, title = " ")+
  geom_hline(yintercept=0, linetype="dashed")+#scale_x_log10()+
  theme_bw()
B.oe<-ggplot(b, aes(D50,Z.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  labs(y="Z-Beta", x="Dispersal (d50)", size=0.25, title=" ")+
  theme_bw()+#scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")


S.cv.oe<-ggplot(b, aes(D50,CV.oe.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  ylab("CV alpha real- CV alpha null")+
  xlab("Dispersal (d50)")+
  labs(title=Pondscapes_names[p])+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()
B.cv.oe<-ggplot(b, aes(D50,CV.oe.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  ylab("CV beta real- CV beta null")+
  xlab("Dispersal (d50)")+
  labs(title=" ")+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()

supp_CV_plot[[p]]<- grid.arrange(S.cv.oe,B.cv.oe, ncol=2)

output_Figure3_data <- output_Figure3_data%>%bind_rows(
  data.frame("Site" = "Giara",
             "Max_RatioAlpha"=max(f.hill.delta(d)),
             "Dist_Max_RatioAlpha"=unique(d[which(f.hill.delta(d)==max(f.hill.delta(d)))]),
             "Min_RatioBeta"=min(B.f.hill.delta(d)),
             "Dist_Min_RatioBeta"=unique(d[which(B.f.hill.delta(d)==min(B.f.hill.delta(d)))]),
             "Max_Z_Alpha"=max(b$Z.S),
             "Dist_Max_Z_Alpha"=b$D50[which(b$Z.S==max(b$Z.S))],
             "Min_Z_Alpha"=min(b$Z.S),
             "Dist_Min_Z_Alpha"=b$D50[which(b$Z.S==min(b$Z.S))],
             "Max_Z_Beta"=max(b$Z.Bet),
             "Dist_Max_Z_Beta"=b$D50[which(b$Z.Bet==max(b$Z.Bet))],
             "Min_Z_Beta"=min(b$Z.Bet),
             "Dist_Min_Z_Beta"=b$D50[which(b$Z.Bet==min(b$Z.Bet))]
  ))

# Final plot assembly
final_plot_v1[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.oe,      x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.oe,      x = 0.85,y = 0,   width = 0.15, height = 1)

final_plot_v2[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.cv.oe,   x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.cv.oe,   x = 0.85,y = 0,   width = 0.15, height = 1)


#Vilanova________________________________________________________________________________________________________####
p=4

a<- Pondscapes_outputs[[p]]
############# Figures Albera
#########################################3
### Plot results: alpha and beta diversity in species` dispersal gradient
### also considering the gradient in local ponds isolation

# Calculation of the minumum and maximum quantiles 
ii.min.g<-which(a[,6]<quantile(a[,6],.05))
ii.max.g<-which(a[,6]>quantile(a[,6],.95))


# Minimum quantile line
d<-a[ii.min.g,1]
S<-a[ii.min.g,2]
gg<-a[ii.min.g,6]
cc<-a[ii.min.g,7]
B<-a[ii.min.g,3]

# Maximum quantile line
d2<-a[ii.max.g,1]
S2<-a[ii.max.g,2]
gg2<-a[ii.max.g,6]
cc2<-a[ii.max.g,7]
B2<-a[ii.max.g,3]

### ALPHA - S
#Function to descrive the minimum curve
M.hill<-nls(S~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=4, Smax=30, q=3, d50=200)   )
summary(M.hill)
ph<-coefficients(M.hill)
f.hill<-function(x) ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3])

#Function to descrive the maximum curve
M.hill.2<-nls(S2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=4, Smax=32, q=3, d50=100)   )
summary(M.hill.2)
ph2<-coefficients(M.hill.2)
f.hill.2<-function(x) ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3])

# Differential delta between the two lines and small plot of it
f.hill.delta<-function(x)(ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3]))/(ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3]))
S.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Alpha ratio",subtitle = "Alpha ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))

# Plot for alpha 
S.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,S)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.2 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(subtitle=Pondscapes_names[p], y="Alpha diversity", x=NULL)+
  
  geom_function(fun=f.hill, col="white",  size=1.8)+ 
  geom_function(fun=f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

S.o.Albera<-S.o.Albera+annotation_custom(ggplotGrob(S.ratio), xmin = log10(400), xmax = log10(3000),  ymin = 1, ymax = 27)

#### BETA
#Function to descrive the minimum curve
B.M.hill<-nls(B~B0+(Bmax-B0)*(d^q)/(d50^q+d^q), start = list(B0=4.396e-01, Bmax=9.179e-01, q=-2.56, d50=250) )
summary(B.M.hill)
B.ph<-coefficients(B.M.hill)
B.f.hill<-function(x) B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3])

#Function to descrive the maximum curve
B.M.hill.2<-nls(B2~B0+(Bmax-B0)*(d2^q)/(d50^q+d2^q), start = list(B0=0, Bmax=10, q=-3, d50=150)   )
summary(B.M.hill.2)
B.ph2<-coefficients(B.M.hill.2)
B.f.hill.2<-function(x) B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3])

# Differential delta between the two lines and small plot of it
B.f.hill.delta<-function(x)(B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3]))/(B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3]))
B.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=B.f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Beta ratio",subtitle = "Beta ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))


# Plot for Beta 
B.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,Be.all)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.02 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(subtitle="", y="Beta diversity", x=NULL)+
  
  geom_function(fun=B.f.hill, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=B.f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

B.o.Albera<-B.o.Albera+annotation_custom(ggplotGrob(B.ratio), xmin = log10(500), xmax = log10(3000),  ymin = .53, ymax = .99)

#Legend extraction of the gradient of Closeness for both alpha and beta (it is the same for both)
legend <- get_legend(ggplot(data.frame(a))+aes(d,Be.all, fill=Gr)+
                       geom_point(shape=21, alpha=0.1,size=3, colour="black")+
                       scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F, name="Degree")+
                       scale_x_log10()+theme_bw())

# Obtention of Z-plots responding to the differential between alpha diversity and beta diversity
#Check function "plotea.null.landscape.D50" for more information
b<-data.frame(Pondscapes__rand_outputs[[p]])
b<-mutate(b, Z.S=(mean.S.real.-mean.S.null.)/sd.S.null.)
b$Z.S[which(b$Z.S==Inf)] <- 0
b<-mutate(b, Z.Bet=(mean.Be.all.real.-mean.Be.all.null.)/sd.Be.all.null.)
b<-mutate(b, CV.oe.S=(sd.S.real./mean.S.real.-sd.S.null./mean.S.null.))
b<-mutate(b, CV.oe.Bet=(sd.Be.all.real./mean.Be.all.real.-sd.Be.all.null./mean.Be.all.null.))
b<-mutate(b, Z.Gama=(Gamma.real-Gamma.null)/sd(Gamma.null))

S.oe<-ggplot(b, aes(D50,Z.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  labs(y="Z-Alpha", x="Dispersal (d50)", size=0.25, title = " ")+
  geom_hline(yintercept=0, linetype="dashed")+scale_x_log10()+
  theme_bw()
B.oe<-ggplot(b, aes(D50,Z.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  labs(y="Z-Beta", x="Dispersal (d50)", size=0.25, title=" ")+
  theme_bw()+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")


S.cv.oe<-ggplot(b, aes(D50,CV.oe.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  ylab("CV alpha real- CV alpha null")+
  xlab("Dispersal (d50)")+
  labs(title=Pondscapes_names[p])+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()
B.cv.oe<-ggplot(b, aes(D50,CV.oe.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  ylab("CV beta real- CV beta null")+
  xlab("Dispersal (d50)")+
  labs(title=" ")+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()

supp_CV_plot[[p]]<- grid.arrange(S.cv.oe,B.cv.oe, ncol=2)

output_Figure3_data <- output_Figure3_data%>%bind_rows(
  data.frame("Site" = "Vila Nova",
             "Max_RatioAlpha"=max(f.hill.delta(d)),
             "Dist_Max_RatioAlpha"=unique(d[which(f.hill.delta(d)==max(f.hill.delta(d)))]),
             "Min_RatioBeta"=min(B.f.hill.delta(d)),
             "Dist_Min_RatioBeta"=unique(d[which(B.f.hill.delta(d)==min(B.f.hill.delta(d)))]),
             "Max_Z_Alpha"=max(b$Z.S),
             "Dist_Max_Z_Alpha"=b$D50[which(b$Z.S==max(b$Z.S))],
             "Min_Z_Alpha"=min(b$Z.S),
             "Dist_Min_Z_Alpha"=b$D50[which(b$Z.S==min(b$Z.S))],
             "Max_Z_Beta"=max(b$Z.Bet),
             "Dist_Max_Z_Beta"=b$D50[which(b$Z.Bet==max(b$Z.Bet))],
             "Min_Z_Beta"=min(b$Z.Bet),
             "Dist_Min_Z_Beta"=b$D50[which(b$Z.Bet==min(b$Z.Bet))]
  ))


# Final plot assembly
final_plot_v1[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.oe,      x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.oe,      x = 0.85,y = 0,   width = 0.15, height = 1)

final_plot_v2[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.cv.oe,   x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.cv.oe,   x = 0.85,y = 0,   width = 0.15, height = 1)



#Guils________________________________________________________________________________________________________####
p=5

a<- Pondscapes_outputs[[p]]
############# Figures Albera
#########################################3
### Plot results: alpha and beta diversity in species` dispersal gradient
### also considering the gradient in local ponds isolation

# Calculation of the minumum and maximum quantiles 
ii.min.g<-which(a[,6]<quantile(a[,6],.05))
ii.max.g<-which(a[,6]>quantile(a[,6],.95))


# Minimum quantile line
d<-a[ii.min.g,1]
S<-a[ii.min.g,2]
gg<-a[ii.min.g,6]
cc<-a[ii.min.g,7]
B<-a[ii.min.g,3]

# Maximum quantile line
d2<-a[ii.max.g,1]
S2<-a[ii.max.g,2]
gg2<-a[ii.max.g,6]
cc2<-a[ii.max.g,7]
B2<-a[ii.max.g,3]

### ALPHA - S
#Function to descrive the minimum curve
M.hill<-nls(S~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=4, Smax=25, q=3, d50=500)   )
summary(M.hill)
ph<-coefficients(M.hill)
f.hill<-function(x) ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3])

#Function to descrive the maximum curve
M.hill.2<-nls(S2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=4, Smax=35, q=1, d50=100)   )
summary(M.hill.2)
ph2<-coefficients(M.hill.2)
f.hill.2<-function(x) ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3])

# Differential delta between the two lines and small plot of it
f.hill.delta<-function(x)(ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3]))/(ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3]))
S.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Alpha ratio",subtitle = "Alpha ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))

# Plot for alpha 
S.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,S)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.2 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(subtitle=Pondscapes_names[p], y="Alpha diversity", x=NULL)+
  
  geom_function(fun=f.hill, col="white",  size=1.8)+ 
  geom_function(fun=f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

S.o.Albera<-S.o.Albera+annotation_custom(ggplotGrob(S.ratio), xmin = log10(600), xmax = log10(3000),  ymin = 0.5, ymax = 22.5)

#### Beta
#Function to descrive the minimum curve
B.M.hill<-nls(B~B0+(Bmax-B0)*(d^q)/(d50^q+d^q), start = list(B0=4.396e-01, Bmax=9.179e-01, q=-2.56, d50=250) )
summary(B.M.hill)
B.ph<-coefficients(B.M.hill)
B.f.hill<-function(x) B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3])

#Function to descrive the maximum curve
B.M.hill.2<-nls(B2~B0+(Bmax-B0)*(d2^q)/(d50^q+d2^q), start = list(B0=0, Bmax=10, q=-1, d50=200)   )
summary(B.M.hill.2)
B.ph2<-coefficients(B.M.hill.2)
B.f.hill.2<-function(x) B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3])

# Differential delta between the two lines and small plot of it
B.f.hill.delta<-function(x)(B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3]))/(B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3]))
B.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=B.f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Beta ratio",subtitle = "Beta ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))


# Plot for beta 
B.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,Be.all)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.02 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(subtitle="", y="Beta diversity", x=NULL)+
  
  geom_function(fun=B.f.hill, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=B.f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

B.o.Albera<-B.o.Albera+annotation_custom(ggplotGrob(B.ratio), xmin = log10(650), xmax = log10(3100),  ymin = .55, ymax = .99)

#Legend extraction of the gradient of Closeness for both alpha and beta (it is the same for both)
legend <- get_legend(ggplot(data.frame(a))+aes(d,Be.all, fill=Gr)+
                       geom_point(shape=21, alpha=0.1,size=3, colour="black")+
                       scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F, name="Degree")+
                       scale_x_log10()+theme_bw())

# Obtention of Z-plots responding to the differential between alpha diversity and beta diversity
#Check function "plotea.null.landscape.D50" for more information
b<-data.frame(Pondscapes__rand_outputs[[p]])
b<-mutate(b, Z.S=(mean.S.real.-mean.S.null.)/sd.S.null.)
b$Z.S[which(b$Z.S==Inf)] <- 0
b<-mutate(b, Z.Bet=(mean.Be.all.real.-mean.Be.all.null.)/sd.Be.all.null.)
b<-mutate(b, CV.oe.S=(sd.S.real./mean.S.real.-sd.S.null./mean.S.null.))
b<-mutate(b, CV.oe.Bet=(sd.Be.all.real./mean.Be.all.real.-sd.Be.all.null./mean.Be.all.null.))
b<-mutate(b, Z.Gama=(Gamma.real-Gamma.null)/sd(Gamma.null))

S.oe<-ggplot(b, aes(D50,Z.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  labs(y="Z-Alpha", x="Dispersal (d50)", size=0.25, title = " ")+
  geom_hline(yintercept=0, linetype="dashed")+scale_x_log10()+
  theme_bw()
B.oe<-ggplot(b, aes(D50,Z.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  labs(y="Z-Beta", x="Dispersal (d50)", size=0.25, title=" ")+
  theme_bw()+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")


S.cv.oe<-ggplot(b, aes(D50,CV.oe.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  ylab("CV alpha real- CV alpha null")+
  xlab("Dispersal (d50)")+
  labs(title=Pondscapes_names[p])+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()
B.cv.oe<-ggplot(b, aes(D50,CV.oe.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  ylab("CV beta real- CV beta null")+
  xlab("Dispersal (d50)")+
  labs(title=" ")+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()

supp_CV_plot[[p]]<- grid.arrange(S.cv.oe,B.cv.oe, ncol=2)

output_Figure3_data <- output_Figure3_data%>%bind_rows(
  data.frame("Site" = "Guils",
             "Max_RatioAlpha"=max(f.hill.delta(d)),
             "Dist_Max_RatioAlpha"=unique(d[which(f.hill.delta(d)==max(f.hill.delta(d)))]),
             "Min_RatioBeta"=min(B.f.hill.delta(d)),
             "Dist_Min_RatioBeta"=unique(d[which(B.f.hill.delta(d)==min(B.f.hill.delta(d)))]),
             "Max_Z_Alpha"=max(b$Z.S),
             "Dist_Max_Z_Alpha"=b$D50[which(b$Z.S==max(b$Z.S))],
             "Min_Z_Alpha"=min(b$Z.S),
             "Dist_Min_Z_Alpha"=b$D50[which(b$Z.S==min(b$Z.S))],
             "Max_Z_Beta"=max(b$Z.Bet),
             "Dist_Max_Z_Beta"=b$D50[which(b$Z.Bet==max(b$Z.Bet))],
             "Min_Z_Beta"=min(b$Z.Bet),
             "Dist_Min_Z_Beta"=b$D50[which(b$Z.Bet==min(b$Z.Bet))]
  ))


# Final plot assembly
final_plot_v1[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.oe,      x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.oe,      x = 0.85,y = 0,   width = 0.15, height = 1)

final_plot_v2[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.cv.oe,   x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.cv.oe,   x = 0.85,y = 0,   width = 0.15, height = 1)


#ROCHA ________________________________________________________________________________________________________####
p=6

a<- Pondscapes_outputs[[p]]
############# Figures Albera
#########################################3
### Plot results: alpha and beta diversity in species` dispersal gradient
### also considering the gradient in local ponds isolation

# Calculation of the minumum and maximum quantiles 
ii.min.g<-which(a[,6]<quantile(a[,6],.05))
ii.max.g<-which(a[,6]>quantile(a[,6],.95))


# Minimum quantile line
d<-a[ii.min.g,1]
S<-a[ii.min.g,2]
gg<-a[ii.min.g,6]
cc<-a[ii.min.g,7]
B<-a[ii.min.g,3]

# Maximum quantile line
d2<-a[ii.max.g,1]
S2<-a[ii.max.g,2]
gg2<-a[ii.max.g,6]
cc2<-a[ii.max.g,7]
B2<-a[ii.max.g,3]

### ALPHA - S
#Function to descrive the minimum curve
M.hill<-nls(S~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=4, Smax=30, q=3, d50=200)   )
summary(M.hill)
ph<-coefficients(M.hill)
f.hill<-function(x) ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3])

#Function to descrive the maximum curve
M.hill.2<-nls(S2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=4, Smax=32, q=3, d50=100)   )
summary(M.hill.2)
ph2<-coefficients(M.hill.2)
f.hill.2<-function(x) ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3])

# Differential delta between the two lines and small plot of it
f.hill.delta<-function(x)(ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3]))/(ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3]))
S.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Alpha ratio",subtitle = "Alpha ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))

# Plot for alpha 
S.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,S)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.2 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  labs(subtitle=Pondscapes_names[p], y="Alpha diversity", x=NULL)+
  
  geom_function(fun=f.hill, col="white",  size=1.8)+ 
  geom_function(fun=f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

S.o.Albera<-S.o.Albera+annotation_custom(ggplotGrob(S.ratio), xmin = log10(300), xmax = log10(2000),  ymin = 2, ymax = 19)

#### Beta
#Function to descrive the minimum curve
B.M.hill<-nls(B~B0+(Bmax-B0)*(d^q)/(d50^q+d^q), start = list(B0=4.396e-01, Bmax=9.179e-01, q=-2.56, d50=250) )
summary(B.M.hill)
B.ph<-coefficients(B.M.hill)
B.f.hill<-function(x) B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3])

#Function to descrive the maximum curve
B.M.hill.2<-nls(B2~B0+(Bmax-B0)*(d2^q)/(d50^q+d2^q), start = list(B0=0, Bmax=10, q=-3, d50=150)   )
summary(B.M.hill.2)
B.ph2<-coefficients(B.M.hill.2)
B.f.hill.2<-function(x) B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3])

# Differential delta between the two lines and small plot of it
B.f.hill.delta<-function(x)(B.ph2[1]+(B.ph2[2]-B.ph2[1])*(x^B.ph2[3])/(B.ph2[4]^B.ph2[3]+x^B.ph2[3]))/(B.ph[1]+(B.ph[2]-B.ph[1])*(x^B.ph[3])/(B.ph[4]^B.ph[3]+x^B.ph[3]))
B.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+
  geom_function(fun=B.f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Beta ratio",subtitle = "Beta ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))


# Plot for beta 
B.o.Albera<-data.frame(a) %>%arrange(desc(d)) %>%ggplot()+aes(d,Be.all)+
  #geom_point(aes(fill=Gr),shape=21, alpha=0.2,size=3, colour="black")+
  geom_jitter(aes(fill=Gr,colour=Gr),height =0.02 ,shape=21, alpha=0.2,size=2)+
  
  scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  scale_color_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F)+
  
  
  labs(subtitle="", y="Beta diversity", x=NULL)+
  
  geom_function(fun=B.f.hill, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill, col="black", linetype="dashed", size=1.5)+
  
  geom_function(fun=B.f.hill.2, col="white",  size=1.8)+ 
  geom_function(fun=B.f.hill.2, col="black", linetype="solid", size=1.5)+
  
  xlab("Dispersal (d50)")+
  
  scale_x_log10()+theme_bw()+theme(legend.position = "none")

B.o.Albera<-B.o.Albera+annotation_custom(ggplotGrob(B.ratio), xmin = log10(500), xmax = log10(2500),  ymin = .5, ymax = .95)

#Legend extraction of the gradient of Closeness for both alpha and beta (it is the same for both)
legend <- get_legend(ggplot(data.frame(a))+aes(d,Be.all, fill=Gr)+
                       geom_point(shape=21, alpha=0.1,size=3, colour="black")+
                       scale_fill_CUNILLERA(palette = "LGTBI", discrete = F, reverse = F, name="Degree")+
                       scale_x_log10()+theme_bw())

# Obtention of Z-plots responding to the differential between alpha diversity and beta diversity
#Check function "plotea.null.landscape.D50" for more information
b<-data.frame(Pondscapes__rand_outputs[[p]])
b<-mutate(b, Z.S=(mean.S.real.-mean.S.null.)/sd.S.null.)
b$Z.S[which(b$Z.S==Inf)] <- 0
b<-mutate(b, Z.Bet=(mean.Be.all.real.-mean.Be.all.null.)/sd.Be.all.null.)
b<-mutate(b, CV.oe.S=(sd.S.real./mean.S.real.-sd.S.null./mean.S.null.))
b<-mutate(b, CV.oe.Bet=(sd.Be.all.real./mean.Be.all.real.-sd.Be.all.null./mean.Be.all.null.))
b<-mutate(b, Z.Gama=(Gamma.real-Gamma.null)/sd(Gamma.null))

S.oe<-ggplot(b, aes(D50,Z.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  labs(y="Z-Alpha", x="Dispersal (d50)", size=0.25, title = " ")+
  geom_hline(yintercept=0, linetype="dashed")+scale_x_log10()+
  theme_bw()
B.oe<-ggplot(b, aes(D50,Z.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  labs(y="Z-Beta", x="Dispersal (d50)", size=0.25, title=" ")+
  theme_bw()+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")


S.cv.oe<-ggplot(b, aes(D50,CV.oe.S)) + geom_point(col="black", fill="grey55", shape=21,size=3)+
  ylab("CV alpha real- CV alpha null")+
  xlab("Dispersal (d50)")+
  labs(title=Pondscapes_names[p])+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()
B.cv.oe<-ggplot(b, aes(D50,CV.oe.Bet)) + geom_point(col="black", fill="grey55", shape=22,size=3)+
  ylab("CV beta real- CV beta null")+
  xlab("Dispersal (d50)")+
  labs(title="")+scale_x_log10()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()

supp_CV_plot[[p]]<- grid.arrange(S.cv.oe,B.cv.oe, ncol=2)

output_Figure3_data <- output_Figure3_data%>%bind_rows(
  data.frame("Site" = "Rocha",
             "Max_RatioAlpha"=max(f.hill.delta(d)),
             "Dist_Max_RatioAlpha"=unique(d[which(f.hill.delta(d)==max(f.hill.delta(d)))]),
             "Min_RatioBeta"=min(B.f.hill.delta(d)),
             "Dist_Min_RatioBeta"=unique(d[which(B.f.hill.delta(d)==min(B.f.hill.delta(d)))]),
             "Max_Z_Alpha"=max(b$Z.S),
             "Dist_Max_Z_Alpha"=b$D50[which(b$Z.S==max(b$Z.S))],
             "Min_Z_Alpha"=min(b$Z.S),
             "Dist_Min_Z_Alpha"=b$D50[which(b$Z.S==min(b$Z.S))],
             "Max_Z_Beta"=max(b$Z.Bet),
             "Dist_Max_Z_Beta"=b$D50[which(b$Z.Bet==max(b$Z.Bet))],
             "Min_Z_Beta"=min(b$Z.Bet),
             "Dist_Min_Z_Beta"=b$D50[which(b$Z.Bet==min(b$Z.Bet))]
  ))


# Final plot assembly
final_plot_v1[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.oe,      x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.oe,      x = 0.85,y = 0,   width = 0.15, height = 1)

final_plot_v2[[p]] <- ggdraw() +
  draw_plot(S.o.Albera,x = 0,   y = 0,   width = 0.3,  height = 1)+
  draw_plot(B.o.Albera,x = 0.3, y = 0,   width = 0.3,  height = 1)+
  draw_plot(legend,    x = 0.6, y = 0.3, width = 0.1,  height = 0.5)+
  draw_plot(S.cv.oe,   x = 0.7, y = 0,   width = 0.15, height = 1)+
  draw_plot(B.cv.oe,   x = 0.85,y = 0,   width = 0.15, height = 1)


# output_Figure3_data printing
save(list ="output_Figure3_data" ,file = "C:/Users/David CM/Dropbox/H2020/RETROMED/Frontiers_Diciembre_2022/Frontiers_Respuestas/output_Figure3_data.RData")

# Image printing
png(filename = "C:/Users/David CM/Dropbox/H2020/RETROMED/Frontiers_Diciembre_2022/Frontiers_Respuestas/Figures_Frontiers/Fig3_PondRetro_Deg_LGTBI.png",
    width =1462*3, height =(length(Pondscapes_outputs)*(320*3)), units = "px",res = 300)
grid.arrange(final_plot_v1[[1]],
             final_plot_v1[[2]],
             final_plot_v1[[3]],
             final_plot_v1[[4]],
             final_plot_v1[[5]],
             final_plot_v1[[6]],
             nrow=6,ncol=1)
dev.off()

png(filename = "C:/Users/David CM/Dropbox/H2020/RETROMED/Frontiers_Diciembre_2022/Frontiers_Respuestas/Figures_Frontiers/Supplementary_Fig3_CV.png",
    width =937*2, height =(length(Pondscapes_outputs)*(398*2)), units = "px",res = 300)
grid.arrange(supp_CV_plot[[1]],
             supp_CV_plot[[2]],
             supp_CV_plot[[3]],
             supp_CV_plot[[4]],
             supp_CV_plot[[5]],
             supp_CV_plot[[6]],
             nrow=6,ncol=1)
dev.off()


# Gammas por cada sistema Supplementary ####

library(tidyverse)

plots_gamma <- list()
for (Pscape in 1:length(Pondscapes_names)) {

plots_gamma[[Pscape]] <- data.frame(Pondscapes_outputs[[Pscape]]) %>% group_by(d) %>% 
    summarise(Mean_S=mean(S), Mean_B=mean(Be.all)) %>% 
    mutate(G=Mean_S*(1+Mean_B)) %>% 
    ggplot(aes(d,G))+
    geom_point(aes(fill=d),shape=21,size=3, colour="black")+
    scale_fill_CUNILLERA(palette = "gespa", discrete = F, reverse = F)+
    labs(subtitle=Pondscapes_names[Pscape], y="Gamma diversity", x=NULL)+
    xlab("Dispersal (d50)")+
    scale_x_log10()+theme_bw()+theme(legend.position = "none")
}

png(filename = "C:/Users/David CM/Dropbox/H2020/RETROMED/Frontiers_Diciembre_2022/Frontiers_Respuestas/Figures_Frontiers/Supplementary_Gamma.png",
    width =400*6, height =250*6, units = "px",res = 300)
gridExtra::grid.arrange(plots_gamma[[1]],plots_gamma[[2]],
                        plots_gamma[[3]],plots_gamma[[4]],
                        plots_gamma[[5]],plots_gamma[[6]], ncol=3)
dev.off()










