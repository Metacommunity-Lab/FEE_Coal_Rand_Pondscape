#########################################################################################################
# This is the general function that run simulations of metacommunity diversity in explicit landscapes
# essentially use the "iteration.model.par" function defined above
#########################################################################################################

iteration.model<-function(replicas, Meta.pool, m.pool, Js, id.module,
                          M.dist, D50, m.max,
                          id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                          Lottery, it, prop.dead.by.it, id.obs){                   
  library(doParallel)                         ### Iterations in parallel
  registerDoParallel(cores=detectCores()-2)   #   Register cores 
  sub.it<-detectCores()-2                     #   Define iterations to run in parallel 
  tandas<-ceiling(replicas/sub.it)            #   define number of cycles to achive defined iteration (replicas)
  out.all<-NULL
  
  for(i in 1:tandas){                         # cycles of iteration runned in parallel
print(c("iteration  ", i, " of ", tandas))    
  out<-iteration.model.par(sub.it=sub.it, Meta.pool=Meta.pool, m.pool=m.pool, Js=Js, id.module=id.module,         # simulation of metacomm 
                                                        M.dist=M.dist, D50=D50, m.max=m.max,
                                                        id.fixed =id.fixed, D50.fixed=D50.fixed, 
                                                        m.max.fixed=m.max.fixed, comm.fixed=comm.fixed,
                                                        Lottery=Lottery, it=it, prop.dead.by.it=prop.dead.by.it, 
                                                        id.obs=id.obs)
    out.all<-rbind(out.all, out)
    }
  out.all
}
##############################################################

#########################################################################################################
# This function run parallel simulations of metacommunity diversity 
# essentially use the "H2020_Coalescent.and.lottery.exp.Kernel.J" function defined above
#########################################################################################################

iteration.model.par<-function(sub.it, Meta.pool, m.pool, Js, id.module,
                          M.dist, D50, m.max,
                          id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                          Lottery, it, prop.dead.by.it, id.obs){                   
#  library(doParallel)                            # if this function is used directly, the parallel lines have to be activated
#  registerDoParallel(cores=detectCores()-2)
#  sub.it<-detectCores()-2
#  tandas<-ceiling(replicas/sub.it)
    
    out<-foreach(i = 1:sub.it, .combine=rbind)%dopar% {
      
      H2020_Coalescent.and.lottery.exp.Kernel.J(Meta.pool=Meta.pool, m.pool=m.pool, Js=Js, id.module=id.module,
                                                M.dist=M.dist, D50=D50, m.max=m.max,
                                                id.fixed =id.fixed, D50.fixed=D50.fixed, 
                                                m.max.fixed=m.max.fixed, comm.fixed=comm.fixed,
                                                Lottery=Lottery, it=it, prop.dead.by.it=prop.dead.by.it, 
                                                id.obs=id.obs)
  }
  out
}



#########################################################################################################
# This function run one simulations of metacommunity diversity 
# It is a general function with several options 
# Meta.pool: Vector of species abundance in the "matacommunity" it is an external vector not affected by local dynamics
# m.pool: migration from the species pool. Is the probability that a reclutant is samples from the "metacommunity" vector
# Js: is the vector of total community size along all local communities (e.g. comunity size proportional to area)
# id.module: vector with local communities membership to spatial modules (externally estimated from landscape structure)
# M.dist: matrix with distances between local communities (could be Euclidian, topological or other)
# m.max: maximum migration between communities. Is migration when distance=0. Typically m.max=1 
#  D50: is the distance at which migration is half m.max
#  id.fixed: are communities that are not involved in the simulation (e.g. river outlet or connection with other metacommunities not simulated)
#  D50.fixed, m.max.fixed: fixed communities may have different dispersal Kernel
#  comm.fixed: vector of species abundances in the fixed community (e.g. Meta.pool)
#  Lottery=c(TRUE,FALSE) If FALSE only coalescent simulation is performed. If TRUE coalescent filling is followed by a lottery dynamic
#  it: number of iterations in the Lottery simulation
#  prop.dead.by.it: fraciton of the local community replaced in each iteration
#  id.obs: Communities for which diversity (alpha, beta) should be returned
#########################################################################################################

H2020_Coalescent.and.lottery.exp.Kernel.J<-function(Meta.pool, m.pool, Js,id.module,
                                                    M.dist, D50, m.max,
                                                    id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                                                    Lottery, it, prop.dead.by.it, id.obs){
  library(vegan)
  #if(dead.by.it>=J)return("dead.by.it cannot be larger than J")
  Meta.pool<-Meta.pool/sum(Meta.pool)               # standardize abundances
  comm.fixed<-comm.fixed/sum(comm.fixed)            # 
  Meta<-NULL                                        # Empty metacommunity
  M.migra<-H2020_migration.matrix.kernel.all(M.dist=M.dist, m.pool=m.pool, D50=D50, m.max=m.max,             # Funciton defined above. It estimates Migration matrix
                                             id.fixed=id.fixed, D50.fixed=D50.fixed, m.max.fixed=m.max.fixed)

  # Coalescent filling start here:
  for(i in 1:ncol(M.migra)){                        # first indiviual in each community randomly selected from the metacommunity
    Meta<-cbind(Meta, rmultinom(1,1,Meta.pool))
  }
  for (ii in 2:max(Js)){                            # local communities incorporate individuals until community size (Js)
    id.j<-which(Js>=ii)                             # identify communities to update (J.temporal<Js)
cat("coalescent construction in J: ", ii," de" ,max(Js),"\n")  # monitoring (automatically off in parallel)
if(length(id.fixed)>0) Meta[,id.fixed]<-comm.fixed*(ii-1)      # scale vector of abundances in the fixed community to the abundance of all other communities
    Pool.neighbor<-(Meta%*%M.migra)         # estimates potential reclutants including immigrants from all communities weighted by local abundances
    if(length(id.j)>1){new<-apply(Pool.neighbor[,id.j],2,born,dead.by.it = 1, M.pool = Meta.pool, m.pool = m.pool)   # random selection of new individuals from reclutants pool 
    Meta[,id.j]<-Meta[,id.j]+new} else {
      Meta[,id.j]<-Meta[,id.j]+born(n = Pool.neighbor[,id.j],dead.by.it = 1, M.pool = Meta.pool, m.pool = m.pool) 
    }                          # upadate communities 
  }
  

  if(Lottery==T){                                        # START LOTTERY ################################################
    dead.by.it<-ceiling(prop.dead.by.it*Js)    # estiamte individual to remove in each iteration and local community
    if(length(id.fixed)>0) dead.by.it[id.fixed]<-0                    # fixed communities are not updated
    max.dead.by.it<-max(dead.by.it)
    if(length(id.fixed)>0)Meta[,id.fixed]<-round(comm.fixed*max(Js),0)               # update abundances of the fixed community to community size
    for(iteration in 1:it){                              # start lottery iterations  
print(c(iteration, " of ", it))
      for(dead in 1:max.dead.by.it) {
        id.no.dead<-which(dead.by.it>=dead)               # identify communities to remove individuals
        Meta[,id.no.dead]<-Meta[,id.no.dead]-apply(Meta[,id.no.dead],2,FUN = change, change=1) # remove individuals along all communities
      }
      Pool.neighbor<-(Meta%*%M.migra)                    # estimates potential reclutant from all communities
      id.comm<-0
      for(reclutants in dead.by.it) {                    # dead.by.it is the vector of number of indiviuals to update in each iteration (a fixed fraction of J)
        id.comm<-id.comm+1
        Meta[,id.comm]<-Meta[,id.comm]+rmultinom(1,reclutants,
                                                 prob = (1-m.pool)*(Pool.neighbor[,id.comm]/sum(Pool.neighbor[,id.comm]))+m.pool*Meta.pool) # random selection of reclutants from neighbours, local or external pool
      }
    }
  }
  
  Meta.com.simulated<-apply(Meta[,id.obs],1,sum)                # vector of species abundances in the simulated metacommuty
  Meta.com.simulated.occ<-ifelse(Meta.com.simulated>0,1,0)      # vector of species that persited at the end of simulations
  Gamma<-sum(Meta.com.simulated.occ)                            # Total richness in the simulated metacommunity (Gamma diversity)
  BB<-as.matrix(vegdist(t(Meta[,id.obs]), method = "jaccard"))  # marix of beta diversity (Jaccard) between local communities
  Bett.all<-apply(BB,2,mean)                                    # average beta diversity 
  Bett.intra<-rep(NA,length(id.module))                         # vectors for beta diveristy between and withing modules
  Bett.inter<-rep(NA,length(id.module))                         #
  for (bb in unique(id.module)){                                # for each module...
    id.i<-which(id.module==bb)                                  # identify communities of this module
    b.intra<-apply(BB[id.i,id.i],2,mean)                        # estiamte beta diversity within
    Bett.intra[id.i]<-b.intra                                   # match beta with community
    b.inter<-apply(BB[-id.i,id.i],2,mean)                       # estiamte inter module beta
    Bett.inter[id.i]<-b.inter                                   # match beta inter with communitiess
  }
  
  Meta.out<-c("m.pool"=m.pool, "Js.max"=max(Js),"Js.min"=min(Js), "D50"=D50, "m.max"=m.max,
          "D50.fixed"=D50.fixed, "m.max.fixed"=m.max.fixed,
          "Lottery"=ifelse(Lottery==TRUE,1,0), "it"=it, "Gamma"=Gamma,
          "S.loc"=apply(ifelse(Meta[,id.obs]>0,1,0),2,sum),
          "B.loc.all"=Bett.all,
          "B.loc.intra.module"=Bett.intra,
          "B.loc.inter.module"=Bett.inter
          )
  Meta.out
}

############################################################################################################################
# Cells at contact in their edges are assumed the zero distance for migration
# m.max: migration between connected cells
# D50 is the distance at which migration decay yo half m.max

H2020_migration.matrix.kernel.all<-function(M.dist, m.pool, D50,m.max, 
                                            id.fixed, D50.fixed, m.max.fixed){
    diag(M.dist)=NA
    M.dist = M.dist-min(M.dist, na.rm=T) # min distance is the distance between neighbour cells. 
                              # connected cells has distance zero and migration m.max
    b = -log(0.5)/D50         # b is estimated | m(D50)=m.max*0.5
    M.migra = m.max*exp(-b*M.dist) 
    if(length(id.fixed)>0){
      b.fixed= -log(0.5)/D50.fixed 
      M.migra[id.fixed,]<- m.max.fixed*exp(-b.fixed*M.dist[id.fixed,]) # migration from outlet
    }                                                          # M.migra is the potential migration between communities, 
  diag(M.migra)<-1                                             # selfrecruitment is considered as 1=m.intra community
  M.migra<-apply(M.migra,2,m_to_1,m.pool)                      # standirize migrations to 1: (m.intra+m.pool+m.neigh=1)
  M.migra
}

############
m_to_1<-function(m, m.pool) (1-m.pool)*m/sum(m) # standarize a vertor of migration to add 1 
                                                # also considering migration from an external pool

born<-function(n, dead.by.it, M.pool, m.pool)rmultinom(1,dead.by.it,(1-m.pool)*(n/sum(n))+m.pool*M.pool)
change<-function(n,change)rmultinom(1,change,n)


######################################################################################################################
# function to resume output of simulations
resume.out<-function(out){
  out2<-list()
  out2[[1]]<-apply(out,2,quantile, 0.5, na.rm=T)    # NA originates if only a single module is involved
  out2[[2]]<-apply(out,2,sd, na.rm=T)
  out2[[3]]<-apply(out,2,quantile, 0.975, na.rm=T)
  out2[[4]]<-apply(out,2,quantile, 0.025, na.rm=T)
  names(out2)<-c("Median", "Standard Deviation", "out.IC.up","out.IC.inf")
  out2
  
}

#
#
#
#
