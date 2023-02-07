##################################################################################################################################
# This function run different NULL MODELs for landscape structure
##################################################################
# a gradient of dispersal abilities is considered (D50: distance at which dispersal is half its maximum)
##################################################################################################################################
# This function run one simulations of metacommunity diversity 
# It is a general function with several options 

################# 
#  n.random: number of random landscape to be generated
#  Meta.pool: Vector of species abundance in the "matacommunity" it is an external vector not affected by local dynamics
#  m.pool: migration from the species pool. Is the probability that a reclutant is samples from the "metacommunity" vector
#  Js: is the vector of total community size along all local communities (e.g. comunity size proportional to area)
#  id.module: vector with local communities membership to spatial modules (externally estimated from landscape structure)
#  D50.min, D50.max: Define the range of D50 values to be xplored (D50 is the distance at which dispersal decay to hal its maximum value)
#. THIS WAS SETTED IN LOG10 SCALE SINCE IS THE SCALE AT WHICH EFFECTS BECOME MORE EVIDENT
#  id.fixed: are communities that are not involved in the simulation (e.g. river outlet or connection with other metacommunities not simulated)
#  D50.fixed, m.max.fixed: fixed communities may have different dispersal Kernel
#  comm.fixed: vector of species abundances in the fixed community (e.g. Meta.pool)
#  Lottery=c(TRUE,FALSE) If FALSE only coalescent simulation is performed. If TRUE coalescent filling is followed by a lottery dynamic
#  it: number of iterations in the Lottery simulation
#  prop.dead.by.it: fraciton of the local community replaced in each iteration
#  id.obs: Communities for which diversity (alpha, beta) should be returned
#  M.X.Y is a two columns matrix with spatial location of communities
#### Effect of Random landscape along dispersal gradient
H2020_Random_Landscape_local_ponds.D50<-function(n.random,Meta.pool, m.pool,  Js, id.module,
                                             D50.min, D50.max, m.max, M.X.Y,
                                             id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                                             Lottery, it, prop.dead.by.it, id.obs){
d50<-10^seq(log10(D50.min), log10(D50.max),,30)
out<-NULL
N<-nrow(M.X.Y)
for(dd in d50){
print(paste("Evaluating disance:  ",dd)  )
# Start estimating diversity pattern in the real landscape  
  out.real.temp<-  iteration.model.par(sub.it=n.random, Meta.pool=Meta.pool, m.pool=m.pool, Js=Js, id.module=id.module, # Function defined in H2020_Lattice_exp_Kernel_J-environment
                                M.dist=as.matrix(dist(M.X.Y)), D50=dd, m.max=m.max,
                                id.fixed=id.fixed, D50.fixed=dd, m.max.fixed=m.max.fixed, comm.fixed=comm.fixed,
                                Lottery=Lottery, it=it, prop.dead.by.it=prop.dead.by.it, id.obs=id.obs)
  out.real.temp<- resume.out(out.real.temp) # Resume results for the real landscape (alpha, beta and gamma diversity)
  S.real<-out.real.temp[[1]][11:(10+N)]
  Be.all.real<-out.real.temp[[1]][(10+N+1):(10+N+N)]
  Gamma.real<-out.real.temp[[1]][10]
  
# Estiamte diversity in randomized landscapes
  out.null.temp<-H2020_Random_Landscape_local_ponds(n.random=n.random, Meta.pool=Meta.pool, m.pool=m.pool, 
                                     Js=Js, id.module=id.module,
                                     D50=dd, m.max=m.max, M.X.Y=M.X.Y,
                                     id.fixed=id.fixed, D50.fixed=D50.fixed, m.max.fixed=m.max.fixed, comm.fixed=comm.fixed,
                                     Lottery=Lottery, it=it, prop.dead.by.it=prop.dead.by.it, id.obs=id.obs)
  out.null.temp<-resume.out(out.null.temp)
  S.null<-out.null.temp[[1]][11:(10+N)]
  Be.all.null<-out.null.temp[[1]][(10+N+1):(10+N+N)]
  Gamma.null<-out.null.temp[[1]][10]
  sd.Gamma.null<-out.null.temp[[2]][10]
  out.t<-c(dd,mean(S.real),sd(S.real),mean(Be.all.real),sd(Be.all.real),Gamma.real,
         mean(S.null),sd(S.null),mean(Be.all.null),sd(Be.all.null), Gamma.null,sd.Gamma.null)
out<-rbind(out,out.t)
}
colnames(out)<-c("D50","mean(S.real)","sd(S.real)","mean(Be.all.real)","sd(Be.all.real)","Gamma.real",
                 "mean(S.null)","sd(S.null)","mean(Be.all.null)","sd(Be.all.null)", "Gamma.null", "sd(Gamma.null)")
out
}
######################
# this function run replicates of the function "iteration.model.par.random" defined below
H2020_Random_Landscape_local_ponds<-function(n.random,Meta.pool, m.pool, Js, id.module,
                                             D50, m.max,M.X.Y,
                                             id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                                             Lottery, it, prop.dead.by.it, id.obs){
  
  library(doParallel)
  registerDoParallel(cores=detectCores()-2)
  sub.it<-detectCores()-1
  tandas<-ceiling(n.random/sub.it)
  out.all<-NULL
  
  for(i in 1:tandas){
    print(c("iteration  ", i, " of ", tandas))    
    out<-iteration.model.par.random(sub.it=sub.it, Meta.pool=Meta.pool, m.pool=m.pool, Js=Js, id.module=id.module,
                             M.X.Y=M.X.Y, D50=D50, m.max=m.max,
                             id.fixed =id.fixed, D50.fixed=D50.fixed, 
                             m.max.fixed=m.max.fixed, comm.fixed=comm.fixed,
                             Lottery=Lottery, it=it, prop.dead.by.it=prop.dead.by.it, 
                             id.obs=id.obs)
    out.all<-rbind(out.all, out)
  }
  out.all
}
  
###########

# This function run in parallel the fucntion "H2020_Coalescent.and.lottery.exp.Kernel.J"
# In each simulation the Distance matrix is randomly generated
iteration.model.par.random<-function(sub.it, Meta.pool, m.pool, Js, id.module,
                              M.X.Y, D50, m.max,
                              id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                              Lottery, it, prop.dead.by.it, id.obs){                   
  #  library(doParallel)
  #  registerDoParallel(cores=detectCores()-2)
  #  sub.it<-detectCores()-2
  #  tandas<-ceiling(replicas/sub.it)
  
  out<-foreach(i = 1:sub.it, .combine=rbind)%dopar% {
    
    H2020_Coalescent.and.lottery.exp.Kernel.J(Meta.pool=Meta.pool, m.pool=m.pool, Js=Js, id.module,
                                              M.dist=as.matrix(dist(cbind(runif(n=nrow(M.X.Y),min = min(M.X.Y[,1]),max = max(M.X.Y[,1])), # M.dist is obtained from randomized locations
                                              runif(n=nrow(M.X.Y),min = min(M.X.Y[,2]),max = max(M.X.Y[,2]))))), 
                                              D50=D50, m.max=m.max,
                                              id.fixed =id.fixed, D50.fixed=D50.fixed, 
                                              m.max.fixed=m.max.fixed, comm.fixed=comm.fixed,
                                              Lottery=Lottery, it=it, prop.dead.by.it=prop.dead.by.it, 
                                              id.obs=id.obs)
  }
  out
}

