rm(list = ls(all.names = TRUE))
#install.packages("Rlab")
library(Rlab)
library(bipartite)
install.packages("devtools")
devtools::install_github("kylebittinger/usedist")
library(usedist)
library("fields")

read.csv("fhWanangChisoTaxonomy for IBM.csv", header=T)->ctfsplot
read.csv("Ant plant list.csv")->ant.plant.list
#subset(ctfsplot,!is.na(ctfsplot$species))->ctfsplot
#subset(ctfsplot,!is.na(ctfsplot$genus))->ctfsplot
subset(ant.plant.list$species,ant.plant.list$ant.plant!="n")->ant.plants
subset(ant.plant.list$species,ant.plant.list$ant.plant=="n")->non.ant.plants
ant.plant<-vector(length=dim(ctfsplot)[1])
for(i in 1:dim(ctfsplot)[1]){
  ant.plant[i]<-ifelse(ctfsplot$Plant.sp[i]%in%ant.plants, "y", "n")
}
ctfs.ant.plants<-cbind(ctfsplot,ant.plant)

ctfs.subset<-subset(ctfs.ant.plants, ctfs.ant.plants$X<50)
ctfs.subset<-subset(ctfs.subset, ctfs.subset$Y<50)

comm <- ctfs.subset

#read.csv("toy data 2.csv")->comm

#########################NEXT STEPS
#make survival conditional on presence of other plants - need to include some negative density dependence from competition with other plants, but also positive density dependence from ant mutualists
#currently just have some neighbour impacts on survival, but nothing species specific
#at what stage to do this? just at the seedling stage? might be simplest... is there some biological basis for tree mortality being density independent beyond a certain age?
#make rules about dispersal and colonisation of new plants by ants - need to include fact that there is sharing of partners
#split plot up into subplots to prevent this from running too slowly at larger scales c.f. youtube video on physics sims
#include non ant-plant species
#

#read.csv("Wanang full plot (with made up ant data).csv")->comm

n.steps=500 # number of time steps to run

car.cap <- 50000

#limits of the plot - currently set up as a torus
x.limits<-c(0,50)
y.limits<-c(0,50)

# function for calculating distance matrix for all individuals
dist.function<-function(input){
  rdist(cbind(input$X,input$Y),compact=FALSE)->dists
}

# survival function to run on a single stem/mutualist combination. take a row as input. Need to make this dependent on nearby other plants and mutualists too.
survival.function<-function(foo,pop.size){ 
  #return(rbern(1,0.85))
  return(rbern(1,(car.cap-pop.size)/car.cap)) #currently just overall carrying capacity of the entire plot - need to change this to relate to density dependence for individual species
}

# new survival function - mortality of 0.9 if another tree within 20 m.
survival.function.2<-function(row.number,distance.matrix){ 
  pairwise.dists <- distance.matrix[,row.number]
  return(ifelse(length(pairwise.dists[pairwise.dists<0.5])>1, rbern(1,0.1), 1))
}

# plant reproduction function, currently uses negative binomial dispersal kernel, same across all species, need to include survival of seedlings dependent on proximity of other plants
reprod.function<-function(tree.row){   
  #for (i in 1:dim(foo)[1]){
  #comm[,1]->tree.row
    rbern(1,0.05)->reprod
    if(reprod==1){
      rexp(1, 1/26)->distance #negative exponential based on mean of median distances of three species from "Is Farther Seed Dispersal Better? Spatial Patterns of Offspring Mortality in Three Rainforest Tree Species with Different Dispersal Abilities"
      runif(1,0,2*pi)->direction
      x.add=distance*sin(direction)
      y.add=distance*cos(direction)
      #foo$X
      #foo$Y
      t(as.matrix(c(tree.row[1],sample(tree.row[2],1),runif(1,0.1,50),(as.numeric(tree.row[4])+x.add),(as.numeric(tree.row[5])+y.add),runif(1,0,1))))->new.t
    }
    if(reprod==0) {new.t<-NA}
    return(new.t)
  }

#preparing vectors for outputs
abundance<-vector(length=n.steps)
richness.lower<-vector(length=n.steps)
richness.upper<-vector(length=n.steps)

#for plotting
par(mfrow=c(2,2))

#run loop to calculate changes in community with time

for (j in 1:n.steps){
  #calculate distances between all trees
  dist.function(comm) -> distances
  
  #run the survival function and remove dead trees
  #pop.size<-dim(comm)[1]
  #surviving<-apply(comm, 1, survival.function, pop.size=pop.size)
  #comm<-comm[surviving==1,]
  
  #local density dependence using survival.function.2
  surviving <- vector(length=dim(comm)[1])
  for (i in 1:dim(comm)[1]){
    survival.function.2(i,distances) -> surviving[i]
  }
  comm<-comm[surviving==1,]
  
  #run the reproduction function and add new trees
  apply(comm, 1, reprod.function)->new.trees
  rm(to.add)
  matrix(unlist(new.trees[!is.na(new.trees)]), ncol = 6, nrow = length(new.trees[!is.na(new.trees)]),byrow=T)->to.add
  
  if(exists("to.add")){
    as.data.frame(to.add)->to.add
    names(comm)->names(to.add)
    rbind(comm,to.add)->comm
    comm <- as.data.frame(comm)
    comm[,c(3,4,5,6)] <- sapply(comm[,c(3,4,5,6)],as.numeric)
    
    #for any dispersal that went off the edge, make it come back on the other side (torus)
    comm$X[comm$X>x.limits[2]]<-comm$X[comm$X>x.limits[2]]%%x.limits[2]
    comm$X[comm$X<x.limits[1]]<-comm$X[comm$X<x.limits[1]]%%x.limits[2]
    comm$Y[comm$Y>y.limits[2]]<-comm$Y[comm$Y>y.limits[2]]%%y.limits[2]
    comm$Y[comm$Y<y.limits[1]]<-comm$Y[comm$Y<y.limits[1]]%%y.limits[2]
    
    #comm$X[comm$X>x.limits[2]]<-comm$X[comm$X>x.limits[2]]-x.limits[2]
    #comm$X[comm$X<x.limits[1]]<-comm$X[comm$X<x.limits[1]]+x.limits[2]
    #comm$Y[comm$Y>y.limits[2]]<-comm$Y[comm$Y>y.limits[2]]-y.limits[2]
    #comm$Y[comm$Y<y.limits[1]]<-comm$Y[comm$Y<y.limits[1]]+y.limits[2]
  }
  #save community summary stats in relevant vectors
  abundance[j]<-dim(comm)[1]
  richness.lower[j]<-length(unique(comm[,1]))
  richness.upper[j]<-length(unique(comm[,2]))
  
  #every ten timesteps, plot the communities
  if(j%%10==0) {
    plot(abundance)
    plot(richness.lower,ylim=c(0,500),xlim=c(0,n.steps),col="red")
    points(richness.upper,col="blue")
    plot(comm$Y~comm$X,ylim=c(0,50),xlim=c(0,50),col=as.factor(comm$Plant.sp))
    plotweb(as.matrix(xtabs(~Plant.sp +Ant.sp, comm, sparse = F)))
    Sys.sleep(0)
  }
  #if there are no individuals left, stop
  if(abundance[j]==0) break
}


#plot(comm$Y~comm$X,ylim=c(-500,500),xlim=c(-500,500),col=as.factor(comm$Plant.sp))
