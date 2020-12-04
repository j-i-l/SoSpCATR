source("network_generator.R")
source("get_meetings.R")
source("SocSpatR.R")

library(foreach)
library(igraph)

#Parameter sets for population and network generation
parsets=list(
list(groups=10,mean.group.size=10,max.group.size=15,d.eff=10,o.dens=10,i.dens=10,sex.eff=NA,m.i.eff=NA,m.o.eff=NA,plot=F),
list(groups=10,mean.group.size=10,max.group.size=15,d.eff=5,o.dens=5,i.dens=10,sex.eff=NA,m.i.eff=NA,m.o.eff=NA,plot=F),
list(groups=10,mean.group.size=10,max.group.size=15,d.eff=0.0001,o.dens=0.000000000001,i.dens=10,sex.eff=NA,m.i.eff=NA,m.o.eff=NA,plot=F)
)

#Parameter sets for simulating interactions
parsets2=c(T,F,F)
parsets3=c(0.8,0.4,0)

#generate a population
pop1=population.generator(groups=10,mean.group.size=10,max.group.size=15,F)

#for each of the 3 parameter sets, example of usage
foreach(i=1:length(parsets))%do%{
	currpar=c(parsets[[i]],population=list(pop1))
	#generate "true" network
	net1=do.call(network.generator,currpar)
	#generate interactions from "true" network
	edges=meetinglist(net1,nmeetings=5000,grouprestrict=parsets2[i],restrict=parsets3[i])

	net=graph_from_data_frame(edges$edgelist,directed=F)

	#get average position for each node
	for(i in V(net)){
		V(net)$x[V(net)==i]=mean(E(net)[.inc(i)]$s.x)
		V(net)$y[V(net)==i]=mean(E(net)[.inc(i)]$s.y)
	}

	#collapse edges
	E(net)$weight <- count_multiple(net)
	net <- simplify(net)

	#get social group membership
	socgroups=cluster_infomap(net)$membership

	#convert to matrix
	locs=as.matrix(data.frame(x=V(net)$x,y=V(net)$y))

	compare_soc_spat(locs,socgroups)
}











