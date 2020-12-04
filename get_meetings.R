meetinglist=function(net1,intn=10,nmeetings=100,grouprestrict=T,restrict=1,plot=F){

require(foreach)
require(igraph)
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }
#size of group = distance from centroid that points can be at
#degree effects propensity to be alone
spat=net1$spat
network=net1$network
network2=graph_from_adjacency_matrix(network,mode="undirected",weighted=T)
inddata=net1$ind_data
inddata$deg=degree(network2)
spat$size=aggregate(indivs~groups,FUN=length,data=inddata)[,2]

scoords=lapply(1:nrow(spat),function(x){
	currx=spat$group.x[x]
	curry=spat$group.y[x]
	#newx=currx+rnorm2(intn,mean=0,spat$size[x])
	#newy=curry+rnorm2(intn,mean=0,spat$size[x])
	newx=currx+rnorm2(intn,mean=0,1)
	newy=curry+rnorm2(intn,mean=0,1)
	data.frame(group=x,s.id=paste(x,1:length(newx)),s.x=newx,s.y=newy)
	
})

allx=c(spat$group.x,do.call(rbind,scoords)$s.x)
ally=c(spat$group.y,do.call(rbind,scoords)$s.y)

if(plot){
	library(randomcoloR)
	n = max(spat$group.id)
	palette = distinctColorPalette(n)

	dev.new()
	plot(spat$group.x,spat$group.y,cex=spat$size,pch=12,col=palette,xlim=c(min(allx)-0.5,max(allx)+0.5),ylim=c(min(ally)-0.5,max(ally)+0.5))

	for(i in 1:length(scoords)){
		ccoords=scoords[[i]]
		points(ccoords$s.x,ccoords$s.y,pch=16,col=palette[i],cex=0.7)
	}

}



scoords2=do.call(rbind,scoords)

maxassoc=max(network)+1

edgelist=foreach(i = 1:nmeetings,.combine=rbind)%do%{
	#randomly choose a first individual
	index1=sample(which(inddata$deg>0),1)
	ind1=inddata[index1,]
	#choose individual to meet based on association
	ind2id=sample(colnames(network),1,prob=(network[row.names(network)==ind1$indivs,]/maxassoc)^2)
	(assoc=(network[row.names(network)==ind1$indivs,]/maxassoc)[colnames(network)==ind2id])
	ind2=inddata[inddata$indiv==ind2id,]
	(samegroup=ind1$groups==ind2$groups)
	if(samegroup&grouprestrict){
		meetplace=scoords2[sample(which(scoords2$group==ind1$groups),1),]
	}else if(ind1$x==ind2$x&ind1$y==ind2$y&grouprestrict){
		meetplace=scoords2[sample(which(scoords2$group==ind1$groups|scoords2$group==ind2$groups),1),]
	}else{
		group1=ind1$groups
		group2=ind2$groups
		x1=ind1$x
		x2=ind2$x
		y1=ind1$y
		y2=ind2$y
	
		halfx=((x1+x2)/2)
		halfy=((y1+y2)/2)

		dists3=as.matrix(dist(data.frame(x=c(halfx,scoords2$s.x),
			y=c(halfy,scoords2$s.y))))
		dists3=dists3[1,2:ncol(dists3)]
		names(dists3)=scoords2$s.id

		dists3=(dists3-min(dists3))/max(dists3-min(dists3))
		dists3=1.001-dists3				
		dists3[dists3<restrict]=0
		dists3=(dists3-min(dists3))/max(dists3-min(dists3))

		#plot((dists3)~dists3)
		dists3=dists3

		if(restrict==0){
			meetplace=scoords2[scoords2$s.id==sample(names(dists3),1),]
		}else{
			meetplace=scoords2[scoords2$s.id==sample(names(dists3),1,prob=dists3),]
		}
		
		if(plot){
			
			
			f=colorRamp(c("white", "blue"),alpha=T)
			#plot(scoords2$s.x[scoords2$group%in%c(group1,group2)],scoords2$s.y[scoords2$group%in%c(group1,group2)],col=rgb(f(dists3),maxColorValue=255),pch=16,cex=0.7)
			plot(scoords2$s.x,scoords2$s.y,col=rgb(f(dists3),maxColorValue=255),pch=16,cex=0.7)

			points(x1,y1,col="red",pch=15)
			points(x2,y2,col="blue",pch=15)
			points(meetplace$s.x,meetplace$s.y,col="purple",pch=16)
		}
		
	}
	data.frame(id1=ind1$indivs,id2=ind2$indivs,group1=ind1$groups,group2=ind2$groups,
		s.id=meetplace$s.id,s.x=meetplace$s.x,s.y=meetplace$s.y,assoc=assoc)
}
return(list(edgelist=edgelist,scoords=scoords2))
}

