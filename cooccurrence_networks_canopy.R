library(cooccur)
library(nlme)
library(vegan)

load('canopy.RData')

data<-net.canopy.dat
data[,4:40]<-1*(data[,4:40]>0)
invert<-colnames(data)[4:28]
algae<-colnames(data)[29:40]

canop<-unique(data$Canopy)
clearing<-unique(data$Clearing)
subdata<-data[,4:40]
subdata<-t(as.matrix(subdata[which(rowSums(subdata)>0),which(colSums(subdata)>0)]))

res<-cooccur(subdata,spp_names=TRUE)

w<-(res$results$obs_cooccur-res$results$exp_cooccur)/res$results$exp_cooccur
#plot(w,res$results$prob_cooccur)


sig_pairs<-res$results[which(res$results$prob_cooccur<0.05 & w>0),10:11]   #significant pairs 
dim(sig_pairs)
net<-c()   #network algae-invertebrates
for (pa in 1:dim(sig_pairs)[1]){
  if ((sig_pairs[pa,1] %in% invert) && (sig_pairs[pa,2] %in% algae)){
    net<-rbind(net,c(as.character(sig_pairs[pa,2]),as.character(sig_pairs[pa,1])))
  }
  if ((sig_pairs[pa,2] %in% invert) && (sig_pairs[pa,1] %in% algae)){
    net<-rbind(net,c(as.character(sig_pairs[pa,1]),as.character(sig_pairs[pa,2])))
  }
}

length(net)

date<-unique(as.character(data$Date))
length(date)

all_res<-list()
sc<-1
for (da in date){
  for (ca in canop){
    for (cl in clearing){
      subdata<-data[which((data$Canopy==ca)*(data$Clearing==cl)*(data$Date==da)==TRUE),4:40]
      spp_pres<-colnames(subdata)[which(colSums(subdata)>0)]
      subnet<-c()
      for (pa in 1:nrow(net)){
        if ((net[pa,1] %in% spp_pres) && (net[pa,2] %in% spp_pres)){
          subnet<-rbind(subnet,net[pa,1:2])
        }
      }
      all_res[[sc]]<-list('canopy'=ca,'clearings'=cl,'date'=da,'net'=subnet)
      sc<-sc+1
      print(c(ca,cl,da))}
  }
}


#every element of the list is now a co-occurrence network for each combination of canopy (4), clearing (2) and date (8)
#total= 64 networks
length(all_res)  
all_res[[1]]$canopy
all_res[[1]]$clearings
all_res[[1]]$date
all_res[[1]]$net

####some analyses on network structure

###null model
#option 1: randomize all
EE<-function(m){ 
  r<-dim(m)[1]
  c<-dim(m)[2]
  n_mat<-matrix(sample(m),r,c)
  while (min(c(rowSums(n_mat),colSums(n_mat)))==0){  
    n_mat<-matrix(sample(m),r,c)}  
  return(n_mat)}

#option 2: keep row totals constant (number of animals per alga) and shuffle row content
FE<-function(m){
  r<-dim(m)[1]
  n_mat<-c()
  for (j in 1:r){n_mat<-rbind(n_mat,sample(m[j,]))}  #shuffle row content
  while (min(c(rowSums(n_mat),colSums(n_mat)))==0){
    n_mat<-c()
    for (j in 1:r){n_mat<-cbind(n_mat,sample(m[j,]))}
  }
  return(n_mat)}

#option 3: keep row and column totals constant
FF<-function(m){   
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}

#all results
res_summary<-c()
for (i in 1:length(all_res)){
  el<-all_res[[i]]$net
  al<-unique(el[,1])
  an<-unique(el[,2])
  m<-matrix(rep(0,length(al)*length(an)),length(al))
  if (min(dim(m))>0){
    for (r in 1:nrow(el)){
      r_id<-which(al==el[r,1])
      c_id<-which(an==el[r,2])
      m[r_id,c_id]<-1
    }
    if (min(dim(m))>1){  
      raw_nodf<-as.numeric(nestednodf(m)$statistic[3])
      n_nodf<-c()
      for (n in 1:1000){
        n_nodf<-c(n_nodf,as.numeric(nestednodf(EE(m))$statistic[3]))
      }
      p_nodf<-sum(n_nodf>raw_nodf)/length(n_nodf)   #p-value 
      z_nodf<-(raw_nodf-mean(n_nodf))/sd(n_nodf)    #z-value
    }
    else {raw_nodf<-NA
    p_nodf<-NA   
    z_nodf<-NA
    }
    fill<-sum(m)/(dim(m)[1]*dim(m)[2])  
    occs<-sum(m)   
    nodes<-sum(dim(m))   
    res_summary<-rbind(res_summary,c(all_res[[i]]$canopy,all_res[[i]]$clearings,all_res[[i]]$date,raw_nodf,p_nodf,z_nodf,fill,occs,nodes))}
  print(i)
}


res_summary<-data.frame(
  'canopy'=as.factor(res_summary[,1]),
  'clearing'=as.factor(res_summary[,2]),
  'date'=as.factor(res_summary[,3]),
  'nodf'=as.numeric(res_summary[,4]), 
  'nodf_p'=as.numeric(res_summary[,4]), 
  'nodf_z'=as.numeric(res_summary[,6]),
  'fill'=as.numeric(res_summary[,7]),   #connectance (n. links/n. all possible links)
  'links'=as.numeric(res_summary[,8]),  #connectivity (n. links)
  'nodes'=as.numeric(res_summary[,9]))  #species richness (n. nodes)


lmeModel<-lme(nodf_z~canopy*clearing,random=~1|date,data=res_summary,na.action=na.omit)
summary(lmeModel)  
plot(ranef(lmeModel))
plot(lmeModel)

lmeModel2<-lme(nodf~canopy*clearing,random=~1|date,data=res_summary,na.action=na.omit)
summary(lmeModel2)   #no significant interaction

lme_fill<-lme(fill~canopy*clearing,random=~1|date,data=res_summary,na.action=na.omit)
summary(lme_fill)  #no significant interaction

lme_links<-lme(links~canopy*clearing,random=~1|date,data=res_summary,na.action=na.omit)
summary(lme_links)  #no significant interaction

lme_nodes<-lme(nodes~canopy*clearing,random=~1|date,data=res_summary,na.action=na.omit)
summary(lme_nodes)  #no significant interaction

pdf("nodf_vs_treatment",height=8,width=8)
par(mar=c(8,5,1,1))
boxplot(nodf_z~canopy*clearing,data=res_summary,cex.names=0.8,las=2,xlab=NULL)
dev.off()
