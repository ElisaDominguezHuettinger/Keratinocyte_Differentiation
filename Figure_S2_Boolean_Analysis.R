rm(list=ls(all=TRUE)) # Clean the workspace

# Set working directory to source file location (should be the same location where the boolean model is stored)
setwd("C:/Users/Elisa/Dropbox/Epidermal_Differentiation/Individual switches/G_switch/Revision PLoS Comp/Final code")


# load the BoolNet library (if it is the first time, please install.package first)
library(BoolNet) 

# Load network
net <- loadNetwork("Boolean_Model_for_Keratinocyte_Differentiation.txt") 

# set basal layer conditions
PKC=0; IL1=0; E6=0; IL4=0; NFkB=0

# define the model with the basal layer conditions (note: we added some additional inputs here that are not discussed in the paper, please simply set to 0)
net_Basal = fixGenes(net, c("PKC", "IL1",  "E6", "IL4", "NFkB"), c(PKC,IL1, E6, IL4, NFkB))

# get the attractors
attr_Basal<- getAttractors(net_Basal)

# plot the attractors
plotAttractors(attr_Basal, onColor="grey0", offColor = "white", drawLegend=FALSE)

# count the number of attractors
Num_Attractors=length(attr_Basal$attractors)
Num_Attractors

# identify the attractor with the high expression of TDM state; it's in the 9th position of the state vector
Num_state_variables=9
position_FLG=9

# calculate the size of the basin of attraction of the high TDM state... by looping over all attractors and checking if TDM state is on, and if so, sum the basin size.
# if the attractor is a cyclic attractor with a fraction of states with the TDM=1, then divide the corresponding basin size 
# by the number of states  
basinSizebasal=c(1:Num_Attractors)
BasinSize_Flg=0
for (att_num in 1:Num_Attractors) {
  basinSizebasal[att_num]=attr_Basal$attractors[[att_num]]$basinSize/2^Num_state_variables
  
  if(sum(getAttractorSequence(attr_Basal,att_num)[,position_FLG])>0)
    {BasinSize_Flg=BasinSize_Flg+ basinSizebasal[att_num]*(sum(getAttractorSequence(attr_Basal,att_num)[,position_FLG])/length(getAttractorSequence(attr_Basal,att_num)[,position_FLG]))}

}

BasinSize_Flg # size of the basin of attraction corresponding to high flg (for cyclic: divided by two)
sum(basinSizebasal)
basinSizebasal*100




############## Calcium challenge experiment - we repeat the same sequence of steps as previously, but this time with PKC=1
# basal layer conditions with high Ca
net_ca_healthy_inputs =fixGenes(net, c("PKC", "IL1",  "E6", "IL4", "NFkB"), c(1,0, 0, 0, 0))

attrCa<- getAttractors(net_ca_healthy_inputs)

plotAttractors(attrCa, onColor="grey0", offColor = "white", drawLegend=FALSE)


BasinSizeCa=c(1:5) #there are 5 attractors in this case
BasinSizeCa_Flg=0

for(AttNumber in 1:5) {
  BasinSizeCa[AttNumber]=attrCa$attractors[[AttNumber]]$basinSize/2^9
  
  if(sum(getAttractorSequence(attrCa,AttNumber)[,9])>0){BasinSizeCa_Flg=BasinSizeCa_Flg+ BasinSizeCa[AttNumber]*(sum(getAttractorSequence(attrCa,AttNumber)[,9])/length(getAttractorSequence(attrCa,AttNumber)[,9]))}
  
  
}


BasinSize_Flg # without calcium

BasinSizeCa_Flg


# compare the sizes of the basins of attraction with TDM=1 with basal or high calcium conditions

barplot(c(BasinSize_Flg,  BasinSizeCa_Flg), main = "Basin size with Flg on",  col = rainbow(2), names.arg = c("Basal", "High Ca"))


## Now we confirm that the fixed point attractors are conserved under asynchronous update


# recall that the attractors with synchronous update regime are:

attr_Basal
attrCa


## first for the basal state
attr_Basal_ASYN<- getAttractors(net_Basal, type="asynchronous")
attr_Basal

plotAttractors(attr_Basal_ASYN, onColor="grey0", offColor = "white", drawLegend=FALSE)
# asynchronous has 4 fixed point attractors

plotAttractors(attr_Basal, onColor="grey0", offColor = "white", drawLegend=FALSE)
# syncronous has 4 fived point attactors and one two state cyclic attractor.

# We want to confirm that the fixed point attractors coincide between the synchronous and asynchronous update regime 
# for this we want to compare each of the 4 attractors of the asyncronous update regime with the 4 attractors of the syncronous, and
# find the attractor number (indexes) such that the atttractors coincide.
# 
MatchingBasal=c(NaN,NaN,NaN,NaN)
for (AttAsyn in 1:4){
  for(AttNumber in 1:4) {
    AA=getAttractorSequence(attr_Basal,AttNumber)
    BB=getAttractorSequence(attr_Basal_ASYN,AttAsyn)
    if( sum(AA-BB)==0){MatchingBasal[AttAsyn]=AttNumber}
  }
}
MatchingBasal


for (ii in 1:4){
  print(ii)
  print(getAttractorSequence(attr_Basal,MatchingBasal[ii]))
  print(getAttractorSequence(attr_Basal_ASYN,ii))
}

MatchingBasal[1]=3

for (ii in 1:4){
  print(ii)
  print(getAttractorSequence(attr_Basal,MatchingBasal[ii]))
  print(getAttractorSequence(attr_Basal_ASYN,ii))
}

## and now for the High calcium state

attrCa_ASYN<- getAttractors(net_ca_healthy_inputs, type="asynchronous")
attrCa_ASYN

plotAttractors(attrCa_ASYN, onColor="grey0", offColor = "white", drawLegend=FALSE)
plotAttractors(attrCa, onColor="grey0", offColor = "white", drawLegend=FALSE)

# 
Matching=c(NaN,NaN,NaN,NaN)
for (AttAsyn in 1:4){
for(AttNumber in 1:4) {
AA=getAttractorSequence(attrCa,AttNumber)
BB=getAttractorSequence(attrCa_ASYN,AttAsyn)
if( sum(AA-BB)==0){Matching[AttAsyn]=AttNumber}
}
}
Matching


for (ii in 1:4){
  print(ii)
  print(getAttractorSequence(attrCa,Matching[ii]))
  print(getAttractorSequence(attrCa_ASYN,ii))
}

