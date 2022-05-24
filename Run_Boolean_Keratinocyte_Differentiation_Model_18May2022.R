
rm(list=ls(all=TRUE)) # Clean the workspace


setwd("C:/Users/Elisa/Dropbox/Epidermal_Differentiation/Individual switches/G_switch/Code G switch")

library(BoolNet) 

# Load network
net <- loadNetwork("Boolean_Model_for_Keratinocyte_Differentiation_19May2022.txt") 


# set basal layer conditions
PKC=0; IL1=0; E6=0; IL4=0; NFkB=0

# define the model with the basal layer conditions
net_Basal = fixGenes(net, c("PKC", "IL1",  "E6", "IL4", "NFkB"), c(PKC,IL1, E6, IL4, NFkB))

# get the attractors
attr_Basal<- getAttractors(net_Basal)
 # plot the attractors

plotAttractors(attr_Basal, mode="graph")

plotAttractors(attr_Basal, onColor="grey0", offColor = "white", drawLegend=FALSE)

Num_Attractors=length(attr_Basal$attractors)
Num_Attractors
Num_state_variables=9
position_FLG=9


basinSizebasal=c(1:Num_Attractors)
BasinSize_Flg=0
for (att_num in 1:Num_Attractors) {
  basinSizebasal[att_num]=attr_Basal$attractors[[att_num]]$basinSize/2^Num_state_variables
  
  if(sum(getAttractorSequence(attr_Basal,att_num)[,position_FLG])>0){BasinSize_Flg=BasinSize_Flg+ basinSizebasal[att_num]*(sum(getAttractorSequence(attr_Basal,att_num)[,position_FLG])/length(getAttractorSequence(attr_Basal,att_num)[,position_FLG]))}

}
BasinSize_Flg
sum(basinSizebasal)
basinSizebasal*100

tt_Basal=getTransitionTable(attr_Basal)


############## Calcium challenge experiment
# basal layer conditions
net_ca_healthy_inputs =fixGenes(net, c("PKC", "IL1",  "E6", "IL4", "NFkB"), c(1,0, 0, 0, 0))

attrCa<- getAttractors(net_ca_healthy_inputs)

plotAttractors(attrCa, mode="graph")
plotAttractors(attrCa, onColor="grey0", offColor = "white", drawLegend=FALSE)


BasinSizeCa=c(1:5)
BasinSizeCa_Flg=0

for(AttNumber in 1:5) {
  BasinSizeCa[AttNumber]=attrCa$attractors[[AttNumber]]$basinSize/2^9
  
  if(sum(getAttractorSequence(attrCa,AttNumber)[,9])>0){BasinSizeCa_Flg=BasinSizeCa_Flg+ BasinSizeCa[AttNumber]*(sum(getAttractorSequence(attrCa,AttNumber)[,9])/length(getAttractorSequence(attrCa,AttNumber)[,9]))}
  
  
}


BasinSize_Flg # without calcium

BasinSizeCa_Flg

barplot(c(BasinSize_Flg,  BasinSizeCa_Flg), main = "Basin size with Flg on",  col = rainbow(2), names.arg = c("Basal", "High Ca"))



BasinSizeCa

sum(BasinSizeCa)

BasinSizeCa=BasinSizeCa*100

#####
AA=getTransitionTable(attrCa)

#### Compute mean time to attractor
# for this we will use two commandos:
# AA$attractorAssignment: for each of the 2^n with n=9 states, it tells us to which attractor that state belongs to;
# ie it is a characterization of the basins of attraction.

# AA$initialState.FLG_AMP[[state]] tells us if the initial state is on or off for that element of the vector

counter=0
average_time=0

for(State in 1:2^9) {
  if (AA$attractorAssignment[[State]]==1 && AA$initialState.FLG_AMP[[State]]==0)
  {counter=counter+1
    average_time=average_time+AA$transitionsToAttractor[[State]]
  }
}

counter/2^9   # fraction of initial conditions in that state
avrage_time=average_time/counter

average_time

## Now we simulate a p53 mutation (HPV))

net_ca_E6 = fixGenes(net, c("PKC", "IL1",  "E6", "IL4", "NFkB"), c(1,0, 1, 0, 0))

attrCaE6<- getAttractors(net_ca_E6)

plotAttractors(attrCaE6, onColor="grey0", offColor = "white", drawLegend=FALSE)


BB=getTransitionTable(attrCaE6)

#####
#path <- getPathToAttractor(net, c(0,1,1,1,1))
#plotSequence(sequence=path)



#### Compute mean time to attractor
counterE6=0
average_timeE6=0

for(State in 1:2^9) {
  if (BB$initialState.FLG_AMP[[State]]==0)
  {counterE6=counterE6+1
  average_timeE6=average_timeE6+BB$transitionsToAttractor[[State]]
  }
}

avrage_timeE6=average_timeE6/counterE6

average_time #without HPV
average_timeE6


barplot(c(average_time,  average_timeE6), main = "Mean time to converge to a Flg on state upon addition of Ca",  col = rainbow(2), names.arg = c("High Ca", "High Ca+HPV"))

