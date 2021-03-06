###install.packages("ape")

#Example of one cycle
nodes = c(1,2,3,4,5) # make the list and call it nodes
nodes # print the list

nodecount = length(nodes) # save the number of nodes in the variable nodecount
tocoalesce = sample(1:nodecount, size=2) # sample 2 different nodes in node list
nodes[tocoalesce[1]] # print the first node sampled
nodes[tocoalesce[2]] # print the second node sampled

coalescencerate = nodecount*(nodecount-1)/2 # calculate the coalescent rate
coalescencetime = rexp(1, rate=coalescencerate) # sample from exponential w. that rate
coalescencetime

nodes <- nodes[-tocoalesce] # remove the two nodes that coalesced
nodes <- c(nodes,2*5-length(nodes)-1) # add the new node
nodes # print the new list

#load script
source("simulatecoalescencetrees.R")

#Run script
library(ape)
par (mfrow=c(2,5))
for (i in c(1:10)){
  print("New Tree")
  yourtree <-simtree(5) # simulate tree with 5 nodes
  ct<-read.tree(text=yourtree);plot(ct,cex=1.5);add.scale.bar(y=1.2,x=0.2,cex = 2,col = "red",lcol="red",lwd=3)
  print(" ")
}

#################################Afternoon
library(coalesceR) #package not avilable, I guess run on server
tree<-sim.tree('hudson',10,100000,100000, 1) # last argument: time change for pop size. BUILD A COALESCENT RELATION. TIME IS IN COALESCENT UNITS = 2N
draw.tree(tree)
mut<-add.mutations(tree, 20*1e-5) # change of mutation rate → we observe copies having groups of common mutations

# we have a coalescent tree + information about a mutation (in color + red dots along the tree for higher rates)
