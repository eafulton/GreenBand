## Network examples
rm(list = ls()) # clear memory
setwd("/Users/ful083/Work/Lenfest_EBFM/Network_&_Criticality_&_Diets/network_analysis/2022_test_curve_shape/")

library(data.table)
library(magrittr)
library(tidyverse)
library(tidygraph)
library(qgraph)
library(visNetwork)
library(intergraph)
library(network)
library(ggraph)
library("qgraph")
library(igraph)

library(ggplot2)
library(AICcmodavg)
library(plyr)
library(stringr)

# Load data
inDir <- "Input_Directory"

inName <- "Input_file_name"
infilecsv <- paste(inDir,inName,".csv",sep="")

use_weights <- 1
dietData <- read.table(infilecsv,header=TRUE, sep=",", row.names = 1)
filecumcsv <- paste(inName,"cumdistrib.csv",sep="")

IDname <- paste(inDir,"Group_IDs.csv",sep="")
localIDs <- read.table(IDname,header=TRUE, sep=",")

##################### CREATE NETWORK STRUCTURE ###########################
# Create edges df
# Calculate Degree in and Degree out
step1 <- reshape2::melt(as.matrix(dietData))
names(step1)[names(step1) == "Var1"] <- "preyName"
names(step1)[names(step1) == "Var2"] <- "predName"
names(step1)[names(step1) == "value"] <- "propDiet"
# Replace names with IDs
step2 <- step1 %>% dplyr::inner_join(localIDs,by=c("preyName" = "Group"))
names(step2)[names(step2) == "ID"] <- "prey_ID"
step2 <- step2 %>% dplyr::inner_join(localIDs,by=c("predName" = "Group"))
names(step2)[names(step2) == "ID"] <- "pred_ID"
nDiet <- step2[step2$propDiet != 0, ]

if (use_weights > 0) {
  nl <- subset(nDiet, select = -c(preyName, predName))
  names(nl)[names(nl) == "prey_ID"] <- "V1"
  names(nl)[names(nl) == "pred_ID"] <- "V2"
  nl$na <- "FALSE"
  nl$V1 <- as.numeric(nl$V1)
  nl$V2 <- as.numeric(nl$V2)
  nl$na <- as.logical.factor(nl$na)
  nl$na <- FALSE
  nl$weight <- nl$propDiet
  nl <- subset(nl, select = -c(propDiet))
} else {
  nl <- subset(nDiet, select = -c(preyName, predName, propDiet))
  names(nl)[names(nl) == "prey_ID"] <- "V1"
  names(nl)[names(nl) == "pred_ID"] <- "V2"
  nl$na <- "FALSE"
  nl$V1 <- as.numeric(nl$V1)
  nl$V2 <- as.numeric(nl$V2)
  nl$na <- as.logical.factor(nl$na)
  nl$na <- FALSE
}

nlv <- localIDs
nlv$na <- "FALSE"
nlv$na <- as.logical.factor(nlv$na)
nlv$na <- FALSE
nlv$vertex.names <- nlv$Group
names(nlv)[names(nlv) == "species_ID"] <- "intergraph_id"
nlv <- subset(nlv, select = -c(Group))

g <- asIgraph(nl, vertices=nlv, directed=TRUE)
set_vertex_attr(g, "label", value=nlv[,2])


##################### IGRAPH WORK UP ###########################
#g <- asIgraph(l$edges, vertices=l$vertexes, directed=TRUE)
#set_vertex_attr(g, "label", value=l$vertexes[,2])

##g <- graph.incidence(dietData, directed=TRUE, mode="in")  ## Don't use this as treating single spp when pred and prey as spearate entities
edges <- get.edgelist(g)   # Can also be extracted using E(g)
nodes <- V(g)$vertex.names    # Name of vertices
E(g)$weight  # Use this to see weight per edge (NULL here as incidence only)

# Adjacency matrix
g[]
# To see only the first line use g[1,]

## Plot the network
#plot(g)

## Use tidygraph instead - create the network (graph) - can also use as_tbl_graph() to convert an igraph or netwprk 
# library network to a tidygraph graph
food_tidy <- as_tbl_graph(g, directed = TRUE)

## Now do real network analysis
centRes <- centrality(food_tidy)
dnonorm <- degree(g, mode = "all", normalized = FALSE)

# Node strength (degree):
centRes$OutDegree # Or InDegree, it's the same in unweighted networks

# Closeness:
centRes$Closeness

# Betweenness:
centRes$Betweenness

# Plotting centrality measures
centralityPlot(food_tidy)
# Compare networks  with this centralityPlot by creating a list of graphs to compare - see Network Analysis in R Cookbook.pdf

#### Other statistical approaches to calculate the indicators
# Number of nodes
length(V(g))

# Find the standalone components - subwebs
clusters(g)   # Output shows its one interconnected cluster (community)
# If it wasn't all one community but you wanted to focus on one community you can suck out just that sub-web using
#subweb<-induced.subgraph(food_tidy, which(clusters(newnet)$membership == x)  where x is the mmerbship number of the subweb you want

# Average path length - can also be done using mean_distance(g)
average.path.length(g)

# Clustering coefficient
transitivity(g)

# Centralization degree score
centralization.degree(g)$centralization

#iGraph degree routine Degree centrality is simplest of the methods, it measures the number of connections between a node and all other nodes. 
d <- degree(g, mode = "in", normalized = TRUE)
write(d, file = "centralityIN.txt", ncolumns = 1, sep = " ")
d <- degree(g, mode = "out", normalized = TRUE)
write(d, file = "centralityOUT.txt", ncolumns = 1, sep = " ")
d <- degree(g, mode = "all", normalized = TRUE)
write(d, file = "centralityALL.txt", ncolumns = 1, sep = " ")
dnonorm <- degree(g, mode = "all", normalized = FALSE)
write(dnonorm, file = "degree.txt", ncolumns = 1, sep = " ")

# Closeness centrality is an evaluation of the proximity of a node to all other
# nodes in a network, not only the nodes to which it is directly connected
c <- closeness(g, mode="in", weights=NA, normalized=T)
write(c, file = "closenessIN.txt", ncolumns = 1, sep = " ")
c <- closeness(g, mode="out", weights=NA, normalized=T)
write(c, file = "closenessOUT.txt", ncolumns = 1, sep = " ")
c <- closeness(g, mode="all", weights=NA, normalized=T)
write(c, file = "closenessALL.txt", ncolumns = 1, sep = " ")

# Betweenness centrality looks for chokepoints in the network
b <- betweenness(g, directed=F, weights=NA, normalized = T)
write(b, file = "betweeness.txt", ncolumns = 1, sep = " ")

# Distance between nodes
bs <- distances(g, v=V(g), to=V(g), weights=NA)

# Distance to specific node (given by node ID number)
distances(g, v=V(g), to=71, weights=NA)

# Look at degree distribution - important for some network structures 
g.degrees <- degree(g)

# Let's count the frequencies of each degree
g.degree.histogram <- as.data.frame(table(g.degrees))

# Need to convert the first column to numbers, otherwise
# the log-log thing will not work
g.degree.histogram[,1] <- as.numeric(g.degree.histogram[,1])

# Now, plot it
ggplot(g.degree.histogram, aes(x = g.degrees, y = Freq)) +
  geom_point() +
  scale_x_continuous("Degree\n(nodes with this amount of connections)",
                     breaks = c(1, 3, 10, 30, 100, 300),
                     trans = "log10") +
  scale_y_continuous("Frequency\n(how many of them)",
                     breaks = c(1, 3, 10, 30, 100, 300, 1000),
                     trans = "log10") +
  ggtitle("Degree Distribution (log-log)") +
  theme_bw()

# Or as a histogram
#hist(degree(g))
ggplot(g.degree.histogram, aes(x = g.degrees)) +
  geom_histogram() +
  scale_x_continuous("Degree\n(nodes with this amount of connections)",
                     breaks = c(1, 3, 10, 30, 100, 300),
                     trans = "log10") +
  scale_y_continuous("Frequency\n(how many of them)",
                     breaks = c(1, 3, 10, 30, 100, 300, 1000),
                     trans = "log10") +
  ggtitle("Degree Distribution (log-log)") +
  theme_bw()
ggplot(g.degree.histogram, aes(x = g.degrees)) +
  geom_histogram(binwidth=1) +
  theme_bw()

# Old method
# Non-cumulative
dd <- degree.distribution(g, mode = "all", cumulative = FALSE)
plot(dd)

# Now cumualtive
dd <- degree.distribution(g,mode = "all", cumulative = TRUE)
# Alternatively plot(degree.distribution(d,mode="out"))
degree_x <- 1:length(dd)

# Fit a power law”
probability <- dd[-1]
nonzero.position <- which(probability != 0)
probability <- probability[nonzero.position]
degree <- 1:max(dnonorm)
degree <- degree[nonzero.position]

# Now  see if there is a power law fit to the distribution
reg <- lm(log(probability) ~ log(degree))
summary(reg)$r.square
cozf <- coef(reg)
plot(probability ~degree)
plot(probability ~ degree, log = "xy")
power.law.fit <- function(x) exp(cozf[[1]]+cozf[[2]]*log(x))
curve(power.law.fit, col = "red", add = T)

# Nicer plot
cumdistrib.df <- data.frame(degree_x, dd)
ggplot(cumdistrib.df, aes(x = degree_x, y = dd)) +
  geom_point() +
  scale_x_log10("Degree") +
  scale_y_log10("Cumulative distribution") +
  #stat_smooth(method="lm")
  #stat_smooth(method = "lm", formula = y ~ poly(x, 5, raw = TRUE), color = 'red')
  theme_bw()

write.csv(cumdistrib.df, file = filecumcsv, row.names = FALSE)

############################ Fitting Power Law Curves ############################

dat <- read.table(filecumcsv,header=TRUE, sep=",")
names(dat)[names(dat) == "degree_x"] <- "x"
names(dat)[names(dat) == "dd"] <- "y"

# Power curve	
# y ~ x ^ (-a)

# Scaled power curve	
# y ~ a * x ^ (-b)

# Power curve with tail	
# y ~ x ^ (-a) * exp(-x / b)

# Scaled power curve with tail	
# y ~ a * x ^ (-b) * exp(-x / c)

# Linear
# y ~ a * x + b

models <- list(lm(y ~ x, data = dat), 
               nls(y ~ (a + b * log(x)), data = dat, start = setNames(coef(lm(y ~ log(x), data = dat)), c("a", "b"))),
               nls(y ~ (a + x ^ b), data = dat, start = list(a = 0, b = -1)), 
               nls(y ~ (a * x ^ b), data = dat, start = list(a = 1, b = -1)),
               nls(y ~ ((x ^ a) * exp(x / b)), data = dat, start = list(a = -1, b = -1)),
               nls(y ~ (a + b * x), data = dat, start = list(a = 1, b = -1))
)



# have a quick look at the visual fit of these models
ggplot(dat, aes(x, y)) + geom_point(size = 5) +
  stat_smooth(method = lm, formula = as.formula(models[[1]]), size = 1, se = FALSE, color = "black") + 
  stat_smooth(method = nls, formula = as.formula(models[[2]]), data = dat, method.args = list(start = setNames(coef(lm(y ~ log(x), data = dat)), c("a", "b"))), size = 1, se = FALSE, color = "green", linetype = 2) +
  stat_smooth(method = nls, formula = as.formula(models[[3]]), data = dat, method.args = list(start = list(a = 0,b = 0)), size = 1, se = FALSE, color = "red", linetype = 2) + 
  stat_smooth(method = nls, formula = as.formula(models[[4]]), data = dat, method.args = list(start = list(a = 0,b = 0)), size = 1, se = FALSE, color = "violet") +
  stat_smooth(method = nls, formula = as.formula(models[[5]]), data = dat, method.args = list(start = list(a = 0,b = 0)), size = 1, se = FALSE, color = "black") +
  stat_smooth(method = nls, formula = as.formula(models[[6]]), data = dat, method.args = list(start = list(a = 0,b = 0)), size = 1, se = FALSE, color = "orange", linetype = 2)

# AIC
ldply(models, function(mod){ data.frame(AICc = AICc(mod), AIC = AIC(mod), model = deparse(formula(mod))) })
models[[5]]
models[[6]]
models[[3]]


# Example output
AICc         AIC                  model
1 -123.242736 -123.874315                  y ~ x
2  -46.786310  -47.417889   y ~ (a + b * log(x))
3   -1.408865   -2.040443          y ~ (a + x^b)
4  -10.569838  -11.201417          y ~ (a * x^b)
5 -129.241956 -129.873535 y ~ ((x^a) * exp(x/b))
6 -123.242736 -123.874315        y ~ (a + b * x)

Nonlinear regression model
model: y ~ ((x^a) * exp(x/b))
data: dat
a        b 
0.2442 -12.8110 
residual sum-of-squares: 0.09678

Number of iterations to convergence: 16 
Achieved convergence tolerance: 9.454e-06

####### Plotting the different curves #######################
# If multiple scenarios are run, you can either merge to a single datarame dat or load them from a compiled file
# Format of teh file is
# Degree  Scenario1  Scenario2  Scenario3
# where degeree increments from 1 to the maximum number of the biggest network
# the values in the Scenario columns are the cumulative_scores per degree for that Scenario.

dat <- read.table("cum_distribs.csv",header=TRUE, sep=",", row.names = 1)

df <- reshape2::melt(as.matrix(dat))
names(df)[names(df) == "Var1"] <- "Degree"
names(df)[names(df) == "Var2"] <- "Scenario"
names(df)[names(df) == "value"] <- "CumScore"
ggplot(df, aes(x=Degree, y=CumScore)) +
  geom_point(aes(colour = factor(Scenario))) +
  geom_line(aes(colour = factor(Scenario))) +
  #scale_color_discrete(name="Scenario", breaks=c("F0", "F_P", "FMSY", "FHalfM"),labels=c("Unfished", "Production scaled", "MMSY", "0.5 * M")) +
  scale_colour_viridis_d("Scenario") +
  theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=20,face="bold"), legend.title=element_text(size=18,face="bold"), legend.text = element_text(size=16))
  
