#outbreaker2 on mphn2
#updated aug10 
setwd("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/mphn_cluster_2/aug9_clean_fastas")
library(outbreaker2)
library(ape)
library(visNetwork)
library(data.table)
library(ggplot2)

#function to extract transmissions post-burnin from OB2 results
get_transmissions = function(x, burnin = burn, min_support = 0.1, labels = NULL, ...){
  if (burnin > max(x$step)) {
    stop("burnin exceeds the number of steps in x")
  }
  x <- x[x$step>burnin,,drop = FALSE]
  
  alpha <- as.matrix(x[,grep("alpha", names(x))])
  colnames(alpha) <- seq_len(ncol(alpha))
  from <- as.vector(alpha)
  to <- as.vector(col(alpha))
  from[is.na(from)] <- 0
  out_dat <- data.frame(xyTable(from,to))
  names(out_dat) <- c("from", "to", "frequency")
  ## Calculate proportion among ancestries
  get_prop <- function(i) {
    ind <- which(out_dat$to == out_dat$to[i])
    out_dat[[3]][i]/sum(out_dat[[3]][ind])
  }
  
  out_dat[,3] <- out_dat[,3] / nrow(x)
  
  out_dat
}

#function to rename transmission data with intuitive labels
rename = function(x){
  if(x==0){
    "MXX"
  }else{
    names(mphn2_dates)[x]
  }
}

#generation time
mean_gen = 2.9
sd_gen = 1.6
shape = (mean_gen^2)/(sd_gen^2)
rate = (mean_gen)/(sd_gen^2)

#incubation time
mean_inc = 3
sd_inc = 1.5
shape_inc = (mean_inc^2)/(sd_inc^2)
rate_inc = (mean_inc)/(sd_inc^2)

#collecting outbreaker2 data
n_iter = 100000 #number of iterations
burn = n_iter*.1 #discard first 10%
mphn2_w_dens = dgamma(1:21, shape=shape, rate=rate)
mphn2_f_dens = dgamma(1:21, shape=shape_inc, rate=rate_inc)
mphn2_dna = read.FASTA("swab_clean.fasta") #fasta
names(mphn2_dna) = c("M01","M03","M02","M04","M05","M06","M07","M08","M09","M10")

#load epi metadata
#earliest of onset dates where available, or collect dates if no symptoms available
mphn2_dates = as.Date(c("27-12-2022", #M01
                        "29-12-2022", #M03
                        "29-12-2022", #M02
                        "23-12-2022", #M04
                        "23-12-2022", #M05
                        "25-12-2022", #M06
                        "24-12-2022", #M07
                        "24-12-2022", #M08
                        "25-12-2022", #M09
                        "23-12-2022"), #M10
                      "%d-%m-%Y") 
names(mphn2_dates) = names(mphn2_dna)
mphn2_ctd = matrix(c("M07","M08", 
                     "M07","M04",
                     "M08","M04", 
                     "M08","M04", 
                     "M07","M01", 
                     "M06","M08",
                     "M03","M04",
                     "M03","M04"),3,2,byrow=TRUE)
mphn2_data = outbreaker_data(dates = mphn2_dates, 
                             dna = mphn2_dna, 
                             w_dens = mphn2_w_dens,
                             f_dens = mphn2_f_dens,
                             ctd = mphn2_ctd)  
mphn2_config = create_config(move_kappa = FALSE, # don't look for missing cases
                             move_pi = FALSE, # don't estimate reporting
                             init_pi = 1,
                             n_iter = n_iter,# set reporting to 1
                             find_import = FALSE, # don't look for additional imported cases
                             init_tree = "star")  # star-like tree as starting point

#run ob2
mphn2_outbreaker = outbreaker(data = mphn2_data, config = mphn2_config)

#get ob2 transmissions
transmissions = get_transmissions(mphn2_outbreaker)
transmissions$from = sapply(transmissions$from, rename)
transmissions$to = sapply(transmissions$to, rename)
names(transmissions) = c("from","to","value")
transmissions = data.table(transmissions)
prob_threshold = 0.25

#filter transmissions based on threshold
plot_transmissions = transmissions[value>=prob_threshold,]

#visualise as graph
nodes = data.frame(id = c("MXX","M10","M09","M08","M07","M06","M05","M04","M03","M02","M01"),
                   label = c("MXX","M10","M09","M08","M07","M06","M05","M04","M03","M02","M01"))
visNetwork(nodes, plot_transmissions, width = "100%") %>% visEdges(arrows ="to")

#table of transmissions and percentiles
transmissions = transmissions[order(value)]
pctl75 = dim(transmissions)[1]*0.75
pctl50 = dim(transmissions)[1]*0.5

#histogram of transmission likelihoods
hist(transmissions$value, breaks=35, main="",xlab="Probability")
abline(v=0.25)
our_pctl = 1-dim(transmissions[value>=0.25])[1]/dim(transmissions)[1]
#add percentiles for comparison
abline(v=0.19,col="red")
abline(v=0.27,col="green")

