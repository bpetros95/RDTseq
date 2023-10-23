#MPHN2 outbreaker2 and distance matrix script, sep28, filtering genomes for only sites with reads in ALL
setwd("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/mphn_cluster_2/sep21_metagenomic_fastas")
library(outbreaker2)
library(ape)
library(visNetwork)
library(data.table)
library(ggplot2)
library(insect)
library(seqinr)
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
mphn2_dna = read.FASTA("mphn2_sep27_mafft_masked_pairwise_aligned.fasta")[-1] #fasta]
names(mphn2_dna) = c("M01","M12","M03","M11","M02","M04","M05","M06","M07","M08","M09","M10")

xbb_dna = mphn2_dna[-c(1,2,4,5)]
ef_dna = mphn2_dna[c(1,2,4,5)]

#mask known problematic sites and sites which have ambiguous reads in any genome
problematic_sites = data.table(read.table("~/Desktop/projects/covid_data_analysis/lhdseq/problematic_sites_sarsCov2.vcf", quote="\""))
problematic_sites_list = problematic_sites[V7=="mask"]$V2

#xbb
xbb_char_seqs = lapply(xbb_dna, function(x){s2c(dna2char(x))})
masked_xbb_char_seqs = lapply(xbb_char_seqs, function(x){return(x[-problematic_sites_list])})
any_genome_ambiguous_xbb = unique(unlist(lapply(masked_xbb_char_seqs,function(x){which(!(x%in%c("A","C","T","G")))})))
masked_xbb_char_seqs = lapply(xbb_char_seqs, function(x){return(x[-any_genome_ambiguous_xbb])})
masked_xbb_strings = lapply(masked_xbb_char_seqs,c2s)
masked_xbb_dna = char2dna(masked_xbb_strings)

#ef
ef_char_seqs = lapply(ef_dna, function(x){s2c(dna2char(x))})
masked_ef_char_seqs = lapply(ef_char_seqs, function(x){return(x[-problematic_sites_list])})
any_genome_ambiguous_ef = unique(unlist(lapply(masked_ef_char_seqs,function(x){which(!(x%in%c("A","C","T","G")))})))
masked_ef_char_seqs = lapply(ef_char_seqs, function(x){return(x[-any_genome_ambiguous_ef])})
masked_ef_strings = lapply(masked_ef_char_seqs,c2s)
masked_ef_dna = char2dna(masked_ef_strings)

#load epi metadata
#earliest of onset dates where available, or collect dates if no symptoms available
xbb_dates = as.Date(c("29-12-2022", #M03
                      "23-12-2022", #M04
                      "23-12-2022", #M05
                      "25-12-2022", #M06
                      "24-12-2022", #M07
                      "24-12-2022", #M08
                      "25-12-2022", #M09
                      "23-12-2022"), #M10
                    "%d-%m-%Y") 
names(xbb_dates) = names(masked_xbb_dna)
xbb_ctd = matrix(c("M07","M08", 
                   "M07","M04",
                   "M08","M04", 
                   "M08","M04",
                   "M06","M08",
                   "M03","M04",
                   "M03","M04"),7,2,byrow=TRUE)
masked_xbb_data = outbreaker_data(dates = xbb_dates, 
                                  dna = masked_xbb_dna, 
                                  w_dens = mphn2_w_dens,
                                  f_dens = mphn2_f_dens,
                                  ctd = xbb_ctd)  

ef_dates = as.Date(c("27-12-2022", #M01
                     "23-12-2022", #M12
                     "22-12-2022", #M11
                     "29-12-2022"), #M02
                   "%d-%m-%Y") 
names(ef_dates) = names(masked_ef_dna)
masked_ef_data = outbreaker_data(dates = ef_dates, 
                                 dna = masked_ef_dna, 
                                 w_dens = mphn2_w_dens,
                                 f_dens = mphn2_f_dens)  
mu = 2e-6
mphn2_config = create_config(move_kappa = FALSE, 
                             move_pi = TRUE,
                             prior_pi=c(1,4),
                             init_mu = mu,
                             move_mu = TRUE,
                             n_iter = n_iter,
                             find_import = TRUE, 
                             init_tree = "star")  

#run ob2
mphn2_outbreaker_xbb = outbreaker(data = masked_xbb_data, config = mphn2_config)
mphn2_outbreaker_ef = outbreaker(data = masked_ef_data, config = mphn2_config)

#get ob2 transmissions
ef_transmissions = get_transmissions(mphn2_outbreaker_ef)
ef_transmissions = data.table(ef_transmissions)
names(ef_transmissions) = c("from","to","value")
ef_names = names(masked_ef_dna)
ef_transmissions$from = as.character(ef_transmissions$from)
ef_transmissions[from!=0]$from = ef_names[as.numeric(ef_transmissions[from!=0]$from)]
ef_transmissions$to = as.character(ef_transmissions$to)
ef_transmissions[to!=0]$to = ef_names[as.numeric(ef_transmissions[to!=0]$to)]

xbb_transmissions = get_transmissions(mphn2_outbreaker_xbb)
xbb_transmissions = data.table(xbb_transmissions)
names(xbb_transmissions) = c("from","to","value")
xbb_names = names(masked_xbb_dna)
xbb_transmissions$from = as.character(xbb_transmissions$from)
xbb_transmissions[from!=0]$from = xbb_names[as.numeric(xbb_transmissions[from!=0]$from)]
xbb_transmissions$to = as.character(xbb_transmissions$to)
xbb_transmissions[to!=0]$to = xbb_names[as.numeric(xbb_transmissions[to!=0]$to)]

transmissions = rbind(ef_transmissions, xbb_transmissions)
transmissions = transmissions[order(value,decreasing=FALSE)]
write.csv(transmissions,"sep28_mphn2_multipleruns_obs_transmission_prob.csv")
pctl75 = dim(transmissions)[1]*0.75
pctl50 = dim(transmissions)[1]*0.5


#histogram of transmission likelihoods
hist(transmissions$value, breaks=35, main="",xlab="Probability")
abline(v=0.25,lwd=2,col="blue")
our_pctl = 1-dim(transmissions[value>=0.25])[1]/dim(transmissions)[1]
#add percentiles for comparison
#abline(v=transmissions[round(pctl75),]$value,col="red") #75th pctl
#abline(v=transmissions[round(pctl50),]$value,col="green") #50th pctl


prob_threshold = 0.25
#plot_transmissions = xbb_transmissions[value>=prob_threshold,]
plot_transmissions = rbind(ef_transmissions[value>=prob_threshold,], xbb_transmissions[value>=prob_threshold])
#visualise as graph
nodes = data.frame(id = c("0","M12","M11","M10","M09","M08","M07","M06","M05","M04","M03","M02","M01"),
                   label = c("0","M12","M11","M10","M09","M08","M07","M06","M05","M04","M03","M02","M01"),
                   color="blue")
visNetwork(nodes, plot_transmissions, width = "100%") %>% visEdges(arrows ="to")



#distance matrix per lineage
ef_dist_mat = dist.dna(masked_ef_dna,model="N",pairwise.deletion=TRUE,as.matrix=TRUE)
xbb_dist_mat = dist.dna(masked_xbb_dna,model="N",pairwise.deletion=TRUE,as.matrix=TRUE)

