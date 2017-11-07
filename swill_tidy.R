# Code associated with the manuscript:
# "Support amongst UK pig farmers and agricultural stakeholders for the use of food losses in animal feed"


# -------------------------- CONTENTS ------------------------------
#
# i. libraries
# 1. Load the survey data and tidy up the data.frame
# 2. Plot S2 Appendix Fig 1.
# 3. Plot S2 Appendix Fig 2.
# 4. Other information about study respondents and survey
# 5. Plot S2 Appendix Fig 3.
# 6. Models of acceptability of different food losses as feed
# ... plot Fig 1
# ... plot S2 Appendix Fig 4
# 7. Plot Fig 2, swill vs conventional feed 
# 8. Plot S2 Appendix Fig 5
# 9. Plot S2 Appendix Fig 6
# 10. Plot Fig 3, support for relegalisation
# 11. Plot S2 Appendix Fig 7, traditional/unnatural prractice
# 12. Factor analysis of respondent values and perceptions of swill
# ... S2 Appendix Figs 8 & 9
# 13. Model of support for legal change
# ... plot Figs 4, 5, & 6
# ... Plot S2 Appendix Figs 10, 11, & 12


# ----------------- i. libraries --------------------


# Load packages used in this analysis.


# data manipulation
library(plyr)
library(dplyr)
library(tidyr)
# plotting
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(RColorBrewer)
# Bayesian modelling
library(rethinking)


# ----- multiplot function ----- #
# This can be used to put multiple ggplot graphs on the same window
# grid.arrange also works well
multiplot <- function(..., plotlist=NULL, file, 
                      cols=2,               # edit number of cols
                      layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# ----- function for model averaging of cumulative link models ----- #
# Create a custom ensemble function for model averaging, because the ensemble function in "rethinking" package
# ... does not handle ordered logistic models (there is an error in the definition of 
# ... link_out, because link.list is an array, not a list for these models)
# ... This arises because we predict the cumulative probability of each Likert score (1-5),
# ... rather than a single dependent variable.
ensemble_custom <- function (..., data, n = 1000, func = WAIC, weights, refresh = 0, 
                             replace = list(), do_link = TRUE, do_sim = TRUE) {
  L <- list(...)
  # check there is more than one model
  if (is.list(L[[1]]) && length(L) == 1) 
    L <- L[[1]]
  mnames <- match.call()
  mnames <- as.character(mnames)[2:(length(L) + 1)]
  # check if weights are explicitly stated in function call
  if (missing(weights)) {
    # if more than one model, calculate the relative weights of each model (sums to 1)
    if (length(L) > 1) {
      use_func <- func
      ictab <- compare(..., func = use_func, refresh = refresh, 
                       n = n, sort = FALSE)
      # ictab <- compare(m6,m7, func = use_func, refresh = refresh, 
      #                  n = n, sort = FALSE)
      rownames(ictab@output) <- mnames
      weights <- ictab@output$weight
    }
    # there is only one model, the weight == 1
    else {
      ictab <- NA
      weights <- 1
    }
  }
  # if weights are NOT explicitly stated in function call 
  # ... then use current weights (with their sum set to 1)
  else {
    ictab <- NA
    weights <- weights/sum(weights)
  }
  # idx is the weights multiplied by number of samples
  idx <- round(weights * n)
  # make lists for link() and sim() outputs
  link.list <- list("empty")
  sim.list <- list("empty")
  # if the data are not specified,  calculate link.list and sim.list
  if (missing(data)) {
    if (do_link == TRUE) 
      link.list <- lapply(L, function(m) link(m, n = n, 
                                              refresh = refresh, replace = replace))
    if (do_sim == TRUE) 
      sim.list <- lapply(L, function(m) sim(m, n = n, refresh = refresh, 
                                            replace = replace))
  }
  # if the data ARE specified,  calculate link.list and sim.list
  else {
    if (do_link == TRUE) 
      link.list <- lapply(L, function(m) link(m, data = data, 
                                              n = n, refresh = refresh, replace = replace))
    if (do_sim == TRUE) 
      sim.list <- lapply(L, function(m) sim(m, data = data, 
                                            n = n, refresh = refresh, replace = replace))
  }
  # start = 1 for all models,
  # idx_end = model weight * 1000
  idx_start <- rep(1, length(idx))
  idx_end <- idx
  # For collections of more than one model:
  # ... not sure
  if (length(L) > 1) 
    # for 2:num(models)
    for (i in 2:length(idx)) {
      
      # idx_start is the minimum of (the idx_start of the previous model + the weighting*1000 of the previous model) or the number of iterations (i.e. the max weighting)
      idx_start[i] <- min(idx_start[i - 1] + 
                            idx[i - 1], 
                          n)
      # idx_start is the minimum of (the idx_start of the current model, + weight*1000 of current model) or the 
      idx_end[i] <- min(idx_start[i] + idx[i] - 1, n)
      if (i == length(idx)) 
        idx_end[i] <- n
    }
  link_out <- link.list[[1]]
  sim_out <- sim.list[[1]]
  for (i in 1:length(idx)) {
    if (idx[i] > 0) {
      idxrange <- idx_start[i]:idx_end[i]
      if (do_link == TRUE) 
        # Nb !!! here I add a second comma, to the ensemble function
        # ... since link_out has three dimensions for ordered categorical modelling 
        link_out[idxrange, ,] <- link.list[[i]][idxrange, ,
                                                ]
      if (do_sim == TRUE) 
        sim_out[idxrange, ] <- sim.list[[i]][idxrange, 
                                             ]
    }
  }
  result <- list(link = link_out, sim = sim_out)
  names(weights) <- mnames
  idxtab <- cbind(idx_start, idx_end)
  rownames(idxtab) <- mnames
  attr(result, "weights") <- weights
  attr(result, "indices") <- idxtab
  attr(result, "ictab") <- ictab
  result
}


# -------------------------- 1. Load the survey data and tidy up the data.frame -----------------------------


# Load the data
# ... Nb you will need to save the data in the [R] working directory (or tell [R] where to find it) 
survey_data <- read.csv("swill_survey_data.csv", header=T)   # Load data
head(survey_data)                                            # Vizualise the data.frame
dim(survey_data)                                             # 165 rows, 88 cols


# Make mini df with the questions ONLY
# ... i.e. the first row
qu.df <- survey_data[1,]


# And make a df without the questions
# ... and the second row, the "ImportId".
no.qu.df <- survey_data[-c(1,2),]
dim(survey_data); dim(no.qu.df) # two rows should be removed


# Add a column to group respondents into pig farmer/other agricultural
pig_farmers <- c("Farmer/farm manager of both pig and poultry", 
                 "Pig farmer/pig farm manager")
no.qu.df$job_group <- ifelse(no.qu.df$Q57 %in% pig_farmers, "Pig farmer", "Other")
unique(no.qu.df$job_group)
rm(pig_farmers) # tidy up


# tidy up
rm(survey_data)


# ------------------------------- 2. Plot S2 Appendix Fig 1. ----------------------------------


# Split by job
# ... Look at number of respondents per job
table(no.qu.df$Q57) # Look at number of respondents per job
no.qu.df$Q57 <- factor(no.qu.df$Q57, 
                       levels=rev(c("Pig farmer/pig farm manager", "Farmer/farm manager of both pig and poultry",
                                    "Agricultural advisor", "Poultry farmer/poultry farm manager","Trader", 
                                "Veterinarian", "Feed processor", "Student", "Retailer", "Food service industry",
                                "Other: involved in the animal industry", "Other: NOT involved in the animal industry")),
                       labels=rev(c("Pig farmer/pig farm manager", "Farmer/farm manager of\nboth pigs and poultry",
                                    "Agricultural advisor", "Poultry farmer or \npoultry farm manager","Trader", 
                                "Veterinarian", "Feed processor", "Student", "Retailer", "Food service industry",
                                "Other: involved in\nthe animal industry", "Other: not involved in\nthe animal industry")))
# no.qu.df$job_group <- factor(no.qu.df$job_group, 
#                              levels=c("Pig farmer", "Other"))
Fig1a <- ggplot(no.qu.df, aes(Q57, fill = job_group)) + geom_bar(stat="count") + coord_flip() + 
  xlab("Profession of respondents") + 
  ylab("Number of respondents") + 
  theme_bw() + 
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=12), 
        legend.text = element_text(size=12), 
        legend.position = "none") +
  scale_fill_manual(values = rev(c("#0072B2","#D55E00")) )
Fig1a


# ... what do Other, non-animal industry do?
# ... NB these six include one former pig farmer and two respondents not giving any details.
no.qu.df[no.qu.df$Q57 =="Other: NOT involved in the animal industry",c("Q57", "Q57_12_TEXT")]
nrow(no.qu.df[no.qu.df$Q57 =="Other: NOT involved in the animal industry",])


# What is the split in ages and gender?
# ... no difference in age profile between pig farmers/other
table(no.qu.df$Q55)        # 142 males, 21 females
summary(no.qu.df$Q56)      # median age bracket is 31-50
Fig1b <- ggplot(no.qu.df, aes(Q56, fill=job_group)) + geom_bar(stat="count", position="dodge") + 
  xlab("Respondent's age") + 
  ylab("Number of respondents") + 
  scale_fill_manual(values = rev(c("#0072B2","#D55E00")),
                    breaks=c("Pig farmer", "Other")) + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.key.size = unit(2, 'lines'),
        legend.position = c(0.2,0.85))
Fig1b


# Produce Fig 1.
multiplot(Fig1a, Fig1b)


# tidy up
rm(Fig1a, Fig1b)


# ----------------------------------- 3. Plot S2 Appendix Fig 2. -----------------------------------


# Load Eurostat data for pig farm size distribution in UK in 2013
eurostat <- read.csv("data/eurostat_pig_farm_size.csv", header=T)
eurostat2 <- ddply(eurostat, .(groups_to_match), summarize, num_farms = sum(num_holdings)) # Group together farm size classes to match survey data categories
eurostat2



# Looking only at pig farmers, what farm sizes are there?
pig.farmer.df <- no.qu.df[no.qu.df$job_group == "Pig farmer",] # make a separate df
nrow(pig.farmer.df) # 82 pig farmers


# What farm sizes in the survey?
# ... Mostly large farms >1000 pigs
# .. Nb seven farm size records are entered incorrectly as dates, rather than ranges
# ... i.e. 1-9 = 9-Jan, and 10-99 = Oct-99. I correct this below.
# ... I also pool the "1000-4999", "5000+" groups into a "1000+" group to match the size classes in the Eurostat data
table(pig.farmer.df$Q58) # There are seven responses erroneously recorded as dates (?)
pig.farmer.df[pig.farmer.df$Q58 %in% c("9-Jan", "Oct-99"),]               # for the erroneous farm size data, the rest of their data looks fine
levels(pig.farmer.df$Q58) <- c(levels(pig.farmer.df$Q58), "1-9", "10-99") # Add new factor levels
pig.farmer.df$Q58[pig.farmer.df$Q58 =="9-Jan"] <- "1-9"
pig.farmer.df$Q58[pig.farmer.df$Q58 =="Oct-99"] <- "10-99"
levels(pig.farmer.df$Q58)
pig.farmer.df$Q58 <- factor(pig.farmer.df$Q58, levels = c("1-9", "10-99","100-199",
                                                          "200-399", "400-999",
                                                          "1000-4999", "5000+"))
pig.farmer.df$groups_to_match <- ifelse(pig.farmer.df$Q58 == "5000+" | pig.farmer.df$Q58 == "1000-4999", "1000+", as.character(pig.farmer.df$Q58)) # Make groups_to_match column, equal to Q58 values, but grouping together all farms >1000 pigs
unique(pig.farmer.df$groups_to_match)              # check that the "1000-4999" & "5000+" values don't appear
head(pig.farmer.df[, c("Q58","groups_to_match")])  # Check the data look OK.


# Sum the number of rows per farm size
survey_farm_sizes <- pig.farmer.df %>% 
  group_by(groups_to_match) %>% 
  count(.)
# survey_farm_sizes <- as.data.frame(survey_farm_sizes)
survey_farm_sizes         # Look at the data
sum(survey_farm_sizes$n)  # the number of responses should sum to 82


# Standardize tables (this isn't really necessary, unless I merge the two df)
colnames(survey_farm_sizes)[2] <- "num_farms" # standardize column name


# individual plots
survey_farm_sizes$groups_to_match <- factor(survey_farm_sizes$groups_to_match, # Reorder farm sizes for plotting
                                            levels= c("1-9", "10-99","100-199",
                                                      "200-399", "400-999",
                                                      "1000+"))
pSurvey <- ggplot(survey_farm_sizes, aes(groups_to_match, num_farms)) + geom_bar(stat="identity", fill = "#0072B2") + 
  theme_bw() + 
  xlab("Number of pigs per farm") + 
  ylab("Number of survey respondents") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14)) + 
  annotate("text", label = "a", x = .75, y = 60, size = 8,fontface =2)
eurostat2$groups_to_match <- factor(survey_farm_sizes$groups_to_match,         # Reorder farm sizes for plotting
                                    levels= c("1-9", "10-99","100-199",
                                              "200-399", "400-999",
                                              "1000+"))
pEuro <- ggplot(eurostat2, aes(groups_to_match, num_farms)) + geom_bar(stat="identity", fill = "#009E73") + 
  theme_bw() + 
  xlab("Number of pigs per farm") + 
  ylab("Number of pig farms in UK") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14)) + 
  annotate("text", label = "b", x = 0.75, y = 5000, size = 8,fontface =2) + 
  ylim(0, 5000)


# Plot comparison
multiplot(pSurvey, pEuro)


# What proportion of pigs are reared on farms of different sizes in the UK?
sum(eurostat[eurostat$groups_to_match %in% c("1-9","10-99"),"num_heads"]) / sum(eurostat$num_heads)
eurostat[eurostat$groups_to_match == "1000+","num_heads"] / sum(eurostat$num_heads) # 84% of pigs reared on farms with > 1000 pigs


# tidy up
rm(pEuro, pSurvey, survey_farm_sizes, 
   eurostat2, pig.farmer.df, eurostat)


# ------------------------ 4. Other information about study respondents and survey ------------------------


# What proportion of the surveys were done on paper vs on tablets?
# 158 electronic, 5 on paper.
table(no.qu.df$survey_administered)


# How long did it take people to complete the electronic survey?
# ... median is 18 minutes.
# ... Nb that many of the longer time periods are not accurate:
# ... it is likely that the survey was registered as started (i.e. the start page opened by one of the survey team),
# ... but it was a number of minutes before they enticed a visitor to complete it 
# >>> To avoid upwardly biasing the estimated length, I therefore quote the median.
summary(no.qu.df$survey_duration_num)


# ------------------------------------- 5. Plot S2 Appendix Fig. 3 ---------------------------------------------------


# For each source of food waste, calculate the mean response per respondent and plot
cols_to_keep <- c("Q1_1", "Q1_2", "Q1_3","Q1_4", "Q1_5",   # list of questions to be inlcuded
                  "Q1_6", "Q1_7", "Q1_8", "Q1_9", "Q1_10",
                  "Q2_1", "Q2_2", "Q2_3","Q2_4", "Q2_5",   # list of questions 
                  "Q2_6", "Q2_7", "Q2_8", "Q2_9", "Q2_10",
                  "Q3_1", "Q3_2", "Q3_3","Q3_4", "Q3_5",   # list of questions 
                  "Q3_6", "Q3_7", "Q3_8", "Q3_9", "Q3_10")
qu.diff.feeds <- qu.df[,cols_to_keep]                      # make a df of the relevant questions
qu.diff.feeds2 <- as.data.frame(t(qu.diff.feeds))          # Transpose so questions are in rows
qu.diff.feeds2$feedstuff <- row.names(qu.diff.feeds2)      # Make a new column with the question identifiers
colnames(qu.diff.feeds2)[1] <- "question"                  # give the first col a name
qu.diff.feeds2$feedstuff_name <- sub(".* - ", "", qu.diff.feeds2$question) # Make a new column with the name of the feedstuff
head(qu.diff.feeds2)                                                       # Look at the data
qu.diff.feeds3 <- qu.diff.feeds2[,c("feedstuff", "feedstuff_name")]        # drop extra column
no.qu.df$uniqueID <- 1:nrow(no.qu.df)                                      # Add in a column to identify each respondent. This is needed later when using spread: https://stackoverflow.com/questions/25960394/unexpected-behavior-with-tidyr
no.qu.diff.feeds_v2 <- gather(no.qu.df[,c("uniqueID",cols_to_keep, "job_group")], feedstuff, acceptability, Q1_1:Q3_10) # convert wide to long
no.qu.diff.feeds_v3 <- join(no.qu.diff.feeds_v2, qu.diff.feeds3, by = "feedstuff", type="inner")
dim(no.qu.diff.feeds_v2); dim(no.qu.diff.feeds_v3)              # check no rows were lost
no.qu.diff.feeds_v3$question_block <- gsub("_.*$", "", no.qu.diff.feeds_v3$feedstuff )
head(no.qu.diff.feeds_v2)


# convert qualitative likert scores to quantitative (1-5 scale)
translate_df <- as.data.frame(                             # Make a df that I use to convert qualitative scores <- quantitative
  cbind(rep(1:5,3),
        c("Very negative", "Negative", "Neither positive nor negative", "Positive", "Very positive",
          "Very uncomfortable", "Uncomfortable", "Neither comfortable nor uncomfortable", "Comfortable", "Very comfortable",
          "Very dissatisfied", "Dissatisfied", "Neither satisfied nor dissatisfied", "Satisfied", "Very satisfied")))
colnames(translate_df) <- c("likert.quant","acceptability")
translate_df                                               # look at the data
no.qu.diff.feeds_v4 <- join(no.qu.diff.feeds_v3, translate_df, by = "acceptability", type = "inner") # join the dfs
dim(no.qu.diff.feeds_v3); dim(no.qu.diff.feeds_v4)             # ensure no rows lost
rm(qu.diff.feeds, qu.diff.feeds2, qu.diff.feeds3,          # tidy up
   no.qu.diff.feeds_v2, no.qu.diff.feeds_v3, translate_df, 
   cols_to_keep)


# Spread df, so one row per respondent
head(no.qu.diff.feeds_v4); dim(no.qu.diff.feeds_v4)
class(no.qu.diff.feeds_v4$likert.quant)      # is the likert scale numeric? No.
no.qu.diff.feeds_v4$likert.quant <- as.numeric(as.character(no.qu.diff.feeds_v4$likert.quant)) # Convert to numeric
class(no.qu.diff.feeds_v4$likert.quant)      # is the likert scale numeric? Yes.
no.qu.diff.feeds_v5 <- as.data.frame( 
  no.qu.diff.feeds_v4 %>%                    # Take the feed acceptability data
    group_by(feedstuff_name) %>%             # for each feedstuff, 
    select(-c(feedstuff, acceptability)) %>% # remove unwanted columns
    spread(question_block,likert.quant))     # make each of the triplicate questions a column
head(no.qu.diff.feeds_v5)


# Calculate the mean per respondent/feed 
# ... This can be done, since there was very high internal consistency (alpha ~ 0.95)
# ... among the three different constructs of feed acceptability
no.qu.diff.feeds_v5$mean_Q123 <- (no.qu.diff.feeds_v5$Q1 + no.qu.diff.feeds_v5$Q2 + no.qu.diff.feeds_v5$Q3) /3
head(no.qu.diff.feeds_v5)


# Calculate the grand mean per source of food losses and job group
overall_acceptability <- no.qu.diff.feeds_v5 %>%
  group_by(job_group, feedstuff_name) %>%
  summarize(N = n(),                       # n() function in dplyr counts the number of rows
            grand_mean = mean(mean_Q123),
            grand_sd = sd(mean_Q123),      # Calculate the standard deviation 
            grand_se = grand_sd / sqrt(N)) # And the standard error of the mean


# Tidy up the names of the feeds, so they're on multiple lines when plotted
overall_acceptability$feedstuff_name <- factor(overall_acceptability$feedstuff_name, 
                                               levels = rev(c("Heat-treated household food leftovers",    
                                                              "Heat-treated, unsold chicken sandwiches from supermarkets",
                                                              "Heat-treated, unsold bacon sandwiches from supermarkets",
                                                              "Misshapen chocolates from chocolate factories",
                                                              "Heat-treated leftovers from a college canteen",
                                                              "Unsold egg sandwiches from supermarkets",
                                                              "Unsold bread from supermarkets", 
                                                              "Unsold confectionary containing porcine gelatine",
                                                              "Biscuit crumbs from biscuit factories",
                                                              "Heat-treated restaurant leftovers")), 
                                               labels = rev(c("Heat-treated household\nfood leftovers", 
                                                              "Heat-treated, unsold chicken\nsandwiches from supermarkets",
                                                              "Heat-treated, unsold bacon\nsandwiches from supermarkets",
                                                              "Misshapen chocolates\nfrom chocolate factories",
                                                              "Heat-treated leftovers\nfrom a college canteen",
                                                              "Unsold egg sandwiches\nfrom supermarkets",
                                                              "Unsold bread\nfrom supermarkets",
                                                              "Unsold confectionary\ncontaining porcine gelatine",
                                                              "Biscuit crumbs from\nbiscuit factories",
                                                              "Heat-treated restaurant\nleftovers")))


# Tidy up the names of the feeds, so they're on multiple lines when plotted
no.qu.diff.feeds_v5$feedstuff_name <- factor(no.qu.diff.feeds_v5$feedstuff_name, 
                                             levels = rev(c("Heat-treated household food leftovers",    
                                                            "Heat-treated, unsold chicken sandwiches from supermarkets",
                                                            "Heat-treated, unsold bacon sandwiches from supermarkets",
                                                            "Misshapen chocolates from chocolate factories",
                                                            "Heat-treated leftovers from a college canteen",
                                                            "Unsold egg sandwiches from supermarkets",
                                                            "Unsold bread from supermarkets", 
                                                            "Unsold confectionary containing porcine gelatine",
                                                            "Biscuit crumbs from biscuit factories",
                                                            "Heat-treated restaurant leftovers")), 
                                             labels = rev(c("Heat-treated household\nfood leftovers", 
                                                            "Heat-treated, unsold chicken\nsandwiches from supermarkets",
                                                            "Heat-treated, unsold bacon\nsandwiches from supermarkets",
                                                            "Misshapen chocolates\nfrom chocolate factories",
                                                            "Heat-treated leftovers\nfrom a college canteen",
                                                            "Unsold egg sandwiches\nfrom supermarkets",
                                                            "Unsold bread\nfrom supermarkets",
                                                            "Unsold confectionary\ncontaining porcine gelatine",
                                                            "Biscuit crumbs from\nbiscuit factories",
                                                            "Heat-treated restaurant\nleftovers")))


# S2App Figure 3: Density plot to show the spread in the data
S2AppFig3 <- ggplot(no.qu.diff.feeds_v5, aes(mean_Q123-1, col= job_group, fill=job_group)) + 
    geom_density(alpha = 0.4) + 
  xlab("Mean acceptability of sources of pig feed\n(scored from 1-5, Very unacceptable - Very acceptable)") + 
  ylab("") + 
    theme_bw() + 
  theme(legend.position = c(0.85, 0.15), 
        legend.title = element_blank(),
        axis.text = element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=14), 
        strip.text = element_text(size=14), 
        strip.background = NULL)  + 
  scale_x_continuous(breaks=0:4, 
                     limits = c(0,4),
                     labels=0:4 + 1) + 
  geom_vline(data = as.data.frame(overall_acceptability), aes(xintercept=grand_mean - 1, col = job_group), size=1, linetype="dashed") + 
  facet_wrap(~feedstuff_name,
             ncol = 3,) +
  scale_fill_manual(values = c("#D55E00","#0072B2"),
                    guide = guide_legend(reverse=TRUE,
                                         direction = "horizontal")) +
  scale_color_manual(values = c("#D55E00","#0072B2"),
                    guide = guide_legend(reverse=TRUE)) +
  guides(col=FALSE)
S2AppFig3


# tidy up
rm(overall_acceptability, no.qu.diff.feeds_v5, S2AppFig3)


# ------------------------------- 6. Models of acceptability of different food losses as feed ----------------------------------


# Compare models including information about the feedstuffs:
# ... i.e. legal status, whether or not it contains ABPs, whehther or not intra-species recycling is possible
# ... with a model without those predictors
# Both models include both respondent and feedstuff random intercepts 
# ... the respondent_ID is included because we have three measures per feedstuff:respondent combination
# ... (from the three different constructs of acceptability)


# Load the rethinking package
library(rethinking)


# Prepare the dataset 
head(no.qu.diff.feeds_v4); dim(no.qu.diff.feeds_v4) # one row per question (i.e. 3 per respondent:feedstuff combination) 


# i) Add in labels for feed sources that are legal/illegal
# ... legal = 1, illegal = 0.
legal_feeds <- c("Biscuit crumbs from biscuit factories", "Misshapen chocolates from chocolate factories",
                 "Unsold bread from supermarkets", "Unsold confectionary containing porcine gelatine",
                 "Unsold egg sandwiches from supermarkets")
no.qu.diff.feeds_v4$legality <- ifelse(no.qu.diff.feeds_v4$feedstuff_name %in% legal_feeds, 1, 0)
rm(legal_feeds)    # tidy up 


# ii) Add in labels for ABP/noABP
# ABP feed = 1, non-ABP feed = 0
ABP_feeds <- c("Unsold confectionary containing porcine gelatine",
               "Unsold egg sandwiches from supermarkets",
               "Heat-treated household food leftovers",                    
               "Heat-treated leftovers from a college canteen",
               "Heat-treated restaurant leftovers",           
               "Heat-treated, unsold bacon sandwiches from supermarkets",
               "Heat-treated, unsold chicken sandwiches from supermarkets")
no.qu.diff.feeds_v4$ABP <- ifelse(no.qu.diff.feeds_v4$feedstuff_name %in% ABP_feeds, 1, 0)
rm(ABP_feeds)        # tidy up


# iii) Add in labels for intra-species recycling/not
# 1 = some intra-species recycling, 0 = no risk of intrascpecies recycling
intra_spp_feeds <- c("Heat-treated household food leftovers",                    
                     "Heat-treated leftovers from a college canteen",
                     "Heat-treated restaurant leftovers",           
                     "Heat-treated, unsold bacon sandwiches from supermarkets",
                     "Unsold confectionary containing porcine gelatine")
no.qu.diff.feeds_v4$intra_spp <- ifelse(no.qu.diff.feeds_v4$feedstuff_name %in% intra_spp_feeds, 1, 0)
rm(intra_spp_feeds)   # tidy up
head(no.qu.diff.feeds_v4)


# Use the CUMULATIVE LINK function (used later for model starting values)
logit <- function(x) log(x/(1-x)) # convenience function


# discrete proportion of each response value
pr_k <- table( no.qu.diff.feeds_v4$likert.quant ) / nrow(no.qu.diff.feeds_v4)
cum_pr_k <- cumsum( pr_k )                  # cumsum converts to cumulative proportions
logit( cum_pr_k )                           # the logit intercept parameter estimate (averaged across all feeds - this is just used for starting values)


# Make a data.frame for fitting models
no.qu.diff.feeds_v4_copy <- no.qu.diff.feeds_v4 # Make a copy since I will change the columns (to make them STAN friendly)
no.qu.diff.feeds_v4_copy$job_group_ID <- ifelse(no.qu.diff.feeds_v4_copy$job_group == "Pig farmer", 1, 0)       # Stan takes only integer group types
no.qu.diff.feeds_v4_copy$feedstuff_ID <- coerce_index(no.qu.diff.feeds_v4_copy$feedstuff_name)   # Stan takes only integer group types
no.qu.diff.feeds_v4_copy$respondent_ID <- coerce_index(no.qu.diff.feeds_v4_copy$uniqueID)   # Stan takes only integer group types
no.qu.diff.feeds_v4stan <- no.qu.diff.feeds_v4_copy[,c("respondent_ID", "job_group_ID","feedstuff_ID","likert.quant", 
                                                       "legality", "ABP", "intra_spp" )]                                 # Keep only cols required for this model
colnames(no.qu.diff.feeds_v4stan)[which(colnames(no.qu.diff.feeds_v4stan) == "likert.quant")] <- "likertQuant" # remove "." from variable names
class(no.qu.diff.feeds_v4stan$likertQuant) # check that these are integers. Nope. 
no.qu.diff.feeds_v4stan$likertQuant <- as.integer(as.character(no.qu.diff.feeds_v4stan$likertQuant)) # convert to integer
class(no.qu.diff.feeds_v4stan$likertQuant) # check that these are integers Yes.
class(no.qu.diff.feeds_v4stan$job_group_ID) # check that these are integers. Nope. 
no.qu.diff.feeds_v4stan$job_group_ID <- as.integer(as.character(no.qu.diff.feeds_v4stan$job_group_ID)) # convert to integer
class(no.qu.diff.feeds_v4stan$job_group_ID) # check that these are integers Yes.
table_feed_codes <- as.data.frame(no.qu.diff.feeds_v4_copy %>%                    # Make a mini table which has the feedstuff ID and name
                                    group_by(feedstuff_name, feedstuff_ID) %>%    # ... this is used later to interpret results
                                    summarize(mean_likert = mean(likert.quant)) ) 
rm(no.qu.diff.feeds_v4_copy)                # tidy up


# Model AC1
# Fit model with respondent,FEED*job, ABP/intra-spp/legality & legality:job interaction
mAC1 <- map2stan(
  alist(
    likertQuant ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    
    # linear models
    phi <- aF[feedstuff_ID] + aR[respondent_ID] +  # summary linear model
      bJF[feedstuff_ID]*job_group_ID + 
      bL*legality + bI*intra_spp + bA*ABP +
      bJL*legality*job_group_ID,   
    
    # adaptive priors
    c(aF,bJF)[feedstuff_ID] ~
      dmvnorm2( c(0, bJ), sigma_F, Rho_F ),
    c(aR)[respondent_ID] ~ dnorm(0, sigma_R),
    
    # Fixed priors
    c(bJ, bL, bI, bJL, bA) ~ dnorm(0,1),                  # priors for slope
    sigma_F ~ dcauchy(0,2),                       # priors for feedstuff_ID prior
    sigma_R ~ dcauchy(0,2),                       # priors for respondent_ID prior
    Rho_F ~ dlkjcorr(4),                          # prior for the correlation matrix
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=no.qu.diff.feeds_v4stan,
  start=list(cutpoints = c(logit(cum_pr_k)[[1]], logit(cum_pr_k)[[2]],     # base starting values on overall mean
                           logit(cum_pr_k)[[3]], logit(cum_pr_k)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Look at m1 chains and coefficients
# Model with varying intercepts for respondent has MUCH lower WAIC
# ... than model without
# plot(mAC1)                  # Check the trace plot. Looks OK.
# precis(mAC1, depth=2)


# model AC2
# Fit model with respondent, FEED*job and legal/ABP/intra-spp
mAC2 <- map2stan(
  alist(
    likertQuant ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    
    # linear models
    phi <- aF[feedstuff_ID] + aR[respondent_ID] +  # summary linear model
      bJF[feedstuff_ID]*job_group_ID + 
      bL*legality + bA*ABP + bI*intra_spp,   
    
    # adaptive priors
    c(aF, bJF)[feedstuff_ID] ~
      dmvnorm2( c(0, bJ), sigma_F, Rho_F ),
    c(aR)[respondent_ID] ~ dnorm(0, sigma_R),
    
    # Fixed priors
    c(bJ,bL,bA,bI) ~ dnorm(0,1),                  # priors for slope
    sigma_F ~ dcauchy(0,2),                       # priors for feedstuff_ID prior
    sigma_R ~ dcauchy(0,2),                       # priors for respondent_ID prior
    Rho_F ~ dlkjcorr(4),                          # prior for the correlation matrix
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=no.qu.diff.feeds_v4stan,
  start=list(cutpoints = c(logit(cum_pr_k)[[1]], logit(cum_pr_k)[[2]],     # base starting values on overall mean
                           logit(cum_pr_k)[[3]], logit(cum_pr_k)[[4]]) ),
  iter = 10000, chains=1, cores = 1)



# Look at m2 chains and coefficients
# plot(mAC2)                  # Check the trace plot.
# precis(mAC2, depth=2)


# model AC3
# Fit model with respondent, FEED*job and ABP/intra-spp
mAC3 <- map2stan(
  alist(
    likertQuant ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    
    # linear models
    phi <- aF[feedstuff_ID] + aR[respondent_ID] +  # summary linear model
      bJF[feedstuff_ID]*job_group_ID + 
      bA*ABP + bL*legality,   
    
    # adaptive priors
    c(aF, bJF)[feedstuff_ID] ~
      dmvnorm2( c(0, bJ), sigma_F, Rho_F ),
    c(aR)[respondent_ID] ~ dnorm(0, sigma_R),
    
    # Fixed priors
    c(bJ,bA,bL) ~ dnorm(0,1),                  # priors for slope
    sigma_F ~ dcauchy(0,2),                       # priors for feedstuff_ID prior
    sigma_R ~ dcauchy(0,2),                       # priors for respondent_ID prior
    Rho_F ~ dlkjcorr(4),                          # prior for the correlation matrix
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=no.qu.diff.feeds_v4stan,
  start=list(cutpoints = c(logit(cum_pr_k)[[1]], logit(cum_pr_k)[[2]],     # base starting values on overall mean
                           logit(cum_pr_k)[[3]], logit(cum_pr_k)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check m3 chains
# plot(mAC3)                  # Check the trace plot.
# precis(mAC3, depth=2)


# model AC4
# Fit model with respondent, FEED*job and legal status/intra-spp
mAC4 <- map2stan(
  alist(
    likertQuant ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    
    # linear models
    phi <- aF[feedstuff_ID] + aR[respondent_ID] +  # summary linear model
      bJF[feedstuff_ID]*job_group_ID + 
      bL*legality + bI*intra_spp,   
    
    # adaptive priors
    c(aF,bJF)[feedstuff_ID] ~
      dmvnorm2( c(0, bJ), sigma_F, Rho_F ),
    c(aR)[respondent_ID] ~ dnorm(0, sigma_R),
    
    # Fixed priors
    c(bJ,bL, bI) ~ dnorm(0,1),                  # priors for slope
    sigma_F ~ dcauchy(0,2),                       # priors for feedstuff_ID prior
    sigma_R ~ dcauchy(0,2),                       # priors for respondent_ID prior
    Rho_F ~ dlkjcorr(4),                          # prior for the correlation matrix
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=no.qu.diff.feeds_v4stan,
  start=list(cutpoints = c(logit(cum_pr_k)[[1]], logit(cum_pr_k)[[2]],     # base starting values on overall mean
                           logit(cum_pr_k)[[3]], logit(cum_pr_k)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Look at chains
# plot(mAC4)                  # Check the trace plot.
# precis(mAC4, depth=2)


# model AC5
# Fit model with respondent, FEED*job and ABP/intra-spp
mAC5 <- map2stan(
  alist(
    likertQuant ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    
    # linear models
    phi <- aF[feedstuff_ID] + aR[respondent_ID] +  # summary linear model
      bJF[feedstuff_ID]*job_group_ID + 
      bA*ABP + bI*intra_spp,   
    
    # adaptive priors
    c(aF, bJF)[feedstuff_ID] ~
      dmvnorm2( c(0, bJ), sigma_F, Rho_F ),
    c(aR)[respondent_ID] ~ dnorm(0, sigma_R),
    
    # Fixed priors
    c(bJ,bA,bI) ~ dnorm(0,1),                  # priors for slope
    sigma_F ~ dcauchy(0,2),                       # priors for feedstuff_ID prior
    sigma_R ~ dcauchy(0,2),                       # priors for respondent_ID prior
    Rho_F ~ dlkjcorr(4),                          # prior for the correlation matrix
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=no.qu.diff.feeds_v4stan,
  start=list(cutpoints = c(logit(cum_pr_k)[[1]], logit(cum_pr_k)[[2]],     # base starting values on overall mean
                           logit(cum_pr_k)[[3]], logit(cum_pr_k)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Look at chains
# plot(mAC5)                  # Check the trace plot.
# precis(mAC5, depth=2)


# Fit model AC6
# ... model without characteristics of feed
mAC6 <- map2stan(
  alist(
    likertQuant ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    
    # linear models
    phi <- aF[feedstuff_ID] + aR[respondent_ID] +  # summary linear model
      bJF[feedstuff_ID]*job_group_ID,   
    
    # adaptive priors
    c(aF, bJF)[feedstuff_ID] ~
      dmvnorm2( c(0, bJ), sigma_F, Rho_F ),
    c(aR)[respondent_ID] ~ dnorm(0, sigma_R),
    
    # Fixed priors
    bJ ~ dnorm(0,1),                              # priors for slope
    sigma_F ~ dcauchy(0,2),                       # priors for feedstuff_ID prior
    sigma_R ~ dcauchy(0,2),                       # priors for respondent_ID prior
    Rho_F ~ dlkjcorr(4),                          # prior for the correlation matrix
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=no.qu.diff.feeds_v4stan,
  start=list(cutpoints = c(logit(cum_pr_k)[[1]], logit(cum_pr_k)[[2]],     # base starting values on overall mean
                           logit(cum_pr_k)[[3]], logit(cum_pr_k)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Look at chains
# Model with varying intercepts for respondent has MUCH lower WAIC
# ... than model without
# plot(mAC6)                  # Check the trace plot. Looks OK.
# precis(mAC6, depth=2)


# Model AC7
# Fit model with respondent, ABP/intra-spp/legality/job & legality:job interaction
# ... i.e. without feed:job interaction
mAC7 <- map2stan(
  alist(
    likertQuant ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    
    # linear models
    phi <- aF[feedstuff_ID] + aR[respondent_ID] +  # summary linear model
      bL*legality + bI*intra_spp + bA*ABP +
      bJ*job_group_ID + 
      bJL*legality*job_group_ID,   
    
    # adaptive priors
    aF[feedstuff_ID]  ~ dnorm(0, sigma_F), 
    aR[respondent_ID] ~ dnorm(0, sigma_R),
    
    # Fixed priors
    c(bJ, bL, bI, bJL, bA) ~ dnorm(0,1),                  # priors for slope
    c(sigma_F,sigma_R) ~ dcauchy(0,2),                       # priors for feedstuff_ID prior
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=no.qu.diff.feeds_v4stan,
  start=list(cutpoints = c(logit(cum_pr_k)[[1]], logit(cum_pr_k)[[2]],     # base starting values on overall mean
                           logit(cum_pr_k)[[3]], logit(cum_pr_k)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Model with varying intercepts for respondent has MUCH lower WAIC
# ... than model without
# plot(mAC7)                  # Check the trace plot. Looks OK.
# precis(mAC7, depth=2)


# Compare models
# ... nb that because the WAIC of models 2-6 are very similar, 
# ... their relative model weights may differ between runs
compare(mAC1, mAC2, mAC3, mAC4, mAC5, mAC6, mAC7)
coeftab(mAC1, mAC2, mAC3, mAC4, mAC5, mAC6, mAC7)


# Plot the coefficients for all models with weightings >0
# ... i.e. not model mAC7 (which has weight of 0.00)
# >>> here I extract the parameter values, confidence intervals, and weights from the models
comparison <- compare(mAC1, mAC2, mAC3, mAC4, mAC5, mAC6, mAC7)
mAC1_res <- precis(mAC1, depth=2, pars=c("bA", "bI", "bL","sigma_F","sigma_R", "bJL"))@output
mAC1_res$coef <- row.names(mAC1_res)    # Extract the coefficient names
colnames(mAC1_res) <- c("Mean","StdDev", "lwr.89", "upr.89", "n_eff", "Rhat", "coef")
mAC1_res$model <- "mAC1"
mAC1_res$weighting <- comparison@output[which(row.names(comparison@output) == "mAC1"),"weight"]
mAC2_res <- precis(mAC2, depth=2, pars=c("bA", "bI", "bL","sigma_F","sigma_R", "bJL"))@output
mAC2_res$coef <- row.names(mAC2_res)    # Extract the coefficient names
colnames(mAC2_res) <- c("Mean","StdDev", "lwr.89", "upr.89", "n_eff", "Rhat", "coef")
mAC2_res$model <- "mAC2"
mAC2_res$weighting <- comparison@output[which(row.names(comparison@output) == "mAC2"),"weight"]
mAC3_res <- precis(mAC3, depth=2, pars=c("bA", "bI", "bL","sigma_F","sigma_R", "bJL"))@output
mAC3_res$coef <- row.names(mAC3_res)    # Extract the coefficient names
colnames(mAC3_res) <- c("Mean","StdDev", "lwr.89", "upr.89", "n_eff", "Rhat", "coef")
mAC3_res$model <- "mAC3"
mAC3_res$weighting <- comparison@output[which(row.names(comparison@output) == "mAC3"),"weight"]
mAC4_res <- precis(mAC4, depth=2, pars=c("bA", "bI", "bL","sigma_F","sigma_R"))@output
mAC4_res$coef <- row.names(mAC4_res)    # Extract the coefficient names
colnames(mAC4_res) <- c("Mean","StdDev", "lwr.89", "upr.89", "n_eff", "Rhat", "coef")
mAC4_res$model <- "mAC4"
mAC4_res$weighting <- comparison@output[which(row.names(comparison@output) == "mAC4"),"weight"]
mAC5_res <- precis(mAC5, depth=2, pars=c("bA", "bI", "bL","sigma_F","sigma_R"))@output
mAC5_res$coef <- row.names(mAC5_res)    # Extract the coefficient names
colnames(mAC5_res) <- c("Mean","StdDev", "lwr.89", "upr.89", "n_eff", "Rhat", "coef")
mAC5_res$model <- "mAC5"
mAC5_res$weighting <- comparison@output[which(row.names(comparison@output) == "mAC5"),"weight"]
mAC6_res <- precis(mAC6, depth=2, pars=c("bA", "bI", "bL","sigma_F","sigma_R"))@output
mAC6_res$coef <- row.names(mAC6_res)    # Extract the coefficient names
colnames(mAC6_res) <- c("Mean","StdDev", "lwr.89", "upr.89", "n_eff", "Rhat", "coef")
mAC6_res$model <- "mAC6"
mAC6_res$weighting <- comparison@output[which(row.names(comparison@output) == "mAC6"),"weight"]
models_res <- rbind( mAC1_res, mAC2_res, mAC3_res, mAC4_res, mAC5_res, mAC6_res)


# Rebael coef for plotting
models_res$coef <- factor(models_res$coef,
                          levels=rev(c("bA", "bI", "bL","bJL", "sigma_F[1]", "sigma_F[2]", "sigma_R")),
                          labels=rev(c("Feeds containing animal by-products", "Feeds with potential for\nintra-species recycling", "Legality (legal feeds)",
                                       "Legal feeds:pig farmer interaction",
                                       "Standard deviation between feeds", "Standard deviation of acceptability of feeds\nbetween pig farmers and others", "Standard deviation between respondents")))


# save model output
# write.table(models_res, "results_tables/20171007_Fig1_output.csv", sep=",", row.names=F)
# models_res <- read.csv("results_tables/20171007_Fig1_output.csv", header=T)


# tidy up
rm(comparison, mAC1_res, mAC2_res, mAC3_res, 
   mAC4_res, mAC5_res, mAC6_res)


# 6i) Plot model output: Fig 1 ---- 
myCols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00")            # Make colour palette
models_res$model <- factor(models_res$model,           # Order the models, in order of their weighting
                                levels=rev(c("mAC1", "mAC2","mAC3","mAC4", "mAC5", "mAC6")),
                                labels=rev(c("AC1", "AC2","AC3","AC4", "AC5", "AC6")))

# Nb reverse the order for coefficients, so sigmas are plotted last
models_res$coef <- factor(models_res$coef,
                          levels=rev(c("Feeds containing animal by-products", "Feeds with potential for\nintra-species recycling", "Legality (legal feeds)",
                                       "Legal feeds:pig farmer interaction",
                                       "Standard deviation between feeds", "Standard deviation of acceptability of feeds\nbetween pig farmers and others", "Standard deviation between respondents")))

Fig1a <- ggplot(models_res, aes(coef, Mean, col = model, size = weighting*100)) + 
  geom_point(position=position_dodge(.45)) + 
  geom_errorbar(aes(ymin=lwr.89, ymax=upr.89, col = model),
                size=1,    # Thinner lines
                width=.2,
                position=position_dodge(.45)
  )  +
  coord_flip() + 
  geom_hline(yintercept = 0) + 
  ylab("Coefficient value") + 
  xlab("") + 
  scale_color_manual(values = myCols, 
                     # name = "Models with\ngreatest weighting", 
                     name = "Model ID", 
                     guide = guide_legend(reverse=TRUE)) +
  # guides(fill=guide_legend(reverse = T)) + 
  scale_radius(name="Model\nweight (%)",
               limits = c(8,42), 
               guide=FALSE)  + # remove size legend
  theme_bw()  +
  theme(legend.position = c(0.12,0.3), 
        plot.margin = unit(c(1,1,1.5,1.2),"cm"), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        legend.text = element_text(size=10), 
        legend.title = element_text(size=12, face="bold"), 
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        panel.border=element_blank(),
        axis.line = element_line(), # axis.line adds in the x and y axis lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  + 
  geom_vline(xintercept = 3.5, linetype = "dashed", col="darkgray")
Fig1a


# plot again, with size legend
Fig1b <- ggplot(models_res, aes(coef, Mean, col = model, size = weighting*100)) + geom_point(position=position_dodge(.45)) + 
  geom_errorbar(aes(ymin=lwr.89, ymax=upr.89, col = model),
                size=1,    # Thinner lines
                width=.2,
                position=position_dodge(.45)
  )  +
  coord_flip() + 
  geom_hline(yintercept = 0) + 
  ylab("Coefficient value") + 
  xlab("") + 
  scale_color_manual(values = myCols, 
                     # name = "Models with\ngreatest weighting", 
                     name = "Model ID", 
                     guide = FALSE) +
  scale_size_continuous(limits = c(8,42),              # standardise the weights displayed in the legend
                        name="Model\nweight (%)")  + # specify name for size legend
  theme_bw()  +
  theme(legend.position = c(0.15,0.435), 
        plot.margin = unit(c(1,1,1.5,1.2),"cm"), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        legend.text = element_text(size=10), 
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        legend.title = element_text(size=12, face="bold"))  + 
  geom_vline(xintercept = 3.5, linetype = "dashed", col="darkgray")


# Step 2) Extract the size legend - leg2
leg <- gtable_filter(ggplot_gtable(ggplot_build(Fig1b)), "guide-box") 


# Step 3) Plot and add second legend
Fig1 <- Fig1a + 
  annotation_custom(grob = leg, 
                    ymin = -1.1, ymax = 0,    # Place legend
                    xmin = 1.18, xmax = 5)      # Nb because of coord_flip() the x and y axes are flipped
grid.newpage()
grid.draw(Fig1)


# Step 4) Add plot title
g = ggplotGrob(Fig1)
g = gtable_add_grob(g, grobTree(textGrob("Lower acceptability", x=0, hjust=0),
                                textGrob("Greater acceptability", x=1, hjust=1)),
                    t=1, l=4)
grid.draw(g)


# tidy up
rm(Fig1a, Fig1b, Fig1, myCols, models_res, g, leg)


# 6ii) Plot S2 Appendix Figure 4, acceptability of feeds with ABPs or intra-spp recycling -----


# Plot effect size for inclusion of animal by-products (ABP) on simulated likert score
# ... to do this, I create a data.frame which is used to generate model output
# ... for 1000 random respondents and randomly selected combinations of other variables 
d_pred_feed_chrc_ABP <- data.frame(                          # Make data.frame for feeds with ABPs
  likertQuant=5,                                             # needed so sim knows number of possible levels
  respondent_ID= sample(no.qu.diff.feeds_v4stan$respondent_ID,100, replace=T),          # Simulated for 100 different respondents
  job_group_ID= sample(c(0,1), 100, replace=T),              # Respondents are either pig farmer or other
  feedstuff_ID= sample(c(2,3,4,5,6,9,10), 100, replace=T) ,  # list of feedstuffs containing ABPs 
  legality=sample(c(0,1), 100, replace=T),                   # Feed is either legal or not
  ABP = 1, 
  intra_spp = rep(0, 100))                                   # To measure only effect size of ABP, this must be the same as in d_pred_feed_chrc_noABP
d_pred_feed_chrc_noABP <- data.frame(                        # Make data.frame for feeds without ABPs
  likertQuant=5,                                             # needed so sim knows number of possible levels
  respondent_ID=sample(no.qu.diff.feeds_v4stan$respondent_ID,100, replace=T), # Simulated for 100 different respondents
  job_group_ID=sample(c(0,1), 100, replace=T),                 # Respondents are either pig farmer or other
  feedstuff_ID=sample(c(1,7,8), 100, replace=T),                 # list of feedstuffs not containing ABPs 
  legality= sample(c(0,1), 100, replace=T),                  # Feed is either legal or not
  ABP = 0, 
  intra_spp = rep(0, 100))                                   # since feed does not contain ABP, it cannot cannot have intra-spp recycling either
d_pred_feed_chrc_ABP_sim <- ensemble_custom(mAC1, mAC2, mAC3, mAC4, mAC5, mAC6, data=d_pred_feed_chrc_ABP)      # simulate data 
d_pred_feed_chrc_ABP_sim_mu <- apply(d_pred_feed_chrc_ABP_sim$sim,2,mean)      # Calculate the mean Likert score per respondent
d_pred_feed_chrc_noABP_sim <- ensemble_custom(mAC1, mAC2, mAC3, mAC4, mAC5, mAC6, data=d_pred_feed_chrc_noABP)
str(d_pred_feed_chrc_noABP_sim)
d_pred_feed_chrc_noABP_sim_mu <- apply(d_pred_feed_chrc_noABP_sim$sim,2,mean)
compare_ABP <- as.data.frame(cbind(d_pred_feed_chrc_noABP_sim_mu, d_pred_feed_chrc_ABP_sim_mu))
colnames(compare_ABP) <- c("no_ABP", "ABP")
compare_ABP_long <- gather(compare_ABP, ABP_category, likert_score, no_ABP:ABP)
compare_ABP_long$ABP_category <- factor(compare_ABP_long$ABP_category, 
                                        levels=c("no_ABP", "ABP"),
                                        labels=c("Feed not containing\nABP","Feed containing\nABP")) 

# Plot effect of ABP inclusion on feed acceptability
compare_ABP_long2 <- compare_ABP_long %>%   # Make a mini-df which contains the means in each group
  group_by(ABP_category) %>%                # ... this is used for plotting
  summarize(mean_per_group = mean(likert_score))
myCols2 <- c("#009E73","#999999")     # Make colour palette
pABP <- 
  ggplot(compare_ABP_long, aes(likert_score, fill=ABP_category)) + 
  geom_density(alpha=0.4) + 
  theme_bw() + 
  xlab("Mean simulated likert score") + 
  ylab("Response density") + 
  theme(axis.title = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12), 
        legend.position=c(0.25,0.85), 
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.key.size = unit(2, 'lines')) + 
    xlim(1,5) + 
  annotate("text", label = "a", x = 1, y = .7, size = 8,fontface =2) +
  geom_vline(data = compare_ABP_long2, aes(xintercept = mean_per_group,  # add in vertica lines at mean
             col = ABP_category), size=1, linetype="dashed") +
  scale_fill_manual(values = myCols2) +
  scale_color_manual(values = myCols2)
pABP
mean(d_pred_feed_chrc_noABP_sim_mu)      # Compare the difference in the mean between feed without ABP
mean(d_pred_feed_chrc_ABP_sim_mu)        # ... and feed with ABP


# tidy up
rm(d_pred_feed_chrc_ABP, d_pred_feed_chrc_noABP,
   d_pred_feed_chrc_ABP_sim, d_pred_feed_chrc_noABP_sim,
   d_pred_feed_chrc_ABP_sim_mu, d_pred_feed_chrc_noABP_sim_mu,
   compare_ABP, compare_ABP_long, compare_ABP_long2)


# And calculate the effect size for the difference between feeds
# ... with potential for intra-species recycling, and those without
# ... to do this, I create a data.frame which is used to generate model output
# ... for 100 random respondents and randomly selected combinations of other variables 
d_pred_feed_chrc_IS <- data.frame(                           # Make data.frame for feeds with intra-spp recycling
  likertQuant=5,                                             # needed so sim knows number of possible levels
  respondent_ID= sample(no.qu.diff.feeds_v4stan$respondent_ID,1000, replace=T),          # Simulated for 100 different respondents
  job_group_ID= sample(c(0,1), 1000, replace=T),              # Respondents are either pig farmer or other
  feedstuff_ID= sample(c(2,3,4,5,6,9,10), 1000, replace=T) ,  # list of feedstuffs containing ABPs 
  legality=sample(c(0,1), 1000, replace=T),                   # Feed is either legal or not
  ABP = 1,                                                   # must contain ABP if potential for intra-spp recycling
  intra_spp = rep(1, 1000))                                   # Feed with intra-spp recycling
d_pred_feed_chrc_noIS <- data.frame(                         # Make data.frame for feeds without ABPs
  likertQuant=5,                                             # needed so sim knows number of possible levels
  respondent_ID=sample(no.qu.diff.feeds_v4stan$respondent_ID,1000, replace=T), # Simulated for 100 different respondents
  job_group_ID=rep(c(0,1), length.out= 1000),                 # Respondents are either pig farmer or other
  feedstuff_ID=rep(c(1,7,8),length.out=1000),                 # list of feedstuffs not containing ABPs 
  legality= sample(c(0,1), 1000, replace=T),                  # Feed is either legal or not
  ABP = 1,                                                   # must contain ABP to measure effect-size of intra-spp recycling
  intra_spp = rep(0, 1000))                                   # No intra-spp recycling
d_pred_feed_chrc_IS_sim <- ensemble_custom(mAC1,mAC2,mAC3, mAC4, mAC5, mAC6, data=d_pred_feed_chrc_IS)      # simulate data 
d_pred_feed_chrc_IS_sim_mu <- apply(d_pred_feed_chrc_IS_sim$sim,2,mean)      # Calculate the mean Likert score per respondent
d_pred_feed_chrc_noIS_sim <- ensemble_custom(mAC1, mAC2, mAC3, mAC4, mAC5, mAC6, data=d_pred_feed_chrc_noIS)
str(d_pred_feed_chrc_noIS_sim)
d_pred_feed_chrc_noIS_sim_mu <- apply(d_pred_feed_chrc_noIS_sim$sim,2,mean)
compare_IS <- as.data.frame(cbind(d_pred_feed_chrc_noIS_sim_mu, d_pred_feed_chrc_IS_sim_mu))
colnames(compare_IS) <- c("no_IS", "IS")
compare_IS_long <- gather(compare_IS, IS_category, likert_score, no_IS:IS)
compare_IS_long$IS_category <- factor(compare_IS_long$IS_category, 
                                        levels=c("no_IS", "IS"),
                                        labels=c("Feeds for which intra-species\nrecycling cannot occur","Feeds for which intra-species\nrecycling can occur")) 
# plot the effect of potential intra-species recycling on acceptability
compare_IS_long2 <- compare_IS_long %>%   # Make a mini-df which contains the means in each group
  group_by(IS_category) %>%                # ... this is used for plotting
  summarize(mean_per_group = mean(likert_score))
pIS <- 
  ggplot(compare_IS_long, aes(likert_score, fill = IS_category)) + 
  geom_density(alpha=0.4) + 
  theme_bw() + 
  xlab("Mean simulated likert score") + 
  ylab("Response density") + 
  theme(axis.title = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12), 
        legend.position=c(0.25,0.85), 
        legend.title = element_blank(),
        legend.key.size = unit(2, 'lines'), 
        legend.text= element_text(size=12)) + 
  scale_fill_manual(values = myCols2) +
  xlim(1,5) + 
  annotate("text", label = "b", x = 1, y = .52, size = 8,fontface =2) +
  geom_vline(data = compare_IS_long2, aes(xintercept = mean_per_group,  # add in vertica lines at mean
                                           col = IS_category), size=1, linetype="dashed") +
  scale_color_manual(values = myCols2)
pIS
mean(d_pred_feed_chrc_noIS_sim_mu)      # Compare the difference in the mean between feed without intra-spp recycling
mean(d_pred_feed_chrc_IS_sim_mu)        # ... and feed with intra-spp recycling


# Plot the ABP and intra-spp recycling plots together
# ... S2 Appendix Fig 4
multiplot(pABP, pIS)
rm(pABP, pIS)    # tidy up


# tidy up
rm(d_pred_feed_chrc_IS, d_pred_feed_chrc_noIS,
   d_pred_feed_chrc_IS_sim, d_pred_feed_chrc_noIS_sim,
   d_pred_feed_chrc_IS_sim_mu, d_pred_feed_chrc_noIS_sim_mu,
   compare_IS, compare_IS_long, compare_IS_long2, myCols2)


# tidy up
rm(myCols, mAC6_res, mAC5_res, mAC1_res, mAC4_res, 
   four_models_res, g, four_models_res2, table_feed_codes)


# tidy up
rm(mAC1, mAC2, mAC3, mAC4, mAC5, mAC6, mAC7, 
   no.qu.diff.feeds_v4stan, 
   no.qu.diff.feeds_v4, 
   logit, pr_k, cum_pr_k)


# ------------------------------- 7. Plot Fig. 2 - swill vs conventional feed -------------------------------



# What are the questions?
qu.df


# Compared with feeding conventional grain- and soybean-based feed, 
# ... heat-treated swill is:
head(no.qu.df)


# Look at all the variants of Q1:
# Q22_6:Q30_1


# Specify the color palette used for these plots
# my.cols <- c("gray", rev(brewer.pal(5, "RdYlGn"))) # make a custom colour palette
my.cols <- c("gray", rev(brewer.pal(5, "Spectral"))) # make a custom colour palette


# Make a theme and scale for Figure 5.
scaleFig5 <- scale_fill_manual(values = my.cols,
                    # direction = -1,
                    guide = guide_legend(reverse=TRUE,
                                         direction = "horizontal",
                                         keyheight = unit(4, units = "mm"),
                                         keywidth = unit(28 / length(labels), units = "mm"),
                                         title.position = 'top',
                                         nrow = 1,
                                         byrow = T,
                                         label.position = "bottom"
                    ))
themeFig5_v2 <-
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size=12),
        axis.text.x = element_blank(), 
        legend.text = element_text(size=10),
        strip.text = element_text(size=14), 
        strip.background = NULL)


# environmental impact
# ... Q22_6
# ... the blanks are "Don't know"
no.qu.dfE <- no.qu.df %>%      # Count the responses per job group
  group_by(job_group) %>% 
  count(Q22_6)
no.qu.dfE
levels(no.qu.dfE$Q22_6) <-  c(levels(no.qu.dfE$Q22_6), "Don't know")  # Convert blanks in data.frame to "Don't know" (which is what they were in the survey)
no.qu.dfE$Q22_6[no.qu.dfE$Q22_6==""] <- "Don't know"
unique(no.qu.dfE$Q22_6)                                               # check factor levels include "Don't know"
no.qu.dfE$Q22_6 <- factor(no.qu.dfE$Q22_6,
                          levels = c("Don't know", "Much less damaging to the environment", "Less damaging to the environment", 
                                     "Neither more nor less damaging to the environment",
                                     "More damaging to the environment", "Much more damaging to the environment"), 
                          labels = c("Don't know", "Much less\ndamaging", "Less damaging", 
                                     "Neither more nor\nless damaging",
                                     "More damaging", "Much more\ndamaging"))
no.qu.dfE$title <- "Swill is __ damaging to the environment than conventional feed"
# no.qu.dfE$job_group <- factor(no.qu.dfE$job_group, 
#                               levels=c("Other","Pig farmer"),
#                               labels=c("Other","Pig\nfarmer"))
pEnv <-  ggplot(no.qu.dfE, aes(job_group, n, fill=Q22_6)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  themeFig5_v2 + 
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 + 
  facet_wrap(~title)


# tidy up
rm(no.qu.dfE)


# Cost
no.qu.dfC <- no.qu.df %>% group_by(job_group) %>% 
  count(Q26_1)
no.qu.dfC
levels(no.qu.dfC$Q26_1) <-  c(levels(no.qu.dfC$Q26_1), "Don't know")
no.qu.dfC$Q26_1[no.qu.dfC$Q26_1==""] <- "Don't know"
unique(no.qu.dfC$Q26_1)
no.qu.dfC$Q26_1 <- factor(no.qu.dfC$Q26_1,
                          levels = c("Don't know",
                                     "Much lower cost",
                                     "Lower cost",
                                     "Costs the same",
                                     "Higher cost",
                                     "Much higher cost"), 
                          labels = c("Don't know",
                                     "Much lower\ncost",
                                     "Lower cost",
                                     "Costs the\nsame",
                                     "Higher cost",
                                     "Much higher\ncost"))
no.qu.dfC$title <- "Swill is __ cost than conventional feed"
# no.qu.dfC$job_group <- factor(no.qu.dfC$job_group, 
#                               levels=c("Other","Pig farmer"),
#                               labels=c("Other","Pig\nfarmer"))
pCost <-  ggplot(no.qu.dfC, aes(job_group, n, fill=Q26_1)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  themeFig5_v2 +
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 + 
  facet_wrap(~title)


# tidy up
rm(no.qu.dfC)  


# Nutrition
# ... Q24_1
no.qu.dfN <- no.qu.df %>%       
  group_by(job_group) %>% 
  count(Q24_1)
no.qu.dfN
levels(no.qu.dfN$Q24_1) <-  c(levels(no.qu.dfN$Q24_1), "Don't know")     # Convert blanks in data.frame to "Don't know"
no.qu.dfN$Q24_1[no.qu.dfN$Q24_1==""] <- "Don't know"
unique(no.qu.dfN$Q24_1)
no.qu.dfN$Q24_1 <- factor(no.qu.dfN$Q24_1,
                          levels = c("Don't know", rev(c( "Much less nutritious", "Less nutritious", 
                                                          "Neither more or less nutritious", "More nutritious",
                                                          "Much more nutritious"))), 
                          labels = c("Don't know", rev(c("Much less\nnutritious", "Less nutritious", 
                                                         "Neither more nor\nless nutritious", "More nutritious",
                                                         "Much more\nnutritious"))))
no.qu.dfN$title <- "Swill is __ nutritious than conventional feed"
# no.qu.dfN$job_group <- factor(no.qu.dfN$job_group, 
#                               levels=c("Other","Pig farmer"),
#                               labels=c("Other","Pig\nfarmer"))

pNut <-  ggplot(no.qu.dfN, aes(job_group, n, fill=Q24_1)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  themeFig5_v2 + 
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  xlab("Swill is __ nutritious\nthan conventional feed") + 
  ylab("Number of respondents") + 
  scaleFig5 + 
  facet_wrap(~title)


# tidy up
rm(no.qu.dfN)  


# Variability of nutritiousness
no.qu.dfV <- no.qu.df %>% group_by(job_group) %>% 
  count(Q25_1)
no.qu.dfV
levels(no.qu.dfV$Q25_1) <-  c(levels(no.qu.dfV$Q25_1), "Don't know")
no.qu.dfV$Q25_1[no.qu.dfV$Q25_1==""] <- "Don't know"
unique(no.qu.dfV$Q25_1)
no.qu.dfV$Q25_1 <- factor(no.qu.dfV$Q25_1,
                          levels = c("Don't know", "Much less variable in nutritional content",
                                     "Less variable in nutritional content",
                                     "Neither more or less variable in nutritional content",
                                     "More variable in nutritional content",
                                     "Much more variable in nutritional content"), 
                          labels = c("Don't know", "Much less\nvariable",
                                     "Less variable",
                                     "Neither more nor\nless variable",
                                     "More variable",
                                     "Much more\nvariable"))
no.qu.dfV$title <- "Swill is __ variable"
# no.qu.dfV$job_group <- factor(no.qu.dfV$job_group, 
#                               levels=c("Other","Pig farmer"),
#                               labels=c("Other","Pig\nfarmer"))

pVar <-  ggplot(no.qu.dfV, aes(job_group, n, fill=Q25_1)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=14),
        axis.title.y = element_blank(), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(size=14), 
        strip.background = NULL) + 
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 +
  ylab("Number of respondents") + 
  facet_wrap(~title)
pVar


# tidy up
rm(no.qu.dfV)


# Disease risk
no.qu.dfD <- no.qu.df %>% group_by(job_group) %>% 
  count(Q27_1)
no.qu.dfD
levels(no.qu.dfD$Q27_1) <-  c(levels(no.qu.dfD$Q27_1), "Don't know")
no.qu.dfD$Q27_1[no.qu.dfD$Q27_1==""] <- "Don't know"
unique(no.qu.dfD$Q27_1)
no.qu.dfD$Q27_1 <- factor(no.qu.dfD$Q27_1,
                          levels = c("Don't know",
                                     "A much lower disease risk",
                                     "A lower disease risk",
                                     "Neither a higher nor lower disease risk",
                                     "A higher disease risk",
                                     "A much higher disease risk"), 
                          labels = c("Don't know",
                                     "A much lower\ndisease risk",
                                     "A lower\ndisease risk",
                                     "Neither a higher\nnor lower risk",
                                     "A higher\ndisease risk",
                                     "A much higher\ndisease risk"))
no.qu.dfD$title <- "Swill is __ disease risk"
pDiz <-  ggplot(no.qu.dfD, aes(job_group, n, fill=Q27_1)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=14),
        axis.title.y = element_blank(), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(size=14), 
        strip.background = NULL) + 
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  xlab("Swill is __ disease risk\nthan conventional feed") + 
  ylab("Number of respondents") + 
  scaleFig5 + 
  facet_wrap(~title)
pDiz


# tidy up
rm(no.qu.dfD)  


# Microbiological safety
no.qu.dfM <- no.qu.df %>% group_by(job_group) %>% 
  count(Q28_1)
no.qu.dfM
levels(no.qu.dfM$Q28_1) <-  c(levels(no.qu.dfM$Q28_1), "Don't know")
no.qu.dfM$Q28_1[no.qu.dfM$Q28_1==""] <- "Don't know"
unique(no.qu.dfM$Q28_1)
no.qu.dfM$Q28_1 <- factor(no.qu.dfM$Q28_1,
                          levels = c("Don't know",
                                     rev(c("Has much lower microbiological safety",
                                           "Has lower microbiological safety",
                                           "Neither higher nor lower microbiological safety",
                                           "Has higher microbiological safety",
                                           "Has much higher microbiological safety"))), 
                          labels = c("Don't know",rev(c(
                            "Much lower\nmicrobiological\nsafety",
                            "Lower\nmicrobiological\nsafety",
                            "Neither higher\nnor lower\nmicro. safety",
                            "Higher\nmicrobiological\nsafety",
                            "Much higher\nmicrobiological\nsafety"))))
no.qu.dfM$title <- "Swill has __ microbiological safety"
pMic <-  ggplot(no.qu.dfM, aes(job_group, n, fill=Q28_1)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  themeFig5_v2 +
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  xlab("Swill has __ microbiological safety\nthan conventional feed") + 
  ylab("Number of respondents") + 
  # scale_fill_brewer(palette="Spectral")
  scaleFig5 + 
  facet_wrap(~title)
pMic


# tidy up
rm(no.qu.dfM)  


# Chemical safety
no.qu.dfS <- no.qu.df %>% group_by(job_group) %>% 
  count(Q29_1)
no.qu.dfS
levels(no.qu.dfS$Q29_1) <-  c(levels(no.qu.dfS$Q29_1), "Don't know")
no.qu.dfS$Q29_1[no.qu.dfS$Q29_1==""] <- "Don't know"
unique(no.qu.dfS$Q29_1)
no.qu.dfS$Q29_1 <- factor(no.qu.dfS$Q29_1,
                          levels = c("Don't know",
                                     rev(c("Has much lower chemical safety",
                                           "Has lower chemical safety",
                                           "Neither higher nor lower chemical safety",
                                           "Has higher chemical safety",
                                           "Has much higher chemical safety"))), 
                          labels = c("Don't know",rev(c(
                            "Much lower\nchemical safety",
                            "Lower\nchemical safety",
                            "Neither higher\nnor lower\nchemical safety",
                            "Higher\nchemical safety",
                            "Much higher\nchemical safety"))))
no.qu.dfS$title <- "Swill has __ chemical safety"
pChm <-  ggplot(no.qu.dfS, aes(job_group, n, fill=Q29_1)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  themeFig5_v2 +
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 + 
  facet_wrap(~title)
pChm


# tidy up
rm(no.qu.dfS)


# Ethicalness
no.qu.dfT <- no.qu.df %>% group_by(job_group) %>% 
  count(Q30_1)
no.qu.dfT
levels(no.qu.dfT$Q30_1) <-  c(levels(no.qu.dfT$Q30_1), "Don't know")
no.qu.dfT$Q30_1[no.qu.dfT$Q30_1==""] <- "Don't know"
unique(no.qu.dfT$Q30_1)
no.qu.dfT$Q30_1 <- factor(no.qu.dfT$Q30_1,
                          levels = c("Don't know",
                                     rev(c("Much less ethical",
                                           "Less ethical",
                                           "Neither more or less ethical",
                                           "More ethical",
                                           "Much more ethical"))), 
                          labels = c("Don't know",rev(c(
                            "Much less\nethical",
                            "Less ethical",
                            "Neither more nor\nless ethical",
                            "More ethical",
                            "Much more\nethical"))))
no.qu.dfT$title <- "Swill is __ ethical"
pEth <-  ggplot(no.qu.dfT, aes(job_group, n, fill=Q30_1)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  themeFig5_v2 +
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  xlab("Swill is __ ethical\nthan conventional feed") + 
  ylab("Number of respondents") + 
  scaleFig5 + 
  facet_wrap(~title)
pEth


# tidy up
rm(no.qu.dfT)

# plot comparisons of impacts of swill and conventional feed
multiplot(pEnv, pCost)
multiplot(pEth, pNut)
multiplot(pChm, pMic)
multiplot(pVar, pDiz)


# tidy up
rm(pEnv, pCost, pEth, pNut, 
   pChm, pMic, pVar, pDiz)



#  ---------------------------- 8. Plot Appendix Fig 6. pig performance with swill or conventional feed ----------------------------


# Use the same scale and theme as for Figure 5


# Growth rates
no.qu.dfG <- no.qu.df %>% group_by(job_group) %>% 
  count(Q31_1)
no.qu.dfG
levels(no.qu.dfG$Q31_1) <-  c(levels(no.qu.dfG$Q31_1), "Don't know")
no.qu.dfG$Q31_1[no.qu.dfG$Q31_1==""] <- "Don't know"
unique(no.qu.dfG$Q31_1)
no.qu.dfG$Q31_1 <- factor(no.qu.dfG$Q31_1,
                          levels = c(rev(c("Much slower growth rates", "Slower growth rates", "Neither faster nor slower growth rates", 
                                           "Faster growth rates", "Much faster growth rates","Don't know"))), 
                          labels = c("Don't know",rev(c("Much slower\ngrowth rates", "Slower\ngrowth rates", "Neither faster\nnor slower", 
                                                        "Faster\ngrowth rates", "Much faster\ngrowth rates"))))
no.qu.dfG$title <- "Pigs fed swill have: __ growth rates"
pGro <-  ggplot(no.qu.dfG, aes(job_group, n, fill=Q31_1)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=14),   # theme includes x-axis label
        axis.title.y = element_blank(), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(size=14), 
        strip.background = NULL) +  
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  xlab("Pigs fed swill have:\n__ growth rates") + 
  ylab("Number of respondents") + 
  scaleFig5 + 
  facet_wrap(~title)


# tidy up
rm(no.qu.dfG)  


# FCR
no.qu.dfF <- no.qu.df %>% group_by(job_group) %>% 
  count(Q32_1)
no.qu.dfF
levels(no.qu.dfF$Q32_1) <-  c(levels(no.qu.dfF$Q32_1), "Don't know")
no.qu.dfF$Q32_1[no.qu.dfF$Q32_1==""] <- "Don't know"
unique(no.qu.dfF$Q32_1)
no.qu.dfF$Q32_1 <- factor(no.qu.dfF$Q32_1,
                          levels = c("Don't know","Much lower feed conversion ratios (more efficient)", "Lower feed conversion ratios (more efficient)", "Has neither higher nor lower feed conversion ratios", 
                                     "Higher feed conversion ratios (less efficient)", "Much higher feed conversion ratios (less efficient)"), 
                          labels = c("Don't know", "Much lower\nFCR", "Lower FCR", "Neither higher nor\nlower FCR", 
                                     "Higher FCR", "Much higher\nFCR"))
no.qu.dfF$title <- "Pigs fed swill have: __ feed conversion ratios"
pFCE <-  ggplot(no.qu.dfF, aes(job_group, n, fill=Q32_1)) + geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=14),   # theme includes x-axis label
        axis.title.y = element_blank(), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(size=14), 
        strip.background = NULL) + 
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  xlab("Pigs fed swill have:\n__ feed conversion ratios") + 
  ylab("Number of respondents") + 
  # scale_fill_brewer(palette="Spectral")
  scaleFig5 + 
  facet_wrap(~title)


# tidy up
rm(no.qu.dfF)  


# Welfare
no.qu.dfW <- no.qu.df %>% group_by(job_group) %>% 
  count(Q33_1)
no.qu.dfW
levels(no.qu.dfW$Q33_1) <-  c(levels(no.qu.dfW$Q33_1), "Don't know")
no.qu.dfW$Q33_1[no.qu.dfW$Q33_1==""] <- "Don't know"
unique(no.qu.dfW$Q33_1)
no.qu.dfW$Q33_1 <- factor(no.qu.dfW$Q33_1,
                          levels = c("Don't know",rev(c( "Much lower welfare", "Lower welfare", "Neither higher nor lower welfare", 
                                                         "Higher welfare", "Much higher welfare"))), 
                          labels = c("Don't know",rev(c( "Much lower\nwelfare", "Lower welfare", "Neither higher\nnor lower\nwelfare", 
                                                         "Higher welfare", "Much higher\nwelfare"))))
no.qu.dfW$title <- "Pigs fed swill have: __ welfare"
pWel <-  ggplot(no.qu.dfW, aes(job_group, n, fill=Q33_1)) + geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_blank(),   # theme includes x-axis label
        axis.title.y = element_blank(),
        axis.text = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(size=14), 
        strip.background = NULL) + 
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  xlab("Pigs fed swill have:\n__ welfare") + 
  ylab("Number of respondents") + 
  scaleFig5 + 
  facet_wrap(~title)
  
  
# tidy up
rm(no.qu.dfW)  


# feed costs
no.qu.dfFC <- no.qu.df %>% group_by(job_group) %>% 
  count(Q34_1)
no.qu.dfFC
levels(no.qu.dfFC$Q34_1) <-  c(levels(no.qu.dfFC$Q34_1), "Don't know", "Much higher feed costs")
no.qu.dfFC$Q34_1[no.qu.dfFC$Q34_1==""] <- "Don't know"
unique(no.qu.dfFC$Q34_1)
dim(no.qu.dfFC)
no.qu.dfFC[11,] <- c("Other", "Much higher feed costs", 0)
no.qu.dfFC[12,] <- c("Pig farmer", "Much higher feed costs", 0)
no.qu.dfFC$Q34_1 <- factor(no.qu.dfFC$Q34_1,
                          levels = c("Don't know", "Much lower feed costs", "Lower feed costs", "Neither higher nor lower feed costs", 
                                     "Higher feed costs","Much higher feed costs"), 
                          labels = c("Don't know", "Much lower\nfeed costs", "Lower\nfeed costs", "Neither higher\nnor lower\nfeed costs", 
                                     "Higher\nfeed costs", "Much higher\nfeed costs"))
no.qu.dfFC$n <- as.numeric(as.character(no.qu.dfFC$n))
no.qu.dfFC$title <- "Pigs fed swill have: __ feed costs"
pCst <-  ggplot(no.qu.dfFC, aes(job_group, n, fill=Q34_1)) + geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_blank(),   # theme includes x-axis label
        axis.title.y = element_blank(),
        axis.text = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(size=14), 
        strip.background = NULL)  + 
  ylim(0,82) +  # set ylim to 82, the number of pig farmers in the sample
  xlab("Pigs fed swill have:\n__ feed costs") + 
  ylab("Number of respondents") + 
  scaleFig5 + 
  facet_wrap(~title)


# tidy up
rm(no.qu.dfFC)


# plot four categories of impacts of will-fed pigs
multiplot(pCst, pWel)
multiplot(pGro, pFCE)


# tidy up
rm(pCst, pGro, pWel , pFCE)


# --------------------------------- 9. Plot S2 Appendix Fig. 6, pork produced on diets of swill or conventional feed ---------------------------------


# Use the same scale and theme as for Figure 5
theme_pork <- 
  theme(legend.title=element_blank(),
      legend.position="bottom",
      axis.title.x=element_blank(),   # theme includes x-axis label
      axis.title.y = element_blank(),
      axis.text = element_text(size=12), 
      legend.text = element_text(size=10),
      strip.text = element_text(size=14), 
      strip.background = NULL)


# Smelling
no.qu.df3 <- no.qu.df %>% group_by(job_group) %>% 
  count(Q43_2)
no.qu.df3
levels(no.qu.df3$Q43_2) <-  c(levels(no.qu.df3$Q43_2), "Don't know")
no.qu.df3$Q43_2[no.qu.df3$Q43_2==""] <- "Don't know"
unique(no.qu.df3$Q43_2)
dim(no.qu.df3)
no.qu.df3$Q43_2 <- factor(no.qu.df3$Q43_2,
                          levels = rev(c("Much worse smelling", "Worse smelling", "Neither better nor worse smelling", 
                                         "Better smelling","Much better smelling", "Don't know")), 
                          labels = rev(c("Much worse\nsmelling", "Worse\nsmelling", "Neither better\nnor worse", 
                                         "Better\nsmelling","Much better\nsmelling", "Don't know")))
no.qu.df3$n <- as.numeric(as.character(no.qu.df3$n))
no.qu.df3$title <- "Pork from swill-fed pigs is: __ smelling"
pPSme <-  ggplot(no.qu.df3, aes(job_group, n, fill=Q43_2)) +
  geom_bar(stat="identity") + coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=14),   # theme includes x-axis label
        axis.title.y = element_blank(), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(size=14), 
        strip.background = NULL) + 
  xlab("Pork from swill-fed pigs\n is: __ smelling") + 
  ylab("Number of respondents") + 
  ylim(0,82) +   # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 + 
  facet_wrap(~title)


# tidy up
rm(no.qu.df3)


# Fattiness
no.qu.df3 <- no.qu.df %>% group_by(job_group) %>% 
  count(Q40_1)
no.qu.df3
levels(no.qu.df3$Q40_1) <-  c(levels(no.qu.df3$Q40_1), "Don't know")
no.qu.df3$Q40_1[no.qu.df3$Q40_1==""] <- "Don't know"
unique(no.qu.df3$Q40_1)
dim(no.qu.df3)
no.qu.df3$Q40_1 <- factor(no.qu.df3$Q40_1,
                          levels = c("Don't know","Much less fatty", "Less fatty", "Neither more nor less fatty", 
                                     "More fatty","Much more fatty"), 
                          labels = c("Don't know", "Much less\nfatty", "Less fatty", "Neither more\nnor less fatty", 
                                     "More fatty","Much more\nfatty"))
no.qu.df3$title <- "Pork from swill-fed pigs is: __ fatty"
no.qu.df3$n <- as.numeric(as.character(no.qu.df3$n))
pPfat <-  ggplot(no.qu.df3, aes(job_group, n, fill=Q40_1)) +
  geom_bar(stat="identity") + coord_flip() + 
  theme_bw() + 
  theme_pork + 
  xlab("Pork from swill-fed pigs\n is: __ fatty") + 
  # ylab("Number of respondents") + 
  ylim(0,82) +   # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 +
  facet_wrap(~title)


# tidy up
rm(no.qu.df3)


# Colour
no.qu.df3 <- no.qu.df %>% group_by(job_group) %>% 
  count(Q41_1)
no.qu.df3
levels(no.qu.df3$Q41_1) <-  c(levels(no.qu.df3$Q41_1), "Don't know", "Much darker in colour", "Much lighter in colour")
no.qu.df3$Q41_1[no.qu.df3$Q41_1==""] <- "Don't know"
unique(no.qu.df3$Q41_1)
dim(no.qu.df3)
no.qu.df3[10,] <- c("Other", "Much darker in colour", 0)        # Add in rows where there were no responses (for plotting purposes)
no.qu.df3[11,] <- c("Pig farmer", "Much darker in colour", 0)   # (for plotting purposes)
no.qu.df3[10,] <- c("Pig farmer", "Much lighter in colour", 0)  # (for plotting purposes)
no.qu.df3$Q41_1 <- factor(no.qu.df3$Q41_1,
                          levels = rev(c("Much darker in colour","Darker in colour", "Neither lighter nor darker in colour", # NB no respondents said it was "Much darker"
                                         "Lighter in colour","Much lighter in colour", "Don't know")), 
                          labels = rev(c("Much darker","Darker", "Neither lighter\nnor darker", # NB no respondents said it was "Much darker"
                                         "Lighter","Much lighter", "Don't know")))
no.qu.df3$n <- as.numeric(as.character(no.qu.df3$n))
no.qu.df3$title <- "Pork from swill-fed pigs is: __ in colour"
pPcol <-  ggplot(no.qu.df3, aes(job_group, n, fill=Q41_1)) +
  geom_bar(stat="identity") + coord_flip() + 
  theme_bw() + 
  theme_pork + 
  theme(axis.title.x=element_blank()) + # to remove the y-axis label
  xlab("Pork from swill-fed pigs\n is: __ in colour") + 
  # ylab("Number of respondents") + 
  ylim(0,82) +   # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 +
  facet_wrap(~title)


# tidy up
rm(no.qu.df3)


# Tastiness
no.qu.df3 <- no.qu.df %>% group_by(job_group) %>% 
  count(Q42_1)
no.qu.df3
levels(no.qu.df3$Q42_1) <-  c(levels(no.qu.df3$Q42_1), "Don't know")
no.qu.df3$Q42_1[no.qu.df3$Q42_1==""] <- "Don't know"
unique(no.qu.df3$Q42_1)
dim(no.qu.df3)
no.qu.df3$Q42_1 <- factor(no.qu.df3$Q42_1,
                          levels = rev(c("Much less tasty", "Less tasty", "Neither more nor less tasty",
                                         "More tasty","Much more tasty", "Don't know")), 
                          labels = rev(c("Much less\ntasty", "Less tasty", "Neither more\nnor less tasty",
                                         "More tasty","Much more\ntasty", "Don't know")))
no.qu.df3$n <- as.numeric(as.character(no.qu.df3$n))
no.qu.df3$title <- "Pork from swill-fed pigs is: __ tasty"
pPtst <-  ggplot(no.qu.df3, aes(job_group, n, fill=Q42_1)) +
  geom_bar(stat="identity") + coord_flip() +
  theme_bw() + 
  theme_pork + 
  theme(axis.title.x=element_blank()) + # to remove the y-axis label
  # xlab("Pork from swill-fed pigs\n is: __ tasty") + 
  # ylab("Number of respondents") + 
  ylim(0,82) +   # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 +
  facet_wrap(~title)


# tidy up
rm(no.qu.df3)


# Marketability
no.qu.df3 <- no.qu.df %>% group_by(job_group) %>% 
  count(Q44_1)
no.qu.df3
levels(no.qu.df3$Q44_1) <-  c(levels(no.qu.df3$Q44_1), "Don't know")
no.qu.df3$Q44_1[no.qu.df3$Q44_1==""] <- "Don't know"
unique(no.qu.df3$Q44_1)
dim(no.qu.df3)
no.qu.df3$Q44_1 <- factor(no.qu.df3$Q44_1,
                          levels=rev(c("Much less marketable", "Less marketable", "Neither more nor less marketable",
                                       "More marketable","Much more marketable", "Don't know")), 
                          labels=rev(c("Much less\nmarketable", "Less\nmarketable", "Neither more\nnor less",
                                       "More\nmarketable","Much more\nmarketable", "Don't know")))

no.qu.df3$n <- as.numeric(as.character(no.qu.df3$n))
no.qu.df3$title <- "Pork from swill-fed pigs is: __ marketable"
pPMrk <-  ggplot(no.qu.df3, aes(job_group, n, fill=Q44_1)) +
  geom_bar(stat="identity") + coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=14),   # theme includes x-axis label
        axis.title.y = element_blank(), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(size=14), 
        strip.background = NULL) + 
  xlab("Pork from swill-fed pigs\n is: __ marketable") + 
  ylab("Number of respondents") + 
  ylim(0,82) +   # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 +
  facet_wrap(~title)

# tidy up
rm(no.qu.df3)


# Profitability
no.qu.df3 <- no.qu.df %>% group_by(job_group) %>% 
  count(Q46_1)
no.qu.df3
levels(no.qu.df3$Q46_1) <-  c(levels(no.qu.df3$Q46_1), "Don't know")
no.qu.df3$Q46_1[no.qu.df3$Q46_1==""] <- "Don't know"
unique(no.qu.df3$Q46_1)
dim(no.qu.df3)
no.qu.df3$Q46_1 <- factor(no.qu.df3$Q46_1,
                          levels=rev(c("Much less profitable", "Less profitable", "Neither more nor less profitable",
                                       "More profitable","Much more profitable", "Don't know")), 
                          labels=rev(c("Much less\nprofitable", "Less\nprofitable", "Neither more\nnor less",
                                       "More profitable","Much more\nprofitable", "Don't know")))

no.qu.df3$n <- as.numeric(as.character(no.qu.df3$n))
no.qu.df3$title <- "Pork from swill-fed pigs is: __ profitable"
pPPrf <-  ggplot(no.qu.df3, aes(job_group, n, fill=Q46_1)) +
  geom_bar(stat="identity") + coord_flip() + 
  theme_bw() + 
  theme_pork + 
  xlab("Pork from swill-fed pigs\n is: __ profitable") + 
  # ylab("Number of respondents") + 
  ylim(0,82) +   # set ylim to 82, the number of pig farmers in the sample
  scaleFig5 +
  facet_wrap(~title)


# tidy up
rm(no.qu.df3)


# multiplot
multiplot(pPPrf, pPfat)
multiplot(pPcol, pPtst)
multiplot(pPSme, pPMrk) # NB remove the second plot when joining images


# tidy up
rm(pPPrf, pPMrk, pPtst, 
   pPcol, pPfat, pPSme)
rm(theme_pork, themeFig5_v2, scaleFig5, my.cols)


#   --------------------------------------- 10. Plot Fig 3, support for relegalisation  ---------------------------------------


# Plot of support for relegalization and willingness to use swill in farm
Q48.cols <- c(rev(brewer.pal(5, "Spectral"))) # make a custom colour palette
scale_Q48 <-   scale_fill_manual(values = Q48.cols,
                                 # direction = -1,
                                 guide = guide_legend(reverse=TRUE,
                                                      direction = "horizontal",
                                                      keyheight = unit(4, units = "mm"),
                                                      keywidth = unit(28 / length(labels), units = "mm"),
                                                      title.position = 'top',
                                                      nrow = 1,
                                                      byrow = T,
                                                      label.position = "bottom"
                                 )) 
cols_to_keep <- c("Q48", "Q51")     # keep data on support for relegalisation and willingness to use swill
unique(no.qu.df$Q48)
unique(no.qu.df$Q51)
levels(no.qu.df$Q51) <- c(levels(no.qu.df$Q51), "Definitely yes")
no.qu.df$Q51[which(no.qu.df$Q51 == "Absolutely")] <-  "Definitely yes"
qu.values <- qu.df[,cols_to_keep]                      # make a df of the relevant questions
qu.values2 <- as.data.frame(t(qu.values))          # Transpose so questions are in rows
qu.values2$question_num <- row.names(qu.values2)      # Make a new column with the question identifiers
colnames(qu.values2)[1] <- "question"                  # give the first col a name
qu.values2$question_short <- sub(".*), ", "", qu.values2$question) # Make a new column with the name of the feedstuff
head(qu.values2)                                                       # Look at the data
qu.values3 <- qu.values2[,c("question_short", "question_num")]        # drop extra column
no.qu.values <- gather(no.qu.df[,c(cols_to_keep, "job_group")], question_num, value, Q48:Q51) # convert wide to long
no.qu.values2 <- join(no.qu.values, qu.values2, by = "question_num", type="inner")
dim(no.qu.values2); dim(no.qu.values)     # check no rows were lost
no.qu.values2$value[which(no.qu.values2$value == "")] <- NA
no.qu.values2$value <- factor(no.qu.values2$value, 
                              levels = rev(c("Definitely not", "Probably not", "Might or might not", "Probably yes", "Definitely yes")),
                              labels= rev(c("Definitely not", "Probably not", "Might or\nmight not", "Probably yes", "Definitely yes")))
no.qu.values2$question_short <- factor(no.qu.values2$question_short, 
                                       levels=c("would you support the re-legalisation of swill?",
                                                "would you consider using swill on your farm?"),
                                       labels=c("Would you support the relegalisation of swill?",
                                                "Would you consider using swill on your farm?"))


# Nb I remove manually the gap for "Other" willingess to use
# ... scales = "free" in facet_wrap would make the bars different sizes
Fig3 <- ggplot(na.omit(no.qu.values2[,c("job_group","value", "question_short" )]), aes(job_group, fill=value)) + geom_bar(stat="count") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=14),   # theme includes x-axis label
        axis.title.y = element_text(size=14), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=12),
        strip.text = element_text(size=14), 
        strip.background = NULL) +
  xlab("") +
  ylab("Number of respondents") + 
  facet_wrap(~question_short,
            ncol = 1) + 
  scale_Q48
Fig3


# tidy up
rm(scale_Q48, Q48.cols)
rm(cols_to_keep, qu.values, qu.values2, qu.values3,
   no.qu.values, no.qu.values2, Fig3) 


#  ---------------------------------- 11. S2 Appendix Fig 7, traditional/unnatural prractice ----------------------------------


# Most people think that swill is NOT unnatural
# ... and that it is a traditional farming practice
Q48.cols <- c(rev(brewer.pal(5, "Spectral"))) # make a custom colour palette
cols_to_keep <- c("Q53_1", "Q53_2")
qu.values <- qu.df[,cols_to_keep]                      # make a df of the relevant questions
qu.values2 <- as.data.frame(t(qu.values))          # Transpose so questions are in rows
qu.values2$question_num <- row.names(qu.values2)      # Make a new column with the question identifiers
colnames(qu.values2)[1] <- "question"                  # give the first col a name
qu.values2$question_short <- sub(".* - ", "", qu.values2$question) # Make a new column with the name of the feedstuff
head(qu.values2)                                                       # Look at the data
qu.values3 <- qu.values2[,c("question_short", "question_num")]        # drop extra column
no.qu.values <- gather(no.qu.df[,c(cols_to_keep, "job_group")], question_num, value, Q53_1:Q53_2) # convert wide to long
no.qu.values2 <- join(no.qu.values, qu.values2, by = "question_num", type="inner")
dim(no.qu.values2); dim(no.qu.values)     # check no rows were lost
head(no.qu.values2)
no.qu.values2$value <- factor(no.qu.values2$value, 
                              levels=rev(c("Definitely not", "Probably not", "Not sure", "Probably yes", "Definitely yes")))
S2AppFig7 <- ggplot(no.qu.values2, aes(job_group, fill=value)) + geom_bar(stat="count") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.title.x=element_text(size=14),   # theme includes x-axis label
        axis.title.y = element_text(size=14), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=12),
        strip.text = element_text(size=14), 
        strip.background = NULL) +
  xlab("") +
  ylab("Number of respondents") + 
  facet_wrap(~question_short, 
             ncol=1) + 
  scale_fill_manual(values = Q48.cols,
                 # direction = -1,
                 guide = guide_legend(reverse=TRUE,
                                      direction = "horizontal",
                                      keyheight = unit(4, units = "mm"),
                                      keywidth = unit(28 / length(labels), units = "mm"),
                                      title.position = 'top',
                                      nrow = 1,
                                      byrow = T,
                                      label.position = "bottom"
                 )) 
S2AppFig7


# tidy up
rm(cols_to_keep, qu.values, qu.values2, qu.values3,
   no.qu.values, no.qu.values2, S2AppFig7, 
   Q48.cols) 



#  ------------------------------------ 12. Factor analysis of respondent values and perceptions of swill ------------------------------------


# i) Tidy up data.frame for Polychoric correlation FA of:
# ii) Values (How important are the following...)
# iii) Perceived impacts of swill


# select out the questions that are important
# ... then convert all responses to numeric
names(no.qu.df)
PCA_cols <- c("Q4_1", "Q4_2", "Q4_4", "Q4_5","Q4_6", "Q4_7", "Q4_8", "Q4_9", "Q4_10","Q4_11", "Q4_12", "Q4_13",   # questions about impact of swill
              "Q54_1","Q54_2","Q54_3","Q54_4","Q54_5","Q54_6","Q54_7","Q54_8","Q54_9","Q54_10","Q54_11","Q54_12") # questions about what is important to farmers
length(PCA_cols) # there are 43 questions included in the PCA
PCA.df <- no.qu.df[,PCA_cols]              
rm(PCA_cols)         # tidy up


# 12i) Plot S2 Appendix Fig 8, respondent values ----


# convert wide to long
qu_values <- c("Q54_1", "Q54_2", "Q54_3", "Q54_4",      # list of relevant questions. 
                  "Q54_5", "Q54_6", "Q54_7", "Q54_8", "Q54_9",
                  "Q54_10", "Q54_11", "Q54_12")
qu.priorities <- qu.df[,qu_values]                   # make a df of the relevant questions
dim(qu.priorities)                                      # should be 12 cols
qu.priorities2 <- as.data.frame(t(qu.priorities))          # Transpose so questions are in rows
qu.priorities2$question_id <- row.names(qu.priorities2)      # Make a new column with the question identifiers
colnames(qu.priorities2)[1] <- "question"                  # give the first col a name
qu.priorities2$impact <- sub(".* - ", "", qu.priorities2$question) # Make a new column with the name of the feedstuff
head(qu.priorities2)                                                       # Look at the data
qu.priorities2 <- qu.priorities2[,c("question_id", "impact")]        # drop extra column
no.qu.priorities <- gather(no.qu.df[,c(qu_values,"job_group")], question_id, acceptability, Q54_1:Q54_12) # convert wide to long
no.qu.priorities2 <- join(no.qu.priorities, qu.priorities2, by = "question_id", type="inner")
dim(no.qu.priorities2); dim(no.qu.priorities)     # check no rows were lost
head(no.qu.priorities2)
no.qu.priorities2$acceptability <- factor(no.qu.priorities2$acceptability, 
                                          levels=rev(c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important")),
                                          labels=rev(c("Not at all\nimportant", "Not important", "Neither important\nnor unimportant", "Important", "Very Important")))
no.qu.priorities2$job_group <- factor(no.qu.priorities2$job_group,
                                      levels=c("Pig farmer", "Other"))
require(RColorBrewer)
Fig15cols <- rev(brewer.pal(5, "Spectral"))
Fig15 <- ggplot(no.qu.priorities2, aes(impact, fill=acceptability)) + 
  geom_bar(stat="count") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        # axis.title.x=element_text(size=14),
        axis.title.y=element_blank(),
        axis.title.x = element_text(size=14), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=12),
        strip.text = element_text(size=14),
        strip.background = NULL) +
  scale_fill_manual(values = Fig15cols,
                    # direction = -1,
                    guide = guide_legend(reverse=T,
                                         direction = "horizontal",
                                         keyheight = unit(4, units = "mm"),
                                         keywidth = unit(36 / length(labels), units = "mm"),
                                         title.position = 'top',
                                         nrow = 1,
                                         byrow = T,
                                         label.position = "bottom"
                    )) +
  ylab("Number of respondents") + 
  facet_wrap(~job_group) 
Fig15

# tidy up
rm(qu_values, qu.priorities, qu.priorities2, 
   no.qu.priorities, no.qu.priorities2, Fig15, Fig15cols)


# 12ii) Plot S2 Appendix Fig 9, Perceived impacts of swill ----


# convert wide to long
impacts_qus <- c("Q4_1", "Q4_2", "Q4_4", "Q4_5",          # list of questions. Nb there is no Q_3
                  "Q4_6", "Q4_7", "Q4_8", "Q4_9", "Q4_10",
                  "Q4_11", "Q4_12", "Q4_13")
qu.swill.impacts <- qu.df[,impacts_qus]                   # make a df of the relevant questions
dim(qu.swill.impacts)                                      # should be 12 cols
qu.swill.impacts2 <- as.data.frame(t(qu.swill.impacts))          # Transpose so questions are in rows
qu.swill.impacts2$question_id <- row.names(qu.swill.impacts2)      # Make a new column with the question identifiers
colnames(qu.swill.impacts2)[1] <- "question"                  # give the first col a name
qu.swill.impacts2$impact <- sub(".* - ", "", qu.swill.impacts2$question) # Make a new column with the name of the feedstuff
head(qu.swill.impacts2)                                                       # Look at the data
qu.swill.impacts2 <- qu.swill.impacts2[,c("question_id", "impact")]        # drop extra column
no.qu.swill.impacts <- gather(no.qu.df[,c(impacts_qus,"job_group")], question_id, acceptability, Q4_1:Q4_13) # convert wide to long
no.qu.swill.impacts2 <- join(no.qu.swill.impacts, qu.swill.impacts2, by = "question_id", type="inner")
dim(no.qu.swill.impacts2); dim(no.qu.swill.impacts)     # check no rows were lost
head(no.qu.swill.impacts2)
no.qu.swill.impacts2$impact <- factor(no.qu.swill.impacts2$impact, levels=c("Lower dependence on foreign protein sources",                                                       
                                                                            "Reduce the environmental impact of food waste disposal",                                           
                                                                            "Reduce the environmental impact of pork production",                                                
                                                                            "Help farms reduce feed costs",                                                                      
                                                                            "Help farmers improve profitability",                                                                
                                                                            "Lower consumer acceptance of pork products",                                                        
                                                                            "Increase the risk of an outbreak of foot-and-mouth disease",                                        
                                                                            "Increase the risk of prion diseases like BSE (mad cow disease) or vCJD (Creutzfeldt-Jacob disease)",
                                                                            "Increase the risk of toxins entering the feed",                                                     
                                                                            "Reduce the traceability of feed production",                                                        
                                                                            "Be an efficient way to use food waste",                                                             
                                                                            "Negatively affect the marketability of pork"),
                                      labels = c("Lower dependence on foreign protein sources",                                                       
                                                 "Reduce the environmental impact of food waste disposal",                                           
                                                 "Reduce the environmental impact of pork production",                                                
                                                 "Help farms reduce feed costs",                                                                      
                                                 "Help farmers improve profitability",                                                                
                                                 "Lower consumer acceptance of pork products",                                                        
                                                 "Increase the risk of an outbreak of foot-and-mouth disease",                                        
                                                 "Increase the risk of prion diseases like BSE (mad cow disease)\nor vCJD (Creutzfeldt-Jacob disease)",
                                                 "Increase the risk of toxins entering the feed",                                                     
                                                 "Reduce the traceability of feed production",                                                        
                                                 "Be an efficient way to use food waste",                                                             
                                                 "Negatively affect the marketability of pork"))
no.qu.swill.impacts2$job_group <- factor(no.qu.swill.impacts2$job_group, 
                                         levels = c("Pig farmer", "Other"))
levels(no.qu.swill.impacts2$acceptability) <- c(levels(no.qu.swill.impacts2$acceptability), "Don't know")   # Don't know responses are recorded as NAs, here I correc this for the purposes of plotting
no.qu.swill.impacts2$acceptability[which(no.qu.swill.impacts2$acceptability == "")] <- "Don't know"         # Don't know responses are recorded as NAs, here I correc this for the purposes of plotting
no.qu.swill.impacts2$acceptability <- factor(no.qu.swill.impacts2$acceptability, 
                                             levels=rev(c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree","Don't know")),
                                             labels=rev(c("Totally\ndisagree", "Disagree", "Neither agree\nnor disagree", "Agree", "Totally\nagree","Don't know")))
my.cols <- c("gray", rev(brewer.pal(5, "Spectral"))) # make a custom colour palette
S2AppFig9 <- ggplot(no.qu.swill.impacts2, aes(impact, fill=acceptability)) + geom_bar(stat="count") + 
  coord_flip() + 
  theme_bw() + 
  xlab("") + 
  facet_wrap(~job_group) + 
  theme(legend.title=element_blank(),
        legend.position="bottom",
        # axis.title.x=element_text(size=14),
        axis.title.x=element_text(size=14),
        axis.title.y = element_blank(), 
        axis.text = element_text(size=12), 
        legend.text = element_text(size=12),
        strip.text = element_text(size=14),
        strip.background = NULL) +
  ylab("Number of respondents") + 
  scale_fill_manual(values = my.cols,
                               # direction = -1,
                               guide = guide_legend(reverse=TRUE,
                                                    direction = "horizontal",
                                                    keyheight = unit(4, units = "mm"),
                                                    keywidth = unit(28 / length(labels), units = "mm"),
                                                    title.position = 'top',
                                                    nrow = 1,
                                                    byrow = T,
                                                    label.position = "bottom"
                               )) 
S2AppFig9


# tidy up
rm(impacts_qus, qu.swill.impacts, qu.swill.impacts2, S2AppFig9, 
   no.qu.swill.impacts, no.qu.swill.impacts2)
rm(my.cols)


# 12iii) Convert all values into a numeric scale ---- 
# ... first re-classify each "Don't know" response to the median value (i.e. neither more nor less)
# Nb that of the 163 respondents, only 35 didn't include any NAs in their answers
# ... to the 43 questions included in the PCA
PCA.df2 <- PCA.df            # make a copy
PCA.df2[PCA.df2==""]<-NA     # reclassify "don't know"/"unsure" to NA
nrow(PCA.df2[complete.cases(PCA.df2),]) # Count the number of respondents without any "Don't knows" in their data: 97/163
# Then convert all values into a numeric scale.
unique(PCA.df2$Q4_1)        # check the factor levels of Q4_1
PCA.df2$Q4_1 <- factor(PCA.df2$Q4_1, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                       labels=c(1:5))
unique(PCA.df2$Q4_2)        # check the factor levels of Q4_2
PCA.df2$Q4_2 <- factor(PCA.df2$Q4_2, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                       labels=c(1:5))
unique(PCA.df2$Q4_4)        # check the factor levels of Q4_4
PCA.df2$Q4_4 <- factor(PCA.df2$Q4_4, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                       labels=c(1:5))
unique(PCA.df2$Q4_5)        # check the factor levels of Q4_5
PCA.df2$Q4_5 <- factor(PCA.df2$Q4_5, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                       labels=c(1:5))
unique(PCA.df2$Q4_6)        # check the factor levels of Q4_6
PCA.df2$Q4_6 <- factor(PCA.df2$Q4_6, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                       labels=c(1:5))
unique(PCA.df2$Q4_7)        # check the factor levels of Q4_7
PCA.df2$Q4_7 <- factor(PCA.df2$Q4_7, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                       labels=c(1:5))
unique(PCA.df2$Q4_8)        # check the factor levels of Q4_8
PCA.df2$Q4_8 <- factor(PCA.df2$Q4_8, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                       labels=c(1:5))
unique(PCA.df2$Q4_9)        # check the factor levels of Q4_9
PCA.df2$Q4_9 <- factor(PCA.df2$Q4_9, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                       labels=c(1:5))
unique(PCA.df2$Q4_10)        # check the factor levels of Q4_10
PCA.df2$Q4_10 <- factor(PCA.df2$Q4_10, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                        labels=c(1:5))
unique(PCA.df2$Q4_11)        # check the factor levels of Q4_11
PCA.df2$Q4_11 <- factor(PCA.df2$Q4_11, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                        labels=c(1:5))
unique(PCA.df2$Q4_12)        # check the factor levels of Q4_12
PCA.df2$Q4_12 <- factor(PCA.df2$Q4_12, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                        labels=c(1:5))
unique(PCA.df2$Q4_13)        # check the factor levels of Q4_13
PCA.df2$Q4_13 <- factor(PCA.df2$Q4_13, levels=c("Totally disagree", "Disagree", "Neither agree nor disagree", "Agree", "Totally agree"),
                        labels=c(1:5))
unique(PCA.df2$Q54_1)        # check the factor levels of Q54_1
PCA.df2$Q54_1 <- factor(PCA.df2$Q54_1, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                        labels=c(1:5))
unique(PCA.df2$Q54_2)        # check the factor levels of Q54_2
PCA.df2$Q54_2 <- factor(PCA.df2$Q54_2, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                        labels=c(1:5))
unique(PCA.df2$Q54_3)        # check the factor levels of Q54_3
PCA.df2$Q54_3 <- factor(PCA.df2$Q54_3, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                        labels=c(1:5))
unique(PCA.df2$Q54_4)        # check the factor levels of Q54_4
PCA.df2$Q54_4 <- factor(PCA.df2$Q54_4, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                        labels=c(1:5))
unique(PCA.df2$Q54_5)        # check the factor levels of Q54_5
PCA.df2$Q54_5 <- factor(PCA.df2$Q54_5, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                        labels=c(1:5))
unique(PCA.df2$Q54_6)        # check the factor levels of Q54_6
PCA.df2$Q54_6 <- factor(PCA.df2$Q54_6, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                        labels=c(1:5))
unique(PCA.df2$Q54_7)        # check the factor levels of Q54_7
PCA.df2$Q54_7 <- factor(PCA.df2$Q54_7, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                        labels=c(1:5))
unique(PCA.df2$Q54_8)        # check the factor levels of Q54_8
PCA.df2$Q54_8 <- factor(PCA.df2$Q54_8, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                        labels=c(1:5))
unique(PCA.df2$Q54_9)        # check the factor levels of Q54_9
PCA.df2$Q54_9 <- factor(PCA.df2$Q54_9, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                        labels=c(1:5))
unique(PCA.df2$Q54_10)        # check the factor levels of Q54_10
PCA.df2$Q54_10 <- factor(PCA.df2$Q54_10, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                         labels=c(1:5))
unique(PCA.df2$Q54_11)        # check the factor levels of Q54_11
PCA.df2$Q54_11 <- factor(PCA.df2$Q54_11, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                         labels=c(1:5))
unique(PCA.df2$Q54_12)        # check the factor levels of Q54_12
PCA.df2$Q54_12 <- factor(PCA.df2$Q54_12, levels=c("Not at all important", "Not important", "Neither important nor unimportant", "Important", "Very Important"),
                         labels=c(1:5))


# Look at the df to check that all values are now numbers
head(PCA.df2)               
class(PCA.df2$Q54_12) # these columns are still classified as factors
PCA.df3 <- sapply(PCA.df2, function(x) as.numeric(as.character(x)) ) # convert all values to numeric
PCA.df3 <- as.data.frame(PCA.df3) # convert back to a data.frame (not matrix)
head(PCA.df3)        # compare these with PCA.df2; they are the same.
class(PCA.df3$Q4_1)  # these columns are now numbers
rm(PCA.df2,PCA.df)   # tidy up


# 12iv) Values (How important are the following...) ---- 


# Select the correct 12 columns
qu.df[,c("Q54_1","Q54_2","Q54_3","Q54_4","Q54_5",                # Look at the questions - check I've selected the correct 12
         "Q54_6","Q54_7","Q54_8","Q54_9","Q54_10",
         "Q54_11","Q54_12")]
values_cols <- PCA.df3[, c("Q54_1","Q54_2","Q54_3","Q54_4",      # Make new df with only these columns
                           "Q54_5","Q54_6","Q54_7","Q54_8",
                           "Q54_9","Q54_10","Q54_11","Q54_12")]
ncol(values_cols)                                                # Should be 12 columns
head(values_cols)


# Check for NAs
# ... replace with the median value for factor analysis
nrow(values_cols) - nrow(values_cols[complete.cases(values_cols),]) # One row has an NA
median(values_cols[complete.cases(values_cols),"Q54_3"])            # the median is 5
values_cols[160,"Q54_3"]                                            # Identify the NA value
values_cols[160,"Q54_3"] <- median(values_cols[complete.cases(values_cols),"Q54_3"]) # ... set it to the median


# Load the "psych" package 
library(psych)


# The Velicer MAP achieves a minimum of 0.03  with  2  factors
# ... NB that an ultra-Heywood case was detected when using VSS to select the number of factors to extract
# ... this occurs where there is communality >1 (https://stats.stackexchange.com/questions/216103/do-heywood-cases-render-efa-cfa-solutions-invalid)
myVSS <- VSS(values_cols, rotate ="varimax",
             title="VSS of farmer values\nrotation=varimax")   
myVSS
rm(myVSS)


# factor analysis with polychloric correlation
# ... Nb I get the warning that:
# 66 cells were adjusted for 0 values using the correction for continuity. Examine your data carefully.
# >>> this is because cells with zero counts are replaced with .5 as a correction for continuity (correct=TRUE).
# ... 2 factors explain 45% of variance.
myFA <- fa(values_cols, cor='poly', 
           nfactors=2, 
           rotate="varimax")
myFA$Structure      # 2 factors explain 45% of variance.
myFA$scores         # these are the principal component scores per respondent
myFA$rotation       # this used the varimax rotation
pairs(myFA$scores,panel=panel.smooth) # check the components are not correlated


# Interpreting the principal components of the varimax rotated data


# PC1 strong loadings:
# >>> PC1: Concerned about perception of the pork industry by consumers and disease control:
# ... strong loadings for Communication with consumers, Traceability, Perception of the pork industry, Labelling of the end product, Food safety & Disease control
PC1LOADINGS <- myFA$loadings[,1]
class(PC1LOADINGS)
PC1LOADINGS[PC1LOADINGS>0.5 | PC1LOADINGS < -0.5] # questions with strong loadings
qu.df[c("Q54_1", "Q54_2", "Q54_3","Q54_4","Q54_5", "Q54_6", "Q54_8")]
rm(PC1LOADINGS)    # tidy up


# PC2 strong loadings:
# >>> PC2: Concerned about financial performance and efficiency
# ... strong loadings for Feed prices, Profitability, Efficient use of resources
PC2LOADINGS <- myFA$loadings[,2]
class(PC2LOADINGS)
PC2LOADINGS[PC2LOADINGS>0.5 | PC2LOADINGS < -0.5] # questions with strong loadings
qu.df[c("Q54_7", "Q54_11", "Q54_12")]
rm(PC2LOADINGS)    # tidy up


# Add the first two factor levels to the dataset
# ... ultimately used in a regression model to understand the support for legal change
no.qu.df$PC1_values <- myFA$scores[,1]
no.qu.df$PC2_values <- myFA$scores[,2]


# tidy up
rm(myFA, values_cols)


# 12v) impacts of swill ----


# Select the correct 12 columns
impact_col_names <- c("Q4_1", "Q4_2", "Q4_4", "Q4_5",    # questions about impact of swill
                      "Q4_6", "Q4_7", "Q4_8", "Q4_9", 
                      "Q4_10","Q4_11", "Q4_12", "Q4_13")
qu.df[,impact_col_names] 
impact_cols <- PCA.df3[,impact_col_names]
head(impact_cols)


# Check for NAs
# ... there are 66 rows with NAs (i.e. Don't know answers)
# ... replace with the median value for factor analysis
nrow(impact_cols) - nrow(impact_cols[complete.cases(impact_cols),])   # 66 rows have an NA
for(i in 1:length(impact_col_names)){     # For loop to replace all NAs with the median value for that question
  impact_cols[is.na(impact_cols[,impact_col_names[i] ]),impact_col_names[i]] <- median(impact_cols[complete.cases(impact_cols),impact_col_names[i]])
}
rm(i, impact_col_names)   # tidy up
nrow(impact_cols) - nrow(impact_cols[complete.cases(impact_cols),])   # Zero rows have an NA


# Look at MAP estimate of number of factors to extract
# ... The Velicer MAP achieves a minimum of 0.03  with  2  factors
myVSS_impacts <- VSS(impact_cols, rotate ="varimax",
                     title="VSS of perceptions of impact of swill\nrotation=varimax")   
myVSS_impacts
rm(myVSS_impacts)


# factor analysis with polychloric correlation
# ... Nb I get the warning that:
# 66 cells were adjusted for 0 values using the correction for continuity. Examine your data carefully.
# >>> this is because cells with zero counts are replaced with .5 as a correction for continuity (correct=TRUE).
# ... 2 factors explain 43% of variance.
myFA_impacts <- fa(impact_cols, cor='poly', 
                   nfactors=2, 
                   rotate="varimax")
myFA_impacts$Structure      # 2 factors explain 43% of variance.
myFA_impacts$scores         # these are the principal component scores per respondent
myFA_impacts$rotation       # this used the varimax rotation
pairs(myFA_impacts$scores,panel=panel.smooth) # check the components are not correlated


# Interpreting the principal components of the varimax rotated data


# PC1 strong loadings:
# >>> PC1: Think that swill would be good for the environment, help farms financially, and reduce trade-deficit
# ... strong loadings for Reduce the environmental impact of food waste disposal, Reduce the environmental impact of pork production, Be an efficient way to use food waste, Help farms reduce feed costs, Help farmers improve profitability, Lower dependence on foreign protein sources). 
PC1LOADINGS <- myFA_impacts$loadings[,1]
PC1LOADINGS[PC1LOADINGS>0.5 | PC1LOADINGS < -0.5] # questions with strong loadings
qu.df[c("Q4_1","Q4_2","Q4_4", "Q4_5","Q4_6",  "Q4_12")]
rm(PC1LOADINGS)    # tidy up


# PC2 strong loadings:
# >>> PC2: Swill would: increase disease risk, and risk being unpalatable to producers
# ... strong loadings for Increase the risk of prion diseases like BSE or vCJD,  Reduce the traceability, Negatively affect the marketability of pork, Increase the risk of toxins entering the feed, Increase the risk of an outbreak of foot-and-mouth disease, Lower consumer acceptance of pork products
PC2LOADINGS <- myFA_impacts$loadings[,2]
PC2LOADINGS[PC2LOADINGS>0.5 | PC2LOADINGS < -0.5] # questions with strong loadings
qu.df[c("Q4_7","Q4_8","Q4_9","Q4_10", "Q4_11", "Q4_13")]
rm(PC2LOADINGS)    # tidy up


# Add the first two factor levels to the dataset
# ... ultimately used in a regression model to understand the support for legal change
no.qu.df$PC1_impacts <- myFA_impacts$scores[,1]
no.qu.df$PC2_impacts <- myFA_impacts$scores[,2]


# tidy up
rm(myFA_impacts, impact_cols)
rm(PCA.df3)


#  ------------------------------------ 13. Model of support for legal change ------------------------------------


# Look at the df
str(no.qu.df)


# !!! Correct the model numbers, so they match the text.
# !!! And make it clear which plots are included in the manuscript


# Remove the cols that are not required for the model
no.qu.df_copy <- no.qu.df # Make a copy since I will change the columns (to make them stan friendly)
no.qu.df_copy$job_group_ID <- ifelse(no.qu.df_copy$job_group == "Pig farmer", 1, 0)       # Stan takes only integer group types
no.qu.df_copy <- no.qu.df_copy[,c("Q53_1","Q53_2","Q48","Q55", "Q56", "Q58", 
                                  "Q37", "Q51","Q59", "Q38", "job_group_ID", 
                                  "PC1_values", "PC2_values", "PC1_impacts", "PC2_impacts")]
dim(no.qu.df_copy)        # 163 respondents, 15 columns


# NB the question numbers correspond to the following questions:
# Q48 == question about support for re-legalisazion
# Q55 == gender
# Q56 == age bracket
# ... Questions specifically for pig farmers:
# Q58 == NUMBER OF PIGS
# Q37 == experience using swill on-farm?
# Q51 == would use swill on your farm?
# Q38 == were affected by FMD outbreak?
# Q59 == Do you use wet or dry feeding?


# Need to convert data to integers, to keep STAN happy.
unique(no.qu.df_copy$Q48)                              # Look at factor levels in the question about support to legal change
no.qu.df_copy$legal_support <- no.qu.df_copy$Q48       # Make a new column in which I convert these to numeric
no.qu.df_copy$legal_support <- factor(no.qu.df_copy$legal_support, 
                                      levels=c("Definitely not","Probably not","Might or might not","Probably yes","Definitely yes"),
                                      labels=c(1:5))
no.qu.df_copy$legal_support <- as.integer(as.character(no.qu.df_copy$legal_support))
head(no.qu.df_copy[,c("Q48","legal_support")])    # Check the original data against the numeric transformation
class(no.qu.df_copy$legal_support)                # Check the data are integers (STAN's preferred class)
unique(no.qu.df_copy$Q55)                         # Look at gender factor levels
no.qu.df_copy$gender <- no.qu.df_copy$Q55         # Make a new column for gender
no.qu.df_copy$gender <- factor(no.qu.df_copy$gender, 
                               levels = c("Male", "Female"),
                               labels = c(0,1))   # Male == 0, Female == 1
no.qu.df_copy$gender <- as.integer(as.character(no.qu.df_copy$gender))
class(no.qu.df_copy$gender)                       # Check the data are integers (STAN's preferred class)
head(no.qu.df_copy[,c("Q55","gender")])           # Check the original data against the numeric transformation

unique(no.qu.df_copy$Q56); class(no.qu.df_copy$Q56)    # Age bracket is currently a factor
no.qu.df_copy$AgeClass <- no.qu.df_copy$Q56            # Make a new column in which I convert these to numeric
no.qu.df_copy$AgeClass <- factor(no.qu.df_copy$AgeClass,         # Therefore convert to an ordered factor
                                 levels = c("0-18","19-30","31-50","50+"),
                                 labels = c(1:4))
no.qu.df_copy$AgeClass <- as.integer(as.character(no.qu.df_copy$AgeClass))
head(no.qu.df_copy[,c("Q56","AgeClass")])    # Check the original data against the numeric transformation
class(no.qu.df_copy$AgeClass)                # Check the data are integers (STAN's preferred class)


# 13a) Fit model of support for legal change for all respondents ----


# Create a smaller dataset for STAN modelling of all respondent's support for legal change
df_legal_change <- no.qu.df_copy[,c("legal_support","gender", "AgeClass", "job_group_ID",
                                    "PC1_values", "PC2_values", "PC1_impacts", "PC2_impacts")]
dim(df_legal_change)        # should be 163 rows, 8 cols


# Fit model


# Check for collinearity between the different factors
# ... they are not correlated
pairs(df_legal_change[,c("PC1_values","PC2_values", "PC1_impacts", "PC2_impacts")], 
      panel=panel.smooth)


# Select starting values
pr_LC <- table( df_legal_change$legal_support ) / nrow(df_legal_change)
cum_pr_LC <- cumsum( pr_LC )                  # cumsum converts to cumulative proportions
logit( cum_pr_LC )                           # the logit intercept parameter estimate (averaged across all feeds - this is just used for starting values)


# Maximal model
# Model with age, job*gender interaction & 4x factor loadings
modelAR1 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aA[AgeClass] +  bJ*job_group_ID + bG*gender +          # phi is my linear model
      bJG*job_group_ID*gender + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aA[AgeClass]  ~ dnorm(0, sigma),              # prior for AgeClass varying intercept
    c(bJ, bG, bJG, bPC1v, bPC2v,bPC1i,bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma ~ dcauchy(0,1),                         # priors for AgeClass variance
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=df_legal_change ,
  start=list(cutpoints = c(logit(cum_pr_LC)[[1]], logit(cum_pr_LC)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC)[[3]], logit(cum_pr_LC)[[4]]) ),
  iter = 8000, chains=1, cores = 1)


# Diagnose chains
# Trace plot looks good, Rhat == 1
plot(modelAR1)                  # Looks OK
precis(modelAR1, depth=2)       # Rhat == 1
logistic(coef(modelAR1)[c(1:4)])        # convert the cutpoints into cumulative probabilities
par(mfrow=c(1,1))
plot(precis(modelAR1, depth=2), 
     main="model AR1") # look at the coefficient estimates, depth = 2 plots all intercepts


# Model with age, job, gender, & 4x factor loadings
modelAR2 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aA[AgeClass] +                          # phi is my linear model
      bJ*job_group_ID + bG*gender +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aA[AgeClass]  ~ dnorm(0, sigma),              # prior for AgeClass varying intercept
    c(bJ, bG, bPC1v, bPC2v,bPC1i,bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma ~ dcauchy(0,1),                         # priors for AgeClass variance
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=df_legal_change ,
  start=list(cutpoints = c(logit(cum_pr_LC)[[1]], logit(cum_pr_LC)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC)[[3]], logit(cum_pr_LC)[[4]]) ),
  iter = 8000, chains=1, cores = 1)


# Diagnose chains
plot(modelAR2)                  # Check the trace plot
precis(modelAR2, depth=2)       # Rhat == 1!


# Model without age, but with job*gender interaction
modelAR3 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- bJ*job_group_ID + bG*gender +          # phi is my linear model
      bJG*job_group_ID*gender + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bJ, bG, bJG, bPC1v, bPC2v,bPC1i,bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=df_legal_change ,
  start=list(cutpoints = c(logit(cum_pr_LC)[[1]], logit(cum_pr_LC)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC)[[3]], logit(cum_pr_LC)[[4]]) ),
  iter = 8000, chains=1, cores = 1)


# Diagnose chains
plot(modelAR3)                  # Looks OK
precis(modelAR3, depth=2)       # Rhat == 1


# model AR4
# Model without age
modelAR4 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- bJ*job_group_ID + bG*gender +          # phi is my linear model
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bJ, bG, bPC1v, bPC2v,bPC1i,bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=df_legal_change ,
  start=list(cutpoints = c(logit(cum_pr_LC)[[1]], logit(cum_pr_LC)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC)[[3]], logit(cum_pr_LC)[[4]]) ),
  iter = 8000, chains=1, cores = 1)


# Diagnose chains
plot(modelAR4)                  # Looks OK
precis(modelAR4, depth=2)       # Rhat == 1

# Model AR5
# Model without job
modelAR5 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aA[AgeClass] +                          # phi is my linear model
      bG*gender +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aA[AgeClass]  ~ dnorm(0, sigma),              # prior for AgeClass varying intercept
    c(bG, bPC1v, bPC2v,bPC1i,bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma ~ dcauchy(0,1),                         # priors for AgeClass variance
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=df_legal_change ,
  start=list(cutpoints = c(logit(cum_pr_LC)[[1]], logit(cum_pr_LC)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC)[[3]], logit(cum_pr_LC)[[4]]) ),
  iter = 8000, chains=1, cores = 1)


# Diagnose chains
plot(modelAR5)
precis(modelAR5, depth=2)       # Rhat == 1


# Model with age group & job & 4x factor loadings
modelAR6 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aA[AgeClass] +  bJ*job_group_ID +      # phi is my linear model
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aA[AgeClass]  ~ dnorm(0, sigma),              # prior for AgeClass varying intercept
    c(bJ, bPC1v, bPC2v,bPC1i,bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma ~ dcauchy(0,1),                         # priors for AgeClass variance
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=df_legal_change ,
  start=list(cutpoints = c(logit(cum_pr_LC)[[1]], logit(cum_pr_LC)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC)[[3]], logit(cum_pr_LC)[[4]]) ),
  iter = 8000, chains=1, cores = 1)


# Diagnose chains
# Trace plot looks good, Rhat == 1
plot(modelAR6)                  # Looks OK
precis(modelAR6, depth=2)       # Rhat == 1


# Model AR7
# Model with Job & 4x factor loadings
modelAR7 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- bJ*job_group_ID   +                    # phi is my linear model
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bJ, bPC1v, bPC2v,bPC1i,bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=df_legal_change ,
  start=list(cutpoints = c(logit(cum_pr_LC)[[1]], logit(cum_pr_LC)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC)[[3]], logit(cum_pr_LC)[[4]]) ),
  iter = 8000, chains=1, cores = 1)


# Diagnose chains
# Trace plot looks good, Rhat == 1
plot(modelAR7)                  # Looks OK
precis(modelAR7, depth=2)       # Rhat == 1


# AR8
# Model without age and job
modelAR8 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                        # phi is my linear model
      bG*gender +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bPC1v, bPC2v,bPC1i,bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_legal_change ,
  start=list(cutpoints = c(logit(cum_pr_LC)[[1]], logit(cum_pr_LC)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC)[[3]], logit(cum_pr_LC)[[4]]) ),
  iter = 8000, chains=1, cores = 1)


# Plot
plot(modelAR8)
precis(modelAR8, depth=2)       # Rhat == 1!


# AR9
# Model with 4x factor loadings
modelAR9 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                        # phi is my linear model
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bPC1v, bPC2v,bPC1i,bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                       # prior for thresholds of cumulative likelihood
  ) ,
  data=df_legal_change ,
  start=list(cutpoints = c(logit(cum_pr_LC)[[1]], logit(cum_pr_LC)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC)[[3]], logit(cum_pr_LC)[[4]]) ),
  iter = 8000, chains=1, cores = 1)


# Diagnose chains
# Trace plot looks good, Rhat == 1
plot(modelAR9)                  # Looks OK
precis(modelAR9, depth=2)       # Rhat == 1


# Compare models
# ... the model with job*gender interaction gets the most weight (0.52)
compare(modelAR1, modelAR2, modelAR3, modelAR4, modelAR5, modelAR6, modelAR7, modelAR8, modelAR9)
coeftab(modelAR1, modelAR2, modelAR3, modelAR4, modelAR5, modelAR6, modelAR7, modelAR8, modelAR9)
plot(coeftab(modelAR3, modelAR8, modelAR1, modelAR9, modelAR5))


# ------------------ 13bi) Plot Fig 4, all Respondents ------------------


# Extract the coefficients for models with greatest weights
mAR3_res <- precis(modelAR3, depth=2, pars=c("bG","bJ", "bJG", "bPC1v", "bPC2v","bPC1i", "bPC2i"))@output
mAR3_res$coef <- row.names(mAR3_res)    # Extract the coefficient names
mAR3_res$model <- "model AR3"
mAR3_res$weighting <- 0.35              # this is from compare() above, exact values may vary between runs
mAR8_res <- precis(modelAR8, depth=2, pars=c("bG","bPC1v", "bPC2v","bPC1i", "bPC2i"))@output
mAR8_res$coef <- row.names(mAR8_res)    # Extract the coefficient names
mAR8_res$model <- "model AR8"
mAR8_res$weighting <- 0.22              # this is from compare() above
mAR1_res <- precis(modelAR1, depth=2, pars=c("aA","bG","bJ", "bJG", "bPC1v", "bPC2v","bPC1i", "bPC2i"))@output
mAR1_res$coef <- row.names(mAR1_res)    # Extract the coefficient names
mAR1_res$model <- "model AR1"
mAR1_res$weighting <- 0.11              # this is from compare() above
mAR9_res <- precis(modelAR9, depth=2, pars=c("bPC1v", "bPC2v","bPC1i", "bPC2i"))@output
mAR9_res$coef <- row.names(mAR9_res)    # Extract the coefficient names
mAR9_res$model <- "model AR9"
mAR9_res$weighting <- 0.10              # this is from compare() above
mAR5_res <- precis(modelAR5, depth=2, pars=c("aA","bG","bPC1v", "bPC2v","bPC1i", "bPC2i"))@output
mAR5_res$coef <- row.names(mAR5_res)    # Extract the coefficient names
mAR5_res$model <- "model AR5"
mAR5_res$weighting <- 0.08              # this is from compare() above
mAR_res <- rbind(mAR3_res,mAR8_res,mAR1_res,mAR9_res, mAR5_res)
head(mAR_res)
mAR_res$coef <- factor(mAR_res$coef, 
                         levels=c("aA[1]", "aA[2]", "aA[3]", "aA[4]","bJ", "bJG", "bG", "bPC2i","bPC1i","bPC2v","bPC1v"), 
                         labels=c("Age: 0-18",
                                  "Age: 19-30",
                                  "Age: 31-50",
                                  "Age: 50+",
                                  "Job (Pig-farmer)",
                                  "Gender*Job interaction",
                                  "Gender (female)", 
                                  "Factor 2, perception of swill:\nSwill would increase disease risk\n& be unpalatable to consumers",
                                  "Factor 1, perception of swill:\nSwill would benefit the environment, help\nfarms financially, & reduce trade-deficits",
                                  "Factor 2, respondent's values:\nRespondents placing high importance on\nresource efficiency & farm finances",
                                  "Factor 1, respondent's values:\nRespondents placing high importance on\ndisease control & the perception of the industry") )
mAR_res$model <- factor(mAR_res$model,
                        levels = c("model AR9","model AR8","model AR5","model AR3","model AR1"),
                        labels = c("AR9","AR8","AR5","AR3","AR1"))
colnames(mAR_res) <- c("Mean","StdDev", "lwr.89", "upr.89", "n_eff", "Rhat", "coef", "model", "weighting")


# tidy up
rm(mAR3_res,mAR8_res,mAR1_res,mAR9_res, mAR5_res)


# save model output
# write.table(mAR_res, "results_tables/Fig4_20171009.csv", sep=",", row.names=F)
# mAR_res <- read.csv("results_tables/Fig4_20171009.csv", header=T)


# Relabel for ceofficients  type 
mAR_res$coef_type <- "Respondent\nCharacteristics"
mAR_res[mAR_res$coef %in% c("Factor 1, respondent's values:\nRespondents placing high importance on\ndisease control & the perception of the industry",
                            "Factor 2, respondent's values:\nRespondents placing high importance on\nresource efficiency & farm finances"),"coef_type"] <- "Factors:\nRespondent values"
mAR_res[mAR_res$coef %in% c("Factor 1, perception of swill:\nSwill would benefit the environment, help\nfarms financially, & reduce trade-deficits",
                            "Factor 2, perception of swill:\nSwill would increase disease risk\n& be unpalatable to consumers"),"coef_type"] <- "Factors:\nPerception of swill"
head(mAR_res)



# if loaded mRes from saved version, need to reorder the coef and model factor levels
# ... and relabel the coefficients (shorter labels)
mAR_res$coef <- factor(mAR_res$coef, 
                       levels=c("Age: 0-18",
                                "Age: 19-30",
                                "Age: 31-50",
                                "Age: 50+",
                                "Job (Pig-farmer)",
                                "Gender*Job interaction",
                                "Gender (female)", 
                                "Factor 2, perception of swill:\nSwill would increase disease risk\n& be unpalatable to consumers",
                                "Factor 1, perception of swill:\nSwill would benefit the environment, help\nfarms financially, & reduce trade-deficits",
                                "Factor 1, respondent's values:\nRespondents placing high importance on\ndisease control & the perception of the industry",
                                "Factor 2, respondent's values:\nRespondents placing high importance on\nresource efficiency & farm finances"),
                       labels=c("Age: 0-18",
                                "Age: 19-30",
                                "Age: 31-50",
                                "Age: 50+",
                                "Job (Pig-farmer)",
                                "Gender:Job interaction",
                                "Gender (female)", 
                                "Swill would increase disease risk &\nbe unpalatable to consumers",
                                "Swill would benefit the environment, help\nfarms financially, & reduce trade-deficits",
                                "Place importance on disease control &\nthe perception of the industry",
                                "Place importance on farm finances &\nresource efficiency"))
mAR_res$model <- factor(mAR_res$model, 
                        levels = c("AR1", "AR3", "AR5", "AR8", "AR9"))


# # Plot Figure 4, not including age class
# # ... model weights shown as different point sizes


# # Step 1) Plot with the weights legend
mAR_res_v2 <- mAR_res[-which(mAR_res$coef %in% c("Age: 0-18", "Age: 19-30", "Age: 31-50", "Age: 50+")),]
myCols5 <- c("#999999", "#E69F00", "#0072B2", "#D55E00", "#009E73")  # Make colour palette for 5 colours
Fig4a <- ggplot(mAR_res_v2, aes( Mean, coef, col = model,size = weighting*100)) +
  # geom_point(position=position_dodge(.9)) +
  stat_sum(position = position_dodgev(height = 0.9)) +
  geom_errorbarh(aes(xmin=lwr.89, xmax=upr.89, col = model),
                size=1,    # Thinner lines
                width=.3,
                position=position_dodgev(height = 0.9)
  )  +
  facet_grid(coef_type ~., scales= "free", switch = "y") +
  # coord_flip() +
  geom_vline(xintercept = 0) +
  xlab("Coefficient value") +
  theme_bw()  +
  theme(legend.position = c(0.9,0.65),
        legend.title = element_text(face="bold"),
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", colour = NA),
        axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12),
        axis.line = element_line(colour="black")) +
  scale_color_manual(values = myCols5,
                     guide=FALSE) +     # removes the colour legend
  scale_size_continuous(limits = c(8,42),              # standardise the weights displayed in the legend
                        name="Model\nweight (%)")


# plot with the colour legend
Fig4b <- ggplot(mAR_res_v2, aes( Mean, coef, col = model,size = weighting*100)) +
  # geom_point(position=position_dodge(.9)) +
  stat_sum(position = position_dodgev(height = 0.9)) +
  geom_errorbarh(aes(xmin=lwr.89, xmax=upr.89, col = model),
                 size=.4,    # Thinner lines
                 width=.3,
                 position=position_dodgev(height = 0.9)
  )  +
  facet_grid(coef_type ~., scales= "free", switch = "y") +
  # coord_flip() +
  geom_vline(xintercept = 0) +
  xlab("Coefficient value") +
  theme_bw()  +
  theme(legend.position = c(0.9,0.65),
        legend.title = element_text(face="bold"),
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        strip.text = element_text(size = 12),
        axis.line = element_line(colour="black")) +
  scale_color_manual(values = myCols5,
                     name = "Model ID",
                     breaks=c("AR1", "AR3", "AR5", "AR8", "AR9")
                     ) +
  scale_size_continuous(limits = c(8,42),              # standardise the weights displayed in the legend
                        name="Model\nweight (%)",
                        guide=FALSE)

# Step 2) Extract the size legend - leg2
leg <- gtable_filter(ggplot_gtable(ggplot_build(Fig4b)), "guide-box")


# Step 3) Plot and add second legend
Fig4 <- Fig4a +
  annotation_custom(grob = leg,
                    xmin = 1.5, xmax = 2,        # Place legend
                    ymin = 5, ymax = 7.72)      # Nb because of coord_flip() the x and y axes are flipped
grid.newpage()
grid.draw(Fig4)


# In this case, I added the two titles in manually
# # Add in two titles
# # ... code from: https://stackoverflow.com/questions/21997715/add-ggplot-annotation-outside-the-panel-or-two-titles
# g = ggplotGrob(Fig4)
# g = gtable_add_grob(g, grobTree(textGrob("Lower support for relegalisation", x=0, hjust=0),
#                                 textGrob("Greater support for relegalisation", x=1, hjust=1)),
#                     t=1, l=4)
# grid.draw(g)


# tidy up
rm(Fig4, Fig4a, Fig4b, g, leg)


# ------------------ 13bii) Plot S2 Appendix Figure 10, including age class ------------------ 
# # Step 1 plot without color legend
S2App_Fig10a <- ggplot(mAR_res, aes(Mean, coef, col = model,size = weighting*100)) + 
  stat_sum(position = position_dodgev(height = 0.9)) +
  geom_errorbarh(aes(xmin=lwr.89, xmax=upr.89, col = model),
                 size=.4,    # Thinner lines
                 width=.3,
                 position=position_dodgev(height = 0.9)
  )  +
  facet_grid(coef_type ~., scales= "free", switch = "y") +
  # coord_flip() +
  geom_vline(xintercept = 0) +
  xlab("Coefficient value") +
  theme_bw()  +
  theme(legend.position = c(0.9,0.65),
        legend.title = element_text(face="bold"),
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        strip.text = element_text(size = 12),
        axis.line = element_line(colour="black")) +
  scale_color_manual(values = myCols5,
                     name = "Model ID",
                     guide = FALSE) +
  scale_size_continuous(limits = c(8,42),              # standardise the weights displayed in the legend
                        name="Model\nweight (%)")


# plot with the colour legend
S2App_Fig10b <- ggplot(mAR_res, aes(coef, Mean, col = model,size = weighting*100)) + 
  stat_sum(position = position_dodgev(height = 0.9)) +
  geom_errorbarh(aes(xmin=lwr.89, xmax=upr.89, col = model),
                 size=.4,    # Thinner lines
                 width=.3,
                 position=position_dodgev(height = 0.9)
  )  +
  facet_grid(coef_type ~., scales= "free", switch = "y") +
  # coord_flip() +
  geom_vline(xintercept = 0) +
  xlab("Coefficient value") +
  theme_bw()  +
  theme(legend.position = c(-1,14.6),           # depending on size of plot window, must fiddle with this value to line up legends
        legend.title = element_text(face="bold"),
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        strip.text = element_text(size = 12),
        axis.line = element_line(colour="black")) +
  scale_color_manual(values = myCols5,
                     name = "Model ID") +
  scale_size_continuous(limits = c(8,42),              # standardise the weights displayed in the legend
                        name="Model\nweight (%)",
                        guide = FALSE)

# Step 2) Extract the size legend - leg2
leg <- gtable_filter(ggplot_gtable(ggplot_build(S2App_Fig10b)), "guide-box")


# Step 3) Plot and add second legend
S2App_Fig10 <- S2App_Fig10a +
  annotation_custom(grob = leg,
                    ymin = 1.5, ymax = 2.4,    # Place legend
                    xmin = 2.45, xmax = 3)      # Nb because of coord_flip() the x and y axes are flipped
grid.newpage()
grid.draw(S2App_Fig10)


# Add titles manually (facet_grid overrides this code)
# ... code from: https://stackoverflow.com/questions/21997715/add-ggplot-annotation-outside-the-panel-or-two-titles
# require(grid)
# require(gtable)
# g = ggplotGrob(S2App_Fig10)
# g = gtable_add_grob(g, grobTree(textGrob("Lower support for relegalisation", x=0, hjust=0),
#                                 textGrob("Greater support for relegalisation", x=1, hjust=1)),
#                     t=1, l=4)
# grid.draw(g)


# tidy up
rm(mAR_res,S2App_Fig10, S2App_Fig10a, S2App_Fig10b, g)


# Check for differences between male and female respondents in the original dataset
# 1/7 female farmers had FMD on farm (14%)
# 21/74 male farmers had FMD on farm (28%)
table(no.qu.df_copy$Q38, no.qu.df_copy$gender)


# tidy up
rm(df_legal_change, cum_pr_LC, pr_LC)
rm(myCols5)
rm(modelAR1, modelAR2, modelAR3, modelAR4, modelAR5, modelAR6, modelAR7, modelAR8, modelAR9)


# 13c) Model of support for legal change amongst pig farmers only ----


# Q58 == NUMBER OF PIGS
# Q37 == experience using swill on-farm?
# Q51 == would use swill on your farm?
# Q38 == were affected by FMD outbreak?
# Q59 == Do you use wet or dry feeding?


# First, make dataset for pig farmers only
no.qu.df_pig_farm <- no.qu.df_copy[no.qu.df_copy$job_group_ID == 1,]
dim(no.qu.df_pig_farm)                        # Should be 82 rows


# Create a smaller dataset for STAN modelling of all respondent's support for legal change
# ... Need to convert data to integers, to keep STAN happy.
no.qu.df_pig_farm$swill_exp <- no.qu.df_pig_farm$Q37           # Make a new column in which I convert these to numeric
no.qu.df_pig_farm$swill_exp <- factor(no.qu.df_pig_farm$swill_exp, 
                                      levels=c("Yes", "No"),
                                      labels=c(1,0))
no.qu.df_pig_farm$swill_exp <- as.integer(as.character(no.qu.df_pig_farm$swill_exp))
head(no.qu.df_pig_farm[,c("Q37","swill_exp")])    # Check the original data against the numeric transformation
class(no.qu.df_pig_farm$swill_exp)                # Check the data are integers (STAN's preferred class)
levels(no.qu.df_pig_farm$Q58) <- c(levels(no.qu.df_pig_farm$Q58), "1-9", "10-99") # Add new factor levels, since some of the sizes are incorrecly interpreted by [R] as dates
no.qu.df_pig_farm$Q58[no.qu.df_pig_farm$Q58 =="9-Jan"] <- "1-9"                   # Correct farm sizes from dates to number of pigs
no.qu.df_pig_farm$Q58[no.qu.df_pig_farm$Q58 =="Oct-99"] <- "10-99"                # Correct farm sizes from dates to number of pigs
no.qu.df_pig_farm$farm_size <- no.qu.df_pig_farm$Q58           # Make a new column in which I convert these to integers
no.qu.df_pig_farm$farm_size <- factor(no.qu.df_pig_farm$farm_size, levels = c("1-9", "10-99","100-199",
                                                                              "200-399", "400-999",
                                                                              "1000-4999", "5000+"),
                                      labels=c(1:7))
no.qu.df_pig_farm$farm_size <- as.integer(as.character(no.qu.df_pig_farm$farm_size))
head(no.qu.df_pig_farm[,c("Q58","farm_size")])    # Check the original data against the numeric transformation
class(no.qu.df_pig_farm$farm_size)                # Check the data are integers (STAN's preferred class)
no.qu.df_pig_farm$will_use_swill <- no.qu.df_pig_farm$Q51           # Make a new column in which I convert these to numeric
no.qu.df_pig_farm$will_use_swill <- factor(no.qu.df_pig_farm$will_use_swill, 
                                           levels=c("Definitely not", "Probably not", "Might or might not", "Probably yes", "Absolutely"),
                                           labels=c(1:5))
no.qu.df_pig_farm$will_use_swill <- as.integer(as.character(no.qu.df_pig_farm$will_use_swill))
head(no.qu.df_pig_farm[,c("Q51","will_use_swill")])    # Check the original data against the numeric transformation
no.qu.df_pig_farm$FMD_exp <- no.qu.df_pig_farm$Q38     # Make a new column in which I convert these to numeric
table(no.qu.df_pig_farm$FMD_exp)                       # There were 3 "Not sure" responses
no.qu.df_pig_farm$FMD_exp[no.qu.df_pig_farm$FMD_exp=="Not sure"] <- "No" # Reclassify Not sure responses as "No"
no.qu.df_pig_farm$FMD_exp <- factor(no.qu.df_pig_farm$FMD_exp, 
                                    levels=c("Yes", "No"),
                                    labels=c(1,0))
no.qu.df_pig_farm$FMD_exp <- as.integer(as.character(no.qu.df_pig_farm$FMD_exp))
head(no.qu.df_pig_farm[,c("Q38","FMD_exp")])             # Check the original data against the numeric transformation
no.qu.df_pig_farm$feed_tech <- no.qu.df_pig_farm$Q59     # Make a new column in which I convert these to integers
no.qu.df_pig_farm$feed_tech <- factor(no.qu.df_pig_farm$feed_tech, 
                                      levels=c("Wet feeding", "Dry feeding"),
                                      labels=c(1,0))
no.qu.df_pig_farm$feed_tech <- as.integer(as.character(no.qu.df_pig_farm$feed_tech))
head(no.qu.df_pig_farm[,c("Q59","feed_tech")])           # Check the original data against the numeric transformation


df_LC_pf <- no.qu.df_pig_farm[,c("legal_support","gender", "AgeClass", 
                                 "PC1_values", "PC2_values", "PC1_impacts", "PC2_impacts",
                                 "farm_size", "will_use_swill","FMD_exp", "swill_exp", "feed_tech")]
dim(df_LC_pf)         # should be 82 rows, 12 cols
rm(no.qu.df_pig_farm) # tidy up


# Fit model for pig farmers support for re-legalisation
# Select starting values
pr_LC_farm <- table( df_LC_pf$legal_support ) / nrow(df_LC_pf)
cum_pr_LC_farm <- cumsum( pr_LC_farm )                  # cumsum converts to cumulative proportions
logit( cum_pr_LC_farm )                           # the logit intercept parameter estimate (averaged across all feeds - this is just used for starting values)


# model FS1
mLC_farm <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aA[AgeClass] + aF[farm_size] +         # phi is my linear model
      bG*gender + bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aA[AgeClass]  ~ dnorm(0, sigma_A),              # prior for AgeClass varying intercept
    aF[farm_size]  ~ dnorm(0, sigma_F),             # prior for farm_size varying intercept
    c(bG, bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma_A ~ dcauchy(0,1),                         # priors for AgeClass variance
    sigma_F ~ dcauchy(0,1),                         # priors for farm_size variance
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Plot
plot(mLC_farm)
precis(mLC_farm, depth=2)               # Some Rhat == 1
logistic(coef(mLC_farm)[c(1:4)])        # convert the cutpoints into cumulative probabilities
par(mfrow=c(1,1))
plot(precis(mLC_farm, depth=2), 
     main="model mLC_farm") # look at the coefficient estimates, depth = 2 plots all intercepts



# Model FS2
# Model without farm size
mLC_farm2 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aA[AgeClass] +          # phi is my linear model
      bG*gender + bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aA[AgeClass]  ~ dnorm(0, sigma_A),              # prior for AgeClass varying intercept
    c(bG, bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma_A ~ dcauchy(0,1),                         # priors for AgeClass variance
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Plot
plot(mLC_farm2)
precis(mLC_farm2, depth=2)               # All Rhat == 1


# Model FS3
# model without AgeClass
mLC_farm3 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aF[farm_size] +         # phi is my linear model
      bG*gender + bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aF[farm_size]  ~ dnorm(0, sigma_F),             # prior for farm_size varying intercept
    c(bG, bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma_F ~ dcauchy(0,1),                         # priors for farm_size variance
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Plot
plot(mLC_farm3)
precis(mLC_farm3, depth=2)               # All Rhat == 1


# Model FS4
# model without AgeClass & farm_size
mLC_farm4 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm4)
precis(mLC_farm4, depth=2)               # All Rhat == 1


# Model FS5
# model without AgeClass, farm_size, gender
mLC_farm5 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm5)
precis(mLC_farm5, depth=2)               # All Rhat == 1


# Model FS6
# model without AgeClass, farm_size, swill_exp
mLC_farm6 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bFMD*FMD_exp + bFT*feed_tech +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm6)
precis(mLC_farm6, depth=2)               # All Rhat == 1


# Model FS7
# model without AgeClass, farm_size, feed_tech
mLC_farm7 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bFMD*FMD_exp +  bG*gender + bS*swill_exp + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bFMD, bS, bG, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm7)
precis(mLC_farm7, depth=2)               # All Rhat == 1


# Model FS8
# model without AgeClass, farm_size, FMD_exp
mLC_farm8 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bS*swill_exp + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bS, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm8)
precis(mLC_farm8, depth=2)               # All Rhat == 1


# Model FS9
# model without AgeClass, farm_size, FMD_exp, swill_exp
mLC_farm9 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm9)
precis(mLC_farm9, depth=2)               # All Rhat == 1


# Model FS10
# model without AgeClass, farm_size, gender, swill_exp
mLC_farm10 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bFMD*FMD_exp + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm10)
precis(mLC_farm10, depth=2)               # All Rhat == 1


# Model FS11
# model without AgeClass, farm_size, gender, FMD_exp
mLC_farm11 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bS*swill_exp + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bS, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm11)
precis(mLC_farm11, depth=2)               # All Rhat == 1


# Model FS12
# model without AgeClass, farm_size, gender, feed_tech
mLC_farm12 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bS*swill_exp + bFMD*FMD_exp + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bS, bFMD, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm12)
precis(mLC_farm12, depth=2)               # All Rhat == 1


# Model FS13
# model without AgeClass, farm_size, swill_exp, feed_tech
mLC_farm13 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bFMD*FMD_exp + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bFMD, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm13)
precis(mLC_farm13, depth=2)               # All Rhat == 1


# Model FS14
# model without AgeClass, farm_size, FMD_exp, feed_tech
mLC_farm14 <- map2stan(
  alist(
    legal_support ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bS*swill_exp + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bS, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_LC_farm)[[1]], logit(cum_pr_LC_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_LC_farm)[[3]], logit(cum_pr_LC_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mLC_farm14)
precis(mLC_farm14, depth=2)               # All Rhat == 1


# Compare models
compare(mLC_farm, mLC_farm2, mLC_farm3, mLC_farm4, 
        mLC_farm5, mLC_farm6, mLC_farm7, mLC_farm8, 
        mLC_farm9, mLC_farm10, mLC_farm11, mLC_farm12, 
        mLC_farm13, mLC_farm14)                    # compare models w/ combinations of swill_exp, gender, farm_size
plot(coeftab(mLC_farm9, mLC_farm10, mLC_farm11))   # compare top three models


# ---------------------------- 13di) Plot Fig 5 ---------------------------- 


# Plot the coefficients of the three models with largest weights
comparison <- compare(mLC_farm, mLC_farm2, mLC_farm3, mLC_farm4,   # extract model weights
                      mLC_farm5, mLC_farm6, mLC_farm7, mLC_farm8, 
                      mLC_farm9, mLC_farm10, mLC_farm11, mLC_farm12, 
                      mLC_farm13, mLC_farm14)
mLC_farm10_res <- precis(mLC_farm10, depth=2, pars=c("bFMD","bFT", "bPC1v", "bPC2v","bPC1i", "bPC2i"))@output
mLC_farm10_res$coef <- row.names(mLC_farm10_res)    # Extract the coefficient names
mLC_farm10_res$model <- "FS10"
mLC_farm10_res$weighting <- comparison@output[which(row.names(comparison@output) == "mLC_farm10"),"weight"]
mLC_farm9_res <- precis(mLC_farm9, depth=2, pars=c("bG","bFT", "bPC1v", "bPC2v","bPC1i", "bPC2i"))@output
mLC_farm9_res$coef <- row.names(mLC_farm9_res)    # Extract the coefficient names
mLC_farm9_res$model <- "FS9"
mLC_farm9_res$weighting <- comparison@output[which(row.names(comparison@output) == "mLC_farm9"),"weight"]
mLC_farm11_res <- precis(mLC_farm11, depth=2, pars=c("bS","bFT", "bPC1v", "bPC2v","bPC1i", "bPC2i"))@output
mLC_farm11_res$coef <- row.names(mLC_farm11_res)    # Extract the coefficient names
mLC_farm11_res$model <- "FS11"
mLC_farm11_res$weighting <- comparison@output[which(row.names(comparison@output) == "mLC_farm11"),"weight"]
mLC_farm91011_res <- rbind(mLC_farm10_res, mLC_farm9_res, mLC_farm11_res)
head(mLC_farm91011_res)
mLC_farm91011_res$coef <- factor(mLC_farm91011_res$coef, 
                                 levels=c("bG", "bS", "bFMD","bFT", "bPC2i","bPC1i","bPC2v","bPC1v"), 
                                 labels=c("Gender (female)", 
                                          "Experience using swill on farm",
                                          "Farm directly affected\nby 2001 FMD outbreak",
                                          "Feeding technique\n(wet feed)",
                                          "Factor 2, perception of swill:\nSwill would increase disease risk\n& be unpalatable to consumers",
                                          "Factor 1, perception of swill:\nSwill would benefit the environment, help\nfarms financially, & reduce trade-deficits",
                                          "Factor 2, respondent's values:\nRespondents placing high importance on\nresource efficiency & farm finances",
                                          "Factor 1, respondent's values:\nRespondents placing high importance on\ndisease control & perception of the industry") )
colnames(mLC_farm91011_res) <- c("Mean","StdDev", "lwr.89", "upr.89", "n_eff", "Rhat", "coef", "model", "weight")
mLC_farm91011_res$model <- factor(mLC_farm91011_res$model,                         # change the order of models, to match their weights
                                  levels = rev(c("FS9", "FS10", "FS11")))

# Save model output
# write.table(mLC_farm91011_res, "results_tables/Fig5_output.csv", sep=",", row.names=F)
# mLC_farm91011_res <- read.csv("results_tables/Fig5_output.csv", header = T)


# Relabel for ceofficients  type 
mLC_farm91011_res$coef_type <- "Respondent & farm\ncharacteristics"
mLC_farm91011_res[mLC_farm91011_res$coef %in% c("Factor 1, respondent's values:\nRespondents placing high importance on\ndisease control & perception of the industry",
                                                "Factor 2, respondent's values:\nRespondents placing high importance on\nresource efficiency & farm finances"),"coef_type"] <- "Factors:\nRespondent values"
mLC_farm91011_res[mLC_farm91011_res$coef %in% c("Factor 1, perception of swill:\nSwill would benefit the environment, help\nfarms financially, & reduce trade-deficits",
                                                "Factor 2, perception of swill:\nSwill would increase disease risk\n& be unpalatable to consumers"),"coef_type"] <- "Factors:\nPerception of swill"
head(mLC_farm91011_res)


# Shorten coefficient titles
mLC_farm91011_res$coef <- factor(mLC_farm91011_res$coef,
                                 levels=c("Gender (female)",
                                          "Experience using swill on farm",
                                          "Farm directly affected\nby 2001 FMD outbreak",
                                          "Feeding technique\n(wet feed)",
                                          "Factor 2, perception of swill:\nSwill would increase disease risk\n& be unpalatable to consumers",
                                          "Factor 1, perception of swill:\nSwill would benefit the environment, help\nfarms financially, & reduce trade-deficits",
                                          "Factor 1, respondent's values:\nRespondents placing high importance on\ndisease control & perception of the industry",
                                          "Factor 2, respondent's values:\nRespondents placing high importance on\nresource efficiency & farm finances"),
                                 labels=c("Gender (female)",
                                          "Have previously used swill on their farm",
                                          "Farm directly affected\nby 2001 FMD outbreak",
                                          "Feeding technique\n(wet feed)",
                                          "Swill would increase disease risk &\nbe unpalatable to consumers",
                                          "Swill would benefit the environment, help\nfarms financially, & reduce trade-deficits",
                                          "Place importance on disease control &\nthe perception of the industry",
                                          "Place importance on farm finances &\nresource efficiency"))


# tidy up
rm(mLC_farm9_res, mLC_farm10_res, mLC_farm11_res)


# plot Fig 5
# ... a bit of fiddling required to specify the location of both legends
# Step 1) plot Model ID legend
myCols <- c("#0072B2","#D55E00", "#009E73")            # Make colour palette
Fig5a <- ggplot(mLC_farm91011_res, aes(Mean, coef, col = model, size = weight*100)) + 
  stat_sum(position = position_dodgev(height = 0.6)) +
  geom_errorbarh(aes(xmin=lwr.89, xmax=upr.89, col = model),
                 size=1,    # Thinner lines
                 width=.3,
                 position=position_dodgev(height = 0.6)
  )  +
  facet_grid(coef_type ~., scales= "free", switch = "y") +
  # coord_flip() +
  geom_vline(xintercept = 0) +
  xlab("Coefficient value") +
  theme_bw()  +
  theme(legend.position = c(0.7,0.35),
        legend.title = element_text(face="bold"),
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", colour = NA),
        axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12),
        axis.line = element_line(colour="black")) +
  scale_color_manual(values = myCols, 
                     name = "Model ID", 
                     guide = guide_legend(reverse=TRUE)) +
  scale_size_continuous(guide = FALSE,                  # remove the size legend
                        limits = c(8,42),              # standardise the weights displayed in the legend
                        name="Model\nweight (%)")
Fig5a

Fig5b <- ggplot(mLC_farm91011_res, aes(Mean, coef, col = model, size = weight*100)) + 
  stat_sum(position = position_dodgev(height = 0.6)) +
  geom_errorbarh(aes(xmin=lwr.89, xmax=upr.89, col = model),
                 size=1,    # Thinner lines
                 width=.3,
                 position=position_dodgev(height = 0.6)
  )  +
  facet_grid(coef_type ~., scales= "free", switch = "y") +
  # coord_flip() +
  geom_vline(xintercept = 0) +
  xlab("Coefficient value") +
  theme_bw()  +
  theme(legend.position = c(0.85,0.35),
        legend.title = element_text(face="bold"),
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", colour = NA),
        axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12),
        axis.line = element_line(colour="black")) +
  scale_size_continuous(limits = c(8,42),              # standardise the weights displayed in the legend
                        name="Model\nweight (%)") +
  scale_color_manual(values = myCols, 
                     name = "Model ID", 
                     guide = FALSE) 
Fig5b  

# Nb I manually placed the legends in this plot, as I couldn't get the second legend to 
# ... align nicely


# # Step 2) Extract the size legend - leg2
# leg <- gtable_filter(ggplot_gtable(ggplot_build(Fig5b)), "guide-box") 
# 
# 
# # Step 3) Plot and add second legend
# Fig5 <- Fig5a + 
#   annotation_custom(grob = leg, 
#                     ymin = 0, ymax = 2,         # Place legend
#                     xmin = 0, xmax = 5.5)       # Nb because of coord_flip() the x and y axes are flipped
# grid.newpage()
# grid.draw(Fig5)


# # Add in titles: 
# # ... code from: https://stackoverflow.com/questions/21997715/add-ggplot-annotation-outside-the-panel-or-two-titles
# g = ggplotGrob(Fig5)
# g = gtable_add_grob(g, grobTree(textGrob("Lower support for relegalisation", x=0, hjust=0),
#                                 textGrob("Greater support for relegalisation", x=1, hjust=1)),
#                     t=1, l=4)
# grid.draw(g)


# tidy up
rm(Fig5, Fig5a, leg, g, myCols, comparison, mLC_farm91011_res)


# Nb there are only seven women in the sample of pig farmers 
nrow(df_LC_pf[df_LC_pf$gender == 1,])


# ---------------- 13dii) Plot S2 Appendix Fig 11, feed tech effect size ---------------- 


# First, simulate the effect size for feed_tech
# >>> this is equivalent to a 0.37 lower score if the respondent uses wet feed.
farmer_feed_tech_pred_dry <- 
  data.frame(
    legal_support=5,                          #### needed so sim knows num levels
    feed_tech = rep(0,1000),
    FMD_exp = sample(c(0,1), 1000, replace=T),   # 
    swill_exp = sample(c(0,1), 1000, replace=T), #
    gender = sample(c(0,1), 1000, replace=T),
    PC1_values = sample(df_LC_pf$PC1_values, 1000, replace=T),  # Randomly draw 200 values from the original dataset
    PC2_values = sample(df_LC_pf$PC2_values, 1000, replace=T),  # Randomly draw 200 values from the original dataset
    PC1_impacts = sample(df_LC_pf$PC1_impacts, 1000, replace=T),  # Randomly draw 200 values from the original dataset
    PC2_impacts = sample(df_LC_pf$PC2_impacts, 1000, replace=T)   # Randomly draw 200 values from the original dataset
  )
farmer_feed_tech_pred_wet <- 
  data.frame(
    legal_support=5,                          #### needed so sim knows num levels
    feed_tech = rep(1,1000),
    FMD_exp = sample(c(0,1), 10, replace=T),   # 
    swill_exp = sample(c(0,1), 1000, replace=T), #
    gender = sample(c(0,1), 1000, replace=T),
    PC1_values = sample(df_LC_pf$PC1_values, 1000, replace=T),  # Randomly draw 200 values from the original dataset
    PC2_values = sample(df_LC_pf$PC2_values, 1000, replace=T),  # Randomly draw 200 values from the original dataset
    PC1_impacts = sample(df_LC_pf$PC1_impacts, 1000, replace=T),  # Randomly draw 200 values from the original dataset
    PC2_impacts = sample(df_LC_pf$PC2_impacts, 1000, replace=T)   # Randomly draw 200 values from the original dataset
  )
farmer_feed_tech_pred_sim_dry <- ensemble_custom(mLC_farm10, mLC_farm9, mLC_farm11, 
                                                 data=farmer_feed_tech_pred_dry)  # Sample posterior for each of the 1000 randomly-generated respondents
farmer_feed_tech_pred_sim_wet <- ensemble_custom(mLC_farm10, mLC_farm9, mLC_farm11, 
                                                 data=farmer_feed_tech_pred_wet)  # Sample posterior for each of the 1000 randomly-generated respondents
farmer_feed_tech_pred_sim_dry_mu_0 <- as.data.frame(apply(farmer_feed_tech_pred_sim_dry$sim,2,mean))
colnames(farmer_feed_tech_pred_sim_dry_mu_0) <- "Mean"
farmer_feed_tech_pred_sim_dry_mu_0$feed_tech <- 0
farmer_feed_tech_pred_sim_wet_mu_1 <- as.data.frame(apply(farmer_feed_tech_pred_sim_wet$sim,2,mean))
colnames(farmer_feed_tech_pred_sim_wet_mu_1) <- "Mean"
farmer_feed_tech_pred_sim_wet_mu_1$feed_tech <- 1
farmer_feed_tech_pred_sim_mu <- rbind(farmer_feed_tech_pred_sim_dry_mu_0, farmer_feed_tech_pred_sim_wet_mu_1)
head(farmer_feed_tech_pred_sim_mu); dim(farmer_feed_tech_pred_sim_mu)
tapply(farmer_feed_tech_pred_sim_mu$Mean, farmer_feed_tech_pred_sim_mu$feed_tech, mean) # Calculate the difference between farmers (not) affected by FMD


# Produce density plot
farmer_feed_tech_pred_sim_mu2 <- farmer_feed_tech_pred_sim_mu %>%   # Make a mini-df which contains the means in each group
  group_by(feed_tech) %>%                # ... this is used for plotting
  summarize(mean_per_group = mean(Mean))
myCols2 <- c("#009E73","#999999")     # Make colour palette
farmer_feed_tech_pred_sim_mu$feed_tech <- factor(farmer_feed_tech_pred_sim_mu$feed_tech, 
                                                 levels = c(0,1),
                                                 labels = c("Farmer using\ndry feed",
                                                            "Farmer using\nwet feed"))
farmer_feed_tech_pred_sim_mu2$feed_tech <- factor(farmer_feed_tech_pred_sim_mu2$feed_tech, 
                                                 levels = c(0,1),
                                                 labels = c("Farmer using\ndry feed",
                                                            "Farmer using\nwet feed"))
FigX <- 
  ggplot(farmer_feed_tech_pred_sim_mu, aes(Mean, fill = feed_tech)) + 
  geom_density(alpha=0.4) + 
  theme_bw() + 
  xlab("Mean simulated response:\nWould you support the re-legalisation of swill?\n(1 = Definitely not, 5 = Definitely yes)") + 
  ylab("Response density") + 
  theme(axis.title = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12), 
        legend.position=c(0.25,0.85), 
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.key.size = unit(2, 'lines')) + 
  xlim(1,5) + 
  annotate("text", label = "a", x = 1, y = .7, size = 8,fontface =2) +
  geom_vline(data = farmer_feed_tech_pred_sim_mu2, aes(xintercept = mean_per_group,  # add in vertica lines at mean
                                           col = feed_tech), size=1, linetype="dashed") +
  scale_fill_manual(values = myCols) +
  scale_color_manual(values = myCols) 
FigX
  

# tidy up after plotting
rm(farmer_feed_tech_pred_sim_dry, farmer_feed_tech_pred_sim_wet,  # remove items used for plotting
   farmer_feed_tech_pred_sim_wet_mu_1, farmer_feed_tech_pred_sim_dry_mu_0, farmer_feed_tech_pred_sim_mu,
   farmer_feed_tech_pred_sim_mu2, farmer_feed_tech_pred_dry, farmer_feed_tech_pred_wet, FigX, myCols2)


# tidy up 
# .... clean working environment of all loaded models
rm(pr_LC_farm, cum_pr_LC_farm)                           # remove starting values
rm(mLC_farm, mLC_farm2, mLC_farm3, mLC_farm4,            # remove models
   mLC_farm5, mLC_farm6, mLC_farm7, mLC_farm8, 
   mLC_farm9, mLC_farm10, mLC_farm11, mLC_farm12, 
   mLC_farm13, mLC_farm14)


# 13e) Model of willingness to use swill on farm ----


# Fit model for willingness to use swill on farm
# Select starting values
pr_WU_farm <- table( df_LC_pf$will_use_swill ) / nrow(df_LC_pf)
cum_pr_WU_farm <- cumsum( pr_WU_farm )        # cumsum converts to cumulative proportions
logit( cum_pr_WU_farm )                       # the logit intercept parameter estimate (averaged across all feeds - this is just used for starting values)


# Model WU1 ("Willingess to Use")
# Use long chain to overcome divergent iterations
mWU <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aA[AgeClass] + aF[farm_size] +         # phi is my linear model
      bG*gender + bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aA[AgeClass]  ~ dnorm(0, sigma_A),              # prior for AgeClass varying intercept
    aF[farm_size]  ~ dnorm(0, sigma_F),             # prior for farm_size varying intercept
    c(bG, bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma_A ~ dcauchy(0,1),                         # priors for AgeClass variance
    sigma_F ~ dcauchy(0,1),                         # priors for farm_size variance
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Plot
plot(mWU)
precis(mWU, depth=2)               # Rhat == 1
logistic(coef(mWU)[c(1:4)])        # convert the cutpoints into cumulative probabilities


# Model WU2
# Model without farm size
mWU2 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aA[AgeClass] +          # phi is my linear model
      bG*gender + bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aA[AgeClass]  ~ dnorm(0, sigma_A),              # prior for AgeClass varying intercept
    c(bG, bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma_A ~ dcauchy(0,1),                         # priors for AgeClass variance
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Plot
plot(mWU2)
precis(mWU2, depth=2)               # Rhat == 1


# Model WU3
# model without AgeClass
mWU3 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,   # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <- aF[farm_size] +         # phi is my linear model
      bG*gender + bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    aF[farm_size]  ~ dnorm(0, sigma_F),             # prior for farm_size varying intercept
    c(bG, bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    sigma_F ~ dcauchy(0,1),                         # priors for farm_size variance
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Plot
plot(mWU3)
precis(mWU3, depth=2)               # All Rhat == 1


# Model WU4
# model without AgeClass & farm_size
mWU4 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU4)
precis(mWU4, depth=2)               # All Rhat == 1


# Model WU5
# model without AgeClass, farm_size, gender
mWU5 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bFMD*FMD_exp + bS*swill_exp + bFT*feed_tech +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bS, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU5)
precis(mWU5, depth=2)               # All Rhat == 1


# Model WU6
# model without AgeClass, farm_size, swill_exp
mWU6 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bFMD*FMD_exp + bFT*feed_tech +  
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU6)
precis(mWU6, depth=2)               # All Rhat == 1


# Model WU7
# model without AgeClass, farm_size, feed_tech
mWU7 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bFMD*FMD_exp +  bG*gender + bS*swill_exp + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bFMD, bS, bG, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU7)
precis(mWU7, depth=2)               # All Rhat == 1


# Model WU8
# model without AgeClass, farm_size, FMD_exp
mWU8 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bS*swill_exp + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bS, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU8)
precis(mWU8, depth=2)               # All Rhat == 1


# compare first few models fit
# ... models 8,7,5 get similar weighting
# ... low weight to models with age class or farm size
compare(mWU8, mWU7, mWU6, mWU5, mWU4, mWU3, mWU2, mWU)
plot(coeftab(mWU8, mWU7, mWU5))


# Model WU9
# model without AgeClass, farm_size, FMD_exp, swill_exp
mWU9 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU9)
precis(mWU9, depth=2)               # All Rhat == 1


# Model WU10
# model without AgeClass, farm_size, gender, swill_exp
mWU10 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bFMD*FMD_exp + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bFMD, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU10)
precis(mWU10, depth=2)               # All Rhat == 1


# Model WU11
# model without AgeClass, farm_size, gender, FMD_exp
mWU11 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bS*swill_exp + bFT*feed_tech + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bS, bFT, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU11)
precis(mWU11, depth=2)               # All Rhat == 1


# Model WU12
# model without AgeClass, farm_size, gender, feed_tech
mWU12 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bS*swill_exp + bFMD*FMD_exp + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bS, bFMD, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU12)
precis(mWU12, depth=2)               # All Rhat == 1


# Model WU13
# model without AgeClass, farm_size, swill_exp, feed_tech
mWU13 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bFMD*FMD_exp + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bFMD, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU13)
precis(mWU13, depth=2)               # All Rhat == 1


# Model WU14
# model without AgeClass, farm_size, FMD_exp, feed_tech
mWU14 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + bS*swill_exp + 
      bPC1v*PC1_values + bPC2v*PC2_values + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG, bS, bPC1v, bPC2v, bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU14)
precis(mWU14, depth=2)               # All Rhat == 1


# Model WU15
# model with only gender & perceptions
mWU15 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bG*gender + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bG,  bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU15)
precis(mWU15, depth=2)               # All Rhat == 1


# Model WU16
# model with only FMD experience & perceptions
mWU16 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bFMD*FMD_exp + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bFMD,  bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU16)
precis(mWU16, depth=2)               # All Rhat == 1


# Model WU17
# model with only swill_exp & perceptions
# i.e. model without AgeClass, farm_size, FMD_exp, feed_tech, gender, values
mWU17 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bS*swill_exp + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bS,  bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)
  

# Check chains
plot(mWU17)
precis(mWU17, depth=2)               # All Rhat == 1


# Model WU18
# model with only swill_exp & perceptions
# i.e. model without AgeClass, farm_size, FMD_exp, feed_tech, gender, values
mWU18 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                           # phi is my linear model
      bFT*feed_tech + 
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bFT,  bPC1i, bPC2i) ~ dnorm(0,10),             # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU18)
precis(mWU18, depth=2)               # All Rhat == 1


# Model WU19
# model without AgeClass, farm_size, FMD_exp, feed_tech, gender, values
mWU19 <- map2stan(
  alist(
    will_use_swill ~ dordlogit(phi, cutpoints) ,     # NB that map2stan prefers an explicit vector of intercept parameters, hence create cutpoints
    phi <-                                          # phi is my linear model
      bPC1i*PC1_impacts + bPC2i*PC2_impacts,
    c(bPC1i, bPC2i) ~ dnorm(0,10),  # priors for slopes
    cutpoints ~ dnorm(0,10)                         # prior for thresholds of cumulative likelihood
  ) ,
  data=df_LC_pf ,
  start=list(cutpoints = c(logit(cum_pr_WU_farm)[[1]], logit(cum_pr_WU_farm)[[2]],     # base starting values on overall mean
                           logit(cum_pr_WU_farm)[[3]], logit(cum_pr_WU_farm)[[4]]) ),
  iter = 10000, chains=1, cores = 1)


# Check chains
plot(mWU19)
precis(mWU19, depth=2)               # All Rhat == 1


# compare models
# ... models 41 and 42 have the highest weights
# ... model 41 has 0.54 of the weight
compare(mWU, mWU2, mWU3, mWU4, 
        mWU5, mWU6, mWU7, mWU8, 
        mWU9, mWU10, mWU11, mWU12, 
        mWU13, mWU14, mWU15, mWU16,
        mWU17, mWU18, mWU19)


# ----------------------- 13f) Plot Fig 6, willingess to use swill ----------------------- 


# Plot the coefficients of the two models with greatest weights 
# ... (cumlatively 60% of model weights)
# ... Step 1) extract model output
mWU17_res <- precis(mWU17, depth=2, pars=c("bS", "bPC1i", "bPC2i"))@output   # extract model output
mWU17_res$coef <- row.names(mWU17_res)    # Extract the coefficient names
mWU17_res$model <- "model WU17"
mWU17_res$weighting <- 0.42               # this comes from the compare() function called above
mWU19_res <- precis(mWU19, depth=2, pars=c("bPC1i", "bPC2i"))@output   # extract model output
mWU19_res$coef <- row.names(mWU19_res)    # Extract the coefficient names
mWU19_res$model <- "model WU19"
mWU19_res$weighting <- 0.18               # this comes from the compare() function called above
mWU_res <- rbind(mWU17_res,mWU19_res)
mWU_res$coef <- factor(mWU_res$coef, 
                         levels=c("bS", "bPC2i","bPC1i"), 
                         labels=c("Have previously used swill on their farm", 
                                  "Factor 2, perception of swill:\nSwill would increase disease risk\n& be unpalatable to consumers",
                                  "Factor 1, perception of swill:\nSwill would benefit the environment, help\nfarms financially, & reduce trade-deficits"))
mWU_res$model <- factor(mWU_res$model, 
                        levels= c("model WU17", "model WU19"),
                        labels= c("WU17","WU19"))
colnames(mWU_res) <- c("Mean","StdDev", "lwr.89", "upr.89", "n_eff", "Rhat", "coef", "model", "weighting")


# save the Fig6 model output
# write.table(mWU_res, "results_tables/Fig6_output.csv", sep=",", row.names=F)
# mWU_res <- read.csv("results_tables/Fig6_output.csv", header=T)


# Relabel for ceofficients  type 
mWU_res$coef_type <- "Respondent & farm\ncharacteristics"
mWU_res[mWU_res$coef %in% c("Factor 1, respondent's values:\nRespondents placing high importance on\ndisease control & perception of the industry",
                                                "Factor 2, respondent's values:\nRespondents placing high importance on\nresource efficiency & farm finances"),"coef_type"] <- "Factors:\nRespondent values"
mWU_res[mWU_res$coef %in% c("Factor 1, perception of swill:\nSwill would benefit the environment, help\nfarms financially, & reduce trade-deficits",
                                                "Factor 2, perception of swill:\nSwill would increase disease risk\n& be unpalatable to consumers"),"coef_type"] <- "Factors:\nPerception of swill"
head(mWU_res)


# Shorten coefficient titles
mWU_res$coef <- factor(mWU_res$coef,
                                 levels=c("Have previously used swill on their farm",
                                          "Factor 2, perception of swill:\nSwill would increase disease risk\n& be unpalatable to consumers",
                                          "Factor 1, perception of swill:\nSwill would benefit the environment, help\nfarms financially, & reduce trade-deficits"),
                                 labels=c("Have previously used swill on their farm",
                                          "Swill would increase disease risk &\nbe unpalatable to consumers",
                                          "Swill would benefit the environment, help\nfarms financially, & reduce trade-deficits"))


# Step 2) Plot with sizes legend
myCols2 <- c("#0072B2", "#D55E00")  # Make colour palette for 2 colours
Fig6a <- ggplot(mWU_res, aes(Mean, coef, col = model,size = weighting*100)) + 
  stat_sum(position = position_dodgev(height = 0.5)) +
  geom_errorbarh(aes(xmin=lwr.89, xmax=upr.89, col = model),
                 size=1,    # Thinner lines
                 width=.3,
                 position=position_dodgev(height = 0.5)
  )  +
  facet_grid(coef_type ~., scales= "free", switch = "y") +
  # coord_flip() +
  geom_vline(xintercept = 0) +
  xlab("Coefficient value") +
  theme_bw()  +
  theme(legend.position = c(0.6,0.55),
        legend.title = element_text(face="bold"),
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", colour = NA),
        axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12),
        axis.line = element_line(colour="black")) +
  scale_color_manual(values = myCols2, 
                     name="Model ID",
                     guide=FALSE) +     # removes the colour legend
  scale_size_continuous(limits = c(8,42),              # standardise the weights displayed in the legend
                name="Model\nweight (%)")   # specify name for size legend

  
# Step 3) Plot with colour legend
Fig6b <- ggplot(mWU_res, aes(Mean, coef,col = model,size = weighting*100)) + 
  stat_sum(position = position_dodgev(height = 0.5)) +
  geom_errorbarh(aes(xmin=lwr.89, xmax=upr.89, col = model),
                 size=1,    # Thinner lines
                 width=.3,
                 position=position_dodgev(height = 0.5)
  )  +
  facet_grid(coef_type ~., scales= "free", switch = "y") +
  # coord_flip() +
  geom_vline(xintercept = 0) +
  xlab("Coefficient value") +
  theme_bw()  +
  theme(legend.position = c(0.6,0.55),
        legend.title = element_text(face="bold"),
        legend.background = element_rect(fill = "grey95"),
        legend.key = element_rect(fill = "grey95"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white", colour = NA),
        axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12),
        axis.line = element_line(colour="black")) +
  scale_color_manual(values = myCols2, 
                     name="Model ID") + 
  scale_size_continuous(limits = c(8,42),              # standardise the weights displayed in the legend
                        name="Model\nweight (%)",
                        guide = FALSE)  # specify name for size legend
Fig6b  


# Added in second legend, title and arrows manually to the plot


# # Step 4) Extract the colour legend - leg
# leg <- gtable_filter(ggplot_gtable(ggplot_build(Fig6b)), "guide-box") 
# 
# 
# # Step 5) Plot and add second legend
# Fig6 <- Fig6a + 
#   annotation_custom(grob = leg, 
#                     ymin = 0, ymax = 3,    # Place legend
#                     xmin = 1.4, xmax = 1.5)      # Nb because of coord_flip() the x and y axes are flipped
# grid.newpage()
# grid.draw(Fig6)
# 
# 
# # Step 6) Add in titles: 
# # ... code from: https://stackoverflow.com/questions/21997715/add-ggplot-annotation-outside-the-panel-or-two-titles
# g = ggplotGrob(Fig6)
# g = gtable_add_grob(g, grobTree(textGrob("Lower willingness to use swill", x=0, hjust=0),
#                                 textGrob("Greater willingness to use swill", x=1, hjust=1)),
#                     t=1, l=4)
# grid.draw(g)


# tidy up
rm(g, Fig6a, Fig6b, Fig6, mWU17_res, mWU19_res, mWU_res, myCols2, leg)



# --------- 13g) Plot S2 Appendix Figure 12, effect of experience feeding swill on willingness to use --------- 


# Calculate effect size of having previous experience of using swill
# ... Use output from model mWU17 (model with highest weight)
WU_pred <- 
  data.frame(
    will_use_swill=5,                          #### needed so sim knows num levels
    swill_exp = rep(c(0,1), each=1000),          # 100 zeros, 100 1s.
    PC1_impacts = sample(df_LC_pf$PC1_impacts, 2000, replace=T),  # Randomly draw 200 values from the original dataset
    PC2_impacts = sample(df_LC_pf$PC2_impacts, 2000, replace=T)   # Randomly draw 200 values from the original dataset
  )
WU_pred_sim <- sim(mWU17,data=WU_pred)  # Sample posterior for each of the 1000 randomly-generated respondents
str(WU_pred_sim); class(WU_pred_sim)
WU_pred_sim_mu_0 <- as.data.frame(apply(WU_pred_sim[,c(1:1000)],2,mean))
colnames(WU_pred_sim_mu_0) <- "Mean"
WU_pred_sim_mu_0$swill_exp <- 0
WU_pred_sim_mu_1 <- as.data.frame(apply(WU_pred_sim[,c(1001:2000)],2,mean))
colnames(WU_pred_sim_mu_1) <- "Mean"
WU_pred_sim_mu_1$swill_exp <- 1
WU_pred_sim_mu <- rbind(WU_pred_sim_mu_0, WU_pred_sim_mu_1)
head(WU_pred_sim_mu); dim(WU_pred_sim_mu)
tapply(WU_pred_sim_mu$Mean, WU_pred_sim_mu$swill_exp, mean) # Calculate the difference between farmers (not) affected by FMD


# plot
WU_pred_sim_mu$swill_exp <- factor(WU_pred_sim_mu$swill_exp, 
                                   levels = c(0,1),
                                   labels = c("Have never used swill\non their farm",
                                              "Have previously used swill\non their farm"))


# plot
compare_swill <- WU_pred_sim_mu %>%   # Make a mini-df which contains the means in each group
  group_by(swill_exp) %>%                # ... this is used for plotting
  summarize(mean_per_group = mean(Mean))
myCols2 <- c("#009E73","#999999")     # Make colour palette
S2App_Fig12<- ggplot(WU_pred_sim_mu, aes(Mean, fill=swill_exp)) + 
  geom_density(alpha=0.4) + 
  theme_bw() + 
  xlim(1,5) + 
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=12),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'lines'),
        legend.text = element_text(size=12) ) + 
  ylab("Response density") + 
  xlab("Mean simulated response:\nWould you consider using swill on your farm?\n(1 = Definitely not, 5 = Definitely yes)") + 
  geom_vline(data = compare_swill, aes(xintercept = mean_per_group,  # add in vertica lines at mean
                                           col = swill_exp), size=1, linetype="dashed") +
  scale_fill_manual(values = myCols2) +
  scale_color_manual(values = myCols2)
S2App_Fig12


# tidy up
rm(WU_pred, WU_pred_sim, S2App_Fig12, compare_swill, myCols2,
   WU_pred_sim_mu_0, WU_pred_sim_mu_1, WU_pred_sim_mu)
rm(mWU,mWU2,mWU3,mWU4,mWU5,mWU6,mWU7,mWU8,mWU9, mWU10, mWU11, mWU12, mWU13, mWU14, mWU15, mWU16, mWU17, mWU18, mWU19)
rm(pr_WU_farm, cum_pr_WU_farm)
rm(df_LC_pf, no.qu.df_copy)


# Clean the global environment
rm(no.qu.df, qu.df, multiplot, ensemble_custom)


# End of script