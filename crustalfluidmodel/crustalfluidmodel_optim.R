# CRUSTALFLUIDMODEL  Calculates the composition of a graphite-saturated C-O-H 
# vapor at given fO2, T, and P conditions.
# 
# Calculations are similar to those of French (1966), Kawagucci et al. (2013), 
# and McDermott (2015), but this model differs by considering graphite and C1-C3
# species together.
# 
# Required packages: CHNOSZ_1.0.8, rootSolve_1.7
# 
# Last modified 19 March 2017 by D. T. Wang (dtw@mit.edu) R version 3.3.3 
# (x86_64) on Windows 7

# Try to optimize:
# Base case = 28 sec for 3 conditions
#             27 sec when expand srch intrvl and extendInt -> NO
#             23 sec when omit stopifnot conditino check
#             23 sec when omit not subsetting outer array logaeq
#   Profiler says spend 50% time in stode().Call (i.e., calling Fortran),
#      and 28% of time in cgl(), the eos solver in CHNOSZ.
#             20 sec when increase atol
# Optimized = 


# Startup tasks -----------------------------------------------------------

library(CHNOSZ)
library(rootSolve)

# Set up CHNOSZ -----------------------------------------------------------

data(thermo)
add.obigt()

T.units("C")    # default = K
E.units("J")    # default = cal
P.units("bar")  # default = bar


# Definitions -------------------------------------------------------------

# Species
# ..................
#  1  graphite  cr
#  2  CO        g
#  3  CO2       g
#  4  CH4       g
#  5  H2        g
#  6  H2O       g
#  7  O2        g
#  8  ethane    g
#  9  propane   g

species = matrix(c("graphite", "CO", "CO2", "CH4", "H2", "H2O", "O2", "ethane", "propane"), 
                 ncol = 1, byrow = TRUE)
states = matrix(c("cr", rep("g", times=8)), ncol = 1, byrow = TRUE)

#species = species[-c(8,9)]
#states = states[-c(8,9)]

# Reactions
# ..................
# (1) C + O2 = CO2
# (2) C + 0.5 O2 = CO
# (3) CO2 + 4 H2 = CH4 + 2 H2O
# (4) H2 + 0.5 O2 = H2O
# (5) 2 CH4 = C2H6 + H2
# (6) CH4 + C2H6 = C3H8 + H2

rxndef <- read.table(
  "ReactionDefinition.csv",
  header = TRUE,
  skip = 1,
  row.names = NULL,
  sep = ","
)

rxndef[is.na(rxndef)] <- 0  # replace blanks (NA) with zeros
rxndef <- rxndef[ , -c(1)]  # drop first column (rxn label)

#rxndef = rxndef[, -c(8,9)]

# Add'l Constraints
# ..................
# [1] f_CO + f_CO2 + ... + f_propane = P


# Define T, P, fO2 conditions ---------------------------------------------

T = seq(300,700,by=200)    # Celsius
P = 500                  # bar
logfO2 = seq(-41,-15,by=30)

conds <- expand.grid(T=T, P=P, logfO2=logfO2)
# source("calcbufferfO2.R") # calculate/save to disk logfO2 buffered by mineral assemblages

T = seq(300,600,by=10)    # Celsius
T = c(0, 50, 100, 150, 200, 225, 250, 275, seq(300,500, by=10), 525, 550, 575, 600, 650, 700)
P = 1000                  # bar
# logfO2 = seq(-41,-15,by=30)

conds <- expand.grid(T=T, P=P, logfO2=logfO2)
source("calcbufferfO2.R") # calculate/save to disk logfO2 buffered by mineral assemblages

conds <- cbind(TP, buffO2$FMQ)    # SAVETIME, just calculate FMQ : T, P, fO2
colnames(conds) <- c("T", "P", "logfO2")

logact_guess_default = rep(c(0), times=length(species))


# Calculate speciation of CCO system and check for graphite stability -----

source("solvers_optim.R")

toc <- tic()  # could also use command from pkg "tictoc"

logaeq <- array(numeric(), c(dim(conds)[1],length(species)))
logaeq[,which(species == "O2")] <- conds$logfO2

logaeq[,which(species == "graphite")] <- 0   # assign loga graphite = 0

for (icond in 1:dim(conds)[1]) {
  itarget = which(species == "CO2")
  
  soln <-
    logactfinder(
      itarget = itarget,
      irxn = 1,
      icond = icond,
      conditions = conds
    )$root
  
  logaeq[icond, itarget] <- soln
}  

for (icond in 1:dim(conds)[1]) {
  itarget = which(species == "CO")
  
  soln <-
    logactfinder(
      itarget = itarget,
      irxn = 2,
      icond = icond,
      conditions = conds
    )$root
  
  logaeq[icond, itarget] <- soln
}

# identify which conditions for graphite is unstable
igases = which(!(states %in% c("liq", "cr", "cr1", "cr2")))
badconds = which(rowSums(10^logaeq[,igases], na.rm = TRUE) > conds$P)   

# conds <- conds[-badconds,]
# logaeq<- logaeq[-badconds,]

for (icond in 1:dim(conds)[1]) {
  itargets = which( !(species %in% c("graphite", "O2", "CO2", "CO")))
  
  if (icond %in% badconds) {   # don't try to solve if Ptot > P, just fill with NaN's
    soln <- NaN
  } else if (!(icond %in% badconds)) {
    soln <-
      equilibriumfinder(
        itgts = itargets,
        irxn = 3:dim(rxndef)[1],
        icond = icond,
        conditions = conds
      )$root
  }
  
  logaeq[icond, itargets] <- soln
  
  loggases = logaeq[icond,igases]
  Ptot = sum(10^loggases, na.rm = TRUE)
  
  cat(c("  logPtot =", round(log10(Ptot),3), "\n"))
  cat(c("   -----------------", "\n"))
}

colnames(logaeq) <- apply(cbind(species,states), 1, paste, collapse = ",")

write.csv(conds, "conds2_FMQ.csv")
write.csv(logaeq,"logaeq2_FMQ.csv", na = "NaN")

toc()

# 20 Mar 2017,  Calculated at 17 temperature and 8 logfO2 intervals (136
# conditions). These calculations took 26 minutes on a Intel Core i7-6600U CPU @
# 2.60GHz with 8GB RAM.