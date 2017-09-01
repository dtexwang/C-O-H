# Provides functions used in CRUSTALFLUIDMODEL:
#   logactfinder(itarget, irxn, icond, conditions)
#
# Call functions from this file using source("solvers.R").  Some variables 
# from CRUSTALFLUIDMODEL are used (i.e., must be global scope).


# TIC/TOC ----------------------------------------------------------------- 
# Timer function for profiling.  http://stackoverflow.com/a/30392068
# 
# > toc <- tic() # start
# > toc()        # end

tic <- function () { 
  now <- proc.time()
  function () { 
    proc.time() - now 
  }
}


# LOGACTFINDER ------------------------------------------------------------

# Finds logact of a species (index itarget)
# given T, P, and logfO2 (in conditions variable).

logactfinder <- function(itarget, irxn, icond, conditions)
{
  cat(c("Finding the equilibrium activity (fugacity) of:", "\n"))
  cat(c("   -----------------", "\n"))
  cat(c("  ", species[itarget], states[itarget], "\n"))
  cat(c("   -----------------", "\n"))
  
  cat(c("        T =", conditions$T[icond], "C", "\n"))
  cat(c("        P =", conditions$P[icond], "bar", "\n"))
  cat(c("   logfO2 =", conditions$logfO2[icond], "\n"))
  cat(c("   -----------------", "\n"))
  
  iO2 = which(species == "O2")
  
  modelfun <- function(logain = 0)
  {
    logact_guess <- logact_guess_default
    logact_guess[itarget] <- logain
    logact_guess[iO2] <- conditions$logfO2[icond]
    
    capture.output(
      result <- subcrt(
        species,
        states,
        as.matrix(rxndef)[irxn,],
        T = conditions$T[icond],
        P = conditions$P[icond],
        logact = logact_guess
      )
    )
    
    stopifnot(dim(result$out[,])[1] == 1)         # check that only one condition was tested
    
    rootfn = result$out$logQ - result$out$logK    # function to zero
    # print(rootfn)
    return(rootfn)
  }
  
  logact_searchintvl = c(-40, 10)   # OPT -10,5 to -40,10
  rt <- uniroot(modelfun, logact_searchintvl, extendInt = "no")  # OPT yes -> no
  
  cat(c("   logact =", round(rt$root, 3), "\n"))
  cat(c("   -----------------", "\n"))
  return(rt)
}


# EQUILIBRUMFINDER --------------------------------------------------------

# Finds logact of all other gaseous species to match the activities already
# solved for.

equilibriumfinder <- function(itgts, irxn, icond, conditions)
{
  cat(c("Finding the equilibrium activities (fugacity) of:", "\n"))
  cat(c("   -----------------", "\n"))
  for (ii in 1:length(itgts)) {
    cat(c("  ", sprintf("%-8s",species[itgts[ii]]), states[itgts[ii]], "\n"))
  }
  cat(c("   -----------------", "\n"))
  
  cat(c("        T =", conditions$T[icond], "C", "\n"))
  cat(c("        P =", conditions$P[icond], "bar", "\n"))
  cat(c("   logfO2 =", conditions$logfO2[icond], "\n"))
  cat(c("   -----------------", "\n"))
  
  ifixed = which((species %in% c("graphite", "O2", "CO2", "CO")))
  
  fixedacts <- logaeq[icond,ifixed]
  
  modelfun <- function(logain)
  {
    logact_guess <- logact_guess_default
    logact_guess[itgts] <- logain
    logact_guess[ifixed] <- fixedacts  # OPT replace with nonsubsetting logaeq # set O2, CO2, CO activities to those already solved
    
    rootfn <- array(numeric(), dim = c(length(irxn)+1))  # the +1 is for Ptot
    
    for (jj in 1:length(irxn)) {
    capture.output(
      result <- subcrt(
        species,
        states,
        as.matrix(rxndef)[irxn[jj],],
        T = conditions$T[icond],
        P = conditions$P[icond],
        logact = logact_guess
      )
    )
    
    # stopifnot(dim(result$out[,])[1] == 1)   # OPT comment          # check that only one condition was tested
    
    rootfn[jj] = result$out$logQ - result$out$logK    # function to zero
    }
    
    loggases = logact_guess[igases]
    Ptot = sum(10^loggases)
    
    rootfn[length(rootfn)] = log10(Ptot) - log10(conditions$P[icond])
    
    return(rootfn)
  }
  guessin = array(0, c(length(itgts)))
  rts <- multiroot(modelfun, guessin, atol = 1e-3)  # OPT increase atol from 1e-8 default to 1e-3
  
  cat(c("   logact =", "\n"))
  for (ii in 1:length(itgts)) {
    cat(c("    ", sprintf("%-4s",species[itgts[ii]]), "=", round(rts$root[ii],3), "\n"))
  }
  cat(c("   -----------------", "\n"))
  return(rts)
}
