######################################################################################################
# helsim01 - R version of the individual-based model in Medley 1989 thesis, and Anderson&Medley 1985
# Original version Graham Medley, July 2013
#
# Results format is a list containing the parameters ($params, last list element) 
# and results from individual realisations (list elements 1,2,3...)
#
# The results from each individual realization is a list. 
# Its elements are:
# lists for each time point specified by the outTimings parameter. 
# $repNo - the realization number (don't have a use for this at the moment). 
#
# The list for each time point in each realization contains: 
# worms: a data.frame of total and female worms at that time for each host.
# hosts: a data.frame with birth dates and death dates for each host
# freeLiving: Reservoir intensity
# time: time at which the data stored in this list element was outputed
# 
#######################################################################################################



# Remove objects
rm(list=ls())

closeAllConnections()
runtime = proc.time()


############################################################################
# WHEN RUNNING SIMULATIONS ON CLUSTER USE THIS BLOCK TO READ IN 
# PARAMETERS FROM COMMAND LINE/BATCH FILE
############################################################################

## command line arguments. 
args <- commandArgs(trailingOnly = TRUE)

outFileStub <- args[1]				# stub for output file name
dataInputPath <- args[2]			# path to parameter files
outPath <- args[3]				# path to where data should be outputted
paramFileName <- args[4]			# name of file with simulation parameters
currentDemogName <- args[5]			# name of demographies used in simulation
seed <- as.numeric(args[6])			# seed for repeated runs

delay <- as.numeric(args[7]) / 10		
print(paste("delay:", delay))


############################################################################
# WHEN RUNNING SIMULATIONS ON DESKTOP USE THIS BLOCK TO SET PARAMETERS 
############################################################################

#outFileStub <- ""
#dataInputPath <- "R:\\ModellingOutput\\NTDmc_STH\\STH-Covid2019\\Ascaris\\Variable-delay\\" 	# replace by your own path
#outPath <- "R:\\ModellingOutput\\NTDmc_STH\\STH-Covid2019\\Ascaris\\Variable-delay\\"		# replace by your own path
#paramFileName <- "AscarisParameters.txt"	
#currentDemogName <- "KenyaKDHS"


############################################################################
# SET WORKING DIRECTORY AND PARAMETER FILE
############################################################################

setwd(dataInputPath)						# function files etc should be in input directory
parameterPath <- paste0(dataInputPath, paramFileName)	# path plus filename of 


############################################################################
# SET UP INITIAL CONDITIONS FOR SIMULATION
############################################################################


# Source functions for stochastic code
source(paste0(dataInputPath, "ParasiteFunctions.R"))
source(paste0(dataInputPath, "helsim_FUNC.R"))


# Set up the output and recording
#logFile <- file(paste0(outPath, outFileStub, "_log.txt"), "w")
#cat(file=logFile,"Log file is open!\n")


# Read in parameters
testExist <- file(parameterPath,"r")
if(!isOpen(testExist)){
	#cat(file=logFile,"Parameter file not found: ",parameterPath,"\n")
	stop("Parameter file not found.")
} else {
	#cat(file=logFile,"Parameter file found: ",parameterPath,"\n")
}
close(testExist)
params <- readParams(fileName=parameterPath, demogName=currentDemogName)


# Configure parameters
params <- configure(params) 


###################################################################
# WHEN RUN ON CLUSTER THESE ARE SET VIA COMMAND LINE (BATCH FILE)
# OTHERWISE SET DIRECTLY FROM PARAMETER FILE

params$delay <- delay

###################################################################

# update chemotimings from default after setting delay via batch file
params$chemoTimings <- updateChemoTimings(params)

params$psi <- getPsi(params)

# Call function to return equilibrium reservoir level for IC. 
equi <- getEquilibrium(params)
params$equiData <- equi


# DEFINE INFECTION INTENSITY THRESHOLDS AND CONVERT TO EGG COUNTS
# the inensity thresholds should really be set in an external configuration file because they vary by species
# ASCARIS: moderate intensity: 5000, high intensity: 50000
# HOOKWORM: moderate intensity: 2000, high intensity: 4000
# TRICHURIS: moderate intensity: 1000, high intensity: 10000
# 24 is a conversion factor from epg (eggs per gram faeces) to egg counts on a microscopy slide of a faecal smear test

moderateIntensityCount <- 5000 / 24
highIntensityCount <- 50000 / 24


# Specify a single realisation/rep as a function
doRealization <- function(params)
{
	# Setup the simulation
	simData <- setupSD(params)	# Set up simulation data for village	

	time <- 0										# Set initial time to 0
	FLlast <- time									# Time at which to update simulation data using the doFreeLive function
	outTimes <- params$outTimings							# Times at which data should be recorded	
	chemoTimes <- params$chemoTimings						# Times at which chemotherapy is given

	nextOutIndex <- which.min(outTimes)						# Determine time when data should be recorded next
	nextOutTime <- outTimes[nextOutIndex]
	ageingInt <- 1/52     								# Check ages every week - should this be in a parameter file?
	nextAgeTime <- ageingInt							# Time at which individuals' age is advanced next
	nextChemoIndex <- which.min(chemoTimes)					# Determin time at which individuals receive next chemotherapy	 
	nextChemoTime <- chemoTimes[nextChemoIndex]				 	

	nextStep <- min(nextOutTime, time+params$maxStep, nextAgeTime, nextChemoTime)	# Determine which event will take place next	

	results <- list() 								# initialise empty list to store results 	
	outCount <- 1									# keep count of outputs

	
	# Run stochastic algorithm	
	while((time<params$maxTime))    
	{
		
		rates <- calcRates(params, simData, time)
		sumRates <- sum(rates)
    
		if(sumRates > 0.001) 
		{
			tstep <- rexp(1, sumRates)	# draw time interval from exponential distribution
		}else
		{
			tstep <- params$maxTime
		}
				
		if((time+tstep)<nextStep) 
		{
			time <- time + tstep			
			simData <- doEvent(rates, params, simData, time)

		}else 
		{ 
			time <- nextStep

		  	## Predetermined event block.
			simData <- doFreeLive(params, simData, nextStep-FLlast)
			FLlast <- nextStep			 
			timeBarrier <- nextStep + 0.001

			# Ageing and death 
			if(timeBarrier>nextAgeTime)
			{ 
		  		simData <- doDeath(params, simData, time) 
				nextAgeTime <- nextAgeTime + ageingInt
			}

			
			# Chemotherapy
			if(timeBarrier>nextChemoTime)
			{
				simData <- doDeath(params, simData, time)				# update age groups/deaths just before chemotherapy
				before <- simData$worms
				simData <- doChemo(params, simData, time)	
				after <- simData$worms		
				chemoTimes[nextChemoIndex] <- params$maxTime + 10		# increase current chemo time beyond maxTime and
				nextChemoIndex <- which.min(chemoTimes)				# choose nextChemoTime
				nextChemoTime <- chemoTimes[nextChemoIndex]		
			}
			
			
	
			# Record
			if(timeBarrier>nextOutTime )  
			{ 	
				currentRes <- list()
				currentRes$worms <- simData$worms
				currentRes$hosts <- simData$demography
				currentRes$freeLiving <- simData$freeLiving
				currentRes$time <- time
				currentRes$prevKKSAC <- getAgeCatSampledPrevByVillage(simData, time, 1, TRUE, 1.0, c(5,14))
				currentRes$prevKKAll <- getAgeCatSampledPrevByVillage(simData, time, 1, TRUE, 1.0, c(0, 100))

				currentRes$prevTrue <- villageTruePrev(simData, TRUE, 1, params)
				currentRes$prevTrue2 <- villageTruePrev2(simData)

				currentRes$meanIntensitySAC <- getMeanInfectionIntensity(simData, time, 1, TRUE, 1.0, c(5,14))
				currentRes$meanIntensityAll <- getMeanInfectionIntensity(simData, time, 1, TRUE, 1.0, c(0, 100))

				currentRes$prevHISAC <- getMediumHeavyPrevalenceByVillage(simData, time, 1, TRUE, 1.0, c(5,14), highIntensityCount)		
				currentRes$prevMHISAC <- getMediumHeavyPrevalenceByVillage(simData, time, 1, TRUE, 1.0, c(5,14), moderateIntensityCount)
				

				results[[outCount]] <- currentRes
				outCount <- outCount+1
				
				outTimes[nextOutIndex] <- params$maxTime+10
				nextOutIndex <- which.min(outTimes)
				nextOutTime <- outTimes[nextOutIndex]
			}		
			
			nextStep <- min(nextOutTime, time+params$maxStep, nextAgeTime)
			
		} ## end of predetermined event block.
	} ## end of while loop. 
	
	
	finalResults <- list(results=results)
	
	return(finalResults)
}


#cat(file=logFile, "Setup complete after ", (proc.time() - runtime)[["elapsed"]], " secs.\n")


## Set up and run multiple instances of doRealization() in parallel

#cat(file=logFile, "Setting up parallel stuff.\n")

# if possible the loops should be implemented so that all iterations are reproducible
require("foreach")
require("doParallel")
require("doRNG")

cores <- detectCores()
cl <- makeCluster(cores)
#cl <- makeCluster(params$nNodes)
registerDoParallel(cl)
#registerDoRNG(seed=params$seed)

#cat(file=logFile,"Starting parallel stuff.\n")

set.seed(params$seed)
#foreachResults <- foreach(i=1:5) %dopar% doRealization(params)
foreachResults <- foreach(i=1:params$numReps) %dorng% doRealization(params)

stopImplicitCluster() # On advice from Wes. 

foreachResults$params <- params


#cat(file=logFile,"Total elapsed time = ",(proc.time() - runtime)[["elapsed"]],"\n")


resultsFile <- paste0(outPath, outFileStub, "results_as_R0_", params$R0, "_k", params$k, "_delay", params$delay, "_seed", params$seed, ".RData")

# we store an R data file with a tree of lists that represents the entire state of the simulation
# in the last time step
# how you can extract a data frame with prevalences from the results is shown in the files
# test_distributions.R and timeSeriesHistograms.R
save(file=resultsFile,results=foreachResults)


#cat(file=logFile,"Log file closing...\n")
#close(logFile)


runtime = proc.time() - runtime
print(runtime)


