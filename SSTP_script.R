########################################################################
########################################################################
##                                                                    ##        
## Spatially shifting temporal points: estimating pooled              ##
## within-time-series variograms for scarce hydrological data          ##
##                                                                    ##
## Avit Kumar Bhowmik                                                 ##
## E-mail: bhowmik@uni-landau.de                                      ##
## Quantitative Landscape Ecology,                                    ## 
## Institute for Environmental Sciences                               ##
## University of Koblenz-Landau                                       ##
## Fortstra?e 7, 76829 Landau in der Pfalz, Germany                   ##
##                                                                    ##
## Pedro Cabral                                                       ##
## E-mail: pcabral@novaims.unl.pt                                     ##
## NOVA IMS, Universidade Nova de Lisboa                              ## 
## 1070-312 Lisboa, Portugal                                          ## 
##                                                                    ##
##                                                                    ##
## Corresponding Author:                                              ##
## Avit Kumar Bhowmik                                                 ##
## Telephone: +49 6341 280-31331                                      ##
## Fax: +49 6341 280-31326                                            ##
##                                                                    ##
########################################################################
########################################################################
##                    Supplementary materials 2                       ##
########################################################################
########################################################################
##                                                                    ##
##                               SSTP                                 ##
##                                                                    ##
##        Estimation of pooled within-time series variograms          ##
##                               with                                 ##
##              spatially  shifting temporal points                   ##
##                                                                    ##
########################################################################
########################################################################

## The free open source software package “R” is required with three packages “spacetime”,
## “intamap” and “gstat” to perform the program (see references in the paper for documentations).
## R software package can be downloaded from http://www.r-project.org/. Once the software is
## downloaded and installed (please follow the installation guide provided on the website) the
## mentioned packages can be installed by executing the following codes in the command console:

install.packages("spacetime")
install.packages("intamap")
install.packages("gstat")

# If it asks for the geographic region while installing, select the region of interest and press Enter.
# The installed packages need to be initiated in R environment by the following commands:

library(spacetime)
library(intamap)
library(gstat)

# A directory needs to be set as the working directory of R and the example data “Data_SSTP.Rdata”,
# which can be temporarilly downloaded from: https://www.dropbox.com/s/eeqkiq54ttwkp0c/Test_Data.Rdata,
# and will be archived permanently in PANAGEA. The data needs to be copied in the directory and loaded
# in the R environment. If a folder called “Sample_Data” in the directory “C:\” (Windows OS) contains
# the example data, the working directory can be set and the data can be loaded in R and checked for
# details by executing the following codes:

setwd("C:/Sample_Data")

# Note for Windows OS users: if you copy the code on a Windows OS, please change the slash “\” to
# back slash “/” while setting the working directory. 

# Linux example
# setwd("/media/storage/projects/SSTP/Data/")
# Mac OS X example
# setwd("/Users/<username>/projects/SSTP/Data/")

# The data will be loaded into the R environment and checked for its componenets
load("SM3_Bhowmik_Cabral.Rdata")
str(Data_SSTP)

# Data_SSTP is the “spacetimedataframe” object where the spatial,
# temporal and attribute information are stored. “spPoints” are the spatial coordinates with
# WGS84 reference, “timepoints” are the temporal points (years) between 1993-2007 and the
# “dataObj” is the attribute (PRCPTOT) at each spPoint and each timePoint. Two separation
# distances need to be defined for spatial shifting of the temporal points. As described
# in the paper, the separations distance between two shifted temporal points is set to
# 1111km ≈ 10 decimal degree as following:
  
sepDist <- 10

# To calculate the smallest- and largest-spatial-lags for a year or within a time series, see
# the "spDists" function. Now the arbitrary spatial coordinates need to be created to shift the
# temporal data. Since 15 years of PRCPTOT data (1993-2007) need to be shifted and assigned to
# different spatial coordinates, the spatial coordinates sets will be created accordingly.
# The spatial coordinates can be obtained from Data_SSTP and new coordinates sets can be created
# by executing the following codes:

allCoords <- Data_SSTP@sp@coords
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(sepDist,32),rep(0,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(0,32),rep(sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(-sepDist,32),rep(0,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(0,32),rep(-sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(sepDist,32),rep(sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(-sepDist,32),rep(sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(sepDist,32),rep(-sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(-sepDist,32),rep(-sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(0,32),rep(2*sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(0,32),rep(-2*sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(sepDist,32),rep(2*sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(-sepDist,32),rep(2*sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(sepDist,32),rep(-2*sepDist,32)))
allCoords <- rbind(allCoords,Data_SSTP@sp@coords+cbind(rep(-sepDist,32),rep(-2*sepDist,32)))

# Now the PRCPTOT data from 1993-2007 need to be assigned to the 15 spatial coordinates sets
# and finally a single “SpatialPointsDataFrame” object can be created with the shifted temporal
# data using the following code. To give the arbitrary coordinates an unique id, their rownames are
# defined as the rownames of the shifted PRCPTOT data. The year (time step) information will
# also be stored. 

allData <- Data_SSTP[,1:15]@data
allData$Years <- rep(1993:2007, each=32)
rownames(allCoords)<-rownames(allData)
allSpPoints <-  SpatialPointsDataFrame(allCoords,allData,proj4string=CRS(proj4string(Data_SSTP)))

# The output object can be plotted to check the temporally observed PRCPTOT values at the
# shifted spatial points, which now corresponds to the temporal points. The output plot can be
# compared to the file “Output 1.tiff”.

spplot(allSpPoints["PRCPTOTWet"], col.regions=bpy.colors(), scales=list(draw=T), colorkey=T)

# For the further steps, the “NA” values need to be removed from the object:

allSpPoints<-allSpPoints[!is.na(allSpPoints@data[[1]]),]

# The next step is to check for anisotropy within the spatially shifted period. This can be done
# by executing the following functions from the “intamap” package in R. Please consult the package
# vignettes for the methodology of anisotropy estimation.

Set  <- 1:15
spacetimedata <- lapply(Set, function(i) {x = Data_SSTP[,i]; x$years = i+1992; x})
spacedata <- do.call(rbind, spacetimedata)
params=NULL
estimateAnisotropy(spacedata[!is.na(spacedata$PRCPTOTWet),], PRCPTOTWet, PRCPTOTWet~1)


# The output “ratio” value is the ratio of the minor and major axis of the ellipse (2.133343) and
# direction is the angle of the anisotropic axis from the normal east to the clockwise directions
# (41.90442). This ratio of the minor and major axis of the ellipse requires transformation before
# input in “gstat” functions as it takes ratio between major and minor axes (A:B) = 1 / 2.133343
# = 0.47 as input. The anisotropy angle remains the same because alpha = anisotropy angle (Ø) from
# normal North to anti-clockwise direction = alternate interior angle of the angle of the anisotropic
# axis from the normal east to the clockwise directions = 41.90. See package vignettes for details.

# The next step is to compute the within-strata variogram applying the functions of “gstat”
# package. While computing the variograms, the smallest and largest-spatial-lags and the
# anisotropy parameters are controlled with the "width" and "cutoff" arguments, respectively
# as described in the paper. Finally, the computed variogram can be plotted. Users may try 
# with other variogram models and parameters.

PRCPTOTWet.dir.SSTP = variogram(PRCPTOTWet ~ 1, allSpPoints, alpha=41.90442, width=27.51, cutoff=550)
plot(PRCPTOTWet.dir.SSTP)
PRCPTOTWetdir.model.SSTP = vgm(psill=5, model="Pow", range=1.93, nugget=185800, anis=c(41.90442, 0.4687479))
PRCPTOTWetdir.fit.SSTP = fit.variogram(PRCPTOTWet.dir.SSTP, PRCPTOTWetdir.model.SSTP, fit.sills = TRUE,
                                       fit.ranges = TRUE, fit.method = 7)
plot(PRCPTOTWet.dir.SSTP, PRCPTOTWetdir.fit.SSTP)


########################################################################
##                               End                                  ##
########################################################################