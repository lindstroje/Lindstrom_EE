library(vegan)
library(MDMR)
library(dplyr)
library(ggplot2)


ORD <- read.csv("0_ORD_juneaugyears.csv")

##
# PERMANOVA #####
##

ORD_Com <- ORD[,4:ncol(ORD)] #"community" information, taking out all columns not species count data
ORD_Plot_Matrix <- as.matrix(ORD_Com)
ORD_Plot_Matrix <- sqrt(ORD_Plot_Matrix)
ORD.dist <- vegdist(ORD_Plot_Matrix, method="bray")


ORD$Month <- as.factor(ORD$Month)
ORD$Plot <- as.factor(ORD$Plot)
ORD$Year <- as.factor(ORD$Year)
str(ORD)

set.seed(123)
ORDadonis <-adonis(ORD.dist~  Plot * Month + Plot *Year , data = ORD, permutations = 999)
ORDadonis

# with year and month added, significant differences based on plot and interaction of plot and month

###### RSC ####
RSC <- read.csv("0_RSC_juneaugyears.csv")

##
# PERMANOVA #####
##

RSC_Com <- RSC[,4:ncol(RSC)] #"community" information, taking out all columns not species count data
RSC_Plot_Matrix <- as.matrix(RSC_Com)
RSC_Plot_Matrix <- sqrt(RSC_Plot_Matrix)
RSC.dist <- vegdist(RSC_Plot_Matrix, method="bray")


RSC$Month <- as.factor(RSC$Month)
RSC$Plot <- as.factor(RSC$Plot)
RSC$Year <- as.factor(RSC$Year)
str(RSC)

set.seed(123)
RSCadonis <-adonis(RSC.dist~  Plot * Month + Plot *Year , data = RSC, 
                   permutations = 999)
RSCadonis


## NMDS for RSC ##
library (vegan)
RSC_PLOT <- read.csv("./RSC_fulldata2.csv")

RSC_Plot_Com <- RSC_PLOT[,4:ncol(RSC_PLOT)] #"community" information, taking out all columns not species count data
RSC_Plot_Matrix <- as.matrix(RSC_Plot_Com)

set.seed(123)
nmdsRSCPlot = metaMDS(RSC_Plot_Matrix, distance = "bray", k=3, wascores = TRUE, trymax = 999) #using Bray-Curtis distance met
nmdsRSCPlot
plot(nmdsRSCPlot)

yarrr_1 <- yarrr::piratepal(palette = "basel") # Color palette to account for various seed treatments

#str(RSC_PLOT)


# -------------------------------------------------
# shapes are year and treatment is color

#shapes <- (pch=c(15,16))[as.factor(RSC_PLOT$Year)]

shapes<- (pch=c(15,0,16,17,2,19))[as.factor(RSC_PLOT$Month)]

mds.fig.RSC <- ordiplot(nmdsRSCPlot, type = "none", ylim = c(-1,1.1), xlim=c(0,0.5))
points(mds.fig.RSC, "sites", pch = shapes , col = yarrr_1[1], cex=1, select = RSC_PLOT$Plot == "A")
points(mds.fig.RSC, "sites", pch = shapes , col = yarrr_1[2], cex=1, select = RSC_PLOT$Plot == "B") 
points(mds.fig.RSC, "sites", pch = shapes , col = yarrr_1[3], cex=1, select = RSC_PLOT$Plot == "C") 
points(mds.fig.RSC, "sites", pch = shapes , col = yarrr_1[4], cex=1, select = RSC_PLOT$Plot == "D") 
points(mds.fig.RSC, "sites", pch = shapes , col = yarrr_1[5], cex=1, select = RSC_PLOT$Plot == "E") 
points(mds.fig.RSC, "sites", pch = shapes , col = yarrr_1[6], cex=1, select = RSC_PLOT$Plot == "Multi")

# Elipses and legend
ordiellipse(nmdsRSCPlot, RSC_PLOT$Year, conf = 0.9, label = TRUE, pch = 2, col = "black")
#title(main = "RSC Plots RSCination Grouped by Treatment")
legend("bottomright", legend = c("June 2020", "June 2021", "July", "Aug 2020", "Aug 2021", "Sept."), pch=shapes, title="Month")
legend("bottomleft", legend = c("A", "B", "C", "D", "E", "Multi"), pch=16, col = yarrr_1[1:6], title = "Seed Treatment", title.col = "black")


###
# RSC SPLIT YEARS
library(dplyr)

# RSC Locations
RSC_locations <- read.csv("RSC_dist.csv")
RSC_locations <-as.matrix(RSC_locations)
RSC_locations[,1] <- as.numeric(RSC_locations[,1])
RSC_locations[,2] <- as.numeric(RSC_locations[,2])
RSC_locations[,3] <- as.numeric(RSC_locations[,3])
RSC_locations[,4] <- as.numeric(RSC_locations[,4])
RSC_locations[,5] <- as.numeric(RSC_locations[,5])
RSC_locations[,6] <- as.numeric(RSC_locations[,6])
RSC_locations_dist <- as.dist(RSC_locations, upper=T)


# RSC June 2020 filter data
RSC_6_2020 <- RSC %>%
  filter(RSC$Year==2020, RSC$Month==6)

RSC_Com <- RSC_6_2020[,4:ncol(RSC_6_2020)] #"community" information, taking out all columns not species count data
RSC_Plot_Matrix <- as.matrix(RSC_Com)
RSC_Plot_Matrix <- sqrt(RSC_Plot_Matrix)
RSC.dist <- vegdist(RSC_Plot_Matrix, method="bray")

# Mantel test to look at seed treatment/site differences

RSC_6_2020_Mantel <- mantel(xdis=RSC.dist, ydis=RSC_locations_dist)
RSC_6_2020_Mantel
##################

RSC_8_2020 <- RSC %>%
  filter(RSC$Year==2020, RSC$Month==8)

RSC_Com <- RSC_8_2020[,4:ncol(RSC_8_2020)] #"community" information, taking out all columns not species count data
RSC_Plot_Matrix <- as.matrix(RSC_Com)
RSC_Plot_Matrix <- sqrt(RSC_Plot_Matrix)
RSC.dist <- vegdist(RSC_Plot_Matrix, method="bray")


RSC_8_2020_Mantel <- mantel(xdis=RSC.dist, ydis = RSC_locations_dist)
RSC_8_2020_Mantel

#### 2021 Data

RSC_6_2021 <- RSC %>%
  filter(RSC$Year==2021, RSC$Month==6)
RSC_Com <- RSC_6_2021[,4:ncol(RSC_6_2021)] #"community" information, taking out all columns not species count data
RSC_Plot_Matrix <- as.matrix(RSC_Com)
RSC_Plot_Matrix <- sqrt(RSC_Plot_Matrix)
RSC.dist <- vegdist(RSC_Plot_Matrix, method="bray")

RSC_6_2021_Mantel <- mantel(xdis=RSC.dist, ydis = RSC_locations_dist)
RSC_6_2021_Mantel



RSC_8_2021 <- RSC %>%
  filter(RSC$Year==2021, RSC$Month==8)
RSC_Com <- RSC_8_2021[,4:ncol(RSC_8_2021)] #"community" information, taking out all columns not species count data
RSC_Plot_Matrix <- as.matrix(RSC_Com)
RSC_Plot_Matrix <- sqrt(RSC_Plot_Matrix)
RSC.dist <- vegdist(RSC_Plot_Matrix, method="bray")

RSC_8_2021_Mantel <- mantel(xdis=RSC.dist, ydis = RSC_locations_dist)
RSC_8_2021_Mantel


RSC_locations <- read.csv("RSC_dist.csv")
RSC_locations <-as.matrix(RSC_locations)
RSC_locations[,1] <- as.numeric(RSC_locations[,1])
RSC_locations[,2] <- as.numeric(RSC_locations[,2])
RSC_locations[,3] <- as.numeric(RSC_locations[,3])
RSC_locations[,4] <- as.numeric(RSC_locations[,4])
RSC_locations[,5] <- as.numeric(RSC_locations[,5])
RSC_locations[,6] <- as.numeric(RSC_locations[,6])
RSC_locations_dist <- as.dist(RSC_locations, upper=T)

help(as.dist)

####### ORDWAY MANTEL TESTS

ORD <- read.csv("0_ORD_juneaugyears.csv")

ORD_locations <- read.csv("ORD_dist.csv")
ORD_locations <-as.matrix(ORD_locations)
ORD_locations[,1] <- as.numeric(ORD_locations[,1])
ORD_locations[,2] <- as.numeric(ORD_locations[,2])
ORD_locations[,3] <- as.numeric(ORD_locations[,3])
ORD_locations[,4] <- as.numeric(ORD_locations[,4])
ORD_locations[,5] <- as.numeric(ORD_locations[,5])
ORD_locations[,6] <- as.numeric(ORD_locations[,6])
ORD_locations_dist <- as.dist(ORD_locations, upper=T)


# ORD June 2020 filter data
ORD_6_2020 <- ORD %>%
  filter(ORD$Year==2020, ORD$Month==6)

ORD_Com <- ORD_6_2020[,4:ncol(ORD_6_2020)] #"community" information, taking out all columns not species count data
ORD_Plot_Matrix <- as.matrix(ORD_Com)
ORD_Plot_Matrix <- sqrt(ORD_Plot_Matrix)
ORD.dist <- vegdist(ORD_Plot_Matrix, method="bray")

# Mantel test to look at seed treatment/site differences

ORD_6_2020_Mantel <- mantel(xdis=ORD.dist, ydis=ORD_locations_dist)
ORD_6_2020_Mantel
##################

ORD_8_2020 <- ORD %>%
  filter(ORD$Year==2020, ORD$Month==8)

ORD_Com <- ORD_8_2020[,4:ncol(ORD_8_2020)] #"community" information, taking out all columns not species count data
ORD_Plot_Matrix <- as.matrix(ORD_Com)
ORD_Plot_Matrix <- sqrt(ORD_Plot_Matrix)
ORD.dist <- vegdist(ORD_Plot_Matrix, method="bray")


ORD_8_2020_Mantel <- mantel(xdis=ORD.dist, ydis = ORD_locations_dist)
ORD_8_2020_Mantel

#### 2021 Data

ORD_6_2021 <- ORD %>%
  filter(ORD$Year==2021, ORD$Month==6)
ORD_Com <- ORD_6_2021[,4:ncol(ORD_6_2021)] #"community" information, taking out all columns not species count data
ORD_Plot_Matrix <- as.matrix(ORD_Com)
ORD_Plot_Matrix <- sqrt(ORD_Plot_Matrix)
ORD.dist <- vegdist(ORD_Plot_Matrix, method="bray")

ORD_6_2021_Mantel <- mantel(xdis=ORD.dist, ydis = ORD_locations_dist)
ORD_6_2021_Mantel



ORD_8_2021 <- ORD %>%
  filter(ORD$Year==2021, ORD$Month==8)
ORD_Com <- ORD_8_2021[,4:ncol(ORD_8_2021)] #"community" information, taking out all columns not species count data
ORD_Plot_Matrix <- as.matrix(ORD_Com)
ORD_Plot_Matrix <- sqrt(ORD_Plot_Matrix)
ORD.dist <- vegdist(ORD_Plot_Matrix, method="bray")

ORD_8_2021_Mantel <- mantel(xdis=ORD.dist, ydis = ORD_locations_dist)
ORD_8_2021_Mantel


ORD_locations <- read.csv("ORD_dist.csv")
ORD_locations <-as.matrix(ORD_locations)
ORD_locations[,1] <- as.numeric(ORD_locations[,1])
ORD_locations[,2] <- as.numeric(ORD_locations[,2])
ORD_locations[,3] <- as.numeric(ORD_locations[,3])
ORD_locations[,4] <- as.numeric(ORD_locations[,4])
ORD_locations[,5] <- as.numeric(ORD_locations[,5])
ORD_locations[,6] <- as.numeric(ORD_locations[,6])
ORD_locations_dist <- as.dist(ORD_locations, upper=T)


#BETADISPR Analysis

ORD <- read.csv("0_ORD_juneaugyears.csv")
ORD_Com <- ORD[,4:ncol(ORD)] #"community" information, taking out all columns not species count data
ORD_Plot_Matrix <- as.matrix(ORD_Com)
ORD.dist <- vegdist(ORD_Plot_Matrix, method="bray")


dispersion <-betadisper(ORD.dist, group = ORD$Plot)
permutest(dispersion) #not significant, but close not significantly different for seed treatment

dispersionyear <-betadisper(ORD.dist, group=ORD$Year)
permutest(dispersionyear) #not significantly different in year to year variation

plot(dispersion, hull = FALSE, ellipse = TRUE)

ordMDS<-metaMDS(ORD_Plot_Matrix, distance = "bray", k=3, trymax = 35, autotransform = TRUE)
ordMDS
stressplot(ordMDS)

##pull points from MDS
NMDS1 <- ordMDS$points[,1] ##also found using scores(birdMDS)
NMDS2 <- ordMDS$points[,2]
ord.plot<-cbind(ORD, NMDS1, NMDS2)


#plot ordination
p<-ggplot(ord.plot, aes(NMDS1, NMDS2, color=ORD$Plot))+
  geom_point(position=position_jitter(.1), shape=3)+##separates overlapping points
  stat_ellipse(type='t',size =1)+ ##draws 95% confidence interval ellipses
  theme_minimal()
p

ORD2<- ORD%>% 
  mutate(plot_year= paste0(Plot, Year))

newbd <-betadisper(ORD.dist, group = ORD2$plot_year)
permutest(newbd, pairwise = TRUE)


#### RSC DISPERSION

RSC <- read.csv("0_RSC_juneaugyears.csv")
RSC_Com <- RSC[,4:ncol(RSC)] #"community" information, taking out all columns not species count data
RSC_Plot_Matrix <- as.matrix(RSC_Com)
RSC.dist <- vegdist(RSC_Plot_Matrix, method="bray")


dispersion <-betadisper(RSC.dist, group = RSC$Plot)
permutest(dispersion) #not significant, but close not significantly different for seed treatment

dispersionyear <-betadisper(RSC.dist, group=RSC$Year)
permutest(dispersionyear) #not significantly different in year to year variation

plot(dispersion, hull = FALSE, ellipse = TRUE)

RSCMDS<-metaMDS(RSC_Plot_Matrix, distance = "bray", k=3, trymax = 35, autotransform = TRUE)
RSCMDS
stressplot(RSCMDS)

RSC2<- RSC%>% 
  mutate(plot_year= paste0(Plot, Year))

newbd <-betadisper(RSC.dist, group = RSC2$plot_year)
permutest(newbd, pairwise = TRUE)