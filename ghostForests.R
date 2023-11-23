rm(list=ls())
options(warn = 1)
setwd("PATH TO THIS FILE")

realData <- read.csv("./rawdataReal data.csv", stringsAsFactors = TRUE)
realData$Ecosystem <- relevel(realData$Ecosystem, ref="Live Mangrove")

library(lme4)
library(performance)
library(usdm)
library(DHARMa)
library(ranger)
library(readxl)
library(zoo)
library(emmeans)

relative.axis.point<-function(proportion, axis){

  if(axis=="x"){
    return(par("usr")[1] + proportion*(par("usr")[2]-par("usr")[1]))
  }
  
  if(axis=="y"){
    return(par("usr")[3] + proportion*(par("usr")[4]-par("usr")[3]))
  }
}

mer.ci <- function(model, newdata, sims, ncores=1, re.form=NA){
  
  if(!is.na(re.form)){re.form = as.formula(re.form)}
  # bootstrap predictions 
  bb <- bootMer(model,
                nsim=sims,
                FUN=function(fit){predict(fit, newdata, re.form=re.form)},
                parallel="multicore",
                ncpus = ncores)
  
  # extract quantiles over simulations to get CIs
  cis <- t(apply(bb$t, 2, function(x){quantile(x, c(0.025, 0.975), na.rm=TRUE)}))
  
  return(data.frame(fit = bb$t0,
                    lower = cis[,1],
                    upper = cis[,2]))
  
}

cols <- c("#6c923f", "#877060", "#048892")

seaComp <- as.data.frame(read_excel("./rawdata/230715_RawData.xlsx", sheet = "Seagrass Composition"))

seaComp$Location <- na.locf(seaComp$Location)
seaComp$Ecosystem <- na.locf(seaComp$Ecosystem)
seaComp$`Transect / Plot` <- na.locf(seaComp$`Transect / Plot`)

seaCompDf <- do.call("rbind", lapply(split(seaComp, 
                                           f= ~ seaComp$Location + seaComp$Ecosystem + seaComp$`Transect / Plot`),
                                     function(x){
                                       
                                     output <- x[1,1:3]
                                      
                                     sp <- x[!is.na(x$`Species Present`),]
                                     
                                     if(nrow(sp)==0){
                                       output$div = 0
                                       output$propZostera = NA
                                       output$propHalophila = NA
                                       output$propCymodocea = NA
                                       return(output)
                                     } else {
                                       output$propZostera = mean(sp$`Species Composition (%)`[grepl("muelleri", sp$`Species Present`)])
                                       output$propHalophila = mean(sp$`Species Composition (%)`[grepl("ovalis", sp$`Species Present`)])
                                       output$propCymodocea = mean(sp$`Species Composition (%)`[grepl("serrulata", sp$`Species Present`)])
                                       output$div = sum(output[,grepl("prop", colnames(output))] > 0, na.rm=TRUE)
                                       return(output)
                                     }  
                                     
                                     }))
seaCompDf$propZostera[is.na(seaCompDf$propZostera) | is.nan(seaCompDf$propZostera)] = 0
seaCompDf$propHalophila[is.na(seaCompDf$propHalophila) | is.nan(seaCompDf$propHalophila)] = 0
seaCompDf$propCymodocea[is.na(seaCompDf$propCymodocea) | is.nan(seaCompDf$propCymodocea)] = 0
seaCompDf$Ecosystem[seaCompDf$Ecosystem == "Seagrass Meadow"] = "Seagrass"

realData <- merge(realData, seaCompDf,
                  by.x=c("Site", "Ecosystem", "Transect"),
                  by.y=c("Location", "Ecosystem", "Transect / Plot"),
                  all.x=TRUE, all.y=FALSE, sort=FALSE)

realData$TurbidityL <- log(realData$Turbidity)

highestWave <- apply(realData[,grepl("Wave", colnames(realData))], 1, which.max)
realData$highestWave <- realData[,grepl("Wave", colnames(realData))][cbind(1:nrow(realData),highestWave)]

# add in mean mangrove dbh

mangComp <- as.data.frame(read_excel("./rawdata/230715_RawData.xlsx", sheet = "Mangrove Composition"))

mangComp$Location <- na.locf(mangComp$Location)
mangComp$Ecosystem <- na.locf(mangComp$Site)
mangComp$`Plot` <- na.locf(mangComp$`Plot`)

# add in forest death ####

ghostDeath <- data.frame(Site = levels(realData$Site),
                         deathYear = c(2017, 2003, 2019,
                                       2021, 2008, 2001))

realData <- merge(realData, ghostDeath,
                  by.x="Site", by.y="Site", sort=FALSE)
realData$deathYear[realData$Ecosystem != "Dead Mangrove"]

# add in seagrass buffer ####

seaBuff <- read.csv("./rawdata/seagrassBuffer.csv", stringsAsFactors = TRUE)

realData <- merge(realData, seaBuff,
                  by.x="Site", by.y="Site",
                  all.x=TRUE, all.y=TRUE, sort=FALSE)

# how does seagrass in the surrounding area reflect seagrass in the sampled
# sites?

grassSub <- realData[realData$Ecosystem=="Seagrass",]

cor(grassSub[,c("Seagrass.Cover", "grassHa100", "grassHa250", "grassProp100", "grassProp250")])

plot(grassSub$grassHa250 ~ grassSub$grassProp250)

# environmental models over habitat types (new for R1) ####

# for models with valid data across all habitat types
allVars <- c("Sea.Temp", "Air.Temp", "Sediment", "TurbidityL", "Seagrass.Cover")
allMs <- lapply(allVars, function(var){
         print(var)
         modelData <- realData
         modelData$resp = modelData[,var]
         lmer(resp ~ Ecosystem + (1|Site), data=modelData)
       })
sapply(lapply(allMs, simulateResiduals), plot)

allContrasts <- lapply(allMs, function(mod){
  modCont = contrast(emmeans(mod, ~ "Ecosystem"), interaction="pairwise")
  modConf = confint(modCont)
  
  data.frame(contrast = summary(modCont)$Ecosystem_pairwise,
             estimate = summary(modConf)$estimate,
             lci = summary(modConf)$lower.CL,
             uci = summary(modConf)$upper.CL,
             p = summary(modCont)$p.value)
  })

allPerformance <- compare_performance(allMs)

mangVars <- c("Canopy.Cover", "Stem.Density", "DBH")
mangMs <- lapply(mangVars, function(var){
                   
                   print(var)
                   modelData <- droplevels(realData[realData$Ecosystem != "Seagrass",])
                   modelData$resp = modelData[,var]
                   lmer(resp ~ Ecosystem + (1|Site), data=modelData)
                   
                 })
sapply(lapply(mangMs, simulateResiduals), plot)

mangContrasts <- lapply(mangMs, function(mod){
  modCont = contrast(emmeans(mod, ~ "Ecosystem"), interaction="pairwise")
  modConf = confint(modCont)
  
  data.frame(contrast = summary(modCont)$Ecosystem_pairwise,
             estimate = summary(modConf)$estimate,
             lci = summary(modConf)$lower.CL,
             uci = summary(modConf)$upper.CL,
             p = summary(modCont)$p.value)
})
mangPerformance <- compare_performance(mangMs)

perfTable <- data.frame(biophys = c(allVars, mangVars),
                        r2m = c(allPerformance$R2_marginal, mangPerformance$R2_marginal),
                        r2c = c(allPerformance$R2_conditional, mangPerformance$R2_conditional))
write.csv(perfTable, "./outputs/biophysPerformance.csv")

contrTable <- do.call("rbind", c(allContrasts, mangContrasts))
contrTable$biophys = c(rep(allVars, each=3), mangVars)

contrSummary <- matrix("-", ncol=length(unique(contrTable$biophys)), 
                       nrow=3, 
                       dimnames=list(unique(contrTable$contrast),
                                     unique(contrTable$biophys)),
                       byrow=TRUE)
contrSummary[c(1:15,16,19, 22)] <- paste0(sprintf("%.2f", contrTable$estimate), " (",
                                      sprintf("%.2f", contrTable$lci), " - ",
                                      sprintf("%.2f", contrTable$uci), ")",
                                      ifelse(contrTable$p < 0.05, "*", ""))

contrSummary <- t(contrSummary)

write.csv(contrSummary, "./outputs/biophysContrast.csv")

# ecosystem comparison plots ####

pdf("./plots/envModelPlot.pdf", height=7, width=5.5)
par(mfcol=c(4,2), mar=c(0,3,0,1), oma=c(2,2,0.5,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
sapply(1:length(c(allVars, mangVars)), function(n){
  
  print(n)
  var <- c(allVars, mangVars)[n]
  tempM <- c(allMs, mangMs)[[n]]
  
  if(n<=length(allVars)){
    pred.df <- data.frame(Ecosystem=levels(realData$Ecosystem))
    modelData <- realData
    pred.df$pos = 1:3
    } else {
    modelData <- droplevels(realData[realData$Ecosystem != "Seagrass",])
    pred.df <- data.frame(Ecosystem=levels(modelData$Ecosystem))
    pred.df$pos <- 1:2
    }
    
  pred.df <- cbind(pred.df, mer.ci(tempM, newdata=pred.df, sims=999, ncores=6))
  
  plot(NULL, type="n", xaxt="n", xlab="", ylab="",
       xlim=c(0.5,3.5), ylim=range(modelData[,var]))
  
  points(y=modelData[,var], x=jitter(c(1:3)[modelData$Ecosystem], amount=0.25), 
         pch=c(25, 24, 21)[modelData$Ecosystem], bg=cols[modelData$Ecosystem],
         cex=1.15, col="white", lwd=0.5)
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, col=rgb(1,1,1,0.5))
  
  segments(x0=pred.df$pos, x1=pred.df$pos, y0=pred.df$lower, y1=pred.df$upper,
           col=cols, lwd=2)
  points(y=pred.df$fit, x=pred.df$pos, pch=c(25,24,21), bg=cols, cex=1.5)
  
  if(n %in% c(4,8)){
    axis(side=1, at=1:3, labels=c("Live\nmangrove", "Ghost\nforest", "Seagrass\n"), mgp=c(3,1,0))
  } else {axis(side=1, at=1:3, labels=NA)}
  box()
  mtext(side=2, line=1.75, las=0,
        text=c(expression("SST ("*degree*"C)"),
               expression("Air temp ("*degree*"C)"),
               expression("Sediment size ("*mu*"m)"),
               "ln(Turbidity (NTU))",
               "Seagrass cover (%)",
               "Canopy cover (%)",
               expression("Stem density (stems plot"^-1*")",
                          "Diameter at breast height (cm)"))[n], cex=0.8)
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.935, "y"),
       labels=paste0("(", LETTERS[c(1,3,5,7,2,4,6,8)][n], ")"), font=2, adj=0)
  
})
dev.off()

# raw environment plots ####

pdf("./plots/rawEnvPlot.pdf", height=7, width=5.5)
par(mfcol=c(4,2), mar=c(0,3,0,1), oma=c(2,2,0.5,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plotPos <- (1:6)[realData$Site] + seq(-0.25,0.25,len=3)[realData$Ecosystem]
sapply(1:length(c(allVars, mangVars)), function(n){
  
  var <- c(allVars, mangVars)[n]
  plot(realData[,var] ~ plotPos, type="n", xaxt="n", xlab="", ylab="")
  rect(xleft=c(1.5,3.5,5.5), xright=c(2.5,4.5,6.5), 
       ytop=par('usr')[4], ybottom=par("usr")[3], border=NA, col="grey90")
  
  points(realData[,var] ~ plotPos, pch=c(25, 24, 21)[realData$Ecosystem], bg=cols[realData$Ecosystem],
         cex=1.15, lwd=0.75)
  if(n %in% c(4,8)){
    axis(side=1, at=1:6, labels=c("AP", "Bm", "Bo", "GB", "SS", "W"))
  } else {axis(side=1, at=1:6, labels=NA)}
  box()
  mtext(side=2, line=1.75, las=0,
        text=c(expression("SST ("*degree*"C)"),
               expression("Air temp ("*degree*"C)"),
               expression("Sediment size ("*mu*"m)"),
               "ln(Turbidity (NTU))",
               "Seagrass cover (%)",
               "Canopy cover (%)",
               expression("Stem density (stems plot"^-1*")",
                          "Diameter at breast height (cm)"))[n], cex=0.8)
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.935, "y"),
       labels=paste0("(", LETTERS[n], ")"), font=2, adj=0)
  
  if(n==5){
    legend(x=relative.axis.point(0, "x"),
           y=relative.axis.point(0.8, "y"),
           pch=c(24,25,21), pt.bg=cols, bg="white",
           legend=c("Live mangrove", "Ghost forest", "Seagrass"))
  }
  })

dev.off()


# models of seagrass presence in mangrove sites ####

mangs <- droplevels(realData[realData$Ecosystem == "Dead Mangrove", ])
colnames(mangs)[grepl("prop", colnames(mangs))] =
  paste0(colnames(mangs)[grepl("prop", colnames(mangs))], "Mang")

seaMean <- tapply(realData$Seagrass.Cover[realData$Ecosystem == "Seagrass"],
                  realData$Site[realData$Ecosystem == "Seagrass"],
                  mean)/100
seaMean <- data.frame(Site = rownames(seaMean),
                      meanGrass = seaMean)
seaMean$grassDiv <- tapply(realData$div[realData$Ecosystem == "Seagrass"],
                           realData$Site[realData$Ecosystem == "Seagrass"],
                           max, na.rm=TRUE)
seaMean$propZosteraGrass <- tapply(realData$propZostera[realData$Ecosystem == "Seagrass"],
                                   realData$Site[realData$Ecosystem == "Seagrass"],
                                   mean, na.rm=TRUE)

mangs <- merge(mangs, seaMean, 
               by.x="Site", by.y="Site",
               all.x=TRUE, all.y=FALSE, sort=FALSE)

# first let's try a model that includes only dispersal variables: Amount of
# seagrass around and how long it's been since death
mangs$timeDead <- 2023 - mangs$deathYear

# factor of seagrass presence, setting 1 to ref level so we predict probability
# of occurence (rather than probability of absence)
mangs$SeagrassFact <- as.factor(mangs$Seagrass)
mangs$SeagrassFact <- relevel(mangs$SeagrassFact, "1")
mangs$Seagrass.Cover <- mangs$Seagrass.Cover * 100

library(MixRF)

model <- ranger(SeagrassFact ~ timeDead + meanGrass + grassHa250 + Canopy.Cover + Rhizophora +
                              Sea.Temp + TurbidityL + Sediment + highestWave + Site,
                data = mangs, probability=TRUE, keep.inbag=TRUE, importance="impurity")
plot(model$variable.importance ~ as.factor(names(model$variable.importance)))

testData <- mangs[,c("timeDead", "meanGrass", "grassHa250", "Canopy.Cover", "Rhizophora",
                     "Sea.Temp", "TurbidityL", "Sediment", "highestWave", "Site")]
tempPred <- predict(model, data=testData, type="se")




pdf("./plots/seagrassRF.pdf", height=5.7, width=11.5)

par(mar=c(0,0,0,0), tcl=-0.25, ps=10, las=1, mgp=c(3,0.5,0))

xBuffs <- c(0.075, 0.99)
xWidths <- abs(diff(xBuffs) / 5)
xLefts <- xBuffs[1] + xWidths * 0:5
split.screen(rbind(c(xLefts[1],xLefts[2],0.575,0.99),
                   c(xLefts[2],xLefts[3],0.575,0.99),
                   c(xLefts[3],xLefts[4],0.575,0.99),
                   c(xLefts[4],xLefts[5],0.575,0.99),
                   c(xLefts[5],xLefts[6],0.575,0.99),
                   
                   c(xLefts[1],xLefts[2],0.08,0.495),
                   c(xLefts[2],xLefts[3],0.08,0.495),
                   c(xLefts[3],xLefts[4],0.08,0.495),
                   c(xLefts[4],xLefts[5],0.08,0.495),
                   c(xLefts[5],xLefts[6],0.08,0.495)))

# partial prob plots

vars <- c("timeDead", "meanGrass", "grassHa250", "Canopy.Cover", "Rhizophora",
          "Sea.Temp", "TurbidityL", "Sediment", "highestWave")
frameMat <- as.data.frame(matrix(data=colMeans(mangs[,vars]), nrow=200, ncol=length(vars), byrow = TRUE))
frameMat[,vars == "Rhizophora"] = 0
frameMat <- cbind(frameMat, "Beachmere")
colnames(frameMat) = c(vars, "Site")

probPreds <- lapply(c(vars, "Site"), function(v){
  
  tempMat <- frameMat
  
  if(v != "Site"){
  tempMat[,v] = seq(min(mangs[,v]), max(mangs[,v]), len=200)  
  } else {
  tempMat <- tempMat[1:length(unique(realData$Site)),]
  tempMat[,v] = unique(realData$Site)
  }
  tempPred <- predict(model, data=tempMat, type="se")
  
  return(cbind(tempMat[,v], 
               fit=tempPred$predictions[,1],
               se=tempPred$se[,1]))
  
})

localImp <- model$variable.importance
localImp <- localImp / sum(localImp)

lapply(1:length(probPreds), function(n){
  print(n)
  
  screen(n)
  x<-probPreds[[n]]
  
  if(!n %in% c(5,10)){
  plot(x[,1:2], type="l", ylim=c(0,0.5), col="red", xlab="", ylab="", xaxt="n", yaxt="n",
       lwd=1.5)
  axis(side=1, mgp=c(3,0.2,0))
  polygon(x=c(x[,1], rev(x[,1])),
          y=c(x[,2] + x[,3], rev(x[,2] - x[,3])),
          border=NA, col=rgb(1,0,0,0.2))
  } else {
    
  uniqueVars <- x[x[,1] %% 1 == 0, ]
  plot(uniqueVars[,1:2], ylim=c(0,0.7), type="n", col="red", xlab="", ylab="", xaxt="n", yaxt="n",
       lwd=1.5)
  
  if(n==5){uniqueVars[,1] = c(0.25, 0.75)} else{
    uniqueVars <- uniqueVars[order(uniqueVars[,1]),]
    uniqueVars[,1] = seq(1.35,5.65,len=6)}
  
  segments(x0=uniqueVars[,1], x1=uniqueVars[,1],
           y0=uniqueVars[,2] + uniqueVars[,3],
           y1=uniqueVars[,2] - uniqueVars[,3], col="red")
  points(uniqueVars[,2] ~ uniqueVars[,1], pch=21, bg="red")
  if(n == 5){
  axis(side=1, at=uniqueVars[,1], labels=c("Absent", "Present"), mgp=c(3,0.2,0))
  } else {
  axis(side=1, at=uniqueVars[,1], labels=c("AP", "Bm", "Bo", "GB", "SS", "W"), mgp=c(3,0.2,0))
  }
}
  
  mtext(side=1, line=1.15,
        text=c("Time since forest death (yrs)", 
               "Mean plot seagrass cover (%)", 
               "Seagrass area (ha: 250m radius)", 
               "Canopy cover (%)",
               "Rhizophora presence",
               expression("SST ("*degree*"C)"),
               "ln(Turbidity (NTU))",
               expression("Sediment size ("*mu*"m)"),
               "Wave height (m)", 
               "Site")[n])
  
  if(n %in% c(1,6)){
    axis(side=2)
    mtext(side=2, line=2, text="P(Seagrass presence)", las=0)
  } else {
    axis(side=2, labels=NA)
  }
  
  text(x=par("usr")[1], y=par("usr")[4]-0.035, pos=4, offset=0.25,
       labels=paste0("(", LETTERS[n], ")"), font=2)
  
  text(x=par("usr")[1], y=par("usr")[4]-0.035, pos=4, offset=1.35,
       labels=paste0("Importance = ", sprintf("%.2f", localImp[n]*100), "%"))
  close.screen(n)
})

close.screen(all.screens=TRUE)
dev.off()

# Model predictions
pdf("./plots/seagrassRFPreds.pdf", height=4.5, width=3)
par(mar=c(0.5,3.5,0.5,0.5), tcl=-0.25, ps=10, las=1, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, xlim=c(0.7,6.3), ylim=c(0,1), xlab="", ylab="", xaxt="n")

uniquePreds <- data.frame(fit = tempPred$predictions[,1],
                          se =   tempPred$se[,1])
uniquePreds$pos <- rep(1:6, each=3) + rep(seq(-0.2,0.2, len=3), 6)

segments(x0=uniquePreds$pos, 
         x1=uniquePreds$pos,
         y0=uniquePreds[,1] - uniquePreds[,2],
         y1=uniquePreds[,1] + uniquePreds[,2],
         col=rgb(1,0,0,1))

points(y=uniquePreds[,1],
       x=uniquePreds$pos, 
       pch=21, bg="red", lwd=1)

siteMeans <- tapply(mangs$Seagrass, factor(mangs$Site, levels=unique(mangs$Site)), mean)
points(y=siteMeans, x=1:6, pch=24, bg=cols[2])

text(x=1:6,
     y=siteMeans, 
     labels=unique(mangs$Site),
     pos=c(2,3,3,3,3,2))
mtext(side=2, line=2, text="P(Seagrass presence in ghost forest)", las=0)

legend(x=par("usr")[1], y=0.5,
       pch=c(24,21), pt.bg=c(cols[2], "red"),
       legend=c("Observed", "Predicted"))

dev.off()

# variable importance plot ####
pdf("./plots/RFimportance.pdf", height=4.5, width=8)
par(mar=c(3,3.5,0.5,3.5), tcl=-0.25, ps=10, las=1, mgp=c(3,0.5,0))
plot(y=model$variable.importance, x=(1:length(model$variable.importance)), pch=16,
     ylim=c(0,2), xaxt="n", xlab="", ylab="")
axis(side=1, at=1:10, labels = c("Time of\nforest death", "Seagrass\ncover", "Seagrass\narea", "Canopy\ncover", "Rhizophora\npresence",
                                "SST\n", "Turbidity\n", "Sediment\nsize", "Wave\nheight", "Site\n"), mgp=c(3,1,0),
     gap.axis=0)
mtext(side=2, line=2, las=0, text="Variable importance (Gini index)")
mtext(side=4, line=2, las=0, text="Variable importance (% of total)")
impRatio = sum(model$variable.importance)
axis(side=4, at=impRatio * seq(0, 0.35, 0.05), labels=seq(0, 0.35, 0.05)*100)

dev.off()

# models of seagrass at sites based on environment ####

# add in env vars? Only Comp 1 is viable in a model, as the rest correlate too strongly
# with seagrass?

cor(grassSub[,c("Seagrass.Cover", "grassHa250", "grassProp250")])

grassM1 <- lm(grassHa250 ~ Seagrass.Cover, data=grassSub)
grassM2 <- lm(grassProp250 ~ Seagrass.Cover, data=grassSub)
plot(simulateResiduals(grassM2))

pdf("./plots/seagrassCover.pdf", height=4.5, width=4.5, useDingbats = FALSE)
par(mar=c(3,3,1,1), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(grassSub$grassProp250 ~ grassSub$grassHa250, type="n", ylim=c(0,1), xaxt="n",
     xlab="",ylab="")
axis(side=1, mgp=c(3,0.2,0))
mtext(side=1,line=1.25, text="Seagrass area (ha: 250m radius)")
mtext(side=2,line=1.75, text="Proportion seagrass (250m radius)", las=0)
points((grassSub$grassProp250 + seq(-0.03,0.03,len=3)) ~ grassSub$grassHa250,
       cex=seq(0.5,3,len=5)[grassCut],
       bg=cols[3], pch=21)
text(x=unique(grassSub$grassHa250),y=unique(grassSub$grassProp250),
     labels=unique(grassSub$Site), pos=c(2,4,4,4,2,4))

legend(x="topright",
       pch=rep(21,5),
       legend=c("Seagrass", rep("", 4)),
       pt.bg=c(rep(cols[3],5)),
       col=c(rep("black", 5)),
       pt.cex=c(seq(0.5,3,len=5)),
       pt.lwd=c(rep(1,5)))
dev.off()

plot(grassSub$Seagrass.Cover ~ grassSub$grassProp250)
plot(grassSub$grassHa250 ~ grassSub$grassProp250)

# does adding in ENV variables make seagrass model better? ####

# do env PCAs correlate with seagrass vars?

grassEnvCor <- cor(mangs[,c("timeDead", "meanGrass", "grassHa250", "grassProp250", "Comp.1", "Comp.2", "Comp.3", "Comp.4")])[1:3,]

modelEnv <- ranger(SeagrassFact ~ timeDead + meanGrass + grassHa250 + grassProp250 + Comp.1 + Comp.2 + Comp.3 + Comp.4, 
                data = mangs, probability=TRUE, keep.inbag=TRUE, importance="impurity")
plot(modelEnv$variable.importance ~ as.factor(names(modelEnv$variable.importance)))
model$prediction.error
modelEnv$prediction.error

modelEnv1 <- ranger(SeagrassFact ~ Comp.1 + Comp.2 + Comp.3 + Comp.4, 
                   data = mangs, probability=TRUE, keep.inbag=TRUE, importance="impurity")
modelEnv1$prediction.error

plot(x=NULL, y=NULL, xlim=c(1,8), ylim=c(0,1), xaxt="n", xlab="", ylab="")
points(y=model$variable.importance / sum(model$variable.importance), 
       x=1:length(model$variable.importance))
points(y=modelEnv$variable.importance / sum(modelEnv$variable.importance), x=1:length(modelEnv$variable.importance), col="red")
axis(side=1, labels=NA)

# what about species compositional data?

modelComp <- ranger(SeagrassFact ~ timeDead + meanGrass + grassHa250 + grassProp250 + 
                      Rhizophora + grassDiv + propZosteraGrass, 
                   data = mangs, probability=TRUE, keep.inbag=TRUE, importance="impurity")
plot(modelComp$variable.importance ~ as.factor(names(modelComp$variable.importance)))
model$prediction.error
modelComp$prediction.error

table(mangs$Site, mangs$grassDiv)

# models of PC axes and variation ####

realData <- cbind(realData, pc$scores[,1:4])

# ecosystem explains 0 so removed as random effect in PC1
pc1mS <- lmer(Comp.1 ~ 1 + (1|Site), data=realData)
pc1mSE <- lmer(Comp.1 ~ 1 + (1|Site) + (1|Ecosystem), data=realData)
performance(pc1mSE)
plot(simulateResiduals(pc1m))

pc2mS <- lmer(Comp.2 ~ 1 + (1|Site), data=realData)
pc2mSE <- lmer(Comp.2 ~ 1 + (1|Site) + (1|Ecosystem), data=realData)
performance(pc2mSE)
plot(simulateResiduals(pc2m))

pc3mS <- lmer(Comp.3 ~ 1 + (1|Site), data=realData)
pc3mSE <- lmer(Comp.3 ~ 1 + (1|Site) + (1|Ecosystem), data=realData)
performance(pc3mSE)
plot(simulateResiduals(pc3m))

pc4mS <- lmer(Comp.4 ~ 1 + (1|Site), data=realData)
pc4mSE <- lmer(Comp.4 ~ 1 + (1|Site) + (1|Ecosystem), data=realData)
performance(pc4mSE)
plot(simulateResiduals(pc4m))

