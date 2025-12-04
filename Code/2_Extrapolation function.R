Surv_analysis1<-function(data){
  data<-data[-(1:4),]
  na_cols <- colSums(is.na(data)) == nrow(data)
  data <- data[, !na_cols]
  colnames(data)<-c("Time","Event","Time","Event","Time","Event","Time","Event","Time","Event","Time","Event")

  ####-----####
  #  Group 1  #
  ####-----####
  df1<-na.omit(data[,3:4])
  df1$Time<-as.numeric(df1$Time)
  df1$Event<-as.numeric(df1$Event)
  

  df2<-na.omit(data[,1:2])
  df2$Time<-as.numeric(df2$Time)
  df2$Event<-as.numeric(df2$Event)
  

  df5<-na.omit(data[,5:6])
  df5$Time<-as.numeric(df5$Time)
  df5$Event<-as.numeric(df5$Event)
  

  ##################################################
  #############   generate  results   ##############
  ##################################################
  # extrapolation 1
  trt1<-as.data.frame(rename(df1,"censrec"="Event","recyrs"="Time"))
  trt1$recyrs<-trt1$recyrs/12
  trt1$rectime<-trt1$recyrs*365
  trt1$rectime2 <- as.integer(trt1$rectime/(365.24/12)) + 1
  lttrt1 <- lifeTable(trt1, timeColumn = "rectime2", eventColumn = "censrec")
  ltHaz_1 <- data.frame(hazKM = lttrt1$Output$hazard, Time = (seq(1:length(lttrt1$Output[,1]))-0.5)/12,
                        AtRisk = lttrt1$Output$atRisk, Events = lttrt1$Output$events)
  ltHaz_1$hazLT = ltHaz_1$Events / (ltHaz_1$AtRisk - ltHaz_1$Events/2)
  ltHaz_1$lnTime <- log(ltHaz_1$Time)
  ltHaz_1$MyId <- 1:dim(ltHaz_1)[1] # Generate id variable
  ltHaz_1$EventsL <- lag(ltHaz_1$Events)
  ltHaz_1$EventsL[1] <- 0
  ltHaz_1$surv <- lttrt1$Output$S
  ltHaz_1$timedelta<-ltHaz_1$Time[2]-ltHaz_1$Time[1]
  in_ipd1<-data.frame("time"=trt1$recyrs,"event"=trt1$censrec)
  
  # extrapolation 2
  trt2<-as.data.frame(rename(df2,"censrec"="Event","recyrs"="Time"))
  trt2$recyrs<-trt2$recyrs/12
  trt2$rectime<-trt2$recyrs*365
  trt2$rectime2 <- as.integer(trt2$rectime/(365.24/12)) + 1
  lttrt2 <- lifeTable(trt2, timeColumn = "rectime2", eventColumn = "censrec")
  ltHaz_2 <- data.frame(hazKM = lttrt2$Output$hazard, Time = (seq(1:length(lttrt2$Output[,1]))-0.5)/12,
                        AtRisk = lttrt2$Output$atRisk, Events = lttrt2$Output$events)
  ltHaz_2$hazLT = ltHaz_2$Events / (ltHaz_2$AtRisk - ltHaz_2$Events/2)
  ltHaz_2$lnTime <- log(ltHaz_2$Time)
  ltHaz_2$MyId <- 1:dim(ltHaz_2)[1] # Generate id variable
  ltHaz_2$EventsL <- lag(ltHaz_2$Events)
  ltHaz_2$EventsL[1] <- 0
  ltHaz_2$surv <- lttrt2$Output$S
  ltHaz_2$timedelta<-ltHaz_2$Time[2]-ltHaz_2$Time[1]
  in_ipd2<-data.frame("time"=trt2$recyrs,"event"=trt2$censrec)
  
  # extrapolation 5
  trt5<-as.data.frame(rename(df5,"censrec"="Event","recyrs"="Time"))
  trt5$recyrs<-trt5$recyrs/12
  trt5$rectime<-trt5$recyrs*365
  trt5$rectime2 <- as.integer(trt5$rectime/(365.24/12)) + 1
  lttrt5 <- lifeTable(trt5, timeColumn = "rectime2", eventColumn = "censrec")
  ltHaz_5 <- data.frame(hazKM = lttrt5$Output$hazard, Time = (seq(1,length(lttrt5$Output[,1]))-0.5)/12,
                        AtRisk = lttrt5$Output$atRisk, Events = lttrt5$Output$events)
  ltHaz_5$hazLT = ltHaz_5$Events / (ltHaz_5$AtRisk - ltHaz_5$Events/2)
  ltHaz_5$lnTime <- log(ltHaz_5$Time)
  ltHaz_5$MyId <- 1:dim(ltHaz_5)[1] # Generate id variable
  ltHaz_5$EventsL <- lag(ltHaz_5$Events)
  ltHaz_5$EventsL[1] <- 0
  ltHaz_5$surv <- lttrt5$Output$S
  ltHaz_5$timedelta<-ltHaz_5$Time[2]-ltHaz_5$Time[1]
  in_ipd5<-data.frame("time"=trt5$recyrs,"event"=trt5$censrec)
  
  #New Time Data 1
  follow_up <- round(ltHaz_1$Time[length(ltHaz_1$Time)]*12,0) #origin study (Month)
  numMod <- 6 # Models considered
  dfHazEst11 <- array(dim=c(numMod, length(ltHaz_5$Time)))
  dfHazEst12 <- array(dim=c(numMod, length(ltHaz_5$Time)))
  Newtime_1 <- data.frame(Time = ltHaz_5$Time, AtRisk = 1)
  Newtime_1$MyId <- 1:dim(Newtime_1)[1]
  Newtime_1$MyId <- ifelse(Newtime_1$MyId > follow_up, follow_up, Newtime_1$MyId)  # Random effects: Using last observed ID for extrapolation
  Newtime_1$EventsL <- 0
  Newtime_1$EventsL[1:follow_up] <- lag(ltHaz_1$Events)
  Newtime_1$EventsL[1] <- 0
  Newtime_1$EventsL <- ifelse(Newtime_1$MyId > follow_up, 0, Newtime_1$EventsL) # AR: Using last observed event count for extrapolation
  Newtime_1$timedelta<-Newtime_1$Time[2]-Newtime_1$Time[1]
  Newtime_1$lnTime<-log(Newtime_1$Time)
  dfGOF1 <- data.frame(matrix(nrow=12, ncol=5))
  colnames(dfGOF1) <- c("Model","LnL","Params","AIC","BIC")
  # Below is constant for when have to derive log-likelihood
  llCons1 <- sum(ltHaz_1$Events*log(ltHaz_1$AtRisk) - log(factorial(ltHaz_1$Events)))
  
  #New Time Data 2
  follow_up <- round(ltHaz_2$Time[length(ltHaz_2$Time)]*12,0) #origin study (Month)
  numMod <- 6 # Models considered
  dfHazEst21 <- array(dim=c(numMod, length(ltHaz_5$Time)))
  dfHazEst22 <- array(dim=c(numMod, length(ltHaz_5$Time)))
  Newtime_2 <- data.frame(Time = ltHaz_5$Time, AtRisk = 1)
  Newtime_2$MyId <- 1:dim(Newtime_2)[1]
  Newtime_2$MyId <- ifelse(Newtime_2$MyId > follow_up, follow_up, Newtime_2$MyId)  # Random effects: Using last observed ID for extrapolation
  Newtime_2$EventsL <- 0
  Newtime_2$EventsL[1:follow_up] <- lag(ltHaz_2$Events)
  Newtime_2$EventsL[1] <- 0
  Newtime_2$EventsL <- ifelse(Newtime_2$MyId > follow_up, 0, Newtime_2$EventsL) # AR: Using last observed event count for extrapolation
  Newtime_2$timedelta<-Newtime_2$Time[2]-Newtime_2$Time[1]
  Newtime_2$lnTime<-log(Newtime_2$Time)
  dfGOF2 <- data.frame(matrix(nrow=12, ncol=5))
  colnames(dfGOF2) <- c("Model","LnL","Params","AIC","BIC")
  # Below is constant for when have to derive log-likelihood
  llCons2 <- sum(ltHaz_2$Events*log(ltHaz_2$AtRisk) - log(factorial(ltHaz_2$Events)))
  
  ###-----------###
  #Extrapolation-1
  ###-----------###
  #----Standard survival model (SSM)----
  temp_GOF<- data.frame(matrix(nrow=7, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  MyDists <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
  for (i in 1:7){
    glmTemp <- flexsurvreg(Surv(time, event) ~ 1, data = in_ipd1, dist = MyDists[[i]])
    temp_GOF[i,1] <- MyDists[[i]]
    temp_GOF[i,2] <- glmTemp$loglik
    temp_GOF[i,3] <- glmTemp$npars
    temp_GOF[i,4] <- glmTemp$AIC
    temp_GOF[i,5] <- log(length(in_ipd1$time))*glmTemp$npars-2*glmTemp$loglik
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  ssm_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  ssm_bic<-temp_GOF$Model[1]
  
  glm_ssm1<-flexsurvreg(Surv(time, event) ~ 1, data = in_ipd1, dist = ssm_aic)
  dfHazEst11[1,] <- summary(glm_ssm1, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[ssm_aic]] <- summary(glm_ssm1, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[1]<-ssm_aic
  dfGOF1$LnL[1]<-sum(ltHaz_1$Events*log(ltHaz_1[[ssm_aic]]) - ltHaz_1[[ssm_aic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[1]<-glm_ssm1$npars
  dfGOF1$AIC[1]<-2*glm_ssm1$npars-2*dfGOF1$LnL[1]
  dfGOF1$BIC[1]<-log(length(in_ipd1$time))*glm_ssm1$npars-2*dfGOF1$LnL[1]
  
  glm_ssm2<-flexsurvreg(Surv(time, event) ~ 1, data = in_ipd1, dist = ssm_bic)
  dfHazEst12[1,] <- summary(glm_ssm2, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[ssm_bic]] <- summary(glm_ssm2, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[2]<-ssm_bic
  dfGOF1$LnL[2]<-sum(ltHaz_1$Events*log(ltHaz_1[[ssm_bic]]) - ltHaz_1[[ssm_bic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[2]<-glm_ssm2$npars
  dfGOF1$AIC[2]<-2*glm_ssm2$npars-2*dfGOF1$LnL[2]
  dfGOF1$BIC[2]<-log(length(in_ipd1$time))*glm_ssm2$npars-2*dfGOF1$LnL[2]
  #-----end ssm -----
  
  #----Fractional polynomial (FP)----
  temp_GOF<- data.frame(matrix(nrow=44, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  #-----FP1 -----
  MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
  for (i in 1:7){
    glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    temp_GOF$Model[i] <- paste(MyPowers[[1]][i],"NA",sep = ",") 
    temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  }
  ### run for 0
  glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  temp_GOF$Model[8] <- paste(0,"NA",sep = ",")
  temp_GOF$LnL[8] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  temp_GOF$Params[8]<-extractAIC(glmTemp)[1]
  temp_GOF$AIC[8] <- extractAIC(glmTemp)[2]
  temp_GOF$BIC[8] <- -2*temp_GOF$LnL[8] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  
  #-----FP2 -----
  myLnL <- array(dim=36)
  myAIC <- array(dim=36)
  MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
  index <- 9
  for (i in 1:7){
    for (j in 1:7){
      if (j > i) {
        glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][j])+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1) 
        temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
        temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
        temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
        temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
        temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
        index <- index + 1
      }
    }
  }
  for (i in 1:7) {
    glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][i]*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)# 
    temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
    temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    index <- index + 1
    
  }
  
  for (i in 1:7) {
    glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)# 
    temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
    temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    index <- index + 1
  }
  
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)# 
  temp_GOF$Model[index] <- paste(0,0,sep = ",")
  temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
  temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
  temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  
  temp_GOF<-arrange(temp_GOF,AIC)
  fp_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  fp_bic<-temp_GOF$Model[1]
  
  FP_pow_1<-as.numeric((strsplit(fp_aic,","))[[1]][1])
  FP_pow_2<-as.numeric((strsplit(fp_aic,","))[[1]][2])
  if(is.na(FP_pow_2)){
    if(FP_pow_1 == 0){
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    } else {
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    }
  } else if (FP_pow_1 == FP_pow_2 & FP_pow_2 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else if (FP_pow_1==0 & FP_pow_2== 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else if (FP_pow_2==0 & FP_pow_1 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else {
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  }
  bestfp<-paste("FP_aic",fp_aic,sep = ":")
  dfHazEst11[2,] <- predict(modFP, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestfp] <- predict(modFP, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[3]<-bestfp
  dfGOF1$LnL[3]<-(extractAIC(modFP)[2] - 2*extractAIC(modFP)[1])*(-0.5)
  dfGOF1$Params[3]<-extractAIC(modFP)[1]
  dfGOF1$AIC[3]<-extractAIC(modFP)[2]
  dfGOF1$BIC[3]<--2*dfGOF1$LnL[3] + log(length(ltHaz_1$Time))*extractAIC(modFP)[1]
  
  FP_pow_1<-as.numeric((strsplit(fp_bic,","))[[1]][1])
  FP_pow_2<-as.numeric((strsplit(fp_bic,","))[[1]][2])
  if(is.na(FP_pow_2)){
    if(FP_pow_1 == 0){
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    } else {
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    }
  } else if (FP_pow_1 == FP_pow_2 & FP_pow_2 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else if (FP_pow_1==0 & FP_pow_2== 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else if (FP_pow_2==0 & FP_pow_1 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else {
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  }
  bestfp<-paste("FP_bic",fp_bic,sep = ":")
  dfHazEst12[2,] <- predict(modFP, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestfp] <- predict(modFP, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[4]<-bestfp
  dfGOF1$LnL[4]<-(extractAIC(modFP)[2] - 2*extractAIC(modFP)[1])*(-0.5)
  dfGOF1$Params[4]<-extractAIC(modFP)[1]
  dfGOF1$AIC[4]<-extractAIC(modFP)[2]
  dfGOF1$BIC[4]<--2*dfGOF1$LnL[4] + log(length(ltHaz_1$Time))*extractAIC(modFP)[1]
  #-----end fp -----
  
  #----Restricted cubic splines (RCS)----
  bc2 <- subset(in_ipd1, event==1)
  temp_GOF<- data.frame(matrix(nrow=5, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  for (i in 1:5){
    fit<-try(glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_1))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- 'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i] <- 'error'
      temp_GOF$BIC[i] <- 'error'
    }
    else{
      gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_1)      
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
      temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
      temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  rcs_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  rcs_bic<-temp_GOF$Model[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=rcs_aic+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+rcs_aic))), length=rcs_aic+2),family=binomial(link=cloglog), data=ltHaz_1)
  bestrcs<-paste("RCS_aic",rcs_aic,sep = ":")
  dfHazEst11[3,] <- predict(glmTemp, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestrcs] <- predict(glmTemp, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[5]<-bestrcs
  dfGOF1$LnL[5]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF1$Params[5]<-extractAIC(glmTemp)[1]
  dfGOF1$AIC[5]<-extractAIC(glmTemp)[2]
  dfGOF1$BIC[5]<--2*dfGOF1$LnL[5] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=rcs_bic+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+rcs_bic))), length=rcs_bic+2),family=binomial(link=cloglog), data=ltHaz_1)
  bestrcs<-paste("RCS_bic",rcs_bic,sep = ":")
  dfHazEst12[3,] <- predict(glmTemp, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestrcs] <- predict(glmTemp, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[6]<-bestrcs
  dfGOF1$LnL[6]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF1$Params[6]<-extractAIC(glmTemp)[1]
  dfGOF1$AIC[6]<-extractAIC(glmTemp)[2]
  dfGOF1$BIC[6]<--2*dfGOF1$LnL[6] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  #-----end rcs -----
  
  #----Royston-Parmar models (RP)----
  temp_GOF<- data.frame(matrix(nrow=18, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  MyScale <- list("hazard","odds","normal")
  for (i in 1:3){
    for (j in 0:5){
      fit<-try(glmTemp <- flexsurvspline(Surv(time, event) ~ 1, data = in_ipd1, k = j, scale = MyScale[[i]]))
      if("try-error" %in% class(fit)) {
        temp_GOF[(i-1)*6+1+j,1] <- "error"
        temp_GOF[(i-1)*6+1+j,2] <- "error"
        temp_GOF[(i-1)*6+1+j,3] <- "error"
        temp_GOF[(i-1)*6+1+j,4] <- "error"
        temp_GOF[(i-1)*6+1+j,5] <- "error"
      }
      else{
        flexsurvspline(Surv(time, event) ~ 1, data = in_ipd1, k = j, scale = MyScale[[i]])
        temp_GOF[(i-1)*6+1+j,1] <- paste(MyScale[[i]],j,sep = ",")
        temp_GOF[(i-1)*6+1+j,2] <- glmTemp$loglik
        temp_GOF[(i-1)*6+1+j,3] <- glmTemp$npars
        temp_GOF[(i-1)*6+1+j,4] <- glmTemp$AIC
        temp_GOF[(i-1)*6+1+j,5] <- log(length(in_ipd1$time))*glmTemp$npars-2*glmTemp$loglik
      }
    }
  }
  temp_GOF<-temp_GOF %>% filter(Model != "error")
  temp_GOF<-arrange(temp_GOF,AIC)
  rp_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  rp_bic<-temp_GOF$Model[1]
  
  temp_scale<-strsplit(rp_aic,",")[[1]][1]
  temp_k<-as.numeric(strsplit(rp_aic,",")[[1]][2])
  glm_rp1<-flexsurvspline(Surv(time, event) ~ 1, data = in_ipd1, k = temp_k, scale = temp_scale)
  dfHazEst11[4,] <- summary(glm_rp1, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[rp_aic]] <- summary(glm_rp1, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[7]<-rp_aic
  dfGOF1$LnL[7]<-sum(ltHaz_1$Events*log(ltHaz_1[[rp_aic]]) - ltHaz_1[[rp_aic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[7]<-glm_rp1$npars
  dfGOF1$AIC[7]<-2*glm_rp1$npars-2*dfGOF1$LnL[7]
  dfGOF1$BIC[7]<-log(length(in_ipd1$time))*glm_rp1$npars-2*dfGOF1$LnL[7]
  
  temp_scale<-strsplit(rp_bic,",")[[1]][1]
  temp_k<-as.numeric(strsplit(rp_bic,",")[[1]][2])
  glm_rp2<-flexsurvspline(Surv(time, event) ~ 1, data = in_ipd1, k = temp_k, scale = temp_scale)
  dfHazEst12[4,] <- summary(glm_rp2, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[rp_bic]] <- summary(glm_rp2, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[8]<-rp_bic
  dfGOF1$LnL[8]<-sum(ltHaz_1$Events*log(ltHaz_1[[rp_bic]]) - ltHaz_1[[rp_bic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[8]<-glm_rp2$npars
  dfGOF1$AIC[8]<-2*glm_rp2$npars-2*dfGOF1$LnL[8]
  dfGOF1$BIC[8]<-log(length(in_ipd1$time))*glm_rp2$npars-2*dfGOF1$LnL[8]
  #-----end rp -----
  
  #----Generalised additive models (GAM)----
  bc2 <- subset(in_ipd1, event==1)
  temp_GOF<- data.frame(matrix(nrow=10, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  
  for (i in 1:5){
    fit<-try(glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time,
                                                                                                                                               seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_1))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- 'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i] <- 'error'
      temp_GOF$BIC[i] <- 'error'
    }
    else{
      gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time,
                                                                                                                                        seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_1)
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
      temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
      temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  gam_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  gam_bic<-temp_GOF$Model[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=gam_aic+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, 
                                                                                                                                          seq(from=0, to=1, by=1/(1+gam_aic))), length=gam_aic+2),family=binomial(link=cloglog), data=ltHaz_1)
  bestgam<-paste("gam_aic",gam_aic,sep = ":")
  dfHazEst11[5,] <- predict(glmTemp, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestgam] <- predict(glmTemp, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[9]<-bestgam
  dfGOF1$LnL[9]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF1$Params[9]<-extractAIC(glmTemp)[1]
  dfGOF1$AIC[9]<-extractAIC(glmTemp)[2]
  dfGOF1$BIC[9]<--2*dfGOF1$LnL[9] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=gam_bic+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, 
                                                                                                                                          seq(from=0, to=1, by=1/(1+gam_bic))), length=gam_bic+2),family=binomial(link=cloglog), data=ltHaz_1)
  bestgam<-paste("gam_bic",gam_bic,sep = ":")
  dfHazEst12[5,] <- predict(glmTemp, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestgam] <- predict(glmTemp, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[10]<-bestgam
  dfGOF1$LnL[10]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF1$Params[10]<-extractAIC(glmTemp)[1]
  dfGOF1$AIC[10]<-extractAIC(glmTemp)[2]
  dfGOF1$BIC[10]<--2*dfGOF1$LnL[10] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  #-----end gam -----
  
  #----Parametric mixture models (PMM)----
  ###----model code----####
  flexsurvmixture <- function(formula, dists, control, sr.control){
    
    n.groups <- length(dists)
    
    dist_list                 <- list()
    dist_list$name            <- paste0(n.groups,"groups: ",paste0(dists,"_", collapse=""))
    dist_list$pars            <- do.call(compile.pars, args = list(dists, n.groups))
    dist_list$location        <- "p1"
    dist_list$transforms      <- do.call(compile.transform, args = list(dists, n.groups))
    dist_list$inv.transforms  <- do.call(compile.inv.transform, args = list(dists, n.groups))
    dist_list$inits           <- function(t){do.call(compile.inits, args = list(t = t, dists, n.groups))}
    
    pfun <- dfun <- list()
    for(i in 1:n.groups){
      pfun[[i]] = get(paste0("p", dists[i]))
      dfun[[i]] = get(paste0("d", dists[i]))
    }
    
    dfns_list = list(
      p = function(q, ...) pmixsurv(pfun, q, n.groups, ...),
      d = function(x, ...) dmixsurv(dfun, x, n.groups, ...),
      mean = function(...) mean_mixsurv(pfun, ...),
      rmst = function(t, ...) rmst_mixsurv(pfun, t, ...) 
    )
    
    optim <- list()
    out <- do.call(
      "flexsurvreg",
      append(
        list(
          formula,
          dist = dist_list,
          dfns = dfns_list,
          control = control, 
          sr.control = sr.control
        ),
        optim
      )
    )
    return(out)
  }
  pmixsurv = function(pfun, q, n.groups, ...) {
    dots <- list(...)
    args <- dots
    args$lower.tail <- F
    args$log.p <- F
    
    theta                     <- c()
    theta[1]                  <- args$p1; args$p1 <- NULL
    if(n.groups > 2){theta[2] <- args$p2; args$p2 <- NULL}
    if(n.groups > 3){theta[3] <- args$p3; args$p3 <- NULL}
    theta <- multinom.p(theta)
    
    args1 <- args; args1[(endsWith(names(args1), "2") | endsWith(names(args1), "3") | endsWith(names(args1), "4"))] <- NULL
    args2 <- args; args2[(endsWith(names(args2), "1") | endsWith(names(args2), "3") | endsWith(names(args2), "4"))] <- NULL
    args3 <- args; args3[(endsWith(names(args3), "1") | endsWith(names(args3), "2") | endsWith(names(args3), "4"))] <- NULL
    args4 <- args; args4[(endsWith(names(args4), "1") | endsWith(names(args4), "2") | endsWith(names(args4), "3"))] <- NULL
    
    names(args1)[endsWith(names(args1), "1")] <- gsub("1", "", names(args1)[endsWith(names(args1), "1")])
    names(args2)[endsWith(names(args2), "2")] <- gsub("2", "", names(args2)[endsWith(names(args2), "2")])
    names(args3)[endsWith(names(args3), "3")] <- gsub("3", "", names(args3)[endsWith(names(args3), "3")])
    names(args4)[endsWith(names(args4), "4")] <- gsub("4", "", names(args4)[endsWith(names(args4), "4")])
    
    out <- (theta[1] * do.call(pfun[[1]], append(list(q), args1)) + 
              theta[2] * do.call(pfun[[2]], append(list(q), args2)) +
              ifelse(n.groups>2, theta[3] * do.call(pfun[[3]], append(list(q), args3)),0) +
              ifelse(n.groups>3, theta[4] * do.call(pfun[[4]], append(list(q), args4)),0)) 
    
    if (is.null(dots$lower.tail) || dots$lower.tail) {
      out <- 1 - out
    }
    if (!is.null(dots$log.p) && dots$log.p) {
      out <- log(out)
    }
    return(out)
  }
  dmixsurv = function(dfun, x, n.groups, ...) {
    dots <- list(...)
    args <- dots
    args$log <- F
    
    theta                     <- c()
    theta[1]                  <- args$p1; args$p1 <- NULL
    if(n.groups > 2){theta[2] <- args$p2; args$p2 <- NULL}
    if(n.groups > 3){theta[3] <- args$p3; args$p3 <- NULL}
    theta <- multinom.p(theta)
    
    args1 <- args; args1[(endsWith(names(args1), "2") | endsWith(names(args1), "3") | endsWith(names(args1), "4"))] <- NULL
    args2 <- args; args2[(endsWith(names(args2), "1") | endsWith(names(args2), "3") | endsWith(names(args2), "4"))] <- NULL
    args3 <- args; args3[(endsWith(names(args3), "1") | endsWith(names(args3), "2") | endsWith(names(args3), "4"))] <- NULL
    args4 <- args; args4[(endsWith(names(args4), "1") | endsWith(names(args4), "2") | endsWith(names(args4), "3"))] <- NULL
    
    names(args1)[endsWith(names(args1), "1")] <- gsub("1", "", names(args1)[endsWith(names(args1), "1")])
    names(args2)[endsWith(names(args2), "2")] <- gsub("2", "", names(args2)[endsWith(names(args2), "2")])
    names(args3)[endsWith(names(args3), "3")] <- gsub("3", "", names(args3)[endsWith(names(args3), "3")])
    names(args4)[endsWith(names(args4), "4")] <- gsub("4", "", names(args4)[endsWith(names(args4), "4")])
    
    
    out <- (theta[1] * do.call(dfun[[1]], append(list(x), args1)) + 
              theta[2] * do.call(dfun[[2]], append(list(x), args2)) +  
              ifelse(n.groups>2,theta[3] * do.call(dfun[[3]], append(list(x), args3)),0) + 
              ifelse(n.groups>3,theta[4] * do.call(dfun[[4]], append(list(x), args4)),0))
    
    if (!is.null(dots$log) && dots$log) {
      out <- log(out)
    }
    return(out)
  }
  rmst_mixsurv = function(pfun, q, t, ...) {
    args <- list(...)
    out <- do.call(
      rmst_generic,
      append(
        list(
          function(q, ...) pmixsurv(pfun, q, ...),
          t = t
        ),
        args
      )
    )
    return(out)
  }
  mean_mixsurv = function(pfun, ...) {
    if(p1 > 0) {
      out <- Inf
    }else {
      args <- list(...)
      out <- do.call(
        rmst_generic,
        append(
          list(
            pfun[[1]],
            t = Inf,
            start = 0
          ),
          args
        )
      )
    }
    return(out)
  }
  multinom.p = function(theta) {
    p <- rep(NA,length(theta) + 1)
    for (i in 1:length(theta)) {
      p[i] <- exp(theta[i])/(1 + sum(exp(theta)))
    }
    p[length(theta) + 1] <- 1 / (1 + sum(exp(theta)))
    return(p)
  }
  #exp
  exp.inits <- function(t) {
    1/mean(t)
  }
  exp.pars <- flexsurv.dists$exp$pars
  exp.t <- list(log)
  exp.it <- list(exp)
  #gamma
  gamma.inits <- function(t) {
    m = mean(t)
    v = var(t)
    c(m^2/v, m/v)
  }
  gamma.pars <- flexsurv.dists$gamma$pars
  gamma.t <- list(log, log)
  gamma.it <- list(exp, exp)
  #gen gamma
  gengamma.inits <- function(t) {
    lt <- log(t[t > 0])
    c(mean(lt), sd(lt), 0)
  }
  gengamma.pars <- flexsurv.dists$gengamma$pars
  gengamma.t <- list(identity, log, identity)
  gengamma.it <- list(identity, exp, identity)
  #gompertz
  gompertz.inits <- function(t) {
    c(0.001, 1/mean(t))
  }
  gompertz.pars <- flexsurv.dists$gompertz$pars
  gompertz.t <- list(identity, log)
  gompertz.it <- list(identity, exp)
  #llogis
  llogis.inits <- function(t) {
    scale <- median(t)
    shape <- 1/log(quantile(t, 0.25)/scale, base = 3)
    if (shape < 0) 
      shape <- 1
    c(shape, scale)
  }
  llogis.pars <- flexsurv.dists$llogis$pars
  llogis.t <- list(log, log)
  llogis.it <- list(exp, exp)
  #lnorm
  lnorm.inits <- function(t) {
    lt <- log(t[t > 0])
    c(mean(lt), sd(lt))
  }
  lnorm.pars <- flexsurv.dists$lnorm$pars
  lnorm.t <- list(identity, log)
  lnorm.it <- list(identity, exp)
  #weibull
  weibull.inits <- function(t) {
    lt <- log(t[t > 0])
    c(1.2/var(lt), exp(mean(lt) + 0.572)) 
  }
  weibull.pars <- flexsurv.dists$weibull$pars
  weibull.t <- list(log, log)
  weibull.it <- list(exp, exp)
  #compile
  compile.inits <- function(t, group.dists, n.groups) {
    if(n.groups >= 1) inits.groups <- rep(0.01, n.groups - 1)
    inits1 <- get(paste0(group.dists[1], ".inits"))
    if(!is.na(group.dists[2])) inits2 <- get(paste0(group.dists[2], ".inits")) else inits2 <- function(t){NULL}
    if(!is.na(group.dists[3])) inits3 <- get(paste0(group.dists[3], ".inits")) else inits3 <- function(t){NULL}
    if(!is.na(group.dists[4])) inits4 <- get(paste0(group.dists[4], ".inits")) else inits4 <- function(t){NULL}
    if(n.groups == 1) return(inits1(t))
    if(n.groups == 2) return(c(inits.groups, inits1(t), inits2(t)))
    if(n.groups == 3) return(c(inits.groups, inits1(t), inits2(t), inits3(t)))
    if(n.groups == 4) return(c(inits.groups, inits1(t), inits2(t), inits3(t), inits4(t)))
    if(n.groups >= 4) stop("too many distributions specified, use 4 groups or fewer")
  }
  compile.pars <- function(group.dists, n.groups) {
    pars.groups = NULL
    if(n.groups >= 1) for(i in 1:(n.groups-1)) pars.groups[i] <- paste0("p",i)
    pars1 <- paste0(get(paste0(group.dists[1], ".pars")),1)
    if(!is.na(group.dists[2])) pars2 <- paste0(get(paste0(group.dists[2], ".pars")),2) else pars2 <- function(){NULL}
    if(!is.na(group.dists[3])) pars3 <- paste0(get(paste0(group.dists[3], ".pars")),3) else pars3 <- function(){NULL}
    if(!is.na(group.dists[4])) pars4 <- paste0(get(paste0(group.dists[4], ".pars")),4) else pars4 <- function(){NULL}
    if(n.groups == 1) return(pars1)
    if(n.groups == 2) return(c(pars.groups, pars1, pars2))
    if(n.groups == 3) return(c(pars.groups, pars1, pars2, pars3))
    if(n.groups == 4) return(c(pars.groups, pars1, pars2, pars3, pars4))
    if(n.groups >= 4) stop("too many distributions specified, use 4 groups or fewer")
  }
  compile.transform <- function(group.dists, n.groups) {
    tgroups = list(NULL)
    if(n.groups >= 1) for(i in 1:(n.groups-1)) tgroups[[i]] <- identity
    t1 <- get(paste0(group.dists[1], ".t"))
    if(!is.na(group.dists[2])) t2 <- get(paste0(group.dists[2], ".t")) else t2 <- function(){NULL}
    if(!is.na(group.dists[3])) t3 <- get(paste0(group.dists[3], ".t")) else t3 <- function(){NULL}
    if(!is.na(group.dists[4])) t4 <- get(paste0(group.dists[4], ".t")) else t4 <- function(){NULL}
    if(n.groups == 1) return(t1)
    if(n.groups == 2) return(c(tgroups, t1, t2))
    if(n.groups == 3) return(c(tgroups, t1, t2, t3))
    if(n.groups == 4) return(c(tgroups, t1, t2, t3, t4))
    if(n.groups >= 4) stop("too many distributions specified, use 4 groups or fewer")
  }
  compile.inv.transform <- function(group.dists, n.groups) {
    itgroups = list(NULL)
    if(n.groups >= 1) for(i in 1:(n.groups-1)) itgroups[[i]] <- identity
    it1 <- get(paste0(group.dists[1], ".it"))
    if(!is.na(group.dists[2])) it2 <- get(paste0(group.dists[2], ".it")) else it2 <- function(){NULL}
    if(!is.na(group.dists[3])) it3 <- get(paste0(group.dists[3], ".it")) else it3 <- function(){NULL}
    if(!is.na(group.dists[4])) it4 <- get(paste0(group.dists[4], ".it")) else it4 <- function(){NULL}
    if(n.groups == 1) return(it1)
    if(n.groups == 2) return(c(itgroups, it1, it2))
    if(n.groups == 3) return(c(itgroups, it1, it2, it3))
    if(n.groups == 4) return(c(itgroups, it1, it2, it3, it4))
    if(n.groups >= 4) stop("too many distributions specified, use 4 groups or fewer")
  }
  ###--------####
  ## simulate results ##
  #newdata
  survgef <- Surv(time=in_ipd1$time , event=in_ipd1$event==1)
  # plot(survgef)
  
  # dists list #
  dists_MM<-as.data.frame(matrix(nrow=28,ncol=2))
  colnames(dists_MM)<-c("dist1",'dist2')
  distlist<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
  
  index1<-1
  for (i in 1:7) {
    for (j in i:7){
      dists_MM$dist1[index1]<- distlist[i]
      dists_MM$dist2[index1]<- distlist[j]
      index1<-index1+1
    }
  }
  
  temp_GOF<- data.frame(matrix(nrow=28, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  
  for (i in 1:28) {
    fit<-try(mixfit <- flexsurvmixture(survgef ~ 1, dists = c(dists_MM$dist1[i],dists_MM$dist2[i]), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8)))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i]<-paste(dists_MM$dist1[i],dists_MM$dist2[i],sep = ",")
      temp_GOF$LnL[i]<-'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i]<-'error'
      temp_GOF$BIC[i]<-'error'
    }
    else{
      flexsurvmixture(survgef ~ 1, dists = c(dists_MM$dist1[i],dists_MM$dist2[i]), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))      
      temp_GOF$Model[i]<-paste(dists_MM$dist1[i],dists_MM$dist2[i],sep = ",")
      temp_GOF$LnL[i]<-mixfit$loglik
      temp_GOF$Params[i]<-mixfit$npars
      temp_GOF$AIC[i]<-(-2*mixfit$loglik+2*mixfit$npars)
      temp_GOF$BIC[i]<-log(length(survgef))*mixfit$npars-2*mixfit$loglik
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  pmm_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  pmm_bic<-temp_GOF$Model[1]
  
  temp_dist1<-strsplit(pmm_aic,",")[[1]][1]
  temp_dist2<-strsplit(pmm_aic,",")[[1]][2]
  glm_pmm1<-flexsurvmixture(survgef ~ 1, dists = c(temp_dist1,temp_dist2), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))
  dfHazEst11[6,] <- summary(glm_pmm1, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[pmm_aic]] <- summary(glm_pmm1, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[11]<-pmm_aic
  dfGOF1$LnL[11]<-sum(ltHaz_1$Events*log(ltHaz_1[[pmm_aic]]) - ltHaz_1[[pmm_aic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[11]<-glm_pmm1$npars
  dfGOF1$AIC[11]<-2*glm_pmm1$npars-2*dfGOF1$LnL[11]
  dfGOF1$BIC[11]<-log(length(in_ipd1$time))*glm_pmm1$npars-2*dfGOF1$LnL[11]
  
  temp_dist1<-strsplit(pmm_bic,",")[[1]][1]
  temp_dist2<-strsplit(pmm_bic,",")[[1]][2]
  glm_pmm2<-flexsurvmixture(survgef ~ 1, dists = c(temp_dist1,temp_dist2), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))
  dfHazEst12[6,] <- summary(glm_pmm2, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[pmm_bic]] <- summary(glm_pmm2, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[12]<-pmm_bic
  dfGOF1$LnL[12]<-sum(ltHaz_1$Events*log(ltHaz_1[[pmm_bic]]) - ltHaz_1[[pmm_bic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[12]<-glm_pmm2$npars
  dfGOF1$AIC[12]<-2*glm_pmm2$npars-2*dfGOF1$LnL[12]
  dfGOF1$BIC[12]<-log(length(in_ipd1$time))*glm_pmm2$npars-2*dfGOF1$LnL[12]
  #-----end pmm -----
  
  ###-----------###
  #Extrapolation-2
  ###-----------###
  #----Standard survival model (SSM)----
  temp_GOF<- data.frame(matrix(nrow=7, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  MyDists <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
  for (i in 1:7){
    glmTemp <- flexsurvreg(Surv(time, event) ~ 1, data = in_ipd2, dist = MyDists[[i]])
    temp_GOF[i,1] <- MyDists[[i]]
    temp_GOF[i,2] <- glmTemp$loglik
    temp_GOF[i,3] <- glmTemp$npars
    temp_GOF[i,4] <- glmTemp$AIC
    temp_GOF[i,5] <- log(length(in_ipd2$time))*glmTemp$npars-2*glmTemp$loglik
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  ssm_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  ssm_bic<-temp_GOF$Model[1]
  
  glm_ssm1<-flexsurvreg(Surv(time, event) ~ 1, data = in_ipd2, dist = ssm_aic)
  dfHazEst21[1,] <- summary(glm_ssm1, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[ssm_aic]] <- summary(glm_ssm1, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[1]<-ssm_aic
  dfGOF2$LnL[1]<-sum(ltHaz_2$Events*log(ltHaz_2[[ssm_aic]]) - ltHaz_2[[ssm_aic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[1]<-glm_ssm1$npars
  dfGOF2$AIC[1]<-2*glm_ssm1$npars-2*dfGOF2$LnL[1]
  dfGOF2$BIC[1]<-log(length(in_ipd2$time))*glm_ssm1$npars-2*dfGOF2$LnL[1]
  
  glm_ssm2<-flexsurvreg(Surv(time, event) ~ 1, data = in_ipd2, dist = ssm_bic)
  dfHazEst22[1,] <- summary(glm_ssm2, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[ssm_bic]] <- summary(glm_ssm2, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[2]<-ssm_bic
  dfGOF2$LnL[2]<-sum(ltHaz_2$Events*log(ltHaz_2[[ssm_bic]]) - ltHaz_2[[ssm_bic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[2]<-glm_ssm2$npars
  dfGOF2$AIC[2]<-2*glm_ssm2$npars-2*dfGOF2$LnL[2]
  dfGOF2$BIC[2]<-log(length(in_ipd2$time))*glm_ssm2$npars-2*dfGOF2$LnL[2]
  #-----end ssm -----
  
  #----Fractional polynomial (FP)----
  temp_GOF<- data.frame(matrix(nrow=44, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  #-----FP1 -----
  MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
  for (i in 1:7){
    glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    temp_GOF$Model[i] <- paste(MyPowers[[1]][i],"NA",sep = ",") 
    temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  }
  ### run for 0
  glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  temp_GOF$Model[8] <- paste(0,"NA",sep = ",")
  temp_GOF$LnL[8] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  temp_GOF$Params[8]<-extractAIC(glmTemp)[1]
  temp_GOF$AIC[8] <- extractAIC(glmTemp)[2]
  temp_GOF$BIC[8] <- -2*temp_GOF$LnL[8] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  
  #-----FP2 -----
  myLnL <- array(dim=36)
  myAIC <- array(dim=36)
  MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
  index <- 9
  for (i in 1:7){
    for (j in 1:7){
      if (j > i) {
        glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][j])+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2) 
        temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
        temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
        temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
        temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
        temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
        index <- index + 1
      }
    }
  }
  for (i in 1:7) {
    glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][i]*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)# 
    temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
    temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
    index <- index + 1
    
  }
  
  for (i in 1:7) {
    glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)# 
    temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
    temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
    index <- index + 1
  }
  
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)# 
  temp_GOF$Model[index] <- paste(0,0,sep = ",")
  temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
  temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
  temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  
  temp_GOF<-arrange(temp_GOF,AIC)
  fp_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  fp_bic<-temp_GOF$Model[1]
  
  FP_pow_1<-as.numeric((strsplit(fp_aic,","))[[1]][1])
  FP_pow_2<-as.numeric((strsplit(fp_aic,","))[[1]][2])
  if(is.na(FP_pow_2)){
    if(FP_pow_1 == 0){
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    } else {
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    }
  } else if (FP_pow_1 == FP_pow_2 & FP_pow_2 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else if (FP_pow_1==0 & FP_pow_2== 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else if (FP_pow_2==0 & FP_pow_1 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else {
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  }
  bestfp<-paste("FP_aic",fp_aic,sep = ":")
  dfHazEst21[2,] <- predict(modFP, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestfp] <- predict(modFP, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[3]<-bestfp
  dfGOF2$LnL[3]<-(extractAIC(modFP)[2] - 2*extractAIC(modFP)[1])*(-0.5)
  dfGOF2$Params[3]<-extractAIC(modFP)[1]
  dfGOF2$AIC[3]<-extractAIC(modFP)[2]
  dfGOF2$BIC[3]<--2*dfGOF2$LnL[3] + log(length(ltHaz_2$Time))*extractAIC(modFP)[1]
  
  FP_pow_1<-as.numeric((strsplit(fp_bic,","))[[1]][1])
  FP_pow_2<-as.numeric((strsplit(fp_bic,","))[[1]][2])
  if(is.na(FP_pow_2)){
    if(FP_pow_1 == 0){
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    } else {
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    }
  } else if (FP_pow_1 == FP_pow_2 & FP_pow_2 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else if (FP_pow_1==0 & FP_pow_2== 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else if (FP_pow_2==0 & FP_pow_1 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else {
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  }
  bestfp<-paste("FP_bic",fp_bic,sep = ":")
  dfHazEst22[2,] <- predict(modFP, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestfp] <- predict(modFP, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[4]<-bestfp
  dfGOF2$LnL[4]<-(extractAIC(modFP)[2] - 2*extractAIC(modFP)[1])*(-0.5)
  dfGOF2$Params[4]<-extractAIC(modFP)[1]
  dfGOF2$AIC[4]<-extractAIC(modFP)[2]
  dfGOF2$BIC[4]<--2*dfGOF2$LnL[4] + log(length(ltHaz_2$Time))*extractAIC(modFP)[1]
  #-----end fp -----
  
  #----Restricted cubic splines (RCS)----
  bc2 <- subset(in_ipd2, event==1)
  temp_GOF<- data.frame(matrix(nrow=5, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  for (i in 1:5){
    fit<-try(glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_2))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- 'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i] <- 'error'
      temp_GOF$BIC[i] <- 'error'
    }
    else{
      gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_2)      
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
      temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
      temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  rcs_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  rcs_bic<-temp_GOF$Model[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=rcs_aic+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+rcs_aic))), length=rcs_aic+2),family=binomial(link=cloglog), data=ltHaz_2)
  bestrcs<-paste("RCS_aic",rcs_aic,sep = ":")
  dfHazEst21[3,] <- predict(glmTemp, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestrcs] <- predict(glmTemp, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[5]<-bestrcs
  dfGOF2$LnL[5]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF2$Params[5]<-extractAIC(glmTemp)[1]
  dfGOF2$AIC[5]<-extractAIC(glmTemp)[2]
  dfGOF2$BIC[5]<--2*dfGOF2$LnL[5] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=rcs_bic+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+rcs_bic))), length=rcs_bic+2),family=binomial(link=cloglog), data=ltHaz_2)
  bestrcs<-paste("RCS_bic",rcs_bic,sep = ":")
  dfHazEst22[3,] <- predict(glmTemp, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestrcs] <- predict(glmTemp, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[6]<-bestrcs
  dfGOF2$LnL[6]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF2$Params[6]<-extractAIC(glmTemp)[1]
  dfGOF2$AIC[6]<-extractAIC(glmTemp)[2]
  dfGOF2$BIC[6]<--2*dfGOF2$LnL[6] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  #-----end rcs -----
  
  #----Royston-Parmar models (RP)----
  temp_GOF<- data.frame(matrix(nrow=18, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  MyScale <- list("hazard","odds","normal")
  for (i in 1:3){
    for (j in 0:5){
      fit<-try(glmTemp <- flexsurvspline(Surv(time, event) ~ 1, data = in_ipd2, k = j, scale = MyScale[[i]]))
      if("try-error" %in% class(fit)) {
        temp_GOF[(i-1)*6+1+j,1] <- "error"
        temp_GOF[(i-1)*6+1+j,2] <- "error"
        temp_GOF[(i-1)*6+1+j,3] <- "error"
        temp_GOF[(i-1)*6+1+j,4] <- "error"
        temp_GOF[(i-1)*6+1+j,5] <- "error"
      }
      else{
        flexsurvspline(Surv(time, event) ~ 1, data = in_ipd2, k = j, scale = MyScale[[i]])
        temp_GOF[(i-1)*6+1+j,1] <- paste(MyScale[[i]],j,sep = ",")
        temp_GOF[(i-1)*6+1+j,2] <- glmTemp$loglik
        temp_GOF[(i-1)*6+1+j,3] <- glmTemp$npars
        temp_GOF[(i-1)*6+1+j,4] <- glmTemp$AIC
        temp_GOF[(i-1)*6+1+j,5] <- log(length(in_ipd2$time))*glmTemp$npars-2*glmTemp$loglik
      }
    }
  }
  temp_GOF<-temp_GOF %>% filter(Model != "error")
  temp_GOF<-arrange(temp_GOF,AIC)
  rp_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  rp_bic<-temp_GOF$Model[1]
  
  temp_scale<-strsplit(rp_aic,",")[[1]][1]
  temp_k<-as.numeric(strsplit(rp_aic,",")[[1]][2])
  glm_rp1<-flexsurvspline(Surv(time, event) ~ 1, data = in_ipd2, k = temp_k, scale = temp_scale)
  dfHazEst21[4,] <- summary(glm_rp1, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[rp_aic]] <- summary(glm_rp1, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[7]<-rp_aic
  dfGOF2$LnL[7]<-sum(ltHaz_2$Events*log(ltHaz_2[[rp_aic]]) - ltHaz_2[[rp_aic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[7]<-glm_rp1$npars
  dfGOF2$AIC[7]<-2*glm_rp1$npars-2*dfGOF2$LnL[7]
  dfGOF2$BIC[7]<-log(length(in_ipd2$time))*glm_rp1$npars-2*dfGOF2$LnL[7]
  
  temp_scale<-strsplit(rp_bic,",")[[1]][1]
  temp_k<-as.numeric(strsplit(rp_bic,",")[[1]][2])
  glm_rp2<-flexsurvspline(Surv(time, event) ~ 1, data = in_ipd2, k = temp_k, scale = temp_scale)
  dfHazEst22[4,] <- summary(glm_rp2, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[rp_bic]] <- summary(glm_rp2, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[8]<-rp_bic
  dfGOF2$LnL[8]<-sum(ltHaz_2$Events*log(ltHaz_2[[rp_bic]]) - ltHaz_2[[rp_bic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[8]<-glm_rp2$npars
  dfGOF2$AIC[8]<-2*glm_rp2$npars-2*dfGOF2$LnL[8]
  dfGOF2$BIC[8]<-log(length(in_ipd2$time))*glm_rp2$npars-2*dfGOF2$LnL[8]
  #-----end rp -----
  
  #----Generalised additive models (GAM)----
  bc2 <- subset(in_ipd2, event==1)
  temp_GOF<- data.frame(matrix(nrow=10, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  for (i in 1:5){
    fit<-try(glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time,
                                                                                                                                               seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_2))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- 'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i] <- 'error'
      temp_GOF$BIC[i] <- 'error'
    }
    else{
      gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time,
                                                                                                                             seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_2)
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
      temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
      temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  gam_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  gam_bic<-temp_GOF$Model[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=gam_aic+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, 
                                                                                                                                          seq(from=0, to=1, by=1/(1+gam_aic))), length=gam_aic+2),family=binomial(link=cloglog), data=ltHaz_2)
  bestgam<-paste("gam_aic",gam_aic,sep = ":")
  dfHazEst21[5,] <- predict(glmTemp, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestgam] <- predict(glmTemp, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[9]<-bestgam
  dfGOF2$LnL[9]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF2$Params[9]<-extractAIC(glmTemp)[1]
  dfGOF2$AIC[9]<-extractAIC(glmTemp)[2]
  dfGOF2$BIC[9]<--2*dfGOF2$LnL[9] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=gam_bic+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, 
                                                                                                                                          seq(from=0, to=1, by=1/(1+gam_bic))), length=gam_bic+2),family=binomial(link=cloglog), data=ltHaz_2)
  bestgam<-paste("gam_bic",gam_bic,sep = ":")
  dfHazEst22[5,] <- predict(glmTemp, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestgam] <- predict(glmTemp, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[10]<-bestgam
  dfGOF2$LnL[10]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF2$Params[10]<-extractAIC(glmTemp)[1]
  dfGOF2$AIC[10]<-extractAIC(glmTemp)[2]
  dfGOF2$BIC[10]<--2*dfGOF2$LnL[10] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  #-----end gam -----
  
  #----Parametric mixture models (PMM)----
  ###--------####
  ## simulate results ##
  #newdata
  survgef <- Surv(time=in_ipd2$time , event=in_ipd2$event==1)
  # plot(survgef)
  
  # dists list #
  dists_MM<-as.data.frame(matrix(nrow=28,ncol=2))
  colnames(dists_MM)<-c("dist1",'dist2')
  distlist<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
  
  index1<-1
  for (i in 1:7) {
    for (j in i:7){
      dists_MM$dist1[index1]<- distlist[i]
      dists_MM$dist2[index1]<- distlist[j]
      index1<-index1+1
    }
  }
  
  temp_GOF<- data.frame(matrix(nrow=28, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  
  for (i in 1:28) {
    fit<-try(mixfit <- flexsurvmixture(survgef ~ 1, dists = c(dists_MM$dist1[i],dists_MM$dist2[i]), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8)))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i]<-paste(dists_MM$dist1[i],dists_MM$dist2[i],sep = ",")
      temp_GOF$LnL[i]<-'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i]<-'error'
      temp_GOF$BIC[i]<-'error'
    }
    else{
      flexsurvmixture(survgef ~ 1, dists = c(dists_MM$dist1[i],dists_MM$dist2[i]), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))      
      temp_GOF$Model[i]<-paste(dists_MM$dist1[i],dists_MM$dist2[i],sep = ",")
      temp_GOF$LnL[i]<-mixfit$loglik
      temp_GOF$Params[i]<-mixfit$npars
      temp_GOF$AIC[i]<-(-2*mixfit$loglik+2*mixfit$npars)
      temp_GOF$BIC[i]<-log(length(survgef))*mixfit$npars-2*mixfit$loglik
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  pmm_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  pmm_bic<-temp_GOF$Model[1]
  
  temp_dist1<-strsplit(pmm_aic,",")[[1]][1]
  temp_dist2<-strsplit(pmm_aic,",")[[1]][2]
  glm_pmm1<-flexsurvmixture(survgef ~ 1, dists = c(temp_dist1,temp_dist2), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))
  dfHazEst21[6,] <- summary(glm_pmm1, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[pmm_aic]] <- summary(glm_pmm1, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[11]<-pmm_aic
  dfGOF2$LnL[11]<-sum(ltHaz_2$Events*log(ltHaz_2[[pmm_aic]]) - ltHaz_2[[pmm_aic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[11]<-glm_pmm1$npars
  dfGOF2$AIC[11]<-2*glm_pmm1$npars-2*dfGOF2$LnL[11]
  dfGOF2$BIC[11]<-log(length(in_ipd2$time))*glm_pmm1$npars-2*dfGOF2$LnL[11]
  
  temp_dist1<-strsplit(pmm_bic,",")[[1]][1]
  temp_dist2<-strsplit(pmm_bic,",")[[1]][2]
  glm_pmm2<-flexsurvmixture(survgef ~ 1, dists = c(temp_dist1,temp_dist2), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))
  dfHazEst22[6,] <- summary(glm_pmm2, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[pmm_bic]] <- summary(glm_pmm2, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[12]<-pmm_bic
  dfGOF2$LnL[12]<-sum(ltHaz_2$Events*log(ltHaz_2[[pmm_bic]]) - ltHaz_2[[pmm_bic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[12]<-glm_pmm2$npars
  dfGOF2$AIC[12]<-2*glm_pmm2$npars-2*dfGOF2$LnL[12]
  dfGOF2$BIC[12]<-log(length(in_ipd2$time))*glm_pmm2$npars-2*dfGOF2$LnL[12]
  #-----end pmm -----
  dfHaz_1aic <-t(dfHazEst11)
  colnames(dfHaz_1aic) <- c("SSM","FP","RCS","RP","GAM","PMM")
  dfHaz_1aic <- cbind(data.frame(time=ltHaz_5$Time), dfHaz_1aic)
  dfSurv_1aic<-as.data.frame(matrix(nrow = length(dfHaz_1aic$time),ncol = 7))
  dfSurv_1aic[[1]]<-dfHaz_1aic$time
  for ( i in 2:7 ) {
    dftemp<-data.frame(dfHaz_1aic$time,dfHaz_1aic[[i]])
    dftemp<-dftemp %>% 
      dplyr::arrange(dftemp[[1]]) %>% 
      dplyr::mutate(cumhaz = cumsum(dftemp[[2]])) %>% 
      dplyr::mutate(survProp = exp(-1*cumhaz))
    dfSurv_1aic[[i]]<-dftemp[[4]]
  }
  colnames(dfSurv_1aic) <- c("Time","SSM","FP","RCS","RP","GAM","PMM")
  
  dfHaz_1bic <-t(dfHazEst12)
  colnames(dfHaz_1bic) <- c("SSM","FP","RCS","RP","GAM","PMM")
  dfHaz_1bic <- cbind(data.frame(time=ltHaz_5$Time), dfHaz_1bic)
  dfSurv_1bic<-as.data.frame(matrix(nrow = length(dfHaz_1bic$time),ncol = 7))
  dfSurv_1bic[[1]]<-dfHaz_1bic$time
  for ( i in 2:7 ) {
    dftemp<-data.frame(dfHaz_1bic$time,dfHaz_1bic[[i]])
    dftemp<-dftemp %>% 
      dplyr::arrange(dftemp[[1]]) %>% 
      dplyr::mutate(cumhaz = cumsum(dftemp[[2]])) %>% 
      dplyr::mutate(survProp = exp(-1*cumhaz))
    dfSurv_1bic[[i]]<-dftemp[[4]]
  }
  colnames(dfSurv_1bic) <- c("Time","SSM","FP","RCS","RP","GAM","PMM")

  dfHaz_2aic <-t(dfHazEst21)
  colnames(dfHaz_2aic) <- c("SSM","FP","RCS","RP","GAM","PMM")
  dfHaz_2aic <- cbind(data.frame(time=ltHaz_5$Time), dfHaz_2aic)
  dfSurv_2aic<-as.data.frame(matrix(nrow = length(dfHaz_2aic$time),ncol = 7))
  dfSurv_2aic[[1]]<-dfHaz_2aic$time
  for ( i in 2:7 ) {
    dftemp<-data.frame(dfHaz_2aic$time,dfHaz_2aic[[i]])
    dftemp<-dftemp %>% 
      dplyr::arrange(dftemp[[1]]) %>% 
      dplyr::mutate(cumhaz = cumsum(dftemp[[2]])) %>% 
      dplyr::mutate(survProp = exp(-1*cumhaz))
    dfSurv_2aic[[i]]<-dftemp[[4]]
  }
  colnames(dfSurv_2aic) <- c("Time","SSM","FP","RCS","RP","GAM","PMM")
  
  dfHaz_2bic <-t(dfHazEst22)
  colnames(dfHaz_2bic) <- c("SSM","FP","RCS","RP","GAM","PMM")
  dfHaz_2bic <- cbind(data.frame(time=ltHaz_5$Time), dfHaz_2bic)
  dfSurv_2bic<-as.data.frame(matrix(nrow = length(dfHaz_2bic$time),ncol = 7))
  dfSurv_2bic[[1]]<-dfHaz_2bic$time
  for ( i in 2:7 ) {
    dftemp<-data.frame(dfHaz_2bic$time,dfHaz_2bic[[i]])
    dftemp<-dftemp %>% 
      dplyr::arrange(dftemp[[1]]) %>% 
      dplyr::mutate(cumhaz = cumsum(dftemp[[2]])) %>% 
      dplyr::mutate(survProp = exp(-1*cumhaz))
    dfSurv_2bic[[i]]<-dftemp[[4]]
  }
  colnames(dfSurv_2bic) <- c("Time","SSM","FP","RCS","RP","GAM","PMM")
  
  output=list(dfSurv_1aic=dfSurv_1aic,
              dfSurv_1bic=dfSurv_1bic,
              dfSurv_2aic=dfSurv_2aic,
              dfSurv_2bic=dfSurv_2bic,
              dfGOF1=dfGOF1,
              dfGOF2=dfGOF2,
              ltHaz_1=ltHaz_1,
              ltHaz_2=ltHaz_2,
              ltHaz_5=ltHaz_5
              )
  return(output)
}


####-------------------------------------------------------------------------------------------------------------------------------------------------------####
####-------------------------------------------------------------------------------------------------------------------------------------------------------####
####-------------------------------------------------------------------------------------------------------------------------------------------------------####


Surv_analysis2<-function(data){
  data<-data[-(1:4),]
  na_cols <- colSums(is.na(data)) == nrow(data)
  data <- data[, !na_cols]
  colnames(data)<-c("Time","Event","Time","Event","Time","Event","Time","Event","Time","Event","Time","Event")
 
  df1<-na.omit(data[,9:10])
  df1$Time<-as.numeric(df1$Time)
  df1$Event<-as.numeric(df1$Event)
  

  df2<-na.omit(data[,7:8])
  df2$Time<-as.numeric(df2$Time)
  df2$Event<-as.numeric(df2$Event)
  

  df5<-na.omit(data[,11:12])
  df5$Time<-as.numeric(df5$Time)
  df5$Event<-as.numeric(df5$Event)
  
  ##################################################
  #############   generate  results   ##############
  ##################################################
  # extrapolation 1
  trt1<-as.data.frame(rename(df1,"censrec"="Event","recyrs"="Time"))
  trt1$recyrs<-trt1$recyrs/12
  trt1$rectime<-trt1$recyrs*365
  trt1$rectime2 <- as.integer(trt1$rectime/(365.24/12)) + 1
  lttrt1 <- lifeTable(trt1, timeColumn = "rectime2", eventColumn = "censrec")
  ltHaz_1 <- data.frame(hazKM = lttrt1$Output$hazard, Time = (seq(1:length(lttrt1$Output[,1]))-0.5)/12,
                        AtRisk = lttrt1$Output$atRisk, Events = lttrt1$Output$events)
  ltHaz_1$hazLT = ltHaz_1$Events / (ltHaz_1$AtRisk - ltHaz_1$Events/2)
  ltHaz_1$lnTime <- log(ltHaz_1$Time)
  ltHaz_1$MyId <- 1:dim(ltHaz_1)[1] # Generate id variable
  ltHaz_1$EventsL <- lag(ltHaz_1$Events)
  ltHaz_1$EventsL[1] <- 0
  ltHaz_1$surv <- lttrt1$Output$S
  ltHaz_1$timedelta<-ltHaz_1$Time[2]-ltHaz_1$Time[1]
  in_ipd1<-data.frame("time"=trt1$recyrs,"event"=trt1$censrec)
  
  # extrapolation 2
  trt2<-as.data.frame(rename(df2,"censrec"="Event","recyrs"="Time"))
  trt2$recyrs<-trt2$recyrs/12
  trt2$rectime<-trt2$recyrs*365
  trt2$rectime2 <- as.integer(trt2$rectime/(365.24/12)) + 1
  lttrt2 <- lifeTable(trt2, timeColumn = "rectime2", eventColumn = "censrec")
  ltHaz_2 <- data.frame(hazKM = lttrt2$Output$hazard, Time = (seq(1:length(lttrt2$Output[,1]))-0.5)/12,
                        AtRisk = lttrt2$Output$atRisk, Events = lttrt2$Output$events)
  ltHaz_2$hazLT = ltHaz_2$Events / (ltHaz_2$AtRisk - ltHaz_2$Events/2)
  ltHaz_2$lnTime <- log(ltHaz_2$Time)
  ltHaz_2$MyId <- 1:dim(ltHaz_2)[1] # Generate id variable
  ltHaz_2$EventsL <- lag(ltHaz_2$Events)
  ltHaz_2$EventsL[1] <- 0
  ltHaz_2$surv <- lttrt2$Output$S
  ltHaz_2$timedelta<-ltHaz_2$Time[2]-ltHaz_2$Time[1]
  in_ipd2<-data.frame("time"=trt2$recyrs,"event"=trt2$censrec)
  
  # extrapolation 5
  trt5<-as.data.frame(rename(df5,"censrec"="Event","recyrs"="Time"))
  trt5$recyrs<-trt5$recyrs/12
  trt5$rectime<-trt5$recyrs*365
  trt5$rectime2 <- as.integer(trt5$rectime/(365.24/12)) + 1
  lttrt5 <- lifeTable(trt5, timeColumn = "rectime2", eventColumn = "censrec")
  ltHaz_5 <- data.frame(hazKM = lttrt5$Output$hazard, Time = (seq(1,length(lttrt5$Output[,1]))-0.5)/12,
                        AtRisk = lttrt5$Output$atRisk, Events = lttrt5$Output$events)
  ltHaz_5$hazLT = ltHaz_5$Events / (ltHaz_5$AtRisk - ltHaz_5$Events/2)
  ltHaz_5$lnTime <- log(ltHaz_5$Time)
  ltHaz_5$MyId <- 1:dim(ltHaz_5)[1] # Generate id variable
  ltHaz_5$EventsL <- lag(ltHaz_5$Events)
  ltHaz_5$EventsL[1] <- 0
  ltHaz_5$surv <- lttrt5$Output$S
  ltHaz_5$timedelta<-ltHaz_5$Time[2]-ltHaz_5$Time[1]
  in_ipd5<-data.frame("time"=trt5$recyrs,"event"=trt5$censrec)
  
  #New Time Data 1
  follow_up <- round(ltHaz_1$Time[length(ltHaz_1$Time)]*12,0) #origin study (Month)
  numMod <- 6 # Models considered
  dfHazEst11 <- array(dim=c(numMod, length(ltHaz_5$Time)))
  dfHazEst12 <- array(dim=c(numMod, length(ltHaz_5$Time)))
  Newtime_1 <- data.frame(Time = ltHaz_5$Time, AtRisk = 1)
  Newtime_1$MyId <- 1:dim(Newtime_1)[1]
  Newtime_1$MyId <- ifelse(Newtime_1$MyId > follow_up, follow_up, Newtime_1$MyId)  # Random effects: Using last observed ID for extrapolation
  Newtime_1$EventsL <- 0
  Newtime_1$EventsL[1:follow_up] <- lag(ltHaz_1$Events)
  Newtime_1$EventsL[1] <- 0
  Newtime_1$EventsL <- ifelse(Newtime_1$MyId > follow_up, 0, Newtime_1$EventsL) # AR: Using last observed event count for extrapolation
  Newtime_1$timedelta<-Newtime_1$Time[2]-Newtime_1$Time[1]
  Newtime_1$lnTime<-log(Newtime_1$Time)
  dfGOF1 <- data.frame(matrix(nrow=12, ncol=5))
  colnames(dfGOF1) <- c("Model","LnL","Params","AIC","BIC")
  # Below is constant for when have to derive log-likelihood
  llCons1 <- sum(ltHaz_1$Events*log(ltHaz_1$AtRisk) - log(factorial(ltHaz_1$Events)))
  
  #New Time Data 2
  follow_up <- round(ltHaz_2$Time[length(ltHaz_2$Time)]*12,0) #origin study (Month)
  numMod <- 6 # Models considered
  dfHazEst21 <- array(dim=c(numMod, length(ltHaz_5$Time)))
  dfHazEst22 <- array(dim=c(numMod, length(ltHaz_5$Time)))
  Newtime_2 <- data.frame(Time = ltHaz_5$Time, AtRisk = 1)
  Newtime_2$MyId <- 1:dim(Newtime_2)[1]
  Newtime_2$MyId <- ifelse(Newtime_2$MyId > follow_up, follow_up, Newtime_2$MyId)  # Random effects: Using last observed ID for extrapolation
  Newtime_2$EventsL <- 0
  Newtime_2$EventsL[1:follow_up] <- lag(ltHaz_2$Events)
  Newtime_2$EventsL[1] <- 0
  Newtime_2$EventsL <- ifelse(Newtime_2$MyId > follow_up, 0, Newtime_2$EventsL) # AR: Using last observed event count for extrapolation
  Newtime_2$timedelta<-Newtime_2$Time[2]-Newtime_2$Time[1]
  Newtime_2$lnTime<-log(Newtime_2$Time)
  dfGOF2 <- data.frame(matrix(nrow=12, ncol=5))
  colnames(dfGOF2) <- c("Model","LnL","Params","AIC","BIC")
  # Below is constant for when have to derive log-likelihood
  llCons2 <- sum(ltHaz_2$Events*log(ltHaz_2$AtRisk) - log(factorial(ltHaz_2$Events)))
  
  ###-----------###
  #Extrapolation-1
  ###-----------###
  #----Standard survival model (SSM)----
  temp_GOF<- data.frame(matrix(nrow=7, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  MyDists <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
  for (i in 1:7){
    glmTemp <- flexsurvreg(Surv(time, event) ~ 1, data = in_ipd1, dist = MyDists[[i]])
    temp_GOF[i,1] <- MyDists[[i]]
    temp_GOF[i,2] <- glmTemp$loglik
    temp_GOF[i,3] <- glmTemp$npars
    temp_GOF[i,4] <- glmTemp$AIC
    temp_GOF[i,5] <- log(length(in_ipd1$time))*glmTemp$npars-2*glmTemp$loglik
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  ssm_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  ssm_bic<-temp_GOF$Model[1]
  
  glm_ssm1<-flexsurvreg(Surv(time, event) ~ 1, data = in_ipd1, dist = ssm_aic)
  dfHazEst11[1,] <- summary(glm_ssm1, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[ssm_aic]] <- summary(glm_ssm1, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[1]<-ssm_aic
  dfGOF1$LnL[1]<-sum(ltHaz_1$Events*log(ltHaz_1[[ssm_aic]]) - ltHaz_1[[ssm_aic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[1]<-glm_ssm1$npars
  dfGOF1$AIC[1]<-2*glm_ssm1$npars-2*dfGOF1$LnL[1]
  dfGOF1$BIC[1]<-log(length(in_ipd1$time))*glm_ssm1$npars-2*dfGOF1$LnL[1]
  
  glm_ssm2<-flexsurvreg(Surv(time, event) ~ 1, data = in_ipd1, dist = ssm_bic)
  dfHazEst12[1,] <- summary(glm_ssm2, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[ssm_bic]] <- summary(glm_ssm2, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[2]<-ssm_bic
  dfGOF1$LnL[2]<-sum(ltHaz_1$Events*log(ltHaz_1[[ssm_bic]]) - ltHaz_1[[ssm_bic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[2]<-glm_ssm2$npars
  dfGOF1$AIC[2]<-2*glm_ssm2$npars-2*dfGOF1$LnL[2]
  dfGOF1$BIC[2]<-log(length(in_ipd1$time))*glm_ssm2$npars-2*dfGOF1$LnL[2]
  #-----end ssm -----
  
  #----Fractional polynomial (FP)----
  temp_GOF<- data.frame(matrix(nrow=44, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  #-----FP1 -----
  MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
  for (i in 1:7){
    glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    temp_GOF$Model[i] <- paste(MyPowers[[1]][i],"NA",sep = ",") 
    temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  }
  ### run for 0
  glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  temp_GOF$Model[8] <- paste(0,"NA",sep = ",")
  temp_GOF$LnL[8] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  temp_GOF$Params[8]<-extractAIC(glmTemp)[1]
  temp_GOF$AIC[8] <- extractAIC(glmTemp)[2]
  temp_GOF$BIC[8] <- -2*temp_GOF$LnL[8] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  
  #-----FP2 -----
  myLnL <- array(dim=36)
  myAIC <- array(dim=36)
  MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
  index <- 9
  for (i in 1:7){
    for (j in 1:7){
      if (j > i) {
        glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][j])+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1) 
        temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
        temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
        temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
        temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
        temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
        index <- index + 1
      }
    }
  }
  for (i in 1:7) {
    glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][i]*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)# 
    temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
    temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    index <- index + 1
    
  }
  
  for (i in 1:7) {
    glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)# 
    temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
    temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    index <- index + 1
  }
  
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)# 
  temp_GOF$Model[index] <- paste(0,0,sep = ",")
  temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
  temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
  temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  
  temp_GOF<-arrange(temp_GOF,AIC)
  fp_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  fp_bic<-temp_GOF$Model[1]
  
  FP_pow_1<-as.numeric((strsplit(fp_aic,","))[[1]][1])
  FP_pow_2<-as.numeric((strsplit(fp_aic,","))[[1]][2])
  if(is.na(FP_pow_2)){
    if(FP_pow_1 == 0){
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    } else {
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    }
  } else if (FP_pow_1 == FP_pow_2 & FP_pow_2 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else if (FP_pow_1==0 & FP_pow_2== 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else if (FP_pow_2==0 & FP_pow_1 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else {
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  }
  bestfp<-paste("FP_aic",fp_aic,sep = ":")
  dfHazEst11[2,] <- predict(modFP, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestfp] <- predict(modFP, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[3]<-bestfp
  dfGOF1$LnL[3]<-(extractAIC(modFP)[2] - 2*extractAIC(modFP)[1])*(-0.5)
  dfGOF1$Params[3]<-extractAIC(modFP)[1]
  dfGOF1$AIC[3]<-extractAIC(modFP)[2]
  dfGOF1$BIC[3]<--2*dfGOF1$LnL[3] + log(length(ltHaz_1$Time))*extractAIC(modFP)[1]
  
  FP_pow_1<-as.numeric((strsplit(fp_bic,","))[[1]][1])
  FP_pow_2<-as.numeric((strsplit(fp_bic,","))[[1]][2])
  if(is.na(FP_pow_2)){
    if(FP_pow_1 == 0){
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    } else {
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
    }
  } else if (FP_pow_1 == FP_pow_2 & FP_pow_2 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else if (FP_pow_1==0 & FP_pow_2== 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else if (FP_pow_2==0 & FP_pow_1 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  } else {
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_1)
  }
  bestfp<-paste("FP_bic",fp_bic,sep = ":")
  dfHazEst12[2,] <- predict(modFP, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestfp] <- predict(modFP, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[4]<-bestfp
  dfGOF1$LnL[4]<-(extractAIC(modFP)[2] - 2*extractAIC(modFP)[1])*(-0.5)
  dfGOF1$Params[4]<-extractAIC(modFP)[1]
  dfGOF1$AIC[4]<-extractAIC(modFP)[2]
  dfGOF1$BIC[4]<--2*dfGOF1$LnL[4] + log(length(ltHaz_1$Time))*extractAIC(modFP)[1]
  #-----end fp -----
  
  #----Restricted cubic splines (RCS)----
  bc2 <- subset(in_ipd1, event==1)
  temp_GOF<- data.frame(matrix(nrow=5, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  for (i in 1:5){
    fit<-try(glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_1))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- 'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i] <- 'error'
      temp_GOF$BIC[i] <- 'error'
    }
    else{
      gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_1)      
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
      temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
      temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  rcs_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  rcs_bic<-temp_GOF$Model[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=rcs_aic+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+rcs_aic))), length=rcs_aic+2),family=binomial(link=cloglog), data=ltHaz_1)
  bestrcs<-paste("RCS_aic",rcs_aic,sep = ":")
  dfHazEst11[3,] <- predict(glmTemp, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestrcs] <- predict(glmTemp, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[5]<-bestrcs
  dfGOF1$LnL[5]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF1$Params[5]<-extractAIC(glmTemp)[1]
  dfGOF1$AIC[5]<-extractAIC(glmTemp)[2]
  dfGOF1$BIC[5]<--2*dfGOF1$LnL[5] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=rcs_bic+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+rcs_bic))), length=rcs_bic+2),family=binomial(link=cloglog), data=ltHaz_1)
  bestrcs<-paste("RCS_bic",rcs_bic,sep = ":")
  dfHazEst12[3,] <- predict(glmTemp, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestrcs] <- predict(glmTemp, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[6]<-bestrcs
  dfGOF1$LnL[6]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF1$Params[6]<-extractAIC(glmTemp)[1]
  dfGOF1$AIC[6]<-extractAIC(glmTemp)[2]
  dfGOF1$BIC[6]<--2*dfGOF1$LnL[6] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  #-----end rcs -----
  
  #----Royston-Parmar models (RP)----
  temp_GOF<- data.frame(matrix(nrow=18, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  MyScale <- list("hazard","odds","normal")
  for (i in 1:3){
    for (j in 0:5){
      fit<-try(glmTemp <- flexsurvspline(Surv(time, event) ~ 1, data = in_ipd1, k = j, scale = MyScale[[i]]))
      if("try-error" %in% class(fit)) {
        temp_GOF[(i-1)*6+1+j,1] <- "error"
        temp_GOF[(i-1)*6+1+j,2] <- "error"
        temp_GOF[(i-1)*6+1+j,3] <- "error"
        temp_GOF[(i-1)*6+1+j,4] <- "error"
        temp_GOF[(i-1)*6+1+j,5] <- "error"
      }
      else{
        flexsurvspline(Surv(time, event) ~ 1, data = in_ipd1, k = j, scale = MyScale[[i]])
        temp_GOF[(i-1)*6+1+j,1] <- paste(MyScale[[i]],j,sep = ",")
        temp_GOF[(i-1)*6+1+j,2] <- glmTemp$loglik
        temp_GOF[(i-1)*6+1+j,3] <- glmTemp$npars
        temp_GOF[(i-1)*6+1+j,4] <- glmTemp$AIC
        temp_GOF[(i-1)*6+1+j,5] <- log(length(in_ipd1$time))*glmTemp$npars-2*glmTemp$loglik
      }
    }
  }
  temp_GOF<-temp_GOF %>% filter(Model != "error")
  temp_GOF<-arrange(temp_GOF,AIC)
  rp_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  rp_bic<-temp_GOF$Model[1]
  
  temp_scale<-strsplit(rp_aic,",")[[1]][1]
  temp_k<-as.numeric(strsplit(rp_aic,",")[[1]][2])
  glm_rp1<-flexsurvspline(Surv(time, event) ~ 1, data = in_ipd1, k = temp_k, scale = temp_scale)
  dfHazEst11[4,] <- summary(glm_rp1, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[rp_aic]] <- summary(glm_rp1, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[7]<-rp_aic
  dfGOF1$LnL[7]<-sum(ltHaz_1$Events*log(ltHaz_1[[rp_aic]]) - ltHaz_1[[rp_aic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[7]<-glm_rp1$npars
  dfGOF1$AIC[7]<-2*glm_rp1$npars-2*dfGOF1$LnL[7]
  dfGOF1$BIC[7]<-log(length(in_ipd1$time))*glm_rp1$npars-2*dfGOF1$LnL[7]
  
  temp_scale<-strsplit(rp_bic,",")[[1]][1]
  temp_k<-as.numeric(strsplit(rp_bic,",")[[1]][2])
  glm_rp2<-flexsurvspline(Surv(time, event) ~ 1, data = in_ipd1, k = temp_k, scale = temp_scale)
  dfHazEst12[4,] <- summary(glm_rp2, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[rp_bic]] <- summary(glm_rp2, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[8]<-rp_bic
  dfGOF1$LnL[8]<-sum(ltHaz_1$Events*log(ltHaz_1[[rp_bic]]) - ltHaz_1[[rp_bic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[8]<-glm_rp2$npars
  dfGOF1$AIC[8]<-2*glm_rp2$npars-2*dfGOF1$LnL[8]
  dfGOF1$BIC[8]<-log(length(in_ipd1$time))*glm_rp2$npars-2*dfGOF1$LnL[8]
  #-----end rp -----
  
  #----Generalised additive models (GAM)----
  bc2 <- subset(in_ipd1, event==1)
  temp_GOF<- data.frame(matrix(nrow=10, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  for (i in 1:5){
    fit<-try(glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time,
                                                                                                                                               seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_1))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- 'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i] <- 'error'
      temp_GOF$BIC[i] <- 'error'
    }
    else{
      gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time,
                                                                                                                             seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_1)
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
      temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
      temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  gam_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  gam_bic<-temp_GOF$Model[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=gam_aic+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, 
                                                                                                                                          seq(from=0, to=1, by=1/(1+gam_aic))), length=gam_aic+2),family=binomial(link=cloglog), data=ltHaz_1)
  bestgam<-paste("gam_aic",gam_aic,sep = ":")
  dfHazEst11[5,] <- predict(glmTemp, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestgam] <- predict(glmTemp, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[9]<-bestgam
  dfGOF1$LnL[9]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF1$Params[9]<-extractAIC(glmTemp)[1]
  dfGOF1$AIC[9]<-extractAIC(glmTemp)[2]
  dfGOF1$BIC[9]<--2*dfGOF1$LnL[9] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=gam_bic+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, 
                                                                                                                                          seq(from=0, to=1, by=1/(1+gam_bic))), length=gam_bic+2),family=binomial(link=cloglog), data=ltHaz_1)
  bestgam<-paste("gam_bic",gam_bic,sep = ":")
  dfHazEst12[5,] <- predict(glmTemp, newdata=Newtime_1, type="response") # Extrapolated
  ltHaz_1[bestgam] <- predict(glmTemp, newdata=ltHaz_1, type="response")  # Within-sample
  dfGOF1$Model[10]<-bestgam
  dfGOF1$LnL[10]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF1$Params[10]<-extractAIC(glmTemp)[1]
  dfGOF1$AIC[10]<-extractAIC(glmTemp)[2]
  dfGOF1$BIC[10]<--2*dfGOF1$LnL[10] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
  #-----end gam -----
  
  #----Parametric mixture models (PMM)----
  ###----model code----####
  flexsurvmixture <- function(formula, dists, control, sr.control){
    
    n.groups <- length(dists)
    
    dist_list                 <- list()
    dist_list$name            <- paste0(n.groups,"groups: ",paste0(dists,"_", collapse=""))
    dist_list$pars            <- do.call(compile.pars, args = list(dists, n.groups))
    dist_list$location        <- "p1"
    dist_list$transforms      <- do.call(compile.transform, args = list(dists, n.groups))
    dist_list$inv.transforms  <- do.call(compile.inv.transform, args = list(dists, n.groups))
    dist_list$inits           <- function(t){do.call(compile.inits, args = list(t = t, dists, n.groups))}
    
    pfun <- dfun <- list()
    for(i in 1:n.groups){
      pfun[[i]] = get(paste0("p", dists[i]))
      dfun[[i]] = get(paste0("d", dists[i]))
    }
    
    dfns_list = list(
      p = function(q, ...) pmixsurv(pfun, q, n.groups, ...),
      d = function(x, ...) dmixsurv(dfun, x, n.groups, ...),
      mean = function(...) mean_mixsurv(pfun, ...),
      rmst = function(t, ...) rmst_mixsurv(pfun, t, ...) 
    )
    
    optim <- list()
    out <- do.call(
      "flexsurvreg",
      append(
        list(
          formula,
          dist = dist_list,
          dfns = dfns_list,
          control = control, 
          sr.control = sr.control
        ),
        optim
      )
    )
    return(out)
  }
  pmixsurv = function(pfun, q, n.groups, ...) {
    dots <- list(...)
    args <- dots
    args$lower.tail <- F
    args$log.p <- F
    
    theta                     <- c()
    theta[1]                  <- args$p1; args$p1 <- NULL
    if(n.groups > 2){theta[2] <- args$p2; args$p2 <- NULL}
    if(n.groups > 3){theta[3] <- args$p3; args$p3 <- NULL}
    theta <- multinom.p(theta)
    
    args1 <- args; args1[(endsWith(names(args1), "2") | endsWith(names(args1), "3") | endsWith(names(args1), "4"))] <- NULL
    args2 <- args; args2[(endsWith(names(args2), "1") | endsWith(names(args2), "3") | endsWith(names(args2), "4"))] <- NULL
    args3 <- args; args3[(endsWith(names(args3), "1") | endsWith(names(args3), "2") | endsWith(names(args3), "4"))] <- NULL
    args4 <- args; args4[(endsWith(names(args4), "1") | endsWith(names(args4), "2") | endsWith(names(args4), "3"))] <- NULL
    
    names(args1)[endsWith(names(args1), "1")] <- gsub("1", "", names(args1)[endsWith(names(args1), "1")])
    names(args2)[endsWith(names(args2), "2")] <- gsub("2", "", names(args2)[endsWith(names(args2), "2")])
    names(args3)[endsWith(names(args3), "3")] <- gsub("3", "", names(args3)[endsWith(names(args3), "3")])
    names(args4)[endsWith(names(args4), "4")] <- gsub("4", "", names(args4)[endsWith(names(args4), "4")])
    
    out <- (theta[1] * do.call(pfun[[1]], append(list(q), args1)) + 
              theta[2] * do.call(pfun[[2]], append(list(q), args2)) +
              ifelse(n.groups>2, theta[3] * do.call(pfun[[3]], append(list(q), args3)),0) +
              ifelse(n.groups>3, theta[4] * do.call(pfun[[4]], append(list(q), args4)),0)) 
    
    if (is.null(dots$lower.tail) || dots$lower.tail) {
      out <- 1 - out
    }
    if (!is.null(dots$log.p) && dots$log.p) {
      out <- log(out)
    }
    return(out)
  }
  dmixsurv = function(dfun, x, n.groups, ...) {
    dots <- list(...)
    args <- dots
    args$log <- F
    
    theta                     <- c()
    theta[1]                  <- args$p1; args$p1 <- NULL
    if(n.groups > 2){theta[2] <- args$p2; args$p2 <- NULL}
    if(n.groups > 3){theta[3] <- args$p3; args$p3 <- NULL}
    theta <- multinom.p(theta)
    
    args1 <- args; args1[(endsWith(names(args1), "2") | endsWith(names(args1), "3") | endsWith(names(args1), "4"))] <- NULL
    args2 <- args; args2[(endsWith(names(args2), "1") | endsWith(names(args2), "3") | endsWith(names(args2), "4"))] <- NULL
    args3 <- args; args3[(endsWith(names(args3), "1") | endsWith(names(args3), "2") | endsWith(names(args3), "4"))] <- NULL
    args4 <- args; args4[(endsWith(names(args4), "1") | endsWith(names(args4), "2") | endsWith(names(args4), "3"))] <- NULL
    
    names(args1)[endsWith(names(args1), "1")] <- gsub("1", "", names(args1)[endsWith(names(args1), "1")])
    names(args2)[endsWith(names(args2), "2")] <- gsub("2", "", names(args2)[endsWith(names(args2), "2")])
    names(args3)[endsWith(names(args3), "3")] <- gsub("3", "", names(args3)[endsWith(names(args3), "3")])
    names(args4)[endsWith(names(args4), "4")] <- gsub("4", "", names(args4)[endsWith(names(args4), "4")])
    
    
    out <- (theta[1] * do.call(dfun[[1]], append(list(x), args1)) + 
              theta[2] * do.call(dfun[[2]], append(list(x), args2)) +  
              ifelse(n.groups>2,theta[3] * do.call(dfun[[3]], append(list(x), args3)),0) + 
              ifelse(n.groups>3,theta[4] * do.call(dfun[[4]], append(list(x), args4)),0))
    
    if (!is.null(dots$log) && dots$log) {
      out <- log(out)
    }
    return(out)
  }
  rmst_mixsurv = function(pfun, q, t, ...) {
    args <- list(...)
    out <- do.call(
      rmst_generic,
      append(
        list(
          function(q, ...) pmixsurv(pfun, q, ...),
          t = t
        ),
        args
      )
    )
    return(out)
  }
  mean_mixsurv = function(pfun, ...) {
    if(p1 > 0) {
      out <- Inf
    }else {
      args <- list(...)
      out <- do.call(
        rmst_generic,
        append(
          list(
            pfun[[1]],
            t = Inf,
            start = 0
          ),
          args
        )
      )
    }
    return(out)
  }
  multinom.p = function(theta) {
    p <- rep(NA,length(theta) + 1)
    for (i in 1:length(theta)) {
      p[i] <- exp(theta[i])/(1 + sum(exp(theta)))
    }
    p[length(theta) + 1] <- 1 / (1 + sum(exp(theta)))
    return(p)
  }
  #exp
  exp.inits <- function(t) {
    1/mean(t)
  }
  exp.pars <- flexsurv.dists$exp$pars
  exp.t <- list(log)
  exp.it <- list(exp)
  #gamma
  gamma.inits <- function(t) {
    m = mean(t)
    v = var(t)
    c(m^2/v, m/v)
  }
  gamma.pars <- flexsurv.dists$gamma$pars
  gamma.t <- list(log, log)
  gamma.it <- list(exp, exp)
  #gen gamma
  gengamma.inits <- function(t) {
    lt <- log(t[t > 0])
    c(mean(lt), sd(lt), 0)
  }
  gengamma.pars <- flexsurv.dists$gengamma$pars
  gengamma.t <- list(identity, log, identity)
  gengamma.it <- list(identity, exp, identity)
  #gompertz
  gompertz.inits <- function(t) {
    c(0.001, 1/mean(t))
  }
  gompertz.pars <- flexsurv.dists$gompertz$pars
  gompertz.t <- list(identity, log)
  gompertz.it <- list(identity, exp)
  #llogis
  llogis.inits <- function(t) {
    scale <- median(t)
    shape <- 1/log(quantile(t, 0.25)/scale, base = 3)
    if (shape < 0) 
      shape <- 1
    c(shape, scale)
  }
  llogis.pars <- flexsurv.dists$llogis$pars
  llogis.t <- list(log, log)
  llogis.it <- list(exp, exp)
  #lnorm
  lnorm.inits <- function(t) {
    lt <- log(t[t > 0])
    c(mean(lt), sd(lt))
  }
  lnorm.pars <- flexsurv.dists$lnorm$pars
  lnorm.t <- list(identity, log)
  lnorm.it <- list(identity, exp)
  #weibull
  weibull.inits <- function(t) {
    lt <- log(t[t > 0])
    c(1.2/var(lt), exp(mean(lt) + 0.572)) 
  }
  weibull.pars <- flexsurv.dists$weibull$pars
  weibull.t <- list(log, log)
  weibull.it <- list(exp, exp)
  #compile
  compile.inits <- function(t, group.dists, n.groups) {
    if(n.groups >= 1) inits.groups <- rep(0.01, n.groups - 1)
    inits1 <- get(paste0(group.dists[1], ".inits"))
    if(!is.na(group.dists[2])) inits2 <- get(paste0(group.dists[2], ".inits")) else inits2 <- function(t){NULL}
    if(!is.na(group.dists[3])) inits3 <- get(paste0(group.dists[3], ".inits")) else inits3 <- function(t){NULL}
    if(!is.na(group.dists[4])) inits4 <- get(paste0(group.dists[4], ".inits")) else inits4 <- function(t){NULL}
    if(n.groups == 1) return(inits1(t))
    if(n.groups == 2) return(c(inits.groups, inits1(t), inits2(t)))
    if(n.groups == 3) return(c(inits.groups, inits1(t), inits2(t), inits3(t)))
    if(n.groups == 4) return(c(inits.groups, inits1(t), inits2(t), inits3(t), inits4(t)))
    if(n.groups >= 4) stop("too many distributions specified, use 4 groups or fewer")
  }
  compile.pars <- function(group.dists, n.groups) {
    pars.groups = NULL
    if(n.groups >= 1) for(i in 1:(n.groups-1)) pars.groups[i] <- paste0("p",i)
    pars1 <- paste0(get(paste0(group.dists[1], ".pars")),1)
    if(!is.na(group.dists[2])) pars2 <- paste0(get(paste0(group.dists[2], ".pars")),2) else pars2 <- function(){NULL}
    if(!is.na(group.dists[3])) pars3 <- paste0(get(paste0(group.dists[3], ".pars")),3) else pars3 <- function(){NULL}
    if(!is.na(group.dists[4])) pars4 <- paste0(get(paste0(group.dists[4], ".pars")),4) else pars4 <- function(){NULL}
    if(n.groups == 1) return(pars1)
    if(n.groups == 2) return(c(pars.groups, pars1, pars2))
    if(n.groups == 3) return(c(pars.groups, pars1, pars2, pars3))
    if(n.groups == 4) return(c(pars.groups, pars1, pars2, pars3, pars4))
    if(n.groups >= 4) stop("too many distributions specified, use 4 groups or fewer")
  }
  compile.transform <- function(group.dists, n.groups) {
    tgroups = list(NULL)
    if(n.groups >= 1) for(i in 1:(n.groups-1)) tgroups[[i]] <- identity
    t1 <- get(paste0(group.dists[1], ".t"))
    if(!is.na(group.dists[2])) t2 <- get(paste0(group.dists[2], ".t")) else t2 <- function(){NULL}
    if(!is.na(group.dists[3])) t3 <- get(paste0(group.dists[3], ".t")) else t3 <- function(){NULL}
    if(!is.na(group.dists[4])) t4 <- get(paste0(group.dists[4], ".t")) else t4 <- function(){NULL}
    if(n.groups == 1) return(t1)
    if(n.groups == 2) return(c(tgroups, t1, t2))
    if(n.groups == 3) return(c(tgroups, t1, t2, t3))
    if(n.groups == 4) return(c(tgroups, t1, t2, t3, t4))
    if(n.groups >= 4) stop("too many distributions specified, use 4 groups or fewer")
  }
  compile.inv.transform <- function(group.dists, n.groups) {
    itgroups = list(NULL)
    if(n.groups >= 1) for(i in 1:(n.groups-1)) itgroups[[i]] <- identity
    it1 <- get(paste0(group.dists[1], ".it"))
    if(!is.na(group.dists[2])) it2 <- get(paste0(group.dists[2], ".it")) else it2 <- function(){NULL}
    if(!is.na(group.dists[3])) it3 <- get(paste0(group.dists[3], ".it")) else it3 <- function(){NULL}
    if(!is.na(group.dists[4])) it4 <- get(paste0(group.dists[4], ".it")) else it4 <- function(){NULL}
    if(n.groups == 1) return(it1)
    if(n.groups == 2) return(c(itgroups, it1, it2))
    if(n.groups == 3) return(c(itgroups, it1, it2, it3))
    if(n.groups == 4) return(c(itgroups, it1, it2, it3, it4))
    if(n.groups >= 4) stop("too many distributions specified, use 4 groups or fewer")
  }
  ###--------####
  ## simulate results ##
  #newdata
  survgef <- Surv(time=in_ipd1$time , event=in_ipd1$event==1)
  # plot(survgef)
  
  # dists list #
  dists_MM<-as.data.frame(matrix(nrow=28,ncol=2))
  colnames(dists_MM)<-c("dist1",'dist2')
  distlist<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
  
  index1<-1
  for (i in 1:7) {
    for (j in i:7){
      dists_MM$dist1[index1]<- distlist[i]
      dists_MM$dist2[index1]<- distlist[j]
      index1<-index1+1
    }
  }
  
  temp_GOF<- data.frame(matrix(nrow=28, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  
  for (i in 1:28) {
    fit<-try(mixfit <- flexsurvmixture(survgef ~ 1, dists = c(dists_MM$dist1[i],dists_MM$dist2[i]), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8)))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i]<-paste(dists_MM$dist1[i],dists_MM$dist2[i],sep = ",")
      temp_GOF$LnL[i]<-'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i]<-'error'
      temp_GOF$BIC[i]<-'error'
    }
    else{
      flexsurvmixture(survgef ~ 1, dists = c(dists_MM$dist1[i],dists_MM$dist2[i]), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))      
      temp_GOF$Model[i]<-paste(dists_MM$dist1[i],dists_MM$dist2[i],sep = ",")
      temp_GOF$LnL[i]<-mixfit$loglik
      temp_GOF$Params[i]<-mixfit$npars
      temp_GOF$AIC[i]<-(-2*mixfit$loglik+2*mixfit$npars)
      temp_GOF$BIC[i]<-log(length(survgef))*mixfit$npars-2*mixfit$loglik
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  pmm_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  pmm_bic<-temp_GOF$Model[1]
  
  temp_dist1<-strsplit(pmm_aic,",")[[1]][1]
  temp_dist2<-strsplit(pmm_aic,",")[[1]][2]
  glm_pmm1<-flexsurvmixture(survgef ~ 1, dists = c(temp_dist1,temp_dist2), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))
  dfHazEst11[6,] <- summary(glm_pmm1, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[pmm_aic]] <- summary(glm_pmm1, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[11]<-pmm_aic
  dfGOF1$LnL[11]<-sum(ltHaz_1$Events*log(ltHaz_1[[pmm_aic]]) - ltHaz_1[[pmm_aic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[11]<-glm_pmm1$npars
  dfGOF1$AIC[11]<-2*glm_pmm1$npars-2*dfGOF1$LnL[11]
  dfGOF1$BIC[11]<-log(length(in_ipd1$time))*glm_pmm1$npars-2*dfGOF1$LnL[11]
  
  temp_dist1<-strsplit(pmm_bic,",")[[1]][1]
  temp_dist2<-strsplit(pmm_bic,",")[[1]][2]
  glm_pmm2<-flexsurvmixture(survgef ~ 1, dists = c(temp_dist1,temp_dist2), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))
  dfHazEst12[6,] <- summary(glm_pmm2, t=Newtime_1$Time, type="hazard")[[1]]$est/12
  ltHaz_1[[pmm_bic]] <- summary(glm_pmm2, t=ltHaz_1$Time, type="hazard")[[1]]$est/12
  dfGOF1$Model[12]<-pmm_bic
  dfGOF1$LnL[12]<-sum(ltHaz_1$Events*log(ltHaz_1[[pmm_bic]]) - ltHaz_1[[pmm_bic]]*ltHaz_1$AtRisk) + llCons1
  dfGOF1$Params[12]<-glm_pmm2$npars
  dfGOF1$AIC[12]<-2*glm_pmm2$npars-2*dfGOF1$LnL[12]
  dfGOF1$BIC[12]<-log(length(in_ipd1$time))*glm_pmm2$npars-2*dfGOF1$LnL[12]
  #-----end pmm -----
  
  ###-----------###
  #Extrapolation-2
  ###-----------###
  #----Standard survival model (SSM)----
  temp_GOF<- data.frame(matrix(nrow=7, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  MyDists <- list("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
  for (i in 1:7){
    glmTemp <- flexsurvreg(Surv(time, event) ~ 1, data = in_ipd2, dist = MyDists[[i]])
    temp_GOF[i,1] <- MyDists[[i]]
    temp_GOF[i,2] <- glmTemp$loglik
    temp_GOF[i,3] <- glmTemp$npars
    temp_GOF[i,4] <- glmTemp$AIC
    temp_GOF[i,5] <- log(length(in_ipd2$time))*glmTemp$npars-2*glmTemp$loglik
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  ssm_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  ssm_bic<-temp_GOF$Model[1]
  
  glm_ssm1<-flexsurvreg(Surv(time, event) ~ 1, data = in_ipd2, dist = ssm_aic)
  dfHazEst21[1,] <- summary(glm_ssm1, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[ssm_aic]] <- summary(glm_ssm1, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[1]<-ssm_aic
  dfGOF2$LnL[1]<-sum(ltHaz_2$Events*log(ltHaz_2[[ssm_aic]]) - ltHaz_2[[ssm_aic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[1]<-glm_ssm1$npars
  dfGOF2$AIC[1]<-2*glm_ssm1$npars-2*dfGOF2$LnL[1]
  dfGOF2$BIC[1]<-log(length(in_ipd2$time))*glm_ssm1$npars-2*dfGOF2$LnL[1]
  
  glm_ssm2<-flexsurvreg(Surv(time, event) ~ 1, data = in_ipd2, dist = ssm_bic)
  dfHazEst22[1,] <- summary(glm_ssm2, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[ssm_bic]] <- summary(glm_ssm2, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[2]<-ssm_bic
  dfGOF2$LnL[2]<-sum(ltHaz_2$Events*log(ltHaz_2[[ssm_bic]]) - ltHaz_2[[ssm_bic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[2]<-glm_ssm2$npars
  dfGOF2$AIC[2]<-2*glm_ssm2$npars-2*dfGOF2$LnL[2]
  dfGOF2$BIC[2]<-log(length(in_ipd2$time))*glm_ssm2$npars-2*dfGOF2$LnL[2]
  #-----end ssm -----
  
  #----Fractional polynomial (FP)----
  temp_GOF<- data.frame(matrix(nrow=44, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  #-----FP1 -----
  MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
  for (i in 1:7){
    glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    temp_GOF$Model[i] <- paste(MyPowers[[1]][i],"NA",sep = ",") 
    temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  }
  ### run for 0
  glmTemp <- glm (cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  temp_GOF$Model[8] <- paste(0,"NA",sep = ",")
  temp_GOF$LnL[8] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  temp_GOF$Params[8]<-extractAIC(glmTemp)[1]
  temp_GOF$AIC[8] <- extractAIC(glmTemp)[2]
  temp_GOF$BIC[8] <- -2*temp_GOF$LnL[8] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  
  #-----FP2 -----
  myLnL <- array(dim=36)
  myAIC <- array(dim=36)
  MyPowers <- list(c(-2,-1,-0.5,0.5,1,2,3))
  index <- 9
  for (i in 1:7){
    for (j in 1:7){
      if (j > i) {
        glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][j])+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2) 
        temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
        temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
        temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
        temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
        temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
        index <- index + 1
      }
    }
  }
  for (i in 1:7) {
    glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^MyPowers[[1]][i]*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)# 
    temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
    temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
    index <- index + 1
    
  }
  
  for (i in 1:7) {
    glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)# 
    temp_GOF$Model[index] <- paste(MyPowers[[1]][i],MyPowers[[1]][j],sep = ",")
    temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
    temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
    temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
    temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
    index <- index + 1
  }
  
  glmTemp <- glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)# 
  temp_GOF$Model[index] <- paste(0,0,sep = ",")
  temp_GOF$LnL[index] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  temp_GOF$Params[index]<-extractAIC(glmTemp)[1]
  temp_GOF$AIC[index] <- extractAIC(glmTemp)[2]
  temp_GOF$BIC[index] <- -2*temp_GOF$LnL[index] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  
  temp_GOF<-arrange(temp_GOF,AIC)
  fp_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  fp_bic<-temp_GOF$Model[1]
  
  FP_pow_1<-as.numeric((strsplit(fp_aic,","))[[1]][1])
  FP_pow_2<-as.numeric((strsplit(fp_aic,","))[[1]][2])
  if(is.na(FP_pow_2)){
    if(FP_pow_1 == 0){
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    } else {
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    }
  } else if (FP_pow_1 == FP_pow_2 & FP_pow_2 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else if (FP_pow_1==0 & FP_pow_2== 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else if (FP_pow_2==0 & FP_pow_1 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else {
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  }
  bestfp<-paste("FP_aic",fp_aic,sep = ":")
  dfHazEst21[2,] <- predict(modFP, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestfp] <- predict(modFP, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[3]<-bestfp
  dfGOF2$LnL[3]<-(extractAIC(modFP)[2] - 2*extractAIC(modFP)[1])*(-0.5)
  dfGOF2$Params[3]<-extractAIC(modFP)[1]
  dfGOF2$AIC[3]<-extractAIC(modFP)[2]
  dfGOF2$BIC[3]<--2*dfGOF2$LnL[3] + log(length(ltHaz_2$Time))*extractAIC(modFP)[1]
  
  FP_pow_1<-as.numeric((strsplit(fp_bic,","))[[1]][1])
  FP_pow_2<-as.numeric((strsplit(fp_bic,","))[[1]][2])
  if(is.na(FP_pow_2)){
    if(FP_pow_1 == 0){
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    } else {
      modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
    }
  } else if (FP_pow_1 == FP_pow_2 & FP_pow_2 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2*lnTime) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else if (FP_pow_1==0 & FP_pow_2== 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^0*lnTime) + I(Time^0*lnTime*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else if (FP_pow_2==0 & FP_pow_1 != 0){
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^MyPowers[[1]][i]) + I(Time^0*lnTime)+ offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  } else {
    modFP<-glm(cbind(Events,AtRisk-Events) ~ I(Time^FP_pow_1) + I(Time^FP_pow_2) + offset(log(timedelta)), family=binomial(link=cloglog), data=ltHaz_2)
  }
  bestfp<-paste("FP_bic",fp_bic,sep = ":")
  dfHazEst22[2,] <- predict(modFP, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestfp] <- predict(modFP, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[4]<-bestfp
  dfGOF2$LnL[4]<-(extractAIC(modFP)[2] - 2*extractAIC(modFP)[1])*(-0.5)
  dfGOF2$Params[4]<-extractAIC(modFP)[1]
  dfGOF2$AIC[4]<-extractAIC(modFP)[2]
  dfGOF2$BIC[4]<--2*dfGOF2$LnL[4] + log(length(ltHaz_2$Time))*extractAIC(modFP)[1]
  #-----end fp -----
  
  #----Restricted cubic splines (RCS)----
  bc2 <- subset(in_ipd2, event==1)
  temp_GOF<- data.frame(matrix(nrow=5, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  for (i in 1:5){
    fit<-try(glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_2))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- 'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i] <- 'error'
      temp_GOF$BIC[i] <- 'error'
    }
    else{
      gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_2)      
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
      temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
      temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  rcs_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  rcs_bic<-temp_GOF$Model[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=rcs_aic+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+rcs_aic))), length=rcs_aic+2),family=binomial(link=cloglog), data=ltHaz_2)
  bestrcs<-paste("RCS_aic",rcs_aic,sep = ":")
  dfHazEst21[3,] <- predict(glmTemp, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestrcs] <- predict(glmTemp, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[5]<-bestrcs
  dfGOF2$LnL[5]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF2$Params[5]<-extractAIC(glmTemp)[1]
  dfGOF2$AIC[5]<-extractAIC(glmTemp)[2]
  dfGOF2$BIC[5]<--2*dfGOF2$LnL[5] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=rcs_bic+2, fx=TRUE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, seq(from=0, to=1, by=1/(1+rcs_bic))), length=rcs_bic+2),family=binomial(link=cloglog), data=ltHaz_2)
  bestrcs<-paste("RCS_bic",rcs_bic,sep = ":")
  dfHazEst22[3,] <- predict(glmTemp, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestrcs] <- predict(glmTemp, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[6]<-bestrcs
  dfGOF2$LnL[6]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF2$Params[6]<-extractAIC(glmTemp)[1]
  dfGOF2$AIC[6]<-extractAIC(glmTemp)[2]
  dfGOF2$BIC[6]<--2*dfGOF2$LnL[6] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  #-----end rcs -----
  
  #----Royston-Parmar models (RP)----
  temp_GOF<- data.frame(matrix(nrow=18, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  MyScale <- list("hazard","odds","normal")
  for (i in 1:3){
    for (j in 0:5){
      fit<-try(glmTemp <- flexsurvspline(Surv(time, event) ~ 1, data = in_ipd2, k = j, scale = MyScale[[i]]))
      if("try-error" %in% class(fit)) {
        temp_GOF[(i-1)*6+1+j,1] <- "error"
        temp_GOF[(i-1)*6+1+j,2] <- "error"
        temp_GOF[(i-1)*6+1+j,3] <- "error"
        temp_GOF[(i-1)*6+1+j,4] <- "error"
        temp_GOF[(i-1)*6+1+j,5] <- "error"
      }
      else{
        flexsurvspline(Surv(time, event) ~ 1, data = in_ipd2, k = j, scale = MyScale[[i]])
        temp_GOF[(i-1)*6+1+j,1] <- paste(MyScale[[i]],j,sep = ",")
        temp_GOF[(i-1)*6+1+j,2] <- glmTemp$loglik
        temp_GOF[(i-1)*6+1+j,3] <- glmTemp$npars
        temp_GOF[(i-1)*6+1+j,4] <- glmTemp$AIC
        temp_GOF[(i-1)*6+1+j,5] <- log(length(in_ipd2$time))*glmTemp$npars-2*glmTemp$loglik
      }
    }
  }
  temp_GOF<-temp_GOF %>% filter(Model != "error")
  temp_GOF<-arrange(temp_GOF,AIC)
  rp_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  rp_bic<-temp_GOF$Model[1]
  
  temp_scale<-strsplit(rp_aic,",")[[1]][1]
  temp_k<-as.numeric(strsplit(rp_aic,",")[[1]][2])
  glm_rp1<-flexsurvspline(Surv(time, event) ~ 1, data = in_ipd2, k = temp_k, scale = temp_scale)
  dfHazEst21[4,] <- summary(glm_rp1, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[rp_aic]] <- summary(glm_rp1, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[7]<-rp_aic
  dfGOF2$LnL[7]<-sum(ltHaz_2$Events*log(ltHaz_2[[rp_aic]]) - ltHaz_2[[rp_aic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[7]<-glm_rp1$npars
  dfGOF2$AIC[7]<-2*glm_rp1$npars-2*dfGOF2$LnL[7]
  dfGOF2$BIC[7]<-log(length(in_ipd2$time))*glm_rp1$npars-2*dfGOF2$LnL[7]
  
  temp_scale<-strsplit(rp_bic,",")[[1]][1]
  temp_k<-as.numeric(strsplit(rp_bic,",")[[1]][2])
  glm_rp2<-flexsurvspline(Surv(time, event) ~ 1, data = in_ipd2, k = temp_k, scale = temp_scale)
  dfHazEst22[4,] <- summary(glm_rp2, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[rp_bic]] <- summary(glm_rp2, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[8]<-rp_bic
  dfGOF2$LnL[8]<-sum(ltHaz_2$Events*log(ltHaz_2[[rp_bic]]) - ltHaz_2[[rp_bic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[8]<-glm_rp2$npars
  dfGOF2$AIC[8]<-2*glm_rp2$npars-2*dfGOF2$LnL[8]
  dfGOF2$BIC[8]<-log(length(in_ipd2$time))*glm_rp2$npars-2*dfGOF2$LnL[8]
  #-----end rp -----
  
  #----Generalised additive models (GAM)----
  bc2 <- subset(in_ipd2, event==1)
  temp_GOF<- data.frame(matrix(nrow=10, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  for (i in 1:5){
    fit<-try(glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time,
                                                                                                                                               seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_2))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- 'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i] <- 'error'
      temp_GOF$BIC[i] <- 'error'
    }
    else{
      gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=i+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time,
                                                                                                                             seq(from=0, to=1, by=1/(1+i))), length=i+2), family=binomial(link=cloglog), data=ltHaz_2)
      temp_GOF$Model[i] <- i
      temp_GOF$LnL[i] <- (extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
      temp_GOF$Params[i]<-extractAIC(glmTemp)[1]
      temp_GOF$AIC[i] <- extractAIC(glmTemp)[2]
      temp_GOF$BIC[i] <- -2*temp_GOF$LnL[i] + log(length(ltHaz_1$Time))*extractAIC(glmTemp)[1]
    }
  }
  
  temp_GOF<-arrange(temp_GOF,AIC)
  gam_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  gam_bic<-temp_GOF$Model[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=gam_aic+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, 
                                                                                                                                          seq(from=0, to=1, by=1/(1+gam_aic))), length=gam_aic+2),family=binomial(link=cloglog), data=ltHaz_2)
  bestgam<-paste("gam_aic",gam_aic,sep = ":")
  dfHazEst21[5,] <- predict(glmTemp, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestgam] <- predict(glmTemp, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[9]<-bestgam
  dfGOF2$LnL[9]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF2$Params[9]<-extractAIC(glmTemp)[1]
  dfGOF2$AIC[9]<-extractAIC(glmTemp)[2]
  dfGOF2$BIC[9]<--2*dfGOF2$LnL[9] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  
  glmTemp <- gam(cbind(Events,AtRisk-Events) ~ s(Time, bs="cr", k=gam_bic+2, fx=FALSE) + offset(log(timedelta)), knots=list(Time=quantile(bc2$time, 
                                                                                                                                          seq(from=0, to=1, by=1/(1+gam_bic))), length=gam_bic+2),family=binomial(link=cloglog), data=ltHaz_2)
  bestgam<-paste("gam_bic",gam_bic,sep = ":")
  dfHazEst22[5,] <- predict(glmTemp, newdata=Newtime_2, type="response") # Extrapolated
  ltHaz_2[bestgam] <- predict(glmTemp, newdata=ltHaz_2, type="response")  # Within-sample
  dfGOF2$Model[10]<-bestgam
  dfGOF2$LnL[10]<-(extractAIC(glmTemp)[2] - 2*extractAIC(glmTemp)[1])*(-0.5)
  dfGOF2$Params[10]<-extractAIC(glmTemp)[1]
  dfGOF2$AIC[10]<-extractAIC(glmTemp)[2]
  dfGOF2$BIC[10]<--2*dfGOF2$LnL[10] + log(length(ltHaz_2$Time))*extractAIC(glmTemp)[1]
  #-----end gam -----
  
  #----Parametric mixture models (PMM)----
  ###--------####
  ## simulate results ##
  #newdata
  survgef <- Surv(time=in_ipd2$time , event=in_ipd2$event==1)
  # plot(survgef)
  
  # dists list #
  dists_MM<-as.data.frame(matrix(nrow=28,ncol=2))
  colnames(dists_MM)<-c("dist1",'dist2')
  distlist<-c("exp","weibull","gamma","lnorm","gompertz","llogis","gengamma")
  
  index1<-1
  for (i in 1:7) {
    for (j in i:7){
      dists_MM$dist1[index1]<- distlist[i]
      dists_MM$dist2[index1]<- distlist[j]
      index1<-index1+1
    }
  }
  
  temp_GOF<- data.frame(matrix(nrow=28, ncol=5))
  colnames(temp_GOF) <- c("Model","LnL","Params","AIC","BIC")
  
  for (i in 1:28) {
    fit<-try(mixfit <- flexsurvmixture(survgef ~ 1, dists = c(dists_MM$dist1[i],dists_MM$dist2[i]), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8)))
    if("try-error" %in% class(fit)) {
      temp_GOF$Model[i]<-paste(dists_MM$dist1[i],dists_MM$dist2[i],sep = ",")
      temp_GOF$LnL[i]<-'error'
      temp_GOF$Params[i]<-'error'
      temp_GOF$AIC[i]<-'error'
      temp_GOF$BIC[i]<-'error'
    }
    else{
      flexsurvmixture(survgef ~ 1, dists = c(dists_MM$dist1[i],dists_MM$dist2[i]), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))      
      temp_GOF$Model[i]<-paste(dists_MM$dist1[i],dists_MM$dist2[i],sep = ",")
      temp_GOF$LnL[i]<-mixfit$loglik
      temp_GOF$Params[i]<-mixfit$npars
      temp_GOF$AIC[i]<-(-2*mixfit$loglik+2*mixfit$npars)
      temp_GOF$BIC[i]<-log(length(survgef))*mixfit$npars-2*mixfit$loglik
    }
  }
  temp_GOF<-arrange(temp_GOF,AIC)
  pmm_aic<-temp_GOF$Model[1]
  temp_GOF<-arrange(temp_GOF,BIC)
  pmm_bic<-temp_GOF$Model[1]
  
  temp_dist1<-strsplit(pmm_aic,",")[[1]][1]
  temp_dist2<-strsplit(pmm_aic,",")[[1]][2]
  glm_pmm1<-flexsurvmixture(survgef ~ 1, dists = c(temp_dist1,temp_dist2), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))
  dfHazEst21[6,] <- summary(glm_pmm1, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[pmm_aic]] <- summary(glm_pmm1, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[11]<-pmm_aic
  dfGOF2$LnL[11]<-sum(ltHaz_2$Events*log(ltHaz_2[[pmm_aic]]) - ltHaz_2[[pmm_aic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[11]<-glm_pmm1$npars
  dfGOF2$AIC[11]<-2*glm_pmm1$npars-2*dfGOF2$LnL[11]
  dfGOF2$BIC[11]<-log(length(in_ipd2$time))*glm_pmm1$npars-2*dfGOF2$LnL[11]
  
  temp_dist1<-strsplit(pmm_bic,",")[[1]][1]
  temp_dist2<-strsplit(pmm_bic,",")[[1]][2]
  glm_pmm2<-flexsurvmixture(survgef ~ 1, dists = c(temp_dist1,temp_dist2), control = list(reltol = 1e-8), sr.control= list(reltol = 1e-8))
  dfHazEst22[6,] <- summary(glm_pmm2, t=Newtime_2$Time, type="hazard")[[1]]$est/12
  ltHaz_2[[pmm_bic]] <- summary(glm_pmm2, t=ltHaz_2$Time, type="hazard")[[1]]$est/12
  dfGOF2$Model[12]<-pmm_bic
  dfGOF2$LnL[12]<-sum(ltHaz_2$Events*log(ltHaz_2[[pmm_bic]]) - ltHaz_2[[pmm_bic]]*ltHaz_2$AtRisk) + llCons2
  dfGOF2$Params[12]<-glm_pmm2$npars
  dfGOF2$AIC[12]<-2*glm_pmm2$npars-2*dfGOF2$LnL[12]
  dfGOF2$BIC[12]<-log(length(in_ipd2$time))*glm_pmm2$npars-2*dfGOF2$LnL[12]
  #-----end pmm -----
  dfHaz_1aic <-t(dfHazEst11)
  colnames(dfHaz_1aic) <- c("SSM","FP","RCS","RP","GAM","PMM")
  dfHaz_1aic <- cbind(data.frame(time=ltHaz_5$Time), dfHaz_1aic)
  dfSurv_1aic<-as.data.frame(matrix(nrow = length(dfHaz_1aic$time),ncol = 7))
  dfSurv_1aic[[1]]<-dfHaz_1aic$time
  for ( i in 2:7 ) {
    dftemp<-data.frame(dfHaz_1aic$time,dfHaz_1aic[[i]])
    dftemp<-dftemp %>% 
      dplyr::arrange(dftemp[[1]]) %>% 
      dplyr::mutate(cumhaz = cumsum(dftemp[[2]])) %>% 
      dplyr::mutate(survProp = exp(-1*cumhaz))
    dfSurv_1aic[[i]]<-dftemp[[4]]
  }
  colnames(dfSurv_1aic) <- c("Time","SSM","FP","RCS","RP","GAM","PMM")
  
  dfHaz_1bic <-t(dfHazEst12)
  colnames(dfHaz_1bic) <- c("SSM","FP","RCS","RP","GAM","PMM")
  dfHaz_1bic <- cbind(data.frame(time=ltHaz_5$Time), dfHaz_1bic)
  dfSurv_1bic<-as.data.frame(matrix(nrow = length(dfHaz_1bic$time),ncol = 7))
  dfSurv_1bic[[1]]<-dfHaz_1bic$time
  for ( i in 2:7 ) {
    dftemp<-data.frame(dfHaz_1bic$time,dfHaz_1bic[[i]])
    dftemp<-dftemp %>% 
      dplyr::arrange(dftemp[[1]]) %>% 
      dplyr::mutate(cumhaz = cumsum(dftemp[[2]])) %>% 
      dplyr::mutate(survProp = exp(-1*cumhaz))
    dfSurv_1bic[[i]]<-dftemp[[4]]
  }
  colnames(dfSurv_1bic) <- c("Time","SSM","FP","RCS","RP","GAM","PMM")
  
  dfHaz_2aic <-t(dfHazEst21)
  colnames(dfHaz_2aic) <- c("SSM","FP","RCS","RP","GAM","PMM")
  dfHaz_2aic <- cbind(data.frame(time=ltHaz_5$Time), dfHaz_2aic)
  dfSurv_2aic<-as.data.frame(matrix(nrow = length(dfHaz_2aic$time),ncol = 7))
  dfSurv_2aic[[1]]<-dfHaz_2aic$time
  for ( i in 2:7 ) {
    dftemp<-data.frame(dfHaz_2aic$time,dfHaz_2aic[[i]])
    dftemp<-dftemp %>% 
      dplyr::arrange(dftemp[[1]]) %>% 
      dplyr::mutate(cumhaz = cumsum(dftemp[[2]])) %>% 
      dplyr::mutate(survProp = exp(-1*cumhaz))
    dfSurv_2aic[[i]]<-dftemp[[4]]
  }
  colnames(dfSurv_2aic) <- c("Time","SSM","FP","RCS","RP","GAM","PMM")
  
  dfHaz_2bic <-t(dfHazEst22)
  colnames(dfHaz_2bic) <- c("SSM","FP","RCS","RP","GAM","PMM")
  dfHaz_2bic <- cbind(data.frame(time=ltHaz_5$Time), dfHaz_2bic)
  dfSurv_2bic<-as.data.frame(matrix(nrow = length(dfHaz_2bic$time),ncol = 7))
  dfSurv_2bic[[1]]<-dfHaz_2bic$time
  for ( i in 2:7 ) {
    dftemp<-data.frame(dfHaz_2bic$time,dfHaz_2bic[[i]])
    dftemp<-dftemp %>% 
      dplyr::arrange(dftemp[[1]]) %>% 
      dplyr::mutate(cumhaz = cumsum(dftemp[[2]])) %>% 
      dplyr::mutate(survProp = exp(-1*cumhaz))
    dfSurv_2bic[[i]]<-dftemp[[4]]
  }
  colnames(dfSurv_2bic) <- c("Time","SSM","FP","RCS","RP","GAM","PMM")
  
  output=list(dfSurv_1aic=dfSurv_1aic,
              dfSurv_1bic=dfSurv_1bic,
              dfSurv_2aic=dfSurv_2aic,
              dfSurv_2bic=dfSurv_2bic,
              dfGOF1=dfGOF1,
              dfGOF2=dfGOF2,
              ltHaz_1=ltHaz_1,
              ltHaz_2=ltHaz_2,
              ltHaz_5=ltHaz_5
  )
  return(output)
}