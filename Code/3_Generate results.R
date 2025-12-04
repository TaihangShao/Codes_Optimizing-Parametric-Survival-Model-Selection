
Arm_res<-function(input_list){
  dfSurv_1aic=input_list$dfSurv_1aic
  dfSurv_1bic=input_list$dfSurv_1bic
  dfSurv_2aic=input_list$dfSurv_2aic
  dfSurv_2bic=input_list$dfSurv_2bic
  dfGOF1=input_list$dfGOF1
  dfGOF2=input_list$dfGOF2
  ltHaz_1=input_list$ltHaz_1
  ltHaz_2=input_list$ltHaz_2
  ltHaz_5=input_list$ltHaz_5
  
  dfSurv_1aic=dfSurv_1aic[-(1:length(ltHaz_1$Time)),]
  dfSurv_1bic=dfSurv_1bic[-(1:length(ltHaz_1$Time)),]
  dfSurv_2aic=dfSurv_2aic[-(1:length(ltHaz_2$Time)),]
  dfSurv_2bic=dfSurv_2bic[-(1:length(ltHaz_2$Time)),]
  ltHaz_51=ltHaz_5[-(1:length(ltHaz_1$Time)),]
  ltHaz_52=ltHaz_5[-(1:length(ltHaz_2$Time)),]
  
  RMST<-as.data.frame(matrix(nrow = 4,ncol = 6))
  colnames(RMST) <- c("SSM","FP","RCS","RP","GAM","PMM")
  rownames(RMST) <- c("1_AIC","1_BIC","2_AIC","2_BIC")


  
  for (i in 1:6) {
    RMST[1,i]<-sum((dfSurv_1aic$Time[2]-dfSurv_1aic$Time[1])*dfSurv_1aic[i+1])
    RMST[2,i]<-sum((dfSurv_1bic$Time[2]-dfSurv_1bic$Time[1])*dfSurv_1bic[i+1])
    RMST[3,i]<-sum((dfSurv_2aic$Time[2]-dfSurv_2aic$Time[1])*dfSurv_2aic[i+1])
    RMST[4,i]<-sum((dfSurv_2bic$Time[2]-dfSurv_2bic$Time[1])*dfSurv_2bic[i+1])
    }
  RMST_KM1<-sum((ltHaz_51$timedelta[1])*ltHaz_51$surv)
  RMST_KM2<-sum((ltHaz_52$timedelta[1])*ltHaz_52$surv)
  
  #Mean RMST
  MRMST<-RMST
  MRMST[1:2,]<-MRMST[1:2,]/length(ltHaz_51$Time)
  MRMST[3:4,]<-MRMST[3:4,]/length(ltHaz_52$Time)
  
  #RMSD 
  RMSD<-RMST
  RMSD[1:2,]<-RMSD[1:2,]-RMST_KM1
  RMSD[3:4,]<-RMSD[3:4,]-RMST_KM2
  
  #MAE (ARMST)
  ARMST<-as.data.frame(matrix(nrow = 4,ncol = 6))
  colnames(ARMST) <- c("SSM","FP","RCS","RP","GAM","PMM")
  rownames(ARMST) <- c("1_AIC","1_BIC","2_AIC","2_BIC")
  
  for (i in 1:6) {
    ARMST[1,i]<-sum(abs((dfSurv_1aic$Time[2]-dfSurv_1aic$Time[1])*dfSurv_1aic[i+1]-(ltHaz_51$timedelta[1])*ltHaz_51$surv))/length(ltHaz_51$Time)
    ARMST[2,i]<-sum(abs((dfSurv_1bic$Time[2]-dfSurv_1bic$Time[1])*dfSurv_1bic[i+1]-(ltHaz_51$timedelta[1])*ltHaz_51$surv))/length(ltHaz_51$Time)
    ARMST[3,i]<-sum(abs((dfSurv_2aic$Time[2]-dfSurv_2aic$Time[1])*dfSurv_2aic[i+1]-(ltHaz_52$timedelta[1])*ltHaz_52$surv))/length(ltHaz_52$Time)
    ARMST[4,i]<-sum(abs((dfSurv_2bic$Time[2]-dfSurv_2bic$Time[1])*dfSurv_2bic[i+1]-(ltHaz_52$timedelta[1])*ltHaz_52$surv))/length(ltHaz_52$Time)
    }
  
  #MAPE
  MAPE<-as.data.frame(matrix(nrow = 4,ncol = 6))
  colnames(MAPE) <- c("SSM","FP","RCS","RP","GAM","PMM")
  rownames(MAPE) <- c("1_AIC","1_BIC","2_AIC","2_BIC")
  
  for (i in 1:6) {
    MAPE[1,i]<-sum(abs((dfSurv_1aic$Time[2]-dfSurv_1aic$Time[1])*dfSurv_1aic[i+1]-(ltHaz_51$timedelta[1])*ltHaz_51$surv)/(ltHaz_51$timedelta[1])*ltHaz_51$surv)/length(ltHaz_51$Time)
    MAPE[2,i]<-sum(abs((dfSurv_1bic$Time[2]-dfSurv_1bic$Time[1])*dfSurv_1bic[i+1]-(ltHaz_51$timedelta[1])*ltHaz_51$surv)/(ltHaz_51$timedelta[1])*ltHaz_51$surv)/length(ltHaz_51$Time)
    MAPE[3,i]<-sum(abs((dfSurv_2aic$Time[2]-dfSurv_2aic$Time[1])*dfSurv_2aic[i+1]-(ltHaz_52$timedelta[1])*ltHaz_52$surv)/(ltHaz_52$timedelta[1])*ltHaz_52$surv)/length(ltHaz_52$Time)
    MAPE[4,i]<-sum(abs((dfSurv_2bic$Time[2]-dfSurv_2bic$Time[1])*dfSurv_2bic[i+1]-(ltHaz_52$timedelta[1])*ltHaz_52$surv)/(ltHaz_52$timedelta[1])*ltHaz_52$surv)/length(ltHaz_52$Time)
  }
  output<-list(RMST=RMST,
               RMST_KM1=RMST_KM1,
               RMST_KM2=RMST_KM2,
               MRMST=MRMST,
               RMSD=RMSD,
               ARMST=ARMST,
               MAPE=MAPE)
  return(output)
  }