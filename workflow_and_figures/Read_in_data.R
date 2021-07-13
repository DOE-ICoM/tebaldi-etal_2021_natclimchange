remove(list=ls())
filedir<-"~/Desktop/ForCasey/"
source(paste0(filedir,"Rcode/allfunctions.r"))

rl1k<-read.csv(filedir%&%"CSV/Kirezci/Kirezci_ESLs.csv")
rl1r<-read.csv(filedir%&%"CSV/Rasmussen/Rasmussen_ESLs.csv")
rl1v<-read.csv(filedir%&%"CSV/Vousdoukas/Vousdoukas_ESLs.csv")



for(whichdecade in seq(2020,2100,by=10)){
proj1k<-array(dim=c(11,8,8086,2))
proj1r<-array(dim=c(11,8,200,2))
proj1v<-array(dim=c(11,8,7472,2))

dimnames(proj1k)<-list(c("ID" ,    "year" ,  "q0.01" , "q0.05",  "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99",  "q0.995", "q0.999"),
                       rasmussenScenarios,NULL,c("rasmussen","bvw"))
dimnames(proj1r)<-list(c("ID" ,    "year" ,  "q0.01" , "q0.05",  "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99",  "q0.995", "q0.999"),
                       rasmussenScenarios,NULL,c("rasmussen","bvw"))
dimnames(proj1v)<-list(c("ID" ,    "year" ,  "q0.01" , "q0.05",  "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99",  "q0.995", "q0.999"),
                       rasmussenScenarios,NULL,c("rasmussen","bvw"))


for(scenario in rasmussenScenarios){
    for(proj in c("rasmussen","bvw")){
        temp<-read.csv(file=filedir%&%"CSV/Kirezci/Kirezci_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv")
        proj1k[,scenario,,proj]<-t(temp)
        temp<-read.csv(file=filedir%&%"CSV/Rasmussen/Rasmussen_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv")
        proj1r[,scenario,,proj]<-t(temp)
        temp<-read.csv(file=filedir%&%"CSV/Vousdoukas/Vousdoukas_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv")
        proj1v[,scenario,,proj]<-t(temp)
    }
}

save(list=c("rl1k","rl1r","rl1v","proj1k","proj1r","proj1v"),file=filedir%&%"Rdatasets/RData_esl_slr_"%&%whichdecade)
}


# now source(filedir%&%"Rcode/probabilisticprojections_and_TWL100TWL1difference_fordatapub.r")
