remove(list=ls())

filedir<-"~/Desktop/ForCasey/"
source(paste0(filedir,"Rcode/allfunctions.r"))

for(whichyear in c(2020,2030,2040,2050,2060,2070,2080,2090,2100)){
    load("Rdatasets/RData_rl_proj_matched_"%&%whichyear)
    print(dimnames(rl1r)[[2]])
    rl1r<-rl1r[,-c(2,3,18)]
    rl1k<-rl1e
    proj1k<-proj1e
    dimnames(proj1k)[[3]]<-dimnames(rl1k)[[1]]
    dimnames(proj1v)[[3]]<-dimnames(rl1v)[[1]]
    dimnames(proj1r)[[3]]<-dimnames(rl1r)[[1]]
    dimnames(proj1k)[[4]][[1]]<-"rasmussen"
    dimnames(proj1v)[[4]][[1]]<-"rasmussen"
    dimnames(proj1r)[[4]][[1]]<-"rasmussen"

    save(list=c("rl1k","rl1v","rl1r","proj1k","proj1v","proj1r"),file=filedir%&%"Rdatasets/RData_esl_slr_"%&%whichyear)
}

for(whichdecade in c(2020,2030,2040,2050,2060,2070,2080,2090,2100)){
#write these out as CSV files
load(filedir%&%"Rdatasets/RData_esl_slr_"%&%whichdecade)
for(scenario in rasmussenScenarios){
    for(proj in c("rasmussen","bvw")){
        temp<-t(proj1k[,scenario,,proj])
        write.csv(temp,file="~/Desktop/ForCasey/CSV/Kirezci_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv",row.names=F)
        temp<-t(proj1r[,scenario,,proj])
        write.csv(temp,file="~/Desktop/ForCasey/CSV/Rasmussen_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv",row.names=F)
        temp<-t(proj1v[,scenario,,proj])
        write.csv(temp,file="~/Desktop/ForCasey/CSV/Vousdoukas_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv",row.names=F)
    }
}
}
