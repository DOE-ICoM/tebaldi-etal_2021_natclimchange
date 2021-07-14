remove(list=ls())
filedir<-"./"
source(paste0(filedir,"Rcode/allfunctions.r"))

whichdecade<-2100
load(filedir%&%"Rdatasets/RData_esl_slr_2100")
load(filedir%&%"Rdatasets/RData_convolutions_fisherinfo_"%&%whichdecade)


eliminate.rasmussen<-numeric(0)
for(i in 1:200)if(all(is.na(rasmussen.shape.scale.TWL.samples[i,,])))eliminate.rasmussen<-c(eliminate.rasmussen,i)
if(length(eliminate.rasmussen)>0)rasmussensubsetID<-seq(length(dimnames(rasmussen.shape.scale.TWL.samples)[[1]]))[-eliminate.rasmussen]

eliminate.kirezci<-numeric(0)
for(i in 1:8086)if(all(is.na(kirezci.shape.scale.TWL.samples[i,,])))eliminate.kirezci<-c(eliminate.kirezci,i)
if(length(eliminate.kirezci)>0)kirezcisubsetID<-seq(length(dimnames(kirezci.shape.scale.TWL.samples)[[1]]))[-eliminate.kirezci]

eliminate.vousdoukas<-numeric(0)
for(i in 1:7472)if(all(is.na(vousdoukas.shape.scale.TWL.samples[i,,])))eliminate.vousdoukas<-c(eliminate.vousdoukas,i)
if(length(eliminate.vousdoukas)>0)vousdoukassubsetID<-seq(length(dimnames(vousdoukas.shape.scale.TWL.samples)[[1]]))[-eliminate.vousdoukas]
if(length(eliminate.vousdoukas)==0)vousdoukassubsetID<-seq(length(dimnames(vousdoukas.shape.scale.TWL.samples)[[1]]))

rmatched<-rl1r[rasmussensubsetID,]
kmatched<-rl1k[kirezcisubsetID,]
vmatched<-rl1v[vousdoukassubsetID,]

rfreq<-rasmussen.freq.probproj[,,rasmussensubsetID,]
kfreq<-kirezci.freq.probproj[,,kirezcisubsetID,]
vfreq<-vousdoukas.freq.probproj[,,vousdoukassubsetID,]

rasmussentokirezci<-find.match(rmatched,kmatched)
rasmussentovousdoukas<-find.match(rmatched,vmatched)
dim(rasmussentokirezci)  #187 2
dim(rasmussentovousdoukas)  #183 2

tt<-table(c(rasmussentokirezci[,1],rasmussentovousdoukas[,1]))
common<-as.numeric(names(tt)[tt==2])
rasmussentokirezciandvousdoukas<-cbind(common,rasmussentokirezci[,2][match(common,rasmussentokirezci[,1])],rasmussentovousdoukas[,2][match(common,rasmussentovousdoukas[,1])])




commonlocs.kirezci<-rasmussentokirezciandvousdoukas[,2]
commonlocs.vousdoukas<-rasmussentokirezciandvousdoukas[,3]
commonlocs.rasmussen<-rasmussentokirezciandvousdoukas[,1]


qq<-"q0.5"
whenone<-array(0,dim=c(length(commonlocs.kirezci),3,2))
dimnames(whenone)<-list(NULL,c("RP1r","RP1k","RP1v"),c("rasmussen","bvw"))


for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locr<-commonlocs.rasmussen[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){

        tempr<-rfreq[qq,,locr,slrm]
        tempk<-kfreq[qq,,lock,slrm]
        tempv<-vfreq[qq,,locv,slrm]


        whenone[ii,"RP1r",slrm]<-seq(length(tempr)-1)[tempr[-1]==1][1]
        whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]
        whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]
    }

}


whenone[is.na(whenone)]<-9


label.100to1.q0.5<-vote3.label.100to1.q0.5<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]
    tt<-table(temp)
    whichtt<-as.numeric(names(tt)[cumsum(tt)>3])
    vote3.label.100to1.q0.5[ii]<-label.100to1.q0.5[ii]<-ifelse(length(whichtt)>0,whichtt,median(temp))

}

###q0.05

qq<-"q0.05"
whenone<-array(0,dim=c(length(commonlocs.kirezci),3,2))
dimnames(whenone)<-list(NULL,c("RP1r","RP1k","RP1v"),c("rasmussen","bvw"))



for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locr<-commonlocs.rasmussen[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){

        tempr<-rfreq[qq,,locr,slrm]
        tempk<-kfreq[qq,,lock,slrm]
        tempv<-vfreq[qq,,locv,slrm]

        whenone[ii,"RP1r",slrm]<-seq(length(tempr)-1)[tempr[-1]==1][1]
        whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]
        whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]

         }
}

whenone[is.na(whenone)]<-9



vote3.label.100to1.q0.05<-label.100to1.q0.05<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]

    vote3.label.100to1.q0.05[ii]<-label.100to1.q0.05[ii]<-min(temp)

}



##q0.95


qq<-"q0.95"
whenone<-array(0,dim=c(length(commonlocs.kirezci),3,2))
dimnames(whenone)<-list(NULL,c("RP1r","RP1k","RP1v"),c("rasmussen","bvw"))


for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locr<-commonlocs.rasmussen[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){

        tempr<-rfreq[qq,,locr,slrm]
        tempk<-kfreq[qq,,lock,slrm]
        tempv<-vfreq[qq,,locv,slrm]

        whenone[ii,"RP1r",slrm]<-seq(length(tempr)-1)[tempr[-1]==1][1]
        whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]
        whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]

         }
}

whenone[is.na(whenone)]<-9



vote3.label.100to1.q0.95<-label.100to1.q0.95<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]

    vote3.label.100to1.q0.95[ii]<-label.100to1.q0.95[ii]<-max(temp)

}

colpal<-brewer.pal(11,"Spectral")[c(1:4,7:11)]
scenariolabels<-c("1.5C","2.0C","2.0C+","2.5C","3.0C","4.0C","5.0C","5.0C+")
names(scenariolabels)<-rasmussenScenarios

gridcoo.commonlocs<-rmatched[,c("longitude","latitude")][commonlocs.rasmussen,]

jpeg("change_in_freq_100to1_q05.jpg",quality=100,height=800,width=1200)
par(fg="gray70")
plot(gridcoo.commonlocs,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Central Estimate")
axis(1,at=seq(-200,200,by=50))
axis(2, at= seq(-100,100,by=20),las=1)
box()
abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[label.100to1.q0.5],cex=1.5)
points(seq(-170,170,length=9),rep(-87,9),pch=19,col=colpal,cex=2)
par(fg=1)
text(seq(-170,170,length=9),rep(-92,9),labels=c(scenariolabels,"none"),cex=1.5)
dev.off()





jpeg("change_in_freq_100to1_q005.jpg",quality=100,height=800,width=1200)
par(fg="gray70")
plot(gridcoo.commonlocs,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Lower Bound")
axis(1,at=seq(-200,200,by=50))
axis(2, at= seq(-100,100,by=20),las=1)
box()
abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[label.100to1.q0.05],cex=1.5)
points(seq(-170,170,length=9),rep(-87,9),pch=19,col=colpal,cex=2)

par(fg=1)
text(seq(-170,170,length=9),rep(-92,9),labels=c(scenariolabels,"none"),cex=1.5)
dev.off()



jpeg("change_in_freq_100to1_q095.jpg",quality=100,height=800,width=1200)
par(fg="gray70")
plot(gridcoo.commonlocs,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Upper Bound")
axis(1,at=seq(-200,200,by=50))
axis(2, at= seq(-100,100,by=20),las=1)
box()
abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[label.100to1.q0.95],cex=1.5)
points(seq(-170,170,length=9),rep(-87,9),pch=19,col=colpal,cex=2)

par(fg=1)
text(seq(-170,170,length=9),rep(-92,9),labels=c(scenariolabels,"none"),cex=1.5)
dev.off()





ll<-length(label.100to1.q0.5)
round(table(label.100to1.q0.5)/ll,dig=4)
#100 to 1, q=0.5
#1     2     3     4     6     7     8     9
#0.5419 0.1117 0.0279 0.0056 0.0112 0.0056 0.0950 0.2011

round(table(label.100to1.q0.05)/ll,dig=4)
#100 to 1, q=0.05

#    1      2      7
#0.9888 0.0056 0.0056

round(table(label.100to1.q0.95)/ll,dig=4)
#100 to 1, q=0.95

#     1      2      3      6      7      8      9
#0.0223 0.0112 0.0503 0.0223 0.0279 0.0615 0.8045




####now only match vousdoukas and kirezci

#only consider rows that have been sampled successfully

vousdoukastokirezci<-find.match(vmatched,kmatched)

commonlocs.kirezci<-vousdoukastokirezci[,2]
commonlocs.vousdoukas<-vousdoukastokirezci[,1]

qq<-"q0.5"
whenone<-array(0,dim=c(length(commonlocs.kirezci),2,2))
dimnames(whenone)<-list(NULL,c("RP1k","RP1v"),c("rasmussen","bvw"))



for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){


        tempk<-kfreq[qq,,lock,slrm]
        tempv<-vfreq[qq,,locv,slrm]


        whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]



        whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]

}
}

whenone[is.na(whenone)]<-9


vote2.label.100to1.q0.5<-label.100to1.q0.5<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]
    tt<-table(temp)
    whichtt<-as.numeric(names(tt)[cumsum(tt)>2])
    vote2.label.100to1.q0.5[ii]<-label.100to1.q0.5[ii]<-ifelse(length(whichtt)>0,whichtt,median(temp))


}

###q0.05

qq<-"q0.05"
whenone<-array(0,dim=c(length(commonlocs.kirezci),2,2))
dimnames(whenone)<-list(NULL,c("RP1k","RP1v"),c("rasmussen","bvw"))



for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]

    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){


                tempk<-kfreq[qq,,lock,slrm]
                tempv<-vfreq[qq,,locv,slrm]


                whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]



                whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]

            }
        }



whenone[is.na(whenone)]<-9



vote2.label.100to1.q0.05<-label.100to1.q0.05<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]

    vote2.label.100to1.q0.05[ii]<-label.100to1.q0.05[ii]<-min(temp)
    }



##q0.95


qq<-"q0.95"
whenone<-array(0,dim=c(length(commonlocs.kirezci),2,2))
dimnames(whenone)<-list(NULL,c("RP1k","RP1v"),c("rasmussen","bvw"))


for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){


                tempk<-kfreq[qq,,lock,slrm]
                tempv<-vfreq[qq,,locv,slrm]


                    whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]
                    whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]



            }
        }


whenone[is.na(whenone)]<-9


label.100to1.q0.95<-vote2.label.100to1.q0.95<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]

    vote2.label.100to1.q0.95[ii]<-label.100to1.q0.95[ii]<-max(temp)

}


gridcoo.commonlocs<-kmatched[,c("longitude","latitude")][commonlocs.kirezci,]

jpeg("evonly_change_in_freq_100to1_q05.jpg",quality=100,height=800,width=1200)
par(fg="gray70")
plot(gridcoo.commonlocs,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Central Estimate")
axis(1,at=seq(-200,200,by=50))
axis(2, at= seq(-100,100,by=20),las=1)
box()
abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[label.100to1.q0.5],cex=1.5)
points(seq(-170,170,length=9),rep(-87,9),pch=19,col=colpal,cex=2)

par(fg=1)
text(seq(-170,170,length=9),rep(-92,9),labels=c(scenariolabels,"none"),cex=1.5)
dev.off()





jpeg("evonly_change_in_freq_100to1_q005.jpg",quality=100,height=800,width=1200)
par(fg="gray70")
plot(gridcoo.commonlocs,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Lower Bound")
axis(1,at=seq(-200,200,by=50))
axis(2, at= seq(-100,100,by=20),las=1)
box()
abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[label.100to1.q0.05],cex=1.5)
points(seq(-170,170,length=9),rep(-87,9),pch=19,col=colpal,cex=2)

par(fg=1)
text(seq(-170,170,length=9),rep(-92,9),labels=c(scenariolabels,"none"),cex=1.5)
dev.off()




jpeg("evonly_change_in_freq_100to1_q095.jpg",quality=100,height=800,width=1200)
par(fg="gray70")
plot(gridcoo.commonlocs,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Upper Bound")
axis(1,at=seq(-200,200,by=50))
axis(2, at= seq(-100,100,by=20),las=1)
box()
abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[label.100to1.q0.95],cex=1.5)
points(seq(-170,170,length=9),rep(-87,9),pch=19,col=colpal,cex=2)

par(fg=1)
text(seq(-170,170,length=9),rep(-92,9),labels=c(scenariolabels,"none"),cex=1.5)
dev.off()



ll<-length(label.100to1.q0.5)
round(table(label.100to1.q0.5)/ll,dig=4)
#100 to 1, q=0.5
#     1      2      3      4      5      6      7      8      9
#0.4333 0.1026 0.0383 0.0043 0.0008 0.0373 0.0319 0.1178 0.2337


round(table(label.100to1.q0.05)/ll,dig=4)
#100 to 1, q=0.05
#     1      2      3      4      6      7      8
# 0.9842 0.0100 0.0019 0.0001 0.0008 0.0023 0.0005


round(table(label.100to1.q0.95)/ll,dig=4)
#100 to 1, q=0.95
#    1      2      3      4      5      6      7      8      9
# 0.0736 0.0154 0.0685 0.0019 0.0003 0.0144 0.0519 0.0943 0.6797


#now source("convolve_alltheway_fordatapub.R")

#now compare distributions after running the full convolution code and producing all3/all2
###one single plot:



jpeg(filedir%&%"pics/Fig1.jpg",quality=100,width=1000,height=1800)
par(mar=c(5,5,5,2)+0.1,mfrow=c(6,2))
temp1<-vote3.label.100to1.q0.5
temp2<-all3.label.100to1.q0.5
hist(temp1,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(a)",col=colpal,cex.main=3)
#main="Majority vote using median RP predictions from 6 methods",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")

hist(temp2,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(b)",col=colpal,cex.main=3)
#    main="Median RP predictions after full convolution",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")

temp1<-vote2.label.100to1.q0.5
temp2<-all2.label.100to1.q0.5
hist(temp1,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(c)",col=colpal,cex.main=3)
#main="Majority vote using median RP predictions from 6 methods",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")

hist(temp2,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(d)",col=colpal,cex.main=3)
#    main="Median RP predictions after full convolution",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")



temp1<-vote3.label.100to1.q0.05
temp2<-all3.label.100to1.q0.05
hist(temp1,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(e)",col=colpal,cex.main=3)
#    main="Minimum vote using 5th percentile of RP predictions from 6 methods",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")

hist(temp2,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(f)",col=colpal,cex.main=3)
#    main="5th percentile RP predictions after full convolution",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")


temp1<-vote2.label.100to1.q0.05
temp2<-all2.label.100to1.q0.05
hist(temp1,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(g)",col=colpal,cex.main=3)
#    main="Minimum vote using 5th percentile of RP predictions from 6 methods",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")

hist(temp2,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(h)",col=colpal,cex.main=3)
#    main="5th percentile RP predictions after full convolution",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")


temp1<-vote3.label.100to1.q0.95
temp2<-all3.label.100to1.q0.95
hist(temp1,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(i)",col=colpal,cex.main=3)
#    main="Maximum vote using 95th percentile RP predictions from 6 methods",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")

hist(temp2,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(j)",col=colpal,cex.main=3)
#    main="95th percentile RP predictions after full convolution",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")

temp1<-vote2.label.100to1.q0.95
temp2<-all2.label.100to1.q0.95
hist(temp1,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(k)",col=colpal,cex.main=3)
#    main="Maximum vote using 95th percentile RP predictions from 6 methods",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")

hist(temp2,breaks=seq(0.5,9.5),prob=T,axes=F,xlab="",ylab=length(temp1)%&%" sites",ylim=c(0,1),xlim=c(0.5,9.5),
     main="(l)",col=colpal,cex.main=3)
#    main="95th percentile RP predictions after full convolution",col=colpal)
axis(1,at=seq(1,9),labels=c(scenariolabels,"none"),las=1)
axis(2,las=2,at=seq(0,1,by=0.1),label=seq(0,100,by=10)%&%"%")

dev.off()





