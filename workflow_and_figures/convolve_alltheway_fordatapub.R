#this needs to be run in the same session after votingsystem_fordatapub.R

qs<-c(0.010, 0.050, 0.167, 0.500, 0.833, 0.950, 0.990, 0.995, 0.999)

vMC<-vousdoukas.shape.scale.TWL.samples[vousdoukassubsetID,,]
kMC<-kirezci.shape.scale.TWL.samples[kirezcisubsetID,,]
rMC<-rasmussen.shape.scale.TWL.samples[rasmussensubsetID,,]

whichdecade<-2100

#or you can run all of them to look at which decade...for(whichdecade in seq(2020,2100,by=10)){
    print(whichdecade)
load(paste0(filedir%&%"Rdatasets/RData_esl_slr_",whichdecade))
subproj1v<-proj1v[,,vousdoukassubsetID,]
subproj1k<-proj1k[,,kirezcisubsetID,]  #get these "subsetID" from convolve_fisherinfo.r
subproj1r<-proj1r[,,rasmussensubsetID,]

commonlocs.kirezci<-rasmussentokirezciandvousdoukas[,2]   #get this from convolve_fisherinfo.r
commonlocs.vousdoukas<-rasmussentokirezciandvousdoukas[,3]
commonlocs.rasmussen<-rasmussentokirezciandvousdoukas[,1]


all3f.int.probproj<-all3f.freq.probproj<-array(dim=c(11,9,length(commonlocs.kirezci)))
dimnames(all3f.int.probproj)<-dimnames(all3f.freq.probproj)<-list(c("mean","sd","q0.01","q0.05" , "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99" , "q0.995", "q0.999"),
                                    c("current",rasmussenScenarios),NULL)

set.seed(123)
for(ii in seq(length(commonlocs.kirezci))){
    lock<-commonlocs.kirezci[ii]
    locr<-commonlocs.rasmussen[ii]
    locv<-commonlocs.vousdoukas[ii]
    print(ii)

    eslv<-vmatched[locv,]    #these "matched" are from convolve_fisherinfo.r
    eslsamplev<-100*vMC[locv,,"TWL100"]  #randomized sample from shape/scale generation
    scalev<-vMC[locv,,"scale"]
    shapev<-vMC[locv,,"shape"]
    rpscurrentv<-RP(eslsamplev/100,shapev,scalev,eslv$thresh,eslv$lambda)

    eslk<-kmatched[lock,]
    eslsamplek<-100*kMC[lock,,"TWL100"]  #randomized sample from shape/scale generation
    scalek<-kMC[lock,,"scale"]
    shapek<-kMC[lock,,"shape"]
    rpscurrentk<-RP(eslsamplek/100,shapek,scalek,eslk$thresh,eslk$lambda)

    eslr<-rmatched[locr,]
    eslsampler<-100*rMC[locr,,"TWL100"]  #randomized sample from a normal approximation to 100-yr event
    scaler<-rMC[locr,,"scale"]
    shaper<-rMC[locr,,"shape"]
    rpscurrentr<-RP(eslsampler/100,shaper,scaler,eslr$thresh,eslr$lambda)


    eslsample<-c(eslsamplev,eslsamplek,eslsampler)

    rpscurrent<-c(rpscurrentv,rpscurrentk,rpscurrentr)

        all3f.int.probproj[,"current",ii]<-c(mean(eslsample),sqrt(var(eslsample)),
                                                           quantile(eslsample,prob=qs-qs[1],na.rm=TRUE))
        all3f.freq.probproj[,"current",ii]<-c(mean(rpscurrent),sqrt(var(rpscurrent)),
                                                            quantile(rpscurrent,prob=qs-qs[1],na.rm=TRUE))
            for(scenario in rasmussenScenarios){

                for(slrmethod in c("rasmussen","bvw")){
                slr.qs<-subproj1v[-c(1,2),scenario,locv,slrmethod]
                if(!all(is.na(slr.qs))){
                    slrsample<-unfold.quantiles(slr.qs)

                    ll<-length(slrsample)
                    slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                    eslsamplev<-eslsamplev[1:ll]


                    sssample<--slrsample+eslsamplev
                    rpssamplev<-RP(sssample/100,shapev[1:ll],scalev[1:ll],eslv$thresh,eslv$lambda)
                    sssamplev<-slrsample+eslsamplev
                    }
                else{
                    rpssamplev<-sssamplev<-rep(NA,ll)
                }


                slr.qs<-subproj1k[-c(1,2),scenario,lock,slrmethod]
                if(!all(is.na(slr.qs))){
                    slrsample<-unfold.quantiles(slr.qs)

                    ll<-length(slrsample)
                    slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                    eslsamplek<-eslsamplek[1:ll]


                    sssample<--slrsample+eslsamplek
                    rpssamplek<-RP(sssample/100,shapek[1:ll],scalek[1:ll],eslk$thresh,eslk$lambda)
                    sssamplek<-slrsample+eslsamplek
                    }
                else{
                    rpssamplek<-sssamplek<-rep(NA,ll)
                }

                slr.qs<-subproj1r[-c(1,2),scenario,locr,slrmethod]
                if(!all(is.na(slr.qs))){
                    slrsample<-unfold.quantiles(slr.qs)

                    ll<-length(slrsample)
                    slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                    eslsampler<-eslsampler[1:ll]


                    sssample<--slrsample+eslsampler
                    rpssampler<-RP(sssample/100,shaper[1:ll],scaler[1:ll],eslr$thresh,eslr$lambda)
                    sssampler<-slrsample+eslsampler
                    }
                else{
                    rpssampler<-sssampler<-rep(NA,ll)
                }

                if(slrmethod=="rasmussen"){
                    rpssample<-c(rpssamplev,rpssamplek,rpssampler)
                sssample<-c(sssamplev,sssamplek,sssampler)}

                else{
                    rpssample<-c(rpssample,rpssamplev,rpssamplek,rpssampler)
                    sssample<-c(sssample,sssamplev,sssamplek,sssampler)}
                }

                all3f.freq.probproj[,scenario,ii]<-c(mean(rpssample),sqrt(var(rpssample)),
                                                              quantile(rpssample,prob=qs-qs[1],na.rm=TRUE))

                all3f.int.probproj[,scenario,ii]<-c(mean(sssample),sqrt(var(sssample)),
                                                              quantile(sssample,prob=qs-qs[1],na.rm=TRUE))
            }
}



###



commonlocs.kirezci<-vousdoukastokirezci[,2]
commonlocs.vousdoukas<-vousdoukastokirezci[,1]


all2f.int.probproj<-all2f.freq.probproj<-array(dim=c(11,9,length(commonlocs.kirezci)))
dimnames(all2f.int.probproj)<-dimnames(all2f.freq.probproj)<-list(c("mean","sd","q0.01","q0.05" , "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99" , "q0.995", "q0.999"),
                                                                c("current",rasmussenScenarios),NULL)

set.seed(123)
for(ii in seq(length(commonlocs.kirezci))){
    lock<-commonlocs.kirezci[ii]
    locv<-commonlocs.vousdoukas[ii]
    print(ii)

    eslv<-vmatched[locv,]
    eslsamplev<-100*vMC[locv,,"TWL100"]  #randomized sample from shape/scale generation
    scalev<-vMC[locv,,"scale"]
    shapev<-vMC[locv,,"shape"]
    rpscurrentv<-RP(eslsamplev/100,shapev,scalev,eslv$thresh,eslv$lambda)

    eslk<-kmatched[lock,]
    eslsamplek<-100*kMC[lock,,"TWL100"]  #randomized sample from shape/scale generation
    scalek<-kMC[lock,,"scale"]
    shapek<-kMC[lock,,"shape"]
    rpscurrentk<-RP(eslsamplek/100,shapek,scalek,eslk$thresh,eslk$lambda)

    eslsample<-c(eslsamplev,eslsamplek)

    rpscurrent<-c(rpscurrentv,rpscurrentk)

    all2f.int.probproj[,"current",ii]<-c(mean(eslsample),sqrt(var(eslsample)),
                                                      quantile(eslsample,prob=qs-qs[1],na.rm=TRUE))
    all2f.freq.probproj[,"current",ii]<-c(mean(rpscurrent),sqrt(var(rpscurrent)),
                                                       quantile(rpscurrent,prob=qs-qs[1],na.rm=TRUE))

        for(scenario in rasmussenScenarios){


            for(slrmethod in c("rasmussen","bvw")){

            slr.qs<-subproj1v[-c(1,2),scenario,locv,slrmethod]
            if(!all(is.na(slr.qs))){
                slrsample<-unfold.quantiles(slr.qs)

                ll<-length(slrsample)
                slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                eslsamplev<-eslsamplev[1:ll]


                sssample<--slrsample+eslsamplev
                rpssamplev<-RP(sssample/100,shapev[1:ll],scalev[1:ll],eslv$thresh,eslv$lambda)
                sssamplev<-slrsample+eslsamplev
            }
            else{
                rpssamplev<-sssamplev<-rep(NA,ll)
            }


            slr.qs<-subproj1k[-c(1,2),scenario,lock,slrmethod]
            if(!all(is.na(slr.qs))){
                slrsample<-unfold.quantiles(slr.qs)

                ll<-length(slrsample)
                slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                eslsamplek<-eslsamplek[1:ll]


                sssample<--slrsample+eslsamplek
                rpssamplek<-RP(sssample/100,shapek[1:ll],scalek[1:ll],eslk$thresh,eslk$lambda)
                sssamplek<-slrsample+eslsamplek
            }
            else{
                rpssamplek<-sssamplek<-rep(NA,ll)
            }


if(slrmethod=="rasmussen"){
            rpssample<-c(rpssamplev,rpssamplek)
            sssample<-c(sssamplev,sssamplek)}
            else{
                rpssample<-c(rpssample,rpssamplev,rpssamplek)
                sssample<-c(sssample,sssamplev,sssamplek)}

}
            all2f.freq.probproj[,scenario,ii]<-c(mean(rpssample),sqrt(var(rpssample)),
                                                          quantile(rpssample,prob=qs-qs[1],na.rm=TRUE))

            all2f.int.probproj[,scenario,ii]<-c(mean(sssample),sqrt(var(sssample)),
                                                         quantile(sssample,prob=qs-qs[1],na.rm=TRUE))
        }
}



save(list=c(objects(pattern="all3f"),objects(pattern="all2f")),file=paste0(filedir,"Rdatasets/RData_fullconvolutions_fisherinfo_",whichdecade))

#}


cols<-brewer.pal(9,"Blues")






grid<-expand.grid(list(x=seq(-150,150,by=50), y=seq(-80,80,by=20)))
alldecades<-seq(2020,2100,by=10)
#for(whichdecade in alldecades){
    print(whichdecade)

    load(paste0(filedir,"Rdatasets/RData_fullconvolutions_fisherinfo_",whichdecade))
whenone<-matrix(0,nrow(rasmussentokirezciandvousdoukas),3)
dimnames(whenone)<-list(NULL,c("q0.05","q0.5","q0.95"))


for(ii in seq(nrow(rasmussentokirezciandvousdoukas))){

for(qq in c("q0.05","q0.5","q0.95")){


        temp<-all3f.freq.probproj[qq,,ii]

        whenone[ii,qq]<-seq(length(temp)-1)[temp[-1]==1][1]


     }
}

whenone[is.na(whenone)]<-9



all3.label.100to1.q0.05<-label.100to1.q0.05<-whenone[,"q0.05"]

all3.label.100to1.q0.5<-label.100to1.q0.5<-whenone[,"q0.5"]

all3.label.100to1.q0.95<-label.100to1.q0.95<-whenone[,"q0.95"]



colpal<-brewer.pal(11,"Spectral")[c(1:4,7:11)]
scenariolabels<-c("1.5C","2.0C","2.0C+","2.5C","3.0C","4.0C","5.0C","5.0C+")

gridcoo.commonlocs<-rmatched[,c("longitude","latitude")][rasmussentokirezciandvousdoukas[,1],]

jpeg(paste0(filedir,"pics/ED_Figure4_panel_a.jpg"),quality=100,height=800,width=1200)
par(fg="gray70")
plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Central Estimate")
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



jpeg(paste0(filedir,"pics/ED_Figure4_panel_c.jpg"),quality=100,height=800,width=1200)
par(fg="gray70")
plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Lower Bound")
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




jpeg(paste0(filedir,"pics/ED_Figure4_panel_e.jpg"),quality=100,height=800,width=1200)
par(fg="gray70")
plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Upper Bound")
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





ll<-length(label.100to1.q0.05)
round(table(label.100to1.q0.5)/ll,dig=4)
#100 to 1, q=0.5
#     1      2      3      4      5      6      7      9

#0.5810 0.1285 0.0838 0.0279 0.0335 0.0894 0.0279 0.0279
round(table(label.100to1.q0.05)/ll,dig=4)
#100 to 1, q=0.05
#     1      2      7

#0.9888 0.0056 0.0056
round(table(label.100to1.q0.95)/ll,dig=4)
#100 to 1, q=0.95
#     1      2      3      4      5      6      7      8      9

#0.1006 0.0279 0.0670 0.0056 0.0223 0.0894 0.0279 0.0503 0.6089


assign(paste0("all3.whenone.",whichdecade),whenone)

####now only match vousdoukas and kirezci




whenone<-matrix(0,nrow(vousdoukastokirezci),3)
dimnames(whenone)<-list(NULL,c("q0.05","q0.5","q0.95"))


for(ii in seq(nrow(vousdoukastokirezci))){

    for(qq in c("q0.05","q0.5","q0.95")){


            temp<-all2f.freq.probproj[qq,,ii]

            whenone[ii,qq]<-seq(length(temp)-1)[temp[-1]==1][1]

    }
    }
whenone[is.na(whenone)]<-9
whenten[is.na(whenten)]<-9

all2.label.100to1.q0.05<-label.100to1.q0.05<-whenone[,"q0.05"]

all2.label.100to1.q0.5<-label.100to1.q0.5<-whenone[,"q0.5"]

all2.label.100to1.q0.95<-label.100to1.q0.95<-whenone[,"q0.95"]



gridcoo.commonlocs<-vmatched[,c("longitude","latitude")][vousdoukastokirezci[,1],]

jpeg(paste0(filedir,"pics/ED_Figure4_panel_b.jpg"),quality=100,height=800,width=1200)
par(fg="gray70")
plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Central Estimate")
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





jpeg(paste0(filedir,"pics/ED_Figure4_panel_d.jpg"),quality=100,height=800,width=1200)
par(fg="gray70")
plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Lower Bound")
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




jpeg(paste0(filedir,"pics/ED_Figure4_panel_f.jpg"),quality=100,height=800,width=1200)
par(fg="gray70")
plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="100yr event becoming annual, \n Upper Bound")
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



ll<-length(label.100to1.q0.05)
round(table(label.100to1.q0.5)/ll,dig=4)
#100 to 1, q=0.5
#    1      2      3      4      5      6      7      8      9

#0.6596 0.0869 0.0527 0.0270 0.0291 0.0625 0.0301 0.0088 0.0433
round(table(label.100to1.q0.05)/ll,dig=4)
#100 to 1, q=0.05
#      1      2      3      4      5      6      7      8

#0.9570 0.0353 0.0043 0.0003 0.0004 0.0015 0.0007 0.0005
round(table(label.100to1.q0.95)/ll,dig=4)
#100 to 1, q=0.95
#    1      2      3      4      5      6      7      8      9

#0.1580 0.0280 0.0361 0.0218 0.0070 0.0625 0.1527 0.0730 0.4608




assign(paste0("all2.whenone.",whichdecade),whenone)


#}





#after running the above for all decades, need to put this all together in datasets that simply point at the first decade by which 100 -> 1

all3.gridcoo.commonlocs<-rmatched[,c("longitude","latitude")][rasmussentokirezciandvousdoukas[,1],]
all3.whenone<-all3.whenten<-array(dim=c(dim(all3.whenone.2020),9))
dimnames(all3.whenone)<-list(NULL,c("q0.05", "q0.5",  "q0.95"),alldecades)

for(yy in alldecades){
    all3.whenone[,,as.character(yy)]<-get(paste0("all3.whenone.",yy))


}




all2.gridcoo.commonlocs<-vmatched[,c("longitude","latitude")][vousdoukastokirezci[,1],]
all2.whenone<-array(dim=c(dim(all2.whenone.2020),9))
dimnames(all2.whenone)<-list(NULL,c("q0.05", "q0.5",  "q0.95"),alldecades)

for(yy in alldecades){
    all2.whenone[,,as.character(yy)]<-get(paste0("all2.whenone.",yy))

}



save(list=c("all2.whenone",
            "all3.whenone","all3.gridcoo.commonlocs","all2.gridcoo.commonlocs"),
     file=filedir%&%"Rdatasets/RData_labels_100to1_fullconvolution_fisherinfo_2020_2100")



####
alldecades<-seq(2020,2100,by=10)


dd<-dim(all3.whenone)
all3.whenone.byscenario<-array(dim=c(dd[1:2],8))
for(ii in 1:dd[1]){
    for(qq in 1:dd[2]){
        temp1<-all3.whenone[ii,qq,]
        if(all(temp1==9))all3.whenone.byscenario[ii,qq,]<-2110
        else{
            yy<-allyears[temp1!=9]
            temp1<-temp1[temp1!=9]
            utemp1<-unique(temp1)
            uyy<-get.early(temp1,yy)
            all3.whenone.byscenario[ii,qq,utemp1]<-uyy}

    }}

dimnames(all3.whenone.byscenario)<-list(NULL,c("q0.05","q0.5","q0.95"),rasmussenScenarios)

temp<-all3.whenone.byscenario
temp[,,1][is.na(temp[,,1])]<-2110
for(ii in 2:8)temp[,,ii][is.na(temp[,,ii])]<-temp[,,ii-1][is.na(temp[,,ii])]
all3.whenone.byscenario<-temp



###plot first decade by which the change in frequency happens, by scenario, for small set of locations
library(RColorBrewer)
colpal<-brewer.pal(11,"Spectral")[c(1:4,6:11)]
decades<-seq(2020,2110,by=10)
decadelabels<-c("2020","2030","2040","2050","2060","2070","2080","2090","2100")
grid<-expand.grid(list(x=seq(-150,150,by=50), y=seq(-80,80,by=20)))
gridcoo.commonlocs<-all3.gridcoo.commonlocs

for(scenario in 1:8){
scenariolabel<-rasmussenScenarios[scenario]
jpeg(paste0(filedir,"pics/FigurelikeS",11+scenario,"_smalllocationset.jpg"),quality=100,height=800,width=1200)
par(fg="gray70")
plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="Decade by which 100yr event becomes annual under , \n Central Estimate")
axis(1,at=seq(-200,200,by=50))
axis(2, at= seq(-100,100,by=20),las=1)
box()
abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[match(all3.whenone.byscenario[,2,scenario],decades)],cex=1.5)

points(seq(-170,170,length=10),rep(-87,10),pch=19,col=colpal,cex=2)
par(fg=1)
text(seq(-170,170,length=10),rep(-92,10),labels=c(decadelabels,"not yet"),cex=1.5)
dev.off()

jpeg(paste0(filedir,"pics/FigurelikeS",11+scenario,"_q0.05_smalllocationset.jpg"),quality=100,height=800,width=1200)
par(fg="gray70")
plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="Decade by which 100yr event becomes annual under , \n Central Estimate")
axis(1,at=seq(-200,200,by=50))
axis(2, at= seq(-100,100,by=20),las=1)
box()
abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[match(all3.whenone.byscenario[,1,scenario],decades)],cex=1.5)

points(seq(-170,170,length=10),rep(-87,10),pch=19,col=colpal,cex=2)
par(fg=1)
text(seq(-170,170,length=10),rep(-92,10),labels=c(decadelabels,"not yet"),cex=1.5)
dev.off()

jpeg(paste0(filedir,"pics/FigurelikeS",11+scenario,"_q0.95_smalllocationset.jpg"),quality=100,height=800,width=1200)
par(fg="gray70")
plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="Decade by which 100yr event becomes annual under , \n Central Estimate")
axis(1,at=seq(-200,200,by=50))
axis(2, at= seq(-100,100,by=20),las=1)
box()
abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[match(all3.whenone.byscenario[,3,scenario],decades)],cex=1.5)

points(seq(-170,170,length=10),rep(-87,10),pch=19,col=colpal,cex=2)
par(fg=1)
text(seq(-170,170,length=10),rep(-92,10),labels=c(decadelabels,"not yet"),cex=1.5)
dev.off()

}


#now for large set of locations (shown in SI)


dd<-dim(all2.whenone)
all2.whenone.byscenario<-array(dim=c(dd[1:2],8))
for(ii in 1:dd[1]){
    for(qq in 1:dd[2]){
        temp1<-all2.whenone[ii,qq,]
        if(all(temp1==9))all2.whenone.byscenario[ii,qq,]<-2110
        else{
            yy<-allyears[temp1!=9]
            temp1<-temp1[temp1!=9]
            utemp1<-unique(temp1)
            uyy<-get.early(temp1,yy)
            all2.whenone.byscenario[ii,qq,utemp1]<-uyy}

    }}

dimnames(all2.whenone.byscenario)<-list(NULL,c("q0.05","q0.5","q0.95"),rasmussenScenarios)

temp<-all2.whenone.byscenario
temp[,,1][is.na(temp[,,1])]<-2110
for(ii in 2:8)temp[,,ii][is.na(temp[,,ii])]<-temp[,,ii-1][is.na(temp[,,ii])]
all2.whenone.byscenario<-temp


###plot first decade, by scenario
library(RColorBrewer)
colpal<-brewer.pal(11,"Spectral")[c(1:4,6:11)]
decades<-seq(2020,2110,by=10)
decadelabels<-c("2020","2030","2040","2050","2060","2070","2080","2090","2100")
grid<-expand.grid(list(x=seq(-150,150,by=50), y=seq(-80,80,by=20)))
gridcoo.commonlocs<-all2.gridcoo.commonlocs

for(scenario in 1:8){
    scenariolabel<-rasmussenScenarios[scenario]
    jpeg(paste0(filedir,"pics/FigureS",11+scenario,".jpg"),quality=100,height=800,width=1200)
    par(fg="gray70")
    plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="Decade by which 100yr event becomes annual under , \n Central Estimate")
    axis(1,at=seq(-200,200,by=50))
    axis(2, at= seq(-100,100,by=20),las=1)
    box()
    abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
    map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
    points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[match(all2.whenone.byscenario[,2,scenario],decades)],cex=1.5)

    points(seq(-170,170,length=10),rep(-87,10),pch=19,col=colpal,cex=2)
    par(fg=1)
    text(seq(-170,170,length=10),rep(-92,10),labels=c(decadelabels,"not by 2100"),cex=1.5)
    dev.off()

    jpeg(paste0(filedir,"pics/FigureS",11+scenario,"_0.05.jpg"),quality=100,height=800,width=1200)
    par(fg="gray70")
    plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="Decade by which 100yr event becomes annual under , \n Central Estimate")
    axis(1,at=seq(-200,200,by=50))
    axis(2, at= seq(-100,100,by=20),las=1)
    box()
    abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
    map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
    points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[match(all2.whenone.byscenario[,1,scenario],decades)],cex=1.5)

    points(seq(-170,170,length=10),rep(-87,10),pch=19,col=colpal,cex=2)
    par(fg=1)
    text(seq(-170,170,length=10),rep(-92,10),labels=c(decadelabels,"not by 2100"),cex=1.5)
    dev.off()

    jpeg(paste0(filedir,"pics/FigureS",11+scenario,"_0.95.jpg"),quality=100,height=800,width=1200)
    par(fg="gray70")
    plot(grid,type="n",las=1,xlab="",ylab="",ylim=c(-90,90), xlim=c(-180,180),axes=F)#,cex.main=2),main="Decade by which 100yr event becomes annual under , \n Central Estimate")
    axis(1,at=seq(-200,200,by=50))
    axis(2, at= seq(-100,100,by=20),las=1)
    box()
    abline(v=seq(-150,150,by=50), h=seq(-80,80,by=20), col="gray78")
    map("world",add=TRUE,interior=F,fill=TRUE,col="gray70")
    points(gridcoo.commonlocs$longitude,gridcoo.commonlocs$latitude,pch=19,col=colpal[match(all2.whenone.byscenario[,3,scenario],decades)],cex=1.5)

    points(seq(-170,170,length=10),rep(-87,10),pch=19,col=colpal,cex=2)
    par(fg=1)
    text(seq(-170,170,length=10),rep(-92,10),labels=c(decadelabels,"not by 2100"),cex=1.5)
    dev.off()

}

cdfs<-array(dim=c(8,3,10))

dimnames(cdfs)<-list(rasmussenScenarios,c("lower","central","upper"),c(decadelabels,"notyet"))



for(i in 1:8){
    tt<-numeric(0)
    for(decade in decades)tt<-c(tt,sum(all2.whenone.byscenario[,2,i]==decade))
    cdfs[i,2,]<-cumsum(tt)/7283
    tt<-numeric(0)
    for(decade in decades)tt<-c(tt,sum(all2.whenone.byscenario[,1,i]==decade))
    cdfs[i,1,]<-cumsum(tt)/7283
    tt<-numeric(0)
    for(decade in decades)tt<-c(tt,sum(all2.whenone.byscenario[,3,i]==decade))
    cdfs[i,3,]<-cumsum(tt)/7283
}


temp<-cdfs[1,,]

for(scenario in rasmussenScenarios){
    temp<-rbind(temp,cdfs[scenario,,])
}
library(xtable)
xtable(temp)
