remove(list=ls())
filedir<-"./"
source(paste0(filedir,"Rcode/allfunctions.r"))

load(filedir%&%"Rdatasets/RData_esl_slr_2100")

rsites.gpd.parameters<-read.delim(filedir%&%"CSV/Rasmussen/gpd_parameters_uhawaii19.tsv")  #need the information about the record length in years
#which for rasmussen varies


allparams.collect<-numeric(0)
allyears<-numeric(0)
for(ii in seq(nrow(rl1r))){
    allparam<-rl1r[ii,]
    id<-allparam$SLRID
    mm<-match(id,rsites.gpd.parameters$PSMSL_ID)
    if(!is.na(mm)&allparam$shape>(-0.5)){  #this is a condition for the Fisher Information matrix application to be valid
        allparam2<-rsites.gpd.parameters[mm,]
        years<-allparam2$GPD_Record_End-allparam2$GPD_Record_Start + 1   #this info to be used in the Fisher Info formula
        allyears<-rbind(allyears,c(id,years))
        allscalesshapes<-simulate.bivariate(10000,allparam$scale,allparam$shape,allparam$lambda,years)
        allTWLs<-RL.N(100,allscalesshapes[,2],allscalesshapes[,1],allparam$thresh,allparam$lambda)
        qqs<-quantile(allTWLs,prob=c(0.05,0.50,0.95))
        allparams.collect<-rbind(allparams.collect,
                         c(allparam[-c(2,3,4,5,16,17,18)],quantile(allscalesshapes[,1],prob=c(0.05,0.50,0.95)),var(allscalesshapes[,1]),quantile(allscalesshapes[,2],prob=c(0.05,0.50,0.95)),var(allscalesshapes[,2]),qqs))

}}

dimnames(allparams.collect)<-list(NULL,c(names(rl1r)[-c(2,3,4,5,16,17,18)],"Scale05","Scale50","Scale95","VScale","Shape05","Shape50","Shape95","VShape","TWL05s","TWLs","TWL95s"))
rasmussen.allparams.collect<-allparams.collect
rasmussen.allyears<-allyears

allparams.collect<-numeric(0)
allyears<-numeric(0)
for(ii in seq(nrow(rl1k))){
    allparam<-rl1k[ii,]
    id<-allparam$SLRID
    if(allparam$shape>(-0.5)){
        years<-35  #this is common to all sites for the model-based ESL estimates
        allyears<-rbind(allyears,c(id,years))
        allscalesshapes<-simulate.bivariate(10000,allparam$scale,allparam$shape,allparam$lambda,years)
        allTWLs<-RL.N(100,allscalesshapes[,2],allscalesshapes[,1],allparam$thresh,allparam$lambda)
        qqs<-quantile(allTWLs,prob=c(0.05,0.50,0.95))
        allparams.collect<-rbind(allparams.collect,
                                 c(allparam[c(1,6,7,8,9)],qqs,mean(allTWLs)))

    }}

dimnames(allparams.collect)<-list(NULL,c(names(rl1k)[c(1,6,7,8,9)],"TWL05s","TWL50s","TWL95s","TWLs"))
kirezci.allparams.collect<-allparams.collect
kirezci.allyears<-allyears

allparams.collect<-numeric(0)
allyears<-numeric(0)
for(ii in seq(nrow(rl1v))){
    allparam<-rl1v[ii,]
    id<-allparam$SLRID
    if(allparam$shape>(-0.5)){
        years<-170    #this is assigned for the estimate to work, even if it is not the length of years used (Vousdoukas et al. use a different approach to the estimation)
        allyears<-rbind(allyears,c(id,years))
        allscalesshapes<-simulate.bivariate(10000,allparam$scale,allparam$shape,allparam$lambda,years)
        allTWLs<-RL.N(100,allscalesshapes[,2],allscalesshapes[,1],allparam$thresh,allparam$lambda)
        qqs<-quantile(allTWLs,prob=c(0.05,0.50,0.95))
        allparams.collect<-rbind(allparams.collect,
                                 c(allparam[c(1,6,7,8,9,10)],qqs,mean(allTWLs),sqrt(var(allTWLs))))

    }}

dimnames(allparams.collect)<-list(NULL,c(names(rl1v)[c(1,6,7,8,9,10)],"TWL05s","TWL50s","TWL95s","TWLs","TWLses"))

vousdoukas.allparams.collect<-allparams.collect
vousdoukas.allyears<-allyears


#rasmussen

nsample<-1000
rasmussen.shape.scale.TWL.samples<-array(dim=c(nrow(rl1r),nsample,4))
for(ii in 1:nrow(rl1r)){
    id<-rl1r$SLRID[ii]
    mm<-match(id,rasmussen.allyears[,1])
    if(!is.na(mm)){
    shape<-rl1r$shape[ii]
    scale<-rl1r$scale[ii]
    lambda<-rl1r$lambda[ii]
    thresh<-rl1r$thresh[ii]
    allscalesshapes<-simulate.bivariate(nsample,scale,shape,lambda,rasmussen.allyears[mm,2])
    allTWL100s<-RL.N(100,allscalesshapes[,2],allscalesshapes[,1],thresh,lambda)
    allTWL1s<-RL.N(1.1,allscalesshapes[,2],allscalesshapes[,1],thresh,lambda)
    rasmussen.shape.scale.TWL.samples[ii,,]<-cbind(allscalesshapes,allTWL100s,allTWL1s)
}}

dimnames(rasmussen.shape.scale.TWL.samples)<-list(rl1r$SLRID,NULL,c("scale","shape","TWL100","TWL1s"))
rl1r.diff.100.1<-rasmussen.shape.scale.TWL.samples[,,3]-rasmussen.shape.scale.TWL.samples[,,4]



rasmussen.int.probproj<-rasmussen.freq.probproj<-array(dim=dim(proj1r)+c(0,1,0,0))
dimnames(rasmussen.int.probproj)<-dimnames(rasmussen.freq.probproj)<-list(c("mean","sd",dimnames(proj1r)[[1]][-c(1,2)]),
                                                            c("current",dimnames(proj1r)[[2]]),dimnames(proj1r)[[3]],
                                                            dimnames(proj1r)[[4]])

qs<-c(0.010,0.050,0.167,0.500,0.833,0.950,0.990,0.995,0.999)
samplesizes<-c(40,117,333,333,117,40,5,4)

set.seed(123)
for(ii in seq(nrow(rl1r))){
    print(ii)
    id<-rl1r$SLRID[ii]
    if(!all(is.na(rasmussen.shape.scale.TWL.samples[as.character(id),,]))){
    eslsample<-rasmussen.shape.scale.TWL.samples[as.character(id),,"TWL100"]
    shapesample<-rasmussen.shape.scale.TWL.samples[as.character(id),,"shape"]
    scalesample<-rasmussen.shape.scale.TWL.samples[as.character(id),,"scale"]
    allparam<-rl1r[ii,]
    thresh<-rep(allparam$thresh,1000)
    lambda<-rep(allparam$lambda,1000)

    rpscurrent<-RP(eslsample,shapesample,scalesample,thresh,lambda)
    eslsample<-eslsample*100

    for(slrmethod in c("rasmussen","bvw")){
        rasmussen.int.probproj[,"current",ii,slrmethod]<-c(mean(eslsample),sqrt(var(eslsample)),
                                                    quantile(eslsample,prob=qs-qs[1]))
        rasmussen.freq.probproj[,"current",ii,slrmethod]<-c(mean(rpscurrent),sqrt(var(rpscurrent)),
                                                     quantile(rpscurrent,prob=qs-qs[1]))

        for(scenario in rasmussenScenarios){

            slr.qs<-proj1r[-c(1,2),scenario,ii,slrmethod]
            if(!all(is.na(slr.qs))){
                slrsample<-numeric(0)
                for(ss in seq(length(samplesizes))){
                    a<-slr.qs[ss]
                    b<-slr.qs[ss+1]
                    uu<-runif(samplesizes[ss],a,b)
                    slrsample<-c(slrsample,uu)}
                ll<-length(slrsample)
                slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                eslsample<-eslsample[1:ll]
                sssample<-slrsample+eslsample
                rasmussen.int.probproj[,scenario,ii,slrmethod]<-c(mean(sssample),sqrt(var(sssample)),
                                                           quantile(sssample,prob=qs-qs[1]))
                sssample<--slrsample+eslsample
                rpssample<-RP(sssample/100,shapesample[1:ll],scalesample[1:ll],thresh[1:ll],lambda[1:ll])

                rasmussen.freq.probproj[,scenario,ii,slrmethod]<-c(mean(rpssample),sqrt(var(rpssample)),
                                                            quantile(rpssample,prob=qs-qs[1],na.rm=TRUE))
            }
        }
        }
}
}


#kirezci

nsample<-1000
kirezci.shape.scale.TWL.samples<-array(dim=c(nrow(rl1k),nsample,4))
for(ii in 1:nrow(rl1k)){
    id<-rl1k$SLRID[ii]
        shape<-rl1k$shape[ii]
        if(shape>-0.5){
        scale<-rl1k$scale[ii]
        lambda<-rl1k$lambda[ii]
        thresh<-rl1k$thresh[ii]
        allscalesshapes<-simulate.bivariate(nsample,scale,shape,lambda,35)
        allTWL100s<-RL.N(100,allscalesshapes[,2],allscalesshapes[,1],thresh,lambda)
        allTWL1s<-RL.N(1.1,allscalesshapes[,2],allscalesshapes[,1],thresh,lambda)
        kirezci.shape.scale.TWL.samples[ii,,]<-cbind(allscalesshapes,allTWL100s,allTWL1s)
    }}

dimnames(kirezci.shape.scale.TWL.samples)<-list(rl1k$SLRID,NULL,c("scale","shape","TWL100","TWL1"))
rl1k.diff.100.1<-kirezci.shape.scale.TWL.samples[,,3]-kirezci.shape.scale.TWL.samples[,,4]


kirezci.int.probproj<-kirezci.freq.probproj<-array(dim=dim(proj1k)+c(0,1,0,0))
dimnames(kirezci.int.probproj)<-dimnames(kirezci.freq.probproj)<-list(c("mean","sd",dimnames(proj1k)[[1]][-c(1,2)]),
                                                              c("current",dimnames(proj1k)[[2]]),dimnames(proj1k)[[3]],
                                                              dimnames(proj1k)[[4]])

qs<-c(0.010,0.050,0.167,0.500,0.833,0.950,0.990,0.995,0.999)
samplesizes<-c(40,117,333,333,117,40,5,4)

set.seed(123)
for(ii in seq(nrow(rl1k))){
    print(ii)
    id<-rl1k$SLRID[ii]
    if(!all(is.na(kirezci.shape.scale.TWL.samples[as.character(id),,]))){
        eslsample<-kirezci.shape.scale.TWL.samples[as.character(id),,"TWL100"]
        shapesample<-kirezci.shape.scale.TWL.samples[as.character(id),,"shape"]
        scalesample<-kirezci.shape.scale.TWL.samples[as.character(id),,"scale"]
        allparam<-rl1k[ii,]
        thresh<-rep(allparam$thresh,1000)
        lambda<-rep(allparam$lambda,1000)

        rpscurrent<-RP(eslsample,shapesample,scalesample,thresh,lambda)
        eslsample<-eslsample*100

        for(slrmethod in c("rasmussen","bvw")){
            kirezci.int.probproj[,"current",ii,slrmethod]<-c(mean(eslsample),sqrt(var(eslsample)),
                                                         quantile(eslsample,prob=qs-qs[1]))
            kirezci.freq.probproj[,"current",ii,slrmethod]<-c(mean(rpscurrent),sqrt(var(rpscurrent)),
                                                          quantile(rpscurrent,prob=qs-qs[1]))

            for(scenario in rasmussenScenarios){

                slr.qs<-proj1k[-c(1,2),scenario,ii,slrmethod]
                if(!all(is.na(slr.qs))){
                    slrsample<-numeric(0)
                    for(ss in seq(length(samplesizes))){
                        a<-slr.qs[ss]
                        b<-slr.qs[ss+1]
                        uu<-runif(samplesizes[ss],a,b)
                        slrsample<-c(slrsample,uu)}
                    ll<-length(slrsample)
                    slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                    eslsample<-eslsample[1:ll]
                    sssample<-slrsample+eslsample
                    kirezci.int.probproj[,scenario,ii,slrmethod]<-c(mean(sssample),sqrt(var(sssample)),
                                                                quantile(sssample,prob=qs-qs[1]))
                    sssample<--slrsample+eslsample
                    rpssample<-RP(sssample/100,shapesample[1:ll],scalesample[1:ll],thresh[1:ll],lambda[1:ll])

                    kirezci.freq.probproj[,scenario,ii,slrmethod]<-c(mean(rpssample),sqrt(var(rpssample)),
                                                                 quantile(rpssample,prob=qs-qs[1],na.rm=TRUE))
                }
            }
        }
    }
}



####vousdoukas

nsample<-1000
vousdoukas.shape.scale.TWL.samples<-array(dim=c(nrow(rl1v),nsample,4))
for(ii in 1:nrow(rl1v)){
    id<-rl1v$SLRID[ii]
        shape<-rl1v$shape[ii]
        if(shape>-0.5){
        scale<-rl1v$scale[ii]
        lambda<-rl1v$lambda[ii]
        thresh<-rl1v$thresh[ii]
        allscalesshapes<-simulate.bivariate(nsample,scale,shape,lambda,170)
        allTWL100s<-RL.N(100,allscalesshapes[,2],allscalesshapes[,1],thresh,lambda)
        allTWL1s<-RL.N(1.1,allscalesshapes[,2],allscalesshapes[,1],thresh,lambda)
        vousdoukas.shape.scale.TWL.samples[ii,,]<-cbind(allscalesshapes,allTWL100s,allTWL1s)
    }}

dimnames(vousdoukas.shape.scale.TWL.samples)<-list(rl1v$SLRID,NULL,c("scale","shape","TWL100","TWL1"))
rl1v.diff.100.1<-vousdoukas.shape.scale.TWL.samples[,,3]-vousdoukas.shape.scale.TWL.samples[,,4]



vousdoukas.int.probproj<-vousdoukas.freq.probproj<-array(dim=dim(proj1v)+c(0,1,0,0))
dimnames(vousdoukas.int.probproj)<-dimnames(vousdoukas.freq.probproj)<-list(c("mean","sd",dimnames(proj1v)[[1]][-c(1,2)]),
                                                              c("current",dimnames(proj1v)[[2]]),dimnames(proj1v)[[3]],
                                                              dimnames(proj1v)[[4]])

qs<-c(0.010,0.050,0.167,0.500,0.833,0.950,0.990,0.995,0.999)
samplesizes<-c(40,117,333,333,117,40,5,4)

set.seed(123)
for(ii in seq(nrow(rl1v))){
    print(ii)
    id<-rl1v$SLRID[ii]
    if(!all(is.na(vousdoukas.shape.scale.TWL.samples[as.character(id),,]))){
        eslsample<-vousdoukas.shape.scale.TWL.samples[as.character(id),,"TWL100"]
        shapesample<-vousdoukas.shape.scale.TWL.samples[as.character(id),,"shape"]
        scalesample<-vousdoukas.shape.scale.TWL.samples[as.character(id),,"scale"]
        allparam<-rl1v[ii,]
        thresh<-rep(allparam$thresh,1000)
        lambda<-rep(allparam$lambda,1000)

        rpscurrent<-RP(eslsample,shapesample,scalesample,thresh,lambda)
        eslsample<-eslsample*100

        for(slrmethod in c("rasmussen","bvw")){
            vousdoukas.int.probproj[,"current",ii,slrmethod]<-c(mean(eslsample),sqrt(var(eslsample)),
                                                         quantile(eslsample,prob=qs-qs[1]))
            vousdoukas.freq.probproj[,"current",ii,slrmethod]<-c(mean(rpscurrent),sqrt(var(rpscurrent)),
                                                          quantile(rpscurrent,prob=qs-qs[1]))

            for(scenario in rasmussenScenarios){

                slr.qs<-proj1v[-c(1,2),scenario,ii,slrmethod]
                if(!all(is.na(slr.qs))){
                    slrsample<-numeric(0)
                    for(ss in seq(length(samplesizes))){
                        a<-slr.qs[ss]
                        b<-slr.qs[ss+1]
                        uu<-runif(samplesizes[ss],a,b)
                        slrsample<-c(slrsample,uu)}
                    ll<-length(slrsample)
                    slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                    eslsample<-eslsample[1:ll]
                    sssample<-slrsample+eslsample
                    vousdoukas.int.probproj[,scenario,ii,slrmethod]<-c(mean(sssample),sqrt(var(sssample)),
                                                                quantile(sssample,prob=qs-qs[1]))
                    sssample<--slrsample+eslsample
                    rpssample<-RP(sssample/100,shapesample[1:ll],scalesample[1:ll],thresh[1:ll],lambda[1:ll])

                    vousdoukas.freq.probproj[,scenario,ii,slrmethod]<-c(mean(rpssample),sqrt(var(rpssample)),
                                                                 quantile(rpssample,prob=qs-qs[1],na.rm=TRUE))
                }
            }
        }
    }
}

save(list=c(objects(pattern="probproj"),objects(pattern="shape.scale.TWL.samples") ),file=filedir%&%"Rdatasets/RData_convolutions_fisherinfo_2100")
#note that by loading a different decade's data the same procedure creates convolutions for earlier times.

#fixing the bug:

kirezci.freq.probproj[,4,,"bvw"]<-kirezci.freq.probproj[,3,,"bvw"]
kirezci.freq.probproj[,9,,"bvw"]<-kirezci.freq.probproj[,8,,"bvw"]
rasmussen.freq.probproj[,4,,"bvw"]<-rasmussen.freq.probproj[,3,,"bvw"]
rasmussen.freq.probproj[,9,,"bvw"]<-rasmussen.freq.probproj[,8,,"bvw"]
vousdoukas.freq.probproj[,4,,"bvw"]<-vousdoukas.freq.probproj[,3,,"bvw"]
vousdoukas.freq.probproj[,9,,"bvw"]<-vousdoukas.freq.probproj[,8,,"bvw"]
kirezci.int.probproj[,4,,"bvw"]<-kirezci.int.probproj[,3,,"bvw"]
kirezci.int.probproj[,9,,"bvw"]<-kirezci.int.probproj[,8,,"bvw"]
rasmussen.int.probproj[,4,,"bvw"]<-rasmussen.int.probproj[,3,,"bvw"]
rasmussen.int.probproj[,9,,"bvw"]<-rasmussen.int.probproj[,8,,"bvw"]
vousdoukas.int.probproj[,4,,"bvw"]<-vousdoukas.int.probproj[,3,,"bvw"]
vousdoukas.int.probproj[,9,,"bvw"]<-vousdoukas.int.probproj[,8,,"bvw"]

save(list=c(objects(pattern="probproj"),objects(pattern="shape.scale.TWL.samples") ),file=filedir%&%"Rdatasets/RData_convolutions_fisherinfo_2100")



mm1<-match(dimnames(rl1r.diff.100.1)[[1]],dimnames(rl1k.diff.100.1)[[1]])
mm2<-match(dimnames(rl1r.diff.100.1)[[1]],dimnames(rl1v.diff.100.1)[[1]])
subsetrasmussen.all3<-seq(length(mm1))[!is.na(mm1)&!is.na(mm2)]
subsetrl1ktorasmussen<-mm1[subsetrasmussen.all3]
subsetrl1vtorasmussen<-mm2[subsetrasmussen.all3]

diffrasmussen<-rl1r.diff.100.1[subsetrasmussen.all3,]
diffkirezci<-rl1k.diff.100.1[subsetrl1ktorasmussen,]
diffvousdoukas<-rl1v.diff.100.1[subsetrl1vtorasmussen,]

TWL100.TWL1.diff.all3<-cbind(diffrasmussen,diffkirezci,diffvousdoukas)
dimnames(TWL100.TWL1.diff.all3)<-list(dimnames(rl1r.diff.100.1)[[1]][subsetrasmussen.all3],NULL)
temp<-rl1r[subsetrasmussen.all3,c("longitude","latitude")]
temp1<-apply(TWL100.TWL1.diff.all3,1,median,na.rm=TRUE)
temp2<-sqrt(apply(TWL100.TWL1.diff.all3,1,var,na.rm=TRUE))
TWL100.TWL1.diff.all3.stats<-as.matrix(cbind(temp,temp1,temp2))
dimnames(TWL100.TWL1.diff.all3.stats)<-list(dimnames(rl1r.diff.100.1)[[1]][subsetrasmussen.all3],c("longitude","latitude","median","stdev"))

mm3<-match(dimnames(rl1k.diff.100.1)[[1]],dimnames(rl1v.diff.100.1)[[1]])
subsetkirezci.all2<-seq(length(mm3))[!is.na(mm3)]
subsetrl1vtokirezci<-mm3[!is.na(mm3)]

diffkirezci<-rl1k.diff.100.1[subsetkirezci.all2,]
diffvousdoukas<-rl1v.diff.100.1[subsetrl1vtokirezci,]

TWL100.TWL1.diff.all2<-cbind(diffkirezci,diffvousdoukas)
dimnames(TWL100.TWL1.diff.all2)<-list(dimnames(rl1k.diff.100.1)[[1]][subsetkirezci.all2],NULL)
temp<-rl1k[subsetkirezci.all2,c("longitude","latitude")]
temp1<-apply(TWL100.TWL1.diff.all2,1,median,na.rm=TRUE)
temp2<-sqrt(apply(TWL100.TWL1.diff.all2,1,var,na.rm=TRUE))
TWL100.TWL1.diff.all2.stats<-as.matrix(cbind(temp,temp1,temp2))
dimnames(TWL100.TWL1.diff.all2.stats)<-list(dimnames(rl1k.diff.100.1)[[1]][subsetkirezci.all2],c("longitude","latitude","median","stdev"))

rm(diffrasmussen,diffkirezci,diffvousdoukas,tempdiff)

save(list=objects(pattern="diff"),file=filedir%&%"Rdatasets/RData_diff_TWL100_TWL1")
