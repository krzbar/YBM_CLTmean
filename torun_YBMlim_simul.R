library(mvSLOUCH)
library(TreeSim)
library(parallel)

##============== user controlled variables =========================
runnum<-"" ## suffix of directories and files, can be used for distinguishing between numerous reruns
numcores<-8 ## number of cores used for the analyses
N<-10000
vn<-c(30,100,1000,10000,20000)
graphical_extension<-"pdf" ## in what format should the graphics be created in
b_use_random_seed<-FALSE ## use previous random seeds to replicate results
b_save_random_seeds<-TRUE ## if random seeds should be saved
##============== end of user controlled variables ==================


fsimul_Xbarlim<-function(N,cl){
    parSapply(cl,sapply(runif(N),function(x){ceiling((1+x)/(1-x))},simplify=TRUE),function(k){
	vExp<-rexp(2*k)
	vLap<-vExp[1:k]-vExp[(k+1):(2*k)]
	sum((1/sqrt(2*(1:k)))*vLap)
    },simplify=TRUE)
}

f_doYBMsimul<-function(n,N,runnum,main_directory,graphic_ext="pdf",b_save_random_seeds=TRUE,b_use_random_seed=FALSE){
    source("to_source_for_simulation_YBMn.R")
    n_directory<-paste0(main_directory,"n_",n,"/")
    dir.create(n_directory, showWarnings = FALSE)
    part_n_directory<-paste0(n_directory,"PartialResults/")
    dir.create(part_n_directory, showWarnings = FALSE)
    

    b_randomseed_file_exists<-FALSE
    Rseed_filename<-paste0(n_directory,"Rseed_BBMsimulstudy_n_",n,runnum,".RData")
    if(b_use_random_seed){##setup the random number generator if replicating results
        if(file.exists(Rseed_filename)){		    
	    load(Rseed_filename)
	    Rseed<-.Random.seed
	    rm(.Random.seed)
	    RNGkind(kind = RNG_kind[1], normal.kind = RNG_kind[2], sample.kind = RNG_kind[3])
	    RNGversion(RNG_version)     
	    assign('.Random.seed', Rseed, pos=.GlobalEnv)
	    b_randomseed_file_exists<-TRUE
	}			
    }else{rexp(1)}

    Rseed<-.Random.seed
    RNG_kind<-RNGkind()
    RNG_version<-getRversion()

    if((!b_use_random_seed)||(!b_randomseed_file_exists)){
	if (b_save_random_seeds){
	    save(.Random.seed,RNG_kind,RNG_version,file=Rseed_filename)	
	    Rseed<-.Random.seed
	}
    }	

    sink(paste0(n_directory,"progress",runnum,".txt"))
    vYBM<-fsimul_XbarYBM(N,n,runnum,part_n_directory,b_save_random_seeds,b_use_random_seed)
    sink()
    
    filename<-paste0(n_directory,"Hist_YBM_n_",n,runnum,".",graphic_ext)
    get(graphic_ext)(filename)
    hist(vYBM,main=paste0("n=",n),breaks=100,col="black",cex.main=2,cex.axis=1.5,freq=FALSE,xlab="",ylab="",xlim=c(-10,10),ylim=c(0,0.4))
    dev.off()    
    res<-list(vYBM=vYBM,Rseed=Rseed,RNG_kind=RNG_kind,RNGversion=RNG_version)
    save(vYBM,Rseed,RNG_kind,RNG_version,file=paste0(n_directory,"Result_YBMsimulstudy_n_",n,runnum,".RData"))
    res
}


main_directory<-paste0("./YBMlim_simulations",runnum,"/")
dir.create(main_directory, showWarnings = FALSE)
if (!exists("cl")){cl <- makeCluster(getOption("cl.cores", numcores),outfile="")}


rexp(1)


b_randomseed_file_exists<-FALSE
Rseed_filename<-paste0(main_directory,"Rseed_YBMlimitsimulstudy",runnum,".RData")
if(b_use_random_seed){##setup the random number generator if replicating results
    if(file.exists(Rseed_filename)){		    
	load(Rseed_filename)
	Rseed<-.Random.seed
	rm(.Random.seed)
	RNGkind(kind = RNG_kind[1], normal.kind = RNG_kind[2], sample.kind = RNG_kind[3])
	RNGversion(RNG_version)     
	assign('.Random.seed', Rseed, pos=.GlobalEnv)
	b_randomseed_file_exists<-TRUE
    }			
}else{rexp(1)}

Rseed<-.Random.seed
RNG_kind<-RNGkind()
RNG_version<-getRversion()

if((!b_use_random_seed)||(!b_randomseed_file_exists)){
    if (b_save_random_seeds){
	save(.Random.seed,RNG_kind,RNG_version,file=Rseed_filename)	
	Rseed<-.Random.seed
    }
}	

vYBMXbarlim<-fsimul_Xbarlim(N,cl)
filename<-paste0(main_directory,"Hist_YBMXbarlim",runnum,".",graphical_extension)
get(graphical_extension)(filename)
hist(vYBMXbarlim,main="Laplace mixture",breaks=100,col="black",cex.main=2,cex.axis=1.5,freq=FALSE,xlab="",ylab="",xlim=c(-10,10),ylim=c(0,0.4))
dev.off()    
save(vYBMXbarlim,Rseed,RNG_kind,RNG_version,file=paste0(main_directory,"Result_YBMlimitsimulstudy",runnum,".RData"))

parSapply(cl,vn,f_doYBMsimul,N=N,runnum=runnum,main_directory=main_directory,graphical_extension,b_save_random_seeds,b_use_random_seed)