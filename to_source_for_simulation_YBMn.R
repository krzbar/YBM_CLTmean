library(TreeSim)
library(mvSLOUCH)

fsimul_XbarYBM<-function(N,n,runnum,part_n_directory,b_save_random_seeds,b_use_random_seed){
    sapply(1:N,function(i,n,runnum,part_n_directory,b_save_random_seeds,b_use_random_seed){
        cat("Doing simulation: ");cat(i);cat(" out of ");cat(N);cat(" for number of tips=");cat(n);cat(" ");cat(format(Sys.time(), "%T %F"));cat("\n");

	b_randomseed_file_exists<-FALSE
        Rseed_filename<-paste0(part_n_directory,"IterationRseed_YBMsimulstudy_n_",n,runnum,"_iter_",i,".RData")
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
	}

        Rseed_iter<-.Random.seed
	RNG_kind<-RNGkind()
	RNG_version<-getRversion()

	if((!b_use_random_seed)||(!b_randomseed_file_exists)){
    	    if (b_save_random_seeds){
        	save(.Random.seed,RNG_kind,RNG_version,file=Rseed_filename)
        	Rseed_iter<-.Random.seed
    	    }
	}

        iter<-i
        apephyltree<-sim.bd.taxa(n,1,1,0)[[1]]
        s2<-1;x0<-0
        ## in Laplace limit we have the root edge! So this effects the limit variance!
        x1<-rnorm(1,mean=0,sd=sqrt(s2*apephyltree$root.edge))
        mData<-simulBMProcPhylTree(apephyltree, matrix(x1,1,1), matrix(s2,1,1), dropInternal = TRUE)
        sample_mean<-mean(mData[,1])
        save(apephyltree,x1,mData,Rseed_iter,iter,sample_mean,file=paste0(part_n_directory,"IterationResult_YBMsimulstudy_n_",n,runnum,"_iter_",i,".RData"))
        sample_mean
    },n=n,runnum=runnum,part_n_directory=part_n_directory,b_save_random_seeds=b_save_random_seeds,b_use_random_seed=b_use_random_seed,simplify=TRUE)        
}

