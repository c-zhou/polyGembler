



.kosambi <- function(r) .25*log((1+2*r)/(1-2*r))
.haldane <- function(r) -.5*log(1-2*r)
.haldane_r <- function(d) .5*(1-exp(-2*d))
.kosambi_r <- function(d) .5*(exp(4*d)-1)/(exp(4*d)+1)

.cm_d <- function(r, mapfn=c("haldane","kosambi")) {
	match.arg(mapfn)
	mapfn = mapfn[1]
    r[r >= .5] = .4999999
    r[r <= -.5] = -.4999999
    if(mapfn=="haldane") {
    	return(-50*log(1-2*r))
    } else {
    	return(25*log((1+2*r)/(1-2*r)))
    }
}

.fit_mds_model <- function(mds_file, ncores=1) {

	ispcs = c(T,T,T,T,F,F,T,T,T,T,F,F,T,T,T,T,F,F)
	ndims = c(2,2,3,3,-1,-1,2,2,3,3,-1,-1,2,2,3,3,-1,-1)
	weightfns = c("lod","lod2","lod","lod2","lod","lod2","lod","lod2","lod","lod2","lod","lod2","lod","lod2","lod","lod2","lod","lod2")
	mapfns = c("kosambi","kosambi","kosambi","kosambi","kosambi","kosambi","haldane","haldane","haldane","haldane","haldane","haldane","none","none","none","none","none","none")
	nmods = 18
	
	cat(paste0("  Fitting ", nmods, " MDS models using ", ncores, " cores.\n"))
	
	registerDoParallel(ncores) ## register doParallel for child process
	maps <- foreach(i=1:nmods) %dopar% {
		ispc = ispcs[i]
		ndim = ndims[i]
		weightfn = weightfns[i]
		mapfn = mapfns[i]
		
		map_i = tryCatch({
					if(ispc) {
						calc.maps.pc(mds_file, ndim=ndim, weightfn=weightfn, mapfn=mapfn)
					} else {
						calc.maps.sphere(mds_file, weightfn=weightfn, mapfn=mapfn)
					}
				}, error = function(cond) {
					NULL
				})
		map_i
	}
	stopImplicitCluster()
	
	map = NULL
	stress = Inf
	
	for(i in 1:nmods) {
		map_i = maps[[i]]
		if(is.null(map_i)) next
		stress_i = if(ispcs[i]) {map_i$smacofsym$stress} else {map_i$smacofsphere$stress}
		
		if(stress_i<stress) {
			stress = stress_i
			map = map_i
			cat(paste0("  **Model ", i, " stress-1 value: ", stress_i, "\n"))
		} else {
			cat(paste0("    Model ", i, " stress-1 value: ", stress_i, "\n"))
		}
	}
	
	return(map)
}

preorder_mds <- function(clus, distanceMat, lodMat, ncores=1) {
	if(length(clus)<3) return(NA)

	wd <- tempdir()
	dir <- getwd()
	setwd(wd)
	on.exit(setwd(dir))
	
	temp_file <- basename(tempfile(tmpdir = wd))
	tmp_file_mds  <- paste(temp_file, ".txt", sep = "")

    dists = distanceMat[clus, clus]
    lods = lodMat[clus, clus]
    n = length(clus)
	
    sink(tmp_file_mds)
    cat(n); cat("\n")
    for(i in 1:n) {
        for(j in 1:n) {
            if(i>=j) next
            cat(i); cat("\t")
            cat(j); cat("\t")
            cat(dists[i,j]); cat("\t")
            cat(lods[i,j]); cat("\n")
        }
    }
    sink()

	map = .fit_mds_model(tmp_file_mds, ncores)
	unlink(tmp_file_mds)
	
	if(is.null(map)) stop("no model fitted.")
	
	return(as.numeric.factor(map$locimap$locus))
}

dist_mds <- function(clus, mds_file, distanceAll, lodAll, indexMat) {
    n = length(clus)
    names = c()
    for(i in 1:n) names=c(names,c(paste0(clus[i],"(+)"), paste0(clus[i],"(-)")))
	dists = matrix(NA, nrow=n*2, ncol=n*2, dimnames=list(names,names))
    lods = matrix(NA, nrow=n*2, ncol=n*2, dimnames=list(names,names))
	selfR = matrix(0, nrow=2,ncol=2)
    selfL = matrix(1000000, nrow=2,ncol=2)
    
	for(i in 1:n) {
        i1 = (i-1)*2+1
        i2 = i*2
        for(j in 1:n) {
            j1 = (j-1)*2+1
            j2 = j*2
            if(i==j) {
                dists[i1:i2,j1:j2] = selfR
				lods[i1:i2,j1:j2] = selfL
            } else if(i<j) {
                k = indexMat[clus[i],clus[j]]
				if(k==-1) return(NULL);
                r = distanceAll[k,]
                dists[i1:i2,j1:j2] = matrix(r,ncol=2,byrow=T)
                dists[j1:j2,i1:i2] = matrix(r,ncol=2,byrow=F)
				l = lodAll[k,]
                lods[i1:i2,j1:j2] = matrix(l,ncol=2,byrow=T)
                lods[j1:j2,i1:i2] = matrix(l,ncol=2,byrow=F)
            }
        }
    }
	
	n = n*2
    sink(mds_file)
    cat(n); cat("\n")
    for(i in 1:n) {
        for(j in 1:n) {
            if(i>=j) next
            cat(i); cat("\t")
            cat(j); cat("\t")
            cat(dists[i,j]); cat("\t")
            cat(lods[i,j]); cat("\n")
        }
    }
    sink()
}

ordering_mds <- function(clus, distanceAll, lodAll, indexMat, fid="", ncores=1) {

    if(length(clus)==0) return(NULL);
    if(length(clus)==1) return(list(order=clus,
                                   oriO=c(paste0(clus,"(+)"), paste0(clus,"(-)")),
                                   cost=0));
	wd <- tempdir()
	dir <- getwd()
	setwd(wd)
	on.exit(setwd(dir))
	
	temp_file <- basename(tempfile(tmpdir = wd))
	## using 'fid' to make sure the tmp files are not overwritten 
	## by each other although highly unlikely
	tmp_file_mds  <- paste(temp_file, fid, sep = "")
	
    dist_mds(clus, tmp_file_mds, distanceAll, lodAll, indexMat)
	
	map = .fit_mds_model(tmp_file_mds, ncores)
	unlink(tmp_file_mds)
	
	if(is.null(map)) return(NULL)
	
    oR = as.numeric.factor(map$locimap$locus)
	n = length(oR)
	trig = rep(T, n/2)
	
	c1 = map$length
	o = rep(NA, n/2)
	oriO = rep(NA, n)
	
	k = 0
	for(i in 1:n) {
		a = oR[i]
		b = ceiling(a/2)
		if(trig[b]) {
			k = k+1
			o[k] = clus[b]
			if(a%%2==1) {
				oriO[k*2-1] = paste0(clus[b],"(+)")
				oriO[k*2]   = paste0(clus[b],"(-)")
			} else {
				oriO[k*2-1] = paste0(clus[b],"(-)")
				oriO[k*2]   = paste0(clus[b],"(+)")
			}
		}
		trig[b] = F
	}

	 return(list(order = o,
                oriO = oriO,
                cost = c1));
}

INF_CM = .cm_d(0.5)

sum_finite <- function(x) {sum(x[is.finite(x)])}

as.numeric.factor <- function(f) {as.numeric(levels(f))[f]}

dist_tsp <- function(clus, distanceAll, indexMat, preorder=NA, nn=1) {
    n = length(clus)
    names = c()
    for(i in 1:n) names=c(names,c(paste0(clus[i],"(+)"),
                                  paste0(clus[i],"(-)")))
    names = c(names," __DUMMY__")
    d = matrix(NA, nrow=n*2+1, ncol=n*2+1, dimnames=list(names,names))
    selfM = matrix(0,nrow=2,ncol=2)
    selfM[1,2] <- selfM[2,1] <--Inf
    for(i in 1:n) {
        i1 = (i-1)*2+1
        i2 = i*2
        for(j in 1:n) {
            j1 = (j-1)*2+1
            j2 = j*2
            if(i==j) {
                d[i1:i2,j1:j2] = selfM
            } else if(i<j) {
                if(indexMat[clus[i],clus[j]]==-1) return(NULL);
                r = distanceAll[indexMat[clus[i],clus[j]],]
                d[i1:i2,j1:j2] = matrix(r,ncol=2,byrow=T)
                d[j1:j2,i1:i2] = matrix(r,ncol=2,byrow=F)
            }
        }
    }
	d = .cm_d(d)
	d[d<0] = -Inf
	MAX <- if(n*2+1<10) {2^16} else {2^31-1}
	
	if(any(is.na(preorder))) {	
		max_d = max(d[is.finite(d)])
		d = d+max_d+2
		max_d = max_d*2+2
		s = floor(log10(MAX/max_d/(n*2+1)))
		d = round(d*10^s)
		d[d==-Inf] = 1
		d[is.na(d)] = 2
		diag(d) = 0
    } else {
		for(i in 1:n) {
			inf = c()
			if(i-nn>1) inf = c(inf, 1:(i-nn-1))
			if(i+nn<n) inf = c(inf, (i+nn+1):n)
			inf = preorder[inf]
			inf = c(2*(inf-1)+1, 2*inf)
			k = preorder[i]
			k = c(2*(k-1)+1, 2*k)
			d[k, inf] = Inf
			d[inf, k] = Inf
		}
		max_d = max(d[is.finite(d)])
		d = d+max_d+2
		max_d = max_d*2+2
		sInf = MAX/(n*2+1)
		s = floor(log10(sInf/max_d))
		d = round(d*10^s)
		d[d==-Inf] = 1
		d[d==Inf] = sInf
		d[is.na(d)] = 2
		diag(d) = 0
	}
	
	d
}

.find_prog <- function(prog) {
	if(!is.null(concorde_path()))
      prog_path <- paste(concorde_path(), .Platform$file.sep, prog, sep ="")
    else prog_path <- prog
	prog_path
}

ordering_tsp <- function(clus, distanceAll, indexMat, method="concorde", preorder=NA, nn=1) {

    if(length(clus)==0) return(NA);
    if(length(clus)==1) return(list(order=clus,
                                   oriO=c(paste0(clus,"(+)"), paste0(clus,"(-)")),
                                   cost=0));
	wd <- tempdir()
	dir <- getwd()
	setwd(wd)
	on.exit(setwd(dir))
	
    d = dist_tsp(clus, distanceAll, indexMat, preorder, nn)
	
	#tour = solve_TSP(TSP(d), method, control = list(precision=0))
	temp_file <- basename(tempfile(tmpdir = wd))
	tmp_file_in  <- paste(temp_file, ".dat", sep = "")
	tmp_file_out <- paste(temp_file, ".sol", sep = "")
	write_TSPLIB(TSP(d), file = tmp_file_in, precision = 0)	
	system2(.find_prog("concorde"),
		args =  paste("-x -o", tmp_file_out, tmp_file_in),
    )
	if(!file.access(tmp_file_out) == 0)
		stop("Solving TSP failed.")
	tour <- scan(tmp_file_out, what = integer(0), quiet = TRUE)
	tour <- tour[-1] + 1L
	unlink(c(tmp_file_in, tmp_file_out, "file*", "Ofile*"))

    if(is.null(tour)) {
    	stop("Solving TSP failed.")
	} else{
		print("##Sovling TSP succeed.")
	}
	
    o = as.integer(tour)
    w = which(o==length(o))
    if(w==1 || w==length(o)) o = o[-w]
    else o = o[c((w+1):length(o),1:(w-1))]
    o = colnames(d)[o]
    c1 = length(tour)
    o1 = gsub("\\(\\+\\)|\\(\\-\\)","",o)
    n = length(o)
    o2 = o1[seq(1,n,2)]
    if(!all(o2==o1[seq(2,n,2)]))
        stop("TSP tour not valid.")
    
    return(list(order = o2,
                oriO = o,
                cost = c1));
}

.simply_write_files <- function(in_RData, in_map, out_file) {
	load(in_RData)
	
	dC = .read_map_file(in_map)$dC
	
	sink(paste0(out_file,".mct"))
	cat("group\tLG"); cat(1); cat("\n")
	cat(scaffs); cat("(+)\t0");cat("\n")
	cat(scaffs); cat("(-)\t");cat(dC);cat("\n")
	sink()
	
	sink(paste0(out_file,".par"))
	cat("-c ");cat(scaffs);cat("\n")
	sink()
	
	sink(paste0(out_file,".log"))
	cat("$contigs/scaffolds\n")
	cat(scaffs[1]); cat("\n");
	
	cat("\n$groups\n")
	cat("1\n")
	
	cat("\n$cm\n")
	cat(dC); cat("\n")
	
	cat("\n$order\n")
	cat("LG");cat(1);cat("\t\t")
	cat(scaffs);cat("\n")
	
	cat("\n$orientation\n")
	cat(scaffs); cat("(+)-")
	cat(scaffs); cat("(-)\n")
	sink()
}

.read_map_file <- function(in_map) {
	cDx = readLines(in_map)
	w = grep("^C",cDx)
	tC = gsub("^C ","",cDx[w])
	w = c(w,length(cDx)+1)
	dC = rep(NA, length(w)-1)
	for(i in 1:length(dC)) {
		rf = as.numeric(strsplit(cDx[w[i]+2],",")[[1]])
		j = w[i]+3
		while( j<w[i+1] ) {
			rf = rbind(rf, as.numeric(strsplit(cDx[j],",")[[1]]))
			j = j+1
		}
		means = rf
		if(!is.null(dim(rf))) means = apply(rf,2,mean)
		dC[i] = sum(.cm_d(means))
	}
	list(dC=dC, tC=tC)
}

nn_joining <- function(in_RData, out_file, nn=2, max_r=.haldane_r(.5)) {
	
    load(in_RData)
    diag(distanceMat) = Inf
	
	clusts = vector("list", length = n)
	for(i in 1:n) {
		d = distanceMat[i,]
        ux = sort(unique(d))
        nb = c()
        for(x in ux) {
			if(x>max_r) break
			nb = c(nb, which(d==x))
            if(length(nb)>=nn) break
        }
		clusts[[i]] = c(i, nb)
	}
	trig = rep(T, n)
	
	for(i in 1:(n-1)) {
		if(!trig[i]) next
		a = clusts[[i]]
		for(j in (i+1):n) {
			if(!trig[j]) next
			b = clusts[[j]]
			if(length(a)>length(b)&&all(b%in%a)) {
				trig[j] = F
			} else if(all(a%in%b)) {
				trig[i] = F
				break
			}
		}
	}
		
    sink(out_file)
    for(i in 1:n) {
		if(!trig[i]) next
		if(length(clusts[[i]])==1) {
			cat("-c ")
            cat(scaffs[clusts[[i]]])
            cat("\n")
            next
		}
			
		sink("/dev/null")
		o = ordering_tsp(clusts[[i]], distanceAll, indexMat)
		sink()
		oR = as.numeric(o$order); oO = o$oriO
        cat("-c ")
        cat(paste(scaffs[oR], collapse=":"))
        sepe = rep(NA, length(oR)-1)
        reverse = rep("false", length(oR))

        for(j in 1:(length(oR)-1)) {
            c_1 = oO[2*j]
            c_2 = oO[2*j+1]
            if(length(grep("(\\+)",c_1))>0) reverse[j]="true"
            sepe[j] =
            if(length(grep("(\\+)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],1]
            } else if(length(grep("(\\+)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],2]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],3]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],4]
            }
        }

        if(length(grep("(\\+)",oO[length(oO)]))>0) reverse[j+1]="true"
        cat(" -s "); cat(paste(sepe, collapse=":"))
        cat(" -r "); cat(paste(reverse, collapse=":"))
        cat("\n")
    }
    sink()
}

linkage_mapping <- function(in_RData, in_map, out_file, max_r=.haldane_r(0.5), make_group=TRUE, ncores=1) {

    load(in_RData)
    
	if(length(scaffs)==1) {
		.simply_write_files(in_RData, in_map, out_file);
		return()
	}
	
	max_d = .cm_d(max_r)
	
	if(make_group) {
		nng = .cm_d(distanceMat)
		diag(nng) = INF_CM
   		nng[nng>max_d] = INF_CM
    	nng = INF_CM-nng
		
		g=graph_from_adjacency_matrix(nng, mode = "undirected", weighted=T, diag=F)
    	clus=membership(cluster_infomap(g))
	
		#### check each cluster to remove chimeric joins
        #### could be misassembly
		maxc = max(clus)
		for(u in 1:maxc) {
			ci = which(clus==u)
			cn = length(ci)
			if(cn>2) {
				## make distance matrix
				dists = matrix(Inf, nrow=cn*2, ncol=cn*2)
				for(i in 1:(cn-1)) {
					i1 = (i-1)*2+1
					i2 = i*2
					for(j in (i+1):cn) {
						j1 = (j-1)*2+1
						j2 = j*2
						k = indexMat[ci[i],ci[j]]
						if(k==-1) stop("genetic mapping exit with errors!!!")
						r = distanceAll[k,]
						dists[i1:i2,j1:j2] = matrix(r,ncol=2,byrow=T)
						dists[j1:j2,i1:i2] = matrix(r,ncol=2,byrow=F)
					}
				}
				
				min_r = apply(dists,1,min)
				chims = which(min_r>max_r)
				if(length(chims)==0) next
				chims = unique(floor(chims/2+.5))
				for(chim in chims) {
					clus[ci[chim]] = max(clus)+1
				}
				print(paste0("#chimeric joins in linkage group ",u,": ",length(chims)))
			}
		}
		rm(nng)
	} else {
		clus = rep(1, length(scaffs))
	}
	
    clusts=list()
    for(i in 1:max(clus)) clusts[[i]] = which(clus==i)
    
	nc = rep(NA,length(clusts))
	for(i in 1:length(nc)) nc[i] = length(clusts[[i]])
	nco = order(nc, decreasing=T)
	
	mc = length(clusts)
	nn = sum(table(clus)>30)
	if(nn==0) nn=1
	
	if(ncores<=nn) {
		nt_p = ncores
		nt_c = rep(1, mc)
	} else {
		nt_p = nn
		nt = floor(ncores/nn)
		nt_c = rep(nt, mc)
		if(ncores>nn*nt) nt_c[1:(ncores-nn*nt)] = nt+1 
	}
	
	registerDoParallel(nt_p) ## register doParallel for parent process
	cat(paste0("####Ordering with MDS using ", ncores, " cores.\n"))
	o <- foreach(i=1:mc) %dopar% {
		ordering_mds(clusts[[nco[i]]], distanceAll, lodAll, indexMat, fid=paste0(".lg",i), nt_c[i])
	}
	stopImplicitCluster()
	
	for(i in 1:mc) {
		if(is.null(o[[i]])) {
			cat(paste0("####Linkage group ",i," MDS ordering failed. Ordering with TSP.\n"))
			o[[i]] = ordering_tsp(clusts[[nco[i]]], distanceAll, indexMat)
		}
	}
	
	mm = .read_map_file(in_map);
	dC = mm$dC
	tC = mm$tC
	
	lgCM = rep(NA,length(o))
    sink(paste0(out_file,".mct"))
    for(i in 1:length(o)) {
        cat("group\t"); cat("LG"); cat(i); cat("\n")
        oR = as.numeric(o[[i]]$order); oO = o[[i]]$oriO
        if(length(oR)<2) {
        	cat(scaffs[oR[1]]); cat("(+)\t0");cat("\n")
        	d = dC[tC==scaffs[oR[1]]]
			cat(scaffs[oR[1]]); cat("(-)\t");cat(d);cat("\n")
			cat("\n")
            next
        }

        d = 0
        nh = nchar(oO[1])
        cat(paste0(scaffs[as.numeric(substr(oO[1],1,nh-3))],substr(oO[1],nh-2,nh)))
        cat("\t"); cat(d); cat("\n")

        d = dC[tC==scaffs[oR[1]]]
        nh = nchar(oO[2])
        cat(paste0(scaffs[as.numeric(substr(oO[2],1,nh-3))],substr(oO[2],nh-2,nh)))
        cat("\t"); cat(d); cat("\n")

        for(j in 1:(length(oR)-1)) {
            c_1 = oO[2*j]
            c_2 = oO[2*j+1]
            dist_j =
            if(length(grep("(\\+)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],1]
            } else if(length(grep("(\\+)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],2]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],3]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],4]
            }

            d = d+.cm_d(dist_j)
            nh = nchar(oO[2*j+1])
            cat(paste0(scaffs[as.numeric(substr(oO[2*j+1],1,nh-3))],substr(oO[2*j+1],nh-2,nh)))
            cat("\t"); cat(d); cat("\n")

            d = d+dC[tC==scaffs[oR[j+1]]]
            nh = nchar(oO[2*j+2])
            cat(paste0(scaffs[as.numeric(substr(oO[2*j+2],1,nh-3))],substr(oO[2*j+2],nh-2,nh)))
            cat("\t"); cat(d); cat("\n")
        }

        cat("\n")
        lgCM[i] = d
    }
    sink()

    sink(paste0(out_file,".par"))
    for(i in 1:length(o)) {
        oR = as.numeric(o[[i]]$order); oO = o[[i]]$oriO
        cat("-c ")
        cat(paste(scaffs[oR], collapse=":"))
        if(length(oR)<2) {
            cat("\n")
            next
        }
        sepe = rep(NA, length(oR)-1)
        reverse = rep("false", length(oR))
        for(j in 1:(length(oR)-1)) {
            c_1 = oO[2*j]
            c_2 = oO[2*j+1]
            if(length(grep("(\\+)",c_1))>0) reverse[j]="true"
            sepe[j] =
            if(length(grep("(\\+)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],1]
            } else if(length(grep("(\\+)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],2]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],3]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                distanceAll[indexMat[oR[j],oR[j+1]],4]
            }
        }
        if(length(grep("(\\+)",oO[length(oO)]))>0) reverse[j]="true"
        cat(" -s "); cat(paste(sepe, collapse=":"))
        cat(" -r "); cat(paste(reverse, collapse=":"))
        cat("\n")
    }
    sink()

    sink(paste0(out_file,".log"))
    cat("$contigs/scaffolds\n")
    for(i in 1:length(scaffs)) {
        cat(scaffs[i]); 
        cat("\n");
    }
    cat("\n$groups\n")
    for(i in seq(1,length(scaffs),20)) {
        for(j in 0:19) {
            if( (i+j)<length(scaffs) ){
                cat(clus[i+j])
                cat(" ")
            }
        }
        cat("\n");
    }
    cat("\n$cm\n")
    cat(paste(lgCM,collapse=","))
    cat("\n")
    cat("\n$order\n")
    for(i in 1:length(o)) {
        cat("LG");cat(i);cat("\t\t")
        cat(paste(scaffs[as.numeric(o[[i]]$order)],collapse="-"));
        cat("\n")
    }
    cat("\n$orientation\n")
    for(i in 1:length(o)) {
        cat("LG");cat(i);cat("\t\t")
        oO = o[[i]]$oriO
        for(j in 1:length(oO)) {
            if(length(oO)>1) {
                nh=nchar(oO[j])
                oO[j]=paste0(scaffs[as.numeric(substr(oO[j],1,nh-3))],substr(oO[j],nh-2,nh))
            } else {
                oO[j]=scaffs[as.numeric(oO[j])]
            }
        }
        cat(paste(oO,collapse="-"));
        cat("\n")
    }
    sink()
	
	return(0)
}

