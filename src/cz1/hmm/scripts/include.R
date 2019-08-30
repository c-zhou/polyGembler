

.kosambi <- function(r) .25*log((1+2*r)/(1-2*r))
.haldane <- function(r) -.5*log(1-2*r)
.haldane_r <- function(d) .5*(1-exp(-2*d))
.kosambi_r <- function(d) .5*(exp(4*d)-1)/(exp(4*d)+1)

.cm_d <- function(r) {
    r[r > .5] = .5
    r[r < -.5] = -.5
    25*log((1+2*r)/(1-2*r))
}

distTSP <- function(all, distanceAll, indexMat) {
    n = length(all)
    names = c()
    for(i in 1:n) names=c(names,c(paste0(all[i],"(+)"),
                                  paste0(all[i],"(-)")))
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
            if(i==j)
                d[i1:i2,j1:j2] = selfM
            else if(i<j) {
                if(indexMat[all[i],all[j]]==-1) return(NULL);
                r = distanceAll[indexMat[all[i],all[j]],]
                d[i1:i2,j1:j2] = matrix(r,ncol=2,byrow=T)
                d[j1:j2,i1:i2] = matrix(r,ncol=2,byrow=F)
            }
        }
    }
    m = max(d, na.rm=T)
    d = d+m+1e-6
    d[d==-Inf] = 1e-12
    d[is.na(d)] = 1e-6
    diag(d) = 0
    list(dMat = d,
         dMax = m,
         dSelf = 1e-12,
         dNa = 1e-6)
}

ordering <- function(all, distanceAll, indexMat, method="concorde") {

    if(length(all)==0) return(NA);
    if(length(all)==1) return(
                              list(order=all,
                                   oriO=all,
                                   cost=0,
                                   error=0));

    d = distTSP(all, distanceAll, indexMat)
    
    if("R.utils" %in% (.packages())) {
    	tour = NULL
    	att = 0
    	timeout = 300
		while(is.null(tour)&&att<12) {
    		tour = withTimeout(solve_TSP(TSP(d$dMat), method), timeout=timeout, onTimeout="silent")
    		att = att+1
    		if(is.null(tour)) warning(paste0("Solving TSP failed in attempt: ", att, ", in time limit: "+timeout+"s."))
    		timeout = timeout+300
    	}
    } else {
    	tour = solve_TSP(TSP(d$dMat), method)
    }
    
    if(is.null(tour)) {
    	stop("Solving TSP failed.")
	} else{
		print("##Sovling TSP succeed.")
	}
	 
    o = as.integer(tour)
    w = which(o==length(o))
    if(w==1 || w==length(o)) o = o[-w]
    else o = o[c((w+1):length(o),1:(w-1))]
    o = colnames(d$dMat)[o]
    c = tour_length(tour)-d$dMax*(length(all)-1)-d$dSelf*length(all)-d$dNa*2
    o2 = unique(gsub("\\(\\+\\)|\\(\\-\\)","",o))

    return(list(order = o2,
                oriO = o,
                dist = d,
                cost = c));
}

errorCount <- function(oo) {
    position_ = c()
    chr_ = c()
    all_splits_ = strsplit(oo,"_")
    for(i in 1:length(all_splits_)) {
        position_ = c(position_,as.numeric(all_splits_[[i]][2]))
        chr_ = c(chr_,all_splits_[[i]][1])
    }
    unique_chr_ = unique(chr_)

    n = length(oo)
    egN = n*(n-1)/2
    egn = egN
    eoN = 0;
    eon = 0;
    for(i in 1:length(unique_chr_)) {
        w = which(chr_==unique_chr_[i])
        e = .error(order(position_[w]))
        eoN = eoN+e[[2]]
        eon = eon+e[[1]]
        egn = egn-length(w)*(length(w)-1)/2
    }

    list(egn=egn, egN=egN, eon=eon, eoN=eoN)
}

.error <- function(oo) {
    n = length(oo)

    N = n*(n-1)/2
    e = 0
    for(i in 1:n)
        for(j in i:n)
            if(i!=j && oo[i]>oo[j]) e=e+1
    if(e>N/2) return(.error(oo[n:1]))
    list(e, N)
}

nmi <- function(cA, cB) {

    if(length(cA)!=length(cB)) return(NA)

    a = unique(cA)
    b = unique(cB)

    I = matrix(nrow=length(a), ncol=length(b))
    for(i in 1:length(a))
        for(j in 1:length(b))
            I[i, j] = length(intersect(
                                       which(cA==a[i]),
                                       which(cB==b[j])
                                       ))

    sA = apply(I,1,sum)
    sB = apply(I,2,sum)

    N = length(cA)

    iAB = 0
    for(i in 1:dim(I)[1])
        for(j in 1:dim(I)[2])
            if(I[i,j]!=0)
                iAB = iAB+I[i,j]*log2(I[i,j]*N/sA[i]/sB[j])

    hA = -sum(sA*log2(sA/N))
    hB = -sum(sB*log2(sB/N))

    2*iAB/(hA+hB)
}

fm <- function(cA, cB, beta=1) {
    if(length(cA)!=length(cB)) return(NA)

    a = unique(cA)
    b = unique(cB)

    rho = 0
    for(i in 1:length(a)) {
        wA = which(cA==a[i])
        for(j in 1:length(b)) {
            wB = which(cB==b[j])
            nAB = length(intersect(wA,wB))
            rho = rho+nAB*(nAB-1)/2
        }
    }

    ca = table(cA)
    nA = sum(ca*(ca-1)/2)
    cb = table(cB)
    nB = sum(cb*(cb-1)/2)

    p = rho/nA
    r = rho/nB

    (1+beta^2)*p*r/(beta^2*p+r)
}


two_point <- function(in_RData, twopoint_file, out_file, fix_rf=NA) {
	
	load(in_RData)
	diag(distanceMat) = Inf
	twopoint=matrix(as.character(unlist(read.table(twopoint_file))),ncol=2)
	n1 = dim(twopoint)[1]
	n2 = dim(twopoint)[2]
	clus=matrix(NA, nrow=n1, ncol=n2)
	for(i in 1:n1)
		for(j in 1:n2)
			clus[i,j] = which(scaffs==twopoint[i,j])
	
	sink(out_file)
	for(i in 1:dim(clus)[1]) {
		if(is.na(clus[i,2])) {
			cat("-c ")
			cat(scaffs[clus[i,1]])
			cat("\n")
			next
		}
		sink("/dev/null")
		o = ordering(clus[i,], distanceAll, indexMat)
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
			
			if(length(grep("(\\+)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
				sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],1]
			} else if(length(grep("(\\+)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
				sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],2]
			} else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
				sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],3]
			} else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
				sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],4]
			}
		}
		
		if(length(grep("(\\+)",oO[length(oO)]))>0) reverse[j+1]="true"
		if(is.na(fix_rf)) {
			cat(" -s "); cat(paste(sepe, collapse=":"))
		} else {
			cat(" -s "); cat(paste(rep(fix_rf, n2-1), collapse=":"))
		}
		cat(" -r "); cat(paste(reverse, collapse=":"))
		cat("\n")
	}
	sink()
}

nearest_neighbour_joining <- function(in_RData, out_file, nn=1, rf_thresh=.kosambi_r(.5)) {
	
    load(in_RData)
    diag(distanceMat) = Inf
    clus=matrix(NA, ncol=nn+1, nrow=0)
    for(i in 1:n) {
        w = sort(distanceMat[i,],index.return=T)
        x = w$x
        ix = w$ix
        ux = unique(x)
        s = 0
        ss = list()
        for(j in ux) {
            s = s+sum(x==j)
            ss[[length(ss)+1]] = which(x==j)
            if(s>=nn) break
        }
                
        p = i
        j = 1
        while(j<length(ss)) {
            p = c(p, ix[ss[[j]]])
            j = j+1
        }
        lp = length(p)
        cp = combs(ix[ss[[length(ss)]]],nn+1-lp)
        n = 0
        for(j in 1:dim(cp)[1]) {
            ng=sort(c(p,cp[j,]))
            if(any(distanceMat[combs(ng,2)]>rf_thresh))
                next
            clus=rbind(clus,ng)
            n = n+1
        }
        if(n==0) clus=rbind(clus,c(p,rep(NA,nn)))
    }
    clus=unique(clus)

    sink(out_file)
    for(i in 1:dim(clus)[1]) {
        if(is.na(clus[i,2])) {
            cat("-c ")
            cat(scaffs[clus[i,1]])
            cat("\n")
            next
        }
        sink("/dev/null")
        o = ordering(clus[i,], distanceAll, indexMat)
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
            
            if(length(grep("(\\+)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],1]
            } else if(length(grep("(\\+)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],2]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],3]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],4]
            }
        }

        if(length(grep("(\\+)",oO[length(oO)]))>0) reverse[j+1]="true"
        cat(" -s "); cat(paste(sepe, collapse=":"))
        cat(" -r "); cat(paste(reverse, collapse=":"))
        cat("\n")
    }
    sink()
}

traverse <- function(adj_matrix) {
    n = dim(adj_matrix)[1]
    ## need to be symmetric
    for(i in 1:(n-1))
        for(j in (i+1):n)
            if(adj_matrix[i,j]>0 || adj_matrix[j,i]>0)
                adj_matrix[i,j]=adj_matrix[j,i]=1

    visited = rep(F, n)
    clans = rep(NA, n)
    cn = 0
    for(i in 1:n) {
        if(visited[i]) next
        cn = cn+1
        clans[i] = cn
        visited[i] = T
        w = which(adj_matrix[i,]>0)
        while(length(w)>0) {
            w0 = w[1]
            w=w[-1]
            if(visited[w0]) next
            clans[w0] = cn
            visited[w0] = T
            w=c(w,which(adj_matrix[w0,]>0))
        }
    }

    clans
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

genetic_linkage_map <- function(in_RData, in_map, out_file, 
                                rf_thresh=.kosambi_r(.5), 
								make_group=TRUE, check=FALSE, ncore=1) {

    load(in_RData)
    
	if(length(scaffs)==1) {
		.simply_write_files(in_RData, in_map, out_file);
		return()
	}
	
    diag(distanceMat) = Inf

	if(make_group) {
    	nng=distanceMat
   		nng[nng>rf_thresh]=1
    	nng=1-nng

    	g=graph_from_adjacency_matrix(nng,mode = "undirected",
                                  	weighted=T, diag=F)
    	clus=membership(cluster_infomap(g))
	} else {
		clus = rep(1, length(scaffs))
	}
	
    all_clusters_=list()
    for(i in 1:max(clus)) all_clusters_[[i]] = which(clus==i)
    sa = matrix(nrow=0, ncol=2)

    for(i in 1:length(all_clusters_)) {
        all = all_clusters_[[i]]
        if(length(all)==1) {
            ass = getSingleAssignment(all, distanceMat)
            r = ass$rf
            if(r>rf_thresh) next
            c = ass$scaffs
            if(length(unique(clus[c]))==1)
                sa = rbind(sa,c(all,clus[c[1]]))
        }
    }
    
    if(dim(sa)[1]>0) {
        for(i in 1:dim(sa)[1]) {
            c = sa[i,1]
            g = sa[i,2]
            all_clusters_[[g]] = c(all_clusters_[[g]],c)
            clus[c] = g
        }
    }

    deleted_clusters_ = setdiff(1:length(all_clusters_), unique(clus))
    if(length(deleted_clusters_) > 0) {
        all_clusters_ = all_clusters_[-deleted_clusters_]
        for(i in 1:length(all_clusters_)) clus[all_clusters_[[i]]]=i
    }

    o = list()
    nc = rep(NA,length(all_clusters_))
    for(i in 1:length(nc)) nc[i] = length(all_clusters_[[i]])
    nco = order(nc,decreasing=T)
    
    if(ncore>1) {
    	registerDoParallel(ncore)
    	o <- foreach (i = 1:length(nco)) %dopar% {
 	 		ordering(all_clusters_[[nco[i]]], distanceAll, indexMat)
		}
    	stopImplicitCluster()
    } else {
    	for(i in 1:length(nco)) {
        	o[[i]] = ordering(all_clusters_[[nco[i]]], distanceAll, indexMat)
    	}
    }
	
	if(check) {
		print("checking linkage groups...")
		while(TRUE) {
			sepeAll = list()
			
			br = TRUE
			for(i in 1:length(o)) {
				oo = o[[i]]
				oR = as.numeric(oo$order)
				oO = o[[i]]$oriO
				sepe = rep(NA, length(oR)-1)
				if(length(sepe)==0) {
					all_clusters_[[length(all_clusters_)+1]] = as.numeric(oo$order)
					next
				}
				for(j in 1:(length(oR)-1)) {
					c_1 = oO[2*j]
					c_2 = oO[2*j+1]
					distance_j =
							if(length(grep("(\\+)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
								distanceAll[indexMat[oR[j],oR[j+1]],1]
							} else if(length(grep("(\\+)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
								distanceAll[indexMat[oR[j],oR[j+1]],2]
							} else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
								distanceAll[indexMat[oR[j],oR[j+1]],3]
							} else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
								distanceAll[indexMat[oR[j],oR[j+1]],4]
							}
					sepe[j] = distance_j
				}
				sepeAll[[length(sepeAll)+1]] = sepe
				br = br&all(sepe<=rf_thresh)
			}
			
			if(br) break
			
			all_clusters_ = list()
			for(i in 1:length(o)) {
				oo = o[[i]]
				oR = as.numeric(oo$order)
				k = c(1,1)
				new_clus = list()
				sepe = sepeAll[[i]]
				for(j in 1:(length(oR)-1)) {
					if(sepe[j]>rf_thresh) {
						a = length(new_clus)+1
						new_clus[[a]] = rep(NA,2)
						new_clus[[a]][1] = k[1]
						new_clus[[a]][2] = k[2]
						k[1] = j+1
						k[2] = j+1
					} else {
						k[2] = k[2]+1
					}
				}
				if(k[1]<=length(oR)) {
					a = length(new_clus)+1
					new_clus[[a]] = rep(NA,2)
					new_clus[[a]][1] = k[1]
					new_clus[[a]][2] = k[2]
				}
				
				oR = as.numeric(oo$order)
				for(j in 1:length(new_clus)) {
					all_clusters_[[length(all_clusters_)+1]] = 
							c(oR[new_clus[[j]][1]:new_clus[[j]][2]])
				}
			}
			
			o = list()
			nc = rep(NA,length(all_clusters_))
			for(i in 1:length(nc)) nc[i] = length(all_clusters_[[i]])
			nco = order(nc,decreasing=T)
			
			if(ncore>1) {
    			registerDoParallel(ncore)
    			o <- foreach (i = 1:length(nco)) %dopar% {
 	 				ordering(all_clusters_[[nco[i]]], distanceAll, indexMat)
				}
    			stopImplicitCluster()
    		} else {
				for(i in 1:length(nco)) {
					o[[i]] = ordering(all_clusters_[[nco[i]]], distanceAll, indexMat)
				}
			}
		}
	}
   
	mm = .read_map_file(in_map);
	dC = mm$dC
	tC = mm$tC
	
	lgCM = rep(NA,length(o))
    sink(paste0(out_file,".mct"))
    for(i in 1:length(o)) {
        cat("group\t"); cat("LG"); cat(i); cat("\n")
        or = as.numeric(o[[i]]$order)

        if(length(or)<2) {
            cat(scaffs[or[1]]); cat("\t0\n\n")
            next
        }
        
        oo = o[[i]]$oriO
        dist = .cm_d(o[[i]]$dist$dMat-o[[i]]$dist$dMax-o[[i]]$dist$dNa)
        for(c in or) {
            f=paste0(c,"(+)")
            r=paste0(c,"(-)")
            d0=dC[tC==scaffs[c]]
            dist[f,r]=dist[r,f]=d0#kosambi_r(d0)
        }

        nh = nchar(oo[1])
        cat(paste0(scaffs[as.numeric(substr(oo[1],1,nh-3))],substr(oo[1],nh-2,nh)))
        cat("\t0\n")
        d = 0
        for(j in 2:length(oo)) {
            nh = nchar(oo[j])
            cat(paste0(scaffs[as.numeric(substr(oo[j],1,nh-3))],substr(oo[j],nh-2,nh)))
            cat("\t")
            d = d+dist[oo[j-1], oo[j]]
            cat(d)
            cat("\n")
        }
    	lgCM[i] = d
        cat("\n")
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
            distance_j =
            if(length(grep("(\\+)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],1]
            } else if(length(grep("(\\+)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],2]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\+)",c_2))>0) {
                sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],3]
            } else if(length(grep("(\\-)",c_1))>0&&length(grep("(\\-)",c_2))>0) {
                sepe[j]=distanceAll[indexMat[oR[j],oR[j+1]],4]
            }
        }
        if(length(grep("(\\+)",oO[length(oO)]))>0) reverse[j]="true"
        cat(" -s "); cat(paste(sepe, collapse=":"))
        cat(" -r "); cat(paste(reverse, collapse=":"))
        # cat(" ## "); cat(paste(oO, collapse="-"))
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
        oo = o[[i]]$oriO
        for(j in 1:length(oo)) {
            if(length(oo)>1) {
                nh=nchar(oo[j])
                oo[j]=paste0(scaffs[as.numeric(substr(oo[j],1,nh-3))],substr(oo[j],nh-2,nh))
            } else {
                oo[j]=scaffs[as.numeric(oo[j])]
            }
        }
        cat(paste(oo,collapse="-"));
        cat("\n")
    }
    sink()

    list(group=clus, order=o)

}

getSingleAssignment <- function(id, distance, exclude=c(),include=c()) {
    a = distance[id,]
    if(length(include)>0) {
        idx = c()
        for(i in include) idx = c(idx, i)
        a[-idx] = Inf
    } else if(length(exclude)>0)
        for(i in exclude) a[i] = Inf
    list(rf=min(a),scaffs=which(a==min(a)))
}

