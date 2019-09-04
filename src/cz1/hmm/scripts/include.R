



.kosambi <- function(r) .25*log((1+2*r)/(1-2*r))
.haldane <- function(r) -.5*log(1-2*r)
.haldane_r <- function(d) .5*(1-exp(-2*d))
.kosambi_r <- function(d) .5*(exp(4*d)-1)/(exp(4*d)+1)

.cm_d <- function(r) {
    r[r >= .5] = .4999999
    r[r <= -.5] = -.4999999
    25*log((1+2*r)/(1-2*r))
}

writeMDS <- function(clus, mds, distanceMat, lodMat) {
	dist = distanceMat[clus, clus]
	lods = lodMat[clus, clus]
	n = length(clus)
	sink(mds)
	cat(n); cat("\n")
	for(i in 1:n) {
		for(j in 1:n) {
			if(i>=j) next
			cat(i); cat("\t")
			cat(j); cat("\t")
			cat(dist[i,j]); cat("\t")
			cat(lods[i,j]); cat("\n")
		}
	}
	sink()
}

INF_CM = .cm_d(0.5)

sum_finite <- function(x) {sum(x[is.finite(x)])}

as.numeric.factor <- function(f) {as.numeric(levels(f))[f]}

distTSP <- function(clus, distanceAll, indexMat, preorder=NULL, nn = 3) {
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
	
	if(is.null(preorder)) {	
		max_d = max(d, na.rm=T)
		d = d+max_d+2
		max_d = max_d*2+2
		s = floor(MAX/max_d/(n*2+1))
		d = round(d*s)
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
		sInf = floor(MAX/(n*2+1))
		s = sInf/max_d-1
		d = round(d*s)
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

ordering <- function(clus, distanceAll, indexMat, method="concorde", preorder=NULL, nn=3) {

    if(length(clus)==0) return(NA);
    if(length(clus)==1) return(
                              list(order=clus,
                                   oriO=clus,
                                   cost=0,
                                   error=0));
	wd <- tempdir()
	dir <- getwd()
	setwd(wd)
	on.exit(setwd(dir))
	
    d = distTSP(clus, distanceAll, indexMat, preorder, nn)
	
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
	unlink(c(tmp_file_in, tmp_file_out))

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
    c1 = tour_length(tour)
    o1 = gsub("\\(\\+\\)|\\(\\-\\)","",o)
    n = length(o)
    o2 = o1[seq(1,n,2)]
    if(!all(o2==o1[seq(2,n,2)]))
        stop("TSP tour not valid.")
    
    return(list(order = o2,
                oriO = o,
                cost = c1));
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
		if(is.na(fix_rf)) {
			cat(" -s "); cat(paste(.kosambi_r(sepe/100), collapse=":"))
		} else {
			cat(" -s "); cat(paste(rep(fix_rf, n2-1), collapse=":"))
		}
		cat(" -r "); cat(paste(reverse, collapse=":"))
		cat("\n")
	}
	sink()
}

nearest_neighbour_joining <- function(in_RData, out_file, nn=1, max_r=.kosambi_r(0.5)) {
	
    load(in_RData)
	
	max_d = .cm_d(max_r)
	distanceMat = .cm_d(distanceMat)
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
            if(any(distanceMat[combs(ng,2)]>max_d))
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
        cat(" -s "); cat(paste(.kosambi_r(sepe/100), collapse=":"))
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

genetic_linkage_map <- function(in_RData, in_map, out_file, max_r=.kosambi_r(0.5), make_group=TRUE, nn=3) {

    load(in_RData)
    
	if(length(scaffs)==1) {
		.simply_write_files(in_RData, in_map, out_file);
		return()
	}
	
	max_d = .cm_d(max_r)
	nng = .cm_d(distanceMat)
	
	if(make_group) {
    	diag(nng) = INF_CM
   		nng[nng>max_d] = INF_CM
    	nng = INF_CM-nng
		
		g=graph_from_adjacency_matrix(nng, mode = "undirected", weighted=T, diag=F)
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
            ass = getSingleAssignment(all, nng)
            d = ass$d
            if(d>max_d) next
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

	rm(nng)

	po = list()
	tmp_file = tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".txt")
	for(i in 1:length(all_clusters_)) {
		clust = all_clusters_[[i]]
		if(length(clust)>3) {
			writeMDS(clust, tmp_file, distanceMat, lodMat)
			maps = list()
			maps[[1]] <- calc.maps.pc(tmp_file, ndim=2, weightfn="lod", mapfn="kosambi")
			maps[[2]] <- calc.maps.pc(tmp_file, ndim=2, weightfn="lod2", mapfn="kosambi")
			maps[[3]] <- calc.maps.pc(tmp_file, ndim=3, weightfn="lod", mapfn="kosambi")
			maps[[4]] <- calc.maps.pc(tmp_file, ndim=3, weightfn="lod2", mapfn="kosambi")
			maps[[5]] <- calc.maps.sphere(tmp_file, weightfn="lod", mapfn="kosambi")
			maps[[6]] <- calc.maps.sphere(tmp_file, weightfn="lod2", mapfn="kosambi")
			lens = rep(NA, 6)
			for(j in 1:6) lens[j] = maps[[j]]$length
			map = maps[[which(lens==min(lens))[1]]]
			po[[i]] = as.numeric.factor(map$locimap$locus)
			unlink(tmp_file)
		} else {
			po[[i]] = c(1:length(clust))
		}
	}
	
	o = list()
    nc = rep(NA,length(all_clusters_))
    for(i in 1:length(nc)) nc[i] = length(all_clusters_[[i]])
    nco = order(nc, decreasing=T)
	for(i in 1:length(nco)) {
    	o[[i]] = ordering(all_clusters_[[nco[i]]], distanceAll, indexMat, preorder=po[[nco[i]]], nn=nn)
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

            d = d+dist_j
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
        cat(" -s "); cat(paste(.kosambi_r(sepe/100), collapse=":"))
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

    list(group=clus, order=o)
}

getSingleAssignment <- function(id, distance, exclude=c(), include=c()) {
    a = distance[id,]
    if(length(include)>0) {
        idx = c()
        for(i in include) idx = c(idx, i)
        a[-idx] = Inf
    } else if(length(exclude)>0)
        for(i in exclude) a[i] = Inf
    list(d=min(a),scaffs=which(a==min(a)))
}

