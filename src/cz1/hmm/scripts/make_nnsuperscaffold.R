.libPaths(c("/shares/common/groups/Group-Coin/c.zhou/old_home/R_lib","/opt/pkg/2015Q2/lib/R/library"))

require(combinat)
require(seriation)
require(argparse)
require(TSP)

library(combinat)
library(seriation)
library(argparse)
library(TSP)

concorde_path("/shares/common/groups/Group-Coin/c.zhou/old_home/sw/TSP_solver")

kosambi <- function(r) .25*log((1+2*r)/(1-2*r))
kosambi_r <- function(d) .5*(exp(4*d)-1)/(exp(4*d)+1)
haldane <- function(r) -.5*log(1-2*r)
haldane_r <- function(d) .5*(1-exp(-2*d))
#cm_d <- function(r) 50*log(1/(1-2*r))
cm_d <- function(r) {
        r[r > .5] = .5
        r[r < -.5] = -.5
        25*log((1+2*r)/(1-2*r))
}

getDistance <- function(a, b) {
    distanceMat[which(contigs==a),which(contigs==b)]
}

geneticmapping <- function(dat_, clus_, map_, contigs_, output="map.txt", m_="map.mct", param="map.par") {

        load(dat_)

    contig1 = read.table(contigs_)
        contigAll = as.character(unlist(contig1[,1]))
        sizes = rep(NA,length(contigs))
        for(i in 1:length(sizes))
                sizes[i]=as.numeric(as.character(contig1[contigs[i]==contigAll,3]))
        clus_in = readLines(clus_)

            clus = rep(NA,length(contigs))
    all_clusters_ = list()

        for(i in 3:length(clus_in)) {
                lg_i = strsplit(clus_in[i],"-")[[1]]
                s=rep(NA,length(lg_i))
                for(j in 1:length(lg_i)) {

                        w = which(lg_i[j]==contigs)
                        clus[w] = i-2
                        s[j] = w
                }
                all_clusters_[[i-2]] = s
        }


        thresh=kosambi_r(.5)
        nng=distanceMat
        nng[nng>thresh]=1
        nng=1-nng
        library(igraph)
        g=graph_from_adjacency_matrix(nng,mode = "undirected",weighted=T, diag=F)

        cl=cluster_infomap(g)
        # cl=cluster_label_prop(g)
        clus = membership(cl)
        all_clusters_=list()
        for(i in 1:max(clus)) all_clusters_[[i]] = which(clus==i)
    sa = matrix(nrow=0, ncol=2)
    for(i in 1:length(all_clusters_)) {
        all = all_clusters_[[i]]
        if(length(all)==1) {
            ass = getSingleAssignment(all, distanceMat)
            print(ass$rf)
            r  =ass$rf
            if(r>thresh) next
            c = ass$contigs
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
    for(i in 1:length(nco)) {
        o[[i]] = ordering(all_clusters_[[nco[i]]], distanceAll, indexMat)
    }
    #all_splits_ = strsplit(contigs,"_")
    #chr_ = c();
    #for(i in 1:length(all_splits_)) chr_ = c(chr_, all_splits_[[i]][1])
    #unique_chr_ = unique(chr_)
    #true_clus = rep(NA, length(chr_))
    #for(i in 1:length(unique_chr_)) true_clus[chr_==unique_chr_[i]] = i
    #eval = list(nmi=nmi(true_clus, clus), fm=fm(true_clus, clus))

    contigCOR = read.table(contigs_)[,c(1,4)]
    cDx = read.table(map_)
    cD = cDx[,c(1,dim(cDx)[2])]
    allRF = strsplit(as.character(unlist(cD[,2])),",")
    dC = rep(NA, length(allRF))
    for(i in 1:length(allRF))
        dC[i] = sum(cm_d(as.numeric(allRF[[i]])))
    tC = gsub("\\*","",cD[,1])

    sink(m_)
    for(i in 1:length(o)) {
        cat("group\t"); cat("LG"); cat(i); cat("\n")
        or = as.numeric(o[[i]]$order)

        if(length(or)<2) {
            cat(contigs[or[1]]); cat("\t0\n\n")
            next
        }

        oo = o[[i]]$oriO
        dist = cm_d(o[[i]]$dist$dMat-o[[i]]$dist$dMax-o[[i]]$dist$dNa)

        for(c in or) {
            f=paste0(c,"(+)")
            r=paste0(c,"(-)")
            d0=dC[tC==contigCOR[contigCOR[,1]==contigs[c],2]]
            dist[f,r]=dist[r,f]=d0#kosambi_r(d0)
        }
        nh = nchar(oo[1])
        cat(paste0(contigs[as.numeric(substr(oo[1],1,nh-3))],substr(oo[1],nh-2,nh)))
        cat("\t0\n")
        d = 0
        for(j in 2:length(oo)) {
            nh = nchar(oo[j])
            cat(paste0(contigs[as.numeric(substr(oo[j],1,nh-3))],substr(oo[j],nh-2,nh)))
            cat("\t")
            # d = d+kosambi(dist[oo[j-1], oo[j]])
            # d = d+cm_d(dist[oo[j-1], oo[j]])
            d = d+dist[oo[j-1], oo[j]]
            cat(d)
            cat("\n")
        }
        cat("\n")
    }
    sink()

    sink(param)
    for(i in 1:length(o)) {
        oR = as.numeric(o[[i]]$order); oO = o[[i]]$oriO
        cat("-c ")
        cat(paste(contigs[oR], collapse=":"))
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
        cat(" -S "); cat(paste(sepe, collapse=":"))
        cat(" -R "); cat(paste(reverse, collapse=":"))
        # cat(" ## "); cat(paste(oO, collapse="-"))
        cat("\n")
    }
    sink()

    sink(output)
    cat("$contigs\n")
    for(i in 1:length(contigs)) {cat(contigs[i]); cat("\n");}
    cat("\n$groups\n")
    for(i in seq(1,length(contigs),20)) {
        for(j in 0:19) {
            if( (i+j)<length(contigs) ){
                cat(clus[i+j])
                cat(" ")
            }
        }
        cat("\n");
    }
    cat("\n$order\n")
    L = rep(0, length(o))
    for(i in 1:length(o)) {
        L[i] = sum(sizes[as.numeric(o[[i]]$order)])
        cat("LG");cat(i);cat("[");cat(L[i]);cat("]");cat("\t\t")
        cat(paste(contigs[as.numeric(o[[i]]$order)],collapse="-"));
        cat("\n")
    }
    cat("\n$orientation\n")
    for(i in 1:length(o)) {
        cat("LG");cat(i);cat("[");cat(L[i]);cat("]");cat("\t\t")
        oo = o[[i]]$oriO
        for(j in 1:length(oo)) {
            if(length(oo)>1) {
                nh=nchar(oo[j])
                oo[j]=paste0(contigs[as.numeric(substr(oo[j],1,nh-3))],substr(oo[j],nh-2,nh))
            } else {
                oo[j]=contigs[as.numeric(oo[j])]
            }
        }
        cat(paste(oo,collapse="-"));
        cat("\n")
    }
    cat("\n$total\n")
    cat(sum(L)); cat("\n")
    sink()

    list(group=clus, order=o, eval=eval)

    save(list=ls(), file=paste0(output,".RData"))
}

getSingleAssignment <- function(id, distance, exclude=c(),include=c()) {
    a = distance[id,]
    if(length(include)>0) {
        idx = c()
        for(i in include) idx = c(idx, i)
        a[-idx] = Inf
    } else if(length(exclude)>0)
        for(i in exclude) a[i] = Inf
    list(rf=min(a),contigs=which(a==min(a)))
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

ordering <- function(all, distanceAll, indexMat) {

    if(length(all)==0) return(NA);
    if(length(all)==1) return(
                              list(order=all,
                                   oriO=all,
                                   cost=0,
                                   error=0));

    d = distTSP(all, distanceAll, indexMat)
    tour = solve_TSP(TSP(d$dMat), method="concorde")
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
                cost = c,
                error=errorCount(o2,contigs)
                ));
}

errorCount <- function(oo, contigs) {
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
    #tmp = rep(0,n)
    #for(i in 1:n) tmp[oo[i]]=i
    #tmp = tmp[ro]

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

parser <- ArgumentParser()
parser$add_argument("--dat",
                    help="Specifies the file RData file to load");
parser$add_argument("--clus",
                    help="Specifies the file contains the linkage groups information");
parser$add_argument("--map",
                    help="Specifies the file contains the contig size estimation");
parser$add_argument("--contig",
                    help="Specifies the file contains the contig information");
parser$add_argument("--o",
                    help="Specifies the output file name");
parser$add_argument("--m",
                    help="Specifies the output genetic map file")
parser$add_argument("--p",
                    help="Specifies the output parameter file for multipoint analysis")

args <- parser$parse_args()
dat_ = args$dat
clus_ = args$clus
map_ = args$map
contig_ = args$contig
o_ = args$o
m_ = args$m
p_ = args$p

gm = geneticmapping(dat_, clus_, map_, contig_, o_, m_, p_)
