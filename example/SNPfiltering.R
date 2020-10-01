#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(updog))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(XNomial))

parser <- ArgumentParser()
parser$add_argument("-i", "--input", required = T,
                    help="Input VCF file")
parser$add_argument("-o", "--output", required = T,
                    help="Output VCF file")
parser$add_argument("-parents", "--parents", default = NULL, type="character",
                    help="Two parental IDs separated with \":\" (default: %(default)s)")
parser$add_argument("-p", "--ploidy", default = 2, type="double",
                    help="Output file (default: %(default)s)")
parser$add_argument("-f", "--minF1", default = 50, type="double",
                    help="Minimum number of F1 progeny (default: %(default)s)")
parser$add_argument("-het", "--minHetD", default = 5, type="double",
                    help="Minimum allele depth for hetrozygotes (default: %(default)s)")
parser$add_argument("-hom", "--minHomD", default = 3, type="double",
                    help="Minimum allele depth for homozygotes (default: %(default)s)")
parser$add_argument("-minMaxGL", "--minMaxGL", default = 0.8, type="double",
                    help="Minimum genotype likelihood for most likely genotype (default: %(default)s)")
parser$add_argument("-pv", "--pv", default = 0.001, type="double",
                    help="P-value threshold for segregation (default: %(default)s)")
parser$add_argument("-minAvgD", "--minAvgD", default = 3, type="double",
                    help="Minimum average allele depth (default: %(default)s)")
parser$add_argument("-maxAvgD", "--maxAvgD", default = 50, type="double",
                    help="Maximum average allele depth (default: %(default)s)")
parser$add_argument("-minF", "--minF", default = 0.05, type="double",
                    help="Minimum minor allele frequency (default: %(default)s)")
parser$add_argument("-maxThreads", "--maxThreads", default = 8, type="double",
                    help="Maximum threads used in parallel (default: %(default)s)")

args <- parser$parse_args()
infl = args$input
outf = args$output
parents = args$parents
ploidy = args$ploidy
minF1 = args$minF1
minHetD = args$minHetD
minHomD = args$minHomD
minMaxGL = args$minMaxGL
pv = args$pv
minAvgD = args$minAvgD
maxAvgD = args$maxAvgD
minF = args$minF
maxThreads = args$maxThreads
minD = min(minHetD, minHomD)

conf6=matrix(c(0.001666667,0.001666667,0.001666667,0.001666667,0.001666667,0.495833333,0.495833333,
               0.001666667,0.001666667,0.001666667,0.001666667,0.199333333,0.594666667,0.199333333,
               0.001666667,0.001666667,0.001666667,0.051083333,0.446416667,0.446416667,0.051083333,
               0.001666667,0.001666667,0.001666667,0.199333333,0.594666667,0.199333333,0.001666667,
               0.001666667,0.001666667,0.001666667,0.495833333,0.495833333,0.001666667,0.001666667,
               0.001666667,0.001666667,0.001666667,0.001666667,0.24875,0.495833333,0.24875,
               0.001666667,0.001666667,0.001666667,0.1005,0.397,0.397,0.1005,
               0.001666667,0.001666667,0.026375,0.24875,0.446416667,0.24875,0.026375,
               0.001666667,0.001666667,0.1005,0.397,0.397,0.1005,0.001666667,
               0.001666667,0.001666667,0.24875,0.495833333,0.24875,0.001666667,0.001666667,
               0.001666667,0.001666667,0.495833333,0.495833333,0.001666667,0.001666667,0.001666667,
               0.001666667,0.001666667,0.0412,0.238866667,0.436533333,0.238866667,0.0412,
               0.001666667,0.01155,0.120266667,0.36735,0.36735,0.120266667,0.01155,
               0.001666667,0.0412,0.238866667,0.436533333,0.238866667,0.0412,0.001666667,
               0.001666667,0.1005,0.397,0.397,0.1005,0.001666667,0.001666667,
               0.001666667,0.199333333,0.594666667,0.199333333,0.001666667,0.001666667,0.001666667,
               0.0041375,0.046141667,0.246279167,0.406883333,0.246279167,0.046141667,0.0041375,
               0.01155,0.120266667,0.36735,0.36735,0.120266667,0.01155,0.001666667,
               0.026375,0.24875,0.446416667,0.24875,0.026375,0.001666667,0.001666667,
               0.051083333,0.446416667,0.446416667,0.051083333,0.001666667,0.001666667,0.001666667,
               0.0412,0.238866667,0.436533333,0.238866667,0.0412,0.001666667,0.001666667,
               0.1005,0.397,0.397,0.1005,0.001666667,0.001666667,0.001666667,
               0.199333333,0.594666667,0.199333333,0.001666667,0.001666667,0.001666667,0.001666667,
               0.24875,0.495833333,0.24875,0.001666667,0.001666667,0.001666667,0.001666667,
               0.495833333,0.495833333,0.001666667,0.001666667,0.001666667,0.001666667,0.001666667), nrow=7)

conf4=matrix(c(0.0025,0.0025,0.0025,0.49625,0.49625,
               0.0025,0.0025,0.167083333,0.660833333,0.167083333,
               0.0025,0.0025,0.49625,0.49625,0.0025,
               0.0025,0.0025,0.249375,0.49625,0.249375,
               0.0025,0.084791667,0.413958333,0.413958333,0.084791667,
               0.0025,0.249375,0.49625,0.249375,0.0025,
               0.0025,0.49625,0.49625,0.0025,0.0025,
               0.029930556,0.221944444,0.49625,0.221944444,0.029930556,
               0.084791667,0.413958333,0.413958333,0.084791667,0.0025,
               0.167083333,0.660833333,0.167083333,0.0025,0.0025,
               0.249375,0.49625,0.249375,0.0025,0.0025,
               0.49625,0.49625,0.0025,0.0025,0.0025),nrow=5)

conf2=matrix(c(0.005,0.4975,0.4975,
               0.25125,0.4975,0.25125,
               0.4975,0.4975,0.005),nrow=3)

process <- function(dat, P1, P2, f1, ploidy=2, minF1=50, minMaxGL=0.8, pv=0.001, minD=3, minAvgD=3, maxAvgD=30, minF=0.05) {
    sizevec = snp$sizevec
    refvec = snp$refvec
    p1size = snp$p1size
    p1ref = snp$p1ref
    p2size = snp$p2size
    p2ref = snp$p2ref
    f1all = snp$f1
    
    uout <- flexdog(refvec=refvec,sizevec=sizevec,ploidy=ploidy,model="f1",p1ref=p1ref,p1size=p1size,p2ref=p2ref,p2size=p2size)

    maxpp = uout$maxpostprob
    f1pass = which(maxpp>=minMaxGL&sizevec>=minD)
    n = length(f1pass)
    if(n<minF1) return(list(pass=F, message=paste0("Genotype call too low (",n,").")))
    d = sum(sizevec[f1pass])
    if(d<n*minAvgD) return(list(pass=F, message="Mean allele depth too small."))
    if(d>n*maxAvgD) return(list(pass=F, message="Mean allele depth too large."))
    geno = ploidy-uout$geno
    gtt = table(geno[f1pass])
    if(sum(gtt>=minF*n)<2) return(list(pass=F, message="Minor allele frequency too small."))
    gcount = rep(0,ploidy+1)
    gcount[as.numeric(names(gtt))+1] = gtt
    nc = dim(conf)[2]
    pvals = rep(0, nc)
    for(i in 1:nc)
        pvals[i] = xmulti(gcount,conf[,i]*n)$pLLR
    maxPv = max(pvals)
    if(maxPv<pv) return(list(pass=F, message=paste0("Distort segregation (p-value=",maxPv,").")))

    shift = 0
    if(!is.null(P1)) shift=shift+1
    if(!is.null(P2)) shift=shift+1

    m_arr = rep(NA, 9+shift+n)
    m_arr[1] = snp$chrom
    m_arr[2] = snp$pos
    m_arr[3] = snp$id
    m_arr[4] = snp$ref
    m_arr[5] = snp$alt
    m_arr[6] = "."
    m_arr[7] = "PASS"
    m_arr[8] = paste0("DP=",sum(sizevec)+p1size+p2size,
                      ";AB=",round(uout$bias,digits=3),
                      ";ER=",round(uout$seq,digits=3),
                      ";OD=",round(uout$od,digits=3),
                      ";NS=",n,
                      ";SG=",paste(gcount,collapse=","),
                      ";PV=",round(maxPv,digits=3))
    m_arr[9] = "GT:AD:DP:GQ:PL"
    
    if(!is.null(P1)) {
        p1pl = phredQualityScore(ploidy-uout$par$p1geno)
        m_arr[10] = paste(c(genoStr(ploidy-uout$par$p1geno), paste(c(p1ref,p1size-p1ref),collapse=","), p1size, p1pl$gq, p1pl$pl), collapse=":")
    }
    if(!is.null(P2)) {
        p2pl = phredQualityScore(ploidy-uout$par$p2geno)
        m_arr[9+shift] = paste(c(genoStr(ploidy-uout$par$p2geno), paste(c(p2ref,p2size-p2ref),collapse=","), p2size, p2pl$gq, p2pl$pl), collapse=":")
    }
    postmat = -10*log10(uout$postmat)
    postmat = round(postmat-rep(apply(postmat,1,min),ploidy+1))
    postmat = postmat[,(ploidy+1):1]
    w = 0
    for(f in f1) {
        w = w+1
        i = which(f1all==f)
        if(length(i)==0||!i%in%f1pass) {
            m_arr[9+shift+w] = "."
        } else {
            f1pl = phredQualityScore(geno[i], postmat[i,])
            m_arr[9+shift+w] = paste(c(genoStr(geno[i]), paste(c(refvec[i],sizevec[i]-refvec[i]),collapse=","), sizevec[i], f1pl$gq, f1pl$pl), collapse=":")
        }
    }

    list(pass=T, message=paste(m_arr, collapse="\t"))
}

alleleDepth <- function(infoStr, indx) {
    ss = strsplit(infoStr, ":")[[1]]
    if(length(ss)<indx)
        return(c(0,0))
    ss = strsplit(ss[indx],",")[[1]]
    as.numeric(ss)
}

parseVCF <- function(vcfFile, P1=NULL, P2=NULL) {
    conn <- file(vcfFile, open="r")
    indx = list()
    dd = 0
    f1 = NULL
    nf1 = NA
    dat = list()
    while(length(ln <- readLines(conn,n=1,warn=FALSE))>0) {
        if(startsWith(ln, '##')) next
        if(startsWith(ln, '#')) {
            ss = strsplit(ln, "\\s+")[[1]]
            for(i in 10:length(ss)) {
                s = strsplit(ss[i],":")[[1]][1]
                if(is.null(indx[[s]])) {
                    indx[[s]] = i
                } else {
                    indx[[s]] = c(indx[[s]], i)
                }
            }
            f1 = setdiff(names(indx),c(P1,P2))
            nf1 = length(f1)
            next
        }
        ss = strsplit(ln, "\\s+")[[1]]
        if(dd==0) {
            ss1 = strsplit(ss[9], ":")[[1]]
            for(i in 1:length(ss1))
                if(ss1[i]=="AD")
                    dd=i
            if(dd==0) stop("no AD field!!!")
        }

        if(grepl(",",ss[5])) next  ##NOT biallelic SNPs
        
        p1ref = 0
        p1alt = 0
        if(!is.null(P1)) {
            for(i in indx[[P1]]) {
                d = alleleDepth(ss[i],dd)
                p1ref = p1ref+d[1]
                p1alt = p1alt+d[2]
            }
        }
        
        p2ref = 0
        p2alt = 0
        if(!is.null(P2)) {
            for(i in indx[[P2]]) {
                d = alleleDepth(ss[i],dd)
                p2ref = p2ref+d[1]
                p2alt = p2alt+d[2]
            }
        }
        
        f1ref = rep(0, nf1)
        f1alt = rep(0, nf1)
        for(i in 1:nf1) {
            f = f1[i]
            for(j in indx[[f]]) {
                d = alleleDepth(ss[j],dd)
                f1ref[i] = f1ref[i]+d[1]
                f1alt[i] = f1alt[i]+d[2]
            }
        }        

        snp = list()
        snp$chrom=ss[1]
        snp$pos=ss[2]
        snp$id=ss[3]
        snp$ref=ss[4]
        snp$alt=ss[5]
        snp$p1ref=p1ref
        snp$p1alt=p1alt
        snp$p2ref=p2ref
        snp$p2alt=p2alt
        snp$f1ref=f1ref
        snp$f1alt=f1alt
        
        dat = c(dat,list(snp))
    }
    close(conn)

    list(f1=f1, dat=dat)
}

phredQualityScore <- function(genotype, qs=NULL) {
    if(is.null(qs)) {
        qs = rep(30, ploidy+1)
        qs[genotype+1] = 0
        gq = 30
    } else {
        stopifnot(qs[genotype+1]==0)
        gq = sort(qs)[2]
    }
    pl = paste(qs, collapse=",")
    list(gq=gq, pl=pl)
}

genoStr <- function(geno) {
    stopifnot(geno>=0&&geno<=ploidy)
    strs = rep(1,ploidy)
    if(geno<ploidy)
        strs[1:(ploidy-geno)] = 0
    paste(strs,collapse="/")
}

writeVcfHeader <- function(conn, P1, P2, f1) {
    writeLines("##fileformat=VCFv4.0", conn)
    writeLines("##Updog=<ID=flexdog,Version=1.1.0,Description=\"Flexible genotyping for polyploids\">", conn)
    writeLines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", conn)
    writeLines("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">", conn)
    writeLines("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">", conn)
    writeLines("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">", conn)
    writeLines("##FORMAT=<ID=PL,Number=.,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">", conn)
    writeLines("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", conn)
    writeLines("##INFO=<ID=AB,Number=1,Type=Float,Description=\"Estimated Allelic Bias\">", conn)
    writeLines("##INFO=<ID=ER,Number=1,Type=Float,Description=\"Estimated Sequencing Error Rate\">", conn)
    writeLines("##INFO=<ID=OD,Number=1,Type=Float,Description=\"Estimated Overdispersion\">", conn)
    writeLines("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number Samples\">", conn)
    writeLines("##INFO=<ID=SG,Number=.,Type=Integer,Description=\"Segregation\">", conn)
    writeLines("##INFO=<ID=PV,Number=1,Type=Float,Description=\"Pvalue\">", conn)
    writeLines(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", P1, P2, f1), collapse="\t"), conn)
}


#### main program
P1=NULL
P2=NULL
if(!is.null(parents)) {
    parents=strsplit(parents,":")[[1]]
    P1 = parents[1]
    P2 = parents[2]
}
if(ploidy==2) {
    conf=conf2
} else if(ploidy==4) {
    conf=conf4
} else if(ploidy==6) {
    conf=conf6
} else {
    stop("segregation configuration not defined.")
}
data = parseVCF(infl,P1,P2)
f1 = data$f1
snps = data$dat

cl <- makeCluster(maxThreads)
registerDoParallel(cl)

M=length(snps)

collected <- foreach (i=1:M, .combine='c', .multicombine=TRUE, .packages=c('updog','XNomial')) %dopar% {
    # if(i%%1000==0) cat(paste0("#markers processed: ",i,"/",M,"\n"))
    
    snp = snps[[i]]
    f1ref=snp$f1ref
    f1alt=snp$f1alt
    sizevec=f1ref+f1alt
    hom=which((f1ref==0|f1alt==0)&sizevec>0)
    het=which(f1ref>0&f1alt>0)
    mm=c(hom[sizevec[hom]>=minHomD],het[sizevec[het]>=minHetD])
    
    snp$refvec=f1ref[mm]
    snp$sizevec=sizevec[mm]
    snp$p1size=snp$p1ref+snp$p1alt
    snp$p2size=snp$p2ref+snp$p2alt
    snp$f1=f1[mm] 

    if(length(mm)<minF1) {
        p <- NULL   
    } else {
        p <- tryCatch({
            process(snp, P1, P2, f1, ploidy, minF1, minMaxGL, pv, minD, minAvgD, maxAvgD, minF)
        }, error = function(cond) {
            cat(paste0("....Data point ",i," error throwed.\n"))
            NULL
        })
    }
    list(p)
}
stopCluster(cl)

vcfConn=file(outf, "w")
writeVcfHeader(vcfConn, P1, P2, f1)
m=0
for(i in 1:M) {
    p = collected[[i]]
    if(is.null(p)){
        m=m+1
        next
    }
    if(p$pass) {
        writeLines(p$message, vcfConn)
        cat(paste(c(i, snps[[i]]$id, "PASS=TRUE"), collapse="\t")); cat("\n")
    } else {
        m=m+1
        cat(paste(c(i, snps[[i]]$id, "PASS=FALSE", p$message), collapse="\t")); cat("\n")
    }
}
close(vcfConn)

cat(paste0("####markers filtered out: ",m,"/",M,"\n"))


