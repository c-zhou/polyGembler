
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-i", "--input", required = T, 
                    help="Input RData file.")
parser$add_argument("-n", "--nn", type="double", default=2,
                    help="Number of nearest neighbours [default %(default)s].")
parser$add_argument("-l", "--lod", type="double", default = 3,
                    help="LOD score threshold for grouping [default %(default)s].")
parser$add_argument("-r", "--rf", type="double", default = 0.38,
                    help="Recombination frequency threshold for grouping [default %(default)s].")
parser$add_argument("-o", "--output", required = T, 
                    help="Output file name.")
parser$add_argument("-c", "--concorde", required = T,
                    help="Concorde executable file path.")
parser$add_argument("-d", "--include",
                    help="External R package directories. Multiple paths are seperated by ':'.")
parser$add_argument("-t", "--tmpdir", required = F,
                    help="Temporary file directory.")
                    
args <- parser$parse_args()

in_RData = args$input
nn = args$nn
min_lod = args$lod
max_rf = args$rf
out_file = args$output
concorde_path = args$concorde
external_lib = args$include
tmpdir = args$tmpdir

if(!is.null(external_lib)) {
    libs = strsplit(external_lib,":")[[1]]
    .libPaths(c(.libPaths(),libs))
}

suppressPackageStartupMessages(library(TSP))

allPackages = rownames(installed.packages())

if(!is.null(tmpdir)) {
	if("unixtools" %in% allPackages) {
        suppressPackageStartupMessages(do.call('library', list("unixtools")))
    	set.tempdir(tmpdir)
    	cat(paste0("Setting TMPDIR ", tmpdir, "\n"))
    } else {
        warning(paste0("Package ",package," is not available. Using system default TMPDIR."))
	}
}

concorde_path(concorde_path)

initial.options <- commandArgs(trailingOnly = FALSE)
file.name <- "--file="
script.name <- sub(file.name, "", initial.options[grep(file.name, initial.options)])
source(paste(sep="/", dirname(script.name), "include.R"))

nn_joining(in_RData, out_file, nn, min_lod, max_rf)

## render all warning messages if has any
warnings()
