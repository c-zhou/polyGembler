
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-i", "--input", required = T,
                    help="Input RData file.")
parser$add_argument("-m", "--map", required = T,
                    help="Input MAP file.")
parser$add_argument("-r", "--distance", type="double", default = 0.381,
                    help="Recombination frequency threshold for grouping [default %(default)s].")
parser$add_argument("-1", "--no-group", action='store_true', 
                    help="Do not build linkage groups, all contigs into one linkage group [default FALSE].")		
parser$add_argument("-x", "--check-chimeric", action='store_true',
                    help="Check chimeric joins during grouping [default FALSE].")
parser$add_argument("-p", "--cores", type="double", default=1,
					help="Concorde executable file path.")			
parser$add_argument("-c", "--concorde", required = T,
                    help="Concorde executable file path.")
parser$add_argument("-l", "--include",
                    help="External R package locations. Multiple paths are seperated by ':'.")
parser$add_argument("-t", "--tmpdir", required = F,
                    help="Temporary file directory.")
parser$add_argument("-o", "--output", required = T,
					help="Output files name prefix.")
                    
args <- parser$parse_args()

in_RData = args$input
in_map = args$map
max_r = args$distance
no_group = args$no_group
check_chimeric = args$check_chimeric
ncores = args$cores
concorde_path = args$concorde
external_lib = args$include
tmpdir = args$tmpdir
out_file = args$output

if(!is.null(external_lib)) {
    libs = strsplit(external_lib,":")[[1]]
    .libPaths(c(.libPaths(),libs))
}

suppressPackageStartupMessages(library(TSP))
suppressPackageStartupMessages(library(MDSMap))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(doParallel))

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

linkage_mapping(in_RData, in_map, out_file, max_r, !no_group, check_chimeric, ncores)

## render all warning messages if has any
warnings()
