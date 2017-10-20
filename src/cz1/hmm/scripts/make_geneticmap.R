
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-i", "--input", required = T,
                    help="Input RData file.")
parser$add_argument("-m", "--map", required = T,
                    help="Input MAP file.")
parser$add_argument("-r", "--rf", default = 0.38, type="double",
                    help="Recombination frequency threshold for grouping [default %(default)s].")
parser$add_argument("-o", "--output", required = T,
                    help="Output files name prefix.")
parser$add_argument("-1", "--group", action='store_false', 
                    help="Build linkage groups, otherwise all contigs into one group.")
parser$add_argument("-2", "--check", action='store_true', 
					help="Check and break linkage groups after ordering if distance greater than the predefined threshold.")			
parser$add_argument("-c", "--concorde", required = T,
                    help="Concorde executable file path.")
parser$add_argument("-l", "--include",
                    help="External R package locations. Multiple paths are seperated by ':'.")

args <- parser$parse_args()

in_RData = args$input
in_map = args$map
rf_thresh = args$rf
out_file = args$output
make_group = args$group
make_check = args$check
concorde_path = args$concorde
external_lib = args$include

if(!is.null(external_lib)) {
    libs = strsplit(external_lib,":")[[1]]
    .libPaths(c(.libPaths(),libs))
}

suppressPackageStartupMessages(library(TSP))
suppressPackageStartupMessages(library(igraph))

concorde_path(concorde_path)

initial.options <- commandArgs(trailingOnly = FALSE)
file.name <- "--file="
script.name <- sub(file.name, "", initial.options[grep(file.name, initial.options)])
source(paste(sep="/", dirname(script.name), "include.R"))

genetic_linkage_map(in_RData, in_map, out_file, rf_thresh, make_group, make_check)

