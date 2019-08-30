
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-i", "--input", required = T, 
                    help="Input RData file.")
parser$add_argument("-n", "--nn", type="double", default=2,
                    help="Number of nearest neighbours [default %(default)s].")
parser$add_argument("-d", "--distance", default = 50, type="double",
                    help="Genetic distance (centimorgan) threshold for grouping [default %(default)s].")
parser$add_argument("-o", "--output", required = T, 
                    help="Output file name.")
parser$add_argument("-c", "--concorde", required = T,
                    help="Concorde executable file path.")
parser$add_argument("-l", "--include",
                    help="External R package locations. Multiple paths are seperated by ':'.")

args <- parser$parse_args()

in_RData = args$input
nn = args$nn
max_d=args$distance
out_file = args$output
concorde_path = args$concorde
external_lib = args$include

if(!is.null(external_lib)) {
    libs = strsplit(external_lib,":")[[1]]
    .libPaths(c(.libPaths(),libs))
}

suppressPackageStartupMessages(library(caTools))
suppressPackageStartupMessages(library(TSP))

concorde_path(concorde_path)

initial.options <- commandArgs(trailingOnly = FALSE)
file.name <- "--file="
script.name <- sub(file.name, "", initial.options[grep(file.name, initial.options)])
source(paste(sep="/", dirname(script.name), "include.R"))

nearest_neighbour_joining(in_RData, out_file, nn-1, max_d)

