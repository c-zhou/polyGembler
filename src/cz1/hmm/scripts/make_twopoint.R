
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-i", "--input", required = T, 
                    help="Input RData file.")
parser$add_argument("-2", "--twopoint", required = T, 
					help="Input two-point cluster file.")
parser$add_argument("-r", "--rf", type="double",
                    help="Fixed recombination frequency [default read from input file].")
parser$add_argument("-o", "--output", required = T, 
                    help="Output file name.")
parser$add_argument("-c", "--concorde", required = T,
                    help="Concorde executable file path.")
parser$add_argument("-l", "--include",
                    help="External R package locations. Multiple paths are seperated by ':'.")

args <- parser$parse_args()

in_RData = args$input
twopoint_file = args$twopoint
fixed_rf = args$rf
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

if(is.null(fixed_rf)) {
	two_point(in_RData, twopoint_file, out_file)
} else {
	two_point(in_RData, twopoint_file, out_file, fixed_rf)
} 

