#!/bin/bash

# envirnment variables
export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

which transrate

asm=$1
ref=$2
outdir=$3
threads=$4




call="transrate --assembly $asm --reference $ref --output 2_transrate_dir/$outdir --threads $threads"

echo "Executing: $call"

eval $call