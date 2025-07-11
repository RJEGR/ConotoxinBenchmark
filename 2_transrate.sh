#!/bin/bash

# envirnment variables
export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

which transrate


forward_fq=$1
reverse_fq=$2
asm=$3
ref=$4
outdir=$5


call="transrate --left $forward_fq --right $reverse_fq --assembly $asm --reference $ref --output 2_transrate_dir/$outdir --threads 20"

# call="transrate --assembly $asm --reference $ref --output 2_transrate_dir/$outdir --threads 20"

echo "Executing: $call"

eval $call

exit
