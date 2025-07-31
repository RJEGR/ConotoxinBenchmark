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


mkdir -p transrate_tmp_dir

call="transrate --left $forward_fq --right $reverse_fq --assembly $asm --reference $ref --output transrate_tmp_dir/$outdir --threads 20"


echo "Executing: $call"

eval $call

find transrate_tmp_dir -name 'contigs.csv' -exec bash -c '
  for src; do
    dst="2_transrate_contigs_dir/${src#transrate_tmp_dir}"
    echo $dst
    mkdir -p "$(dirname "$dst")"
    cp "$src" "$dst"

  done
' bash {} +

exit
