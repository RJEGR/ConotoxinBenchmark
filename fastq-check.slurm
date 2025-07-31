EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/fastq-pair-master/build
export PATH=$PATH:$EXPORT


forward=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/art/200x/I3_superfamily//I3_superfamily_PE1.fq
reverse=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/art/200x/I3_superfamily//I3_superfamily_PE2.fq

fastq_pair $forward $reverse



find oases_auto/transrate_transcripts -name 'contigs.csv' -exec bash -c '
  for src; do
    dst="2_transrate_contigs_dir/${src#oases_auto/transrate_transcripts}"
    echo $dst
    mkdir -p "$(dirname "$dst")"
    cp "$src" "$dst"

  done
' bash {} +

du -sh