#!/bin/sh

export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

which transrate

run_transrate() {
    
    REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs

    mkdir -p transrate_tmp_dir
    
    find */*_FASTA_DIR -maxdepth 1 -type f -name "*.fa" | while read -r asm; do
        
        echo "Processing FASTA : $asm"

        # Substring to extract the sf prefix

        FAMILYPREFIX=$(echo $asm | awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /_FASTA_DIR/) print $i}')


        echo ${FAMILYPREFIX%_FASTA_DIR}

        ref=$REFDIR/${FAMILYPREFIX%_FASTA_DIR}.fasta
        
        
        TRANSRATE_DIR=${asm%.fa}_dir

        call="transrate --assembly $asm --reference $ref --output transrate_tmp_dir/$TRANSRATE_DIR --threads 20"

        echo "Executing: $call"

        eval $call

        echo "Finished processing FASTA : $asm"
    
    done



}

run_transrate


mkdir -p transrate_contigs_dir

find transrate_tmp_dir -name 'contigs.csv' -exec bash -c '
        for src; do
              dst="transrate_contigs_dir/${src#transrate_tmp_dir}"
              echo $dst
              mkdir -p "$(dirname "$dst")"
              cp "$src" "$dst"

        done
' bash {} +

# rm -rf transrate_tmp_dir

exit