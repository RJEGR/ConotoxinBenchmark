#!/bin/bash

#SBATCH --job-name=transrate
#SBATCH --output=transrate-%j.log
#SBATCH --error=transrate-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=6-00:00:00
#SBATCH -p cicese

# 1) Configuración del entorno
export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

# 2) Parámetros
threads=24
REF_DIR="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/inputs"
ASM_DIR="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/art/ensambles"
TRANSRATE_ROOT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/fernando_pub/reads_artificiales/art/transrate"

mkdir -p "$TRANSRATE_ROOT"

# 3) Para cada referencia…
find "$REF_DIR" -maxdepth 1 -type f -name '*.fasta' | while read -r ref; do
  sample=$(basename "${ref%.*}")

  # 4) Busca sus ensamblajes correspondientes
  for asm in "$ASM_DIR"/*_"${sample}"_trinity_ensamble.Trinity.fasta; do
    # comprueba que exista
    [[ -f "$asm" ]] || continue

    # extrae la cobertura (el primer campo antes de '_')
    coverage=$(basename "$asm" | cut -d_ -f1)

    # 5) Crea un directorio de salida para este análisis
    outdir="$TRANSRATE_ROOT/${coverage}_${sample}_transrate"
    mkdir -p "$outdir"

    echo "=== Transrate: ${coverage}_${sample} ==="

    # 6) Ejecuta Transrate
    transrate \
      --assembly "$asm" \
      --reference "$ref" \
      --threads "$threads" \
      --output "$outdir"
  done
done

echo "✔ Todos los análisis de Transrate han finalizado."
