# Creating ground truth dataset

Abinitio we started wrangling conoServer records from hosted [website](https://www.conoserver.org/download/conoserver_nucleic.xml.gz). You can figure out how to do it yourself using r (and r package xml, biostrings, rsample and tidyverse) to chose quality data and power and sensitivity analysis. We provided the example script `ConoServerDB.R` to help you. Then using it as ground truth dataset for the simulated data set using the following steps.

![alt text](https://github.com/RJEGR/ConotoxinBenchmark/blob/main/Figures/Analysis_schema.png)

A flowchart illustrating the analysis schema for creating a ground truth dataset and simulating RNAseq data. The diagram shows sequential steps including data wrangling from the ConoServer database, running simulation scripts, generating depth files, evaluating assembly precision, and performing subsampling analysis. Each step is visually connected with arrows, and the chart includes labels such as Step 1 Simulate RNAseq, Step 2 Assembly Evaluation, Step 3 Subsampling Analysis, and Step 4 Kmer Analysis. The environment is a clean digital workspace with a logical progression of tasks. The tone is informative and methodical.


### Step 1.
Use the pipeline script `Simulate_rnaseq.sh` to create a simulated data set. With flag `-d` Define the folder where reference fasta file used to simulated data are allocated (`vfolds_resampling_dir`)

```bash
Simulate_rnaseq.sh -d vfolds_resampling_dir
# sbatch simulate_rnaseq.sh -d vfolds_resampling_dir
```

This step will produce a folder with the simulated data set, allocated in subfolder of the prefix `${prefix%.fasta}`, and Manifest files in txt extension. Next time you decide use any simulated PE100 RNAseq data you need to use Manifest file to call the files, so please do not delete or move the structure of the folders where the simulated datasets were allocated.

### Step 1.2

Additonal output generated from the step 1 are depth files where coverage for every simulated PE100 RNAseq data. Use Effort.R example script to summarise depth files and evaluate the sampling effort

### Step 2.

One of the null-hypothesis to test relly whether or not the assembler method have an effect in precision to assembly conotoxins. Thereby, this step evaluate the precision of the assembly method using the simulated data set. To do that, we used the script `Assemblers.sh` to run the analysis using the simulated data set. This step requieres a configuration file where each line specifies an individual tool and its run-command in the format 'tool=command'. See the `run.config` file for an example. Also 
 
```bash
./Assemblers.sh -s all -c run.config
#sbatch run_assembly_batches.sh
```



### Step 3
We also test whether or not sample complexity (as linear function of the number of reads) as effect in precision to assembly conotoxins. To do that, we used the script `Subsampling.sh` to run the analysis using individual simulated (-s all) data set and different subset of proportions (-p).

```bash
for P in $(seq 0.1 0.1 1); do
     ./Subsampling.sh -s all -p $P
done

```

#### Step 3.2

Additonal output generated from the step 3 are &&&. Use Effort.R example script to summarise depth files and evaluate the sampling effort



### Step 4

```bash
kmer_analysis.sh
```

### Step 5
Conotoxin classification
Select baseline dataset
Evals ConoSorter, blast, hhm, hhsuite,

### Step 6
Using real-data 
Downloading step
Assembly (spades, trinity, plass)
Evals metrics (busco, Nx)
Annotate


```bash

find . -name "*_1.fastq" | while read -r fastq_file; do

output_bs=$(basename "$fastq_file" _1.fastq)

relative_path="$PWD/${fasta_file##*/}" # Extract the filename from the path
absolute_path="$(cd "$(dirname "$relative_path")" && pwd)/$(basename "$relative_path")"

echo "Absolute path: $absolute_path"

forward_file=$absolute_path/${output_bs}_1.fastq
reverse_file=$absolute_path/${output_bs}_2.fastq

echo "Running artificial step for $forward_file and $reverse_file."

echo -e "$output_bs\t$forward_file\t$reverse_file\t$reference_file" >> "${output_bs}_samples.txt"

done


sample_count=$(find . -maxdepth 1 -name "*_samples.txt" | wc -l)

if [ "$sample_count" -gt 5 ]; then
     find . -maxdepth 1 -name "*_samples.txt" > batches.log
     split -l 5 batches.log batch_

     for batch in batch_*; do
          subfolder="${batch}_dir"
          mkdir -p "$subfolder"
          while read -r sample_file; do
               cp -r Assemblers.sh ./run_assemblers_dir run.suitable.config "$sample_file" "$subfolder/"
          done < "$batch"
     done

     rm batches.log batch_*
fi

```

line_count=$(find . -maxdepth 1 -name "*_samples.txt" | wc -l)

if [ "$line_count" -gt 5 ]; then
     find . -maxdepth 1 -name "*_samples.txt" | split -l 5 - batch_

     for batch in batch_*; do
          subfolder="${batch%.*}_dir"
          mkdir -p "$subfolder"
          while read -r line; do
               echo "$line" > "$subfolder/$(echo "$line" | cut -f1).txt"
          done < "$batch"
     done
fi

```

```bash

run_assembler_batches() {
     local txt_count

     txt_count=$(find . -maxdepth 1 -name "*_samples.txt" | wc -l)
     
     if [ "$txt_count" -le 5 ]; then
          
          echo "sbatch ./Assemblers.sh -s all -c run.suitable.config"
     
     else
          # Split txt files into batches of 5 and run Assemblers.sh for each batch
          find . -maxdepth 1 -name "*_samples.txt" | split -l 5 - batch_
          
          for batch in batch_*; do
               files=$(cat "$batch" | tr '\n' ',' | sed 's/,$//')

               echo "sbatch ./Assemblers.sh -s "$files" -c run.suitable.config"
          done

          rm batch_*
     fi
}

run_assembly_batches.sh

```
