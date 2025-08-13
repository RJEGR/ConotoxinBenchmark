# Creating ground truth dataset

Abinitio we started wrangling conoServer records from hosted [website](https://www.conoserver.org/download/conoserver_nucleic.xml.gz). You can figure out how to do it yourself using r (and r package xml, biostrings, rsample and tidyverse) to chose quality data and power and sensitivity analysis. We provided the example script `ConoServerDB.R` to help you. Then using it as ground truth dataset for the simulated data set using the following steps.

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

One of the null-hypothesis to test relly whether or not the assembler method have an effect in precision to assembly conotoxins. Thereby, this step evaluate the precision of the assembly method using the simulated data set. To do that, we used the script `Sampling_analysis.boostrap.slurm` to run the analysis using the simulated data set. This step requieres a configuration file where each line specifies an individual tool and its run-command in the format 'tool=command'. See the `run.config` file for an example. Also 
 
```bash
Assemblers.sh -s all -c run.config
#sbatch run_assembly_batches.sh
```



### Step 3
We also test whether or not sample complexity (as linear function of the number of reads) as effect in precision to assembly conotoxins. To do that, we used the script `Subsampling.sh` to run the analysis using individual simulated data set. For increment power and precision analysis we did this step 4-folds to replicate the analysis.

```bash
./Subsampling.sh -s $manifest -o boostrap_dir_${i}
# sbatch Sampling_analysis.boostrap.slurm

```

### Step 4

```bash
kmer_analysis.sh
```