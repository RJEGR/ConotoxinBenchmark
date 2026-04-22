# Benchmarking the Accuracy of De Novo Assembly Conotoxins: A Transcriptomic Approach to Drug Discovery



![Alt text](https://github.com/RJEGR/ConotoxinBenchmark/blob/main/Figures/Figure1.png)
*Figure 1: The database schema describes the workflow used to generate a supplementary database after analyzing RNA-seq data.*

# Authors

Díaz-Castillo Brizuela-Rodríguez, Licea-Navarro and Gómez-Reyes (2026)

# Description

<div align="justify">
Motivation: Transcriptome data is set to play an increasingly important role in the development of venomics projects, which serves as a source to validate toxin genes to large-scale proteomics. Nevertheless, some limitation in the current proteotranscriptomic workflows relies on the ambiguous matches between both, transcriptome and mass spectrometry data. De Novo transcriptome assembly has been the technique used to obtain a catalog of candidate toxins. The accuracy of the technique to recover full-length toxins is influenced by several factors related to the empirical error from the sequencing platform and the biological composition of the venom RNA sample. While transcriptomic data ideally enhance the annotation of venom proteomes, the evaluation of the sensitivity and precision of assemblers for transcriptome is a key step toward establishing standards for bias identification and ensuring efficient workflow in conotoxin proteo-transcriptomic research.

Result: Here, we present a comparative study in which 12 open-source software  designed to assemble transcriptomes in a reference-free fashion were uniformly tested on a large set of simulated and real RNA-seq datasets spanning different groups of conotoxin members to different gene superfamilies and frameworks. We built > 100 assemblies and evaluated their performance on a combination of reference-free/reference-based precision, and sensitivity metrics. We demonstrate that assembly bias and performance in conotoxins are influenced by three factors: the depth of the sequence data, the sequence similarity between conotoxins (Low Entropy), and the choice or configuration of the assembly software, particularly the k-mer size parameter. The adjustment of software configuration may decrease chimerisms and collapse true family gene variations of conotoxins in the transcriptomic space. The reliance solely on transcriptomic analysis for the discovery of conotoxins results in a substantial underestimation of their actual diversity. 

Conclusion: The evaluation presented here serves both to assess the current status of the problem and to identify the most promising approaches to ensure further progress. We recommend integrating proteomic data and, preferably, acquiring reference genomes, which may significantly improve the matching of true toxins.

</div>

