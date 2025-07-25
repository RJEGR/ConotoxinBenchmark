#!/bin/bash
#SBATCH --job-name=metaanalysis
#SBATCH -N 1
#SBATCH --mem=120GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

# Check for --help flag


while getopts "s:c:h" opt; do
    case $opt in
        s) MANIFEST_FILE="$OPTARG";;
        c) CONFIG_FILE="$OPTARG";;
        h)
            echo "Usage: $(basename $0) -s <Manifest> -c <Config>"
            echo
            echo "Arguments:"
            echo "  -s <Manifest>: A text file containing sample information. The first column should contain sample names, and the second column should contain file paths."
            echo "  -c <Config>:   A configuration file where each line specifies a tool and its command in the format 'tool=command'."
            echo
            echo "Description:"
            echo "  This script performs meta-analysis by concatenating reads from the Manifest file, configuring tools from the Config file, and executing the specified commands."
            exit 0
            ;;
        ?)
            echo "Invalid option: -$OPTARG" >&2
            echo "Run '$(basename $0) -h' for usage information."
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

Manifest=$MANIFEST_FILE 
CONFIG=$CONFIG_FILE 

if [[ -z "$Manifest" || -z "$CONFIG" ]]; then
    echo "Error: Missing required arguments. Both -s (Manifest) and -c (Config) are required."
    echo "Run '$(basename $0) -h' for usage information."
    exit 1
fi

single_reads=${Manifest%.txt}.concat.fq


call=`awk '{print $2}' $Manifest | tr "\n" " "`

cat $call > $single_reads


CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE

#./Run_trinity.sh $single_reads $CPU $MEM
#./Run_spades.sh $single_reads $CPU $MEM

chmod +x run_assemblers_dir/*sh

while read -r line; do 
    # Skip lines starting with '#'
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    [[ -z "${line// }" ]] && continue

    tool=$(echo $line | awk -F "=" '{print $1}')
    
    bs=${single_reads%_*}

    call=$(echo $line | awk -F "=" '{print $2}') 

    chkp_file="1_${tool}_${bs}.chkp"   
   
    if [ ! -f "$chkp_file" ]; then

        OUTDIR=${bs}_${tool}_dir
        
        mkdir -p $OUTDIR

        call=$(eval echo $call $single_reads $CPU $MEM $OUTDIR) # Expands variables in the command to fix issue with quotes         

        echo "Executing: $call"

        eval $call

        touch $chkp_file

    else
        echo "${chkp_file} checkpoint already exists."
    fi
done < $CONFIG


# End of script


exit
