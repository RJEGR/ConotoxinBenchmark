
#!/bin/bash
#SBATCH --job-name=transrate_coverage
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -o out.%J
#SBATCH -e err.%J

mkdir -p chkp_dir

for i in $(ls *fa); do
    
    echo "Processing data: $i"
    
    manifest=${i%_*}.txt

    chkp_file="${i%.fa}.chkp" 

    if [ ! -f chkp_dir/"$chkp_file" ]; then

    call="./Accuracy.sh -s $manifest -d $i -m reference"
    
    echo $call
    
    eval $call
        
    touch chkp_dir/"$chkp_file"

    else
        echo "${chkp_file} checkpoint already exists in chkp_dir."          
    
    continue
    fi

done

exit