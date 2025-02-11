root="$(git rev-parse --show-toplevel)"
source "${root}/gut_contents.env"

mkdir -p "${root}/output_data/slurm_logs/"

for run in "${runs[@]}"; do
    for primer in "${primers[@]}"; do
        if [ "$primer" == "primerset1" ]; then
            echo "Trimming: $run with $primer"
            echo sbatch --job-name "${run}_${primer}_denoise" --output "${root}/output_data/slurm_logs/${run}_${primer}_denoise.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",ftrunc="${ps1_trunc[0]}",rtrunc="${ps1_trunc[1]}" "${root}/code/denoise.sb"
            sbatch --job-name "${run}_${primer}_denoise" --output "${root}/output_data/slurm_logs/${run}_${primer}_denoise.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",ftrunc="${ps1_trunc[0]}",rtrunc="${ps1_trunc[1]}" "${root}/code/denoise.sb"

        else
            echo "Trimming: $run with $primer"
            echo sbatch --job-name "${run}_${primer}_denoise" --output "${root}/output_data/slurm_logs/${run}_${primer}_denoise.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",ftrunc="${ps2_trunc[0]}",rtrunc="${ps2_trunc[1]}" "${root}/code/denoise.sb"
            sbatch --job-name "${run}_${primer}_denoise" --output "${root}/output_data/slurm_logs/${run}_${primer}_denoise.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",ftrunc="${ps2_trunc[0]}",rtrunc="${ps2_trunc[1]}" "${root}/code/denoise.sb"

        fi
    done
done