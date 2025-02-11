root="$(git rev-parse --show-toplevel)"
source "${root}/gut_contents.env"

mkdir -p "${root}/output_data/slurm_logs/"

for run in "${runs[@]}"; do
    for primer in "${primers[@]}"; do
        if [ "$primer" == "primerset1" ]; then
            echo "Trimming: $run with $primer"
            echo sbatch --job-name "${run}_${primer}_trim_primers" --output "${root}/output_data/slurm_logs/${run}_${primer}_trim_primers.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",fprimer="${ps1_primers[0]}",rprimer="${ps1_primers[1]}" "${root}/code/trim_primers.sb"
            sbatch --job-name "${run}_${primer}_trim_primers" --output "${root}/output_data/slurm_logs/${run}_${primer}_trim_primers.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",fprimer="${ps1_primers[0]}",rprimer="${ps1_primers[1]}" "${root}/code/trim_primers.sb"
        else
            echo "Trimming: $run with $primer"
            echo sbatch --job-name "${run}_${primer}_trim_primers" --output "${root}/output_data/slurm_logs/${run}_${primer}_trim_primers.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",fprimer="${ps2_primers[0]}",rprimer="${ps2_primers[1]}" "${root}/code/trim_primers.sb"
            sbatch --job-name "${run}_${primer}_trim_primers" --output "${root}/output_data/slurm_logs/${run}_${primer}_trim_primers.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",fprimer="${ps2_primers[0]}",rprimer="${ps2_primers[1]}" "${root}/code/trim_primers.sb"
        fi
    done
done