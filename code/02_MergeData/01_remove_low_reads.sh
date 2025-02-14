root="$(git rev-parse --show-toplevel)"
source "${root}/gut_contents.env"

mkdir -p "${root}/output_data/slurm_logs/"

for run in "${runs_to_filter[@]}"; do
    for primer in "${primers[@]}"; do
        if [ "$primer" == "primerset1" ]; then
            echo "Trimming: $run with $primer"
            echo sbatch --job-name "${run}_${primer}_remove_low_reads" --output "${root}/output_data/slurm_logs/${run}_${primer}_remove_low_reads.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",ps_lo_thresh="${ps1_lo_thresh}" "${root}/code/02_MergeData/remove_low_reads.sb"
            sbatch --job-name "${run}_${primer}_remove_low_reads" --output "${root}/output_data/slurm_logs/${run}_${primer}_remove_low_reads.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",ps_lo_thresh="${ps1_lo_thresh}" "${root}/code/02_MergeData/remove_low_reads.sb"

        else
            echo "Trimming: $run with $primer"
            echo sbatch --job-name "${run}_${primer}_remove_low_reads" --output "${root}/output_data/slurm_logs/${run}_${primer}_remove_low_reads.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",ps_lo_thresh="${ps2_lo_thresh}", "${root}/code/02_MergeData/remove_low_reads.sb"
            sbatch --job-name "${run}_${primer}_remove_low_reads" --output "${root}/output_data/slurm_logs/${run}_${primer}_remove_low_reads.slurm.log" --export=root="${root}",primerset="${primer}",run="${run}",ps_lo_thresh="${ps2_lo_thresh}", "${root}/code/02_MergeData/remove_low_reads.sb"
        fi
    done
done