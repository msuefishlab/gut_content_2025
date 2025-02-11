root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

for run in "${runs[@]}"; do
    for primer in "${primers[@]}"; do
        echo "Processing: $run with $primer"
        echo sbatch --job-name ${run}_${primer}_inspect_denoise --output ${root}"/output_data/slurm_logs/"${run}"_"${primer}"_inspect_denoise.slurm.log" --export=root=${root},primerset=${primer},run=${run} ${root}/code/inspect_denoise.sb
        sbatch --job-name ${run}_${primer}_inspect_denoise --output ${root}"/output_data/slurm_logs/"${run}"_"${primer}"_inspect_denoise.slurm.log" --export=root=${root},primerset=${primer},run=${run} ${root}/code/inspect_denoise.sb

    done
done