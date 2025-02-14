root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

for run in "${runs[@]}"; do
    for primer in "${primers[@]}"; do
        echo "Processing: $run with $primer"
        echo sbatch --job-name ${run}_${primer}_import_qiime --output ${root}"/output_data/slurm_logs/"${run}"_"${primer}"_import.slurm.log" --export=root=${root},primerset=${primer},run=${run} ${root}/code/01_Filter_And_Denoise/load_qiime_data.sb
        sbatch --job-name ${run}_${primer}_import_qiime --output ${root}"/output_data/slurm_logs/"${run}"_"${primer}"_import.slurm.log" --export=root=${root},primerset=${primer},run=${run} ${root}/code/01_Filter_And_Denoise/load_qiime_data.sb
    done
done