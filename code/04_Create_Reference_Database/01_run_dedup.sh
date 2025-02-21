root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

echo sbatch --job-name BOLD_dedup --output ${root}"/output_data/slurm_logs/BOLD_dedup.slurm.log" --export=root=${root} ${root}/code/04_Create_Reference_Database/dedup.sb
sbatch --job-name BOLD_dedup --output ${root}"/output_data/slurm_logs/BOLD_dedup.slurm.log" --export=root=${root} ${root}/code/04_Create_Reference_Database/dedup.sb
