root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

echo sbatch --job-name BOLD_create_consensus --output ${root}"/output_data/slurm_logs/BOLD_create_consensus.slurm.log" --export=root=${root} ${root}/code/04_Create_Reference_Database/create_consensus.sb
sbatch --job-name BOLD_create_consensus --output ${root}"/output_data/slurm_logs/BOLD_create_consensus.slurm.log" --export=root=${root} ${root}/code/04_Create_Reference_Database/create_consensus.sb
