root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

echo sbatch --job-name COInr_Train --output ${root}"/output_data/slurm_logs/COInr_Train.slurm.log" --export=root=${root} ${root}/code/04_Reference_Database/train.sb
sbatch --job-name COInr_Train --output ${root}"/output_data/slurm_logs/COInr_Train.slurm.log" --export=root=${root} ${root}/code/04_Reference_Database/train.sb
