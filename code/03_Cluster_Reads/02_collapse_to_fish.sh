root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

for primer in "${primers[@]}"; do
	echo "Processing: $primer"
	echo sbatch --job-name ${primer}_collapse --output ${root}"/output_data/slurm_logs/"${primer}"_collapse.slurm.log" --export=root=${root},primer=${primer} ${root}/code/03_Cluster_Reads/collapse.sb
	sbatch --job-name ${primer}_collapse --output ${root}"/output_data/slurm_logs/"${primer}"_collapse.slurm.log" --export=root=${root},primer=${primer} ${root}/code/03_Cluster_Reads/collapse.sb
done