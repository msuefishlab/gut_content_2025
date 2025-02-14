root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

for primer in "${primers[@]}"; do
	echo "Processing: $primer"
	echo sbatch --job-name ${primer}_cluster --output ${root}"/output_data/slurm_logs/"${primer}"_cluster.slurm.log" --export=root=${root},primer=${primer} ${root}/code/03_Cluster_Reads/cluster.sb
	sbatch --job-name ${primer}_cluster --output ${root}"/output_data/slurm_logs/"${primer}"_cluster.slurm.log" --export=root=${root},primer=${primer} ${root}/code/03_Cluster_Reads/cluster.sb
done