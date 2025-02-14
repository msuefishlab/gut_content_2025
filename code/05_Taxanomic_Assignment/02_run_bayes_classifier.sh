root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

for primer in "${primers[@]}"; do
	echo "Processing: $primer"
	echo sbatch --job-name ${primer}_nbClassifier --output ${root}"/output_data/slurm_logs/"${primer}"_nbClassifier.slurm.log" --export=root=${root},primer=${primer} ${root}/code/05_Taxanomic_Assignment/nbClassifier.sb
	sbatch --job-name ${primer}_nbClassifier --output ${root}"/output_data/slurm_logs/"${primer}"_nbClassifier.slurm.log" --export=root=${root},primer=${primer} ${root}/code/05_Taxanomic_Assignment/nbClassifier.sb
done