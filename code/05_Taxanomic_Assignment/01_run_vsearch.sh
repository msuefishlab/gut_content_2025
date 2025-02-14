root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"

mkdir -p ${root}/output_data/slurm_logs/

for primer in "${primers[@]}"; do
	echo "Processing: $primer"
	echo sbatch --job-name ${primer}_vsearch --output ${root}"/output_data/slurm_logs/"${primer}"_vsearch.slurm.log" --export=root=${root},primer=${primer} ${root}/code/05_Taxanomic_Assignment/vsearch.sb
	sbatch --job-name ${primer}_vsearch --output ${root}"/output_data/slurm_logs/"${primer}"_vsearch.slurm.log" --export=root=${root},primer=${primer} ${root}/code/05_Taxanomic_Assignment/vsearch.sb

done