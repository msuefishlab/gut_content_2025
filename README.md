	cd images/
    singularity build --remote qiime2.sif qiime2_mod.def 


	singularity build qiime2019.sif docker://quay.io/qiime2/core:2019.1


	For R computations:

		    docker build --platform linux/amd64 -t r4guts .
    		docker tag gwas_tools jasongallant/r4guts
    		docker push jasongallant/r4guts