# convert fast5 to pod5, recall bases with dorado

pip3 install pod5 --break-system-packages

fast5 ~/Desktop/fast5files/*.fast5 --output output_pod5s/ --one-to-one ~/Desktop/fast5files/

samtools fasta ~/Desktop/recalled_sup.bam > ~/Desktop/recalled_sup.fa