# Use the official Miniconda3 image as the base
FROM continuumio/miniconda3:latest

# Create the conda environment "mkcoinr" with Python 3.9 and install packages
RUN conda create -n mkcoinr python=3.9 -y && \
	conda install -n mkcoinr -c bioconda blast -y && \
	conda install -n mkcoinr -c bioconda vsearch -y && \
	/opt/conda/envs/mkcoinr/bin/python -m pip install cutadapt nsdpy

# Set the PATH so that executables in the mkcoinr environment are used by default
ENV PATH /opt/conda/envs/mkcoinr/bin:$PATH

# Default command (adjust as needed)
CMD ["python"]
