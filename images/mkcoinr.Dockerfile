# Use the official Miniconda3 image as the base
FROM continuumio/miniconda3:latest

# Create the conda environment "mkcoinr" with Python 3.9 and install packages
RUN conda create -n mkcoinr python=3.9 -y && \
	conda install -n mkcoinr -c bioconda blast -y && \
	conda install -n mkcoinr -c bioconda vsearch -y && \
	/opt/conda/envs/mkcoinr/bin/python -m pip install cutadapt nsdpy

# Install Git
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*

# Clone the repository into /opt/mkCOInr
RUN git clone https://github.com/meglecz/mkCOInr.git /opt/mkCOInr

# Add the conda environment's bin directory and the repository's scripts directory to PATH
ENV PATH /opt/mkCOInr/scripts:/opt/conda/envs/mkcoinr/bin:$PATH

# Default command (adjust as needed)
CMD ["python"]
