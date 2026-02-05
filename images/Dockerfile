FROM rocker/verse:4

USER root
RUN chown -R rstudio:rstudio /home/rstudio/

# Install taxize and devtools (and their dependencies)
RUN install2.r --error --deps TRUE \
	taxize \
	devtools \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# Install bold from the archived tarball URL
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/bold/bold_1.3.0.tar.gz', repos = NULL, type = 'source')"

# Bioconductor packages (use Rscript for logic)
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('Biostrings', ask=FALSE, update=FALSE)"


RUN apt-get clean && apt-get -y update && apt-get install -y --no-install-recommends libxt6 imagemagick
