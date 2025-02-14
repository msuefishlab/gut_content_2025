FROM rocker/verse:4

USER root
RUN chown -R rstudio:rstudio /home/rstudio/

RUN install2.r --error --deps TRUE \
	bold \
	taxize \
	devtools \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN apt-get clean && apt-get -y update && apt-get install -y --no-install-recommends libxt6 imagemagick