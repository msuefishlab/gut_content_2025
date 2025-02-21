#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"gut_contents.env"


IMAGE=${IMAGE:-jasongallant/r4guts:latest}

PORT=${RSTUDIO_PORT:-8787}
export PASSWORD=$(openssl rand -base64 15)

HPC_ENV="local"


function get_port {
    # lsof doesn't return open ports for system services, so we use netstat
    # until ! lsof -i -P -n | grep -qc ':'${PORT}' (LISTEN)';
    
    until ! netstat -ln | grep "  LISTEN  " | grep -iEo  ":[0-9]+" | cut -d: -f2 | grep -wqc ${PORT};
    do
        ((PORT++))
        echo "Checking port: ${PORT}"
    done
    echo "Got one !"
}

IMAGE_SLASHED=$(echo $IMAGE | sed 's/:/\//g')
RSTUDIO_HOME=${root}/images/.rstudio-rocker/${IMAGE_SLASHED}/session
RSTUDIO_TMP=${root}/images/.rstudio-rocker/${IMAGE_SLASHED}/tmp
# RSITELIB=${HOME}/.rstudio-rocker/${IMAGE_SLASHED}/site-library
R_LIBS_USER=${root}/images/.rstudio-rocker/${IMAGE_SLASHED}
#mkdir -p ${HOME}/.rstudio
mkdir -p ${RSTUDIO_HOME}
#mkdir -p ${RSITELIB}
mkdir -p ${R_LIBS_USER}
mkdir -p ${RSTUDIO_TMP}
mkdir -p ${RSTUDIO_TMP}/var/run


echo
echo "Finding an available port ..."
get_port

LOCALPORT=${PORT}
# LOCALPORT=8787
PUBLIC_IP=$(curl https://checkip.amazonaws.com)

echo "On you local machine, open an SSH tunnel like:"
# echo "  ssh -N -L ${LOCALPORT}:localhost:${PORT} ${USER}@m3-bio1.erc.monash.edu.au"
echo "  ssh -N -L ${LOCALPORT}:localhost:${PORT} ${USER}@$(hostname -f)"
echo "  or"
echo "  ssh -N -L ${LOCALPORT}:localhost:${PORT} ${USER}@${PUBLIC_IP}"

# For smux/srun/sbatch jobs, route via the login node to a the compute node where rserver runs - not working for me
# echo "  ssh -N -L ${LOCALPORT}:${HOSTNAME}:${PORT} ${USER}@m3.massive.org.au"
echo
echo "Point your web browser at http://localhost:${LOCALPORT}"
echo
echo "Login to RStudio with:"
echo "  username: ${USER}"
echo "  password: ${PASSWORD}"
echo
echo "Protip: You can choose your version of R from any of the tags listed here: https://hub.docker.com/r/rocker/rstudio/tags"
echo "        and set the environment variable IMAGE, eg"
echo "        IMAGE=rocker/rstudio:4.1.1 $(basename "$0")"
echo
echo "Starting RStudio Server (R version from image ${IMAGE})"

# Set some locales to suppress warnings
LC_CTYPE="C"
LC_TIME="C"
LC_MONETARY="C"
LC_PAPER="C"
LC_MEASUREMENT="C"

SINGULARITYENV_PASSWORD="${PASSWORD}" \
singularity exec --bind ${HOME}:/home/rstudio \
					--bind ${RSTUDIO_HOME}:${HOME}/.rstudio \
					--bind ${R_LIBS_USER}:${R_LIBS_USER} \
					--bind ${RSTUDIO_TMP}:/tmp \
					--bind=${RSTUDIO_TMP}/var:/var/lib/rstudio-server \
					--bind=${RSTUDIO_TMP}/var/run:/var/run/rstudio-server \
					--env R_LIBS_USER=${R_LIBS_USER} \
					${rimage} \
					rserver --auth-none=0 --auth-pam-helper-path=pam-helper --www-port=${PORT} --server-user $(whoami)
					# --bind ${RSITELIB}:/usr/local/lib/R/site-library \

printf 'rserver exited' 1>&2