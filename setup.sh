#!/bin/bash
# ----------------------------------------------------------------------------------------------------
# Title:       setup.sh
# Description: Setup conda environment for Facets and install all the python and R libraries.
# Dependency:  conda>=4.7; environment.yaml (in the same directory as this script)
# Usage:       ./setup.sh <optional_environment_name>
# ----------------------------------------------------------------------------------------------------

#############
# Functions #
#############

HOSTNAME=$(hostname -s)

function printe () {
	# print error message
        printf "\n$HOSTNAME::$(date +'%T')::ERROR $@\n"
}

function printi () {
	# print info
        printf "\n$HOSTNAME::$(date +'%T')::INFO $@\n"
}

function cleanup() {
	# Housekeeping function that will run at the end regardless of 
	# exit status. Upon successful setup, deactivate environment.
	# If setup is unsuccessful at any stage, remove the environment.
        EXITCODE=$?
        if [[ $EXITCODE != 0 ]] && [[ -e ${FACETS_ENV_PATH} ]]; then
		source deactivate
                printe "Installation unsuccessful. Removing conda environment: ${FACETS_ENV}"
                conda env remove -n ${FACETS_ENV}
	fi
        exit $EXITCODE
}
trap cleanup EXIT


###########################
# Setup Conda Environment #
###########################

# Check for environment name
FACETS_ENV=$1
[[ ! -z $FACETS_ENV ]] || {
        printi "No preferred environment name given. Conda environment name will be set to 'FACETS'";
        FACETS_ENV="FACETS";
        }

# set library paths
[[ ! -z "${LIBRARY_PATH}" ]] || {
        printi "LIBRARY_PATH not set. Will be set to /usr/lib64";
        #export LIBRARY_PATH="/usr/lib64";
        }

[[ ! -z "${LD_LIBRARY_PATH}" ]] || {
        printi "LD_LIBRARY_PATH not set. Will be set to /usr/lib64";
        #export LD_LIRBARY_PATH="/usr/lib64";
        }

# Check for conda installation
type -p conda > /dev/null || {
        printe "Conda installation is required! Visit https://www.anaconda.com/distribution/";
        exit 1;
        }

# Set pythonpath
export PYTHONPATH="" # to avoid conflicts with system python libraries


[[ -e "${PWD}/facets_environment.yaml" ]] || {
        printe "${PWD}/facets_environment.yaml not found."; 
        exit 1;
        }

# get conda path
CONDA=$(type -p conda | sed "s/conda is //")
[[ -e $CONDA ]] || {
        printe "Cannot find conda binary. Make sure conda is added to your PATH variable.";
        echo $CONDA
        exit 1;
        }

# ensure that no other env with the same name exist
FACETS_ENV_PATH=$(echo $CONDA | sed "s/bin\/conda/envs\/${FACETS_ENV}/")
[[ ! -e ${FACETS_ENV_PATH} ]] || {
        printi "${FACETS_ENV_PATH} already exist.";
        #exit 0;
        }

printi "Creating conda environment: $FACETS_ENV"
conda env create --name $FACETS_ENV --file=$PWD/facets_environment.yaml

EXITCODE=$?
[[ $EXITCODE == 0 ]] || exit $EXITCODE;

# TODO:
# Use conda pre-built R libraries. Since some of the libraries are in dev, lets install it using R for now.
# Activate environment
printi "Activating ${FACETS_ENV}"
source activate ${FACETS_ENV}
[[ $EXITCODE == 0 ]] || {
        printe "Cannot activate ${FACETS_ENV}.";
        exit $EXITCODE;
        }

RSCRIPT=$(echo $(which python) | sed "s/python/Rscript/")

[[ -e "$RSCRIPT" ]] || {
        printe "Cannot locate Rscript in the environment $FACETS_ENV. Was your environment setup properly?";
        exit 1;
        }

printi "Installing R packages using devtools..."

$RSCRIPT -e 'devtools::install_github("rptashkin/facets2n", ref= "master", force=T); devtools::install_github("taylor-lab/facets-suite", ref = "feature/facets2N");'

EXITCODE=$?
[[ $EXITCODE == 0 ]] || {
        printe "Error during R package installation from github.";
        exit $EXITCODE;
        }

printi "Environment setup complete! \^.^/"
conda deactivate
