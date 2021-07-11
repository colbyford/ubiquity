FROM mcr.microsoft.com/azureml/base:openmpi3.1.2-ubuntu18.04

RUN conda install -c r -y pip=20.1.1 openssl=1.1.1c r-base rpy2 && \
    conda install -c conda-forge -y mscorefonts && \
    conda clean -ay && \
    pip install --no-cache-dir azureml-defaults
    
ENV TAR="/bin/tar"
RUN R -e "install.packages(c('remotes', 'reticulate', 'optparse', 'here', 'azuremlsdk'), repos = 'https://cloud.r-project.org/')"
RUN R -e "remotes::install_github('https://github.com/colbyford/ubiquity')"


## Install Optional Dependencies
RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev

RUN R -e "install.packages(c('survival', 'tidyverse'), repos = 'https://cloud.r-project.org/')"