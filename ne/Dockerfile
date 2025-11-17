FROM rocker/r2u

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential cmake \
  zlib1g pkg-config git zlib1g-dev 

RUN R -e "install.packages('terra')"

RUN R -e "install.packages('BiocManager')"

RUN R -e "BiocManager::install('SNPRelate')"

RUN R -e "install.packages('dartR.base')"

RUN R -e "install.packages('dartR.popgen')"

RUN R -e "install.packages('tidyverse')"

RUN R -e "install.packages('hierfstat')"


