FROM rocker/tidyverse:3.6.1

LABEL maintainer="Philip Jonsson <philip.jonsson@gmail.com>"
ENV HTSLIB_VERSION 1.5
ENV HTSTOOLS_VERSION 0.1.1
ENV FACETSSUITE_VERSION Rpackagev2
ENV SNP_PILEUP /usr/bin/snp-pileup

# Install necessary stuff
RUN apt-get update && apt-get install -y g++ tar bzip2 libbz2-dev liblzma-dev

# Download and install htslib
RUN cd /tmp \
    && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar xvjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure \
    && make \
    && make install \
    && cp libhts.so* /usr/lib

# Download snp-pileup repo, install snp-pileup executable
RUN cd /tmp \
    && wget https://github.com/mskcc/htstools/archive/snp_pileup_${HTSTOOLS_VERSION}.tar.gz \
    && tar xzvf snp_pileup_${HTSTOOLS_VERSION}.tar.gz \
    && cd htstools-snp_pileup_${HTSTOOLS_VERSION} \
    && g++ -std=c++11 snp-pileup.cpp -lhts -o snp-pileup \
    && mv snp-pileup /usr/bin

# Clean up tmp
RUN rm -rf /tmp/*

# Add Facets Suite to the container
RUN mkdir /facets-suite
ADD . /facets-suite
RUN cd /facets-suite && \
    Rscript -e "devtools::install()"
ENV PATH=/facets-suite:$PATH

# install Facets
RUN wget https://github.com/mskcc/facets/archive/v0.5.14.zip -O facets_v0.5.14.zip && \
    unzip facets_v0.5.14.zip && \
    rm -f facets_v0.5.14.zip && \
    cd facets-0.5.14 && \
    Rscript -e "devtools::install()"
