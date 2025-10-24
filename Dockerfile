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

# Download facetsSuite repo
RUN wget -O facets-suite-${FACETSSUITE_VERSION}.zip https://github.com/mskcc/facets-suite/archive/${FACETSSUITE_VERSION}.zip \
    && unzip facets-suite-${FACETSSUITE_VERSION}.zip \
    && mv facets-suite-${FACETSSUITE_VERSION}/*-wrapper.R /usr/bin \
    && chmod +x /usr/bin/*wrapper.R

# Install package
RUN cd facets-suite-${FACETSSUITE_VERSION} \
    && Rscript -e "devtools::install()"

# Always run Rscript as vanilla
RUN export Rscript="Rscript --vanilla"
