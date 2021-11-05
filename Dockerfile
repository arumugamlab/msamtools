################## BASE IMAGE ######################
FROM alpine:3.14.2

################## METADATA ######################
LABEL base_image="alpine:3.14"
LABEL version="1"
LABEL software="msamtools"
LABEL software.version="1.0.1"
LABEL about.summary="microbiome-related extension to samtools"
LABEL about.home="https://github.com/arumugamlab/msamtools"
LABEL about.documentation="https://github.com/arumugamlab/msamtools"
LABEL about.license_file="https://github.com/arumugamlab/msamtools"
LABEL about.license="SPDX:MIT"
LABEL about.tags="Metagenomics, Genomics"
LABEL extra.identifiers.biotools="msamtools"
LABEL org.opencontainers.image.authors="arumugam@sund.ku.dk"

################## MAINTAINER ######################
MAINTAINER Mani Arumugam <arumugam@sund.ku.dk>

################## INSTALLATION ######################

ENV MSAM_VERSION 1.0.1

RUN apk --no-cache update \
    && apk --no-cache upgrade \
    && apk add --no-cache gcc libc-dev wget zlib-dev make bash \
    && cd /tmp \
    && wget http://arumugamlab.sund.ku.dk/Tutorials/msamtools-$MSAM_VERSION.tar.gz -O - | tar xfz - \
    && cd msamtools-$MSAM_VERSION/ \
    && ./configure --prefix=/usr \
    && make install \
    && /usr/bin/install -c deps/samtools/samtools-1.9/samtools /usr/bin \
    && rm -rf /tmp/msamtools-$MSAM_VERSION \
    && apk del gcc libc-dev wget zlib-dev make

ENTRYPOINT ["msamtools"]
