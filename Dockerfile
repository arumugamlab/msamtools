################## BASE IMAGE ######################
FROM alpine:3.23

################## METADATA ######################
LABEL base_image="alpine:3.23"
LABEL version="1"
LABEL software="msamtools"
LABEL software.version="1.1.3"
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

ENV MSAM_VERSION 1.1.3

RUN apk --no-cache update \
    && apk --no-cache upgrade \
    && apk add --no-cache argtable2 \
    && apk add --no-cache --virtual .build_deps gcc libc-dev wget zlib-dev make bash argtable2-dev \
    && cd /tmp \
    && wget https://github.com/arumugamlab/msamtools/releases/download/$MSAM_VERSION/msamtools-$MSAM_VERSION.tar.gz -O - | tar xfz - \
    && cd msamtools-$MSAM_VERSION/ \
    && ./configure --prefix=/usr \
    && make install \
    && /usr/bin/install -c deps/samtools/samtools-1.9/samtools /usr/bin \
    && cd /tmp \
    && rm -rf msamtools-$MSAM_VERSION \
    && apk del .build_deps

ENTRYPOINT ["msamtools"]
