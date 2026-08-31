################## BUILDING NOTES ######################
#
# Build from the published GitHub release (default):
#   docker build -t msamtools:1.2.0 .
#
# Build from a local release tarball for testing:
#   Place msamtools-<version>.tar.gz in the build context, then run:
#   docker build --build-arg MSAM_SOURCE=local --build-arg MSAM_VERSION=<version> -t msamtools:test .
#
# MSAM_SOURCE may be "remote" (default) or "local".

################## BUILD ARGUMENTS ######################
ARG MSAM_VERSION=1.2.0
ARG MSAM_SOURCE=remote

################## REMOTE SOURCE ######################
FROM alpine:3.23 AS source-remote

ARG MSAM_VERSION

RUN apk add --no-cache wget \
    && wget -q \
        "https://github.com/arumugamlab/msamtools/releases/download/${MSAM_VERSION}/msamtools-${MSAM_VERSION}.tar.gz" \
        -O /msamtools.tar.gz


################## LOCAL SOURCE ######################
FROM scratch AS source-local

ARG MSAM_VERSION

COPY msamtools-${MSAM_VERSION}.tar.gz /msamtools.tar.gz


################## SELECT SOURCE ######################
ARG MSAM_SOURCE
FROM source-${MSAM_SOURCE} AS source


################## BASE IMAGE ######################
FROM alpine:3.23

ARG MSAM_VERSION

################## METADATA ######################
LABEL base_image="alpine:3.23"
LABEL version="1"
LABEL software="msamtools"
LABEL software.version="${MSAM_VERSION}"
LABEL about.summary="microbiome-related extension to samtools"
LABEL about.home="https://github.com/arumugamlab/msamtools"
LABEL about.documentation="https://github.com/arumugamlab/msamtools"
LABEL about.license_file="https://github.com/arumugamlab/msamtools/blob/master/LICENSE"
LABEL about.license="SPDX:MIT"
LABEL about.tags="Metagenomics, Genomics"
LABEL extra.identifiers.biotools="msamtools"
LABEL org.opencontainers.image.authors="arumugam@sund.ku.dk"

################## INSTALLATION ######################

ENV MSAM_VERSION=${MSAM_VERSION}

COPY --from=source /msamtools.tar.gz /tmp/msamtools.tar.gz

RUN apk --no-cache update \
    && apk --no-cache upgrade \
    && apk add --no-cache argtable2 \
    && apk add --no-cache --virtual .build_deps \
        gcc \
        libc-dev \
        zlib-dev \
        make \
        bash \
        wget \
        bzip2 \
        argtable2-dev \
    && cd /tmp \
    && tar xfz msamtools.tar.gz \
    && cd "msamtools-${MSAM_VERSION}" \
    && HTSLIB_VERSION="$(awk '$1 == "htslib_version" {print $3}' versions.txt)" \
    && test -n "${HTSLIB_VERSION}" \
    && ./configure --prefix=/usr \
    && make install \
    && cd /tmp \
    && wget -q \
        "https://github.com/samtools/samtools/releases/download/${HTSLIB_VERSION}/samtools-${HTSLIB_VERSION}.tar.bz2" \
        -O samtools.tar.bz2 \
    && tar xjf samtools.tar.bz2 \
    && cd "samtools-${HTSLIB_VERSION}" \
    && ./configure \
        --without-curses \
        --disable-bz2 \
        --disable-lzma \
        --disable-libcurl \
        --without-libdeflate \
        --disable-lzma \
        --disable-s3 \
        --disable-gcs \
        --disable-plugins \
        --disable-ref-cache \
    && make samtools \
    && /usr/bin/install -c samtools /usr/bin/samtools \
    && /usr/bin/install -c -m0644 LICENSE /usr/share/licenses/msamtools/samtools-LICENSE \
    && cd /tmp \
    && rm -rf \
        "msamtools-${MSAM_VERSION}" \
        "samtools-${HTSLIB_VERSION}" \
        msamtools.tar.gz \
        samtools.tar.bz2 \
    && apk del .build_deps

ENTRYPOINT ["msamtools"]
