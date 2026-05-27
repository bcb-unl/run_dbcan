FROM mambaorg/micromamba:2-ubuntu22.04 AS build

USER $MAMBA_USER

WORKDIR /tmp

# Conda: Python + DIAMOND (not bioconda dbcan — avoids channel lag vs GitHub releases)
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.deps.yml /tmp/environment.deps.yml
RUN micromamba env create -f /tmp/environment.deps.yml \
  && micromamba clean --all --yes \
  && rm /tmp/environment.deps.yml

ENV PATH=/opt/conda/envs/run_dbcan/bin:$PATH
ENV MAMBA_DEFAULT_ENV=run_dbcan

# Install dbcan from the same tree as the Release tag (passed as BUILD_VERSION for setuptools-scm)
ARG BUILD_VERSION=0.0.0.dev0
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${BUILD_VERSION}

COPY --chown=$MAMBA_USER:$MAMBA_USER . /src
WORKDIR /src

RUN pip install --no-cache-dir hatchling "hatch-vcs" "setuptools-scm>=7,<8" \
  && pip install --no-cache-dir .

FROM mambaorg/micromamba:2-ubuntu22.04 AS final

ARG BUILD_VERSION=0.0.0.dev0

LABEL base_image="micromamba:2-ubuntu22.04"
LABEL version="1"
LABEL software="dbcan and dbcan utils"
LABEL software.version="${BUILD_VERSION}"
LABEL about.summary="Docker image for dbcan and dbcan utils"
LABEL about.tags="annotation, metagenomics, dbcan"
LABEL about.license="MIT"

COPY --from=build /opt/conda/envs/run_dbcan /opt/conda/envs/run_dbcan

ENV ENV_NAME=run_dbcan
ENV PATH=/opt/conda/envs/run_dbcan/bin:$PATH
ENV MAMBA_DEFAULT_ENV=run_dbcan

WORKDIR /app

CMD ["run_dbcan", "--help"]

MAINTAINER Xinpeng Zhang <xzhang55@huskers.unl.edu>, Haidong Yi <hyi@stjude.org>
