FROM continuumio/miniconda3@sha256:456e3196bf3ffb13fee7c9216db4b18b5e6f4d37090b31df3e0309926e98cfe2
LABEL authors="Marcos Cámara Donoso, Christina Chatzipantsiou, Athanasios Kousathanas" \
      description="Docker image containing all software requirements for the lifebit-ai/traits pipeline"

RUN apt-get --allow-releaseinfo-change update
RUN apt-get update && \
    apt-get install -y \
              build-essential \
              wget \
              unzip \
              autoconf \
              zlib1g-dev \
              libbz2-dev \
              liblzma-dev \
              libcurl4-gnutls-dev \
              libssl-dev \
              libgsl0-dev \
              libxt-dev \
              procps

COPY environment.yml /
RUN conda env create python=2.7 -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/traits/bin:$PATH

USER root

WORKDIR /data/

CMD ["bash"]