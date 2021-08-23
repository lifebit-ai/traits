FROM continuumio/miniconda3@sha256:456e3196bf3ffb13fee7c9216db4b18b5e6f4d37090b31df3e0309926e98cfe2

RUN apt-get update --allow-releaseinfo-change
RUN apt-get install -y \
              wget \
              procps


COPY environment.yml /
RUN conda env create python=2.7 -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/traits/bin:$PATH


USER root

WORKDIR /data/

CMD ["bash"]