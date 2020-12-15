FROM nfcore/base:1.10.2
LABEL authors="Marcos Cámara Donoso, Christina Chatzipantsiou, Athanasios Kousathanas" \
      description="Docker image containing all software requirements for the lifebit-ai/traits pipeline"

RUN apt-get update && \
    apt-get install -y \
              build-essential \
              git \
              unzip \
              autoconf \
              zlib1g-dev \
              libbz2-dev \
              liblzma-dev \
              libcurl4-gnutls-dev \
              libssl-dev \
              libgsl0-dev \
              libperl-dev \
              libxt-dev \
              speedtest-cli \
              procps

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/traits/bin:$PATH

RUN git clone https://github.com/bulik/ldsc.git && \
    cd ldsc

# Correct the shebang of the files to point to python2
RUN for i in `ls /ldsc/*py` ; do sed -i 's/python/python2/g'  $i; done
RUN for i in `ls /ldsc/ldscore/*py | grep -v __init__` ; do sed -i '1 i #!/usr/bin/env python2' $i ; done

RUN mkdir /opt/bin
COPY bin/* /opt/bin/

RUN find /ldsc/ -type f -iname "*.py" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.R" -exec chmod +x {} \;

RUN touch .Rprofile
RUN touch .Renviron

ENV PATH="$PATH:/ldsc/"
ENV PATH="$PATH:/opt/bin/"

USER root

WORKDIR /data/

CMD ["bash"]