FROM nfcore/base
LABEL authors="Sam Minot, Ph.D." \
      description="Docker image containing all requirements for nf-core/aligngenomes pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-aligngenomes-1.0dev/bin:$PATH
