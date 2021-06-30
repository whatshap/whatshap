FROM ubuntu:16.04

# These commands copy your files into the specified directory in the image
# and set that as the working location
COPY . /app
WORKDIR /app

# Install fpa through conda
RUN apt-get update && \
    apt-get install -y curl make g++ zlib1g-dev && \
    curl -LO https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda install -c bioconda fpa minimap2
RUN cd BMEAN && ./install.sh && cd .. && make
ENV PATH=/app/:${PATH}

LABEL Name=CONSENT
