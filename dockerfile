# Use a minimal base Ubuntu image
FROM ubuntu:20.04

# Avoid timezone prompt during packages installation
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary packages
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    git \
    libgomp1 \
    locales \
    && rm -rf /var/lib/apt/lists/*

# Set up environment
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

# Install Miniconda and Mamba
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda install -c conda-forge mamba -y && \
    /opt/conda/bin/conda clean -ya

# Create environment and install packages
RUN mamba create -n switch_finder_env python && \
    /bin/bash -c "source activate switch_finder_env && \
    mamba install -y -c conda-forge pandas numpy scikit-learn statsmodels && \
    mamba install -y -c bioconda viennarna"

# Copy necessary files
COPY programs /programs

# Set environment variables for scripts
ENV RNAstructure_path=/programs/RNAstructure
ENV RNApathfinder_path=/programs/RNApathfinder
ENV DATAPATH=/programs/RNAstructure/data_tables

WORKDIR switchfinder

CMD [ "/bin/bash" ]

RUN echo "source activate switch_finder_env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
#ENTRYPOINT ["python", "app.py"]