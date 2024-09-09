FROM ghcr.io/apptainer/apptainer:latest

ARG DEBIAN_FRONTEND=noninteractive
 
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config \
    git-all \
    wget \
    curl \
    libbz2-dev \
    zlib1g-dev \
    python3-dev \
    libffi-dev \
    r-base \
    r-base-dev 

#rust install
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y && \
    . /root/.cargo/env && \
    cargo build --release

#julia install 
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.5-linux-x86_64.tar.gz --directory-prefix=/ && \
    tar zxvf julia-1.10.5-linux-x86_64.tar.gz && \
    export PATH="$PATH:/julia-1.10.5/bin" && \
    julia -e 'import Pkg; Pkg.add(["Glob", "CSV", "DataFrames", "CodecZlib", "ArgParse"])'

#cargo install proseg
RUN cargo install proseg

ENV NUMBA_CACHE_DIR=/work/numba_cache
ENV MPLCONFIGDIR=/work/mpl_cache


#RUN wget https://dl.fedoraproject.org/pub/epel/8/Everything/x86_64/Packages/s/singularity-3.8.0-1.el8.x86_64.rpm && \
#    alien -d singularity-3.8.0-1.el8.x86_64.rpm && \
#    apt-get install ./singularity_3.8.0-2_amd64.deb
#RUN curl -L https://dl.google.com/go/go1.17.5.linux-amd64.tar.gz -o go1.17.5.linux-amd64.tar.gz
#RUN tar -C /usr/local -xzf go1.17.5.linux-amd64.tar.gz
#RUN rm go1.17.5.linux-amd64.tar.gz
#RUN chmod -R 777 /cmdstan/*
#RUN R -e "library(cmdstanr);cmdstanr::set_cmdstan_path(path = list.dirs('/cmdstan')[[2]]);cpp_options <- list('CXX' = 'clang++','CXXFLAGS+'= '-march=native',PRECOMPILED_HEADERS = FALSE);rebuild_cmdstan()"
#RUN git clone https://github.com/stan-dev/cmdstan.git -b v2.32.2 /home/cmdstan.github --recursive 
