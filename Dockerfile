    # Use an official Conda image as a base
    FROM continuumio/miniconda3

    # Set environment variables
    ENV PATH /opt/conda/envs/TIPP/bin:$PATH

    # Install build-essential and g++
    RUN apt-get update && apt-get install -y build-essential g++ zlib1g-dev && rm -rf /var/lib/apt/lists/*

    # Create the TIPP environment and install required packages
    RUN conda init bash

    RUN bash -c "source /root/.bashrc && conda create --override-channels -c conda-forge -n TIPP -y" && \
    bash -c "source /root/.bashrc && conda activate TIPP && \
    conda install --override-channels -c conda-forge -c bioconda -c r \
       boost=1.85.0 \
       diamond \
       minimap2 \
       spoa \
       mcl \
       r-pheatmap \
       r-igraph \
       bioconductor-biostrings \
       r-stringdist \
       r-ggplot2 \
       trf \
       graphaligner \
       python=3.9 \
       flye -y && \
    pip install tiara==1.0.3 && \
    conda clean --all --yes"

    # Clone the TIPP repository and run the install script
    RUN git clone https://github.com/Wenfei-Xian/TIPP.git && \
    cd TIPP && \
    bash install.sh && \
    CURRENT_PATH=$(pwd)/src && \
    echo "export PATH=$CURRENT_PATH:\$PATH" >> ~/.bashrc && \
    echo "export PATH=$CURRENT_PATH/seqtk:\$PATH" >> ~/.bashrc && \
    echo "export PATH=$CURRENT_PATH/kmc3/bin:\$PATH" >> ~/.bashrc

    # Prepare commands
    RUN chmod +x /TIPP/src/*.pl && cp /TIPP/src/*.pl /usr/local/bin && \
    cp /TIPP/src/*.r /usr/local/bin && \
    cp /TIPP/src/readskmercount /usr/local/bin/ && \
    cp -r /TIPP/src/kmc3 /usr/local/bin/ && \
    cp -r /TIPP/src/seqtk /usr/local/bin/

    # Entry point
    CMD ["/bin/bash", "-c", "source /root/.bashrc && conda activate TIPP && bash"]
