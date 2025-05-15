# Use a base image with both Python and R
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy environment files
COPY snakemake/envs/*.yml /app/envs/

# Create conda environment
RUN conda env create -f envs/instrain.yml

# Install R and required packages
RUN conda install -c conda-forge r-base r-essentials r-renv

# Copy the pipeline code
COPY . /app/

# Set up R environment
RUN Rscript -e "renv::init()" && \
    Rscript -e "renv::restore()"

# Set environment variables
ENV PATH="/opt/conda/envs/instrain/bin:${PATH}"

# Set the entrypoint
ENTRYPOINT ["snakemake"]
CMD ["--help"] 