# Build stage first: create env with conda, then use conda pack to reduce size of env
# See https://pythonspeed.com/articles/conda-docker-image-size/
FROM python:3.8.13-slim AS build

# Set working directory
WORKDIR /app

# Not sure why this is necessary
ENV PATH="/root/miniconda3/bin:${PATH}"

# Expose port for Streamlit
EXPOSE 8501

# Install git (to clone project), wget (to download miniconda), qnd libgl1 (for cadquery) - maybe the latter can be removed here
RUN apt-get update && apt-get install -y git wget libgl1-mesa-glx && rm -rf /var/lib/apt/lists/*

# Download and install miniconda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# Clone project into container
RUN git clone -b streamlit https://github.com/danieljuschus/pemfc-aircraft-sizing

# Create conda venv
COPY environment.yml .
RUN conda env create -f environment.yml

# Install conda-pack to venv
RUN conda install -c conda-forge conda-pack

# Use conda-pack to create a standalone enviornment with smaller size in /venv:
RUN conda-pack -n myenv -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# Fix up paths
RUN /venv/bin/conda-unpack

# The runtime-stage image; we can use Debian as the base image since the Conda env also includes Python for us.
FROM debian:buster AS runtime

# Install libgl1 for cadquery
RUN apt-get update && apt-get install -y libgl1-mesa-glx

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

# Make run commands use venv
SHELL ["/bin/bash", "--login", "-c"]

# Copy all necessary files and set working dir
COPY app/gui.py app/
COPY app/models app/models/
WORKDIR /app

# Run GUI
ENTRYPOINT source /venv/bin/activate && streamlit run gui.py