FROM python:3.8.13-slim AS build

# Set working directory
WORKDIR /app

ENV PATH="/root/miniconda3/bin:${PATH}"
#ARG PATH="/root/miniconda3/bin:${PATH}"

# Expose port for Streamlit
EXPOSE 8501

# Install git (to clone project), wget (to download miniconda), qnd libgl1 (for cadquery)
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

RUN conda install -c conda-forge conda-pack

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n myenv -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack


# The runtime-stage image; we can use Debian as the
# base image since the Conda env also includes Python
# for us.
FROM debian:buster AS runtime

RUN apt-get update && apt-get install -y libgl1-mesa-glx

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

# Make run commands use venv
#SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]
#RUN echo "export PATH="/root/miniconda3/bin:$PATH"" >>~/.bashrc \
#    && /bin/bash -c "source ~/.bashrc" \
#RUN echo "conda activate myenv" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Run GUI
COPY app/gui.py app/
COPY app/test.py app/
#COPY entrypoint.sh .
COPY app/models app/models/

WORKDIR /app

ENTRYPOINT source /venv/bin/activate && streamlit run gui.py
#ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "myenv", "streamlit", "run", "test.py"]
#ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "myenv", "streamlit", "run", "gui.py", "--server.port=8501", "--server.address=0.0.0.0"]
#RUN chmod +x entrypoint.sh
#ENTRYPOINT ["./entrypoint.sh"]