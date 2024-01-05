# Build stage first: create env with conda, then use conda pack to reduce size of env
# See https://pythonspeed.com/articles/conda-docker-image-size/
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install git (to clone project) and other stuff which Streamlit apparently needs
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    && rm -rf /var/lib/apt/lists/*

# Clone project into container
RUN git clone -b streamlit https://github.com/danieljuschus/pemfc-aircraft-sizing

# Workdir into cloned repo
WORKDIR /app/pemfc-aircraft-sizing

# Install Python requirements
RUN pip install --no-cache-dir -r requirements_gui.txt

# Expose port for Streamlit
EXPOSE 8501

# Run GUI
WORKDIR /app/pemfc-aircraft-sizing/app
ENTRYPOINT ["streamlit", "run", "gui.py"]