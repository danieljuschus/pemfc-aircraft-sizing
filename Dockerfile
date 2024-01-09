FROM python:3.11-slim

WORKDIR /app

# Install stuff which Streamlit apparently needs
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    libgl1-mesa-glx libxrender1\
    xvfb\
    && rm -rf /var/lib/apt/lists/*

# Copy entire repository into workdir
COPY . .

# Install Python requirements
RUN pip install --no-cache-dir -r requirements_gui.txt

# Expose port for Streamlit
EXPOSE 8501

# Run GUI
WORKDIR /app/app
ENTRYPOINT ["streamlit", "run", "gui.py"]
