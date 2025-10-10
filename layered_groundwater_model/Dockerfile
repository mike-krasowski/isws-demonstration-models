FROM python:3.13

#set the working directory inside of container ??????
WORKDIR /usr/local

# Copy your FloPy script and any other model files into the container
COPY ../ ./isws

# Install FloPy using pip
RUN pip install -r ./isws/config/requirements.txt


# Download the MODFLOW 6 executables (example version)
ARG MF_URL='https://github.com/MODFLOW-ORG/executables/releases/download/23.0/linux.zip'


# Download the executable and place it in the system's PATH
RUN apt-get update && apt-get install -y wget unzip ffmpeg && \
    wget -O /tmp/linux.zip $MF_URL && \
    unzip /tmp/linux.zip -d /tmp/linux && \
    cp /tmp/linux/mf6 /usr/local/bin/mf6 && \
    cp /tmp/linux/mp7 /usr/local/bin/mp7 && \
    chmod +x /usr/local/bin/mf6 && \
    chmod +x /usr/local/bin/mp7 && \
    rm -rf /tmp/linux*


# Make port 8889 available to the world outside this container
EXPOSE 85

# Run Jupyter Notebook when the container launches
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=85", "--no-browser", "--allow-root"]

# CMD ["python", "main_lx.py"]
