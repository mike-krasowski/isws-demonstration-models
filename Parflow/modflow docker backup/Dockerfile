FROM python:3.13

#set the working directory inside of container ??????
WORKDIR /usr/local

# Copy your FloPy script and any other model files into the container
COPY .. ./isws

# Install FloPy using pip
# RUN pip install jupyter flopy numpy matplotlib
RUN pip install -r ./isws/config/requirements.txt


# Download the MODFLOW 6 executable (example version)
ARG MF6_VERSION="mf6.6.3_linux"
ARG MF6_URL='https://github.com/MODFLOW-ORG/modflow6/releases/download/6.6.3/mf6.6.3_linux.zip'

# Download the executable and place it in the system's PATH
RUN apt-get update && apt-get install -y wget unzip && \
    wget -O /tmp/mf6.zip $MF6_URL && \
    unzip /tmp/mf6.zip -d /tmp/mf6 && \
    cp /tmp/mf6/$MF6_VERSION/bin/mf6 /usr/local/bin/mf6 && \
    chmod +x /usr/local/bin/mf6 && \
    rm -rf /tmp/mf6*


# Make the executable runable
# RUN chmod +x /usr/local/bin/linux/mf6

# Make port 8889 available to the world outside this container
EXPOSE 80

# # Add the MODFLOW 6 executable to the system's PATH
# ENV PATH="{C:/Users/esrag/Documents/GitHub/isws-demonstration-models/layered_groundwater_model/config}"

# # This could be running a Python script, starting a Jupyter notebook, etc.
# Run Jupyter Notebook when the container launches


CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=80", "--no-browser", "--allow-root"]

# CMD ["python", "main.py"]

