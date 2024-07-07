FROM python:3.9.19-slim-bullseye

RUN pip install requests
RUN apt update
RUN apt install -y --no-install-recommends git
RUN git clone https://github.com/lculibrk/vep_project
WORKDIR vep_project
RUN pip install ./

