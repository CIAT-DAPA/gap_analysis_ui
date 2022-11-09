FROM rocker/shiny:4.1.0

RUN rm -rf /srv/shiny-server/*

WORKDIR /srv/shiny-server/

RUN apt-get clean \
    && apt-get update \
    && apt-get install -y \
        default-jdk 
    #&& Rscript /srv/shiny-server/imports.R

COPY ./src/ ./

# docker build -t stevensotelo/gap_analysis_ui:latest .
# docker tag xx stevensotelo/gap_analysis_ui:latest

# docker run --name gap_ui --rm -p 3838:3838 stevensotelo/gap_analysis_ui:latest