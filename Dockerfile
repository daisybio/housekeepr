FROM rocker/shiny:4.0.3
MAINTAINER Markus List <markus.list@wzw.tum.de>

#install system packages
RUN apt-get update && apt-get install -y libxml2-dev libssl-dev libcurl4-openssl-dev libmariadbclient-dev libpq-dev libv8-dev liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev

#watchtower
LABEL com.centurylinklabs.watchtower.stop-signal="SIGKILL"

#copy shiny app to work dir
RUN mkdir /srv/housekeepr
WORKDIR /srv/housekeepr
COPY . .

#update shiny server conf and configure it to run housekeepr in single app mode
RUN sed -i 's/site_dir \/srv\/shiny-server;/app_dir \/srv\/housekeepr;/g' /etc/shiny-server/shiny-server.conf

#make annotationhub cache dir
RUN mkdir -p –m 777 /root/.cache/AnnotationHub

#install R packages via renv
ENV RENV_VERSION 0.12.2
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org')); \ 
          remotes::install_github('rstudio/renv@${RENV_VERSION}'); \
          renv::restore(); \
          library(AnnotationHub); \
          setAnnotationHubOption('ASK', FALSE); \
          ah <- AnnotationHub()"

#fix bug where app greys out on SSL connection
RUN echo 'sanitize_errors off;disable_protocols xdr-streaming xhr-streaming iframe-eventsource iframe-htmlfile;' >> /etc/shiny-server/shiny-server.conf

# make shiny owner of housekeepr-dir 
RUN chown -R shiny ./

EXPOSE 3838
