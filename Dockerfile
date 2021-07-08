FROM rocker/shiny:4.0.3
MAINTAINER Markus List <markus.list@wzw.tum.de>

#install system packages
RUN apt-get update && apt-get install -y libxml2-dev libssl-dev libcurl4-openssl-dev libmariadbclient-dev libpq-dev libv8-dev liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev

#watchtower
#LABEL com.centurylinklabs.watchtower.stop-signal="SIGKILL"

#set work dir
RUN mkdir /srv/housekeepr
WORKDIR /srv/housekeepr

#install R packages via renv
ENV RENV_VERSION 0.12.3
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY renv.lock renv.lock
RUN R -e "renv::restore()"

# make shiny owner of housekeepr-dir 
RUN chown -R shiny:shiny /srv

# copy shiny app
COPY ./*.R ./
COPY ./init_data.RData .
COPY ./www ./www

# download AnnotationHub to speed up loading time
USER shiny
RUN R -e "library(AnnotationHub); \
          setAnnotationHubOption('CACHE', '~/.myHub'); \
          ah <- AnnotationHub(ask = FALSE)"
          
USER root


EXPOSE 3838
