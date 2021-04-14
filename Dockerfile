# Get shiny server plus tidyverse image
FROM rocker/shiny-verse:4.0.4

# Install system libraries
RUN apt-get update && apt-get install -y \
    sudo \
    mono-devel \ 
    mono-complete \
    xdg-utils
  
# Install R packages 
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('andrewsali/shinycssloaders')"
RUN R -e "install.packages('magrittr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plotly', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('hexbin', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('protViz', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('RSQLite', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('scales', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinythemes', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('rlang', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('testthat', repos='http://cran.rstudio.com/')"

# Install rawDiag
RUN R -e "install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawDiag_0.0.41.tar.gz', repo=NULL)"

RUN R -e "devtools::install_github('andrewsali/plotlyBars')"

# Copy the app and files to the image
COPY app.R /srv/shiny-server/
COPY files /srv/shiny-server/files

# Select port
EXPOSE 3838

# Run app
CMD ["R", "-e", "library(shiny); source('/srv/shiny-server/app.R');"]