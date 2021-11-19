FROM rocker/shiny
MAINTAINER Direct All Questions Directly To Devin Schweppe "dkschwep@uw.edu"

RUN R -e "install.packages(c('shiny', 'rmarkdown'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('shinydashboard','RSQLite','stringr','stringi','dplyr','readr','data.table','purrr','blob','shinybusy'), repos='https://cran.rstudio.com/')"


#Stage Shiny files, etc.
RUN mkdir /root/Repos/
COPY . /root/Repos/MSPtoDB/

EXPOSE 3838

CMD ["Rscript", "/root/Repos/MSPtoDB/shiny/app.R"]
