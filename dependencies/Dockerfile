FROM rocker/shiny
MAINTAINER Direct All Questions Directly To Devin Schweppe "dkschwep@uw.edu"


RUN R -e "install.packages(c('shiny','rmarkdown','shinydashboard','RSQLite','stringr','stringi','dplyr','readr','data.table','purrr','blob','shinybusy','BiocParallel','compiler','vroom','shinyWidgets'), repos='https://cran.rstudio.com/')"


RUN R -e "install.packages(c('shiny','rmarkdown','shinydashboard','RSQLite','stringr','stringi','dplyr','readr','data.table','purrr','blob','shinybusy'), repos='https://cran.rstudio.com/')"

COPY ./shiny-app/* /srv/shiny-server/
CMD ["/usr/bin/shiny-server"]
