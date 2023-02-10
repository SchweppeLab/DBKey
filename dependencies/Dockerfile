FROM rocker/shiny
MAINTAINER Direct All Questions Directly To Devin Schweppe "dkschwep@uw.edu"

RUN R -e "install.packages(c('shiny','rmarkdown','shinydashboard','RSQLite','stringr','stringi','dplyr','readr','data.table','purrr','blob','shinybusy','BiocManager','compiler','vroom','shinyWidgets','shinyjs'), repos='https://cran.rstudio.com/')"

RUN R -e "BiocManager::install('BiocParallel')"
