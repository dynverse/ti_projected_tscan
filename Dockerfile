FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_cran("TSCAN")'

LABEL version 0.1.5

ADD . /code

ENTRYPOINT Rscript /code/run.R
