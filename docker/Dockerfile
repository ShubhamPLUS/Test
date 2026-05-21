## ProteoForge — batch analysis image
## Base: Bioconductor release with R 4.4+
FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="ProteoForge"
LABEL org.opencontainers.image.description="Quantitative proteomics analysis platform"
LABEL org.opencontainers.image.source="https://github.com/proteoforge/proteoforge"
LABEL org.opencontainers.image.licenses="MIT"

## System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      libfontconfig1-dev \
      libcairo2-dev \
      libxt-dev \
      pandoc \
      quarto \
    && rm -rf /var/lib/apt/lists/*

## renv for reproducible package restore
RUN Rscript -e "install.packages('renv', repos='https://cloud.r-project.org')"

## Copy lockfile and restore packages (cached layer)
WORKDIR /app
COPY renv.lock ./
RUN Rscript -e "renv::restore(prompt=FALSE)"

## Copy package source and install
COPY . .
RUN Rscript -e "devtools::install('.', dependencies=TRUE, upgrade='never')"

## Entry point: CLI script
ENTRYPOINT ["Rscript", "/app/inst/scripts/proteoforge_cli.R"]
CMD ["--help"]
