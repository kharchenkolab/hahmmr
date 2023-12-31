% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflows.R
\name{get_bulk}
\alias{get_bulk}
\title{Produce combined bulk expression and allele profile}
\usage{
get_bulk(
  count_mat,
  lambdas_ref,
  df_allele,
  gtf,
  genetic_map = NULL,
  min_depth = 0,
  nu = 1,
  verbose = TRUE
)
}
\arguments{
\item{count_mat}{matrix Gene expression counts}

\item{lambdas_ref}{matrix Reference expression profiles}

\item{df_allele}{dataframe Allele counts}

\item{gtf}{dataframe Transcript gtf}

\item{genetic_map}{dataframe Genetic map}

\item{min_depth}{integer Minimum coverage to filter SNPs}

\item{nu}{numeric Phase switch rate}

\item{verbose}{logical Whether to print progress}
}
\value{
dataframe Pseudobulk gene expression and allele profile
}
\description{
Produce combined bulk expression and allele profile
}
\examples{
bulk_example = get_bulk(
    count_mat = gene_counts_example,
    lambdas_ref = ref_hca,
    df_allele = df_allele_example,
    gtf = gtf_hg38)
}
