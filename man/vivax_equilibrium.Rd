% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{vivax_equilibrium}
\alias{vivax_equilibrium}
\title{Vivax equilibrium solution with biting heterogeneity}
\usage{
vivax_equilibrium(age, ft, EIR, p, v_eq = "full")
}
\arguments{
\item{age}{vector of age groups, in units of years}

\item{ft}{proportion of clinical cases effectively treated}

\item{EIR}{EIR for adults, in units of infectious bites per person per year}

\item{p}{vector of model parameters}
}
\description{
Returns the vivax equilibrium states for the model of White et al.,
  (2018).
}
