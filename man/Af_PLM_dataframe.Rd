% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Af_PLM_dataframe.R
\name{Af_PLM_dataframe}
\alias{Af_PLM_dataframe}
\title{Function to create a dataframe of the Protein Language Model probabilities and ranks of the mutations along the edges of B cell lineage trees.}
\usage{
Af_PLM_dataframe(AntibodyForests_object, sequence.name, path_to_probabilities)
}
\arguments{
\item{AntibodyForests_object}{AntibodyForests-object, output from Af_build()}

\item{sequence.name}{character, name of the sequence column in the AntibodyForests object (example VDJ_sequence_aa_trimmed)}

\item{path_to_probabilities}{character, path to the folder containing probability matrices for all sequences. Probability matrices should be in CSV format and the filename should include sampleID_clonotypeID_nodeNR, matching the AntibodyForests-object.}
}
\value{
A dataframe with the sample, clonotype, node numbers, number of substitutions, mean substitution rank and mean substitution probability
}
\description{
Function to create a dataframe of the Protein Language Model probabilities and ranks of the mutations along the edges of B cell lineage trees.
}
\examples{
\dontrun{
PLM_dataframe <- Af_PLM_dataframe(AntibodyForests_object = AntibodyForests::small_af,
                 sequence.name = "VDJ_sequence_aa_trimmed",
                 path_to_probabilities = "/directory/ProbabilityMatrix")
}
}
