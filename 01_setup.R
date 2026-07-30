# Cannabis Use and Everyday Functioning Among People Living With HIV
# Script: 01_setup.R
# Purpose: 01 Setup
#
# Remove participant identifiers, local usernames, credentials,
# restricted output, and machine-specific file paths before publication.

rm(list = ls())

data_dir <- "data"
table_dir <- "output/tables"
figure_dir <- "output/figures"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)

# Add library() calls below.
# Save final session information with:
# writeLines(capture.output(sessionInfo()), "session-info.txt")
