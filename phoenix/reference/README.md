# WinNonlin reference outputs

Two kinds of file live here:

**Fetched (gitignored).** `run_all.R` downloads the published Phoenix
WinNonlin 6.3/7.0 output for `datasets::Theoph` from the NonCompart validation
repository. Not vendored, so provenance stays with the source.

**Yours (commit these).** Export "Final Parameters Pivoted" from your own
Phoenix NCA object following `../PHOENIX_SETUP.md` and save it here with any
filename not matching `winnonlin_theoph_*.csv`. `run_all.R` prefers a local
export over the fetched reference.

When you commit your own export, record alongside it:

- Phoenix version and build (Help -> About)
- Licence tier and enabled modules
- Calculation settings as actually set
- SHA-256 of the input file

Without that, the comparison is unreproducible a year from now.
