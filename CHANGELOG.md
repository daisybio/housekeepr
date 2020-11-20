# 1.0-rc13
- bugfix: restoring from URL was broken due to recent changes

# 1.0-rc12
- fixes #3

# 1.0-rc11
- fixes #1

# 1.0-rc10
- related to #1
- fixed #2
- restructuring order of observers and renderers

# 1.0-rc9
- bugfix: selection of genes in table wouldn't correctly update charts

# 1.0-rc8
- improved responsiveness of results visualization

# 1.0-rc7
- improved responsiveness by showing boxes with spinner

# 1.0-rc6
- bugfix: multiple ensembl ids with same symbols led to an error
- duplicated symbols are now unificated by appending ensembl ids

# 1.0-rc5
- bugfix: selection in table sometimes failed to update plots correctly
- improved height of charts to avoid overlaps

# 1.0-rc4
- rounding numbers in table
- improve general responsiveness of app

# 1.0-rc3
- added presentation mode which is solely for showing results and status, but cannot start new analyses
- making session status more robust using database
- removed rownames of table; default ordering by rank
- various bug fixes

# 1.0-rc2
- showing result url while analysis is running
- showing time stamp of when result was produced
- if "most recent ensembl" is chosen, storing resolved ensembl releases instead for reproducibility

# 1.0-rc1
- first release candidate
