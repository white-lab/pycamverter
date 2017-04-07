# Change Log

## 0.3.2 (2017-04-06)

Features

  - Added auto-maybe for peptides with best ranking by MASCOT.

## 0.3.1 (2017-04-05)

Features

  - Added functionality to reprocess individual scans without PTM combination
    limits.

Bugfixes

  - Fixed main process hanging when a subprocess encountered an error.

## 0.3.0 (2017-03-24)

Features

  - Added support for ProteomeDiscoverer 2.1.
  - Several optimizations to get working on pS/T data sets.
  - Added pypy support.
  - Added support for runs with multiple raw files.
  - Limit generated sequences to 10 combinations per peptide to reduce pS/T
    file sizes.
  - Reduced memory usage by using python generators.
  - Use SQLite to store exported data. +This change is not backwards compatible
    with version of CAMV / pycamverter < 0.3.0+.

Bugfixes

  - Cleanup temp data files after runs.
  - Fixed py2.7 errors.

## 0.2.0 (2017-03-17)

Features

  - Added support for using ProteomeDiscoverer search files, alongside MASCOT.
  - Started annotating export files with generator's pycamverter version.
  - Support csv scan lists.

Bugfixes

  - Fixed issue with peptideData being blank when ambiguous peptides were
    present.
  - Fixed issues with viewing peptides with similar types of modifications
  - Merge proteins from ambiguous peptide matches

## 0.1.2 (2017-03-10)

Features

  - Upped ProteoWizard version to 3.0.10577
  - Added C13 isotope labeling
  - Improved speed of spectra matching

## 0.1.1 (2017-03-07)

Features

  - Use python multiprocessing to cut data processing time in half.

## 0.1.0 (2017-03-06)

Features

  - Initial release, spun off from pyproteome.
