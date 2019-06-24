# Change Log

## 0.15.1 (2019-06-24)

Features

  - Updated ProteoWizard to latest version.

Bug fixes

  - Updated regex matcher for protein accession / description.

## 0.15.0 (2019-01-28)

Features

  - Added filtering by minimum ion score.

## 0.14.0 (2019-01-28)

Bug fixes

  - Fixed calculating b ion mass for unlabeled peptides.

## 0.13.9 (2019-01-28)

Bug fixes

  - Show output of msconvert.exe on error.

## 0.13.8 (2019-01-25)

Bug fixes

  - Don't error on unlabeled runs.

## 0.13.7 (2019-01-22)

Bug fixes

  - Fixed bug in phosphoRS module.

## 0.13.6 (2019-01-18)

Bug fixes

  - Fixed very slow scan processing on PD2.2 search files.


## 0.13.5 (2019-01-17)

Bug fixes

  - Fixed phosphoRS PD2.2 bug.

## 0.13.4 (2019-01-17)

Bug fixes

- Fixed phosphoRS for PD1.4, PD2.2.

## 0.13.3 (2019-01-17)

Bug fixes

  - Fixed os.path.dirname error on startup.

## 0.13.2 (2019-01-10)

Bug fixes

  - Set pycamv.log path to same folder as .msf file.

## 0.13.1 (2019-01-04)

Bug fixes

  - Fixed errors in building PyCAMVerter.exe.

## 0.13.0 (2019-01-03)

Features

  - Added support for ProteomeDiscoverer 2.2.

## 0.12.0 (2018-04-30)

Features

  - Added support for phosphoRS phosphorylation assignments.

## 0.11.2 (2018-04-04)

Bug fixes

  - Fixed win32api error in deployments.

## 0.11.1 (2018-02-07)

Bug fixes

  - Fix crash on dynamic N/C-term modifications.
  - Fix crash on empty workflow elements.

## 0.11.0 (2017-11-10)

Features

  - Add protein sequence, uniprot accessions, peptide offsets to output.

## 0.10.0 (2017-09-07)

Features

  - Add fragment intensity info to database, allowing quick export of
    quantitation data.

## 0.9.0 (2017-08-30)

Bug fixes

  - Fixed bug in calculating the mass of parent ion fragments with C13 additions.
  - Only estimate max C13 number using ions within precursor isolation window.

## 0.8.3 (2017-07-19)

Features

  - Use low priority threads during processing to avoid impacting other tasks.

## 0.8.2 (2017-06-08)

Bug fixes

  - Fixed crash when X amino acids included modifications other than "Mapping*".

## 0.8.1 (2017-05-25)

Bug fixes

  - Fixed crash when no peaks are detected in the precursor window.

## 0.8.0 (2017-05-24)

Features

  - Determine C13 fragments based on the number of C13 peaks seen within the
    isolation window.

## 0.7.0 (2017-05-17)

Features

  - Added support for TMT11plex

## 0.6.1 (2017-05-15)

Bug fixes

  - Fixed some cases where intermediate processing files were not deleted when
    a processing run was interrupted or errored.
  - Fixed several type errors.

## 0.6.0 (2017-05-12)

Features

  - Added the ability to limit CPU processes via command line arguments.

Bug fixes

  - Allow for mapping of "X" in a peptide sequence to any amino acid
    (not just S/N).
  - Fixed bug in determining precursor peptide charge in rare cases.
  - Fixed crash when multiple proteins exist with the same description.

## 0.5.0 (2017-05-04)

Bug fixes

  - Allow for alternative scans for quant scan

## 0.4.0 (2017-05-02)

Features

  - Added the ability to import validation data from CAMV-Matlab sessions.

## 0.3.3 (2017-04-11)

Features

  - Improved performance of scan lookups with SQL indices.

Bug fixes

  - Fixed crash on MASCOT xml search data.

## 0.3.2 (2017-04-06)

Features

  - Added auto-maybe for peptides with best ranking by MASCOT.

## 0.3.1 (2017-04-05)

Features

  - Added functionality to reprocess individual scans without PTM combination
    limits.

Bug fixes

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

Bug fixes

  - Cleanup temp data files after runs.
  - Fixed py2.7 errors.

## 0.2.0 (2017-03-17)

Features

  - Added support for using ProteomeDiscoverer search files, alongside MASCOT.
  - Started annotating export files with generator's pycamverter version.
  - Support csv scan lists.

Bug fixes

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
