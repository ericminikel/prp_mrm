This repository holds the code and data for the following manuscript:

Minikel EV & Kuhn E (equal contribution), Cocco AR, Vallabh SM, Hartigan CR, Reidenbach AG, Safar JG, Raymond GJ, McCarthy MD, O’Keefe R, Llorens F, Zerr I, Capellari S, Parchi P, Schreiber SL, Carr SA. **Domain-specific quantification of prion protein in cerebrospinal fluid by targeted mass spectrometry.** _Mol Cell Proteomics_. 2019 Sep 26; [PMID: 31558565](https://pubmed.ncbi.nlm.nih.gov/31558565/). [Full text](https://www.mcponline.org/content/18/12/2388).

In this manuscript we describe prion protein (PrP) multiple reaction monitoring (MRM) &mdash; a targeted mass spec assay &mdash; and what we learned from applying this assay to clinical samples of cerebrospinal fluid from prion disease and non-prion disease patients.

In this repository you can:

+ View the [data as exported from Skyline](/data/skyline/), the [data after processing in R](/data/processed/), and the [script that does the processing](/src/process_data.R).
+ Download a [summary table of the MRM data for clinical samples](/data/summary/clinical_sample_mrm_results.tsv) and a [summary of the pilot LC-MS/MS DDA data](https://github.com/ericminikel/prp_mrm/blob/master/data/dda/csf_dda_sumbyintensity_190331.tsv).
+ Reproduce display items and statistics from the manuscript by running [src/mrm_figures.R](src/mrm_figures.R). This will write several figures and tables to the [display_items/script_generated/](display_items/script_generated/) directory; editable versions of display items created manually in PowerPoint or Excel are provided in [display_items/manual/](display_items/manual/).

Note that this manuscript has changed in revision since the bioRxiv posting. The current `master` branch reflects the revisions. For the latest code corresponding to the bioRxiv version, see the [repo as of June 10](https://github.com/ericminikel/prp_mrm/tree/d90c639b21aa2a2268598aaa5ca9e5b110b66f7a).

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

