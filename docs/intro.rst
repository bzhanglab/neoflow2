Introduction
============

Neoflow is a streamlined computational workflow that integrates WES and MS/MS 
proteomics data for neoantigen prioritization to facilitate cancer immunotherapy.
It includes four modules: 

* variant annotation and customized database construction;
* variant peptide identification including MS/MS searching, FDR estimation, PepQuery validation;
* human leukocyte antigen (HLA) typing; and 
* MHC-binding prediction and neoantigen prioritization. 

The four modules are streamlined using `Nextflow <http://www.nextflow.io/>`_ and Docker. Neoflow supports both label free and iTRAQ/TMT data. 

The processes for each of the four modules are included in the following Nextflow files:

* variant annotation and customized database construction: ``db_construction.nf``;
* variant peptide identification including MS/MS searching, FDR estimation, PepQuery validation: ``msms.nf``;
* human leukocyte antigen (HLA) typing: ``hla_typing.nf``;
* MHC-binding prediction and neoantigen prioritization: ``neoantigen.nf``.
