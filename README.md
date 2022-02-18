# FORALL

Guidelines for the R project
Data and analysis results are stored in ExpressionSet objects, figures and tables.

The proteins eSet is generated in init.R.

Output is stored in the output directory.

Several scripts (analysis*.R) perform functions on the proteins object. The scripts read in a proteins object and their output is stored in respective subfolders of output. This can be tables, figures or modified proteins objects.

Accessory functions are stored in functions and sourced in analysis scripts.

Each figure is built within a figure_*.R script. 


Publication-ready figure components are stored in output/figures.



Sample data objects are provided in the repository for illustration and convenience, these files may not be complete

To run the workflow
Copy raw data to ./data, e.g. symbols_table.txt (Nextflow output).


Copy meta data files to ./meta.


init.R


analysis*.R


figure_*.R


Find results in ./output.
