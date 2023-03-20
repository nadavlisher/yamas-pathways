# YaMAS (YOLO Microbiome Analysis System)

YaMAS is a package designed to easily download DNA datasets from the NCBI SRA website. It is developed by the YOLO lab team, and is designed to be simple, efficient, and easy to use for non-programmer users.

## Installation

To install YaMAS, you can use pip:

```
pip install yamas
```

## Dependencies
Before proceeding with the installation of YaMAS, please make sure that all dependencies are fulfilled. In case any of the dependencies are missing, the program will not run as expected. Please refer to the installation instructions and ensure that all requirements are met before proceeding.
- YaMAS should be downloaded in a [qiime2](https://docs.qiime2.org/2023.2/) enviorment.
- [SRA-toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) package should be downloaded in the enviorment.
- [Entrez](https://anaconda.org/bioconda/entrez-direct) package should be downloaded in the enviorment.
- Exporting a project requires a downloaded [classifier file](https://data.qiime2.org/2022.8/common/gg-13-8-99-nb-classifier.qza).

## Getting Started

YaMAS provides an easy-to-use interface in the terminal.

To download the visualization file for multiple projects and save each file to a separate folder, use the following command:
```
yamas --download PRJEB12345
```
To export an OTU (Operational Taxonomic Unit), taxonomy, and phylogeny tree for a single project, use the following command:
```
yamas --export project_path trim trunc classifier_file threads
```
Arguments:
- project_path: path to the project directory (created by YaMAS in the previous step).
- classifier_file: path to the trained classifier file. 
- trim & trunc: choose graph edges. 
- threads: specifies the number of threads to use for parallel processing, which can speed up the export process (default is 12).

