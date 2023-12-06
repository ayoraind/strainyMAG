## Workflow to resolve strains from long-read assemblies using Strainy.
### Usage

```

=================================================================================
 STRAIN RESOLUTION OF ASSEMBLIES USING STRAINY: TAPIR Pipeline version 1.0dev
=================================================================================
 The typical command for running the pipeline is as follows:
        nextflow run main.nf --reads "PathToReadFile(s)" --output_dir "PathToOutputDir" 

        Mandatory arguments:
         --reads                        Query fastqz file of sequences you wish to supply as input (e.g., "/MIGE/01_DATA/01_FASTQ/T055-8-*.fastq.gz")
         --output_dir                   Output directory (e.g., "/MIGE/01_DATA/03_ASSEMBLY")
         
        Optional arguments:
         --help                         This usage statement.
         --version                      Version statement

```


## Introduction
This pipeline attempts to resolve strains from long-read assemblies. This Nextflow pipeline was adapted from the [Strainy github page](https://github.com/katerinakazantseva/stRainy/) and the [metagenomeStrainy_ONT_pipeline github page](https://github.com/katerinakazantseva/MetagenomeStrainy_ONT_pipeline/).  
The CheckM2 database directory (has to take the name 'CheckM2_database') is expected to be present in the same directory hosting the main.nf file. Because github does not accept files of size greater than 1 GB, a typical CheckM2_database is not present in this github page. A typical CheckM2 database can be downloaded [here](https://zenodo.org/records/5571251).

## Sample command
An example of a command to run this pipeline is:

```
nextflow run main.nf --reads "Sample_files/*.fastq.gz" --output_dir "test2" 
```


## Word of Note
This is an ongoing project at the Microbial Genome Analysis Group, Institute for Infection Prevention and Hospital Epidemiology, �niversit�tsklinikum, Freiburg. The project is funded by BMBF, Germany, and is led by [Dr. Sandra Reuter](https://www.uniklinik-freiburg.de/institute-for-infection-prevention-and-control/microbial-genome-analysis.html).


## Authors and acknowledgment
The TAPIR (Track Acquisition of Pathogens In Real-time) team.
