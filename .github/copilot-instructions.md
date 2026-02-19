## Core Context
This is a repo for managing viral datasets. It pulls from genbank with IDs, and is currently run through the nextflow script vgtk-init.nf

It is designed to be flexible on viruses, allowing segemented and not, the main viruses for now are RABV, flu, and HCV

For all new code, add tests to /home3/oml4h/RABV-gTK/tests/unit
## data info
It needs running in the conda env vgtk, which can be set up with the environment.yml file in the repo.


## ongoing issues
Optionally including gisaid data
Allowing updates rather than pulling all data from genbank each time and re-processing massive datasets

## behaviour preferences

ask clarifying questions on architecture etc. 