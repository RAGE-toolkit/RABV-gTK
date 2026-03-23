# Building RABV-gDB

Building **RABV-gDB** is a straightforward process. All required prerequisite files are located in the `generic/rabv` directory. If you want to modify the existing reference set, follow the steps below.

## 1. Editing the Reference File

The reference file is located at:

`generic/rabv/ref_list.txt`

This file contains the accession numbers used as the reference database for the BLAST step. These accessions act as the core reference set against which all other sequences are mapped.

To add new reference sequences in the future, simply append the relevant accession numbers to `ref_list.txt`.

## 2. Building the Reference Alignment and Clade Assignments

The reference alignment is generated using a diverse set of rabies virus sequences selected from the following studies. Major and minor clade assignments are also based on these publications:

- **Troupin et al. (2016)**  
  *Large-Scale Phylogenomic Analysis Reveals the Complex Evolutionary History of Rabies Virus in Multiple Carnivore Hosts*  
  *PLoS Pathogens*

- **Kuzmin et al. (2012)**  
  *Molecular Inferences Suggest Multiple Host Shifts of Rabies Viruses from Bats to Mesocarnivores in Arizona during 2001–2009*  
  *PLoS Pathogens*

In addition, reference sequences for the **AM5** and **AM3c** minor clades were included from the following studies:

- **Brunker et al. (2025)**  
  *Genomic characterization of a dog-mediated rabies outbreak in El Pedregal, Arequipa, Peru*  
  *PLoS Neglected Tropical Diseases*  
  PMID: 40043048

- **Caraballo et al. (2021)**  
  *A Novel Terrestrial Rabies Virus Lineage Occurring in South America: Origin, Diversification, and Evidence of Contact between Wild and Domestic Cycles*  
  *Viruses*  
  PMID: 34960753

The selected sequences were aligned using **MAFFT** with default parameters.

## 3. Building the Reference Tree

The reference tree is constructed from the reference alignment sequences using **IQ-TREE** with default parameters.

We also evaluated the available reference sequences to identify the best-fitting substitution model before finalizing the tree.

## 4. Creating the Gene Information Table

The gene information table is located at:

`generic/rabv/Tables/gene_info.csv`

This table is used in the web version of gDB:

http://gdb-web-dev.cvr.gla.ac.uk/

It provides gene and protein annotation details for visualization in the web interface. All genes/proteins from query sequences are annotated based on this table.

If CDS start or end coordinates are modified, the same changes must also be updated in `gene_info.csv` to ensure accurate visualization.
