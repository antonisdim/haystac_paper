HAYSTAC: High-AccuracY and Scalable Taxonomic Assignment of MetagenomiC data 
===

Introduction 
------------

This is the official repository for the HAYSTAC paper. 

Alignment based metagenomics
----------------------------

`haystack` is an easy to use pipeline for metagenomic identifications 

You can easily:

1. Construct a database 
2. Prepare your sample for analysis, including trimming sequencing adapters, collapsing paired end reads and also downloading data from the SRA
3. Perform species identification 
4. Perform a chemical damage pattern analysis 

Performance tests
-----------------

All the code, setup script, run scripts and config files for the performance tests can be found under the performance_test directory.

You will need to clone the repository and then sequentially run the setup and run scripts from inside their directory.

Simulation datasets 
-------------------

For the simulation datasets you can fetch the database that was used in the apper and then follow the code examples. For the simulation we use gargammel (REF), so you will need to install that program first.

The relevant scripts can be found under the simulating_datasets directory. 

Analysis  
--------

For every analysis we have performed in this paper, we include the commands used in scripts under the analysis directory. 

License 
-------

MIT

Citation 
--------

Happening soon 
