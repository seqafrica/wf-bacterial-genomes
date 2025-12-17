<!---This section of documentation typically contains a list of things the workflow can perform also any other intro.--->

This workflow is used to produce long-read de novo assemblies of bacterial and fungal genomes, or to perform reference-based variant calling.

The workflow also performs analysis of the assemblies, such as species identification, antimicrobial resistance (AMR) analysis, and sequence typing through an optional `--isolates` mode.

In brief, this workflow will perform the following: 

+ De novo (or reference-based) assembly of bacterial and fungal genomes 
+ Variant calling (`--reference` mode only)
+ Annotation of regions of interest within the assembly
+ Identification of plasmid regions within the assembly
+ MLST and genomic similarity based species identification (`--isolates` mode only)
+ Identify genes and SNVs associated with AMR (`--isolates` mode only)

<figure>
<img src="docs/images/wf-bac-genomes.drawio.png" alt="wf-bacterial-genomes overview schematic."/>
<figcaption>Schematic depicting wf-bacterial-genomes workflow.</figcaption>
</figure>
