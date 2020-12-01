# warbler_inversions


### This is a collection of workflows and scripts used in "Inversions maintain differences between migratory phenotypes of a songbird"


#### General workflows are available in the *workflows* folder and have been divided into:

##### Genome assembly and curation: Assemblies.txt
##### Repeat and gene annotation: Annotations.txt
##### Mapping and variant calling of resequencing data: Mapping_variant_calling.txt
##### Analyses of bionano and linked read data: Bionano_linked read data
##### Dating analyses: Dating.txt

#### Customized scripts and configuration files can be found in the *scripts* folder


#### Estimate coverage from barcoded molecules

##### Extract alignments from the focal scaffolds 
```
samtools view -h P6352_101.sorted.bam Scaffold68 | samtools sort -tBX - > P6352_101.Scaffold68.bc.bam
samtools view -h P6352_102.sorted.bam Scaffold68 | samtools sort -tBX - > P6352_102.Scaffold68.bc.bam

samtools view -h P6352_101.sorted.bam Scaffold156 | samtools sort -tBX - > P6352_101.Scaffold156.bc.bam
samtools view -h P6352_102.sorted.bam Scaffold156 | samtools sort -tBX - > P6352_102.Scaffold156.bc.bam

samtools view -h P6352_101.sorted.bam Scaffold29:55700000- | samtools sort -tBX - > P6352_101.Scaffold29.bc.bam
samtools view -h P6352_102.sorted.bam Scaffold29:55700000- | samtools sort -tBX - > P6352_102.Scaffold29.bc.bam

```
