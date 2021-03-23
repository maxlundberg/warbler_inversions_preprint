# Breakpoint analyses


## Whole-genome alignments

Align the southern assembly to the northern assembly using `SatsumaSynteny`

```
southern_genome="ww_southern_10x_bionano.filt.new.fasta"
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta"

SatsumaSynteny -q $southern_genome -t $genome -n 19 -o $genome.southern_ww.final.satsuma > $genome.southern_ww.final.satsuma.stdout
```

Also align the northern assembly to the genome of the collared flycatcher

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta"
flycatcher_genome="FicAlb.1.5.genome.ncbi.fasta"

SatsumaSynteny -q $genome -t $flycatcher_genome -n 19 -o $genome.flycatcher.satsuma > $genome.flycatcher.satsuma.stdout
```


## Optical maps


Map the southern optical map to the northern genome

```
perl Solve3.2.2_08222018/Pipeline/08222018/fa2cmap_multi_color.pl -i ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -e BspQI 1 -o ww_pacbio_bionano

python Solve3.2.2_08222018/Pipeline/08222018/runCharacterize.py -t Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -r ww_pacbio_bionano/ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt_BSPQI_0kb_0labels.cmap -q exp_refineFinal1_contigs.cmap -n19 -a optArguments_human.xml

python Solve3.2.2_08222018/Pipeline/08222018/runSV.py -t Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -r ww_pacbio_bionano/ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt_BSPQI_0kb_0labels.cmap -q exp_refineFinal1_contigs.cmap -T 19 -j 19 -a optArguments_human.xml -e alignref/exp_refineFinal1_contigs.err -o southern_bionano_vs_ref_genome_sv
 
cd southern_bionano_vs_ref_genome_sv
python smap_to_vcf_v2.py -s exp_refineFinal1_contigs.smap -o southern_bionano_vs_northern_ref
 
```

Map the northern optical map to the southern genome

```
perl Solve3.2.2_08222018/Pipeline/08222018/fa2cmap_multi_color.pl -i ww_southern_10x_bionano.filt.new.fasta -e BspQI 1 -o ww_southern_10x_bionano.filt  


python Solve3.2.2_08222018/Pipeline/08222018/runCharacterize.py -t Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -r ww_southern_10x_bionano.filt/ww_southern_10x_bionano.filt_BSPQI_0kb_0labels.cmap -q exp_refineFinal1_contigs.cmap -n19 -a optArguments_human.xml

python Solve3.2.2_08222018/Pipeline/08222018/runSV.py -t Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -r ww_southern_10x_bionano.filt/ww_southern_10x_bionano.filt_BSPQI_0kb_0labels.cmap -q exp_refineFinal1_contigs.cmap -T 19 -j 19 -a optArguments_human.xml -e alignref/exp_refineFinal1_contigs.err -o northern_bionano_vs_southern_genome_sv
 
cd northern_bionano_vs_southern_genome_sv

python smap_to_vcf_v2.py -s exp_refineFinal1_contigs.smap -o northern_bionano_vs_southern_ref

```

Also hybridize the southern optical map to the northern genome in order to visualize the differences close to breakpoints. This time don't perform any cuts in the reference
assembly

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta"

perl Solve3.2.2_08222018/HybridScaffold/08222018/hybridScaffold.pl -f -B 2 -N 1 -r Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -n $genome -b exp_refineFinal1_contigs.cmap -c hybridScaffold_config_solve_aggressive.xml -o $genome.bionano.southern.solve_config.aggressive.no_cut.output

```


## Map 10x chromium data to detect structural variants

Here we will use the `longranger wgs` pipeline to check for structural variants in the divergent regions, in particular inversion breakpoints

Map the southern and northern linked reads to the northern assembly

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta"
genome_prefix=$(echo $genome | sed 's/.fasta//')

longranger mkref $genome

java -jar $PICARD_HOME/picard.jar CreateSequenceDictionary R=./refdata-$genome_prefix/fasta/genome.fa O=./refdata-$genome_prefix/fasta/genome.dict

#Set sex to female as it has two copies of the sex chromosome (Z), as would be the case in the female of mammals
longranger wgs --reference=refdata-$genome_prefix --id=P6352_101_uppmax_longranger_wgs --fastqs=./northern_barcoded_reads --sex=female
longranger wgs --reference=refdata-$genome_prefix --id=P6352_101_uppmax_longranger_wgs_southern_genome --fastqs=./southern_barcoded_reads --sex=female
```

Map southern and northern linked reads to the southern assembly. In this case we will extract the 500 largest scaffolds from the assembly

```
genome="ww_southern_10x_bionano.filt.500_largest_scaffolds.fasta"
genome_prefix=$(echo $genome | sed 's/.fasta//')

longranger mkref $genome

java -jar $PICARD_HOME/picard.jar CreateSequenceDictionary R=./refdata-$genome_prefix/fasta/genome.fa O=./refdata-$genome_prefix/fasta/genome.dict

#Set sex to female as it has two copies of the sex chromosome (Z), as would be the case in the female of mammals
longranger wgs --reference=refdata-$genome_prefix --id=P6352_101_uppmax_longranger_wgs_southern_genome --fastqs=./ --sex=female
longranger wgs --reference=refdata-$genome_prefix --id=P6352_102_uppmax_longranger_wgs_southern_genome --fastqs=Sample_P6352_102/ --sex=female
```


## Calculate coverage of linked reads 

Map linked reads of the northern and southern sample to the northern genome

```
bwa index ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta

bwa mem ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -p northern_barcoded.fastq.gz -t 19 -C | samtools view -bS - > P6352_101.bam
samtools sort -@ 20 P6352_101.bam > P6352_101.sorted.bam
samtools index P6352_101.sorted.bam

bwa mem ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -p southern_barcoded.fastq.gz -t 19 -C | samtools view -bS - > P6352_102.bam
samtools sort -@ 20 P6352_102.bam > P6352_102.sorted.bam
samtools index P6352_102.sorted.bam
```

Map linked reads of the northern and southern sample to the southern genome

```
bwa index ww_southern_10x_bionano.filt.new.fasta

bwa mem ww_southern_10x_bionano.filt.new.fasta -p northern_barcoded.fastq.gz -t 19 -C | samtools view -bS - > P6352_101.southern_genome.bam
samtools sort -@ 20 P6352_101.southern_genome.bam > P6352_101.southern_genome.sorted.bam
samtools index P6352_101.southern_genome.sorted.bam

bwa mem ww_southern_10x_bionano.filt.new.fasta -p southern_barcoded.fastq.gz -t 19 -C | samtools view -bS - > P6352_102.southern_genome.bam
samtools sort -@ 20 P6352_102.southern_genome.bam > P6352_102.southern_genome.sorted.bam
samtools index P6352_102.southern_genome.sorted.bam
```

Extract alignments from the scaffolds in the divergent region on chromosome 5 (Scaffold68 in the northern genome and Scaffold36 in the southern genome)
and in the divergent region on chromosome 1 in the southern genome (Scaffold11).

```
#Northern genome
samtools view -h P6352_101.sorted.bam Scaffold68 | samtools sort -tBX - > P6352_101.Scaffold68.bc.bam
samtools view -h P6352_102.sorted.bam Scaffold68 | samtools sort -tBX - > P6352_102.Scaffold68.bc.bam

#Southern genome 
samtools view -h P6352_101.southern_genome.sorted.bam Scaffold11 | samtools sort -tBX - > P6352_101.southern_genome.Scaffold11.bc.bam
samtools view -h P6352_102.southern_genome.sorted.bam Scaffold11 | samtools sort -tBX - > P6352_102.southern_genome.Scaffold11.bc.bam

samtools view -h P6352_101.southern_genome.sorted.bam Scaffold36 | samtools sort -tBX - > P6352_101.southern_genome.Scaffold36.bc.bam
samtools view -h P6352_102.southern_genome.sorted.bam Scaffold36 | samtools sort -tBX - > P6352_102.southern_genome.Scaffold36.bc.bam
```

Extract molecule positions (bed file) for each sample and region using `tigmint`

```
molecule_size=10000

bam_prefixes=$(ls *bc.bam | sed 's/\.bam//')

for prefix in $bam_prefixes
do
tigmint-molecule $prefix.bam -s $molecule_size -q 1  | sort -k2,2n -k3,3n > $prefix.molecules.bed
done
```

Quantify the number of molecules overlapping non-overlapping 1 kb windows using `BEDTools`

```
cat ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.1kb_windows.bed | grep "Scaffold68\b" > Scaffold68.1kb_windows.bed
cat ww_southern_10x_bionano.filt.new.fasta.1kb_windows.bed | grep "Scaffold11\b" > Scaffold11.1kb_windows.bed
cat ww_southern_10x_bionano.filt.new.fasta.1kb_windows.bed | grep "Scaffold36\b" > Scaffold36.1kb_windows.bed


bedtools coverage -a Scaffold68.1kb_windows.bed -b P6352_101.Scaffold68.bc.molecules.bed > P6352_101.Scaffold68.molecules.1kb_coverage.bed
bedtools coverage -a Scaffold68.1kb_windows.bed -b P6352_102.Scaffold68.bc.molecules.bed > P6352_102.Scaffold68.molecules.1kb_coverage.bed

bedtools coverage -a Scaffold11.1kb_windows.bed -b P6352_101.southern_genome.Scaffold11.bc.molecules.bed > P6352_101.Scaffold11.molecules.1kb_coverage.bed
bedtools coverage -a Scaffold11.1kb_windows.bed -b P6352_102.southern_genome.Scaffold11.bc.molecules.bed > P6352_102.Scaffold11.molecules.1kb_coverage.bed

bedtools coverage -a Scaffold36.1kb_windows.bed -b P6352_101.southern_genome.Scaffold36.bc.molecules.bed > P6352_101.Scaffold36.molecules.1kb_coverage.bed
bedtools coverage -a Scaffold36.1kb_windows.bed -b P6352_102.southern_genome.Scaffold36.bc.molecules.bed > P6352_102.Scaffold36.molecules.1kb_coverage.bed
```





