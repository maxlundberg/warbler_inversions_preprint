# Northern genome assembly


## Optical map assembly

Bionano molecules were assembled *de novo* in Irysview using the `pipelineCL.py` script. 

```
pipelineCL.py -U -d -T 228 -j 228 -N 4 -i 5 -a optArguments_human.xml -w -t /home/bionano/tools/ -l output -b Molecules.bnx -C clusterArguments.xml
```

The contigs from the first assembly were used to determine noise parameters for a second assembly

```
pipelineCL.py -U -d -T 228 -j 228 -N 4 -i 5 -a optArguments_human.xml -w -t /home/bionano/tools/ -l output -b Molecules.bnx -r BN001_1_1_merged.cmap -y -C clusterArguments.xml
```


## Pacbio long-read assembly

The pacbio data was first assembled *de novo* using `HGAP` and unzipped using `Falcon Unzip` using a coverage of 40x prereads


We used 10x chromium Illumina reads from the same sample to polish the raw assembly further

First, trimmed and barcoded reads were extracted from the raw chromium reads using `longranger basic` with default settings. The barcoded reads were mapped to
the assembly using `bwa mem`, converted into bam files using `samtools` and duplicates were marked with ´picardtools´.

As a reference for mapping, the primary contigs and haplotigs were combined into a single fasta file (ps_036_unzip_40x.primary_and_haplotigs.fasta). 

```
genome="ps_036_unzip_40x.primary_and_haplotigs.fasta"
bwa index $genome
bwa mem $genome -p barcoded.fastq.gz -t 19 -C | samtools view -bS - > $genome.P6352_101.bam
samtools sort -@ 20 $genome.P6352_101.bam > $genome.P6352_101.sorted.bam
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=false I=$genome.P6352_101.sorted.bam O=$genome.P6352_101.sorted.dedup.bam M=$genome.P6352_101.sorted.bam.duplicatedata.txt
samtools index $genome.P6352_101.sorted.dedup.bam
```

Next, `Pilon` was used together with the mapped reads to polish the assembly

```
genome="ps_036_unzip_40x.primary_and_haplotigs.fasta"
bam_file="ps_036_unzip_40x.primary_and_haplotigs.fasta.P6352_101.sorted.dedup.bam"
java -Xmx1000G -jar pilon.jar --threads 20 --genome $genome --frags $bam_file --fix indels --output $genome --changes --diploid
```






## Hybridize Illumina-polished primary contigs with bionano data and create supercontigs

The primary contigs was hybridized with the optical map (bionano) from the same sample to cut missassembled contigs and to identify contigs that overlap with each other.
Before performing further scaffolding of the genome, we will attempt join as many of the overlapping contigs as possible and create supercontigs.

```
genome="ps_036_unzip_40x.primary_and_haplotigs.pilon.primary_contigs.fasta"
perl Solve3.2.2_08222018/HybridScaffold/08222018/hybridScaffold.pl -f -B 2 -N 2 -r Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -n $genome -b exp_refineFinal1_contigs.cmap -c hybridScaffold_config_solve_aggressive.xml -o $genome.bionano.solve_config.aggressive.output
```

Extract the cut primary contigs and rename them

```
cat ps_036_unzip_40x.primary_and_haplotigs.pilon.primary_contigs.fasta.cut.fasta | sed 's/\:/_/g' > ps_036_unzip_40x.primary_and_haplotigs.pilon.primary_contigs.fasta.cut.renamed.fasta
```

Also extract information about the contigs that are found in each of the superscaffolds in bionano hybrid assembly and the estimated gap sizes between

```
cat exp_refineFinal1_contigs_bppAdjust_cmap_ps_036_unzip_40x_primary_and_haplotigs_pilon_primary_contigs_fasta_NGScontigs_HYBRID_SCAFFOLD.agp | sed 's/:/_/g' | get_supercontigs_from_bionano_scaffolding.pl > supercontigs_pacbio_unzipped.list
cat exp_refineFinal1_contigs_bppAdjust_cmap_ps_036_unzip_40x_primary_and_haplotigs_pilon_primary_contigs_fasta_NGScontigs_HYBRID_SCAFFOLD.gap | sed 's/:/_/g' | get_gaps_in_supercontigs_from_bionano_scaffolding.pl supercontigs_pacbio_unzipped.list > supercontigs_pacbio_unzipped.est_gap_sizes.out
```


Extract contigs that will be used for each supercontig. The sequences will be imported in to `GAP5`. 

```
mkdir supercontigs_unzip

genome="ps_036_unzip_40x.primary_and_haplotigs.pilon.primary_contigs.fasta.cut.renamed.fasta"
samtools faidx $genome

cat supercontigs_pacbio_unzipped.list | cut -f1,5 | while read id contig
do
samtools faidx $genome $contig >> supercontigs_unzip/sequences.supercontig.$id.fasta
done

#Use tg_index to be able to import sequences into gap5
cd supercontigs_unzip

fastas=$(ls *.fasta)

for fasta in $fastas
do
tg_index $fasta -o $fasta.database
done
```

Add the supercontigs that could be created in GAP5 to the assembly and remove the contigs that were used to create them. Use scripts from
`kentUtils` to filter the assembly

```
cat supercontigs_unzip/Supercontig*.fasta > supercontigs.fasta

faSomeRecords $genome -exclude contigs_in_supercontigs.txt ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.minus_supercontig_contigs.fasta

cat supercontigs.fasta ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.minus_supercontig_contigs.fasta > ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.fasta
```



## Purge haplotigs


In this step we will remove contigs that could be nested within larger contigs. This approach may remove too many of the repeat-rich scaffolds and also remove contigs
that are preferentially hybridized with the bionano map. As a reference for this analysis, we will therefore use the bionano hybrid assembly including (cut) contigs that have not been
incorporated into superscaffolds


Combine superscaffolds and unmapped cut contigs, and use `minimap2` to map pacbio subreads to the assembly

```
cat exp_refineFinal1_contigs_bppAdjust_cmap_ps_036_unzip_40x_primary_and_haplotigs_pilon_primary_contigs_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta exp_refineFinal1_contigs_bppAdjust_cmap_ps_036_unzip_4
0x_primary_and_haplotigs_pilon_primary_contigs_fasta_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta | sed 's/_obj//;s/:/_/;s/[|]/-/g' > ps_036_unzip_40x.primary_and_haplotigs.pilon.bionano.renamed.fasta
samtools faidx ps_036_unzip_40x.primary_and_haplotigs.pilon.bionano.renamed.fasta

genome="ps_036_unzip_40x.primary_and_haplotigs.pilon.bionano.renamed.fasta"
reads="ps_036_subreads.fastq.gz"
minimap2 -a -x map-pb -t 19 $genome $reads | samtools view -bS -F 256 - > $genome.pacbio_reads.bam
samtools sort -@ 20 $genome.pacbio_reads.bam > $genome.pacbio_reads.sorted.bam
samtools index $genome.pacbio_reads.sorted.bam
```

Run purge haplotigs and remove redundant sequences

```
purge_haplotigs readhist -genome ps_036_unzip_40x.primary_and_haplotigs.pilon.bionano.renamed.fasta -bam ps_036_unzip_40x.primary_and_haplotigs.pilon.bionano.renamed.fasta.pacbio_reads.sorted.bam
purge_haplotigs contigcov -in ps_036_unzip_40x.primary_and_haplotigs.pilon.bionano.renamed.fasta.pacbio_reads.sorted.bam.gencov -low 10 -high 85 -mid 34 -junk 105
purge_haplotigs purge -genome ps_036_unzip_40x.primary_and_haplotigs.pilon.bionano.renamed.fasta -c coverage_stats.csv -dotplots -bam ps_036_unzip_40x.primary_and_haplotigs.pilon.bionano.renamed.fasta.pacbio_reads.sorted.bam

#Find the ID of all contigs not included in supercontigs that are considered repeats or haplotigs
cat curated.reassignments.tsv | grep -v "#" | awk '{if($6!="KEEP" && $1!~/Super/){print $1}}' | sed 's/-/\|/g' > haplotigs_repeats_remove.list

#To this list we will also add three contigs included in short superscaffolds that are classified as haplotigs and map to other superscaffolds to a very high degree
and have a high degree of haploid coverage. Also remove six very short contigs.
cat haplotigs_repeats_remove.list short_contigs.list | sort | uniq > contigs_to_remove.list
faSomeRecords ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.fasta -exclude purge_haplotigs_unzip_bionano/contigs_to_remove.list ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.fasta
```



## Scaffold with 10x chromium data

First make the reference compatible with arcs by only retaining the number in the id

```
cat ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.fasta | perl -e 'my %contigs; while(my $input=<STDIN>){
	if($input=~/^>(\d+)/){$seq_id=$1;if(exists($contigs{$seq_id})){$contigs{$seq_id}++;print ">$seq_id","0",$contigs{$seq_id},"\n"}
	else{$contigs{$seq_id}=1;print ">$seq_id","\n"}}
	else{print $input}}' > 	ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta
```

Map the 10x chromium reads to the assembly using `bwa mem`

```
genome="ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta"

bwa index $genome
bwa mem $genome -p barcoded.fastq.gz -t 19 -C | samtools view -bS - > $genome.P6352_101.bam
samtools sort -@ 20 -n $genome.P6352_101.bam > $genome.P6352_101.ns.bam
```

Use `arcs`and `LINKS v1.8.6` to perform the scaffolding

```
genome="ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta"
echo "$genome.P6352_101.ns.bam" > bamfiles.txt

arcs -f $genome -a bamfiles.txt --dist_est

graph=$genome.scaff_s98_c5_l0_d0_e30000_r0.05_original.gv

python makeTSVfile.py $graph $genome.scaff_s98_c5_l0_d0_e30000_r0.05.tigpair_checkpoint.tsv $genome

touch empty.fof

~/links_v1.8.6/LINKS -f ../$genome -s empty.fof -b $genome.scaff_s98_c5_l0_d0_e30000_r0.05 -l 5 -a 0.3

#Rename the output files
mv *.scaffolds.fa ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.fasta
mv *r0.05.assembly_correspondence.tsv ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.assembly_correspondence.tsv

#Change the scaffold names to facilitate downstream analyses
cat ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.fasta | perl -ne 'if($_=~/^>([^,]+)/){print ">$1\n"}else{print $_}' > ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.fasta 
```



## Scaffold with bionano data

Here we will perform a second round of bionano scaffolding

```
genome="ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.fasta"
perl Solve3.2.2_08222018/HybridScaffold/08222018/hybridScaffold.pl -f -B 2 -N 2 -r Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -n $genome -b exp_refineFinal1_contigs.cmap -c hybridScaffold_config_solve_aggressive.xml -o $genome.bionano.solve_config.aggressive.output
```

There are some scaffolds that overlap each other in the superscaffolds and these could be extracted and joined in `GAP5`. First, rename the cut scaffolds.

```
cat ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.fasta.cut.fasta | sed 's/\:/_/g' > ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.cut.fasta 
```

Extract information about the scaffolds that are found in each of the superscaffolds in bionano hybrid assembly and the estimated gap sizes between them

```
exp_refineFinal1_contigs_bppAdjust_cmap_ps_036_unzip_40x_primary_contigs_pilon_fasta_cut_renamed_with_supercontigs_new_filt2_num_fasta_arcs_l5_r0_3_renamed_fasta_NGScontigs_HYBRID_SCAFFOLD.agp | sed 's/:/_/g' | get_supercontigs_from_bionano_scaffolding.pl > supercontigs_pacbio_unzipped.list
exp_refineFinal1_contigs_bppAdjust_cmap_ps_036_unzip_40x_primary_contigs_pilon_fasta_cut_renamed_with_supercontigs_new_filt2_num_fasta_arcs_l5_r0_3_renamed_fasta_NGScontigs_HYBRID_SCAFFOLD.gap | sed 's/:/_/g' | get_gaps_in_supercontigs_from_bionano_scaffolding.pl supercontigs_pacbio_unzipped.list > supercontigs_pacbio_unzipped.est_gap_sizes.out
```


Extract scaffolds that are found in superscaffolds where there are overlaps

```
mkdir merged_scaffolds_bionano_final

genome="ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.cut.fasta"


cat supercontigs_pacbio_unzipped.list | cut -f1,5 | while read id contig
do
samtools faidx $genome $contig >> merged_scaffolds_bionano_final/sequences.superscaffold.$id.fasta
done

cd merged_scaffolds_bionano_final/

fastas=$(ls *.fasta)

for fasta in $fastas
do
tg_index $fasta -o $fasta.database
done
```


Merge the joined superscaffolds with the rest of the scaffolds

```
cat merged_scaffolds_bionano_final/Superscaffold.*.fasta | sed 's/\:/_/g;s/-//g' > Superscaffolds.pacbio.arcs.bionano.fasta

faSomeRecords bionano_40x_unzip_supercontigs_filt_arcs/$genome -exclude scaffolds_to_remove_pacbio_unzip_arcs_cut.list $genome.minus_superscaffold_contigs.fasta

cat Superscaffolds.pacbio.arcs.bionano.fasta $genome.minus_superscaffold_contigs.fasta > $genome.new.fasta
```


Scaffold the genome once more with bionano data, this time without allowing cuts to the scaffolds

```
genome="ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.cut.fasta.new.fasta"

perl Solve3.2.2_08222018/HybridScaffold/08222018/hybridScaffold.pl -f -B 2 -N 1 -r Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -n $genome -b exp_refineFinal1_contigs.cmap -c hybridScaffold_config_solve_aggressive.xml -o $genome.bionano.solve_config.aggressive.no_cut.output
```


Remove a very short scaffold and also remove two chromosome 1 scaffold where there is little support for cutting. The chromosome 1 scaffolds will be replaced by the original uncut sequence from the
previous scaffolding step (arcs)

```
cat *HYBRID_SCAFFOLD.fasta *HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta > ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.cut.fasta.new.bionano.fasta

faSomeRecords ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.cut.fasta.new.bionano.fasta -exclude scaffolds.to.remove.list ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.cut.fasta.new.bionano.filt.fasta

#extract scaffold27 from the original arcs file
samtools faidx ../../ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.out/ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.fasta scaffold27 > scaffold27.fasta


cat scaffold27.fasta ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.cut.fasta.new.bionano.filt.fasta > ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.cut.fasta.new.bionano.filt.new.fasta 

```

Close gaps in the assembly using `PBJelly`. The pipeline will be run with default settings except for the support stage:

```
Jelly.py support pbjelly.xml -x "--spanOnly --capturedOnly"
```



## Second round of Illumina polishing

We performed a second round of Illumina polishing. However, as the number of sequences in the assembly was below 500, we used the `longranger align` pipeline to align reads
in a barcode-aware fashion 

```
genome="jelly.out.fasta"
genome_prefix=$(echo $genome | sed 's/.fasta$//')

longranger mkref $genome

longranger align --reference=./refdata-$genome_prefix --id=P6352_101_uppmax_longranger_basic --fastqs=./
```

Use `Pilon`to polish the genome

```
genome="jelly.out.fasta"
bam_file="possorted_bam.bam"

java -Xmx1000G -jar pilon.jar --threads 20 --genome $genome --frags $bam_file --fix indels --output $genome --changes --diploid
```



## Synteny check

We aligned the genome against the genomes of zebra finch and chicken downloaded from Ensembl using `SatsumaSynteny`

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.fasta"

chicken_genome="Gallus_gallus.GRCg6a.dna.toplevel.new.fa"
zf_genome="Taeniopygia_guttata.taeGut3.2.4.dna.toplevel_final.fa"

SatsumaSynteny -q $genome -t $chicken_genome -n 19 -o $genome.chicken_new.satsuma > $genome.chicken_new.satsuma.stdout

SatsumaSynteny -q $genome -t $zf_genome -n 19 -o $genome.zf.satsuma > $genome.zf.satsuma.stdout

cp -r $genome.zf.satsuma $outdir
cp $genome.zf.satsuma.stdout $outdir
```

One scaffold was cut as it showed a chromosome fusion not present in the other bird species



## Assembly of mitochondrial genome

The mitochondrial sequence was missing in the original assembly. For this purpose we added the mitochondrion from the Lundberg et al (2017) assembly
to the new assembly, mapped linked reads to it and created a new consensus sequence from the reads

```
genome="ps_036_unzip_40x.primary_contigs.pilon.fasta.cut.renamed.with_supercontigs.new.filt2.num.fasta.arcs.l5.r0.3.renamed.cut.fasta.new.bionano.filt.new.pilon.fasta"

cat $genome ww_mt_final_old_assembly.fasta > $genome.with_old_mtdna.fasta

bwa index $genome.with_old_mtdna.fasta
bwa mem $genome.with_old_mtdna.fasta -p barcoded.fastq.gz -t 19 -C | samtools view -bS - > $genome.with_old_mtdna.fasta.P6352_101.bam
samtools sort -@ 20 $genome.with_old_mtdna.fasta.P6352_101.bam > $genome.with_old_mtdna.fasta.P6352_101.sorted.bam
samtools index $genome.with_old_mtdna.fasta.P6352_101.sorted.bam
```

From the aligned reads call variants on the mitochondrial scaffold. The variants will be filtered based on quality and excluded variants located in high coverage intervals
(12800-13200,14100-15800). 

```
freebayes --ploidy 1 --bam mitochondrial_alignments.bam -f ../ww_mt_final_old_assembly.fasta  > mito_raw.vcf
vcftools --vcf mito_raw.vcf --exclude-bed problem_intervals_mt.bed --stdout --recode --recode-INFO-all  --minQ 30 > mito_filt.vcf

bcftools consensus --fasta-ref ../ww_mt_final_old_assembly.fasta mito_filt.vcf.gz | sed 's/ww_mt_final/MT/' > mt_scaffold_new.fasta
```




# Southern genome assembly


## Linked read *de novo* assembly

Use `Supernova` to assemble 10x chromium reads *de novo*

```
supernova run --id=P6352_102_supernova_version2 --fastqs=Sample_P6352_102/ --indices=SI-GA-C2 --lanes 3
supernova mkoutput --asmdir=P6352_102_supernova_version2/outs/assembly/ --outprefix=P6352_102_supernova_version2 --style=pseudohap
``` 

## Scaffold with bionano data

Bionano molecules were assembled *de novo* using the same workflow as for the northern sample (see above). 

Next, the supernova assembly was hybridized with the bionano contigs

```
genome="P6352_102_supernova_version2.fasta"
 
perl Solve3.2.2_08222018/HybridScaffold/08222018/hybridScaffold.pl -f -B 2 -N 2 -r Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -n $genome -b exp_refineFinal1_contigs.cmap -c hybridScaffold_config_solve_aggressive.xml -o $genome.bionano.solve_config.aggressive.output

```

From the raw hybrid assembly we removed scaffolds that were only comprised of ambiguous nucleotides, were entire duplicates of other scaffolds or contained mostly
contaminant sequences. We also replaced two scaffolds that appeared incorrectly split by the bionano data with the corresponding uncut scaffold in the supernova assembly. 

