# Resequencing data


## Trim reads

Trim raw reads using `trimmomatic`

```
java -jar trimmomatic.jar PE -threads 16 -basein ${sample}_R1_001.fastq.gz -baseout ${sample}_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
```

## Map reads to assembly

Map reads for each sample using `bwa mem` with a list of samples and read prefixes. Read duplicates will be removed using `picardtools`.

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta"
sample_list="sample_info_phyllocopus_new_genome.txt"

bwa index $genome

cat $sample_list | while read name prefix
do
bwa mem $genome -M -t 19 -R "@RG\tID:${name}\tLB:$name\tSM:$name\tPL:Illumina" ${prefix}_[12]P.fastq.gz | samtools view -bS - > $name.bam
samtools sort -@ 20 $name.bam > $name.sorted.bam
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$name.sorted.bam O=$name.sorted.nodup.bam M=$name.sorted.bam.duplicatedata.txt
samtools index $name.sorted.nodup.bam
cp $name.sorted.nodup.bam* $out_dir
cp $name.sorted.bam.duplicatedata.txt $out_dir
done
```


## Call variants

Call variants using `freebayes`. Use `GNU parallel` to make the variant calling more efficient

```
#create a size-sorted list of scaffolds and remove the mitochondrial scaffold from the list
cat ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta.fai | grep -v MT | sort -k2,2nr | cut -f1 > ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta.scaffolds.list

ls *.bam > bamfiles.list

parallel --jobs 20 "freebayes --bam-list bamfiles.list -f ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta --region {} | gzip -c > ./vcf_files/{}.vcf.gz; echo {} processed" :::: ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta.scaffolds.list

vcf-concat vcf_files/*.vcf.gz | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.raw.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.raw.vcf.gz
```


## Filter variants

Use a combination of `vcflib` and `vcftools` to filter the raw set of variants


First, filter variants based on quality, strand support and read placement

```
vcffilter -f 'QUAL > 30 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0' freebayes_old_reseq_plus_new_samples_pacbio_bionano.raw.vcf.gz | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt1.vcf.gz&
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt1.vcf.gz
```

Next, filter based on minimum coverage and remove sites that are missing a maximum of four genotypes in each population

```
vcffilter -g "DP > 4" freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt1.vcf.gz | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.vcf.gz

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.vcf.gz --keep ~/southern_samples_reseq_new.list --max-missing-count 4 --removed-sites --stdout > removed_sites_southern.txt #5454894 sites
vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.vcf.gz --keep ~/northern_samples_reseq_new.list --max-missing-count 4 --removed-sites --stdout > removed_sites_northern.txt #2328213 sites

cat removed_sites_southern.txt removed_sites_northern.txt | grep -v CHROM | sort -k1,1 -k2,2n | uniq > removed_sites_combined.txt  #5735154 sites

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.vcf.gz --exclude-positions removed_sites_combined.txt --recode --recode-INFO-all --stdout | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz&
```

Remove sites with a mean coverage twice that of the median mean coverage of all sites and monomorphic sites

```
vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz --site-mean-depth --stdout  --not-chr MT > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.mean_depths.out&

#Get median in R
depths=read.delim("freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.mean_depths.out")
median(depths[,3])
[1] 15.1818
#

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz --max-meanDP 30 --removed-sites --stdout > excessive_coverage_positions.out

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz --freq --stdout | sed '1d' | grep ":1" > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.monomorphic_pos.out

cat excessive_coverage_positions.out freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.monomorphic_pos.out | grep -v CHROM | cut -f1-2 | sort -k1,1 -k2,2n | uniq > excessive_coverage_monomorphic_positions.out 

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz --exclude-positions excessive_coverage_monomorphic_positions.out --recode --recode-INFO-all --stdout | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.vcf.gz
```

Decompose variants and remove variants overlapping annotated repeats

```
vcfallelicprimitives -kg freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.vcf.gz | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.vcf.gz

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.vcf.gz --recode --recode-INFO-all --exclude-bed ww_pacbio_bionano_final.repeats.corrected.new.bed --stdout | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.vcf.gz
```


For several analyses later on, such as fst calculations, we will need to fix how a small number of missing genotypes have been recoded in vcftools

```
zcat freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.vcf.gz | perl -ne 'if($_=~/\t\.:/){$_=~s/\t\.:/\t\.\/\.:/g} print $_' | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz
```


## Calculate FST

To calculate FST and extract highly differentiated positions we will use `vcftools`

```
fst_thresholds="0.5 0.6 0.7 0.8 0.9 1"
for fst_threshold in $fst_thresholds
do
export fst_threshold
zcat freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz | vcftools --vcf - --weir-fst-pop ~/southern_samples_reseq_new.list  --weir-fst-pop ~/northern_samples_reseq_new.list --stdout | sed '1d' | perl -lane 'if($F[2]>=$ENV{fst_threshold}){print $F[0],"\t",$F[1]-1,"\t",$F[1],"\t",$F[2]}' > fst_${fst_threshold}_variants.bed
done
```

Also calculate FST for bi-allelic SNPs in non-overlapping 10 kb windows

```
vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz --remove-indels --max-alleles 2 --weir-fst-pop ~/southern_samples_reseq_new.list --weir-fst-pop ~/northern_samples_reseq_new.list --fst-window-size 10000 --stdout > fst_10kb_window_biallelic_SNPs.new.out
```


## Annotate highly differentiated variants in the divergent regions

Extract variants that have a FST of at least 0.7 between homozygotes of northern and southern haplotypes in each of the divergent regions
For this purpose use `vcftools` and `bcftools` to efficiently retrieve variants from each region 

```
mkdir highly_diff_variants_div_region


#Chr1
bcftools view freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz Scaffold156 | bgzip -c > Scaffold156.vcf.gz

bcftools view Scaffold156.vcf.gz | vcftools --vcf - --weir-fst-pop ~/chr1_southern.txt  --weir-fst-pop ~/chr1_northern.txt --stdout | perl -lane 'if($F[2]>=0.7){print $F[0],"\t",$F[1],"\t",$F[2]}' > highly_diff_variants_div_region/chr1_highly_diff_positions.fst.out
cat highly_diff_variants_div_region/chr1_highly_diff_positions.fst.out | cut -f1-2 > highly_diff_variants_div_region/chr1_highly_diff_positions.list

bcftools view Scaffold156.vcf.gz | vcftools --vcf - --positions highly_diff_variants_div_region/chr1_highly_diff_positions.list --recode --recode-INFO-all --stdout > highly_diff_variants_div_region/chr1_highly_differentiated_variants.vcf



#Chr3
bcftools view freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz Scaffold29 | bgzip -c > Scaffold29.vcf.gz

bcftools view Scaffold29.vcf.gz Scaffold29:55740000- | vcftools --vcf - --weir-fst-pop ~/chr3_southern.txt  --weir-fst-pop ~/chr3_northern.txt --stdout | perl -lane 'if($F[2]>=0.7){print $F[0],"\t",$F[1],"\t",$F[2]}' > highly_diff_variants_div_region/chr3_highly_diff_positions.fst.out
cat highly_diff_variants_div_region/chr3_highly_diff_positions.fst.out | cut -f1-2 > highly_diff_variants_div_region/chr3_highly_diff_positions.list

bcftools view Scaffold29.vcf.gz Scaffold29:55740000- | vcftools --vcf - --positions highly_diff_variants_div_region/chr3_highly_diff_positions.list --recode --recode-INFO-all --stdout > highly_diff_variants_div_region/chr3_highly_differentiated_variants.vcf


#Chr5
bcftools view freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz Scaffold68 | bgzip -c > Scaffold68.vcf.gz

bcftools view Scaffold68.vcf.gz | vcftools --vcf - --weir-fst-pop ~/chr5_southern.txt  --weir-fst-pop ~/chr5_northern.txt --stdout | perl -lane 'if($F[2]>=0.7){print $F[0],"\t",$F[1],"\t",$F[2]}' > highly_diff_variants_div_region/chr5_highly_diff_positions.fst.out
cat highly_diff_variants_div_region/chr5_highly_diff_positions.fst.out | cut -f1-2 > highly_diff_variants_div_region/chr5_highly_diff_positions.list

bcftools view Scaffold68.vcf.gz | vcftools --vcf - --positions highly_diff_variants_div_region/chr5_highly_diff_positions.list --recode --recode-INFO-all --stdout > highly_diff_variants_div_region/chr5_highly_differentiated_variants.vcf

```

Use `snpeff` and `snpsift` to annotate the variants. As for annotation we will use manually curated genes from the divergent regions 

```
mkdir ww
cd genomes

mv ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta ww.fa

cd ../ww
mv ww_webapollo_annotations_20200221.gff3.gz genes.gff

java -jar snpEff.jar build -gff3 -v ww 


java -jar snpEff.jar ann ww chr1_highly_differentiated_variants.vcf  > chr1_highly_differentiated_variants.annotated.vcf
mv snpEff_genes.txt chr1_snpEff_genes.txt
mv snpEff_summary.html chr1_snpEff_summary.html

java -jar snpEff.jar ann ww chr3_highly_differentiated_variants.vcf  > chr3_highly_differentiated_variants.annotated.vcf
mv snpEff_genes.txt chr3_snpEff_genes.txt
mv snpEff_summary.html chr3_snpEff_summary.html

java -jar snpEff.jar ann ww chr5_highly_differentiated_variants.vcf  > chr5_highly_differentiated_variants.annotated.vcf
mv snpEff_genes.txt chr5_snpEff_genes.txt
mv snpEff_summary.html chr5_snpEff_summary.html


cat chr1_highly_differentiated_variants.annotated.vcf | snpEff/scripts/vcfEffOnePerLine.pl | java -jar SnpSift.jar extractFields - CHROM POS REF ALT ANN[*].GENE ANN[*].EFFECT ANN[*].IMPACT ANN[*].ERRORS | egrep "MODERATE|HIGH" | uniq > moderate_high_effect_variants.list 
cat chr3_highly_differentiated_variants.annotated.vcf | snpEff/scripts/vcfEffOnePerLine.pl | java -jar SnpSift.jar extractFields - CHROM POS REF ALT ANN[*].GENE ANN[*].EFFECT ANN[*].IMPACT ANN[*].ERRORS | egrep "MODERATE|HIGH" | uniq  >> moderate_high_effect_variants.list 
cat chr5_highly_differentiated_variants.annotated.vcf | snpEff/scripts/vcfEffOnePerLine.pl | java -jar SnpSift.jar extractFields - CHROM POS REF ALT ANN[*].GENE ANN[*].EFFECT ANN[*].IMPACT ANN[*].ERRORS | egrep "MODERATE|HIGH" | uniq  >> moderate_high_effect_variants.list 

```


## Call structural variants in the resequenced samples

First run `delly`. To speed things up we will run delly in parallel using `gnu parallel`

```
ls *.bam | sed 's/.bam//' > samples.list


#DELETIONS
parallel --jobs 20 'delly call -t DEL -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.del.round1.bcf {}.bam' :::: samples.list
delly merge -t DEL -o del.sites.bcf *.del.round1.bcf
parallel --jobs 20 'delly call -t DEL -v del.sites.bcf -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.del.round2.bcf {}.bam' :::: samples.list
bcftools merge -m id -O b -o ww.merged.del.bcf *.del.round2.bcf

#DUPLICATIONS
parallel --jobs 20 'delly call -t DUP -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.dup.round1.bcf {}.bam' :::: samples.list
delly merge -t DUP -o dup.sites.bcf *.dup.round1.bcf
parallel --jobs 20 'delly call -t DUP -v dup.sites.bcf -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.dup.round2.bcf {}.bam' :::: samples.list
bcftools merge -m id -O b -o ww.merged.dup.bcf *.dup.round2.bcf

#INV
parallel --jobs 20 'delly call -t INV -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.inv.round1.bcf {}.bam' :::: samples.list
delly merge -t INV -o inv.sites.bcf *.inv.round1.bcf
parallel --jobs 20 'delly call -t INV -v inv.sites.bcf -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.inv.round2.bcf {}.bam' :::: samples.list
bcftools merge -m id -O b -o ww.merged.inv.bcf *.inv.round2.bcf

#INS
parallel --jobs 20 'delly call -t INS -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.ins.round1.bcf {}.bam' :::: samples.list
delly merge -t INS -o ins.sites.bcf *.ins.round1.bcf
parallel --jobs 20 'delly call -t INS -v ins.sites.bcf -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.ins.round2.bcf {}.bam' :::: samples.list
bcftools merge -m id -O b -o ww.merged.ins.bcf *.ins.round2.bcf

#BND
parallel --jobs 20 'delly call -t BND -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.bnd.round1.bcf {}.bam' :::: samples.list
delly merge -t BND -o bnd.sites.bcf *.bnd.round1.bcf
parallel --jobs 20 'delly call -t BND -v bnd.sites.bcf -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -o {}.bnd.round2.bcf {}.bam' :::: samples.list
bcftools merge -m id -O b -o ww.merged.bnd.bcf *.bnd.round2.bcf

types="bnd del dup ins inv"

for type in $types; do bcftools index ww.merged.$type.bcf; done
```

Filter the raw ouput

```
delly filter -m 1 -t DEL -f germline -o ww.merged.del.filt.bcf ww.merged.del.bcf
delly filter -m 1 -t DUP -f germline -o ww.merged.dup.filt.bcf ww.merged.dup.bcf
delly filter -m 1 -t INV -f germline -o ww.merged.inv.filt.bcf ww.merged.inv.bcf
delly filter -m 1 -t INS -f germline -o ww.merged.ins.filt.bcf ww.merged.ins.bcf
delly filter -m 1 -t BND -f germline -o ww.merged.bnd.filt.bcf ww.merged.bnd.bcf
```

Calculate FST between northern and southern homozygotes in the divergent regions. This will only be done for deletions and insertions, as these are the
only categories containing highly differentiated variants of higher quality

```
#Deletions
bcftools view ww.merged.del.filt.bcf | vcftools --vcf - --weir-fst-pop ~/chr1_southern.txt --weir-fst-pop ~/chr1_northern.txt --stdout | sed '1d' | grep -v "nan" | awk '{if($1=="Scaffold156" && $3>=0.7){print $0}}' > high_diff_del.delly.out
bcftools view ww.merged.del.filt.bcf | vcftools --vcf - --weir-fst-pop ~/chr3_southern.txt --weir-fst-pop ~/chr3_northern.txt --stdout | sed '1d' | grep -v "nan" | awk '{if($1=="Scaffold29" && $3>=0.7){print $0}}' >> high_diff_del.delly.out
bcftools view ww.merged.del.filt.bcf | vcftools --vcf - --weir-fst-pop ~/chr5_southern.txt --weir-fst-pop ~/chr5_northern.txt --stdout | sed '1d' | grep -v "nan" | awk '{if($1=="Scaffold68" && $3>=0.7){print $0}}' >> high_diff_del.delly.out

cat high_diff_del.delly.out  | cut -f1-2 | sort -k1,1 -k2,2n > high_diff_del.delly.pos.list

bcftools view ww.merged.del.filt.bcf | vcftools --vcf - --positions high_diff_del.delly.pos.list --recode --recode-INFO-all --stdout > high_diff_del.delly.vcf

cat high_diff_del.delly.vcf | grep -v "#" | perl -lane '$F[7]=~/END=(\d+)/; print $F[0],"\t",$F[1]-1,"\t",$1,"\t","Delly_DEL_",$1-$F[1]+1,"bp"' | sort -k1,1 -k2,2n > high_diff_del.delly.bed


#Insertions
bcftools view ww.merged.ins.filt.bcf | vcftools --vcf - --weir-fst-pop ~/chr1_southern.txt --weir-fst-pop ~/chr1_northern.txt --stdout | sed '1d' | grep -v "nan" | awk '{if($1=="Scaffold156" && $3>=0.7){print $0}}' > high_diff_ins.delly.out
bcftools view ww.merged.ins.filt.bcf | vcftools --vcf - --weir-fst-pop ~/chr3_southern.txt --weir-fst-pop ~/chr3_northern.txt --stdout | sed '1d' | grep -v "nan" | awk '{if($1=="Scaffold29" && $3>=0.7){print $0}}'  >> high_diff_ins.delly.out
bcftools view ww.merged.ins.filt.bcf | vcftools --vcf - --weir-fst-pop ~/chr5_southern.txt --weir-fst-pop ~/chr5_northern.txt --stdout | sed '1d' | grep -v "nan" | awk '{if($1=="Scaffold68" && $3>=0.7){print $0}}'  >> high_diff_ins.delly.out

cat high_diff_ins.delly.out  | cut -f1-2 | sort -k1,1 -k2,2n > high_diff_ins.delly.pos.list

bcftools view ww.merged.ins.filt.bcf | vcftools --vcf - --positions high_diff_ins.delly.pos.list --recode --recode-INFO-all --stdout > high_diff_ins.delly.vcf

cat high_diff_ins.delly.vcf | grep -v "#" | perl -lane '$F[7]=~/END=(\d+);.*INSLEN=(\d+)/; print $F[0],"\t",$F[1]-1,"\t",$1,"\t","Delly_INS_",$2,"bp"' | sort -k1,1 -k2,2n > high_diff_ins.delly.bed
```


Also check SVs with `MANTA`

```
python configManta.py --bam=1A05.sorted.nodup.bam --bam=UK06.sorted.nodup.bam --bam=1N12.sorted.nodup.bam --bam=1M13.sorted.nodup.bam --bam=0G03.sorted.dedup.bam --bam=0G04.sorted.dedup.bam --bam=0G10.sorted.dedup.bam --bam=0J01.sorted.dedup.bam --bam=1K05.sorted.dedup.bam --bam=1K10.sorted.dedup.bam --bam=1L17.sorted.dedup.bam --bam=1L19.sorted.dedup.bam --bam=1L20.sorted.dedup.bam --bam=1M08.sorted.dedup.bam --bam=1O01.sorted.dedup.bam --bam=1O04.sorted.dedup.bam --bam=1O06.sorted.dedup.bam --bam=1P02.sorted.dedup.bam --bam=3K06.sorted.dedup.bam --bam=7A12.sorted.dedup.bam --bam=96A01.sorted.dedup.bam --bam=96B07.sorted.dedup.bam --referenceFasta ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta

MantaWorkflow/runWorkflow.py -j 20 -m local
```

Filter variants based on the "PASS" flag
```
vcftools --gzvcf diploidSV.vcf.gz --remove-filtered-all --recode --recode-INFO-all  --stdout | gzip -c > diploidSV.pass.vcf.gz

```

Calculate FST using `vcftools`and extract SVs with FST>=0.7 between northern and southern homozygotes in each of the regions  


```
zcat diploidSV.pass.vcf.gz | vcftools --vcf - --weir-fst-pop ~/chr1_southern.txt --weir-fst-pop ~/chr1_northern.txt --stdout | sed '1d' | grep Scaffold156 | awk '{if($3>=0.7 && $3!~/nan/){print $0}}' > manta_highly_diff_svs.out
zcat diploidSV.pass.vcf.gz | vcftools --vcf - --weir-fst-pop ~/chr3_southern.txt --weir-fst-pop ~/chr3_northern.txt --stdout | sed '1d' | grep Scaffold29 | awk '{if($3>=0.7 && $3!~/nan/){print $0}}' >> manta_highly_diff_svs.out
zcat diploidSV.pass.vcf.gz | vcftools --vcf - --weir-fst-pop ~/chr5_southern.txt --weir-fst-pop ~/chr5_northern.txt --stdout | sed '1d' | grep Scaffold68 | awk '{if($3>=0.7 && $3!~/nan/){print $0}}' >> manta_highly_diff_svs.out

cat manta_highly_diff_svs.out | cut -f1-2 | sort -k1,1 -k2,2n > manta_highly_diff_svs.pos.list

vcftools --gzvcf diploidSV.pass.vcf.gz --positions manta_highly_diff_svs.pos.list --stdout --recode --recode-INFO-all > manta_highly_diff_variants.vcf

cat manta_highly_diff_variants.vcf | grep -v "#" | perl -lane '$F[7]=~/END=(\d+);SVTYPE=([^;]+);SVLEN=-*(\d+)/; print $F[0],"\t",$F[1]-1,"\t",$1,"\t","Manta_",$2,"_",$3,"bp"' | sort -k1,1 -k2,2n > manta_highly_diff_variants.bed 
```


Check distance to the closest genes between the manta and delly SV annotations

```
cat high_diff_del.delly.bed high_diff_ins.delly.bed manta_highly_diff_variants.bed | sort -k1,1 -k2,2n > highly_diff_svs_delly_manta_merged.bed

bedtools closest -a highly_diff_svs_delly_manta_merged.bed -b ww_webapollo_annotations_20200221.gene_intervals.bed -d > highly_diff_svs_delly_manta_merged.overlap.genes.bed 
```

Also check overlap between SVs and specific gene features. For this purpose, we will extract exon intervals and also "strict intron" intervals, i.e. gene intervals not
overlapping with exons

```
zcat ww_webapollo_annotations_20200221.gff3.gz | perl -lane  '{if($F[2] eq "exon"){print $F[0],"\t",$F[3]-1,"\t",$F[4],"\t",$F[8]}}' | sort -k1,1 -k2,2n > ww_webapollo_annotations_20200221.exons.bed
bedtools subtract -a ww_webapollo_annotations_20200221.gene_intervals.bed -b ww_webapollo_annotations_20200221.exons.bed > ww_webapollo_annotations_20200221.introns.bed

bedtools intersect -a highly_diff_svs_delly_manta_merged.bed -b annotation/ww_webapollo_annotations_20200221.exons.bed -wo | awk '{print $0,"\t","exon (",$7-$6,"bp)"}' | sed 's/_INS_/\tINS\t/;s/_DEL_/\tDEL\t/' > highly_diff_svs_delly_manta_merged.exon_overlap.out
bedtools intersect -a highly_diff_svs_delly_manta_merged.bed -b annotation/ww_webapollo_annotations_20200221.introns.bed -wo | awk '{print $0,"\t","intron (",$7-$6,"bp)"}' | sed 's/_INS_/\tINS\t/;s/_DEL_/\tDEL\t/' > highly_diff_svs_delly_manta_merged.intron_overlap.out
```


## Get coverage for resequenced samples

Get coverage in 10 kb windows for the resequencing samples using `BEDTools`

```
bed_file="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.10kb_windows.bed"
bedtools multicov -p -q 1 -bed $bed_file -bams 1A05.sorted.nodup.bam UK06.sorted.nodup.bam 0G03.sorted.dedup.bam 0G04.sorted.dedup.bam 0G10.sorted.dedup.bam 0J01.sorted.dedup.bam 1P02.sorted.dedup.bam 3K06.sorted.dedup.bam 7A12.sorted.dedup.bam 96A01.sorted.dedup.bam 96B07.sorted.dedup.bam 1M13.sorted.nodup.bam 1N12.sorted.nodup.bam 1K05.sorted.dedup.bam 1K10.sorted.dedup.bam 1L17.sorted.dedup.bam 1L19.sorted.dedup.bam 1L20.sorted.dedup.bam 1M08.sorted.dedup.bam 1O01.sorted.dedup.bam 1O04.sorted.dedup.bam 1O06.sorted.dedup.bam > ww_read_counts.10kb.windows.out 
```


## Get divergent region haplotypes for resequenced samples

To determine the divergent region genotypes (southern/northern haplotypes) for each resequenced sample we extracted genotypes for SNPs present on the SNP array
in Lundberg et al. 2017. 

First map the SNP array probes (final_probes.fasta) to the northern assembly using `gmap`

```
gmap_build -D ./ -d ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta

gmap -d ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt -D ./ -t5 -n1 -f coords final_probes.fasta > final_probes.new_genome.gmap.coords.out
```

Extract the position of the targeted SNPs from the alignments of the probes e

```
cat final_probes.new_genome.gmap.coords.out | extract_ref_pos_snps_array_gmap_coords.pl > final_probes.new_genome.gmap.coords.pos.info.out
cat final_probes.new_genome.gmap.coords.pos.info.out | grep -v "#" | cut -f3 > mapped_snps.txt 
```

Filter the original SNP array data file using `plink` as in Lundberg et al 2017 and only select SNPs with probes that map in the new genome

```
plink --file 2014-01_All_12_Plates_mod  --autosome-num 40 --maf 0.01 --geno 0.05 --mind 0.05 --exclude loci_to_remove.txt --extract mapped_snps.txt --remove ind_to_remove.txt --recode --out trimmed_new --allow-extra-chr --update-map /snps_plink_update_position.txt
```


Adjust the major allele to the reference allele, flip genotypes of SNPs found on the reverse strand and convert the data into a vcf file

```
plink --file trimmed_new --recode --allow-extra-chr --autosome-num 40 --a2-allele final_probes.new_genome.gmap.coords.pos.info.out 4 3 '#' --out trimmed_new_ref_flipped

cat trimmed_new_ref_flipped.log | perl -ne 'if($_=~/variant (\S+)\.$/){print "$1\n"}' > snps_to_flip.list

plink --file trimmed_new --recode --flip snps_to_flip.list --allow-extra-chr --autosome-num 40 --a2-allele final_probes.new_genome.gmap.coords.pos.info.out 4 3 '#' --out trimmed_new_ref_flipped_2

plink --file trimmed_new_ref_flipped_2 --allow-extra-chr --autosome-num 40 --out snparray_data_filt --update-chr final_probes.new_genome.gmap.coords.pos.info.out 1 3 '#' --update-map final_probes.new_genome.gmap.coords.pos.info.out 2 3 '#' --make-bed

plink --bfile snparray_data_filt --recode vcf --allow-extra-chr --out snparray_data_filt

cat snparray_data_filt.vcf | perl -ne 'if($_=~/CHROM/){$_=~s/\d+\_//g}; print $_' > snparray_data_filt.new.vcf
```

Genotype the resequenced samples and the two linked read samples for these variants using `freebayes` and use `vcflib` and `vcftools` to decompose and select the correct variants.
We will also remove eight SNPs that have more than one alternative allele and convert the data into plink format.

```
cat snparray_data_filt.new.vcf | grep -v "#" | perl -lane 'if($F[0] eq "Scaffold156" || $F[0] eq "Scaffold29" && $F[1]>55759500 || $F[0] eq "Scaffold68"){print $F[0],"\t",$F[1]-1,"\t",$F[1]}' > snparray_data_filt.new.pos.div.region.bed 

cat snparray_data_filt.new.pos.div.region.bed | cut -f1,3 > snparray_data_filt.new.pos.div.region.out

freebayes -f ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta --report-monomorphic --targets snparray_data_filt.new.pos.div.region.bed -L ~/bam_files_reseq_data.list > snp_array_snps_div_regions_reseq_10x_samples.vcf

vcfallelicprimitives -kg snp_array_snps_div_regions_reseq_10x_samples.vcf | vcftools --vcf - --positions snparray_data_filt.new.pos.div.region.out --recode --stdout | sed 's/_102/-102/;s/_101_uppmax_longranger_wgs/-101/' > snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.vcf

cat snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.vcf | cut -f1-5 | grep -v "#" | grep "," | cut -f1-2 > snp_positions_remove.list

vcftools --vcf snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.vcf  --recode --stdout --exclude-positions snp_positions_remove.list | add_snp_id_from_other_vcf.pl --ref_vcf_file snparray_data_filt.new.vcf > snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.vcf

plink --vcf snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.vcf --recode --allow-extra-chr --out snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new
```

Combine the data from the SNParray with that of the resequenced samples

```
cat snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.vcf  | grep -v "#" | cut -f3 > snps_div_region_reseq_chromium.list

plink --bfile snparray_data_filt --extract snps_div_region_reseq_chromium.list --merge snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.ped  snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.map --recode --out merged_data --allow-extra-chr

plink --noweb --file merged_data --make-bed --out merged_bin --allow-extra-chr
```

Use `InvClust` in R to get inversion genotypes in the three divergent regions 


```
plinkdata=read.plink("merged_bin.bed","merged_bin.bim","merged_bin.fam")
plink_genos=plinkdata$genotypes
plink_annot=plinkdata$map

genos=as.data.frame(array(NA,c(nrow(plink_genos),4)))
genos[,1]=paste(plinkdata$fam[,2])

genos2=as.data.frame(array(NA,c(nrow(plink_genos),5)))
genos2[,1]=paste(plinkdata$fam[,2])


#chr1
roi=data.frame(chr="Scaffold156",LBP=1,RBP=12000000,reg= "inv1")
invcall=invClust(roi=roi,wh=1,geno=plink_genos,annot=plink_annot,dim=2)
values=invcall[[2]]$y[,1:2]
samples_plink=plinkdata$fam[,2]
write.table(cbind(samples_plink,values),"invclust_mds_coordinates_chrom1.txt",sep="\t",row.names=F,col.names=F,quote=F)

genos[,2]=paste(invGenotypes(invcall))


#chr3
roi=data.frame(chr="Scaffold29",LBP=55E6,RBP=80E6,reg= "inv3")
invcall=invClust(roi=roi,wh=1,geno=plink_genos,annot=plink_annot,dim=2)
values=invcall[[2]]$y[,1:2]
samples_plink=plinkdata$fam[,2]
write.table(cbind(samples_plink,values),"invclust_mds_coordinates_chrom3.txt",sep="\t",row.names=F,col.names=F,quote=F)

genos[,3]=paste(invGenotypes(invcall))

#chr5
roi=data.frame(chr="Scaffold68",LBP=1,RBP=6000000,reg= "inv5")
invcall=invClust(roi=roi,wh=1,geno=plink_genos,annot=plink_annot,dim=2)
values=invcall[[2]]$y[,1:2]
samples_plink=plinkdata$fam[,2]
write.table(cbind(samples_plink,values),"invclust_mds_coordinates_chrom5.txt",sep="\t",row.names=F,col.names=F,quote=F)


genos[,4]=paste(invGenotypes(invcall))

write.table(genos,file="invclust_genotypes_snparray_resequencing_new_samples.txt",sep="\t",row.names=F,col.names=F,quote=F)
```









