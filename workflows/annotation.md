# Annotation


## Repeatmask genome

First use `repeatmodeler` for *de novo* identification of repeats

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.fasta"

BuildDatabase -name ww_pacbio_bionano_final $genome
RepeatModeler -engine ncbi -pa 7 -database ww_pacbio_bionano_final | tee ww_pacbio_bionano_final_repeatmodeler.out
```

Next use `Repeatmasker` to mask the genome. For this part we will download the repbase repeat release from 20181026 and the Dfam consensus 20171107 library.
Combine this with the output from repeatmodeler (consensi.fa.classified).

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.fasta"

RepeatMasker/util/queryRepeatDatabase.pl -species aves > aves_repeats.fasta 

cat aves_repeats.fasta consensi.fa.classified > combined_aves_rmodeler.lib.fa 

RepeatMasker -pa 7 -s -lib combined_aves_rmodeler.lib.fa -dir $genome

cat ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.fasta.out | perl -lane 'next if($.<4);$,="\t";print $F[4],$F[5]-1,$F[6],$F[9]' | sort -k1,1 -k2,2n > ww_pacbio_bionano_final.repeats.bed
bgzip -c ww_pacbio_bionano_final.repeats.bed > ww_pacbio_bionano_final.repeats.bed.gz
tabix -p bed ww_pacbio_bionano_final.repeats.bed.gz
```

Softmask the genome with `BEDTools`

```
bedtools maskfasta -fi ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.fasta -bed ww_pacbio_bionano_final.repeats.bed -soft -fo ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta
samtools faidx ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta
```


## Map protein data

Download protein data from ensembl for chicken, zebra finch and flycatcher, as well as bird proteins from uniprot with transcript evidence, protein evidence or that are manually curated. Divide
this dataset into batches of 100 sequences using `fasta-splitter` (http://kirill-kryukov.com/study/tools/fasta-splitter/) to facilitate parallelization


```
mkdir protein_parts_2019

protein_files="Gallus_gallus_v6a_proteins_ensembl.fasta Ficedula_albicollis_1.4_proteins_ensembl.fasta ebra_finch_ensembl_proteins.fasta aves_swissprot_feb_2019.fasta aves_uniprot_protein_evidence_feb_2019.fasta aves_uniprot_transcript_evidence_feb_2019.fasta" 

for file in $protein_files
do
fasta-splitter.pl --part-size 100 --measure count $file --out-dir protein_parts_2019
done
```

Map batches of proteins to the softmasked assembly using `exonerate` and speed up the process using `gnu-parallel`

```
parallel -j 20 'exonerate -t ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta  -q ./protein_parts_2019/{}.fasta --model protein2genome --bestn 1 --showtargetgff --fsmmemory 5000 --showvulgar no --showalignment no --showquerygff no --ryo "AveragePercentIdentity: %pi\n" --softmasktarget > ./exonerate_gff/{}.gff' :::: protein_parts.list 

cat exonerate_gff/aves_*.gff > exonerate_aves_uniprot_combined.gff
cat exonerate_gff/Ficedula_albicollis_1.4_proteins_ensembl.part-*.gff > exonerate_flycatcher_combined.gff
cat exonerate_gff/Gallus_gallus_v6a_proteins_ensembl.part-*.gff > exonerate_chicken_combined.gff
cat exonerate_gff/zebra_finch_ensembl_proteins.part-*.gff > exonerate_zebra_finch_combined.gff

cat exonerate*_combined.gff | egrep -v "Hostname|Average" > exonerate_for_augustus.gff
```

Extract hints for augustus using scripts from `BRAKER`

```
align2hints.pl  --in=exonerate_for_augustus.gff  --out=exonerate_for_augustus.hints.new.gff --prg=exonerate --genome ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta

cat exonerate_for_augustus.hints.new.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl > exonerate_hints.nr.new.gff
```


## RNAseq data

Use `gsnap` to map the RNAseq data to the softmasked assembly. This data has previously been trimmed using `trimgalore` with default settings

```
gmap_build -db=ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm -D ./ ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta

database="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm"

prefixes=$(ls *.gz | cut -b1-29 | sort | uniq)

for prefix in $prefixes
do
gsnap -d $database -D ./ --gunzip --novelsplicing=1 --nofails --format=sam  -t 19 --read-group-name $prefix $prefix*.gz | samtools view -bS - > $prefix.bam 
samtools sort -@ 19 $prefix.bam > $prefix.sorted.bam
samtools index $prefix.sorted.bam
done

ls *.sorted.bam > bamfiles.list
samtools merge -@ 19 -b bamfiles.list merged_ww_rnaseq.bam
samtools index merged_ww_rnaseq.bam
```


Create hints for augustus using auxiliary augustus scripts, `samtools` and `strand_cov` (https://github.com/pmenzel/stranded-coverage)

```
samtools sort -n merged_ww_rnaseq.bam > merged_ww_rnaseq.ns.bam
filterBam --uniq --paired --pairwiseAlignment --in merged_ww_rnaseq.ns.bam --out rnaseq.ns.filt.bam
samtools sort rnaseq.ns.filt.bam  > rnaseq.filt.sorted.bam
bam2hints --intronsonly --in=rnaseq.filt.sorted.bam --out=intron_hints_rnaseq.gff

strand_cov -o rnaseq.filt.sorted rnaseq.filt.sorted.bam

cat rnaseq.filt.sorted.plus.wig | wig2hints.pl width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep  --UCSC=stranded.track --radius=4.5 --pri=4 --strand="+" > exonpart_hints_rnaseq_plus.gff
cat rnaseq.filt.sorted.minus.wig | wig2hints.pl width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep  --UCSC=stranded.track --radius=4.5 --pri=4 --strand="-" > exonpart_hints_rnaseq_minus.gff
cat exonpart_hints_rnaseq_plus.gff exonpart_hints_rnaseq_minus.gff > exonpart_hints_rnaseq_both_strands.gff
```

Combine RNAseq hints and protein hints and divide them into separate files for each scaffold to make gene prediction more efficient

```
cat intron_hints_rnaseq.gff exonpart_hints_rnaseq_both_strands.gff exonerate_hints.nr.new.gff | sort -k1,1 -k4,4n> rnaseq_exonerate_protein.hints.new.gff


mkdir augustus_rnaseq_protein_hints_new

cat scaffolds.list | while read scaffold
do
touch augustus_rnaseq_protein_hints_new/$scaffold.hints.gff
done

cat rnaseq_exonerate_protein.hints.new.gff | perl -lane '$,="\t"; open(OUT,">>augustus_rnaseq_protein_hints_new/$F[0].hints.gff");print OUT @F; close(OUT)'
```


## Generate gene models 

Use `augustus` to predict genes. For the prediction we will use specific-specific parameters that originate from training based on transcript data mapped
to a short-read genome. For details about the training see the augustus training workflow in the scripts folder. 

```
mkdir scaffold_fasta

cat scaffolds.list | while read scaffold
do
samtools faidx ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta $scaffold > scaffold_fasta/$scaffold.fasta
done

mkdir augustus_rnaseq_intron_exonpart_exonerate_new_gffs
parallel -j 20 "augustus --alternatives-from-evidence=true --softmasking=true --UTR=on --gff3=on --hintsfile=augustus_rnaseq_protein_hints_new/{}.hints.gff --species=ww_old_genome_webapollo_1200 --allow_hinted_splicesites=atac --extrinsicCfgFile=augustus_config/extrinsic/extrinsic.M.RM.E.W.P.own.12.cfg scaffold_fasta/{}.fasta > augustus_rnaseq_intron_exonpart_exonerate_new_gffs/{}.gff3" :::: scaffolds.list

cat augustus_rnaseq_intron_exonpart_exonerate_new_gffs/*gff3 | join_aug_pred.pl > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.new.gff3
```


## Assemble transcripts

Also assemble transcripts from the RNAseq data using `Hisat2` and `Stringtie`

```
hisat2-build ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.hisat2

prefixes=$(ls *.gz | cut -b1-29 | sort | uniq)

for prefix in $prefixes
do
hisat2 -x $genome_database -1 ${prefix}_1_val_1.fq.gz -2 ${prefix}_2_val_2.fq.gz -p 19 --rna-strandness RF --dta | samtools view -bS - > $prefix.bam
samtools sort -@ 19 $prefix.bam > $prefix.sorted.bam
samtools index $prefix.sorted.bam
done

ls *.sorted.bam > bamfiles.list
samtools merge -@ 19 -b bamfiles.list merged_ww_rnaseq.hisat2.new.bam
samtools index merged_ww_rnaseq.hisat2.new.bam

stringtie --rf -p 4 merged_ww_rnaseq.hisat2.new.bam > stringtie.new.gtf
gffread -E stringtie.new.gtf -o- > stringtie.new.gff3
gffread -w transcripts.stringtie.new.fasta -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.sm.fasta stringtie.new.gtf
```


## Synteny-transferred models

Transfer chicken annotations to the willow warbler genome by first aligning the warbler and chicken genome with `satsuma` and then using `kraken`

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.fasta"
chicken_genome="Gallus_gallus.GRCg6a.dna.toplevel.new.fa"
SatsumaSynteny -q $genome -t $chicken_genome -n 19 -o $genome.chicken_new.satsuma

RunKraken -c kraken_chicken_ww_pacbio_bionano.config -s Gallus_gallus.GRCg6a.95.gtf -S chicken -T warbler -o kraken_chicken.gtf
CleanKrakenFiles -i kraken_chicken.gtf > kraken_chicken.filt.gtf
```


## Functional annotation of gene models

Collect info from the augustus predictions

```
cat augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.new.gff3 | collect_info_from_augustus_gff.pl > summary_augustus.out
```

Extract amino acid sequence for the longest isoform and perform reciprocal blast searches against chicken protein (longest translations)

```
getAnnoFasta.pl augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.new.gff3
mv augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.new3.aa augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.fasta

cat augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.fasta | extract_longest_aug_transcript.pl > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.fasta 

makeblastdb -in Gallus_gallus_v6a_proteins_ensembl.longest_transcript.fasta -out Gallus_gallus_v6a_proteins_ensembl.longest_transcript.fasta -title Gallus_gallus_v6a_proteins_ensembl.longest_transcript.fasta -dbtype prot

blastp -query augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.fasta -db ../Gallus_gallus_v6a_proteins_ensembl.longest_transcript.fasta -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads 7 -evalue 1e-5 > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.chicken_longest_transcript.blast.out
 
makeblastdb -in augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.fasta -out augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.fasta -title augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.fasta -dbtype prot

blastp -query ../Gallus_gallus_v6a_proteins_ensembl.longest_transcript.fasta -db augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.fasta -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads 7 -evalue 1e-5 > Gallus_gallus_v6a_proteins_ensembl.longest_transcript.augustus_longest_transcript.blast.out

parse_reciprocal_blast_data.pl  --query_vs_target augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.chicken_longest_transcript.blast.out --target_vs_query Gallus_gallus_v6a_proteins_ensembl.longest_transcript.augustus_longest_transcript.blast.out > reciprocal_blast_summary.out
```


Intersect the CDS of the augustus annotations with the CDS synteny-transferred kraken annotations

```
cat augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.new.gff3 | awk '{if($3=="CDS"){print $0}}' > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.cds.gff3

cat kraken_chicken.filt.corrected.gtf | awk '{if($3=="CDS"){print $0}}' > kraken_chicken.filt.corrected.cds.gtf

bedtools intersect -wao -s -a augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.cds.gff3 -b kraken_chicken.filt.corrected.cds.gtf >  augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.cds.chicken_kraken.overlap.out

summarize_gene_models_synteny_intersect.pl --intersection_file augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.cds.chicken_kraken.overlap.out > summary_overlap_augustus_kraken.out
```


Also blast against vertebrate swissprot proteins (N=85,546)

```
makeblastdb -in swissprot_vertebrates_20190627.fasta -title swissprot_vertebrates_20190627.fasta  -out swissprot_vertebrates_20190627.fasta -dbtype prot
blastp -query augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.fasta -db swissprot_vertebrates_20190627.fasta -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads 7 -evalue 1e-5 > augustus_longest_transcript_swissprot_vertebrates_blast.out
```

Identify protein domains using `Interproscan`. First split the data into 20 batches using `fasta-splitter`

```
mkdir aug_proteins_split interpro_augustus_proteins

fasta-splitter.pl --part-size 1000 --measure count augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.fasta --out-dir aug_proteins_split

ls aug_proteins_split/ | sed 's/.fasta//' > protein_parts_augustus.list

parallel -j 20 'interproscan.sh -i aug_proteins_split/{}.fasta -dp -pa -appl Pfam,PANTHER,SMART,CDD --goterms --iprlookup --output-dir ./interpro_augustus_proteins' :::: protein_parts_augustus.list 

cat interpro_augustus_proteins/augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.part-*.tsv > interpro_augustus_proteins.combined.tsv

summarize_interpro_output.pl --interpro_output_file interpro_augustus_proteins.combined.tsv > interpro_augustus_proteins.combined.summary.out
```

Summarize the functional annotation and select only gene models that has a hit to a protein or contains a protein domain

```
summarize_annotation_data.pl --gene_info summary_augustus.out --blast_summary_file reciprocal_blast_summary.out --kraken_intersect_file summary_overlap_augustus_kraken.out --swissprot_blast augustus_longest_transcript_swissprot_vertebrates_blast.out --swissprot_info_file swissprot_vertebrates_20190627_info.txt --interpro_summary interpro_augustus_proteins.combined.summary.out > augustus_annotation_combined.info.txt

change_gene_names_aug_models_gff3.pl --aug_gff augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.new.gff3 --name_list old_names_new_names_augustus.txt > ww_augustus_new_names.filt.gff3
```

