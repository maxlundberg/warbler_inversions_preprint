## Call variants

Use `freebayes` to call variants in the four high coverage resequenced samples

```
freebayes-v1.3.1 -f ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -b 1A05.sorted.nodup.bam -b UK06.sorted.nodup.bam -b 1L19.sorted.dedup.bam -b 1M13.sorted.nodup.bam -b DW83.sorted.nodup.bam -v pacbio.vcf
gzip -c pacbio.vcf > pacbio.vcf.gz
```


## Run gIMble

Preprocess the data. The file `bam_symlinks` provides a list of the bam files used for variant calling

```
gIMble preprocess -f ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -v pacbio.vcf.gz -b bam_symlinks -g 2 -m 8 -M 0.75 -t 8
```


Next, set up zarr stores

```
gIMble setup -g chromosome_1.gimble.genomefile -b chromosome_1.gimble.bed -s chromosome_1.gimble.samples -v gimble.vcf.gz -o chromosome_1.64
gIMble setup -g chromosome_3.gimble.genomefile -b chromosome_3.gimble.bed -s chromosome_3.gimble.samples -v gimble.vcf.gz -o chromosome_3.64
gIMble setup -g chromosome_5.gimble.genomefile -b chromosome_5.gimble.bed -s chromosome_5.gimble.samples -v gimble.vcf.gz -o chromosome_5.64
```


Get blocks

```
gIMble blocks -z chromosome_1.64.z -l 64 --force
gIMble blocks -z chromosome_3.64.z -l 64 --force
gIMble blocks -z chromosome_5.64.z -l 64 --force
```


Set up models

```
gIMble model -p 2 -s A,B -n 1,1 -j A,B
gIMble model -p 2 -s A,B -n 1,1 -j A,B -m 'A>B'
gIMble model -p 2 -s A,B -n 1,1 -j A,B -m 'B>A'
```


For each divergent region, optimise each model

```
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -b -t 16,1 -i 1000
```


Simulate SI data for each divergent region

```
gIMble simulate -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.sims_1_SI.ini -b 72432 -r 100 -t 8 -l optimised_SI
gIMble simulate -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.sims_3_SI.ini -b 83537 -r 100 -t 8 -l optimised_SI
gIMble simulate -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.sims_5_SI.ini -b 28351 -r 100 -t 8 -l optimised_SI
```


For each simulated SI dataset, optimise under an SI and the best fitting (IM) model

```
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
```


For the chromosome 3 region, simulate and optimise IM2 data to obtain 95% CIs

```
gIMble simulate -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.sims_3_SI.ini -b 83537 -r 100 -t 8 -l optimised_IM
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_IM
```