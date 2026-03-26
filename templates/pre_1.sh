#!/bin/bash
paste meta_CADD_head p > header_meta

#cp header.meta meta_CADDALL.txt
#cp header.meta meta_CADD15.txt
#cp header.meta meta_CADD20.txt
(zgrep '^#' f5.vcf.gz && zgrep -v '^#' f5.vcf.gz | sort -u -k1,1 -k2,2 -k4,4 -k5,5 -k6,6) | bgzip > f5_dedup.vcf.gz
bcftools view -G f5_dedup.vcf.gz --threads $task.cpus -Ov -o p1.vcf

grep -v '#' p1.vcf >f6
cut -f 1-8 f6 >p1
cut -f 1-2,4-5 p1 >c1
cut -f 1-2,4-5 p1 | sed 's/\t/\_/g' >c1a
#sed '/#/d' p1.vcf >p1


##variant info
#cut -f 1-2,4-5 p1 > c1
#cut -f 1-2,4-5 p1 | sed 's/\t/\_/g' >c1a

#align the order of alt allele as appears in c1

cut -f 4 c1 >alt
cut -f 8 p1 >c
#cut -f 3 -d';' c | sed 's/CSQ\=//g' | \
#    sed 's/-A/a/g' |\
#    sed 's/-C/c/g' |\
#    sed 's/-G/g/g' |\
#    sed 's/-T/t/g' >csq
awk -F";" '{for (i=1;i<=NF;i++) if ($i ~/CSQ\=/) print$i}' p1 | sed 's/CSQ\=//g' | sed 's/-A/a/g' | sed 's/-C/c/g' | sed 's/-G/g/g' | sed 's/-T/t/g'  > csq

paste p1 csq | awk '$5 !~/,/' | cut -f 1-7,9 >p1_s
awk '$5 ~/,/' p1 >p1_m

paste p1 csq | awk '$5 ~/,/' > p1_order_1
paste p1 csq | awk '$5 ~/,/' |while read i; do echo $i | cut -f 5 -d' ' | sed 's/\,/\n/g' > j; echo '*' > k; echo $i | cut -f 9 -d' ' | sed 's/\,/\n/g' >> k ; cat j |while read l; do grep -w "$l" k; done > x1; c1=$(wc -l x1 | cut -f 1 -d' '); c2=$(wc -l j | cut -f 1 -d' '); if [ "$c1" -eq "$c2" ]; then paste -sd',' x1; else paste -sd',' j; fi; done > order

paste p1_m order |awk '$9 !~/\|/' |cut -f 5 > alt_re
paste p1_m order |awk '$9 !~/\|/' |cut -f 1-8 > p1_re
#paste p1_m order |awk '$9 !~/,/' |cut -f 8 | cut -f 3 -d';' | cut -f 1 > csq_re
paste p1_m order |awk '$9 !~/\|/' |cut -f 8 | awk -F";" '{for (i=1;i<=NF;i++) if ($i ~/CSQ\=/) print$i}' |sed 's/CSQ\=//g' > csq_re
paste p1_m order |awk '$9 ~/\|/' |cut -f 1-7,9 > p1_1

##repeat
paste alt_re csq_re |sed 's/CSQ\=//g' |while read i; do echo $i | cut -f 1 -d' ' | sed 's/\,/\n/g' | awk '{if (length($1)==1) print"--"; else print$i}' > j; echo '*' >k; echo $i | cut -f 2 -d' ' | sed 's/\,/\n/g' >> k ; cat j |while read l; do m=${l:1};grep -w "$m" k; done > x1 ; c1=$(wc -l x1 | cut -f 1 -d' '); c2=$(wc -l j | cut -f 1 -d' '); if [ "$c1" -eq "$c2" ]; then paste -sd',' x1; else paste -sd',' j; fi; done > order_re



paste p1_re order_re |awk '$9 ~/,/' |cut -f 1-7,9 > p1_2

cat p1_s p1_1 p1_2 |\
        sort -k1,1 -k2,2n |\
            awk -F"\t" '{print$1"_"$2"_"$4"_"$5,$6,$7,$8}'  >p1_order
awk 'NR==FNR{a[$1]=$0; next} {print a[$1]}' p1_order c1a >p1_u


#awk 'NR==FNR{x++} END{ if(x!=FNR){print"mismatch ERROR on ma ordering"} }' p1 p1_u
##tba: n of alleles? suppose n.alt=10 atm
cut -f 4 -d' ' p1_u|awk -F"," '{OFS=FS}{for (i=1;i<=NF;i++) if ($i ~/\*/) $i="*|*|*|||||||"}1' > c_u

#cut -f 8 p1_u >c_u


##allele funtional consequence
cut -f 2 -d'|' c_u  >c2
cut -f 1-8 f6 >> f61.vcf
##gene with ensemblID; Note: there are 806 x-genes crossing chunks

#cut -f 3-4 -d'|' c_u|sed 's/|/_/g' >c3
#sort -u c3 > gene.lst
bedtools intersect \
        -wao \
        -a p1.vcf \
        -b p50.bed IBD.bed |\
        cut -f 1-5,13 >p1.bed
datamash -g 1,2,3,4,5 collapse 6 <p1.bed |\
    cut -f 6 >c3

#cut -f 3-4 -d'|' c_u|sed 's/|/_/g' >c3
perl -ne 'print join("\n", split(/\,/,$_));print("\n")' c3 |sort -u |grep -E 'locus|ENSG'>gene.lst
#perl -ne 'print join("\n", split(/\,/,$_));print("\n")' c3 |sort -u |grep -E 'locus'>gene.lst

##AF
cut -f 3 -d';' c_u | awk -F"|" '{OFS="\t"}{if ($5>0) print$6,$15,$24,$33,$42,$51,$60,$69,$78,$87; else print$8,$17,$26,$35,$44,$53,$62,$71,$80,$89}' >c4
#| sed 's/AF\=//g' >c4a
#awk -F, -v OFS=, 'NR==FNR{if(max<10)max=10;next};
#                           {NF=10}1' c4a{,} | sed 's/\,/\t/g' >c4



##raw_score_all
cut -f 3 -d';' c_u |awk -F"|" '{OFS="\t"}{print$9,$18,$27,$36,$45,$54,$63,$72,$81,$90}' >c5

##phred_score >=15, which set smaller scores as 0
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<1.387112){$i="";}}1' c5 >c5a

##phred_score >=20
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<2.097252){$i="";}}1' c5 >c5b
# Step 1: Extract the 3rd field (semicolon-separated) and select specific subfields
###perl -F';' -ane 'print join("\t", (split(/\|/, $F[2]))[8,17,26,35,44,53,62,71,80,89]), "\n"' c_u > c5

# Step 2: Apply Phred score filters
###perl -F'\t' -ane 'for(@F){$_="" if $_ ne "" && $_ < 1.387112}; print join("\t", @F), "\n"' c5 > c5a
###perl -F'\t' -ane 'for(@F){$_="" if $_ ne "" && $_ < 2.097252}; print join("\t", @F), "\n"' c5 > c5b

##genotype
##zgrep -v '#' f5.vcf.gz | cut -f 10- | awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++) $i=substr($i,1,3)}1' >c6
#zcat f5.vcf.gz|grep -v '#' | cut -f 10- |awk '
#BEGIN { OFS="\t" }
#{
#    for (i = 1; i <= NF; i++) {
#        $i = substr($i, 1, 3)
#    }
#    print
#}' > c6

#rm alt_re order*

#rm k j
