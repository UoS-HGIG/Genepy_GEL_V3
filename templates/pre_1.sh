#!/bin/bash

## pre_1.sh with ONLY the CADD and allele-frequency extraction changed.
## Everything else (dedup, allele re-ordering, gene assignment) is the original.
##
## What is different:
##   * CSQ field positions come from the ##INFO=<ID=CSQ ...> header at run time
##     (the old code derived gS/gE/eS/eE with awk rules that kept overwriting,
##     so gS landed on the LAST gnomADg_AF* field and the gE range reached into
##     gnomADg_AC / gnomADg_AN -- allele COUNTS competing to be a frequency)
##   * values are read by splitting each CSQ block on "|", so nothing depends
##     on a fixed number of VEP fields
##   * max AF = maximum over ALL gnomADg ancestral-group AF fields; when that
##     is 0 or absent, the maximum over the gnomADe fields is used instead
##   * an allele with no AF in either source records 0; columns past the last
##     ALT allele stay empty

## ---- settings for the new extraction ---------------------------------------
CADD_FIELD=${CADD_FIELD:-CADD_RAW}      # exact CSQ field name of the raw score
GNOMADG_PREFIX=${GNOMADG_PREFIX:-gnomADg}
GNOMADE_PREFIX=${GNOMADE_PREFIX:-gnomADe}
NMAX=${NMAX:-10}                        # number of AF / CADD columns
AF_MISSING=${AF_MISSING:-0}             # allele exists but neither source has an AF
AF_PAD=${AF_PAD:-}                      # column with no allele at all

paste meta_CADD_head p > header_meta

#cp header.meta meta_CADDALL.txt
#cp header.meta meta_CADD15.txt
#cp header.meta meta_CADD20.txt
(zgrep '^#' f5.vcf.gz && zgrep -v '^#' f5.vcf.gz | sort -u -k1,1 -k2,2 -k4,4 -k5,5 -k6,6) | bgzip > f5_dedup.vcf.gz
bcftools view -G f5_dedup.vcf.gz --threads ${THREADS:-1} -Ov -o p1.vcf

grep -v '#' p1.vcf >f6
cut -f 1-8 f6 >p1
cut -f 1-2,4-5 p1 >c1
cut -f 1-2,4-5 p1 | sed 's/\t/\_/g' >c1a

##variant info

#align the order of alt allele as appears in c1

cut -f 4 c1 >alt
cut -f 8 p1 >c
awk -F";" '{for (i=1;i<=NF;i++) if ($i ~/CSQ\=/) print$i}' p1 | sed 's/CSQ\=//g' | sed 's/-A/a/g' | sed 's/-C/c/g' | sed 's/-G/g/g' | sed 's/-T/t/g'  > csq

paste p1 csq | awk '$5 !~/,/' | cut -f 1-7,9 >p1_s
awk '$5 ~/,/' p1 >p1_m

paste p1 csq | awk '$5 ~/,/' > p1_order_1
paste p1 csq | awk '$5 ~/,/' |while read i; do echo $i | cut -f 5 -d' ' | sed 's/\,/\n/g' > j; echo '*' > k; echo $i | cut -f 9 -d' ' | sed 's/\,/\n/g' >> k ; cat j |while read l; do grep -w "^$l" k; done > x1; c1=$(wc -l x1 | cut -f 1 -d' '); c2=$(wc -l j | cut -f 1 -d' '); if [ "$c1" -eq "$c2" ]; then paste -sd',' x1; else paste -sd',' j; fi; done > order

paste p1_m order |awk '$9 !~/\|/' |cut -f 5 > alt_re
paste p1_m order |awk '$9 !~/\|/' |cut -f 1-8 > p1_re
paste p1_m order |awk '$9 !~/\|/' |cut -f 8 | awk -F";" '{for (i=1;i<=NF;i++) if ($i ~/CSQ\=/) print$i}' |sed 's/CSQ\=//g' > csq_re
paste p1_m order |awk '$9 ~/\|/' |cut -f 1-7,9 > p1_1

##repeat
paste alt_re csq_re |sed 's/CSQ\=//g' |while read i; do echo $i | cut -f 1 -d' ' | sed 's/\,/\n/g' | awk '{if (length($1)==1) print"--"; else print$i}' > j; echo '*' >k; echo $i | cut -f 2 -d' ' | sed 's/\,/\n/g' >> k ; cat j |while read l; do m=${l:1};grep -w "^$m" k; done > x1 ; c1=$(wc -l x1 | cut -f 1 -d' '); c2=$(wc -l j | cut -f 1 -d' '); if [ "$c1" -eq "$c2" ]; then paste -sd',' x1; else paste -sd',' j; fi; done > order_re

paste p1_re order_re |awk '$9 ~/,/' |cut -f 1-7,9 > p1_2

cat p1_s p1_1 p1_2 |\
        sort -k1,1 -k2,2n |\
            awk -F"\t" '{print$1"_"$2"_"$4"_"$5,$6,$7,$8}'  >p1_order
awk 'NR==FNR{a[$1]=$0; next} {print a[$1]}' p1_order c1a >p1_u


###############################################################################
## CSQ FIELD INDEX -- derived from the header instead of the gS/gE/eS/eE rules
###############################################################################
csq_f=$(zgrep -m1 "^##INFO=<ID=CSQ" f5_dedup.vcf.gz | sed -E 's/.*Format: *//; s/">.*//; s/"$//')
if [ -z "${csq_f}" ]; then
    echo "ERROR: no '##INFO=<ID=CSQ ... Format: ...' line in f5_dedup.vcf.gz" >&2
    exit 1
fi

## NF_CSQ = fields per block, CADD_I = position of the raw score,
## GG_LIST / GE_LIST = comma separated positions of the gnomADg / gnomADe
## frequency fields. A field qualifies only if its name contains AF and
## contains neither AC nor AN, so allele counts cannot be read as frequencies.
read -r NF_CSQ CADD_I GG_LIST GE_LIST < <(
  echo "${csq_f}" | awk -F'|' \
      -v gg="${GNOMADG_PREFIX}" -v ge="${GNOMADE_PREFIX}" -v cadd="${CADD_FIELD}" '
  {
      ci = 0; gl = ""; el = ""
      for (i = 1; i <= NF; i++) {
          if ($i == cadd) ci = i
          if ($i !~ /AF/ || $i ~ /AC/ || $i ~ /AN/) continue
          if (index($i, gg) == 1) gl = gl (gl == "" ? "" : ",") i
          if (index($i, ge) == 1) el = el (el == "" ? "" : ",") i
      }
      print NF, ci, (gl == "" ? "NA" : gl), (el == "" ? "NA" : el)
  }')

[ "${CADD_I}" -gt 0 ] 2>/dev/null || {
    echo "ERROR: CSQ field '${CADD_FIELD}' not found in the header" >&2; exit 1; }
[ "${GG_LIST}" != "NA" ] || echo "WARNING: no ${GNOMADG_PREFIX} AF fields in the header" >&2
[ "${GE_LIST}" != "NA" ] || echo "WARNING: no ${GNOMADE_PREFIX} AF fields in the header" >&2

## a record of what was picked, for checking after any re-annotation
echo "${csq_f}" | awk -F'|' -v ci="${CADD_I}" -v gl="${GG_LIST}" -v el="${GE_LIST}" '
{
    split(gl, g, ","); for (i in g) R[g[i]] = "gnomADg_max_AF"
    split(el, e, ","); for (i in e) R[e[i]] = (R[e[i]] == "" ? "" : R[e[i]] ",") "gnomADe_max_AF"
    R[ci] = (R[ci] == "" ? "" : R[ci] ",") "CADD_raw"
    print "index\tname\trole"
    for (i = 1; i <= NF; i++) print i "\t" $i "\t" (i in R ? R[i] : "")
}' > csq_field_index.tsv
echo "CSQ fields: ${NF_CSQ} | ${CADD_FIELD} at ${CADD_I} | ${GNOMADG_PREFIX} AF at ${GG_LIST} | ${GNOMADE_PREFIX} AF at ${GE_LIST}" >&2

## placeholder for the * allele, widened to the real number of CSQ fields
star_blk=$(awk -v n="${NF_CSQ}" 'BEGIN{s="*|*|*"; for(i=4;i<=n;i++) s=s"|"; print s}')

##tba: n of alleles? suppose n.alt=10 atm
## NOTE: $i ~/*/ in the original is an unanchored quantifier; /\*/ is meant.
cut -f 4 -d' ' p1_u | awk -F"," -v OFS="," -v sb="${star_blk}" \
    '{for (i=1;i<=NF;i++) if ($i ~ /\*/) $i=sb}1' > c_u

## The extraction reads CSQ block k as ALT allele k, which is what the
## re-ordering above is for. That holds only while there is exactly ONE CSQ
## block per ALT allele; if VEP emitted one block per TRANSCRIPT, or the
paste alt c_u | awk -F"\t" -v nmax="${NMAX}" '
{
    n++
    na = split($1, A, ",")
    nb = split($2, B, ",")
    if (na != nb) { mm++; if (mm <= 5) bad = bad sprintf("  row %d: %d ALT vs %d CSQ blocks\n", NR, na, nb) }
    if (na > nmax) tr++
}
END {
    printf("variants: %d | ALT/CSQ block mismatches: %d | with >%d ALT alleles (extra ignored): %d\n",
           n, mm, nmax, tr) > "/dev/stderr"
    if (mm > 0) {
        printf("WARNING: %d rows where the CSQ block count differs from the ALT count;\n", mm) > "/dev/stderr"
        printf("         on those rows column k is NOT allele k. First few:\n%s", bad) > "/dev/stderr"
    }
}'


##allele funtional consequence
cut -f 2 -d'|' c_u  >c2
echo "##fileformat=VCFv4.2" > f61.vcf
cut -f 1-8 f6 >> f61.vcf
##gene with ensemblID; Note: there are 806 x-genes crossing chunks

bedtools intersect \
        -wao \
        -a f61.vcf \
        -b p50.bed |\
        cut -f 1-5,12 >p1.bed
datamash -g 1,2,3,4,5 collapse 6 <p1.bed |\
    cut -f 6 >c3

perl -ne 'print join("\n", split(/\,/,$_));print("\n")' c3 |sort -u |grep -E 'ENSG'>gene.lst


###############################################################################
## AF -- one column per ALT allele, positions taken from the field index
##
## MAXIMUM AF over all gnomADg ancestral groups; when that maximum is 0 or
## absent, the maximum over the gnomADe groups is used instead. An allele with
## no AF in either source records ${AF_MISSING}; a column with no allele at all
## stays "${AF_PAD}".
###############################################################################
awk -F"," -v OFS="\t" -v gl="${GG_LIST}" -v el="${GE_LIST}" \
    -v nmax="${NMAX}" -v miss="${AF_MISSING}" -v pad="${AF_PAD}" '
BEGIN { NUMRE = "^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$" }
function isnum(s) { return (s != "" && s != "." && s ~ NUMRE) }
# maximum over the positions listed in lst for the block already split into F
function blockmax(lst,   n, a, i, v, best) {
    if (lst == "NA") return ""
    n = split(lst, a, ",");  best = ""
    for (i = 1; i <= n; i++) {
        v = F[a[i]]
        if (isnum(v) && (best == "" || v + 0 > best + 0)) best = v
    }
    return best
}
{
    line = ""
    for (k = 1; k <= nmax; k++) {
        val = pad
        if (k <= NF) {
            split($k, F, "|")
            g = blockmax(gl)
            if (isnum(g) && g + 0 > 0) val = g
            else {
                e = blockmax(el)
                val = isnum(e) ? e : (isnum(g) ? g : miss)
            }
        }
        line = line (k == 1 ? "" : OFS) val
    }
    print line
}' c_u > c4

###############################################################################
## raw_score_all -- CADD raw score at the position given by the field index
###############################################################################
awk -F"," -v OFS="\t" -v ci="${CADD_I}" -v nmax="${NMAX}" '
BEGIN { NUMRE = "^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$" }
function isnum(s) { return (s != "" && s != "." && s ~ NUMRE) }
{
    line = ""
    for (k = 1; k <= nmax; k++) {
        val = ""
        if (k <= NF) { split($k, F, "|"); if (isnum(F[ci])) val = F[ci] }
        line = line (k == 1 ? "" : OFS) val
    }
    print line
}' c_u > c5

##phred_score >=15, which set smaller scores as 0
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<1.387112){$i="";}}1' c5 >c5a

##phred_score >=20
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<2.097252){$i="";}}1' c5 >c5b

##genotype
##zgrep -v '#' f5.vcf.gz | cut -f 10- | awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++) $i=substr($i,1,3)}1' >c6

#rm alt_re order*
#rm k j
