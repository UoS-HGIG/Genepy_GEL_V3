#!/bin/bash
###zcat f5.vcf.gz | awk 'BEGIN { OFS="\t" } !/^#/ { for (i = 10; i <= NF; i++) { $i = substr($i, 1, 3) } print }' > c6
###awk 'BEGIN { OFS="\t" } !/^#/ { for (i = 10; i <= NF; i++) { $i = substr($i, 1, 3) } print }' f5.vcf > c6
> meta_CADDALL.txt
> meta_CADD15.txt
> meta_CADD20.txt
region=$1
echo "$region"
mkdir -p "metafiles15_${region}"
mkdir -p "metafiles20_${region}"
mkdir -p "metafilesALL_${region}"
touch "metafilesALL_${region}/1.txt"
touch "metafiles15_${region}/1.txt"
touch "metafiles20_${region}/1.txt"


cp header_meta meta_CADDALL.txt
cp header_meta meta_CADD15.txt
cp header_meta meta_CADD20.txt

#zgrep -v "#" f5.vcf.gz | perl -F'\t' -lane 'print join("\t", map { substr($_,0,3) } @F[9..$#F])' > c6
#zgrep -v '#' f5.vcf.gz | cut -f 10- | awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++) $i=substr($i,1,3)}1' >c6
zgrep -v '#' f5_dedup.vcf.gz | cut -f10- | perl -F'\t' -lane 'print join("\t", map { substr($_,0,3) } @F)' > c6

##merge;
files=(c1 c2 c3 c4 c5 c6)

# Get the number of lines in the first file
ref_lines=$(wc -l < "${files[0]}")

aligned=true

for file in "${files[@]}"; do
    lines=$(wc -l < "$file")
    if [[ "$lines" -ne "$ref_lines" ]]; then
        aligned=false
        break
    fi
done

if $aligned; then
    echo "All C* files have the same number of lines."
else
    echo "C* Files do not have the same number of lines."
fi
paste c1 c2 c3 c4 c5 c6 >> meta_CADDALL.txt
paste c1 c2 c3 c4 c5a c6 >> meta_CADD15.txt
paste c1 c2 c3 c4 c5b c6 >> meta_CADD20.txt
mv c6 metafilesALL_${region}/.
empty_count=$(awk -F'\t' 'NR > 1 && ($17 == "" || $17 == ".")  && $4 != "*"' meta_CADDALL.txt | wc -l)

if [[ "$empty_count" -eq 0 ]]; then
    echo "MetaCADD_ALL >>>> All rows have values in the SCORE1 (column 17)."
else
    echo "MetaCADD_ALL >>>> Found $empty_count rows with missing values in SCORE1 (column 17)."
fi
##CADD_all#
while read gene;
    do
        
        cp header_meta "metafilesALL_${region}/${gene}_CADDALL.meta"
        grep -a -w "$gene" meta_CADDALL.txt >> "metafilesALL_${region}/${gene}_CADDALL.meta"
        ###bgzip ${gene}_CADDALL.meta
    done < gene.lst

##CADD_15
while read gene;
    do
        cp header_meta "metafiles15_${region}/${gene}_CADD15.meta"
        grep -a -w "$gene" meta_CADD15.txt >> "metafiles15_${region}/${gene}_CADD15.meta"
        ###bgzip ${gene}_CADD15.meta
    done < gene.lst

##CADD_20
while read gene;
    do
        cp header_meta "metafiles20_${region}/${gene}_CADD20.meta"
        grep -a -w "$gene" meta_CADD20.txt >> "metafiles20_${region}/${gene}_CADD20.meta"
        ###bgzip ${gene}_CADD20.meta
   done < gene.lst
