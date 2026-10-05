#!/bin/bash

## pre_1.sh -- VCF (VEP + gnomAD + CADD) -> per-allele annotation columns
##
## Requires VEP run with ALLELE_NUM, CADD_RAW and the gnomAD
## frequency fields are all located by NAME in the ##INFO=<ID=CSQ ...> header at
## run time, so nothing depends on a fixed CSQ layout.
##
## ALLELE FREQUENCY (GNOMAD_v4):
##   * the max is taken over SELECTED ancestry groups only. ami, asj, fin and
##     mid are left out, along with oth/remaining: a founder variant that is
##     common in a bottlenecked population is not thereby common generally, and
##     treating it as such filters out exactly the founder pathogenic alleles
##     worth finding. gnomAD excludes the same groups from its own grpmax.
##   * gnomAD EXOMES first; the genome fields are consulted only when the exome
##     value is zero or absent
##
##
## ONLY THE FIRST ${NMAX} ALT ALLELES ARE KEPT. Any beyond that are ignored --
## the meta files have a fixed ${NMAX} AF columns and ${NMAX} CADD columns, and
## the count of truncated variants is reported to stderr.

## ---- settings ---------------------------------------------------------------
CADD_FIELD=${CADD_FIELD:-CADD_RAW}
GNOMADG_PREFIX=${GNOMADG_PREFIX:-gnomADg}
GNOMADE_PREFIX=${GNOMADE_PREFIX:-gnomADe}
NMAX=${NMAX:-10}                        # ALT alleles kept; the rest are dropped
AF_MISSING=${AF_MISSING:-0}             # allele present, no value in either source
AF_PAD=${AF_PAD:-}                      # column with no allele at all

## ---- frequency selection ----------------------------------------------------
## The AF token is matched as a whole underscore-token, so this works for either
## naming convention without being told which:
##   gnomADe_AFR_AF    VEP --af_gnomade / --af_gnomadg
##   gnomADe_AF_afr    --custom straight off the gnomAD v4 VCF
AF_METRIC=${AF_METRIC:-AF}

## Groups the max is taken over.
AF_GROUPS=${AF_GROUPS:-afr,amr,eas,nfe,sas}

## Groups deliberately left out. gnomAD reports an AF for every one of them, so
## each name here is doing real work.
EXCLUDE_GROUPS=${EXCLUDE_GROUPS:-ami,asj,fin,mid,oth,remaining}

## Tokens that disqualify a field whatever else it contains: sex-split and
## subset frequencies are not ancestry groups, and AF_raw is pre-filter.
## Without this, "gnomADe_AF_nfe_female" would be read as an nfe frequency.
EXCLUDE_TOKENS=${EXCLUDE_TOKENS:-male,female,xx,xy,raw,controls,non,cancer,neuro,topmed}

## The dataset-wide ungrouped field (gnomADe_AF) is NOT used by default: it
## pools every population including the excluded ones, which defeats the point.
## 1 = use it as a last resort when no per-group field exists.
ALLOW_OVERALL=${ALLOW_OVERALL:-0}

paste meta_CADD_head p > header_meta

(zgrep '^#' f5.vcf.gz && zgrep -v '^#' f5.vcf.gz | sort -u -k1,1 -k2,2 -k4,4 -k5,5 -k6,6) | bgzip > f5_dedup.vcf.gz
bcftools view -G f5_dedup.vcf.gz --threads ${THREADS:-1} -Ov -o p1.vcf

grep -v '#' p1.vcf >f6
cut -f 1-8 f6 >p1
cut -f 1-2,4-5 p1 >c1
cut -f 4 c1 >alt

## ---- CSQ field index: ALLELE_NUM, CADD and the frequency fields, by name ----
## A frequency field qualifies only if ALL of these hold, on underscore-separated
## TOKENS rather than substrings:
##   * one token is exactly AF_METRIC,
##   * one token is one of AF_GROUPS,
##   * no token is in EXCLUDE_GROUPS or EXCLUDE_TOKENS.
## Position-independent, so gnomADe_AFR_AF and gnomADe_AF_afr both match.
##
## Token equality is what makes this safe. "gnomADe_AC_afr" contains the
## substring "af" twice and any /AF/ regex reads it as a frequency -- that is the
## bug that made pre_1.sh report AF 0 for a variant whose true max AF was 0.22.
## It carries no "AF" token, so it is rejected here.
csq_f=$(grep -m1 "^##INFO=<ID=CSQ" p1.vcf | sed -E 's/.*Format: *//; s/">.*//; s/"$//')
[ -n "${csq_f}" ] || { echo "ERROR: no CSQ Format line in p1.vcf" >&2; exit 1; }

METRIC="${AF_METRIC}"

read -r NF_CSQ AN_I CADD_I GG_LIST GE_LIST GG_N GE_N < <(
  echo "${csq_f}" | awk -F'|' -v gg="${GNOMADG_PREFIX}" -v ge="${GNOMADE_PREFIX}" \
                             -v cadd="${CADD_FIELD}" -v metric="${METRIC}" \
                             -v groups="${AF_GROUPS}" -v exg="${EXCLUDE_GROUPS}" \
                             -v ext="${EXCLUDE_TOKENS}" -v allowall="${ALLOW_OVERALL}" '
  function pfx(name, p) { return (index(tolower(name), tolower(p)) == 1) }
  function norm(s) { gsub(/[ _]/, "_", s); return tolower(s) }
  # returns 1 keep, 0 reject; sets WHY
  function qualifies(name,   n, A, i, t, hasm, hasg) {
      WHY = ""
      n = split(tolower(name), A, "_")
      for (i = 1; i <= n; i++) {
          t = A[i]
          if (t in BADT) { WHY = "dropped:subset_token_" t; return 0 }
          if (t in EXCG) { WHY = "dropped:excluded_group_" t; return 0 }
          if (t == tolower(metric)) hasm = 1
          if (t in KEEP)            hasg = 1
      }
      if (!hasm) { WHY = "dropped:other_metric"; return 0 }
      if (!hasg) {
          if (allowall + 0 == 1) return 1          # ungrouped, used as last resort
          WHY = "dropped:ungrouped_pools_excluded_populations"
          return 0
      }
      return 1
  }
  {
      ng = split(groups, G, ","); for (i = 1; i <= ng; i++) KEEP[tolower(G[i])] = 1
      ne = split(exg,    E, ","); for (i = 1; i <= ne; i++) EXCG[tolower(E[i])] = 1
      nt = split(ext,    T, ","); for (i = 1; i <= nt; i++) BADT[tolower(T[i])] = 1
      for (i = 1; i <= NF; i++) {
          ## Space and underscore are treated as the same character when
          ## matching the CADD field: the plugin is seen writing both
          ## "CADD_RAW" and "CADD RAW" depending on version, and an exact
          ## match on the wrong one aborts the run.
          if ($i == "ALLELE_NUM")        an = i
          if (norm($i) == norm(cadd))    ci = i
          if (!pfx($i, gg) && !pfx($i, ge)) continue
          if (!qualifies($i)) continue
          if (pfx($i, gg)) { gl = gl (gl == "" ? "" : ",") i; gn = gn (gn == "" ? "" : ",") $i }
          if (pfx($i, ge)) { el = el (el == "" ? "" : ",") i; en = en (en == "" ? "" : ",") $i }
      }
      print NF, an+0, ci+0, (gl == "" ? "NA" : gl), (el == "" ? "NA" : el),
            (gn == "" ? "NA" : gn), (en == "" ? "NA" : en)
  }')

[ "${AN_I}"   -gt 0 ] || { echo "ERROR: ALLELE_NUM not in the CSQ header -- re-run VEP with --allele_number" >&2; exit 1; }
[ "${CADD_I}" -gt 0 ] || { echo "ERROR: CSQ field '${CADD_FIELD}' not found" >&2; exit 1; }

## Hard failure, not a warning, when neither source yields a usable field. A
## warning here would emit a column of ${AF_MISSING} for every variant and the
## pipeline would look like it had worked.
if [ "${GG_LIST}" = "NA" ] && [ "${GE_LIST}" = "NA" ]; then
    {
      echo "ERROR: no '${METRIC}' field found for any of the groups: ${AF_GROUPS}"
      echo "       gnomAD-prefixed fields that ARE in this CSQ header:"
      echo "${csq_f}" | tr '|' '\n' | grep -iE "^(${GNOMADG_PREFIX}|${GNOMADE_PREFIX})" | sed 's/^/         /'
      echo "       Expected per-group fields such as gnomADe_AFR_AF or gnomADe_AF_afr."
      echo "       If only an ungrouped gnomADe_AF is listed, re-run VEP so the"
      echo "       per-population frequencies are included, or set ALLOW_OVERALL=1"
      echo "       to use the ungrouped AF -- which pools ${EXCLUDE_GROUPS}."
    } >&2
    exit 1
fi
[ "${GE_LIST}" != "NA" ] || echo "WARNING: no ${GNOMADE_PREFIX} ${METRIC} group fields; genomes only" >&2
[ "${GG_LIST}" != "NA" ] || echo "WARNING: no ${GNOMADG_PREFIX} ${METRIC} group fields; exomes only" >&2

## Audit table. Fields that look like a frequency but were dropped are labelled
## with WHY, so the group exclusions are verifiable rather than implicit.
echo "${csq_f}" | awk -F'|' -v an="${AN_I}" -v ci="${CADD_I}" \
                           -v gl="${GG_LIST}" -v el="${GE_LIST}" \
                           -v metric="${METRIC}" -v groups="${AF_GROUPS}" \
                           -v exg="${EXCLUDE_GROUPS}" -v ext="${EXCLUDE_TOKENS}" \
                           -v allowall="${ALLOW_OVERALL}" \
                           -v gg="${GNOMADG_PREFIX}" -v ge="${GNOMADE_PREFIX}" '
function norm(s) { gsub(/[ _]/, "_", s); return tolower(s) }
function qualifies(name,   n, A, i, t, hasm, hasg) {
    WHY = ""
    n = split(tolower(name), A, "_")
    for (i = 1; i <= n; i++) {
        t = A[i]
        if (t in BADT) { WHY = "dropped:subset_token_" t; return 0 }
        if (t in EXCG) { WHY = "dropped:excluded_group_" t; return 0 }
        if (t == tolower(metric)) hasm = 1
        if (t in KEEP)            hasg = 1
    }
    if (!hasm) { WHY = "dropped:other_metric"; return 0 }
    if (!hasg) {
        if (allowall + 0 == 1) return 1
        WHY = "dropped:ungrouped_pools_excluded_populations"
        return 0
    }
    return 1
}
{
    ng = split(groups, G, ","); for (i = 1; i <= ng; i++) KEEP[tolower(G[i])] = 1
    ne = split(exg,    E, ","); for (i = 1; i <= ne; i++) EXCG[tolower(E[i])] = 1
    nt = split(ext,    T, ","); for (i = 1; i <= nt; i++) BADT[tolower(T[i])] = 1
    split(gl, a, ","); for (i in a) R[a[i]] = "gnomADg_max_" metric
    split(el, b, ","); for (i in b) R[b[i]] = (R[b[i]] == "" ? "" : R[b[i]] ",") "gnomADe_max_" metric
    R[ci] = "CADD_raw"; R[an] = "allele_number"
    print "index\tname\trole"
    for (i = 1; i <= NF; i++) {
        role = (i in R) ? R[i] : ""
        if (role == "" && (index(tolower($i), tolower(gg)) == 1 ||
                           index(tolower($i), tolower(ge)) == 1)) {
            qualifies($i); role = WHY
        }
        ## Frequency fields that are NOT gnomAD-prefixed are labelled too.
        ## VEP MAX_AF is the maximum across EVERY population it knows, the
        ## excluded ones included, so it is exactly the field a well-meaning
        ## reader would reach for and must be visibly rejected, not silent.
        else if (role == "" && qualifies($i) == 0 && WHY != "dropped:other_metric")
            role = "dropped:not_gnomAD_" substr(WHY, 9)
        print i "\t" $i "\t" role
    }
}' > csq_field_index.tsv

echo "CSQ fields ${NF_CSQ} | ALLELE_NUM ${AN_I} | ${CADD_FIELD} ${CADD_I}" >&2
echo "${METRIC} over groups ${AF_GROUPS} (excluded: ${EXCLUDE_GROUPS})" >&2
echo "  ${GNOMADE_PREFIX} [first]  ${GE_N}" >&2
echo "  ${GNOMADG_PREFIX} [fallback] ${GG_N}" >&2

## ---- allele alignment: one pass, capped at NMAX alleles ---------------------
## c_u holds the CSQ blocks re-ordered to ALT order, in the same row order as
## p1. Alleles past NMAX are dropped. A slot with no annotation (the * spanning
## deletion, or an allele VEP did not annotate) gets a placeholder of the right
## field width so the "|" count stays correct.
SB=$(awk -v n="${NF_CSQ}" 'BEGIN{s="*|*|*"; for(i=4;i<=n;i++) s=s"|"; print s}')
EB=$(awk -v n="${NF_CSQ}" 'BEGIN{s="";      for(i=2;i<=n;i++) s=s"|"; print s}')

awk -F'\t' -v AN="${AN_I}" -v NMAX="${NMAX}" -v SB="${SB}" -v EB="${EB}" '
/^#/ { next }
{
    nalt = split($5, ALT, ",")
    keep = (nalt > NMAX) ? NMAX : nalt
    if (nalt > NMAX) NOVER++
    if (nalt > 1)    NMULTI++
    NVAR++

    csq = ""
    ni = split($8, I, ";")
    for (i = 1; i <= ni; i++) if (substr(I[i],1,4) == "CSQ=") { csq = substr(I[i],5); break }

    delete OUT; delete NPER
    nb = (csq == "") ? 0 : split(csq, B, ",")
    for (i = 1; i <= nb; i++) {
        split(B[i], F, "|")
        k = F[AN] + 0
        if (k < 1 || k > keep) continue
        if (++NPER[k] == 1) OUT[k] = B[i]        # first block per allele
    }

    line = ""
    for (k = 1; k <= keep; k++) {
        if (k in OUT)           blk = OUT[k]
        else if (ALT[k] == "*") { blk = SB; NSTAR++ }
        else                    { blk = EB; NMISS++ }
        if (NPER[k] > 1) NMULTIBLK++
        line = line (k == 1 ? "" : ",") blk
    }
    print line
}
END {
    printf("variants %d | multi-allelic %d | * alleles %d\n", NVAR, NMULTI, NSTAR) > "/dev/stderr"
    if (NOVER > 0)
        printf("NOTE: %d variant(s) have >%d ALT alleles; alleles %d+ were IGNORED\n",
               NOVER, NMAX, NMAX+1) > "/dev/stderr"
    if (NMISS > 0)
        printf("WARNING: %d ALT allele(s) had no CSQ block; placeholder emitted\n", NMISS) > "/dev/stderr"
    if (NMULTIBLK > 0)
        printf("NOTE: %d allele(s) had >1 CSQ block (per transcript); first kept -- use VEP --pick\n",
               NMULTIBLK) > "/dev/stderr"
}' p1.vcf > c_u

## ---- consequence -----------------------------------------------------------
cut -f 2 -d'|' c_u > c2

## ---- gene with ensemblID; note: 806 x-genes cross chunks --------------------
echo "##fileformat=VCFv4.2" > f61.vcf
cut -f 1-8 f6 >> f61.vcf
bedtools intersect -wao -a f61.vcf -b p50.bed | grep -v '^#' | cut -f 1-5,12 > p1.bed
datamash -g 1,2,3,4,5 collapse 6 < p1.bed | cut -f 6 > c3
perl -ne 'print join("\n", split(/\,/,$_));print("\n")' c3 | sort -u | grep -E 'ENSG' > gene.lst

## ---- frequency: max over selected groups, EXOMES first ---------------------
## Order is exomes then genomes. A zero means "not seen in that dataset", so a
## zero exome value falls through to the genomes; if both are zero or absent the
## result is AF_MISSING.
awk -F"," -v OFS="\t" -v gl="${GG_LIST}" -v el="${GE_LIST}" \
    -v nmax="${NMAX}" -v miss="${AF_MISSING}" -v pad="${AF_PAD}" '
BEGIN { NUMRE = "^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$" }
function isnum(s) { return (s != "" && s != "." && s ~ NUMRE) }
function blockmax(lst,   n, a, i, v, best) {
    if (lst == "NA") return ""
    n = split(lst, a, ","); best = ""
    for (i = 1; i <= n; i++) { v = F[a[i]]
        if (isnum(v) && (best == "" || v + 0 > best + 0)) best = v }
    return best
}
function usable(v) { return (v != "" && v + 0 > 0) }
{
    line = ""
    for (k = 1; k <= nmax; k++) {
        val = pad
        if (k <= NF) {
            split($k, F, "|")
            e = blockmax(el)                      # gnomAD exomes first
            if (usable(e)) { val = e; NE++ }
            else {
                g = blockmax(gl)                  # then genomes
                if (usable(g))    { val = g; NG++ }
                else if (e != "") { val = e; NE++ }   # exome reported 0
                else if (g != "") { val = g; NG++ }   # genome reported 0
                else              { val = miss; NMISS++ }
            }
        }
        line = line (k == 1 ? "" : OFS) val
    }
    print line
}
END {
    printf("frequency source: %d allele(s) from exomes, %d from genomes, %d with no value (set to %s)\n",
           NE, NG, NMISS, miss) > "/dev/stderr"
}' c_u > c4

## ---- CADD raw --------------------------------------------------------------
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

##phred >=15
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<1.387112){$i="";}}1' c5 >c5a
##phred >=20
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<2.097252){$i="";}}1' c5 >c5b

##genotype
##zgrep -v '#' f5.vcf.gz | cut -f 10- | awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++) $i=substr($i,1,3)}1' >c6
