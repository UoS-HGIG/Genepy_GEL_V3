# GENEPY_GEL_AGGV3

A Nextflow (DSL2) pipeline that computes [GenePy](https://github.com/UoS-HGIG/GenePy-2) per-gene
pathogenicity scores from the Genomics England **aggV3** aggregated VCFs, inside the
GEL Research Environment.

GenePy collapses the variants a person carries in a gene into a single per-gene score,
combining a deleteriousness metric (CADD), an allele frequency, and the individual's
genotype. This pipeline produces the per-gene input matrices and runs the scoring, shard
by shard, for one chromosome per invocation.

---

## Pipeline overview

One invocation processes **one chromosome**, driven by the aggV3 shard/subshard layout.
Each subshard VCF flows through eight processes:

| # | Process | Module | What it does |
|---|---|---|---|
| 1 | `CADD_score` | `modules/CADD.nf` | Strips genotypes (`bcftools view -G`) to make a sites-only `p1.vcf`, removes the `chr` prefix, and runs CADD 1.6 **on indels only** — SNV scores come from the precomputed `whole_genome_SNVs.tsv.gz` at the VEP step. Emits `wes_<subshard>.raw.tsv.gz`. |
| 2 | `site_qc_cadd15` | `modules/site_qc_cadd15.nf` | Applies aggV3 site QC (`MEDIAN_DP>=8`, `MEDIAN_GQ>=10`, `MISSINGNESS_RATE<=0.12`; the `_XX`/`_XY` variants on chrX), restricts to positions with CADD phred ≥ 15, and restores the site-QC file to its multi-allelic representation via `INFO/SOURCE_RECORD` so it can be intersected against `p1`. |
| 3 | `VEP_score` | `modules/VEP.nf` | VEP with the CADD plugin, and `--af_gnomade --af_gnomadg`. The CSQ block is then annotated back onto the genotype-bearing VCF with `bcftools annotate -c '+CSQ'`. |
| 4 | `Pre_processing_1` | `modules/Pre_pr1.nf` | Genotype QC. Adds `FORMAT/AB` (`LAD[:1]/DP`), then sets to missing any genotype failing `FT=PASS` and `DP>=8 & AB>=0.15`, keeping hom-ref calls regardless. |
| 5 | `Pre_processing_2` | `modules/Pre_pr2.nf` | Runs `templates/pre_1.sh`: parses the CSQ block into the per-allele annotation columns `c1`–`c5` (locus, consequence, gene, allele frequency, CADD raw), plus `gene.lst`. |
| 6 | `Pre_processing_3` | `modules/Pre_pr3.nf` | Runs `templates/pre_2.sh`: splits into per-gene `.meta` matrices at two CADD thresholds, `metafiles15_*` (phred ≥ 15) and `metafiles20_*` (phred ≥ 20). |
| 7 | `Reatt_Genes` | `modules/Gene_reattach.nf` | A gene straddling a shard boundary appears in more than one subshard. This concatenates those fragments into `dup15/` and `dup20/` so the gene is scored once, whole. |
| 8 | `Genepy_score` | `modules/Genepy.nf` | Floors zero/empty frequency columns to `3.98e-6`, then runs `templates/genepy.py` with the karyotype file. Emits `*.meta.txt` per gene. |

### Where the annotation logic lives

Steps 5 and 6 are thin Nextflow wrappers; the real work is in **`templates/pre_1.sh`** and
**`templates/pre_2.sh`**. `pre_1.sh` is where CADD and allele frequency are extracted from
the CSQ block, and it is the file to read (and change) when the annotation definition moves.

---

## Requirements

- Nextflow (DSL2)
- Containers, pulled from quay.io:
  - `quay.io/parsboy1987/cadd1.6_ens95:v1` — CADD, bcftools, the pre-processing shell steps
  - `quay.io/parsboy1987/vep115_bcftools:latest` — VEP
  - `quay.io/parsboy1987/genepy:v1` — Python/R for scoring
- GEL RE resources referenced by default in `nextflow.config`:
  - VEP cache at `/nas/weka.gel.zone/pgen_public_data_resources/vep_resources/VEP111/`
  - CADD 1.6 data at `/tools/aws-workspace-apps/CADD/1.6/CADD-scripts/data/GRCh38_v1.6/`

---

## Parameters

Required (no default — must be supplied):

| Parameter | Meaning |
|---|---|
| `--chr` | Chromosome to process, e.g. `chr22`. Selects shards via `multiallelic_shards_bed`. |
| `--shard_path` | Root of the aggV3 shard tree. Globbed as `<shard_path>/shard-*/subshard-*/dragen.vcf.gz`. |
| `--multiallelic_shards_bed` | Shard manifest. Column 1 is the chromosome, columns 5 and 6 the shard and subshard numbers. |
| `--base_site_qc` | Root of the site-QC tree: `<base_site_qc>/shard-N/subshard-M/dragen.gel.siteqc.vcf.gz`. |
| `--annotations_cadd` | CADD 1.6 annotation directory, symlinked into the container. |
| `--kary` | Karyotype / sex file passed to `genepy.py`. |
| `--plugin3` | Additional CADD ≥ 15 region BED, concatenated with the per-subshard regions. |
| `--genomad_indx1`, `--genomad_indx2` | Staged into the VEP task. See *Known issues*. |

Defaulted (override only if you need to):

| Parameter | Default |
|---|---|
| `--outDir` | `Results` |
| `--gene_code_bed` | `templates/gencode.v45.annotation_2.bed` |
| `--header_meta` | `header.meta_org` |
| `--genepy_py` | `templates/genepy.py` |
| `--templates` | `templates` |
| `--plugin1` | `whole_genome_SNVs.tsv.gz` (CADD SNVs) |
| `--plugin2` | `gnomad.genomes.r3.0.indel.tsv.gz` (CADD indels) |
| `--vep_plugins` | `templates/Plugins` |
| `--homos_vep` | VEP 111 cache path |

---

## Usage

```bash
nextflow run main.nf \
  --chr chr22 \
  --shard_path            /path/to/aggV3/shards \
  --multiallelic_shards_bed /path/to/multiallelic_shards.bed \
  --base_site_qc          /path/to/siteqc \
  --annotations_cadd      /tools/aws-workspace-apps/CADD/1.6/CADD-scripts/data/annotations/GRCh38_v1.6 \
  --plugin3               /path/to/cadd15_regions.bed \
  --kary                  /path/to/karyotype.txt \
  --outDir                Results_chr22 \
  -resume
```

`-resume` is worth using by default: the per-subshard tasks are expensive and independent.

### Scope of a run

`--chr` selects which shards are eligible, via the chromosome column of
`multiallelic_shards_bed`. On top of that, `main.nf` carries an explicit filter that narrows
the run further:

```groovy
chrx = Channel.fromPath(shard_path_pattern, checkIfExists: true)
.filter { vcf_file ->
    vcf_file.parent.parent.name == "shard-97" &&
    vcf_file.parent.name in ["subshard-15", "subshard-16"]
}
```


---

## Outputs

```
<outDir>/
  <shard>/<subshard>/           per-subshard intermediates (c1..c5, f3/f5/f6, VEP VCFs, meta_CADD*.txt)
  <chr>/
    metafiles15_*/              per-gene .meta matrices, CADD phred >= 15
    metafiles20_*/              per-gene .meta matrices, CADD phred >= 20
    dup15/, dup20/              genes reassembled across shard boundaries
    <chr>_<cadd>_dup.lst        which genes needed reassembly
    Genepy_score/15/*.meta.txt  GenePy scores
    Genepy_score/20/*.meta.txt
  <shard_dir_name>/Report.html, Timeline.html
```

### `.meta` column layout

Each row is a variant; each row carries the per-allele annotation followed by genotypes.
Columns 7–16 are the ten allele-frequency slots (one per ALT allele, up to ten; a site with
fewer alleles leaves the remainder blank) and columns 17–26 the matching CADD raw scores.
`Genepy_score` relies on this fixed offset when it floors the frequency columns, so the
column count in `pre_1.sh` and the `for (i=7;i<=16;i++)` loop in `modules/Genepy.nf` must
stay in step.

---

## Resources

Per-task, from `nextflow.config`, with memory and CPU doubling on each retry:

| Label | CPUs | Memory | maxForks |
|---|---|---|---|
| `CADD_score` | 8 | 64 GB | 20 |
| `site_qc_cadd15` | 8 | 64 GB | 30 |
| `VEP_score` | 8 | 64 GB | 20 |
| `Pre_processing_1` | 8 | 48 GB | 20 |
| `Pre_processing_2` | 8 | 16 GB | 20 |
| `Pre_processing_3` | 8 | 16 GB | 20 |
| `Reatt_Genes` | 8 | 16 GB | 10 |
| `Genepy_score` | 8 | 32 GB | 20 |

Global ceiling: `maxForks = 20`, `time = 24.h`, `maxRetries = 3`.

`Report.html` and `Timeline.html` are written per run and are the place to get real CPU-hour
figures for costing a larger run.

---

## Verifying the annotation fields after a run

`pre_1.sh` locates CADD, `ALLELE_NUM` and the gnomAD frequency fields **by name** from the
`##INFO=<ID=CSQ ...Format: ...>` header rather than by fixed column position, and writes
`csq_field_index.tsv` listing every CSQ field with the role it was given — including which
gnomAD fields were rejected and why.

Read that file after the first run. The frequency step takes the maximum across selected
genetic ancestry groups (`afr, amr, eas, nfe, sas`), deliberately excluding bottlenecked and
founder populations (`ami, asj, fin, mid`) and the `oth`/`remaining` bin, so that a founder
variant common in one such group is not treated as common generally. That requires
**per-population** gnomAD fields in the CSQ block; if VEP emitted only the ungrouped
`gnomADe_AF`/`gnomADg_AF`, the step fails loudly and prints the gnomAD fields it did find.

Note also that the frequencies come from whatever gnomAD release is bundled in the VEP cache
in use (`VEP111`), which is **not** the same as the gnomAD version behind aggV3 itself. Worth
stating explicitly in any write-up of the scores.

