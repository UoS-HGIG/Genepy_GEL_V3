/*
    Module: Make_gene_ann
    Builds the gene-assignment BED used by Pre_processing_2 (Genecode_Ann) from the
    gene-boundary BED6 (gene_bed): each gene's core plus its CRE, one row per gene.
    CRE = +/- cre_flank bp around the most upstream TSS (gene start on +, gene end on -),
    the same definition genepy.py uses for scoring, so assignment and scoring agree.
    Runs only when params.Genecode_Ann is left empty.
*/
process Make_gene_ann {
    tag "core+CRE BED (flank ${flank} bp)"
    executor 'local'

    input:
    path(gene_bed)
    val(flank)

    output:
    path("gene_core_CRE.bed")

    shell:
    """
    awk -v F=${flank} 'BEGIN{OFS="\\t"} {if(\$6=="+"){t=\$2; s=t-F; e=(\$3>t+F+1?\$3:t+F+1)} else {t=\$3-1; s=(\$2<t-F?\$2:t-F); e=\$3+F} if(s<0)s=0; print \$1,s,e,\$4,\$5,\$6}' ${gene_bed} > gene_core_CRE.bed
    """
}
