#!/usr/bin/env nextflow


phenotype_filename = params.in.replace(".tsv","")

process split_phenotypes {
    
    executor 'local'
    
    publishDir "analysis/${phenotype_filename}/phenotypes/", mode: 'copy'

    input:
        file 'infile.tsv' from Channel.fromPath(params.in)
    output:
        file '*.tsv' into phenotypes


    """
    #!/usr/bin/env python
    import re
    from slugify import slugify
    with open('infile.tsv', 'r') as f:
        lines = [re.split("[\\t|,]", x) for x in f.read().splitlines()]
        vars = lines[0][1:]
        for n, v in enumerate(vars):
            v_slug = slugify(v)
            file_path = "${phenotype_filename}-" + slugify(v) + ".tsv"
            with open(file_path, 'w') as p_out:
                p_out.write("strain\\t" + slugify(v) + "\\n")
                for l in lines[1:]:
                    if l[n+1] != "NA":
                        p_out.write(l[0] +"\\t" + l[n+1] + "\\n")
    """

}

phenotypes_split = phenotypes.flatten()


process perform_mapping {

    publishDir "analysis/${phenotype_filename}/mapping_Rdata/", mode: 'copy', pattern: '*.Rda'
    
    tag { input_tsv }
    
    input:
        file input_tsv from phenotypes_split

    output:
        file '*-mapping.tsv' into phenotype_mappings
        file '*-mapping.Rda' into phenotype_mappings_R
        file '*-genes.Rda' into phenotype_fine_mappings_R

    """
    #!/usr/bin/env Rscript --vanilla

    library(readr)
    library(cegwas)
    library(dplyr)
    library(data.table)

    df <- readr::read_tsv("${input_tsv}")

    phenotype_name <- names(df)[[2]]

    pheno <- process_pheno(df)
    mapping_df <- gwas_mappings(pheno,  cores = 1, mapping_snp_set = F)
    p_mapping_df <- process_mappings(mapping_df, phenotype_df = pheno, CI_size = 50, snp_grouping = 200)
    genes <- process_correlations(variant_correlation(p_mapping_df, condition_trait = F))

    readr::write_tsv(p_mapping_df, paste0(phenotype_name, '-mapping.tsv'))
    save(p_mapping_df, file = paste0(phenotype_name, '-mapping.Rda'))
    save(genes, file = paste0(phenotype_name, '-genes.Rda'))
    
    """
}

process concatenate_mappings {

    publishDir "analysis/${phenotype_filename}/mapping/", mode: 'copy'
    
    input:
        val(input_mapping) from phenotype_mappings.toSortedList()
    output:
        file("${phenotype_filename}-mappings.tsv")

    """
        cat ${input_mapping.join(" ")} > ${phenotype_filename}-mappings.tsv
    """

}
