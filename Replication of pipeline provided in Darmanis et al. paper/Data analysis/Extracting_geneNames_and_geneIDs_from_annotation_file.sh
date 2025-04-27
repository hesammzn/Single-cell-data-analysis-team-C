awk '$3 == "gene" { 
    match($0, /gene_id "[^"]+"/, a); 
    match($0, /gene_name "[^"]+"/, b); 
    if (a[0] && b[0]) 
        print a[0] "\t" b[0] 
}' gencode.v19.chr_patch_hapl_scaff.annotation.gtf | sed 's/gene_id "//;s/"//;s/gene_name "//' | sort -u > gene_id_gene_name.tsv
