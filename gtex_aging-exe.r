suppressMessages(library(dplyr))

gtex_rds = function(
    dir  = NULL, # Directory name
    name = NULL  # Name of gene list
) {
    source('src/pdtime.r'); t0=Sys.time()
    dir.create(dir)
    gtex = readRDS('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz.rds')
    gtex %>% dim %>% print
    pdtime(t0,2) %>% message
    
    # Save filtered gene_tpm data as RDS file
    paste0('Save filtered gene_tpm data as RDS file\n') %>% cat
    fctr = read.delim(paste0(name,'.tsv'))
    paste0('  Factors\t=',) %>% cat; nrow(fctr) %>% print
    if(nrow(fctr)<50) { paste(fctr[,2],collapse=', ') %>% message }
    Ensgids = unlist(strsplit(as.character(gtex$Name),'[.]'))[c(TRUE,FALSE)]
    paste0('  Total Ensgids n\t=',) %>% cat; length(Ensgids) %>% print
    gtex_factors = fctr$ENSGid
    paste0('  Matched Ensgids n\t=',) %>% cat; length(gtex_factors) %>% print
    Ensgids_which = which(Ensgids %in% gtex_factors); #Ensgids_which %>% print
    gtex_factors = gtex[Ensgids_which,]
    paste0('  dim of data\t=') %>% cat; dim(gtex_factors) %>% print
    f_name = paste0(dir,'/gtex_',name,'.rds')
    saveRDS(gtex_factors,f_name)
    paste0('Write RDS file: ',f_name,'\n') %>% cat
    pdtime(t0,1) %>% message
}

gtex_aging = function(
    dir  = NULL,
    name = NULL,
    k    = NULL # temp variable for debug
) {
    source('src/pdtime.r'); t0=Sys.time()
    gtex_factors = readRDS(paste0(dir,'/gtex_',name,'.rds'))
    paste0('Read factor RDS file\t=') %>% cat; dim(gtex_factors) %>% print
    gtex_colids  = read.delim('gtex_SAMPID_SUBJID_SMTSD.tsv')
    paste0('Read Gtex meta file\t=') %>% cat; dim(gtex_colids) %>% print
    
    source('src/gtex_violin.r')
    factor = gtex_factors$Description
    n = length(factor)
    for(i in k:k) { # Debug: Why only one plot can draw?
        paste0('**',i,'/',n,', Draw plot for ',factor[i],'... ') %>% cat
        gtex_violin(
            gtex_factors = gtex_factors,
            gtex_colids  = gtex_colids,
            factor       = factor[i],
            organ_rank   = 1, # 1,12,30
            dir          = dir,
            fig          = 'png'
        )
        'done\n\n' %>% cat
    }
}