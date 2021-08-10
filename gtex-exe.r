#!/usr/bin/env Rscript

suppressMessages(library(dplyr))

### Functions for GTEx plots.ipynb ###
gtex_rds = function(
    dir  = NULL, # Directory name
    name = NULL  # Name of gene list
) {
    source('src/pdtime.r'); t0=Sys.time()
    ifelse(!dir.exists(dir), dir.create(dir),'')
    gtex = readRDS('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz.rds')
    gtex %>% dim %>% print
    pdtime(t0,2) %>% message
    
    # Save filtered gene_tpm data as RDS file
    paste0('Save filtered gene_tpm data as RDS file\n') %>% cat
    fctr = read.delim(paste0(name,'.tsv'))
    paste0('  Factors\t=') %>% cat; nrow(fctr) %>% print
    if(nrow(fctr)<50) { paste(fctr[,2],collapse=', ') %>% message }
    Ensgids = unlist(strsplit(as.character(gtex$Name),'[.]'))[c(TRUE,FALSE)]
    paste0('  Total Ensgids n\t=') %>% cat; length(Ensgids) %>% print
    gtex_factors = fctr$ENSGid
    paste0('  Matched Ensgids n\t=') %>% cat; length(gtex_factors) %>% print
    Ensgids_which = which(Ensgids %in% gtex_factors); #Ensgids_which %>% print
    gtex_factors = gtex[Ensgids_which,]
    paste0('  dim of data\t\t=') %>% cat; dim(gtex_factors) %>% print
    f_name = paste0(dir,'/gtex_',name,'.rds')
    saveRDS(gtex_factors,f_name)
    paste0('Write RDS file: ',f_name,'\n') %>% cat
    pdtime(t0,1) %>% message
}

gtex_aging_test = function(
    dir  = NULL,
    name = NULL,
    k    = NULL # organ rank
) {
    gtex_factors = readRDS(paste0(dir,'/gtex_',name,'.rds'))
    paste0('Read factor RDS file\t=') %>% cat; dim(gtex_factors) %>% print
    gtex_colids  = read.delim('gtex_SAMPID_SUBJID_SMTSD_v8.tsv')
    paste0('Read Gtex meta file\t=') %>% cat; dim(gtex_colids) %>% print
    
    source('src/gtex_violin.r')
    factor = gtex_factors$Description
    if(!is.null(k)) { k = k
    } else k = 1
    n = length(factor)
    for(i in 1:n) { # Debug: Why only one plot can draw?
        paste0('**',i,'/',n,', Draw plot for ',factor[i],'... \n') %>% cat
        gtex_violin(
            gtex_factors = gtex_factors,
            gtex_colids  = gtex_colids,
            factor       = factor[i],
            organ_rank   = k, # 1,12,30
            dir          = dir,
            fig          = 'png'
        )
        'done\n\n' %>% cat
    }
}
### End: GTEx plots.ipynb ###


### Functions for GTEx heatmap.ipynb ###
age_over = function(
    f_gene = NULL,
    f_age  = NULL,
    out    = NULL
) {
    # Read files
    paste0('Read, ',f_gene,' = ') %>% cat
    geneset = read.delim(f_gene)
    dim(geneset) %>% print

    paste0('Read, ',f_age,' = ') %>% cat
    age_li = readRDS(f_age)
    length(age_li) %>% print

    # Draw venn diagram
    paste0('\n  Extract gene lists.. ') %>% cat
    ensg_geneset = geneset$Ensgid %>% unique
    age_df = data.table::rbindlist(age_li)
    ensg_gtex = age_df$genes %>% unique
    ensg_li = list(ensg_geneset,ensg_gtex)
    names(ensg_li) = c(tools::file_path_sans_ext(f_gene %>% basename),'Age-correlated genes_Bonf0.05')
    'done\n' %>% cat

    source('src/venn_analysis.r')
    o = venn_fig(ensg_li,out)


    # Aggregate input geneset table by Ensgid
    paste0('\n  Extract Ensgids = ') %>% cat
    ensg_uniq = geneset$Ensgid %>% unique
    n = length(ensg_uniq)
    print(n)

    # Overlap input geneset to GTEx Age DB
    paste0('  Overlap the genes with GTEx Age DB:\n\n') %>% cat
    n = length(age_li)
    age_geneset_li = lapply(c(1:n),function(i) {
        paste0('  ',i,' ',age_li[[i]]$Tissue[1],'... ') %>% cat
        age_sub = subset(age_li[[i]], genes %in% geneset$Ensgid)
        if(nrow(age_sub)>0) {
            dim(age_sub) %>% print
            return(age_sub)
        } else {
            'NULL\n' %>% cat
            return(NULL)
        }
    })

    # Remove NULL element
    paste0('  Removing NULL from GTEx = ') %>% cat
    age_geneset_li = Filter(Negate(is.null),age_geneset_li)
    length(age_geneset_li) %>% print

    # Merge annoataions
    paste0('  Merge annotations = ') %>% cat
    age_geneset_df = data.table::rbindlist(age_geneset_li)
    colnames(age_geneset_df) = c('Tissue','Ensgid','Age.Coef','AveExpr','t','P.Value','Bonf','B')
    age_geneset_ann = merge(age_geneset_df,geneset,by='Ensgid',all.x=T)
    dim(age_geneset_ann) %>% print

    # Save as TSV file
    f_name1 = paste0(out,'/age_bonf0.05_ann.tsv')
    write.table(age_geneset_ann,f_name1,quote=F,row.names=F,sep='\t')
    paste0('Write TSV: ',f_name1,'\n') %>% cat

    # Save as RDS file
    f_name2 = paste0(out,'/age_bonf0.05_ann.rds')
    saveRDS(age_geneset_ann,f_name2)
    paste0('Write RDS: ',f_name2,'\n') %>% cat
}


age_hist = function(
    f_aging = NULL,
    out     = NULL
) {
    suppressMessages(library(reshape))
    suppressMessages(library(ggplot2))
    suppressMessages(library(plyr))

    # Read file
    paste0('Read, ',f_aging,' = ') %>% cat
    age_geneset_df = readRDS(f_aging)
    m = ncol(age_geneset_df)
    col_nm = colnames(age_geneset_df)
    dim(age_geneset_df) %>% print
    paste0('  Gene n = ') %>% cat; unique(age_geneset_df$Ensgid) %>% length %>% print


    # Preparing histogram dataset
    paste0('  Preparing histogram dataset... ') %>% cat
    ## Count tissue genes
    paste0('count tissue.. ') %>% cat
    tb_df = as.character(age_geneset_df$Tissue) %>% table %>% data.frame
    colnames(tb_df) = c('Tissues','Gene_n')
    length(tb_df$Tissues) %>% print

    ## descending order by Gene_n
    tb_df$Tissues = factor(tb_df$Tissues,levels=tb_df$Tissues[order(-tb_df$Gene_n)])

    ## Up-regulated genes
    age_up = subset(age_geneset_df,Age.Coef>0)
    up_age_n = table(age_up$Tissue) %>% data.frame
    colnames(up_age_n) = c('Tissues','Age_up')
    df1 = merge(tb_df,up_age_n,by='Tissues',all.x=T)

    ## Down-regulated genes
    age_dn = subset(age_geneset_df,Age.Coef<0)
    dn_age_n = table(age_dn$Tissue) %>% data.frame
    colnames(dn_age_n) = c('Tissues','Age_down')
    df2 = merge(df1,dn_age_n,by='Tissues',all.x=T)

    ## Reshape the table
    df3 = melt(df2[,c(1,3,4)],id='Tissues')
    colnames(df3) = c('Tissues','Direction','Gene_n')
    df3$Direction = factor(df3$Direction,levels=c('Age_down','Age_up'))
    df = ddply(df3,c('Tissues'),transform,label_y=cumsum(Gene_n))


    # Draw histogram
    p=ggplot(df,aes(x=Tissues,y=Gene_n,fill=Direction))+theme_bw()+
        geom_bar(stat='identity')+
        scale_fill_manual(values=c('cyan3','coral1'))+
        geom_text(aes(y=label_y,label=Gene_n),
                hjust=.3,vjust=.5,angle=45)+
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
    print(p)
    f_name = paste0(out,'/hist_bonf0.05.png')
    ggsave(f_name,p,width=6,height=6,units='in',dpi=300)
    paste0('  Write PNG: ',f_name,'\n\n') %>% cat
}


age_heat = function(
    f_aging = NULL,
    out     = NULL,
    stat    = 'bonf0.05',
    over_n  = 1,       # Tissue count criteria
    wh      = c(9,7) # Heatmap size
) {
    suppressMessages(library(tidyr))
    suppressMessages(library(ComplexHeatmap))
    suppressMessages(library(circlize))
    suppressMessages(library(dendextend))

    # 1. Read file
    paste0('Read, ',f_aging,' = ') %>% cat
    pro_age = readRDS(f_aging)
    dim(pro_age) %>% print

    ## Extract gene expression slope values
    organ_df = data.frame(
        Ensgid = pro_age$Ensgid,
        Organs = pro_age$Tissue
    )
    gene_organ = table(organ_df$Ensgid) %>% data.frame
    paste0('  Tissue common/specific genes:\n') %>% cat
    table_gene_organ   = gene_organ[order(-gene_organ$Freq),]
    table_gene_organ_n = table(table_gene_organ$Freq) %>% data.frame
    colnames(table_gene_organ_n) = c('Organ_n','Gene_n')
    tb = table_gene_organ_n[order(table_gene_organ_n$Organ_n,decreasing=T),]
    print(tb)

    ## Write TSV file
    f_name = paste0(out,'/age_heatmap_',stat,'.tsv')
    write.table(tb,f_name,sep='\t',row.names=F,quote=F)
    paste0('  Write TSV: ',f_name,'\n') %>% cat


    # 2. Count tissue number by each gene
    paste0('\nCount Gene-Tissue numbers = ') %>% cat
    organ_n_df = table(pro_age$Ensgid) %>% data.frame
    colnames(organ_n_df) = c('Ensgid','Organ_n')
    organ_ann = pro_age %>% left_join(organ_n_df,by='Ensgid')
    dim(organ_ann) %>% print
    paste0('  Gene N = ',unique(organ_ann$Ensgid) %>% length,'\n') %>% cat


    # 3. Unique organ numbers
    organs = unique(pro_age$Tissue)
    paste0('\nOrgans_n = ',length(organs),'\n') %>% cat

    ## Paste as Gene.name_Ensgid
    paste0('  Paste Gene name and Engsid... ') %>% cat
    organ_ann$sybl_ensg = paste0(
        organ_ann$Symbol,' ',
        organ_ann$Ensgid)
    paste0('done\n') %>% cat

    ## Spread Age.Coef table by tissue
    paste0('  Spread Tissue-Age.Coef = ') %>% cat
    col_nm = c('Tissue','Age.Coef','sybl_ensg','Organ_n')
    organ_ann_df = organ_ann %>% dplyr::select(col_nm %>% all_of) %>% spread(Tissue, Age.Coef)
    organ_ann_df = organ_ann_df[order(-organ_ann_df$Organ_n),]
    dim(organ_ann_df) %>% print

    ## Add aggregated organ list
    paste0('  Add organ list.. ') %>% cat
    organ_ann_df$organ_list = apply(organ_ann_df,1,function(row) {
        row1 = row[c(3:length(row))]
        row1 = row1[!is.na(row1)]
        return(paste0(names(row1),collapse=', '))
    }) %>% unlist
    paste0('done\n') %>% cat
    
    ## Row sort by directions
    col_nm = colnames(organ_ann_df)
    n = ncol(organ_ann_df)
    for(i in n:3) {
        organ_ann_df = organ_ann_df[order(organ_ann_df %>% dplyr::select(col_nm[i]),decreasing=T),]
    }
    organ_ann_df = organ_ann_df[order(organ_ann_df$Organ_n,decreasing=T),]

    ## Save as RDS file
    f_name = paste0(out,'/heatmap_',stat,'_beta.rds')
    saveRDS(organ_ann_df,f_name)
    paste0('  Write RDS: ',f_name,'\n') %>% cat


    # 4. Prepare heatmap matrix
    paste0('\nPrepare to matrix... \n') %>% cat

    ## Row reorder by directions
    '  reorder.. row.. ' %>% cat
    over_n_df = subset(organ_ann_df,Organ_n>=over_n)

    ## Column reorder by beta hits
    'col.. ' %>% cat
    m1 = ncol(over_n_df)-1
    col_nm1 = colnames(over_n_df)
    over_subcol1 = over_n_df %>% dplyr::select(col_nm1[c(3:m1)] %>% all_of)
    col_na = apply(over_subcol1,2,function(x) { sum(is.na(x)) })
    col_nm2 = c(colnames(over_n_df)[1:2], names(col_na)[order(col_na)])
    over_n_df = over_n_df %>% relocate(col_nm2 %>% all_of)

    ## Row reorder by tissues
    'tissue.. ' %>% cat
    col_nm3 = colnames(over_n_df)
    for(i in m1:3) {
        over_n_df = over_n_df[order(
            over_n_df %>% dplyr::select(col_nm3[i] %>% all_of),
            decreasing=T),]
    }

    ## Row reorder by directions
    '\n  direction.. ' %>% cat
    over_subcol2 = over_n_df %>% dplyr::select(col_nm1[c(3:m1)])
    row_which = apply(over_subcol2,1,function(row) {
        up_n = which(row>0) %>% length
        dn_n = which(row<0) %>% length
        if(up_n>0 & dn_n==0)       return(1)
        else if(up_n==0 & dn_n>0)  return(2)
        else if(up_n>0  & dn_n>0)  return(3)
        else if(up_n==0 & dn_n==0) return(NA)
    }) %>% unlist
    over_n_df = over_n_df[order(row_which),]
        
    ## Row split
    'split.. ' %>% cat
    organ_split = factor(over_n_df$Organ_n)
    organ_split_levels = levels(organ_split)
    k = length(organ_split_levels)
    organ_split = factor(
        organ_split,
        levels = organ_split_levels[k:1]
    )
        
    ## Convert df to matrix
    'convert.. ' %>% cat
    m = ncol(over_n_df)-1
    col_nm = colnames(over_n_df)[3:m]
    mat1 = as.matrix(over_n_df %>% dplyr::select(all_of(col_nm)))
    mat1 = mat1[,colSums(is.na(mat1))<nrow(mat1)] # remove NA only columns
    row.names(mat1) = over_n_df$sybl_ensg
    'done\n' %>% cat


    # 5. Draw heatmap
    paste0('  heatmap dim = ') %>% cat
    dim(mat1) %>% print
    f_name = paste0(out,'/beta_organ_',stat,'_over_',over_n,'.png')
    mat_max = max(mat1,na.rm=T); mat_min = min(mat1,na.rm=T)
    col.rg = c(mat_min,0,mat_max)
    cell.cols = colorRamp2(col.rg,c('deepskyblue','white','indianred1'))
    h = Heatmap(
        mat1,
        col=cell.cols,
        name='Beta',
        cluster_columns=F,
        cluster_rows=F,
        split=organ_split,
        row_title='Aging genes',
        column_title='Tissue',
        row_names_side='left',
        rect_gp=gpar(col='black'),
        column_names_max_height=unit(3,"in"),
        row_names_max_width = unit(3,"in")
    )
    png(f_name, width=wh[1], height=wh[2], units='in', res=150)
    print(h)
    #draw(h,heatmap_legend_side='right')
    #dev.copy(png,f_name,width=wh[1],height=wh[2],units='in',res=150)
    paste0('  Draw PNG file: ',f_name,'\n\n') %>% cat
    #dev.off()
    graphics.off() # killing all devices
}

### Function for /src/gtex_violin.r ###
age_violin = function(
    pro_age, gene_exp, meta,
    ensgid     = NULL,
    organ_rank = NULL,
    out        = NULL
) {
    suppressMessages(library(ggplot2))
    suppressMessages(library(Hmisc))

    # 1. Filter & Prepare data
    '\nFilter & Prepare GTEx transcriptome data.. ' %>% cat
    ensgids = lapply(gene_exp$Name, function(x) { return(strsplit(x,"\\.")[[1]][1]) }) %>% unlist
    row_which = which(ensgids==ensgid)
    gtex_rownm = gene_exp[row_which,2] #%>% paste0(collapse=' ')
    gtex_values = gene_exp[row_which,3:ncol(gene_exp)] %>% unlist # vector
    gtex_values_colids = colnames(gene_exp)[3:ncol(gene_exp)]
    'done\n' %>% cat

    'Merging meta data.. ' %>% cat
    gtex_values_coldf = merge(data.frame(SAMPID=gtex_values_colids),meta,by='SAMPID')
    gtex_t = data.frame(exp=gtex_values,SAMPID=gtex_values_colids)
    colnames(gtex_t) = c(gtex_rownm,'SAMPID')
    gtex_values_coldf_ = merge(gtex_values_coldf,gtex_t,by='SAMPID')
    dim(gtex_values_coldf_) %>% print

    'Convert gender and extract data.. ' %>% cat
    gtex_values_coldf_$SEX[gtex_values_coldf_$SEX==1] <- 'Male'
    gtex_values_coldf_$SEX[gtex_values_coldf_$SEX==2] <- 'Female'
    gtex_values_coldf_ = gtex_values_coldf_[order(gtex_values_coldf_$SMTSD_SYM),]
    'done\n' %>% cat

    # 2. Scale gene expressions
    '\nScale gene expression values.. ' %>% cat
    pro_age_sub = subset(pro_age, Ensgid==ensgid)
    organs = pro_age_sub$Tissue
    organ_rank = organ_rank %>% as.numeric
    dat_organ = subset(gtex_values_coldf_, SMTSD_SYM %in% organs[organ_rank])
    dat_rows = dat_organ[,c(1,3:4,6,10)]
    dat_val  = log10(dat_organ[,11]+1)
    dat_s  = scale(as.numeric(dat_val))
    dat = cbind(dat_rows,value.s=dat_s,value=dat_val)
    dim(dat) %>% print

    # 3. Draw violin plot
    paste0('Draw plots: ',organs[organ_rank],'\n') %>% cat
    dat_ = subset(dat, SMTSD_SYM %in% organs[organ_rank])
    title = paste0(
        gtex_rownm,' (',organs[organ_rank], '), Age.Coef = ',round(pro_age$Age.Coef[organ_rank],3),
        ', Bonf = ',pro_age$Bonf[organ_rank]%>%formatC(format='e',digits=1)
    )
    
    f_name1 = paste0(out,'/',gtex_rownm,'_',organs[organ_rank],'_scatter.png')
    p1=ggplot(dat_,aes(x=AGE_db,y=value.s,colour=SEX))+theme_bw()+
        geom_point()+
        geom_smooth(method=lm)+
        labs(title=title,y='Gene expression')
    ggsave(f_name1,p1,width=6,height=5,units='in',dpi=300)

    paste0('  ',organs[organ_rank],', Draw plot 1: ',f_name1,'\n') %>% cat

    p2=ggplot(dat_,aes(x=AGE,y=value.s,colour=SEX))+theme_bw()+
        geom_violin(trim=F)+
        stat_summary(fun.data='mean_sdl',fun.args=list(mult=1),
                     geom='crossbar',width=0,position=position_dodge(0.9))+
        stat_summary(fun=median,fun.min=median,fun.max=median,
                     geom='crossbar',width=0.3,position=position_dodge(0.9))+
        labs(title=title,y='Gene expression')
    f_name2 = paste0(out,'/',gtex_rownm,'_',organs[organ_rank],'_violin.png')
    ggsave(f_name2,p2,width=8,height=5,units='in',dpi=300)
    paste0('  ',organs[organ_rank],', Draw plot 2: ',f_name2,'\n') %>% cat
}

age_violin_multi = function(
    f_aging      = NULL,
    ensgid_multi = NULL,
    organ_rank   = NULL,
    out          = NULL
) {
    # 1. Read files
    paste0('Read, ',f_aging,' = ') %>% cat
    pro_age = readRDS(f_aging)
    pro_age = pro_age[order(pro_age$Bonf),]
    dim(pro_age) %>% print
    f_gene_exp = 'gtex_age_gene_expression_bonf0.05.rds'
    paste0('Read, ',f_gene_exp,' = ') %>% cat
    gene_exp = readRDS(f_gene_exp)
    dim(gene_exp) %>% print
    f_meta = 'gtex_SAMPID_SUBJID_SMTSD_v8.tsv'
    paste0('Read, ',f_meta,' = ') %>% cat
    meta = read.delim(f_meta)
    dim(meta) %>% print

    # 2. Run age_violin with loop
    ensgid_vec = strsplit(ensgid_multi,',')[[1]]
    n = length(ensgid_vec)
    for(i in 1:n) {
        age_violin(pro_age,gene_exp,meta,ensgid_vec[i],organ_rank,out)
    }
}
### End: /src/gtex_violin.rs###


# Parsing args
## Ref: https://stackoverflow.com/questions/3433603/parsing-command-line-arguments-in-r-scripts/3434100#3434100
## Document: https://cran.r-project.org/web/packages/argparser/argparser.pdf
suppressMessages(library(argparser))

## Functions
p = arg_parser('gtex-exe.r | 2020-08-21')
p = add_argument(p,'--example', flag=T, help='See command examples')
p = add_argument(p,'--age',  flag=T, help='Function for overlapping input gene list to GTEx age DB')
p = add_argument(p,'--hist', flag=T, help='Function for drawing histogram by tissues')
p = add_argument(p,'--heat', flag=T, help='Function for drawing heatmap of aging genes')
p = add_argument(p,'--viol', flag=T, help='Function for drawing violin plot of an aging gene')

## Arguments
p = add_argument(p,'--gene',    help='Input TXT file as gene list with format: <Symbol> <Name> <Ensgid>')
p = add_argument(p,'--geneage', help='Input RDS file from --age function')
p = add_argument(p,'--out',     help='Out dir path')
p = add_argument(p,'--dbage',   help='RDS file path of GTEx Age DB')
p = add_argument(p,'--ensgid',  help='Input a Ensg ID for drawing violin plot')
p = add_argument(p,'--organrank',help='Input a organ rank number to draw violin plot')

## Command examples
argv = parse_args(p)
if(argv$exam==TRUE) {'
Command Examples:

    Rscript gtex-exe.r --help
    Rscript gtex-exe.r --age --gene OUT_DIR/input_22.tsv --dbage gtex_age_list_bonf_0.05.rds --out OUT_DIR
    Rscript gtex-exe.r --hist --geneage OUT_DIR/age_bonf0.05_ann.rds --out OUT_DIR
    Rscript gtex-exe.r --heat --geneage OUT_DIR/age_bonf0.05_ann.rds --out OUT_DIR
    Rscript gtex-exe.r --viol --geneage OUT_DIR/age_bonf0.05_ann.rds --ensgid ENSG00000099810,ENSG00000137145 --organarnk 1 --out OUT_DIR
' %>% cat
}

# Run functions by parsed args
source('src/pdtime.r'); t0=Sys.time()
if(argv$age==TRUE) {
    '\n** Function age_over start **\n' %>% cat
    age_over(argv$gene,argv$dbage,argv$out)
} else if(argv$hist==TRUE) {
    '\n** Function age_hist start **\n' %>% cat
    age_hist(argv$geneage,argv$out)
} else if(argv$heat==TRUE) {
    '\n** Function age_heat start **\n' %>% cat
    age_heat(argv$geneage,argv$out)
} else if(argv$viol==TRUE) {
    '\n** Function age_violin_multi start **\n' %>% cat
    age_violin_multi(argv$geneage,argv$ensgid,argv$organrank,argv$out)
}
pdtime(t0,1) %>% message
