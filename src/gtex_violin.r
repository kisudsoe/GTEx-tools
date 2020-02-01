suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
gtex_violin = function(
    gtex_factors = NULL,
    gtex_colids  = NULL,
    factor       = NULL,
    organ_rank   = 1,
    dir          = NULL,
    fig          = NULL # c(png, svg)
) {
    # 1. Prepare data
    gtex_values = gtex_factors[,3:ncol(gtex_factors)]
    rownames(gtex_values) = gtex_factors$Description
    gtex_values_colids = colnames(gtex_values)
    
    gtex_values_coldf = merge(data.frame(SAMPID=gtex_values_colids),gtex_colids,by='SAMPID')
    gtex_t            = cbind(t(gtex_values),SAMPID=gtex_values_colids)
    colnames(gtex_t)  = c(as.character(gtex_factors$Description),'SAMPID')
    gtex_values_coldf_= merge(gtex_values_coldf,gtex_t,by='SAMPID')
    
    gtex_values_coldf_$SEX[gtex_values_coldf_$SEX==1] <- 'Male'
    gtex_values_coldf_$SEX[gtex_values_coldf_$SEX==2] <- 'Female'
    gtex_values_coldf_ = gtex_values_coldf_[order(gtex_values_coldf_$SMTSD_SYM),]
    organs = levels(gtex_values_coldf_$SMTSD_SYM)
    
    dat_rows   = gtex_values_coldf_[,c(1,3:5,10)] # Extract gtex_colids
    col_which1 = which(colnames(gtex_values_coldf_)==factor)
    dat_val    = gtex_values_coldf_[,col_which1]
    
    # 2. Calculating aging corrolations r and its p-val by organs and sex
    dat_s = NULL; cor_df = NULL
    for(i in 1:length(organs)) {
        #print(paste0(i,'_',organs[i]))
        dat_organ = subset(gtex_values_coldf_,SMTSD_SYM %in% organs[i])
        dat_age_  = dat_organ$AGE_n
        col_which2= which(colnames(dat_organ)==factor)
        dat_val_  = dat_organ[,col_which2] #<- ALDH1A1, ALDH1A2, ALDH1A3
        dat_s_    = scale(as.numeric(unlist(dat_val_)))
        dat_s     = c(dat_s,dat_s_)
        cor_test  = cor.test(as.numeric(dat_age_),as.numeric(dat_val_))
		
        dat_org_m = subset(dat_organ,SEX=='Male')
        col_which3= which(colnames(dat_org_m)==factor)
        if(nrow(dat_org_m)==0) {
            cor_r_male=NA
            cor_p_male=NA
        } else {
        	cor_test_m = cor.test(
				as.numeric(dat_org_m$AGE_n),
				as.numeric(dat_org_m[,col_which3])) #<- ALDH1A1, ALDH1A2, ALDH1A3
			cor_r_male = cor_test_m$estimate
			cor_p_male = cor_test_m$p.value
        }

        dat_org_f = subset(dat_organ,SEX=='Female')
        if(nrow(dat_org_f)==0) {
            cor_r_female=NA
            cor_p_female=NA
        } else {
        	cor_test_f = cor.test(
				as.numeric(dat_org_f$AGE_n),
				as.numeric(dat_org_f[,col_which3])) #<- ALDH1A1, ALDH1A2, ALDH1A3
			cor_r_female = cor_test_f$estimate
			cor_p_female = cor_test_f$p.value
        }
        cor_result = data.frame(organ=organs[i],
                                cor_r=cor_test$estimate,
                                cor_p=cor_test$p.value,
                                cor_r_male=cor_r_male,
                                cor_p_male=cor_p_male,
                                cor_r_female=cor_r_female,
                                cor_p_female=cor_p_female)
        if(i==1) { cor_df = cor_result
        } else {
            cor_   = cor_result
            cor_df = rbind(cor_df,cor_)
        }
    }
    dat = cbind(dat_rows,value.s=dat_s,value=dat_val)
    cor_gene = cor_df[order(cor_df$cor_p),]
    f_name1 = paste0(dir,'/gene_',factor,'_cor.tsv')
    write.table(cor_gene,f_name1,sep='\t',row.names=F)
    paste0('Write file: ',f_name1,'\n') %>% cat
    
    # 3. Draw violin plot by best organ
    dat_= subset(dat, SMTSD_SYM %in% cor_gene$organ[organ_rank]) # filter by best organ
    title1 = paste0(factor,' (',cor_gene$organ[organ_rank],
                    '), r = ',round(cor_gene$cor_r[organ_rank],3),
                    ', p = ',cor_gene$cor_p[organ_rank]%>%formatC(format='e',digits=1))
    p=ggplot(dat_,aes(x=AGE,y=value.s))+theme_bw()+
        geom_violin(trim=F)+
        stat_summary(fun.data=mean_sdl,fun.args=list(mult=1),geom='crossbar',width=0)+
        #stat_summary(fun.data='mean_sdl',funargs=list(mult=1),geom='pointrange')+
        stat_summary(fun.y=median,fun.ymin=median,fun.ymax=median,geom='crossbar',width=0.4)+
        #geom_dotplot(binaxis='y',stackdir='center',fill='gray',
        #             stackratio=1,dotsize=0.5)+
        labs(title=title1,y='Gene expression')
        #labs(title='GDF11, r= -0.001, p= 0.8619',y='Gene expression')
    print(p)
    f_name2 = paste0(dir,'/gene_',factor,'_',cor_gene$organ[organ_rank],'_plot.',fig)
    if(fig=='png') { dev.copy(png,f_name2,width=7,height=5,units='in',res=300)
    } else dev.copy(svg,f_name2,width=7,height=5)
    dev.off()
    paste0('Draw plot1: ',f_name2,'\n') %>% cat
    
    table(dat_$SEX) %>% print
    title2 = paste0(factor,' (',cor_gene$organ[organ_rank],
					'); Male(r =',round(cor_gene$cor_r_male[organ_rank],3),
                    ', p =',cor_gene$cor_p_male[organ_rank]%>%formatC(format='e',digits=1),
                    '), Female (r =',round(cor_gene$cor_r_female[organ_rank],3),
                    ',p =',cor_gene$cor_p_female[organ_rank]%>%formatC(format='e',digits=1),')')
    p=ggplot(dat_,aes(x=AGE,y=value.s,colour=SEX))+#theme_bw()+
        geom_violin(trim=F)+
        #geom_boxplot(aes(colour=SEX))+
        #geom_dotplot(aes(fill=SEX),binaxis='y',stackdir='center',
        #             stackratio=1,dotsize=0.5,position=position_dodge(0.9))+
        stat_summary(fun.data='mean_sdl',fun.args=list(mult=1),
                     geom='crossbar',width=0,position=position_dodge(0.9))+
        #stat_summary(aes(colour=SEX),fun.data='mean_sdl',funargs=list(mult=1),
        #             geom='pointrange',position=position_dodge(0.9))+
        stat_summary(fun.y=median,fun.ymin=median,fun.ymax=median,
                     geom='crossbar',width=0.4,position=position_dodge(0.9))+
        labs(title=title2,y='Gene expression')
    print(p)
    f_name3 = paste0(dir,'/gene_',factor,'_',cor_gene$organ[organ_rank],'_sex_plot.',fig)
    if(fig=='png') { dev.copy(png,f_name3,width=9,height=5,units='in',res=300)
    } else dev.copy(svg,f_name3,width=9,height=5)#,units='in',res=300)
    while(!is.null(dev.list())) dev.off()
    paste0('Draw plot2: ',f_name3,'\n') %>% cat
}