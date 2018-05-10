require(psych)

anova_cal <- function(input_path, conditions, selective_dic, output_path){
    
    dir.create(output_path, showWarnings = FALSE)
    
    aov_results <- matrix(, nrow = 0, ncol = 4, dimnames = list(NULL, c('region','contrast','F','p')))
    for (con in conditions){
        fs <- list.files(input_path, pattern = con, full.names = TRUE)
        
        
        for (region in selective_dic[[con]]){
            print(region)
            all_data = data.frame()
            for (f in fs){
                data = read.csv(f, header = TRUE)
                contrast <- strsplit(f, '/')[[1]]
                contrast <- paste(strsplit(contrast[length(contrast)], '_')[[1]][1:3], collapse = '_')
                tmp <- subset(data, select = c('Family_ID', region, 'twin', 'condition'))
                
                colnames(tmp) <- c('Family_ID', 'region', 'twin', 'condition')
                fit <- aov(region ~ twin*condition+Error(Family_ID/condition), tmp)
                
                results <- summary(fit)
                results_sort <- list(region,
                                     contrast,
                                     results$`Error: Family_ID:condition`[[1]]$`F value`[2],
                                     results$`Error: Family_ID:condition`[[1]]$`Pr(>F)`[2])
                aov_results <- rbind(aov_results, results_sort)
                all_data <- rbind(all_data, tmp)
            }
            fit <- aov(region ~ twin*condition+Error(Family_ID/condition), all_data)
            
            results <- summary(fit)
            results_sort <- list(region, 
                                 paste(con, 'vs_3', sep = '_'),
                                 results$`Error: Family_ID:condition`[[1]]$`F value`[2],
                                 results$`Error: Family_ID:condition`[[1]]$`Pr(>F)`[2])
            aov_results <- rbind(aov_results, results_sort)
        }
        
    }
    rownames(aov_results) <- list(1:nrow(aov_results))[[1]]
    write.csv(aov_results, paste(output_path,'aov_results_detail.csv'))
}

selective_dic <- list(body = c('body','LOS','MTG','ITG','OTS'), 
                      face = c('face','OFA','FFA'), 
                      place = c('place','TOS','PPA','RSC'), 
                      tool = c('tool','LO','ITGobject','pFs'))

conditions <- c('body', 'face', 'place', 'tool')



input_path <- '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_pattern/gss+50twin/for_anova/'
output_path <- '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_pattern/gss+50twin/anova_results/'

anova_cal(input_path, conditions, selective_dic, output_path)




