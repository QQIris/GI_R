require(psych)

icc_cal <- function(input_path, conditions, colid, output_path){
    for (con in conditions){
        mzfile <- list.files(input_path, pattern = paste(con, '.+mztwin', sep = ''), full.names = TRUE)
        dzfile <- list.files(input_path, pattern = paste(con, '.+dztwin', sep = ''), full.names = TRUE)
        
        mzall <- read.csv(mzfile, header = TRUE)
        dzall <- read.csv(dzfile, header = TRUE)
        
        iccresults <- matrix(,length(colid),4,dimnames=list(colid,c('iccmz','p','iccdz','p')))
        
        for (n in 1:length(colid)){
            print(colid[n])
            rois   <- paste(colid[n],c(rep(1,1),rep(2,1)),sep="")
            mzdata <- mzall[rois]
            dzdata <- dzall[rois]
            
            mzresults <- ICC(na.omit(mzdata), missing=TRUE)
            dzresults <- ICC(na.omit(dzdata), missing=TRUE)
            iccresults[n,] <- c(mzresults$results$ICC[2], mzresults$results$p[2],
                                dzresults$results$ICC[2], dzresults$results$p[2])
            
          
        dir.create(output_path)
        write.csv(iccresults, paste(output_path, con, '_icc_resutls.csv', sep = ''))
        }
    }
    
}

colid <- c('body','LOS','MTG','ITG','OTS','face','OFA','FFA','place','TOS','PPA','RSC','tool','LO','ITGobject','pFs')
condition <- c('body', 'face', 'place', 'tool')

input_path <- '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_icc/regressoutHM+gss/rm_outliers/thr2/twin_form/'
output_path <- '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_icc/regressoutHM+gss/rm_outliers/thr2/results/'
icc_cal(input_path, condition, colid, output_path)
