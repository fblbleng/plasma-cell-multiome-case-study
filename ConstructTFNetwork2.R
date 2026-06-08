# ============================================================
# ConstructTFNetwork2.R
#
# A patched copy of hdWGCNA's ConstructTFNetwork(), required by
# 19_hdwgcna_tf_regulatory_network.R.
#
# Why this exists. The stock hdWGCNA ConstructTFNetwork was
# written against xgboost 1.x. xgboost 2.0.0 changed two things
# this function depends on:
#   (a) the cross-validation model accessor moved from
#       xgb$models[[i]] to xgb$cv_predict$models[[i]]
#   (b) the CV-predict callback was renamed from cb.cv.predict
#       to xgb.cb.cv.predict
# On a server with xgboost >= 2.0.0 the stock function errors.
# This version detects the xgboost major version and branches
# accordingly.
#
# It also subsets the motif table on the 'motif_name' column,
# which is reliably present in GetMotifs() output, rather than a
# 'gene_name' column that some hdWGCNA versions do not return.
#
# Usage:
#   source("ConstructTFNetwork2.R")
#   seu <- ConstructTFNetwork2(seu, model_params = model_params)
#
# Dependencies: hdWGCNA (for Get/Set accessors), xgboost, dplyr.
# Requires stats::cor to be the resolved cor() (cor is masked by
# WGCNA and IRanges); set this before calling:
#   conflicted::conflicts_prefer(stats::cor)
# ============================================================

ConstructTFNetwork2 <- function(
        seurat_obj,
        model_params,
        nfold = 5,
        callbacks = NULL,
        wgcna_name = NULL
){

    if(is.null(wgcna_name)){ wgcna_name <- seurat_obj@misc$active_wgcna }
    CheckWGCNAName(seurat_obj, wgcna_name)

    # get the motif information from the Seurat object
    motif_matrix <- GetMotifMatrix(seurat_obj)
    motif_df <- GetMotifs(seurat_obj)

    if(is.null(motif_df)){
        stop("Motif info not found in the seurat_obj. Please run MotifScan first.")
    }

    if(! 'motif_name' %in% colnames(motif_df)){
        stop('motif_name column missing in motif table (GetMotifs(seurat_obj)).')
    }

    # detect xgboost 2.x (changed model accessor and callback names)
    check_xgboost_new <- packageVersion('xgboost') >= '2.0.0'

    if(is.null(callbacks)){
        if(check_xgboost_new){
            callbacks <- list(xgboost::xgb.cb.cv.predict(save_models = TRUE))
        } else {
            callbacks <- list(xgboost::cb.cv.predict(save_models = TRUE))
        }
    }

    # subset the motif_df to motifs whose TF is in the Seurat object
    motif_df <- subset(motif_df, motif_name %in% rownames(seurat_obj))

    # expression matrix
    datExpr <- as.matrix(GetDatExpr(seurat_obj, wgcna_name = wgcna_name))
    genes_use <- colnames(datExpr)
    genes_use <- genes_use[genes_use %in% rownames(motif_matrix)]

    importance_df <- data.frame()
    eval_df <- data.frame()

    pb <- utils::txtProgressBar(min = 0, max = length(genes_use),
                                style = 3, width = 50, char = "=")
    counter <- 1

    for(cur_gene in genes_use){

        setTxtProgressBar(pb, counter)

        # TFs whose motif is near this gene
        cur_tfs <- names(which(motif_matrix[cur_gene, ]))
        cur_tfs <- subset(motif_df, motif_ID %in% cur_tfs) %>% .$motif_name %>% unique
        cur_tfs <- cur_tfs[cur_tfs %in% genes_use]

        if(cur_gene %in% cur_tfs){
            cur_tfs <- cur_tfs[cur_tfs != cur_gene]
        }
        x_vars <- datExpr[, cur_tfs]
        y_var <- as.numeric(datExpr[, cur_gene])

        if(length(cur_tfs) < 2){
            print(paste0('Not enough putative TFs, skipping ', cur_gene))
            next
        }

        tf_cor <- as.numeric(cor(x = as.matrix(x_vars), y = y_var))
        names(tf_cor) <- cur_tfs

        if(all(y_var == 0)){
            print(paste0('skipping ', cur_gene))
            next
        }

        xgb <- xgboost::xgb.cv(
            params = model_params,
            data = xgboost::xgb.DMatrix(x_vars, label = y_var),
            nrounds = 100,
            showsd = FALSE,
            nfold = nfold,
            callbacks = callbacks,
            verbose = FALSE
        )

        xgb_eval <- as.data.frame(xgb$evaluation_log)
        xgb_eval$variable <- cur_gene

        # average importance across folds (xgboost 2.x model accessor)
        importance <- Reduce('+', lapply(1:nfold, function(i){
            if(check_xgboost_new){
                cur_model <- xgb$cv_predict$models[[i]]
            } else {
                cur_model <- xgb$models[[i]]
            }
            cur_imp <- xgboost::xgb.importance(
                feature_names = colnames(x_vars), model = cur_model)
            ix <- match(colnames(x_vars), as.character(cur_imp$Feature))
            cur_imp <- as.matrix(cur_imp[ix, -1])
            cur_imp[is.na(cur_imp)] <- 0
            cur_imp
        })) / nfold
        importance <- as.data.frame(importance)

        importance$tf <- colnames(x_vars)
        importance$gene <- cur_gene
        importance$Cor <- as.numeric(tf_cor)

        importance <- importance %>%
            dplyr::select(c(tf, gene, Gain, Cover, Frequency, Cor))
        importance <- dplyr::arrange(importance, -Gain)

        importance_df <- rbind(importance_df, importance)
        eval_df <- rbind(eval_df, xgb_eval)
        counter <- counter + 1
    }
    close(pb)

    seurat_obj <- SetTFNetwork(seurat_obj, importance_df, wgcna_name = wgcna_name)
    seurat_obj <- SetTFEval(seurat_obj, eval_df, wgcna_name = wgcna_name)
    seurat_obj
}
