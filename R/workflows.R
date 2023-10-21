#' @import dplyr
#' @import stringr
#' @import glue
#' @import ggplot2
#' @importFrom methods is as
#' @import patchwork
#' @importFrom grDevices colorRampPalette
#' @importFrom data.table as.data.table
#' @importFrom stats dnbinom na.omit optim pnorm setNames start
#' @useDynLib hahmmr
NULL


########################### Preprocessing ############################


#' Produce combined bulk expression and allele profile
#'
#' @param count_mat matrix Gene expression counts
#' @param lambdas_ref matrix Reference expression profiles
#' @param df_allele dataframe Allele counts
#' @param gtf dataframe Transcript gtf
#' @param min_depth integer Minimum coverage to filter SNPs
#' @param nu numeric Phase switch rate
#' @param genetic_map dataframe Genetic map
#' @param verbose logical Whether to print progress
#' @return dataframe Pseudobulk gene expression and allele profile
#' @export
get_bulk = function(count_mat, lambdas_ref, df_allele, gtf, genetic_map = NULL, min_depth = 0, nu = 1, verbose = TRUE) {

    count_mat = check_matrix(count_mat)

    fit = fit_ref_sse(rowSums(count_mat), lambdas_ref, gtf)

    exp_bulk = get_exp_bulk(
            count_mat,
            fit$lambdas_bar,
            gtf,
            verbose = verbose
        ) %>%
        filter((logFC < 5 & logFC > -5) | Y_obs == 0) %>%
        mutate(sse = fit$sse)

    if (verbose) {
        message(paste0('Fitted weights: ', paste0(signif(fit$w, 2), collapse = ',')))
    }

    allele_bulk = get_allele_bulk(
        df_allele,
        gtf,
        genetic_map,
        nu = nu,
        min_depth = min_depth)
            
    bulk = combine_bulk(
        allele_bulk = allele_bulk,
        exp_bulk = exp_bulk)

    if (any(duplicated(bulk$snp_id))) {
        stop('duplicated SNPs, check genotypes')
    }

    # doesn't work with 0s in the ref
    # TODO: should just make gene NA so we don't miss SNPs in that gene
    bulk = bulk %>% mutate(
        gene = ifelse(lambda_ref == 0, NA, gene),
        Y_obs = ifelse(lambda_ref == 0, NA, Y_obs),
        lambda_ref = ifelse(lambda_ref == 0, NA, lambda_ref)
    )

    # reannotate lost genes after joining
    gene_dict = allele_bulk %>% filter(!is.na(gene)) %>% 
        {setNames(.$gene, .$snp_id)} 

    bulk = bulk %>% 
        mutate(gene = ifelse(snp_id %in% names(gene_dict), gene_dict[as.character(snp_id)], gene)) %>% 
        select(-any_of(c('gene_start', 'gene_end', 'gene_length'))) %>%
        left_join(
            gtf %>% select(gene, gene_start, gene_end),
            by = 'gene'
        ) %>% 
        mutate(gene_index = as.integer(factor(gene, unique(gene[gene != '' | is.na(gene)]))))

    bulk = bulk %>%
        mutate(CHROM = as.character(CHROM)) %>%
        mutate(CHROM = ifelse(CHROM == 'X', 23, CHROM)) %>%
        mutate(CHROM = factor(as.integer(CHROM)))

    return(bulk)
}

#' Check the format of a count matrix
#' @param count_mat matrix Count matrix
#' @return matrix Count matrix
#' @keywords internal
check_matrix = function(count_mat) {

    # Make sure that the count matrix is of type integer
    if (!is.numeric(count_mat)) {
        msg = "The parameter 'count_mat' should be of type 'integer'. Please fix."
        stop(msg)
    } else if (all(count_mat != as.integer(count_mat))) {
        msg = "The parameter 'count_mat' should be of type 'integer'. Please fix."
        stop(msg)
    } else if (any(duplicated(rownames(count_mat)))) {
        msg = "Please remove duplicated genes in count matrix"
        stop(msg)
    }

    return(count_mat)
}

#' Aggregate into pseudobulk alelle profile
#' @param df_allele dataframe Single-cell allele counts
#' @param gtf dataframe Transcript gtf
#' @param genetic_map dataframe Genetic map
#' @param nu numeric Phase switch rate
#' @param min_depth integer Minimum coverage to filter SNPs
#' @return dataframe Pseudobulk allele profile
#' @export
get_allele_bulk = function(df_allele, gtf, genetic_map = NULL, nu = 0.5, min_depth = 0) {

    overlap_transcript = GenomicRanges::findOverlaps(
            df_allele %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$POS,
                    end = .$POS)
            )},
            gtf %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$gene_start,
                    end = .$gene_end)
            )}
        ) %>%
        as.data.frame() %>%
        setNames(c('snp_index', 'gene_index')) %>%
        left_join(
            df_allele %>% mutate(snp_index = 1:n()) %>%
                select(snp_index, snp_id),
            by = c('snp_index')
        ) %>%
        left_join(
            gtf %>% mutate(gene_index = 1:n()),
            by = c('gene_index')
        ) %>%
        arrange(snp_index, gene) %>%
        distinct(snp_index, `.keep_all` = T)

    df_allele = df_allele %>%
        select(-any_of(c('gene', 'gene_start', 'gene_end'))) %>% 
        left_join(
            overlap_transcript %>% select(snp_id, gene, gene_start, gene_end),
            by = c('snp_id')
        )

    df_allele = df_allele %>%
        mutate(AR = AD/DP) %>%
        filter(GT %in% c('1|0', '0|1')) %>%
        arrange(CHROM, POS) %>%
        group_by(CHROM) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ungroup() %>%
        filter(DP >= min_depth) %>%
        mutate(
            pBAF = ifelse(GT == '1|0', AR, 1-AR),
            pAD = ifelse(GT == '1|0', AD, DP - AD),
            # R = ifelse(pBAF == AR, -1, 1),
            R = ifelse(GT == '1|0', -1, 1)
        ) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM))) %>%
        filter(!(CHROM == 6 & POS < 33480577 & POS > 28510120)) %>%
        arrange(CHROM, POS) %>%
        annot_cm(genetic_map) %>%
        group_by(CHROM) %>%
        filter(n() > 1) %>%
        mutate(
            inter_snp_cm = c(NA, cM[2:length(cM)] - cM[1:(length(cM)-1)]),
            p_s = switch_prob_cm(inter_snp_cm, nu = nu)
        ) %>%
        ungroup() %>%
        mutate(gene = ifelse(gene == '', NA, gene))

    df_allele = df_allele %>% mutate(logFC = 0) %>% filter(DP > 0)

    return(df_allele)
}

#' Aggregate into bulk expression profile
#' @param count_mat dgCMatrix Gene expression counts
#' @param lambdas_ref matrix Reference expression profiles
#' @param gtf dataframe Transcript gtf
#' @return dataframe Pseudobulk gene expression profile
#' @keywords internal
get_exp_bulk = function(count_mat, lambdas_ref, gtf, verbose = FALSE) {

    depth_obs = sum(count_mat)
    
    mut_expressed = filter_genes(count_mat, lambdas_ref, gtf)
    count_mat = count_mat[mut_expressed,,drop=FALSE]
    lambdas_ref = lambdas_ref[mut_expressed]

    bulk_obs = count_mat %>%
        rowSums() %>%
        data.frame() %>%
        setNames('Y_obs') %>%
        tibble::rownames_to_column('gene') %>%
        mutate(lambda_obs = (Y_obs/depth_obs)) %>%
        mutate(lambda_ref = lambdas_ref[gene]) %>%
        mutate(d_obs = depth_obs) %>%
        left_join(gtf, by = "gene") 
    
    # annotate using GTF
    bulk_obs = bulk_obs %>%
        mutate(gene = droplevels(factor(gene, gtf$gene))) %>%
        mutate(gene_index = as.integer(gene)) %>%
        arrange(gene) %>%
        mutate(gene = as.character(gene)) %>%
        mutate(CHROM = factor(CHROM)) %>%
        mutate(
            logFC = log2(lambda_obs) - log2(lambda_ref),
            lnFC = log(lambda_obs) - log(lambda_ref),
            logFC = ifelse(is.infinite(logFC), NA, logFC),
            lnFC = ifelse(is.infinite(lnFC), NA, lnFC)
        )
    
    return(bulk_obs)
}


#' filter for mutually expressed genes
#' @param count_mat dgCMatrix Gene expression counts
#' @param lambdas_ref named numeric vector A reference expression profile
#' @param gtf dataframe Transcript gtf
#' @return vector Genes that are kept after filtering
#' @keywords internal
filter_genes = function(count_mat, lambdas_ref, gtf, verbose = FALSE) {

    genes_keep = gtf$gene %>% 
        intersect(rownames(count_mat)) %>%
        intersect(names(lambdas_ref))

    # https://www.sciencedirect.com/science/article/pii/S1357272520301990
    genes_exclude = gtf %>%
        filter(CHROM == 6 & gene_start < 33480577 & gene_end > 28510120) %>%
        pull(gene)

    genes_keep = genes_keep[!genes_keep %in% genes_exclude]

    count_mat = count_mat[genes_keep,,drop=FALSE]
    lambdas_ref = lambdas_ref[genes_keep]
    lambdas_obs = rowSums(count_mat)/sum(count_mat)

    min_both = 2

    mut_expressed = ((lambdas_ref * 1e6 > min_both & lambdas_obs * 1e6 > min_both) |
        (lambdas_ref > mean(lambdas_ref[lambdas_ref != 0])) |
        (lambdas_obs > mean(lambdas_obs[lambdas_obs != 0]))) &
        (lambdas_ref > 0)

    retained = names(mut_expressed)[mut_expressed]

    if (verbose) {
        message(glue('number of genes left: {length(retained)}'))
    }

    return(retained)
}

#' Annotate genetic distance between markers
#' @param bulk dataframe Pseudobulk profile
#' @param genetic_map dataframe Genetic map
#' @return dataframe Annotated pseudobulk profile
#' @keywords internal
annot_cm = function(bulk, genetic_map) {

    if ('cM' %in% colnames(bulk)) {
        return(bulk)
    } else {
        if (is.null(genetic_map)) {
            stop('Genetic map needs to be provided if cM is not in annotated')
        }
    }

    bulk = bulk %>% ungroup()
    
    marker_map = GenomicRanges::findOverlaps(
            bulk %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$POS,
                    end = .$POS)
            )},
            genetic_map %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$start,
                    end = .$end)
            )}
        ) %>%
        as.data.frame() %>%
        setNames(c('marker_index', 'map_index')) %>%
        left_join(
            bulk %>% mutate(marker_index = 1:n()) %>%
                select(marker_index, snp_id),
            by = c('marker_index')
        ) %>%
        left_join(
            genetic_map %>% mutate(map_index = 1:n()),
            by = c('map_index')
        ) %>%
        arrange(marker_index, -start) %>%
        distinct(marker_index, `.keep_all` = T) %>%
        select(snp_id, cM)

    bulk = bulk %>% left_join(marker_map, by = 'snp_id') %>%
        filter(!is.na(cM))

    return(bulk)
    
}

#' predict phase switch probablity as a function of genetic distance
#' @param d numeric vector Genetic distance in cM
#' @param nu numeric Phase switch rate
#' @param min_p numeric Minimum phase switch probability 
#' @return numeric vector Phase switch probability
#' @keywords internal
switch_prob_cm = function(d, nu = 1, min_p = 1e-10) {

    if (nu == 0) {
        p = rep(0, length(d))
    } else {
        p = (1-exp(-2*nu*d))/2
        p = pmax(p, min_p)
    }

    p = ifelse(is.na(d), 0, p)

    return(p)
}

#' Fit a reference profile from multiple references using constrained least square
#' @param Y_obs vector 
#' @param lambdas_ref matrix 
#' @param gtf dataframe 
#' @return fitted expression profile
#' @keywords internal
fit_ref_sse = function(Y_obs, lambdas_ref, gtf, min_lambda = 2e-6, verbose = FALSE) {

    if (any(duplicated(rownames(lambdas_ref)))) {
        stop('Duplicated genes in lambdas_ref')
    }

    if (length(dim(lambdas_ref)) == 1 | is.null(dim(lambdas_ref))) {
        return(list('w' = 1, 'lambdas_bar' = lambdas_ref))
    }
    
    Y_obs = Y_obs[Y_obs > 0]

    # take the union of expressed genes across cell type
    genes_common = gtf$gene %>% 
        intersect(names(Y_obs)) %>%
        intersect(rownames(lambdas_ref)[rowMeans(lambdas_ref) > min_lambda])

    if (verbose) {
        message(glue('{length(genes_common)} genes common in reference and observation'))
    }

    Y_obs = Y_obs[genes_common]
    lambdas_obs = Y_obs/sum(Y_obs)
    lambdas_ref = lambdas_ref[genes_common,,drop=FALSE]

    n_ref = ncol(lambdas_ref)
    
    fit = optim(
        fn = function(w) {
            w = w/sum(w)
            sum(log(lambdas_obs/as.vector(lambdas_ref %*% w))^2)
        },
        method = 'L-BFGS-B',
        par = rep(1/n_ref, n_ref),
        lower = rep(1e-6, n_ref)
    )

    w = fit$par
    w = w/sum(w)
    w = setNames(w, colnames(lambdas_ref))

    lambdas_bar = lambdas_ref %*% w %>% {setNames(as.vector(.), rownames(.))}

    return(list('w' = w, 'lambdas_bar' = lambdas_bar, 'sse' = fit$value/length(Y_obs)))
}

#' Combine allele and expression pseudobulks
#' @param allele_bulk dataframe Bulk allele profile
#' @param exp_bulk dataframe Bulk expression profile
#' @param gtf dataframe Transcript gtf
#' @return dataframe Pseudobulk allele and expression profile
#' @keywords internal
combine_bulk = function(allele_bulk, exp_bulk, gtf) {

    bulk = allele_bulk %>% 
        select(-any_of(c("gene_start", "gene_end", "logFC"))) %>%
        full_join(
            exp_bulk,
            by = c("CHROM", "gene")
        ) %>%
        mutate(
            snp_id = ifelse(is.na(snp_id), gene, snp_id),
            POS = ifelse(is.na(POS), gene_start, POS),
            # phase switch is forbidden if not heteroSNP
            p_s = ifelse(is.na(p_s), 0, p_s)
        ) %>%
        arrange(CHROM, POS) %>%
        filter(!(CHROM == 6 & POS < 33480577 & POS > 28510120)) %>%
        group_by(CHROM) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ungroup()

    # get rid of duplicate gene expression values
    bulk = bulk %>% 
        group_by(CHROM, gene) %>%
        mutate(
            Y_obs = ifelse(
                !is.na(gene) & n() > 1,
                c(unique(Y_obs), rep(NA, n()-1)), Y_obs
            )
        ) %>%
        ungroup() %>% 
        mutate(
            lambda_obs = Y_obs/d_obs,
            logFC = log2(lambda_obs/lambda_ref),
            lnFC = log(lambda_obs/lambda_ref)
        ) %>%
        mutate_at(
            c('logFC', 'lnFC'),
            function(x) ifelse(is.infinite(x), NA, x)
        )
    
    return(bulk)
    
}

########################### Wrkflows for CNV detection ############################

#' Check if columns are present in dataframe
#' @param df dataframe Dataframe 
#' @param cols character vector Column names
#' @return dataframe Dataframe 
#' @keywords internal
check_cols = function(df, cols) {

    if (!all(cols %in% colnames(df))) {
        stop(glue('Missing columns: {paste0(setdiff(cols, colnames(df)), collapse = ", ")}'))
    }

    return(df)

}

#' Analyze allele profile
#' @param bulk dataframe Bulk allele profile
#' @param t numeric Transition probability
#' @param theta_min numeric Minimum allele fraction
#' @param gamma numeric Overdispersion parameter
#' @param nu numeric Phase switch rate
#' @param r numeric Alternative allele count bias
#' @param hmm character HMM model to use (S3 or S5)
#' @param fit_theta logical Whether to fit theta_min
#' @param fit_gamma logical Whether to fit gamma
#' @param theta_start numeric Starting value for theta_min
#' @param verbose logical Whether to print progress
#' @return dataframe Bulk allele profile with CNV states
#' @export
analyze_allele = function(bulk, t = 1e-5, theta_min = 0.08, gamma = 20, nu = 0.5, r = 0.015, 
    hmm = 'S5', fit_theta = FALSE, fit_gamma = FALSE, theta_start = 0.05, verbose = TRUE) {

    required_cols = c('CHROM', 'POS', 'REF', 'ALT', 'snp_id', 'snp_index', 'GT', 'pAD', 'DP', 'R', 'gene', 'inter_snp_cm', 'p_s')

    bulk = check_cols(bulk, required_cols)

    bulk = bulk %>% mutate(pBAF = pAD/DP) %>% filter(!is.na(pBAF))

    # update transition probablity
    bulk = bulk %>% 
        mutate(p_s = switch_prob_cm(inter_snp_cm, nu = UQ(nu)))
    
    if (fit_gamma) {
        gamma = bulk %>% 
            filter(!is.na(AD)) %>%
            {fit_gamma(.$AD, .$DP, r = r)}

        if (verbose) {
            message(glue('Fitted gamma: {gamma}'))
        }
    }

    if (fit_theta) {
        if (verbose) {message('Fitting theta_min ..')}
        bulk = bulk %>%
            group_by(CHROM) %>%
            mutate(
                approx_theta_post_s3(pAD, DP, R, p_s, t = t, gamma = gamma, r = r, start = theta_start)
            ) %>%
            ungroup()
    } else {
        bulk$theta_min = theta_min
    }

    if (hmm == 'S3') {
        run_hmm = run_allele_hmm_s3
    } else if (hmm == 'S5') {
        run_hmm = run_allele_hmm_s5
    } else {
        stop('hmm must be S3 or S5')
    }
    
    bulk = bulk %>%
        group_by(CHROM) %>%
        mutate(
            state = run_hmm(pAD = pAD, DP = DP, R = R, p_s = p_s, theta_min = unique(theta_min), gamma = UQ(gamma), r = r),
            cnv_state = str_remove(state, '_up|_down')
        ) %>%
        mutate(
            haplo_theta_min = case_when(
                str_detect(state, 'up') ~ 'major',
                str_detect(state, 'down') ~ 'minor',
                TRUE ~ ifelse(pBAF > 0.5, 'major', 'minor')
            ),
            major_count = ifelse(haplo_theta_min == 'major', pAD, DP - pAD),
            minor_count = DP - major_count
        ) %>%
        ungroup() %>%
        annot_segs() %>%
        smooth_segs() %>% 
        annot_segs() %>%
        mutate(seg_size = seg_end - seg_start)
    
    segs_retest = bulk %>%
        filter(cnv_state != 'neu') %>%
        group_by(seg) %>%
        summarise(
            theta_hat = theta_hat_seg(sum(major_count), sum(minor_count)),
            approx_theta_post_s2(pAD, DP, R, p_s, gamma = UQ(gamma), start = unique(theta_hat)),
            LLR = calc_allele_LLR(pAD, DP, R, p_s, theta_mle, gamma = UQ(gamma), r = r)
        )

    bulk = bulk %>% 
        select(-any_of(c('LLR', 'theta_hat', 'theta_mle', 'theta_sigma'))) %>%
        left_join(segs_retest, by = 'seg')

    # get posterior haplotype counts
    bulk = bulk %>% mutate(cnv_state_post = cnv_state) %>%
        classify_alleles() 
    
    return(bulk)
}

#' Analyze allele and expression profile
#' @param bulk dataframe Bulk allele and expression profile
#' @param t numeric Transition probability
#' @param theta_min numeric Minimum allele fraction
#' @param gamma numeric Overdispersion parameter
#' @param logphi_min numeric Minimum log2 fold change
#' @param hmm character HMM model to use (S7 or S15)
#' @param min_genes integer Minimum number of genes per segment
#' @param exclude_neu logical Whether to exclude neutral segments in retest
#' @param nu numeric Phase switch rate
#' @param r numeric Alternative allele count bias
#' @param fit_theta logical Whether to fit theta_min
#' @param fit_gamma logical Whether to fit gamma
#' @param theta_start numeric Starting value for theta_min
#' @param verbose logical Whether to print progress
#' @return dataframe Bulk allele and expression profile with CNV states
#' @export
analyze_joint = function(
    bulk, t = 1e-5, gamma = 20, theta_min = 0.08, logphi_min = 0.25, hmm = 'S15',
    nu = 1, min_genes = 10, r = 0.015, theta_start = 0.05, exclude_neu = TRUE, fit_gamma = FALSE, fit_theta = FALSE, verbose = TRUE
) {

    required_cols = c('CHROM', 'POS', 'REF', 'ALT', 'snp_id', 'snp_index', 'GT', 'AD', 'pAD', 'DP', 'R', 'Y_obs', 'logFC', 
        'd_obs', 'lambda_ref', 'gene', 'inter_snp_cm', 'p_s')

    bulk = check_cols(bulk, required_cols)

    bulk = bulk %>% mutate(pBAF = pAD/DP)

    # update transition probablity
    bulk = bulk %>% 
        filter(DP > 0 | is.na(DP)) %>%
        mutate(p_s = switch_prob_cm(inter_snp_cm, nu = UQ(nu)))

    bulk = find_diploid(bulk, gamma = gamma, t = t, theta_min = theta_min, r = r)

    # fit expression baseline
    fit = bulk %>%
        filter(!is.na(Y_obs)) %>%
        filter(logFC < 8 & logFC > -8) %>%
        filter(diploid) %>%
        {fit_lnpois_cpp(.$Y_obs, .$lambda_ref, unique(.$d_obs))}
        
    bulk = bulk %>% mutate(mu = fit[1], sig = fit[2])

    if (verbose) {
        message(glue('Fitted sigma: {fit[2]}'))
    }

    # fit allele overdispersion
    if (fit_gamma) {
        gamma = bulk %>% 
            filter(diploid) %>%
            filter(!is.na(AD)) %>%
            {fit_gamma(.$AD, .$DP, r = r)}

        if (verbose) {
            message(glue('Fitted gamma: {gamma}'))
        }
    }

    if (fit_theta) {
        if (verbose) {message('Fitting theta_min ..')}
        bulk = bulk %>%
            group_by(CHROM) %>%
            mutate(
                approx_theta_post_s3(pAD[!is.na(pAD)], DP[!is.na(pAD)], R[!is.na(pAD)], p_s[!is.na(pAD)], 
                    t = t, gamma = gamma, r = r, start = theta_start)
            ) %>%
            ungroup()
    } else {
        bulk$theta_min = theta_min
    }

    if (hmm == 'S7') {
        run_hmm = run_joint_hmm_s7
    } else if (hmm == 'S15') {
        run_hmm = run_joint_hmm_s15
    } else {
        stop('hmm must be either S7 or S15')
    }
    
    # run joint HMM
    bulk = bulk %>% 
        group_by(CHROM) %>%
        mutate(state = 
            run_hmm(
                pAD = pAD,
                DP = DP, 
                R = R,
                p_s = p_s,
                Y_obs = Y_obs, 
                lambda_ref = lambda_ref, 
                d_total = na.omit(unique(d_obs)),
                phi_amp = 2^(logphi_min),
                phi_del = 2^(-logphi_min),
                mu = mu,
                sig = sig,
                t = t,
                gamma = UQ(gamma),
                theta_min = unique(theta_min),
                r = r
            )
        ) %>%
        mutate(cnv_state = str_remove(state, '_down|_up')) %>%
        annot_segs(var = 'cnv_state') %>%
        smooth_segs(min_genes = min_genes) %>%
        mutate(cnv = ifelse(cnv_state == 'neu', 0, 1)) %>%
        annot_segs(var = 'cnv') %>%
        mutate(state = ifelse(cnv_state == 'neu', 'neu', state)) %>% 
        mutate(cnv = ifelse(cnv_state == 'neu', 0, 1))

    # annotate phi MLE
    bulk = bulk %>%
        group_by(seg) %>%
        mutate(
            approx_phi_post(
                Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)),
                mu = mu[!is.na(Y_obs)],
                sig = sig[!is.na(Y_obs)]
            ),
            p_amp = 1-pnorm(1, mean = phi_mle, sd = phi_mle_sig),
            p_del = 1-p_amp
        ) %>%
        ungroup()
    
    # annotate theta MLE
    segs_retest = bulk %>%
        filter((cnv_state != 'neu') | (!exclude_neu)) %>%
        group_by(seg) %>%
        summarise(
            approx_theta_post_s2(pAD[!is.na(pAD)], DP[!is.na(pAD)], R[!is.na(pAD)], p_s[!is.na(pAD)], gamma = UQ(gamma), r = r, start = 0.1),
            LLR_y = calc_allele_LLR(pAD[!is.na(pAD)], DP[!is.na(pAD)], R[!is.na(pAD)], p_s[!is.na(pAD)], theta_mle, gamma = UQ(gamma), r = r),
            LLR_x = calc_exp_LLR(
                Y_obs[!is.na(Y_obs)], 
                lambda_ref[!is.na(Y_obs)], 
                unique(na.omit(d_obs)), 
                unique(phi_mle), 
                mu = mu[!is.na(Y_obs)], 
                sig = sig[!is.na(Y_obs)]),
            LLR = LLR_x + LLR_y
        )

    bulk = bulk %>% 
        select(-any_of(c('LLR', 'theta_hat', 
            'theta_mle', 'theta_mle_sig', 
            'theta_map', 'theta_map_sig', 
            'LLR_y', 'LLR_x'))) %>%
        left_join(segs_retest, by = 'seg')

    bulk = bulk %>% mutate(
            cnv_state_post = case_when(
                cnv_state == 'neu' ~ 'neu',
                p_amp > 0.999 ~ 'amp',
                p_del > 0.999 ~ 'del',
                TRUE ~ 'loh' 
            ),
            state_post = state
        )

    bulk = bulk %>% classify_alleles()

    # store these info here
    bulk$nu = nu 
    bulk$gamma = gamma

    return(bulk)
}

#' Annotate copy number segments after HMM decoding 
#' @param bulk dataframe Pseudobulk profile
#' @return a pseudobulk dataframe
#' @keywords internal
annot_segs = function(bulk, var = 'cnv_state') {

    bulk = bulk %>% 
            group_by(CHROM) %>%
            arrange(CHROM, snp_index) %>%
            mutate(boundary = c(0, get(var)[2:length(get(var))] != get(var)[1:(length(get(var))-1)])) %>%
            group_by(CHROM) %>%
            mutate(seg = paste0(CHROM, generate_postfix(cumsum(boundary)+1))) %>%
            arrange(CHROM) %>%
            mutate(seg = factor(seg, unique(seg))) %>%
            ungroup() %>%
            group_by(seg) %>%
            mutate(
                seg_start = min(POS),
                seg_end = max(POS),
                seg_start_index = min(snp_index),
                seg_end_index = max(snp_index),
                n_genes = length(na.omit(unique(gene))),
                n_snps = sum(!is.na(pAD))
            ) %>%
            ungroup()

    return(bulk)
}


#' Generate alphabetical postfixes
#' @param n vector of integers
#' @return vector of alphabetical postfixes
#' @keywords internal
generate_postfix <- function(n) {

    if (any(is.na(n))) {
        stop("Segment number cannot contain NA")
    }

    alphabet <- letters

    postfixes <- sapply(n, function(i) {
        postfix <- character(0)
        while (i > 0) {
            remainder <- (i - 1) %% 26
            i <- (i - 1) %/% 26
            postfix <- c(alphabet[remainder + 1], postfix)
        }
        paste(postfix, collapse = "")
    })

    return(postfixes)
}

#' Smooth the segments after HMM decoding 
#' @param bulk dataframe Pseudobulk profile
#' @param min_genes integer Minimum number of genes to call a segment
#' @return dataframe Pseudobulk profile
#' @keywords internal
smooth_segs = function(bulk, min_genes = 10) {

    bulk = bulk %>% group_by(seg) %>%
        mutate(
            cnv_state = ifelse(n_genes <= min_genes, NA, cnv_state)
        ) %>%
        ungroup() %>%
        group_by(CHROM) %>%
        mutate(cnv_state = zoo::na.locf(cnv_state, fromLast = FALSE, na.rm=FALSE)) %>%
        mutate(cnv_state = zoo::na.locf(cnv_state, fromLast = TRUE, na.rm=FALSE)) %>%
        ungroup()

    chrom_na = bulk %>% group_by(CHROM) %>% summarise(all_na = all(is.na(cnv_state)))

    if (any(chrom_na$all_na)) {
        chroms_na = paste0(chrom_na$CHROM[chrom_na$all_na], collapse = ',')
        msg = glue("No segments containing more than {min_genes} genes for CHROM {chroms_na}.")
        stop(msg)
    }

    return(bulk)
}

#' Estimate of imbalance level theta in a segment
#' @param major_count vector of major allele count
#' @param minor_count vector of minor allele count
#' @return estimate of theta
#' @keywords internal
theta_hat_seg = function(major_count, minor_count) {
    major_total = sum(major_count)
    minor_total = sum(minor_count)
    MAF = major_total/(major_total + minor_total)
    return(MAF - 0.5)
}

find_diploid = function(bulk, gamma = 20, r = 0.015, theta_min = 0.08, t = 1e-5, debug = FALSE, verbose = TRUE) {

    # define imbalanced regions
    bulk = bulk %>% 
        group_by(CHROM) %>%
        mutate(
            state = run_allele_hmm_s3(pAD, DP, R, p_s, theta_min = UQ(theta_min), gamma = UQ(gamma), r = r),
            cnv_state = str_remove(state, '_down|_up')
        ) %>% 
        ungroup()
    
    bulk = bulk %>% mutate(diploid = cnv_state == 'neu')

    return(bulk)
}


#' Laplace approximation of the posterior of expression fold change phi
#' @param Y_obs numeric vector Gene expression counts
#' @param lambda_ref numeric vector Reference expression levels
#' @param d numeric Total library size
#' @param alpha numeric Shape parameter of the gamma distribution
#' @param beta numeric Rate parameter of the gamma distribution
#' @param mu numeric Mean of the normal distribution
#' @param sig numeric Standard deviation of the normal distribution
#' @param lower numeric Lower bound of phi
#' @param upper numeric Upper bound of phi
#' @param start numeric Starting value of phi
#' @return numeric MLE of phi and its standard deviation
#' @keywords internal
approx_phi_post = function(Y_obs, lambda_ref, d, alpha = NULL, beta = NULL, mu = NULL, sig = NULL, lower = 0.2, upper = 10, start = 1) {
    
    if (length(Y_obs) == 0) {
        return(tibble('phi_mle' = 1, 'phi_sigma' = 0))
    }
    
    start = max(min(1, upper), lower)

    l = function(phi) {l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi = phi)}

    fit = optim(
        start,
        function(phi) {
            -l(phi)
        },
        method = 'L-BFGS-B',
        lower = lower,
        upper = upper,
        hessian = TRUE
    )

    mean = fit$par
    sd = sqrt(as.numeric(1/(fit$hessian)))

    if (is.na(sd)) {
        mean = 1
        sd = 0
    }

    return(tibble('phi_mle' = mean, 'phi_mle_sig' = sd))
}

#' Find optimal theta of a chromosome using forward-backward; uses a 3-state allele HMM
#' @keywords internal
approx_theta_post_s3 = function(pAD, DP, R, p_s, t, upper = 0.45, start = 0.05, gamma = 20, r = 0.015) {

    gamma = unique(gamma)

    if (length(gamma) > 1) {
        stop('gamma has to be a single value')
    }
    
    if (length(pAD) <= 10) {
        return(tibble('theta_mle' = 0, 'theta_sigma' = 0))
    }

    theta_grid = seq(start, 0.2, 0.01)

    liks = sapply(
        theta_grid,
        function(theta) {
            calc_allele_lik_s3(pAD, DP, R, p_s, t, abs(theta), gamma = gamma, r = r)
        }
    )

    theta_start = theta_grid[which.max(liks)]
    
    fit = optim(
        theta_start, 
        function(theta) {-calc_allele_lik_s3(pAD, DP, R, p_s, t, abs(theta), gamma = gamma, r = r)},
        method = 'L-BFGS-B',
        lower = -upper,
        upper = upper,
        hessian = FALSE
    )
        
    return(tibble('theta_min' = abs(fit$par)))
}


#' Find optimal theta of a CNV segment using forward-backward; uses a 2-state allele HMM
#' @keywords internal
approx_theta_post_s2 = function(pAD, DP, R, p_s, upper = 0.499, start = 0.25, gamma = 20, r = 0.015) {

    gamma = unique(gamma)

    if (length(gamma) > 1) {
        stop('gamma has to be a single value')
    }
    
    if (length(pAD) <= 10) {
        return(tibble('theta_mle' = 0, 'theta_sigma' = 0))
    }

    fit = optim(
        start, 
        function(theta) {-calc_allele_lik_s2(pAD, DP, R, p_s, abs(theta), gamma = gamma, r = r)},
        method = 'L-BFGS-B',
        lower = -upper,
        upper = upper,
        hessian = TRUE
    )
    
    mu = abs(fit$par)
    sigma = sqrt(as.numeric(1/(fit$hessian)))

    if (is.na(sigma)) {
        sigma = 0
    }
    
    return(tibble('theta_mle' = mu, 'theta_mle_sig' = sigma))
}


#' Calculate LLR for an allele HMM
#' @param pAD numeric vector Phased allele depth
#' @param DP numeric vector Total allele depth
#' @param R numeric vector Allelic bias direction
#' @param p_s numeric vector Phase switch probabilities
#' @param theta_mle numeric MLE of imbalance level theta (alternative hypothesis)
#' @param theta_0 numeric Imbalance level in the null hypothesis 
#' @param gamma numeric Dispersion parameter for the Beta-Binomial allele model
#' @return numeric Log-likelihood ratio
#' @keywords internal
calc_allele_LLR = function(pAD, DP, R, p_s, theta_mle, theta_0 = 0, gamma = 20, r = 0.015) {
    if (length(pAD) <= 1) {
        return(0)
    }
    l_1 = calc_allele_lik_s2(pAD, DP, R, p_s, theta = theta_mle, gamma = gamma, r = r) 
    AD = ifelse(R == -1, pAD, DP-pAD)
    l_0 = l_bbinom(AD, DP, gamma * (0.5 - r), gamma * (0.5 + r))
    return(l_1 - l_0)
}


#' Calculate LLR for an expression HMM
#' @param Y_obs numeric vector Gene expression counts
#' @param lambda_ref numeric vector Reference expression levels 
#' @param d numeric vector Total library size
#' @param phi_mle numeric MLE of expression fold change phi (alternative hypothesis)
#' @param mu numeric Mean parameter for the PLN expression model
#' @param sig numeric Dispersion parameter for the PLN expression model
#' @param alpha numeric Hyperparameter for the gamma poisson model (not used)
#' @param beta numeric Hyperparameter for the gamma poisson model (not used)
#' @return numeric Log-likelihood ratio
#' @keywords internal
calc_exp_LLR = function(Y_obs, lambda_ref, d, phi_mle, mu = NULL, sig = NULL, alpha = NULL, beta = NULL) {
    if (length(Y_obs) == 0) {
        return(0)
    }

    l_1 = l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi = phi_mle)
    l_0 = l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi = 1)
    
    return(l_1 - l_0)
}

#' fit gamma maximum likelihood
#' @param AD numeric vector Variant allele depth
#' @param DP numeric vector Total allele depth
#' @param r numeric Degree of allele bias
#' @return a fit
#' @keywords internal
fit_gamma = function(AD, DP, r = 0.015, start = 20) {
    fit = optim(
        par = start,
        fn = function(gamma) {
            -l_bbinom(AD, DP, gamma * (0.5 - r), gamma * (0.5 + r))
        },
        method = 'L-BFGS-B',
        lower = 1e-04,
        upper = 1e4
    )
    gamma = fit$par
    return(gamma)
}


#' classify alleles using viterbi and forward-backward
#' @param bulk dataframe Pesudobulk profile
#' @return dataframe Pesudobulk profile
#' @keywords internal
classify_alleles = function(bulk) {

    allele_bulk = bulk %>%
        filter(!cnv_state_post %in% c('neu')) %>%
        filter(!is.na(pAD)) %>% 
        group_by(CHROM, seg) %>%
        filter(n() > 1)

    if (nrow(allele_bulk) == 0) {
        return(bulk)
    }
    
    allele_post = allele_bulk %>%
        group_by(CHROM, seg) %>%
        mutate(
            p_up = forward_back_allele(get_allele_hmm_s2(pAD, DP, R, p_s, theta = unique(theta_mle), gamma = 20))[,1],
            haplo_post = case_when(
                p_up >= 0.5 & GT == '1|0' ~ 'major',
                p_up >= 0.5 & GT == '0|1' ~ 'minor',
                p_up < 0.5 & GT == '1|0' ~ 'minor',
                p_up < 0.5 & GT == '0|1' ~ 'major'
            )
        ) %>%
        ungroup() %>%
        select(snp_id, p_up, haplo_post)

    bulk = bulk %>% 
            select(-any_of(colnames(allele_post)[!colnames(allele_post) %in% c('snp_id')])) %>%
            left_join(
                allele_post,
                by = c('snp_id')
            )

    return(bulk)
}


## for tidyverse (magrittr & dplyr) functions 
if (getRversion() >= "2.15.1"){
  utils::globalVariables(c(".", "AD", "AD_all", "ALT", "AR", "CHROM", "DP", "DP_all", "FILTER", 
    "GT", "ID", "LLR", "LLR_sample","LLR_x", "LLR_y", "L_x_a", "L_x_d", "L_x_n", 
    "L_y_a", "L_y_d", "L_y_n", "MAF", "OTH", "OTH_all", "POS",
    "QUAL", "REF", "TP", "UQ", "Y_obs", "acen_hg19", "acen_hg38", "annot", 
    "branch", "cM", "chrom_sizes_hg19", "chrom_sizes_hg38", "clone",
    "clone_opt", "clone_size", "cluster", "cnv_state", "cnv_state_expand", "cnv_state_map",
    "cnv_state_post", "cnv_states", "compartment", "component", "cost", "d_obs", "diploid",
    "down", "edges", "end_x", "end_y", "exp_rollmean", "expected_colnames", "extract",
    "frac_overlap_x", "frac_overlap_y", "from", "from_label", "from_node", "gaps_hg19",
    "gaps_hg38", "gene", "gene_end", "gene_index", "gene_length", "gene_snps", "gene_start",
    "group", "groupOTU", "gtf_hg19", "gtf_hg38", "gtf_mm10", "haplo", "haplo_naive", "haplo_post",
    "haplo_theta_min", "het", "hom_alt", "i", "inter_snp_cm", "inter_snp_dist", 
     "j", "keep", "theta_hat_roll", "phi_mle_sig", 'R',
    "l31_x", "l31_y", "l_clone", "l_clone_x", "l_clone_y", "label", "lambda_obs", "lambda_ref", 
    "len_overlap", "len_x", "len_y", "lnFC", "lnFC_i", "lnFC_j", "lnFC_max_i", "lnFC_max_j",
    "logBF", "logBF_x", "logBF_y", "logFC", "loh", "major", "major_count", "marker_index", 
    "minor", "minor_count", "mu", "n_chrom_snp", "n_genes", "n_mut", 
    "n_snps", "n_states", "name", "node", "node_phylo", "nodes", "p", "pAD", "pBAF", "p_0",
    "p_1", "p_amp", "p_bamp", "p_bdel", "p_cnv", "p_del", "p_loh", "p_max", "p_n", "p_neu", "p_s", "p_up",
    "phi_mle", "phi_mle_roll", "phi_sigma", "pnorm.range", "potential_missing_columns",
    "precision", "prior_amp", "prior_bamp", "prior_bdel", "prior_clone", "prior_del",
    "prior_loh", "s", "seg", "seg_cons", "seg_end", "seg_end_index", "seg_label", "seg_length",
    "seg_start", "seg_start_index", "segs_consensus", "seqnames", "set_colnames", "sig",
    "site", "size", "snp_id", "snp_index", "snp_rate", "start_x", "start_y", "state", "state_post",
    "superclone", "theta_hat", "theta_level", "theta_mle", "theta_sigma", "to", "to_label",
    "to_node", "total", "value", "variable", "vcf_meta", "width", "x", "y"))
}