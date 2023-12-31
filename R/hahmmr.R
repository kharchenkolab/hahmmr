#########################################################################

# HaHMMR: Haplotype-aware Hidden Markov Models for CNV detection from RNA

#########################################################################


############ allele HMMs ############

#' Make a 2-state allele HMM - no transitions to netural state
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param p_s numeric vector Phase switch probabilities
#' @param theta numeric Haplotype imbalance
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @return HMM object
#' @keywords internal
get_allele_hmm_s2 = function(pAD, DP, R, p_s, theta, gamma = 20, r = 0.015) {

    states = c("theta_up", "theta_down")

    N = length(p_s)

    Pi = sapply(p_s, function(p_s) {c(1 - p_s, p_s, p_s, 1 - p_s)}) %>% 
        array(dim = c(2, 2, N))
    
    prior = c(0.5, 0.5)
    theta_states = outer(R*r, c(0.5 + theta, 0.5 - theta), FUN = "+")
    # theta_states = outer((1+2*R*r), c(0.5 + theta, 0.5 - theta), FUN = "*")
    theta_states = pmax(pmin(theta_states, 1),0)
    alpha_states = theta_states * gamma
    beta_states = (1 - theta_states) * gamma

    hmm = list(
        x = pAD, 
        logPi = log(Pi),
        delta = prior, 
        alpha = alpha_states, 
        beta = beta_states,
        d = DP,
        N = N,
        M = 2,
        K = 1,
        states = states
    )

    return(hmm)
}

#' Make a 3-state allele HMM - allow transitions to netural state
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param p_s numeric vector Phase switch probabilities
#' @param theta_min numeric Haplotype imbalance
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @param t numeric Transition probability between haplotype states
#' @param r numeric Variant mapping bias
#' @return HMM object
#' @keywords internal
get_allele_hmm_s3 = function(pAD, DP, R, p_s, t, theta_min, gamma = 20, r = 0.015) {

    gamma = unique(gamma)

    # states
    states = c("neu", "theta_up", "theta_down")

    N = length(pAD)
    M = length(states)

    # transition matrices
    calc_trans_mat_s3 = function(p_s, t) {
        matrix(
            c(1-t, t/2, t/2, 
             t, (1-t)*(1-p_s), (1-t)*p_s, 
             t, (1-t)*p_s, (1-t)*(1-p_s)),
            ncol = 3,
            byrow = TRUE
        )
    }

    Pi = sapply(p_s, function(p_s) {calc_trans_mat_s3(p_s, t)}) %>% 
        array(dim = c(M, M, N))

    theta_states = outer(R*r, c(0.5, 0.5 + theta_min, 0.5 - theta_min), FUN = "+")
    # theta_states = outer((1+2*R*r), c(0.5, 0.5 + theta_min, 0.5 - theta_min), FUN = "*")
    theta_states = pmax(pmin(theta_states, 1),0)
    alpha_states = theta_states * gamma
    beta_states = (1 - theta_states) * gamma
            
    hmm = list(
        x = pAD, 
        logPi = log(Pi), 
        delta = c(1-t, t/2, t/2), 
        alpha = alpha_states,
        beta = beta_states,
        d = DP,
        N = N,
        M = M,
        p_s = p_s,
        states = states
    )
    
    return(hmm)
}


#' Viterbi algorithm for allele HMM
#' @param hmm HMM object; expect variables x (allele depth), d (total depth),
#' logPi (log transition prob matrix), delta (prior for each state), 
#' alpha (alpha for each state), beta (beta for each state), 
#' states (states), p_s (phase switch probs)
#' @return character vector Decoded states
#' @keywords internal
viterbi_allele <- function(hmm) {

    N <- hmm$N
    M <- hmm$M
    nu <- matrix(NA, nrow = N, ncol = M)
    z <- rep(NA, N)

    logprob = sapply(1:M, function(m) {

        l_x = dbbinom(x = hmm$x, size = hmm$d, alpha = hmm$alpha[,m], beta = hmm$beta[,m], log = TRUE)

        l_x[is.na(l_x)] = 0

        return(l_x)

    })

    z = viterbi_compute(log(hmm$delta), logprob, hmm$logPi, N, M, nu, z)
        
    return(z)
}

#' Forward-backward algorithm for allele HMM
#' @param hmm HMM object; expect variables x (allele depth), d (total depth),
#' logPi (log transition prob matrix), delta (prior for each state), 
#' alpha (alpha for each state), beta (beta for each state), 
#' states (states), p_s (phase switch probs)
#' @return numeric matrix; posterior probabilities
#' @examples
#' forward_back_allele(pre_likelihood_hmm)
#' @export
forward_back_allele = function(hmm) {

    # case of one-data point
    if (hmm$N == 1) {
        return(NA)
    }

    N <- hmm$N
    M <- hmm$M
        
    logprob = sapply(1:M, function(m) {

        l_x = dbbinom(x = hmm$x, size = hmm$d, alpha = hmm$alpha[,m], beta = hmm$beta[,m], log = TRUE)

        l_x[is.na(l_x)] = 0

        return(l_x)

    })
        
    logphi <- log(hmm$delta)
        
    marginals = forward_backward_compute(logphi, logprob, hmm$logPi, N, M)

    colnames(marginals) = hmm$states

    return(marginals)
}

#' Only compute total log likelihood from an allele HMM
#' @param hmm HMM object; expect variables x (allele depth), d (total depth),
#' logPi (log transition prob matrix), delta (prior for each state), 
#' alpha (alpha for each state), beta (beta for each state), 
#' states (states), p_s (phase switch probs)
#' @return numeric; total log likelihood
#' @examples
#' likelihood_allele(pre_likelihood_hmm)
#' @export
likelihood_allele = function(hmm) {
        
    N <- hmm$N
    M <- hmm$M
        
    logprob = sapply(1:M, function(m) {

        l_x = dbbinom(x = hmm$x, size = hmm$d, alpha = hmm$alpha[,m], beta = hmm$beta[,m], log = TRUE)

        l_x[is.na(l_x)] = 0

        return(l_x)

    })
        
    logphi <- log(hmm$delta)
        
    LL <- likelihood_compute(logphi, logprob, hmm$logPi, N, M)

    return(LL)
}

#' Run a 5-state allele-only HMM - two theta levels
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param p_s numeric vector Phase switch probabilities
#' @param t numeric Transition probability between copy number states
#' @param theta_min numeric Minimum haplotype frequency deviation threshold
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @param prior numeric vector Prior probabilities for each state
#' @param ... Additional parameters
#' @return character vector Decoded states
#' @examples
#' with(bulk_example, {
#'     run_allele_hmm_s5(pAD = pAD, DP = DP, R = R, p_s = p_s, theta_min = 0.08, gamma = 30)
#' })
#' @export
run_allele_hmm_s5 = function(pAD, DP, p_s, t = 1e-5, theta_min = 0.08, gamma = 20, prior = NULL, ...) {

    gamma = unique(gamma)

    if (length(gamma) > 1) {
        stop('More than one gamma parameter')
    }

    # states
    states = c("neu", "theta_1_up", "theta_1_down", "theta_2_up", "theta_2_down")

    N = length(pAD)
    M = 5

    # transition matrices
    calc_trans_mat_s5 = function(p_s, t) {
        matrix(
            c(1-t, t/4, t/4, t/4, t/4, 
             t/2, (1-t)*(1-p_s), (1-t)*p_s, t/4, t/4,
             t/2, (1-t)*p_s, (1-t)*(1-p_s), t/4, t/4,
             t/2, t/4, t/4, (1-t)*(1-p_s), (1-t)*p_s,
             t/2, t/4, t/4, (1-t)*p_s, (1-t)*(1-p_s)),
            ncol = 5,
            byrow = TRUE
        )
    }

    Pi = sapply(p_s, function(p_s) {calc_trans_mat_s5(p_s, t)}) %>% 
        array(dim = c(M, M, N))

    # intitial probabilities
    if (is.null(prior)) {
        prior = rep(1/5, 5)
    }

    theta_1 = theta_min
    theta_2 = 0.4
    alphas = gamma * c(0.5, 0.5 + theta_1, 0.5 - theta_1, 0.5 + theta_2, 0.5 - theta_2)
    betas = gamma * c(0.5, 0.5 - theta_1, 0.5 + theta_1, 0.5 - theta_2, 0.5 + theta_2)
            
    hmm = list(
        x = pAD, 
        logPi = log(Pi), 
        delta = prior, 
        alpha = matrix(rep(alphas, N), ncol = M, byrow = TRUE),
        beta = matrix(rep(betas, N), ncol = M, byrow = TRUE),
        d = DP,
        N = N,
        M = M,
        states = states
    )
    
    mpc = states[viterbi_allele(hmm)]
    
    return(mpc)
}

#' Run a 3-state allele HMM
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param R numeric vector Variant mapping bias direction
#' @param p_s numeric vector Phase switch probabilities
#' @param t numeric Transition probability between copy number states
#' @param theta_min numeric Minimum haplotype frequency deviation threshold
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @return character vector Decoded states
#' @keywords internal
run_allele_hmm_s3 = function(pAD, DP, R, p_s, t = 1e-5, theta_min = 0.08, gamma = 20, r = 0.015) {
    hmm = get_allele_hmm_s3(pAD, DP, R, p_s, t = t, theta_min = theta_min, gamma = gamma, r = r)
    mpc = hmm$states[viterbi_allele(hmm)]
    return(mpc)
}

#' Calculate allele likelihoods for 2-state allele HMM
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param R numeric vector Variant mapping bias direction
#' @param p_s numeric vector Phase switch probabilities
#' @param theta numeric Haplotype imbalance
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @return numeric vector Allele likelihoods
#' @keywords internal          
calc_allele_lik_s2 = function (pAD, DP, R, p_s, theta, gamma = 20, r = 0.015) {
    hmm = get_allele_hmm_s2(pAD, DP, R, p_s, theta, gamma = gamma, r = r)
    LL = likelihood_allele(hmm)
    return(LL)
}

#' Calculate allele likelihoods for 3-state allele HMM
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param p_s numeric vector Phase switch probabilities
#' @param theta numeric Haplotype imbalance
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @keywords internal       
calc_allele_lik_s3 = function(pAD, DP, R, p_s, t, theta, gamma = 20, r = 0.015) {
    hmm = get_allele_hmm_s3(pAD, DP, R, p_s, t, theta, gamma = gamma, r = r)
    LL = likelihood_allele(hmm)
    return(LL)
}

############ Joint HMM for CNV detection ############

#' Run 7-state joint HMM on a pseudobulk profile
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param R numeric vector Variant mapping bias direction
#' @param p_s numeric vector Phase switch probabilities
#' @param theta_min numeric Minimum haplotype imbalance threshold
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @param Y_obs numeric vector Observed gene counts
#' @param lambda_ref numeric vector Reference expression rates
#' @param d_total integer Total library size for expression counts
#' @param phi_del numeric Expected fold change for deletion
#' @param phi_amp numeric Expected fold change for amplification
#' @param mu numeric Global expression bias
#' @param sig numeric Global expression variance
#' @param t numeric Transition probability between copy number states
#' @param r numeric Variant mapping bias
#' @param debug logical Whether to print debug messages
#' @return character vector Decoded states
#' @keywords internal
run_joint_hmm_s7 = function(
    pAD, DP, R, p_s, Y_obs, lambda_ref, d_total, 
    theta_min = 0.08, phi_del = 2^(-0.25), phi_amp = 2^(0.25), 
    t = 1e-5, mu = 0, sig = 1, gamma = 20, r = 0.015,
    debug = FALSE
) {

    # states
    states = c("1" = "neu", 
        "2" = "del_up", "3" = "del_down",
        "4" = "loh_up", "5" = "loh_down", 
        "6" = "amp_up", "7" = "amp_down")

    states_cn = str_remove(states, '_up|_down')
    states_phase = str_extract(states, 'up|down')

    # relative abundance of states
    w = c('neu' = 1, 'del' = 1, 'amp' = 1, 'loh' = 1)
        
    prior = sapply(1:length(states), function(to){
            get_trans_probs_s7(
                t = t, p_s = 0, w,
                cn_from = 'neu', phase_from = NA,
                cn_to = states_cn[to], phase_to = states_phase[to])
        })

    states_index = 1:length(states)
    
    # transition matrices
    As = calc_trans_mat_s7(t, p_s, w, states_cn, states_phase)

    # display(As[,,10])

    theta_u = 0.5 + theta_min
    theta_d = 0.5 - theta_min

    # parameters for each state
    theta_states = c(0.5, rep(c(theta_u, theta_d), 3))
    theta_states = outer(R*r, theta_states, FUN = "+")
    # theta_states = outer((1+2*R*r), theta_states, FUN = "*")
    alpha_states = theta_states * gamma
    beta_states = (1 - theta_states) * gamma
    
    phi_states = c(1, rep(phi_del, 2), rep(1, 2), rep(phi_amp, 2))

    N = length(Y_obs)

    if (length(mu) == 1) {
        mu = rep(mu, N)
        sig = rep(sig, N)
    }

    if (length(d_total) == 1) {
        d_total = rep(d_total, N)
    }

    hmm = list(
        x = pAD, 
        d = DP,
        y = Y_obs,
        l = d_total,
        lambda = lambda_ref,
        mu = mu,
        sig = sig,
        logPi = log(As), 
        phi = phi_states,
        delta = prior, 
        alpha = alpha_states,
        beta = beta_states,
        states = states,
        p_s = p_s
    )

    MPC = states[viterbi_joint(hmm)]
        
    return(MPC)
}


#' Calculate the transition matrix for 7-state joint HMM
#' @param t numeric; CNV state transition probability
#' @param p_s numeric; phase switch probability
#' @param w numeric; relative abundance of states
#' @param states_cn character; CNV states
#' @param states_phase character; haplotype phase states
#' @return array; transition matrix
#' @keywords internal
calc_trans_mat_s7 = function(t, p_s, w, states_cn, states_phase) {

    sapply(1:length(states_cn), function(from) {
        sapply(1:length(states_cn), function(to) {
            get_trans_probs_s7(t, p_s, w, states_cn[from], states_phase[from], states_cn[to], states_phase[to])
        }) %>% t
    }) %>% t %>%
    array(dim = c(length(states_cn), length(states_cn), length(p_s)))

}

#' Helper function to calculate transition porbabilities for 7-state joint HMM
#' cn/phase are sclars, only p_s is vectorized
#' @param t numeric; CNV state transition probability
#' @param p_s numeric; phase switch probability
#' @param w numeric; relative abundance of states
#' @param cn_from character; originating CNV state
#' @param phase_from character; originating haplotype phase state
#' @param cn_to character; destination CNV state
#' @param phase_to character; destination haplotype phase state
#' @return numeric; transition probability
#' @keywords internal
get_trans_probs_s7 = function(t, p_s, w, cn_from, phase_from, cn_to, phase_to) {

    if (cn_from == 'neu' & cn_to == 'neu') {
        p_s = rep(0.5, length(p_s))
    }

    if (cn_from == cn_to) {
        if (is.na(phase_from) & is.na(phase_to)) {
            p = 1-t
            p = rep(p, length(p_s))
        } else if (phase_from == phase_to) {
            p = (1-t) * (1-p_s)
        } else {
            p = (1-t) * p_s
        }
    } else {
        if (cn_from != 'neu') {
            w[c('del', 'amp', 'loh')] = 0
        }
        p = t * w[[cn_to]]/sum(w[names(w)!=cn_from])
        if (!is.na(phase_to)) {
            p = p/2
        }
        p = rep(p, length(p_s))
    }
    
    return(p)
}

#' Viterbi algorithm for joint HMM, can only handle one set of observations
#' @param hmm HMM object; expect variables x (allele depth), d (total depth), y (expression count), l (cell total library size),
#' lambda (reference expression rate), mu (global expression mean), sig (global expression standard deviation), 
#' logPi (log transition prob matrix), phi (expression fold change for each state), delta (prior for each state), 
#' alpha (alpha for each state), beta (beta for each state), states (states), p_s (phase switch probs)
#' @return character vector Decoded states
#' @keywords internal
viterbi_joint = function(hmm) {

    N <- length(hmm$x)
    M <- nrow(hmm$logPi[,,1])
    nu <- matrix(NA, nrow = N, ncol = M)
    z <- rep(NA, N)
    
    logprob = sapply(1:M, function(m) {

        l_x = dbbinom(x = hmm$x, size = hmm$d, alpha = hmm$alpha[,m], beta = hmm$beta[,m], log = TRUE)

        l_x[is.na(l_x)] = 0

        if (!is.null(hmm$y)) {
            l_y = rep(0,N)
            valid = !is.na(hmm$y)

            l_y[valid] = dpoilog(
                    x = hmm$y[valid],
                    sig = hmm$sig[valid],
                    mu = hmm$mu[valid] + log(hmm$phi[m] * hmm$l[valid] * hmm$lambda[valid]),
                    log = TRUE
                )
        } else {
            l_y = 0
        }

        return(l_x + l_y)

    })

    z = viterbi_compute(log(hmm$delta), logprob, hmm$logPi, N, M, nu, z)
        
    return(z)
}

#' Run 15-state joint HMM on a pseudobulk profile
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param p_s numeric vector Phase switch probabilities
#' @param theta_min numeric Minimum haplotype imbalance threshold
#' @param theta_neu numeric Haplotype imbalance threshold for neutral state
#' @param bal_cnv logical Whether to include balanced CNV states
#' @param r numeric Variant mapping bias
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @param Y_obs numeric vector Observed gene counts
#' @param lambda_ref numeric vector Reference expression rates
#' @param d_total integer Total library size for expression counts
#' @param phi_del numeric Expected fold change for deletion
#' @param phi_amp numeric Expected fold change for amplification
#' @param phi_bamp numeric Expected fold change for balanced amplification
#' @param phi_bdel numeric Expected fold change for balanced deletion
#' @param mu numeric Global expression bias
#' @param sig numeric Global expression variance
#' @param t numeric Transition probability between copy number states
#' @param prior numeric vector Prior probabilities for each state
#' @param exp_only logical Whether to only use expression data
#' @param allele_only logical Whether to only use allele data
#' @param classify_allele logical Whether to classify allele states
#' @param debug logical Whether to print debug messages
#' @param ... Additional parameters
#' @return character vector Decoded states
#' @examples
#' with(bulk_example, {
#'     run_joint_hmm_s15(pAD = pAD, DP = DP, p_s = p_s, Y_obs = Y_obs, lambda_ref = lambda_ref, 
#'     d_total = na.omit(unique(d_obs)), mu = mu, sig = sig, t = 1e-5, gamma = 30, theta_min = 0.08)
#' })
#' @export
run_joint_hmm_s15 = function(
    pAD, DP, p_s, Y_obs = 0, lambda_ref = 0, d_total = 0, theta_min = 0.08, theta_neu = 0,
    bal_cnv = TRUE, phi_del = 2^(-0.25), phi_amp = 2^(0.25), phi_bamp = phi_amp, phi_bdel = phi_del, 
    mu = 0, sig = 1, r = 0.015,
    t = 1e-5, gamma = 18,
    prior = NULL, exp_only = FALSE, allele_only = FALSE,
    classify_allele = FALSE, debug = FALSE, ...
) {

    # states
    states = c(
        "1" = "neu", "2" = "del_1_up", "3" = "del_1_down", "4" = "del_2_up", "5" = "del_2_down",
        "6" = "loh_1_up", "7" = "loh_1_down", "8" = "loh_2_up", "9" = "loh_2_down", 
        "10" = "amp_1_up", "11" = "amp_1_down", "12" = "amp_2_up", "13" = "amp_2_down", 
        "14" = "bamp", "15" = "bdel"
    )

    states_cn = str_remove(states, '_up|_down')
    states_phase = str_extract(states, 'up|down')

    # relative abundance of states
    w = c('neu' = 1, 'del_1' = 1, 'del_2' = 1e-10, 'loh_1' = 1, 'loh_2' = 1e-10, 'amp_1' = 1, 'amp_2' = 1e-10, 'bamp' = 1e-4, 'bdel' = 1e-10)
        
    # intitial probabilities
    if (is.null(prior)) {
        # encourage CNV from telomeres
        prior = sapply(1:length(states), function(to){
                get_trans_probs_s15(
                    t = min(t * 100, 1), p_s = 0, w,
                    cn_from = 'neu', phase_from = NA,
                    cn_to = states_cn[to], phase_to = states_phase[to])
            })
    }

    # to do: renormalize the probabilities after deleting states
    states_index = 1:length(states)

    if (!bal_cnv) {
        states_index = 1:13
    }
        
    if (exp_only) {
        pAD = rep(NA, length(pAD))
        p_s = rep(0, length(p_s))
    }
    
    if (allele_only) {
        states_index = c(1, 6:9)

        Y_obs = rep(NA, length(Y_obs))
    }

    if (classify_allele) {
        states_index = c(6,7)
    }
    
    # transition matrices
    As = calc_trans_mat_s15(t, p_s, w, states_cn, states_phase)

    theta_u_1 = 0.5 + theta_min
    theta_d_1 = 0.5 - theta_min

    theta_u_2 = 0.9
    theta_d_2 = 0.1

    theta_u_neu = 0.5 + theta_neu
    theta_d_neu = 0.5 - theta_neu

    # parameters for each state
    alpha_states = gamma * c(theta_u_neu, rep(c(theta_u_1, theta_d_1, theta_u_2, theta_d_2), 3), theta_u_neu, theta_u_neu)
    beta_states = gamma * c(theta_d_neu, rep(c(theta_d_1, theta_u_1, theta_d_2, theta_u_2), 3), theta_d_neu, theta_d_neu)
    phi_states = c(1, rep(phi_del, 2), rep(0.5, 2), rep(1, 4), rep(phi_amp, 2), rep(2.5, 2), phi_bamp, phi_bdel)
    
    # subset for relevant states
    prior = prior[states_index]
    As = As[states_index, states_index,]
    alpha_states = alpha_states[states_index]
    beta_states = beta_states[states_index]
    phi_states = phi_states[states_index]
    states = states[states_index] %>% setNames(1:length(.))
                
    N = length(Y_obs)
    M = length(states)

    if (length(mu) == 1) {
        mu = rep(mu, N)
        sig = rep(sig, N)
    }

    if (length(d_total) == 1) {
        d_total = rep(d_total, N)
    }

    hmm = list(
        x = pAD, 
        d = DP,
        y = Y_obs,
        l = d_total,
        lambda = lambda_ref,
        mu = mu,
        sig = sig,
        logPi = log(As), 
        phi = phi_states,
        delta = prior, 
        alpha = matrix(rep(alpha_states, N), ncol = M, byrow = TRUE),
        beta = matrix(rep(beta_states, N), ncol = M, byrow = TRUE),
        states = states,
        p_s = p_s
    )

    MPC = states[viterbi_joint(hmm)]
        
    return(MPC)
}


#' Calculate the transition matrix for 15-state joint HMM
#' @param t numeric; CNV state transition probability
#' @param p_s numeric; phase switch probability
#' @param w numeric; relative abundance of states
#' @param states_cn character; CNV states
#' @param states_phase character; haplotype phase states
#' @return array; transition matrix
#' @keywords internal
calc_trans_mat_s15 = function(t, p_s, w, states_cn, states_phase) {

    sapply(1:length(states_cn), function(from) {
        sapply(1:length(states_cn), function(to) {
            get_trans_probs_s15(t, p_s, w, states_cn[from], states_phase[from], states_cn[to], states_phase[to])
        }) %>% t
    }) %>% t %>%
    array(dim = c(length(states_cn), length(states_cn), length(p_s)))

}




#' Helper function to calculate transition porbabilities for 15-state joint HMM
#' cn/phase are sclars, only p_s is vectorized
#' @param t numeric; CNV state transition probability
#' @param p_s numeric; phase switch probability
#' @param w numeric; relative abundance of states
#' @param cn_from character; originating CNV state
#' @param phase_from character; originating haplotype phase state
#' @param cn_to character; destination CNV state
#' @param phase_to character; destination haplotype phase state
#' @return numeric; transition probability
#' @keywords internal
get_trans_probs_s15 = function(t, p_s, w, cn_from, phase_from, cn_to, phase_to) {

    if (cn_from == 'neu' & cn_to == 'neu') {
        p_s = rep(0.5, length(p_s))
    }

    if (cn_from == cn_to) {
        if (is.na(phase_from) & is.na(phase_to)) {
            p = 1-t
            p = rep(p, length(p_s))
        } else if (phase_from == phase_to) {
            p = (1-t) * (1-p_s)
        } else {
            p = (1-t) * p_s
        }
    } else {
        p = t * w[[cn_to]]/sum(w[names(w)!=cn_from])
        if (!is.na(phase_to)) {
            p = p/2
        }
        p = rep(p, length(p_s))
    }
    
    return(p)
}

#' Generalized viterbi algorithm for joint HMM, can handle multiple sets of observations
#' @param hmm HMM object; expect variables x (allele depth), d (total depth), y (expression count), l (cell total library size),
#' lambda (reference expression rate), mu (global expression mean), sig (global expression standard deviation), 
#' logPi (log transition prob matrix), phi (expression fold change for each state), delta (prior for each state), 
#' alpha (alpha for each state), beta (beta for each state), states (states), p_s (phase switch probs)
#' @return character vector; decoded states
#' @keywords internal
viterbi_joint_mat <- function(hmm) {

    N <- nrow(hmm$x)
    M <- nrow(hmm$logPi[,,1])
    K <- ncol(hmm$x)
    nu <- matrix(NA, nrow = N, ncol = M)
    z <- rep(NA, N)
    
    logprob = sapply(1:M, function(m) {

        sapply(
            1:K,
            function(k) {

                l_x = dbbinom(x = hmm$x[,k], size = hmm$d[,k], alpha = hmm$alpha[m,k], beta = hmm$beta[m,k], log = TRUE)

                l_x[is.na(l_x)] = 0

                if (!is.null(hmm$y)) {
                    l_y = rep(0,N)
                    valid = !is.na(hmm$y[,k])

                    l_y[valid] = dpoilog(
                            x = hmm$y[valid,k],
                            sig = hmm$sig[valid,k],
                            mu = hmm$mu[valid,k] + log(hmm$phi[m,k] * hmm$l[valid,k] * hmm$lambda[valid,k]),
                            log = TRUE
                        )
                } else {
                    l_y = 0
                }

                return(l_x + l_y)

            }
        ) %>% rowSums()

    })

    z = viterbi_compute(log(hmm$delta), logprob, hmm$logPi, N, M, nu, z)
        
    return(z)
}