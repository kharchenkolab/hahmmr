########################### Visualizations ############################
# default color palette
## pal = RColorBrewer::brewer.pal(n = 8, 'Set1')
pal = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
getPalette = colorRampPalette(pal)

#' @keywords internal
cnv_colors = c("neu" = "gray",
        "neu_up" = "darkgray", "neu_down" = "gray",
        "del_up" = "royalblue", "del_down" = "darkblue",
        "loh_up" = "darkgreen", "loh_down" = "olivedrab4",
        "amp_up" = "red", "amp_down" = "tomato3",
        "del_1_up" = "royalblue", "del_1_down" = "darkblue",
        "loh_1_up" = "darkgreen", "loh_1_down" = "olivedrab4",
        "amp_1_up" = "red", "amp_1_down" = "tomato3",
        "del_2_up" = "royalblue", "del_2_down" = "darkblue",
        "loh_2_up" = "darkgreen", "loh_2_down" = "olivedrab4",
        "amp_2_up" = "red", "amp_2_down" = "tomato3",
        "del_up_1" = "royalblue", "del_down_1" = "darkblue",
        "loh_up_1" = "darkgreen", "loh_down_1" = "olivedrab4",
        "amp_up_1" = "red", "amp_down_1" = "tomato3",
        "del_up_2" = "royalblue", "del_down_2" = "darkblue",
        "loh_up_2" = "darkgreen", "loh_down_2" = "olivedrab4",
        "amp_up_2" = "red", "amp_down_2" = "tomato3",
        "bamp" = "salmon", "bdel" = "skyblue",
        "amp" = "tomato3", "loh" = "olivedrab4", "del" = "royalblue",
        "theta_up" = "darkgreen", "theta_down" = "olivedrab4",
        "theta_1_up" = "darkgreen", "theta_1_down" = "olivedrab4",
        "theta_2_up" = "darkgreen", "theta_2_down" = "olivedrab4",
        "theta_up_1" = "darkgreen", "theta_down_1" = "olivedrab4",
        "theta_up_2" = "darkgreen", "theta_down_2" = "olivedrab4",
        '0|1' = 'red', '1|0' = 'blue','major' = '#66C2A5', 'minor' = '#FC8D62')

#' @keywords internal
cnv_labels = names(cnv_colors) %>%
    stringr::str_remove_all('_') %>%
    stringr::str_to_upper() %>%
    stringr::str_replace('UP', '(major)') %>%
    stringr::str_replace('DOWN', '(minor)') %>%
    stringr::str_replace('LOH', 'CNLoH') %>%
    setNames(names(cnv_colors))


#' Plot a pseudobulk HMM profile
#'
#' @param bulk dataframe Pseudobulk profile
#' @param use_pos logical Use marker position instead of index as x coordinate
#' @param allele_only logical Only plot alleles
#' @param min_LLR numeric LLR threshold for event filtering
#' @param min_depth numeric Minimum coverage depth for a SNP to be plotted
#' @param exp_limit numeric Expression logFC axis limit
#' @param phi_mle logical Whether to plot estimates of segmental expression fold change
#' @param theta_roll logical Whether to plot rolling estimates of allele imbalance
#' @param dot_size numeric Size of marker dots
#' @param dot_alpha numeric Transparency of the marker dots
#' @param legend logical Whether to show legend
#' @param exclude_gap logical Whether to mark gap regions and centromeres
#' @param genome character Genome build, either 'hg38' or 'hg19'
#' @param text_size numeric Size of text in the plot
#' @param raster logical Whether to raster images
#' @return ggplot Plot of pseudobulk HMM profile
#' @examples
#' p = plot_psbulk(bulk_example)
#' @export
plot_psbulk = function(
        bulk, use_pos = TRUE, allele_only = FALSE, min_LLR = 5, min_depth = 8, exp_limit = 2,
        phi_mle = TRUE, theta_roll = FALSE, dot_size = 0.8, dot_alpha = 0.5, legend = TRUE,
        exclude_gap = TRUE, genome = 'hg38', text_size = 10, raster = FALSE
    ) {

    if (!all(c('state_post', 'cnv_state_post') %in% colnames(bulk))) {
        bulk = bulk %>%
            mutate(
                state_post = state,
                cnv_state_post = cnv_state
            )
    }

    # filter events by LLR
    if (min_LLR != 0) {
        bulk = bulk %>% mutate(
            LLR = ifelse(is.na(LLR), 0, LLR),
            cnv_state_post = ifelse(LLR < min_LLR, 'neu', cnv_state_post),
            state_post = ifelse(LLR < min_LLR, 'neu', state_post)
        )
    }

    # mark clonal LOH
    if ('loh' %in% colnames(bulk)) {
        bulk = bulk %>% mutate(state_post = ifelse(loh, 'del', state_post))
    }

    if (use_pos) {
        marker = 'POS'
        marker_label = 'Genomic position'
    } else {
        marker = 'snp_index'
        marker_label = 'SNP index'
    }

    # fix retest states
    bulk = bulk %>%
        mutate(
            theta_level = ifelse(str_detect(state_post, '_2'), 2, 1),
            state_post = ifelse(
                cnv_state_post %in% c('amp', 'loh', 'del'),
                ifelse(p_up > 0.5, paste0(cnv_state_post, '_', theta_level, '_', 'up'), paste0(cnv_state_post, '_', theta_level, '_', 'down')),
                state_post
        ))

    # correct for baseline bias
    if (!allele_only) {
        bulk = bulk %>% mutate(logFC = logFC - mu)
    } else {
        bulk = bulk %>% mutate(logFC = NA)
    }

    D = bulk %>%
        mutate(logFC = ifelse(logFC > exp_limit | logFC < -exp_limit, NA, logFC)) %>%
        mutate(pBAF = ifelse(DP >= min_depth, pBAF, NA)) %>%
        mutate(pHF = pBAF) %>%
        as.data.table %>%
        data.table::melt(measure.vars = c('logFC', 'pHF'))

    if (allele_only) {
        D = D %>% filter(variable == 'pHF')
    }

    p = ggplot(
            D,
            aes(x = get(marker), y = value, color = state_post),
            na.rm=TRUE
        )

    if (use_pos & exclude_gap) {

        if (genome == 'hg38') {
            gaps = gaps_hg38 %>% filter(end - start > 1e+06)
            acen = acen_hg38
        } else if (genome == 'hg19') {
            gaps = gaps_hg19 %>% filter(end - start > 1e+06)
            acen = acen_hg19
        } else if (genome == 'mm10') {
            gaps = data.frame(CHROM = 1, start = 1, end = 1)
            acen = data.frame()
        } else {
            stop("Genome version must hg38, hg19 or mm10")
        }

        segs_exclude = rbind(gaps, acen) %>%
            mutate(CHROM = factor(as.integer(CHROM))) %>%
            rename(seg_start = start, seg_end = end) %>%
            filter(CHROM %in% bulk$CHROM)

        if (nrow(segs_exclude) > 0) {
            p = p + geom_rect(inherit.aes = FALSE, data = segs_exclude,
                aes(xmin = seg_start, xmax = seg_end, ymin = -Inf, ymax = Inf),
                fill = "gray95")
        }
    }

    legend_breaks = c("neu", "loh_up", "loh_down", "del_up", "del_down", "amp_up", "amp_down", "bamp", "bdel")

    p = p + geom_point(
            aes(shape = str_detect(state_post, '_2'), alpha = str_detect(state_post, '_2')),
            size = dot_size,
            na.rm = TRUE
        ) +
        geom_hline(
            data = data.frame(y = c(0,1), variable = 'pHF'),
            aes(yintercept = y),
            size = 0, alpha = 0
        ) +
        suppressWarnings(scale_alpha_discrete(range = c(dot_alpha, 1))) +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15)) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, 'mm'),
            panel.spacing.y = unit(1, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.title = element_text(size = text_size),
            strip.text = element_text(size = text_size),
            axis.title = element_text(size = text_size),
            legend.text = element_text(size = text_size),
            plot.margin = margin(t = 1, r = 0, b = 1, l = 0, 'cm')
        ) +
        facet_grid(variable ~ CHROM, scales = 'free', space = 'free_x') +
        # scale_x_continuous(expand = expansion(add = 5)) +
        scale_color_manual(
            values = cnv_colors,
            limits = names(cnv_colors),
            breaks = legend_breaks,
            labels = cnv_labels[legend_breaks],
            na.translate = FALSE
        ) +
        guides(
            color = guide_legend(title = "CNV state", override.aes = aes(size = 3), ncol = 1),
            fill = 'none', alpha = 'none', shape = 'none'
        ) +
        xlab(marker) +
        ylab('')

    if (!allele_only) {
        p = p + geom_hline(
                data = data.frame(y = c(-exp_limit, exp_limit), variable = 'logFC'),
                aes(yintercept = y),
                size = 0, alpha = 0)
    }

    if (!legend) {
        p = p + guides(color = 'none', fill = 'none', alpha = 'none', shape = 'none')
    }

    if (phi_mle & (!allele_only)) {
        segs = bulk %>%
            distinct(CHROM, seg, seg_start, seg_start_index, seg_end, seg_end_index, phi_mle) %>%
            mutate(variable = 'logFC') %>%
            filter(log2(phi_mle) < exp_limit)

        if (use_pos) {
            start = 'seg_start'
            end = 'seg_end'
        } else {
            start = 'seg_start_index'
            end = 'seg_end_index'
        }

        p = p + geom_segment(
            inherit.aes = FALSE,
            data = segs,
            aes(x = get(start), xend = get(end), y = log2(phi_mle), yend = log2(phi_mle)),
            color = 'darkred',
            size = 0.5
        ) +
        geom_hline(data = data.frame(variable = 'logFC'), aes(yintercept = 0), color = 'gray30', linetype = 'dashed')
    } else if (!allele_only) {
        p = p + geom_line(
            inherit.aes = FALSE,
            data = bulk %>% mutate(variable = 'logFC') %>% filter(log2(phi_mle_roll) < exp_limit),
            aes(x = get(marker), y = log2(phi_mle_roll), group = '1'),
            color = 'darkred',
            size = 0.35
        ) +
        geom_hline(data = data.frame(variable = 'logFC'), aes(yintercept = 0), color = 'gray30', linetype = 'dashed')
    }

    if (theta_roll) {
        p = p +
            geom_line(
                inherit.aes = FALSE,
                data = D %>% mutate(variable = 'pHF'),
                aes(x = snp_index, y = 0.5 - theta_hat_roll, color = paste0(cnv_state_post, '_down')),
                # color = 'black',
                size = 0.35
            ) +
            geom_line(
                inherit.aes = FALSE,
                data = D %>% mutate(variable = 'pHF'),
                aes(x = snp_index, y = 0.5 + theta_hat_roll, color = paste0(cnv_state_post, '_up')),
                # color = 'gray',
                size = 0.35
            )
    }

    p = p + xlab(marker_label)

    if (raster) {
        p = ggrastr::rasterize(p, layers = 'Point', dpi = 300)
    }

    return(p)
}

#' Plot a group of pseudobulk HMM profiles
#'
#' @param bulks dataframe Pseudobulk profiles annotated with "sample" column
#' @param ncol integer Number of columns
#' @param title logical Whether to add titles to individual plots
#' @param title_size numeric Size of titles
#' @param ... additional parameters passed to plot_psbulk()
#' @return a ggplot object
#' @examples
#' p = plot_bulks(bulk_example)
#' @export
plot_bulks = function(
    bulks, ..., ncol = 1, title = TRUE, title_size = 8
    ) {

    if (!'sample' %in% colnames(bulks)) {
        bulks$sample = 1
    }

    plot_list = bulks %>%
        split(.$sample) %>%
        lapply(
            function(bulk) {

                sample = unique(bulk$sample)
                n_cells = unique(bulk$n_cells)

                p = plot_psbulk(
                        bulk, ...
                    ) +
                    theme(
                        title = element_text(size = title_size),
                        axis.text.x = element_blank(),
                        axis.title = element_blank(),
                        plot.margin = margin(t = 0, r = 0, b = 0.25, l = 0, 'cm')
                    )

                if (title) {
                    if (is.null(n_cells)) {
                        title_text = sample
                    } else {
                        title_text = glue('{sample} (n={n_cells})')
                    }
                    p = p + ggtitle(title_text)
                }

                return(p)
            }
        )

    panel = wrap_plots(plot_list, ncol = ncol, guides = 'collect')

    return(panel)
}
