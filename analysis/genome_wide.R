# qrsh -l h_rt=48:00:00,highp,h_data=32G
# module load gcc/9.3.0
# module load R/4.0.2
library('dplyr')
library('tidyr')
library('nimble')
library('ggplot2')
library('cowplot')
library('ggrepel')
library('scales')
theme_set(theme_cowplot())
library('waterbear')

raw_counts = read.table('../data/WillScreen.count.txt', header = TRUE, sep = '\t')

ordering = c('1', '2', '3', '4')

c_name = grep('^(A|B)[0-9]', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)

sample_mapping = mutate(sample_mapping, sample = sub('[0-9]*$', '', c_name))
sample_mapping = mutate(sample_mapping, bin = sub('(A|B)', '', c_name))

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

devtools::install('~/Dropbox/postdoc/waterbear')

counts_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

# debugonce(wb_make_object)
wo = wb_make_object(counts_array, gene_mapping, bin_size_prior = c(0.15, 0.35, 0.35, 0.15),
  control_guide_regex = NA)

# debugonce(wb_em_start)
wo = wb_em_start(wo)

source('../../waterbear/R/nimble.R')

# ~1hr
system.time({
n_model = nimbleModel(sample_specific_dispersion_model_no_controls,
  data = wo$data, constants = wo$const, inits = wo$init)
})

nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)

system.time({
n_configuration = configureMCMC(n_model)
sampler_control = list(order = wo$const$order)
# n_configuration$removeSamplers('gene_inclusion')
# n_configuration$addSampler(target = 'gene_inclusion', type = 'rank_RW',
#   control = sampler_control)
})

n_configuration$addMonitors(
  c('gene_inclusion',
    'total_shift',
    'guide_shift',
    'gene_shift',
    'dispersion',
    'sample_dispersion',
    'psi'
    ))

system.time({
n_mcmc_build = buildMCMC(n_configuration)
})
system.time({
C_n_model = compileNimble(n_model)
})
system.time({
C_n_mcmc = compileNimble(n_mcmc_build, project = n_model, showCompilerOutput = TRUE)
})

n_samples = 4000
n_burnin = 2000
n_chains = 4
cur_seed = 326

system.time({
samples = runMCMC(
  C_n_mcmc, niter = n_samples, nburnin = n_burnin,
  nchains = n_chains,
  setSeed = cur_seed,
  summary = TRUE)
})

saveRDS(samples, 'samples.rds')

# run parallel

cl = makeCluster(4)

run_mcmc = function(seed, wo, n_samples, n_burnin) {
  library('nimble')
  # source('/u/project/hjp/hjp/waterbear/nimble.R')
  ddirchmulti <- nimbleFunction(
    run = function(x = double(1), alpha = double(1), size = double(0),
      log = integer(0, default = 0)) {
      returnType(double(0))
      logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) -
        sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) +
          size)
      if(log) return(logProb)
      else return(exp(logProb))
    })

  rdirchmulti <- nimbleFunction(
    run = function(n = integer(0), alpha = double(1), size = double(0)) {
      returnType(double(1))
      if(n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
      p <- rdirch(1, alpha)
      return(rmulti(1, size = size, prob = p))
    })

  dirichlet_to_normal_bins = nimbleFunction(
    run = function(alpha = double(1)) {
      returnType(double(1))
      a0 = sum(alpha)
      p = alpha / a0
      qp = numeric(length(p))
      qp[1] = p[1]
      for (i in 2:length(p)) {
        qp[i] = qp[i - 1] + p[i]
      }
      cutoff = numeric(length(qp) - 1)
      for (i in 1:(length(qp) - 1)) {
        cutoff[i] = probit(qp[i])
      }
      return(cutoff)
    })

  cutoff_to_probability = nimbleFunction(
    run = function(cutoff = double(1), offset = double(0)) {
      returnType(double(1))
      p = numeric(length(cutoff) + 1)
      p[1] = phi(cutoff[1] - offset)
      total = p[1]
      for (i in 2:length(cutoff)) {
        p[i] = phi(cutoff[i] - offset) - phi(cutoff[i - 1] - offset)
        total = total + p[i]
      }
      p[length(cutoff) + 1] = 1 - total
      return(p)
    })

  sample_specific_dispersion_model_no_controls = nimbleCode({
    dispersion ~ dexp(1.0 / dispersion_prior_mean)
    for (n in 1:N) {
      bin_alpha[n, 1:N_bins] ~ ddirch(bin_alpha_prior[n, 1:N_bins])
      cutoffs[n, 1:N_cutoffs] <- dirichlet_to_normal_bins(bin_alpha[n, 1:N_bins])
      sample_dispersion[n] ~ dexp(1.0 / dispersion)
    }
    psi ~ dbeta(10, 10)
    sigma_gene ~ dgamma(1, 0.10)
    sigma_guide ~ dgamma(1, 0.10)
    for (gene in 1:N_genes) {
      gene_inclusion[gene] ~ dbern(psi)
      gene_shift[gene] ~ dnorm(0, sd = sigma_gene)
    }
    for (g in 1:N_guides) {
      guide_shift[g] ~ dnorm(gene_shift[guide_to_gene[g]], sd = sigma_guide)
      total_shift[g] <- gene_inclusion[guide_to_gene[g]] * guide_shift[g]
      for (n in 1:N) {
        q[n, g, 1:N_bins] <- cutoff_to_probability(cutoffs[n, 1:N_cutoffs], total_shift[g])
        x[n, guide_data_index[g], 1:N_bins] ~ ddirchmulti(
          sample_dispersion[n] * q[n, g, 1:N_bins],
          x_total[n, guide_data_index[g]])
      }
    }
  })

  n_model = nimbleModel(sample_specific_dispersion_model_no_controls,
    data = wo$data, constants = wo$const, inits = wo$init)

  # nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
  # nimbleOptions(MCMCsaveHistory = TRUE)

  n_configuration = configureMCMC(n_model)
  sampler_control = list(order = wo$const$order)
  # n_configuration$removeSamplers('gene_inclusion')
  # n_configuration$addSampler(target = 'gene_inclusion', type = 'rank_RW',
  #   control = sampler_control)

  n_configuration$addMonitors(
    c('gene_inclusion',
      'total_shift',
      'guide_shift',
      'gene_shift',
      'dispersion',
      'sample_dispersion',
      'psi'
      ))

  n_mcmc_build = buildMCMC(n_configuration)
  C_n_model = compileNimble(n_model)
  C_n_mcmc = compileNimble(n_mcmc_build, project = n_model, showCompilerOutput = TRUE)

  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    setSeed = seed,
    nchains = 1)

  return(samples)
}

# chain_output <- parLapply(cl = this_cluster, X = 1:4,
#                           fun = run_MCMC_allcode,
#                           data = myData)
all_chains = parLapply(cl = cl, X = 1:4, wo = wo, n_samples = 10000, n_burnin = 5000,
  fun = run_mcmc)

tmp = run_mcmc(1, wo = wo, n_samples = 10000, n_burnin = 5000)

all_chains = parLapply(cl = cl, X = 1:4, wo = wo, n_samples = 10000, n_burnin = 5000,
  fun = run_mcmc)

# end parallel execution

########################################################################
# analysis of of samples (local)
########################################################################

samples = readRDS('./samples.rds')

debugonce(lfsr)

wb_lfsr = lfsr(samples)

wb_lfsr = mutate(wb_lfsr,
  mapping = as.integer(sub('\\]', '', sub('gene_shift\\[', '', mapping_name))))



# wb_lfsr = unlist(wb_lfsr)
# data.frame(mapping = 1:length())

wb_samples = list(samples = samples)

debugonce(waterbear:::wb_recode)
gi = waterbear:::wb_recode(samples, wo)

sample_summary = samples$summary$all.chains
rn = rownames(sample_summary)
sample_summary[grepl('bin_alpha', rn), ]
sample_summary[grepl('dispersion', rn), ]
sample_summary[grepl('psi', rn), ]

gi %>% arrange(desc(Mean))

tmp = dplyr::select(gi, gene, inclusion = Mean)
# gs = waterbear:::wb_recode(samples, wo, 'gene_shift')
gs = wb_recode(samples, wo, 'gene_shift')
# gs = dplyr::select(gs, gene, gene_effect = Mean)

gs = left_join(gs, select(wb_lfsr, -mapping_name), by = 'mapping')

interesting_genes = c('B2m',
  'AU040320',
  'H2-Q7',
  # 'Nlrc5',
  'Gpaa1',
  'Gpr108',
  'H2-Q6')
gs = left_join(gs, data.frame(gene = interesting_genes, label = interesting_genes),
  by = 'gene')

p = ggplot(gs, aes(Mean, -log10(lfsr)))
p = p + geom_point(aes(color = factor(lfsr < 0.10), alpha = 1 - sqrt(lfsr)))
p = p + geom_text_repel(aes(label = label),
  min.segment.length = 0,
  arrow = arrow(length = unit(0.015, "npc")),
  nudge_x = .1,
  point.padding = 1
  )
p = p + scale_color_manual(values = c('gray', '#D3837A'))
p = p + xlab('Gene effect size')
# p = p + ylab('Gene effect size')
p = p + theme(legend.position="none")
save_plot('no_nlrc5_lfsr.png', p, base_asp = 4 / 3)
# p

p = ggplot(gs, aes(Mean, -log10(lfsr)))
p = p + geom_point(aes(color = factor(lfsr < 0.10), alpha = 1 - sqrt(lfsr)))
# p = p + geom_text_repel(aes(label = label),
#   min.segment.length = 0,
#   arrow = arrow(length = unit(0.015, "npc")),
#   nudge_x = .1,
#   point.padding = 1
#   )
p = p + scale_color_manual(values = c('gray', '#D3837A'))
p = p + xlab('Gene effect size')
# p = p + ylab('Gene effect size')
p = p + theme(legend.position="none")
save_plot('lfsr_no_labels.png', p, base_asp = 4 / 3)
# p

interesting_genes2 = c(
  'Pign',
  'Pigs',
  'Pigv',
  'Pigyl',
  'Piga'
)

gs = left_join(gs, data.frame(gene = interesting_genes2, label2 = interesting_genes2),
  by = 'gene')

p = ggplot(gs, aes(Mean, -log10(lfsr)))
p = p + geom_point(aes(color = factor(lfsr < 0.10), alpha = 1 - sqrt(lfsr)))
p = p + geom_text_repel(aes(label = label2),
  min.segment.length = 0,
  arrow = arrow(length = unit(0.015, "npc")),
  nudge_x = .1,
  point.padding = 1
  )
p = p + scale_color_manual(values = c('gray', '#D3837A'))
p = p + xlab('Gene effect size')
# p = p + ylab('Gene effect size')
p = p + theme(legend.position="none")
save_plot('lfsr_pig.png', p, base_asp = 4 / 3)
# p

interesting_genes = c('B2m', 'AAVR', 'H2-Q7', 'Gpaa1', 'Gpr108', 'H2-Q6')

write.table(dplyr::select(gs, -c(label, mapping)), sep = '\t', 'lfsr.tsv',
  quote = FALSE, row.names = FALSE)

inner_join(gs, data.frame(gene = interesting_genes), by = 'gene')

filter(gs, gene == '')


wb_results = inner_join(gs, tmp, by = 'gene')
wb_results = dplyr::mutate(wb_results, fdr = 1 - inclusion)

mageck_results = read.table(
  '~/data/aav/2021-08-24_mouseAAV-CRISPR-screen/MaGeck_tests/A1B1vsA4B4.gene_summary.txt',
  header = TRUE, stringsAsFactors = FALSE)
mageck_results = dplyr::mutate(mageck_results, mageck_fdr = pmin(neg.fdr, pos.fdr),
  mageck_lfc = pos.lfc)

all_results = dplyr::left_join(wb_results,
  dplyr::select(mageck_results, id, mageck_fdr, mageck_lfc), by = c('gene' = 'id'))

p = ggplot(all_results, aes(gene_effect, -mageck_lfc, color = factor(fdr < 0.05)))
p = p + geom_point(alpha = 0.10)
p

# I can already tell you that we independently validated B2M / H2Q7 / GPR108 and AU040320

dplyr::filter(all_results, grepl('b2m', gene, ignore.case = TRUE))
dplyr::filter(all_results, grepl('h2-q7', gene, ignore.case = TRUE))
dplyr::filter(all_results, grepl('^h2', gene, ignore.case = TRUE))
dplyr::filter(all_results, grepl('gpr108', gene, ignore.case = TRUE))
dplyr::filter(all_results, grepl('au040320', gene, ignore.case = TRUE))



all_results = dplyr::mutate(all_results,
  significant = (function(w, m, w_alpha, m_alpha) {
    res = vector('character', length(m))
    m[is.na(m)] = 1
    w[is.na(w)] = 1
    for (i in 1:length(m)) {
      if (w[i] < w_alpha && m[i] > m_alpha) {
        res[i] = 'waterbear'
      } else if (m[i] < m_alpha && w[i] > w_alpha){
        res[i] = 'MAGeCK'
      } else if (w[i] < w_alpha && m[i] < m_alpha) {
        res[i] = 'both'
      } else {
        res[i] = 'neither'
      }
    }
    res
  })(fdr, mageck_fdr, 0.01, 0.05))

dir.create('../img')

p = ggplot(all_results, aes(gene_effect, -mageck_lfc, color = significant))
p = p + geom_point(alpha = 0.10)
save_plot('../img/mageck_waterbear_all.png', p, base_width = 10, base_height = 8)

p = ggplot(dplyr::filter(all_results, significant != 'neither'), aes(gene_effect, -mageck_lfc, color = significant))
p = p + geom_point(alpha = 0.40)
save_plot('../img/mageck_waterbear_sig.png', p, base_width = 10, base_height = 8)

# add-in annotations by william
william_summary = read.csv('../data/A1B1vsA4B4.gene_summary.csv', header = TRUE)
william_summary = dplyr::select(william_summary, gene = id, comment = Comment)

all_results = left_join(all_results, william_summary, by = 'gene')
all_results = dplyr::arrange(all_results, desc(comment))


write.csv(all_results, '../results/all_results.csv')

dir.create('../results')

all_results = dplyr::arrange(all_results, desc(abs(gene_effect)))
write(all_results$gene[all_results$fdr < 0.01],
  '../results/go_waterbear.txt')

write(all_results$gene[all_results$mageck_fdr < 0.05],
  '../results/go_mageck.txt')
write(all_results$gene, '../results/go_background.txt')

# getting the guide shift


# debugonce(waterbear:::wb_recode)
wb_recode_guide = function (wb_samples, wo, extract_regex)
{
    if (!is.null(wb_samples$summary$all.chains)) {
        s = wb_samples$summary$all.chains
    }
    gi = data.frame(s, mapping = rownames(s))
    gi = dplyr::filter(gi, grepl(extract_regex, mapping))
    gi = dplyr::mutate(gi, mapping = sub(paste0(extract_regex,
        "\\["), "", mapping))
    gi = dplyr::mutate(gi, mapping = sub("\\]", "", mapping))
    gi = dplyr::mutate(gi, mapping = as.integer(mapping))
    gm = dplyr::distinct(dplyr::select(wo$test_guide_names, gene,
        i))
    gi = inner_join(gi, gm, by = c('mapping' = "i"))
    gi
}
guide_effect = wb_recode_guide(samples, wo, 'guide_shift')

nlrc5 = dplyr::filter(guide_effect, gene == 'Nlrc5' | gene == 'B2m')

p = ggplot(nlrc5, aes(factor(mapping), Mean, color = factor(mapping)))
p = p + geom_point()
p = p + geom_errorbar(aes(ymin = X95.CI_low, ymax = X95.CI_upp))
p = p + facet_wrap(~gene)
p


guide_effect = dplyr::select(gs, gene, gene_effect = Mean)

# manhattan plot

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

keys = wb_results$gene[100:105]

tmp = select(txdb, keys = keys, columns = c('CDSCHROM', 'CDSSTART'), keytype = 'GENEID')

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

annot<-getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), mart=ensembl)


inner_join(annot, data.frame(mgi_symbol = keys), by = 'mgi_symbol')

tmp = inner_join(wb_results, annot, by = c('gene' = 'mgi_symbol'))
tmp$chromosome_name %>% unique

tmp = dplyr::mutate(tmp,
  chromosome_name = factor(chromosome_name, levels = c(as.character(1:19), 'X', 'Y')))


# adapted from: https://www.r-graph-gallery.com/101_Manhattan_plot.html

manhattan = group_by(tmp, chromosome_name) %>%
  summarize(chr_len = max(end_position)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  left_join(tmp, ., by=c("chromosome_name"="chromosome_name")) %>%
  arrange(chromosome_name, start_position) %>%
  mutate( BPcum=start_position+tot)


axisdf = manhattan %>% group_by(chromosome_name) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

manhattan_significant = dplyr::mutate(manhattan, high_pip = inclusion > 0.90)
manhattan_significant = dplyr::filter(manhattan_significant, high_pip)

p = ggplot(manhattan, aes(BPcum, inclusion ))
p = p + geom_point( aes(color=as.factor(chromosome_name)), alpha=0.8, size=1.3)
p = p + scale_color_manual(values = rep(c("grey", "skyblue"), length(unique(manhattan$chromosome_name))))
p = p + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
p = p + scale_x_continuous( label = axisdf$chromosome_name, breaks= axisdf$center )
# p = p + pseudo_log_trans(sigma = 1, base = exp(1))
p

p = ggplot(manhattan, aes(BPcum, gene_effect ))
p = p + geom_point( aes(color=as.factor(chromosome_name)), alpha=0.8, size=1.3)
p = p + scale_color_manual(values = rep(c("grey", "skyblue"),
    length(unique(manhattan$chromosome_name))))
p = p + geom_point(size = 1.5, color = 'orange', data = manhattan_significant)
p = p + geom_hline(yintercept = 0, linetype = 3)
p = p + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
p = p + scale_x_continuous( label = axisdf$chromosome_name, breaks= axisdf$center )
# p = p + pseudo_log_trans(sigma = 1, base = exp(1))
p = p + geom_text_repel(aes(label = highlight_gene))
p = p + theme(legend.position = 'none')
p = p + xlab('')
p = p + ylab('gene effect size')
save_plot(p, file = '../img/manhattan_gene_effect.png', base_width = 14, base_height = 10)
p

# p = ggplot(manhattan, aes(BPcum, gene_effect ))
p = ggplot(manhattan, aes(BPcum, inclusion ))
p = p + geom_point( aes(color=as.factor(chromosome_name)), alpha=0.8, size=1.3)
p = p + scale_color_manual(values = rep(c("grey", "skyblue"), length(unique(manhattan$chromosome_name))))
p = p + geom_point(size = 1.5, color = 'orange', data = manhattan_significant)
p = p + geom_hline(yintercept = 0.9, linetype = 3)
p = p + scale_y_continuous(trans=scales::tanh())
p = p + scale_x_continuous( label = axisdf$chromosome_name, breaks= axisdf$center )
# p = p + pseudo_log_trans(sigma = 1, base = exp(1))
p = p + geom_text_repel(aes(label = highlight_gene))
p = p + theme(legend.position = 'none')
p = p + ylab('posterior inclusion probability')
p = p + xlab('')
save_plot(p, file = '../img/manhattan.png', base_width = 14, base_height = 10)

genes_of_interest = base::scan(text = 'B2m
AU040320
H2-Q7
Il2rg
Gpaa1
Gpr108
H2-Q6
H2-D1', sep = '\n', what = character())

genes_of_interest = data.frame(gene = genes_of_interest, interest = rep(TRUE, length(genes_of_interest)))

wb_results = left_join(wb_results, genes_of_interest, by = 'gene')

wb_results = dplyr::mutate(wb_results,
  highlight_gene = ifelse(!is.na(interest), gene, interest))

cbb = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

p = ggplot(wb_results, aes(gene_effect, log10(inclusion)))
p = p + geom_point(aes(color = as.factor(inclusion > 0.90)))
p = p + geom_text_repel(aes(label = highlight_gene))
p = p + theme(legend.position = 'none')
p = p + scale_color_manual(values = cbb[c(2, 1)])
p = p + xlab('gene effect size')
p = p + ylab('log10(posterior inclusion probability)')
save_plot(p, file = '../img/volcano.png', base_height = 10)

p = ggplot(wb_results, aes(gene_effect, inclusion))
p = p + geom_point(aes(color = as.factor(inclusion > 0.90)))
p = p + geom_text_repel(aes(label = highlight_gene))
p = p + theme(legend.position = 'none')
p = p + scale_color_manual(values = cbb[c(2, 1)])
p = p + xlab('gene effect size')
p = p + ylab('posterior inclusion probability')
save_plot(p, file = '../img/volcano_linear.png', base_height = 10)

p = ggplot(wb_results, aes(gene_effect, -log10(1 - inclusion)))
p = p + geom_point(aes(color = as.factor(inclusion > 0.90)))
p = p + geom_text_repel(aes(label = highlight_gene))
p = p + theme(legend.position = 'none')
p = p + scale_color_manual(values = cbb[c(2, 1)])
p = p + xlab('gene effect size')
p = p + ylab('-log10(1 - posterior inclusion probability)')
save_plot(p, file = '../img/volcano_minus_log10.png', base_height = 10)
