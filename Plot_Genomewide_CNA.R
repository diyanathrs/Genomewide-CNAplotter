# Generating co-occurance figs for AML genome-wide CNA paper by Jim Allan's group, NuCancer
library(ggplot2)
library(dplyr)
library(BSgenome)
library(data.table) 
library(GenomicRanges)
library(karyoploteR)
library(plyranges)
library(regioneR)
library(locuszoomr)

################
## New plot ###
################

# load data
files.2.plt <- list.files('co_occur/','.txt')

lapply(files.2.plt[2:8], plot_cooccr.focal2, wh.len=35)
lapply(files.2.plt[10:26], plot_cooccr.focal2, wh.len=35)

plot_cooccr.focal2(files.2.plt[15], wh.len = 35, var.peak = 30) #use 25-45 - issues with 4


plot_cooccr.cna<- function(file, wh.len, var.peak) {
  dat <- read.delim(paste0('co_occur/',file), header = T)
  # remove smaller loh
  range(dat$Length)
  loh <- dat %>% dplyr::filter(Event %in% c('Allelic Imbalance', 'LOH')) %>% dplyr::filter(Length > 10e6)
  #min.length <- 0.5e6
  #dat <- subset(dat, dat$Length > min.length)
  dat <- subset(dat, !dat$Event %in% c('Allelic Imbalance', 'LOH'))
  # sort events
  dat$Event <- factor(dat$Event, levels = c('CN Loss','CN Gain', 'High Copy Gain','Homozygous Copy Loss'))
  dat <- dat[order(dat$Event),]
  table(dat$Event)
  table(loh$Event)
  range(loh$Length)
  range(dat$Length)
  
  name <- gsub(x = file,' ','_')
  name <- gsub(x = name,'.txt','')
  
  #to Granges
  gr.all <- toGRanges(dat$Chromosome.Region)
  gr.all$sample <- dat$Sample
  gr.all$event <- dat$Event
  #create vector for cols
  gr.all$color <- rep('#2476e0',length(gr.all))
  gr.all$color[gr.all$event %in% c('CN Loss')] <- '#db3d3d'
  gr.all$color[gr.all$event=='Homozygous Copy Loss'] <- '#820707' # 820707
  gr.all$color[gr.all$event=='High Copy Gain'] <- '#092073'
  
  if (nrow(loh) > 0) {
    gr.loh <- toGRanges(loh$Chromosome.Region)
    gr.loh$sample <- loh$Sample
    gr.loh$event <- gr.loh$Event}
  length(unique(gr.all$sample))
  # remove event for testing
  # gr dummy for white 0 line
  chr.lens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
  # Create GRanges covering each chromosome
  gr.dummy <- GRanges(seqnames = names(chr.lens[1:22]), ranges = IRanges(start = 1, end = chr.lens[1:22]))
  gr.dummy <- gr.dummy+1e6
  
  table(gr.all$event)  
  sm.lst <- unique(gr.all$sample)
  
  # divide gr.all to small and large 
  # expand <2mb regions for 5mb
  # keep only recurrent focals - do this to gr.all
  # select peaks within 70% loss.mx/gain.mx
  gr.small <- gr.all[width(gr.all) < 6e6]
  gr.large <- gr.all[width(gr.all) > 6e6]
  
  # test <- data.frame(rec.filtered)
  peaks.loss <- count_overlaps(disjoin(gr.all), gr.all[gr.all$event %in% c('CN Loss','Homozygous Copy Loss')])
  peaks.gain <- count_overlaps(disjoin(gr.all), gr.all[gr.all$event %in% c('CN Gain', 'High Copy Gain')])
  if (missing(var.peak)) {
    peak <- min(max(peaks.loss), max(peaks.gain))
    #add variable peak arg
  } else {
    peak <- var.peak }
  loss.rec <- disjoin(gr.all)[peaks.loss >= peak]
  gain.rec <- disjoin(gr.all)[peaks.gain >= peak]
  
  # new method - create disjs count overlaps, discard low peaks then reduce
  # range.loss <- count_overlaps(disjoin(gr.all), gr.all[gr.all$event %in% c('CN Loss','Homozygous Copy Loss')])
  # peaks.loss <- count_overlaps(range.loss, gr.all[gr.all$event %in% c('CN Loss','Homozygous Copy Loss')])
  # loss.rec <- range.loss[peaks.loss >= max(peaks.loss)*0.2]
  # range.gain <- GenomicRanges::reduce(gr.all[gr.all$event %in% c('CN Gain', 'High Copy Gain')])
  # peaks.gain <- count_overlaps(range.gain, gr.all[gr.all$event %in% c('CN Gain', 'High Copy Gain')])
  # gain.rec <- range.gain[peaks.gain >= max(peaks.gain)*0.2]
  
  # combine both loss.rec and gain rec then filter by comb ranges > gr.small > flank gr.small
  #rec.filtered <- c(loss.rec, gain.rec)
  #rec.filtered <- resize(rec.filtered, width = 1e6,  fix="center")
  rec.filtered <- GenomicRanges::reduce(c(loss.rec, gain.rec))
  gr.small.fil <- filter_by_overlaps(gr.small, rec.filtered)
  # merge fil'ed small rec cnas back to filtered gr.all (gr.l
  gr.all <- c(gr.large, gr.small.fil) # this is for density plot so no need to flank
  #test <- data.frame(gr.all)
  # flank the size of small rec
  range(width(gr.small.fil))
  gr.small.fil <- flank(gr.small.fil, width = 3e6, both = T)
  # remove
  #gr.small.fil <- gr.small.fil[!gr.small.fil$event=='CN Gain']
  
  # for the ymax have to get the highest peak for either gain or loss.
  # make disjoints and count overlaps
  loss.mx <- max(count_overlaps(disjoin(gr.all), gr.all[gr.all$event=='CN Loss']))
  gain.mx <- max(count_overlaps(disjoin(gr.all), gr.all[gr.all$event=='CN Gain']))
  hloss.mx <- max(count_overlaps(disjoin(gr.all), gr.all[gr.all$event=='Homozygous Copy Loss']))
  hgain.mx <- max(count_overlaps(disjoin(gr.all), gr.all[gr.all$event=='High Copy Gain']))
  # take max
  y.max <- max(c(loss.mx, gain.mx))
  
  # change plot params based on sample length
  pp <- getDefaultPlotParams(plot.type=3) # 60 and 120
  pp$topmargin <- 10
  pp$bottommargin <- 5
  pp$data1height <- min(length(sm.lst)*0.9, 22) # only have a max of 25
  pp$data2height <- max(length(sm.lst)*1.4, 20) #
  pp$data1inmargin <- 0.05
  pp$data2inmargin <- 0.05
  pp$leftmargin <- 0.055
  pp$rightmargin <- 0.02
  cyto.col <- getCytobandColors(color.schema = 'only.centromeres')
  cyto.col[10] <- '#000000'
  
  #pdf(file = paste0('co_occr_out/',name,'.pdf'),width = 11, height = max(length(sm.lst)*0.1,4))
  pdf(file = paste0('co_occr_out/',name,'_2.pdf'),width = 11, height = 3+length(sm.lst)*0.1)
  
  kp <- plotKaryotype(plot.type=3, plot.params = pp, chromosomes = c("autosomal"), 
                      labels.plotter = NULL, ideogram.plotter = NULL)
  
  kpAddCytobandsAsLine(kp, lwd = 12, color.table = cyto.col)
  #kpAddChromosomeSeparators(kp, lty = 2)
  kpAddChromosomeNames(kp, chr.names = paste0(c(1:22)),yoffset = 1, cex=0.8)
  kpDataBackground(kp, data.panel = 1, color = '#ededed')
  kpDataBackground(kp, data.panel = 2, color = '#ededed')
  kpAddMainTitle(kp, name, cex=0.7)
  # kpAddBaseNumbers(kp)
  # only keep focal smalls in gr.all
  try(kpPlotCoverage(kp, gr.all[gr.all$event %in% c('CN Gain')], r0 = 0.52, r1=1, col = '#2476e0', show.0.cov = T, 
                     ymax = y.max), silent = T)
  try(kpPlotCoverage(kp, gr.all[gr.all$event %in% c('High Copy Gain')], r0 = 0.52, r1=1, col = '#092073', 
                     show.0.cov = T, ymax = y.max), silent = T)
  kpPlotCoverage(kp, gr.all[gr.all$event == 'CN Loss'], r0 = 0.48, r1=0, col = '#db3d3d',show.0.cov = T, 
                 ymax = y.max)
  try(kpPlotCoverage(kp, gr.all[gr.all$event == 'Homozygous Copy Loss'], r0 = 0.48, r1=0, col = '#820707',
                     show.0.cov = T, ymax = y.max), silent = T)
  kpPlotCoverage(kp, gr.dummy, r0 = 0.49, r1=0, col = '#ededed', show.0.cov = T, ymax = wh.len) # if len is 20 =30, 50=40 
  kpPlotCoverage(kp, gr.dummy, r0 = 0.51, r1=1, col = '#ededed', show.0.cov = T, ymax = wh.len)
  if (nrow(loh) > 0) {
    kpPlotRegions(kp, gr.loh, r0 = 0.49, r1=0.51, col = '#dba825' , data.panel=1, avoid.overlapping = F, num.layers = 1)}
  
  bar.height <- 0 # only have border for larger alterations
  gr.high <- gr.all[gr.all$event %in% c('High Copy Gain', 'Homozygous Copy Loss')]
  
  for (sm in seq_along(sm.lst)) {
    r0 <- bar.height 
    r1 <- bar.height + 1/length(sm.lst)
    # kpPlotRegions(kp, gr.all[gr.all$sample==sm.lst[sm]], col = gr.all[gr.all$sample==sm.lst[sm]]$color, 
    #               r0 = r0, r1 = r1, data.panel = 2, num.layers = 1,avoid.overlapping = F, 
    #               border = 'gray') #border = '#ededed'
    if (nrow(loh) > 0) {
      kpPlotRegions(kp, gr.loh[gr.loh$sample==sm.lst[sm]], col = '#dba825', border = '#ededed',
                    r0 = r0, r1 = r1, data.panel = 2, avoid.overlapping = F) }
    kpPlotRegions(kp, gr.small.fil[gr.small.fil$sample==sm.lst[sm]], col = gr.small.fil[gr.small.fil$sample==sm.lst[sm]]$color, 
                  r0 = r0, r1 = r1, data.panel = 2, avoid.overlapping = F, border = '#ededed')
    kpPlotRegions(kp, gr.large[gr.large$sample==sm.lst[sm]], col = gr.large[gr.large$sample==sm.lst[sm]]$color, 
                  r0 = r0, r1 = r1, data.panel = 2, avoid.overlapping = F, border = '#ededed') 
    # kpPlotRegions(kp, gr.high[gr.high$sample==sm.lst[sm]], col = gr.high[gr.high$sample==sm.lst[sm]]$color, 
    #               r0 = r0, r1 = r1, data.panel = 2, avoid.overlapping = F, border = '#ededed') 
    # check if cell line
    if (sm.lst[sm] == 'GSM5397277_4HF_SNP6_25') {
      kpAddLabels(kp, labels="MV4-11", r0=r0, r1=r1, data.panel = 2, cex=0.5) #GSM888153 HL-60
    } else if (sm.lst[sm] == 'GSM888549 OCI-AML3') {
      kpAddLabels(kp, labels="OCI-AML3", r0=r0, r1=r1, data.panel = 2, cex=0.5)
    } else if (sm.lst[sm] == 'GSM888153 HL-60') {
      kpAddLabels(kp, labels="HL-60", r0=r0, r1=r1, data.panel = 2, cex=0.5)
    }  else {
      kpAddLabels(kp, labels=paste("Sample",sm), r0=r0, r1=r1, data.panel = 2, cex=0.5)
    }
    bar.height <- r1
  }
  
  name
  dev.off()
}


