setwd('~/d/sci/src/prp_mrm')
options(stringsAsFactors=FALSE)
library(sqldf)
library(reshape2)
library(minpack.lm)
library(beeswarm)

percent = function(proportion,digits=2,format='fg') {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format=format),"%",sep="") ) )
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

expand_range = function(x,by=0.5) {
  xmin = min(x,na.rm=T) - by
  xmax = max(x,na.rm=T) + by
  return (c(xmin, xmax))
}

imgmode = 'png'
# imgmode = 'pdf'
imgsave = get(imgmode) # get the function that saves as that type of image
if (imgmode == 'png') {
  resx = 600 # multiplier for image width, height, and resolution
} else if (imgmode == 'pdf') {
  resx = 1
}


peptides = read.table('data/ref/peptides.tsv', sep='\t', header=T, comment.char='', quote='')
hupeps = subset(peptides, human=='yes')

codons = read.table('data/ref/codons.tsv', sep='\t', header=T, comment.char='', quote='')
elisa_data = read.table('~/d/sci/src/csf_prp_quantification/data/processed/csf.tsv',sep='\t',header=T)
clinsamp = read.table('data/samples/samples.tsv', sep='\t', header=T, comment.char='', quote='')
pparms = read.table('data/ref/mrm_pparms.tsv', sep='\t', header=T, comment.char='', quote='')

# get ELISA details including plate #
if (exists('elisa')) {
  rm(elisa)
}
plates = read.table('~/d/sci/src/csf_prp_quantification/data/elisa/meta/plates.tsv',sep='\t',header=TRUE)
plates = plates[-which(plates$plateno=='6_alt'),]
plates$plateno = as.integer(plates$plate)
for (plate in plates$plateno) {
  filename = paste('~/d/sci/src/csf_prp_quantification/data/elisa/processed/',formatC(plate,width=2,flag='0'),'_summary.tsv',sep='')
  if (file.exists(filename)) {
    temp = read.table(filename,sep='\t',header=TRUE,quote='',comment.char='')
    temp$plate = plate
    if (all(is.na(temp$flag))) {
      temp$flag = '' # when all are blank, it reads them as NA which is annoying, so fix it here.
    }
    if (exists('elisa')) {
      elisa = rbind(elisa,temp)
    } else {
      elisa = temp
    }
  }
}
clinsamp$elisa_plate = elisa$plate[match(clinsamp$sample, elisa$sample)]


# load all MRM data into one table
runs = read.table('data/meta/runs.tsv',sep='\t',header=T,quote='',comment.char='')
if (exists(mrm_data)) {
  rm(mrm_data)
}
for (i in 1:nrow(runs)) {
  run = runs$expt_no[i]
  meta = read.table(paste('data/meta/',runs$metafile[i],sep=''), sep='\t',header=T,quote='',comment.char='')
  proc = read.table(paste('data/processed/',runs$metafile[i],sep=''), sep='\t',header=T,quote='',comment.char='')
  proc$description = meta$description[match(proc$sample,meta$skyline_id)]
  if ('sample' %in% colnames(meta)) {
    proc$indiv_id = meta$sample[match(proc$sample,meta$skyline_id)]
  } else {
    proc$indiv_id = ''
  }
  proc$label = meta$label[match(proc$sample,meta$skyline_id)]
  proc$run = run
  if (i == 1) {
    mrm_data = proc
  } else {
    mrm_data = rbind(mrm_data, proc)
  }
}
# fix weird artifact - for files where 'flag' field is always empty it comes in as NA instead of ''
mrm_data$flag[is.na(mrm_data$flag)] = ''
nrow(mrm_data)

# cast everything to numerics so that sqlite can handle division
mrm_data$light_area = as.numeric(mrm_data$light_area)
mrm_data$n15_area = as.numeric(mrm_data$n15_area)
mrm_data$heavy_area = as.numeric(mrm_data$heavy_area)

# check how many entries have which modified sequence (how often is met-ox peptide used)
table(mrm_data$modified_sequence)

# check which samples used met-ox peptides
mrm_data[mrm_data$flag=='met-ox',c('run','label')]

# create ratio columns
colnames(mrm_data)
mrm_data$nh = as.numeric(mrm_data$n15_area) / as.numeric(mrm_data$heavy_area)
mrm_data$lh = as.numeric(mrm_data$light_area) / as.numeric(mrm_data$heavy_area)
mrm_data$ln = as.numeric(mrm_data$light_area) / as.numeric(mrm_data$n15_area)

# here: process clinical data
# - joining in ELISA and other covariates
# - joining in plotting parameters (x, color, etc.)

# prepare summary tables for clinical samples
clin_all = subset(mrm_data, run %in% 16:20 & indiv_id != '9-IPC')
# check sample size
length(unique(clin_all$indiv_id))

# clin will be a table with one row per individual and peptide (summarized across 2 technical replicates)
clin1 = sqldf("
select   run, indiv_id, peptide,
         avg(case when flag='' then light_area else null end) light,
         avg(case when flag='' then n15_area else null end) n15,
         avg(case when flag='' then heavy_area else null end) heavy,
         avg(case when flag='' then ln else null end) ln_mean,
         stdev(case when flag='' then ln else null end) ln_sd, 
         avg(case when flag='' then nh else null end) nh_mean,
         stdev(case when flag='' then nh else null end) nh_sd, 
         avg(case when flag='' then lh else null end) lh_mean,
         stdev(case when flag='' then lh else null end) lh_sd, 
         sum(case when flag='' then 1 else 0 end) n
from     clin_all
where    indiv_id in (select sample from clinsamp)
and      peptide in (select peptide from peptides where human = 'yes')
group by 1, 2, 3
;")

# note that 12 replicates, including both replicates of just 1 sample, are met-ox'ed
# therefore for clin_all the sum of the "n" column is 648 instead of 660
# and for clin1, the nrow() is 329 instead of 330
# to allow the rank/partition trick later on to work, we need to have a null row
# present for the sample with missing VVEQ data
# clin_all$indiv_id[clin_all$flag=='met-ox'] # IZ_408 is the sample where both were oxidized

clin = sqldf("
select   e.sample indiv_id, e.prp_ngml elisa_ngml, e.prp_flag elisa_flag,
         e.hb_ngml, e.hb_flag, e.dc_mgml, e.dc_flag,
         e.prion_category, e.sample_source source, e.elisa_plate,
         c.run, c.peptide, c.light, c.n15, c.heavy, c.ln_mean, c.ln_sd,
         c.nh_mean, c.nh_sd, c.lh_mean, c.lh_sd, c.n
from     clinsamp e, clin1 c
where    e.sample = c.indiv_id -- use this for joining
;")

clin$ncorder = peptides$ncorder[match(clin$peptide, peptides$peptide)]

# this should be 330
nrow(clin)

clin$xpep = peptides$ncorder[match(clin$peptide, peptides$peptide)]
clin$color_pep = peptides$color[match(clin$peptide, peptides$peptide)]
clin$xdiag = pparms$x[match(clin$prion_category, pparms$prion_category)]
clin$color_diag = pparms$color[match(clin$prion_category, pparms$prion_category)]

clin$xall = clin$xpep*4 + clin$xdiag


# get the N of each type of clinical sample
if (file.exists('data_received/csf-sample-details.txt')) {
  sampdeets = read.table('data_received/csf-sample-details.txt',sep='\t',header=T,quote='',comment.char='')
  sampdeets$indiv_id[!grepl('-',sampdeets$indiv_id)] = paste('IZ',sampdeets$indiv_id[!grepl('-',sampdeets$indiv_id)],sep='_')
  sampdeets = subset(sampdeets, indiv_id %in% clin$indiv_id)
  table(clin$prion_category)/6
  table(sampdeets$diag_subcat)
  table(sampdeets$diagnosis[sampdeets$diag_subcat=='other'])
}


# prepare summary tables for selectivity samples in analytical validation stages
selectivity_samples = read.table('data/meta/selectivity_samples.tsv',sep='\t',header=T)

selectivity_data = sqldf("
                         select   s.species, m.*
                         from     mrm_data m, selectivity_samples s
                         where    m.run = s.run and m.label = s.sample
                         ;")

selectivity_data$expectedness = ''
for (species in unique(selectivity_samples$species)) {
  selectivity_data$expectedness[selectivity_data$species==species & selectivity_data$peptide %in% peptides$peptide[peptides[,species]=='yes']] = 'sequence-matched'
  selectivity_data$expectedness[selectivity_data$species==species & selectivity_data$peptide %in% peptides$peptide[peptides[,species]=='no']] = 'not matched'
}

# remove selectivity data for VVEQMCITQYEK from expts below #8 because we didn't yet have that peptide at the time
selectivity_data = selectivity_data[!(selectivity_data$peptide=='VVEQMCITQYEK' & selectivity_data$run < 8),]

# get ready to plot selectivity
selectivity_legend = data.frame(label=c('sequence-matched','not matched'), col=c('#FF7216','#C7C7C7'))
selectivity_tranches = data.frame(species=c('human','cyno','rat','mouse'),tranche=0:-3)
selectivity_data$tranche = selectivity_tranches$tranche[match(selectivity_data$species, selectivity_tranches$species)]
tranche_size = nrow(peptides) + 1
selectivity_data$allorder = peptides$allorder[match(selectivity_data$peptide,peptides$peptide)]
selectivity_data$y = tranche_size * selectivity_data$tranche - peptides$allorder[match(selectivity_data$peptide,peptides$peptide)]


# Table S1. Performance of all peptides in analytical validation stages
# Rows: 9 peptides + met-ox forms (group by peptide, modified_sequence or flag)
# Cols: retention time mean +- sd, ion, count(*) n, tech replicate cv, mean L:H ratio

selectivity_data$id_wo_replicate = gsub('-[123]','',selectivity_data$id)
selectivity_data$id_wo_replicate[selectivity_data$id %in% c('G1','G2')] = 'G'
selectivity_data$id_wo_replicate[selectivity_data$id %in% c('IPC1','IPC2')] = 'IPC'

# compute technical duplicate mean CV for each peptide based on selectivity data
techcv1 = sqldf("
                select   peptide, run, id_wo_replicate, avg(lh) mean_lh, stdev(lh) sd_lh, count(*) n
                from     selectivity_data s
                where    expectedness = 'sequence-matched'
                group by 1, 2, 3
                having   count(*) > 1
                order by 1, 2, 3
                ;")
techcv1$cv = techcv1$sd_lh/techcv1$mean_lh

mean_lhrs1 = sqldf("
                   select   peptide, run, id_wo_replicate, avg(lh) mean_lh
                   from     selectivity_data s
                   where    expectedness = 'sequence-matched'
                   group by 1, 2, 3
                   order by 1, 2, 3
                   ;")

mean_lhrs2 = sqldf("
                   select   peptide, avg(mean_lh) mean_mean_lh, count(*) n_samp_lh
                   from     mean_lhrs1
                   group by 1
                   order by 1
                   ;")

techcv2 = sqldf("
                select   p.allorder, t.peptide, avg(t.cv) mean_cv, count(*) n_samp, sum(n) n_total_rep
                from     techcv1 t, peptides p
                where    t.peptide = p.peptide
                group by 1
                order by 1
                ;")

techcv2$mean_mean_lh = mean_lhrs2$mean_mean_lh[match(techcv2$peptide, mean_lhrs2$peptide)]
techcv2$n_samp_lh = mean_lhrs2$n_samp_lh[match(techcv2$peptide, mean_lhrs2$peptide)]

techcv2$mean_mean_lh = signif(techcv2$mean_mean_lh, digits=2)

techcv2


table_s5_prep = sqldf("
                      select   allorder, peptide, modified_sequence, ion, count(*) n_rt, avg(retention_time) mean_rt, stdev(retention_time) sd_rt
                      from     selectivity_data
                      where    expectedness == 'sequence-matched'
                      and      run != 6
                      group by 1, 2, 3, 4
                      order by 1, 2, 3, 4
                      ;")

table_s5_prep$mean_lh = techcv2$mean_mean_lh[match(table_s5_prep$peptide, techcv2$peptide)]
table_s5_prep$n_samp_lh = techcv2$n_samp_lh[match(table_s5_prep$peptide, techcv2$peptide)]
table_s5_prep$mean_cv = techcv2$mean_cv[match(table_s5_prep$peptide, techcv2$peptide)]
table_s5_prep$n_samp = techcv2$n_samp[match(table_s5_prep$peptide, techcv2$peptide)]

table_s5_prep$retention_time = paste(formatC(table_s5_prep$mean_rt, format='f', digits=1), '±', formatC(table_s5_prep$sd_rt, format='f', digits=1), ' (N=', table_s5_prep$n_rt, ')', sep='')
table_s5_prep$lh = paste(formatC(table_s5_prep$mean_lh, format='f', digits=1), ' (N=', table_s5_prep$n_samp_lh, ')', sep='')
table_s5_prep$cv = paste(percent(table_s5_prep$mean_cv), ' (N=', table_s5_prep$n_samp, ')', sep='')

table_s5 = table_s5_prep[,c('allorder','peptide','ion','retention_time','lh','cv')]
table_s5

write.table(table_s5, 'display_items/script_generated/table-s5.tsv',sep='\t',col.names=T,row.names=F,quote=F)

# stats for table S1 legend
length(unique(paste(selectivity_data$run, selectivity_data$label)))
length(unique(paste(selectivity_data$run, selectivity_data$id_wo_replicate)))
sqldf("select species, count(*) n from (select species, id_wo_replicate from selectivity_data group by 1,2) group by 1;")

# calculate concentration of heavy peptide
# 100 fmol per 30 uL sample
moles = 100 * 1e-15
vol = 30 * 1e-6
molarity = moles/vol # mol/L
approx_mw = 22000 # g/mol
masspervol = molarity * approx_mw # g/L
ngml = masspervol * (1e9) / (1e3) # convert g/L to ng/ML


# Figure 5. Application of the PrP MRM assay for preclinical species of interest
# A: grid of expected peptide results based on sequence match
# B (was S5B): Selectivity data collapsed by peptide
# C (was S5E): Dilution linearity of mouse peptides in brain homogenate
imgsave(paste('display_items/script_generated/figure-5.',imgmode,sep=''),width=6.5*resx,height=3*resx,res=resx)
layout_matrix = matrix(c(1,1,1,2,2,3,3),nrow=1,byrow=T)
layout(layout_matrix)

y_axis_data = sqldf("
                    select   y, peptide, expectedness
                    from     selectivity_data s
                    group by 1, 2, 3
                    order by 1, 2, 3
                    ;")
y_axis_data$color = selectivity_legend$col[match(y_axis_data$expectedness,selectivity_legend$label)]
selectivity_data$color = selectivity_legend$col[match(selectivity_data$expectedness,selectivity_legend$label)]
ylims = c(min(y_axis_data$y)-.5, max(y_axis_data$y)+.5) # for one species: c(min(peptides$allorder)-.5, max(peptides$allorder)+.5)
selectivity_data$y_collapsed = -peptides$allorder[match(selectivity_data$peptide,peptides$peptide)]
collapsed_ylims=c(-9.5,-.5)

par(mar=c(5,10,5,2))
xlims = c(0.5,3.0)
species = c('human','cynomolgus','mouse','rat') # colnames(peptides)[4:7]
colnames(peptides)[4:7] = species
species_x = c(1,1.5,2,2.5)
plot(NA,NA,xlim=xlims,ylim=collapsed_ylims,ann=FALSE,axes=FALSE,xaxs='i',yaxs='i')
mtext(side=2,at=-peptides$allorder,text=peptides$peptide,col=peptides$color,font=1,las=2,cex=0.8)
species_matrix = as.matrix(peptides[,species] == 'yes')
mtext(side=1, line=0.5, text='sequence expected')
par(xpd=T)
points(x=c(1,2), y=rep(min(collapsed_ylims),2)-1.5, pch=c(15,4), lwd=c(NA,2), col=selectivity_legend$col, cex=1.1)
text(x=c(1,2), y=rep(min(collapsed_ylims),2)-1.5, pos=4, col=selectivity_legend$col, labels=c('yes','no'), cex=1.1)
text(x=species_x, y=rep(max(collapsed_ylims),4)+0.5, labels=species, adj=c(0,0), srt=45)
par(xpd=F)
for (row in 1:nrow(species_matrix)) {
  for (col in 1:ncol(species_matrix)) {
    if (species_matrix[row,col]) {
      points(species_x[col], -row, pch=15, col=selectivity_legend$col[selectivity_legend$label=='sequence-matched'])
    } else {
      points(species_x[col], -row, pch=4, lwd=2, col=selectivity_legend$col[selectivity_legend$label=='not matched'])
    }
  }
}
mtext('A', side=3, cex=2, adj = -0.9, line = 0.3)

par(mar=c(5,1,5,2))
xlims = c(floor(min(log10(selectivity_data$lh))), ceiling(max(log10(selectivity_data$lh))))
xats = min(xlims):max(xlims)
plot(NA,NA,xlim=xlims,ylim=collapsed_ylims,ann=FALSE,axes=FALSE,xaxs='i',yaxs='i')
abline(v=xats, lwd=0.25)
axis(side=1,at=xats)
mtext(side=1, line=2.5, text='observed abundance')
mtext(side=1, line=3.5, text='log10(L:H)', cex=0.8)
points(x=log10(selectivity_data$lh), y=selectivity_data$y_collapsed, col=selectivity_data$color, pch=20)
mtext('B', side=3, cex=2, adj = -0.2, line = 0.3)

mocurve = subset(mrm_data, run == 9)
mocurve$fracwt = as.numeric(gsub('x.*','',mocurve$description))
cv_step1 = sqldf("
                 select   peptide, fracwt, avg(lh) mean_lh, stdev(lh) sd_lh
                 from     mocurve
                 where    fracwt > 0
                 group by 1, 2
                 order by 1, 2
                 ;")
cv_step2 = sqldf("
                 select   c.peptide, avg(c.sd_lh / c.mean_lh) mean_cv
                 from     cv_step1 c, peptides p
                 where    c.peptide = p.peptide
                 and      p.mouse = 'yes'
                 group by 1
                 order by 1
                 ;")
cv_step2
#         peptide    mean_cv
# 1     ESQAYYDGR 0.08301260
# 2    GENFTETDVK 0.14291497
# 3 RPKPGGWNTGGSR 0.37349903
# 4  VVEQMCVTQYQK 0.06612986
# 5   YPGQGSPGGNR 0.07887789
best_mo_peptides = subset(cv_step2, mean_cv < .10)
best_mo_peptides$col = peptides$color[match(best_mo_peptides$peptide,peptides$peptide)]
mocurve = subset(mocurve, peptide %in% best_mo_peptides$peptide)
mocurve$col = best_mo_peptides$col[match(mocurve$peptide, best_mo_peptides$peptide)]

mo_rel = sqldf("
               select   mall.peptide, mall.col, mall.fracwt, mall.lh, m100.mean_lh, mall.lh / m100.mean_lh rel
               from     mocurve mall, (select peptide, avg(lh) mean_lh from mocurve where fracwt == 1 group by 1 order by 1) m100
               where    mall.peptide = m100.peptide
               order by 1
               ;")

best_mo_peptides$intercept = 0
best_mo_peptides$slope = 0
best_mo_peptides$adjr2 = 0
for (i in 1:nrow(best_mo_peptides)) {
  m = lm(rel ~ fracwt, data=subset(mo_rel, peptide==best_mo_peptides$peptide[i]))
  best_mo_peptides$intercept[i] = coefficients(m)["(Intercept)"]
  best_mo_peptides$slope[i] = coefficients(m)["fracwt"]
  best_mo_peptides$adjr2[i] = summary(m)$adj.r.squared
}
percent(range(best_mo_peptides$adjr2),digits=3)

par(mar=c(5,3,5,3))
plot(NA, NA, xlim=c(-0.05,1.05), ylim=c(0,1.1), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=unique(mo_rel$fracwt), labels=NA, lwd=1, lwd.ticks=1, tck=-0.02)
axis(side=1, at=c(0,.5,1), labels=c('100%\nKO','50/50\nmix','100%\nWT'), lwd=0, lwd.ticks=1, tck=-0.05, cex.axis=0.8)
axis(side=2, at=(0:4)/4, labels=percent((0:4)/4), lwd=1, lwd.ticks=1, las=2)
mtext(side=1, line=2.5, text='BH mixture proportion')
mtext(side=2, line=3, text='normalized L:H ratio')
set.seed(1)
points(jitter(mo_rel$fracwt,amount=0.025), mo_rel$rel, col=mo_rel$col, pch=20)
for (i in 1:nrow(best_mo_peptides)) {
  abline(a=best_mo_peptides$intercept[i], b=best_mo_peptides$slope[i], lwd=1, col=best_mo_peptides$col[i])
}
legend('topleft', best_mo_peptides$peptide, pch=20, lwd=1, col=best_mo_peptides$col, text.col=best_mo_peptides$col, text.font=1, cex=0.65, bty='n')
mtext('C', side=3, cex=2, adj = -0.2, line = 0.3)

dev.off()






# Figure S5. QC & partial analytical validation of the PrP MRM assay.
# A-C (was S7A, B, E) : peptide abundance by day
# D (was S10A): Correlation of peptides with one another
# E (was S5D): Linearity of hi- into low human CSF sample

imgsave(paste('display_items/script_generated/figure-s05.',imgmode,sep=''),width=6.5*resx,height=6.5*resx,res=resx)

layout_matrix = matrix(c(1,1,2,2,3,3,4,4,4,4,4,4,5,5,5,6,6,6),nrow=3,byrow=T)
layout(layout_matrix, heights=c(1,.5,1))

# this figure requires access to raw Skyline data
# load all MRM data into one table
runs = read.table('data/meta/runs.tsv',sep='\t',header=T,quote='',comment.char='')
if (exists('mrm_raw')) {
  rm(mrm_raw)
}
for (batch in 1:5) {
  rowno = which(runs$expt_no == batch + 15)
  rawdat = read.table(paste('data/skyline/',runs$datafile[rowno],sep=''), sep=',',header=T,quote='',comment.char='')
  rawdat$batch = batch
  if (exists('mrm_raw')) {
    mrm_raw = rbind(mrm_raw, rawdat)
  } else {
    mrm_raw = rawdat
  }
}
colnames(mrm_raw) = gsub('[^a-z0-9_]','_',tolower(colnames(mrm_raw)))
mrm_raw = subset(mrm_raw, peptide_sequence %in% peptides$peptide[peptides$human=='yes'])

colnames(mrm_raw)
mrm_raw$light_area = as.numeric(mrm_raw$light_area)
mrm_raw$x15n_area = as.numeric(mrm_raw$x15n_area)
mrm_raw$heavy_area = as.numeric(mrm_raw$heavy_area)

pep_batch_totals = sqldf("
                         select   batch,
                         peptide_sequence,
                         peptide_modified_sequence,
                         sum(light_area) l,
                         sum(x15n_area) n15,
                         sum(heavy_area) h,
                         sum(x15n_area/heavy_area) nh,
                         sum(light_area/x15n_area) ln,
                         sum(light_area/heavy_area) lh
                         from     mrm_raw
                         group by 1, 2, 3
                         order by 1, 2, 3
                         ;")

batch_totals = sqldf("
                     select   batch,
                     sum(light_area) l,
                     sum(x15n_area) n15,
                     sum(heavy_area) h,
                     sum(x15n_area/heavy_area) nh,
                     sum(light_area/x15n_area) ln,
                     sum(light_area/heavy_area) lh
                     from     mrm_raw
                     group by 1
                     order by 1
                     ;")

pep_batch = sqldf("
                  select   b.batch, p.ncorder,p.peptide, pb.peptide_modified_sequence, p.color,
                  pb.l/b.l lx,
                  pb.n15/b.n15 nx,
                  pb.h/b.h hx,
                  pb.nh/b.nh nhx,
                  pb.ln/b.ln lnx,
                  pb.lh/b.lh lhx
                  from     peptides p, pep_batch_totals pb, batch_totals b
                  where    p.peptide = pb.peptide_sequence
                  and      pb.batch = b.batch
                  order by 1, 2
                  ;")

plot_params = sqldf("
                    select   m.peptide_modified_sequence, p.peptide, p.ncorder, p.color
                    from     mrm_raw m, peptides p
                    where    m.peptide_sequence = p.peptide
                    group by 1, 2
                    order by 3
                    ;")
plot_params$order_w_metox  = 1:7
plot_params$angle = NA
plot_params$density = NA
plot_params$angle[plot_params$peptide_modified_sequence=='VVEQM[+16]C[+57]ITQYER'] = 45
plot_params$density[plot_params$peptide_modified_sequence=='VVEQM[+16]C[+57]ITQYER'] = 40


multiplot_params = data.frame(var=c('lx','nx','lnx'),
                              disp=c('light','15N','light:15N'))
par(mar=c(4,5,4,2))
for (i in 1:nrow(multiplot_params)) {
  var = multiplot_params$var[i]
  x = dcast(pep_batch, formula = peptide_modified_sequence ~ batch, value.var=var)
  x$order_w_metox = plot_params$order_w_metox[match(x$peptide_modified_sequence, plot_params$peptide_modified_sequence)]
  x = x[order(x$order_w_metox),]
  mat = as.matrix(x[,as.character(1:5)])
  barplot(mat, col=plot_params$color, density=plot_params$density, angle=plot_params$angle, yaxt='n', border=NA)
  axis(side=2, at=(0:4)/4, labels=percent((0:4)/4), lwd=1, lwd.ticks=1, las=2)
  mtext(side=1, line=2.5, text='day number')
  mtext(side=2, line=3.0, text='% of total')
  mtext(side=3, line=1, text=multiplot_params$disp[i])
  mtext(LETTERS[i], side=3, cex=2, adj = 0,0, line = 1.5)
}

# wide empty plot to hold master legend
par(mar=c(1,2,1,2))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', axes=F, ann=F)
legend(x=0.0, y=1, cex=1.6, legend=plot_params$peptide_modified_sequence[1:4], fill=plot_params$color[1:4], text.col=plot_params$color[1:4], text.font=2, density=plot_params$density[1:4], angle=plot_params$angle[1:4], bty='n')
legend(x=0.4, y=1, cex=1.6, legend=plot_params$peptide_modified_sequence[5:7], fill=plot_params$color[5:7], text.col=plot_params$color[5:7], text.font=2, density=plot_params$density[5:7], angle=plot_params$angle[5:7], bty='n')


scatter_data = sqldf("
                     select   refpep.peptide refpeptide,
                     other.peptide otherpeptide,
                     other.color_pep color,
                     refpep.ln_mean refln,
                     refpep.lh_mean reflh,
                     other.ln_mean othln,
                     other.lh_mean othlh
                     from     clin refpep, clin other
                     where    refpep.indiv_id = other.indiv_id
                     and      refpep.peptide = 'VVEQMCITQYER'
                     and      other.peptide != 'VVEQMCITQYER'
                     ;")

lims = c(0,55)
xats = 0:60
xbigs = (0:6)*10
par(mar=c(4,4,4,2))
plot(NA, NA, xlim=lims, ylim=lims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=xats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=xbigs, labels=xbigs, lwd=0, lwd.ticks=1, tck=-0.050)
axis(side=2, at=xats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025, las=2)
axis(side=2, at=xbigs, labels=xbigs, lwd=0, lwd.ticks=1, tck=-0.050, las=2)
mtext(side=1, line=2.0, text='reference peptide L:15N ratio')
mtext(side=2, line=2.0, text='other peptide L:15N ratio')
mtext(side=3, line=0, text='L:15N')
points(x=scatter_data$refln, y=scatter_data$othln, col=scatter_data$color, pch=20, cex=0.75)
legend('topleft', legend=hupeps$peptide[c(1,2,3,4,6)], pch=20, col=hupeps$color[c(1,2,3,4,6)], text.col=hupeps$color[c(1,2,3,4,6)], text.font=2, bty='n', cex=0.8)
mtext('D', side=3, cex=2, adj = 0,0, line = 1.5)


expt11 = subset(mrm_data, run==11 & peptide %in% peptides$peptide[peptides$human=='yes'])
hilo = subset(expt11, grepl('^A',label))
hilo$proportion_a = as.numeric(gsub('-[12]','',gsub('A','',hilo$label)))/100
hilo$col = peptides$color[match(hilo$peptide, peptides$peptide)]

hilo_rel = sqldf("
                 select   hall.peptide, hall.proportion_a, hall.lh, h100.mean_lh, hall.lh / h100.mean_lh rel
                 from     hilo hall, (select peptide, avg(lh) mean_lh from hilo where proportion_a == 1 group by 1 order by 1) h100
                 where    hall.peptide = h100.peptide
                 order by 1
                 ;")

hilo_peptides = peptides[peptides$human=='yes',c('ncorder','peptide','color')]
hilo_peptides$intercept = 0
hilo_peptides$slope = 0
hilo_peptides$adjr2 = 0
for (i in 1:nrow(hilo_peptides)) {
  m = lm(rel ~ proportion_a, data=subset(hilo_rel, peptide==hilo_peptides$peptide[i]))
  hilo_peptides$intercept[i] = coefficients(m)["(Intercept)"]
  hilo_peptides$slope[i] = coefficients(m)["proportion_a"]
  hilo_peptides$adjr2[i] = summary(m)$adj.r.squared
}
percent(range(hilo_peptides$adjr2),digits=3)

par(mar=c(4,4,4,2))
plot(NA, NA, xlim=c(-0.05,1.05), ylim=c(0,1.05), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=c(0,.25,.5,.75,1), labels=c('100% low','','50/50 mix','','100% high'), lwd=1, lwd.ticks=1)
axis(side=2, at=(0:4)/4, labels=percent((0:4)/4), lwd=1, lwd.ticks=1, las=2)
mtext(side=1, line=2.5, text='mixture proportion of human CSF')
mtext(side=2, line=3, text='normalized L:H ratio')
set.seed(1) # for jitter
points(jitter(hilo_rel$proportion_a,amount=0.025), hilo_rel$rel, col=hilo$col, pch=20)
for (i in 1:nrow(hilo_peptides)) {
  abline(a=hilo_peptides$intercept[i], b=hilo_peptides$slope[i], lwd=1, col=hilo_peptides$color[i])
}
par(xpd=T)
legend(x=-0.05, y=1.1, peptides$peptide[peptides$human=='yes'],col=peptides$color[peptides$human=='yes'],text.font=2,text.col=peptides$color[peptides$human=='yes'],pch=20,lwd=1,cex=0.8,bty='n')
par(xpd=F)
mtext('E', side=3, cex=2, adj = -0.2, line = 0.3)


dev.off()








imgsave(paste('display_items/script_generated/figure-2.',imgmode,sep=''),width=6.5*resx,height=8*resx,res=resx)

layout_matrix = as.matrix(c(1,2,3,4),nrow=4,byrow=T)
layout(layout_matrix)

par(mar=c(2,4,3,1))

indsmry = sqldf("
                select   indiv_id, avg(light) l, avg(n15) n, avg(heavy) h, avg(ln_mean) ln, avg(nh_mean) nh, avg(lh_mean) lh
                from     clin
                group by 1 order by 1
                ;")

ratios = sqldf("
               select   c.indiv_id, h.ncorder, h.peptide, h.color, c.light/i.l lx, c.n15/i.n nx, c.heavy/i.h hx,
               c.ln_mean / i.ln lnx, c.nh_mean / i.nh nhx, c.lh_mean / i.lh lhx
               from     clin c, indsmry i, hupeps h
               where    c.indiv_id = i.indiv_id
               and      c.peptide = h.peptide
               order by 1, 2, 3
               ;")

xats = rep(1:9,10) * 10^(rep(0:9,each=9))
xbigs = 10^(0:9)
xbigs_labs = paste('1e',0:9,sep='')

paneladj = 0.0

range(ratios$lx, na.rm=T)
range(ratios$nx, na.rm=T)
range(ratios$lnx, na.rm=T)

ylims = c(3e3,3e7)
xlims = c(0.5, 6.5)

plot(NA, NA, xlim=xlims, ylim=ylims, log='y', axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=xats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=xbigs, labels=xbigs_labs, lwd=0, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=2.5, text='peak area')
for (indiv in unique(clin$indiv_id)) {
  subs = subset(clin, indiv_id==indiv)
  subs = subs[with(subs, order(ncorder)),]
  points(subs$ncorder, subs$light, type='l', lwd=0.25, col='#CCCCCC')
  points(subs$ncorder, subs$light, pch=20, col=subs$color_pep)
}
mtext(side=3, line=0, text='light', cex=1, font=1)
mtext('A', side=3, cex=2, adj = paneladj, line = 0.5)

ylims = c(3e3, 3e6)
plot(NA, NA, xlim=xlims, ylim=ylims, log='y', axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=xats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=xbigs, labels=xbigs_labs, lwd=0, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=2.5, text='peak area')
for (indiv in unique(clin$indiv_id)) {
  subs = subset(clin, indiv_id==indiv)
  subs = subs[with(subs, order(ncorder)),]
  points(subs$ncorder, subs$n15, type='l', lwd=0.25, col='#CCCCCC')
  points(subs$ncorder, subs$n15, pch=20, col=subs$color_pep)
}
mtext(side=3, line=0, text=expression(''^'15'*'N'), cex=1, font=1)
mtext('B', side=3, cex=2, adj = paneladj, line = 0.5)

ylims = c(0.1, 100)
r_ats = rep(1:9,10) * 10^(rep(-3:6,each=9))
r_bigs = 10^(-3:6)
r_bigs_labs = formatC(r_bigs, digits=0, format='fg', big.mark=',')
plot(NA, NA, xlim=xlims, ylim=ylims, log='y', axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=r_ats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=r_bigs, labels=r_bigs_labs, lwd=0, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=2, text='peak area ratio')
for (indiv in unique(clin$indiv_id)) {
  subs = subset(clin, indiv_id==indiv)
  subs = subs[with(subs, order(ncorder)),]
  points(subs$ncorder, subs$ln_mean, type='l', lwd=0.25, col='#CCCCCC')
  points(subs$ncorder, subs$ln_mean, pch=20, col=subs$color_pep)
}
mtext(side=3, line=0, text=expression('light:'^'15'*'N'), cex=1, font=1)
mtext('C', side=3, cex=2, adj = paneladj, line = 0.5)

# empty panel for bottom legend material
par(mar=c(2,4,0,1))
plot(NA, NA, xlim=xlims, ylim=c(0,1), axes=F, ann=F, xaxs='i')
par(xpd=T)
mtext(side=3, line=0.5, at=hupeps$ncorder, text=hupeps$peptide, col=hupeps$color, font=2, cex=0.6)
par(xpd=F)
# arrows to indicate N vs. C terminal
nctermcol = '#5C5C5C'
par(xpd=T)
yarr = 0.95
xprop = .15
arrows(x0=min(xlims), x1=min(xlims) + xprop*max(xlims), y0=yarr, code=1, length=0.1, angle=30, lwd=2, col=nctermcol)
text(x=min(xlims) + xprop*max(xlims), y=yarr, pos=4, labels='N terminus', col=nctermcol, font=2)
arrows(x0=max(xlims), x1=max(xlims) + min(xlims) - xprop*max(xlims), y0=yarr, code=1, length=0.1, angle=30, lwd=2, col=nctermcol)
text(x=max(xlims) + min(xlims) - xprop*max(xlims), y=yarr, pos=2, labels='C terminus', col=nctermcol, font=2)
par(xpd=F)

dev.off()



# Table 1. Quantification of peptides in human CSF
# Rows: 6 human peptides (asterisk on VVEQ count met-ox replicates)
# Cols: H, N, L, L:N, implied mean ng/mL, intra-run CV (all samples), inter-run CV (IPC)
table1_prep = sqldf("
select   peptide, 
         avg(ln_sd / ln_mean) mean_intrarun_cv
from     clin
group by 1
;")

table1_interindiv = sqldf("
select   peptide, 
         stdev(ln_mean) / avg(ln_mean) interindiv_cv 
from     clin
group by 1
;")


table1_interrun = sqldf("
select   p.ncorder, ipc.peptide, stdev(ipc.mean_ln)/avg(ipc.mean_ln) mean_interrun_cv
from     (select peptide, run, avg(ln) mean_ln from mrm_data where run in (16,17,18,19,20) and indiv_id=='9-IPC' group by 1, 2) ipc,
         peptides p
where    p.peptide = ipc.peptide and p.human = 'yes'
group by 1, 2
order by 1, 2
;")

table1_all = sqldf("
select   p.ncorder, 
         p.codons, 
         p.peptide, 
         t1a.mean_intrarun_cv, 
         t1b.mean_interrun_cv,
         t1c.interindiv_cv
from     table1_prep t1a, table1_interrun t1b, table1_interindiv t1c, peptides p
where    t1a.peptide = t1b.peptide and t1a.peptide = p.peptide and t1c.peptide = p.peptide
order by 1
;")

table1_disp = table1_all[,c('codons','peptide','mean_intrarun_cv','mean_interrun_cv','interindiv_cv')]
table1_disp$mean_ln = formatC(table1_disp$mean_ln, format='f', digits=1)
table1_disp$implied_ngml = formatC(table1_disp$implied_ngml, format='f', digits=0)
table1_disp$mean_intrarun_cv = gsub(' ','',paste(formatC(table1_disp$mean_intrarun_cv*100, digits=0, format='f'),"%",sep="") )
table1_disp$mean_interrun_cv = gsub(' ','',paste(formatC(table1_disp$mean_interrun_cv*100, digits=0, format='f'),"%",sep="") )
table1_disp$interindiv_cv = gsub(' ','',paste(formatC(table1_disp$interindiv_cv*100, digits=0, format='f'),"%",sep="") )
table1_disp

write.table(table1_disp, 'display_items/script_generated/table-1.tsv',sep='\t',col.names=T,row.names=F,quote=F)



#### Begin Table S4
table_s4_prep = sqldf("
                      select   peptide, run-15 batch, avg(nh_sd/nh_mean) mean_cv, avg(nh_mean) mean_nh
                      from     clin
                      where    nh_mean is not null and nh_sd is not null
                      group by 1, 2
                      order by 1, 2
                      ;")
table_s4_prep
table_s4_left = dcast(table_s4_prep, peptide ~ batch, value.var='mean_cv')
table_s4_right = dcast(table_s4_prep, peptide ~ batch, value.var='mean_nh')

for (i in 2:6) {
  table_s4_left[,i] = gsub(' ','',paste(formatC(table_s4_left[,i]*100, digits=1, format='f'),"%",sep="") )
}

for (i in 2:6) {
  table_s4_right[,i] = formatC(table_s4_right[,i], format='f', digits=2)
}

table_s4_joined = sqldf("
select   p.ncorder, p.peptide, l.*, r.*
from     peptides p, table_s4_left l, table_s4_right r
where    p.peptide = l.peptide and p.peptide = r.peptide
;")
table_s4_joined = table_s4_joined[,c(-3, -9)] # delete the extraneous "peptide" cols
table_s4_joined = table_s4_joined[,-1]

colnames(table_s4_joined) = c('peptide','cv1','cv2','cv3','cv4','cv5','mean1','mean2','mean3','mean4','mean5')

write.table(table_s4_joined, 'display_items/script_generated/table-s4.tsv', sep='\t', col.names=T, row.names=F, quote=F)
#### -- End Table S4





# calculate response factors and perform normalization
expt15 = subset(mrm_data, run == 15 & peptide %in% peptides$peptide[peptides$human=='yes'])
eparms = data.frame(label=LETTERS[1:12],spike=rep(c(0,2.4,24,240),each=3),x=rep(0:3,each=3))
expt15$spike = eparms$spike[match(expt15$label, eparms$label)]
expt15$x = eparms$x[match(expt15$label, eparms$label)]
expt15$pep_col = peptides$color[match(expt15$peptide, peptides$peptide)]
xlabs = sqldf("select x, spike from eparms group by 1, 2 order by 1, 2;")
xbreak = 0.5
xlims = c(-0.5, 3.5)

expt15$nl = 1/expt15$ln # 15N:L ratio is reciprocal of the L:15N ratio
true_n = 24.2 # true 15N spike concentration in clinical samples = 24.2 ng/mL
for (i in 1:nrow(hupeps)) { # for each peptide
  # use dose response experiment data (expt15) for normalization
  m = lm(spike ~ nl + 0, data = subset(expt15, spike > 0 & peptide==hupeps$peptide[i])) # linear model: spiked 15N concentration (2.4, 24, 240 ng/mL) as a function of 15N:L ratio
  hupeps$slope[i] = m$coefficients['nl'] # extract the slope for these curves for each peptide. range from 447 for VVEQ to 39 for RPKP
}
hupeps$rf = max(hupeps$slope) / hupeps$slope # set response factor equal to the max slope (the one for VVEQ...) divided by each peptide's slope
clin$rf = hupeps$rf[match(clin$peptide, hupeps$peptide)] # match the response factors into the clin table
clin$raw_ngml = clin$ln_mean * true_n # calculate the "raw" (non-normalized) ng/mL concentration as the L:15N ratio times the spike concentration (24.2 ng/mL)
clin$norm_ngml = clin$ln_mean * clin$rf * true_n # calculate the normalized ng/mL concentration as the raw concentration times the response factor

# note: some people found it non-intuitive that the linear model is spike ~ nl + 0 and not spike ~ ln + 0
# in other words, why is the curve fit to the 15N:L ratio and not the L:15N ratio?
# the answer is that if you do spike ~ ln + 0, the L:15N ratio is highest for the highest-responding 
# peptide, so the slope, which is fitted as the coefficient of ln, is lowest for the highest-responding peptide
# conversely, if you fit spike ~ nl + 0, the 15N:L ratio is lowest for the highest-responding peptide,
# so its slope is the highest. the above model reflects this, giving a slope of 447 ng/mL for VVEQ and
# only 39 ng/mL for RPKP. thus the slopes reflect how well the peptides "respond" in terms of light peak area
# also consider that you want to fix the intercept at 0, in other words, you want to specify 
# that the model should assume that the 15N:L ratio is 0 when the 15N spike concentration is 0
# if you model based on L:15N ratio, you can't specify this condition.


table_s3_prep = sqldf("
select   h.peptide, 
         cast(round(avg(light)/1000000,2) as text)||'±'||cast(round(stdev(light)/1000000,2) as text) meansd_light,
         cast(round(avg(n15)/1000000,2) as text)||'±'||cast(round(stdev(n15)/1000000,2) as text) meansd_n15,
         cast(round(avg(light/n15),1) as text)||'±'||cast(round(stdev(light/n15),1) as text) meansd_ln,
         avg(raw_ngml) raw_ngml,
         avg(h.rf) respfact,
         avg(norm_ngml) norm_ngml
from     clin c, hupeps h
where    c.peptide = h.peptide
group by 1
order by h.ncorder
;")

# average normalized PrP concentration across all peptides
mean(table_s3_prep$norm_ngml)

table_s3_disp = table_s3_prep
table_s3_disp$raw_ngml = round(table_s3_disp$raw_ngml, 0)
table_s3_disp$respfact = round(table_s3_disp$respfact, 1)
table_s3_disp$norm_ngml = round(table_s3_disp$norm_ngml, 0)

table_s3_disp

write.table(table_s3_disp, 'display_items/script_generated/table-s3.tsv',sep='\t',col.names=T,row.names=F,quote=F)




#### Begin Figure S6 - depiction of normalization scheme
imgsave(paste('display_items/script_generated/figure-s06.',imgmode,sep=''),width=6.5*resx,height=6.5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4),nrow=2,byrow=T)
layout(layout_matrix)

paneladj = -0.1

par(mar=c(4,5,3,1))

spike_lims = c(0.1,2.5)
ratio_lims = c(.001,10)

plot(NA, NA, xlim=ratio_lims, ylim=spike_lims, log='x', ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=r_ats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=r_bigs, labels=r_bigs_labs, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=log10(xlabs$spike), labels=xlabs$spike, las=2)
mtext(side=2, line=2.5, text=expression(''^'15'*'N spike (ng/mL)'))
mtext(side=1, line=2.5, text=expression(''^'15'*'N:light peak area ratio'))
for (i in 1:nrow(hupeps)) {
  points(x=r_ats, y=log10(r_ats*hupeps$slope[i]), type='l', col=hupeps$color[i])
}
points(expt15$nl, jitter(log10(expt15$spike),0.2), pch=20, col=expt15$pep_col)
mtext(side=3, line=0, text='raw response curves', cex=1, font=1)
mtext('A', side=3, cex=2, adj = paneladj, line = 0.5)


expt15$rf = hupeps$rf[match(expt15$peptide, hupeps$peptide)]
plot(NA, NA, xlim=ratio_lims, ylim=spike_lims, log='x', ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=r_ats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=r_bigs, labels=r_bigs_labs, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=log10(xlabs$spike), labels=xlabs$spike, las=2)
mtext(side=2, line=2.5, text=expression(''^'15'*'N spike (ng/mL)'))
mtext(side=1, line=2.5, text=expression(''^'15'*'N:light peak area ratio'))
for (i in 1:nrow(hupeps)) {
  points(x=r_ats, y=log10(r_ats*hupeps$slope[i]*hupeps$rf[i]), type='l', lwd=7-i, col=hupeps$color[i])
}
points(expt15$nl/expt15$rf, jitter(log10(expt15$spike),0.2), pch=20, col=expt15$pep_col)
mtext(side=3, line=0, text='normalized response curves', cex=1, font=1)
mtext('B', side=3, cex=2, adj = paneladj, line = 0.5)

par(mar=c(7,5,3,1))

ylims = c(0,2500)
xlims = c(0.5, 6.5)

yats = (0:5)*500

plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=yats, lwd=1, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=3, text='raw [PrP] (ng/mL)')
for (indiv in unique(clin$indiv_id)) {
  subs = subset(clin, indiv_id==indiv)
  subs = subs[with(subs, order(ncorder)),]
  points(subs$ncorder, subs$raw_ngml, type='l', lwd=0.25, col='#CCCCCC')
  points(subs$ncorder, subs$raw_ngml, pch=20, col=subs$color_pep)
}
mtext(side=3, line=0, text='raw concentration', cex=1, font=1)
mtext(side=1, line=0.5, at=hupeps$ncorder, text=hupeps$peptide, font=2, las=2, col=hupeps$color, cex=0.6)
mtext('C', side=3, cex=2, adj = paneladj, line = 0.5)



plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=yats, lwd=1, lwd.ticks=1, tck=-0.05, las=2)
mtext(side=2, line=3, text='normalized [PrP] (ng/mL)')
for (indiv in unique(clin$indiv_id)) {
  subs = subset(clin, indiv_id==indiv)
  subs = subs[with(subs, order(ncorder)),]
  points(subs$ncorder, subs$norm_ngml, type='l', lwd=0.25, col='#CCCCCC')
  points(subs$ncorder, subs$norm_ngml, pch=20, col=subs$color_pep)
}
mtext(side=3, line=0, text='normalized concentration', cex=1, font=1)
mtext(side=1, line=0.5, at=hupeps$ncorder, text=hupeps$peptide, font=2, las=2, col=hupeps$color, cex=0.6)
mtext('D', side=3, cex=2, adj = paneladj, line = 0.5)


dev.off() #### End Figure S6





### Begin Table S4 - CV and L:15N by quartile

# Darn, the R sqldf implementation of sqlite does not support rank/partition!
# sqldf("
# select   peptide, indiv_id, ln_mean, rank () over (partition by peptide order by ln_mean asc) lnrank
# from     clin
# order by 1, 4
# ;")
# Do this manually instead:
# sort by peptide and L:N ratio
clin = clin[with(clin, order(peptide, ln_mean)),]
clin$ln_rank = rep(1:55, 6)
clin$ln_rank[is.na(clin$light)] = NA # remove the missing VVEQ sample-peptide combination
# wow, that was actually way easier than rank/partition
# now check my work:
for (currpeptide in peptides$peptide[peptides$human=='yes']) {
  subs = subset(clin, peptide==currpeptide)
  cor = cor.test(subs$ln_mean, subs$ln_rank, method='spearman')
  if (cor$estimate < .99999) {
    cat(paste(currpeptide,' has spearman correlation ',cor$estimate,'\n'))
    flush.console()
    break
  }
}
clin$percentile = (clin$ln_rank-1)/(nrow(clin)/6)
clin$quartile = floor(clin$percentile * 4)/4

table_s4_prep = sqldf("
                      select   peptide, quartile, avg(ln_sd/ln_mean) mean_cv, avg(ln_mean) mean_ln
                      from     clin
                      where    ln_mean is not null and ln_sd is not null
                      group by 1, 2
                      order by 1, 2
                      ;")
table_s4_prep

table_s4_left = dcast(table_s4_prep, peptide ~ quartile, value.var='mean_cv')
table_s4_right = dcast(table_s4_prep, peptide ~ quartile, value.var='mean_ln')

for (i in 2:5) {
  table_s4_left[,i] = gsub(' ','',paste(formatC(table_s4_left[,i]*100, digits=1, format='f'),"%",sep="") )
}

for (i in 2:5) {
  table_s4_right[,i] = formatC(table_s4_right[,i], format='f', digits=1)
}

table_s4_joined = sqldf("
                        select   p.ncorder, p.peptide, l.*, r.*
                        from     peptides p, table_s4_left l, table_s4_right r
                        where    p.peptide = l.peptide and p.peptide = r.peptide
                        ;")
table_s4_joined = table_s4_joined[,c(-3, -8)] # delete the extraneous "peptide" cols
table_s4_joined = table_s4_joined[,-1] # delete the ncorder col

colnames(table_s4_joined) = c('peptide','cv0_24','cv25_49','cv50_74','cv75_100','mean0_24','mean25_49','mean50_74','mean75_100')

write.table(table_s4_joined, 'display_items/script_generated/table-s4.tsv', sep='\t', col.names=T, row.names=F, quote=F)
#### -- End Table S4







# Figure 3. All PrP peptides are decreased in the CSF of prion disease patients.
# A. 6 peptides
# B. ELISA

imgsave(paste('display_items/script_generated/figure-3.',imgmode,sep=''),width=6.5*resx,height=2.5*resx,res=resx)

layout_matrix = matrix(c(1,1,1,2),nrow=1,byrow=T)
layout(layout_matrix)

par(mar=c(4,5,5,4))

fig2_alpha = 0.65
ysublab = -.2083

peplabs = sqldf("
                select   peptide, (ncorder)*4 + 2 xmid
                from     peptides
                where    human = 'yes'
                order by 2 asc
                ;")


xlims = expand_range(clin$xall)
ylims = c(60, 3000)
yats = c(10*(1:9), 100*(1:9), 1000*(1:9))
ybigs = c(30, 100, 300, 1000)

par(mar=c(4,5,5,1))
plot(NA, NA, xlim=xlims, ylim=ylims, log='y', xaxs='i', yaxs='i', ann=F, axes=F)

for (xval in unique(clin$xall)) {
  rows = clin$xall == xval & !is.na(clin$norm_ngml)
  clin$xbee[rows] = clin$xall[rows] + beeswarm(log10(norm_ngml) ~ xall, data=clin[rows,], method='square', corral='random', corralWidth=0.5, do.plot=F)$x - 1
}

axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=peplabs$xmid, labels=peplabs$peptide, lwd=0, lwd.ticks=0, cex.axis=0.65, font=2, line=-0.75)
#mtext(side=1, line=3, text='MRM results by peptide')
axis(side=2, at=yats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=ybigs, lwd=0, lwd.ticks=1, tck=-0.05, las=2)
abline(h=1, lwd=1, lty=3)
mtext(side=2, line=3.5, text='MRM CSF [PrP] (ng/mL)')
par(xpd=T)
points(clin$xbee, clin$norm_ngml, pch=19, col=alpha(clin$color_diag,fig2_alpha), lwd=0)
par(xpd=F)

# summarize by diagnosis and by peptide
diag_pep = sqldf("
                 select   m.peptide, m.prion_category, m.xall, avg(norm_ngml) mean, stdev(norm_ngml) sd, count(*) n
                 from     clin m
                 group by 1, 2, 3
                 order by 3 asc
                 ;")
diag_pep$l95 = diag_pep$mean - 1.96*diag_pep$sd/sqrt(diag_pep$n)
diag_pep$u95 = diag_pep$mean + 1.96*diag_pep$sd/sqrt(diag_pep$n)
arrows(x0=diag_pep$xall, y0=diag_pep$l95, y1=diag_pep$u95, code=3, length=0.05, angle=90, col='black')
points(x=diag_pep$xall, y=diag_pep$mean, pch=15, col='black')

par(xpd=T)
legend(x=max(xlims)*.7, y=max(ylims)*5, legend=pparms$prion_category, col=pparms$color, text.col=pparms$color, pch=19, text.font=2, bty='n')
par(xpd=F)

# arrows to indicate N vs. C terminal
nctermcol = '#5C5C5C'
par(xpd=T)
yarr = ysublab * max(ylims)
xprop = .15
arrows(x0=min(xlims), x1=min(xlims) + xprop*(max(xlims)-min(xlims)), y0=yarr, code=1, length=0.1, angle=30, lwd=2, col=nctermcol)
text(x=min(xlims) + xprop*(max(xlims)-min(xlims)), y=yarr, pos=4, labels='N terminus', col=nctermcol, font=2)
arrows(x0=max(xlims), x1=max(xlims) - xprop*(max(xlims)-min(xlims)), y0=yarr, code=1, length=0.1, angle=30, lwd=2, col=nctermcol)
text(x=max(xlims) - xprop*(max(xlims)-min(xlims)), y=yarr, pos=2, labels='C terminus', col=nctermcol, font=2)
par(xpd=F)

mtext('A', side=3, cex=2, adj = -0.1, line = 0.8)



# preparing for plotting of ELISA data
clinsamp$x = pparms$x[match(clinsamp$prion_category, pparms$prion_category)]
clinsamp$color = alpha(pparms$color[match(clinsamp$prion_category, pparms$prion_category)],fig2_alpha)
for (prion_category in pparms$prion_category) {
  rows = clinsamp$prion_category == prion_category
  clinsamp$xbee[rows] = clinsamp$x[rows] + beeswarm(prp_ngml ~ x, data=clinsamp[rows,], do.plot=F)$x - 1
}
csmry = sqldf("
              select   prion_category, x, color, avg(prp_ngml) mean, stdev(prp_ngml) sd, count(*) n
              from     clinsamp
              group by 1, 2, 3
              order by 2 asc
              ;")
csmry$l95 = csmry$mean - 1.96*csmry$sd/sqrt(csmry$n)
csmry$u95 = csmry$mean + 1.96*csmry$sd/sqrt(csmry$n)



xlims = c(0.5, 3.5)
ylims = c(12, 1700)
par(mar=c(4,5,5,3))
plot(NA, NA, xlim=xlims, ylim=ylims, log='y', xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)

#redo beeswarming because size of graphics device matters
for (prion_category in pparms$prion_category) {
  rows = clinsamp$prion_category == prion_category
  clinsamp$xbee[rows] = clinsamp$x[rows] + beeswarm(log10(prp_ngml) ~ x, data=clinsamp[rows,], method='square', corral='random', corralWidth=0.5, do.plot=F)$x - 1
}


points(clinsamp$xbee, clinsamp$prp_ngml, pch=19, col=clinsamp$color, lwd=0)
arrows(x0=csmry$x, y0=csmry$l95, y1=csmry$u95, code=3, length=0.05, angle=90, col='black')
points(csmry$x, csmry$mean, pch=15, col='black')
axis(side=1, at=c(0,4), labels=NA, lwd=1, lwd.ticks=0)
#axis(side=1, at=pparms$x, labels=pparms$prion_category, lwd=0, lwd.ticks=0)
axis(side=1, at=2, labels='BetaPrion Human ELISA', lwd=0, lwd.ticks=0, font=2, cex.axis=.65, line=-0.75)
ylab = ysublab*max(ylims)
par(xpd=T)
text(x=2, y=ylab, labels='full-length PrP')
par(xpd=F)
axis(side=2, at=yats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=ybigs, lwd=0, lwd.ticks=1, tck=-0.05, las=2)
#axis(side=2, at=(0:5)*100, lwd=1, lwd.ticks=1, las=2)
abline(h=10, lwd=2, lty=3, col='red')
# mtext(side=4, at=10, las=2, col='red', line=0.5, text='LLQ')
mtext(side=2, line=3, text='ELISA CSF [PrP] (ng/mL)')
mtext('B', side=3, cex=2, adj = -0.1, line = 0.8)


dev.off()





imgsave(paste('display_items/script_generated/figure-4.',imgmode,sep=''),width=6.5*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(c(1,2,2),nrow=1,byrow=T)
layout(layout_matrix)

# check correlation between all peptides and each other and ELISA

# old color ramp from 0 to 1
# color_ramp = c('#FFFFFF','#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081','#000000')
# corparams = data.frame(floor=(0:10)/10, color=color_ramp)

# color ramp from -1 to 1
color_ramp = c("#230007","#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a", "#fc9272", "#fcbba1", "#fee0d2", "#fff5f0",'#FFFFFF','#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081','#000000')
corparams = data.frame(floor=(-10:10)/10, color=color_ramp)

pep_elisa = peptides[peptides$human=='yes',c('ncorder','peptide')]
for (i in 1:nrow(pep_elisa)) {
  subs = subset(clin, peptide == pep_elisa$peptide[i])
  cor = cor.test(subs$norm_ngml, subs$elisa_ngml, method='spearman', alternative='two.sided')
  pep_elisa$r[i] = signif(cor$estimate,digits=2)
  pep_elisa$p[i] = signif(cor$p.value,digits=2)
}
pep_elisa$y = 7-pep_elisa$ncorder
pep_elisa$x = 1
pep_elisa$floor = floor(pep_elisa$r*10)/10
pep_elisa$col = corparams$color[match(pep_elisa$floor, corparams$floor)]

xlims = c(0.5,1.5)
ylims = c(0.5,6.5)
par(mar=c(10,12,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=pep_elisa$y, labels=pep_elisa$peptide, las=2, lwd=0, lwd.ticks=0)
axis(side=1, at=mean(xlims), labels='ELISA', lwd=0, lwd.ticks=0, las=2)
points(pep_elisa$x, pep_elisa$y, pch=15, col=pep_elisa$col, cex=5)
# text(pep_elisa$x, pep_elisa$y, labels=formatC(pep_elisa$r,format='f',digits=2), cex=0.8, font=2, col='white')

mtext('A', side=3, cex=2, adj = 0.0, line = 0.5)

crosspeps = sqldf("
                  select   p1.ncorder x, 7-p2.ncorder y, p1.peptide pep1, p2.peptide pep2
                  from     peptides p1, peptides p2
                  where    p1.human = 'yes' and p2.human = 'yes'
                  and      p1.ncorder <= p2.ncorder
                  order by p1.ncorder asc, p2.ncorder asc
                  ;")
crosspeps$r = 0.0
crosspeps$p = 0.0
for (i in 1:nrow(crosspeps)) {
  c1 = clin[clin$peptide==crosspeps$pep1[i],c('indiv_id','norm_ngml')]
  c2 = clin[clin$peptide==crosspeps$pep2[i],c('indiv_id','norm_ngml')]
  colnames(c1) = c('indiv_id','pep1_ngml')
  c1$pep2_ngml = c2$norm_ngml[match(c1$indiv_id, c2$indiv_id)]
  cor = cor.test(c1$pep1_ngml, c1$pep2_ngml, method='spearman', alternative='two.sided')
  crosspeps$r[i] = signif(cor$estimate,digits=2)
  crosspeps$p[i] = signif(cor$p.value,digits=2)
}

crosspeps$floor = floor(crosspeps$r*10)/10
crosspeps$col = corparams$color[match(crosspeps$floor, corparams$floor)]

par(mar=c(10,2,3,1))
lims = c(0.5, 6.5)
plot(NA, NA, xlim=c(min(lims),9.5), ylim=lims, pch=15, col=crosspeps$col, cex=5, axes=F, ann=F, xaxs='i', yaxs='i')
#axis(side=2, at=7-peptides$ncorder[peptides$human=='yes'], labels=peptides$peptide[peptides$human=='yes'], las=2, lwd=0, lwd.ticks=0)
axis(side=1, at=peptides$ncorder[peptides$human=='yes'], labels=peptides$peptide[peptides$human=='yes'], las=2, lwd=0, lwd.ticks=0)
points(crosspeps$x, crosspeps$y, pch=15, col=crosspeps$col, cex=5)
# text(crosspeps$x, crosspeps$y, labels=formatC(crosspeps$r,format='f',digits=2), cex=0.8, font=2, col='white')

# draw rectangles for inter- and intra-domain
# N-N
rect(xleft=.6,xright=1.4,ybottom=4.6,ytop=5.4,col=NA,lwd=2)
segments(x0=1.5, y0=5.5, x1=2, y1=6, lwd=1.5)
text(x=2,y=6, pos=4,labels='N-N intradomain')
# N-C
rect(xleft=.6,xright=2.4,ybottom=.6,ytop=4.4,col=NA,lwd=2)
segments(x0=2.5, y0=4.5, x1=3, y1=5, lwd=1.5)
text(x=3,y=5, pos=4,labels='interdomain')
# C-C
polygon(x=c(2.6,2.6,3.4,3.4,4.4,4.4,5.4,5.4,2.6), y=c(.6,3.4,3.4,2.4,2.4,1.4,1.4,0.6,0.6),lwd=2)
segments(x0=3.5, y0=3.5, x1=4, y1=4, lwd=1.5)
text(x=4,y=4, pos=4,labels='C-C intradomain')

scale_xleft = 6.2
scale_xright = 8.8
scale_factor = (scale_xright - scale_xleft) / 21
scale_ybot = 5.1
scale_ytop = 5.4
par(xpd=T)
rect(xleft=scale_xleft,xright=scale_xright,ybottom=scale_ybot,ytop=scale_ytop,lwd=2.5)
rect(xleft=scale_xleft+(0:20)*scale_factor, xright=scale_xleft+(1:21)*scale_factor, ybottom=rep(scale_ybot,5), ytop=rep(scale_ytop,5), border=NA, col=corparams$color)
text(x=(scale_xright + scale_xleft)/2, y=scale_ytop, pos=3, labels="Spearman's correlation", cex=1.0)
text(x=c(scale_xleft,mean(c(scale_xleft,scale_xright)),scale_xright),y=c(scale_ybot,scale_ybot,scale_ybot),pos=1,labels=c('-1','0','1'),cex=1.0,font=3)
par(xpd=F)

mtext('B', side=3, cex=2, adj = 0.0, line = 0.5)

dev.off() # -- end Figure 4






# Stats for results section:

# Create a composite ng/mL call based on all (or maybe just Steve's favorite 4?) peptides
# lm(mrm_ngml ~ total_protein)
# lm(elisa_ngml ~ mrm_ngml)
# lm(elisa_ngml ~ mrm_ngml + total_protein)
# lm(mrm_ngml ~ hemoglobin)
# t.test(0 transfers chaps vs. 0 transfers neat), or something similar, from expt13


clin_smry = sqldf("
select   indiv_id, avg(norm_ngml) mrm_prp, avg(elisa_ngml) elisa_prp, avg(dc_mgml) total_protein, avg(hb_ngml) hemoglobin
from     clin
group by 1
order by 1
;")
mean(clin_smry$mrm_prp)
mean(clin_smry$elisa_prp)

cor.test(clin_smry$mrm_prp, clin_smry$elisa_prp, method='spearman', alternative='two.sided')

cor.test(clin_smry$mrm_prp, clin_smry$hemoglobin, method='spearman', alternative='two.sided')

m = lm(mrm_prp ~ hemoglobin, data=clin_smry)
summary(m)

m = lm(mrm_prp ~ total_protein, data=clin_smry)
summary(m)

m = lm(elisa_prp ~ mrm_prp, data=clin_smry)
summary(m)

m = lm(elisa_prp ~ total_protein, data=clin_smry)
summary(m)


m = lm(elisa_prp ~ mrm_prp + total_protein, data=clin_smry)
summary(m)



deterg = subset(mrm_data, run == 13 & peptide %in% peptides$peptide[peptides$human=='yes'])
deterg$chaps = ''
deterg$chaps[grepl('^C',deterg$id)] = '0.03% CHAPS'
deterg$chaps[grepl('^N',deterg$id)] = 'neat'
deterg$transfers = as.integer(gsub('[AB]','',gsub('[CN]','',deterg$id)))
deterg$replicate = substr(deterg$id,3,3)

deterg_baseline = sqldf("
                        select   peptide, avg(lh) baseline_lh, count(*) n
                        from     deterg
                        where    chaps = '0.03% CHAPS'
                        and      transfers = 0
                        and      peptide in (select peptide from peptides where human = 'yes')
                        group by 1
                        order by 1
                        ;")
deterg_baseline 

deterg$baseline = deterg_baseline$baseline_lh[match(deterg$peptide, deterg_baseline$peptide)]

deterg$rel = deterg$lh / deterg$baseline


m = lm(rel ~ chaps, data=subset(deterg, transfers==0))
summary(m)
percent(1/(1+coefficients(m)[2]) - 1) # % by which recovery is increased in presence of chaps
summary(m)$coefficients['chapsneat','Pr(>|t|)'] # p value for this difference. 
# double check that lm() in R with a categorical variable is just equivalent to ANOVA
anova = aov(rel ~ chaps, data=subset(deterg, transfers==0))
summary(anova)[[1]][["Pr(>F)"]]







# finally, prepare a table with individual clinical sample data for GitHub
processed_data_table = sqldf("
select   peptide, indiv_id, prion_category diagnostic_category, light, n15, ln_mean ln_ratio, norm_ngml mrm_ngml
from     clin
order by ncorder, indiv_id
;")

write.table(processed_data_table, 'data/summary/clinical_sample_mrm_results.tsv',sep='\t',na = '',row.names=F,col.names=T,quote=F)
