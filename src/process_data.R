setwd('~/d/sci/src/prp_mrm')
options(stringsAsFactors=FALSE)
library(sqldf)
library(reshape2)

peptides = read.table('data/ref/peptides.tsv', sep='\t', header=T, quote='', comment.char='')
metox = read.table('data/ref/metox.tsv', sep='\t', header=T, quote='', comment.char='')

runs = read.table('data/meta/runs.tsv', sep='\t', header=T)

for (i in 1:nrow(runs)) {
  
  if (runs$datafile[i] == '') {
    next # this allows me to have runs.tsv contain expts I've done but haven't yet received skyline files for
  }
  
  meta = read.table(paste('data/meta/',runs$metafile[i],sep=''), sep='\t', header=T, quote='', comment.char='')
  
  rawdata = read.table(paste('data/skyline/',runs$datafile[i],sep=''), sep=',', header=T, quote='', comment.char='')
  colnames(rawdata) = gsub('[^a-z0-9_]+','_',tolower(colnames(rawdata)))
  
  # subset to useful columns and rename them
  # if we have 15N in this experiment
  if (any(grepl('x15n', colnames(rawdata)))) {
    subdata = rawdata[,c('replicate_name','peptide_sequence', 'peptide_modified_sequence', 'peptide_retention_time', 'fragment_ion', 'light_area', 'heavy_area','x15n_area')]
    colnames(subdata) = c('sample', 'peptide', 'modified_sequence', 'retention_time', 'ion', 'light_area', 'heavy_area','n15_area')
  } else {
    subdata = rawdata[,c('replicate_name','peptide_sequence', 'peptide_modified_sequence', 'peptide_retention_time', 'fragment_ion', 'light_area', 'heavy_area')]
    colnames(subdata) = c('sample', 'peptide', 'modified_sequence', 'retention_time', 'ion', 'light_area', 'heavy_area')
    subdata$n15_area = 0
  }
  
  subdata$heavy_area = as.numeric(subdata$heavy_area)
  subdata$n15_area = as.numeric(subdata$n15_area)
  subdata$light_area = as.numeric(subdata$light_area)
  
  # for possible met-ox peptides, choose the more abundant one and flag if it is met-ox
  # only do this if all three versions - heavy, light, and 15N are more abundant for met-ox
  # note use of >= for n15_area comparison, in case this is later applied to cases where all are 0
  metox_compare = sqldf("
  select   red.sample, red.peptide, red.heavy_area red_heavy_area, ox.heavy_area ox_heavy_area, red.n15_area red_n15_area, ox.n15_area ox_n15_area, red.light_area red_light_area, ox.light_area ox_light_area,
           case when ox.light_area > red.light_area and ox.heavy_area > red.heavy_area and ox.n15_area >= red.n15_area then 'met-ox' else '' end flag
           -- red.n15_area red_n15_area, ox.n15_area ox_n15_area
  from     subdata red, subdata ox
  where    red.sample = ox.sample
  and      red.peptide = ox.peptide
  and      red.ion = ox.ion
  and      red.ion in (select best_transition from metox)
  and      red.modified_sequence in (select modified_sequence from metox where flag = '')
  and      ox.modified_sequence in (select modified_sequence from metox where flag = 'met-ox')
  order by 1
  ;")
  
  # subset to best transitions. for now, keep dups for peptides with a met-ox form
  data = sqldf("
  select   s.*
  from     subdata s, peptides p
  where    s.peptide = p.peptide
  and      s.ion = p.best_transition
  ;")
  
  # now go through and choose the reduced or oxidized form as appropriate
  data = sqldf("
  select   d.*, c.flag flag, m.flag required_flag
  from     data d left outer join metox_compare c
  on       d.sample = c.sample and d.peptide = c.peptide
  left outer join metox m
  on       d.peptide = m.peptide and d.modified_sequence = m.modified_sequence
  ;")
  rows_to_delete = which(data$flag != data$required_flag)
  if (length(rows_to_delete) > 0) {
    data = data[-rows_to_delete,] # remove those rows
  }
  data = data[,-which(colnames(data)=='required_flag')] # remove the duplicate flag column
  data$flag[is.na(data$flag)] = ''
  
  data$id = meta$label[match(data$sample, meta$skyline_id)]
  data = data[!is.na(data$id),] # remove things not in the metadata file

  # data$lhr = data$light_area / data$heavy_area
  # # we don't really care about >3 digits to the right of the decimal, so round them off
  # data$lhr = round(data$lhr, digits=5)
  # 
  # data$n15_area[data$n15_area == '#N/A'] = NA
  # data$n15_area = as.numeric(data$n15_area)
  # data$n15hr = data$n15_area / data$heavy_area
  # data$n15hr = round(data$n15hr, digits=5)
  
  roundmean = function(x) round(mean(x), digits=5)
  
  #widetable = dcast(data[,c('peptide','id','lhr')], peptide ~ id, value.var = 'lhr', fun.aggregate = roundmean)
  #write.table(format(widetable, ndigits=3, nsmall=3, scientific=F), paste('data/processed/', runs$metafile[i], sep=''), sep='\t', row.names=F, col.names=T, quote=F)
  
  write.table(data, paste('data/processed/', runs$metafile[i], sep=''), sep='\t', row.names=F, col.names=T, quote=F)
}

