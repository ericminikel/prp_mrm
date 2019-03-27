options(stringsAsFactors=FALSE)
setwd('~/d/sci/src/mrm_prp/')

# read in list of samples from Parchi & Zerr cohorts, with data from previous ELISA study
samples = read.table('data/samples/samples.tsv',sep='\t',header=T)

# set random seed so this script produces the *same* random result every time
set.seed(1)

# each batch can be 11 samples (in duplicate = 22 replicates, + 2 replicats of interplate control = full 24)
batchno = c(rep(1:5,each=11),rep(6,each=8))
samples$batchno = sample(batchno,replace=FALSE)

# confirm that the two cohorts & three diagnostic categories are distributed throughout
table(samples[,c('batchno','sample_source')])
table(samples[,c('batchno','prion_category')])

# write out a minimal table sorted by batchno to put into plans for each lab day
batches = samples[with(samples,order(batchno, sample)), c('batchno','sample')]

write.table(batches,'data/samples/batches.tsv',sep='\t',quote=F,row.names=F,col.names=T)

# name duplicates of each, A1-X1, etc. and insert an interplate control at a random
# order within each batch
deids = data.frame(batch=c(rep(1:5,each=24),rep(6,each=18)),
                   label=as.character(NA),
                   sample=as.character(NA))
deids$label = paste(rep(LETTERS[1:24],6)[1:138], deids$batch, sep='')
ipc_order = c(sample(1:12, size=5, replace=TRUE), sample(1:9, size=1, replace=TRUE))
ipc_indices = rep(ipc_order*2, each=2) + rep(c(-1,0),6)
ipc_absolute_indices = (rep(1:6,each=2)-1)*24 + ipc_indices
deids$sample[ipc_absolute_indices] = '9-IPC'
# now insert the sample IDs back in in order, around the IPCs
j = 1
for (i in 1:nrow(batches)) {
  while (!is.na(deids$sample[j])) {
    j = j + 2
  } 
  deids$sample[j:(j+1)] = batches$sample[i]
}

write.table(deids, 'data/samples/deids.tsv',sep='\t',quote=F,row.names=F,col.names=T)

# write out individual metadata files
last_expt = 15
for (batchno in 1:6) {
  expt = last_expt + batchno
  meta = data.frame(prep_date='TBD',
                    label=deids$label[deids$batch==batchno],
                    skyline_id='TBD',
                    sample=deids$sample[deids$batch==batchno],
                    description='')
  write.table(meta,paste('data/meta/expt',expt,'.tsv',sep=''),sep='\t',quote=F,row.names=F,col.names=T)
}

# note for anyone reading this script trying to reconcile the 6 batches here vs. 5 in data/
# we originally had N=63 samples divided into batches of 11, 11, 11, 11, 11, and 8. the 
# final batch (what should have been expt21) failed - no peptides were recovered (light nor heavy nor 15N),
# apparently due to a desalt solvent issue. we didn't have any additional aliquots of the final 8 samples
# suitable for repeating the experiment, so we proceeded to analyze just the 55 that had been run successfully
