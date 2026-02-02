suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(genomation))

option_list = list(
  make_option('--meta', action='store', type='character',
              help='Filename for metadata defining experiment'),
  make_option('--windowed', action='store_true', default=FALSE,
              help='Look at 1kb windows if selected'),
  make_option('--covariates', action='store_true', default='FALSE', 
              help='Use one of the preconfigured covariates based on base'),
  make_option('--base', action='store', type='character',
              help='Base name for generated output files')
)

args <- commandArgs(TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args)

## set base for file names
if ( opt$windowed ) {
    base <- paste(opt$base, 'windowed', sep='-')
} else {
    base <- opt$base
}

## are we using covariates
if ( opt$covariates ) {
   base <- paste(base, 'with-covariates', sep='-')
}

meta <- read.table( file = opt$meta, header=TRUE, sep="\t", stringsAsFactors=FALSE,)

## read the files to a methylRawList object: myobj
preobj=methRead(as.list(meta$file),
                sample.id=as.list(meta$sample),
                assembly="GRCm39",
                treatment=meta$experiment,
                context="CpG",
                mincov = 15
                )

# windowed if they asked for it
if ( opt$windowed ) {
    myobj <-tileMethylCounts(preobj,win.size=1000,step.size=1000)
} else {
    myobj <- preobj
}

meth=unite(myobj, destrand=TRUE)

## print QC pdf
pdf(paste(base, 'qc.pdf', sep='-'), width=8.5, height=11)

walk(myobj, function(x) getMethylationStats(x,plot=TRUE,both.strands=F))
walk(myobj, function(x) getCoverageStats(x,plot=TRUE,both.strands=F))

clusterSamples(meth, dist="correlation", method="ward.D", plot=TRUE)


getCorrelation(meth,plot=TRUE)

PCASamples(meth)

dev.off()

## We may have covariates
if ( opt$covariates ) {

    if ( opt$base == 'HDM-v-PBS' ) {

        covariates = data.frame(sex=as.factor(meta$sex),
                                genotype=as.factor(meta$genotype))
        
    } else if ( opt$base == 'ESR1-v-WT' || opt$base == 'ESR2-v-WT' ) {
        
        covariates = data.frame(sex=as.factor(meta$sex),
                                treatment=as.factor(meta$treatment))

    } else if ( opt$base == 'ESR1_Male-v-ESR2_Male' || opt$base == 'ESR1_Female-v-ESR2_Female' ) {
        covariates = data.frame(treatment=as.factor(meta$treatment))

    } else if ( opt$base == 'ESR1_Male-v-ESR1_Female' || opt$base == 'ESR2_Male-v-ESR2_Female' ) {
        covariates = data.frame(treatment=as.factor(meta$treatment))

    } else {
        stop("Could not find covariates for this base")
    }
    
    myDiff=calculateDiffMeth(meth, covariates=covariates)
    
} else {
    myDiff=calculateDiffMeth(meth)
}



myDiff10p=getMethylDiff(myDiff,difference=10,qvalue=0.05)

write.table(myDiff10p, paste(base, 'diff-met-out.tsv', sep='-'),
            quote=F, sep='\t', na='NA', col.names=NA)

write.table(myDiff,paste(base, 'all-out.tsv', sep='-'), 
            quote=F, sep='\t', na='NA', col.names=NA)

gene.obj=readTranscriptFeatures('GRCm39-features.bed')

annot <- annotateWithGeneParts(as(myDiff10p,"GRanges"),gene.obj, intersect.chr=T)

write.table(annot@dist.to.TSS, paste(base, 'annotation.tsv', sep='-'), 
            quote=F, sep='\t', na='NA', col.names=NA)
write.table(annot@members, paste(base, 'members.tsv', sep='-'), 
            quote=F, sep='\t', na='NA', col.names=NA)

## For windowed, write the by-sample map
if ( opt$windowed ) {
  write.table(meth,paste(base, 'by-sample.tsv', sep='-'), quote=F, sep='\t', na='NA')
} 

save.image(paste(base, '.RData', sep=''))





