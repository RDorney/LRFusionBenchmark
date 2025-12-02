#run featurecounts
library(Rsubread)
library(tidyr)
bamfiles <- list.files("/bioinformatics/ryley/Library_Benchmark", pattern = "*.sorted.bam$", full.names = TRUE) 
transcriptomebam <- grepl(pattern = "transcriptome", bamfiles)
bamfiles <- bamfiles[!transcriptomebam]

featureCountsAnalysis <- featureCounts(files=bamfiles,
              annot.inbuilt="hg38",
              useMetaFeatures=TRUE,
              allowMultiOverlap = TRUE,
              countMultiMappingReads=TRUE,
              genome="/opt/JAFFA-version-2.3/hg38.fa",
              countChimericFragments=TRUE,
              nthreads=4)

#featureCounts(files=bamfiles,
#              annot.ext=~/reference_files/gencode.v47.annotation.gtf,
#              isGTFAnnotationFile=TRUE,
#              useMetaFeatures=TRUE,
#              allowMultiOverlap = TRUE,
#              countMultiMappingReads=TRUE,
#              genome="/opt/JAFFA-version-2.3/hg38.fa",
#              countChimericFragments=TRUE,
#              nthreads=6)
