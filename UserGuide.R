#Example
library(meripDeep)

GENE_ANNO_GTF = system.file("extdata", "example.gtf", package="exomePeak2")

f1 = system.file("extdata", "IP1.bam", package="exomePeak2")
f2 = system.file("extdata", "IP2.bam", package="exomePeak2")
f3 = system.file("extdata", "IP3.bam", package="exomePeak2")
f4 = system.file("extdata", "IP4.bam", package="exomePeak2")
IP_BAM = c(f1,f2,f3,f4)
f1 = system.file("extdata", "Input1.bam", package="exomePeak2")
f2 = system.file("extdata", "Input2.bam", package="exomePeak2")
f3 = system.file("extdata", "Input3.bam", package="exomePeak2")
INPUT_BAM = c(f1,f2,f3)

f1 = system.file("extdata", "treated_IP1.bam", package="exomePeak2")
TREATED_IP_BAM = c(f1)
f1 = system.file("extdata", "treated_Input1.bam", package="exomePeak2")
TREATED_INPUT_BAM = c(f1)

txdb = makeTxDbFromGFF(GENE_ANNO_GTF)

#Peak Calling with GFF
library(BSgenome.Hsapiens.UCSC.hg19)
Peaks <- peakCalling(bam_IP = IP_BAM,
                     bam_input = INPUT_BAM,
                     txdb = txdb,
                     genome = BSgenome.Hsapiens.UCSC.hg19)

#Peak Calling with TxDb package
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
Peaks <- peakCalling(bam_IP = IP_BAM,
                     bam_input = INPUT_BAM,
                     txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                     genome = BSgenome.Hsapiens.UCSC.hg19)

#Access to metadata columns
mcols(Peaks)

#Plot transcript topology
#list_grl should be a list of GRanges / GRangesList
#names of the list will be the label
plotTopology(list_grl = Peaks, txdb = txdb, savePrefix = "topology")
