library(Biostrings)
library(dplyr)

# Input sequence files are from database searches of Los Alamos HIV-1 (all subtypes)
# sequence databases where the ranges in the file names were used to retrieve sequences
# where the ranges correspond to Los Alamos alignments to HXB2. 

# U5 HMM alignment
#-------------------------------------------------------------------------------

o <- readDNAStringSet('search_allSub_9569-9719.fasta')
r <- as.character(reverseComplement(o))

# Exclude sequences that do not have CA near their ends (within 5 NT of end). 
i <- stringr::str_locate(r, 'TG')
i <- i[,1]-1
z <- i > 5
o <- o[!z]
i <- i[!z]

# Los Alamos searches were for 150 NT. Limit sequences to 100 NT.
u5 <- unique(subseq(o, width(o)-100-i, width(o)-i))

# Use CD-HIT-EST to cluster sequences and retrieve sequences in the top cluster
# so that we do not use sequences that are too dissimilar to build the MSA for the HMM.
s <- unlist(lapply(strsplit(names(u5), '\\.'), '[[', 1))
d <- tibble(id = gsub('\\.' , ';', names(u5)), subtype = s, seq = as.character(u5))
x <- DNAStringSet(d$seq)
names(x) <- d$id

# Cluster at 80% seq id. CD-HIT-EST will not go lower than 80%.
writeXStringSet(x, 'u5.fasta')
system('~/software/cd-hit_4.8.1/cd-hit-est -i u5.fasta -o u5.cdhit -c 0.8 -U 10 -gap -6 -sc 1 -d 0')

# Extract read ids from largest cluster.
f <- readLines('u5.cdhit.clstr')
f <- f[2:which(grepl('>Cluster 1', f))[1]-1]
x <- x[names(x) %in% sub('>', '', unlist(stringr::str_extract(f, '>([^\\.]+)')))]

# Exclude reads that cause alignments to go awry (based on inspection of MSA without exclusions).
# https://alignmentviewer.org

exclude <- c("A1D;UG;2007;191955;MW006055",
             "B;US;2019;5120-WK26-PL12-D21;OM205091",
             "B;US;2019;5111-WK26-PL7-D14;OM204179",
             "B;US;2019;5111-WK26-PL13-I10;OM204010",
             "C;IN;2000;Hy49;JQ268702")

x <- x[! names(x) %in% exclude]

# Build MSA.
writeXStringSet(x, 'u5_select.fasta')
system('~/software/mafft/bin/mafft --maxiterate 1000 --globalpair u5_select.fasta > u5_select.mafft')

# Make sure all alignments end on CA.
z <- readDNAStringSet('u5_select.mafft')
table(subseq(z, width(z)-1, width(z)))
writeXStringSet(z, 'U5.aln')


# U3 HMM alignment
#-------------------------------------------------------------------------------

o <- readDNAStringSet('search_allSub_1-150.fasta')
i <- stringr::str_locate(as.character(o), 'TG')
i <- i[,1]-1
z <- i > 5
o <- o[!z]
i <- i[!z]

# Los Alamos searches were for 150 NT. Limit sequences to 100 NT.
u3 <- unique(subseq(o, 1+i, 100+i))
s <- unlist(lapply(strsplit(names(u3), '\\.'), '[[', 1))
d <- tibble(id = gsub('\\.' , ';', names(u3)), subtype = s, seq = as.character(u3))
x <- DNAStringSet(d$seq)
names(x) <- d$id

# Cluster at 80% seq id.
writeXStringSet(x, 'u3.fasta')
system('~/software/cd-hit_4.8.1/cd-hit-est -i u3.fasta -o u3.cdhit -c 0.8 -U 10 -gap -6 -sc 1 -d 0')

f <- readLines('u3.cdhit.clstr')
f <- f[2:which(grepl('>Cluster 1', f))[1]-1]
x <- x[names(x) %in% sub('>', '', unlist(stringr::str_extract(f, '>([^\\.]+)')))]


# Exclude reads that cause alignments to go awry (based on inspection of MSA without exclusions).
exclude <- c("01_AE;CN;2006;FJ054;DQ859180",
             "B;SE;1996;18_sw_a;AF063179",
             "B;SE;1996;18_sw_c;AF063180",
             "D;UG;1991;UG270;AB485651",
             "02_AG;JP;-;YM0178;AB195666")

x <- x[! names(x) %in% exclude]

# Build MSA.
writeXStringSet(x, 'u3_select.fasta')
system('~/software/mafft/bin/mafft --maxiterate 1000 --globalpair u3_select.fasta > u3_select.mafft')

# Make sure all alignments start with TG.
z <- readDNAStringSet('u3_select.mafft')
table(subseq(z, 1, 3))

writeXStringSet(z, 'U3.aln')
writeXStringSet(reverseComplement(z), 'U3_rc.aln')


# Build HMMs.
#-------------------------------------------------------------------------------

system('~/software/hmmer-3.3.2/bin/hmmbuild HIV1_1-100_U5.hmm U5.aln')
system('~/software/hmmer-3.3.2/bin/hmmbuild HIV1_1-100_U3.hmm U3.aln')
system('~/software/hmmer-3.3.2/bin/hmmbuild HIV1_1-100_U3_RC.hmm U3_rc.aln')
