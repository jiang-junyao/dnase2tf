# dnase2tf
It's hard to install original version of dnase2tf in windows system, so i put dnase2tf on github to install it handily by: 
```
devtools::install_github('jiang-junyao/dnase2tf')
```
For more detail of this package, please see: https://sourceforge.net/projects/dnase2tfr/

Parameter of dnase2tf:
# dnase2tf(datafilepath, hotspotfilename, mapfiledir, outputfilepath, ...)
# 	finds footprint candidates from the DNase-I seq data (given by
# 	datafilepath) on the user specified regions of interest (given by
# 	hotspotfilename).  Mappability file location (mapfiledir) and the
# 	filename for the output file are also required.
# 
# Arguments
# 
# Datafilepath: the directory path+ the file name prefix of the DNaseI
# 	sequence reads.  The input files are split corresponding to each
# 	chromosome and shares the common "prefixes" and 'suffixes' in their file
# 	names.   Those names differ between chromosomes.    These files are
# 	assumed to be in the format 'DATA NAME(Prefix)' + '_chr#' + '.txt'.   For
# 	example, typical data files would be named as  "DHS_Dex_chr1.txt?,
# 	'DHS_Dex_chr2.txt', etc and located in the same folder.  In the data
# 	files, each row represents the genomic location of a mapped DNA sequence
# 	(chromosome, start, stop, strand).   Tabs separate fields and coordinates
# 	are 1-based.
# 
# Hotspotfilename:  the directory path of the region of interest in BED
# 	file format (0-based).
# 
# Mapfiledir: the directory path of the mappability file of the reference
# 	genome.   The program assumes that mappability files named by chromosomes
# 	as 'chr1b.txt', 'chr2b.txt', etc.
# 
# Outputfilepath: the directory path + the file name prefix of the output.
# 	This program outputs in both BED and BEDGraph file formats with various
# 	FDR (False Discovery Ratio) starting from 0% to 100%.
# 
# Options
# 
# assemseqdir : the directory path of the assembly sequence files (chr1.fa, chr2.fa, etc..)
# 	which can be downloaded from the following links:
# 		http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz  (hg19)
# 		http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz  (mm9)
# 	If omitted, dinucleotide bias correction is disabled in the computation.
# 
#  
# dftfilename: the dinucleotide bias percentage (DFT) file.  A DFT file is
# 	a tab-deliminated text file containing expected and observed proportions
# 	dinucleotide DNAs measure from the data.  The source code of the DFT file
# 	generation program is also provided with this package and the user may
# 	compile and run the program.  This option is ignored if dinucleotide bias
# 	correction is not used.
# 
# numworker: the number of CPU cores to run concurrently. 
# 	To prevent the problem caused by out of memory error, try
# 	the small number first.  By default, 2 workers are assigned for the
# 	computation.
#   
# maxw: The maximum width of footprint candidates.  '30' by default.
# 
# minw: The minimum width of footprint candidates.  '6' by default.
#
# z_threshold: z-score threshold of initial footprint candidate selection.
# '-2' is set to a default value.
