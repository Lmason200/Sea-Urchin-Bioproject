The first things to do once your Jupyter environment is open is to connect to the data store:
Click on the Discovery Environment icon on the left (blue swooshy thing). The parameters for the IRODS Setup Config should be mostly filled out. You just need to enter your cyverse user name and password and hit submit. A directory listing of your files in the Data Store should appear below.
Next, we have to initialize. Click on Terminal in the launcher to open a terminal window.
Type iinit
Enter your password again (Note: no text will appear as you type) and press Enter.

> import qiime2 as q2

##Trim Primers
(May need to be entered as one text line)

! qiime cutadapt trim-paired --i-demultiplexed-sequences paired-end-demux.qza --o-trimmed-sequences paired-end-trimmed-demux.qza --p-cores 16 --p-front-f AATGATACGGCGACCACCGAGATCTACACTATGGTAATTGTGTGCCAGCMGCCGCGGTAA --p-front-r CAAGAGAAGACGGCATACGAGATNNNNNNAGTCAGTCAGCCGGACTACHVGGGTWTCTAAT --p-minimum-length 215 --p-discard-untrimmed --verbose 

## Create Visualization File

! qiime demux summarize \
--i-data trimmed-sequences paired-end-trimmed-demux.qza
--o-visualization trimmed-sequences paired-end-trimmed-demux.qzv

## Load Visualization File 

> q2.Visualization.load ("paired-end-trimmed-demux.qzv")

## Using DADA 2 to Clean Up Samples

! qiime dada2 denoise-single \
--i-demultiplexed-seqs paired-end-trimmed-demux.qza \
--p-trunc-len 180 \
--p-n-threads 0 \
--output-dir work/DADA2_denoising_output \
--verbose
#### Reverse reads were poor quality so only foward were run

## Visualizing Results

####DADA2 STATS

! qiime metadata tabulate \
--m-input-file work/DADA2_denoising_output/denoising_stats.qza \
--o-visualization work/DADA2_denoising_output/denoising_stats.qzv

> q2.Visualization.load("work/DADA2_denoising_output/denoising_stats.qzv

####SEQUENCES

! qiime feature-table tabulate-seqs \
--i-data work/DADA2_denoising_output/representative_sequences.qza \
--o-visualization work/DADA2_denoising_output/rep_seqs.qzv 

> q2.Visualization.load("work/DADA2_denoising_output/rep_seqs.qzv")


####FREQUENCY

! qiime feature-table summarize \
--i-table work/DADA2_denoising_output/table.qza \
--o-visualization work/DADA2_denoising_output/table.qzv 

> q2.Visualization.load("work/DADA2_denoising_output/table.qzv")


## Assigning Taxonomy
#### Download the SILVA classifier
! wget https://data.qiime2.org/2019.10/common/silva-132-99-515-806-nb-classifier.qza
! mv silva-132-99-515-806-nb-classifier.qza work/

#### Classify the Sequences

! qiime feature-classifier classify-sklearn \
--i-classifier work/silva-132-99-515-806-nb-classifier.qza \
--i-reads work/DADA2_denoising_output/representative_sequences.qza \
--output-dir work/classified_sequences \
--verbose

#### Visualize the Classified Sequences

! qiime metadata tabulate \
--m-input-file work/classified_sequences/classification.qza \
--o-visualization work/classified_sequences/taxonomy.qzv
  
> q2.Visualization.load("work/classified_sequences/taxonomy.qzv")

## Making Phylogenetic Tree

##### Make new folder for phylogeny 

> ! mkdir work/phylogeny

! qiime alignment mafft \
--i-sequences work/DADA2_denoising_output/representative_sequences.qza \
--o-alignment work/phylogeny/aligned-rep-seqs.qza

! qiime alignment mask \
--i-alignment work/phylogeny/aligned-rep-seqs.qza \
--o-masked-alignment work/phylogeny/masked-aligned-rep-seqs.qza
> > 
! qiime phylogeny fasttree \
--i-alignment work/phylogeny/masked-aligned-rep-seqs.qza \
--o-tree work/phylogeny/fasttree-tree.qza

! qiime phylogeny midpoint-root \
--i-tree work/phylogeny/fasttree-tree.qza \
--o-rooted-tree work/phylogeny/fasttree-tree-rooted.qza


## Exporting The Files

> ! mkdir export

> ! qiime tools export \
>   --input-path work/DADA2_denoising_output/table.qza \
>   --output-path work/export/table

#### Convert the BIOM Files to tsv to Be Used in R

> ! biom convert \
> -i work/export/table/feature-table.biom \
> -o work/export/table/table.tsv --to-tsv

#### Exporting the Rep Sequence files (fasta)

! qiime tools export \
  --input-path work/DADA2_denoising_output/representative_sequences.qza \
  --output-path work/export/rep-seqs.fasta
  
 ! qiime tools export \
  --input-path work/classified_sequences/taxonomy.qzv \
  --output-path export/taxonomy/TSVTaxonomyDirectoryFormat
  
  ! qiime tools export \
  --input-path work/phylogeny/fasttree-tree-rooted.qza \
  --output-path work/export/exported-tree
  
! iput -r work/export/








###Summary

After exporting the files, the next step was to make them into interpretable data. We did this in R-Studio. In R-Studio, the data was made into tables and certain factors of interest were compared based on the variations in data. 

Time Zero- Time 26 
Time 26- Time 30
Gut- Feces
Feces- Seawater


For the Time Zero- Time 26 data, the presence of large amounts of Gammaproteobacteria and Betaproteobacteria is noteworthy, as well as Rhodospirilla and the expected Vibrionaceae at Time 26. At Time Zero, there is more diversity.

For the Time 26- Time 30 data, there was a lot less Gammaproteobacteria and more diversity such as Flavobacteria, Rhizobiales and Kiritimatella at Time 30. 

For the comparison between Gut and Feces, there is a great deal of Gammaproteobacteria in the Feces, in particular, Cellvibrionales and Caulobacterales. The Gut had the majority Bacteroidia with also some Clostridia. 

For the comparison between Feces and Seawater, there are numbers of Gammaproteobacteria and Alphaproteobacteria as well as Clostridia and Bacilli in the Seawater. However in the Feces, it is a majority Alphaproteobacteria, Gammaproteobacteria, and Bacteroidia. 


Overall, there was an generally an increase in Alpha and Gammaproteobacteria when the temperature increased. Between the tissue types there is an huge difference in bacterial diversity when comparing the Feces to other tissue types. 