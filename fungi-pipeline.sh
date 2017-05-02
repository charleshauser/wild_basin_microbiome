# split libraries
split_libraries_fastq.py -v -q 0 --max_bad_run_length 250 --min_per_read_length_fraction 0.001 --store_demultiplexed_fastq -m $WB_fungi/012815CHits-mapping2.txt -i $WB_fungi/012815CHits-full.fastq -o sl_all2 -b $WB_fungi/processed_seqs/barcodes.fastq --barcode_type hamming_8 --phred_offset 33

# move into created directory
cd sl_all2

# get fastq stats
usearch7 -fastq_filter seqs.fastq -log seqs.stats.log.txt

# dereplicate and sort
usearch7 -derep_fulllength seqs.fna -output seqs_not_trimmed.filt.derep.fasta -sizeout

# remove singletons
usearch7 -sortbysize seqs_not_trimmed.filt.derep.fasta -output seqs_not_trimmed.filt.derep.mc2.fasta -minsize 2

# run ITSx direct
ITSx --reset T --cpu 8 --only_full T --save_regions all -t F -i seqs_not_trimmed.filt.derep.mc2.fasta -o ITSx.output2

# cluster OTUs
usearch7 -cluster_otus seqs.filt.derep.mc2.itsonly.fasta -otus seqs.filt.derep.mc2.itsonly.repset.fasta

# chimera UNITE
usearch7 -uchime_ref seqs.filt.derep.mc2.itsonly.repset.fasta -db $QIIME_DIR/unite/ITS1_ITS2_datasets/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS1.fasta -strand plus -minh 0.5 -nonchimeras seqs.filt.derep.mc2.itsonly.repset.nochimeras.fasta -chimeras seqs.filt.derep.mc2.itsonly.repset.chimeras.fasta -uchimealns seqs.filt.derep.mc2.itsonly.repset.chimeraalns.txt -threads 8

# relable representative sequence for each OTU with OTU label
python $BINF/bin/drive4_py/fasta_number.py seqs.filt.derep.mc2.itsonly.repset.nochimeras.fasta OTU_ > seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUs.fasta

# map original seqs to OTU representative sequence with 97% similarity
usearch7 -usearch_global seqs.fna -db seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUs.fasta -strand plus -id 0.97 -uc otu.map.uc -threads 8

# make tab-delimited OTU table
python $BINF/bin/drive5_py/uc2otutab_mod.py otu.map.uc > seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.txt

# convert tab-delimited OTU table to biom format
biom convert --table-type="OTU table" --to-json -i seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.txt -o seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.biom

# assign taxonomy using BLAST
parallel_assign_taxonomy_blast.py -v -O 4 -t $QIIME_DIR/unite/sh_taxonomy_qiime_ver7_dynamic_01.08.2015.txt -r $QIIME_DIR/unite/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta -i seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUs.fasta -o blast_assigned_taxonomy

# add taxonomy to biom table
biom add-metadata -i seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.biom -o seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.blasttax.biom --observation-metadata-fp blast_assigned_taxonomy_ver2/seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUs_tax_assignments.txt --observation-header OTUID,taxonomy --sc-separated taxonomy

# remove OTUs lacking "k__Fungi" taxonomy label from biom table
filter_taxa_from_otu_table.py -i seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.blasttax.biom -o seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.blasttax.onlyfungi.biom -p k__Fungi

# convert biom table with taxonomy to tab-delimited table
biom convert -i seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.blasttax.onlyfungi.biom -o seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.blasttax.onlyfungi.txt --to-tsv --header-key taxonomy --output-metadata-id ConsensusLineage

# move back up to original directory
cd ..

# summarize biom taxonomy table
biom summarize_table -i sl_all2/seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.blasttax.biom -o biome.table.summary.txt

# run core diversity analyses to produce rarefaction curve, histogram plots, and 3D PCoA plots
core_diversity_analyses.py -i sl_all2/seqs.filt.derep.mc2.itsonly.repset.notpadded.nochimeras.OTUTable.blasttax.biom -o cdfungi2/ -m $WB_fungi/012815CHits-mapping2.txt -e 28941 --nonphylogenetic_diversity

# get 2D PCoA plots using bray-curtis matrix made in core diversity analysis
make_2d_plots.py -i cdfungi2/bdiv_even28941/bray_curtis_pc.txt -m $WB_fungi/012815CHits-mapping2.txt -o 2d_plots/