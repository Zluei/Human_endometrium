pyscenic grn \
	--num_workers 50 \
	--output adj.str_all_count.tsv \
	--method grnboost2 \
	str_all_count.loom ~/software/SCENIC/hs_hgnc_curated_tfs.txt

pyscenic ctx adj.str_all_count.tsv \
	/home/zhangleilei/software/SCENIC/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
	--annotations_fname ~/software/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
	--expression_mtx_fname str_all_count.loom --mode "dask_multiprocessing" \
	--output str_all_count_regulon.csv \
	--num_workers 50 \
	--mask_dropouts
pyscenic aucell str_all_count.loom str_all_count_regulon.csv \
	--output str_all_count_SCENIC.loom \
	--num_workers 50
