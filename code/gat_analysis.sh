# Genomic enrichment


module unload python/2.7.5
module load autoload python/3.5.2


#annotated files
        #broad_intersection_q05_merged_Final.bed
        #regular_intersection_q05_merged_Final.csv         


#regular
#simple repeats
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/simpleRepeat_hg19_newOrder.txt --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_simpleRepeats.log > regular_simpleRepeats_gat_enrich.bed


#repeat families
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/rmsk_hg19_sorted_newOrder.txt --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_rmsk.log > regular_rmsk_gat_enrich.bed


#Histone marks
#gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/data_drip/wgEncodeSydhHistoneMcf7H3k09me3UcdPk_sorted.txt --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_SydhHistoneMcf7H3k09me3.log > regular_SydhHistoneMcf7H3k09me3_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/data_drip/wgEncodeSydhHistoneMcf7H3k27acUcdPk.bed_1  --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_SydhHistoneMcf7H3k27ac.log > regular_SydhHistoneMcf7H3k27ac_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/data_drip/wgEncodeSydhHistoneMcf7H3k36me3bUcdPk.bed_1 --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_SydhHistoneMcf7H3k36me3b.log > regular_SydhHistoneMcf7H3k36me3b_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/data_drip/wgEncodeSydhHistoneU2osH3k36me3bUcdPk.bed_1 --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_SydhHistoneU2osH3k36me3b.log > regular_SydhHistoneU2osH3k36me3b_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/data_drip/wgEncodeSydhHistoneU2osH3k9me3UcdPk.bed_1 --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_SydhHistoneU2osH3k9me3.log > regular_SydhHistoneU2osH3k9me3_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/data_drip/wgEncodeSydhHistoneMcf7H3k27me3bUcdPk.bed_1 --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_SydhHistoneMcf7H3k27me3.log > regular_SydhHistoneMcf7H3k27me3_gat_enrich.bed


#Chromatine states
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/data_drip/wgEncodeBroadHmmK562HMM.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=.log > regular_HmmK562HMM_gat_enrich.bed


#G4
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_K_PDS_minus_hits_intersect.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_K_PDS_minus_gat_enrich..log> regular_Na_K_PDS_minus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_K_plus_hits_intersect.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_K_plus.log > regular_Na_K_plus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_PDS_minus_hits_intersect.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_PDS_minus.log  > regular_Na_PDS_minus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_PDS_plus_hits_intersect.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_PDS_plus.log  > regular_Na_PDS_plus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_K_minus_hits_intersect.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_K_minus.log > regular_Na_K_minus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_K_PDS_plus_hits_intersect.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_K_PDS.log > regular_Na_K_PDS_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_K_PDS_minus_hits_intersect_sorted.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_K_PDS_minus.log > regular_Na_K_PDS_minus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_K_plus_hits_intersect_sorted.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_K_plus.log > regular_Na_K_plus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_PDS_minus_hits_intersect_sorted.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_PDS_minus.log > regular_Na_PDS_minus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_PDS_plus_hits_intersect_sorted.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_PDS_plus.log > regular_Na_PDS_plus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_K_minus_hits_intersect_sorted.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_K_minus.log  > regular_Na_K_minus_gat_enrich.bed
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/G4-nonB_DB/GSE63874_Na_K_PDS_plus_hits_intersect_sorted.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_Na_K_PDS.log > regular_Na_K_PDS_plus_gat_enrich.bed


#Cpg islands
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/CpGIsland/cpgIslandExt_1.txt --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_cpgIsland.log > regular_cpgIsland_gat_enrich.bed


#Microsatellite
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/data_drip/Microsatellite.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_microsatellite.log > regular_microsatellite_gat_enrich.bed


#SegmentalDups
gat-run.py --segments=regular_intersection_q05_merged_Final.csv --annotations=../../data_annotation/data_drip/SegmentalDups.bed --workspace=../../data_annotation/hg19_chrs_coords.bed --ignore-segment-tracks --num-samples=1000 --log=regular_SegmentalDups.log > regular_SegmentalDups_gat_enrich.bed





