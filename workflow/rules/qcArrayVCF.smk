# # TO DO

# Evaluate GenTrain_Score and Cluster_Sep of final array VCF

# Zero out SNPs with GenTrain < 0.67.
# Manually inspect from 0.67 to 0.7 (the clusterSep)

# Zero out SNPs with cluster separation <= 0.425.
# Manually inspect from 0.4 to 0.45.

# Provide scatter plot of Gen Train and ClusterSep scores for ALL
# Provide scatter plot of Gen Train and ClusterSep scores for those to manually inspect - HTML Interactive
# Provide scatter plot of Gen Train and ClusterSep scores for those to keep GenTrain > 0.7 and ClusterSep > 0.4


# Provide plot for each snp in list of manual inspection thresholds
# install R script in Container from freeseek
	# sudo apt install r-cran-optparse r-cran-ggplot2 r-cran-data.table r-cran-gridextra
	# /bin/rm -f $HOME/bin/gtc2vcf_plot.R
	# wget -P $HOME/bin https://raw.githubusercontent.com/freeseek/gtc2vcf/master/gtc2vcf_plot.R
	# chmod a+x $HOME/bin/gtc2vcf_plot.R
# Plot variant for illumina
	# gtc2vcf_plot.R \
	# --illumina \
	# --vcf input.vcf \
	# --chrom 11 \
	# --pos 66328095 \
	# --png rs1815739.png


rule get_GenTrain_and_ClusterSep:
	input:
		vcf4="data/snp_array/vcf/4-extracted-dbSNP-rsids.vcf",
	output:
		QCvalues="data/snp_array/vcf/qcArray-GenTrain-and-ClusterSep-Scores.txt",
	# benchmark:
	# 	repeat("benchmarks/get_GenTrain_and_ClusterSep.array", 3),
	shell:
		"""
		bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/GenTrain_Score\t%INFO/Cluster_Sep\n' {input.vcf4} > data/snp_array/vcf/temp.GenTrainClusterSep.noHead.txt
		echo -e "CHROM\tPOS\tID\tREF\tALT\tGenTrain_Score\tCluster_Sep" | cat - data/snp_array/vcf/temp.GenTrainClusterSep.noHead.txt > {output.QCvalues}
		rm data/snp_array/vcf/temp.GenTrainClusterSep.noHead.txt
		"""

rule GenTrain_ClusterSep_QC:
	input:
		QCscores="data/snp_array/vcf/qcArray-GenTrain-and-ClusterSep-Scores.txt",
	output:
		fullScatterPlot="data/snp_array/qc/all-qc-scores-scatterplot.png",
		#- file of snps >= 0.7 GenTrain
		#- file of snps >= 0.45 ClusterSep
		keeperDF="data/snp_array/qc/SNPs-above-GenTrain-and-ClusterSep-Threshold-DataFrame.txt",
		keeperSNPs="data/snp_array/qc/SNPs-above-GenTrain-and-ClusterSep-Threshold-rsidLIST.txt",
		#- 0.67 =< file of snps < 0.7 GenTrain for manual checks
		#	- variant plot for each snp
		GenTrain_Manual_checksDF="data/snp_array/qc/manualChecks/GenTrain-SNPs-to-Check-DataFrame.txt",
		GenTrain_Manual_checksSNPs="data/snp_array/qc/manualChecks/GenTrain-SNPs-to-Check-rsidLIST.txt",
		#- 0.4 =< file of snps < 0.45 GenTrain for manual checks
		#	- variant plot for each snp
		ClusterSep_Manual_checksDF="data/snp_array/qc/manualChecks/ClusterSep-SNPs-to-Check-DataFrame.txt",
		ClusterSep_Manual_checksSNPs="data/snp_array/qc/manualChecks/ClusterSep-SNPs-to-Check-rsidLIST.txt",
		#- file of snps to EXCLUDE - 
		#- 0.45 GenTrain for manual checks
		#	- variant plot for each snp
		GenTrain_ToDropDF="data/snp_array/qc/toDrop/GenTrain-SNPs-to-DROP-DataFrame.txt",
		GenTrain_ToDropSNPs="data/snp_array/qc/toDrop/GenTrain-SNPs-to-DROP-rsidLIST.txt",
		ClusterSep_ToDropDF="data/snp_array/qc/toDrop/ClusterSep-SNPs-to-DROP-DataFrame.txt",
		ClusterSep_ToDropSNPs="data/snp_array/qc/toDrop/ClusterSep-SNPs-to-DROP-rsidLIST.txt",
	conda: "../../env/GenTrain_ClusterSep_QC.yaml",
	# benchmark:
	# 	repeat("benchmarks/GenTrain_ClusterSep_QC.array", 3),
	params:
		workdirPath=config['workdirPath'],
		qcDir="data/snp_array/qc/manualChecks/",
		dropDir="data/snp_array/qc/toDrop/",
		GenTrainUpperThreshold=config['QCarray']['GenTrainUpperThreshold'],
		GenTrainLowerThreshold=config['QCarray']['GenTrainLowerThreshold'],
		ClusterSepUpperThreshold=config['QCarray']['ClusterSepUpperThreshold'],
		ClusterSepLowerThreshold=config['QCarray']['ClusterSepLowerThreshold'],
	# benchmark:
	# 	repeat("benchmarks/GenTrain_ClusterSep_QC.array", 3),
	script:
		"../scripts/GenTrain-ClusterSep-QC.py"

# rule scatterplot:
# 	input:
# 		keeperDF="data/snp_array/qc/SNPs-above-GenTrain-and-ClusterSep-Threshold-DataFrame.txt",
# 	output:
# 		keeperScatterplot

rule remove_snps_under_thresholds:
	input:
		vcf4="data/snp_array/vcf/4-extracted-dbSNP-rsids.vcf",
		GenTrain_ToDropSNPs="data/snp_array/qc/toDrop/GenTrain-SNPs-to-DROP-rsidLIST.txt",
		ClusterSep_ToDropSNPs="data/snp_array/qc/toDrop/ClusterSep-SNPs-to-DROP-rsidLIST.txt",
	output:
		vcf5="data/snp_array/vcf/5-passing-QC-rsids.vcf",
		vcf5gz="data/snp_array/vcf/5-passing-QC-rsids.vcf.gz",
	# benchmark:
	# 	repeat("benchmarks/remove_snps_under_thresholds.array", 3),
	script:
		"../scripts/remove-non-passing-qc-snps.py"

rule clean_snpArray_VCF_index:
	input:
		vcf5gz="data/snp_array/vcf/5-passing-QC-rsids.vcf.gz",
	output:
		vcf5TBI="data/snp_array/vcf/5-passing-QC-rsids.vcf.gz.tbi",
	resources:
		mem_mb=1500,
		runtime="00:30:00",
	# benchmark:
	# 	repeat("benchmarks/clean_snpArray_VCF_index.array", 1)
	shell:
		"""
		bcftools index --tbi {input.vcf5gz} -o {output.vcf5TBI}
		"""
