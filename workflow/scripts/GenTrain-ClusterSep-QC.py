from os import path
from snakemake.shell import shell
import pandas as pd
import matplotlib.pyplot as plt

# Extract arguments.
# extra = snakemake.params.get("extra", "")
wd = snakemake.params.get("workdirPath", "")
qcDir = snakemake.params.get("qcDir", "")
dropDir = snakemake.params.get("dropDir", "")
GTUpper = snakemake.params.get("GenTrainUpperThreshold", "")
GTLower = snakemake.params.get("GenTrainLowerThreshold", "")
CSUpper = snakemake.params.get("ClusterSepUpperThreshold", "")
CSLower = snakemake.params.get("ClusterSepLowerThreshold", "")


log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "mkdir -p {wd}{qcDir}"
)

shell(
    "mkdir -p {wd}{dropDir}"
)

df = pd.read_table(wd + snakemake.input.QCscores, header =0)
x = df['GenTrain_Score'].to_numpy()
y = df['Cluster_Sep'].to_numpy()

plt.scatter(x, y, color='black', s=0.2)
plt.xlabel('GenTrain_Score')
plt.ylabel('Cluster_Sep')
# plt.axhline(y=0.45, xmin= 0.7, xmax=1, color='r', linestyle='-')
# plt.axhline(y=0.4, xmin=0.67, xmax=1, color='b', linestyle='-')
# plt.axvline(x=0.67, ymin=0.4, ymax=1, color='b', linestyle='-')
# plt.axvline(x=0.7, ymin=0.45, ymax=1, color='r', linestyle='-')
plt.axhline(y=0.45, color='r', linestyle='-')
plt.axhline(y=0.4, color='b', linestyle='-')
plt.axvline(x=0.67, color='b', linestyle='-')
plt.axvline(x=0.7, color='r', linestyle='-')
plt.savefig(wd + snakemake.output.fullScatterPlot)


filterDF = df[(df['GenTrain_Score'].astype(float) >= GTUpper) & (df['Cluster_Sep'].astype(float) >= CSUpper)]
filterDF.to_csv(wd + snakemake.output.keeperDF, sep='	', header=True, index=False, na_rep='N/A')

shell(
    "awk '{{print $3}}' {wd}{snakemake.output.keeperDF} > {wd}{snakemake.output.keeperSNPs}"
)


filter2DF = df[(df['GenTrain_Score'].astype(float) < GTUpper) & (df['GenTrain_Score'].astype(float) >= GTLower)]
filter2DF.to_csv(wd + snakemake.output.GenTrain_Manual_checksDF, sep='	', header=True, index=False, na_rep='N/A')

shell(
    "awk '{{print $3}}' {wd}{snakemake.output.GenTrain_Manual_checksDF} > {wd}{snakemake.output.GenTrain_Manual_checksSNPs}"
)


filter3DF = df[(df['Cluster_Sep'].astype(float) < CSUpper) & (df['Cluster_Sep'].astype(float) >= CSLower)]
filter3DF.to_csv(wd + snakemake.output.ClusterSep_Manual_checksDF, sep='	', header=True, index=False, na_rep='N/A')

shell(
    "awk '{{print $3}}' {wd}{snakemake.output.ClusterSep_Manual_checksDF} > {wd}{snakemake.output.ClusterSep_Manual_checksSNPs}"
)


#### TO DROP #####

filter4DF = df[(df['GenTrain_Score'].astype(float) < GTLower)]
filter4DF.to_csv(wd + snakemake.output.GenTrain_ToDropDF, sep='	', header=True, index=False, na_rep='N/A')

shell(
    "awk '{{print $3}}' {wd}{snakemake.output.GenTrain_ToDropDF} > {wd}{snakemake.output.GenTrain_ToDropSNPs}"
)


filter5DF = df[(df['Cluster_Sep'].astype(float) < CSLower)]
filter5DF.to_csv(wd + snakemake.output.ClusterSep_ToDropDF, sep='	', header=True, index=False, na_rep='N/A')

shell(
    "awk '{{print $3}}' {wd}{snakemake.output.ClusterSep_ToDropDF} > {wd}{snakemake.output.ClusterSep_ToDropSNPs}"
)