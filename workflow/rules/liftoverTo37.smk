NAMES=config['LiftoverSub']['filename']

# remove chr string 
# awk '{ gsub(/chr/,"", $1); print }' OFS="\t" /homer/herrick/EUR-face-ancestry/2022/2-reference-IDs/data/13-Last-overlap-EUR1kgenhgdp.vcf > /homer/herrick/EUR-face-ancestry/2022/2-reference-IDs/data/15-Last-overlap-EUR1kgenhgdp-noChrString.vcf

#add chr string
# awk '$1="chr"$1' OFS='\t' 7ABCD-chr-pos-from-1000genomes-38.txt > 7ABCD-with-chr-added-Regions-file-for-HGDP-overlap-extract.txt


rule query_positions:
    input:
        # startVCF="data/liftover/{name}.vcf.gz",
        startVCF=expand("data/liftover/{name}.vcf", name=NAMES),
    output:
        # positionList="data/liftover/{name}.peekPositionList.txt",
        rsidList="data/liftover/rsidList.{name}.txt",
    # script:
    #     "../scripts/query-positions.py"
    shell:
        """
        awk '{{print $3}}' {input.startVCF} > {output.rsidList}
        """

rule get_ids_from_dbsnp_file:
    input:
        # positionList="data/liftover/{name}.peekPositionList.txt",
        rsidList="data/liftover/rsidList.{name}.txt",
        tmpoFile="dbSNP/tempFile.to.remove",
        # get_reference_assembly_version
    output:
        # dbsnpGuideFile="data/liftover/Peek{name}.peekPositionList.txt",
        dbsnpExtractedIDsFile="data/liftover/dbSNP-IDs-from-{name}-rsIDList.txt",
    params:
        dbsnpDir="dbSNP/",
        dbsnpFile=config['dbSNP']['file'],
    shell:
        """
        zgrep -w -f {input.rsidList} {params.dbsnpDir}{params.dbsnpFile} > {output.dbsnpExtractedIDsFile}
        """

rule create_guide_file:
    input:
        dbsnpExtractedIDsFile="data/liftover/dbSNP-IDs-from-{name}-rsIDList.txt",
    output:
        guideFile="data/liftover/GuideFile-for-{name}-liftover.txt",
    script:
        "../scripts/create-guide-file-37.py"


rule ensure_chr_string_in_vcf:
    input:
        startVCF="data/liftover/{name}.vcf",
    output:
        startVCFwChr="data/liftover/{name}.chrString.vcf"
    script:
        "../scripts/ensure-no-chr-string-in-start-vcf-37.py"


rule final_liftover:
    input:
        guideFile="data/liftover/GuideFile-for-{name}-liftover.txt",
        startVCF="data/liftover/{name}.vcf",
    output:
        liftedVCF="data/liftover/{version}/liftedOver.{name}.vcf",
    params:
        refAssemblyVersion=config['LiftoverSub']['desiredVersion'],
    shell:
        """
        mkdir -p data/liftover/{params.refAssemblyVersion}
        awk 'NR==FNR {A[$1,$3] = $2; next} ($1,$3) in A{$2=A[$1,$3]}1' OFS='\t' {input.guideFile} {input.startVCF} > data/liftover/{params.refAssemblyVersion}/{output.liftedVCF}
        """

