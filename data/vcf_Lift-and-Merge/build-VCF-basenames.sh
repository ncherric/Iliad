#!/bin/bash

mv ../../config/mergeTheseVCFs.txt ../../config/mergeTheseVCFs-previousProject.txt
ls *.{vcf,vcf.gz} | sed "s/.vcf.gz//" | sed "s/.vcf//" > ../../config/mergeTheseVCFs.txt