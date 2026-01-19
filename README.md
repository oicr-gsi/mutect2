# mutect2

muTect2 is a tool from GATK suite. muTect2 calls somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations. The caller uses a Bayesian somatic genotyping model. This wrapper workflow accepts either paired tumor/normal inputs or tumor-only inputs in BAM format. It also accepts optional Panel Of Normals and other parameters. The workflow outputs a filtered vcf along with index and metrics files.

## Overview

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools/archive/0.1.19.tar.gz)


## Usage

### Cromwell
```
java -jar cromwell.jar run mutect2.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumorBam`|File|Input tumor file (bam or sam).
`tumorBai`|File|Index for tumorBam
`reference`|String|the reference genome for input sample
`outputFileNamePrefix`|String|prefix of output file


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`normalBam`|File?|None|Input normal file (bam or sam).
`normalBai`|File?|None|Index for noramlBam
`intervalFile`|String?|None|One or more genomic intervals over which to operate
`intervalsToParallelizeBy`|String?|None|Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4)
`pon`|File?|None|panel of normal
`ponIdx`|File?|None|index of pon


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitStringToArray.lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`splitStringToArray.memory`|Int|1|Memory allocated to job (in GB)
`splitStringToArray.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`getChrCoefficient.memory`|Int|1|Memory allocated for this job
`getChrCoefficient.timeout`|Int|1|Hours before task timeout
`runMutect2.mutectTag`|String|"mutect2"|version tag for mutect
`runMutect2.mutect2ExtraArgs`|String?|None|placehoulder for extra arguments
`runMutect2.threads`|Int|4|Requested CPU threads
`runMutect2.memory`|Int|32|Memory allocated to job (in GB).
`runMutect2.minMemory`|Int|6|A minimum amount of memory allocated to the task, overrides the scaled RAM setting
`runMutect2.timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mergeVCFs.memory`|Int|4|Memory allocated for job
`mergeVCFs.timeout`|Int|12|Hours before task timeout
`mergeStats.memory`|Int|4|Memory allocated for job
`mergeStats.timeout`|Int|5|Hours before task timeout
`filter.filterExtraArgs`|String?|None|placehoulder for extra arguments
`filter.memory`|Int|16|Memory allocated for job
`filter.timeout`|Int|12|Hours before task timeout
`mutect2Metrics.modules`|String|"python/3.10.6 pandas/2.1.3"|Names and versions of modules
`mutect2Metrics.memory`|Int|4|Memory allocated for job
`mutect2Metrics.timeout`|Int|2|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`filteredVcfFile`|File|the filtered vcf file|vidarr_label: filteredVcfFile
`filteredVcfIndex`|File|Index of filtered vcf file|vidarr_label: filteredVcfIndex
`mergedUnfilteredStats`|File|Stats for merged unfiltered files|vidarr_label: mergedUnfilteredStats
`filteringStats`|File|Stats for filtering process|vidarr_label: filteringStats
`metricsFile`|File|Summary metrics derived from the filtered Mutect2 VCF.|


./commands.txt found, printing out the content...
## Commands
  
 This section lists command(s) run by mutect2 workflow
 
 
 * Running mutect2
 
 Mutect requires alignment files (.bam) as an input. Also, various other resources are needed. The results produced by the workflow
 are variant calls, provisioned as a .vcf file.
 
 ### Split interval string by comma
 
 ```
     echo "~{intervalsToParallelizeBy}" | tr '~{lineSeparator}' '\n'
 ```
 
 ### Getting Chromosome Scaling Coefficient
 
 ```
     CHROM=$(echo ~{region} | sed 's/:.*//')
     LARGEST=$(grep SN:chr ~{refDict} | cut -f 3 | sed 's/LN://' | sort -n | tail -n 1)
     grep -w SN:$CHROM ~{refDict} | cut -f 3 | sed 's/.*://' | awk -v largest_chr=$LARGEST '{print int(($1/largest_chr + 0.1) * 10)/10}'
 ```
 
 ### Mutect2 (GATK)
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx1g" GetSampleName -R ~{refFasta} -I ~{tumorBam} -O tumor_name.txt -encode
     tumor_command_line="-I ~{tumorBam} -tumor `cat tumor_name.txt`"
 
     cp ~{refFai} .
     cp ~{refDict} .
 
     if [ -f "~{normalBam}" ]; then
       gatk --java-options "-Xmx1g" GetSampleName -R ~{refFasta} -I ~{normalBam} -O normal_name.txt -encode
       normal_command_line="-I ~{normalBam} -normal `cat normal_name.txt`"
     else
       normal_command_line=""
     fi
 
     if [ -f "~{intervalFile}" ]; then
       if ~{intervalsProvided} ; then
         intervals_command_line="-L ~{sep=" -L " intervals} -L ~{intervalFile} -isr INTERSECTION"
       else
         intervals_command_line="-L ~{intervalFile}"
       fi
     else
       if ~{intervalsProvided} ; then
         intervals_command_line="-L ~{sep=" -L " intervals} "
       fi
     fi
 
     if [[ ! -z "~{gnomad}" ]] ; then
       germline_resource_line="--germline-resource ~{gnomad}"
     else 
       germline_resource_line=""
     fi
 
     gatk --java-options "-Xmx~{memory-8}g" Mutect2 \
     -R ~{refFasta} \
     $tumor_command_line \
     $normal_command_line \
     $germline_resource_line \
     ~{"-pon " + pon} \
     $intervals_command_line \
     -O "~{outputVcf}" \
     ~{mutect2ExtraArgs}
 ```
 
 ### MergeVcfs (GATK)
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{memory-3}g" MergeVcfs \
     -I ~{sep=" -I " vcfs} \
     -O ~{outputName}
 ```
 
 ### MergeMutectStats (GATK)
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{memory-3}g" MergeMutectStats \
     -stats ~{sep=" -stats " stats} \
     -O ~{outputStats}
 ```
 
 ### FilterMutectCalls (GATK)
 
 ```
     set -euo pipefail
 
     cp ~{refFai} .
     cp ~{refDict} .
 
     gatk --java-options "-Xmx~{memory-4}g" FilterMutectCalls \
     -V ~{unfilteredVcf} \
     -R ~{refFasta} \
     -O ~{filteredVcfName} \
     ~{"-stats " + mutectStats} \
     --filtering-stats ~{filteredVcfName}.stats \
     ~{filterExtraArgs}
 
     bgzip -c ~{filteredVcfName} > ~{filteredVcfName}.gz
     bgzip -c ~{unfilteredVcf} > ~{unfilteredVcfName}.gz
 
     gatk --java-options "-Xmx~{memory-5}g" IndexFeatureFile -I ~{filteredVcfName}.gz
     gatk --java-options "-Xmx~{memory-5}g" IndexFeatureFile -I ~{unfilteredVcfName}.gz
 ```
 
 ### Pull metrics from filtered vcf files 
 ```
 set -euo pipefail
 
   python3 - <<CODE
   import pandas as pd
   import numpy as np
   import sys
 
   NUM_DP = 2
   FILTER = 6
   REF = 3
   ALT = 4
 
   vcf_file = "~{filteredVcf}"
   out_tsv = "~{metricsTsvName}"
 
   results = {
       "num_calls": 0,
       "num_PASS": 0,
       "num_SNPs": 0,
       "num_multi_SNPs": 0,
       "num_indels": 0,
       "titv_ratio": np.nan,
       "num_MNPs": 0,
       "pct_ti": 0,
       "pct_tv": 0,
   }
 
   total_ti = 0
   total_tv = 0
   chunksize = 5 * 10**4
 
   try:
     df_iter = pd.read_csv(
       vcf_file,
       sep="\t",
       comment="#",
       header=None,
       compression="gzip",
       chunksize=chunksize,
     )
   except pd.errors.EmptyDataError:
     pd.DataFrame(results, index=[0]).to_csv(out_tsv, sep="\t", index=False)
     sys.exit(0)
 
   for sub_df in df_iter:
     results["num_calls"] += len(sub_df)
 
     pass_df = sub_df[sub_df[FILTER] == "PASS"]
     results["num_PASS"] += len(pass_df)
 
     if len(pass_df) == 0:
         continue
 
     snps = pass_df[
         (pass_df[REF].str.len() == 1) &
         (pass_df[ALT].str.len() == 1)
     ]
     results["num_SNPs"] += len(snps)
 
     results["num_multi_SNPs"] += len(
         pass_df[
             (pass_df[REF].str.len() == 1) &
             (pass_df[ALT].str.contains(","))
         ]
     )
 
     results["num_MNPs"] += len(
         pass_df[
             (pass_df[REF].str.len() > 1) &
             (pass_df[REF].str.len() == pass_df[ALT].str.len())
         ]
     )
 
     total_ti += len(
         snps[(snps[REF] + snps[ALT]).isin(["AG", "GA", "CT", "TC"])]
     )
     total_tv += len(
         snps[(snps[REF] + snps[ALT]).isin(
             ["AC", "AT", "CA", "CG", "GT", "GC", "TA", "TG"]
         )]
     )
 
   results["num_indels"] = (
     results["num_PASS"]
     - results["num_SNPs"]
     - results["num_multi_SNPs"]
     - results["num_MNPs"]
   )
 
   if results["num_PASS"] > 0:
     results["pct_ti"] = round(total_ti / results["num_PASS"], NUM_DP)
     results["pct_tv"] = round(total_tv / results["num_PASS"], NUM_DP)
 
   results["titv_ratio"] = (
     np.nan if total_tv == 0 else round(total_ti / total_tv, NUM_DP)
   )
 
   pd.DataFrame(results, index=[0]).to_csv(out_tsv, sep="\t", index=False)
   CODE
 ``` ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_