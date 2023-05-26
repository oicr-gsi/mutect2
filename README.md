# mutect2

Somatic short variant analysis.

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
`gatk`|String|gatk version to be used
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
`gnomad`|File?|None|Genome Aggregation Database
`gnomadIdx`|File?|None|Index of gnomad


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitStringToArray.lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`splitStringToArray.memory`|Int|1|Memory allocated to job (in GB)
`splitStringToArray.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`runMutect2.mutectTag`|String|"mutect2"|version tag for mutect
`runMutect2.mutect2ExtraArgs`|String?|None|placehoulder for extra arguments
`runMutect2.threads`|Int|4|Requested CPU threads
`runMutect2.memory`|Int|32|Memory allocated to job (in GB).
`runMutect2.timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mergeVCFs.memory`|Int|4|Memory allocated for job
`mergeVCFs.timeout`|Int|12|Hours before task timeout
`mergeStats.memory`|Int|4|Memory allocated for job
`mergeStats.timeout`|Int|5|Hours before task timeout
`filter.filterExtraArgs`|String?|None|placehoulder for extra arguments
`filter.memory`|Int|16|Memory allocated for job
`filter.timeout`|Int|12|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`filteredVcfFile`|File|the filtered vcf file
`filteredVcfIndex`|File|Index of filtered vcf file
`mergedUnfilteredStats`|File|Stats for merged unfiltered files
`filteringStats`|File|Stats for filtering process


## Commands
 
 This section lists command(s) run by mutect2 workflow
 
  
  ```
      echo "~{intervalsToParallelizeBy}" | tr '~{lineSeparator}' '\n'
 ```
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
  
      gatk --java-options "-Xmx~{memory-8}g" Mutect2 \
      -R ~{refFasta} \
      $tumor_command_line \
      $normal_command_line \
      ~{"--germline-resource " + gnomad} \
      ~{"-pon " + pon} \
      $intervals_command_line \
      -O "~{outputVcf}" \
      ~{mutect2ExtraArgs}
    ```
  ```
      set -euo pipefail
  
      gatk --java-options "-Xmx~{memory-3}g" MergeVcfs \
      -I ~{sep=" -I " vcfs} \
      -O ~{outputName}
    ```
  ```
      set -euo pipefail
  
      gatk --java-options "-Xmx~{memory-3}g" MergeMutectStats \
      -stats ~{sep=" -stats " stats} \
      -O ~{outputStats}
    ```
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
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
