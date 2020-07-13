# mutect2

Somatic short variant analysis.

## Dependencies

* [gatk 4.1.1.0](https://software.broadinstitute.org/gatk/download/index)
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


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`normalBam`|File?|None|Input normal file (bam or sam).
`bedFile`|File?|None|Interval file to split for scattering.
`numChunk`|Int?|None|Number of files to split the interval file into.
`outputTumorNamePrefix`|String|basename(tumorBam,'.bam')|Output prefix, either input file basename or custom string.
`outputNormalNamePrefix`|String?|basename(select_first([normalBam, ""]),'.bam')|Output prefix, either input file basename or custom string.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitIntervals.modules`|String|"gatk/4.1.2.0 hg19/p13"|Environment module names and version to load (space separated) before command execution
`splitIntervals.gatk`|String|"$GATK_ROOT/gatk"|GATK to use
`splitIntervals.refFasta`|String|"$HG19_ROOT/hg19_random.fa"|Path to the reference fasta
`splitIntervals.memory`|Int|32|Memory allocated for job
`splitIntervals.timeout`|Int|72|Hours before task timeout
`tumorOnlyModeBed.modules`|String|"gatk/4.1.2.0 hg19/p13 samtools/1.9"|
`tumorOnlyModeBed.refFasta`|String|"$HG19_ROOT/hg19_random.fa"|
`tumorOnlyModeBed.gatk`|String|"$GATK_ROOT/gatk"|
`tumorOnlyModeBed.mutectTag`|String|"mutect2_gatk"|
`tumorOnlyModeBed.threads`|Int|8|
`tumorOnlyModeBed.memory`|Int|32|
`tumorOnlyModeBed.timeout`|Int|96|
`tumorMatchedNormalBed.modules`|String|"gatk/4.1.2.0 hg19/p13 samtools/1.9"|
`tumorMatchedNormalBed.gatk`|String|"$GATK_ROOT/gatk"|
`tumorMatchedNormalBed.refFasta`|String|"$HG19_ROOT/hg19_random.fa"|
`tumorMatchedNormalBed.mutectTag`|String|"mutect2_gatk"|
`tumorMatchedNormalBed.threads`|Int|8|
`tumorMatchedNormalBed.memory`|Int|32|
`tumorMatchedNormalBed.timeout`|Int|96|
`vcfMerge.modules`|String|"gatk/4.1.2.0 hg19/p13"|Environment module names and version to load (space separated) before command execution
`vcfMerge.gatk`|String|"$GATK_ROOT/gatk"|GATK to use
`vcfMerge.refFasta`|String|"$HG19_ROOT/hg19_random.fa"|
`vcfMerge.memory`|Int|32|Memory allocated for job
`vcfMerge.timeout`|Int|72|Hours before task timeout
`tumorOnlyMode.modules`|String|"gatk/4.1.2.0 hg19/p13 samtools/1.9"|Environment module names and version to load (space separated) before command execution.
`tumorOnlyMode.refFasta`|String|"$HG19_ROOT/hg19_random.fa"|Path to the reference fasta.
`tumorOnlyMode.gatk`|String|"$GATK_ROOT/gatk"|GATK to use.
`tumorOnlyMode.mutectTag`|String|"mutect2_gatk"|Metric tag is used as a file extension for output.
`tumorOnlyMode.interval`|File?|None|An interval or list of intervals to include in analysis.
`tumorOnlyMode.threads`|Int|8|Requested CPU threads
`tumorOnlyMode.memory`|Int|32|Memory allocated for job.
`tumorOnlyMode.timeout`|Int|96|Hours before task timeout
`tumorMatchedNormal.modules`|String|"gatk/4.1.2.0 hg19/p13 samtools/1.9"|Environment module names and version to load (space separated) before command execution
`tumorMatchedNormal.interval`|File?|None|An interval or list of intervals to include in analysis.
`tumorMatchedNormal.gatk`|String|"$GATK_ROOT/gatk"|GATK to use
`tumorMatchedNormal.refFasta`|String|"$HG19_ROOT/hg19_random.fa"|Path to the reference fasta
`tumorMatchedNormal.mutectTag`|String|"mutect2_gatk"|Metric tag is used as a file extension for output
`tumorMatchedNormal.threads`|Int|8|Requested CPU threads
`tumorMatchedNormal.memory`|Int|32|Memory allocated for job
`tumorMatchedNormal.timeout`|Int|96|Hours before task timeout
`vcfProcess.modules`|String|"samtools/1.9"|
`vcfProcess.memory`|Int|32|
`vcfProcess.timeout`|Int|72|


### Outputs

Output | Type | Description
---|---|---
`unfilteredVcfFile`|File|None
`unfilteredVcfIndex`|File|None
`filteredVcfFile`|File|None
`filteredVcfIndex`|File|None

## Rationale explained 

gatk SplitIntervals was used for splitting the interval file for scattering because it balances the genomic space interrogated by each of the chunks, versus just using an equal number of intervals in each.

## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with wdl_doc_gen (https://github.com/oicr-gsi/wdl_doc_gen/)_
