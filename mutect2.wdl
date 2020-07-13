version 1.0

workflow mutect2 {
  input {
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
    File? intervalFile
    String? intervalsToParallelizeBy
    File? pon
    File? ponIdx
    File? gnomad
    File? gnomadIdx
  }

  parameter_meta {
    tumorBam: "Input tumor file (bam or sam)."
    normalBam: "Input normal file (bam or sam)."
  }

  meta {
    author: "Angie Mosquera, Alexander Fortuna"
    email: "amosquera@oicr.on.ca, afortuna@oicr.on.ca"
    description: "Somatic short variant analysis."
    dependencies: [
    {
      name: "gatk/4.1.1.0",
      url: "https://software.broadinstitute.org/gatk/download/index"
    },
    {
      name: "samtools/1.9",
      url: "https://github.com/samtools/samtools/archive/0.1.19.tar.gz"
    }]
  }

  call splitStringToArray {
    input:
      intervalsToParallelizeBy = intervalsToParallelizeBy
  }

  String outputBasename = basename(tumorBam, '.bam')
  Boolean intervalsProvided = if (defined(intervalsToParallelizeBy)) then true else false

  scatter(subintervals in splitStringToArray.out) {
    call runMutect2 {
      input:
        intervals = subintervals,
        intervalsProvided = intervalsProvided,
        intervalFile = intervalFile,
        tumorBam = tumorBam,
        tumorBai = tumorBai,
        normalBam = normalBam,
        normalBai = normalBai,
        pon = pon,
        ponIdx = ponIdx,
        gnomad = gnomad,
        gnomadIdx = gnomadIdx,
        outputBasename = outputBasename
    }
  }

  Array[File] unfilteredVcfs = runMutect2.unfilteredVcf
  Array[File] unfilteredVcfIndices = runMutect2.unfilteredVcfIdx
  Array[File] unfilteredStats = runMutect2.stats

  call mergeVCFs {
    input:
      vcfs = unfilteredVcfs,
      vcfIndices = unfilteredVcfIndices
  }

  call mergeStats {
    input:
      stats = unfilteredStats
  }

  call filter {
    input:
      intervalFile = intervalFile,
      unfilteredVcf = mergeVCFs.mergedVcf,
      unfilteredVcfIdx = mergeVCFs.mergedVcfIdx,
      mutectStats = mergeStats.mergedStats
  }


  output {
    File unfilteredVcfFile = filter.unfilteredVcfGz
    File unfilteredVcfIndex = filter.unfilteredVcfTbi
    File filteredVcfFile = filter.filteredVcfGz
    File filteredVcfIndex = filter.filteredVcfTbi
    File mergedUnfilteredStats = mergeStats.mergedStats
    File filteringStats = filter.filteringStats
  }
}

task splitStringToArray {
  input {
    String? intervalsToParallelizeBy
    String lineSeparator = ","
    Int memory = 1
    Int timeout = 1
    String modules = ""
  }

  command <<<
    echo "~{intervalsToParallelizeBy}" | tr '~{lineSeparator}' '\n'
  >>>

  output {
    Array[Array[String]] out = read_tsv(stdout())
  }
}

task runMutect2 {
  input {
    String modules = "gatk/4.1.6.0 hg19/p13"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String refFai = "$HG19_ROOT/hg19_random.fa.fai"
    String refDict = "$HG19_ROOT/hg19_random.dict"
    String mutectTag = "mutect2"
    File? intervalFile
    Array[String]? intervals
    Boolean intervalsProvided
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
    File? pon
    File? ponIdx
    File? gnomad
    File? gnomadIdx
    String? mutect2ExtraArgs
    String outputBasename
    Int threads = 4
    Int memory = 32
    Int timeout = 24
  }

  String outputVcf = if (defined(normalBam)) then outputBasename + "." + mutectTag + ".vcf" else outputBasename + "." + mutectTag + ".tumor_only.vcf"
  String outputVcfIdx = outputVcf + ".idx"
  String outputStats = outputVcf + ".stats"

  command <<<
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
  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputVcf}"
    File unfilteredVcfIdx = "~{outputVcfIdx}"
    File stats = "~{outputStats}"
  }
}

task mergeVCFs {
  input {
    String modules = "gatk/4.1.6.0 hg19/p13"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    Array[File] vcfs
    Array[File] vcfIndices
    Int memory = 4
    Int timeout = 12
  }

  parameter_meta {
    modules: "Environment module names and version to load (space separated) before command execution"
    vcfs: "Vcf's from scatter to merge together"
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
  }

  meta {
    output_meta: {
      mergedVcf: "Merged vcf, unfiltered.",
      mergedVcfIdx: "Merged vcf index, unfiltered."
    }
  }

  String outputName = basename(vcfs[0])

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{memory-3}g" MergeVcfs \
    -I ~{sep=" -I " vcfs} \
    -O ~{outputName}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedVcf = "~{outputName}"
    File mergedVcfIdx = "~{outputName}.idx"
  }
}

task mergeStats {
  input {
    String modules = "gatk/4.1.6.0"
    Array[File]+ stats
    Int memory = 4
    Int timeout = 5
  }

  String outputStats = basename(stats[0])

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{memory-3}g" MergeMutectStats \
    -stats ~{sep=" -stats " stats} \
    -O ~{outputStats}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedStats = "~{outputStats}"
  }
}

task filter {
  input {
    String modules = "gatk/4.1.6.0 hg19/p13 samtools/1.9"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String refFai = "$HG19_ROOT/hg19_random.fa.fai"
    String refDict = "$HG19_ROOT/hg19_random.dict"
    File? intervalFile
    File unfilteredVcf
    File unfilteredVcfIdx
    File mutectStats
    String? filterExtraArgs
    Int memory = 16
    Int timeout = 12
  }

  String unfilteredVcfName = basename(unfilteredVcf)
  String filteredVcfName = basename(unfilteredVcf, ".vcf") + ".filtered.vcf"

  command <<<
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
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcfGz = "~{unfilteredVcfName}.gz"
    File unfilteredVcfTbi = "~{unfilteredVcfName}.gz.tbi"
    File filteredVcfGz = "~{filteredVcfName}.gz"
    File filteredVcfTbi = "~{filteredVcfName}.gz.tbi"
    File filteringStats = "~{filteredVcfName}.stats"
  }
}
