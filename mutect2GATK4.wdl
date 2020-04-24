version 1.0

workflow mutect2GATK4 {
  input {
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
    File? intervals
    Int? scatterCount
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

  String outputBasename = basename(tumorBam, '.bam')

  call splitIntervals {
    input:
      intervals = intervals,
      scatterCount = scatterCount
  }

  scatter(subintervals in splitIntervals.intervalFiles) {
    call runMutect2 {
      input:
        intervals = subintervals,
        tumorBam = tumorBam,
        tumorBai = tumorBai,
        normalBam = normalBam,
        normalBai = normalBai,
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
      intervals = intervals,
      unfilteredVcf = mergeVCFs.mergedVcf,
      unfilteredVcfIdx = mergeVCFs.mergedVcfIdx,
      mutectStats = mergeStats.mergedStats
  }


  output {
    File unfilteredVcfFile = filter.unfilteredVcfGz
    File unfilteredVcfIndex = mergeVCFs.mergedVcfIdx
    File filteredVcfFile = filter.filteredVcfGz
    File filteredVcfIndex = filter.filteredVcfIdx
    File mergedUnfilteredStats = mergeStats.mergedStats
    File filteringStats = filter.filteringStats
  }
}

task splitIntervals {
  input {
    String modules = "gatk/4.1.6.0 hg19/p13"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String refFai = "$HG19_ROOT/hg19_random.fa.fai"
    String refDict = "$HG19_ROOT/hg19_random.dict"
    File? intervals
    Int? scatterCount
    String? splitIntervalsExtraArgs
    Int memory = 16
    Int timeout = 72
  }

  parameter_meta {
    modules: "Environment module names and version to load (space separated) before command execution"
    refFasta: "Path to the reference fasta"
    intervals: "Interval file to split for scattering"
    scatterCount: "Number of files to split the interval file into"
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
  }

  meta {
    output_meta: {
      intervals: "Interval files that were split to be used for scattering."
    }
  }

  command <<<
    mkdir interval-files

    cp ~{refFai} .
    cp ~{refDict} .

    gatk --java-options "-Xmx~{memory}g" SplitIntervals \
    -R ~{refFasta} \
    ~{"-L " + intervals} \
    --scatter-count ~{scatterCount} \
    -O interval-files \
    ~{splitIntervalsExtraArgs}

    cp interval-files/*.interval_list .
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Array[File] intervalFiles = glob("*.interval_list")
  }
}

task runMutect2 {
  input {
    String modules = "gatk/4.1.6.0 hg19/p13 samtools/1.9"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String refFai = "$HG19_ROOT/hg19_random.fa.fai"
    String refDict = "$HG19_ROOT/hg19_random.dict"
    String mutectTag = "mutect2_gatk"
    File? intervals
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
    String? mutect2ExtraArgs
    String outputBasename
    Int threads = 8
    Int memory = 32
    Int timeout = 96
  }

  String outputVcf = outputBasename + "." + mutectTag + ".vcf"
  String outputVcfIdx = outputVcf + ".idx"
  String outputStats = outputVcf + ".stats"

  command <<<
    tumor_name=$(samtools view -H ~{tumorBam} | grep '^@RG'| sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
    tumor_command_line="-I ~{tumorBam} -tumor ${tumor_name}"

    cp ~{refFai} .
    cp ~{refDict} .

    if [ -f ~{normalBam} ]; then
      normal_name=$(samtools view -H ~{normalBam} | grep '^@RG'| sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
      normal_command_line="-I ~{normalBam} -normal ${normal_name}"
    fi

    gatk --java-options "-Xmx~{memory-6}g" Mutect2 \
    -R ~{refFasta} \
    $tumor_command_line \
    $normal_command_line \
    ~{"-L " + intervals} \
    -O "~{outputVcf}"
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
    Int memory = 8
    Int timeout = 72
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

  String outputName = basename(vcfs[0], ".vcf")

  command <<<
    gatk --java-options "-Xmx~{memory}g" MergeVcfs \
    -I ~{sep=" -I " vcfs} \
    -O ~{outputName}.vcf
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
    Int timeout = 72
  }

  String outputStats = basename(stats[0])
  command <<<
    gatk --java-options "-Xmx~{memory}g" MergeMutectStats \
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
    File? intervals
    File unfilteredVcf
    File unfilteredVcfIdx
    File mutectStats
    String? filterExtraArgs
    Int memory = 4
    Int timeout = 72
  }

  String unfilteredVcfName = basename(unfilteredVcf)
  String filteredVcfName = basename(unfilteredVcf, ".vcf") + ".filtered.vcf"

  command <<<
    cp ~{refFai} .
    cp ~{refDict} .

    gatk --java-options "-Xmx~{memory}g" FilterMutectCalls -V \
    -R ~{refFasta} \
    -O ~{filteredVcfName} \
    ~{"-stats " + mutectStats} \
    --filtering-stats ~{filteredVcfName}.stats \
    ~{filterExtraArgs}

    bgzip -c ~{filteredVcfName} > ~{filteredVcfName}.gz
    bgzip -c ~{unfilteredVcfName} > ~{unfilteredVcfName}.gz
  >>>

  output {
    File unfilteredVcfGz = "~{unfilteredVcfName}.gz"
    File filteredVcfGz = "~{filteredVcfName}.gz"
    File filteredVcfIdx = "~{filteredVcfName}.idx"
    File filteringStats = "~{filteredVcfName}.stats"
  }
}
