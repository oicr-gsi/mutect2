version 1.0

workflow mutect2GATK4 {
  input {
    File tumorBam
    File? normalBam
    File? bedFile
    Int? numChunk
    String outputTumorNamePrefix = basename(tumorBam, '.bam')
    String? outputNormalNamePrefix = basename(select_first([normalBam, ""]), '.bam')
  }

  parameter_meta {
    tumorBam: "Input tumor file (bam or sam)."
    normalBam: "Input normal file (bam or sam)."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
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

  # Interval file provided, perform scatter/gather
  if(defined(bedFile)) {

    call splitIntervals {
      input:
        bedFile = bedFile,
        numChunk = numChunk
    }

    Array[File] intervals = splitIntervals.intervals

    scatter(interval in intervals) {
      if (!defined(normalBam)) {
        call tumorOnlyMode as tumorOnlyModeBed {
          input:
            tumorBam = tumorBam,
            interval = interval,
            outputTumorNamePrefix = outputTumorNamePrefix
        }
      }

      if (defined(normalBam)) {
        call tumorMatchedNormal as tumorMatchedNormalBed {
          input:
            tumorBam = tumorBam,
            normalBam = normalBam,
            interval = interval,
            outputNormalNamePrefix = outputNormalNamePrefix
        }
      }

      File? outputVcf = select_first([tumorOnlyMode.vcfFile, tumorMatchedNormal.vcfFile])
      File? outputVcfStats = select_first([tumorOnlyMode.vcfStats, tumorMatchedNormal.vcfStats])
    }

    call vcfMerge {
      input:
        vcfs = outputVcf,
        vcfStats = outputVcfStats
    }
  }

  # Interval file not provided, call Mutect2 without -L
  if(!defined(bedFile)) {

    if (!defined(normalBam)) {
      call tumorOnlyMode {
        input:
          tumorBam = tumorBam,
          outputTumorNamePrefix = outputTumorNamePrefix
      }
    }

    if (defined(normalBam)) {
      call tumorMatchedNormal {
        input:
          tumorBam = tumorBam,
          normalBam = normalBam,
          outputNormalNamePrefix = outputNormalNamePrefix
      }
    }

    File? unfilteredVcf = select_first([tumorOnlyMode.unfilteredVcf, tumorMatchedNormal.unfilteredVcf])
    File? filteredVcf = select_first([tumorOnlyMode.filteredVcf, tumorMatchedNormal.filteredVcf])
  }

  call vcfProcess {
    input:
      unfilteredVcf = select_first([vcfMerge.unfilteredVcf, unfilteredVcf]),
      filteredVcf = select_first([vcfMerge.filteredVcf, filteredVcf])
  }

  output {
    File unfilteredVcfFile = vcfProcess.unfilteredVcfGz
    File unfilteredVcfIndex = vcfProcess.unfilteredVcfIndex
    File filteredVcfFile = vcfProcess.filteredVcfGz
    File filteredVcfIndex = vcfProcess.filteredVcfIndex
  }
}

task splitIntervals {
  input {
    String modules = "gatk/4.1.2.0 hg19/p13"
    String gatk = "$GATK_ROOT/gatk"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    File? bedFile
    Int? numChunk
    Int memory = 32
    Int timeout = 72
  }

  parameter_meta {
    modules: "Environment module names and version to load (space separated) before command execution"
    gatk: "GATK to use"
    refFasta: "Path to the reference fasta"
    bedFile: "Interval file to split for scattering"
    numChunk: "Number of files to split the interval file into"
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

    ~{gatk} SplitIntervals \
    -R ~{refFasta} \
    -L ~{bedFile}
    --scatter-count ~{numChunk}
    -O interval-files

    cd interval-files
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Array[File] intervals = glob("*.intervals")
  }
}

task tumorOnlyMode {
  input {
    String modules = "gatk/4.1.2.0 hg19/p13 samtools/1.9"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String gatk = "$GATK_ROOT/gatk"
    String mutectTag = "mutect2_gatk"
    File tumorBam
    File? interval
    String outputTumorNamePrefix
    Int threads = 8
    Int memory = 32
    Int timeout = 96
  }

  parameter_meta {
      gatk: "GATK to use."
      tumorBam: "Input tumor file (bam or sam)."
      interval: "An interval or list of intervals to include in analysis."
      refFasta: "Path to the reference fasta."
      mutectTag: "Metric tag is used as a file extension for output."
      outputTumorNamePrefix: "Output prefix, either input file basename or custom string."
      memory: "Memory allocated for job."
      threads: "Requested CPU threads"
      timeout: "Hours before task timeout"
      modules: "Environment module names and version to load (space separated) before command execution."
    }

  meta {
    output_meta : {
      vcfFile: "VCF file output of BMPP (bam file) for tumor only mode."
    }
  }

  String vcf = "~{outputTumorNamePrefix}.~{mutectTag}.tumor_only.vcf"
  String vcfStats_ = "~{outputTumorNamePrefix}.~{mutectTag}.tumor_only.vcf.stats"

  command <<<
    tumor_name=$(samtools view -H ~{tumorBam} | grep '@RG')

    if [ -f "~{interval}" ]; then
      ~{gatk} -Xmx~{memory-6}G Mutect2 \
      -R ~{refFasta} \
      -I ~{tumorBam} \
      -tumor $tumor_name \
      -L ~{interval}
      -O ~{vcf}
    else
      ~{gatk} -Xmx~{memory-6}G Mutect2 \
      -R ~{refFasta} \
      -I ~{tumorBam} \
      -tumor $tumor_name \
      -O ~{outputTumorNamePrefix}.~{mutectTag}.tumor_only.unfiltered.vcf

      ~{gatk} FilterMutectCalls \
      -R ~{refFasta} \
      -V ~{outputTumorNamePrefix}.~{mutectTag}.tumor_only.unfiltered.vcf \
      -O ~{outputTumorNamePrefix}.~{mutectTag}.tumor_only.filtered.vcf
    fi
  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File? vcfFile = "~{vcf}"
    File? vcfStats = "~{vcfStats_}"
    File? unfilteredVcf = "~{outputTumorNamePrefix}.~{mutectTag}.tumor_only.unfiltered.vcf"
    File? filteredVcf = "~{outputTumorNamePrefix}.~{mutectTag}.tumor_only.filtered.vcf"
  }
}

task tumorMatchedNormal {
  input {
    String modules = "gatk/4.1.2.0 hg19/p13 samtools/1.9"
    File tumorBam
    File? normalBam
    File? interval
    String? outputNormalNamePrefix
    String gatk = "$GATK_ROOT/gatk"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String mutectTag = "mutect2_gatk"
    Int threads = 8
    Int memory = 32
    Int timeout = 96
  }

  parameter_meta {
    gatk: "GATK to use"
    tumorBam: "Input tumor file (bam or sam)"
    normalBam: "Input normal file (bam or sam)"
    interval: "An interval or list of intervals to include in analysis."
    refFasta: "Path to the reference fasta"
    mutectTag: "Metric tag is used as a file extension for output"
    outputNormalNamePrefix: "Output prefix, either input file basename or custom string"
    memory: "Memory allocated for job"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
    modules: "Environment module names and version to load (space separated) before command execution"
  }

  meta {
    output_meta : {
      vcfFile: "VCF file output of BMPP (bam file) for tumour with a matched normal."
    }
  }

  String vcf = "~{outputNormalNamePrefix}.~{mutectTag}.vcf"
  String vcfStats_ = "~{outputNormalNamePrefix}.~{mutectTag}.vcf.stats"

  command <<<
    tumor_name=$(samtools view -H ~{tumorBam} | grep '@RG')
    normal_name=$(samtools view -H ~{normalBam} | grep '@RG')

    if [ -f "~{interval}" ]; then
      ~{gatk} -Xmx~{memory-6}G Mutect2 \
      -R ~{refFasta} \
      -I ~{tumorBam} \
      -I ~{normalBam} \
      -tumor $tumor_name \
      -normal $normal_name \
      -L ~{interval}
      -O ~{vcf}
    else
      ~{gatk} -Xmx~{memory-6}G Mutect2 \
      -R ~{refFasta} \
      -I ~{tumorBam} \
      -I ~{normalBam} \
      -tumor $tumor_name \
      -normal $normal_name \
      -O ~{outputNormalNamePrefix}.~{mutectTag}.unfiltered.vcf

      ~{gatk} FilterMutectCalls \
      -R ~{refFasta} \
      -V ~{outputNormalNamePrefix}.~{mutectTag}.unfiltered.vcf \
      -O ~{outputNormalNamePrefix}.~{mutectTag}.filtered.vcf
    fi
  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File? vcfFile = "~{vcf}"
    File? vcfStats = "~{vcfStats_}"
    File? unfilteredVcf = "~{outputNormalNamePrefix}.~{mutectTag}.unfiltered.vcf"
    File? filteredVcf = "~{outputNormalNamePrefix}.~{mutectTag}.filtered.vcf"
  }
}

task vcfMerge {
  input {
    String modules = "gatk/4.1.2.0 hg19/p13"
    String gatk = "$GATK_ROOT/gatk"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    Array[File?] vcfs
    Array[File?] vcfStats
    Int memory = 32
    Int timeout = 72
  }

  parameter_meta {
    modules: "Environment module names and version to load (space separated) before command execution"
    gatk: "GATK to use"
    vcfs: "Vcf's from scatter to merge together"
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
  }

  meta {
    output_meta: {
      unfilteredVcf: "Merged vcf, unfiltered.",
      unfilteredVcfIndex: "Merged vcf index, unfiltered.",
      filteredVcf: "Merged vcf, processed through gatk FilterMutectCalls.",
      filteredVcfIndex: "Merged vcf index, processed through gatk FilterMutectCalls."
    }
  }

  String outputPrefix = basename(select_first([vcfs[0], ""]), ".vcf")

  command <<<
    ~{gatk} GatherVcfs \
    -I ~{sep="-I " vcfs}
    -O ~{outputPrefix}.unfiltered.vcf

    ~{gatk} MergeMutectStats \
    -stats ~{sep="-stats " vcfStats} \
    -O merged.stats

    ~{gatk} FilterMutectCalls \
    -R ~{refFasta} \
    -stats merged.stats \
    -V ~{outputPrefix}.unfiltered.vcf \
    -O ~{outputPrefix}.filtered.vcf
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputPrefix}.unfiltered.vcf"
    File filteredVcf = "~{outputPrefix}.filtered.vcf"
  }
}

task vcfProcess {
  input {
    String modules = "samtools/1.9"
    File unfilteredVcf
    File filteredVcf
    Int memory = 32
    Int timeout = 72
  }

  String unfilteredOutputPrefix = basename(unfilteredVcf, ".vcf")
  String filteredOutputPrefix = basename(filteredVcf, ".vcf")

  command <<<
    bgzip -c ~{unfilteredOutputPrefix}.vcf > ~{unfilteredOutputPrefix}.vcf.gz
    bgzip -c ~{filteredOutputPrefix}.vcf > ~{filteredOutputPrefix}.vcf.gz

    tabix -p vcf ~{unfilteredOutputPrefix}.vcf.gz
    tabix -p vcf ~{filteredOutputPrefix}.vcf.gz
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcfGz = "~{unfilteredOutputPrefix}.vcf.gz"
    File unfilteredVcfIndex = "~{unfilteredOutputPrefix}.vcf.gz.tbi"
    File filteredVcfGz = "~{filteredOutputPrefix}.vcf.gz"
    File filteredVcfIndex = "~{filteredOutputPrefix}.vcf.gz.tbi"
  }
}
