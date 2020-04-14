version 1.0

workflow mutect2GATK4 {
  input {
    File tumorBam
    File? normalBam
    File bedFile
    Int numChunk
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
    description: "Workflow to run MuTect2GATK4"
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

  call splitIntervals {
    input:
      bedFile = bedFile,
      numChunk = numChunk
  }

  Array[File] intervals = splitIntervals.intervals

  scatter(interval in intervals) {
    if (!defined(normalBam)) {
      call tumorOnlyMode {
        input:
          tumorBam = tumorBam,
          interval = interval,
          outputTumorNamePrefix = outputTumorNamePrefix
      }
    }

    if (defined(normalBam)) {
      call tumorMatchedNormal {
        input:
          tumorBam = tumorBam,
          normalBam = normalBam,
          interval = interval,
          outputNormalNamePrefix = outputNormalNamePrefix
      }
    }

    File outputVcf = select_first([tumorOnlyMode.vcfFile, tumorMatchedNormal.vcfFile])
  }

  call vcfMerge {
    input:
      vcfs = outputVcf,
      outputTumorNamePrefix = outputTumorNamePrefix
  }

  output {
    File outputMergedVcf = vcfMerge.outputMergedVcf
  }
}

task splitIntervals {
  input {
    String modules = "gatk/4.1.2.0 hg19/p13"
    String gatk = "$GATK_ROOT/gatk"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    File bedFile
    Int numChunk
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
    File interval
    String outputTumorNamePrefix
    Int memory = 18
    Int timeout = 72
  }

  parameter_meta {
      gatk: "GATK to use."
      tumorBam: "Input tumor file (bam or sam)."
      interval: "An interval or list of intervals to include in analysis."
      refFasta: "Path to the reference fasta."
      mutectTag: "Metric tag is used as a file extension for output."
      outputTumorNamePrefix: "Output prefix, either input file basename or custom string."
      memory: "Memory allocated for job."
      timeout: "Hours before task timeout"
      modules: "Environment module names and version to load (space separated) before command execution."
    }

    meta {
      output_meta : {
        vcfFile: "VCF file output of BMPP (bam file) for tumor only mode."
      }
    }

  command <<<
    tumor_name=$(samtools view -H ~{tumorBam} | grep '@RG')

    ~{gatk} -Xmx~{memory-6}G Mutect2 \
    -R ~{refFasta} \
    -I ~{tumorBam} \
    -tumor $tumor_name \
    -L ~{interval}
    -O ~{outputTumorNamePrefix}.~{mutectTag}.vcf
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File vcfFile = "~{outputTumorNamePrefix}.~{mutectTag}.vcf"
  }
}

task tumorMatchedNormal {
  input {
    String modules = "gatk/4.1.2.0 hg19/p13 samtools/1.9"
    File tumorBam
    File? normalBam
    File interval
    String? outputNormalNamePrefix
    String gatk = "$GATK_ROOT/gatk"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String mutectTag = "mutect2_gatk"
    Int memory = 18
    Int timeout = 72
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
    timeout: "Hours before task timeout"
    modules: "Environment module names and version to load (space separated) before command execution"
  }

  meta {
    output_meta : {
      vcfFile: "VCF file output of BMPP (bam file) for tumour with a matched normal."
    }
  }

  command <<<
    tumor_name=$(samtools view -H ~{tumorBam} | grep '@RG')
    normal_name=$(samtools view -H ~{normalBam} | grep '@RG')

    ~{gatk} -Xmx~{memory-6}G Mutect2 \
    -R ~{refFasta} \
    -I ~{tumorBam} \
    -I ~{normalBam} \
    -normal $normal_name \
    -L ~{interval}
    -O ~{outputNormalNamePrefix}.~{mutectTag}.vcf
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File vcfFile = "~{outputNormalNamePrefix}.~{mutectTag}.vcf"
  }
}

task vcfMerge {
  input {
    String modules = "gatk/4.1.2.0"
    String gatk = "$GATK_ROOT/gatk"
    Array[File] vcfs
    String outputTumorNamePrefix
    Int memory = 32
    Int timeout = 72
  }

  parameter_meta {
    modules: "Environment module names and version to load (space separated) before command execution"
    gatk: "GATK to use"
    vcfs: "Vcf's from scatter to merge together"
    outputTumorNamePrefix: "Output prefix, either input file basename or custom string"
    memory: "Memory allocated for job"
    timeout: "Hours before task timeout"
  }

  meta {
    output_meta: {
      outputMergedVcf: "Merged vcf"
    }
  }

  command <<<
    ~{gatk} GatherVcfs \
    -I ~{sep="-I " vcfs}
    -O ~{outputTumorNamePrefix}.merged.vcf
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File outputMergedVcf = "~{outputTumorNamePrefix}.merged.vcf"
  }
}

