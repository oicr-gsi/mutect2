version 1.0

workflow Mutect2GATK4 {
  input {
    File tumor_bam
    File normal_bam
    String outputTumorNamePrefix = basename(tumor_bam, '.bam')
    String outputNormalNamePrefix = basename(normal_bam, '.bam')
  }

  call runMuTect2GATK4 {
    input:
    tumor_bam = tumor_bam,
    normal_bam = normal_bam,
    outputPrefix = outputTumorNamePrefix
  }

  output {
    File outputMuTect2GATK4 = runMuTect2GATK4.outputMuTect2GATK4
  }

  parameter_meta {
  tumor_bam: "Input tumor file (bam or sam)."
  normal_bam: "Input normal file (bam or sam)."
  outputFileNamePrefix: "Output prefix to prefix output file names with."
}
  meta {
    author: "Alexander Fortuna"
    email: "alexander.fortuna@oicr.on.ca"
    description: "Workflow to run MuTect2GATK4"
    dependencies: [{
      name: "gatk/4.1.1.0",
      url: "https://software.broadinstitute.org/gatk/download/index"
    }]
  }
}

task runMuTect2GATK4 {
  input {
    File tumor_bam
    File normal_bam
    String gatk = "$GATK_ROOT/gatk"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String outputPrefix = "OUTPUT"
    String mutectTag = "mutect2_gatk"
    Int jobMemory = 18
    String modules = "gatk/4.1.2.0 hg19/p13 samtools/1.9"
  }

  parameter_meta {
  gatk: "gatk to use"
  tumor_bam: "Input tumor file (bam or sam)"
  normal_bam: "Input normal file (bam or sam)"
  refFasta: "Path to the reference fasta"
  mutectTag: "metric tag is used as a file extension for output"
  filter: "Picard filter to use"
  outputPrefix: "Output prefix, either input file basename or custom string"
  jobMemory: "memory allocated for Job"
  modules: "Environment module names and version to load (space separated) before command execution"
}

meta {
  output_meta : {
    outputMuTect2GATK4: "Workflow that takes output of BMPP (bam file) and calls SNPs, indels."
  }
}

command <<<

tumor_name=$(samtools view -H ~{tumor_bam} | grep '@RG')
normal_name=$(samtools view -H ~{normal_bam} | grep '@RG')

  ~{gatk} -Xmx~{jobMemory-6}G Mutect2 \
  -R ~{refFasta} \
  -I ~{tumor_bam} \
  -I ~{normal_bam} \
  -tumor $tumor_name \
  -normal $normal_name \
  -O "~{outputTumorNamePrefix}.~{mutectTag}.vcf"
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File outputMuTect2GATK4 = "~{outputTumorNamePrefix}.~{mutectTag}.vcf"
  }
}
