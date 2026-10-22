# TRANA: Usage

## Introduction

TRANA is a Nextflow pipeline for taxonomic profiling of 16S rRNA gene sequences with EMU. FastQC performs read quality control. For Oxford Nanopore reads, NanoPlot provides additional quality reports, optional Porechop_ABI processing trims adapters, and Filtlong filters reads by length and quality. For short reads, Cutadapt trims primers unless disabled. EMU estimates taxonomic abundances, and MultiQC summarises quality-control results.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with a header row. The first three columns are `sample`, `fastq_1`, and `fastq_2`; add the optional fourth column `neg_control` to associate samples with negative controls. Existing three-column samplesheets remain valid.

```bash
--input '[path to samplesheet file]'
```

### Input file example

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline detects whether a sample is single- or paired-end from the presence of `fastq_2`. Leave `fastq_2` empty for single-end reads, including Nanopore reads. Use the column names shown below. Sample names and file paths must not contain spaces.

Set `--seqtype` for the entire run: use separate runs for Nanopore and short-read data. Single- and paired-end short-read samples can share a samplesheet when running with `--seqtype sr`.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column        | Description                                                                                                                                          |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`      | Custom sample name without spaces. Use the same name for multiple sequencing libraries/runs from the same sample.                                    |
| `fastq_1`     | Path to single-end reads or the first file of paired-end reads. Must be gzipped with the extension `.fastq.gz` or `.fq.gz`.                          |
| `fastq_2`     | Path to the second file of paired-end reads. Must be gzipped with the extension `.fastq.gz` or `.fq.gz`. Leave empty for single-end reads.           |
| `neg_control` | Optional fourth column. Exact `sample` name of the associated negative control, without spaces. Leave empty for negative-control samples.            |
| `FW_primer`   | Optional forward primer, mapped to `meta.fw_primer`. A nonempty value overrides `--FW_primer`; a blank or omitted column uses the parameter instead. |
| `RV_primer`   | Optional reverse primer, mapped to `meta.rv_primer`. A nonempty value overrides `--RV_primer`; a blank or omitted column uses the parameter instead. |

An [example samplesheet](../assets/examples/samplesheet_medium.csv) has been provided with the pipeline.

### Per-sample negative controls

Place `neg_control` after `fastq_1` and `fastq_2`. Include each negative control as its own sample row with read files, but leave its `neg_control` cell empty. Each named control must exactly match a `sample` entry in the same run's samplesheet; the pipeline rejects missing controls before emitting reads to downstream processes. Controls can appear anywhere in the samplesheet, and several samples can share a control. For samples sequenced more than once, use the same control assignment on every row. The pipeline rejects different controls or a mix of assigned and blank controls for the same sample. Repeated rows with all assignments blank or missing remain valid.

```csv
sample,fastq_1,fastq_2,neg_control
SAMPLE_01,/data/sample_01_R1.fastq.gz,/data/sample_01_R2.fastq.gz,NEG_01
SAMPLE_02,/data/sample_02.fastq.gz,,NEG_01
NEG_01,/data/neg_01.fastq.gz,,
```

The control name identifies another sample, not a file path. Blank entries do not assign a control; samplesheets without the column are still accepted.

## Running the pipeline

Run from a local checkout with your samplesheet, EMU database directory, and Krona taxonomy file:

```bash
nextflow run main.nf \
  --input /data/samplesheet.csv \
  --outdir results_emu \
  --db /data/emu_database \
  --krona_taxonomy_tab /data/krona/taxonomy.tab \
  --seqtype map-ont \
  -profile singularity
```

The EMU database directory must contain `species_taxid.fasta` and `taxonomy.tsv`. Use `--run_krona false` if Krona plots are not needed; in that case, omit `--krona_taxonomy_tab`.

This command uses the `singularity` profile. Replace it with `docker` to run with Docker. Do not add the `test` profile to a production run, since it supplies test inputs and a test database.

To run the bundled test configuration from a local checkout:

```bash
nextflow run main.nf -profile test,docker --outdir results_test
```

### Sequencing type

| `--seqtype` | Input                             | Processing                                                                                 |
| ----------- | --------------------------------- | ------------------------------------------------------------------------------------------ |
| `map-ont`   | Nanopore reads (default)          | Optional adapter trimming, Filtlong filtering, and EMU with the minimap2 `map-ont` preset. |
| `lr:hq`     | High-accuracy Nanopore reads      | The same preprocessing as `map-ont`, with the minimap2 `lr:hq` preset.                     |
| `sr`        | Single- or paired-end short reads | Cutadapt primer trimming and EMU with the minimap2 `sr` preset.                            |

Use `lr:hq` for high-accuracy ONT reads with error levels at or below approximately 1%. See [this preset comparison](https://mbhall88.github.io/post/minimap2-lrhq-preset-testing/) and the [minimap2 release notes](https://github.com/lh3/minimap2/releases/tag/v2.27). Other values, including `map-pb`, are not currently accepted by the workflow.

### Long-read preprocessing

| Parameter                                | Default | Purpose                                                                                              |
| ---------------------------------------- | ------- | ---------------------------------------------------------------------------------------------------- |
| `--adapter_trimming`                     | `false` | Enable ONT adapter trimming with Porechop_ABI.                                                       |
| `--quality_filtering`                    | `true`  | Enable Filtlong length and quality filtering. Set to `false` to bypass filtering.                    |
| `--longread_qc_qualityfilter_minlength`  | `1200`  | Minimum read length in bases.                                                                        |
| `--longread_qc_qualityfilter_maxlength`  | `1800`  | Maximum read length in bases.                                                                        |
| `--longread_qc_qualityfilter_min_mean_q` | `94`    | Filtlong minimum mean quality threshold, expressed as percentage identity rather than a Phred score. |
| `--save_preprocessed_reads`              | `false` | Publish preprocessed reads.                                                                          |

Use `--downsample_n_reads N` to subsample after preprocessing for either long or short reads. Set `--save_downsampled_reads true` to publish the subsampled reads.

### Short-read primer trimming

For short reads, set `--seqtype sr` and provide primers through `--FW_primer` and `--RV_primer`, or through the corresponding samplesheet columns. Nonempty samplesheet values override the parameter for that sample; blank values fall back to the parameter.

```bash
nextflow run main.nf \
  --input /data/shortread_samplesheet.csv \
  --outdir results_shortreads \
  --db /data/emu_database \
  --seqtype sr \
  --FW_primer ACTCCTACGGGAGGCAGCAG \
  --RV_primer GGACTACHVGGGTWTCTAAT \
  --run_krona false \
  -profile docker
```

Replace these example primers with those used for your assay. For already primer-trimmed reads, use `--skip_cutadapt true`. Cutadapt defaults to a minimum overlap of 3 bases (`--cutadapt_min_overlap`) and a maximum error rate of 0.1 (`--cutadapt_max_error_rate`). Use `--retain_untrimmed true` to retain reads without primer matches, and `--save_intermediates true` to publish trimmed reads.

### Results and working files

See the [output documentation](output.md) for the reports and abundance tables produced by the pipeline.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow.log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When running `nextflow run genomic-medicine-sweden/TRANA`, Nextflow downloads and caches the pipeline from GitHub. Update that cached copy with:

```bash
nextflow pull genomic-medicine-sweden/TRANA
```

This does not update a local checkout used with `nextflow run main.nf`; update that checkout separately with Git. For reproducible runs from GitHub, select a release or commit with `-r <revision>`.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command).

#### memory issues

If a process runs out of memory, provide a custom configuration with `-c` and adjust that process's resource request as shown below. Check any `process.resourceLimits` settings in your selected profiles or site configuration, since they can cap CPU, memory, and time requests.

#### Advanced option on process level

The custom config below can be provided to the pipeline via the [`-c`](#-c) parameter.

```nextflow
process {
    withName: 'GMS_TRANA:TRANA:EMU_ABUNDANCE' {
        memory = 50.GB
    }
}
```

> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).
