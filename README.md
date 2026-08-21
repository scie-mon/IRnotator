# IRnotator

IRnotator is a Nextflow DSL2 workflow for the identification, topology-based classification, manual curation, and consolidation of ionotropic receptor (IR) gene annotations.

The recommended workflow accepts a genome FASTA together with one or more GFF annotation files. It derives normalised protein sequences, detects candidate IR proteins with HMMER, obtains DeepTMHMM topology predictions, supports manual review, and produces curated candidate sets and merged IR isoform annotations.

## Contents

- [Installation and prerequisites](#installation-and-prerequisites)
- [Quick start via the controller](#quick-start-via-the-controller)
- [Direct Nextflow usage](#direct-nextflow-usage)
- [Input routes and DeepTMHMM backends](#input-routes-and-deeptmhmm-backends)
- [Manual-review cycle](#manual-review-cycle)
- [Parameters](#parameters)
- [Outputs, reruns, and result salvage](#outputs-reruns-and-result-salvage)
- [HPC and container configuration](#hpc-and-container-configuration)

## Installation and prerequisites

### Requirements

- [Nextflow](https://www.nextflow.io/)
- Java, as required by the installed Nextflow release
- Git
- Docker, or Apptainer/Singularity on HPC systems
- A local DeepTMHMM distribution only when using `--deeptmhmm_mode container`
- A BioLib Python environment and valid BioLib credentials only when using `--deeptmhmm_mode biolib`

IRnotator runs input normalisation, HMMER, and container-mode DeepTMHMM in containers. The selected container runtime must be able to retrieve or access the configured images.

> [!NOTE]
> Local DeepTMHMM installation requires a separately obtained licensed distribution. IRnotator can instead reuse externally generated DeepTMHMM results; see [Outputs, reruns, and result salvage](#outputs-reruns-and-result-salvage).

### Obtain IRnotator

```bash
git clone git@github.com:scie-mon/IRnotator.git
cd IRnotator
```

### Install the controller command

IRnotator provides `bin/irnotator.py` as its recommended controller. Install a user-local symlink so that it can be called as `irnotator` from any directory:

```bash
mkdir -p "$HOME/.local/bin"
chmod +x bin/irnotator.py
ln -sfn "$(pwd)/bin/irnotator.py" "$HOME/.local/bin/irnotator"
```

Ensure that `$HOME/.local/bin` is on your `PATH`. If not:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Add this `export` command to `~/.bashrc` or `~/.zshrc` to make the setting persistent. Verify the installation:

```bash
irnotator --help
```

For a shared HPC installation, use a group-managed directory already on users' `PATH` rather than `$HOME/.local/bin`.

### Local DeepTMHMM setup

IRnotator supports local DeepTMHMM execution through `--deeptmhmm_mode container`. The licensed DeepTMHMM distribution is not bundled with IRnotator.

1. Obtain the DeepTMHMM standalone package through the [official request form](https://forms.cloud.microsoft/Pages/ResponsePage.aspx?id=rURVGECW5U69G3ew50KVeHoV9Bi7sdNBnLmAyRfuhiJUM1lRRVY0MExIMEU1V0UxN1ZDVkFFUzlVQi4u).

2. Extract the received archive into a persistent location readable from every execution node:

   ```bash
   mkdir -p /path/to/software
   unzip DeepTMHMM-v1.0.zip -d /path/to/software
   ```

3. Confirm that the extracted directory contains either `predict.py` directly, or `deeptmhmm/predict.py`. These are the two layouts supported by IRnotator's local runner.

4. Configure the path in `IRnotator/nextflow.config` or a site-specific Nextflow config:

   ```groovy
   params {
       deeptmhmm_mode = 'container'
       deeptmhmm_dir  = '/absolute/path/to/DeepTMHMM'
   }
   ```

5. Supply a configuration file at launch time if not using `IRnotator/nextflow.config`:

   ```bash
   irnotator -c local_deeptmhmm.config ...
   ```

   Alternatively, supply the path directly:

   ```bash
   irnotator --deeptmhmm_dir '/absolute/path/to/DeepTMHMM' ...
   ```

The default configuration uses the `interpro/deeptmhmm:1.0` container image. With Docker, the image is normally pulled automatically on first use. It can also be pre-fetched with:

```bash
docker pull interpro/deeptmhmm:1.0
```

> [!IMPORTANT]
> When using an HPC executor, `deeptmhmm_dir` must be readable from compute nodes. Use a shared filesystem path or configure the required filesystem bind mounts in your site-specific Nextflow configuration.

### GPU support

IRnotator can request GPU access for local container-mode DeepTMHMM with:

```bash
irnotator \
  --deeptmhmm_mode container \
  --deeptmhmm_gpu true \
  ...
```

For Docker, `--deeptmhmm_gpu true` adds the container option `--gpus all`.
The host must have a supported NVIDIA GPU setup, compatible drivers, and a
container runtime configured for GPU access. On HPC systems, also request a GPU
through the scheduler in the site-specific Nextflow configuration.

> [!NOTE]
> IRnotator runs DeepTMHMM independently for each candidate sequence. GPU
> acceleration therefore does not necessarily improve total runtime: per-job
> container startup, model loading, filesystem I/O, and scheduler overhead can
> dominate for many short single-sequence tasks.

A capable multicore CPU allocation can be as fast as, or faster than, GPU mode
for this execution pattern.

Use GPU mode when local benchmarks demonstrate a meaningful improvement for the
target dataset, or when CPU capacity is constrained and compatible GPU
resources are readily available.

### Run the bundled test dataset

IRnotator includes a small genome-and-annotation test dataset in
`test/multi_gff/`. The repository's `nextflow.config` is pre-configured to use
this dataset, allowing an installation test without providing external input
files.

Run the test through the controller:

```bash
irnotator --outdir test_out
```

Or run the workflow directly with Nextflow:

```bash
nextflow run main.nf -profile docker --outdir test_out
```

#### Run the test without DeepTMHMM

The test dataset also includes precomputed topology predictions in:

```text
test/multi_gff/deeptmhmm_predicts/
```

These predictions are **not used automatically**. If no local DeepTMHMM
installation or BioLib backend is available, configure this directory as a
salvage location:

```groovy
params {
    deeptmhmm_salvage_paths = [
        "${projectDir}/test/multi_gff/deeptmhmm_predicts"
    ]
}
```

> [!NOTE]
> `deeptmhmm_predicts/` must be supplied explicitly because it is not part of
> the default salvage-path configuration.

## Inputs

IRnotator requires one protein-source input mode, an HMM profile directory, and
optionally one or more directories containing reusable DeepTMHMM results.

### Protein-source input modes

Choose **exactly one** of the following modes.

| Mode | Required parameters | Description |
|---|---|---|
| Genome + annotation | `--genome_fasta`, `--annot_gff` | Recommended mode. IRnotator derives and normalises protein sequences from a genome FASTA and one or more annotation GFF files. This mode supports downstream IR isoform merging. |
| Protein FASTA | `--proteins_faa` | Alternative mode for an already available protein dataset. Genome/GFF-specific normalisation and final GFF merging are not available in this mode. |

For genome + annotation mode, `--annot_gff` may specify one or more GFF files,
for example:

```bash
--genome_fasta genome.fna \
--annot_gff 'annotations/*.gff3'
```

The GFF protein identifier attribute is controlled by
`--gff_protein_attribute` and defaults to `protein_id`.

### IR HMM profiles

IRnotator searches protein sequences with HMMER against all HMM profiles in
`--hmm_dir`.

The repository provides the canonical IR profile set in `IRnotator/hmms/` which are used by default:

```text
hmms/
├── PF00060.hmm
└── PF10613.hmm
```

Alternatively, provide a directory containing custom HMM profile files

```bash
--hmm_dir /path/to/custom_hmms
```

 or add custom HMMs to `hmms/`.

### DeepTMHMM result cache and salvage paths

IRnotator maintains `deeptmhmm_dump/` as its local DeepTMHMM result cache. At
the start of each run, the pipeline sanitizes this directory and first attempts
to reuse compatible cached results from it.

All valid DeepTMHMM results—whether generated during the current run or
salvaged from an external location—are added to `deeptmhmm_dump/`. This allows
subsequent runs to reuse previously obtained predictions before submitting new
DeepTMHMM jobs.

Additional external result locations can be supplied with
`deeptmhmm_salvage_paths`:

```groovy
params {
    deeptmhmm_salvage_paths = [
        "${projectDir}/deeptmhmm_dump",
        "/path/to/previous_deeptmhmm_results",
        "/path/to/additional_deeptmhmm_results"
    ]
}
```

Keep `${projectDir}/deeptmhmm_dump` as the first entry in
`deeptmhmm_salvage_paths` to preserve the default cache-first behaviour. If
`deeptmhmm_salvage_paths` is overridden without this entry, IRnotator will not
search its local `deeptmhmm_dump/` cache before checking external salvage
locations or running new DeepTMHMM jobs.

Only sequences for which no compatible cached or salvaged result is found are
submitted to DeepTMHMM when `--run_deeptmhmm true`. To prevent new topology
runs entirely, set:

```bash
--run_deeptmhmm false
```

If no local DeepTMHMM installation, BioLib execution, or compatible salvaged
result is available, IRnotator writes the affected candidate sequences to:

```text
<outdir>/deeptmhmm/unresolved_deeptmhmm.faa
```

Run DeepTMHMM externally on this FASTA file and add the resulting valid
per-sequence result directories to `deeptmhmm_dump/` or another configured
salvage path. Then rerun IRnotator with the same inputs and resume state. The
pipeline will salvage the newly available results and continue from the
topology-processing stage rather than reprocessing completed upstream tasks.


## Quick start via the controller

The `irnotator` controller is the recommended interface for
genome-and-annotation runs. It forwards Nextflow options and IRnotator
parameters to `main.nf`, selects the Docker profile unless another profile is
specified, and manages pipeline resumption across manual-review rounds.

After each workflow round that requires curation, the controller serves the
local manual-review application from `<outdir>/manual_review/`. It manages the
review session, accepts browser-submitted or externally exported `review.tsv`
decisions, validates and records the review state, and launches subsequent
Nextflow rounds until the review workflow is complete.

### Minimal genome-and-annotation run

```bash
irnotator \
  --genome_fasta /path/to/genome.fna \
  --annot_gff '/path/to/annotations/*.gff3' \
  --deeptmhmm_mode container \
  --deeptmhmm_dir /path/to/DeepTMHMM \
  --outdir results
```

Parameters supplied on the command line override values in `nextflow.config`.

## Direct Nextflow usage

```bash
nextflow run main.nf \
  -profile docker \
  -c irnotator.local.config \
  --genome_fasta /path/to/genome.fna \
  --annot_gff '/path/to/annotations/*.gff3' \
  --hmm_dir /path/to/ir_hmms \
  --outdir results \
  -resume
```

Unlike the controller, direct Nextflow execution leaves `-resume` under user control. Direct execution prepares the same manual-review files, but does not serve the local review interface or orchestrate successive review rounds. Supply `--review_tsv` explicitly, or place it at `<outdir>/manual_review/review.tsv` before a resumed run.

## Input routes and DeepTMHMM backends

| Route | Required parameters | Use case |
|---|---|---|
| Genome + annotation | `--genome_fasta`, `--annot_gff` | Recommended route; supports final IR isoform merging. |
| Protein FASTA | `--proteins_faa` | Alternative route for an existing protein set. |

### DeepTMHMM backends

| Backend | Set with | Requirements |
|---|---|---|
| Local container | `--deeptmhmm_mode container` | Licensed standalone distribution and `--deeptmhmm_dir` |
| BioLib | `--deeptmhmm_mode biolib` | BioLib Python environment and credentials |
| Reuse existing results | `--run_deeptmhmm false` with `--deeptmhmm_salvage_paths` | Compatible completed DeepTMHMM results |

## Manual-review cycle

1. Run IRnotator through `irnotator`.
2. Review the package under `<outdir>/manual_review/` in the controller-served browser interface.
3. Submit or export `review.tsv` as instructed by the controller.
4. The controller applies decisions and resumes the workflow for the next review round.

For direct Nextflow execution, use `--review_tsv /path/to/review.tsv` or place the completed file at `<outdir>/manual_review/review.tsv`. Then, re-run the pipeline with `-resume`.

## Parameters

CLI `--parameter value` settings override `nextflow.config` values.

| Parameter | Default | Description |
|---|---:|---|
| `proteins_faa` | `null` | Input protein FASTA |
| `genome_fasta` | repository test path | Genome FASTA |
| `annot_gff` | repository test glob | One or more GFF annotation files |
| `gff_protein_attribute` | `protein_id` | GFF protein attribute |
| `translation_table` | `1` | Translation table |
| `allow_internal_stops` | `false` | Permit internal stops during normalisation |
| `isoform_overlap_fraction` | `0.60` | Isoform-merge overlap fraction |
| `hmm_dir` | `test/hmms` | IR HMM directory |
| `hmm_evalue` | `1e-5` | HMMER E-value threshold |
| `outdir` | `null` | Output directory |
| `review_tsv` | `null` | Completed review TSV |
| `run_deeptmhmm` | `true` | Run topology prediction |
| `deeptmhmm_mode` | `container` | `container` or `biolib` |
| `deeptmhmm_dir` | `null` | Local DeepTMHMM distribution |
| `deeptmhmm_model` | `DTU/DeepTMHMM:1.0.24` | BioLib application identifier passed to the cloud runner; used only when `--deeptmhmm_mode biolib`. |
| `deeptmhmm_salvage_paths` | `["${projectDir}/deeptmhmm_dump"]` | Reusable results directories |
| `deeptmhmm_keep_going` | `false` | Continue after individual topology failures |
| `deeptmhmm_gpu` | `false` | Request GPU support when `--deeptmhmm_mode container`  |
| `normalizer_container` | `quay.io/biocontainers/biopython:1.84` | Normalisation image |
| `hmmer_container` | `biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1` | HMMER image |
| `deeptmhmm_container` | `interpro/deeptmhmm:1.0` | DeepTMHMM image |
| `deeptmhmm_biolib_container` | `null` | Optional BioLib container |

## Outputs, reruns, and result salvage

Results are published under `--outdir`:

| Location | Contents |
|---|---|
| `input/` | Normalised input records |
| `hmmsearch/` | Per-profile HMMER outputs |
| `hmm_hits/` | Collected hits and candidate sequences |
| `deeptmhmm/` | Finalised topology results and unresolved sequences |
| `classification/` | Candidate summaries, classifications, decisions, and passed-protein FASTA |
| `manual_review/` | Review package and controller files |
| `merged_gff/` | Merged IR annotations for genome + GFF input |
| `echo/` | Echoed input records, when produced |

With the controller, rerun the same command as instructed by the controller. With direct Nextflow, rerun with `-resume`. Retain the work directory to reuse cache entries. Unresolved topology sequences are written to `<outdir>/deeptmhmm/unresolved_deeptmhmm.faa`.

## HPC and container configuration

Configure cluster-specific executor settings in a site config rather than editing the repository configuration:

```groovy
profiles {
  hpc {
    process.executor = 'slurm'
    process.queue = 'standard'
    apptainer.enabled = true
    docker.enabled = false
    apptainer.autoMounts = true
  }
}

params {
  hmm_dir = '/shared/reference/irnotator/hmms'
  deeptmhmm_mode = 'container'
  deeptmhmm_dir = '/shared/software/DeepTMHMM'
}
```

Adapt scheduler directives, container cache, filesystem binds, and resource limits to local policy. Inputs, work directory, output directory, HMM profiles, and `deeptmhmm_dir` must be readable from compute nodes.
