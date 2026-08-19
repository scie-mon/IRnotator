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

Ensure that `$HOME/.local/bin` is on your `PATH`:

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

## Quick start via the controller

The `irnotator` controller is the recommended interface for genome-and-annotation runs. It forwards Nextflow options and IRnotator parameters to `main.nf`, selects the Docker profile unless another profile is specified, and manages pipeline resumption across manual-review rounds.

### Minimal genome-and-annotation run

```bash
irnotator \
  --genome_fasta /path/to/genome.fna \
  --annot_gff '/path/to/annotations/*.gff3' \
  --hmm_dir /path/to/ir_hmms \
  --deeptmhmm_mode container \
  --deeptmhmm_dir /path/to/DeepTMHMM \
  --outdir results
```

Parameters supplied on the command line override values in `nextflow.config`. Always provide explicit `--genome_fasta`, `--annot_gff`, `--hmm_dir`, and `--outdir` values rather than relying on development defaults in the repository configuration.

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

For direct Nextflow execution, use `--review_tsv /path/to/review.tsv` or place the completed file at `<outdir>/manual_review/review.tsv`.

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
| `outdir` | `wrapper_test2` | Output directory |
| `review_tsv` | `null` | Completed review TSV |
| `run_deeptmhmm` | `true` | Run topology prediction |
| `deeptmhmm_mode` | `container` | `container` or `biolib` |
| `deeptmhmm_dir` | repository local path | Local DeepTMHMM distribution |
| `deeptmhmm_model` | `DTU/DeepTMHMM:1.0.24` | BioLib model identifier |
| `deeptmhmm_salvage_paths` | configured list | Reusable results directories |
| `deeptmhmm_keep_going` | `false` | Continue after individual topology failures |
| `deeptmhmm_gpu` | `false` | Request GPU support |
| `normalizer_container` | `quay.io/biocontainers/biopython:1.84` | Normalisation image |
| `hmmer_container` | `biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1` | HMMER image |
| `deeptmhmm_container` | `interpro/deeptmhmm:1.0` | DeepTMHMM image |
| `deeptmhmm_biolib_container` | `null` | Optional BioLib container |

> [!WARNING]
> `hmm_evalue`, `translation_table`, `allow_internal_stops`, and `isoform_overlap_fraction` can affect biological results. Select and report values appropriate to the study.

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
