# Third-party notices

IRnotator is distributed under the GNU General Public License, version 3.0 or later (`GPL-3.0-or-later`). See [LICENSE](LICENSE).

IRnotator orchestrates third-party software and may distribute or use third-party reference data. The inclusion or use of a component listed below does not change that component's licence, copyright, attribution requirements, or terms of use. Users are responsible for complying with the applicable terms of every dependency and reference resource used in their runs.

## Bundled reference data

| Component | Role in IRnotator | Distributed in this repository? | Licence / terms | Provenance |
|---|---|---:|---|---|
| Pfam profile HMM `PF00060` | Reference profile HMM used for IR candidate identification | Yes | CC0 1.0 | Obtain the exact Pfam release/version and source URL from the file-acquisition record before making a release. |
| Pfam profile HMM `PF10613` | Reference profile HMM used for IR candidate identification | Yes | CC0 1.0 | Obtain the exact Pfam release/version and source URL from the file-acquisition record before making a release. |

Pfam is freely available under the Creative Commons Zero 1.0 Universal (`CC0 1.0`) dedication. Where practical, retain the original source, release version, accession, and retrieval date for each bundled HMM.

## Runtime dependencies

| Component | Role in IRnotator | Distributed in this repository? | Licence / terms |
|---|---|---:|---|
| Nextflow | Workflow runtime | No | Apache License 2.0. See the [Nextflow licence information](https://www.nextflow.io/about-us.html). |
| HMMER | Profile-HMM search runtime | No | External dependency. Users must obtain and use HMMER under its applicable licence and distribution terms. |
| DeepTMHMM | Protein topology prediction | No software, model weights, or service credentials are distributed by IRnotator | External dependency. Users must obtain access to DeepTMHMM and comply with the applicable DeepTMHMM and/or BioLib terms. IRnotator does not grant any rights to DeepTMHMM software, model weights, hosted services, or output quotas. |

## Containers and execution environments

IRnotator may execute tools through OCI/Docker-compatible containers or Apptainer/Singularity-compatible container images. Container images are not ordinarily redistributed as part of this repository unless explicitly stated otherwise.

Each image, its operating-system layers, and the software installed within it may have separate copyright notices and licence obligations. Users and deployers are responsible for recording the image name, immutable digest or version tag, registry source, and applicable licence information for the images used in a run.

## Generated files and user data

Genome assemblies, annotations, protein sequences, run outputs, DeepTMHMM result caches, and other user-provided or generated data are not relicensed by IRnotator. Their ownership, access restrictions, redistribution rights, and citation requirements remain the responsibility of the data provider and user.

## Updating these notices

Update this document before each public release when dependencies, container images, bundled reference data, or their versions change. In particular, replace the Pfam provenance placeholders above with the exact release/version, retrieval date, and source URL used to create the distributed HMM files.
