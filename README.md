# `lifebit-ai/traits`

**Workflow for collecting genetic traits information (heritability, genetic correlation) from GWAS summary statistics**.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/) [![Docker](https://img.shields.io/docker/automated/lifebit-ai/traits.svg)](https://hub.docker.com/r/lifebit-ai/traits)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run lifebit-ai/traits -profile binary_h2
    ```

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run lifebit-ai/traits --input_file '<path to file>' --post_analysis 'heritability' --hapmap3_snplist '<path to file>' --ld_scores_tar_bz2 "<path to tar.bz2 file with LD scores>"
    ```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The lifebit-ai/traits pipeline comes with documentation about the pipeline which you can read at [https://lifebit-ai/traits/docs](https://lifebit-ai/traits/docs) or find in the [`docs/` directory](docs).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

## Credits

`lifebit-ai/traits` was originally written by Marcos Cámara Donoso, Christina Chatzipantsiou, Athanasios Kousathanas.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

> **NOTE**: This pipeline was created using the nf-core template.  For further information or help with nf-core pipelines, you can get in touch with the core developers and community on [Slack](https://nfcore.slack.com/channels/lifebit-ai/traits) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  lifebit-ai/traits for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
