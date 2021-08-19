# `lifebit-ai/traits`

**Workflow for collecting genetic traits information (heritability, genetic correlation) from GWAS summary statistics**.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/) [![Docker](https://img.shields.io/docker/automated/lifebit-ai/traits.svg)](https://hub.docker.com/r/lifebit-ai/traits)


## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run main.nf -profile binary_h2
    ```

4. Start running your own analysis!

    ```bash
    nextflow run main.nf --input_gwas_statistics '<path to file>' --post_analysis 'heritability' --hapmap3_snplist '<path to file>' --ld_scores_tar_bz2 "<path to tar.bz2 file with LD scores>"
    ```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The lifebit-ai/traits pipeline comes with documentation about the pipeline which you can read at [https://lifebit-ai/traits/docs](https://lifebit-ai/traits/docs) or find in the [`docs/` directory](docs).


## Credits

`lifebit-ai/traits` was originally written by Marcos Cámara Donoso, Christina Chatzipantsiou, Athanasios Kousathanas.

