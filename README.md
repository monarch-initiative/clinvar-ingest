# clinvar-ingest

| [Documentation](https://monarch-initiative.github.io/clinvar-ingest) |

modular ingest for clinvar variants. Makes variant nodes and then creates associations based on the relative information (phenotype, disease, gene, pathogenicity)

## Ingest files and how nodes and edges are generated
Two files downloaded from clinvar are leveraged in this ingest
- clinvar.vcf which contains a single line per clinvar variant with each variant's associated terms reported in the INFO column and grouped by which submission record(s) they originated from. https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38  (hg38 genome version)
- submission_summary.txt contains a single line per variant record. These records contain in depth information about the variant in question that we can leverage in the ingest process. Multiple records often exist per one variant. https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

Variant nodes (SequenceVariant)
- SequenceVariant nodes are created from clinvar variants that are deemed Pathogenic or Likely Pathogenic. Additionally, only variants that have a clinvar review status of 3 or more stars (4 maximum) will be included. This subset corresponds to the most credible set of variants that are pathogenic within the clinvar dataset.

Variant → Disease edges (VariantToDiseaseAssociation)
 - Disease ids are derived from the ReportedPhenotypeInfo column within the submission_summary.txt file. This column consists of medgen ids that we then map to a mondo id. Alternatively, if a mondo id cannot be found, then the SubmittedPhenotypeInfo column will be used instead. If neither column maps to a mondo id then no edge will be made.
 - Predicates are derived from the ClinicalSignificance column within the submission_summary.txt file. Currently only Pathogenic and Likely pathogenic are included as “causes” and “associated_with_increased_likelihood_of” respectively.

Variant → Phenotype edges (VariantToPhenotypicFeatureAssociation)
- These edges are only created if a Variant → Disease edge can be made. The phenotype terms themselves are derived from the INFO column of the clinvar.vcf file. The information within this column is reported as groups of terms that map back to an individual record(s) from the submission_summary.txt file. Any group of terms that contains the disease id from the Variant → Disease edge will have Variant → Phenotype edges created for all reported Human phenotype ontology terms reported within the group.
- The predicate for these edges in "contributes_to"

Variant → Gene edges (VariantToGeneAssociation)
- These edges are created only if a Variant → Disease edge can be made. Gene symbols are derived from the INFO column within the clinvar.vcf file and gene symbols are mapped to ncbi genes.
- The predicate is_sequence_variant_of is used. Sequence ontology terms are also reported within the INFO column pertaining to the variants "molecular consequence" (MC subfield within INFO). These terms are recored in the "type" slot for the SequenceVariant node that is created (https://biolink.github.io/biolink-model/type/)



## Requirements

- Python >= 3.10
- [Poetry](https://python-poetry.org/docs/#installation)

## Setting Up a New Project

Upon creating a new project from the `cookiecutter-monarch-ingest` template, you can install and test the project:

```bash
cd clinvar-ingest
make install
make test
```

There are a few additional steps to complete before the project is ready for use.

#### GitHub Repository

1. Create a new repository on GitHub.
1. Enable GitHub Actions to read and write to the repository (required to deploy the project to GitHub Pages).
   - in GitHub, go to Settings -> Action -> General -> Workflow permissions and choose read and write permissions
1. Initialize the local repository and push the code to GitHub. For example:

   ```bash
   cd clinvar-ingest
   git init
   git remote add origin https://github.com/<username>/<repository>.git
   git add -A && git commit -m "Initial commit"
   git push -u origin main
   ```

#### Transform Code and Configuration

1. Edit the `download.yaml`, `transform.py`, `transform.yaml`, and `metadata.yaml` files to suit your needs.
   - For more information, see the [Koza documentation](https://koza.monarchinitiative.org) and [kghub-downloader](https://github.com/monarch-initiative/kghub-downloader).
1. Add any additional dependencies to the `pyproject.toml` file.
1. Adjust the contents of the `tests` directory to test the functionality of your transform.

#### Documentation

1. Update this `README.md` file with any additional information about the project.
1. Add any appropriate documentation to the `docs` directory.

> **Note:** After the GitHub Actions for deploying documentation runs, the documentation will be automatically deployed to GitHub Pages.  
> However, you will need to go to the repository settings and set the GitHub Pages source to the `gh-pages` branch, using the `/docs` directory.

#### GitHub Actions

This project is set up with several GitHub Actions workflows.  
You should not need to modify these workflows unless you want to change the behavior.  
The workflows are located in the `.github/workflows` directory:

- `test.yaml`: Run the pytest suite.
- `create-release.yaml`: Create a new release once a week, or manually.
- `deploy-docs.yaml`: Deploy the documentation to GitHub Pages (on pushes to main).
- `update-docs.yaml`: After a release, update the documentation with node/edge reports.


Once you have completed these steps, you can remove this section from the `README.md` file.

## Installation

```bash
cd clinvar-ingest
make install
# or
poetry install
```

> **Note** that the `make install` command is just a convenience wrapper around `poetry install`.

Once installed, you can check that everything is working as expected:

```bash
# Run the pytest suite
make test
# Download the data and run the Koza transform
make download
make run
```

## Usage

This project is set up with a Makefile for common tasks.  
To see available options:

```bash
make help
```

### Download and Transform

Download the data for the clinvar_ingest transform:

```bash
poetry run clinvar_ingest download
```

To run the Koza transform for clinvar-ingest:

```bash
poetry run clinvar_ingest transform
```

To see available options:

```bash
poetry run clinvar_ingest download --help
# or
poetry run clinvar_ingest transform --help
```

### Testing

To run the test suite:

```bash
make test
```

---

> This project was generated using [monarch-initiative/cookiecutter-monarch-ingest](https://github.com/monarch-initiative/cookiecutter-monarch-ingest).  
> Keep this project up to date using cruft by occasionally running in the project directory:
>
> ```bash
> cruft update
> ```
>
> For more information, see the [cruft documentation](https://cruft.github.io/cruft/#updating-a-project)
