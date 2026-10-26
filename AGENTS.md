# Instructions for AI Agents

This repository contains the Quarto book **Mathematical & Statistical Literature Review of Multispecies Stock Assessment Models**, configured in `_quarto.yml` and rendered to HTML in `docs/`.

## Reference Papers & Companion Private Repository

Because several cited publications are not open access, downloaded PDFs and reference software source snippets must **never** be committed to this public repository.

- **Local path**: The downloaded PDFs and extracted software source files (`software_docs/`) live in the companion private repository located at the sibling directory `../multispecies_stock_assessment_papers` (`/home/gustav/Git/multispecies_stock_assessment_papers`).
- **GitHub repository**: `gustavdelius/multispecies_stock_assessment_papers` (private). If the sibling directory is not present locally, clone it with:
  ```bash
  gh repo clone gustavdelius/multispecies_stock_assessment_papers ../multispecies_stock_assessment_papers
  ```
- When verifying mathematical formulations, equations, likelihoods, or software implementation details, read the primary PDFs and `software_docs/` files in `../multispecies_stock_assessment_papers`.

## Citations & Bibliography

- All bibliography entries and DOIs are maintained in `references.bib` and rendered in `references.qmd`.
- Always use bracketed Quarto citation syntax (`[@key]` or `[@key1; @key2]`) for all citations rather than unbracketed `@key` or linking to local PDF files.

## Building the Book

After editing any `.qmd` chapter or `_quarto.yml`, render the HTML book into `docs/`:

```bash
quarto render
```
