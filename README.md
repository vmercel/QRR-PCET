# Path-spread-resolved QTST for proton transfer (QRR)

Manuscript prepared for the *Journal of Chemical Physics* with the official **AIP REVTeX 4.2** class (`aip`, `jcp`).

The AIP sample that this class options come from is `docs/aipsamp.tex` (REVTeX `aipsamp.tex`, AIP Publishing).

```latex
\documentclass[aip, jcp, amsmath, amssymb, reprint]{revtex4-2}
```

For the AIP Peer X-Press review copy, switch `reprint` to `preprint` (12 pt, single column, double spaced).

## Paper

- PDF: [`manuscript/manuscript.pdf`](manuscript/manuscript.pdf)
- Source: [`manuscript/manuscript.tex`](manuscript/manuscript.tex)
- Bibliography: [`references/references.bib`](references/references.bib) (`aipnum4-2.bst`)
- Supplementary: [`supplementary/supplementary.pdf`](supplementary/supplementary.pdf)

## Reproduce

```bash
python -m pip install -r requirements.txt
python run_all.py                          # experiments + figures + LaTeX
python run_all.py --skip_experiments       # figures + compile only
```

Figures are generated from `results/*.json`. Do not hand-edit plot values.

## Zenodo

After the GitHub repository is public:

1. Sign in at [zenodo.org](https://zenodo.org) (GitHub login is easiest).
2. GitHub → Settings → Applications → Zenodo, grant access to this repo.
3. On Zenodo, *Enable* the repository. A new tag/release on GitHub mints a DOI.

`.zenodo.json` in this directory supplies the deposit metadata.

## Citation

See `CITATION.cff`.
