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

Key result files:

- `results/harmonic_benchmark.json` — exact and PIMC $\langle\chi_p\rangle$
- `results/eckart_rates.json` — analytic Eckart rates and KIEs
- `results/eckart_pimc.json` — umbrella-PIMC centroid QTST (QRR at $\kappa=1$)
- `results/adw_landscapes.json` — ADW WHAM PMFs and $p^\ddagger(\chi_p)$

Deuterium centroid QTST at $V_0=12\,k_{\mathrm{B}}T$ is flagged `reliable: false` and is omitted from the rate comparison.

## Archive

This version is archived on Zenodo:

- Record: https://zenodo.org/records/22718950
- Version DOI: [https://doi.org/10.5281/zenodo.22718950](https://doi.org/10.5281/zenodo.22718950)
- Concept DOI (always the latest): [https://doi.org/10.5281/zenodo.22718949](https://doi.org/10.5281/zenodo.22718949)
- GitHub tag: [`v1.0.1`](https://github.com/vmercel/QRR-PCET/tree/v1.0.1)

## Citation

Vubangsi, M., Al-Turjman, F., & Tchoffo, M. (2026). *Path-spread-resolved quantum transition-state theory for proton transfer* (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.22718950

See also `CITATION.cff`.
