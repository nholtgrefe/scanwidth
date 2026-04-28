# Documentation

This folder contains a minimal Sphinx documentation site.

## Local build

```bash
pip install -e .[docs]
sphinx-build -b html docs/source docs/build/html
```

Open `docs/build/html/index.html` in a browser.

The docs configuration assumes the package is installed (for example via
`pip install -e .[docs]`) and does not mutate `sys.path` in `conf.py`.

## GitHub Pages

A workflow is provided at `.github/workflows/docs.yml` to deploy docs to
GitHub Pages on pushes to `main`.

In repository settings, ensure Pages is configured to use **GitHub Actions**
as the source.
