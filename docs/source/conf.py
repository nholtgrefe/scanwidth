"""Sphinx configuration for scanwidth documentation."""

from __future__ import annotations

project = "scanwidth"
author = "Niels Holtgrefe"
copyright = "2026, Niels Holtgrefe"
release = "0.2.0.dev"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": False,
}

templates_path = ["_templates"]
exclude_patterns: list[str] = []

html_theme = "pydata_sphinx_theme"
html_title = "scanwidth"
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/nholtgrefe/scanwidth",
            "icon": "fa-brands fa-github",
        },
    ],
}
html_static_path = ["_static"]
