# STAB Knowledge Base
This is a knowledge base for the Space-Time Analysis Bayes (STAB) group led by Dr. Jon Wakefield at UW/UBC and Dr. Zehang Richard Li at UCSC.

## Site structure

This site is built with [Jekyll](https://jekyllrb.com) and the
[Just the Docs](https://just-the-docs.com) theme, which gives a textbook-style
layout with a collapsible sidebar of chapters/subtopics. Content lives in
plain Markdown files — no build tooling is required to add or edit a page.

Current outline (chapters map to top-level folders/files, subtopics to files
inside them):

1. **Survey Theory** — finite population sampling, randomness, weights,
   direct estimates, variance fix
2. **Survey Setup & Notation** — DHS, MICS, other
3. **SAE & Bayesian Background**
4. **INLA & STAN**
5. **Models** — area-level, unit-level, cluster-level (each with theory,
   when we use them, limitations, code)
6. **Appendix** — group members, dictionary, resources

To add content, just edit the relevant `.md` file — the front matter
(`title`, `parent`, `nav_order`, etc.) controls where it sits in the sidebar.
To add a brand-new page, copy an existing `.md` file's front matter as a
template.

## Publishing on GitHub Pages

1. Push this repo to GitHub (as `stab-knowledge-base` under this org/user,
   to match `baseurl` in [`_config.yml`](_config.yml) — update `baseurl` if
   the repo name or hosting path differs).
2. In the repo's **Settings → Pages**, set **Source** to "Deploy from a
   branch", branch `main`, folder `/ (root)`.
3. GitHub will build the site automatically on every push; it's typically
   live at `https://<org-or-user>.github.io/stab-knowledge-base/` within a
   few minutes.

## Local preview (optional)

With Ruby installed:

```bash
bundle install
bundle exec jekyll serve
```

Then open `http://localhost:4000`.
