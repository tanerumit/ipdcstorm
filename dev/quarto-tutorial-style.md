---
name: quarto-tutorial-style
description: Authoring conventions for Quarto (.qmd) tutorial vignettes in this project. Captures YAML front matter, section numbering, code/output presentation, callouts, cross-references, and reference-list style. Apply when creating a new tutorial, promoting an .md narrative into a .qmd vignette, or auditing existing tutorials for consistency.
version: 1.0.0
source: vignettes/tutorial_1_setup.qmd
---

# Quarto tutorial style

Machine-readable style reference. Each section defines **what** the pattern is and **how** to apply it. Anchor patterns are given as literal snippets; substitute placeholders spelled `<like-this>`.

## 1. When to apply

Apply this style when:

- writing a new `.qmd` vignette under `vignettes/`;
- translating a standalone narrative `.md` file into a vignette (preserve body text verbatim, add the YAML block below, adapt code fences);
- auditing an existing vignette before a docs release.

Do **not** apply to:

- README source (`README.Rmd` uses `github_document` output, not the HTML format below);
- pkgdown reference pages (auto-generated from roxygen).

## 2. YAML front matter

Use this exact block as the starting point for every tutorial. Change only `title` (and `bibliography:` when citations are needed).

```yaml
---
title: "<Tutorial N> - <Short Title>"
subtitle: "Hurricane Wind Hazard Model for the Caribbean"
author: ""
date: today
format:
  html:
    theme: cosmo
    toc: true
    toc-depth: 3
    toc-location: left
    number-sections: false
    code-fold: false
    code-tools: true
    code-copy: true
    highlight-style: github
    smooth-scroll: true
    self-contained: true
execute:
  eval: false
  echo: true
  warning: false
  message: false
editor:
  markdown:
    wrap: sentence
---
```

Rules for each field:

- `title`: sentence case, `-` (hyphen-space-hyphen) between tutorial label and short title. Example: `"Tutorial 4 - Climate Stress-Test Experiment"`.
- `subtitle`: fixed project subtitle. Do not change between tutorials.
- `author`: intentionally empty string. A CI/CD step can inject Git commit authors later.
- `date: today`: render-time date; avoids stale timestamps.
- `number-sections: false`: section numbers are **written explicitly in the heading text** (see §3). Turning this on would cause double numbering.
- `execute.eval: false`: static code. See §4.
- `wrap: sentence`: one sentence per source line. Keeps Git diffs small.

When references are needed, add:

```yaml
bibliography: references.bib
link-citations: true
```

and create / extend `vignettes/references.bib`. Keys use the `surnameYEARshortword` pattern (e.g., `holland1980`, `knaff2015`).

## 3. Heading hierarchy and numbering

Use one `# H1` for the tutorial-scope page title (mirrors the YAML `title`, can be phrased differently).

Subsequent headings are **manually numbered** at every level, with a period after each numeric component. Depth is capped at H4.

```markdown
# Tutorial N - <Title>

## 1. <Section>

### 1.1. <Subsection>

#### 1.1.1. <Sub-subsection>

## 2. <Section>
```

Rules:

- Always put a period after every numeric component: `1.`, `2.1.`, `4.2.5.` — never `1`, `2.1`, `4.2.5`.
- One space after the final period before the title text.
- Title casing matches the sibling tutorials (sentence case with proper nouns capitalised).
- Never rely on Quarto auto-numbering. `number-sections: false` is mandatory when this style is applied.
- Reserve `## 8. References` (or the next unused top-level number) for the reference list — see §9.

## 4. Code blocks

Two distinct rendering intents. Pick per-chunk.

### 4.1. Static R (most common)

For R code that the tutorial wants to **show** but not execute at render time — lets you embed pre-computed output without CI/environment flakiness:

````markdown
``` r
out_base <- run_hazard_model(cfg, targets, seed = 42L)
```

    #> Hazard configuration (preset: "default")
    #>   Study period  : 1970 - present
    #>   Simulation    : 2,000 synthetic years
````

- Opening fence is ``` ```r ``` — no curly braces (static, not a chunk).
- Output block directly under the code block uses **4-space indentation** and begins each line with `#> ` (hash, greater-than, space). This matches `knitr::opts_chunk$set(comment = "#>")`.
- No blank line between the code fence and the indented output.
- Put at least one blank line between the output block and the next paragraph.

### 4.2. Executable R (rare — only for chunks that must run)

When the code must run at render time (live plots, auto-updating tables):

````markdown
```{r}
#| label: <unique-label>
#| fig-width: 9
#| fig-height: 6
<code>
```
````

- Every executable chunk has a unique `#| label:` for the Quarto execution cache.
- Override `eval: true` at the chunk level only when needed.

### 4.3. Non-R code

Other languages use the bare language tag, always static:

````markdown
``` bash
devtools::install_github("tanerumit/ipdcstorm")
```
````

## 5. Callouts (blockquote style)

Use plain markdown blockquotes for tips, warnings, and inspection prompts. **Do not** use Quarto's `::: callout-*` fences in tutorials covered by this style — the blockquote idiom is consistent across Tutorial 1's body.

```markdown
> **<Short title>**
>
> <Body paragraph. Kept to one sentence per line per YAML wrap rule.>
```

Rules:

- Bold-first-line, then a blank `>` line, then body.
- Body may contain nested markdown (tables, lists, inline code).
- Keep callouts short; if the content exceeds 8 lines, promote it to a numbered subsection.

Exception: in tutorials that *do* use callouts (e.g., Tutorial 4's `::: callout-note`), follow Quarto's syntax there. But Tutorial 1 defines the default.

## 6. Tables

Use pipe tables. Column headers separated by pipes, alignment row with dashes.

```markdown
| Column | Description |
|----|----|
| `name` | Short identifier. |
| `lambda` | Annual rate, \(1 / \mathrm{yr}\). |
```

Rules:

- Left alignment everywhere unless numeric columns need right alignment — then use `|---:|` for the numeric columns.
- Header cells in title case; body cells in sentence case ending with a period.
- Inline code (`` `colname` ``) for column / variable / function names.
- Tables embedded inside a blockquote callout are supported (preserve the leading `> ` on every line).

## 7. Cross-references to the package

Link to pkgdown reference pages via:

```markdown
[`function_name()`](https://tanerumit.github.io/ipdcstorm/reference/function_name.md)
```

Rules:

- Link text is inline-code with parentheses — `` `function_name()` `` — so readers see it as an R function call.
- URL uses the pkgdown site domain and the `.md` suffix (pkgdown rewrites `.md` → `.html` automatically).
- Use this form throughout prose; do not link the same function more than twice per section.

Cross-references within the tutorial use explicit section numbers in text plus a `#sec-<slug>` anchor when needed:

```markdown
See [Section 5.3](#sec-background-wind) for the full workflow.
```

Define the anchor at the target heading:

```markdown
### 5.3. Background Wind {#sec-background-wind}
```

## 8. Inline notation

- **Code identifiers**: inline backticks — `` `N_SIM` ``, `` `out_base$rates` ``, `` `"stationary"` ``.
- **Units**: non-breaking or literal spaces are acceptable. Use `kt` (knots), `km`, `m/s`, `hPa`. Don't mix unit systems within one table.
- **Inequality / ranges**: Unicode ≥, ≤, — (em-dash), – (en-dash) for number ranges. Ascii fallback (`>=`, `<=`, `-`) is allowed but be consistent within a file.
- **Math**: Pandoc math delimiters `$...$` for inline, `$$...$$` for display. Keep complex formulas out of tables; describe them in running prose instead.

## 9. References section

When the tutorial cites named work in the body, close with a References section. Use APA 7 formatting manually — do **not** rely on Pandoc citation rendering unless `bibliography:` is declared in YAML (§2).

Template:

```markdown
## <N>. References

<Author, A. A., Author, B. B., & Author, C. C.> (<Year>). <Article title in sentence case>. *<Journal Name>*, *<Volume>*(<Issue>), <pp1>–<pp2>. https://doi.org/<doi>

<Author, A. A.> (<Year>). *<Book or dataset or report title in title case>* (<edition / number>) [<Descriptor, e.g., Data set>]. <Publisher>. https://<url>
```

Rules:

- Section heading uses the next unused top-level number (e.g., `## 8. References` in Tutorial 1, `## 9. References` if a tutorial has an earlier section).
- Entries ordered alphabetically by first-author surname. Same first author → order by year ascending.
- Author list uses `Last, F. M.` with `&` before the final author.
- Article titles: sentence case. Journal, book, database titles: title case.
- Italics on journal names and volume numbers. Issue in `(parentheses)`, not italic.
- Data-set entries carry `[Data set]` per APA 7 §10.9. Report entries include the report number in `(WMO-No. 1194)` style.
- DOI rendered as a clickable `https://doi.org/...` URL. Prefer DOI over publisher URL when both exist.
- One blank line between entries; no bullet markers (the list is styled as paragraphs so Pandoc can CSS-apply hanging indents).

## 10. File layout checklist

A compliant tutorial has, in order:

1. YAML front matter (§2).
2. `# H1` title matching or paraphrasing `title:`.
3. One or more numbered `## H2` sections, each optionally containing `###` / `####` subsections.
4. Mixture of static code (§4.1), output blocks, callouts (§5), tables (§6), prose paragraphs with pkgdown links (§7).
5. Optional `## N. References` in APA 7 (§9).

## 11. Common mistakes to avoid

- Setting `number-sections: true` while also writing manual numbering — this causes Quarto to render `1 1. Overview`.
- Using ``` ```{r} ``` for static code intended to *not* execute — the chunk will try to run and may fail in CI.
- Mixing em-dashes (—) and hyphens (-) inconsistently in headings. Choose one and stick with it per file.
- Copying output blocks from `R` console verbatim without the `#> ` prefix (Quarto renders them as regular code blocks without the comment indicator).
- Linking pkgdown references with `.html` instead of `.md`. Pkgdown rewrites either, but `.md` is the local convention.
- Placing a bibliographic entry without a DOI. Always include a DOI or a stable publisher URL so readers can verify.

## 12. Source and precedence

This style is extracted from `vignettes/tutorial_1_setup.qmd` as the canonical reference. Tutorials 2–4 predate some of these conventions and may diverge in minor details (e.g., Tutorial 4 uses `::: callout-note` blocks and a formal `bibliography:` entry). When a conflict arises between this file and an existing tutorial, **this file wins** for new authoring; existing tutorials are retrofitted incrementally.

Last verified against `tutorial_1_setup.qmd` on package audit.
