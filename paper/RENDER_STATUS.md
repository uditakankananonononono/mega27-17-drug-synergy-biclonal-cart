# Paper render status (Lane A)

- PDF is rendered by `.github/workflows/render_paper.yml` (xelatex, 3 passes) and committed back by CI. Local pdflatex lacks packages; do not render locally.
- The workflow is path-filtered to `paper/main.tex` and the workflow file. Empty commits do NOT trigger it. To retrigger, touch `paper/main.tex` (e.g. a trailing comment line).
- Anomaly 2026-10-07: the run for bc8ada9 rendered successfully but its "Commit rendered PDF" step failed. Cause unknown (no log access from the build box). Retriggering via a main.tex touch (9a47e4e) succeeded and CI committed ff3897f. If a render run shows a failed commit step, retrigger once, then check job steps via the Actions API.
- Current state: PDF at ff3897f, 54pp, matches main.tex. Pages 1, 10-11 (P4-CC1 and abstract clause) inspected visually.
