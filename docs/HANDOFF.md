# Handoff — Chapter 3–5 report review follow-up

Written 2026-09-17 on the laptop, to continue on another machine. Start a
Claude Code session in the repo root and ask it to read this file.

## Where things stand

A supervisor-style review of the Chapter 3, 4 and 5 reports (the Persian
`docs/ch*_report_fa.tex` / `.pdf` and the English `docs/ch3_report.html`,
`docs/ch4_report.html`) raised these issues. Status of each:

| # | Issue | Status |
|---|---|---|
| — | Formula derivations as a Persian appendix (پ) in all three reports | **Done**, pushed |
| — | Text-only fixes (wrong numbers, disclosures, wording) in all reports | **Done**, pushed (`373b9e1`) |
| 3 | Torque boxes sized from the perturbation they face (oracle knowledge) | **Code done** (`ec104e7`); reruns in progress; reports **not yet updated** |
| 2 | Runs counted as walking although the foot lifts off / slips | **Code done** (Ch4 `265ccef`, Ch3 `cca8dd4`, untested); reruns in progress; reports **not yet updated** |
| 1 | The model is 74 kg; published RABBIT is ~32 kg | Open — needs a decision |
| 4 | Three ε values (Ch3 draws conclusions at 0.5; some Ch4 diagnostics at 0.35) | Open — needs reruns at ε = 0.20 |
| 5 | Chapter 4 L1 figure panels labelled "Case II/III" in the reverse of the robust table's numbering | Open — relabel and redraw |
| 5 | Chapter 5 x0 pole-admissibility margins not reported | Open — compute and add |

## Decisions already made (do not re-ask)

- **Torque box rule** (`p.box.rule = 'rating'`, `Chapter4/ch4_params.m`): one
  actuator rating for every controller, scale and load — **556 N·m** (1.25 × the
  heaviest feedforward in the design envelope: scales ≤ 1.5, hip loads ≤ 70 kg,
  on posture_195); **731 N·m** for Case IV (scale 3). `'thesis'` reproduces the
  old per-case rules; old result files resolve to it (`ch4_box_rule`). The L1
  thesis form (box on μ1 only) keeps `p.l1.u_max`.
- **Physical validity** (`Chapter3/Analysis/ch3_validity.m`, used by
  `ch4_validity`): a sample is invalid when **Fz ≤ 0** (lift-off) or
  **|Fx|/Fz > 0.4** (`p.limits.mu_s`, slip). A run is **valid up to the step
  containing its first invalid sample**. `ch3_step` / `ch4_step` record the held
  torque and the contact force (Ch4: TRUE model) at every solver point.
- **Leak / footstrike options** restored in Ch4 (`p.l1.alpha_leak`,
  `p.l1.impact_estimate` = carry | continuous | fold), off by default, so the
  mitigation table's rows can be rerun.

## What to run next

1. Copy from the laptop into `Results/` (git-ignored, not on GitHub):
   `ch4_result_2026-09-17_17-27-44.mat` (robust preset, rating rule) and
   `ch4_result_2026-09-17_17-47-55.mat` (Case IV), plus optionally the figure
   folders `ch4_robust_2026-09-17_17-27-44/`, `ch4_case4_2026-09-17_17-47-55/`.
   Without them the driver simply reruns those two presets (~20 min).
2. In MATLAB at the repo root: `startup`, then `ch4_rerun_all`
   (`Chapter4/Analysis/reruns/`). It resumes and skips finished presets. Stages:
   `robust case4 l1 load ch3tests ch3table long_robust long_sweep long_nine
   long_oos long_summary`. On macOS/Linux, `Chapter4/Analysis/reruns/run_reruns.sh`
   runs one stage per MATLAB session with retries.
3. **Check `ch3tests` first**: the Chapter 3 validity code (`cca8dd4`) has not
   been run. The collocation order test is a *known, intentional* failure;
   anything else failing is new.
4. Outputs: `Results/ch4_result_*.mat` (presets), `Results/reruns/ch3/table.log`,
   `Results/reruns/ch4/progress.log` and `summary.log`.

Expect roughly 2–3 hours for everything after the presets.

## Results so far (rating rule + validity, 25 steps)

Robust preset (`17-27-44`), "valid" = steps before first lift-off/slip:

| Controller | ×1 valid, max‖η‖ | ×1.5 valid, max‖η‖ | ×0.7 valid, max‖η‖ |
|---|---|---|---|
| clfqp | 1, 0.25 | 0, 16.6 | 0, 14.5 |
| clfqp_con | 25, 0.29 | 25, 14.7 | 0, 15.3 |
| rclfqp_con | 25, **2.85** | 25, 4.24 | 5, 4.34 |

Case IV (`17-47-55`, box 731): clfqp falls in step 4 (2 valid), clfqp_con falls
in step 3 (2 valid), rclfqp_con walks 25 steps, all valid, max‖η‖ 1.89.

Notable: validity changes verdicts (unconstrained clfqp slips beyond 0.4 even
with a perfect model); the larger box changes the robust law itself at ×1
(max‖η‖ 1.13 → 2.85, peak torque 195 → 371 N·m).

## Report updates to make when the reruns finish

Chapter 4 (`docs/ch4_report_fa.tex` + rebuild PDF, and `docs/ch4_report.html`):
- the box rule paragraphs (robust section, L1 section, load study) — replace
  the per-case rules and the disclosures added in `373b9e1` with the rating;
- tables `tab:robust`, `tab:case4`, `tab:l1`, `tab:load`, `tab:loads` (box
  column), `tab:track`, `tab:rate`, `tab:plateau`, `tab:mitig`, `tab:window`,
  `tab:oos`, `tab:claims`: add a valid-steps entry, update numbers;
- the "negative normal force" paragraph → now scored by validity;
- abstract, conclusion, and every number quoted from those tables;
- `Chapter4/Analysis/ch4_doc_figures.m`: repin the result stamps, redraw
  figures; `docs/CH4_UNCERTAINTY.md` §5a.
Chapter 3 (`docs/ch3_report_fa.tex`, `docs/ch3_report.html`): the controller
table `tab:ctrl` from `Results/reruns/ch3/table.log`, with valid steps.

## Conventions that bite

- **MATLAB**: never run two `-batch` sessions at once (they hang). On the Mac,
  `-batch` aborts intermittently in libcurl — retry; long jobs checkpoint.
  stdout is buffered, so long runs log to a file.
- **Persian PDFs**: build with `latexmk -xelatex <file>.tex` in `docs/`. Do not
  use `\appendix` — with this xepersian preamble it loops forever; the
  appendices number sections by hand (`\renewcommand{\thesection}{پ}`). If a
  build is killed, delete the `.aux` and `.out` before rebuilding.
- **Figures** are drawn at the width the report prints them, text ≥ 9 pt
  (`ch3_doc_figures`, `ch4_doc_figures`, `ch5_doc_figures`).
- **Commits** end with a `Co-Authored-By` line; keep another session's
  uncommitted `Chapter6/` work out of these commits.
