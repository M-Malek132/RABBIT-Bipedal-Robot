# Handoff — Chapter 3–5 report review, and the Chapter 3 gait campaign

Two pieces of work, both finished, both recorded here rather than pending:

1. **The Chapter 3–5 report review follow-up** (2026-09-17/18, commits
   `ccce627`, `a910b7c`, `517e5d2` and the test-stamp commit after them).
2. **The Chapter 3 gait campaign** (2026-09-18 → 09-23): what this 74 kg model
   can walk at, and within what actuator budget. Jump to
   "Chapter 3 gait campaign" for the headline — a verified gait at the
   **1.2 m/s design speed**, and a torque floor of **162.5 N·m** against a
   declared limit of 120 — at which box a gait also meets the **15 N·s**
   impulse limit (2026-09-26).

Then a third, **finished** (2026-09-25/26): the follow-ups this file listed as
open. Report edits, the thesis-row rerun, the moving-x₀ study and the impulse
march are all done and written into the reports. Jump to **"Follow-ups of
2026-09-25"**.

## Where things stand

A supervisor-style review of the Chapter 3, 4 and 5 reports (the Persian
`docs/ch*_report_fa.tex` / `.pdf` and the English `docs/ch3_report.html`,
`docs/ch4_report.html`) raised these issues. Status of each:

| # | Issue | Status |
|---|---|---|
| — | Formula derivations as a Persian appendix (پ) in all three reports | **Done**, pushed |
| — | Text-only fixes (wrong numbers, disclosures, wording) in all reports | **Done**, pushed (`373b9e1`) |
| 3 | Torque boxes sized from the perturbation they face (oracle knowledge) | **Done** — one 556/731 N·m rating; rerun and written up (`a910b7c`) |
| 2 | Runs counted as walking although the foot lifts off / slips | **Done** — validity scored, tests pass, reports rewritten (`a910b7c`) |
| 1 | The model is 74 kg; published RABBIT is ~32 kg | **Done** — kept at 74 kg, disclosed in Ch3 and Ch4 (Ch5 has no robot) |
| 4 | Three ε values (Ch3 draws conclusions at 0.5; some Ch4 diagnostics at 0.35) | **Done** — Ch3 table and all three Ch4 diagnostics rerun at ε = 0.20 |
| 5 | Chapter 4 L1 figure panels labelled "Case II/III" in the reverse of the robust table's numbering | **Done** — `ch4_plot_uncertainty` numbers panels by mass scale; figures redrawn |
| 5 | Chapter 5 x0 pole-admissibility margins not reported | **Done** — margins added to `ch5_report_fa.tex` and `ch5_report.html` |

All seven review items are closed. What was left was ordinary follow-up, not
review debt; each now has a driver (see "Follow-ups of 2026-09-25"):

- ~~The Chapter 5 runs start at rest, so the limitation is never exercised~~ —
  **done**: `ch5_moving_x0` ran all four configurations from a moving `x0`;
  results in Ch5 `sec:moving` / caveat 4 (summary below).
- ~~The `l1` thesis-form row of `tab:l1` comes from the old run~~ — **done**:
  `ch4_eps020_diagnostics('thesisrow')` reproduces it to the digit and adds the
  valid counts (25 / 0 / 0 at ×1 / ×0.7 / ×1.5); both reports updated.
- ~~`docs/ch4_report.html` figures are re-embedded by hand~~ — **done**: each
  `<img>` names its file (`data-fig="figures/…png"`) and
  `python3 docs/embed_figures.py docs/ch4_report.html` rebuilds them
  (`--check` reports stale copies without writing). The ch3 and ch5 HTML
  images match no file in `docs/figures/`, so they are not tagged.

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
- **L1 box** (2026-09-17): the l1 preset runs at the 556 N·m rating on the
  TOTAL applied torque (`p.l1.constrain_applied = true`, already the default),
  not the thesis form's μ1-only box. The report must say the formulation
  differs from the thesis. The "l1 torque box 65 Nm … raising it to 244 Nm"
  warning from `ch4_load_gait` is harmless here: that `p.l1.u_max` is unused
  under the rating.
- **Workflow**: update the Chapter 3/4 reports and commit only after ALL
  reruns have finished, in one pass.

## How to reproduce the runs behind the reports

Everything below has been run; this is the recipe, not a to-do list. `Results/`
is git-ignored, so on a fresh machine the runs have to be redone (~3 hours) or
the folder copied across.

1. One MATLAB session at a time. In MATLAB at the repo root: `startup`, then
   `ch4_rerun_all` (`Chapter4/Analysis/reruns/`). It resumes and skips finished
   presets. Stages: `robust case4 l1 load ch3tests ch3table long_robust
   long_sweep long_nine long_oos long_summary`. On macOS/Linux,
   `Chapter4/Analysis/reruns/run_reruns.sh` runs one stage per MATLAB session
   with retries.
2. **ε = 0.20 diagnostics** (issue 4): `ch3_table_rerun(0.20)` →
   `Results/reruns/ch3_eps020/table.log`, then `ch4_eps020_diagnostics`
   (`Chapter4/Analysis/reruns/`) → the κ table, the robust growth ablation and
   the thesis predictor's θ̂ overshoot, in `Results/reruns/ch4_eps020/`.
3. Outputs: `Results/ch4_result_*.mat` (presets, the four the reports cite are
   `2026-09-17_17-27-44` robust, `_17-47-55` case4, `_20-33-53` l1,
   `_20-36-27` load), `Results/reruns/ch3/table.log`,
   `Results/reruns/ch4/progress.log` and `summary.log`.
4. Figures: `ch4_doc_figures` reads the four stamps pinned at its head — repin
   them when the results change. Then `python3 docs/embed_figures.py
   docs/ch4_report.html` re-embeds them in the English report.
5. Test suites: `ch3_test_all` (83 checks, 1 known failure, 5.3 s) and
   `ch4_test_all` (54 checks, 0 failures, 5.1 s), both on R2023b, 2026-09-18.
   The collocation order test is a *known, intentional* failure; anything else
   failing is new.

## What the reruns found

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

### Rating reruns (2026-09-17 evening)

All stages done, none failed. Full numbers: `Results/reruns/ch4/summary.log`,
`Results/reruns/ch3/table.log`, `Results/ch4_result_2026-09-17_20-33-53.mat`
(l1) and the load result after it. What validity changes:

- **l1 preset** (box 556 on total torque): `l1_con` valid 25/25 at ×1 and ×1.5
  but **0 at ×0.7** (slip, max μ 1.58); unconstrained `l1` valid 0 at ×0.7 and
  ×1.5 (Fz down to −226 N). clfqp valid 1 at ×1.
- **Robust plateau** (rclfqp_con, 60 steps): ×1 valid **34**/60 (slip in 35),
  ×0.7 valid **5**, ×1.5 valid 58, Case IV valid 120/120.
- **Nine long runs** (default nrm075, l1_con): s100 and s150 valid 60/60;
  s070 valid 0; every random-load sequence and kg46 lose validity within
  2–22 steps (slips, Fz to −560 N); s150 with `l1` valid 0. The other variants
  (leak, cap, fold, continuous, dt 0.5 ms) behave the same; dt 0.5 ms is the
  best on random loads (rnd2 53, rnd3 56, rnd5 60).
- **Predictor-rate sweep** (×1.5): `l1_con` valid 25/25 at every rate except
  a = 800/nrm 0 (falls in 16) and a = 700/nrm 0.75 (falls in 21); `l1` always 0.
- **Ch3 table** (ε 0.5): PD 4/4 valid; clfqp 1; clfqp_con 0, 0, 1 at 60/80/100%.

### ε = 0.20 reruns finished (2026-09-17, issue 4)

`Results/reruns/ch3_eps020/table.log`, `Results/reruns/ch4_eps020/`:

- **Ch3 table**: falls come one step later than at ε 0.5 (clfqp 3, clfqp_con
  1/2/4) but the VALID counts are unchanged (4, 1, 0, 0, 1). The Chapter 3
  conclusion survives at 0.20.
- **κ table**: chatter picture as at 0.35 (κ 0: 417–558 N·m per sample; κ ≥
  0.66: 1–2 N·m). Validity sharpens it: κ 0 → 0 valid everywhere, κ 0.33 → 0
  valid at ×0.7, κ 0.66 → 1 of 3 at ×0.7, κ ≥ 0.99 → 3 of 3 everywhere.
- **Growth**: peak V still grows ×4700 (4.9e-5 → 0.229) over 25 steps, but at
  0.20 all 25 steps are VALID — it no longer fells Case I in step 21 as at
  0.35. Survives norows/nobox/D2=0; vanishes with D1=0 (max‖η‖ 2.85 → 0.26).
  Vanishes without the layer too, but that run is 0 valid (chatter).
- **Predictor**: thesis θ̂ overshoots the true θ **7.2×** at ×0.7 and **2.0×**
  at ×1.5 (3–5× at 0.35); plant predictor 1.29× / 1.01×. Scale 1 ratios are
  meaningless (true θ ≈ 0) and now report NaN.

## Chapter 3 gait campaign (2026-09-18 → 09-23)

Separate from the report review above: a campaign to find out what this 74 kg
model can actually walk at, and within what actuator budget. All the gaits live
in `Results/reruns/` (git-ignored); the drivers are in
`Chapter3/Analysis/reruns/`.

### Speed: the design target is reachable, 0.5 m/s is not

`ch3_params` asks for **1.2 m/s** (`p.v_des`, `enforce_nec1 = true`), but every
stored gait carries `enforce_nec1 = 0`, so the equality was never imposed and
`posture_195` walks at **1.5628 m/s** — 30% over target, with no margin
anywhere.

`ch3_speed_ladder` imposes the equality and marches down. Verified gaits now
exist at 1.5628, 1.50, 1.45, 1.40, 1.35, 1.30, 1.25, **1.20** and 1.15 m/s
(`Results/reruns/speed_ladder/`, `gait_v1200.mat` is the design-speed one).

- **The landing-on-the-cone problem was a speed artefact.** At 1.5628 the
  impact impulse sits exactly on the μ = 0.4 cone (row 13 = +5.6e-9); at 1.20
  the landing friction demand is **0.161**, slack by 5.15.
- **The march walls at 1.15 m/s.** Every step to 1.10, down to 0.013, refuses
  in ~15 s. Torque, stance friction and the Fz floor are all active at once;
  no single relaxation releases it, and the closest pair (μ_s 0.6 **with** the
  Fz floor at 10 N) only reaches 1.10. So sub-1.15 walking is bought by
  weakening the contact model.
- Cold-start routes to low speed all failed: `ch3_stage3_from_scratch` (N=41)
  fails at its bare stage; N=81 with NEC1 off runs the speed away to 1.14 with
  peak torque 2335 Nm; N=81 anchored was still oscillating at 1e-1 after three
  hours when the session was killed. `ch3_speed_march`'s seed
  (`ch3_gait_full_constrained`) is **stale** under the 74 kg model and that
  campaign cannot run as written.

### Torque: the floor is 162.5 N·m, not the declared 120

`ch3_torque_march` cuts the box from the 1.2 m/s gait with the speed free.
Landed rungs, all verified trajectories:

| box | speed | max\|c\| | verify | impulse |
|---|---|---|---|---|
| 180 | 1.2815 | 2.8e-6 | 3.2e-5 | 18.77 Ns |
| 170 | 1.2857 | 6.3e-9 | 5.1e-5 | 18.88 Ns |
| 165 | 1.2935 | 2.2e-6 | 6.3e-5 | 18.87 Ns |
| **162.5** | **1.2921** | **7.3e-11** | **7.1e-5** | 19.07 Ns |
| 160 | — | 6.2e-4, 9.2e-4 | — | refused twice |

- **162.5 is a floor, not a stopping point** (`ch3_torque_wall`): five attempts
  at 160 all refuse, and relaxing the contact limits makes the torque row
  **worse** (μ_s → 0.5 gives 1.1e-2, Fz floor → 25 N gives 1.3e-2) because the
  optimizer spends the new slack on the contact, not on torque. The barrier is
  the trajectory.
- **Speed runs backwards on this branch.** The gait got *faster* as the box
  tightened (1.2815 → 1.2921 m/s). Pinning the speed back to 1.20 collapses the
  solve to a spurious point (max\|c\| 3.0). **1.2 m/s and 120 N·m cannot be
  pursued together by continuation.**
- **~~Half realizable only~~ — closed 2026-09-26.** These gaits carried
  18.8–19.1 N·s against a declared 15, with that gate off. `ch3_impulse_march`
  (box held at 162.5, speed free, 1 N·s rungs) landed 18, 17, 16 and **15 N·s**.
  The last is clean: max c −2.1e-9, verify 1.2e-4 against a 1e-3 tolerance, at
  **1.3932 m/s**, T 0.3144 s, L 0.438 m, min Fz 72.3, μ 0.3995, torque and
  impulse active. File: `Results/reruns/impulse_march/u162_I1500.mat`, on the
  user's PC. Speed ran up again (1.2921 → 1.3932), and the Fz floor went
  inactive from 16 N·s.

### Stride: not a design variable in this transcription (a negative result)

If slow walking on a long stride is what exhausts the contact, shorten the
stride. That route is closed, and it cost about seven hours to establish, so it
is written down rather than left to be rediscovered.

Row 19 was added to `ch3_col_constraints` for the test
(`L_step <= p.limits.step_len_max`, gated, **off by default**), plus
`p.cost_normalize` and `p.cost_scale` in `ch3_col_cost` (both default to the
old behaviour). `ch3_test_all` is green with them: 82 pass, 0 fail, 1 xfail.

**The stride will not come in by even one millimetre.** Every run: the
violation stays exactly constant, the equalities hold, the stride moves ~1e-5 m
and the step size collapses to 1e-8 with first-order optimality ~1e6.

| varied | tried | result |
|---|---|---|
| seed gait | the 1.2 m/s gait (μ, Fz, torque all active) and `posture_195` (friction slack 0.155) | identical stall — "active corner" does not explain it |
| rung size | ceilings 0.40 and 0.43 m; cuts of 1 mm and 0.1 mm | smaller violation did not help |
| objective | per-distance, and torque² alone (`cost_normalize`) | not the obstacle |
| objective scale | ×1e-4 (`cost_scale`), so feasibility dominates the merit | moved 1.3e-5 m — identical |
| problem scaling | `ScaleProblem` off | diverges, spurious, 36 off a rollout |

So step length is an **output** of the dynamics, periodicity and virtual
constraints here — not something to trade. This is a claim about this
formulation and this solver, not a proof that no short-stride gait exists; a
different parametrisation (shorter phase sweep, another posture branch, a cold
solve built around a short stride) might find one. It is not reachable by
continuation from the gaits that exist.

**Four hypotheses were proposed and all four were wrong**: that slow walking
costs torque (falsified — 57 N·m of headroom went unused); that the active
contact corner blocked the stride (falsified by the slack-friction seed); that
the per-distance objective opposed it (falsified by `cost_normalize`); and that
the merit-function balance starved it (falsified by `cost_scale`). Treat any
fifth mechanism with suspicion unless it is tested first.

**Consequence for walking slower:** 1.15 m/s is the floor for this model at
μ = 0.4 with a 50 N normal-force floor. Below it needs a weaker contact model
(μ = 0.6 with a 10 N floor reached 1.10) or a genuinely different gait.

### What this changed in the reports

Committed in `9132285` and `8d0c6b9`. The claim "no verified gait below 195 N·m
has been found" appeared three times in Chapter 3 and is now wrong; all three
are rewritten (`sec:torquefloor` carries the rung table and the wall probe).
Chapter 4's `sec:box120` said the remaining work moved to Chapter 3 — that work
is done and the answer is negative, so the gap to 120 narrows from 1.63× to
1.35× and does not close.

`sec:stridefixed` in Chapter 3 and the matching HTML section carry the negative
result above. `sec:budget` and `sec:box120` have now been ported into
`docs/ch4_report.html`, which previously had none of the 120 N·m material.

Both items this paragraph listed as open are closed (2026-09-25): the PDFs
were rebuilt, and `tab:claims` now scores every claim at 120 N·m as well (see
the follow-ups).

## Follow-ups of 2026-09-25

### Done — report edits, both languages

- **Budget caveat in `tab:claims`** (`ch4_report_fa.tex`, `ch4_report.html`).
  A fourth column scores every claim at the declared 120 N·m, from
  `Results/box120/box120.log`. Remarks 4.6 and 4.7 and §4.2.4 fail there. "L₁
  performs like the CLF-QP" holds only in that both boxed laws fall in step 2,
  and Case IV and the load were never run at 120 N·m ("not run"). The caption
  and the HTML plate note say the other two columns are at 556/731 N·m.
- **Ch4 conclusion.** The Persian one opens with the 120 N·m reading, and no
  longer calls "a gait optimized under a 120 N·m box" the remaining work:
  Chapter 3 tried, and stopped at 162.5. The English "limitations" list gets
  the budget item it lacked.
- **Ch5 tolerance** (`ch5_report_fa.tex` admissibility paragraph,
  `ch5_report.html` table and note). The reports read the first rung,
  −ḣ(0) ≤ p₁h(0) (3.78 / 3.60 m/s), as the tolerable closing speed. It is
  only a ceiling. For carts starting together at speed s, every higher
  derivative of h is zero, so y_i ≥ 0 ⇔ s ≤ h(0)/Σ_{j≤i} 1/p_j, and the five
  exactly checkable rungs allow only **0.894 / 0.852 m/s**. The pendulum's
  15.0 / 12.5 are gone: its x₀ is the top of the end effector's circle, so
  ḣ(0) ≡ 0 for any rotation.
- **HTML figures**: `docs/embed_figures.py`, and `data-fig` on the eight ch4
  images (see the loose ends above).

### Runs done, 2026-09-25 evening, and what they changed

- **Thesis row** (`Results/reruns/ch4_eps020/summary.log`): `l1_con` in the
  thesis form reproduces the 2026-09-13 row exactly (max‖η‖ 0.381 / 14.474 /
  18.704, min Fz 49 / −1403 / −978 N, falls in step 14 at ×1.5) and is valid
  25 / 0 / 0 at ×1 / ×0.7 / ×1.5. The applied torque reaches 199 / 361 / 538 N·m
  against its 243.8 N·m box on μ₁. `tab:l1` (.tex, via a four-field `\fallc`) and
  the HTML row, tag and note are updated; "old box rule" is gone.
- **Moving x₀** (`Results/reruns/ch5_moving_x0/moving_x0.log`). Spring-mass:
  the bisection matches the closed form to 1e-16, and rung 5 binds (s_adm
  0.894 / 0.852 m/s). h stays ≥ 0 through 1.5× s_adm and goes negative at 2×
  (−0.162 / −0.159) and 3×, so the corollary is conservative by 1.5–2×. Poles
  scaled by 1.05c keep h ≥ 0 at 1.1, 1.5 and 3×, but at 3× the peak force goes
  from 84 / 80 N to 413 010 / 71 066 N. Pendulum: s_1 does not exist, rung 3
  binds (s_adm 3.73 / 3.41 rad/s, s_top 3.06 / 2.78). Only 1.87 rad/s at −1.0 m
  completes; 17 other runs, admissible and pole-scaled ones included, go
  non-finite within 0.13–1.17 s with the QP infeasible. **Cause (above s_top):**
  the barrier row's coefficient (k/J_m)∇p_yᵀD⁻¹ is exactly zero at the upright
  start (∇p_y = 0), so h does not have relative degree 4 at x₀. Harmless from
  rest (row slack, y₄ = 4950); fatal once y₄(x₀) < 0. The 1.70 rad/s run at −0.5 m
  fails at 1.17 s for a reason not located. The folded-arm zero of the same
  coefficient (θ₁ = 0, θ₂ = ±π; nominal runs fold to −3.10) is a candidate,
  unconfirmed: `one_run` saves no trajectory. Written into Ch5 `sec:moving`
  (Persian) and caveat 4 (English), with the "inherits the full relative degree"
  and "Theorem 5.1 holds without exception" sentences qualified and the
  limitation in both conclusions updated.

### The three runs: how they were specified, and item 3's results

All three items below are done (item 3's results are at its end). One MATLAB
session at a time. This Mac does **not** hold the
Chapter 3 campaign or most of the 09-17 reruns: `Results/reruns/speed_ladder/`,
`torque_march/`, `stride_march/` and `ch4_eps020/` are absent, as are the
`20-33-53` (l1) and `20-36-27` (load) results. They live on the machine that
produced them.

1. **`ch4_eps020_diagnostics('thesisrow')`**, then `('summary')`. Three
   25-step runs, a few minutes, any machine. Writes
   `Results/reruns/ch4_eps020/thesisrow_s*.mat`, one line each in
   `progress.log`, and THESISROW lines in `summary.log`. **Then**: rewrite
   `tab:l1`'s italic row with the valid count (`\rc` → `\rcv` in the .tex;
   "C, as §4.2 writes it" in the HTML), and drop "old box rule" from the
   caption and the plate note. What stays true: the thesis form's box is
   `p.l1.u_max` = 243.8 N·m, on μ₁ only, at every scale.
2. **`ch5_moving_x0()`**. About 40 closed-loop runs, perhaps half an hour,
   any machine. Log: `Results/reruns/ch5_moving_x0/moving_x0.log`, ending
   `CH5_MOVING_X0_DONE`. Read first: the spring-mass line "bisection agrees
   to …" should be about 1e-12. If it is not, `ch5_ecbf_admissible` and the
   closed form disagree, and that is the finding. **Then**, in both Ch5
   admissibility paragraphs, add for each configuration: s_adm, the multiple
   of it at which h first goes negative (if any), and whether poles scaled
   by c keep h ≥ 0. The pendulum's measured s_adm replaces "can break only
   through the centripetal terms".
3. **`ch3_impulse_march()`**. Needs
   `Results/reruns/torque_march/gait_u162.mat`, so run it on the campaign
   machine or copy that folder here first; without it the script stops with
   a clear error. Log: `Results/reruns/impulse_march/impulse_march.log`,
   ending in a `VERDICT` line and `IMPULSE_MARCH_DONE`. **Then**: Ch3
   `sec:torquefloor` and Ch4 `sec:box120` (both languages) say the gaits are
   "half realizable at best". Replace that with the verdict: either a fully
   realizable 162.5 N·m gait, or where the march stalls (a continuation
   result, not a floor).

   **First run (2026-09-25, user's PC): inconclusive, and its VERDICT line is
   wrong.** The first rung, 19.00 N·s (a 0.4% cut from the 19.07 seed), ended
   "converged to an infeasible point": exitflag −2, max c 1.05e-2, 1430 s. The
   violation was in neither the torque nor the impulse row; mu printed 0.400
   and min Fz 50.0. But fmincon's own log shows it meeting every limit on the
   way (feasibility 2.2e-8 at iteration 19, 2.0e-7 at 25, 1.1e-6 at 85), with
   first-order optimality at 1e8–1e9: the seed sits on a corner, and SQP chased
   the torque cost off the feasible set. The march then quit without bisecting,
   because (19.07 − 19.00)/2 is below the old 0.25 N·s MIN_STEP, and wrote "the
   impulse comes down only to 19.07 Ns". That line is not a finding.
   The driver now keeps the best iterate that meets every limit in true
   units, stops once 30 iterations bring no better one, and uses it if the
   final iterate misses (still mesh-verified). It also cuts 1 N·s rungs
   (18 → 15), bisects to 0.1 N·s, logs the worst constraint row, and words a
   stall as a continuation result.

   **Rerun (2026-09-25 23:17 → 09-26 01:47): done — a gait meets both.**
   | cap | outcome | v m/s | min Fz | max c | verify |
   |---|---|---|---|---|---|
   | 18 | converged (exitflag 2) | 1.3162 | 50.0 | 1.4e-10 | 9.6e-5 |
   | 17 | final missed; iterate 104 kept | 1.3460 | 50.0 | 6.3e-6 (friction) | 1.3e-4 |
   | 16 | converged | 1.3707 | 67.1 | −2.3e-13 | 1.4e-4 |
   | **15** | converged, clean | **1.3932** | 72.3 | −2.1e-9 | 1.2e-4 |

   The 17 N·s rung is only a waypoint: `ch3_col_check_limits` flags its
   6.3e-6 friction excess, while the march's landing test (max c ≤ 1e-4) let
   it through. The 16 and 15 N·s gaits meet every limit on their own. Written
   into Ch3 `sec:torquefloor` (new `tab:impulsemarch`, both languages), the
   Persian Ch3 conclusion, and Ch4 `sec:box120` (both languages).

**PDFs.** Rebuilt and committed twice (`e22a457`, `d14953a`). The impulse
write-up changed `ch3_report_fa.tex` and `ch4_report_fa.tex` again, so those
two need one more latexmk build.

## Conventions that bite

- **MATLAB**: never run two `-batch` sessions at once (they hang). On the Mac,
  `-batch` aborts intermittently in libcurl — retry; long jobs checkpoint.
  stdout is buffered, so long runs log to a file. `diary` **appends**: delete
  the log first or a rerun doubles every count read back from it.
- **Persian PDFs**: build with `latexmk -f -xelatex <file>.tex` in `docs/`.
  **The `-f` is required**: the xepersian preamble raises real font errors
  ("B Nazanin" lacks U+066A/U+066B) and without `-f` latexmk discards the build
  and writes no PDF at all. It exits **12** even on a good build, so check the
  page count and `Latex failed to resolve` in the log, not the exit code
  (2026-09-18: ch3 10 pages, ch4 13, ch5 10, zero unresolved). Do not use
  `\appendix` — with this preamble it loops forever; the appendices number
  sections by hand (`\renewcommand{\thesection}{پ}`). If a build is killed,
  delete the `.aux` and `.out` before rebuilding.
- **Figures** are drawn at the width the report prints them, text ≥ 9 pt
  (`ch3_doc_figures`, `ch4_doc_figures`, `ch5_doc_figures`).
- **Commits** end with a `Co-Authored-By` line; keep another session's
  uncommitted `Chapter6/` work out of these commits.
