# Infection dynamics for ENSCI 3120

An R/Shiny teaching app exported with **Shinylive** and hosted on **GitHub Pages**. Each student's browser runs its own R session. GitHub serves static files; there is no shared Shiny server allocating one R process per student.

## Deploy with GitHub Pages

1. Extract this ZIP. Open the `sir-shiny-app` folder. It contains `app/`, `tests/`, and a hidden `.github/` folder. Upload the contents, not the ZIP, to your repository root.
2. Create an empty GitHub repository, for example `SIRShiny`. Do not initialize it with a README if using the commands below.
3. In that repository choose **Settings → Pages → Build and deployment → Source → GitHub Actions**.
4. Open Git Bash in the extracted `sir-shiny-app` folder and run:

```bash
git init -b main
git add .
git commit -m "Add browser-based infection dynamics teaching app"
git remote add origin https://github.com/loeckete/SIRShiny.git
git push -u origin main
```

Replace the remote URL if you use another repository name. These commands assume a new, empty repository; do not force-push over an existing project. Git tracks `.github` even if your file explorer hides it.

5. Open **Actions → Test and deploy SIR app**. Model and Shiny tests must pass before export and deployment. No shinyapps.io account, deployment token, or paid Shiny server is needed.
6. After the workflow is green, open the URL from the deployment job. For the example repository above it will be `https://loeckete.github.io/SIRShiny/`. **That is an example destination, not an already published app.**
7. Complete the browser checklist below before using it in class. Share the URL early so students can load R before the activity begins.

If your repository already contains a README, clone it first, copy this package's contents into that checkout, then `git add .`, commit, and push. If Actions does not start, check the default branch is `main` and `.github/workflows/deploy.yml` is present at the repository root.

## What students can do

- SIR with lasting immunity, SIS with return to susceptibility, SEIR with an exposed stage, or SIRD with additional disease mortality.
- Change initial populations and rates, record a prediction, and press **Run simulation**.
- View host counts, entries to and exits from the infectious group, and a susceptible–infectious phase plot.
- Keep one completed run as a dashed comparison. Both runs' settings appear below the graph.
- Export results to CSV, including all settings used for the displayed run.
- Open a threshold explanation after making a prediction.

The initial page deliberately has no simulated answer. Editing controls does not change existing results until Run is pressed. A status message identifies changed settings and validation failures. Predictions remain in the current session only and are not submitted to Top Hat or saved remotely.

## Thursday's 20-minute investigation

| Scenario | Initial S / I / R | beta | gamma | m | Duration |
|---|---|---|---|---|---|---|
| Baseline | 99 / 1 / 0 | 0.004 | 0.20 | 0 | 60 days |
| Lower transmission | 99 / 1 / 0 | 0.0015 | 0.20 | 0 | 60 days |
| Faster recovery | 99 / 1 / 0 | 0.004 | 0.40 | 0 | 60 days |

Each preset resets all controls, including model and demography. Click **Load scenario settings**, make a prediction, then click **Run simulation**. Loading is separate from running so students can commit to a prediction.

- **0–2 minutes:** Open app. Reconnect S, I, and R to the reading.
- **2–7:** Predict the baseline, run, and explain the peak. Keep the completed run for comparison.
- **7–12:** Load Lower transmission. Explain which change should occur before running. Compare with baseline.
- **12–17:** Answer the existing short-answer Top Hat question: does a decrease in infectious hosts mean transmission stopped? Use entries and exits to justify the answer.
- **17–20:** Discuss a reasonable explanation. Faster finishers test Faster recovery.

Expected: baseline initial growth is 0.196 hosts/day; susceptible threshold is 50 hosts. Lower transmission raises the threshold to approximately 133.3 hosts and I declines immediately. Faster recovery raises it to 100 hosts and I also declines immediately. New infection can remain positive during decline. These are hypothetical model settings, not disease-specific estimates.

## Model details and notation

See the app's **Model and assumptions** tab for the equations. This implements the four model structures on the EcoEvoApps page, with an original base-R solver and interface. It does not depend on the `ecoevoapps` R package or copy its code, images, or branding.

Transmission is density dependent: `beta*S*I`, with counts in a fixed study area. It is not `beta*S*I/N`. Time is in days. Background birth rate equals background death rate `m`; `v` is the fraction of newborns entering the immune group. When `m=0`, changing `v` has no effect. Use initial R to represent immunity already present at time zero. D counts disease deaths only, starts at zero, and is excluded from the living population N.

The reading's `m` is total loss from the infected group. Here `gamma` is recovery, `m` is background mortality, and `mu` is extra disease mortality. For the default SIR model, background mortality is zero, so the reading's loss rate equals `gamma`.

For SIR/SIS/SIRD with I>0, I increases when `beta*S > gamma+m+mu` (mu is zero except in SIRD). For SEIR, immediate I growth depends instead on `a*E > (gamma+m)*I`. The SEIR transmission-chain threshold shown in the app accounts for survival through E and must not be confused with the immediate I-growth condition.

Continuous deterministic states can be fractional. Very small nonzero I is not a stochastic extinction prediction. Births balance background deaths; without disease mortality, the living total stays fixed. With SIRD, living hosts plus accumulated disease deaths stays fixed.

## Run locally or export manually

Install R, then run from the repository root:

```r
install.packages(c("shiny", "shinylive"))
shiny::runApp("app")
```

To test and export:

```bash
Rscript tests/model-tests.R
Rscript tests/shiny-tests.R
```

```r
shinylive::export("app", "site")
httpuv::runStaticServer("site")
```

Use an HTTP server, not a double-clicked `index.html`. Only `app/` goes into the exported application; tests and instructor documentation are excluded. Runtime dependencies are Shiny and base R. The base-R RK4 solver uses small substeps bounded by the fastest process rate. Model tests check conservation, nonnegativity, an analytical solution, step refinement, and independent SciPy DOP853 reference values.

## Browser checklist before class

- Open the deployed URL in a fresh desktop browser session. Wait for the R runtime and controls.
- Run Baseline; confirm the host and flow curves appear. Save the comparison.
- Load Lower transmission; confirm old plots remain with a changed-settings message until Run. Then check dashed baseline and declining current I.
- Try Faster recovery and all four model choices. Check the phase and equations tabs.
- Download CSV; confirm its settings match the displayed run, including when controls have changed without running.
- Enter a negative initial population; confirm an informative error preserves the prior successful run.
- Test on at least one student-style laptop and phone on classroom Wi-Fi. Mobile plots may require landscape orientation.

**Scale:** computation happens on each device, so there is no shared R-server session limit. This package has not been load-tested with 400 simultaneous clients. Initial runtime downloads, hosting bandwidth, browser compatibility, and campus Wi-Fi can still affect access. Ask students to open it before class and keep the slide-deck fallback graphs available. Reloading loses session predictions and comparisons; CSV downloads preserve numerical results.

## Verification status

The build environment had no R installation and could not fetch one. Independent reference trajectories and a Python reproduction of the numerical scheme were checked locally. Native R/Shiny execution, Shinylive export, and browser behavior are not yet verified here. The included GitHub workflow runs R model and Shiny server tests before publishing, and records R/package versions in `build-info.txt`. A green workflow still needs the browser checks above.

## Sources and attribution

- EcoEvoApps team (Gaurav Kandlikar, Madeline Cowen, Kenji Hayashi, Rosa McGuire, Xinyi Yan, Marcel Vaz), *Infectious disease dynamics (compartment/SIR models)*: https://ecoevoapps.gitlab.io/docs/articles/sir-dynamics.html
- Posit, R Shinylive documentation: https://posit-dev.github.io/r-shinylive/
- GitHub Pages deployment action: https://github.com/actions/deploy-pages

The delivery architecture follows the Shinylive + GitHub Actions + GitHub Pages setup used in your previous app package. The live population-growth repository could not be retrieved in this environment, so this is not claimed to be an exact copy of its current workflow.
