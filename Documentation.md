# Documentation

This repository contains an R implementation of the **corrected
Welch–Satterthwaite equation** for the effective degrees of freedom of a
weighted sum of independent variance components.

The implementation follows von Davier (2026), *A Corrected Welch
Satterthwaite Equation. And: What You Always Wanted to Know About Kish's
Effective Sample but Were Afraid to Ask* (arXiv:2602.20912).

## Contents

- [The corrected equation](#the-corrected-equation)
- [Repository files](#repository-files)
- [Calling the R functions](#calling-the-r-functions)
  - [Setup](#setup)
  - [Two variance components: `cws2()`](#two-variance-components-cws2)
  - [Two or more components: `cws()`](#two-or-more-components-cws)
  - [Applying `cws()` row-wise](#applying-cws-row-wise)
  - [Running the simulation script](#running-the-simulation-script)
  - [Input rules and edge cases](#input-rules-and-edge-cases)
- [Related arXiv papers](#related-arxiv-papers)
- [Classical references](#classical-references)
- [Committing the documentation to GitHub](#committing-the-documentation-to-github)

## The corrected equation

Let

$$
S^{2}=\sum_{k=1}^{K}w_{k}S_{k}^{2}
$$

be a weighted sum of independent variance components $S_{1}^{2},\dots,S_{K}^{2}$,
where $w_{k}>0$ and

$$
\frac{\nu_{k}S_{k}^{2}}{\sigma_{k}^{2}}\sim\chi^{2}_{\nu_{k}}.
$$

The corrected Welch–Satterthwaite estimator is

$$
\nu^{*}_{vD}
=
\frac{\left(\sum_{k=1}^{K}w_{k}S_{k}^{2}\right)^{2}}
{\sum_{k=1}^{K}w_{k}^{2}\frac{\left(S_{k}^{2}\right)^{2}}{\nu_{k}+2}}
-2.
$$

Equivalently, writing $S_{k}^{4}=(S_{k}^{2})^{2}$,

$$
\nu^{*}_{vD}
=
\frac{\left(\sum_{k=1}^{K}w_{k}S_{k}^{2}\right)^{2}}
{\sum_{k=1}^{K}\frac{w_{k}^{2}S_{k}^{4}}{\nu_{k}+2}}
-2.
$$

The original Welch–Satterthwaite estimator omits the $+2$ in the
denominator and the final $-2$:

$$
\nu^{*}_{WS}
=
\frac{\left(\sum_{k=1}^{K}w_{k}S_{k}^{2}\right)^{2}}
{\sum_{k=1}^{K}\frac{w_{k}^{2}S_{k}^{4}}{\nu_{k}}}.
$$

Boardman's (1974) "unbiased Satterthwaite" estimator has the same
denominator as the corrected estimator but does **not** subtract 2:

$$
\nu^{*}_{Boardman}
=
\frac{\left(\sum_{k=1}^{K}w_{k}S_{k}^{2}\right)^{2}}
{\sum_{k=1}^{K}\frac{w_{k}^{2}S_{k}^{4}}{\nu_{k}+2}}
=
\nu^{*}_{vD}+2.
$$

When all variance components are equal, $S_{k}^{2}\equiv\sigma_{E}^{2}$,
and the component degrees of freedom are ignored, the original
Satterthwaite structure reduces to Kish's (1965) effective sample size:

$$
n_{\mathrm{eff}}
=
\frac{\left(\sum_{k=1}^{K}w_{k}\right)^{2}}
{\sum_{k=1}^{K}w_{k}^{2}}.
$$

This last expression is a theoretical special case; the R functions
below implement the corrected equation for estimated variance
components with finite $\nu_{k}$.

### Direct R translation

For vectors `s2`, `nu`, and `w`, the corrected estimator is exactly:

```r
sum(w * s2)^2 / sum(w^2 * s2^2 / (nu + 2)) - 2
```

For $K=2$, `cws2(s2_1, s2_2, nu_1, nu_2, w1, w2)` computes:

```r
s2 <- w1 * s2_1 + w2 * s2_2
s2^2 / (w1^2 * s2_1^2 / (nu_1 + 2) +
        w2^2 * s2_2^2 / (nu_2 + 2)) - 2
```

## Repository files

| File | Purpose |
| --- | --- |
| `corrected_welch_satterthwaite.R` | Defines `cws2()` for two components and `cws()` for $K\ge 1$ components. |
| `simulate_corrected_ws.R` | Simulation study comparing the corrected estimator to the target $K\nu$. |
| `short_satter_07152026.tex` | LaTeX source of the derivation and paper. |
| `LICENSE` | CC0 1.0 Universal public-domain dedication. |
| `Documentation.md` | This file. |

## Calling the R functions

### Setup

The functions use only base R and have no package dependencies. From a
terminal, change to the repository directory and start R:

```bash
cd /home/matthias/python/Corrected-Welch-Satterthwaite
R
```

Then source the function file:

```r
source("corrected_welch_satterthwaite.R")
```

You can also run a single expression directly from the shell:

```bash
Rscript -e 'source("corrected_welch_satterthwaite.R"); cws2(2, 3, 10, 20)'
```

From R, you can source the file by its full path if your working
directory is elsewhere:

```r
source("/home/matthias/python/Corrected-Welch-Satterthwaite/corrected_welch_satterthwaite.R")
```

### Two variance components: `cws2()`

```r
cws2(s2_1, s2_2, nu_1, nu_2, w1 = 1, w2 = 1)
```

| Argument | Description |
| --- | --- |
| `s2_1`, `s2_2` | The two estimated variance components $S_{1}^{2}$, $S_{2}^{2}$. These are **variances**, not standard deviations. |
| `nu_1`, `nu_2` | The degrees of freedom of the two components. Both must be positive. |
| `w1`, `w2` | Positive weights. The default is equal weighting, `w1 = w2 = 1`. |

The return value is a single numeric value: the corrected effective
degrees of freedom.

#### Equal weights

```r
cws2(2, 3, 10, 20)
# [1] 31.67347
```

#### Unequal weights

```r
cws2(4, 9, 5, 8, w1 = 2, w2 = 3)
# [1] 12.93122
```

This is the value produced by the formula:

```r
s2  <- 2 * 4 + 3 * 9                 # 35
den <- 2^2 * 4^2 / (5 + 2) +
       3^2 * 9^2 / (8 + 2)          # 82.04286
s2^2 / den - 2                       # 12.93122
```

#### Using standard deviations

If you have standard deviations `sd1` and `sd2`, square them first:

```r
sd1 <- sqrt(2)
sd2 <- sqrt(3)
cws2(sd1^2, sd2^2, 10, 20)
# [1] 31.67347
```

#### Relationship to `cws()`

For two components, `cws2()` and `cws()` give the same result when the
arguments are supplied consistently:

```r
cws2(2, 3, 10, 20)
# [1] 31.67347
cws(c(2, 3), c(10, 20))
# [1] 31.67347
```

### Two or more components: `cws()`

```r
cws(s2, nu, w = rep(1, length(s2)))
```

| Argument | Description |
| --- | --- |
| `s2` | Numeric vector of variance components $S_{1}^{2},\dots,S_{K}^{2}$, all $\ge 0$. |
| `nu` | Numeric vector of component degrees of freedom, same length as `s2`, all $>0$. |
| `w` | Numeric vector of positive weights, same length as `s2`. Defaults to all ones. |

The return value is a single numeric value:

$$
\nu^{*}_{vD}
=
\frac{\left(\sum_{k=1}^{K}w_{k}S_{k}^{2}\right)^{2}}
{\sum_{k=1}^{K}w_{k}^{2}\left(S_{k}^{2}\right)^{2}/(\nu_{k}+2)}
-2.
$$

#### Basic call

```r
cws(c(2, 3), c(10, 20))
# [1] 31.67347
```

#### Unequal weights

```r
cws(c(2, 3), c(10, 20), w = c(1, 2))
# [1] 30.49231
```

Weights are relative: multiplying every weight by the same positive
constant does not change the result.

```r
cws(c(2, 3), c(10, 20), w = c(1, 2))
# [1] 30.49231
cws(c(2, 3), c(10, 20), w = c(2, 4))
# [1] 30.49231
```

#### Three components

```r
cws(c(1, 2, 3), c(5, 10, 15), w = c(1, 2, 3))
# [1] 29.40575
```

#### Corrected versus original Welch–Satterthwaite

```r
s2 <- c(1, 1)
nu <- c(1, 2)

cws(s2, nu)                    # corrected: 4.857143
# [1] 4.857143

sum(s2)^2 / sum(s2^2 / nu)     # original WS: 2.666667
# [1] 2.666667

cws(s2, nu) + 2                # Boardman: 6.857143
# [1] 6.857143
```

The original formula gives 2.67, while the corrected formula gives
4.86. The gap illustrates the effect of the $+2$ correction when the
component degrees of freedom are small; as the $\nu_k$ grow, the two
formulas become close in relative terms.

#### A zero variance component

If all but one component are zero, the corrected estimator returns the
degrees of freedom of the remaining component:

```r
cws(c(0, 3), c(10, 20))
# [1] 20
```

This matches the behavior described in the paper: zero components drop
out of both sums, and the result is the degrees of freedom of the
component that remains. The smallest possible value of the estimator is
therefore the smallest component degrees of freedom.

### Applying `cws()` row-wise

`cws()` expects a vector, so for a matrix of variance components use
`apply()`:

```r
set.seed(1)
m <- matrix(rchisq(30, df = 5) / 5, nrow = 10, ncol = 3)
apply(m, 1, function(x) cws(x, nu = rep(5, 3)))
#  [1] 14.51439 14.15320 12.55504 17.61801 15.76236
#  [6] 16.21664 ...
```

Here each row contains three simulated variance components, each based
on $\nu=5$ degrees of freedom.

A typical simulation cell can be written as:

```r
simulate_one <- function(K, nu, nsim = 10000) {
  u <- matrix(rchisq(nsim * K, df = nu), nrow = nsim, ncol = K)
  s2 <- u / nu
  df_est <- apply(s2, 1, function(x) cws(x, nu = rep(nu, K)))
  mean(df_est)
}

simulate_one(K = 4, nu = 2)
```

### Running the simulation script

The repository includes a complete simulation study. Run it from the
repository directory with:

```bash
Rscript simulate_corrected_ws.R
```

The script:

- locates `corrected_welch_satterthwaite.R` relative to its own path,
  so it also works when called from another directory;
- sets the seed to `20260715`;
- uses `K = 2, 4, 6, 8, 10, 20, 30, 40`;
- uses $\nu = 1, 2, 3, 4, 5, 10$;
- runs `nsim = 10000` replications per cell;
- prints a data frame with the target degrees of freedom, mean estimate,
  bias, relative bias, Monte Carlo standard error, and RMSE.

The first rows of the output are:

```text
  K nu target_df mean_df   bias relative_bias   mc_se   rmse
  2  1         2   2.248 0.2479      0.123939 0.01051  1.080
  4  1         4   4.566 0.5663      0.141567 0.01900  1.982
  6  1         6   6.745 0.7453      0.124212 0.02534  2.642
  8  1         8   8.929 0.9291      0.116142 0.03043  3.181
 10  1        10  11.102 1.1019      0.110186 0.03575  3.741
 20  1        20  21.438 1.4376      0.071880 0.05529  5.712
```

The last rows are:

```text
  K nu target_df mean_df   bias relative_bias   mc_se   rmse
 10 10       100 100.690 0.6901      0.006901 0.07023  7.057
 20 10       200 200.600 0.5995      0.002998 0.10788 10.804
 30 10       300 300.313 0.3133      0.001044 0.13538 13.541
 40 10       400 400.846 0.8459      0.002115 0.15586 15.608
```

The full run takes only a few seconds on a typical laptop.

To save the output to a file:

```bash
Rscript simulate_corrected_ws.R > simulation_results.txt
```

If you do not want to commit generated output, add
`simulation_results.txt` to `.gitignore`.

### Input rules and edge cases

- `s2_1`, `s2_2`, `nu_1`, `nu_2`, `w1`, and `w2` must each have length 1.
- `s2` and `nu` must have the same length, and `w` must have that same
  length as well.
- Variances must satisfy `s2 >= 0`; degrees of freedom and weights must
  satisfy `nu > 0` and `w > 0`.
- Inputs must be numeric. Invalid inputs trigger an error via
  `stopifnot()`.
- If the denominator is zero, the function returns `NA_real_`. This can
  happen, for example, when every variance component is zero.
- Use finite degrees of freedom. An infinite `nu` is not rejected by the
  input checks but makes the corresponding denominator term zero.
- The functions return the degrees of freedom, not a $p$-value or a
  confidence interval.
- For a sample variance based on $n$ observations, use $\nu=n-1$.
  For jackknife or balanced repeated replication components, the
  component degrees of freedom are often $\nu_k=1$.

## Related arXiv papers

The following arXiv papers are directly related to the correction and
its applications. The first entry is the primary reference for this
repository.

| arXiv | Title | Year | Relevance |
| --- | --- | --- | --- |
| [2602.20912](https://arxiv.org/abs/2602.20912) | *A Corrected Welch Satterthwaite Equation. And: What You Always Wanted to Know About Kish's Effective Sample but Were Afraid to Ask* — M. von Davier | 2026 | Main reference: derivation of the corrected equation, Kish connection, simulation study, jackknife, multiple imputation, and Welch test applications. |
| [2409.14606](https://arxiv.org/abs/2409.14606) | *A Modified Satterthwaite (1941,1946) Effective Degrees of Freedom Approximation* — M. von Davier | 2024 | Earlier modification for small component degrees of freedom; theoretical basis for the correction. |
| [2503.22080](https://arxiv.org/abs/2503.22080) | *An Improved Satterthwaite Effective Degrees of Freedom Correction for Weighted Syntheses of Variance* — M. von Davier | 2025 | Weighted extension of the correction, including Rubin's total variance in multiple imputation. |
| [2603.11705](https://arxiv.org/abs/2603.11705) | *Effective Degrees of Freedom for Balanced Repeated Replication and Paired Jackknife Variance Estimates: A Unified Approach via Stratum Contrasts* — M. von Davier | 2026 | Derives degrees of freedom for BRR and jackknife variance estimators and connects them to the Welch–Satterthwaite approximation. |
| [2609.17880](https://arxiv.org/abs/2609.17880) | *Explicit Gamma-Stein Bounds for Satterthwaite's Approximation* — G. Bailly, F. Rapin, Y. Swan, R. von Sachs | 2026 | Theoretical error bounds for the Satterthwaite approximation of weighted sums of chi-square and gamma variables. |
| [1808.08986](https://arxiv.org/abs/1808.08986) | *The Behrens–Fisher Problem with Covariates and Baseline Adjustments* — C. Cao, M. Pauly, F. Konietschke | 2018 | Uses the Welch–Satterthwaite $t$-test and generalizes it to covariate adjustment. |
| [2212.13372](https://arxiv.org/abs/2212.13372) | *Two-sample Behrens–Fisher Problems for High-Dimensional Data: A Normal Reference F-Type Test* — T. Zhu, P. Wang, J.-T. Zhang | 2022 | Uses the Welch–Satterthwaite chi-square approximation in a high-dimensional testing setting. |

### Suggested reading order

1. **arXiv:2602.20912** — the corrected equation used in this repository.
2. **arXiv:2503.22080** — the weighted extension and multiple-imputation
   application.
3. **arXiv:2409.14606** — the earlier small-$\nu$ modification.
4. **arXiv:2603.11705** — degrees of freedom for BRR and jackknife
   variance estimates.
5. **arXiv:2609.17880** — theoretical approximation guarantees.

### BibTeX for the main paper

```bibtex
@misc{vonDavier2026cws,
  author       = {von Davier, Matthias},
  title        = {A Corrected {Welch} {Satterthwaite} Equation. And: What
                  You Always Wanted to Know About {Kish}'s Effective
                  Sample but Were Afraid to Ask},
  year         = {2026},
  eprint       = {2602.20912},
  archivePrefix = {arXiv},
  primaryClass = {stat.AP},
  doi          = {10.48550/arXiv.2602.20912},
  url          = {https://arxiv.org/abs/2602.20912}
}
```

## Classical references

- Boardman, T. J. (1974). Confidence Intervals for Variance Components
  – A Comparative Monte Carlo Study. *Biometrics*, 30(2), 251–262.
  doi:[10.2307/2529647](https://doi.org/10.2307/2529647)
- Johnson, E. G., & Rust, K. F. (1992). Population Inferences and
  Variance Estimation for NAEP Data. *Journal of Educational Statistics*,
  17(2), 175–190. doi:[10.2307/1165168](https://doi.org/10.2307/1165168)
- Kish, L. (1965). *Survey Sampling*. John Wiley & Sons.
- Kish, L. (1992). Weighting for Unequal Pi. *Journal of Official
  Statistics*, 8(2), 183–200.
- Lipsitz, S., Parzen, M., & Zhao, L. P. (2002). A Degrees-of-Freedom
  Approximation in Multiple Imputation. *Journal of Statistical
  Computation and Simulation*, 72(4), 309–318.
  doi:[10.1080/00949650212848](https://doi.org/10.1080/00949650212848)
- Rubin, D. B., & Schenker, N. (1986). Multiple Imputation for Interval
  Estimation From Simple Random Samples With Ignorable Nonresponse.
  *Journal of the American Statistical Association*, 81(394), 366–374.
  doi:[10.2307/2289225](https://doi.org/10.2307/2289225)
- Satterthwaite, F. E. (1941). Synthesis of Variance. *Psychometrika*,
  6(5), 309–316. doi:[10.1007/BF02288586](https://doi.org/10.1007/BF02288586)
- Satterthwaite, F. E. (1946). An Approximate Distribution of Estimates
  of Variance Components. *Biometrics Bulletin*, 2(6), 110–114.
  doi:[10.2307/3002019](https://doi.org/10.2307/3002019)
- von Davier, M. (2025). An Improved Satterthwaite (1941, 1946)
  Effective df Approximation. *Journal of Educational and Behavioral
  Statistics*. doi:[10.3102/10769986241309329](https://doi.org/10.3102/10769986241309329)
- Welch, B. L. (1947). The Generalization of "Student's" Problem When
  Several Different Population Variances Are Involved. *Biometrika*,
  34(1/2), 28–35. doi:[10.2307/2332510](https://doi.org/10.2307/2332510)

## Committing the documentation to GitHub

The repository is already a Git clone with a configured `origin`
remote. From the repository directory:

```bash
cd /home/matthias/python/Corrected-Welch-Satterthwaite
```

### 1. Check what changed

```bash
git status
```

`Documentation.md` should appear as an untracked file. The LaTeX source
`short_satter_07152026.tex` may also appear as untracked if it has not
been committed yet.

### 2. Stage the documentation

```bash
git add Documentation.md
```

If you also want to commit the LaTeX source in the same commit:

```bash
git add Documentation.md short_satter_07152026.tex
```

Review the staged changes:

```bash
git diff --cached
```

### 3. Commit

```bash
git commit -m "Add documentation for corrected Welch-Satterthwaite functions"
```

### 4. Push to GitHub

```bash
git push origin main
```

### One-line version

```bash
git add Documentation.md && git commit -m "Add documentation for corrected Welch-Satterthwaite functions" && git push origin main
```

### If Git asks who you are

Configure your name and email once:

```bash
git config --global user.name "Your Name"
git config --global user.email "you@example.com"
```

### If the push is rejected

If GitHub reports that the remote contains work you do not have locally,
rebase on the remote branch and push again:

```bash
git pull --rebase origin main
git push origin main
```

### Authentication

If `git push` asks for a username and password, use a GitHub personal
access token (PAT) as the password, or authenticate with the GitHub CLI:

```bash
gh auth login
git push origin main
```

Prefer a credential helper, the GitHub CLI, or an SSH remote over putting
a token directly in the remote URL. If a token has been embedded in a
remote URL, treat it as exposed and rotate it in GitHub settings. You can
list the configured remote names without printing their URLs with:

```bash
git remote
```

If you need to inspect the URL, `git remote get-url origin` prints it;
avoid sharing that output when a token is embedded in the URL.

### Verify the result

After pushing, open the repository on GitHub and confirm that
`Documentation.md` appears on the `main` branch. You can also check
locally:

```bash
git log --oneline -3
git status
```
