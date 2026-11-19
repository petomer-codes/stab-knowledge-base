---
title: Direct Estimators
layout: default
parent: Survey Theory
nav_order: 1
---

# Introduction

In many applications, the target variable is the areas’ population mean,
$p_i$, for $i=1,\ldots,m$ areas. For instance, the under-five child
mortality rate (U5MR). Under a finite population framework (see *\[link
to finite population section\]* for more information), this is a fixed
quantity: if we observed all of the children, we could calculate it
directly, but since we have only a survey sample, we need to estimate
it.

A direct estimator is an estimator derived directly from the survey,
without any model assumption. This approach is called “design-based”,
where the randomness comes from which units are selected. There are many
versions of direct estimators; here we focus on the two most commonly
used.

# Estimators

Let $U$ be a finite population of size $N$, partitioned into $m$
disjoint subpopulations $U = U_1 \cup \cdots \cup U_m$, where
$|U_i|=N_i$ and $\sum_{i=1}^m N_i = N$. Each unit $j \in U$ has an
associated response $y_{ij}$, and the quantity of interest for area $i$
is the area mean
$$p_i = \frac{1}{N_i}\sum_{j \in U_i} y_{ij}, \quad i = 1, \ldots, m.$$

Let $S \subset U$ be a sample of $n$ units, partitioned into
$S = S_1 \cup \cdots \cup S_m$, with $|S_i|=n_i$. Each unit $j$ is
included in the sample with probability $\pi_{ij} = P(j \in S_i)$ and
has a corresponding sampling weight $w_{ij} = 1/\pi_{ij}$. Under this
framework, $y_{ij}$ is treated as fixed and randomness enters only
through the sampling mechanism.

## Horvitz–Thompson estimator

The Horvitz–Thompson (HT) estimator (Horvitz and Thompson 1952) is given
by:
$$\hat{p}_i^{\text{HT}} = \frac{\sum_{j \in S_i} w_{ij} y_{ij}}{N_i}, \quad i = 1, \ldots, m.$$

We can view the multiplication by weights as a correction for unequal
representation in the sample: a sampled unit $j$ with inclusion
probability $\pi_{ij}$ “stands in” for $w_{ij} = 1/\pi_{ij}$ units of
the population. For example, if a unit represents $\pi_{ij}=1/10$ of the
population, then $w_{ij}=10$, and its response is counted ten times so
that it accounts for the units that were not sampled. A more detailed
discussion of how the weights are constructed is given in the *\[weights
section\]*.

A limitation of the HT estimator is that it normalizes the estimator of
the population total, $\sum_{j \in S_i} w_{ij} y_{ij}$, by the known
population size $N_i$, which implicitly assumes that
$N_i = \sum_{j \in S_i} w_{ij}$. This equality is not guaranteed: the
realized sum of weights $\sum_{j \in S_i} w_{ij}$ is itself a random
quantity. Dividing by $N_i$ therefore ignores that randomness, which can
make the estimator unstable. This may happen when the inclusion
probabilities are unevenly assigned within an area, or when the number
of sampled units $n_i$ in area $i$ is random.

## Hájek estimator

The Hájek estimator (Hájek 1971) addresses this by normalizing the
estimated population total with the estimated population size,
$\hat{N}_i=\sum_{j\in S_i}w_{ij}$, in place of the known $N_i$:
$$\hat{p}_i^{\text{Hájek}} = \frac{\sum_{j \in S_i} w_{ij} y_{ij}}{\hat N_i}, \quad i = 1, \ldots, m.$$

Since the Hájek estimator is used throughout, we denote it simply by
$\hat{p}_i^{w}$, i.e., $\hat{p}_i^{\text{Hájek}} = \hat{p}_i^{w}$.

# Expectation and Variance

## HT estimator

### Expectation

Let $I_{ij} = \mathbf{1}[j \in S_i]$, indicating whether unit $j$ is
present in the sample. For any unit, the expectation of its indicator is
exactly the inclusion probability:
$$\pi_{ij} = \mathbb{E}(I_{ij}) = \Pr(j \in S_i).$$

Taking the expectation of the estimator, we can show that it is
unbiased: $$\begin{aligned}
  \mathbb{E}\left(\hat{p}_i^{\text{HT}}\right)
  &= \mathbb{E}\left(\frac{1}{N_i}\sum_{j \in S_i} \frac{y_{ij}}{\pi_{ij}}\right)
   = \mathbb{E}\left(\frac{1}{N_i}\sum_{j \in U_i} \frac{y_{ij}}{\pi_{ij}}\, I_{ij}\right) \\
  &= \frac{1}{N_i}\sum_{j \in U_i} \frac{y_{ij}}{\pi_{ij}}\, \mathbb{E}(I_{ij})
   = \frac{1}{N_i}\sum_{j \in U_i} \frac{y_{ij}}{\pi_{ij}}\, \pi_{ij} \\
  &= \frac{1}{N_i}\sum_{j \in U_i} y_{ij} = p_i.
\end{aligned}$$

### Variance

Let $\pi_{i,jk} = \mathbb{E}(I_{ij}I_{ik}) = \Pr(j,k \in S_i)$ denote
the joint inclusion probability of units $j$ and $k$ in area $i$, with
$\pi_{i,jj} = \pi_{ij}$. Then
$\mathrm{Cov}(I_{ij}, I_{ik}) = \pi_{i,jk} - \pi_{ij}\pi_{ik}$, and
since $y_{ij}$ and $\pi_{ij}$ are fixed, $$\begin{aligned}
  \mathrm{Var}\left(\hat{p}_i^{\text{HT}}\right)
  &= \frac{1}{N_i^2}\,
     \mathrm{Var}\left(\sum_{j \in U_i} \frac{y_{ij}}{\pi_{ij}}\, I_{ij}\right)
   = \frac{1}{N_i^2} \sum_{j \in U_i}\sum_{k \in U_i}
     \frac{y_{ij}}{\pi_{ij}}\frac{y_{ik}}{\pi_{ik}}
     \left(\pi_{i,jk} - \pi_{ij}\pi_{ik}\right) \\
  &= \frac{1}{N_i^2} \left[
     \sum_{j \in U_i} \frac{1-\pi_{ij}}{\pi_{ij}}\, y_{ij}^2
     + \sum_{j \in U_i}\sum_{k \neq j}
     \frac{\pi_{i,jk} - \pi_{ij}\pi_{ik}}{\pi_{ij}\pi_{ik}}\, y_{ij}y_{ik}
     \right].
\end{aligned}$$

Since this variance depends on the values $y_{ij}$ of all units in
$U_i$, including those not in the sample, it cannot be computed from the
data. An unbiased estimator (Rao and Molina 2015, sec. 2.3.1) is
obtained by: $$\hat{V}\left(\hat{p}_i^{\text{HT}}\right)
  = -\frac{1}{2N_i^2}\sum_{j\in S_i}\sum_{k\in S_i,j\not=k}\frac{\pi_{i,jk}-\pi_{ij}\pi_{ik}}{\pi_{i,jk}}\left(\frac{y_{ij}}{\pi_{ij}}-\frac{y_{ik}}{\pi_{ik}}\right)^2.$$

## Hájek estimator

The Hájek estimator is a ratio of two estimated totals,
$\hat{p}_i^{w} = \hat{T}_i / \hat{N}_i$, with
$\hat{T}_i = \sum_{j \in S_i} w_{ij} y_{ij}$ and
$\hat{N}_i = \sum_{j \in S_i} w_{ij}$.

### Expectation

*(to be added)*

### Variance

A Taylor linearization variance estimator (Rao and Molina 2015, sec.
2.4.1) of $\hat{p}_i^{w}$ is obtained from the HT variance estimator
above by replacing $y_{ij}$ with the residual
$\tilde e_{ij} = y_{ij} - \hat{p}_i^{w}$ and $N_i$ with $\hat N_i$:
$$\hat{V}\left(\hat{p}_i^{w}\right)
  = -\frac{1}{2\hat{N}_i^{2}}\sum_{j\in S_i}\sum_{k\in S_i,\,k\neq j}
    \frac{\pi_{i,jk}-\pi_{ij}\pi_{ik}}{\pi_{i,jk}}
    \left(\frac{\tilde e_{ij}}{\pi_{ij}}-\frac{\tilde e_{ik}}{\pi_{ik}}\right)^2 .$$

# Design Consistency

An estimator is design consistent if, as the area sample size goes to
infinity, it converges in probability to the true area mean under the
sampling design alone, without assuming any model.

Under the regularity condition stated in the *\[finite population
section\]*, both the HT and Hájek estimators are design consistent,
since $$\begin{aligned}
    \mathrm{bias}\left(\hat{p}_i^{\text{HT}}\right) = 0, \qquad & \mathrm{Var}\left(\hat{p}_i^{\text{HT}}\right) = O(N_i^{-1}) \\
    \mathrm{bias}\left(\hat{p}_i^{w}\right) = O(N_i^{-1}), \qquad & \mathrm{Var}\left(\hat{p}_i^{w}\right) = O(N_i^{-1}),
\end{aligned}$$ and hence both are design consistent.

# Limitations

In small area estimation, many areas may have a small or nonexistent
sample size. Consequently, the variance of $\hat{p}_i^w$ inflates,
resulting in an unstable estimator. An extreme case arises when
$\hat{p}_i^w=0$ or $1$, or when it is undefined because $n_i=0$ (see the
*\[variance fix section\]* for more detail). Thus, we need a model that
can address these issues, which leads us to the model-based approach in
the *\[model-based section\]*.

# Code

## Overview

The code below computes direct (design-based) estimates of the neonatal
mortality rate (NMR) for Nigeria at the admin-1 level using the DHS 2024
survey, in two ways:

1.  **`surveyPrev::directEST()`**, the pipeline function used in the
    main project.
2.  **`srvyr` / `survey`**, the usual way of computing design-based
    estimates in R, as an alternative to `surveyPrev`.

`survey_mean()` computes both the Hájek estimator and its variance. The
HT estimator is not computed here because it needs the known area size
$N_i$, which DHS does not provide.

``` r
library(srvyr)
library(dplyr)
library(ggplot2)
library(surveyPrev)
```

## 1. Data

### 1.1 Settings

``` r
indicator_label <- "nmr"
year            <- 2024
country         <- "Nigeria"
cc              <- "NGA"   # ISO3 code,
```

### 1.2 DHS survey data

`getDHSindicator()` returns one row per birth, with the cluster,
stratum, sampling weight, and the 0/1 neonatal death indicator
(`value`).

``` r
dhsData <- getDHSdata(country = country, indicator = "CM_ECMR_C_NNR", year = year)
data    <- getDHSindicator(dhsData, indicator = "CM_ECMR_C_NNR", FUN = surveyPrev::NMR)
```

``` r
head(data)
```

    ##   cluster householdID         v022         v023       v024 weight strata value
    ## 1       1           2 sokoto urban sokoto urban north west 491963  urban     0
    ## 2       1           2 sokoto urban sokoto urban north west 491963  urban     0
    ## 3       1           2 sokoto urban sokoto urban north west 491963  urban     0
    ## 4       1          10 sokoto urban sokoto urban north west 491963  urban     0
    ## 5       1          10 sokoto urban sokoto urban north west 491963  urban     0
    ## 6       1          13 sokoto urban sokoto urban north west 491963  urban     0

### 1.3 Geography and cluster-to-area mapping

`clusterInfo()` assigns each DHS cluster (via its GPS location) to an
admin area.

``` r
geo       <- getDHSgeo(country = country, year = year)
poly.adm1 <- sf::st_as_sf(geodata::gadm(country = cc, level = 1, path = tempdir()))
poly.adm2 <- sf::st_as_sf(geodata::gadm(country = cc, level = 2, path = tempdir()))
poly.adm2$admin2.name.full <- paste0(poly.adm2$NAME_1, "_", poly.adm2$NAME_2)

cluster.info <- clusterInfo(geo = geo, poly.adm1 = poly.adm1, poly.adm2 = poly.adm2,
                            by.adm1 = "NAME_1", by.adm2 = "NAME_2")
```

## 2. Direct estimates with `surveyPrev`

We work at the admin-1 level, where each area has enough clusters and
events that the variance fix for degenerate areas is not needed (see the
[variance fix section](#var-fix)), so `var.fix = FALSE`.

``` r
direct_estimates <- directEST(data = data, cluster.info = cluster.info,
                              admin = 1, var.fix = FALSE)

head(direct_estimates$res.admin1)
```

    ##   admin1.name direct.est   direct.var direct.logit.est direct.logit.var direct.logit.prec   direct.se
    ## 1        Abia 0.03587090 1.244342e-04        -3.291299       0.10403634          9.612026 0.011155009
    ## 2     Adamawa 0.04248869 3.300865e-05        -3.115100       0.01994313         50.142583 0.005745316
    ## 3   Akwa Ibom 0.03878761 1.116606e-04        -3.210095       0.08032949         12.448728 0.010566956
    ## 4     Anambra 0.03168507 4.387149e-05        -3.419712       0.04660578         21.456568 0.006623556
    ## 5      Bauchi 0.04779752 7.504430e-05        -2.991804       0.03622834         27.602704 0.008662811
    ## 6     Bayelsa 0.03848557 1.030346e-04        -3.218226       0.07524470         13.289973 0.010150598
    ##   direct.lower direct.upper        cv
    ## 1   0.01938883   0.06542918 0.3109766
    ## 2   0.03255001   0.05528858 0.1352199
    ## 3   0.02262983   0.06570651 0.2724312
    ## 4   0.02098301   0.04758028 0.2090434
    ## 5   0.03341196   0.06794144 0.1812398
    ## 6   0.02284620   0.06412827 0.2637508

## 3. Direct estimates with `srvyr`

### 3.1 Attach admin-1 names

Each birth is linked to its admin-1 area through its cluster. Births in
clusters that could not be assigned to an area are dropped.

``` r
dat <- left_join(data,
                 cluster.info$data[, c("cluster", "admin1.name")],
                 by = "cluster") %>%
  filter(!is.na(admin1.name))
```

### 3.2 Survey design

Clusters (PSUs) nested within strata, with sampling weights. The design
is built on the whole sample, and areas are handled as domains with
`group_by()`.

``` r
# "adjust" centers strata with a single cluster at the grand mean instead of failing
options(survey.lonely.psu = "adjust")

des <- dat %>%
  as_survey_design(ids = cluster, strata = strata, weights = weight, nest = TRUE)
```

### 3.3 Estimates per admin-1 area

`survey_mean()` returns the Hájek estimator and its variance in each
area.

``` r
hajek <- des %>%
  group_by(admin1.name) %>%
  summarise(p_hajek = survey_mean(value, vartype = "var")) %>%
  rename(var_hajek = p_hajek_var)

head(hajek)
```

    ## # A tibble: 6 × 3
    ##   admin1.name p_hajek var_hajek
    ##   <chr>         <dbl>     <dbl>
    ## 1 Abia         0.0359 0.000122 
    ## 2 Adamawa      0.0425 0.0000314
    ## 3 Akwa Ibom    0.0388 0.000108 
    ## 4 Anambra      0.0317 0.0000419
    ## 5 Bauchi       0.0478 0.0000714
    ## 6 Bayelsa      0.0385 0.0000982

## 4. Comparison

``` r
sp <- direct_estimates$res.admin1 %>%
  transmute(admin1.name,
            p_sp   = direct.est,
            var_sp = direct.var)

cmp <- inner_join(hajek, sp, by = "admin1.name")
```

``` r
ggplot(cmp, aes(x = p_sp, y = p_hajek)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point() +
  labs(x = "surveyPrev (directEST)", y = "srvyr (survey_mean)",
       title = "Admin-1 NMR direct estimates") +
  theme_minimal()
```

![](direct-estimators_files/figure-gfm/compare-estimates-1.png)<!-- -->

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-hajek1971discussion" class="csl-entry">

Hájek, Jaroslav. 1971. “Discussion of, ‘an Essay on the Logical
Foundations of Survey Sampling, Part i,’ by D. Basu.” In *Foundations of
Statistical Inference*, edited by V. P. Godambe and D. A. Sprott. Holt,
Rinehart; Winston.

</div>

<div id="ref-HT" class="csl-entry">

Horvitz, D. G., and D. J. Thompson. 1952. “A Generalization of Sampling
Without Replacement from a Finite Universe.” *Journal of the American
Statistical Association* 47 (260): 663–85.
<http://www.jstor.org/stable/2280784>.

</div>

<div id="ref-smebook" class="csl-entry">

Rao, J. N. K., and I. Molina. 2015. *Direct Domain Estimation*. John
Wiley & Sons, Ltd.
https://doi.org/<https://doi.org/10.1002/9781118735855.ch2>.

</div>

</div>
