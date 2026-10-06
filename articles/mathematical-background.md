# Mathematical Background

## Introduction

Equivalence testing first appeared in the field of pharmacokinetics
Metzler (1974), and finds its most common application in the context of
generic drug manufacturing Senn (2021) where the aim is to determine
whether a parameter like the mean variation of treatment response
between a reference and a generic drug falls within a predetermined
equivalence range that we denote $`(-c,c)`$, indicating that drugs are
equivalent. Deviations within this region would be considered
insignificant as they are too small to reflect the dissimilarity of the
therapeutic effects between the compared treatments. In contrast to
standard hypothesis testing for equality of means, where the null
hypothesis assumes that both means are equal, and the alternative
assumes they are not, equivalence testing reverses the roles of the
hypothesis formulations and considers a region rather than a point.
Specifically, it defines the alternative as the equivalence region
within which the parameter of interest must lie for the treatments to be
considered equivalent, and the null hypothesis as the opposite. This
paradigm puts the burden of proof on equivalence, rather than
non-equality, and emphasizes the importance of assessing the similarity
of treatments in addition to their differences.

## Mathematical setup

More formally, the hypotheses of interest are defined as

``` math
\begin{equation*}
    H_0: \; \theta \not\in \Theta_1 \quad vs. \quad H_1: \; \theta \in \Theta_1 := (\theta_L, \, \theta_U),
\end{equation*}
```

where $`\Theta_1=(\theta_L,\theta_U)`$ denotes the range of equivalence
whose limits are known constants, and are usually symmetrical so we
define $`c:=\theta_U=-\theta_L`$.

Existing approaches, such as the state-of-the-art Two One-Sided Tests
(TOST) Schuirmann (1987), have demonstrated conservativeness and
decreased power, particularly for highly variable responses. This can be
demonstrated by computing the power function which would be evaluated at
different points of the parameter space. In particular, evaluating it at
the equivalence bounds yields the size and shows that it is strictly
bounded by the significance level $`\alpha`$. To see this more clearly,
consider a canonical form for the average equivalence problem in the
univariate framework, which consists of two independent random variables
$`\widehat{\theta}`$ and $`\widehat{\sigma}_{\nu}`$ with distributions

``` math
\begin{equation*}
    \widehat{\theta} \sim N\left(\theta, \sigma_{\nu}^2\right), \quad \text{and} \quad \nu\frac{\widehat{\sigma}_{\nu}^2}{\sigma_{\nu}^2} \sim \chi^2_{\nu},
\end{equation*}
```
where $`\theta`$ is the parameter of interest, $`\sigma_{\nu}`$ is the
standard error and $`\nu`$ are the degrees of freedom.

The TOST is based on two test statistics

``` math
\begin{equation*}
T_L := \frac{\widehat{\theta}+c}{\widehat{\sigma}_{\nu}}\sim t_\nu\left(\frac{\theta+c}{\sigma_\nu}\right), \quad\text{and}\quad T_U := \frac{\widehat{\theta}-c}{\widehat{\sigma}_{\nu}}\sim t_\nu\left(\frac{\theta-c}{\sigma_\nu}\right).
\end{equation*}
```
$`\text{H}_{0}`$ is rejected in favor of equivalence if both tests
simultaneously reject their marginal null hypotheses, i.e.,

``` math
\begin{equation*}
     T_L \geq t_{1-\alpha,\nu}, \quad \text{ and} \quad T_U \leq -t_{1-\alpha,\nu}.
\end{equation*}
```

By rearranging the terms, we can define a rejection region for the TOST:
``` math
\begin{equation*}
    C_1 := \left\{\widehat{\theta} \in \mathbb{R},\, \widehat{\sigma}_{\nu} \in \mathbb{R}_+ \,\Big\vert \,  \vert\widehat{\theta}\vert  \leq  c - t_{1-\alpha,\nu} \hat{\sigma}_\nu  \right\}.
\end{equation*}
```

Given $`\alpha`$, $`\theta`$, $`\sigma_{\nu}`$, $`\nu`$ and $`c`$, the
power function of the TOST corresponds to the probability of rejecting
$`H_0`$, i.e., the integral of the joint density of
$`(\hat{\theta},\hat{\sigma}_\nu)`$ over the rejection area $`C_1`$ (see
also Phillips (1990)), that is
``` math
\begin{equation*}
    \begin{aligned}
        p(\alpha, \theta, \sigma_{\nu}, \nu, c) :=& \Pr\left(T_L \geq t_{1-\alpha,\nu}\; \text{{ and}} \;  T_U \leq -t_{1-\alpha,\nu}\,  \big|\;  \alpha, \theta, \sigma_{\nu}, \nu, c \right)\\
        =&\int_0^\infty I(\widehat\sigma_\nu t_{1-\alpha,\nu} < c)\left\{\Phi\left(\frac{c+t_{\alpha,\nu} \widehat\sigma_\nu-\theta}{\sigma_\nu}\right)
    -\Phi\left(-\frac{c+ t_{\alpha,\nu} \widehat\sigma_\nu+\theta}{\sigma_\nu}
    \right) \right\} \\
        \phantom{=}&\times f_W ( \widehat\sigma_\nu \lvert \sigma_\nu,\nu)d\widehat\sigma_\nu.
    \end{aligned}
\end{equation*}
```
Noting that the vector $`(T_L,T_U)`$ has a bivariate non-central
$`t`$-distribution with non-centrality parameters
$`\frac{\theta-c}{\sigma_\nu}`$ and $`\frac{\theta+c}{\sigma_\nu}`$
respectively, we can express the power function in terms of Owen’s
$`Q`$-function Owen (1965):
``` math
\begin{equation*}
    \begin{aligned}
        p(\alpha, \theta, \sigma_{\nu}, \nu, c) :=& Q_{\nu}\left(-t_{1-\alpha,\nu},\, \frac{\theta - c}{\sigma_{\nu}},\, \frac{c \sqrt{\nu}}{\sigma_{\nu} t_{1-\alpha,\nu}} \right) - Q_{\nu}\left(t_{1-\alpha,\nu},\, \frac{\theta + c}{\sigma_{\nu}},\, \frac{c \sqrt{\nu}}{\sigma_{\nu} t_{1-\alpha,\nu}} \right).
    \end{aligned}
\end{equation*}
```

## Conservativeness of the TOST

This formulation allows to demonstrate that, for $`\sigma_{\nu}>0`$, the
TOST is not size-$`\alpha`$ in finite samples as
``` math
\begin{equation*}
    \begin{aligned}
    \omega_\nu(\alpha, c, \sigma_{\nu}) &:= \sup_{\theta \, \in \, \Theta_0}\; p(\alpha, \theta, \sigma_{\nu}, \nu, c) \\
    &< Q_{\nu}\left(-t_{1-\alpha,\nu},\, 0,\, \frac{c \sqrt{\nu}}{\sigma_{\nu} t_{1-\alpha,\nu}} \right) < \Pr\left(T_{\nu} \leq -t_{1-\alpha,\nu}\right) = \alpha,
    \end{aligned}
\end{equation*}
```
where $`\Theta_0:=\mathbb{R}\setminus(-c,c)`$ corresponds to parameter
space under the null.

## Finite sample corrections to the TOST

As $`\omega_\nu(\gamma, \delta, \sigma_{\nu})`$ is continuously
differentiable and strictly increasing in $`\gamma`$ and $`\delta`$,
solving the following matching paradigms represents a finite sample
correction to the TOST and would ensure size-$`\alpha`$ tests:
``` math
\alpha^*:=\text{argzero}_{\gamma \in [\alpha, 0.5)} \;\big[ \omega_\nu(\gamma, c, \sigma_{\nu}) - \alpha\big],
```
``` math
c^*:=\text{argzero}_{\delta \in [c, \infty)} \;\big[ \omega_\nu(\alpha, \delta, \sigma_{\nu}) - \alpha\big].
```

The conditions under which the solutions are singletons and provide a
uniformly more powerful inference are studied in length in Boulaguiem et
al. (2023) and compared to existing methods.

## References

Boulaguiem, Y., J. Quartier, M. Lapteva, et al. 2023. “Finite Sample
Adjustments for Average Equivalence Testing.” *Statistics in Medicine*,
2023–03.

Metzler, C. M. 1974. “Bioavailability – a Problem in Equivalence.”
*Journal of Pharmaceutical Sciences* 30: 309–17.

Owen, D. B. 1965. “A special case of a bivariate non-central
$`t`$-distribution.” *Biometrika* 52: 437–46.

Phillips, K. 1990. “Power of the Two One-Sided Tests Procedure in
Bioequivalence.” *Journal of Pharmacokinetics and Biopharmaceutics* 18:
137–44.

Schuirmann, D. J. 1987. “A Comparison of the Two One-Sided Tests
Procedure and the Power Approach for Assessing the Equivalence of
Average Bioavailability.” *Journal of Pharmacokinetics and
Biopharmaceutics* 15 (6): 657–80.

Senn, S. 2021. *Statistical Issues in Drug Development, 3rd Edition*.
Wiley.
