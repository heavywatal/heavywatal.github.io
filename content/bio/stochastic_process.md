+++
title = 'Stochastic Process'
tags = ["genetics", "math"]
[menu.main]
  parent = "bio"
+++

## Geometric distribution 幾何分布

確率 *p* のBernoulli試行で *k* 回目に初めて当たる確率。
まず外し続ける必要がある。

<div>$$\begin{split}
\mathrm{Prob}[X = k] &= p(1 - p)^{k - 1} \\
\mathrm E[X] &= \frac 1 p \\
\mathrm{Var}[X] &= \frac {1 - p} p
\end{split}$$</div>

## Negative binomial distribution 負の二項分布

確率 *p* のBernoulli試行で *r* 回目の当たりが *k* 回目に出る確率。
幾何分布を *r* 回畳み込んだもの。

<div>$$\begin{split}
\mathrm{Prob}[X = k] &= {k - 1 \choose r - 1} p^r (1 - p)^{k - r} \\
\mathrm E[X] &= \frac r p \\
\mathrm{Var}[X] &= \frac {r(1 - p)} {p^2}
\end{split}$$</div>

## Exponential distribution 指数分布

幾何分布の連続時間バージョン。
時間あたり $\lambda$ 回起こるPoisson過程で時間 $t$ まで起こらない確率 $Q(t)$

<div>$$\begin{split}
Q(t) &= \lim_{\mathrm dt \to 0} (1 - \lambda \mathrm dt) ^ {t/ \mathrm d t}\\
     &= \lim_{n \to \infty} (1 - \frac {\lambda t} n) ^ n\\
     &= e ^ {-\lambda t}
\end{split}$$</div>

{{%div class="note"%}}
ネイピア数

<div>$$\begin{split}
\lim_{x \to \infty} \left(1 + \frac 1 x \right)^x &= e\\
\lim_{x \to \infty} \left(1 + \frac a x \right)^x
   &= \lim_{y \to \infty} \left(1 + \frac 1 y \right)^{ay}\\
   &= e^a
\end{split}$$</div>
{{%/div%}}

{{%div class="note"%}}
別の考え方

<div>$$\begin{split}
Q(t + \mathrm dt) &= Q(t) (1 - \lambda \mathrm dt)\\
Q(t + \mathrm dt) - Q(t) &= -\lambda Q(t) \mathrm dt\\
\frac {\mathrm dQ(t)} {\mathrm dt} &= -\lambda Q(t)\\
\frac 1 {Q(t)} \mathrm dQ(t) &= -\lambda \mathrm dt\\
\ln Q(t) &= -\lambda t\\
Q(t) &= e ^ {-\lambda t}
\end{split}$$</div>
{{%/div%}}

$t$ らへんで初めて起こる確率は、
$t$ まで起こらず次の瞬間に起こる確率
$Q(t) \lambda \mathrm dt = \lambda e ^ {-\lambda t} \mathrm dt$
であり、それを全区間たすと1になる

<div>$$\begin{split}
\int _0^\infty \lambda e ^ {-\lambda t} \mathrm dt
   = -e ^ {\lambda t} \Big | _0^\infty
   = 1
\end{split}$$</div>

よって、最初の事象が起こるまでの待ち時間 $t$ の分布、
すなわち確率密度関数 (PDF: probability density function) は

<div>$$\begin{split}
f(t; \lambda) &= \lambda e ^ {-\lambda t} \\
\mathrm E[t] &= \frac 1 \lambda \\
\mathrm{Var}[t] &= \frac 1 {\lambda^2}
\end{split}$$</div>

また、$t$ までに最低1回は起こる確率、
すなわち累積分布関数 (CDF: cumulative distribution function) は

<div>$$\begin{split}
F(t; \lambda)
   &= \int _0^t f(x, \lambda) \mathrm dx
   &= \int _0^t \lambda e ^ {-\lambda x} \mathrm dx\\
   &= -e ^ {-\lambda x} \Big | _0^t\\
   &= 1 - e ^ {-\lambda t}
\end{split}$$</div>

これは起こらない確率 $Q(t)$ を1から引いたものに等しい。

## Gamma distribution

パラメータの取り方によって書き方がいくつかある。
shape *k* (&gt;0), scale *θ* (&gt;0) を使った一般的な表し方では、
平均待ち時間 *θ* の事象が *k* 回起こるまでの待ち時間の分布と見なせる。

<div>$$\begin{split}
f(x; k, \theta) &= \frac {x^{k-1} \exp[-\frac x \theta]} {\theta^k \Gamma(k)} \\
\mathrm E[x] &= k\theta \\
\mathrm{Var}[x] &= k\theta^2 \\
\end{split}$$</div>

rate parameter $\lambda = \frac 1 \theta$ を使った形だと導出もイメージしやすい。
時間当たり *λ* 回起こるPoisson過程で *k* 回起こるまでの待ち時間の分布、
すなわち指数分布を *k* 個畳み込んだ分布。
cf. [/lectures/wakeley-2-2]({{< relref "lectures/wakeley-2-2.md" >}})

<div>$$\begin{split}
f(t; k, \lambda) &= \lambda e^{-\lambda t}
                \frac {(\lambda t)^{k-1}}
                      {\Gamma(k)} \\
\mathrm E[t] &= \frac k \lambda \\
\mathrm{Var}[t] &= \frac k {\lambda^2} \\
\end{split}$$</div>

内部的に *k* ステップからなる過程が時間 $\mu = \frac k \lambda = k\theta$ で完了する。

<div>$$\begin{split}
f(t; k, \mu) &= e^{-\left(\frac k \mu t \right)}
                \frac {k^k t^{k-1}}
                      {\mu^k \Gamma(k)} \\
\mathrm E[t] &= \mu \\
\mathrm{Var}[t] &= \frac {\mu^2} k \\
\end{split}$$</div>

## Weibull distribution

指数関数がmemerylessな待ち時間なのに対して、
こちらはある時間まで事象が起こらなかったという記憶ありの待ち時間

<div>$$\begin{split}
f(t; k, \lambda) &= \frac k \lambda \left(\frac t \lambda \right)^{k - 1}
                    \exp[-\left(\frac t \lambda \right)^k] \\
F(t; k, \lambda) &= 1 - \exp[-\left(\frac t \lambda \right)^k]
\end{split}$$</div>

-   Scale parameter: $\lambda$
-   Shape parameter: *k*
