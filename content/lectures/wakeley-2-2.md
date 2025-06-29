+++
title = 'Wakeley輪読会 2章2節'
tags = ["genetics", "book"]
[menu.main]
  parent = "lectures"
+++

Book
:   [Coalescent Theory --- An Introduction](https://www.amazon.co.jp/dp/0974707759?&linkCode=ll1&tag=heavywatal-22&linkId=a938e9af55ab8c42f36704afa4824abb)

Author
:   [John Wakeley](http://www.oeb.harvard.edu/faculty/wakeley/wakeleylab.htm)

Publisher
:   [Roberts & Company](http://www.roberts-publishers.com/authors/wakeley-john/coalescent-theory.html)

Errata
:   [PDF](http://www.oeb.harvard.edu/faculty/wakeley/John/Reprints/CorrectionsWakeleyFall2010.pdf)

輪読担当
:   岩嵜航

日程
:   2015-03-20

## 2. Probability Theory

### 2.1 Fundamentals

#### 2.1.1 Events, Probabilities, and Random Variables

#### 2.1.2 Probability Distributions

### 2.2 Poisson Processes

The backbone of the neutral coalescent:

Poisson distribution
:   the number of events that occur over a fixed period of time

Exponential distribution
:   waiting time until a first event occurs

------------------------------------------------------------------------

The poisson process is a counting process.

<div>\[\begin{aligned}
P[K(t) = 1]   &= \lambda t + o(t) \\
P[K(t) \ge 2] &= o(t)
\end{aligned}\]</div>

where
:   $K(t)$ is the number of observed events before $t$. $K(0) = 0$\
    $\lambda$ is the rate of occurrence per unit time\
    $o(t)$ goes to faster than $t$

These implies, within a sufficiently short period of time, $\delta t$,
:   $P[K(\delta t) = 1] \approx \lambda \delta t$\
    $P[K(\delta t) \ge 2]$ is negligible (i.e., two events don't occur at the same time)

The number of events over time 0 to $t$
(or from arbitrary starting time $s$ to $s + t$) is Poisson distributed,
for $k = 0, 1, 2, ...$,

<div>\[\begin{aligned}
P[K(t) = k]          \;&=\; \frac {(\lambda t)^k} {k!} e^{-\lambda t} \\
P[K(t+s) - K(s) = k] \;&=\; P[K(t) = k]
\end{aligned}\]</div>

(Eq. 2.53)
Waiting time to the first event is exponentially distributed,

<div>\[
f_T(t) = \lambda e^{-\lambda t}
\]</div>

*memoryless* (**無記憶性**)
:   The number of coin-tosses required to observe the next heads is independent of previous results.\
    The waiting times between successive events are *i.i.d.* (independent and identically distributed).

------------------------------------------------------------------------

(Eq. 2.54)
The waiting time $W$ until the $n$ th event
(= the sum of $n$ *i.i.d.* waiting times)
can be derived by $n - 1$ successive convolutions:

<div>\[\begin{aligned}
f_{W,1}(t) &= f_T(t) = \lambda e^{-\lambda t} \\
f_{W,2}(t) &= \int _0^t \lambda e^{-\lambda (t-t')} f _{W,1}(t') dt' \\
           &= \lambda ^2 e^{-\lambda t} \left[t'\right]_0^t \\
           &= t\lambda^2 e^{-\lambda t} \\
f_{W,3}(t) &= \int _0^t \lambda e^{-\lambda (t-t')} f _{W,2}(t') dt' \\
           &= \frac {t^2 \lambda^3 e^{-\lambda t}} {2!} \\
           &\;\vdots \\
f_{W,n}(t) &= \lambda e^{-\lambda t} \frac{(\lambda t)^{n-1}} {(n - 1)!}
\end{aligned}\]</div>

This is the **gamma distribution**.

(Eq. 2.55)
The mean and the variance are

<div>\[\begin{aligned}
\operatorname{E}[W] &= \sum^n \operatorname{E}[T] \\
              &= \sum^n \frac 1 \lambda \\
              &= \frac n \lambda \\
\operatorname{Var}[W] &= \sum^n \operatorname{Var}[T] \\
              &= \sum^n \frac 1 {\lambda^2} \\
              &= \frac n {\lambda^2}
\end{aligned}\]</div>

These hold even when $n$ is not an integer
if we replace the factorial $(n - 1)!$ with *gamma function*,
$\Gamma(n) = \int^\infty_0 x^{n-1} e^{-x} dx$.
(Eq. 2.57)

------------------------------------------------------------------------

The coalescent considers the events
that have a very small probability of occurring in any single generation:

-   coalescence (common ancestor event)
-   mutation
-   migration

They will each form a Poisson process.

#### 2.2.1 Poisson Process Results for the Coalescent

##### The Sum of Independent Poissons

Two independent Poisson random variables:
:   $X_1$ with occurrence rate $\lambda _1$\
    $X_2$ with occurrence rate $\lambda _2$

(Eq. 2.58)
The distribution of $Y = X_1 + X_2$ can be obtained by convolution:

<div>\[\begin{aligned}
P[Y=k] \;&=\; \sum _{i=0}^k P[X_1=i] \; P[X_2=k-i] \\
         &=\; \sum _{i=0}^k \frac {\lambda _1^i} {i!} e^{-\lambda _1}
                            \frac {\lambda _2^{k-i}} {(k-i)!} e^{-\lambda _2} \\
         &=\; e^{-\lambda _1} e^{-\lambda _2}
              \sum _{i=0}^k \frac {\lambda _1^i} {i!}
                            \frac {\lambda _2^{k-i}} {(k-i)!} \frac {k!}{k!} \\
         &=\; \frac {e^{-(\lambda _1 + \lambda _2)}} {k!}
              \sum _{i=0}^k {k \choose i} \lambda _1^i \lambda _2^{k-i} \\
         &=\; \frac {e^{-(\lambda _1 + \lambda _2)}} {k!} (\lambda _1 + \lambda _2)^k \\
         &=\; \text{Poisson distribution with occurrence rate } \lambda _1 + \lambda _2
\end{aligned}\]</div>

The sum of independent Poisson processes is another Poisson process.

##### The Probability that the First Event Is of a Particular Type

(Eq 2.60)
The probability that $X_1$ is observed before $X_2$
is given simply by the relative rate of the event
(i.e., as a fraction of the total rate):

<div>\[\begin{aligned}
P[T_1 < T_2]
   &= \int _0^\infty P[T_2>t]\; f _{T_1}(t) dt \\
   &= \int _0^\infty \left(e^{-\lambda _2 t} \right) _\text{Eq. 2.59}\;
                     \lambda _1 e^{-\lambda _1 t} dt \\
   &= \lambda _1 \int _0^\infty e^{-(\lambda _1 + \lambda _2) t} dt \\
   &= \lambda _1 \left[-\frac {e^{-(\lambda _1 + \lambda _2)t}}
                             {\lambda _1 + \lambda _2} \right]_0^\infty \\
   &= \frac {\lambda _1} {\lambda _1 + \lambda _2},
\end{aligned}\]</div>

using (Eq 2.59)

<div>\[\begin{aligned}
P[T>t] \;&=\; \int _t^\infty \lambda e^{-\lambda t} dt \\
         &=\; \left[-e^{-\lambda t} \right]_t^\infty \\
         &=\; e^{-\lambda t}.
\end{aligned}\]</div>

##### The Time to the First Event among Independent Poissons

(Eq. 2.61)
The distribution of $T = \min(T_1, T_2)$

<div>\[\begin{aligned}
P[T>t] \;&=\; P[\min(T_1, T_2) > t] \\
         &=\; P[T_1 > t \;\cap\; T_2 > t] \\
         &=\; P[T_1 > t]\; P[T_2 > t] \\
         &=\; e^{-\lambda _1 t} e^{-\lambda _2 t} \\
         &=\; e^{-(\lambda _1 + \lambda _2) t}
\end{aligned}\]</div>

Therefore, $f_ {\min(T_1, T_2)}(t) = (\lambda _1 + \lambda _2) e^{-(\lambda _2 + \lambda _1)t}$

There is a one-to-one correspondence between cumulative distribution $P[T \le t]$ and probability densities $f_T(t)$

##### The Number of Events Required to See a Particular Outcome

$X_2, X_2, X_2, ..., X_2, \boldsymbol{X_1}$, ...

(Eq. 2.62)
How many $X_2$ (e.g., mutation events) occur
before $X_1$ (e.g., common ancestor event)?

<div>\[\begin{aligned}
P[K=k] \;&=\; P[\text{First }X_1\text{ occurs at }K\text{th trial}] \\
         &=\; P[X_2\text{ occurs }k - 1\text{ times first, then }X_1\text{ occurs}] \\
         &=\; \left(\frac {\lambda _2} {\lambda _1 + \lambda _2} \right)^{k-1}
                    \frac {\lambda _1} {\lambda _1 + \lambda _2}
\end{aligned}\]</div>

This is the *geometric distribution* with the rate
$p = \frac {\lambda _1} {\lambda _1 + \lambda _2}$.

The process is just like a series of *Bernoulli trials* with probability of success
$p = \frac {\lambda _1} {\lambda _1 + \lambda _2}$.

##### Tying All This Together: A Filtered Poisson Process

Reinterpret $f_T(t)$ with the sum rule and product rule
(see Eq. 2.7, 2.8).

<div>\[\begin{aligned}
f_T(t)\; &=\; \sum _{k=1}^\infty f_T(t \mid K=k)\; P[K=k] \\
         &=\; \sum _{k=1}^\infty (\text{Eq. }2.54) (\text{Eq. }2.62) \\
         &=\; \sum _{k=1}^\infty
              (\lambda _1 + \lambda _2) e^{-(\lambda _1 + \lambda _2)t}
              \frac {\{(\lambda _1 + \lambda _2)t\}^{k-1}} {(k-1)!}\;
              \left(\frac {\lambda _2} {\lambda _1 + \lambda _2} \right)^{k-1}
                    \frac {\lambda _1} {\lambda _1 + \lambda _2} \\
         &=\; \lambda _1 e^{-(\lambda _1 + \lambda _2)t}\;
              \sum _{k=1}^\infty \frac {(\lambda _2 t)^{k-1}} {(k-1)!} \\
         &=\; \lambda _1 e^{-(\lambda _1 + \lambda _2)t}\;
              \sum _{k=0}^\infty \frac {(\lambda _2 t)^k} {k!} \\
         &=\; \lambda _1 e^{-(\lambda _1 + \lambda _2)t}\; e^{\lambda _2 t} \\
         &=\; \lambda _1 e^{-\lambda _1 t}
\end{aligned}\]</div>

> [!note]
> Taylor series of $e^x$
>
> <div>\[\begin{aligned}
> e^x \;&=\; 1 + \frac x 1 + \frac {x^2} {2!} + \frac {x^3} {3!} + \cdots \\
>       &=\; \sum _{k=0}^\infty \frac {x^k} {k!}
> \end{aligned}\]</div>

Filtered Poisson process = Poisson process with rate $\lambda p\; (=\lambda _1)$
:   -   total occurrence rate $\lambda = \lambda _1 + \lambda _2$
    -   acceptance rate (proportion of the focal event)
        $p = \frac {\lambda _1} {\lambda _1 + \lambda _2}$

#### 2.2.2 Convolutions of Exponential Distributions

The sum of the waiting times of $n$ events:

If the rate does not change with each event ($\lambda _i = \lambda \text{ for all } i$),
:   → gamma-distributed (Eq. 2.54)

If the rate changes with each event ($\lambda _i \neq \lambda _j \text{ for } i \neq j$),
:   (e.g., in the study of genealogies, the rate of coalescence changes every time a coalescent event occurs)\
    → obtained by successive convolution of exponential distribution (Eq. 2.63, 2.64)

<div>\[\begin{aligned}
f_{T_1 + T_2}(t)
   &= \int _0^t f_{T_1}(s)\; f_{T_2}(t-s)ds \\
   &= \int _0^t \lambda _1 e^{-\lambda _1 s}\; \lambda _2 e^{-\lambda _2 (t-s)}ds \\
   &= \lambda _1 \lambda _2 e^{-\lambda _2 t} \int _0^t e^{-(\lambda _1 -\lambda _2)s}ds \\
   &= \lambda _1 \lambda _2 e^{-\lambda _2 t}
      \left[-\frac 1{\lambda _1 - \lambda _2} e^{-(\lambda _1 -\lambda _2)s} \right]_0^t \\
   &= \frac {\lambda _1} {\lambda _1 - \lambda _2} \lambda _2 e^{-\lambda _2 t}
      \left(1 - e^{-(\lambda _1 -\lambda _2)t} \right) \\
   &= \frac {\lambda _1} {\lambda _1 - \lambda _2} \lambda _2 e^{-\lambda _2 t} +
      \frac {\lambda _2} {\lambda _2 - \lambda _1} \lambda _1 e^{-\lambda _1 t} \\
   &= \frac {\lambda _1} {\lambda _1 - \lambda _2} f _{T_2}(t) +
      \frac {\lambda _2} {\lambda _2 - \lambda _1} f _{T_1}(t) \\
  (&= \text{weighted sum of the original distributions})\\[1ex]
f_{T_1 + T_2 + T_3}(t)
   &= \int _0^t f_{T_1 + T_2}(s)\; f_{T_3}(t-s)ds \\
   &= \frac {\lambda _2} {\lambda _2 - \lambda _1}
      \frac {\lambda _3} {\lambda _3 - \lambda _1} \lambda _1 e^{-\lambda _1 t} +
      \frac {\lambda _3} {\lambda _3 - \lambda _2}
      \frac {\lambda _1} {\lambda _1 - \lambda _2} \lambda _2 e^{-\lambda _2 t} +
      \frac {\lambda _1} {\lambda _1 - \lambda _3}
      \frac {\lambda _2} {\lambda _2 - \lambda _3} \lambda _3 e^{-\lambda _3 t} \\[1ex]
f_{T_1 + T_2 + T_3 + T_4}(t) &= \cdots \\
   &\;\vdots \\
f_{\sum _{i=1}^n T_i}(t)
   &= \sum_{i=1}^n \lambda _i e^{-\lambda _i t}
      \prod _{j=1,\;j \neq i} \frac {\lambda _j} {\lambda _j - \lambda _i}
\end{aligned}\]</div>

We will use this to obtain the distribution of the total waiting time
to the MRCA of the entire samples (Eq. 3.27)
