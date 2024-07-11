+++
title = 'PRML輪読会 2章前半'
tags = ["math", "book"]
[menu.main]
  parent = "lectures"
+++

Author
:   Christopher M. Bishop

Book
: [Pattern Recognition and Machine Learning](https://www.amazon.co.jp/dp/0387310738?&linkCode=ll1&tag=heavywatal-22&linkId=0dacd1cec1bcc3d73dc0a9f27d158183)
: [パターン認識と機械学習 上](https://www.amazon.co.jp/dp/4621061224?&linkCode=ll1&tag=heavywatal-22&linkId=07bf2a676e8cf9d3f62a5ae847fa4962)
: [パターン認識と機械学習 下](https://www.amazon.co.jp/dp/4621061240?&linkCode=ll1&tag=heavywatal-22&linkId=c376a10b53faed6da45c3591d5dbc61a)

Publisher
:   [Springer](http://www.springer.com/computer/image+processing/book/978-0-387-31073-2)

Materials
:   <http://research.microsoft.com/en-us/um/people/cmbishop/prml/>

輪読担当
:   岩嵜航

日程
:   2014-03-10

## 2. Probability Distributions

*density estimation*: **密度推定**
:   与えられた確率変数xのセットから確率分布p(x)をモデル化すること。

*parametric* method: **パラメトリック法**
:   平均や分散といった少数の変数によって規定される分布を仮定して当てはめる。
    頻度主義的には、尤度関数など特定のクライテリアに従って最適化して値を決める。
    ベイズ主義的には、事前分布と観察データから事後分布を得る。

**i.i.d. 独立同分布** あるいは **独立同一分布**
:   **independent and identically distributed** の略。
    確率変数のセットが同じ分布に従っていて独立なとき。

2.1–2.2 離散型確率変数

2.3 *Gaussian distribution*: **正規分布**

2.4 *exponential family*: **指数型分布族**

2.5 *nonparametric* method: **ノンパラメトリック法**

### 2.1 Binary Variables

**二値変数**
$x \in \{0,1\}$

コインの表裏(head/tail)のように、2つの値を取りうる確率変数。
e.g. 表なら *x* = 1、裏なら *x* = 0。

------------------------------------------------------------------------

*Bernoulli* distribution **ベルヌーイ分布**

表の出る確率が *μ* であるコインを1回投げて、
表が出る確率と裏が出る確率はそれぞれ

<div>\[\begin{aligned}
p(x = 1 \mid \mu) &= \mu\\
p(x = 0 \mid \mu) &= 1 - \mu
\end{aligned}\]</div>

*x* の確率分布として書いて平均と分散 (Exercise 2.1) を求めると

<div>\[\begin{aligned}
\text{Bern}(x \mid \mu) &= \mu ^ x (1 - \mu) ^ {1 - x} \\
\operatorname{E}[x] &= 0(1 - \mu) + 1\mu = \mu\\
\operatorname{var}[x] &= \operatorname{E}[(x - \mu)^2] = \mu^2 (1 - \mu) + (1 - \mu)^2 \mu = \mu (1 - \mu)
\end{aligned}\]</div>

パラメータ *μ* の下で *N* 回投げたデータセット $D = \{x_1, ..., x_N\}$
が得られる確率、すなわち尤度は

<div>\[
p(D \mid \mu) = \prod_{n = 1}^N {p(x_n \mid \mu)}
           = \prod_{n = 1}^N {\mu^{x_n} (1 - \mu)^{1 - x_n}}
\]</div>

これを最大化する *μ* を求めるため、まず対数を取って

<div>\[\begin{aligned}
\ln p(D \mid \mu) &= \sum_{n = 1}^N {\ln p(x_n \mid \mu)}\\
               &= \sum_{n = 1}^N {\{x_n \ln \mu + (1 - x_n) \ln (1 - \mu)\}}\\
               &= \ln \mu \sum_{n = 1}^N x_n + \ln (1 - \mu) \sum_{n = 1}^N (1 - x_n)\\
               &= \{\ln \mu - \ln(1 - \mu)\} \sum_{n = 1}^N x_n + N \ln (1 - \mu)\\
\end{aligned}\]</div>

<div class="note">

*sufficient statistic*: **十分統計量**

ここで対数尤度は個々の $x_n$ によらず、
総和 $\sum_n {x_n}$ だけに依存している。
そんな感じのやつを十分統計量と呼ぶが、ここでは詳しく触れない。
</div>

*μ* で微分したものが0になるように

<div>\[\begin{aligned}
\frac{d}{d\mu} \ln p(D \mid \mu)
   = \{\frac{1}{\mu} + \frac{1}{1 - \mu}\} \sum_{n = 1}^N x_n - \frac{N}{1 - \mu}
   &= 0\\
- \frac{N}{1 - \mu} &= - \frac{1}{\mu (1 - \mu)} \sum_{n = 1}^N x_n\\
\mu_\text{ML} &= \frac{1}{N} \sum_{n = 1}^N {x_n}\\
              &= \frac{m}{N}
\end{aligned}\]</div>

結局のところ標本平均と同じ。
コイン3回投げて3回表だと最尤は *μ* = 1 だが、これはどう考えてもover-fittingである。
だいたい半々で表裏が出るっしょ、という事前情報で補正していきたい。

------------------------------------------------------------------------

*Binomial* distribution **二項分布**

確率 *μ* で表の出るコインを *N* 回投げて表が出る回数 *m* の確率分布

<div>\[
\text{Bin}(m \mid N, \mu) = \binom{N}{m} \mu^m (1 - \mu)^{N - m}
\]</div>

<div class="note">

[Figure 2.1](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.1.png)

<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.1.png" alt="Figure 2.1" width="300px">

例えば *N* = 10, *μ* = 0.25 のときの *m* の頻度分布
</div>

1回1回の観察は独立なベルヌーイ試行であり、
そういうときは $\operatorname{E}[x + z] = \operatorname{E}[x] + \operatorname{E}[z]$ かつ
$\operatorname{var}[x + z] = \operatorname{var}[x] + \operatorname{var}[z]$ が成り立つので
(Exercise 1.10)、平均と分散は (Exercise 2.4)

<div>\[\begin{aligned}
\operatorname{E}[m] &= \operatorname{E}[\sum_{n=1}^N x_n] = \sum_{n=1}^N \operatorname{E}[x_n] = \sum_{n=1}^N \mu = N \mu \\
\operatorname{var}[m] &= \operatorname{var}[\sum_{n=1}^N x_n] = \sum_{n=1}^N \operatorname{var}[x_n]
              = \sum_{n=1}^N \mu (1 - \mu) = N \mu (1 - \mu)
\end{aligned}\]</div>

#### 2.1.1 beta distribution

3回続けて表が出たからといって μ=1 なわけがない、
という事前情報をうまく取り入れて過学習を避けたい。
事前分布 *p(μ)* を導入しよう。

<div>\[
\text{posterior} \propto \text{prior} \times \text{likelihood}
\]</div>

尤度が $\mu^x (1 - \mu)^{1 - x}$ という形なので、
事前分布も *μ* と *1 - μ* の累乗にしておくと、
事後分布と事前分布のの関数形が同じになって (**共役性**, *conjugacy*)
いろいろ便利 (後述)。

<div>\[
\text{Beta}(\mu \mid a, b) =
   \frac{\Gamma(a + b)}{\Gamma(a) \Gamma(b)} \mu^{a-1} (1 - \mu)^{b-1}
\]</div>

[ガンマ関数](http://ja.wikipedia.org/wiki/%E3%82%AC%E3%83%B3%E3%83%9E%E9%96%A2%E6%95%B0)
のところは(積分して1になるように)正規化するための係数で、
形を決めるのは後ろの部分。
*a* と *b* はパラメータ *μ* の分布を決めるための **超パラメータ** (*hyperparameter*)。

平均と分散は (Exercise 2.6)

<div>\[\begin{aligned}
\operatorname{E}[\mu] &= \int_0^1 \mu \text{Beta}(\mu \mid a,b)d\mu\\
              &= \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \int_0^1 \mu^{a} (1-\mu)^{b-1} d\mu\\
              &= \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \frac{\Gamma(a+1)\Gamma(b)}{\Gamma(a+b+1)}\\
              &= \frac{\Gamma(a+b)}{\Gamma(a)} \frac{a\Gamma(a)}{(a+b)\Gamma(a+b)}\\
              &= \frac{a}{a + b}\\
\operatorname{var}[\mu] &= \int_0^1 \mu^2 \text{Beta}(\mu \mid a,b)d\mu - \operatorname{E}[\mu]^2\\
                &= \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \frac{\Gamma(a+2)\Gamma(b)}{\Gamma(a+b+2)} - \operatorname{E}[\mu]^2\\
                &= \frac{\Gamma(a+b)}{\Gamma(a)} \frac{a(a+1)\Gamma(a)}{(a+b)(a+b+1)\Gamma(a+b)} - \operatorname{E}[\mu]^2\\
                &= \frac{a(a+1)}{(a+b)(a+b+1)} - (\frac{a}{a+b})^2\\
                &= \frac{ab}{(a + b)^2 (a + b + 1)}
\end{aligned}\]</div>

<div class="note">

[Figure 2.2a](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.2a.png), [Figure 2.2b](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.2b.png), [Figure 2.2c](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.2c.png), [Figure 2.2d](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.2d.png)

<p>
<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.2a.png" alt="Figure 2.2a" width="240px">
<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.2b.png" alt="Figure 2.2b" width="240px">
<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.2c.png" alt="Figure 2.2c" width="240px">
<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.2d.png" alt="Figure 2.2d" width="240px">
</p>
ベータ分布がどんな形になるか、いろんな *a*, *b* でプロットしてみた
</div>

これを二項分布の尤度関数(2.9)と掛け算して得られる事後分布の比例関係は

<div>\[
p(\mu \mid m, N, a, b) \propto \mu^{m + a - 1} (1 - \mu)^{N - m + b - 1}
\]</div>

ベータ分布と同じ形をしているので、同じように正規化できる:

<div>\[\begin{aligned}
p(\mu \mid m, N, a, b) &= \frac{\Gamma(N + a + b)}{\Gamma(m + a) \Gamma(N - m + b)}
                          \mu^{m + a - 1} (1 - \mu)^{N - m + b - 1}\\
                       &= \text{Beta}(\mu \mid m + a, N - m + b)
\end{aligned}\]</div>

この事後分布は、新しいデータを加えて再評価するときに事前分布として使える。

<div class="note">
Figure 2.3

<p>
<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.3a.png" alt="Figure 2.3a" width="240px">
<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.3b.png" alt="Figure 2.3b" width="240px">
<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.3c.png" alt="Figure 2.3c" width="240px">
</p>

*a* = 2, *b* = 2 の事前分布に([Figure 2.3a](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.3a.png))、
1回投げて表が出たというデータが追加されると([Figure 2.3b](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.3b.png))、
事後分布は *a* = 3, *b* = 2 のベータ分布になる([Figure 2.3c](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.3c.png))。
</div>

表の観測数 *m* 回、裏の観測数 *N - m* 回が追加されたら、
事前分布の *a* と *b* をその分それぞれ増やすだけで事後分布が得られる。
したがって超パラメータ *a* と *b* はその時点での表と裏の **有効観測数**
(*effective number of observations*) と解釈することができる。

*sequential* approach **逐次学習**
:   観測毎に事後分布を更新していけるのでリアルタイムで利用できる。
    処理し終わったデータは捨ててもよい。

------------------------------------------------------------------------

コインの真の特性 *μ* を知るのが目的ではなく、
手元のデータを元に次のコイントスの結果を予測したいとしたら、
興味があるのは $p(\mu)$ ではなくて $p(x = 1 \mid D)$

<div>\[\begin{aligned}
p(x = 1 \mid D) = \int_0^1 p(x = 1 \mid \mu) p(\mu \mid D) d\mu
                = \int_0^1 \mu p(\mu \mid D) d\mu
                = \operatorname{E}[\mu \mid D]
                = \frac{m + a}{N + a + b}
\end{aligned}\]</div>

ベータ分布の平均値の式の *a* と *b* に観察数を足しただけ。
事後分布の平均値は常に最尤推定値と事前平均の間にあり、観察数を上げていくと最尤推定の $\frac{m}{N}$ に近づく。
分散は小さく(事後分布はシャープに)なっていく。

------------------------------------------------------------------------

データを生成する分布の上で事後平均を平均すると事前平均になる (Exercise 2.8):

<div>\[\begin{aligned}
\operatorname{E}_D[\operatorname{E}_\theta[\theta \mid D]]
   &\equiv \int\left\{\int \theta p(\theta \mid D) d\theta\right\}p(D)dD\\
   &= \int\int \theta p(\theta,D) d\theta dD\\
   &= \int\left\{\int p(\theta,D) dD \right\} \theta d\theta\\
   &= \int p(\theta)\theta d\theta \equiv \operatorname{E}_\theta[\theta]
\end{aligned}\]</div>

事後分散の平均と事後平均の分散を足すと事前分散になる:

<div>\[\begin{aligned}
\operatorname{E}_D[\operatorname{var}_\theta[\theta \mid D]] + \operatorname{var}_D[\operatorname{E}_\theta[\theta \mid D]]
   &= \operatorname{E}_D[\operatorname{E}_\theta[\theta^2 \mid D] - \operatorname{E}_\theta[\theta \mid D]^2]
      + \operatorname{E}_D[\operatorname{E}_\theta[\theta \mid D]^2] - \operatorname{E}_D[\operatorname{E}_\theta[\theta \mid D]]^2\\
   &= \operatorname{E}_D[\operatorname{E}_\theta[\theta^2 \mid D]] - \operatorname{E}_D[\operatorname{E}_\theta[\theta \mid D]]^2\\
   &= \operatorname{E}_\theta[\theta^2] - \operatorname{E}_\theta[\theta]^2\\
   &= \operatorname{var}_\theta[\theta]
\end{aligned}\]</div>

データによって例外はあるが、平均的には「事前分散 &gt; 事後分散」となる。

### 2.2 Multinomial Variables

**多値変数** $x \in \{1, 2, ..., K\}$
:   サイコロの目のように、3つ以上の値を取りうる確率変数。 e.g. $x \in \{1,2,3,4,5,6\}$

**1-of-K 符号化法**
:   長さ *K* のベクトルのうち $x_k$ だけが1で、そのほかが0。
    例えばサイコロで3が出たらその観察値の表記は

<div>\[
\vec{x} = (0, 0, 1, 0, 0, 0)
\]</div>

------------------------------------------------------------------------

確率 $\mu_k$ で $x_k = 1$ になるとすると、
サイコロを1回振るときの *x* の分布は以下のように表せる。

<div>\[\begin{aligned}
p(\vec{x} \mid \vec{\mu})
  &= \prod _{k = 1}^K {\mu _k^{x _k}}\\
\sum _{\vec{x}} p(\vec{x} \mid \vec{\mu})
  &= \sum _{k=1}^K {\mu _k} = 1\\
\operatorname{E}[\vec{x} \mid \vec{\mu}]
  &= \sum _{\vec{x}} p(\vec{x} \mid \vec{\mu}) \vec{x}
   = \vec{\mu}
\end{aligned}\]</div>

------------------------------------------------------------------------

サイコロを *N* 回振った観察データ *D* に対応する尤度関数は

<div>\[\begin{aligned}
p(D \mid \vec{\mu}) = \prod_{n=1}^N \prod_{k=1}^K \mu_k^{x_{nk}}
                    = \prod_{k=1}^K \mu_k^{\sum_n x_{nk}}
                    = \prod_{k=1}^K \mu_k^{m_k}
\end{aligned}\]</div>

$m_k$ は *N* 回のうち *k* が出た回数。
出る順番や他の出目にはよらず、総和だけでよい、つまりこれも **十分統計量** の例。

$\mu_k$ の和が1になるという拘束条件の下で対数尤度を最大化する
$\mu_k$ を求めるには下記のようにラグランジュ未定乗数法を用いる。

<div>\[
\mu_k^{ML} = - \frac{m_k}{\lambda} = \frac{m_k}{N}
\]</div>

結局、観察総数 *N* のうちその目が出た数の割合が最尤推定値。

<div class="note">

ラグランジュ未定係数法 (Appendix E)

極値を求めたい関数と拘束条件をそれぞれ *f*, *g* で表すと

<div>\[\begin{aligned}
df &= \frac{\partial f}{\partial \mu_1} + ... + \frac{\partial f}{\partial \mu_K} = 0\\
dg &= \frac{\partial g}{\partial \mu_1} + ... + \frac{\partial g}{\partial \mu_K} = 0
\end{aligned}\]</div>
<div>\[
\frac{\partial f}{\partial \mu_k} + \lambda \frac{\partial g}{\partial \mu_k} = 0
\]</div>

て感じで拘束条件のない連立方程式に置き換えられる。
今回の例では

<div>\[\begin{aligned}
f(\vec{\mu}) = \sum_{k=1}^K m_k \ln \mu_k;\;
g(\vec{\mu}) = \sum_{k=1}^K \mu_k - 1
\end{aligned}\]</div>
<div>\[\begin{aligned}
\frac{\partial f}{\partial \mu_k} + \lambda \frac{\partial g}{\partial \mu_k} =
\frac{m_k}{\mu_k} + \lambda &= 0\\
m_k + \lambda \mu_k &= 0\;\therefore \mu_k^{ML} = - \frac{m_k}{\lambda}\\
\sum_k(m_k + \lambda \mu_k) &= 0\\
N + \lambda &= 0\\
\lambda &= -N
\end{aligned}\]</div>
</div>

------------------------------------------------------------------------

**多項分布** (*Multinomial distribution*)

観察総数とパラメータを条件とした、それぞれの出目の同時分布

<div>\[\begin{aligned}
\text{Mult}(m_1, ..., m_K \mid \vec{\mu}, N)
= \binom{N}{m_1, ..., m_K} \prod_{k=1}^{K} \mu_k^{m_k}
= \frac{N!}{m_1! ... m_K!} \prod_{k=1}^{K} \mu_k^{m_k}
\end{aligned}\]</div>

#### 2.2.1 Dirichlet distribution

多項分布の共役事前分布を考えよう。

式の形からするとおそらく $\mu_k$ の累乗にすればいいはず。
ということで、積分して1になるように正規化してみる。

<div>\[\begin{aligned}
\text{Dir}(\vec{\mu} \mid \vec{\alpha})
= \frac{\Gamma(\sum_k{\alpha_k})}{\Gamma(\alpha_1)...\Gamma(\alpha_K)} \prod_{k=1}^K \mu_k^{\alpha_k - 1}
\end{aligned}\]</div>

<div class="note">

[Figure 2.4](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.4.png)

<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.4.png" alt="Figure 2.4" width="200px">

$0 \le \mu_k \le 1$ かつ $\sum_k \mu_k = 1$
という制約下での $K$ 変数のディリクレ分布は $K – 1$ 次元の
**単体** (*simplex*) になる。
</div>

<div class="note">

[Figure 2.5a](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.5a.png), [Figure 2.5b](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.5b.png), [Figure 2.5c](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.5c.png)

<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.5a.png" alt="Figure 2.5a" width="200px">
<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.5b.png" alt="Figure 2.5n" width="200px">
<img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure2.5c.png" alt="Figure 2.5c" width="200px">

いろんな $α$ でのディリクレ分布。
simplexの面が水平軸方向に、密度が垂直軸になっている。
</div>

事後分布はこれと尤度の掛け算に比例する (2.40)。
それを積分して1になるよう正規化する (2.41)。

<div>\[\begin{aligned}
\text{posterior} &\propto \text{prior} \times \text{likelihood}\\
p(\vec{\mu} \mid D, \vec{\alpha})
   &\propto \text{Dir}(\vec{\mu} \mid \vec{\alpha}) \text{Mult}(D \mid \vec{\mu})\\
   &\propto \prod_{k=1}^K {\mu_k^{\alpha_k + m_k - 1}}\\
p(\vec{\mu} \mid D, \vec{\alpha})
&= \frac{\Gamma(\sum_k{\alpha_k} + N)}{\Gamma(\alpha_1 + m_1) ... \Gamma(\alpha_K + m_K)}
   \prod_{k=1}^K {\mu_k^{\alpha_k + m_k - 1}}\\
&= \text{Dir}(\vec{\mu} \mid \vec{\alpha} + \vec{m})
\end{aligned}\]</div>

確かに事後分布もディリクレ分布の形をしている。
*K* = 2 にすると二項分布・ベータ分布の話と一致。
逆に言うと、ディリクレ分布はベータ分布を一般化した多変量ベータ分布と見なせる。
超パラメータ $\alpha_k$ はサイコロで *k* が出た有効観察数のように解釈できる。

<div class="note">

**Johann Peter Gustav Lejeune Dirichlet** (1805–1859)

名前は 'le jeune de Richelet (リシュレから来た若者)' に由来。
最初の論文でフェルマーの最終定理の部分的な証明をして一躍有名に。
作曲家メンデルスゾーンの妹と結婚した。
</div>
