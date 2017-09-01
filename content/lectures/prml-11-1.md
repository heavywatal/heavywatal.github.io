+++
title = 'PRML輪読会 11章1節'
tags = ["math", "book"]
[menu.main]
  parent = "lectures"
+++

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/0387310738/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/612j5Uo43eL._SX180_.jpg" alt="Pattern Recognition and Machine Learning (Information Science and Statistics)" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4621061224/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41O0QFyTHJL._SX160_.jpg" alt="パターン認識と機械学習 上" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4621061240/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/418MuoJetFL._SX160_.jpg" alt="パターン認識と機械学習 下 (ベイズ理論による統計的予測)" /></a>

Author
:   Christopher M. Bishop

Book
:   [Pattern Recognition and Machine Learning](http://www.amazon.co.jp/exec/obidos/ASIN/4621061224/heavywatal-22/)

Publisher
:   [Springer](http://www.springer.com/computer/image+processing/book/978-0-387-31073-2)

Materials
:   <http://research.microsoft.com/en-us/um/people/cmbishop/prml/>

輪読担当
:   岩嵜航

日程
:   2014-12-08

## 11. Sampling Methods

10章では決定論的な近似を見てきた。
この章ではサンプリングを伴う *Monte Carlo* 法を取り扱う。

{{%div class="note"%}}
モンテカルロ法の由来

スタニスワフ・ウラムがソリテアの成功率を考えてた時に思いついて、
同僚のジョン・フォン・ノイマンが計算機上での実用まで持ってったらしい。
モナコ公国のモンテカルロ地区に国営カジノがあって、
ウラムの叔父がそこで負けて親戚から借金したことにちなんで
同僚のニコラス・メトロポリスが命名したらしい。
{{%/div%}}

**目標:**
変数 $\mathbf z$ の分布 $\color{red}{p(\mathbf z)}$ を考えた上で、
ある関数 $\color{blue}{f(\mathbf z)}$ の値がどうなるか予測したい。

{{%div class="note"%}}
[Figure 11.1](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure11.1.png)

<p><img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure11.1.png"/ alt="Figure 11.1" width="300px"></p>
{{%/div%}}

$\color{blue}{f(\textbf{z})}$ の期待値は (**式11.1**)

<div>\[\begin{split}
\mathbb E[f] = \int f(\mathbf z) p(\mathbf z) \mathrm d \mathbf z
\end{split}\]</div>

みたいな感じで表せるが、だいたいは複雑過ぎて解析的に解けないので、そういうときどうしようかという話。

$\color{red}{p(\mathbf z)}$ の分布から $L$ 個サンプリングしてきた
$\mathbf{z}_l$ をそれぞれ $f$ に放り込んで平均を取ってみよう (**式11.2**)。

<div>\[\begin{split}
\hat f = \frac 1 L \sum^L_{l=1} f(\mathbf z_l)
\end{split}\]</div>

その推定値の期待値は (**Exercise 11.1**)

<div>\[\begin{split}
\mathbb E[\hat f] &= \mathbb E \left[\frac 1 L \sum^L_{l=1} f(\mathbf z_l) \right] \\
                  &= \frac 1 L \sum^L_{l=1} \int f(\mathbf z_l) p(\mathbf z_l) \mathrm d \mathbf z_l \\
                  &= \frac 1 L \sum^L_{l=1} \mathbb E[f] \\
                  &= \mathbb E[f]
\end{split}\]</div>

で真の期待値と同じになる。
この推定値の分散は (**Exercise 11.1**, **式 11.3**)

<div>\[\begin{split}
\mathrm{var}[\hat f] &= \mathrm{var} \left[\frac 1 L \sum^L_{l=1} f(\mathbf z_l) \right] \\
                     &= \frac 1 {L^2} \sum^L_{l=1} \mathrm{var}[f(\mathbf z)] \\
                     &= \frac 1 L \mathrm{var}[f] \\
                     &= \frac 1 L \mathbb E[(f - \mathbb E[f])^2]
\end{split}\]</div>

となる。注意すべき点としては:

-   推定精度が次元数によらない
-   基本的には $L$ をそんなに大きく取らなくても(10とか20くらいで)よさそう
-   ただし、サンプルが独立じゃない場合にはその辺を加味した有効サンプル数が十分になるように多めに取るべし
-   $\color{red}{p(\mathbf z)}$ が大きいところで $\color{blue}{f(\mathbf z)}$ がゼロに近くなるような場合、少確率で出てくる大きな値に推定値が引っ張られることがあるので比較的多めに取るべし

------------------------------------------------------------------------

$\color{red}{p(\mathbf z)}$ が実は $p(z_1, z_2, ..., z_M)$ という同時確率だということを思い出そう。
$z_i$ がそれぞれ独立な分布から出てくる場合はいいとして、そうじゃない場合はどうしたらいいか？

{{%div class="note"%}}
[Figure 8.2](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure8.2.png)

<p><img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure8.2.png"/ alt="Figure 8.2" width="200px"></p>
変数の因果関係がこのような閉路なし有向グラフで表せる場合、同時確率は **式 8.4**
$p(x_1)p(x_2)p(x_3)p(x_4 \mid x_1,x_2,x_3)p(x_5 \mid x_1,x_3)p(x_6 \mid x_4)p(x_7 \mid x_4,x_5)$
のように条件付き確率の積で表せる。
{{%/div%}}

依存関係の親となるほうから順に条件付き確率で生成 (*ancestral sampling* **伝承サンプリング**)
していくことにすると、同時確率は一般的に (**式 11.4**)

<div>\[\begin{split}
p(\mathbf z) = \prod_{i=1}^M p(\mathbf z_i \mid \mathrm{pa}_i)
\end{split}\]</div>

というふうに書ける。
変数の一部が観測可能な場合は *logic sampling*
(セクション11.1.4で登場する **重点サンプリング** *importance sampling* の特殊ケース)
が使える。

因果が分からなくて無向グラフで表されるような場合には
$z_1$ から $z_M$ まで一周するだけでは求まらず、
ギブズサンプリング (Gibbs sampling) のような計算量のかかる手法が必要になる。

### 11.1. Basic Sampling Algorithms

コンピュータ上でサンプリングを行うときに真の乱数を使うことは稀で、だいたいは適当なシードから決定論的な過程で擬似乱数を生成することになる。
擬似乱数の質も問題になったりするけどこの本では詳しく扱わない。
いい感じで $(0,1)$ の一様乱数が生成できるものとして進める。

{{%div class="note"%}}
Unix/Linux系OSが提供する乱数

ハードウェア的なノイズから生成した真の乱数は `/dev/random` から読み出せるが、
いくつも生成しようとするとノイズが溜まるまで待たされることになるのであまり使わない。
待ち時間無しにそれなりの擬似乱数を作ってくれるデバイスとして `/dev/urandom`
があるが、ここから毎回読み出すのもコストが高いので、シード生成にのみ使う。
{{%/div%}}

{{%div class="note"%}}
[/cxx/random]({{< relref "cxx/random.md" >}})
{{%/div%}}

#### 11.1.1 Standard distributions

変数 $z$ が $(0,1)$ の一様乱数だとして、
適当な関数をかけて $y = f(z)$ とするとその分布は (**式 11.5**)

<div>\[\begin{split}
p(y) = p(z) \left| \frac {\mathrm dz} {\mathrm dy} \right|
\end{split}\]</div>

となる。
変換後の乱数 $y$ が任意の形の分布 $\color{red}{p(y)}$ に従うようにするにはどうしたらいいか。
$\color{red}{p(y)}$ の不定積分を (**式 11.6**)

<div>\[\begin{split}
z = h(y) \equiv \int _{-\infty}^y p(\hat y) \mathrm d\hat y
\end{split}\]</div>

のように $\color{blue}{h(y)}$ として定義してみると **図 11.2** のような関係になる。

{{%div class="note"%}}
[Figure 11.2](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure11.2.png)

<p><img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure11.2.png"/ alt="Figure 11.2" width="300px"></p>
縦軸を $z$ として青い線を逆関数の目線で見てみると、$y$ が中央付近に来るような $z$ の区間はすごく短いが、$y$ が左から3分の1くらいのところに来るような $z$ の区間はかなり長い。
{{%/div%}}

その不定積分の逆関数を一様乱数にかけて $y = h^{-1}(z)$ とすれば欲しかった分布の乱数が出てくる！

例えば指数分布だと (**式 11.7**)

<div>\[\begin{split}
p(y) &= \lambda \exp(-\lambda y) \\
z = h(y) &= \int_0^y \lambda \exp(-\lambda \hat y) \mathrm d \hat y \\
         &= \left[-\exp(-\lambda \hat y) \right]_0^y \\
         &= 1 - \exp(-\lambda y) \\
\exp(-\lambda y) &= 1 - z \\
     -\lambda y  &= \ln(1 - z) \\
               y &= -\frac {\ln(1 - z)} \lambda
\end{split}\]</div>

となるので、$y = -\lambda^{-1} \ln(1 - z)$ とすれば $y$ は指数分布に従う乱数となる。

------------------------------------------------------------------------

別の例としてコーシー分布も同じように変換できる (**式 11.8**, **Exercise 11.3**)

<div>\[\begin{split}
p(y) &= \frac 1 \pi \frac 1 {1 + y^2} \\
z = h(y) &= \int_{-\infty}^y \frac 1 \pi \frac 1 {1 + \hat y^2} \mathrm d \hat y \\
         &= \frac 1 \pi \left[\arctan(\hat y) \right]_{-\infty}^y \\
         &= \frac 1 \pi \left(\arctan(y) + \frac \pi 2 \right) \\
         &= \frac {\arctan(y)} \pi + \frac 1 2 \\
\arctan(y) &= \pi(z - \frac 1 2) \\
         y &= \tan\left[\pi(z - \frac 1 2)\right]
\end{split}\]</div>

------------------------------------------------------------------------

多変量の場合はヤコビアンを使えばよい

<div>\[\begin{split}
p(y_1, ..., y_M) = p(z_1, ..., z_M) \left| \frac {\partial (z_1, ..., z_M)}
                                                {\partial (y_1, ..., y_M)} \right|
\end{split}\]</div>

例として2系統の独立な正規乱数を生成する *Box-Muller* 法を見てみる。
まず $(-1,1)$ の一様乱数をふたつ $z_1, z_2$ として取ってきて、
$z_1^2 + z_2^2 \leq 1$ を満たさなければ捨てる。
これは下図の円の中に収まる一様乱数だけ取ってくることに相当する。

{{%div class="note"%}}
[Figure 11.3](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure11.2.png)

<p><img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure11.3.png"/ alt="Figure 11.3" width="200px"></p>
{{%/div%}}

$r^2 = z_1^2 + z_2^2$ として (**式 11.10**, **式 11.11**)

<div>\[\begin{split}
y_1 &= z_1 \left(\frac {-2\ln r^2} {r^2}\right)^{1/2}\\
y_2 &= z_2 \left(\frac {-2\ln r^2} {r^2}\right)^{1/2}
\end{split}\]</div>

のように変換すると $y_1$ と $y_2$ の同時分布は

<div>\[\begin{split}
p(y_1, y_2) &= p(z_1, z_2) \left| \frac{\partial(z_1, z_2)} {\partial(y_1, y_2)} \right|\\
            &= \left[\frac 1 {\sqrt{2\pi}} \exp(-\frac {y_1^2} 2) \right]
               \left[\frac 1 {\sqrt{2\pi}} \exp(-\frac {y_2^2} 2) \right]
\end{split}\]</div>

のように表され、それぞれ独立な標準正規乱数になっていることがわかる。
平均と分散を変えたければ、$y = \mu + \sigma z$ のように標準偏差をかけて平均を足せばよい。

{{%div class="note"%}}
C++11 の `std::normal_distribution` や GSL の `gsl_ran_gaussian` でも使われている。
円に収まらないものを棄却する方法ではなく、三角関数を使ってそのまま用いる方法が取られる。
{{%/div%}}

多変量の場合も同様に $\mathbf y = \mathbf \mu + \mathbf{Lz}$ として動かせる。
ただし共分散は $\mathbf \Sigma = \mathbf{LL}^\mathrm T$ として **コレスキー分解** (*Cholesky decomposition*)する。
これは対称行列に特化したLU分解で、$\mathbf L$ は下三角行列になる。
変換後の平均と分散を確かめてみる (**Excersize 11.5**)

<div>\[\begin{split}
\mathbb E[\mathbf y] &= \mathbb E[\mathbf \mu + \mathbf{Lz}] = \mathbf \mu + \mathbf 0 \\
\mathrm{cov}[\mathbf y]
   &= \mathbb E\left[(\mathbf y - \mathbb E[\mathbf y])(\mathbf y - \mathbb E[\mathbf y])^\mathrm T \right] \\
   &= \mathbb E[(\mathbf \mu + \mathbf{Lz} - \mathbf \mu)(\mathbf \mu + \mathbf{Lz} - \mathbf \mu)^\mathrm T] \\
   &= \mathbb E[\mathbf{Lz}(\mathbf{Lz})^\mathrm T] \\
   &= \mathbf{LL}^\mathrm T = \mathbf \Sigma\\
\end{split}\]</div>

ただし一様乱数 $\mathbf z$ については $\mathbb E[\mathbf z] = \mathbf 0$ かつ $\mathbb E[\mathbf{zz}^\mathrm T] = \mathbf I$ 。

------------------------------------------------------------------------

ここで説明したような手法が使えるのは、不定積分の逆関数が簡単に得られるような場合だけ。
より一般的に使える *rejection sampling* と *importance sampling* について、この先で見ていく。
