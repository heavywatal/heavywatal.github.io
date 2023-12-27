+++
title = 'PRML輪読会 3章4節'
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
:   2014-06-30

## 3. Linear Models For Regression

3.1 Linear Basis Function Models: 八島さん、関口さん

3.2 The Bias-Variance Decomposition: チャッキーさん

3.3 Bayesian Linear Regression: 佐伯さん、永田さん

### 3.4 Bayesian Model Comparison

最尤推定における過学習の問題 → 点推定じゃなくて周辺化することで回避しよう

-   訓練データだけでモデルを比較できる (確認データ不要)
-   すべてのデータを訓練に使うことができる (cross-validation不要)
-   複雑性のパラメータも含めて同時に決められる e.g. *relevance vector machine* (Chapter 7)

モデルの不確実さを確率で表し、加法定理・乗法定理を駆使して評価しよう

<div>\[\begin{aligned}
p(X)   &= \sum^Y p(X,Y) \\
p(X,Y) &= p(Y \mid X) P(X)
\end{aligned}\]</div>

<div class="note">

変数

新しい入力: $\mathbf x$\
それに対する出力(予測したい): $t$\
モデルの中のパラメータ: $\mathbf w$\
観察(トレーニング)データ: $\mathcal D$\
*L* 個のモデル: $\mathcal M_1, ..., \mathcal M_L$
</div>

------------------------------------------------------------------------

モデルの事後分布 $p(\mathcal M _i \mid \mathcal D)$ は、

-   $p(\mathcal M _i)$: どのモデルがアリかなという好み(事前分布)と、
-   $p(\mathcal D \mid \mathcal M _i)$:
    そのモデルの下での観察データの出やすさ
    (*model evidence*; *marginal likelihood* **周辺尤度**)

の積に比例する (**式 3.66**)。

<div>\[
p(\mathcal M _i \mid \mathcal D)
   \propto p(\mathcal M _i) p(\mathcal D \mid \mathcal M _i)
\]</div>

これを評価したいんだけど、
モデルの事前分布なんてだいたい分からないので、重要なのは後者のevidence。

<div class="note">

*Bayes factor* **ベイズ因子**

モデル $\mathcal M _j$ に対する $\mathcal M _i$ のevidence比
$\frac {p(\mathcal D \mid \mathcal M _i)} {p(\mathcal D \mid \mathcal M _j)}$
</div>

------------------------------------------------------------------------

**Mixture distribution**

モデルの事後分布が分かれば予測分布 *predictive distribution*
(新しい $\mathbf x$ に対して $t$ がどんな値となるか)
も加法定理と乗法定理より導かれる (**式 3.67**)

<div>\[\begin{aligned}
p(t \mid \mathbf x, \mathcal D)
   &= \sum _{i=1} ^L p(t, \mathcal M _i \mid \mathbf x, \mathcal D) \\
   &= \sum _{i=1} ^L {p(t \mid \mathbf x, \mathcal M _i, \mathcal D) p(\mathcal M _i \mid \mathcal D)}
\end{aligned}\]</div>

これは、それぞれのモデルでの予測分布(入力に対してどういう出力になりそうか)を
事後分布(どのモデルっぽいか)で重み付けした平均した、混合分布。

例えば L=2 でモデルの片方の予測が $t = a$ らへんの鋭いピーク、
もう片方のモデルの予測が $t = b$ らへんの鋭いピークだった場合、
混合分布の予測はその中点 $t = (a + b) / 2$ にピークを持つのではなく、二山になってしまう。

------------------------------------------------------------------------

**Model selection**

パラメータセット $w$ を持つモデル $\mathcal M_i$
のevidenceをまた加法定理と乗法定理でばらしてみると (**式 3.68**)

<div>\[\begin{aligned}
p(\mathcal D \mid \mathcal M _i)
   &= \int p(\mathcal D, \mathbf w \mid \mathcal M _i) \mathrm d \mathbf w \\
   &= \int p(\mathcal D \mid \mathbf w, \mathcal M _i) p(\mathbf w \mid \mathcal M _i) \mathrm d \mathbf w
\end{aligned}\]</div>

パラメータセットの尤度をその確率分布で重み付けして積分したもの、
ってことで周辺尤度と呼ばれるのが納得できる。
また、そのモデルからデータセットが生成される確率
(ただしパラメータは事前分布からランダムに取ったもの) とも理解できる。
この $p(\mathbf w \mid \mathcal M_i)$
はモデルで想定してる何らかの事前分布ってことでいいのかな？

<div class="note">

積分の中身からすると、パラメータの事後分布を求める式の正規化項になる (**式 3.69**)

<div>\[
p(\mathbf w \mid \mathcal D, \mathcal M _i)
   = \frac {p(\mathcal D \mid \mathbf w, \mathcal M _i) p(\mathbf w \mid \mathcal M _i)}
           {p(\mathcal D \mid \mathcal M _i)}
\]</div>
</div>

あるひとつのパラメータ $w$ を持つモデルを考える。

<div class="note">

[Figure 3.12](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure3.12.png) 近似

<p><img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure3.12.png" alt="Figure 3.12" width="300px"></p>
パラメータ $w$ の事前分布(青)と、それよりシャープな事後分布(赤)。
MAP推定値らへんで長方形に分布してるものとして近似。
</div>

[Figure 3.12](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure3.12.png) のように近似すると式3.68の積分をただの掛け算で書き変えられる
(モデル依存の表記を省略, **式 3.70, 3.71**)。

<div>\[\begin{aligned}
p(\mathcal D)
   &= \int p(\mathcal D \mid w) p (w) \mathrm dw \\
   &\simeq \frac 1 {\Delta w _\text{prior}} \int p(\mathcal D \mid w) \mathrm dw \\
   &\simeq \frac 1 {\Delta w _\text{prior}}
           p(\mathcal D \mid w _\text{MAP}) \Delta w _\text{posterior} \\
   &= p(\mathcal D \mid w _\text{MAP}) \frac {\Delta w _\text{posterior}}
                                             {\Delta w _\text{prior}} \\
\ln p(\mathcal D)
   &\simeq \ln p(\mathcal D \mid w _\text{MAP})
         + \ln \left( \frac {\Delta w _\text{posterior}} {\Delta w _\text{prior}} \right)
\end{aligned}\]</div>

第一項は一番いいパラメータの当てはまりの良さ、
第二項はモデルの複雑性によるペナルティ
(事後分布の幅が狭くなるほど大きな負になる)。

$M$ 個のパラメータを持つモデルを考える。
事前分布と事後分布の幅の比が全てのパラメータで等しいとすると (**式 3.72**)

<div>\[\begin{aligned}
p(\mathcal D)
   &= p(\mathcal D \mid w _\text{MAP}) \left(\frac {\Delta w _\text{posterior}}
                                                   {\Delta w _\text{prior}} \right)^M \\
\ln p(\mathcal D)
   &\simeq \ln p(\mathcal D \mid w _\text{MAP})
         + M \ln \left( \frac {\Delta w _\text{posterior}} {\Delta w _\text{prior}} \right)
\end{aligned}\]</div>

パラメータが増える(モデルの複雑性が増す)ごとに第一項は大きくなっていくかもしれないが、
第二項のペナルティも大きな負になっていく。
**中程度が良さそう → 過学習しない！**

<div class="note">

長方形じゃなくてもっとちゃんとしたGaussian近似をSection 4.4.1で
</div>

<div class="note">

[Figure 3.13](http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure3.13.png) **どうして中程度の複雑性のモデルが好まれるか**

<p><img src="http://research.microsoft.com/en-us/um/people/cmbishop/prml/prmlfigs-png/Figure3.13.png" alt="Figure 3.13" width="300px"></p>
横軸はデータセットが取りうる値を1次元で表現。
モデルの複雑性を $\mathcal M _1 < \mathcal M _2 <  \mathcal M _3$ とする。

シンプルなモデル $\mathcal M _1$ は生成(説明)できるデータの範囲が狭く
(いろいろパラメータを変えても似通ったデータセットしか出てこない)、
複雑なモデル $\mathcal M _3$ はいろんなデータを生成できるがそれぞれの重みは低い。
特定のデータセット $\mathcal D _0$ に対しては中程度の複雑さを持つモデル
$\mathcal M _2$ が一番大きいevidenceを持つことになる。
</div>

------------------------------------------------------------------------

**Expected Bayes factor**

$\mathcal M_1$ が真のモデルだとする。
ベイズ因子は個々のデータで見ると
正しくない $\mathcal M_2$ とかで大きくなる場合もあるが、
真の分布の上でを平均すると (**式 3.73**)

<div>\[
\int p(\mathcal D \mid \mathcal M _1)
    \ln \frac {p(\mathcal D \mid \mathcal M _1)}
              {p(\mathcal D \mid \mathcal M _2)} \mathrm d \mathcal D
\]</div>

で *Kullback-Leibler divergence* (Section 1.6.1 **式 1.113**)
と同じ形になり（対数の中身と符号を入れ替え）、
常に正（ただし2つの分布が等しい場合は0）の値をとることが分かっているので、
**平均的には正しいモデルのベイズ因子が大きくなり、好まれる。**
ただし、データを生成する真の分布が *L* 個のモデルの中に含まれてれば、の話。

------------------------------------------------------------------------

**まとめ**

-   Bayesian frameworkでは過学習を避け、訓練データだけでモデル比較できる
-   でもモデルの形に関する仮定は必要で、それが正しくないと誤った結論を導きうる
-   結論は事前分布の特性にかなり依存
    -   非正則事前分布では正規化定数が定義できないためevidenceを定義できない
    -   じゃあ正則事前分布の極限(e.g. 分散∞の正規分布)をとればいいかというと、
        それではevidenceが0に収束してしまう
    -   先に2つのモデルのevidence比を取ってから極限をとるといいかも
-   実際の応用では独立なテストデータを評価用に取っとくのが賢明 (←え、結局？)
