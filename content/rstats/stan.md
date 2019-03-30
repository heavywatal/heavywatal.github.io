+++
date = 2016-06-23T13:06:06+09:00
tags = ["r", "c++"]
title = "Stan"
subtitle = "高速MCMCでパラメータ推定"

[menu.main]
  parent = "rstats"
  weight = 1
+++

数あるMCMCアルゴリズムの中でも効率的なHMC(Hybrid/Hamiltonian Monte Carlo)を用いてベイズ推定を行うツール。
[Pythonやコマンドラインなどいろんな形で利用可能](http://mc-stan.org/interfaces/)だが、
とりあえずRで[RStan](http://mc-stan.org/interfaces/rstan.html)を使ってみる。

http://mc-stan.org/

## インストール

Rから`install.packages('rstan')`で一発。
jagsと違ってstan本体も同時に入れてくれる。
[RStan-Getting-Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
を見ると、時代や環境によってはいろいろ難しいかったのかも。

標準的な(Xcode Command Line Tools とか build-essential的な)開発環境はどっちみち必要。

## 基本的な流れ

1. rstanを読み込む
   ```r
   library(rstan)
   rstan_options(auto_write = TRUE)
   options(mc.cores = parallel::detectCores())
   ```

1. 名前付きlistとしてデータを用意する。
   e.g., 平均10、標準偏差3の正規乱数。
   ```r
   .data = list(x=rnorm(10000, 10, 3))
   .data$n_obs = length(.data$x)
   ```

1.  Stan言語でモデルを記述する。
    別ファイルにしてもいいし、下記のようにR文字列でもいい。
    e.g., 与えられたデータが正規分布から取れてきたとすると、
    その平均と標準偏差はどれくらいだったか？
    ```r
    .stan_code = '
    data {
      int n_obs;
      real[n_obs] x;
    }

    parameters {
      real mu;
      real<lower=0> sigma;
    }

    model {
      x ~ normal(mu, sigma);
    }'
    ```

1. モデルをC++に変換してコンパイルする。
   ファイルから読み込んだ場合は中間ファイル`*.rda`がキャッシュされる。
   ```r
   .model = rstan::stan_model(model_code=.stan_code)
   # or
   .model = rstan::stan_model(file='model.stan')
   ```

1. コンパイル済みモデルを使ってMCMCサンプリング
   ```r
   .fit = rstan::sampling(.model, data=.data, iter=10000, chains=3)
   ```

1. 結果を見てみる
   ```r
   print(.fit)
   summary(.fit)
   plot(.fit)
   pairs(.fit)
   rstan::traceplot(.fit)
   rstan::stan_trace(.fit)
   rstan::stan_hist(.fit)
   rstan::stan_dens(.fit)
   ```

## Stan文法

http://mc-stan.org/documentation/
PDFしか無くて残念

### ブロック

コード内に登場できるブロックは7種類で、順番はこの通りでなければならない。

`functions {...}`
: 関数を定義できる。

`data {...}`
: Rから受け取る定数の宣言。

`transformed data {...}`
: 定数の宣言と代入。
  決め打ちのハイパーパラメータとか。
  決定論的な変換のみ可能。

`parameters {...}`
: サンプリングされる変数の宣言。

`transformed parameters {...}`
: 変数の宣言と代入。
  モデルで使いやすい形にパラメータを変形しておくとか？

`model {...}`
: 唯一の必須ブロック。
  サンプルされないローカル変数を宣言してもよいが、制約をかけることはできない。

`generated quantities {...}`
: サンプリング後の値を使って好きなことをするとこ？
  `normal_rng()`などによる乱数生成が許される唯一のブロック。
  rstanならここを使わずRで結果を受け取ってからどうにかするほうが簡単？

### モデリング

あるパラメータにおけるlog probabilityと近傍での傾きを計算し、
それらを元に次の値にジャンプする、という操作が繰り返される。
modelブロック内で暗黙的に定義されている `target` 変数に対して
`+=` 演算子で対数確率をどんどん加算していく。
(昔は隠れ変数`lp__`や`increment_log_prob()`などを使ってた。)

サンプリング文(sampling statement)はそれを簡単に記述するためのショートカット。
名前とは裏腹に、確率分布からのサンプリングが行われるわけではないので紛らわしい。
例えば以下の表現はほぼ等価。
(定数の扱い方がうまいとかでサンプリング文のほうが効率的らしいけど)

```stan
x ~ normal(0.0, 1.0);
target += normal_lpdf(x | 0.0, 1.0);
target += -0.5 * square(x);
```

確率分布としての正規化はうまいことやっといてくれるから気にしなくていいらしい
(が、`T[,]`によるtruncated distributionではこうやって調整する、
とかいう記述もあるので、そのへんはまだよく分からない)。

<div>\[\begin{aligned}
\log p(x) &\propto -\frac {x^2} 2 \\
     p(x) &\propto \exp \left(- \frac {x^2} 2 \right)
\end{aligned}\]</div>

名のある確率分布はだいたい関数として用意されている。
形のバリエーションとしては:

- 確率密度関数: `*_lpdf(y | ...)`, `*_lpmf(y | ...)`
- 累積分布関数: `*_cdf(y | ...)`, `*_lcdf(y | ...)`
- 相補累積分布関数: `*_lccdf(y | ...)`
- 乱数生成: `*_rng(...)`

(対数版のsuffixは昔は `_cdf_log()`, `_ccdf_log()` という形だった)


### 型

整数(`int`)、実数(`real`)、実数ベクトル(`vector`, `row_vector`)、実数行列(`matrix`)。
内部的に `Eigen::Vector` や `Eigen::Matrix` が使われているので、
可能な限り`for`文よりも行列演算を使うように心がける。
配列(array)は `std::vector` で実装されていて、
整数配列や行列配列など何でも作れるが、行列演算はできない。

宣言時に上限下限を設定できる (constrained integer/real)。

bool型は無くて基本的に整数の1/0。分岐ではnon-zeroがtrue扱い。

```stan
int i;
int v[42];
real x;
real x[42];
int<lower=1,upper=6> dice;

vector[3] v;
row_vector[3] r;
matrix[3, 3] m;

x * v  // vector[3]
r * v  // real
v * r  // matrix[3, 3]
m * v  // vector[3]
m * m  // matrix[3, 3]
m[1]   // row_vector[3]
```

そのほかの特殊な制約つきの型

- `simplex`: 合計が1になる非負実数ベクトル
- `unit_vector`: 二乗和が1になる実数ベクトル
- `ordered`, `positive_ordered`:
  昇順実数ベクトル。降順にしたければ `transformed parameters` ブロックで。
- `cov_matrix`, `corr_matrix`, `cholesky_factor_cov`, `cholesky_factor_corr`


### Tips

条件分岐するときはなるべく`if`文を避けて三項演算子やステップ関数を使うべし、
という言語が多いけどStanでは逆に`if`文を素直に書くほうが良いらしい。
`if_else()`では真値でも両方の引数が評価されちゃうし、
`step()` や `int_step()` からの掛け算は遅いのだとか。

代入演算子は普通に `=` イコール。(昔は `<-` 矢印だった)

対数尤度の値を確認したいときは `print("log_prob: ", target())`

## 可視化

https://www.rdocumentation.org/packages/rstan/topics/Plots

```r
stan_plot()
stan_trace()
stan_scat()
stan_hist()
stan_dens()
stan_ac()

# S3 method
pairs()
print()
```

`stanfit` クラスのmethodとして `plot()` や `traceplot()` が定義されているが、
いくつかのチェックとともに `stan_plot()` 系の関数を呼び出すだけで大きな違いは無さそう。


## トラブル対処

### StanHeaders version is ahead of rstan version

Stanのヘッダーライブラリとrstanは別々のパッケージで提供されていて、
Stan更新への追従にタイムラグがあるらしい。
こんなん開発者側でどうにかして欲しいけど、
とりあえず古い `StanHeaders` を入れてしのぐしかない。
https://github.com/stan-dev/rstan/wiki/RStan-Transition-Periods

```r
install.packages("https://cran.r-project.org/src/contrib/Archive/StanHeaders/StanHeaders_2.9.0.tar.gz", repos=NULL, type='source')
```

https://cran.r-project.org/src/contrib/Archive/StanHeaders/

### 最新版をGitHubからインストール

リポジトリの構造が標準とはちょっと違う
```r
remotes::install_github('stan-dev/rstan', ref='develop', subdir='rstan/rstan')
```

## 関連書籍

<a href="https://www.amazon.co.jp/Stan%E3%81%A8R%E3%81%A7%E3%83%99%E3%82%A4%E3%82%BA%E7%B5%B1%E8%A8%88%E3%83%A2%E3%83%87%E3%83%AA%E3%83%B3%E3%82%B0-Wonderful-R-%E6%9D%BE%E6%B5%A6-%E5%81%A5%E5%A4%AA%E9%83%8E/dp/4320112423/ref=as_li_ss_il?_encoding=UTF8&psc=1&refRID=N4GQH6P8QHKYX0E91BVK&linkCode=li3&tag=heavywatal-22&linkId=6ffd2c8744eace95b257403231a29102" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4320112423&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4320112423" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/gp/product/1482253445/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=9367946532fbf4369b6161aedd15c4a6" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1482253445&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1482253445" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
