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
[RやPythonなどいろんなインターフェイスで利用可能](https://mc-stan.org/interfaces/)。
[RStan](https://mc-stan.org/users/interfaces/rstan.html),
[PyStan](https://mc-stan.org/users/interfaces/pystan.html)
が長らく使われてきたが、
[CmdStanR](https://mc-stan.org/cmdstanr/),
[CmdStanPy](https://mc-stan.org/cmdstanpy/)
への移行が進んできている。

https://mc-stan.org/


## インストール

RやPythonのパッケージを入れてから、それ越しにCmdStan本体を入れる。

```r
install.packages("cmdstanr", repos = c(stan = "https://stan-dev.r-universe.dev"))
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan()
cmdstan_path()
cmdstan_version()
```

```py
%pip3 install cmdstanpy
import cmdstanpy
cmdstanpy.install_cmdstan()
cmdstanpy.cmdstan_path()
cmdstanpy.cmdstan_version()
cmdstanpy.show_versions()
```


## 基本的な流れ

1. cmdstanrを読み込む
   ```r
   library(cmdstanr)
   ```
1. 名前付きlistとしてデータを用意する。
   e.g., 平均10のポアソン乱数。
   ```r
   sample_size = 1000L
   mydata = list(N = sample_size, x = rpois(sample_size, 10))
   ```
1.  Stan言語でモデルを記述する。
    RStanには文字列で渡せたがCmdStanPy, CmdStanRは別ファイル必須。
    e.g., 与えられたデータがポアソン分布から取れてきたとすると、
    その平均はどれくらいだったか？
    ```stan
    data {
      int<lower=0> N;
      array[N] int<lower=0> x;
    }

    parameters {
      real<lower=0> lambda;
    }

    model {
      x ~ poisson(lambda);
    }
    ```
1. モデルをC++に変換してコンパイルする。
   中間ファイルは `*.hpp`。
   ```r
   model = cmdstan_model("model.stan")
   ```
   <https://mc-stan.org/cmdstanr/reference/cmdstan_model.html>
1. コンパイル済みモデルを使ってMCMCサンプリング
   ```r
   fit = model$sample(mydata)
   ```
   <https://mc-stan.org/cmdstanr/reference/model-method-sample.html>
1.  結果を見てみる
    ```r
    print(fit)
    fit$summary()
    fit$cmdstan_summary()
    fit$cmdstan_diagnose()
    fit$sampler_diagnostics()
    fit$diagnostic_summary()
    fit$metadata()
    ```
    <https://mc-stan.org/cmdstanr/reference/CmdStanMCMC.html>
1.  MCMCサンプルを使う。
    ```r
    draws_df = fit$draws(format = "df")
    draws = fit$draws()
    params = names(model$variables()$parameters)
    bayesplot::mcmc_acf_bar(draws, pars = params)
    bayesplot::mcmc_trace(draws, pars = params)
    bayesplot::mcmc_hist(draws, pars = params)
    bayesplot::mcmc_combo(draws, pars = params)
    rhat = bayesplot::rhat(fit)
    neff = bayesplot::neff_ratio(fit)
    bayesplot::mcmc_rhat(rhat)
    bayesplot::mcmc_neff(neff)
    ```
    <https://mc-stan.org/cmdstanr/reference/fit-method-draws.html>


## Stan文法

https://mc-stan.org/documentation/

### ブロック

コード内に登場できるブロックは7種類で、省略可能だが順番はこの通りでなければならない。

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
: モデルで使いやすい形に値を変換しておくとか。
  ここに書いた変数もサンプリングされる。

`model {...}`
: 唯一の必須ブロック。
  サンプルされないローカル変数を宣言してもよいが、制約をかけることはできない。

`generated quantities {...}`
: `normal_rng()`などによる乱数生成が許される唯一のブロック。
  観察値の確信区間とかを

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

- スカラーは整数(`int`) or 実数(`real`)。
- 実数ベクトル(`vector`, `row_vector`)と実数行列(`matrix`)は
  `Eigen::Vector` や `Eigen::Matrix` で実装されているので効率的に行列演算を行える。
- 配列(`array`)は `std::vector` で実装されていて、
  整数配列や行列配列など何でも作れるが、行列演算はできないので生の`for`ループが必要。
  [`int v[3]` のように宣言する書き方は非推奨になった。](https://mc-stan.org/docs/reference-manual/brackets-array-syntax.html)
- 宣言時に上限下限を設定できる (constrained integer/real)。
- bool型は無くて基本的に整数の1/0。分岐ではnon-zeroがtrue扱い。

```stan
int i;
real x;
array[42] int a;
array[42] real y;
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

正規分布のsigmaやポアソン分布のlambdaの値域をちゃんと `real<lower=0>` に絞っているのに
`Scale parameter is 0, but must be positive!`
と怒られることがある。
実害はないけどどうしても警告を消したい場合は
[違うseedを使うとか `step_size = 0.1` のように歩幅を狭めるとかで対処できる](https://discourse.mc-stan.org/t/scale-parameter-is-0-but-must-be-0-can-i-do-anything-to-deal-with-this/19453)。

## `library(bayesplot)`

<https://mc-stan.org/bayesplot/>


## `library(posterior)`

<https://mc-stan.org/posterior/>


## `library(rstanarm)`

<https://mc-stan.org/rstanarm/>

R標準のGLMのような使い心地でStanを動かせるようにするパッケージ。

- formulaでモデルを立てられる。
- data.frameを渡せる。
- パラメータ調整やコンパイルの済んだ部品を組むような形なので試行錯誤が早い。
- ただしcmdstanrではなくrstanを使う。


## `library(brms)`

<https://paul-buerkner.github.io/brms/>

rstanarmと同様にformula形式でStanを動かせるようにする。
相違点:

- [stan-dev](https://github.com/stan-dev)チームの一員ではない。
- コンパイル済みの部品を使わずStanコードを生成する。そのぶん遅いが柔軟。
- CmdStanRをbackendとして使える。


## 関連書籍

<a href="https://www.amazon.co.jp/dp/4320112423/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=a55a7a3616d8d8a6d4516c5a26bca46f&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4320112423&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4320112423" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4065165369/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=75cdd63ab18e3b59eaeac9e628b27ce8&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4065165369&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4065165369" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
