+++
title = 'RプログラミングTips'
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -97
+++

## tidyverse

<a href="https://tidyverse.org/">
<img src="https://tidyverse.tidyverse.org/logo.png" align="right" width="120" height="139">
</a>

解析も作図も**整然データ (tidy data)**を用意するところから始まる。

- わかりやすいスライド: [整然データってなに？ --- @f_nishihara](https://speakerdeck.com/fnshr/zheng-ran-detatutenani)
- 詳しい解説: [整然データとは何か --- @f_nishihara](http://id.fnshr.info/2017/01/09/tidy-data-intro/)
- 原著: [Tidy Data --- @hadley](https://dx.doi.org/10.18637/jss.v059.i10)

>   *tidy datasets are all alike but every messy dataset is messy in its own way*\
>   --- *Hadley Wickham*

[tidyverse](https://www.tidyverse.org/)
はそういう思想に基いて互いに連携するようデザインされたパッケージ群で、
R標準の関数よりも遥かに分かりやすく安全で高機能なものを提供してくれている。

<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?s=books&ie=UTF8&qid=1508340700&sr=1-3&linkCode=li1&tag=heavywatal-22&linkId=c21926e67c143d462d77fc8995212796" target="_blank"><img align="right" border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL110_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img align="right" src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li1&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%A7%E3%81%AF%E3%81%98%E3%82%81%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9-Hadley-Wickham/dp/487311814X/ref=as_li_ss_il?ie=UTF8&qid=1508340144&sr=8-1&keywords=r&linkCode=li1&tag=heavywatal-22&linkId=f08ff14b91c2edb1a0b323a3f9460e34" target="_blank"><img align="right" border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL110_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img align="right" src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li1&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />

- グラフ描画には [ggplot2]({{< relref "ggplot2.md" >}})
- data.frame内の計算・要約・抽出には [dplyr]({{< relref "dplyr.md" >}})
- data.frameの変形・ネストには [tidyr]({{< relref "tidyr.md" >}})
- リストなどに対するループ処理には [purrr]({{< relref "purrr.md" >}})
- 文字列処理には [stringr]({{< relref "stringr.md" >}})
- data.frame <=> CSV/TSV の読み書きには [readr]({{< relref "readr.md" >}})
- list <=> json の読み書きには [jsonlite](https://cran.r-project.org/web/packages/jsonlite/)

[日本語版](https://amzn.to/2yyFRKt)でも[英語版](https://amzn.to/2tbRmVc)でも[公開オンライン版](https://r4ds.had.co.nz/)でもいいのでとにかく **R for Data Science (r4ds)** を読むのが一番。

## 確率分布

<https://cran.r-project.org/web/views/Distributions.html>

<https://en.wikibooks.org/wiki/R_Programming/Probability_Distributions>

### 関数の種類

`d___(x, ...)`
:   確率密度関数 (PDF)。 $P[X = x]$

`p___(q, ..., lower.tail=TRUE, log.p=FALSE)`
:   累積分布関数 (CDF)。
    デフォルトでは左から`q`までの積分 $P[X \leq q]$ 。
    `lower.tail=FALSE` とすると`q`より右の積分(相補CDF) $P[X > q]$ 。
:    第一引数やパラメータ`...`はvector処理されるが、
     後半の論理値引数はvectorを与えても先頭のみが使われる。

    {{%div class="warning"%}}
離散分布の場合は境界の値を含むか含まないかでバー1本分の差が出るので注意。
ホントは `lower.tail=FALSE` でも境界を含むべきだと思うんだけど。
    {{%/div%}}

    {{%div class="note"%}}
`1 - pnorm(...)`や`log(pnorm(...))`のほうが直感的に分かりやすいので、
`lower.tail=FALSE`や`log.p=TRUE`は不要なようにも思われるが、
これらの引数で内部処理させたほうが浮動小数点型の限界付近での計算が正確。
```r
# complementary
1 - pnorm(10, 0, 1)                # 0
pnorm(10, 0, 1, lower.tail=FALSE)  # 7.619853e-24

# log
log(pnorm(10, 0, 1))               # 0
pnorm(10, 0, 1, log.p=TRUE)        # -7.619853e-24
```
    {{%/div%}}

`q___(p, ..., lower.tail=TRUE, log.p=FALSE)`
:   累積分布関数の逆関数。
    左からの累積確率が `p` となる確率変数の値。
    `p___()` の逆関数。
    `lower.tail=FALSE` とすると右から。

`r___(n, ...)`
:   乱数を `n` 個生成する

```r
> dnorm(c(0, 1.96))
[1] 0.39894228 0.05844094
> pnorm(c(0, 1.96))
[1] 0.5000000 0.9750021
> qnorm(c(0.5, 0.975))
[1] 0.000000 1.959964
> rnorm(4)
[1] -1.77327259  0.95713346  0.27941121  0.08387267
```

### 分布の種類

離散

```r
_binom(size, prob)
_geom(prob)
_hyper(m, n, k)
_nbinom(size, prob, mu)
_pois(lambda)
_signrank(n)
_wilcox(m, n)
```

連続

```r
_beta(shape1, shape2)
_cauchy(location=0, scale=1)
_chisq(df)
_exp(rate=1)
_f(df1, df2)
_gamma(shape, rate=1, scale=1/rate)
_lnorm(meanlog=0, sdlog=1)
_logis(location=0, scale=1)
_norm(mean=0, sd=1)
_t(df)
_unif(min=0, max=1)
_weibull(shape, scale=1)
```

## 関数を作る

基本

```r
my_add = function(x, y=1) {
   x + y
}

my_add(13, 29)  # 42
my_add(3)       # 4   (2つめの引数を省略すると1)
```

### 引数の選択肢を指定・チェック `match.arg()`

省略すると先頭の要素が採用される

```r
great_scott = function(name=c('Marty', 'Emmett', 'Biff')) {
    name = match.arg(name)
    print(sprintf('I am %s.', name))
}
great_scott('Marty')     # OK
great_scott('DeLorean')  # Error
great_scott()            # OK, Marty
```

### 引数の有無をチェック

困ったことに、引数不足で呼び出した瞬間にはエラーを出してくれず、
関数の中でその引数の中身を参照しようとしたところで初めてエラーが出て止まる

```r
my_func = function(x) {
    print('WTF! This line is printed anyway.')
    print(x)
}
my_func()
```

関数の中で `missing()` を使えばエラーを出さずに有無をチェックできる。
これを利用して引数省略可能な関数を作ることもできるが、
それなら素直にデフォルト引数を使ったほうがシンプルだし、
利用者は引数を見るだけで意図を読み取れる

```r
good_func = function(x, y=x) {
    x + y
}

bad_func = function(x, y) {
    if (missing(y)) {y = x}
    x + y
}

good_func(3)  # 6
bad_func(3)   # 6
```

### 可変長引数 (`...`)

最後の引数をピリオド３つにし、`list(...)` か `c(...)` で受け取る

```r
func = function(x, y, ...){
    for (i in c(...)) {
        cat(i, "\n")
    }
}
```

### 引数を中身ではなく名前として評価

```r
value = 42
fun1 = function(x) x
fun2 = function(x) substitute(x)
fun3 = function(x) deparse(substitute(x))

> fun1(value)
[1] 42
> fun1(quote(value))
value
> fun2(value)
value
> fun3(value)
[1] "value"
```

詳しくは non-standard expression (NSE) で調べる:

- https://adv-r.hadley.nz/meta
- https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
- https://rlang.tidyverse.org/articles/tidy-evaluation.html


## 最適化・高速化

### 最新の R を使う

-   2.11.0 64bit化
-   2.14.0 すべての標準パッケージがバイトコンパイルされている
-   2.15.? data.frameの操作が劇的に速くなった

### ベクトル化されてる関数・演算子を使う

Rのナマの`for`ループは絶望的に遅い。
四則演算のみならず様々な処理がベクトルのまま行えるようになっているので、そっちを使う。

```r
vec = seq_len(60)
vec + vec
vec * 2

ifelse(vec %% 3,
       ifelse(vec %% 5, vec, 'buzz'),
       ifelse(vec %% 5, 'fizz', 'fizzbuzz'))
```

### list, data.frame, matrix

-   [tidyverse](https://tidyverse.tidyverse.org/) パッケージ群
    ([dplyr]({{< relref "dplyr.md" >}}), [purrr]({{< relref "purrr.md" >}}),
    [tidyr]({{< relref "tidyr.md" >}})など) を介して操作すると楽チンかつ高速。
-   基本的にlistやdata.frameは遅いので、
    matrixで済むもの(数値のみの表など)はmatrixで。

### 並列化

See [foreach #parallel]({{< relref "foreach.md#parallel" >}})

### イテレータでメモリ節約

See [foreach #iterators]({{< relref "foreach.md#iterators" >}})

### ボトルネックを知る

```r
> Rprof()           # start profiling
> some_hard_work()
> Rprof(NULL)       # end profiling
> summaryRprof()
```

### 実行時間の計測

```r
> system.time(測定したいコマンド)
```


## 関連書籍

統計解析やグラフ描画ではなく、
プログラミング言語としてのRを学びたいときに:

<a href="https://www.amazon.co.jp/RStudio%E3%81%A7%E3%81%AF%E3%81%98%E3%82%81%E3%82%8BR%E3%83%97%E3%83%AD%E3%82%B0%E3%83%A9%E3%83%9F%E3%83%B3%E3%82%B0%E5%85%A5%E9%96%80-Garrett-Grolemund/dp/4873117151/ref=as_li_ss_il?ie=UTF8&qid=1489745137&sr=8-1&keywords=rstudio&linkCode=li3&tag=heavywatal-22&linkId=59bf9fd2a28700d591c1bb951fcf6bd4" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117151&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117151" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/%E3%82%A2%E3%83%BC%E3%83%88%E3%83%BB%E3%82%AA%E3%83%96%E3%83%BBR%E3%83%97%E3%83%AD%E3%82%B0%E3%83%A9%E3%83%9F%E3%83%B3%E3%82%B0-Norman-Matloff/dp/4873115795/ref=as_li_ss_il?ie=UTF8&qid=1485613704&sr=8-1&keywords=%E3%82%A2%E3%83%BC%E3%83%88%E3%82%AA%E3%83%96r&linkCode=li3&tag=heavywatal-22&linkId=b322dc8f7f7dc364e861086b0d53a10b" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873115795&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873115795" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
