+++
title = 'RプログラミングTips'
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -97
+++

## tidyverse

<a href="https://tidyverse.org/">
<img src="/_img/hex-stickers/tidyverse.webp" style="float: right;" width="120" height="139">
</a>

解析も作図も **整然データ (tidy data)** を用意するところから始まる。

- [初学者向け講義資料2025](/slides/tohoku2025r/3-structure1.html)
- わかりやすいスライド: [整然データってなに？ by @f_nisihara](https://speakerdeck.com/fnshr/zheng-ran-detatutenani)
- 詳しい解説: [整然データとは何か by @f_nisihara](https://id.fnshr.info/2017/01/09/tidy-data-intro/)
- 原著: [Tidy Data by @hadley](https://dx.doi.org/10.18637/jss.v059.i10)

<figure>

>   *tidy datasets are all alike but every messy dataset is messy in its own way*
<figcaption>— Hadley Wickham</figcaption>
</figure>

[tidyverse](https://www.tidyverse.org/)
はそういう思想に基いて互いに連携するようデザインされたパッケージ群で、
R標準の関数よりも遥かに分かりやすく安全で高機能なものを提供してくれている。

- グラフ描画には [ggplot2]({{< relref "ggplot2.md" >}})
- data.frame内の計算・要約・抽出には [dplyr]({{< relref "dplyr.md" >}})
- data.frameの変形・ネストには [tidyr]({{< relref "tidyr.md" >}})
- リストなどに対するループ処理には [purrr]({{< relref "purrr.md" >}})
- 文字列処理には [stringr]({{< relref "stringr.md" >}})
- data.frame <=> CSV/TSV の読み書きには [readr]({{< relref "readr.md" >}})
- list <=> json の読み書きには [jsonlite](https://cran.r-project.org/web/packages/jsonlite/)

[日本語版](https://amzn.to/2yyFRKt)でも[英語版](https://amzn.to/2tbRmVc)でも[公開オンライン版](https://r4ds.hadley.nz/)でもいいのでとにかく **R for Data Science (r4ds)** を読むのが一番。
即戦力が欲しい場合は、新しく日本語で書かれた「[RユーザのためのRStudio[実践]入門−tidyverseによるモダンな分析フローの世界](https://amzn.to/2u0hmTs)」が取っつきやすい。

## 確率分布

<https://cran.r-project.org/web/views/Distributions.html>

<https://en.wikibooks.org/wiki/R_Programming/Probability_Distributions>

### 関数の種類

`d___(x, ...)`
:   確率密度関数 (PDF)。 $P[X = x]$

`p___(q, ..., lower.tail = TRUE, log.p = FALSE)`
:   累積分布関数 (CDF)。
    デフォルトでは左から`q`までの積分 $P[X \leq q]$ 。
    `lower.tail = FALSE` とすると`q`より右の積分(相補CDF) $P[X > q]$ 。
:    第一引数やパラメータ`...`はvector処理されるが、
     後半の論理値引数はvectorを与えても先頭のみが使われる。

:   離散分布の場合は境界の値を含むか含まないかでバー1本分の差が出るので注意。
    ホントは `lower.tail = FALSE` でも境界を含むべきだと思うんだけど。

:   `1 - pnorm(...)`や`log(pnorm(...))`のほうが直感的に分かりやすいので、
    `lower.tail = FALSE`や`log.p = TRUE`は不要なようにも思われるが、
    これらの引数で内部処理させたほうが浮動小数点型の限界付近での計算が正確。
    ```r
    # complementary
    1 - pnorm(10, 0, 1)                # 0
    pnorm(10, 0, 1, lower.tail = FALSE)  # 7.619853e-24

    # log
    log(pnorm(10, 0, 1))               # 0
    pnorm(10, 0, 1, log.p = TRUE)        # -7.619853e-24
    ```

`q___(p, ..., lower.tail = TRUE, log.p = FALSE)`
:   累積分布関数の逆関数。
    左からの累積確率が `p` となる確率変数の値。
    `p___()` の逆関数。
    `lower.tail = FALSE` とすると右から。

`r___(n, ...)`
:   乱数を `n` 個生成する

```r
dnorm(c(0, 1.96))
## [1] 0.39894228 0.05844094
pnorm(c(0, 1.96))
## [1] 0.5000000 0.9750021
qnorm(c(0.5, 0.975))
## [1] 0.000000 1.959964
rnorm(4)
## [1] -1.77327259  0.95713346  0.27941121  0.08387267
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
_cauchy(location = 0, scale = 1)
_chisq(df)
_exp(rate = 1)
_f(df1, df2)
_gamma(shape, rate = 1, scale = 1 / rate)
_lnorm(meanlog = 0, sdlog = 1)
_logis(location = 0, scale = 1)
_norm(mean = 0, sd = 1)
_t(df)
_unif(min = 0, max = 1)
_weibull(shape, scale = 1)
```

## 関数を作る

基本

```r
my_add = function(x, y = 1) {
   x + y
}

my_add(13, 29)  # 42
my_add(3)       # 4   (2つめの引数を省略すると1)
```

### 引数の選択肢を指定・チェック `match.arg()`

省略すると先頭の要素が採用される

```r
great_scott = function(name = c("Marty", "Emmett", "Biff")) {
    name = match.arg(name)
    print(sprintf("I am %s.", name))
}
great_scott("Marty")     # OK
great_scott("DeLorean")  # Error
great_scott()            # OK, Marty
```

### 引数の有無をチェック

困ったことに、引数不足で呼び出した瞬間にはエラーを出してくれず、
関数の中でその引数の中身を参照しようとしたところで初めてエラーが出て止まる

```r
my_func = function(x) {
    print("WTF! This line is printed anyway.")
    print(x)
}
my_func()
```

関数の中で `missing()` を使えばエラーを出さずに有無をチェックできる。
これを利用して引数省略可能な関数を作ることもできるが、
それなら素直にデフォルト引数を使ったほうがシンプルだし、
利用者は引数を見るだけで意図を読み取れる

```r
good_func = function(x, y = x) {
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

fun1(value)
## [1] 42
fun1(quote(value))
## value
fun2(value)
## value
fun3(value)
## [1] "value"
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

Rの数値や文字列は1個の値でもベクトル扱い。
同じ長さ(または長さ1)の相手との計算はベクトルのまま行える。
自分で `for` 文を書くよりこれを利用するほうが楽チンで高速。

```r
x = c(1, 2, 3)  # 長さ3の数値ベクトル
x + x           # 同じ長さ同士の計算
# [1] 2 4 6

y = 42          # 長さ1の数値ベクトル
x + y           # 長さ3 + 長さ1 = 長さ3 (それぞれ足し算)
# [1] 43 44 45

# 自分でforを書く悪い例。
# コードが無駄に長くなる上に実行速度も遅い。
z = c(0, 0, 0)
for (i in seq_len(3)) {
  z[i] = x[i] + y
}
```

四則演算のみならず様々な処理がベクトルのまま行えるようになっている。
e.g., `log()`, `sqrt()`, `ifelse()`, `paste()`


### 並列化

See [parallel]({{< relref "parallel.md" >}})

### イテレータでメモリ節約

See [parallel #iterators]({{< relref "parallel.md#iterators" >}})

### ボトルネックを知る

```r
Rprof()           # start profiling
some_hard_work()
Rprof(NULL)       # end profiling
summaryRprof()
```

### 実行時間の計測・比較

```r
r_for = function(n) {
  s = 0; for (i in seq_len(n)) {s = s + 1 / i}; s
}
r_vec = function(n) sum(1 / seq_len(n))

Rcpp::cppFunction("double rcpp(int n) {
  double s = 0; for (int i = 1; i <= n; ++i) {s += 1.0 / i;} return s;
}")  # Compilation takes a few seconds here

n = 1000000L
system.time(r_for(n))
system.time(r_vec(n))
system.time(rcpp(n))
```

計測にはr-libチームによる[bench](https://bench.r-lib.org/)が便利。
時間だけでなくメモリや実行結果までチェックした上、可視化までお世話してくれる:
```r
df = bench::mark(r_for(n), r_vec(n), rcpp(n))
df
#    expression          min         mean       median          max  itr/sec     mem_alloc  n_gc n_itr   total_time   result     memory                                             time       gc
#        <char> <bench_time> <bench_time> <bench_time> <bench_time>    <num> <bench_bytes> <num> <int> <bench_time>   <list>     <list>                                           <list>   <list>
# 1:   r_for(n)      30.49ms      31.14ms      31.11ms      32.82ms  32.1107        3.89MB     0    17        529ms 14.39273 <Rprofmem> 1: 30.7ms,30.7ms,30.8ms,31.9ms,31.1ms,31.2ms,... <tbl_df>
# 2:   r_vec(n)       2.25ms       2.85ms       2.74ms       7.87ms 350.2935       11.44MB    35    81        231ms 14.39273 <Rprofmem> 2:  7.87ms,6.67ms,6.73ms,7.65ms,2.81ms,2.8ms,... <tbl_df>
# 3:    rcpp(n)          1ms       1.02ms          1ms       1.79ms 977.6923        2.49KB     0   489        500ms 14.39273 <Rprofmem> 3: 1.79ms,1.05ms,1.04ms,1.04ms,1.03ms,1.03ms,... <tbl_df>
plot(df)
```

ちょっとした比較には
[rbenchmark](https://cran.r-project.org/package=rbenchmark)
の表示がシンプルで見やすい:
```r
rbenchmark::benchmark(r_for(n), r_vec(n), rcpp(n))[,1:4]
#       test replications elapsed relative
# 1 r_for(n)          100   3.968   29.835
# 2 r_vec(n)          100   0.473    3.556
# 3  rcpp(n)          100   0.133    1.000
```

もうちょっと詳しく見たい場合は
[microbenchmark](https://github.com/joshuaulrich/microbenchmark/)
の出力も扱いやすい:
```r
df = microbenchmark::microbenchmark(r_for(n), r_vec(n), rcpp(n), times = 100L)
df
# Unit: milliseconds
#      expr       min        lq      mean    median        uq       max neval
#  r_for(n) 35.716909 36.131130 37.406180 36.791894 38.330989 42.605700   100
#  r_vec(n)  2.630786  3.251540  4.285222  3.412421  5.573207  8.180405   100
#   rcpp(n)  1.216385  1.222994  1.270488  1.231115  1.304536  1.515532   100
str(df)
# Classes ‘microbenchmark’ and 'data.frame':      300 obs. of  2 variables:
#  $ expr: Factor w/ 3 levels "r_for(n)","r_vec(n)",..: 2 1 1 1 2 2 1 1 1 2 ...
#  $ time: num  4722392 42474783 38625124 38819035 3888012 ...
ggplot(df) + aes(expr, time) +
  stat_summary(fun = mean, geom = "bar") +
  geom_jitter(height = 0, alpha = 0.5)
```

## Links

- <https://adv-r.hadley.nz/>
