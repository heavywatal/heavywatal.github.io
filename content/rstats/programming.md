+++
title = 'RプログラミングTips'
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -90
+++

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

## 1列でもmatrixのまま

matrixから添字で1行だけ取ると、自動でvectorにしてくれる。
大抵はそれでいいんだけど、そうしてほしくないときは `drop=FALSE`

```r
> m = matrix(1:8, 2, 4)
> m
     [,1] [,2] [,3] [,4]
[1,]    1    3    5    7
[2,]    2    4    6    8
> m[1, ]
[1] 1 3 5 7
> m[1, , drop=FALSE]
     [,1] [,2] [,3] [,4]
[1,]    1    3    5    7
```

## 文字列をRコマンドとして実行

ショートカットとしてこのような関数を用意しておいて

```r
eval_parse = function(...){
    eval.parent(parse(text=paste0(...)))
}
```

変数や関数を規則的に定義するとか

```r
> v = seq_len(4)
> eval_parse('x', v, ' = ', v ^ 2)
> c(x1, x2, x3, x4)
[1]  1  4  9 16
```

## `attach()` するべからず

データフレーム `some_data` のある要素 `some_item` にアクセスしたいとき、
普通はダラーを使って `some_data$some_item` のようにする。
`attach(some_data)` すると、データフレーム内のすべての要素が独立した変数みたいになって、
`some_item` だけでアクセスできるようになる (`detach()` するまで)。
でもそうすると、ワークスペースは変数でいっぱいになってしまうし、
名前の衝突による不具合が起こる危険性がある。
そこで、範囲限定で一時的に `attach()` 状態を作り出す `with()` を使う。
安全で楽チン。

```r
half = iris$Sepal.Length / 2

with(iris, {
    ## このブレースの中では
    ## attach(iris)されたのと同じように、
    ## 列の名前でダイレクトにアクセスできる。

    half = Sepal.Length / 2
})
```

## 最適化・高速化

### 最新の R を使う

-   2.11.0 64bit化
-   2.14.0 すべての標準パッケージがバイトコンパイルされている
-   2.15.? data.frameの操作が劇的に速くなった

### ベクトル化されてる関数・演算子を使う

ナマの `for` 文や `if` 文を避けるのと同義

```r
> vec = seq_len(1000000)
> ifelse(vec %% 3,
+        ifelse(vec %% 5, vec, 'buzz'),
+        ifelse(vec %% 5, 'fizz', 'fizzbuzz'))
```

### list, data.frame, matrix

-   [dplyr]({{< relref "dplyr.md" >}}), [purrr]({{< relref "purrr.md" >}}),
    [tidyr]({{< relref "tidyr.md" >}}) を介して操作すると楽チンかつ高速。
-   基本的にlistやdata.frameは遅いので、
    matrixで済むもの(数値のみの表など)はmatrixで。

### 並列化

<http://cran.r-project.org/web/views/HighPerformanceComputing.html>

`snow` とか `multicore` が使われてきたが、
バージョン2.14から `parallel` が標準ライブラリに入った。

<http://www.rdocumentation.org/packages/parallel>

CPUコア数を取得する

```r
> library(parallel)
> parallel::detectCores()
[1] 4
```

[plyr]({{< relref "plyr.md" >}}) にも `apply()` 的な処理を簡単に並列化する機能がある。

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

### バイトコンパイル

より機械語に近い、構文解析を済ませた形まで関数を変換させておく (2.13以降)

```r
library(compiler)
cfun = compiler::cmpfun(my_func)
```

------------------------------------------------------------------------

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/1593273843/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/511RD-C-K%2BL._SX160_.jpg" alt="The Art of R Programming: A Tour of Statistical Software Design" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4873115795/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51sGGydPkYL._SX160_.jpg" alt="アート・オブ・Rプログラミング" /></a>
