+++
title = 'parallel'
subtitle = "並行処理 in R"
aliases = ["foreach.html"]
tags = ["r", "concurrent"]
[menu.main]
  parent = "rstats"
  weight = -55
+++



## parallel

- https://cran.r-project.org/web/views/HighPerformanceComputing.html
- https://stat.ethz.ch/R-manual/R-patched/library/parallel/html/00Index.html

Rの並列化では `snow` や `multicore` が使われてきたが、
バージョン2.14からそれらを統合した `parallel` が標準ライブラリに入った。


``` r
options(mc.cores = parallel::detectCores(logical = FALSE))

timeit = function(expr, digits = 1) {
  start = Sys.time()
  expr
  diff = Sys.time() - start
  as.numeric(round(diff, digits = digits))
}

slow_square = function(x) {
  Sys.sleep(0.3)
  x * x
}

lapply(seq_len(4L), slow_square) |> timeit()
purrr::map(seq_len(4L), slow_square) |> timeit()
parallel::mclapply(seq_len(4L), slow_square) |> timeit()
```

```
[1] 1.2
[1] 1.2
[1] 0.3
```

[`mclapply()`](https://stat.ethz.ch/R-manual/R-patched/library/parallel/html/mclapply.html)
は `lapply()` のお手軽並列化バージョン。
UNIX系OSのforkに依存するためWindows不可。

```r
mclapply(X, FUN, ...,
         mc.preschedule = TRUE, mc.set.seed = TRUE,
         mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
         mc.cleanup = TRUE, mc.allow.recursive = TRUE, affinity.list = NULL)
```

[`purrr::map()`]({{< relref "purrr.md" >}})のように無名関数を渡せる
[ラッパー関数 `mcmap()` を書いてみた](https://github.com/heavywatal/rwtl/blob/master/R/parallel.R)。
そのへんをもっとしっかりやった
[future](https://github.com/HenrikBengtsson/future),
[future.apply](https://github.com/HenrikBengtsson/future.apply),
[furrr](https://github.com/DavisVaughan/furrr)
を使っていくのが良さそう。

もっと細かくいろいろ制御したい場合は後述の `foreach` を介して使う。


### `makeCluster()`

`spec`
: いくつのworkerを立ち上げるか。
  物理コア数を取得するには `parallel::detectCores(logical = FALSE)`

`type = "PSOCK"`
: デフォルト。高コストだけどだいたいどの環境でも使える。
  マルチCPUのサーバーで並列化したい場合はこれ。
  `foreach()` で使う場合 `.export=` や `.packages=` の指定が重要。

`type = "FORK"`
: 4コア1CPUとかの普通のデスクトップマシンで気楽に並列化したいならこれ。
  低コストだし `.export=` や `.packages=` を指定せず `foreach()` できる。
  Windowsでは使えないらしいけど。

`outfile = ""`
: `print()`や`message()`などの出力先を標準に戻す。
  デフォルトでは`/dev/null`に捨てられてしまう。


### 乱数と再現性

並列化せず親プロセスで実行すれば当然再現性あり:

``` r
rint = function(..., n = 65535L, size = 3L) sample.int(n, size)

set.seed(19937)
lapply(seq_len(4L), rint) |> simplify2array()
set.seed(19937)
lapply(seq_len(4L), rint) |> simplify2array()
```

```
      [,1]  [,2]  [,3]  [,4]
[1,] 49270 65407  8721 23119
[2,] 19176 60768 65517 62677
[3,] 29429 29809  5210  9959
      [,1]  [,2]  [,3]  [,4]
[1,] 49270 65407  8721 23119
[2,] 19176 60768 65517 62677
[3,] 29429 29809  5210  9959
```

デフォルト `mc.set.seed = TRUE` で並列化すると子プロセスがてんでにシードを設定して再現性なし:

``` r
set.seed(19937)
parallel::mclapply(seq_len(4L), rint) |> simplify2array()
set.seed(19937)
parallel::mclapply(seq_len(4L), rint) |> simplify2array()
```

```
      [,1]  [,2]  [,3]  [,4]
[1,] 62234  3262 12411 62023
[2,] 62793  4436 53581 20557
[3,]  8474 55511  8828 14334
      [,1]  [,2]  [,3]  [,4]
[1,] 44555 57701 37457 48634
[2,] 35817  3149 40988 23286
[3,] 10319 34748 35695 30145
```

`mc.set.seed = FALSE` では親プロセスのシードが参照されて再現性こそあるものの、
いくつかの子プロセス(`mc.cores`の数ずつ？)がセットで同じ乱数列を出してくる:

``` r
invisible(runif(1))
set.seed(19937)
parallel::mclapply(seq_len(6L), rint, mc.set.seed = FALSE, mc.cores = 2L) |> simplify2array()
```

```
      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 49270 49270 65407 65407  8721  8721
[2,] 19176 19176 60768 60768 65517 65517
[3,] 29429 29429 29809 29809  5210  5210
```

`RNGkind("L'Ecuyer-CMRG")` を設定する、
もしくは `set.seed(19937, "L'Ecuyer-CMRG")` のようにシード設定すると
`mc.set.seed = TRUE` の挙動が変わる。
再現性は親プロセスの乱数生成器の状態依存となり、
なおかつ子プロセスがそれぞれ異なる乱数列を出してくれるようになる。
親プロセスの乱数生成器の状態が `mclapply()` や子プロセスによって更新されないことに注意:

``` r
RNGkind("L'Ecuyer-CMRG")
set.seed(19937)
parallel::mclapply(seq_len(4L), rint) |> simplify2array()
parallel::mclapply(seq_len(4L), rint) |> simplify2array()
invisible(runif(1))
parallel::mclapply(seq_len(4L), rint) |> simplify2array()
set.seed(19937)
parallel::mclapply(seq_len(4L), rint) |> simplify2array()
```

```
      [,1]  [,2]  [,3]  [,4]
[1,] 33309 19174 22386 46754
[2,] 38519 53475 45877 24809
[3,] 39904 60478 42220 45717
      [,1]  [,2]  [,3]  [,4]
[1,] 33309 19174 22386 46754
[2,] 38519 53475 45877 24809
[3,] 39904 60478 42220 45717
      [,1]  [,2]  [,3]  [,4]
[1,] 17241 10459 58255 18031
[2,] 39348 54919 19631 64554
[3,] 65371 12315 16237  4515
      [,1]  [,2]  [,3]  [,4]
[1,] 33309 19174 22386 46754
[2,] 38519 53475 45877 24809
[3,] 39904 60478 42220 45717
```

## foreach

<https://CRAN.R-project.org/package=foreach>

カウンター無しでループを書けるようにするパッケージ。
普段は [`purrr::map()`]({{< relref "purrr.md" >}}) とかのほうが使いやすいけど、
前述の `mclapply()` よりも並列化を細かく制御したい場面ではお世話になる。
この場合、橋渡しライブラリ
[`doParallel`](https://CRAN.R-project.org/package=doParallel)
の助けを借りつつクラスタの生成や破棄なども自前でやることになる。


``` r
library(foreach)
library(doParallel)
cluster = makeCluster(getOption("mc.cores", 2L), type = "FORK")
registerDoParallel(cluster)

foreach(x = seq_len(4L), .combine = c) %do% {
  slow_square(x)
} |> timeit()
```

```
[1] 1.2
```

``` r
foreach(x = seq_len(4L)) %dopar% {
  slow_square(x)
} |> timeit()
```

```
[1] 0.3
```

``` r
stopCluster(cluster)
```

### `foreach()`

`.combine (list)`
: 型が既知でvectorが欲しい場合に `c` にするなど

`.multicombine (FALSE)`
: 結果が出る度に二値関数で結合していくか、まとめてか。
  `.combine=c` や `cbind` を指定すると暗黙`TRUE`。

`.maxcombine (100)`
: まとめる場合の最大個数

`.export (NULL)`
: 並列化する場合、処理ブロック内に持ち込みたいオブジェクト

`.packages (NULL)`
: 並列化する場合、名前空間省略で使いたいパッケージ

`.inorder (TRUE)`
: 並列化する場合、順序を保持したいか

`.init`, `.final`, `.noexport`, `.verbose`


## iterators

<https://CRAN.R-project.org/package=iterators>

大抵はメモリを一気に確保してしまう方が速いが、
データがRAMを超えるほど大きいときはそうも言ってられない。
最大要求メモリを減らしたり、
並列`foreach`のノード間通信を減らすためには`iterators`を利用する。

`nextElem(it, ...)`
: イテレータを進めて値を得る。
  デバッグ時にとりあえず全部見たいときは `as.list(it)` が便利。

`icount(n)`
: イテレータ版 `seq_len()`

`icountn(vn)`
: 自然数限定イテレータ版 `expand.grid()`

`iter(obj, by = c("column", "row"))`
: イテレータ版 `purrrlyr::by_row()` のようなもので、
  並列`foreach`で各ノードに巨大data.frameを送りたくない場合に有用。
  data.frame以外も適用可。
  e.g., `iter(diamonds, by = "row")`

`isplit(x, f, drop = FALSE, ...)`
: イテレータ版 `purrrlyr::slice_rows()` のようなもので、
  `f`は列名じゃなくてfactor。
  data.frame以外も適用可。
  e.g., `isplit(diamonds, diamonds$cut)`

`iread.table(file, ..., verbose = FALSE)`, `ireadLines(con, n = 1, ...)`
: ファイルを1行ずつ読み込む

`irbinom(..., count)`, `irnbinom()`, `irnorm()`, `irpois()`, `irunif()`
: 乱数

`idiv(n, ..., chunks, chunkSize)`
: 整数`n`をいい感じに振り分ける。
