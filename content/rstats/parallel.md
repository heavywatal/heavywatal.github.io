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

lapply(seq(1, 4), slow_square) |> timeit()
```

```
[1] 1.2
```

``` r
purrr::map(seq(1, 4), slow_square) |> timeit()
```

```
[1] 1.2
```

``` r
parallel::mclapply(seq(1, 4), slow_square) |> timeit()
```

```
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

foreach(x = seq(1, 4), .combine = c) %do% {
  slow_square(x)
} |> timeit()
```

```
[1] 1.2
```

``` r
foreach(x = seq(1, 4)) %dopar% {
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
