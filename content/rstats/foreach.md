+++
title = 'foreach/parallel'
tags = ["r", "concurrent"]
[menu.main]
  parent = "rstats"
  weight = -55
+++

## foreach

カウンター無しでループを書けるようにするパッケージ。
普段は [`purrr::map()`]({{< relref "purrr.md" >}}) とかのほうが使いやすいけど、
[後述の並列化](#parallel) の場面ではお世話になる。

```r
foreach (mu = seq_len(8), .combine=c) %do% {
    rnorm(1, mu)
}
```

https://CRAN.R-project.org/package=foreach

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


## parallel

- https://cran.r-project.org/web/views/HighPerformanceComputing.html
- https://stat.ethz.ch/R-manual/R-patched/library/parallel/html/00Index.html

Rの並列化では `snow` や `multicore` が使われてきたが、
バージョン2.14からそれらを統合した `parallel` が標準ライブラリに入った。
直接触るのは難しいので、`foreach`(とその橋渡しライブラリ
[`doParallel`](https://CRAN.R-project.org/package=doParallel))
を介して使う。

```r
library(doParallel)
cores = detectCores(logical=FALSE)
cluster = makeCluster(cores, 'FORK')
registerDoParallel(cluster)
foreach (mu = seq_len(8), .combine=c) %dopar% {
    rnorm(1, mu)
}
stopCluster(cluster)
```

クラスタの生成や破棄などをぜーんぶ自動でやってもらいたい場合は
[`hoxo-m/pforeach`](https://github.com/hoxo-m/pforeach)。

[`purrr::map()`]({{< relref "purrr.md" >}})的な使い勝手の
[`map_par()` を書いてみた](https://github.com/heavywatal/rwtl/blob/master/R/parallel.R)。

### `makeCluster()`

`spec`
: いくつのworkerを立ち上げるか。
  CPUコア数を取得するには `parallel::detectCores(logical=FALSE)`

`type='PSOCK'`
: デフォルト。高コストだけどだいたいどの環境でも使える。
  マルチCPUのサーバーで並列化したい場合はこれ。
  `foreach()` で使う場合 `.export=` や `.packages=` の指定が重要。

`type='FORK'`
: 4コア1CPUとかの普通のデスクトップマシンで気楽に並列化したいならこれ。
  低コストだし `.export=` や `.packages=` を指定せず `foreach()` できる。
  Windowsでは使えないらしいけど。

`outfile=''`
: `print()`や`message()`などの出力先を標準に戻す。
  デフォルトでは`/dev/null`に捨てられてしまう。


## iterators

https://CRAN.R-project.org/package=iterators

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

`iter(obj, by=c('column', 'row'))`
: イテレータ版 `purrr::by_row()` のようなもので、
  並列`foreach`で各ノードに巨大data.frameを送りたくない場合に有用。
  data.frame以外も適用可。
  e.g., `iter(iris, by='row')`

`isplit(x, f, drop=FALSE, ...)`
: イテレータ版 `purrr::slice_rows()` のようなもので、
  `f`は列名じゃなくてfactor。
  data.frame以外も適用可。
  e.g., `isplit(iris, iris$Species)`

`iread.table(file, ..., verbose=FALSE)`, `ireadLines(con, n=1, ...)`
: ファイルを1行ずつ読み込む

`irbinom(..., count)`, `irnbinom()`, `irnorm()`, `irpois()`, `irunif()`
: 乱数

`idiv(n, ..., chunks, chunkSize)`
: 整数`n`をいい感じに振り分ける。


## 関連書籍

<a href="https://www.amazon.co.jp/Parallel-R-Q-Ethan-McCallum-ebook/dp/B005Z29QT4/ref=as_li_ss_il?ie=UTF8&qid=1477817549&sr=8-1&keywords=r+parallel&linkCode=li3&tag=heavywatal-22&linkId=6250767635d9095cf1e3a6bace362a83" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=B005Z29QT4&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=B005Z29QT4" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%AB%E3%82%88%E3%82%8B%E3%83%8F%E3%82%A4%E3%83%91%E3%83%95%E3%82%A9%E3%83%BC%E3%83%9E%E3%83%B3%E3%82%B9%E3%82%B3%E3%83%B3%E3%83%94%E3%83%A5%E3%83%BC%E3%83%86%E3%82%A3%E3%83%B3%E3%82%B0-%E7%A6%8F%E5%B3%B6-%E7%9C%9F%E5%A4%AA%E6%9C%97/dp/4883379353/ref=as_li_ss_il?ie=UTF8&qid=1477817567&sr=8-1&keywords=r+%E4%B8%A6%E5%88%97&linkCode=li3&tag=heavywatal-22&linkId=7f12f8273082a42443269e875a2c282b" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4883379353&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4883379353" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
