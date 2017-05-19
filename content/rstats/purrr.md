+++
date = "2016-05-27T18:31:56+09:00"
tags = ["r", "tidyverse"]
title = "purrr"
subtitle = "ループ処理やapply系関数の決定版"
[menu.main]
  parent = "rstats"
  weight = -73
+++

標準のapply族や[plyr]({{< relref "plyr.md" >}})の関数に引導を渡すべく再設計されたようなパッケージ。
[dplyr]({{< relref "dplyr.md" >}}) や [tidyr]({{< relref "tidyr.md" >}}) と組み合わせて真価を発揮する。
purerのように読むらしい。
いまのところ並列化する機能はないので、そういうときは
[foreachとparallel]({{< relref "foreach.md" >}}) を使う。

[tidyverse](https://github.com/tidyverse/tidyverse) に含まれているので、
`install.packages('tidyverse')` で一括インストール、
`library(tidyverse)` で一括ロード。

パイプ演算子 `%>%` については[dplyr]({{< relref "dplyr.md" >}})を参照。

https://github.com/tidyverse/purrr

## list, vector操作

`purrr::map(.x, .f, ...)`
: listの各要素に関数を適用し、listで返す。
  `base::lapply()`や`plyr::llply()`の改良版。
  型の決まったvectorを返す`base::sapply()`や`base::vapply()`の改良版として
  `map_lgl()`, `map_chr()`, `map_int()`, `map_dbl()` がある。
: `purrr::walk(.x, .f, ...)`は、
  `map()`同様に関数を適用しつつ元の値をそのままinvisible返しする亜種。
: `purrr::lmap(.x, .f, ...)`は、
  `.x[[i]]`ではなく`.x[i]`を参照する点で`map()`と異なる亜種。

`purrr::map_df(.x, .f, ..., .id=NULL)`
: listの各要素に関数を適用し、結果を
  `dplyr::bind_rows()`で結合したdata.frameとして返す。
  `plyr::ldply()`の改良版として。

`purrr::map2(.x, .y, .f, ...)`
: 2変数バージョン。`map2_chr()`などの型限定vector版も提供されている。
  3変数以上渡したいときはlistかdata.frameにまとめて下の`pmap()`を使う。

`purrr::pmap(.l, .f, ...)`
: listの中身をparallelに処理するmap。
  `.f`の引数はlistの要素名と一致させる。
  e.g., `pmap(list(a=1:3, b=4:6), function(a, b) {a * b})` 。
  名前無し2要素までなら`map2()`などと同様に`.x`や`.y`を使ったformulaが利用可能。
  e.g., `pmap(list(1:3, 4:6), ~ .x * .y)` 。
: data.frameの正体はlist of columnsなので、
  そのまま`.l`として渡して`plyr::mlply()`的に使える。
  同様に `purrr::pmap_df()` は `plyr::mdply()` として使える。

`purrr::map_if(.x, .p, .f, ...)`
: `.p`が`TRUE`になる要素のみ`.f()`を適用し、残りはそのまま出力。
  `.p`はlogical vectorでもいいし、
  `.x[[i]]`を受け取るpredicate関数でもよい。
: 番号か名前で選ぶには`purrr::map_at(.x, .at, .f, ...)`

`purrr::reduce(.x, .f, ..., .init)`
: 二変数関数を順々に適用して1つの値を返す。
  C++でいう`std::accumulate()`。

`purrr::accumulate(.x, .f, ..., .init)`
: 二変数関数を順々に適用し、過程も含めてvectorで返す。
  例えば `accumulate(1:3, sum)` の結果は `1 3 6` 。
  C++でいう`std::partial_sum()`。

`purrr::invoke(.f, .x=NULL, ..., .env=NULL)`
: 関数に渡す引数があらかじめlistにまとまってるときに使う`do.call()`の改良版。

`purrr::invoke_map(.f, .x=list(NULL), ..., .env=NULL)`
: 関数listを順々に実行してlistで返す。
  引数`.x`は同じ長さのlist of listsか、list of a listをリサイクル。
: e.g., `invoke_map(list(runif, rnorm), list(c(n=3, 0, 1)))`
: 返り値がすべて同じ型で、それぞれ長さ1だった場合、
  特定の型のvectorとして返せる亜種がある:
  `invoke_map_chr()`, `invoke_map_dbl()`, `invoke_map_int()`, `invoke_map_lgl()`
: data.frameを繋げて返す`invoke_map_df()`も。

`purrr::flatten(.x)`
: 階層性のあるlistを一段階解消する。
  階層なしlistをvector化するには明示的に型指定できる
  `flatten_lgl()`, `flatten_int()`, `flatten_dbl()`, `flatten_chr()`, `flatten_df()`
  のほうが標準の`unlist()`よりも安心。

`purrr::keep(.x, .p, ...)`, `discard()`, `compact()`
: listやvectorの要素を `.p` に応じて取捨選択。
  `.p` に関数を渡した場合の挙動は
  `.x[.p(.x)]` じゃなくて `.x[map_lgl(.x, .p, ...)]` となることに注意。

`purrr::split_by(.x, .f, ...)`, `order_by()`, `sort_by()`
: listやvectorを `.f` に応じて切ったり並べ替えたり。


## 無名関数

ラムダ式とも言う。
チルダで始まる `~ a + b` のようなものはRではformulaという形式になるが、
`purrr`のmap系関数はそれを引数`.f`として受け取ることができる。
formula内部では、第一引数を`.x`または`.`として、第二引数を`.y`として参照する。

```r
# with named function
ord = function(x) {strtoi(charToRaw(x), 16L)}
letters %>>% map_int(ord)

# with unnamed function
letters %>>% map_int(function(x) {strtoi(charToRaw(x), 16L)})

# with formula
letters %>>% map_int(~ strtoi(charToRaw(.x), 16L))
```

`purrr::as_function(.f, ...)`
: 上記の機能を担う重要な関数で、`map()`内部でも使われている。
  formulaを`function (.x, .y, . = .x)`という無名関数に変換する。

`purrr::partial(...f, ..., .env, .lazy, .first)`
: 引数を部分的に埋めてある関数を作る。C++でいう `std::bind()`


## その他の便利関数

`purrr::list_along(x)`
: 引数と同じ長さの空listを作る。

`purrr::set_names(x, nm=x)`
: 標準の`setNames(x=nm, nm)`は第二引数のほうが省略不可という気持ち悪い定義だったが、
  この改良版ではその心配が解消されている。
  長さや型のチェックもしてくれる。

`purrr::transpose(.l)`
: 行列転置関数`t()`のlist版。
  例えば、pair of lists <=> list of pairs。
: data.frameに適用するとlist of rowsが得られる。


## `purrrlyr`

data.frame を引数にとるものは purrr 0.2.2.1 から切り離され、
[purrrlyr](https://github.com/hadley/purrrlyr) に移動された。
**これらは今のところdeprecatedではないが、
近いうちにそうなるので早くほかのアプローチに移行せよ** 、とのこと。
https://github.com/hadley/purrrlyr/blob/master/NEWS.md

[`tidyr`]({{< relref "tidyr.md" >}}) でネストして、
[`purrr`]({{< relref "purrr.md" >}}) でその list of data.frames に処理を施し、
[`dplyr`]({{< relref "dplyr.md" >}}) でその変更を元の data.frame に適用する、
というのがtidyverse流のモダンなやり方らしい。

```r
## OLD
iris %>%
  purrrlyr::slice_rows('Species') %>%
  purrrlyr::by_slice(head, .collate='rows')

## NEW
iris %>%
  tidyr::nest(-Species) %>%
  dplyr::mutate(data= purrr::map(data, head)) %>%
  tidyr::unnest()
```


`purrrlyr::dmap(.d, .f, ...)`
: data.frameかgrouped_dfを受け取り、列ごとに`.f`を適用してdata.frameを返す。
  対象列を選べる`dmap_if()`と`dmap_at()`もある。
  全列で長さが揃っていれば怒られないので
  `dplyr::mutate_all()`的にも`dplyr::summarise_all()`的にも使える。
  ただし`.f`に渡せるのは単一の関数のみで`dplyr::funs()`は使えない。
  パッケージ作者は `dplyr` の `mutate_*()` や `summarise_*()` の利用を推奨。

`purrrlyr::slice_rows(.d, .cols=NULL)`
: 指定した列でグループ化してgrouped_dfを返す。
  `dplyr::group_by_(.dots=.cols)` と同じ。

`purrrlyr::by_slice(.d, ..f, ..., .collate=c('list', 'rows', 'cols'), .to='.out', .labels=TRUE)`
: grouped_dfを受け取ってグループごとに関数を適用する。
  `dplyr::do()` とほぼ同じ役割で、一長一短。
  こちらは出力形式をより柔軟に指定できるが、
  中の関数からgrouping variableを参照できないという弱点を持つ。

`purrrlyr::by_row(.d, ..f, ..., .collate=c('list', 'rows', 'cols'), .to='.out', .labels=TRUE)`
: data.frame 1行ごとに関数を適用する。
  `dplyr::rowwise() %>% dplyr::do()`的な処理を一撃で書ける。
: `.to`: 結果listの列名
: `.labels`: 元の`.d`の列をラベルとして結果に残すか
: `.collate`: 結果列の展開方法(下記例)。
  とりあえず`list`にしておいて後で`unnest()`するのが無難か。
: `purrrlyr::invoke_rows()`はかなり似ているが、
  データと関数の順序が逆になっている点と、
  `.f`が`.d`の列名で引数を取るという点で異なる。
  (`by_row()`では1行のdata.frameとして受け取って`.$col`のように参照する)

`purrrlyr::invoke_rows(.f, .d, ..., .collate, .to, .labels)`
: 第一引数が関数になってるところは`invoke()`っぽくて、
  結果の収納方法を`.collate`などで調整できるところは`by_row()`っぽい。


```r
> iris[1:3, 1:2] %>% by_row(~data_frame(x=0:1, y=c('a', 'b')), .collate='list')
Source: local data frame [3 x 3]

  Sepal.Length Sepal.Width           .out
         (dbl)       (dbl)          (chr)
1          5.1         3.5 <tbl_df [2,2]>
2          4.9         3.0 <tbl_df [2,2]>
3          4.7         3.2 <tbl_df [2,2]>
> iris[1:3, 1:2] %>% by_row(~data_frame(x=0:1, y=c('a', 'b')), .collate='rows')
Source: local data frame [6 x 5]

  Sepal.Length Sepal.Width  .row     x     y
         (dbl)       (dbl) (int) (int) (chr)
1          5.1         3.5     1     0     a
2          5.1         3.5     1     1     b
3          4.9         3.0     2     0     a
4          4.9         3.0     2     1     b
5          4.7         3.2     3     0     a
6          4.7         3.2     3     1     b
> iris[1:3, 1:2] %>% by_row(~data_frame(x=0:1, y=c('a', 'b')), .collate='cols')
Source: local data frame [3 x 6]

  Sepal.Length Sepal.Width    x1    x2    y1    y2
         (dbl)       (dbl) (int) (int) (chr) (chr)
1          5.1         3.5     0     1     a     b
2          4.9         3.0     0     1     a     b
3          4.7         3.2     0     1     a     b
```


## 関連書籍

<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?_encoding=UTF8&qid=1485613345&sr=1-1-catcorr&linkCode=li3&tag=heavywatal-22&linkId=6133a1fd9babbf590e304ee9f670fa4a" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
