+++
date = "2016-05-27T18:31:56+09:00"
tags = ["r", "hadley"]
title = "purrr"
subtitle = "apply系関数の究極形"
[menu.main]
  parent = "rstats"
  weight = -73
+++

https://github.com/hadley/purrr

標準のapply族や[plyr]({{< relref "plyr.md" >}})の関数に引導を渡すべく再設計されたようなパッケージ。
[dplyr]({{< relref "dplyr.md" >}}) や [tidyr]({{< relref "tidyr.md" >}}) と組み合わせて真価を発揮する。
purerのように読むらしい。

Rの中で `install.packages('purrr')` としてインストールし、
使う前に `library(purrr)` で読み込む。

パイプ演算子 `%>%` については[dplyr]({{< relref "dplyr.md" >}})を参照。

## list, vector操作

`purrr::map(.x, .f, ...)`
: listの各要素に関数を適用し、listで返す。
  型の決まったvectorを返す亜種として
  `map_lgl()`, `map_chr()`, `map_int()`, `map_dbl()`, `map_df()`
  がある。

`purrr::map_if(.x, .p, .f, ...)`
: `.p`が`TRUE`になる要素のみ`.f()`を適用し、残りはそのまま出力。
  `.p`はlogical vectorでもいいし、
  `.x[[i]]`を受け取るpredicate関数でもよい。

`purrr::map2(.x, .y, .f, ...)`
: 2変数バージョン。`map2_chr()` や `map3()` なども提供されている。

`purrr::walk(.x, .f, ...)`
: mapと同じように関数を適用するが、元の値をそのままinvisible返し。

`purrr::lmap(.x, .f, ...)`, `lmap_if()`, `lmap_at()`
: listの各要素に関数を適用し、listで返す。
  `.x[[i]]`ではなく`.x[i]`を参照する点で`map()`と異なる。

`purrr::reduce()`
: 二変数関数を順々に適用して1つの値を返す。
  C++でいう`std::accumulate()`。

`purrr::flatten(.x)`
: 階層性のあるlistを一段階解消する。
  階層なしlistをvector化するには明示的に型指定できる
  `flatten_lgl()`, `flatten_int()`, `flatten_dbl()`, `flatten_chr()`, `flatten_df()`
  のほうが標準の`unlist()`よりも安心。

`purrr::keep(.x, .p, ...)`, `discard()`, `compact()`
: listやvectorの要素を `.p` に応じて取捨選択。
  わざわざ関数にするほどでもない気もするが、パイプラインで使いやすい形。

`purrr::split_by(.x, .f, ...)`, `order_by()`, `sort_by()`
: listやvectorを `.f` に応じて切ったり並べ替えたり。

## data.frame操作

`purrr::slice_rows(.d, .cols=NULL)`
: 指定した列でグループ化して `grouped_df` を返す。
  `dplyr::group_by_(.dots=.cols)` と同じ。

`purrr::by_slice(.d, ..f, ..., .collate=c('list', 'rows', 'cols'), .to='.out', .labels=TRUE)`
: `grouped_df` を受け取ってグループごとに関数を適用する。
  `dplyr::do()` の改良版。

`purrr::by_row(.d, ..f, ..., .collate=c('list', 'rows', 'cols'), .to='.out', .labels=TRUE)`
: 1行ごとに関数を適用する。
  `dplyr::rowwise() %>% dplyr::do()`的な処理を一撃で書ける。
: `.collate`: 結果列の展開方法。
  とりあえず`list`にしておいて後で`unnest()`するのが無難か。
: `.to`: 結果listの列名
: `.labels`: 元の`.d`の列をラベルとして結果に残すか

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

## 関数

無名関数
: ラムダ式とも言う。
  チルダとドットで簡単に書ける。

```r
ord = function(x) {strtoi(charToRaw(x), 16L)}
letters %>>% map_int(ord)
letters %>>% map_int(function(x) {strtoi(charToRaw(x), 16L)})
letters %>>% map_int(~ strtoi(charToRaw(.), 16L))
```

`purrr::partial()`
: 引数を部分的に埋めてある関数を作る。C++でいう `std::bind()`
