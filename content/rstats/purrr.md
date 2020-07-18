+++
date = 2016-05-27T18:31:56+09:00
tags = ["r", "tidyverse"]
title = "purrr"
subtitle = "ループ処理やapply系関数の決定版"
[menu.main]
  parent = "rstats"
  weight = -73
+++

<a href="https://purrr.tidyverse.org/">
<img src="http://hexb.in/vector/purrr.svg" align="right" width="120" height="139">
</a>

forループやlistの処理などをより簡潔に書けるようにしてくれるパッケージ。
標準のapply系関数よりも覚えやすく読みやすい。
[dplyr]({{< relref "dplyr.md" >}}) や [tidyr]({{< relref "tidyr.md" >}}) と組み合わせて使う。
いまのところ並列化する機能はないので、
それに関しては[foreach/parallel]({{< relref "foreach.md" >}})ページを参照。

[tidyverse](https://tidyverse.tidyverse.org/) に含まれているので、
`install.packages("tidyverse")` で一括インストール、
`library(tidyverse)` で一括ロード。

## list, vector操作

### 各要素に関数を適用するapply系関数

```r
v = list(1, 2L, "3")
check_class = function(x) {paste0(x, " is ", class(x))}

# 自分でfor文を書くと結構大変
results = vector("list", length(v))
for (i in seq_along(v)) {
  results[[i]] = check_class(v[[i]])
}

# 1行で簡潔に記述でき、意図も明確
results = lapply(v, check_class)
results = purrr::map(v, check_class)
results = purrr::map(v, function(x) {paste0(x, " is ", class(x))})
results = purrr::map(v, ~ {paste0(.x, " is ", class(.x))})
```

`purrr::map(.x, .f, ...)`
: list/vector `.x` の各要素に関数 `.f` を適用した結果をlistに詰めて返す。
  `base::lapply()`とほぼ同義。
: `.f` にformulaや数値などを渡すと[関数に変換した上で処理してくれる。(後述)](#無名関数)

`purrr::map_lgl()`, `map_int()`, `map_dbl()`, `map_chr()`
: 型の決まったvectorを返すmap亜種。
  `base::sapply()`や`base::vapply()`よりも覚えやすく読みやすい。

`purrr::map_dfr(.x, .f, ..., .id = NULL)`
: 結果を`dplyr::bind_rows()`で結合したdata.frameとして返すmap亜種。
  例えば、同じ形式のCSVファイルを一気に読んで結合、みたいなときに便利:

```r
files = fs::dir_ls("path/to/data/", glob = "*.csv")
combined_df = purrr::map_dfr(files, readr::read_csv)
```

`purrr::map2(.x, .y, .f, ...)`
: 2変数バージョン。`map2_int()`などの型指定vector版もある。
  3変数以上渡したいときはlistかdata.frameにまとめて次の`pmap()`を使う。

`purrr::pmap(.l, .f, ...)`
: listの中身をparallelに処理するmap。
  関数 `.f`の引数はlistの要素名と一致させるか `...` で受け流す必要がある。
  e.g., `pmap(list(a = 1:3, b = 4:6), function(a, b) {a * b})` 。
: data.frameの正体はlist of columnsなので、
  そのまま`.l`として渡せる。
: `pmap_int` など出力型指定の亜種もある。

`purrr::map_if(.x, .p, .f, ...)`
: `.p`が`TRUE`になる要素のみ`.f()`を適用し、残りはそのまま出力。
  `.p`はlogical vectorでもいいし、
  `.x[[i]]`を受け取るpredicate関数でもよい。
: 番号か名前で選ぶには`purrr::map_at(.x, .at, .f, ...)`

`purrr::walk(.x, .f, ...)`
: `map()`同様に関数を適用しつつ元の値をそのままinvisible返しする亜種。

`purrr::lmap(.x, .f, ...)`
: `.x[[i]]`ではなく`.x[i]`を参照する亜種。

`purrr::imap(.x, .f, ...)`
:  名前や整数インデックスを第二引数で受け取れる亜種。
  `iwalk()` などの派生もある。

`purrr::modify(.x, .f, ...)`
: 入力と同じ型で出力する亜種。
  つまりdata.frameを入れたらlistじゃなくてdata.frameが出てくる。
  e.g., `diamonds %>% modify_if(is.numeric, round)`

`purrr::reduce(.x, .f, ..., .init)`
: 二変数関数を順々に適用して1つの値を返す。
  C++でいう`std::accumulate()`。
  例えば <code>reduce(1:3, \`+\`)</code> の結果は6。

`purrr::accumulate(.x, .f, ..., .init)`
: 二変数関数を順々に適用し、過程も含めてvectorで返す。
  C++でいう`std::partial_sum()`。
  例えば `accumulate(1:3, sum)` の結果は `1 3 6` 。


### list作成・変形・解体

`purrr::list_along(x)`
: `x` と同じ長さの空listを作る `vector("list", length(x))` のショートカット。

`purrr::flatten(.x)`
: 階層性のあるlistを一段階解消する。
  階層なしlistをvector化するには明示的に型指定できる
  `flatten_lgl()`, `flatten_int()`, `flatten_dbl()`, `flatten_chr()`, `flatten_dfr()`
  のほうが標準の`unlist()`よりも安心。

`purrr::keep(.x, .p, ...)`, `discard()`, `compact()`
: listやvectorの要素を `.p` に応じて取捨選択。
  `.p` に関数を渡した場合の挙動は
  `.x[.p(.x)]` じゃなくて `.x[map_lgl(.x, .p, ...)]` となることに注意。

`purrr::pluck(.x, ..., .default = NULL)`
: オブジェクト `.x` 内の要素を引っ張り出す `[[` の強力版。
  `...` には整数、文字列、関数、listで複数指定できる。
  例えば `accessor(x[[1]])$foo` だと読む順が左右に振られるが、
  `pluck(x, 1, accessor, "foo")` だと左から右に読み流せる。

`purrr::cross2(.x, .y, .filter = NULL)`
: listの各要素の組み合わせを作る。
  `.filter` に渡した関数が `TRUE` となるものは除外される。
  名前付きlistを渡す `purrr::cross()` や `purrr::cross_df()` のほうが便利かも。
  vectorなら [`tidyr::crossing()` とか `tidyr::expand()`]({{< relref "tidyr.md" >}}) が使える。

`purrr::transpose(.l)`
: 行列転置関数`t()`のlist版。
  例えば、pair of lists <=> list of pairs。
: data.frameに適用するとlist of rowsが得られる。


### その他

`purrr::invoke(.f, .x = NULL, ..., .env = NULL)`
: list `.x` の中身を引数として関数 `.f` を呼び出す。
: 関数に渡す引数があらかじめlistにまとまってるときに使う`do.call()`の改良版。

```r
params = list(n = 6L, size = 10L, replace = TRUE)
purrr::invoke(sample.int, params)
```

`purrr::invoke_map(.f, .x = list(NULL), ..., .env = NULL)`
: 関数listを順々に実行してlistで返す。
  引数`.x`は同じ長さのlist of listsか、list of a listをリサイクル。
: e.g., `invoke_map(list(runif, rnorm), list(c(n = 3, 0, 1)))`

`purrr::has_element(.x, .y)`
: list `.x` は要素 `.y` を持っている。

`purrr::set_names(x, nm = x)`
: 標準の`setNames(x = nm, nm)`は第二引数のほうが省略不可という気持ち悪い定義だったが、
  この改良版ではその心配が解消されている。
  長さや型のチェックもしてくれる。


## 無名関数

apply/map系関数は、名前のついた関数だけでなく、その場で定義された無名関数も受け取れる。
ごく短い関数や一度しか使わない関数に名前をつけずに済むので便利。
さらにpurrrのmap系関数はformula
(チルダで始まる `~ x + y` のようなもの)
や数値を受け取って関数として処理してくれる。

```r
# named function
ord = function(x) {strtoi(charToRaw(x), 16L)}
map_int(letters, ord)

# unnamed function
map_int(letters, function(x) {strtoi(charToRaw(x), 16L)})

# formula
map_int(letters, ~ strtoi(charToRaw(.x), 16L))

# integer/character
li = list(lower = letters, upper = LETTERS)
map_chr(li, 3L)
map_chr(li, function(x) {x[[3L]]})
```

formula内部では、第一引数を`.x`または`.`として、第二引数を`.y`として参照する。
`..1`, `..2`, `..3` のような形で三つめ以降も参照できる。

`purrr::as_mapper(.f, ...)`
: `map()` 内部で関数への変換機能を担っている関数。
: formulaを受け取ると `function (.x, .y, . = .x)` のような関数に変換する。
: 数値や文字列を受け取ると `[[` による抽出関数に変換する。
  参照先が存在しない場合の値はmap関数の `.default` 引数で指定できる。

`purrr::partial(...f, ..., .env, .lazy, .first)`
: 引数を部分的に埋めてある関数を作る。C++でいう `std::bind()`


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

パイプ演算子 `%>%` については[dplyr]({{< relref "dplyr.md" >}})を参照。

```r
## OLD
diamonds %>%
  purrrlyr::slice_rows("cut") %>%
  purrrlyr::by_slice(head, .collate = "rows")

## NEW
diamonds %>%
  tidyr::nest(-cut) %>%
  dplyr::mutate(data = purrr::map(data, head)) %>%
  tidyr::unnest(data)
```


`purrrlyr::dmap(.d, .f, ...)`
: data.frameかgrouped_dfを受け取り、列ごとに`.f`を適用してdata.frameを返す。
  対象列を選べる`dmap_if()`と`dmap_at()`もある。
  全列で長さが揃っていれば怒られないので
  `dplyr::mutate_all()`的にも`dplyr::summarise_all()`的にも使える。
  ただし`.f`に渡せるのは単一の関数のみで`dplyr::funs()`は使えない。
  パッケージ作者は `dplyr` の `mutate_*()` や `summarise_*()` の利用を推奨。

`purrrlyr::slice_rows(.d, .cols = NULL)`
: 指定した列でグループ化してgrouped_dfを返す。
  `dplyr::group_by_(.dots = .cols)` と同じ。

`purrrlyr::by_slice(.d, ..f, ..., .collate = c("list", "rows", "cols"), .to = ".out", .labels = TRUE)`
: grouped_dfを受け取ってグループごとに関数を適用する。
  `dplyr::do()` とほぼ同じ役割で、一長一短。
  こちらは出力形式をより柔軟に指定できるが、
  中の関数からgrouping variableを参照できないという弱点を持つ。

`purrrlyr::by_row(.d, ..f, ..., .collate = c("list", "rows", "cols"), .to = ".out", .labels = TRUE)`
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


## 関連書籍

<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?s=books&ie=UTF8&qid=1508340700&sr=1-3&linkCode=li3&tag=heavywatal-22&linkId=6a53371cc80e2d6d7fc50fae5d8b862d" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%A7%E3%81%AF%E3%81%98%E3%82%81%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9-Hadley-Wickham/dp/487311814X/ref=as_li_ss_il?ie=UTF8&qid=1508340144&sr=8-1&keywords=r&linkCode=li3&tag=heavywatal-22&linkId=4137d3d3351f8ccab5a93cefdc28fdec" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
