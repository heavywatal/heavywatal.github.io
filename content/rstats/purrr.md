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
<img src="/_img/hex-stickers/purrr.webp" align="right" width="120" height="139">
</a>

forループやlistの処理などをより簡潔に書けるようにしてくれるパッケージ。
標準のapply系関数よりも覚えやすく読みやすい。
[dplyr]({{< relref "dplyr.md" >}}) や [tidyr]({{< relref "tidyr.md" >}}) と組み合わせて使う。
いまのところ並列化する機能はないので、
それに関しては[parallel]({{< relref "parallel.md" >}})ページを参照。

[tidyverse](https://tidyverse.tidyverse.org/) に含まれているので、
`install.packages("tidyverse")` で一括インストール、
`library(tidyverse)` で一括ロード。

## list, vector操作

### 各要素に関数を適用するapply系関数

```r
library(conflicted)
library(tidyverse)

v = list(1, 2L, "3")
check_class = function(x) {paste0(x, " is ", class(x))}

# 自分でfor文を書くと結構大変
results = vector("list", length(v))
for (i in seq_along(v)) {
  results[[i]] = check_class(v[[i]])
}

# 1行で簡潔に記述でき、意図も明確
results = lapply(v, check_class)
results = v |> purrr::map(check_class)
results = v |> purrr::map(function(x) {paste0(x, " is ", class(x))})
results = v |> purrr::map(\(x) {paste0(x, " is ", class(x))})
results = v |> purrr::map(~ {paste0(.x, " is ", class(.x))})
```

`purrr::map(.x, .f, ...)`
: list/vector `.x` の各要素に関数 `.f` を適用した結果をlistに詰めて返す。
  `base::lapply()`とほぼ同義。
: `.f` にformulaや数値などを渡すと[関数に変換した上で処理してくれる。(後述)](#関数-f-として渡せるもの)

`purrr::map_vec(.x, .f, ..., .ptype = NULL)`
: listではなくvectorを返すmap亜種。
: 型は推定してもらえるが `.ptype = integer()` のように指定も可能。
: `base::vapply()` と違って型省略可能。
  `base::sapply()` と違って常にvectorを返す。
: 関数名で型指定する `purrr::map_lgl()`, `map_int()`, `map_dbl()`, `map_chr()` もある。

`purrr::map2(.x, .y, .f, ...)`
: 2変数バージョン。
  3変数以上渡したいときはlistかdata.frameにまとめて次の`pmap()`を使う。
: ほかの `map_*()` 亜種にも同様に提供されている。

`purrr::pmap(.l, .f, ...)`
: listの中身をparallelに処理するmap。
  関数 `.f` の引数はlistの要素名と一致させるか `...` で受け流す必要がある。
  e.g., `pmap(list(a = 1:3, b = 4:6), function(a, b) {a * b})` 。
: data.frameの正体はlist of columnsなのでそのまま`.l`として渡せる。
: ほかの `map_*()` 亜種にも同様に提供されている。

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
  e.g., `diamonds |> modify_if(is.numeric, round)`

`purrr::modify_tree(x, ..., leaf = identity, is_node = NULL, pre = identity, post = identity)`
: listを再帰的に巡りながら関数を適用する。
  最末端のleaf/sinkノードだけに適用するなら `leaf`,
  内側も含めて全ノードに適用するなら `pre` / `post` を使う
  (前者は下りながら `leaf` 適用前、後者は `leaf` 適用後に上りながら)。


### list要素の抽出・変更

`purrr::pluck(.x, ..., .default = NULL)`
: オブジェクト `.x` 内の要素を引っ張り出す `[[` の強力版。
  `...` には整数、文字列、関数、listで複数指定できる。
  例えば `accessor(x[[1]])$foo` だと読む順が左右に振られるが、
  `pluck(x, 1, accessor, "foo")` だと左から右に読み流せる。
: 存在を問うだけの　`pluck_exists()` もある。
: `pluck()<-` も提供されているので代入も可能。
  この用途には
  `purrr::modify_in(.x, .where, .f, ...)` や
  `purrr::assign_in(x, where, value)` もある。

`purrr::chuck(.x ,...)`
: 存在しない場合にエラー終了してくれる `pluck()` 亜種。

`purrr::list_assign(.x, ..., .is_node = NULL)`
: listに対して `dplyr::mutate()` するような感じ。
  変更対象の要素がlistとして存在していても普通に上書き。

`purrr::list_modify(.x, ..., .is_node = NULL)`
: `list_assign()` と似てるが、変更対象の既存listで言及されなかった要素を変更しない、という点で異なる
: e.g., `list(a = list(x = 1)) |> list_modify(a = list(y = 2))` で要素 `x` がそのまま保持されて要素 `a` は長さ2になる。

`purrr::list_merge(.x, ..., .is_node = NULL)`
: `list_modify()` と似てるが、変更対象の既存要素に上書きせずappendする。
: e.g., `list(a = list(x = 1)) |> list_modify(a = list(x = 2))` で要素 `x` が長さ2になる。

`purrr::keep(.x, .p, ...)`, `discard()`, `compact()`
: listやvectorの要素を `.p` に応じて取捨選択。
  `.p` に関数を渡した場合の挙動は
  `.x[.p(.x)]` じゃなくて `.x[map_lgl(.x, .p, ...)]` となることに注意。
: 名前を指定する亜種 `keep_at()`, `discard_at()` もある。

`purrr::some(.x, .p, ...)`, `purrr::every()`, `purrr::none()`
: `.p(.x[[i]])` が {少なくともひとつ `TRUE`, すべて `TRUE`, すべて `FALSE`} なら `TRUE` を返す。

`purrr::has_element(.x, .y)`
: list `.x` は要素 `.y` を持っている。 `some(.x, identical, .y)`


### list変形・解体

`purrr::list_flatten(x, ..., name_spec, name_repair)`
: 階層性のあるlistを浅いほうから一段階解消する。
  結果がすべて整数とかでも勝手にvector化せず、常にlistを返す。
  `unlist(x, recursive = FALSE) |> as.list()` のようなもの。

`purrr::list_simplify(x, ..., strict = TRUE, ptype = NULL)`
: listを一段階解消して同じ長さのvectorにする。
  入力と出力の対応が保たれるので `dplyr::mutate()` の中とかでも使いやすい。
  入力listの要素はすべて互換性のある型かつ長さ1である必要がある。
: `ptype = integer()` のように明示的に型指定できる。
: 型が合わなくてvector化できないようならlistでいいから出力して、というときは
  `strict = FALSE`

`purrr::list_c(x, ..., ptype = NULL)`
: listを一段階解消して要素を連結し、vectorにする。
  `purrr::list_simplify()` とは異なり、対応関係や長さは気にせずとにかく連結する。
  `unlist(x, recursive = FALSE) |> as.vector()` を安全にしたようなもの。
: `ptype = integer()` のように明示的に型指定できる。

`purrr::list_rbind(x, ..., names_to = rlang::zap(), ptype = NULL)`
: list of data.frames を `rbind()` して返す。
  例えば、同じ形式のCSVファイルを一気に読んで結合、みたいなときに便利:
  ```r
  files = fs::dir_ls("path/to/data/", glob = "*.csv")
  combined_df = files |> purrr::map(readr::read_csv) |> purrr::list_rbind()
  ```
: `purrr::list_cbind()` もある。

`purrr::list_transpose(x, ..., template = NULL, simplify = NA, ptype = NULL, default = NULL)`
: 行列転置関数`t()`のlist版。
: 例えば、pair of lists <=> list of pairs。
  data.frameをlistとして渡すとlist of rowsが得られる。


### その他

`purrr::reduce(.x, .f, ..., .init)`
: 二変数関数を順々に適用して1つの値を返す。
  C++でいう`std::accumulate()`。
  例えば ``reduce(1:3, `+`)`` の結果は6。

`purrr::accumulate(.x, .f, ..., .init)`
: 二変数関数を順々に適用し、過程も含めてvectorで返す。
  C++でいう`std::partial_sum()`。
  例えば ``accumulate(1:3, `+`)`` の結果は `1 3 6` 。

`purrr::set_names(x, nm = x)`
: 標準の`setNames(x = nm, nm)`は第二引数のほうが省略不可という気持ち悪い定義だが、
  この改良版ではその心配が解消されている。
  長さや型のチェックもしてくれる。


## 関数 `.f` として渡せるもの

apply/map系関数は、名前のついた関数だけでなく、その場で定義された無名関数も受け取れる。
ごく短い関数や一度しか使わない関数には名前をつけないほうが楽ちん。
R 4.1 からはバックスラッシュを使った短縮表記 `\()` が便利。
purrrのmap系関数はチルダ `~` を使ったformulaを受け取って関数として処理してくれる。

```r
# named function
ord = function(x) {strtoi(charToRaw(x), 16L)}
map_int(letters, ord)

# unnamed function
map_int(letters, function(x) {strtoi(charToRaw(x), 16L)})
map_int(letters, \(x) strtoi(charToRaw(x), 16L))

# formula
map_int(letters, ~ strtoi(charToRaw(.x), 16L))

# integer/character
li = list(lower = letters, upper = LETTERS)
map_chr(li, 3L)
map_chr(li, \(x) x[[3L]])
```

formula内部では、第一引数を`.x`または`.`として、第二引数を`.y`として参照する。
`..1`, `..2`, `..3` のような形で三つめ以降も参照できる。

`purrr::as_mapper(.f, ...)`
: `map()` 内部で関数への変換機能を担っている関数。
: formulaを受け取ると `function(.x, .y, . = .x)` のような関数に変換する。
: 数値や文字列を受け取ると `[[` による抽出関数に変換する。
  参照先が存在しない場合の値はmap関数の `.default` 引数で指定できる。

`purrr::partial(...f, ..., .env, .lazy, .first)`
: 引数を部分的に埋めてある関数を作る。C++でいう `std::bind()`


## deprecated/superseded

`purrr::map_dfr(.x, .f, ..., .id = NULL)`, `map_dfc()`
: 入力と出力が一対一対応しないということでmapファミリーから外され、
  `purrr::map() |> purrr::list_rbind()` に取って代わられた。

`purrr::flatten(.x)`
: `purrr::list_flatten()`, `purrr::list_simplify()`, `purrr::list_c()` に取って代わられた。
  `flatten_lgl()`, `flatten_int()`, `flatten_dbl()`, `flatten_chr()`, `flatten_dfr()`
  も同様。

`purrr::invoke(.f, .x = NULL, ..., .env = NULL)`
: `rlang::exec()` に取って代わられた。
: list `.x` の中身を引数として関数 `.f` を呼び出す。
: 関数に渡す引数があらかじめlistにまとまってるときに使う`do.call()`の改良版。

`purrr::invoke_map(.f, .x = list(NULL), ..., .env = NULL)`
: `purrr::map(.f, rlang::exec, ...)` に取って代わられた。
: 関数listを順々に実行してlistで返す。
  引数`.x`は同じ長さのlist of listsか、list of a listをリサイクル。
: e.g., `invoke_map(list(runif, rnorm), list(c(n = 3, 0, 1)))`

`purrr::list_along(x)`
: `rep_along(x, list())` に取って代わられた。
: `x` と同じ長さの空listを作る `vector("list", length(x))` のショートカット。

`purrr::cross2(.x, .y, .filter = NULL)`
: [`tidyr::crossing()` とか `tidyr::expand()`]({{< relref "tidyr.md" >}}) のほうが推奨。
: listの各要素の組み合わせを作る。
  `.filter` に渡した関数が `TRUE` となるものは除外される。
  名前付きlistを渡す `purrr::cross()` や `purrr::cross_df()` もある。

`purrr::transpose(.l)`
: `purrr::list_transpose()` に取って代わられた。



## 関連書籍

<a href="https://www.amazon.co.jp/dp/1492097403?&linkCode=li3&tag=heavywatal-22&linkId=163b4c2d2d4f43d197e985a033d397c1&language=ja_JP&ref_=as_li_ss_il" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1492097403&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=1492097403" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/487311814X?&linkCode=li3&tag=heavywatal-22&linkId=a289b1f9dbb4f189b4209b374662d6f7&language=ja_JP&ref_=as_li_ss_il" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
