+++
title = 'tidyr'
subtitle = "シンプルなデータ変形ツール"
tags = ["r", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -74
+++

<a href="https://tidyr.tidyverse.org/">
<img src="https://raw.githubusercontent.com/rstudio/hex-stickers/master/SVG/tidyr.svg" align="right" width="120" height="139">
</a>

data.frameを縦長・横広・入れ子に変形・整形するためのツール。
[dplyr]({{< relref "dplyr.md" >}}) や [purrr]({{< relref "purrr.md" >}})
と一緒に使うとよい。
[reshape2]({{< relref "reshape2.md" >}}) を置き換えるべく再設計された改良版。

[tidyverse](https://tidyverse.tidyverse.org/) に含まれているので、
`install.packages("tidyverse")` で一括インストール、
`library(tidyverse)` で一括ロード。

-   <https://r4ds.had.co.nz/tidy-data.html>
-   <https://github.com/tidyverse/tidyr>
-   `vignette("tidy-data")`
-   `demo(package = "tidyr")`
-   https://speakerdeck.com/yutannihilation/tidyr-pivot

パイプ演算子 `%>%` については[dplyr]({{< relref "dplyr.md" >}})を参照。


## Pivoting: 縦長 ↔ 横広

https://tidyr.tidyverse.org/articles/pivot.html

### `tidyr::pivot_longer()` で縦長にする

複数列にまたがっていた値を1列にまとめ、元の列名をその横に添えることで、
data.frameを横広(wide-format)から縦長(long-format)に変形する。
`reshape2::melt()`, `tidyr::gather()` の改良版。

`tidyr::pivot_longer(data, cols, names_to = "name", ..., values_to = "value", ...)`

`cols`
: 動かしたい値が含まれている列。
  コロンで範囲指定、文字列、tidyselect関数なども使える。
  動かさない列を `!` で反転指定するのほうが楽なことも多い。

`names_to`
: 元々列名だったものを入れる列の名前

`values_to`
: 値の移動先の列名

```r
anscombe %>% tibble::rowid_to_column("id")
#>    id x1 x2 x3 x4    y1   y2    y3    y4
#> 1   1 10 10 10  8  8.04 9.14  7.46  6.58
#> 2   2  8  8  8  8  6.95 8.14  6.77  5.76
#> 3   3 13 13 13  8  7.58 8.74 12.74  7.71
#> 4   4  9  9  9  8  8.81 8.77  7.11  8.84
#> 5   5 11 11 11  8  8.33 9.26  7.81  8.47
#> 6   6 14 14 14  8  9.96 8.10  8.84  7.04
#> 7   7  6  6  6  8  7.24 6.13  6.08  5.25
#> 8   8  4  4  4 19  4.26 3.10  5.39 12.50
#> 9   9 12 12 12  8 10.84 9.13  8.15  5.56
#> 10 10  7  7  7  8  4.82 7.26  6.42  7.91
#> 11 11  5  5  5  8  5.68 4.74  5.73  6.89

anscombe_long = anscombe %>% tibble::rowid_to_column("id") %>%
  pivot_longer(!id, names_to = "namae", values_to = "atai") %>%
  print()
#> # tbl_df [88 x 3]
#>       id namae  atai
#>    <int> <chr> <dbl>
#>  1     1    x1 10.00
#>  2     1    x2 10.00
#>  3     1    x3 10.00
#>  4     1    x4  8.00
#> --
#> 85    11    y1  5.68
#> 86    11    y2  4.74
#> 87    11    y3  5.73
#> 88    11    y4  6.89

# anscombe %>% gather("namae", "atai")
```

### `tidyr::pivot_wider()` で横広にする

1列にまとまっていた値を、別の変数に応じて複数の列に並べ直すことで、
data.frameを縦長(long-format)から横広(wide-format)に変形する。
`reshape2::dcast()`, `tidyr::spread()` の改良版。

`tidyr::pivot_wider(data, id_cols = NULL, names_from = name, ..., values_from = value, values_fill = NULL, values_fn = NULL)`

`id_cols`
: ここで指定した列のユニークな組み合わせが変形後にそれぞれ1行になる。
  `!`で反転指定、`:`で範囲指定、文字列、tidyselect関数なども使える。
  デフォルトでは `names_from` と `values_from` で指定されなかった列すべて。

`names_from`
: 新しく列名になる列。"name" という列名なら省略可能。

`values_from`
: 動かしたい値が入っている列。"value" という列名なら省略可能。

`values_fill`
: 存在しない組み合わせのセルを埋める値。
  列によって値を変えたい場合は名前付きリストで渡す。

`values_fn`
: `id_cols` の組み合わせが一意に定まらず複数のvalueを1セルに詰め込む場合の処理関数。
  デフォルトでは警告とともに `list()` が使われる。


```r
anscombe_long %>%
  pivot_wider(names_from = namae, values_from = atai) %>%
  dplyr::select(!id)
#> # tbl_df [11 x 8]
#>       x1    x2    x3    x4    y1    y2    y3    y4
#>    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1    10    10    10     8  8.04  9.14  7.46  6.58
#>  2     8     8     8     8  6.95  8.14  6.77  5.76
#>  3    13    13    13     8  7.58  8.74 12.74  7.71
#>  4     9     9     9     8  8.81  8.77  7.11  8.84
#> --
#>  8     4     4     4    19  4.26  3.10  5.39 12.50
#>  9    12    12    12     8 10.84  9.13  8.15  5.56
#> 10     7     7     7     8  4.82  7.26  6.42  7.91
#> 11     5     5     5     8  5.68  4.74  5.73  6.89

# anscombe_long %>% spread(namae, atai) %>% dplyr::select(!id)
```

カテゴリカル変数を指示変数(ダミー変数)に変換するのにも使える:

```r
pg = PlantGrowth %>% dplyr::slice(c(1, 2, 11, 12, 21, 22)) %>% print()
#   weight group
# 1   4.17  ctrl
# 2   5.58  ctrl
# 3   4.81  trt1
# 4   4.17  trt1
# 5   6.31  trt2
# 6   5.12  trt2
pg %>% tibble::rowid_to_column("id") %>%
  dplyr::mutate(name = group, value = 1L) %>%
  tidyr::pivot_wider(values_fill = 0L) %>%
  dplyr::select(!c(id, ctrl))
#   weight group  trt1  trt2
#    <dbl> <fct> <int> <int>
# 1   4.17  ctrl     0     0
# 2   5.58  ctrl     0     0
# 3   4.81  trt1     1     0
# 4   4.17  trt1     1     0
# 5   6.31  trt2     0     1
# 6   5.12  trt2     0     1
```


### `tidyr::pivot_*` 関数のもっと高度なオプション

`names_sep` や `names_pattern` を指定して
`names_to`, `name_from` に複数の値を渡すと
`tidyr::separate()` / `tidyr::unite()` 的な操作も同時にやってしまえる:

```r
anscombe %>% tibble::rowid_to_column("id") %>%
  tidyr::pivot_longer(!id, names_to = c("axis", "group"), names_sep = 1L) %>%
  print() %>%
  pivot_wider(id, names_from = c(axis, group), names_sep = "_")
#> # tbl_df [88 x 4]
#>       id  axis group value
#>    <int> <chr> <chr> <dbl>
#>  1     1     x     1 10.00
#>  2     1     x     2 10.00
#>  3     1     x     3 10.00
#>  4     1     x     4  8.00
#> --
#> 85    11     y     1  5.68
#> 86    11     y     2  4.74
#> 87    11     y     3  5.73
#> 88    11     y     4  6.89
#> # tbl_df [11 x 9]
#>       id   x_1   x_2   x_3   x_4   y_1   y_2   y_3   y_4
#>    <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1     1    10    10    10     8  8.04  9.14  7.46  6.58
#>  2     2     8     8     8     8  6.95  8.14  6.77  5.76
#>  3     3    13    13    13     8  7.58  8.74 12.74  7.71
#>  4     4     9     9     9     8  8.81  8.77  7.11  8.84
#> --
#>  8     8     4     4     4    19  4.26  3.10  5.39 12.50
#>  9     9    12    12    12     8 10.84  9.13  8.15  5.56
#> 10    10     7     7     7     8  4.82  7.26  6.42  7.91
#> 11    11     5     5     5     8  5.68  4.74  5.73  6.89
```

```r
VADeaths
#>       Rural Male Rural Female Urban Male Urban Female
#> 50-54       11.7          8.7       15.4          8.4
#> 55-59       18.1         11.7       24.3         13.6
#> 60-64       26.9         20.3       37.0         19.3
#> 65-69       41.0         30.9       54.6         35.1
#> 70-74       66.0         54.3       71.1         50.0

VADeaths %>%
  as.data.frame() %>%
  rownames_to_column("age") %>%
  pivot_longer(!age, names_to = c("region", "sex"), names_sep = " ", values_to = "death")
#> # tbl_df [20 x 4]
#>      age region    sex death
#>    <chr>  <chr>  <chr> <dbl>
#>  1 50-54  Rural   Male  11.7
#>  2 50-54  Rural Female   8.7
#>  3 50-54  Urban   Male  15.4
#>  4 50-54  Urban Female   8.4
#> --
#> 17 70-74  Rural   Male  66.0
#> 18 70-74  Rural Female  54.3
#> 19 70-74  Urban   Male  71.1
#> 20 70-74  Urban Female  50.0
```

`names_transform` に関数を指定すると、
列名だったものに適用される。
例えば型変換に使える:

```r
anscombe %>%
  tibble::rowid_to_column("id") %>%
  tidyr::pivot_longer(!id,
    names_to = c("axis", "group"),
    names_sep = 1L,
    names_transform = list(group = as.integer)) %>%
  tidyr::pivot_wider(c(id, group), names_from = axis) %>%
  dplyr::select(!id) %>%
  dplyr::arrange(group)
#> # tbl_df [44 x 3]
#>    group     x     y
#>    <int> <dbl> <dbl>
#>  1     1    10  8.04
#>  2     1     8  6.95
#>  3     1    13  7.58
#>  4     1     9  8.81
#> --
#> 41     4    19 12.50
#> 42     4     8  5.56
#> 43     4     8  7.91
#> 44     4     8  6.89
```

`names_prefix` を使えば、列名の頭に共通して付いてた文字を消せる:

```r
anscombe %>%
  dplyr::select(starts_with("x")) %>%
  tidyr::pivot_longer(everything(), names_prefix = "x")
#> # tbl_df [44 x 2]
#>     name value
#>    <chr> <dbl>
#>  1     1    10
#>  2     2    10
#>  3     3    10
#>  4     4     8
#> --
#> 41     1     5
#> 42     2     5
#> 43     3     5
#> 44     4     8
```

`names_to` に `".value"` という特殊な値を渡すことで、
旧列名から新しい列名が作られ、複数列への縦長変形を同時にできる。
また `names_to = c("name", NA)` のように不要な列を捨てることもできる。

```r
tidy_anscombe = anscombe %>%
  tidyr::pivot_longer(                # 縦長に変形したい
    everything(),                     # すべての列について
    names_to = c(".value", "group"),  # x, yを列名に、1, 2, 3をgroup列に
    names_sep = 1L,                   # 切る位置
    names_transform = list(group = as.integer)) %>%   # 型変換
  dplyr::arrange(group) %>%           # グループごとに並べる
  print()                             # ggplotしたい形！
#> # tbl_df [44 x 3]
#>    group     x     y
#>    <int> <dbl> <dbl>
#>  1     1    10  8.04
#>  2     1     8  6.95
#>  3     1    13  7.58
#>  4     1     9  8.81
#> --
#> 41     4    19 12.50
#> 42     4     8  5.56
#> 43     4     8  7.91
#> 44     4     8  6.89
```

See https://speakerdeck.com/yutannihilation/tidyr-pivot?slide=67 for details.


## Nested data.frame --- 入れ子構造

https://tidyr.tidyverse.org/articles/nest.html

### `tidyr::nest(data, ..., .names_sep = NULL)`

data.frameをネストして(入れ子にして)、list of data.frames のカラムを作る。
内側のdata.frameに押し込むカラムを `...` に指定するか、
外側に残すカラムを `!` で反転指定する。

```r
diamonds %>% nest(NEW_COLUMN = !cut) %>% unnest()
#> # tbl_df [5 x 2]
#>         cut           NEW_COLUMN
#>       <ord>               <list>
#> 1     Ideal <tbl_df [21551 x 9]>
#> 2   Premium <tbl_df [13791 x 9]>
#> 3      Good  <tbl_df [4906 x 9]>
#> 4 Very Good <tbl_df [12082 x 9]>
#> 5      Fair  <tbl_df [1610 x 9]>

# equivalent to
diamonds %>% nest(NEW_COLUMN = c(carat, color:z))
diamonds %>% dplyr::group_nest(cut, .key = "NEW_COLUMN")
diamonds %>% dplyr::nest_by(cut, .key = "NEW_COLUMN")
```

なんでもかんでもフラットなdata.frameにして
[dplyr]({{< relref "dplyr.md" >}})を駆使する時代は終わり、
ネストしておいて[purrr]({{< relref "purrr.md" >}})を適用するのが
tidyverse時代のクールなやり方らしい。

cf. [Hadley Wickham: Managing many models with R (YouTube)](https://www.youtube.com/watch?v=rz3_FDVt9eg)


### `tidyr::unnest(data, cols, ...)`

ネストされたdata.frameを展開してフラットにする。
list of data.framesだけでなく、list of vectorsとかでもよい。

ネストされた列が複数ある場合に曖昧なコードにならないよう
`cols` を明示的に指定することが求められる。

```r
diamonds %>%
  nest(NEW_COLUMN = !cut) %>%
  unnest(NEW_COLUMN)
```

## その他の便利関数

### `tidyr::separate()`

文字列カラムを任意のセパレータで複数カラムに分割。
`tidyr::unite()` の逆。
`reshape2::colsplit()` に相当。

`tidyr::separate(data, col, into, sep = "[^[:alnum:]]", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn", ...)`

`col`
:   切り分けたい列の名前

`into`
:   切り分けたあとの新しい列名を文字列ベクタで

`sep = "[^[:alnum:]]"`
:   セパレータを正規表現で。デフォルトはあらゆる非アルファベット。
:   整数を渡すと位置で切れる。例えば `A4` を `1L` で切ると `A` と `4` に。

`remove = TRUE`
:   切り分ける前の列を取り除くかどうか

`convert = FALSE`
:   切り分け後の値の型変換を試みるか

`extra = "warn"`
:   列数が揃わないときにどうするか: `warn`, `drop`, `merge`

`fill = "warn"`
:   足りない場合にどっち側をNAで埋めるか: `warn`, `right`, `left`。
    つまり、文字を左詰めにするには`right`が正解(紛らわしい)。

```r
va_deaths = VADeaths %>% as.data.frame() %>% tibble::rownames_to_column("class") %>% print()
#>   class Rural Male Rural Female Urban Male Urban Female
#> 1 50-54       11.7          8.7       15.4          8.4
#> 2 55-59       18.1         11.7       24.3         13.6
#> 3 60-64       26.9         20.3       37.0         19.3
#> 4 65-69       41.0         30.9       54.6         35.1
#> 5 70-74       66.0         54.3       71.1         50.0

va_deaths %>%
  tidyr::separate(class, c("lbound", "ubound"), "-", convert = TRUE)
#>   lbound ubound Rural Male Rural Female Urban Male Urban Female
#> 1     50     54       11.7          8.7       15.4          8.4
#> 2     55     59       18.1         11.7       24.3         13.6
#> 3     60     64       26.9         20.3       37.0         19.3
#> 4     65     69       41.0         30.9       54.6         35.1
#> 5     70     74       66.0         54.3       71.1         50.0
```

行方向に分割する `tidyr::separate_rows(data, ..., sep, convert)` もある。

`tidyr::extract(data, col, into, regex, ...)`
を使えば正規表現でもっと細かく指定できる。

名前の似てる `tidyr::extract_numeric(x)` は
文字列から数字部分をnumericとして抜き出す関数だったが今はdeprecatedなので、
新しい[`readr::parse_number()`]({{< relref "readr.md" >}})を使うべし。


### `tidyr::unite(data, col, ..., sep = "_", remove = TRUE, na.rm = FALSE)`

複数カラムを結合して1列にする。
`tidyr::separate()` の逆。

`paste()` とか `stringr::str_c()` でも似たようなことができるけど
`na.rm = TRUE` の挙動が欲しいときに便利。

```r
df = tibble(x = c("x", "x", NA), y = c("y", NA, "y"))

df %>% tidyr::unite(z, c(x, y), sep = "_", remove = FALSE)
#> # tbl_df [3 x 3]
#>       z     x     y
#>   <chr> <chr> <chr>
#> 1   x_y     x     y
#> 2  x_NA     x  <NA>
#> 3  NA_y  <NA>     y
df %>% tidyr::unite(z, c(x, y), sep = "_", remove = FALSE, na.rm = TRUE)
#> # tbl_df [3 x 3]
#>       z     x     y
#>   <chr> <chr> <chr>
#> 1   x_y     x     y
#> 2     x     x  <NA>
#> 3     y  <NA>     y
df %>% dplyr::mutate(z = stringr::str_c(x, y, sep = "_"))
#> # tbl_df [3 x 3]
#>       x     y     z
#>   <chr> <chr> <chr>
#> 1     x     y   x_y
#> 2     x  <NA>  <NA>
#> 3  <NA>     y  <NA>
df %>% dplyr::mutate(z = dplyr::coalesce(x, y))
#> # tbl_df [3 x 3]
#>       x     y     z
#>   <chr> <chr> <chr>
#> 1     x     y     x
#> 2     x  <NA>     x
#> 3  <NA>     y     y
```


### `tidyr::complete(data, ..., fill = list())`

指定した列の全ての組み合わせが登場するように、
指定しなかった列に欠損値`NA`(あるいは任意の値)を補完した行を挿入する。

```r
df %>% complete(key1, key2, fill = list(val1 = 0, val2 = "-"))
```

### `tidyr::expand(data, ...)`

指定した列の全ての組み合わせが登場するような新しいdata.frameを作る。
全ての列を指定すれば`complete()`と同じ効果だが、
指定しなかった列が消えるという点では異なる。

`crossing(...)`はvectorを引数に取る亜種で、
tibble版`expand.grid(...)`のようなもの。

`nesting(...)`は存在するユニークな組み合わせのみ残す、
`nest(data, ...) %>% dplyr::select(!data)`のショートカット。
この結果は`expand()`や`complete()`の引数としても使える。

数値vectorの補完には`full_seq(x, period, tol = 1e-6)`が便利。


### `tidyr::drop_na(data, ...)`

`complete()`の逆。
指定した列に`NA`が含まれてる行を削除する。
何も指定しなければ標準の `data[complete.cases(data),]` と同じ。

### `tidyr::replace_na()`

欠損値 `NA` を好きな値で置き換える。
これまでは `mutate(x = ifelse(is.na(x), 0, x))` のようにしてたところを

```r
df %>% replace_na(list(x = 0, y = "unknown"))
```

逆に、特定の値を`NA`にしたい場合は
[`dplyr::na_if()`]({{< relref "dplyr.md" >}})


### `tidyr::fill()`

`NA` を、その列の直前の `NA` でない値で埋める。
えくせるでセルの結合とかやってしまって、
最初のセルにしか値が無いような場合に使うのかな？


## 関連書籍

<a href="https://www.amazon.co.jp/dp/4297121700?&linkCode=li3&tag=heavywatal-22&linkId=77762a4d0080a840ec5d94df9c0c5ceb&language=ja_JP&ref_=as_li_ss_il" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4297121700&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4297121700" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?s=books&ie=UTF8&qid=1508340700&sr=1-3&linkCode=li3&tag=heavywatal-22&linkId=6a53371cc80e2d6d7fc50fae5d8b862d" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%A7%E3%81%AF%E3%81%98%E3%82%81%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9-Hadley-Wickham/dp/487311814X/ref=as_li_ss_il?ie=UTF8&qid=1508340144&sr=8-1&keywords=r&linkCode=li3&tag=heavywatal-22&linkId=4137d3d3351f8ccab5a93cefdc28fdec" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
