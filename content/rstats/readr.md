+++
title = 'readr'
subtitle = "高速で柔軟なテーブル読み込み"
tags = ["r", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -65
+++

<a href="https://readr.tidyverse.org/">
<img src="/_img/hex-stickers/readr.webp" align="right" width="120" height="139">
</a>

タブ区切りテキストやCSVファイルを読み込んでdata.frameにするツール。
`.gz` や `.xz` などの圧縮ファイルも透過的に読み書き可能。
標準でも `read.table()` や `read.csv()` があるけど、それらと比べて

-   場合により数倍高速・省メモリ
-   列の名前や型を指定しやすい
-   指定した列だけ読み込むこともできる
-   生data.frameではなく安全な [tibble](#tibble) として返してくれる
-   空白行を勝手にスキップする (1.2から `skip_empty_rows = TRUE`)
-   勝手に列名を変更<del>しない</del> する (2.0から `name_repair = "unique"`)
-   <del>`stringsAsFactors = FALSE` とイチイチ書かなくて文字列を読める</del>
    R 4.0 から標準関数もこの挙動。

[tidyverse](https://tidyverse.tidyverse.org/) に含まれているので、
`install.packages("tidyverse")` で一括インストール、
`library(tidyverse)` で一括ロード。
例えば:
```r
library(conflicted)
library(tidyverse)
write_tsv(diamonds, "diamonds.tsv.gz")
read_tsv("diamonds.tsv.gz")
```

最近ちょっと表示がおせっかい過ぎるので `~/.Rprofile` で設定を直す:

```r
options(
  readr.num_columns = 0L,
  readr.show_col_types = FALSE,
  readr.show_progress = FALSE
)
```

-   <https://r4ds.hadley.nz/data-import.html>
-   <https://cran.r-project.org/package=readr>

## 主な関数

### ファイル読み込み

```r
read_delim(file, delim,
  quote = '"',
  escape_backslash = FALSE,
  escape_double = TRUE,
  col_names = TRUE,
  col_types = NULL,
  col_select = NULL,
  id = NULL,
  locale = default_locale(),
  na = c("", "NA"),
  quoted_na = TRUE,
  comment = "",
  trim_ws = FALSE,
  skip = 0,
  n_max = Inf,
  guess_max = min(1000, n_max),
  name_repair = "unique",
  num_threads = readr_threads(),
  progress = show_progress(),
  show_col_types = should_show_types(),
  skip_empty_rows = TRUE,
  lazy = TRUE)
```

**`read_csv(...)`** や **`read_tsv(...)`**
は区切り文字 `delim =` 指定済みのショートカット。

`read_table(...)`
:   連続する空白文字をひとつの区切りと見なして処理

`read_fwf(file, col_positions, ...)`
:   fixed width file. 第二引数の指定方法は

    1.  `fwf_empty(infile, skip = 0, col_names = NULL)` で自動推定
    1.  `fwf_widths(widths, col_names = NULL)` で幅指定
    1.  `fwf_positions(start, end, col_names = NULL)` で開始・終了位置指定

`read_lines(file, skip = 0, n_max = -1L, ...)`, `read_lines_raw(...)`
:   1行を1要素とした文字列ベクタとして読み込む

`read_file(file)`
:   ファイルの内容まるごと文字列で返す


### ファイル書き出し

```r
write_delim(x, path,
  delim = " ",
  na = "NA",
  append = FALSE,
  col_names = !append)
```

**`write_csv(...)`** や **`write_tsv(...)`**
は区切り文字 `delim =` 指定済みのショートカット。

`write_lines(x, path, na = "NA", append = FALSE)`
: vectorやlistを1行ずつ書き出す。

`write_file(x, path, append = FALSE)`
: 文字列をそのまま書き出す。


### 文字列から別の型へ

`parse_number(x, na = c("", "NA"), locale = default_locale())`
:   文字列で最初に登場する数値を抜き出す。
    あるカラムでは末尾に単位が付いちゃってる、みたいな状況でのみ使う。
    それ以外の複雑な判断はしないので、例えば `"6e23"` は単に6になる。

`parse_double(x, ...)`, `parse_integer(x, ...)`
:   文字列をdouble/intの数値として解釈して返す。
    `"6e23"` のような指数形式も大丈夫。
    異物が混じっていた場合は警告してくれる。
    (標準の`as.integer()`とかは黙って小数点以下を切り捨てたりする)

`parse_logical(x, ...)`
:   1/0/T/F/TRUE/FALSE を大文字小文字問わずlogicalに変換。

`parse_factor(x, levels, ordered = FALSE, ...)`

`parse_date(x, format = "", ...)`,
`parse_datetime(x, format = "", ...)`,
`parse_time(x, format = "", ...)`


## 列の型を指定する

https://cran.r-project.org/web/packages/readr/vignettes/column-types.html

基本的には何も指定しなくても数値などを認識していい感じに設定してくれる。
標準の `read.csv()` などと違って暗黙のfactor変換はしない。
整数と実数は区別せずnumeric型で読む(1.2から)。

明示的に型を指定したい場合は `col_types` 引数に `cols()` 関数の値を渡す。
文字列で `"ccdi_"` のように省略することも可能。

```r
read_csv("mydata.csv", col_types="ccdi_")

colsp = cols(length=col_double(), count="i", .default="c")
read_csv("mydata.csv", col_types=colsp)
```

- `[c] col_character()`: 文字列
- `[i] col_integer()`: 整数
- `[d] col_double()`: 実数
- `[l] col_logical()`: TRUE or FALSE
- `[D] col_date(format = "")`: 日付
- `[t] col_time(format = "")`: 時間
- `[T] col_datetime(format = "")`: 日付
- `[n] col_number()`: 数字以外の文字が含まれていても無視して数字として返す
- `[?] col_guess()`: 推測
- `[_] col_skip()`: 列を読まない
- `col_factor(levels, ordered)`: factor

指定した列だけ読むには `cols(..., .default = col_skip())`
とするか `cols_only(...)` を使う。


## Excelファイルを読み込む

<a href="https://readxl.tidyverse.org/">
<img src="/_img/hex-stickers/readxl.webp" align="right" width="120" height="139">
</a>

https://github.com/tidyverse/readxl

自分のデータは絶対にExcel形式ではなくCSVやTSV形式で保存すべきだが、
人から受け取ったファイルや論文のサプリデータがExcelだったら仕方がない。
`readxl` というパッケージを利用すれば、
一旦Officeで開いてCSVに変換するという手間なしで直接Rに読み込める。

Rの中から `install.packages("readxl")` でインストールし、
使う前に `library(readxl)` でパッケージを読み込む。

`excel_sheets(path)`
:   ファイルに含まれるシートの名前を取得

`read_excel(path, sheet = 1, col_names = TRUE, col_types = NULL, na = "", skip = 0)`
:   `.xls` と `xlsx` のどちらの形式でも読める。
    `sheet` は番号でも名前でもいい。
    それ以降の引数については `readr` の関数と同じ。


## 最新版をソースからインストールする

https://github.com/tidyverse/readr

```
remotes::install_github("tidyverse/readr")
# ...
ld: library not found for -lintl
```

見つからないと言われてる `libintl.*` はgettextの一部なのでそいつを入れる。
パス指定で楽をするため、keg-onlyだけど無理矢理シムリンクを張る
(ホントは良くないかも)。

```sh
brew install gettext
brew link gettext --force
```

Homebrewを`/usr/local/`以外に入れている場合は、
それをRに見つけさせるため `~/.R/Makevars` にオプションを書く。

```makefile
LDFLAGS = -L${HOME}/.homebrew/lib
```

再びRで `install_github("tidyverse/readr")` を試みる。




## tibble

<a href="https://tibble.tidyverse.org/">
<img src="/_img/hex-stickers/tibble.webp" align="right" width="120" height="139">
</a>

- https://r4ds.had.co.nz/tibbles.html
- https://github.com/tidyverse/tibble

`tbl_df` クラスが付与された改良版data.frameのことを**tibble**と呼ぶ。
[readr]({{< relref "readr.md" >}}) で読み込んだデータもこの形式になる。

```r
tbl_mtcars = as_tibble(mtcars)
class(tbl_mtcars)
## [1] "tbl_df"     "tbl"        "data.frame"
class(mtcars)
## [1] "data.frame"
```

生のdata.frameとの違いは:

-   巨大なデータをうっかり`print()`しても画面を埋め尽くさない。
    (逆に全体を見たい場合は工夫が必要。後述)
-   列名の部分一致で良しとしない。
    例えば `mtcars$m` は黙ってvectorを返してしまうが、
    `tbl_mtcars$m` は警告つき `NULL` 。
-   型に一貫性があり、勝手に`drop = TRUE`しない。
    例えば `mtcars[,"mpg"]` はvectorになってしまうが、
    `tbl_mtcars[,"mpg"]` はtibbleのまま。
    vectorが欲しい場合は二重四角括弧 `tbl_mtcars[["mpg"]]`。
-   行の名前は使わない。
    行の名前として保持されている情報を使いたい場合は
    `rownames_to_column()` とか `rowid_to_column()`
    で独立した列にしておく必要がある。

新しいtibble 1.4以降では
[pillar](https://github.com/r-lib/pillar/)
というパッケージが有効数字や欠損値などの表示形式を勝手にイジるようになってしまった。
見やすくない上に遅いので私は `registerS3method()` で上書きしている。


### 関数

[`tibble::tibble(...)`](https://tibble.tidyverse.org/reference/tibble.html)
:   tibbleを新規作成。ちょっと昔までは `dplyr::data_frame()` だった。
:   `base::data.frame()` と違ってバグが混入しにくくて便利:
    -   勝手に型変換しない (`stringsAsFactors = FALSE`が基本)
    -   勝手に列名を変えない
    -   長さ1の変数以外はリサイクルしない
    -   引数の評価がlazyに行われるので前の列を利用して後の列を作ったりできる
    -   `tbl_df` クラスを付加
    -   ただし1.4以降のバージョンでは表示が遅くて見にくい

[`tibble::as_tibble(x)`](https://tibble.tidyverse.org/reference/as_tibble.html)
:   既存のdata.frameやmatrixをtibbleに変換。
    ちょっと昔までは `dplyr::tbl_df()` とか `dplyr::as_data_frame()` だった。
    v2.0からは列名がちゃんとついてないとエラーになる。

[`tibble::new_tibble(x, ..., nrow, class = NULL)`](https://tibble.tidyverse.org/reference/new_tibble.html)
:   tibbleのサブクラスを扱う開発者向け `as_tibble()` 。
    検証なし、`nrow` 必須の代わりに高速。
    クラスを先頭に追加できるのも便利。

[`tibble::enframe(x, name = "name", value = "value")`](https://tibble.tidyverse.org/reference/enframe.html)
:   名前付きvectorとかlistを2列のtibbleに変換する。
    `tibble::deframe(x)` はその逆。
    `c(a = 1, b = 2) |> enframe() |> deframe()`

[`tibble::add_row(.data, ..., .before = NULL, .after = NULL)`](https://tibble.tidyverse.org/reference/add_row.html)
:   既存のtibbleに新しいデータを1行追加する。

[`tibble::rownames_to_column(df, var = "rowname")`](https://tibble.tidyverse.org/reference/rownames.html)
:   行の名前をcharacter型で1列目の変数にする。`dplyr::add_rownames()`の後継。
:   `tibble::rowid_to_column(df, var = "rowid")` はそれを整数で。
:   `tibble::column_to_rownames(df, var = "rowname")` はその逆。
:   `tibble::remove_rownames(df)` は消すだけ。

`tibble::glimpse(.data, width = NULL)`
:   データの中身をざっと見る。
    `print()` とか `str()` のようなもの。

`pillar::type_sum(x)`
:   オブジェクトの型

`pillar::obj_sum(x)`
:   `type_sum`とサイズ e.g., `"data.frame [150 x 5]"`


### 設定

表示される行数や幅を調節する項目には以下のようなものがある。
`~/.Rprofile` に書いておけば起動時に勝手に設定される。

```r
height = 30L  # for example
width = 160L
options(
  pillar.neg = FALSE,
  pillar.subtle = FALSE,
  pillar.print_max = height,
  pillar.print_min = height,
  pillar.width = width,
  width = min(width, 10000L)
)
```

- https://github.com/r-lib/pillar/blob/main/R/options.R
- https://github.com/tidyverse/tibble/blob/main/R/options.R


### 大きいtibbleの全体を表示する

普通に `print()` すると大きいtibbleの全体が見えない。
`print()` 関数にオプションを指定したり `utils::page()` を利用したりする必要がある。
RStudioやVSCodeを使っている場合は `View()` でスプレッドシートのように閲覧可能。

```r
# pillar:::print.tbl()にオプションを渡す方式
diamonds |> print(n = Inf, width = Inf)
diamonds |> page("print", n = Inf, width = Inf)

# 標準data.frame用のprint()を呼ぶ方式
diamonds |> base::print.data.frame(max = .Machine$integer.max)
```

オプションをいちいち設定しなくて済むように
[`max_print()`](https://github.com/heavywatal/rwtl/blob/master/R/print.R),
[`less()`](https://github.com/heavywatal/rwtl/blob/master/R/pipe.R)
のような関数を定義しておくのもよい。

`getOption("max.print")` の初期値は99999だがRStudioでは勝手に1000まで下げられる。


## 関連書籍

<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?s=books&ie=UTF8&qid=1508340700&sr=1-3&linkCode=li3&tag=heavywatal-22&linkId=6a53371cc80e2d6d7fc50fae5d8b862d" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%A7%E3%81%AF%E3%81%98%E3%82%81%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9-Hadley-Wickham/dp/487311814X/ref=as_li_ss_il?ie=UTF8&qid=1508340144&sr=8-1&keywords=r&linkCode=li3&tag=heavywatal-22&linkId=4137d3d3351f8ccab5a93cefdc28fdec" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
