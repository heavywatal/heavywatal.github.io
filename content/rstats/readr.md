+++
title = 'readr'
subtitle = "高速で柔軟なテーブル読み込みツール"
tags = ["r", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -65
+++

-   http://r4ds.had.co.nz/data-import.html
-   https://cran.r-project.org/package=readr

タブ区切りテキストやCSVファイルを読み込んでdata.frameにするツール。
圧縮された `***.tsv.gz` なども自動的に展開して読んでくれる。
標準でも `read.table()` や `read.csv()` があるけど、それらと違って

-   場合により10倍ほど高速
-   文字列を勝手にfactor扱いしたりしないので
    `stringsAsFactors=FALSE` とイチイチ書かなくて済む
-   勝手に列名を変更しない
-   列の名前や型を指定しやすい
-   指定した列だけ読み込むこともできる
-   生data.frameではなく安全な [tibble](#tibble) として返してくれる

[tidyverse](https://github.com/tidyverse/tidyverse) に含まれているので、
`install.packages('tidyverse')` で一括インストール、
`library(tidyverse)` で一括ロード。


## 主な関数

### ファイル読み込み

`read_csv(file, col_names=TRUE, col_types=NULL, locale=default_locale(), na=c('', 'NA'), comment='', trim_ws=TRUE, skip=0, n_max=-1, progress=interactive())`
:   カンマ区切りテキスト

`read_tsv(...)`
:   タブ区切りテキスト

`read_table(...)`
:   連続する空白文字をひとつの区切りと見なして処理

`read_fwf(file, col_positions, ...)`
:   fixed width file. 第二引数の指定方法は

    1.  `fwf_empty(infile, skip=0, col_names=NULL)` で自動推定
    2.  `fwf_widths(widths, col_names=NULL)` で幅指定
    3.  `fwf_positions(start, end, col_names=NULL)` で開始・終了位置指定

`read_delim(file, delim, quote='"', escape_backslash=FALSE, ...)`
:   区切り文字など自分で細かく指定

`read_lines(file, skip=0, n_max=-1L, ...)`, `read_lines_raw(...)`
:   1行を1要素とした文字列ベクタを返す

`read_file(file)`
:   ファイルの内容まるごと文字列で返す

書き出し用の関数 `write_***()` も一応付いているが、
まだ圧縮ファイルを書き出せないので微妙。


### 文字列から別の型へ

`parse_number(x, na=c('', 'NA'), locale=default_locale())`
:   文字列で最初に登場する数値を抜き出す。
    あるカラムでは末尾に単位が付いちゃってる、みたいな状況でのみ使う。
    それ以外の複雑な判断はしないので、例えば `'6e23'` は単に6になる。

`parse_double(x, ...)`, `parse_int(x, ...)`
:   文字列をdouble/intの数値として解釈して返す。
    `'6e23'` のような指数形式も大丈夫。
    異物が混じっていた場合は警告してくれる。

`parse_logical(x, ...)`
:   1/0/T/F/TRUE/FALSE を大文字小文字問わずlogicalに変換。

`parse_factor(x, levels, ordered=FALSE, ...)`

`parse_date(x, format='', ...)`, `parse_datetime(x, format='', ...)`, `parse_time(x, format='', ...)`


## 列の型を指定する

https://cran.r-project.org/web/packages/readr/vignettes/column-types.html

基本的には何も指定しなくても数値や文字列を認識していい感じに設定してくれる。
標準の `read.csv()` などと違って暗黙のfactor変換はしない。
明示的に型を指定したい場合は `col_types` 引数に `cols()` 関数の値を渡す。
文字列で `'ccdi_'` のように省略することも可能。

```r
read_csv('mydata.csv', col_types='ccdi_')

colsp = cols(length=col_double(), count='i', .default='c')
read_csv('mydata.csv', col_types=colsp)
```

- `[c] col_character()`: 文字列
- `[i] col_integer()`: 整数
- `[d] col_double()`: 実数
- `[l] col_logical()`: TRUE or FALSE
- `[D] col_date(format='')`: 日付
- `[t] col_time(format='')`: 時間
- `[T] col_datetime(format='')`: 日付
- `[n] col_number()`: 数字以外の文字が含まれていても無視して数字として返す
- `[?] col_guess()`: 推測
- `[_] col_skip()`: 列を読まない
- `col_factor(levels, ordered)`: factor

指定した列だけ読むには `cols(..., .default=col_skip())`
とするか `cols_only(...)` を使う。


## Excelファイルを読み込む

https://github.com/hadley/readxl

自分のデータは絶対にExcel形式ではなくCSVやTSV形式で保存すべきだが、
人から受け取ったファイルや論文のサプリデータがExcelだったら仕方がない。
`readxl` というパッケージを利用すれば、
一旦Officeで開いてCSVに変換するという手間なしで直接Rに読み込める。

Rの中から `install.packages('readxl')` でインストールし、
使う前に `library(readxl)` でパッケージを読み込む。

`excel_sheets(path)`
:   ファイルに含まれるシートの名前を取得

`read_excel(path, sheet=1, col_names=TRUE, col_types=NULL, na='', skip=0)`
:   `.xls` と `xlsx` のどちらの形式でも読める。
    `sheet` は番号でも名前でもいい。
    それ以降の引数については `readr` の関数と同じ。


## 最新版をソースからインストールする

https://github.com/tidyverse/readr

```r
devtools::install_github('tidyverse/readr')
# ...
ld: library not found for -lintl
```

見つからないと言われてる `libintl.*` はgettextの一部なのでそいつを入れる。
パス指定で楽をするため、keg-onlyだけど無理矢理シムリンクを張る
(ホントは良くないかも)。

```sh
% brew install gettext
% brew link gettext --force
```

Homebrewを`/usr/local/`以外に入れている場合は、
それをRに見つけさせるため `~/.R/Makevars` にオプションを書く。

```
LDFLAGS = -L${HOME}/.homebrew/lib
```

再びRで `install_github('tidyverse/readr')` を試みる。




## tibble

http://r4ds.had.co.nz/tibbles.html

`tbl_df` クラスが付与された改良版data.frameのことを**tibble**と呼ぶ。
もともとは [dplyr]({{< relref "dplyr.md" >}}) パッケージで扱っていたが、独立パッケージになった。
[readr]({{< relref "readr.md" >}}) で読み込んだデータもこの形式になる。

```r
> tbl_iris = as_tibble(iris)
> class(tbl_iris)
[1] "tbl_df"     "tbl"        "data.frame"
> class(iris)
[1] "data.frame"
```

生のdata.frameとの違いは:

-   巨大なデータをうっかり`print()`しても画面を埋め尽くさない。
    (逆に全体を見たい場合は工夫が必要。後述)
-   列名の部分一致で良しとしない。
    例えば `iris$Spec` は黙ってvectorを返してしまうが、
    `tbl_iris$Spec` は警告つき `NULL` 。
-   型に一貫性があり、勝手に`drop=TRUE`しない。
    例えば `iris[,'Species']` はvectorになってしまうが、
    `tbl_iris[,'Species']` はtibbleのまま。

### 関数

`tibble::tibble(...)`
:   tibbleを新規作成。ちょっと昔までは `dplyr::data_frame()` だった。
:   `base::data.frame()` と違ってバグが混入しにくくて便利:
    -   勝手に型変換しない (`stringsAsFactors=FALSE`が基本)
    -   勝手に列名を変えない
    -   長さ1の変数以外はリサイクルしない
    -   引数の評価がlazyに行われるので前の列を利用して後の列を作ったりできる
    -   `tbl_df` クラスを付加

`tibble::as_tibble(x)`
:   既存のdata.frameやmatrixをtibbleに変換。
    ちょっと昔までは `dplyr::tbl_df()` とか `dplyr::as_data_frame()` だった。

`tibble::add_row(.data, ...)`
:   既存のtibbleに新しいデータを1行追加する。

`tibble::rownames_to_column(df, var='rowname')`
:   行の名前(無ければ1からの整数)を1列目の変数する`dplyr::add_rownames()`の改良版。
:   `tibble::column_to_rownames(df, var='rowname')` はその逆。
:   `tibble::remove_rownames(df)` は消すだけ。

`tibble::glimpse(.data, width=NULL)`
:   データの中身をざっと見る。
    `print()` とか `str()` のようなもの。

`tibble::type_sum(x)`
:   オブジェクトの型

`tibble::obj_sum(x)`
:   `type_sum`とサイズ e.g., `"data.frame [150 x 5]"`


### 設定

`.Rprofile` に以下のように書くことで、
大きいtibbleの`print()`で表示される行数を調節することができる。
デフォルトはどちらも`10L`。

```r
# Maximum number of rows to print(tbl_df)
options(tibble.print_max=30L)

# Number of rows to print(tbl_df) if exceeded the maximum
options(tibble.print_min=30L)
```

### PAGERで全体を表示する

巨大テーブルをコンパクトに表示してくれるのは助かるけど、逆に全体を見渡したい場合もある。
一時的に最大表示行数を拡大してPAGERに流し込むような関数を書いておくと良い。
普通のdata.frameを見るのにも便利。

```r
less = function(x, method=c('print', 'dput'), max.print=getOption('max.print'), ...) {
    opts = options(max.print=max.print,
            tibble.print_max=max.print,
            tibble.print_min=max.print)
    on.exit(options(opts))
    utils::page(x, match.arg(method), ...)
}

less(iris)
less(as_tibble(iris))
```


