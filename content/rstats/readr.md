+++
title = 'readr'
subtitle = "高速で柔軟なテーブル読み込みツール"
[menu.main]
  parent = "rstats"
+++

<https://github.com/hadley/readr>

タブ区切りテキストやCSVファイルを読み込んで `data.frame` にするツール。
標準でも `read.table()` や `read.csv()` があるけど、それらと違って

-   場合により10倍ほど高速
-   文字列を勝手にfactor扱いしたりしないので
    `stringsAsFactors=FALSE` とイチイチ書かなくて済む
-   勝手に列名を変更しない
-   列の名前や型を指定しやすい
-   返ってくるクラスが `c('tbl_df', 'tbl', 'data.frame')`

Rの中から `install.packages('readr')` でインストールし、
使う前に `library(readr)` でパッケージを読み込む。

## 主な関数

`read_csv(file, col_names=TRUE, col_types=NULL, na='NA', skip=0, ...)`
:   カンマ区切りテキスト

`read_tsv(file, col_names=TRUE, col_types=NULL, na='NA', skip=0, ...)`
:   タブ区切りテキスト

`read_table(file, col_names=TRUE, col_types=NULL, na='NA', skip=0, n_max=-1)`
:   連続する空白文字をひとつの区切りと見なして処理

`read_fwf(file, col_positions, col_types=NULL, ...)`
:   fixed width file. 第二引数の指定方法は

    1.  `fwf_empty(infile, skip=0, col_names=NULL)` で自動推定
    2.  `fwf_widths(widths, col_names=NULL)` で幅指定
    3.  `fwf_positions(start, end, col_names=NULL)` で開始・終了位置指定

`read_delim(file, delim, quote='"', ...)`
:   区切り文字など自分で細かく指定

`read_lines(file, n_max=-1L)`
:   1行を1要素とした文字列ベクタを返す

`read_file(file)`
:   ファイルの内容まるごと文字列で返す

## 列の型を指定する

基本的には何も指定しなくても数値や文字列を認識していい感じに設定してくれる。
標準の `read.csv()` などと違って暗黙のfactor変換もしない。
明示的に型を指定したい場合は `col_types` 引数に名前付きリストを渡す。
文字列を使って `'ccdi_'` のようにも指定できる。

`[_] col_skip()`
:   列を読まない

`[i] col_integer()`
:   整数

`[d] col_double()`
:   実数

`[n] col_numeric()`
:   数字以外の文字が含まれていても無視して数字として返す

`[l] col_logical()`
:   TRUE or FALSE

`[c] col_character()`
:   文字列

`col_factor(levels, ordered)`
:   factor

`col_date(format='')`, `col_datetime(format='', tz='UTC')`
:   日付

```r
read_csv('mydata.csv', col_types='ccdi_')

read_csv('mydata.csv', col_types=list(
   length=col_double(),
   count=col_integer(),
   misc=col_skip())
```

## Excelファイルを読み込む

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
