+++
title = 'tidyr'
subtitle = "シンプルなデータ変形ツール"
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -75
+++

-   <https://github.com/hadley/tidyr>
-   <http://rpubs.com/m_dev/tidyr-intro-and-demos>
-   <http://blog.rstudio.org/2014/07/22/introducing-tidyr/>
-   `vignette("tidy-data")`
-   `demo(package="tidyr")`

[dplyr]({{< relref "dplyr.md" >}}) と共に `data.frame` を整形するためのツール。
[reshape2]({{< relref "reshape2.md" >}}) を置き換えるべく再設計された改良版。

Rの中で `install.packages('tidyr')` としてインストールし、
使う前に `library(tidyr)` で読み込む。

下記のコード例で使うデータ

```r
iris %>% head(3) %>% add_rownames('id')

Source: local data frame [3 x 6]

     id Sepal.Length Sepal.Width Petal.Length Petal.Width Species
  (chr)        (dbl)       (dbl)        (dbl)       (dbl)  (fctr)
1     1          5.1         3.5          1.4         0.2  setosa
2     2          4.9         3.0          1.4         0.2  setosa
3     3          4.7         3.2          1.3         0.2  setosa
```

## `tidyr::gather()`

複数列にまたがっていた値を、カテゴリ変数と値の2列に変換することで、
wide-formatな `data.frame` をlong-formatに変形する。
`reshape2::melt()` に相当。

`tidyr::gather(data, key, value, ..., na.rm=FALSE, convert=FALSE)`

`data`
:   `%>%` 越しに渡す

`key`, `value`
:   出力結果で使う列名

`...`
:   動かす列名を指定。コロンによる範囲指定、マイナスによる除外指定も可。

e.g., `Species` 以外の列について、
元の列名を `kagi` 、値を `atai` に格納した縦長の表に変形

```r
iris %>% head(3) %>% add_rownames('id') %>%
    gather(kagi, atai, -id, -Species)

Source: local data frame [12 x 4]

      id Species         kagi  atai
   (chr)  (fctr)       (fctr) (dbl)
1      1  setosa Sepal.Length   5.1
2      2  setosa Sepal.Length   4.9
3      3  setosa Sepal.Length   4.7
4      1  setosa  Sepal.Width   3.5
5      2  setosa  Sepal.Width   3.0
6      3  setosa  Sepal.Width   3.2
7      1  setosa Petal.Length   1.4
8      2  setosa Petal.Length   1.4
9      3  setosa Petal.Length   1.3
10     1  setosa  Petal.Width   0.2
11     2  setosa  Petal.Width   0.2
12     3  setosa  Petal.Width   0.2
```

## `tidyr::spread()`

`tidyr::gather()` の逆で、
long-formatの `data.frame` をwide-formatに変形する。
`reshape2::dcast()` に相当。
IDとなるような列がないと `Error: Duplicate identifiers` と怒られる。

`tidyr::spread(data, key, value, fill=NA, convert=FALSE, drop=TRUE)`

`data`
:   `%>%` 越しに渡す

`key`
:   ここに指定したカテゴリ変数のぶんだけ新しい列が作られる

`value`
:   値が入ってる列

`fill=NA`
:   該当する組み合わせの値が存在しない場合に何で埋めるか

`convert=FALSE`
:   -

`drop=TRUE`
:   該当する組み合わせの行が存在しない場合に欠落させるか

e.g., `kagi` 内の文字列を新たな列名として横長の表に変形して `atai` を移す

```r
iris %>% head(3) %>% add_rownames('id') %>%
    gather(kagi, atai, -id, -Species) %>%
    spread(kagi, atai)

Source: local data frame [3 x 6]

     id Species Sepal.Length Sepal.Width Petal.Length Petal.Width
  (chr)  (fctr)        (dbl)       (dbl)        (dbl)       (dbl)
1     1  setosa          5.1         3.5          1.4         0.2
2     2  setosa          4.9         3.0          1.4         0.2
3     3  setosa          4.7         3.2          1.3         0.2
```

## `tidyr::separate()`

`reshape2::colsplit()` に相当。
`character` 列を任意のセパレータで複数の列に分割。

`tidyr::separate(data, col, into, sep='[^[:alnum:]]', remove=TRUE, convert=FALSE, extra='warn', fill='warn')`

`data`
:   `%>%` 越しに渡す

`col`
:   切り分けたい列の名前

`into`
:   切り分けたあとの新しい列名を文字列ベクタで

`sep='[^[:alnum:]]'`
:   セパレータを正規表現で。デフォルトはあらゆる非アルファベット。

`remove=TRUE`
:   切り分ける前の列を取り除くかどうか

`convert=FALSE`
:   -

`extra='warn'`
:   列数が揃わないときにどうするか: `warn`, `drop`, `merge`

`fill='warn'`
:   足りない場合にどっちから埋めるか: `warn`, `right`, `left`

`kagi` 列を `part`, `axis` という2列に分割

```r
iris %>% head(3) %>% add_rownames('id') %>%
   gather(kagi, atai, -id, -Species) %>%
   separate(kagi, c('part', 'axis'))

Source: local data frame [12 x 5]

      id Species  part   axis  atai
   (chr)  (fctr) (chr)  (chr) (dbl)
1      1  setosa Sepal Length   5.1
2      2  setosa Sepal Length   4.9
3      3  setosa Sepal Length   4.7
4      1  setosa Sepal  Width   3.5
5      2  setosa Sepal  Width   3.0
6      3  setosa Sepal  Width   3.2
7      1  setosa Petal Length   1.4
8      2  setosa Petal Length   1.4
9      3  setosa Petal Length   1.3
10     1  setosa Petal  Width   0.2
11     2  setosa Petal  Width   0.2
12     3  setosa Petal  Width   0.2
```

逆をやるのが `tidyr::unite(data, col, ..., sep='_', remove=TRUE)` 。

`tidyr::extract(data, col, into, regex, ...)`
を使えば正規表現でもっと細かく指定できる。

名前の似てる `tidyr::extract_numeric(x)` は
文字列から数字部分を抜き出して `numeric` で返す関数。

## `tidyr::nest()`, `tidyr::unnest()`

ネストされた `data.frame` とは、
`vector` ではなく `list` の列を持っていて、
ひとつのセルに複数の値を保持した状態のものを指す。

`nest()` で `data.frame` を圧縮し、
`unnest()` で `list` を展開してフラットにする

```r
iris %>% nest(-Species)

Source: local data frame [3 x 5]
Groups: <by row>

     Species Sepal.Length Sepal.Width Petal.Length Petal.Width
      (fctr)        (chr)       (chr)        (chr)       (chr)
1     setosa    <dbl[50]>   <dbl[50]>    <dbl[50]>   <dbl[50]>
2 versicolor    <dbl[50]>   <dbl[50]>    <dbl[50]>   <dbl[50]>
3  virginica    <dbl[50]>   <dbl[50]>    <dbl[50]>   <dbl[50]>
```

## `tidyr::expand()`

`base::expand.grid()` のラッパー。
指定した列の全ての組み合わせが登場するように、欠損値 `NA` の入った行を挿入する。

## `tidyr::replace_na()`

欠損値 `NA` を好きな値で置き換える。
これまでは `mutate(x= ifelse(is.na(x), 0, x))` のようにしてたところを

```r
df %>% replace_na(list(x=0, y='unknown'))
```

## `tidyr::complete()`

指定した列の全ての組み合わせが登場するように、
他の列を `NA` などにして補完する

```r
df %>% complete(col1, col2, fill=list(col1=0, col2='-'))
```

## `tidyr::fill()`

`NA` を、その列の直前の `NA` でない値で埋める。
えくせるでセルの結合とかやってしまって、
最初のセルにしか値が無いような場合に使うのかな？
