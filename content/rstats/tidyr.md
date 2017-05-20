+++
title = 'tidyr'
subtitle = "シンプルなデータ変形ツール"
tags = ["r", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -70
+++

data.frameを縦長・横長・入れ子に変形・整形するためのツール。
[dplyr]({{< relref "dplyr.md" >}}) や [purrr]({{< relref "purrr.md" >}})
と一緒に使うとよい。
[reshape2]({{< relref "reshape2.md" >}}) を置き換えるべく再設計された改良版。

[tidyverse](https://github.com/tidyverse/tidyverse) に含まれているので、
`install.packages('tidyverse')` で一括インストール、
`library(tidyverse)` で一括ロード。

-   <http://r4ds.had.co.nz/tidy-data.html>
-   <https://github.com/tidyverse/tidyr>
-   `vignette("tidy-data")`
-   `demo(package="tidyr")`

下記のコード例で使うデータ

```r
> iris %>% head(3L) %>% rownames_to_column('id')

  id Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1  1          5.1         3.5          1.4         0.2  setosa
2  2          4.9         3.0          1.4         0.2  setosa
3  3          4.7         3.2          1.3         0.2  setosa
```

パイプ演算子 `%>%` については[dplyr]({{< relref "dplyr.md" >}})を参照。


## `tidyr::gather()` で縦長にする

複数列にまたがっていた値を、カテゴリ変数と値の2列に変換することで、
横長(wide-format)のdata.frameを縦長(long-format)に変形する。
[ggplot2]({{< relref "ggplot2.md" >}})で使いやすいのはこの縦長。
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
> iris %>% head(3L) %>% rownames_to_column('id') %>%
    gather(kagi, atai, -id, -Species)

   id Species         kagi atai
1   1  setosa Sepal.Length  5.1
2   2  setosa Sepal.Length  4.9
3   3  setosa Sepal.Length  4.7
4   1  setosa  Sepal.Width  3.5
5   2  setosa  Sepal.Width  3.0
6   3  setosa  Sepal.Width  3.2
7   1  setosa Petal.Length  1.4
8   2  setosa Petal.Length  1.4
9   3  setosa Petal.Length  1.3
10  1  setosa  Petal.Width  0.2
11  2  setosa  Petal.Width  0.2
12  3  setosa  Petal.Width  0.2
```

## `tidyr::spread()` で横長にする

`tidyr::gather()` の逆で、縦長のdata.frameを横長に変形する。
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
> iris %>% head(3L) %>% rownames_to_column('id') %>%
    gather(kagi, atai, -id, -Species) %>%
    spread(kagi, atai)

  id Species Petal.Length Petal.Width Sepal.Length Sepal.Width
1  1  setosa          1.4         0.2          5.1         3.5
2  2  setosa          1.4         0.2          4.9         3.0
3  3  setosa          1.3         0.2          4.7         3.2
```

## Nested data.frame --- 入れ子構造

### `tidyr::nest(data, ..., .key=data)`

data.frameをネストして(入れ子にして)、list of data.frames のカラムを作る。
内側のdata.frameに押し込むカラムを `...` に指定するか、
外側に残すカラムをマイナス指定する。

```r
iris %>% nest(-Species, .key=NEW_COLUMN)
# A tibble: 3 × 2
     Species        NEW_COLUMN
      <fctr>            <list>
1     setosa <tibble [50 × 4]>
2 versicolor <tibble [50 × 4]>
3  virginica <tibble [50 × 4]>

# equivalent to
iris %>% dplyr::group_by(Species) %>% nest()
iris %>% nest(matches('Length$|Width$'))
```

なんでもかんでもフラットなdata.frameにして
[dplyr]({{< relref "dplyr.md" >}})を駆使する時代は終わり、
ネストしておいて[purrr]({{< relref "purrr.md" >}})を適用するのが
tidyverse時代のクールなやり方らしい。

cf. [Hadley Wickham: Managing many models with R (YouTube)](https://www.youtube.com/watch?v=rz3_FDVt9eg)


### `tidyr::unnest(data, ..., .drop=NA, id=NULL, .sep=NULL)`

ネストされたdata.frameを展開してフラットにする。
list of data.framesだけでなく、list of vectorsとかでもよい。


## その他の便利関数

### `tidyr::separate()`

文字列カラムを任意のセパレータで複数カラムに分割。
`reshape2::colsplit()` に相当。

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
:   足りない場合にどっち側をNAで埋めるか: `warn`, `right`, `left`。
    つまり、文字を左詰めにするには`right`が正解(紛らわしい)。

`kagi` 列を `part`, `axis` という2列に分割

```r
> iris %>% head(3L) %>% rownames_to_column('id') %>%
    gather(kagi, atai, -id, -Species) %>%
    separate(kagi, c('part', 'axis'))

   id Species  part   axis atai
1   1  setosa Sepal Length  5.1
2   2  setosa Sepal Length  4.9
3   3  setosa Sepal Length  4.7
4   1  setosa Sepal  Width  3.5
5   2  setosa Sepal  Width  3.0
6   3  setosa Sepal  Width  3.2
7   1  setosa Petal Length  1.4
8   2  setosa Petal Length  1.4
9   3  setosa Petal Length  1.3
10  1  setosa Petal  Width  0.2
11  2  setosa Petal  Width  0.2
12  3  setosa Petal  Width  0.2
```

逆をやるのが `tidyr::unite(data, col, ..., sep='_', remove=TRUE)` 。

`tidyr::extract(data, col, into, regex, ...)`
を使えば正規表現でもっと細かく指定できる。

名前の似てる `tidyr::extract_numeric(x)` は
文字列から数字部分をnumericとして抜き出す関数だったが今はdeprecatedなので、
新しい[`readr::parse_number()`]({{< relref "readr.md" >}})を使うべし。

### `tidyr::complete(data, ..., fill=list())`

指定した列の全ての組み合わせが登場するように、
指定しなかった列に欠損値`NA`(あるいは任意の値)を補完した行を挿入する。

```r
df %>% complete(key1, key2, fill=list(val1=0, val2='-'))
```

### `tidyr::expand(data, ...)`

指定した列の全ての組み合わせが登場するような新しいdata.frameを作る。
全ての列を指定すれば`complete()`と同じ効果だが、
指定しなかった列が消えるという点では異なる。

`crossing(...)`はvectorを引数に取る亜種で、
tibble版`expand.grid(...)`のようなもの。

`nesting(...)`は存在するユニークな組み合わせのみ残す、
`nest(data, ...) %>% dplyr::select(-data)`のショートカット。
この結果は`expand()`や`complete()`の引数としても使える。

数値vectorの補完には`full_seq(x, period)`が便利。


### `tidyr::drop_na(data, ...)`

`complete()`の逆。
指定した列に`NA`が含まれてる行を削除する。
何も指定しなければ標準の `data[complete.cases(data),]` と同じ。

### `tidyr::replace_na()`

欠損値 `NA` を好きな値で置き換える。
これまでは `mutate(x= ifelse(is.na(x), 0, x))` のようにしてたところを

```r
df %>% replace_na(list(x=0, y='unknown'))
```

逆に、特定の値を`NA`にしたい場合は
[`dplyr::na_if()`]({{< relref "dplyr.md" >}})


### `tidyr::fill()`

`NA` を、その列の直前の `NA` でない値で埋める。
えくせるでセルの結合とかやってしまって、
最初のセルにしか値が無いような場合に使うのかな？


## 関連書籍

<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?_encoding=UTF8&qid=1485613345&sr=1-1-catcorr&linkCode=li3&tag=heavywatal-22&linkId=6133a1fd9babbf590e304ee9f670fa4a" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
