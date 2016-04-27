+++
title = 'dplyr'
subtitle = "高速data.frame処理"
tags = ["r", "hadley"]
[menu.main]
  parent = "rstats"
  weight = -75
+++

-   <https://github.com/hadley/dplyr>
-   <http://cran.r-project.org/web/packages/dplyr/vignettes/introduction.html>
-   <http://blog.rstudio.org/2014/01/17/introducing-dplyr/>
-   <http://www.rstudio.com/resources/cheatsheets/>

データフレームに特化した版 [plyr]({{< relref "plyr.md" >}}) 。
[tidyr]({{< relref "tidyr.md" >}}) と一緒に使うとよい。

`install.packages('dplyr')` でインストールし、
`library(dplyr)` で読み込んでから使う。

## 関数の連結

`plyr::ddply()` 的な処理をもっと柔軟に、見やすく書ける

```r
## plyr
plyr::ddply(plyr::mutate(subset(iris, Species!='setosa', select=-c(Sepal.Width, Sepal.Length)), petal_area=Petal.Length*Petal.Width*0.5), .(Species), numcolwise(mean))

## dplyr
iris %>%
   dplyr::filter(Species != 'setosa') %>%
   dplyr::select(-starts_with('Sepal')) %>%
   dplyr::mutate(petal_area=Petal.Length * Petal.Width * 0.5) %>%
   dplyr::group_by(Species) %>%
   dplyr::summarise_each(funs(mean))

## どちらも結果は
     Species Petal.Length Petal.Width petal_area
1 versicolor        4.260       1.326     2.8602
2  virginica        5.552       2.026     5.6481
```

`%>%`
:   左の値を右の関数に第一引数として渡す。
    `.data %>% func(arg1, arg2)` は `func(.data, arg1, arg2)` になる。
    処理する順に書けるので、次々と関数を適用していくときでも読みやすい。

    {{%div class="note"%}}
現状 `dplyr` は `magrittr` パッケージのものを採用しているが、
`pipeR` パッケージの `%>>%` のほうが柔軟で高速。
<http://renkun.me/pipeR-tutorial/>
    {{%/div%}}

`dplyr::group_by(.data, col1, col2, ..., add=FALSE)`
:   グループごとに区切って次の処理に渡す。 e.g. `summarise()`, `tally()`, `do()` など

`dplyr::rowwise(.data)`
:   行ごとに区切って次の処理に渡す

`dplyr::do(.data, ...)`
:   グループごとに処理する。
    `{}` 内に長い処理を書いてもいいし、関数に渡してもよい。
    グループごとに切りだされた部分は `.` で参照できる。
    出力が `data.frame` じゃないと
    `Error: Results are not data frames at positions: 1`
    のように怒られるが、
    `do(dummy=func(.))` のように名前付きにすると
    `data.frame` に入らないような型でも大丈夫になる。

## コア関数 (verb)

`dplyr::filter(.data, ...)`
:   列の値で行を絞る。`base::subset()` と似たようなもの

    ```r
    > iris %>% dplyr::filter(Sepal.Length<6, Sepal.Width>4)
      Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    1          5.7         4.4          1.5         0.4  setosa
    2          5.2         4.1          1.5         0.1  setosa
    3          5.5         4.2          1.4         0.2  setosa
    ```

`dplyr::slice(.data, ...)`
:   行数を指定して行を絞る。
    `` `[`(i,) `` の代わりに

    ```r
    > iris %>% dplyr::slice(2:4)
      Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    1          4.9         3.0          1.4         0.2  setosa
    2          4.7         3.2          1.3         0.2  setosa
    3          4.6         3.1          1.5         0.2  setosa
    ```

`dplyr::select(.data, ...)`
:   使う列を絞る。複数指定, 範囲指定、負の指定、パターン指定が可能\
    `starts_with(x, ignore.case=FALSE)`\
    `ends_width(x, ignore.case=FALSE)`\
    `contains(x, ignore.case=FALSE)`\
    `matches(x, ignore.case=FALSE)`\
    `num_range('x', 1:5, width=2)`

    ```r
    iris %>% dplyr::select(-(Sepal.Width:Petal.Length))
    iris %>% dplyr::select(ends_with('Length'))
    ```

    `.data[, j, drop=TRUE]` のように1列分をベクタで得たいときは二重角括弧か、
    `pipeR` の `%>>%` を通して括弧に流す。

    ```r
    iris %>% `[[`('Species')
    iris %>>% (Species)
    ```

`dplyr::rename(.data, ...)`
:   列の改名

    ```r
    > iris %>% dplyr::rename(sp=Species)
      Sepal.Length Sepal.Width Petal.Length Petal.Width     sp
    1          5.1         3.5          1.4         0.2 setosa
    2          4.9         3.0          1.4         0.2 setosa
    3          4.7         3.2          1.3         0.2 setosa
    ```

`dplyr::mutate(.data, ...)`
:   既存の列を使って新しい列を作る。
    `base::transform()` とほとんど同じだがそれよりも高速で高機能

    ```r
    # modify existing column
    iris %>% dplyr::mutate(Sepal.Length = log(Sepal.Length))

    # create new column
    iris %>% dplyr::mutate(ln_sepal_length = log(Sepal.Length))
    ```

`dplyr::transmute(.data, ...)`
:   指定した列以外を保持しない版の `mutate()`

`dplyr::distinct(.data, ...)`
:   指定した列に関してユニークな行のみ返す

    ```r
    > iris %>% distinct(Species)
      Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
    1          5.1         3.5          1.4         0.2     setosa
    2          7.0         3.2          4.7         1.4 versicolor
    3          6.3         3.3          6.0         2.5  virginica
    ```

`dplyr::arrange(.data, column1, column2, ...)`
:   指定した列の昇順で `data.frame` の行を並べ替える。
    `arrange(desc(column))` で降順になる。
    `order()` を使うよりもタイピングの繰り返しが少ないし直感的

    ```r
    .data[with(.data, order(col_a, col_b)), ]
    # is equivalent to
    .data %>% dplyr::arrange(col_a, col_b)
    ```

`dplyr::summarise(.data, ...)`
:   列に対する関数をグループごとに適用して1行にしたものを `rbind()` してまとめる

    ```r
    > iris %>% group_by(Species) %>% summarise(minpl=min(Petal.Length), maxsw=max(Sepal.Width))
         Species minpl maxsw
    1     setosa   1.0   4.4
    2 versicolor   3.0   3.4
    3  virginica   4.5   3.8
    ```

`dplyr::summarise_each(.data, funs, ...)`, `dplyr::mutate_each(.data, funs, ...)`
:   複数の列に複数の関数を適用する(`dplyr::funs()` と共に)

    ```r
    iris %>% group_by(Species) %>% summarise_each(funs(min, mean, max), ends_with("Width"))
    ```

## その他の関数

`dplyr::data_frame(...)`
:   より便利でバグの混入しにくい方法で `data.frame` を作る。

    -   勝手に型変換しない (`stringsAsFactors=FALSE`)
    -   勝手に列名を変えない
    -   長さ1の変数以外はリサイクルしない
    -   引数の評価がlazyに行われるので前の列を利用して後の列を作ったりできる
    -   `tbl_df` クラスを付加
    -   ただし `matrix` からの直接変換はサポートしなそう
        <https://github.com/hadley/dplyr/issues/324>

`dplyr::***_join(x, y, by=NULL, copy=FALSE)`
:   `by` で指定した列がマッチするように行を合わせて `cbind()`

    `full_join()`: `x` と `y` の全ての行を保持。\
    `inner_join()`: `x` と `y` の `by` がマッチする行のみ\
    `left_join()`: `x` の全ての行を保持。`y` に複数マッチする行があったらすべて保持。\
    `right_join()`: `y` の全ての行を保持。`x` に複数マッチする行があったらすべて保持。\
    `semi_join()`: `x` の全ての行を保持。`y` に複数マッチする行があっても元の `x` の行だけ保持。\
    `anti_join()`: `y` にマッチしない `x` の行のみ。

    列名が異なる場合は `by` を名前付きベクタにすればよい

`dplyr::bind_rows(...)`, `dplyr::bind_cols(...)`
:   標準の `rbind()`, `cbind()` より効率よく `data.frame` を結合。
    引数は個別でもリストでもよい。

`dplyr::add_rownames(x, var='rowname')`
:   1列目に行名として1からの整数を振る。

`dplyr::ntile(x, n)`
:   数値ベクトル `x` を順位によって `n` 個のクラスに均等分け

`dplyr::n_distinct(x)`
:   高速で簡潔な `length(unique(x))`

`dplyr::last(x, order_by=NULL, default=default_missing(x))`
:   最後の要素にアクセス。
    `x[length(x)]` と書かなくて済む！

`dplyr::glimpse(.data, width=getOption('width'))`
:   データの中身をざっと見る。
    `print()` とか `str()` のようなもの。

`dplyr::lead(x n=1, default=NA, order_by=NULL)`, `dplyr::lag(...)`
:   `x` の中身を `n` だけずらして `default` で埋める。
    `lead()` は前に、`lag()` は後ろにずらす

    ```r
    > lag(seq_len(5), 2)
    [1] NA NA 1 2 3
    ```

`dplyr::top_n(.data, n, wt=NULL)`
:   `.data %>% arrange(wt) %>% head(n)` を一撃で、グループごとに。

`dplyr::tally(x, wt, sort=FALSE)`
:   `summarise(x, n=n())` のショートカット。
    `wt` にカラムを指定して重み付けすることもできる。

`dplyr::count(x, ..., wt=NULL, sort=FALSE)`
:   `group_by(...) %>% tally()` のショートカット。

`dplyr::between(x, left, right)`
:   `left <= x & x <= right` に相当
