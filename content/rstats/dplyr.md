+++
title = 'dplyr'
subtitle = "高速data.frame処理"
tags = ["r", "hadley"]
[menu.main]
  parent = "rstats"
  weight = -75
+++

-   <https://github.com/hadley/dplyr>
-   <https://cran.r-project.org/web/packages/dplyr/vignettes/introduction.html>
-   <https://blog.rstudio.org/2014/01/17/introducing-dplyr/>
-   <https://www.rstudio.com/resources/cheatsheets/>

data.frameに対して抽出(select, filter)、部分的変更(mutate)、要約(summarise)、ソート(arrange)などの処理を施すためのパッケージ。
前作 [plyr]({{< relref "plyr.md" >}}) のうちdata.frameに関する部分を抜き出して強化したパッケージ。
[purrr]({{< relref "purrr.md" >}}) や [tidyr]({{< relref "tidyr.md" >}}) と一緒に使うとよい。

Rの中で `install.packages('dplyr')` としてインストールし、
使う前に `library(dplyr)` で読み込む。

## 関数の連結

一時変数を作ったり、関数を何重にも重ねて書いたりすることなく、
適用する順に処理を記述することができる。

```r
## dplyr
iris %>%
   dplyr::filter(Species != 'setosa') %>%
   dplyr::select(-dplyr::starts_with('Sepal')) %>%
   dplyr::mutate(petal_area=Petal.Length * Petal.Width * 0.5) %>%
   dplyr::group_by(Species) %>%
   dplyr::summarise_each(funs(mean))

## plyr で書くと読みにくい
plyr::ddply(plyr::mutate(subset(iris, Species!='setosa', select=-c(Sepal.Width, Sepal.Length)), petal_area=Petal.Length*Petal.Width*0.5), .(Species), numcolwise(mean))

## どちらも結果は
     Species Petal.Length Petal.Width petal_area
1 versicolor        4.260       1.326     2.8602
2  virginica        5.552       2.026     5.6481
```

`%>%`
:   左の値を右の関数に第一引数として渡す。
    `.data %>% func(arg1, arg2)` は `func(.data, arg1, arg2)` と等価等価。
    処理する順に書けるので、次々と関数を適用していくときでも読みやすい。

    {{%div class="note"%}}
現状 `dplyr` は `magrittr` パッケージのものを採用しているが、
`pipeR` パッケージの `%>>%` のほうが柔軟で高速。
<http://renkun.me/pipeR-tutorial/>
    {{%/div%}}

`dplyr::group_by(.data, col1, col2, ..., add=FALSE)`
:   グループごとに区切って次の処理に渡す。 e.g. `summarise()`, `tally()`, `do()` など

`dplyr::rowwise(.data)`
:   行ごとに区切って次の処理に渡す。
    今後この方法は非推奨となるっぽいので、
    [`purrr`]({{< relref "purrr.md" >}})`::by_row()`を使ったほうが良い。

`dplyr::do(.data, ...)`
:   グループごとに処理する。
    `{}` 内に長い処理を書いてもいいし、関数に渡してもよい。
    グループごとに切りだされた部分は `.` で参照できる。
    出力がdata.frameじゃないと
    `Error: Results are not data frames at positions: 1`
    のように怒られるが、
    `do(dummy=func(.))` のように名前付きにすると
    data.frameに入らないような型でも大丈夫になる。
    今後は[`purrr`]({{< relref "purrr.md" >}})`::by_slice()`を使ったほうが良さそう。

## コア関数 (verb)

`dplyr::select(.data, ...)`
:   使う列を絞る。複数指定, 範囲指定、負の指定、パターン指定が可能\
    `starts_with(x, ignore.case=FALSE)`\
    `ends_width(x, ignore.case=FALSE)`\
    `contains(x, ignore.case=FALSE)`\
    `matches(x, ignore.case=FALSE)`\
    `num_range('x', 1:5, width=2)`

    ```r
    iris %>% dplyr::select(-(Sepal.Width:Petal.Length))
    iris %>% dplyr::select(dplyr::ends_with('Length'))
    ```

    `.data[, j, drop=TRUE]` のように1列分をベクタで得たいときは二重角括弧か、
    `pipeR` の `%>>%` を通して括弧に流す。
    ```r
    iris %>% `[[`('Species')
    iris %>>% (Species)
    ```

`dplyr::filter(.data, ...)`
:   列の値で行を絞る。`base::subset()` と似たようなもの
    ```r
    > iris %>% dplyr::filter(Sepal.Length<6, Sepal.Width>4)
      Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    1          5.7         4.4          1.5         0.4  setosa
    2          5.2         4.1          1.5         0.1  setosa
    3          5.5         4.2          1.4         0.2  setosa
    ```

`dplyr::distinct(.data, ..., .keep_all=FALSE)`
:   指定した列に関してユニークな行のみ返す。
    `filter(!duplicated(.[, ...]))` をよりスマートに。
    指定しなかった列も残すには `.keep_all=TRUE` が必要。
    ```r
    > iris %>% dplyr::distinct(Species)
      Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
    1          5.1         3.5          1.4         0.2     setosa
    2          7.0         3.2          4.7         1.4 versicolor
    3          6.3         3.3          6.0         2.5  virginica
    ```

`dplyr::slice(.data, ...)`
:   行番号を指定して行を絞る。
    `` `[`(i,) `` の代わりに
    ```r
    > iris %>% dplyr::slice(2:4)
      Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    1          4.9         3.0          1.4         0.2  setosa
    2          4.7         3.2          1.3         0.2  setosa
    3          4.6         3.1          1.5         0.2  setosa
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

`dplyr::rename(.data, ...)`
:   列の改名。
    `mutate()`と同じようなイメージで `new=old` と指定。
    ```r
    > iris %>% dplyr::rename(sp=Species)
      Sepal.Length Sepal.Width Petal.Length Petal.Width     sp
    1          5.1         3.5          1.4         0.2 setosa
    2          4.9         3.0          1.4         0.2 setosa
    3          4.7         3.2          1.3         0.2 setosa
    ```

`dplyr::arrange(.data, column1, column2, ...)`
:   指定した列の昇順でdata.frameの行を並べ替える。
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
    > iris %>% dplyr::group_by(Species) %>%
        dplyr::summarise(minpl=min(Petal.Length), maxsw=max(Sepal.Width))
         Species minpl maxsw
    1     setosa   1.0   4.4
    2 versicolor   3.0   3.4
    3  virginica   4.5   3.8
    ```

`dplyr::summarise_all(.data, .funs, ...)`, `dplyr::mutate_all(.data, .funs, ...)`
:   グループ化カラム以外の全てのカラムに関数を適用する。
    ```r
    iris %>% dplyr::group_by(Species) %>%
        dplyr::summarise_all(funs(min, mean, max))
    ```

`dplyr::summarise_at(.data, .cols, .funs, ...)`, `dplyr::mutate_all(.data, .cols, .funs, ...)`
:   select補助関数を使って指定したカラムに関数を適用する。
    `***_each()` は非推奨になった。
    ```r
    iris %>% dplyr::group_by(Species) %>%
        dplyr::summarise_at(vars(dplyr::ends_with("Width")), funs(min, mean, max))
    ```

`dplyr::summarise_if(.data, .predicate, .funs, ...)`, `dplyr::mutate_if(.data, .predicate, .funs, ...)`
:   `.predicate`がTRUEになるカラムだけに関数を適用する。
    ```r
    iris %>% dplyr::group_by(Species) %>%
        dplyr::summarise_if(is.numeric, funs(min, mean, max))
    ```

## その他の関数

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
:   標準の `rbind()`, `cbind()` より効率よくdata.frameを結合。
    引数は個別でもリストでもよい。
    そのほかにも標準の集合関数を置き換えるものが提供されている:
    `intersect()`, `union()`, `union_all()`, `setdiff()`, `setequal()`

`dplyr::if_else(condition, true, false, missing=NULL)`
:   標準の`ifelse()`よりも型に厳しく、高速らしい。
    NAのときにどうするかを指定できるのも大変良い。

`dplyr::coalesce(x, ...)`
:   最初のvectorでNAだったとこは次のvectorのやつを採用、
    という`ifelse(!is.na(x), x, y)`的な処理をする。
    基本的には同じ長さのvectorを渡すが、
    2つめに長さ1のを渡して`tidyr::replace_na()`的に使うのも便利。

`dplyr::na_if(x, y)`
:   `x[x == y] = NA; x` のショートカット

`dplyr::recode(.x, ..., .default=NULL, .missing=NULL)`
:    vectorの値を変更する。e.g.,
     `recode(letters, a='A!', c='C!')`

`dplyr::add_rownames(x, var='rowname')`
:   1列目に行名として1からの整数を振る。

`dplyr::group_indices(.data, ...)`
:   `grouped_df` ではなくグループIDとして1からの整数を返す版 `group_by()`

`dplyr::ntile(x, n)`
:   数値ベクトル `x` を順位によって `n` 個のクラスに均等分け

`dplyr::n_distinct(x)`
:   高速で簡潔な `length(unique(x))`

`dplyr::last(x, order_by=NULL, default=default_missing(x))`
:   最後の要素にアクセス。
    `x[length(x)]` や `tail(x, 1)` よりも楽チンだが、
    安全性重視の `dplyr::nth()` を内部で呼び出すため遅い。

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
:   `left <= x & x <= right` のショートカット。

`dplyr::near(x, y, tol=.Machine$double.eps^0.5)`
:   `abs(x - y) < tol` のショートカット。

`dplyr::case_when(...)`
:   `if {} else if {} else if {} ...` のショートカット。


## tibble

`tbl_df` クラスが付与された改良版data.frameのことを**tibble**と呼ぶ。
もともとはdplyrパッケージで扱っていたが、tibbleパッケージとして独立した。
[readr]({{< relref "readr.md" >}}) で読み込んだデータもこの形式になる。

```r
> tbl_iris = as_tibble(iris)
> class(tbl_iris)
[1] "tbl_df"     "tbl"        "data.frame"
> class(iris)
[1] "data.frame"
```

生のdata.frameとの違いは:

-   うっかり巨大なデータを`print()`しても画面を埋め尽くさない
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

`tibble::glimpse(.data, width=NULL)`
:   データの中身をざっと見る。
    `print()` とか `str()` のようなもの。

`tibble::type_sum(x)`
:   オブジェクトの型

`tibble::obj_sum(x)`
:   `type_sum`とサイズ e.g., `"data.frame [150 x 5]"`

`tibble::remove_rownames(df)` \
`tibble::rownames_to_column(df, var='rowname')` \
`tibble::column_to_rownames(df, var='rowname')`

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
