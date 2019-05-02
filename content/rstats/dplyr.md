+++
title = 'dplyr'
subtitle = "高速data.frame処理"
tags = ["r", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -75
+++

<a href="https://dplyr.tidyverse.org/">
<img src="https://dplyr.tidyverse.org/logo.png" align="right" width="120" height="139">
</a>

data.frameに対して抽出(select, filter)、部分的変更(mutate)、要約(summarise)、ソート(arrange)などの処理を施すためのパッケージ。
前作 [plyr]({{< relref "plyr.md" >}}) のうちdata.frameに関する部分が強化されている。
[purrr]({{< relref "purrr.md" >}}) や [tidyr]({{< relref "tidyr.md" >}}) と一緒に使うとよい。

[tidyverse](https://tidyverse.tidyverse.org/) に含まれているので、
`install.packages("tidyverse")` で一括インストール、
`library(tidyverse)` で一括ロード。

-   <http://r4ds.had.co.nz/transform.html>
-   <https://github.com/tidyverse/dplyr>

## 関数の連結 %>%

<a href="https://magrittr.tidyverse.org/">
<img src="https://magrittr.tidyverse.org/logo.png" align="right" width="120" height="139">
</a>

dplyrではなく[magrittr](https://magrittr.tidyverse.org/)の機能。

`x %>% f(a, b)` は `f(x, a, b)` と等価。
左の値 `x` を第一引数として右の関数 `f()` に渡す。
一時変数を作ったり、関数を何重にも重ねたりすることなく、
適用する順に次々と処理を記述することができるようになる。
慣れれば書きやすく読みやすい。

```r
library(tidyverse)

## with piping
iris %>%
  dplyr::filter(Species != "setosa") %>%
  dplyr::select(-dplyr::starts_with("Sepal")) %>%
  dplyr::mutate(petal_area=Petal.Length * Petal.Width * 0.5) %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise_all(funs(mean))

## with a temporary variable
x = iris
x = dplyr::filter(x, Species != "setosa")
x = dplyr::select(x, -dplyr::starts_with("Sepal"))
x = dplyr::mutate(x, petal_area=Petal.Length * Petal.Width * 0.5)
x = dplyr::group_by(x, Species)
x = dplyr::summarise_all(x, funs(mean))

## with nested functions
dplyr::summarise_all(
  dplyr::group_by(
    dplyr::mutate(
      dplyr::select(
        dplyr::filter(iris, Species != "setosa"),
        -dplyr::starts_with("Sepal")),
      petal_area=Petal.Length * Petal.Width * 0.5),
    Species),
  funs(mean))

## result
     Species Petal.Length Petal.Width petal_area
1 versicolor        4.260       1.326     2.8602
2  virginica        5.552       2.026     5.6481
```

現状では `magrittr` パッケージの `%>%` が広く採用されているが、
`pipeR` パッケージの `%>>%` のほうが高速らしい。
<https://renkun-ken.github.io/pipeR-tutorial/>


## 抽出・絞り込み

### 列

`dplyr::select(.data, ...)`
:   列を絞る。複数指定、範囲指定、負の指定が可能。
    [select helper](https://dplyr.tidyverse.org/reference/select_helpers.html)
    によるパターン指定も便利。
    残るのが1列だけでも勝手にvectorにはならずdata.frameのまま。

    ```r
    iris %>% dplyr::select(Petal.Width, Species)
    iris %>% dplyr::select("Petal.Width", "Species")
    iris %>% dplyr::select(c("Petal.Width", "Species"))
    iris %>% dplyr::select(4:5)
    iris %>% dplyr::select(-c(1:3))
    iris %>% dplyr::select(-(Sepal.Length:Petal.Length))
    iris %>% dplyr::select(matches("^Petal\\.Width$|^Species$"))
    ```

:   文字列変数で指定しようとすると意図が曖昧になるので、
    [unquoting](https://dplyr.tidyverse.org/articles/programming.html#unquoting)
    やpronounで明確に:
    ```r
    Sepal.Length = c("Petal.Width", "Species")
    iris %>% dplyr::select(Sepal.Length)       # ambiguous!
    iris %>% dplyr::select(.data$Sepal.Length) # pronoun => Sepal.Length
    iris %>% dplyr::select(!!Sepal.Length)     # unquote => Petal.Width, Species
    iris %>% dplyr::select(!!!rlang::syms(Sepal.Length))  # Petal.Width, Species
    ```
    これらの指定方法は `rename()` や `pull()` でも有効。
    一方、文字列を受け取れない `distinct()` や `group_by()`
    などの関数には普通のunquoteは通用しない。
    最後の例のように `rlang::syms()` でシンボル化して
    [unquote-splicing](https://dplyr.tidyverse.org/articles/programming.html#unquote-splicing)
    して渡す必要がある。
    ```r
    columns = c("Petal.Width", "Species")
    iris %>% distinct(!!as.name(columns[1L]))
    iris %>% distinct(!!!rlang::syms(columns))
    ```
    詳しくは[宇宙本](https://amzn.to/2u0hmTs)第3章のコラム
    「selectのセマンティクスとmutateのセマンティクス、tidyeval」を参照。

:   列の値によって選べる亜種 `dplyr::select_if(.tbl, .predicate, ...)` もある。


`dplyr::rename(.data, ...)`
:   列の改名。
    `mutate()`と同じようなイメージで `new=old` と指定。
    ```r
    iris %>% dplyr::rename(sp=Species) %>% head(2)
      Sepal.Length Sepal.Width Petal.Length Petal.Width     sp
    1          5.1         3.5          1.4         0.2 setosa
    2          4.9         3.0          1.4         0.2 setosa
    ```
    変数に入った文字列を使う場合も`mutate()`と同様にunquotingで:
    ```r
    old_name = "Species"
    new_name = toupper(old_name)
    iris %>% dplyr::rename(!!new_name := !!old_name) %>% head(2)
      Sepal.Length Sepal.Width Petal.Length Petal.Width SPECIES
    1          5.1         3.5          1.4         0.2  setosa
    2          4.9         3.0          1.4         0.2  setosa
    ```
    名前付きベクターと
    [unquote-splicing](https://dplyr.tidyverse.org/articles/programming.html#unquote-splicing)
    を使えば一括指定できる:
    ```r
    named_vec = setNames(names(iris), LETTERS[1:5])
    iris %>% dplyr::rename(!!!named_vec) %>% head(2L)
        A   B   C   D      E
    1 5.1 3.5 1.4 0.2 setosa
    2 4.9 3.0 1.4 0.2 setosa
    ```

:   リネーム関数を渡せる亜種:<br>
    `rename_all(.tbl, .funs = list(), ...)`<br>
    `rename_at(.tbl, .vars, .funs = list(), ...)`<br>
    `rename_if(.tbl, .predicate, .funs = list(), ...)`<br>

    ```r
    pattern = c("^Sepal" = "Gaku", "^Petal" = "Kaben")
    iris %>% head() %>%
      dplyr::rename_all(stringr::str_replace_all, pattern) %>%
      print()
    ```

`dplyr::pull(.data, var=-1)`
:   指定した1列をvector(またはlist)としてdata.frameから抜き出す。

    ```r
    iris %>% head() %>% dplyr::pull(Species)
    iris %>% head() %>% dplyr::pull("Species")
    iris %>% head() %>% dplyr::pull(5)
    iris %>% head() %>% dplyr::pull(-1)
    iris %>% head() %>% `[[`("Species")
    iris %>% head() %>% {.[["Species"]]}
    iris %>% head() %>% {.$Species}
    {iris %>% head()}$Species
    ```


### 行

`dplyr::filter(.data, ...)`
:   条件を満たす行だけを返す。`base::subset()` と似たようなもの。

    ```r
    iris %>% dplyr::filter(Sepal.Length<6, Sepal.Width>4)
      Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    1          5.7         4.4          1.5         0.4  setosa
    2          5.2         4.1          1.5         0.1  setosa
    3          5.5         4.2          1.4         0.2  setosa
    ```

    <div class="warning">{{<markdownify>}}
評価結果が `NA` となる行は除去される。
特に不等号を使うときやや直感に反するので要注意。
e.g., `filter(gene != "TP53")`
{{</markdownify>}}</div>

:   複数列で評価する亜種:<br>
    `filter_all(.tbl, .vars_predicate)`<br>
    `filter_if(.tbl, .predicate, .vars_predicate)`<br>
    `filter_at(.tbl, .vars, .vars_predicate)`<br>


`dplyr::distinct(.data, ..., .keep_all=FALSE)`
:   指定した列に関してユニークな行のみ返す。
    `base::unique.data.frame()` よりも高速で、
    `filter(!duplicated(.[, ...]))` よりスマートで柔軟。
    指定しなかった列を残すには `.keep_all=TRUE` とする。

    ```r
    iris %>% dplyr::distinct(Species)
         Species
    1     setosa
    2 versicolor
    3  virginica
    ```

`dplyr::slice(.data, ...)`
:   行番号を指定して行を絞る。
    `` `[`(i,) `` の代わりに。

    ```r
    iris %>% dplyr::slice(2:4)
      Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    1          4.9         3.0          1.4         0.2  setosa
    2          4.7         3.2          1.3         0.2  setosa
    3          4.6         3.1          1.5         0.2  setosa
    ```

`dplyr::sample_n(tbl, size, replace=FALSE, weight=NULL)`
:   指定した行数だけランダムサンプルする。
    割合指定の `sample_frac()` もある。


## 列の変更・追加

`dplyr::mutate(.data, ...)`
:   既存の列を変更したり、新しい列を作ったり。
    `base::transform()` の改良版。
    ```r
    # modify existing column
    iris %>% dplyr::mutate(Sepal.Length = log(Sepal.Length))

    # create new column
    iris %>% dplyr::mutate(ln_sepal_length = log(Sepal.Length))
    ```

    変数に入った文字列を列名として使いたい場合は `!!` と
    [unquoting用の代入演算子 `:=`](https://dplyr.tidyverse.org/articles/programming.html#setting-variable-names)
    を使う:
    ```r
    y = "new_column"
    x = "Sepal.Length"

    # unquoting only right hand side
    iris %>% dplyr::mutate(y = log(!!as.name(x))) %>% head(2)
      Sepal.Length Sepal.Width Petal.Length Petal.Width Species        y
    1          5.1         3.5          1.4         0.2  setosa 1.629241
    2          4.9         3.0          1.4         0.2  setosa 1.589235

    # unquoting both sides
    iris %>% dplyr::mutate(!!y := log(!!as.name(x))) %>% head(2)
      Sepal.Length Sepal.Width Petal.Length Petal.Width Species new_column
    1          5.1         3.5          1.4         0.2  setosa   1.629241
    2          4.9         3.0          1.4         0.2  setosa   1.589235
    ```
:   複数列を変更する亜種:<br>
    `mutate_all(.tbl, .funs, ...)`<br>
    `mutate_at(.tbl, .vars, .funs, ..., .cols = NULL)`<br>
    `mutate_each(.tbl, .funs, ...)`<br>
    `mutate_if(.tbl, .predicate, .funs, ...)`<br>

`dplyr::transmute(.data, ...)`
:   指定した列以外を保持しない版 `mutate()` 。
    言い換えると、列の中身の変更もできる版 `select()` 。
:   複数列を変更する亜種:<br>
    `transmute_all(.tbl, .funs, ...)`<br>
    `transmute_at(.tbl, .vars, .funs, ..., .cols = NULL)`<br>
    `transmute_if(.tbl, .predicate, .funs, ...)`<br>


## data.frameの要約・集計・整列

`dplyr::summarise(.data, ...)`
:   指定した列に関数を適用して1行のdata.frameにまとめる。
    グループ化されていたらグループごとに適用して `bind_rows()` する。

    ```r
    iris %>% dplyr::group_by(Species) %>%
        dplyr::summarise(minpl=min(Petal.Length), maxsw=max(Sepal.Width))
         Species minpl maxsw
    1     setosa   1.0   4.4
    2 versicolor   3.0   3.4
    3  virginica   4.5   3.8
    ```

`dplyr::summarise_all(.data, .funs, ...)`
:   全てのカラムに関数を適用する。
    `***_each()` は非推奨になった。

    ```r
    iris %>% dplyr::group_by(Species) %>%
        dplyr::summarise_all(funs(min, max))
    ```

`dplyr::summarise_at(.data, .cols, .funs, ...)`
:   select補助関数を使って指定したカラムに関数を適用する。

    ```r
    iris %>% dplyr::group_by(Species) %>%
        dplyr::summarise_at(vars(dplyr::ends_with("Width")), funs(min, max))
    ```

`dplyr::summarise_if(.data, .predicate, .funs, ...)`
:   `.predicate`がTRUEになるカラムだけに関数を適用する。

    ```r
    iris %>% dplyr::group_by(Species) %>%
        dplyr::summarise_if(is.numeric, funs(min, max))
    ```

`dplyr::tally(x, wt, sort=FALSE)`
:   `summarise(x, n=n())` のショートカット。
    `wt` にカラムを指定して重み付けすることもできる。
:   `dplyr::add_tally()` は元の形を維持したままカウント列を追加。

`dplyr::count(x, ..., wt=NULL, sort=FALSE)`
:   `group_by(...) %>% tally()` のショートカット。
:   `dplyr::add_count()` は元の形を維持したままカウント列を追加。

`dplyr::arrange(.data, column1, column2, ...)`
:   指定した列の昇順でdata.frameの行を並べ替える。
    `arrange(desc(column))` で降順になる。
    `order()` を使うよりもタイピングの繰り返しが少ないし直感的
    ```r
    .data[order(.data$col_a, .data$col_b),]
    # is equivalent to
    .data %>% dplyr::arrange(col_a, col_b)
    ```

`dplyr::top_n(.data, n, wt)`
:   `wt` で指定した列(省略時は右端の列)の降順で `n` 行だけ返す。
    境界にタイがある場合はすべて保持されて `n` 行以上の出力になる。
    元の順番が保持されるので、ソートされた結果が欲しい場合は事後に
    `arrange()` を適用することになる。
    順序の方向が `arrange()` と逆であることに注意。
    昇順で取りたいときは `n` にマイナス指定か `wt` に `desc(X)` (昇順なのに！)。


## data.frameを結合

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


## その他の関数

主に`mutate()`や`filter()`を補助するもの

`dplyr::if_else(condition, true, false, missing=NULL)`
:   標準の`ifelse()`よりも型に厳しく、高速らしい。
    NAのときにどうするかを指定できるのも大変良い。
    ネストせずにスッキリ書ける `dplyr::case_when()` も便利。

`dplyr::coalesce(x, ...)`
:   最初のvectorでNAだったとこは次のvectorのやつを採用、
    という`ifelse(!is.na(x), x, y)`的な処理をする。
    基本的には同じ長さのvectorを渡すが、
    2つめに長さ1のを渡して`tidyr::replace_na()`的に使うのも便利。

`dplyr::na_if(x, y)`
:   `x[x == y] = NA; x` のショートカット

`dplyr::recode(.x, ..., .default=NULL, .missing=NULL)`
:   vectorの値を変更する。e.g.,
    `recode(letters, a="A!", c="C!")`

`dplyr::row_number(x)`
:   `rank(x, ties.method = "first", na.last = "keep")` のショートカット。
    グループ毎に連番を振るのに便利。

`dplyr::ntile(x, n)`
:   数値ベクトル `x` を順位によって `n` 個のクラスに均等分け

`dplyr::n_distinct(x)`
:   高速で簡潔な `length(unique(x))`

`dplyr::n_groups(x)`
:   グループ数。

`dplyr::last(x, order_by=NULL, default=default_missing(x))`
:   最後の要素にアクセス。
    `x[length(x)]` や `tail(x, 1)` よりも楽チンだが、
    安全性重視の `dplyr::nth()` を内部で呼び出すため遅い。

`dplyr::lead(x n=1, default=NA, order_by=NULL)`, `dplyr::lag(...)`
:   `x` の中身を `n` だけずらして `default` で埋める。
    `lead()` は前に、`lag()` は後ろにずらす

    ```r
    lag(seq_len(5), 2)
    ## [1] NA NA 1 2 3
    ```

`dplyr::between(x, left, right)`
:   `left <= x & x <= right` のショートカット。

`dplyr::near(x, y, tol=.Machine$double.eps^0.5)`
:   `abs(x - y) < tol` のショートカット。

`dplyr::case_when(...)`
:   `if {} else if {} else if {} ...` のショートカット。


## グループ化

[`tidyr`]({{< relref "tidyr.md" >}}) でネストして、
[`purrr`]({{< relref "purrr.md" >}}) でその list of data.frames に処理を施し、
[`dplyr`]({{< relref "dplyr.md" >}}) でその変更を元の data.frame に適用する、
というのがtidyverse流のモダンなやり方らしい。

```r
iris %>%
  dplyr::group_nest(Species) %>%
  dplyr::mutate(data= purrr::map(data, head)) %>%
  tidyr::unnest()
```

`dplyr::group_by(.data, ..., add=FALSE)`
:   グループごとに区切って次の処理に渡す。
    e.g. `summarise()`, `tally()` など

`dplyr::group_data(.data)`
:   グループ情報を参照:
    ```r
    iris %>% group_by(Species) %>% group_data()
    #       Species                       .rows
    #        <fctr>                      <list>
    # 1:     setosa             1,2,3,4,5,6,...
    # 2: versicolor       51,52,53,54,55,56,...
    # 3:  virginica 101,102,103,104,105,106,...
    ```

    左側のキー列だけ欲しければ `dplyr::group_keys()` 、<br>
    右端の行番号だけ欲しければ `dplyr::group_rows()` 。

`dplyr::group_nest(.tbl, ..., .key = "data", keep = FALSE)`
:   入れ子 data.frame を作る。
    `group_by(...) %>% tidyr::nest()` のショートカット。

`dplyr::group_split(.tbl, ..., keep = FALSE)`
:   list of data.frames に分割する。

`dplyr::group_indices(.data, ...)`
:   `grouped_df` ではなくグループIDとして1からの整数列を返す版 `group_by()`

`dplyr::rowwise(.data)`
:   行ごとに区切って次の処理に渡す。

`dplyr::do(.data, ...)`
:   **廃止予定？** グループごとに処理する。
    `{}` 内に長い処理を書いてもいいし、関数に渡してもよい。
    グループごとに切りだされた部分は `.` で参照できる。
    出力がdata.frameじゃないと
    `Error: Results are not data frames at positions: 1`
    のように怒られるが、
    `do(dummy=func(.))` のように名前付きにすると
    data.frameに入らないような型でも大丈夫になる。

    ```r
    iris %>%
      dplyr::group_by(Species) %>%
      dplyr::do(head(.))
    ```


## 関連書籍

<a href="https://www.amazon.co.jp/dp/4774198536/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=f0acaf09c5bcbd85ee22c534caede9d1" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4774198536&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4774198536" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?s=books&ie=UTF8&qid=1508340700&sr=1-3&linkCode=li3&tag=heavywatal-22&linkId=6a53371cc80e2d6d7fc50fae5d8b862d" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%A7%E3%81%AF%E3%81%98%E3%82%81%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9-Hadley-Wickham/dp/487311814X/ref=as_li_ss_il?ie=UTF8&qid=1508340144&sr=8-1&keywords=r&linkCode=li3&tag=heavywatal-22&linkId=4137d3d3351f8ccab5a93cefdc28fdec" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
