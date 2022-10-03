+++
title = 'dplyr'
subtitle = "高速data.frame処理"
tags = ["r", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -75
+++

<a href="https://dplyr.tidyverse.org/">
<img src="/_img/hex-stickers/dplyr.webp" align="right" width="120" height="139">
</a>

data.frameに対して抽出(select, filter)、部分的変更(mutate)、要約(summarize)、ソート(arrange)などの処理を施すためのパッケージ。
前作 [plyr]({{< relref "plyr.md" >}}) のうちdata.frameに関する部分が強化されている。
[purrr]({{< relref "purrr.md" >}}) や [tidyr]({{< relref "tidyr.md" >}}) と一緒に使うとよい。

[tidyverse](https://tidyverse.tidyverse.org/) に含まれているので、
`install.packages("tidyverse")` で一括インストール、
`library(tidyverse)` で一括ロード。

-   <https://r4ds.had.co.nz/transform.html>
-   <https://github.com/tidyverse/dplyr>

## 関数の連結 %>%

<a href="https://magrittr.tidyverse.org/">
<img src="/_img/hex-stickers/pipe.webp" align="right" width="120" height="139">
</a>

dplyrではなく[magrittr](https://magrittr.tidyverse.org/)の機能。
R 4.1 から標準でも `|>` として使えるようになった。

`x %>% f(a, b)` は `f(x, a, b)` と等価。
左の値 `x` を第一引数として右の関数 `f()` に渡す。
一時変数を作ったり、関数を何重にも重ねたりすることなく、
適用する順に次々と処理を記述することができるようになる。
慣れれば書きやすく読みやすい。

```r
library(tidyverse)

## with piping
result = diamonds %>%              # 生データから出発して
  select(carat, cut, price) %>%    # 列を抽出して
  filter(carat > 1) %>%            # 行を抽出して
  group_by(cut) %>%                # グループ化して
  summarize(mean(price)) %>%       # 平均を計算
  print()                          # 表示してみる

## with a temporary variable
result = select(diamonds, carat, cut, price) # 列を抽出して
result = filter(result, carat > 1)           # 行を抽出して
result = group_by(result, cut)               # グループ化して
result = summarize(result, mean(price))      # 平均を計算

## with nested functions
result = summarize(                    # 平均を計算
    group_by(                            # グループ化して
      filter(                              # 行を抽出して
        select(diamonds, carat, cut, price), # 列を抽出して
        carat > 1),                        # 行を抽出して
      cut),                              # グループ化して
    mean(price))                       # 平均を計算

## result
        cut mean(price)
      <ord>       <dbl>
1      Fair    7177.856
2      Good    7753.601
3 Very Good    8340.549
4   Premium    8487.249
5     Ideal    8674.227
```

現状では `magrittr` パッケージの `%>%` が広く採用されているが、
`pipeR` パッケージの `%>>%` のほうが高速らしい。
<https://renkun-ken.github.io/pipeR-tutorial/>


## 抽出・絞り込み

### 列

[`dplyr::select(.data, ...)`](https://dplyr.tidyverse.org/reference/select.html)
:   列を絞る。`:`範囲指定、`!`負の指定が可能。
    [selection helpers](https://tidyselect.r-lib.org/reference/language.html)
    によるパターン指定も便利。
    複数の条件を組み合わせるには `&` (AND), `|` (OR) で。
    残るのが1列だけでも勝手にvectorにはならずdata.frameのまま。

    ```r
    diamonds %>% dplyr::select(1, 2, 7)
    diamonds %>% dplyr::select(carat, cut, price)
    diamonds %>% dplyr::select(c("carat", "cut", "price"))
    diamonds %>% dplyr::select(!c(carat, cut, price))
    diamonds %>% dplyr::select(starts_with("c"))
    diamonds %>% dplyr::select(where(is.numeric))
    # diamonds %>% dplyr::select(-carat, -cut, -price)
    ```

    カンマ区切りはORの意味で働くはずだがマイナス指定
    `dplyr::select(-carat, -cut, -price)`
    の場合だけANDのように働く(つまり3列とも抜ける)ので特殊。
    (だからマニュアルから消えたのかな？)

:   文字列ベクタで指定しようとすると意図が曖昧になるので、
    [ヘルパー関数 `all_of()`, `any_of()`](https://tidyselect.r-lib.org/reference/all_of.html)を挟む。あるいは
    [{{embrace}}](https://rlang.r-lib.org/reference/embrace-operator.html),
    [!!unquoting](https://adv-r.hadley.nz/quasiquotation.html),
    [pronoun$](https://rlang.r-lib.org/reference/dot-data.html)を使う:
    ```r
    clarity = c("carat", "cut", "price")
    diamonds %>% dplyr::select(clarity)         # ambiguous!
    diamonds %>% dplyr::select(.data$clarity)   # clarity
    diamonds %>% dplyr::select(all_of(clarity)) # carat, cut, price
    diamonds %>% dplyr::select(any_of(clarity)) # carat, cut, price
    diamonds %>% dplyr::select({{clarity}})     # carat, cut, price
    diamonds %>% dplyr::select(!!clarity)       # carat, cut, price
    diamonds %>% dplyr::select(!!!rlang::syms(clarity))  # carat, cut, price
    ```
    これらの指定方法は `rename()` や `pull()` でも有効。

    一方、文字列を受け取れない `distinct()` や `group_by()`
    などの関数には普通のunquoteは通用しない。
    最後の例のように `rlang::data_syms()` でシンボルのリストを作って
    [!!!unquote-splicing](https://rlang.r-lib.org/reference/splice-operator.html)
    して渡す必要がある。
    ```r
    columns = c("cut", "color")
    diamonds %>% distinct(!!as.name(columns[1L]))
    diamonds %>% distinct(!!!rlang::data_syms(columns))
    ```
    詳しくは [rlangパッケージのドキュメント](https://rlang.r-lib.org/) か
    [宇宙船本](https://amzn.to/3KpiXq0)第3章のコラム
    「selectのセマンティクスとmutateのセマンティクス」を参照。


`dplyr::rename(.data, ...)`
:   列の改名。
    `mutate()`と同じようなイメージで `new = old` と指定。
    ```r
    diamonds %>% dplyr::rename(SIZE = carat)
    #    SIZE       cut color clarity depth table price     x     y     z
    # 1  0.23     Ideal     E     SI2  61.5    55   326  3.95  3.98  2.43
    # 2  0.21   Premium     E     SI1  59.8    61   326  3.89  3.84  2.31
    ```
:   変数に入った文字列を使う場合も`mutate()`と同様に {{embrace}} や !!unquote で:
    ```r
    old_name = "carat"
    new_name = toupper(old_name)
    diamonds %>% dplyr::rename({{new_name}} := {{old_name}})
    # tbl_df [53940 x 10]
    #   CARAT       cut color clarity depth table price     x     y     z
    # 1  0.23     Ideal     E     SI2  61.5    55   326  3.95  3.98  2.43
    # 2  0.21   Premium     E     SI1  59.8    61   326  3.89  3.84  2.31
    ```
    名前付きベクターと
    [!!!unquote-splicing](https://rlang.r-lib.org/reference/splice-operator.html)
    を使えば一括指定できる:
    ```r
    named_vec = setNames(names(diamonds), LETTERS[seq(1, 10)])
    diamonds %>% dplyr::rename(!!!named_vec)
    #       A         B     C     D     E     F     G     H     I     J
    #   <dbl>     <ord> <ord> <ord> <dbl> <dbl> <int> <dbl> <dbl> <dbl>
    # 1  0.23     Ideal     E   SI2  61.5    55   326  3.95  3.98  2.43
    # 2  0.21   Premium     E   SI1  59.8    61   326  3.89  3.84  2.31
    ```
:   `rename_with(.tbl, .fn, .cols = everything(), ...)`
    はリネーム関数を渡せる亜種:
    ```r
    diamonds %>% dplyr::rename_with(toupper, everything())
    ```

`dplyr::pull(.data, var = -1)`
:   指定した1列をvector(またはlist)としてdata.frameから抜き出す。

    ```r
    diamonds %>% head() %>% dplyr::pull(price)
    diamonds %>% head() %>% dplyr::pull("price")
    diamonds %>% head() %>% dplyr::pull(7)
    diamonds %>% head() %>% dplyr::pull(-4)
    diamonds %>% head() %>% `[[`("price")
    diamonds %>% head() %>% {.[["price"]]}
    ```


### 行

`dplyr::filter(.data, ...)`
:   条件を満たす行だけを返す。`base::subset()` と似たようなもの。

    ```r
    diamonds %>% dplyr::filter(carat > 3 & price < 10000)
    #   carat     cut color clarity depth table price     x     y     z
    # 1  3.01 Premium     I      I1  62.7    58  8040  9.10  8.97  5.67
    # 2  3.11    Fair     J      I1  65.9    57  9823  9.15  9.02  5.98
    # 3  3.01 Premium     F      I1  62.2    56  9925  9.24  9.13  5.73
    ```

    {{<div class="warning">}}
評価結果が `NA` となる行は除去される。
特に不等号を使うときやや直感に反するので要注意。
e.g., `filter(gene != "TP53")`
{{</div>}}

:   複数列で条件指定するには `if_any()`, `if_all()` が使える。

    ```r
    diamonds %>% dplyr::filter(if_any(c(x, y, z), ~ .x > 20))
    #   carat       cut color clarity depth table price     x     y     z
    #   <dbl>     <ord> <ord>   <ord> <dbl> <dbl> <int> <dbl> <dbl> <dbl>
    # 1  2.00   Premium     H     SI2  58.9  57.0 12210  8.09 58.90  8.06
    # 2  0.51 Very Good     E     VS1  61.8  54.7  1970  5.12  5.15 31.80
    # 3  0.51     Ideal     E     VS1  61.8  55.0  2075  5.15 31.80  5.12
    diamonds %>% dplyr::filter(if_all(where(is.numeric), ~ .x > 4))
    #   carat     cut color clarity depth table price     x     y     z
    #   <dbl>   <ord> <ord>   <ord> <dbl> <dbl> <int> <dbl> <dbl> <dbl>
    # 1  4.01 Premium     I      I1  61.0    61 15223 10.14 10.10  6.17
    # 2  4.01 Premium     J      I1  62.5    62 15223 10.02  9.94  6.24
    # 3  4.13    Fair     H      I1  64.8    61 17329 10.00  9.85  6.43
    # 4  5.01    Fair     J      I1  65.5    59 18018 10.74 10.54  6.98
    # 5  4.50    Fair     J      I1  65.8    58 18531 10.23 10.16  6.72
    ```


`dplyr::distinct(.data, ..., .keep_all = FALSE)`
:   指定した列に関してユニークな行のみ返す。
    `base::unique.data.frame()` よりも高速で、
    `filter(!duplicated(.[, ...]))` よりスマートで柔軟。
    指定しなかった列を残すには `.keep_all = TRUE` とする。

    ```r
    diamonds %>% dplyr::distinct(cut)
    #         cut
    # 1     Ideal
    # 2   Premium
    # 3      Good
    # 4 Very Good
    # 5      Fair
    ```

`dplyr::slice(.data, ...)`
:   行番号を指定して行を絞る。
    `` `[`(i,) `` の代わりに。

    ```r
    diamonds %>% dplyr::slice(1, 2, 3)
    #   carat     cut color clarity depth table price     x     y     z
    # 1  0.23   Ideal     E     SI2  61.5    55   326  3.95  3.98  2.43
    # 2  0.21 Premium     E     SI1  59.8    61   326  3.89  3.84  2.31
    # 3  0.23    Good     E     VS1  56.9    65   327  4.05  4.07  2.31
    ```

`dplyr::slice_head(.data, ..., n, prop)`, `slice_tail()`
:   先頭・末尾の行を抽出。

`dplyr::slice_sample(.data, ..., n, prop, weight_by = NULL, replace = FALSE)`
:   指定した行数・割合だけランダムサンプルする。
:   `sample_n()` と `sample_frac()` は非推奨。


## 列の変更・追加

`dplyr::mutate(.data, ...)`
:   既存の列を変更したり、新しい列を作ったり。
    `base::transform()` の改良版。
    ```r
    # modify existing column
    diamonds %>% dplyr::mutate(price = price * 107.54)

    # create new column
    diamonds %>% dplyr::mutate(gram = 0.2 * carat)
    ```

    変数に入った文字列を列名として使いたい場合は上記 `select()` のときと同様、
    {{embrace}} や !!unquote を使う。
    左辺(代入先)の名前も文字列を使って表現したい場合は
    [walrus(セイウチ)演算子 `:=`](https://www.tidyverse.org/blog/2020/02/glue-strings-and-tidy-eval/#custom-result-names)
    を使う:
    ```r
    A = "carat"
    B = "gram"

    # unquoting only right hand side
    diamonds %>% dplyr::mutate(B = 0.2 * !!as.name(A))
    #   carat     cut color clarity depth table price     x     y     z     B
    # 1  0.23   Ideal     E     SI2  61.5    55   326  3.95  3.98  2.43 0.046
    # 2  0.21 Premium     E     SI1  59.8    61   326  3.89  3.84  2.31 0.042

    # unquoting both sides
    diamonds %>% dplyr::mutate({{B}} := 0.2 * !!as.name(A))
    #   carat     cut color clarity depth table price     x     y     z  gram
    # 1  0.23   Ideal     E     SI2  61.5    55   326  3.95  3.98  2.43 0.046
    # 2  0.21 Premium     E     SI1  59.8    61   326  3.89  3.84  2.31 0.042
    ```

`dplyr::transmute(.data, ...)`
:   指定した列以外を保持しない版 `mutate()` 。
    言い換えると、列の中身の変更もできる版 `select()` 。


## data.frameの要約・集計・整列

`dplyr::summarize(.data, ..., .group = NULL)`
:   指定した列に関数を適用して1行のdata.frameにまとめる。
    グループ化されていたらグループごとに適用して `bind_rows()` する。

    ```r
    diamonds %>%
      group_by(cut) %>%
      summarize(avg_carat = mean(carat), max_price = max(price))
    #         cut avg_carat max_price
    #       <ord>     <dbl>     <int>
    # 1      Fair 1.0461366     18574
    # 2      Good 0.8491847     18788
    # 3 Very Good 0.8063814     18818
    # 4   Premium 0.8919549     18823
    # 5     Ideal 0.7028370     18806
    ```

:   関数の結果が長さ1じゃなくても大丈夫(v1.0.0):

    ```r
    diamonds %>% dplyr::summarize(range(carat), range(price))
    #   range(carat) range(price)
    #          <dbl>        <int>
    # 1         0.20          326
    # 2         5.01        18823
    ```

:   各グループの結果がtibbleの場合、結合・展開してくれる:

    ```r
    diamonds %>% dplyr::nest_by(cut) %>%
      dplyr::summarize(head(data, 2L))
    ```

    これを利用して複数ファイルを一気読みできる:

    ```r
    path = fs::dir_ls("path/to/data", glob = "*.tsv")
    df = tibble(path) %>%
      dplyr::rowwise() %>%
      dplyr::summarize(readr::read_tsv(path))
    ```

    けどまあ `purrr::map_dfr(path, read_tsv)` のほうが読みやすい気がする。

:   複数カラムに関数を適用するには `across()` を使う:

    ```r
    median.ordered = function(x, na.rm = FALSE) {
      levels(x)[median(as.integer(x), na.rm = na.rm)]
    }

    # summarize_at
    diamonds %>% group_by(cut) %>%
      summarize(across(starts_with("c"), median))

    # summarize_if
    diamonds %>% group_by(cut) %>%
      summarize(across(where(is.numeric), mean, na.rm = TRUE))

    # summarize_all
    diamonds %>% group_by(cut) %>%
      summarize(across(everything(), median, na.rm = TRUE))
    ```

    `*_at()`, `*_if()`, `*_all()`, `*_each()` は非推奨。

`dplyr::tally(x, wt, sort = FALSE)`
:   `summarize(x, n = n())` のショートカット。
    `wt` にカラムを指定して重み付けすることもできる。
:   `dplyr::add_tally()` は元の形を維持したままカウント列を追加。

`dplyr::count(x, ..., wt = NULL, sort = FALSE)`
:   `group_by(...) %>% tally()` のショートカット。
:   `dplyr::add_count()` は元の形を維持したままカウント列を追加。

`dplyr::arrange(.data, column1, column2, ...)`
:   指定した列の昇順でdata.frameの行を並べ替える。
    `arrange(desc(column))` で降順になる。
    `order()` を使うよりもタイピングの繰り返しが少ないし直感的
    ```r
    mtcars[order(mtcars$cyl, mtcars$disp), ]
    # is equivalent to
    mtcars %>% dplyr::arrange(cyl, disp)
    ```

`dplyr::relocate(.data, ..., .before = NULL, .after = NULL)`
:   指定した列を左端(もしくは `.before`/`.after` 別の列)に移動する。

`dplyr::slice_max(.data, order_by, ..., n, prop, with_ties = TRUE)`
:   `order_by, ...` で指定した列の降順で `n` 行だけ返す。
    境界のタイを保持する場合は `n` 行以上の出力になる。
    昇順で取りたいときはマイナス指定か `slice_min()`。
:   `top_n()` は非推奨。


## data.frameを結合

`dplyr::***_join(x, y, by = NULL, copy = FALSE)`
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

`dplyr::if_else(condition, true, false, missing = NULL)`
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

`dplyr::recode(.x, ..., .default = NULL, .missing = NULL)`
:   vectorの値を変更する。e.g.,
    `recode(letters, a = "A!", c = "C!")`

`dplyr::row_number(x)`
:   `rank(x, ties.method = "first", na.last = "keep")` のショートカット。
    グループ毎に連番を振るのに便利。

`dplyr::ntile(x, n)`
:   数値ベクトル `x` を順位によって `n` 個のクラスに均等分け

`dplyr::n_distinct(x)`
:   高速で簡潔な `length(unique(x))`

`dplyr::n_groups(x)`
:   グループ数。

`dplyr::last(x, order_by = NULL, default = default_missing(x))`
:   最後の要素にアクセス。
    `x[length(x)]` や `tail(x, 1)` よりも楽チンだが、
    安全性重視の `dplyr::nth()` を内部で呼び出すため遅い。

`dplyr::lead(x n = 1, default = NA, order_by = NULL)`, `dplyr::lag(...)`
:   `x` の中身を `n` だけずらして `default` で埋める。
    `lead()` は前に、`lag()` は後ろにずらす

    ```r
    lag(seq_len(5), 2)
    ## [1] NA NA 1 2 3
    ```

`dplyr::between(x, left, right)`
:   `left <= x & x <= right` のショートカット。

`dplyr::near(x, y, tol = .Machine$double.eps^0.5)`
:   `abs(x - y) < tol` のショートカット。

`dplyr::case_when(...)`
:   `if {} else if {} else if {} ...` のショートカット。


## グループ化

[`tidyr`]({{< relref "tidyr.md" >}}) でネストして、
[`purrr`]({{< relref "purrr.md" >}}) でその list of data.frames に処理を施し、
[`dplyr`]({{< relref "dplyr.md" >}}) でその変更を元の data.frame に適用する、
というのがtidyverse流のモダンなやり方らしい。
それをもうちょいスマートにやる `group_modify()` がv0.8.1で導入された。

```r
diamonds %>%
  dplyr::group_nest(cut) %>%
  dplyr::mutate(data = purrr::map(data, head, n = 2L)) %>%
  tidyr::unnest()

# since dplyr 0.8.1
diamonds %>%
  dplyr::group_by(cut) %>%
  dplyr::group_modify(~ head(.x, 2L))
```

`dplyr::group_by(.data, ..., add = FALSE, .drop = group_by_drop_default(.data))`
:   グループごとに区切って次の処理に渡す。
    e.g. `summarize()`, `slice()`, `tally()`, `group_modify()` など
:   `.drop = FALSE` とすると行数ゼロになるグループも捨てずに保持できる。

`dplyr::group_data(.data)`
:   グループ情報を参照:
    ```r
    diamonds %>% dplyr::group_by(cut) %>% dplyr::group_data()
    #         cut           .rows
    # 1      Fair    <int [1610]>
    # 2      Good    <int [4906]>
    # 3 Very Good   <int [12082]>
    # 4   Premium   <int [13791]>
    # 5     Ideal   <int [21551]>
    ```

    左側のキー列だけ欲しければ `dplyr::group_keys()` 、<br>
    左端の行番号だけ欲しければ `dplyr::group_rows()` 。

`dplyr::group_nest(.tbl, ..., .key = "data", keep = FALSE)`
:   入れ子 data.frame を作る。
    `group_by(...) %>% tidyr::nest()` のショートカット。

`dplyr::group_split(.tbl, ..., keep = FALSE)`
:   list of data.frames に分割する。
    `.tbl %>% split(.$group)` と同等だが、
    ドットダラーを使わくて済むし複数列をキーにするのも簡単。

`dplyr::group_indices(.data, ...)`
:   `grouped_df` ではなくグループIDとして1からの整数列を返す版 `group_by()`

`dplyr::group_modify(.tbl, .f, ...)`
:   グループごとに `.f` を適用して再結合したdata.frameを返す。
:   `.f` は2つの引数をとる関数。
    1つめは区切られたdata.frame、
    2つめはグループ化のキー(となった列を含む1行のdata.frame)。
:   引数が1つだったり、違うものを2つめに受け取る関数はそのままじゃ渡せないので、
    [purrrでよく見る無名関数]({{< relref "purrr.md#無名関数" >}})
    に包んで渡す。
    このとき2つの引数はそれぞれ `.x`, `.y` として参照できる。

    ```r
    diamonds %>% dplyr::group_by(cut) %>% dplyr::group_modify(~ head(.x, 2L))
    #          cut carat color clarity depth table price     x     y     z
    #  1      Fair  0.22     E     VS2  65.1    61   337  3.87  3.78  2.49
    #  2      Fair  0.86     E     SI2  55.1    69  2757  6.45  6.33  3.52
    #  3      Good  0.23     E     VS1  56.9    65   327  4.05  4.07  2.31
    #  4      Good  0.31     J     SI2  63.3    58   335  4.34  4.35  2.75
    #  5 Very Good  0.24     J    VVS2  62.8    57   336  3.94  3.96  2.48
    #  6 Very Good  0.24     I    VVS1  62.3    57   336  3.95  3.98  2.47
    #  7   Premium  0.21     E     SI1  59.8    61   326  3.89  3.84  2.31
    #  8   Premium  0.29     I     VS2  62.4    58   334  4.20  4.23  2.63
    #  9     Ideal  0.23     E     SI2  61.5    55   326  3.95  3.98  2.43
    # 10     Ideal  0.23     J     VS1  62.8    56   340  3.93  3.90  2.46
    ```

:   `group_map()` は結果を `bind_rows()` せずlistとして返す亜種。
    `group_walk()` は `.f` 適用前の `.tbl` を返す亜種。

`dplyr::rowwise(data, ...)`
:   受け取ったdataに `rowwise_df` クラスを付与して返す。
    これは1行ごとにグループ化された `grouped_df` のようなもので、
    mutate などを適用すると列全体ではなく1行ごとに関数に渡される。
    一旦非推奨となったがv1.0.0で蘇った。
:   ループ処理が重いので、計算内容と行数のバランスに注意。
    例えば `dplyr::c_across()` の併用で横方向の合計や平均を簡潔に書けるが、
    ベクトル演算を使う場合に比べてかなり遅い。

    ```r
    diamonds %>%
      dplyr::rowwise() %>%
      dplyr::mutate(mean = mean(c_across(x:z)))         # slow

    diamonds %>% dplyr::mutate(mean = (x + y + z) / 3)  # fast
    ```

`dplyr::do(.data, ...)`
:   非推奨。
    代わりに `group_modify()` とかを使う。

    ```r
    diamonds %>% dplyr::group_by(cut) %>% dplyr::do(head(., 2L))
    ```


## matrix, array

data.frame を主眼とする dplyr では matrix や array を扱わない。
一時期存在していた `tbl_cube` 関連の機能は
[cubelyr](https://github.com/hadley/cubelyr) に隔離された。

`as.tbl_cube(x, dim_names, met_names, ...)`
:   matrix/arrayからdata.frameの一歩手前に変換する。
    [`reshape2::melt`]({{< relref "reshape2.md" >}}) の改良版。
    これの結果に `tibble::as_tibble()` を適用するとわかりやすい。
:   ただし `dimnames(x)` が空ではダメで、長さの正しい名前付きlistになっている必要がある。

```r
# deprecated
iris3 %>%
  reshape2::melt() %>%
  tibble::as_tibble() %>%
  dplyr::rename(obs = Var1, metrics = Var2, species = Var3)

# new
x = iris3
dimnames(x)[[1L]] = seq_len(dim(iris3)[[1L]])
names(dimnames(x)) = c("obs", "metrics", "species")
cubelyr::as.tbl_cube(x, met_name = "value") %>% as_tibble()
```


## 関連書籍

<a href="https://www.amazon.co.jp/dp/4297121700?&linkCode=li3&tag=heavywatal-22&linkId=77762a4d0080a840ec5d94df9c0c5ceb&language=ja_JP&ref_=as_li_ss_il" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4297121700&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4297121700" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?s=books&ie=UTF8&qid=1508340700&sr=1-3&linkCode=li3&tag=heavywatal-22&linkId=6a53371cc80e2d6d7fc50fae5d8b862d" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%A7%E3%81%AF%E3%81%98%E3%82%81%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9-Hadley-Wickham/dp/487311814X/ref=as_li_ss_il?ie=UTF8&qid=1508340144&sr=8-1&keywords=r&linkCode=li3&tag=heavywatal-22&linkId=4137d3d3351f8ccab5a93cefdc28fdec" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
