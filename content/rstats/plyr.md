+++
title = 'plyr'
subtitle = "データ分割-関数適用-再結合を効率的に"
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = 99
+++

Split, Apply, Combine
:   特定の条件でデータを分割し、それぞれに関数を適用し、再びそれを統合する。
    R標準 `apply` 系の関数よりも直感的な使い方ができ、処理も高速。

Rの中から `install.packages('plyr')` でインストールし、
使う前に `library(plyr)` で読み込む。

{{%div class="warning"%}}
data.frame 処理には [dplyr]({{< relref "dplyr.md" >}})、
list, vector 処理には [purrr]({{< relref "purrr.md" >}})
がある今となっては、このパッケージを使うべき場面はもう無さそう。
{{%/div%}}

-   <http://plyr.had.co.nz/>
-   <http://www.jstatsoft.org/v40/i01>
-   <https://www.rdocumentation.org/packages/plyr>
-   <https://cran.r-project.org/web/packages/plyr/index.html>

## コア関数

    from                           to array  da.fr  list   nothing
    array                             aaply  adply  alply  a_ply
    data.frame                        daply  ddply  dlply  d_ply
    list or vector                    laply  ldply  llply  l_ply
    (Replicates evaluation)           raply  rdply  rlply  r_ply
    (Call a multi-argument function)  maply  mdply  mlply  m_ply

e.g. 複数ファイルを読み込んでひとつのdata.frameにまとめる。
listからdata.frameを作るので `ldply()`

```r
> filenames = list.files(pattern='\\.csv$')
> large_table = ldply(filenames, read.csv)
```

e.g. data.frameについてある列の値でグループ化し、
グループ毎に数値の列の平均を取る。
data.frameからdata.frameを作るので `ddply()`

```r
> ddply(iris, .(Species), numcolwise(mean))
     Species Sepal.Length Sepal.Width Petal.Length Petal.Width
1     setosa        5.006       3.428        1.462       0.246
2 versicolor        5.936       2.770        4.260       1.326
3  virginica        6.588       2.974        5.552       2.026
```

## ヘルパー関数

`plyr::join(x, y, by=NULL, type="left", match="all")`
:   `by` で指定した列の値が等しいものを同じ行として、いい感じに `cbind()`。
    複数行がマッチした場合のデフォルトの挙動は `base::merge()` と同じく
    `match="all"` だが `match="first"` も指定できて、そちらは高速らしい。

    `type=`
    :   `"inner"`: `x` と `y` の `by` がマッチする行のみ\
        `"left"`: `x` の全ての行を保持\
        `"right"`: `y` の全ての行を保持\
        `"full"`: `"left"` の結果の下に、`y` の残りの行を追加

`plyr::join_all(dfs, by=NULL, type="left", match="all")`
:   listに入った複数のdata.frameを再帰的に `join()` する。

`plyr::rename(x, replace)`
:   data.frame列名などを **部分的に** 変更

    ```r
    # replace引数には名前付き文字列vectorを与える
    # 古い名前が名前、新しい名前が値
    plyr::rename(.data, c(col_a = "alpha", col_b = "beta"))
    ```

`plyr::count(.data, vars=NULL, wt_var=NULL)`
:   data.frameのなかで `vars` 列に関してユニークな行数をカウント。
    重み付けに使う列を `wt_var` に指定できる。

`plyr::colwise(.fun, .cols=true, ...)`,
:   関数を列ごとに適用するものに変換する。
    例えば `colwise(mean)(.data)` は `colMeans(.data)` とほぼ同義。
    関数で使えない型が含まれている行の結果には `NA` が入る。
    `numcolwise(.fun, ...)` と `catcolwisw(.fun, ...)`
    はそれぞれ数値の行、カテゴリ変数の行だけに適用する関数を返してくれる。

`plyr::each(func1, func2, ...)`
:   同じ引数に対して複数の関数を並列に作用させる。 e.g. `each(min, max, mean)(1:10)`, `each(head, tail)(.data, n=10)`

`plyr::splat(func)`
:   ひとつのリストや文字列ベクタでまとめて引数を受け取れるような関数に変換する。
    `do.call()` はlistしか取らないがこちらは名前付きベクタも可

    ```r
    > params = c(by=2, length=4)
    > splat(seq)(params)
    [1] 1 3 5 7
    > do.call(seq, as.list(params))
    [1] 1 3 5 7
    ```

## 並列化

`doMC` 越しに `foreach` をバックエンドとして使用する

```r
install.packages("doMC")
library(doMC)
doMC::registerDoMC(parallel::detectCores())

.data = plyr::ldply(lst, func, .parallel=TRUE)
```
