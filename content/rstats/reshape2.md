+++
title = 'reshape2'
subtitle = "柔軟なデータ変形ツール"
tags = ["r", "hadley"]
[menu.main]
  parent = "rstats"
  weight = 99
+++

{{%div class="warning"%}}
`reshape2` はもう古いので、このページの内容も更新しない。
同じ作者が新しく設計しなおした
[tidyr]({{< relref "tidyr.md" >}}) + [dplyr]({{< relref "dplyr.md" >}})
のほうがより高速で洗練されているのでそちらを使おう。
{{%/div%}}

-   <http://had.co.nz/reshape/>
-   <http://cran.r-project.org/web/packages/reshape2/>
-   <http://www.rdocumentation.org/packages/reshape2>
-   <http://seananderson.ca/2013/10/19/reshape.html>

## `melt()`

`data.frame` の複数列の値を、カテゴリ変数1列と値1列の組に変換する。
これにより、変換する列数の分だけ `data.frame` が縦長(long-format)になる。
やや冗長性は増すが、[ggplot2]({{< relref "ggplot2.md" >}}) での作図などさまざまな操作がしやすくなる。

{{%div class="note"%}}
これじゃなくて[tidyr]({{< relref "tidyr.md" >}})の`gather()`を使おう。
{{%/div%}}

```r
reshape2::melt(data, id.vars, measure.vars,
               variable.name="variable", value.name="value",
               na.rm=FALSE, factorsAsStrings=TRUE, ...)
```

`data`
:   `data.frame`

`id.vars`
:   そのまま列として維持したい列名を文字列で指定。
    何も指定しなければ `measure.vars` 以外のすべて。

`measure.vars`
:   列名を `variable` に、値を `value` に分解したい列名を文字列で指定。
    何も指定しなければ `id.vars` 以外のすべて。

`variable.name="variable"`
:   meltされた列名を格納する新しい列の名前。

`value.name="value"`
:   meltされた値を格納する新しい列の名前。

`na.rm=FALSE`
:   `NA` が含まれる行を取り除くかどうか。

### Example

1. ライブラリを読み込んでサンプルデータを見てみる(wide-format)
    ```r
    > library(reshape2)
    > head(reshape2::french_fries)
       time treatment subject rep potato buttery grassy rancid painty
    61    1         1       3   1    2.9     0.0    0.0    0.0    5.5
    25    1         1       3   2   14.0     0.0    0.0    1.1    0.0
    62    1         1      10   1   11.0     6.4    0.0    0.0    0.0
    26    1         1      10   2    9.9     5.9    2.9    2.2    0.0
    63    1         1      15   1    1.2     0.1    0.0    1.1    5.1
    27    1         1      15   2    8.8     3.0    3.6    1.5    2.3
    ```

2. データをlong-formatに整形
   ```r
   > molten = reshape2::melt(reshape2::french_fries,
                             id.vars=c("time", "treatment", "subject", "rep"),
                             variable.name="flavor", na.rm=TRUE)
   > head(molten)
     time treatment subject rep flavor value
   1    1         1       3   1 potato   2.9
   2    1         1       3   2 potato  14.0
   3    1         1      10   1 potato  11.0
   4    1         1      10   2 potato   9.9
   5    1         1      15   1 potato   1.2
   6    1         1      15   2 potato   8.8
   ```

3.  [ggplot2]({{< relref "ggplot2.md" >}}) で作図

    ```r
    > library(ggplot2)
    > gp = ggplot(molten, aes(x=time, y=value, colour=treatment, shape=as.factor(rep)))
    > gp = gp + geom_point(alpha=0.3)
    > gp = gp + geom_smooth(aes(group=treatment), method=loess, se=FALSE)
    > gp = gp + facet_grid(flavor ~ subject)
    > gp
    ```

## `dcast()`, `acast()`

カテゴリ変数を含む `data.frame` を `melt()` と逆方向に
(long-formatからwide-formatへ)整形する。
2次元までの `data.frame` が欲しければ `dcast()` 、
3次元以上の `Array` が欲しければ `acast()` を使う。

{{%div class="note"%}}
これじゃなくて[tidyr]({{< relref "tidyr.md" >}})の`spread()`を使おう。
`fun.aggregate`のように関数をグループごとに適用したい場合は
[dplyr]({{< relref "dplyr.md" >}})の`group_by()`と`summarise()`を使う。
{{%/div%}}

```r
reshape2::dcast(data, formula, fun.aggregate=NULL, ...,
          margins=NULL, subset=NULL, fill=NULL, drop=TRUE,
          value.var=guess_value(data))
```

`data`
:   `melt()` されたような形でカテゴリ変数を含む `data.frame`

`formula`
:   `x_var ~ y_var ~ z_var ~ ...` のような形で出力形式を指定

`fun.aggregate=NULL`
:   `mean` や `sum` など、整形後に同じマスに来る複数の値に適用する関数。
    デフォルトでは `length` が働いて要素数が得られる。

`...`
:   aggregate関数への引数を渡せる

`margins=NULL`
:   列全体の平均や行全体の和などを追加するかどうか

`subset=NULL`
:   適用範囲を限定する e.g., `subset=.(variable=="length")`

`fill=NULL`

`drop=TRUE`

`value.var=guess_value(data)`

### Example

データは上の `melt()` の例で作った `molten`。

`fun.aggregate` を省略すると `length` が適用されて要素数が分かる

```r
> reshape2::acast(molten, treatment ~ flavor)
Aggregation function missing: defaulting to length
  potato buttery grassy rancid painty
1    232     231    232    232    232
2    232     230    232    232    231
3    231     231    231    231    231
```

グループごとの平均値を `data.frame` で

```r
> reshape2::dcast(molten, treatment ~ flavor, mean)
  treatment   potato  buttery    grassy   rancid   painty
1         1 6.887931 1.780087 0.6491379 4.065517 2.583621
2         2 7.001724 1.973913 0.6629310 3.624569 2.455844
3         3 6.967965 1.717749 0.6805195 3.866667 2.525541
```

足し算すると辞書式に並ぶ

```r
> reshape2::acast(molten, treatment ~ flavor + rep, mean)
  potato_1 potato_2 buttery_1 buttery_2  grassy_1  grassy_2 rancid_1 rancid_2 painty_1 painty_2
1 6.772414 7.003448  1.797391  1.762931 0.4456897 0.8525862 4.283621 3.847414 2.727586 2.439655
2 7.158621 6.844828  1.989474  1.958621 0.6905172 0.6353448 3.712069 3.537069 2.315517 2.597391
3 6.937391 6.998276  1.805217  1.631034 0.5895652 0.7706897 3.752174 3.980172 2.038261 3.008621
```

チルダで繋ぐと1次元増える

```r
> reshape2::acast(molten, treatment ~ flavor ~ rep, mean)
, , 1

    potato  buttery    grassy   rancid   painty
1 6.772414 1.797391 0.4456897 4.283621 2.727586
2 7.158621 1.989474 0.6905172 3.712069 2.315517
3 6.937391 1.805217 0.5895652 3.752174 2.038261

, , 2

    potato  buttery    grassy   rancid   painty
1 7.003448 1.762931 0.8525862 3.847414 2.439655
2 6.844828 1.958621 0.6353448 3.537069 2.597391
3 6.998276 1.631034 0.7706897 3.980172 3.008621
```
