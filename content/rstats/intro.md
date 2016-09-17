+++
title = 'R自学自習の基礎知識'
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -99
+++

## インターネットで調べる

-   [R-Tips](http://cse.naro.affrc.go.jp/takezawa/r-tips/r.html)
-   [RjpWiki](http://www.okada.jp.org/RWiki/)
-   <http://seekr.jp/> --- Googleカスタム検索 for R
-   <http://www.rseek.org/> --- Google Custom Search for R
-   <https://www.r-bloggers.com/>
-   <https://cran.r-project.org/web/views/>
-   <http://www.rdocumentation.org/>

## スクリプトを保存

R のコンソールにコマンドを打ち込むと、即座に結果が返ってくる。
一度きりで済むならそのように対話的な方法が分かりやすいかもしれないが、
別のデータにも使い回すとか、タイプミスでやり直すとか、
同じような処理を何度か繰り返す場合には(大抵はそうなる)、
一連の処理をまとめてスクリプト(テキストファイル)に書き出しておくとよい。
そこからのコピペでいろいろ楽ができるし、
あれどうやったんだっけと忘れる心配も減る(無くなりはしない)。
ファイルの拡張子は `.txt` でも何でもいいが `.R` にしとくと分かりやすい。

## パッケージ

便利な関数やサンプルデータなどをひとまとめにしたもの。
R開発チームが公式に作ってるものから、
ユーザーが自分用に作ってアップロードしたものまで、さまざまある。
自分でやろうとしてることは既に誰かがやってくれてる可能性が高いので、
車輪の再発明をする前に、まずは既存のパッケージを調べるべし。

### Standard Packages

https://stat.ethz.ch/R-manual/R-devel/doc/html/packages.html

R の標準機能。
何もしなくても使用可能な状態になっているので、
パッケージであることはあまり意識しなくてもいい。

`base`
: `c()`, `data.frame()`, `sum()` などホントに基本的なもの

`graphics`, `grDevices`, `grid`
: `plot()` などグラフ描画関連

`stats`
: `anova()`, `glm()`, `t.test()` など統計解析関連

`utils`
: `help()`, `install.packages()`, `read.table()` など

ほかに
`compiler`, `datasets`, `methods`, `parallel`,
`splines`, `stats4`, `tcltk`, `tools`

### Recommended Packages

R と一緒にインストールされるが、
使用する前に `library()` で呼び出しておく必要があるパッケージ。

例えば `MASS` に入ってる `stepAIC()` を使うには

```r
> stepAIC(model)
Error: could not find function "stepAIC"
> library(MASS)
> stepAIC(model)  # OK
```

ほかに
`boot`, `class`, `cluster`, `codetools`,
`foreign`, `KernSmooth`, `lattice`, `Matrix`,
`mgcv`, `nlme`, `nnet`, `rpart`, `spatial`, `survival`

### Contributed Packages

数千ものパッケージが有志により開発され、CRANにまとめて公開されている。

特にHadley Wickhamらによる
[tidyverse](https://github.com/hadley/tidyverse) パッケージ群 (例えば
[ggplot2]({{< relref "ggplot2.md" >}}),
[dplyr]({{< relref "dplyr.md" >}}),
[purrr]({{< relref "purrr.md" >}}),
[tidyr]({{< relref "tidyr.md" >}}),
[readr]({{< relref "readr.md" >}}),
[stringr]({{< relref "stringr.md" >}})など)
はどんな解析にも有用で、標準になってもいいくらい便利。
Rの中から下記のようなコマンドで一括インストール・読み込みできる。

```r
> install.packages('tidyverse')
> library(tidyverse)
```

そのほか
<https://cran.r-project.org/web/views/>
で用途別に紹介されている。


## 作業ディレクトリ

R ではどこかのディレクトリ (= フォルダ) に身をおいて作業する。
`read.table()` でファイルを開くときなどに
R がファイルを探すのはこの作業ディレクトリである。
`No such file or directory` と怒られる場合は、
作業ディレクトリとファイルの場所が合っていないかったり、
ファイル名のタイプミスだったりすることが多い。

`getwd()`, `setwd()`
:   現在地の working directory を get/set する関数

    ```r
    > getwd()
    [1] /Users/watal

    > setwd("~/Desktop")

    > getwd()
    [1] /Users/watal/Desktop
    ```

`list.files(path = ".", ...)` または `dir(path = ".", ...)`
:   ディレクトリ内のファイルを列挙する関数。
    path を省略するとワーキングディレクトリが対象となる

    ```r
    > read.table("mydata.txt", header=TRUE)
    Warning in file(file, "rt") :
      cannot open file 'mydata.txt': No such file or directory
    Error in file(file, "rt") : cannot open the connection
    > list.files()
    [1] mydate.txt
    # ファイル名が微妙に違う!!
    ```

## R で調べる

`help(topic, package=NULL, ...)`
:   ヘルプを表示する。
    知りたいオブジェクトの先頭にクエスチョンマークを付けるという方法もある。
    以下の2つは等価

    ```r
    > help(sum)
    > ?sum
    ```

`help.start()`
:   組み込みヘルプをウェブブラウザで開く。
    インストール済みパッケージのヘルプも Packages のリンクから見られる。

`attributes(obj)`, `str(obj)`
:   オブジェクトの属性や構造を調べる。
    `$names` に列挙されてる属性はダラー `$` を挟んで取り出せる

    ```r
    > x = rnorm(100)
    > y = rnorm(100)
    > mymodel = lm(y ~ x)
    > attributes(mymodel)
    $names
     [1] "coefficients"  "residuals"     "effects"       "rank"
     [5] "fitted.values" "assign"        "qr"            "df.residual"
     [9] "xlevels"       "call"          "terms"         "model"

    $class
    [1] "lm"
    >
    > mymodel$coefficients
    (Intercept)           x
     -0.1494812   0.1218096
    >
    > str(mymodel)
    List of 12
     $ coefficients : Named num [1:2] -0.171 -0.198
      ..- attr(*, "names")= chr [1:2] "(Intercept)" "x"
     $ residuals    : Named num [1:100] -0.5143 -1.3148 0.1954 1.6039 -0.0875 ...
      ..- attr(*, "names")= chr [1:100] "1" "2" "3" "4" ...
     $ effects      : Named num [1:100] 1.2027 2.0571 0.1179 1.8384 -0.0946 ...
      ..- attr(*, "names")= chr [1:100] "(Intercept)" "x" "" "" ...
    ```

## データ読み込み

See [readr]({{< relref "readr.md" >}})

## データ処理・整形

See [dplyr]({{< relref "dplyr.md" >}}), [purrr]({{< relref "purrr.md" >}}), and [tidyr]({{< relref "tidyr.md" >}})

## グラフ作図

2D: See [ggplot2]({{< relref "ggplot2.md" >}})

3D: See [rgl]({{< relref "rgl.md" >}})

## 参考図書

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4320018575/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41dIXM4egfL._SX160_.jpg" alt="統計学:Rを用いた入門書" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/400006973X/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51wEV2xJA1L._SX160_.jpg" alt="データ解析のための統計モデリング入門――一般化線形モデル・階層ベイズモデル・MCMC (確率と情報の科学)" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4274067831/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41PFcUnw4aL._SX160_.jpg" alt="The R Tips―データ解析環境Rの基本技・グラフィックス活用集" /></a>
