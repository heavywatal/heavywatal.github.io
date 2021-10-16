+++
title = 'R自学自習の基礎知識'
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -99
+++

## 初学者向け講義資料

- 2021-10 [Rによるデータ前処理実習](/slides/tmd2021/) 東京医科歯科大学 データ関連人材育成プログラム
    - 入門1: 前処理とは。Rを使うメリット。Rの基本。
    - 入門2: データ可視化の重要性と方法。
    - データ構造の処理1: 抽出、集約など。
    - データ構造の処理2: 結合、変形など。
    - データ内容の処理: 数値、文字列、日時など。
- [Rを用いたデータ解析の基礎と応用](https://comicalcommet.github.io/r-training-2021/)
   (石川由希 2021 名古屋大学):<br>
  初心者に寄り添ってさらに噛み砕いた説明。
  [**よくあるエラー集**](https://comicalcommet.github.io/r-training-2021/R_training_2021_7.html)が特に重宝。


## R環境のインストール

R本体
: コマンドを解釈して実行するコア部分
: よく使われる関数なども標準パッケージとして同梱
: https://cran.rstudio.com/ からダウンロードしてインストール

RStudio Desktop
: Rをより快適に使うための総合開発環境(IDE)
: 必須ではないけど、結構みんな使ってるらしい
: https://rstudio.com/ からダウンロードしてインストール

## スクリプトを保存

R のコンソールにコマンドを打ち込むと、即座に結果が返ってくる。
一度きりで済む短い処理ならそのように対話的な方法が分かりやすいかもしれないが、
別のデータにも使い回すとか、タイプミスでやり直すとか、
同じような処理を何度か繰り返す場合には(大抵はそうなる)、
一連の処理をまとめてスクリプト(テキストファイル)に書き出しておくとよい。
ファイルの拡張子は `.txt` でも何でもいいが `.R` にすることが多い。

[R Markdown](https://rmarkdown.rstudio.com/)
を使えばコードと解析結果の図表を同時に保存・プレゼンすることができる。


## パッケージ

便利な関数やサンプルデータなどをひとまとめにしたもの。
R開発チームが公式に作ってるものから、
ユーザーが自分用に作ってアップロードしたものまで、さまざまある。
自分でやろうとしてることは既に誰かがやってくれてる可能性が高いので、
車輪の再発明をする前に、まずは既存のパッケージを調べるべし。
プログラミングは、イチから「書く」のではなく、
既存のパーツを利用して「組む」感覚のほうが近い。


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
stepAIC(model)
## Error: could not find function "stepAIC"
library(MASS)
stepAIC(model)  # OK
```

ほかに
`boot`, `class`, `cluster`, `codetools`,
`foreign`, `KernSmooth`, `lattice`, `Matrix`,
`mgcv`, `nlme`, `nnet`, `rpart`, `spatial`, `survival`

### Contributed Packages

<a href="https://tidyverse.tidyverse.org/">
<img src="https://tidyverse.tidyverse.org/logo.png" align="right" width="120" height="139">
</a>

数千ものパッケージが有志により開発され、CRANにまとめて公開されている。

例えば [rstudio.com/products/rpackages](https://www.rstudio.com/products/rpackages/)
で紹介されているもの、特にHadley Wickhamらによる
[tidyverse](https://www.tidyverse.org/) パッケージ群
([ggplot2]({{< relref "ggplot2.md" >}}),
[dplyr]({{< relref "dplyr.md" >}}),
[tidyr]({{< relref "tidyr.md" >}}),
[purrr]({{< relref "purrr.md" >}}),
[readr]({{< relref "readr.md" >}}),
[stringr]({{< relref "stringr.md" >}})など)
はどんな解析にも有用で、標準になってもいいくらい便利。
というより、tidyverseなき裸のRは使う気が起きない。
Rの中から下記のようなコマンドで一括インストール・読み込みできる。

```r
install.packages("tidyverse")  # 一度やればOK
library(tidyverse)             # 読み込みはRを起動するたびに必要
update.packages()              # たまには更新しよう
```

そのほか
<https://cran.r-project.org/web/views/>
で用途別に紹介されている。

パッケージを作るには [devtools]({{< relref "devtools.md" >}}) を使う。


## 作業ディレクトリ

Rではどこかのディレクトリ(=フォルダ)に身をおいて作業する。
[データファイルを読み込む]({{< relref "readr.md" >}})ときなどに
Rがファイルを探すのはこの作業ディレクトリである。
`No such file or directory` と怒られる場合は、
作業ディレクトリとファイルの場所が合っていないかったり、
ファイル名のタイプミスだったりすることが多い。

作業ディレクトリは解析の途中でホイホイ移動せず、
プロジェクトの最上階などに固定しておいて、
そこからの相対パスでファイルを読み書きするのが分かりやすくて安全。

`getwd()`, `setwd()`
:   現在地 **w**orking **d**irectory をget/setする関数。

    ```r
    getwd()
    ## [1] /Users/watal
    setwd("~/Desktop")
    getwd()
    ## [1] /Users/watal/Desktop
    ```

`list.files(path = ".", ...)` または `dir(path = ".", ...)`
:   ディレクトリ内のファイルを列挙する関数。
    path を省略するとワーキングディレクトリが対象となる

    ```r
    read.table("mydata.txt", header = TRUE)
    ## Warning in file(file, "rt") :
    ##   cannot open file 'mydata.txt': No such file or directory
    ## Error in file(file, "rt") : cannot open the connection
    list.files()
    ## [1] mydate.txt   ファイル名が微妙に違う!!
    ```

## R で調べる

`help(topic, package = NULL, ...)`
:   ヘルプを表示する。
    知りたいオブジェクトの先頭にクエスチョンマークを付けるという方法もある。
    以下の2つは等価

    ```r
    help(sum)
    ?sum
    ```

`help.start()`
:   組み込みヘルプをブラウザで開く。
    インストール済みパッケージのヘルプも Packages のリンクから見られる。

`attributes(obj)`, `str(obj)`
:   オブジェクトの属性や構造を調べる。
    `$names` に列挙されてる属性はダラー `$` を挟んで取り出せる

    ```r
    x = rnorm(100)
    y = rnorm(100)
    mymodel = lm(y ~ x)
    attributes(mymodel)
    ## $names
    ##  [1] "coefficients"  "residuals"     "effects"       "rank"
    ##  [5] "fitted.values" "assign"        "qr"            "df.residual"
    ##  [9] "xlevels"       "call"          "terms"         "model"
    ## $class
    ## [1] "lm"
    mymodel$coefficients
    ## (Intercept)           x
    ##  -0.1494812   0.1218096
    str(mymodel)
    ## List of 12
    ##  $ coefficients : Named num [1:2] -0.171 -0.198
    ##   ..- attr(*, "names")= chr [1:2] "(Intercept)" "x"
    ##  $ residuals    : Named num [1:100] -0.5143 -1.3148 0.1954 1.6039 -0.0875 ...
    ##   ..- attr(*, "names")= chr [1:100] "1" "2" "3" "4" ...
    ##  $ effects      : Named num [1:100] 1.2027 2.0571 0.1179 1.8384 -0.0946 ...
    ##   ..- attr(*, "names")= chr [1:100] "(Intercept)" "x" "" "" ...
    ```

## インターネットで調べる・尋ねる

- [r-wakalang](https://github.com/tokyor/r-wakalang) --- Slack上の日本語コミュニティ
- <https://r4ds.had.co.nz/> --- R for Data Science (体系的に学びたい人はぜひ通読を)
- <https://www.rdocumentation.org/> --- 各パッケージの公式ドキュメントを閲覧
- [Stack Overflow](https://stackoverflow.com/questions/tagged/r)

## データ読み込み

See [readr]({{< relref "readr.md" >}}).

## データ処理・整形

See [dplyr]({{< relref "dplyr.md" >}}), [tidyr]({{< relref "tidyr.md" >}}), and [purrr]({{< relref "purrr.md" >}}).

## グラフ作図

2D: See [ggplot2]({{< relref "ggplot2.md" >}}).

3D: See [rgl]({{< relref "rgl.md" >}}).

## 関連書籍

Rのモダンな使い方、考え方を学ぶにはr4ds。
まずは[公開オンライン版](https://r4ds.had.co.nz/)を読んでみて。<br>
<a href="https://www.amazon.co.jp/dp/1491910399/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=ebfaf931addcf48f288a9fb8dfebd0c9" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/487311814X/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=c508bbce8036b379868bc9628363583f" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />

即戦力が欲しい場合は、新しく日本語で書かれた宇宙本が取っつきやすい:<br>
<a href="https://www.amazon.co.jp/dp/4297121700?&linkCode=li3&tag=heavywatal-22&linkId=77762a4d0080a840ec5d94df9c0c5ceb&language=ja_JP&ref_=as_li_ss_il" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4297121700&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4297121700" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />

統計モデルの基礎は緑本で:<br>
<a href="https://www.amazon.co.jp/dp/400006973X/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=eea6b2ec8a0fbef9df8c163f92cc6ced" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=400006973X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=400006973X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
