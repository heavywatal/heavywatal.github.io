+++
title = 'devtools'
subtitle = "Rパッケージ作成支援"
tags = ["r", "tidyverse", "package"]
[menu.main]
  parent = "rstats"
  weight = -50
+++

自分で書いた関数が多くなってきたら、まとめてパッケージを作るとよい。
少しだけ面倒だが、以下のようなメリットがある。

-   普通に変数や関数を定義するとワークスペースが名前でいっぱいになってしまうが、
    パッケージ内で定義されている変数や関数は `ls()` で出てこないのでスッキリ
-   既存のオブジェクトと名前が衝突するような場合でも、
    `mypackage::func1` のように名前空間を明示的に指定して呼び出せる

CRANに上げる程ではないとしても、
GitHubに公開しておけば誰でも使えるようになるので、
共同研究者と解析環境を共有したり、
ひとつの論文のワークフローを置いておいたり、いろいろ使い道はある。

## Rパッケージ

- <http://r-pkgs.had.co.nz/>
- [How to develop good R packages (for open science)](http://www.masalmon.eu/2017/12/11/goodrpackages/)

### 最低限の作成手順

1.  開発支援パッケージをインストールする:
    `install.packages(c("devtools", "usethis"))`

1.  [usethis](http://usethis.r-lib.org/) の関数をいくつか使って骨組みを作る:

    ```r
    usethis::create_package("hello")
    usethis::use_mit_license()
    usethis::use_roxygen_md()
    usethis::use_package_doc()
    ```

1.  `devtools::check()` で様子を見てみる。
    LaTeX 関連で怒られたら足りないパッケージを入れる:
    `sudo tlmgr install inconsolata helvetic`

1.  ローカルgitリポジトリを作って最初のコミットをする:
    `usethis::use_git()`

1.  GitHubかどこかに空のリポジトリを作る。
    パッケージと同名でなくてもよい。

1.  そのリモートリポジトリにプッシュ:

    ```sh
    git remote add origin https::/github.com/heavywatal/rhello.git
    git push -u origin master
    ```

1.  とりあえず誰でもインストール可能なパッケージができたはず:
    ```r
    devtools::install_github('heavywatal/rhello')
    ```


### 構造

```sh
DESCRIPTION  # 一番大事
NAMESPACE    # 見せるオブジェクトを列挙
README.md    # 全体の説明を簡単に
R/           # Rソースコード
data/        # サンプルデータなど
exec/        # 実行ファイル
inst/        # CITATION
man/         # ヘルプファイル.Rd
src/         # C++ソースコード
tests/
vignettes/
```

### [`DESCRIPTION`](http://r-pkgs.had.co.nz/description.html)

-   どうでも良さそうなファイル名とは裏腹に、ちゃんと書かないと動かない。
-   始めはusethisとかに生成してもらい、他のパッケージを参考にしつつ修正していく。
-   `Imports` に列挙したものは依存パッケージとして同時にインストールされる。
    しかし `library()` 時に同時に読み込まれるわけではない。
-   `Title` はピリオドを含まない一文。
    `Description` はピリオドを含む文章。
-   ライセンスを別ファイルにする場合は `License: file LICENSE` と書く
-   `Authors@R` のとこは後でRで評価されるので変な形。

### [`NAMESPACE`](http://r-pkgs.had.co.nz/namespace.html)

-   [後述のroxygen2](#roxygen2)がソースコードから自動生成するので**直接触らない**。
-   ここで `export()` された関数だけがユーザーから見える。
-   外部パッケージから `importFrom(package, function)`
    された関数はattachされてパッケージ内で利用可能になる。
    `import(package)` で全ての関数をまとめて処理できるけど名前の衝突が怖い。

### [`R/` ソースコード](http://r-pkgs.had.co.nz/r.html)

-   `NAMESPACE` や `man/*.Rd` を自動生成してもらえるように
    [後述のroxygen](#roxygen2)形式でコメントを書く
-   ファイルの数や名前は何でもいいので、開発者が分かりやすいようにしとく。
    例えば、似たような関数をひとつのファイルにまとめ、
    同名の`@rdname`を指定して`man/`と同じ構造にするとか。
-   `library()` や `require()` を書かない。
    必要なパッケージは `DESCRIPTION` の `Imports` に書き、
    `名前空間::関数()` のようにして使う。
    どうしても名前空間を省略したいものだけ `@importFrom` を指定する。
-   パッケージ読み込み時などに実行される特殊関数
    `.onLoad(...)`, `.onAttach(...)`, `.onUnload(...)`
    は慣習的に`zzz.R`というファイルに記述する。

### [`src/` C++ソースコード](http://r-pkgs.had.co.nz/src.html)

-   [Rcpp](https://cran.r-project.org/web/packages/Rcpp)
    が型変換などをスムーズにしてくれる:
    ```c++
    // [[Rcpp::plugins(cpp14)]]
    #include <Rcpp.h>

    //' First example
    //' @param args string vector
    //' @export
    // [[Rcpp::export]]
    Rcpp::IntegerVector len(const std::vector<std::string>& args) {
        return {static_cast<int>(args.size())};
    }
    ```
-   とりあえず `usethis::use_rcpp()` で設定を整える。
    `DESCRIPTION` や `src/.gitignore` などが書き換えられる。
-   `R/***-package.R` など適当なとこに `@useDynLib` の設定を追加:

    ```r
    #' @useDynLib hello, .registration = TRUE
    #' @importFrom Rcpp sourceCpp
    #' @keywords internal
    "_PACKAGE"
    ```
    `@importFrom Rcpp sourceCpp` を省くと、パッケージ利用時に
    `'enterRNGScope' not provided by package 'Rcpp'`
    のようなエラーが出る場合がある
    (明示的に `library(Rcpp)` するなどして既にRcppロード済みの環境では動く)。
-   同じところに `.onUnload` も定義しておく:
    ```r
    .onUnload = function(libpath) {
      library.dynam.unload("hello", libpath)
    }
    ```
    すると `unloadNamespace("hello")` したときに共有ライブラリもちゃんと外れるようになる。
    ちなみに `devtools::unload()` はこれを省略してもちゃんとリロードしてくれる。
-   外部ライブラリのリンクに関する設定など、
    開発者側で指定すべきビルドオプションは `src/Makevars` に指定:
    ```
    CXX_STD=CXX14
    PKG_CPPFLAGS=-DSTRICT_R_HEADERS -I/usr/local/include
    PKG_LIBS=-L/usr/local/lib -Wl,-rpath,/usr/local/lib -lthankyou
    ```
    `STRICT_R_HEADERS` を定義しておくことで余計なマクロ定義を防げる。
    `configure` や [CMake]({{< relref "cmake.md" >}}) を使って
    `src/Makevars.in` から生成する手もある。

-   ユーザ側で指定すべきオプションがある場合は `~/.R/Makevars` に:
    ```
    LLVM_LOC=/usr/local/opt/llvm
    CC=$(LLVM_LOC)/bin/clang
    CXX=$(LLVM_LOC)/bin/clang++
    # CFLAGS=-g -Wall -O2 -march=native -mtune=native
    # CXXFLAGS=$(CFLAGS)
    # CPPFLAGS=-I${HOME}/local/include -I${HOME}/.linuxbrew/include
    # LDFLAGS=-L${HOME}/local/lib -L${HOME}/.linuxbrew/lib
    ```
    例えばMPI依存パッケージをmacOSでビルドでしようとすると
    `clang: error: unsupported option '-fopenmp'`
    と怒られるので `brew install llvm`
    で別のコンパイラを入れて上記のように指定する。

その他の注意点

- `std::abort()` や `std::exit()` は呼び出したRセッションまで殺してしまう。
  例外は投げっぱなしで拾わなくても大丈夫で、
  `std::exception`の派生クラスなら`what()`まで表示してもらえる。
- グローバル変数やクラスのstaticメンバは `dyn.unload()` されるまで生き続ける。
  `parallel::mclapply()` とかでフォークした先での変更は子同士にも親にも影響しない。

### [`vignettes/`](http://r-pkgs.had.co.nz/vignettes.html)

個々の関数の使用例はRソースファイルの `@examples` に書くとして、
複数の関数を組み合わせてどう使うかとか、
パッケージ全体の使い方とかを説明するのが`vignettes/`の役割。
Rmarkdown形式で書いてHTMLやPDFで表示できるので表現力豊か。

`usethis::use_vignette('hello')` で雛形を作ってもらうのが楽。

`pandoc` と `pandoc-citeproc` が必要なので
[Homebrew]({{< relref "homebrew.md" >}}) とかでインストールしておく。

`check()`がデフォルトで`vignette=TRUE`かつ処理がやや重いので、
毎回その重さを受け入れるか、わざわざ`FALSE`指定するかというのは悩みどころ。


### [`inst/`](http://r-pkgs.had.co.nz/inst.html)

ここに入ってるものはインストール先でトップディレクトリに移される謎仕様。

論文で引用されることを想定している場合は `inst/CITATION` を作る。
`citation('ggplot2')` のように参照できるようになる。

### [`demo/`](http://r-pkgs.had.co.nz/misc.html)

vignettesに取って代わられた古い機能。
ソースコード`*.R`を置いておくと`demo()`関数で呼び出せるというだけ。

`check()`でソースコードの中身は実行されないが、
`demo/00Index` というファイルの整合性が取れてないと警告される。
「デモの名前 + スペース3つ以上かタブ + 適当な説明」という形式。
ファイル名から拡張子を取り除いたものがデモの名前になる。
```
mydemo1    Description of demo 1
mydemo2    Description of demo 2
```

人に見せたり`tests/`に入れたりするほどのもんでもないし、
毎回`check()`されなくてもいいけど、
試しに書いたコードを一応取っとく、くらいの用途にはいいかな？



## `devtools`

<a href="https://devtools.r-lib.org/">
<img src="https://devtools.r-lib.org/reference/figures/logo.svg" align="right" width="120" height="139">
</a>

骨組みを作るとこからCRANにデプロイするとこまでお世話してくれる。

### 主な関数

`document(pkg='.', ...)`
:   `roxygen2` を呼び出してソースコードから
    `NAMESPACE` や `man/*.Rd` を自動生成する。

`check(pkg='.', document=TRUE, cleanup=TRUE, cran=TRUE, check_version=FALSE, ...)`
:   パッケージとしての整合性を確認。
    ついでに`document()`は実行できるけど`spell_check()`はできないので手動で。

`test(pkg='.', filter=NULL, ...)`
:   `testthat` を呼び出して `test/` 以下のテストコードを実行する

`install(pkg='.', reload=TRUE, quick=FALSE, local=TRUE, ...)`
:   ローカルにあるディレクトリからインストール

`install_github(repo, username=NULL, ref='master', subdir=NULL, ...)`
:   GitHubリポジトリからインストール

`unload(pkg='.')`
:   `datach('package:XXX')` とか `unloadNamespace(XXX)`
    よりもちゃんとまっさらにパッケージを外す。

`load_all(pkg='.', reset=TRUE, recompile=FALSE, export_all=TRUE, quiet=FALSE)`
:   `install()` せずファイルから直接 `library()` する。
    ロード済みでもまず `unload()` が呼ばれるので安心。

`clean_dll(pkg='.')`
:   `src/` 以下に生成される `.o`, `.so` を消す。
    普段は触る必要ないが、たまにこれが必要な不具合に出くわす。


### 設定

項目の説明を読む

```r
?devtools
```

例えば `.Rprofile` に

```r
options(devtools.desc.author='Watal M. Iwasaki <user@example.com> [aut, cre]')
options(devtools.desc.license='MIT')
```

## `roxygen2`

<a href="https://CRAN.R-project.org/package=roxygen2">
<img src="https://raw.githubusercontent.com/klutometis/roxygen/master/man/figures/logo.png" align="right">
</a>

Rソースコードのコメントから`NAMESPACE`とヘルプ(`man/*.Rd`)を自動生成する。

- <https://cran.r-project.org/web/packages/roxygen2/>
- <https://github.com/klutometis/roxygen>
- <http://r-pkgs.had.co.nz/man.html>
- <http://kbroman.org/pkg_primer/pages/docs.html>
- <https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html>

`roxygen2::roxygenise(package.dir='.', ..., clean=FALSE)`
を直接呼んでもよいが、
基本的には `devtools::document()` を使って間接的に利用する。

### 使い方

```r
#' A simple function to add 1
#' @param x A numeric vector
#' @export
#' @examples
#' increment(42)
increment = function(x) {x + 1}
```

-   `#' ` から始まる行がroxygenコメントとして扱われる。
-   タグは `@` で始まる。
    `@` そのものを入力したいときは重ねて `@@` とする。
-   1行目には必ずタイトルを書く。
    タグ `@title` は省略可能。
    とりあえずコピペして全部同じタイトルにしようとすると、
    関数ごとにユニークなタイトルをつけろと怒られる。
-   2段落目は `@description` 扱いされる。
    省略するとタイトルが流用される。
-   空行だけでは切れ目として扱われないので `NULL` などを置いたりする。
-   dplyrなどで列名を直に指定すると
    `undefined global variables`
    という警告が出るので
    `dplyr::filter(iris, .data$Species == 'setosa')`
    のように pronoun を使って抑える。
    どこかに `#' @importFrom rlang .data` を書いておく。
    ただし `$` によるアクセスがやや遅いので、
    `group_by()` などで多数のグループを処理するときなど気になる場合は
    `!!as.name("carat")` のようにしたほうが速い。
-   Markdownを使うためには
    `Roxygen: list(markdown = TRUE)` を `DESCRIPTION` に加える。
-   `"_PACKAGE"` という文字列の上に書かれたブロックは、
    パッケージそのものに対応するヘルプ `*-package.Rd` になる。
    `@useDynLib` や `@importFrom` など全体に関わる設定はここで。

    ```r
    #' Example package to say hello
    #' @aliases NULL hello-package
    #' @useDynLib hello, .registration = TRUE
    #' @importFrom Rcpp sourceCpp
    #' @importFrom magrittr %>%
    #' @importFrom rlang .data
    #' @keywords internal
    "_PACKAGE"
    ```


### タグ

`@import pkg1, pkg2, ...`
:   `NAMESPACE` で `import()` するパッケージを指定。
    名前の衝突が怖いので基本的に使わない。

`@importFrom pkg func`
:   `NAMESPACE` で `importFrom()` するパッケージと関数を指定。
    e.g., `@importFrom magrittr %>%`

`@export`
:   `NAMESPACE` で `export()` する関数を指定。

`@param arg1 description...`
:   関数の引数。型や役割の説明を書く。

`@inheritParams package::function`
:   未`@param`引数の記述を別の関数から継承する。

`@return description`
:   関数の返り値。

`@examples code...`
:   例となるコードを記述する。
    exportしないやつには書いてはいけない。
    `\dontrun{}` に入れるとチェックから除外される。
    単数形の `@example` は外部ファイルのパスを受け取る

`@rdname basename`
:   `man/`に書き出すRdファイルの名前。
    複数の関数で同じものを指定すればひとつのヘルプにまとめられる。
    このとき`@param`などは共有されるので、
    同じ名前のものはどこかで1度だけ記述する。
    逆に言えば、中身が違うものに同じ名前をつけてはいけないし、
    引数や機能がほとんど重ならない関数をまとめると分かりにくくなる。

`@docType`, `@name`
:   パッケージやデータを記述するのに必要だったが今では不要っぽい。


## 関連書籍

<a href="https://www.amazon.co.jp/Packages-Organize-Test-Document-Share-ebook/dp/B00VAYCHL0/ref=as_li_ss_il?ie=UTF8&qid=1477817362&sr=8-1&keywords=r+package&linkCode=li3&tag=heavywatal-22&linkId=c760a9619f6d31273f668389360b687d" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=B00VAYCHL0&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=B00VAYCHL0" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%83%91%E3%83%83%E3%82%B1%E3%83%BC%E3%82%B8%E9%96%8B%E7%99%BA%E5%85%A5%E9%96%80-_%E3%83%86%E3%82%B9%E3%83%88%E3%80%81%E6%96%87%E6%9B%B8%E5%8C%96%E3%80%81%E3%82%B3%E3%83%BC%E3%83%89%E5%85%B1%E6%9C%89%E3%81%AE%E6%89%8B%E6%B3%95%E3%82%92%E5%AD%A6%E3%81%B6-Hadley-Wickham/dp/4873117593/ref=as_li_ss_il?ie=UTF8&qid=1477817362&sr=8-6&keywords=r+package&linkCode=li3&tag=heavywatal-22&linkId=f322b7223574d738364d6675a58a9de9" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117593&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117593" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
