+++
title = 'devtools'
subtitle = "Rパッケージ作成支援"
tags = ["r", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -50
+++

<http://r-pkgs.had.co.nz/>

自分で書いた関数が多くなってきたら、まとめてパッケージを作るとよい。
少しだけ面倒だが、以下のようなメリットがある。

-   普通に変数や関数を定義するとワークスペースが名前でいっぱいになってしまうが、
    パッケージ内で定義されている変数や関数は `ls()` で出てこないのでスッキリ
-   既存のオブジェクトと名前が衝突するような場合でも、
    `mypackage::func1` のように名前空間を明示的に指定して呼び出せる
-   バイトコンパイルされる分だけ関数がちょっと速くなるかも

CRANに上げる程ではないとしても、
GitHubに公開しておけば誰でも使えるようになるので、
共同研究者と解析環境を共有したり、
ひとつの論文のワークフローを置いておいたり、いろいろ使い道はある。

## Rパッケージ

### 最低限の作成手順

1.  Rを起動して `devtools` をインストールする

    ```r
    install.packages('devtools')
    library(devtools)
    ```

2.  `devtools::create()` で骨組みを作る

    ```r
    setwd('~/tmp/')
    library(devtools)
    pkgname = 'namaespace'
    title_desc = 'Utility to create dummy R packages for namespace'
    github_repo = sprintf('https://github.com/heavywatal/%s', pkgname)
    devtools::create(pkgname, description=list(
        Package=pkgname,
        Title=title_desc,
        Description=paste0(title_desc, '.'),
        `Authors@R`="person('Watal M.', 'Iwasaki', email='user@example.com', role=c('aut', 'cre'))",
        License='MIT',
        Suggests='dplyr, ggplot2, purrr, readr, tidyr',
        Imports='devtools, stringr',
        URL=github_repo,
        BugReports=paste0(github_repo, '/issues'))
    ```

3.  `devtools::check(pkgname)` で様子を見てみる

    {{%div class="note"%}}
LaTeX 関連で怒られたら足りないパッケージを入れる:
`sudo tlmgr install inconsolata helvetic`
    {{%/div%}}

4.  GitHubに空のリポジトリを作る
    -   パッケージと同名でなくてもよい
    -   `README` や `LICENSE` は後から作る

5.  コミットしてプッシュ

    ```sh
    % git init
    % git add --all
    % git commit -m "first commit"
    % git remote add origin git@github.com:heavywatal/namaespace.git
    % git push -u origin master
    ```

6.  とりあえず誰でもインストール可能なパッケージができたはず

    ```r
    install_github('heavywatal/namaespace')
    ```

7.  あとは `R/` にソースコードを置いたり `README.md` を書いたり

### 構造

```sh
DESCRIPTION  # 一番大事
NAMESPACE    # 見せるオブジェクトを列挙
README.md    # 全体の説明を簡単に
R/           # Rソースコード
data/        # サンプルデータなど
man/         # ヘルプファイル.Rd
src/         # C++
tests/
vignettes/
```

### `DESCRIPTION`

-   どうでも良さそうなファイル名とは裏腹に、ちゃんと書かないと動かない
-   始めは `devtools::create()` で自動生成し、それから修正していく
-   `Imports` に列挙したパッケージは同時にインストールされる。
    しかし `library()` 時に同時に読み込まれるわけではない。
-   `Title` はピリオドを含まず、Descriptionはピリオドを含むようにする
-   ライセンスを別ファイルにする場合は `License: file LICENSE` と書く
-   `Authors@R` のとこは後でRで評価されるので変な形

### `NAMESPACE`

-   `roxygen2` (下記) がソースコードから自動生成するので**直接触らない**
-   ここで `export()` された関数だけがユーザーから見える
-   ここで `import()` された関数はパッケージ内でattachされた状態になるが、
    そうしないで毎回 `名前空間::` 越しにアクセスしたほうがよい。
    `magrittr::%>%` とか `pipeR::%>>%` のような演算子は仕方がないので
    `importFrom()` で個別指定する。

### Rソースコード

-   `R/` 以下に配置
-   `NAMESPACE` や `man/*.Rd` を自動生成してもらえるように
    [後述のroxygen](#roxygen2)形式でコメントを書く
-   ファイルの数や名前は何でもいいので、開発者が分かりやすいようにしとく。
    例えば、似たような関数をひとつのファイルにまとめ、同名の`@rdname`を指定し、
    `man/`と同じ構造にするのがちょうどいいのではなかろうか。
-   `library()` や `require()` を書かない。
    必要なパッケージは `DESCRIPTION` の `Imports` に書き、
    `名前空間::関数()` のようにして使う。
    どうしても名前空間を省略したいものだけ `@importFrom` を指定する。
-   パッケージ読み込み時などに実行される特殊関数
    `.onLoad(...)`, `.onAttach(...)`, `.onUnload(...)`
    は慣習的に`zzz.R`というファイルに記述する。

### C++ソースコード

-   `src/` 以下に配置
-   [Rcpp](https://cran.r-project.org/web/packages/Rcpp)
    が型変換などをスムーズにしてくれる。
-   `Rcpp.h`のインクルードといくつかのコメントが必要
    ```c++
    // [[Rcpp::plugins(cpp14)]]
    #include <Rcpp.h>

    //' First example
    //' @param args string vector
    //' @return string
    //' @export
    // [[Rcpp::export]]
    std::string cxx_func(Rcpp::CharacterVector args=Rcpp::CharacterVector::create()) {
        auto vs_args = Rcpp::as<std::vector<std::string>>(args);
        return vs_args;
    }
    ```

-   開発者側で指定すべきビルドオプションは `src/Makevars` に指定
    ```
    CXX_STD = CXX11
    PKG_CPPFLAGS = -isystem ${HOME}/local/include
    PKG_LIBS = -lmy_great_lib
    ```

-   ユーザ側で指定すべきオプションは `~/.R/Makevars` に
    ```
    CFLAGS = -g -Wall -O2 -march=native -mtune=native
    CXXFLAGS = $(CFLAGS)
    CXX1XFLAGS = $(CFLAGS)
    CXX1XSTD = -std=c++14
    LDFLAGS = -L${HOME}/local/lib -L${HOME}/.homebrew/lib
    ```

-   `DESCRIPTION` にいくらか書き足す
    ```
    LinkingTo: Rcpp
    SystemRequirements: C++14
    ```

-   `@docType package` のブロックに `@useDynLib` でパッケージ名を指定する
    ```R
    #' @docType package
    #' @useDynLib namaespace
    NULL
    ```

その他の注意点

- `std::abort()` や `std::exit()` は呼び出したRセッションまで殺してしまう。
  例外は投げっぱなしで拾わなくても大丈夫で、
  `std::exception`の派生クラスなら`what()`まで表示してもらえる。

### vignettes

<http://r-pkgs.had.co.nz/vignettes.html>

個々の関数の使用例はRソースファイルの `@examples` に書くとして、
複数の関数を組み合わせてどう使うかとか、
パッケージ全体の使い方とかを説明するのが`vignettes/`の役割。
Rmarkdown形式で書いてHTMLやPDFで表示できるので表現力豊か。

`build()`や`check()`がデフォルトで`vignette=TRUE`かつ処理がやや重いので、
毎回その重さを受け入れるか、わざわざ`FALSE`指定するかというのは悩みどころ。

`pandoc` と `pandoc-citeproc` が必要なので
[Homebrew]({{< relref "homebrew.md" >}}) とかでインストールしておく。

### demo

<http://r-pkgs.had.co.nz/misc.html>

vignettesに取って代わられた古い機能。
ソースコード`*.R`を置いておくと`demo()`関数で呼び出せるというだけ。

`check()`でソースコードの中身は参照されないが、
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

<https://github.com/hadley/devtools>

骨組みを作るとこからCRANにデプロイするとこまでお世話してくれる。

### 関数

`create(path, description=getOption('devtools.desc'), check=FALSE, rstudio=TRUE)`
:   まっさらな状態から骨組みを作る。
    `path` が既に存在している場合は先に進めないので、
    やり直すときは `system('rm -rf path')` などとして一旦消す必要がある。

`document(pkg='.', ...)`
:   `roxygen2` を呼び出してソースコードから
    `NAMESPACE` や `man/*.Rd` を自動生成する

`check(pkg='.', document=TRUE, cleanup=TRUE, cran=TRUE, check_version=FALSE, ...)`
:   パッケージとしての整合性を確認。
    `document()`も呼ばれる。

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

Rソースコードのコメントから`NAMESPACE`とヘルプ(`man/*.Rd`)を自動生成する。

- <http://cran.r-project.org/web/packages/roxygen2/>
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
#' @export
#' @param x A numeric vector
#' @return A numeric vector
#' @examples
#' increment(42)
increment = function(x) {x + 1}
```

-   `#'` から始まる行がroxygenコメントとして扱われる。
-   タグは `@` で始まる。
    `@` そのものを入力したいときは重ねて `@@` とする。
-   1行目には必ずタイトルを書く。
    タグ `@title` は省略可能。
    とりあえずコピペして全部同じタイトルにしようとすると、
    関数ごとにユニークなタイトルをつけろと怒られる。
-   2段落目は `@description` 扱いされる。
    省略するとタイトルが流用される。
-   空行だけでは切れ目として扱ってくれないので、
    データやパッケージのドキュメントを書く場合は `NULL` を入れる。

### タグ

`@import pkg1, pkg2, ...`
:   `NAMESPACE` で `import()` するパッケージを指定。
    名前の衝突が怖いので基本的に使わない。

`@importFrom pkg func`
:   `NAMESPACE` で `importFrom()` するパッケージを指定。
    e.g., `@importFrom pipeR %>>%`

`@export`
:   `NAMESPACE` で `export()` するパッケージを指定。

`@param arg1 description...`
:   関数の引数。型や役割の説明を書く。

`@inheritParams package::function`
:   未`@param`引数の記述を別の関数から継承する。
    継承の連鎖はしないので、親関数で`@param`定義されてないものはダメ。

`@return description`
:   関数の返り値。

`@examples code...`
:   例となるコードを記述する。
    exportしないやつには書いてはいけない。
    単数形の `@example` は外部ファイルのパスを受け取る

`@rdname basename`
:   `man/`に書き出すRdファイルの名前。
    複数の関数で同じものを指定すればひとつのヘルプにまとめられる。
    このとき`@param`などは共有されるので、
    同じ名前のものはどこかで1度だけ記述する。
    逆に言えば、中身が違うものに同じ名前をつけてはいけないし、
    引数や機能がほとんど重ならない関数をまとめると分かりにくくなる。

`@docType`
:   関数やクラスには不要だが `data` か `package` の場合はここで指定。

`@name name`
:   パッケージやデータはこのタグで明示的に設定する必要がある。
