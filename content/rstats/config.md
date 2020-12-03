+++
title = 'R環境設定'
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -95
+++

https://cran.r-project.org/manuals.html

## インストール

- [R本体](https://cran.rstudio.com/)
- [RStudio (任意)](https://www.rstudio.com/products/rstudio/download/)

Macなら [Homebrew]({{< relref "homebrew.md" >}}) で
`brew install --cask r rstudio` のように入れるのが楽チン。
Caskじゃない `brew install r` のほうだとバイナリ版パッケージが使えなくて毎回ソースからビルドさせられるので大変。

https://cran.r-project.org/doc/manuals/R-admin.html


## 起動オプション

ワークスペースの自動保存や自動復帰は危険なので切っておく。
R.app や RStudio から使う場合はメニューから環境設定みたいなやつを開く。
シェルから使う場合は例えば以下のようなエイリアスを設定する。

```sh
alias r='R --quiet --no-save --no-restore-data'
```

詳しくは `R --help` または
https://cran.r-project.org/doc/manuals/R-intro.html#Invoking-R


## 環境変数

https://cran.r-project.org/doc/manuals/R-admin.html#Environment-variable-index

https://stat.ethz.ch/R-manual/R-patched/library/base/html/EnvVar.html

```r
Sys.getenv()
```


## .Renviron

<https://stat.ethz.ch/R-manual/R-patched/library/base/html/Startup.html>

R起動時に読み込まれ、環境変数を設定するファイル。
Rスクリプトではなく、シェルスクリプトっぽい代入式で書く。
例 (<https://github.com/heavywatal/dotfiles/blob/master/.R/.Renviron>):

```sh
R_USER=${HOME}/.R
R_LIBS_USER=${R_USER}/library
R_ENVIRON_USER=${R_USER}/.Renviron
R_PROFILE_USER=${R_USER}/.Rprofile
R_HISTFILE=${R_USER}/.Rhistory
R_HISTSIZE=65535
LANG=C
LC_CTYPE=en_US.UTF-8
```

探される・読み込まれる順序はだいたい以下のとおり:

1. `$R_ENVIRON`
1. `$R_HOME/etc/Renviron.site`
1. `$R_ENVIRON_USER`
1. `./.Renviron`
1. `~/.Renviron`

読み込ませたくないときは `--no-environ` オプション。


## .Rprofile

https://cran.r-project.org/doc/manuals/R-intro.html#Customizing-the-environment

R起動時に読み込まれるファイル。
中身はRスクリプトなので、パッケージの読み込みや関数の定義など、Rでできることは何でもできるはず。
例: <https://github.com/heavywatal/dotfiles/blob/master/.R/.Rprofile>

`.First()` と `.Last()` はそれぞれ起動時と終了時に実行される関数。
これらが原因で `R CMD` やパッケージ関連の操作が失敗することもあるので、
普通の対話環境でのみ有効になるよう `if (interactive())` などとしておいたほうが安心。

読み込まれる順序はだいたい以下のとおり。
`.Renviron` のほうが先に読み込まれるので、
上記のように `R_PROFILE_USER` を定義しておいて、そこに置いとけば読み込まれる。

1. `$R_PROFILE`
1. `$R_HOME/etc/Rprofile.site`
1. `$R_PROFILE_USER`
1. `./.Rprofile`
1. `~/.Rprofile`

読み込ませたくないときは `--no-init-file` オプション。


### `options()`

<https://stat.ethz.ch/R-manual/R-patched/library/base/html/options.html>

`?options` で項目の一覧を見られる。

スクリプトの結果が設定依存で変わっては困るので、
そういう本質的なものはいじらずに、
表示関連の項目だけ変えるに留めたほうがよい。
例えば `options(stringsAsFactors=FALSE)` とやってしまいたくなるが、
それを設定した環境とそうでない環境で `read.table()` の結果が変わるのは危険。
関数の引数として毎回指定するか、
[readr]({{< relref "readr.md" >}})の関数を使うべし。

`warn=1`
: 警告レベルの設定。
  デフォルト(`warn=0`)では、警告があっても計算は滞り無く進行し、
  最後に "There were 50 or more warnings (use warnings() to see the first 50)"
  などと軽く表示されるだけなので、
  見落としたりして後々大変なバグ取り作業に発展する恐れがある。
  警告が発生するごとに警告文を表示する(`warn=1`)か、
  エラー扱いにして計算をストップするようにしておく(`warn=2`)ことでそれを回避すべし。
  ちなみに負数だと警告無視。

`warnPartialMatchAttr`, `warnPartialMatchDollar`
: listやdata.frameなどの要素を抜き出すとき、
  対象が一意に定まる範囲で変数名の省略が許されてしまう (e.g., `mtcars$m`)。
  これは危険なので、せめて警告がでるように設定する。
  [tibble]({{< relref "readr.md#tibble" >}}) を使うほうがより安全。

`warnPartialMatchArgs`
: 関数の引数名の省略に関する警告。
  自分のコーディングに関しては`TRUE`にしておきたいけど、
  結構いろんなパッケージが警告を発してうるさいので仕方なく`FALSE`。

`showWarnCalls=TRUE`, `showErrorCalls=TRUE`
: 警告やエラーの出処をたどって表示する。
  Pythonほどわかりやすくないが、ちょっとはマシになる。


## ライブラリの管理

https://cran.r-project.org/doc/manuals/R-admin.html#Add_002don-packages

### `~/.R/Makevars`

パッケージをソースコードからビルドするときの設定。
大概はビルド済みのものを入れるので不要だけど、たまに必要になる。

例えば、素のmacOSでOpenMP依存のパッケージをビルドしようとすると
`clang: error: unsupported option '-fopenmp'`
などと怒られる。
[Homebrew]({{< relref "homebrew.md" >}})で
`gcc` か `llvm` をインストールし、
そっちを使ってビルドするように `~/.R/Makevars` で指定する:

```make
# gcc
CC=/usr/local/bin/gcc-8
CXX=/usr/local/bin/g++-8

# llvm
CC=/usr/local/opt/llvm/bin/clang
CXX=/usr/local/opt/llvm/bin/clang++
```

`CPPFLAGS`, `CFLAGS`, `CXXFLAGS`, `LDFLAGS` なども設定可能。

https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Using-Makevars
