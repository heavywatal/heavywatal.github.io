+++
title = 'R環境設定'
tags = ["r"]
weight = -95
+++

<https://cran.r-project.org/manuals.html>

## インストール

1.  OSのソフトウェア・アップデートを基本的に全て適用して再起動。
1.  ファイル名の末尾(`.pdf` とか `.png` とか)の[拡張子を常時表示する](https://duckduckgo.com/?q=拡張子+表示)ようにOSを設定。
1.  <https://cloud.r-project.org/>
    から最新版の **R本体** をダウンロードしてインストール。
    OK連打のデフォルト設定で。
    古いものが既に入っている場合は念のため削除してから。
    - <iconify-icon inline icon="bi:windows"></iconify-icon>
      [Windows → base](https://cloud.r-project.org/bin/windows/base) → `R-*-win.exe`
    - <iconify-icon inline icon="bi:apple"></iconify-icon>
      [Mac](https://cloud.r-project.org/bin/macosx/)
      → `R-*-arm64.pkg` (Apple Silicon) or `R-*-x86_64.pkg` (Intel)
1.  <https://posit.co/download/rstudio-desktop/>
    から最新版の **RStudio** をダウンロードしてインストール。
    古いものが既に入っている場合は念のため削除してから。
1.  RStudioを起動。初回は「開発元の不明なアプリ」を許可するような操作が必要かも。

> [!TIP]
Macなら [Homebrew]({{< relref "homebrew.md" >}}) で
`brew install r-app rstudio` のように入れるのが楽チン。
Caskじゃない `brew install r` のほうだとバイナリ版パッケージが使えなくて毎回ソースからビルドさせられるので大変。

<https://cran.r-project.org/doc/manuals/R-admin.html>


### What’s New?

<https://cran.r-project.org/doc/manuals/r-release/NEWS.html>

- 4.6
  - `%notin%` operator
  - C++20 default
- 4.5
  - [`penguins`](https://stat.ethz.ch/R-manual/R-patched/library/datasets/html/penguins.html)
  - `install.packages()` in parallel
  - `base::use()`
  - C23 default
- 4.4
  - `%||%` operator
- 4.3
  - C17 and C++17 default; C++23 support
  - extraction with the placeholder `_`
- 4.2
  - placeholder `_` for a named argument
- 4.1
  - shorthand function `\(x) x + 1`
  - native pipe operator `|>`
  - C++14 default
- 4.0
  - `StringsAsFactors = FALSE` by default
  - color palettes: R4, Okabe-Ito, etc.
  - raw character strings `r"(...)"`
  - `tools::R_user_dir()`
  - C++20 support
- 3.6
  - C++11 default


## 起動オプション

ワークスペースの自動保存や自動復帰は危険なので切っておく。
R.app や RStudio から使う場合はメニューから環境設定みたいなやつを開いて設定。
シェルから使う場合は例えば以下のようなエイリアスを設定する。

```sh
alias r='R --quiet --no-save --no-restore-data'
```

詳しくは `R --help` または
https://cran.r-project.org/doc/manuals/R-intro.html#Invoking-R


## パッケージのサーチパス

[`?.libPaths`](https://stat.ethz.ch/R-manual/R-patched/library/base/html/libPaths.html)

`.libPaths()` でパッケージのインストール先候補一覧を取得できる。
`install.packages(pkgs, lib, ...)`
の `lib = ` オプションを指定しない場合にこれらが参照される。
また `library()` によるパッケージ読み込みもこれらのパスから。

`.libPaths("newpath")` のように任意のパスを追加することもできるが、
`.Renviron` ファイルなどで環境変数から自動的に設定するほうが便利。
環境変数には優先順位があるので例えば次のように使い分けられる:

- `R_LIBS`: プロジェクトごとの一時的な設定
- `R_LIBS_USER`: ユーザーが常に使いたい設定。空の場合の規定値はOSによって異なる:
    - <iconify-icon inline icon="bi:ubuntu"></iconify-icon>
      Linux: `~/R/%p-library/%v`
    - <iconify-icon inline icon="bi:apple"></iconify-icon>
      Mac: `~/Library/R/%a/%v/library`
    - <iconify-icon inline icon="bi:windows"></iconify-icon>
      Windows: `${LOCALAPPDATA}/R/win-library/%v`\
      (`%LOCALAPPDATA%` は大概 `C:\Users\${username}\AppData\Local`)
- `R_LIBS_SITE`: 管理者が全ユーザーに使わせたい設定
    - R内から `.Library.site` で参照可能
    - 空の場合は `$R_HOME/site-library` になる

ここで設定するパスには `%v` といった記号でRのバージョン情報などを含めることも可能。
古いRでインストールしたパッケージを新しいRで使おうとすると
`package ‘***’ was installed before R x.y.0: please re-install it`
などと怒られるので、バージョン番号を入れておいたほうがいい。
`path.expand()` も適用されるので `~/.R/library` のようなチルダも展開される。

> [!WARNING]
設定がファイルから読み込まれても、
**当該ディレクトリが存在しないと認識されず自動生成もされない**。

Rと一緒についてくる標準パッケージのインストール先は `.Library` で参照可能。
何も設定しないで使うとほかのパッケージもそこに入ってしまう場合があってあんまりよろしくない。


## 環境変数

https://cran.r-project.org/doc/manuals/R-admin.html#Environment-variable-index

[`?"environment variables"`](https://stat.ethz.ch/R-manual/R-patched/library/base/html/EnvVar.html)

```r
Sys.getenv()
```


## .Renviron

[`?Startup`](https://stat.ethz.ch/R-manual/R-patched/library/base/html/Startup.html)

R起動時に読み込まれ、環境変数を設定するファイル。
Rスクリプトではなく、シェルスクリプトっぽい代入式で書く。
例 (<https://github.com/heavywatal/dotfiles/blob/master/.R/.Renviron>):

```sh
R_LIBS_USER=${HOME}/.R/library/%v
R_ENVIRON_USER=${HOME}/.R/.Renviron
R_PROFILE_USER=${HOME}/.R/.Rprofile
R_HISTFILE=${HOME}/.R/.Rhistory
_R_CHECK_SYSTEM_CLOCK_=FALSE
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
普通の対話環境でのみ有効になるよう `if (interactive())` で包んでおいたほうが安心。

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

[`?options`](https://stat.ethz.ch/R-manual/R-patched/library/base/html/options.html)
で項目の一覧を見られる。

ほかの人とスクリプトをやり取りする場合など、
実行結果が設定依存で変わっては困るので、
そういう本質的なものはいじらずに、
表示関連の項目だけ変えるに留めたほうがよい。

`warn=1`
: 警告レベルの設定。
  デフォルト(`warn=0`)では、警告があっても計算は滞り無く進行し、
  最後に "There were 50 or more warnings (use warnings() to see the first 50)"
  などと軽く表示されるだけなので、
  見落としたりして後々大変なバグ取り作業に発展する恐れがある。
  警告が発生するごとに警告文を表示する(`warn=1`)か、
  エラー扱いにして計算をストップするようにしておく(`warn=2`)ことでそれを回避できる。
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

`defaultPackages`
: 起動時(`.First()` 実行よりは後)に自動で読み込むパッケージを指定する。
: 環境変数 `R_DEFAULT_PACKAGES` からも変更可能。
: デフォルトは datasets, utils, grDevices, graphics, stats, methods.
: `conflicted` と `tidyverse` を加えて横着したいところだが、
  そうすると肝心の衝突チェックが機能しなくなる。
  前者だけを加えた上で次のようにフックを設定すれば自動読み込みでチェック有効:
  ```r
  setHook(packageEvent("conflicted", "attach"), \(...) library(tidyverse))
  ```
  ちなみに `conflicted` は knitr chunk 内や `withr::local_package()` では動かない。
  ["conflicted is designed specifically for use in interactive sessions"](https://github.com/r-lib/conflicted/issues/88#issuecomment-1445383091)
  とのこと。


## Advanced

https://cran.r-project.org/doc/manuals/R-admin.html#Add_002don-packages


### 開発者ツール

C, C++, Fortran で書かれたソースコードをビルドするためのツール。
普通はビルド済みのバイナリ版パッケージをインストールするだけので不要。
[Stan]({{< relref "stan.md" >}})でモデルを書くとか、パッケージの開発最新版を使いたいとか、
そういうときに必要になる。

- <iconify-icon inline icon="bi:apple"></iconify-icon>
  Mac:
  - [**Command Line Tools**](https://duckduckgo.com/?q=command+line+tools):
    ターミナルで `xcode-select --install` を実行。
    Xcode環境は不要。
  - [**`gfortran-*.pkg`**](https://mac.r-project.org/tools/) を落として入れる。
    ほかの方法で入れたものを使うのは難しい。
- <iconify-icon inline icon="bi:windows"></iconify-icon>
  Windows: [**Rtools**](https://cloud.r-project.org/bin/windows/Rtools/)
  (R本体のバージョンに合わせる)


### `~/.R/Makevars`

パッケージをソースコードからビルドするときの設定。
不具合の元なので基本的には作らない。

<https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Using-Makevars>
