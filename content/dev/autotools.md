+++
title = 'autoconf, automake'
tags = ["package"]
[menu.main]
  parent = "dev"
+++

環境に合わせて [Makefile]({{< relref "make.md" >}}) を作る仕組み。
以下の3つのツールを組み合わせて使う。

- <https://www.gnu.org/software/autoconf/>
- <https://www.gnu.org/software/automake/>
- <https://www.gnu.org/software/libtool/>

[CMake]({{< relref "cmake.md" >}}) のほうがモダンで高速。


## Commands

`autoscan`
:   指定したディレクトリ(指定しなければカレント)のソースコードを読んで
    `configure.ac` の雛形となる `configure.scan` を作る。
    既存の `configure.ac` もチェックするらしいが、その内容をすべて
    `configure.scan` に反映してくれるわけではなさそうなので
    そのまま上書きしてはダメっぽい。

`aclocal`
:   `configure.ac` を読んで `aclocal.m4` を作る

`automake`
:   `Makefile.am` と `configure.ac` から `Makefile.in` を作る

`autoconf`
:   `configure.ac` と `aclocal.m4` から `configure` を作る

`autoreconf`
:   上記のツールをいい感じに繰り返し呼び出して各種ファイルを更新

## 大まかな流れ

1.  `configure.scan` の雛形を自動生成し、
    `configure.ac` に名前変更:

        autoscan
        mv configure.scan configure.ac

1.  `configure.ac` を適宜編集
1.  `Makefile.am` を作る
1.  その2つのファイルから、自動的にその他のファイルを生成:

        autoreconf --install

1.  できあがった `configure` ファイルを試してみる:

        ./configure --help
        ./configure
        make

1.  `configure.ac` や `Makefile.am` を変更したら
    `autoreconf` で反映させる、を繰り返す

## `configure.ac`

<https://www.gnu.org/software/autoconf/manual/html_node/>

-   `configure.ac` の基本構造:
    <https://www.gnu.org/software/autoconf/manual/html_node/Autoconf-Input-Layout.html>
-   標準マクロ:
    <https://www.gnu.org/software/autoconf/manual/html_node/Autoconf-Macro-Index.html>
-   M4マクロ:
    <https://www.gnu.org/software/autoconf/manual/html_node/M4-Macro-Index.html>

Gitのタグをバージョン番号として取り込む:

    AC_INIT([MyApp], m4_esyscmd([git describe --tags | tr -d '\n']))

## `Makefile.am`

<https://www.gnu.org/software/automake/manual/html_node/>

-   マクロ:
    <https://www.gnu.org/software/automake/manual/html_node/Macro-Index.html>
-   変数:
    <https://www.gnu.org/software/automake/manual/html_node/Variable-Index.html>

インストールするファイルと場所を指定する変数:

    bin_PROBRAMS = beer
    bin_SCRIPTS = beer.sh
    lib_LIBRARIES = libbeer.a
    include_HEADERS = beer.h

ビルドするのに必要な情報をターゲットごとに指定する変数。
`target_ARGNAME` のような形をとる:

    beer_SOURCES = main.cpp
    beer_CPPFLAGS = -DNDEBUG
    beer_CXXFLAGS = -O3
    libbeer_a_SOURCES = lib.cpp

`Makefile` 全体に関わる変数。
ただし上記のターゲット特異的変数に上書きされる:

    AM_CPPFLAGS = -Wall -Wextra
    AM_CXXFLAGS = -O2

ユーザーが指定する `CPPFLAGS` は
`beer_CPPFLAGS` や `AM_CPPFLAGS` を上書きせず、
後ろに並べて使用される。

<https://www.gnu.org/software/automake/manual/html_node/Flag-Variables-Ordering.html>
