+++
title = 'Boost'
subtitle = "ほぼ標準C++ライブラリ"
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

-   <http://www.boost.org/>
-   <http://www.boost.org/doc/libs/release/>

## Installation

### ソースから

<http://boostjp.github.io/howtobuild.html>

<http://www.boost.org/doc/libs/release/more/getting_started/unix-variants.html>

普通の `configure` と `make` じゃないので混乱するけど、まあまあ相同な手順。
gccとclangの両方から使える統一ライブラリを作るのは難しいらしいので、
それぞれのコンパイラで別々にビルドしてインストールする。

1.  <http://sourceforge.net/projects/boost/> から最新ソースを入手して展開。
    とりあえず `boost-jam` とか `boost-build` とかは無視して `boost` 本体のみで結構:

        % wget -O- http://downloads.sourceforge.net/boost/boost_1_61_0.tar.bz2 | tar xj
        % cd boost_1_61_0/

2.  ヘルプを見る:

        % ./bootstrap.sh --help

3.  ビルドすべきライブラリを考える:

        % ./bootstrap.sh --show-libraries

4.  適当なオプションを与えて `bootstrap.sh` を実行:

    ```sh
    % ./bootstrap.sh --without-icu --with-libraries=filesystem,graph,iostreams,program_options,serialization,system,test
    ```

    `b2` がビルドされ、
    `b2` に渡すオプションが書かれた `project-config.jam` が生成される。

5.  ヘルプを見る:

        % ./b2 --help

6. `~/user-config.jam` にツールセットを定義:

        using gcc : 14 : g++-6 : <cxxflags>-std=c++14 -stdlib=libstdc++ <linkflags>-stdlib=libstdc++ ;
        using clang : 14 : clang++ : <cxxflags>-std=c++14 -stdlib=libc++ <linkflags>-stdlib=libc++ ;

7.  ツールセットを指定してビルド:

        % ZLIB_SOURCE=${HOME}/tmp/build/zlib-1.2.8 ./b2 -j2 toolset=gcc-14 link=static runtime-link=shared threading=multi variant=release --layout=system --stagedir=stage_gcc stage 2>&1 | tee stage.log
        % ./b2 -j2 toolset=clang-14 link=static runtime-link=shared threading=multi variant=release --layout=system --stagedir=stage_clang stage 2>&1 | tee stage.log

    {{%div class="note"%}}
新しめのHomebrew gccで使うとzlibのリンクでエラーが出るっぽい。
[zlib公式](http://zlib.net/)からソースを落として展開し、
[`ZLIB_SOURCE`をフルパスで指定する。]
(http://www.boost.org/doc/libs/1_61_0/libs/iostreams/doc/installation.html)
{{%/div%}}


8.  古いやつがあれば消しておく:

        % mv ~/local/boost* ~/local/include/boost ~/.Trash

9.  手動でインストール:

        % mv stage_gcc ~/local/boost-gcc
        % mv stage_clang ~/local/boost-clang
        % rsync -au boost ~/local/include


### パッケージマネージャで簡単インストール

Macなら [Homebrew]({{< relref "mac/homebrew.md" >}}) でもインストールできる:

    % brew install boost --c++11 --without-single

ただし `--layout=tagged` になっているため、
リンクするときは末尾に `-mt` が必要になる。
基本的にはclangからのみ利用可能。

Ubuntuなら [ppa:boost-latest/ppa](https://launchpad.net/~boost-latest/+archive/ppa)
リポジトリを加えて `libboost*-dev` を適当にインストール:

    % sudo add-apt-repository ppa:boost-latest/ppa
    % sudo apt-get update
    % apt-cache search libboost.*-all-dev
