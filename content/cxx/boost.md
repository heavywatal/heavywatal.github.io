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

- http://www.boost.org/doc/libs/release/more/getting_started/unix-variants.html
- http://www.boost.org/build/
- https://boostjp.github.io/howtobuild.html

普通の `configure` と `make` じゃないので混乱するけど、まあまあ相同な手順。
gccとclangの両方から使える統一ライブラリを作るのは難しいらしいので、
それぞれのコンパイラで別々にビルドしてインストールする。

1.  <https://sourceforge.net/projects/boost/> から最新ソースを入手して展開。
    とりあえず `boost-jam` とか `boost-build` とかは無視して `boost` 本体のみで結構:
    ```
    % wget -O- https://downloads.sourceforge.net/boost/boost_1_62_0.tar.bz2 | tar xj
    % cd boost_1_62_0/
    ```

2.  ヘルプを見る `./bootstrap.sh --help`

3.  ビルドすべきライブラリを考える `./bootstrap.sh --show-libraries`

4.  適当なオプションを与えて `bootstrap.sh` を実行:
    ```sh
    % ./bootstrap.sh --without-icu --with-libraries=coroutine2,filesystem,graph,iostreams,program_options,serialization,system,test
    ```

    `b2` がビルドされ、
    `b2` に渡すオプションが書かれた `project-config.jam` が生成される。

5.  ヘルプを見る `./b2 --help`

6. `~/user-config.jam` に [ツールセットを定義]
    (http://www.boost.org/build/doc/html/bbv2/reference/tools.html)。
    `darwin`はMac-gcc用:
    ```
    using gcc : 14 : g++-6 : <compileflags>-fPIC <cxxflags>-std=c++14 ;
    using darwin : 14 : g++-6 : <compileflags>-fPIC <cxxflags>-std=c++14 ;
    using clang : 14 : clang++ : <compileflags>-fPIC <cxxflags>-std=c++14 -stdlib=libc++ <linkflags>-stdlib=libc++ ;
    ```

7.  システム標準zlibをリンクしようとしてエラーになることがあるので、
    [zlib公式](http://zlib.net/)からソースを落として展開し、
    [一緒にビルドされるように]
    (http://www.boost.org/doc/libs/release/libs/iostreams/doc/installation.html)
    `ZLIB_SOURCE`をフルパス指定する。
    ```
    % wget -O- http://zlib.net/zlib-1.2.8.tar.gz | tar xz -C ${HOME}/tmp/build
    % export ZLIB_SOURCE=${HOME}/tmp/build/zlib-1.2.8
    ```

8.  ツールセットを指定してビルド:
    ```
    % ./b2 -j2 toolset=gcc-14 link=static,shared runtime-link=shared threading=multi variant=release --layout=tagged --build-dir=../b2gcc --stagedir=stage/gcc stage
    % ./b2 -j2 toolset=darwin-14 link=static,shared runtime-link=shared threading=multi variant=release --layout=tagged --build-dir=../b2gcc --stagedir=stage/gcc stage
    % ./b2 -j2 toolset=clang-14 link=static,shared runtime-link=shared threading=multi variant=release --layout=tagged --build-dir=../b2clang --stagedir=stage/clang stage
    ```

9.  手動でインストール:
    ```
    % rsync -auv stage/gcc/ ~/local/boost-gcc
    % rsync -auv stage/clang/ ~/local/boost-clang
    % rsync -au boost ~/local/include
    ```

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
