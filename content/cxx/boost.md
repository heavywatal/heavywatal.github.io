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

### パッケージマネージャで

[Homebrew]({{< relref "mac/homebrew.md" >}})/Linuxbrew
で最新版を簡単にインストールできる。
オプションは適当に:

    % brew install boost --c++11 --without-single

`--layout=tagged` でビルドされるため、
リンクするときは末尾に `-mt` が必要になる。


### ソースから

- http://www.boost.org/doc/libs/release/more/getting_started/unix-variants.html
- http://www.boost.org/build/
- https://boostjp.github.io/howtobuild.html

1.  <http://www.boost.org/users/download/> から最新ソースを入手して展開。
    ```
    % wget -O- https://downloads.sourceforge.net/boost/boost_1_63_0.tar.bz2 | tar xj
    % cd boost_1_63_0/
    ```

1.  ビルドすべきライブラリを考える `./bootstrap.sh --show-libraries`

1.  適当なオプションを与えて `bootstrap.sh` を実行:
    ```sh
    % ./bootstrap.sh --help
    % ./bootstrap.sh --without-icu --with-libraries=coroutine2,filesystem,graph,iostreams,program_options,serialization,system,test
    ```
    設定が `project-config.jam` に書き出され、
    `b2` がビルドされる。 `./b2 --help`

1. `~/user-config.jam` に [ツールセットを定義]
    (http://www.boost.org/build/doc/html/bbv2/reference/tools.html)。
    `darwin`はMac-gcc用:
    ```
    using gcc : 14 : g++-6 : <compileflags>-fPIC <cxxflags>-std=c++14 ;
    using darwin : 14 : g++-6 : <compileflags>-fPIC <cxxflags>-std=c++14 ;
    using clang : 14 : clang++ : <compileflags>-fPIC <cxxflags>-std=c++14 -stdlib=libc++ <linkflags>-stdlib=libc++ ;
    ```
    gccとclangの両方から使える統一ライブラリを作るのは難しいらしいので、
    それぞれのコンパイラで別々にビルドしてインストールする。

1.  システム標準zlibをリンクしようとしてエラーになるような場合は、
    [zlib公式](http://zlib.net/)からソースを落として展開し、
    [一緒にビルドされるように]
    (http://www.boost.org/doc/libs/release/libs/iostreams/doc/installation.html)
    `ZLIB_SOURCE`をフルパス指定する。
    ```sh
    % wget -O- http://zlib.net/zlib-1.2.8.tar.gz | tar xz -C ${HOME}/tmp/build
    % export ZLIB_SOURCE=${HOME}/tmp/build/zlib-1.2.8
    ```

1.  ツールセットを指定してビルド:
    ```sh
    % ./b2 -j2 toolset=gcc-14 link=static,shared runtime-link=shared threading=multi variant=release --layout=tagged --build-dir=../b2gcc --stagedir=stage/gcc stage
    % ./b2 -j2 toolset=darwin-14 link=static,shared runtime-link=shared threading=multi variant=release --layout=tagged --build-dir=../b2gcc --stagedir=stage/gcc stage
    % ./b2 -j2 toolset=clang-14 link=static,shared runtime-link=shared threading=multi variant=release --layout=tagged --build-dir=../b2clang --stagedir=stage/clang stage
    ```

1.  prefixを指定してインストール:
    ```sh
    % ./b2 -j2 toolset=gcc-14 link=static,shared runtime-link=shared threading=multi variant=release --layout=tagged --build-dir=../b2gcc --stagedir=stage/gcc --prefix=${HOME}/local install
    ```
    あるいは手動でインストール:
    ```sh
    % rsync -auv stage/gcc/ ~/local/boost-gcc
    % rsync -auv stage/clang/ ~/local/boost-clang
    % rsync -auv boost ~/local/include
    ```
