+++
title = 'Boost'
subtitle = "ほぼ標準C++ライブラリ"
[menu.main]
  parent = "cxx"
+++

-   <http://www.boost.org/>
-   <http://www.boost.org/doc/libs/release/>

## Installation

<http://boostjp.github.io/howtobuild.html>

<http://www.boost.org/doc/libs/release/more/getting_started/unix-variants.html>

普通の configure と `make` じゃないので混乱するけど、まあまあ相同な手順。

1.  <http://sourceforge.net/projects/boost/> から最新ソースを入手して展開。
    とりあえず `boost-jam` とか `boost-build` とかは無視して `boost` 本体のみで結構:

        % wget -O- http://downloads.sourceforge.net/boost/boost_1_59_0.tar.bz2 | tar xj
        % cd boost_1_59_0/

2.  ヘルプを見る:

        % ./bootstrap.sh --help

3.  ビルドすべきライブラリを考える:

        % ./bootstrap.sh --show-libraries

4.  適当なオプションを与えて bootstrap.sh を実行

    -   b2 がビルドされる
    -   b2 に渡すオプションが書かれた
        `project-config.jam` が生成される

    <!-- -->

        % ./bootstrap.sh --without-icu --with-libraries=filesystem,graph,iostreams,program_options,serialization,system,test

5.  ヘルプを見る:

        % ./b2 --help

6.  コンパイラやオプションを指定してビルド:

        % ./b2 -j2 toolset=gcc-5 cxxflags="-std=c++11" link=static runtime-link=shared threading=multi variant=release --layout=system --stagedir=stage_gcc stage 2>&1 | tee stage.log
        % ./b2 -j2 toolset=clang cxxflags="-std=c++11 -stdlib=libc++" linkflags="-stdlib=libc++" link=static runtime-link=shared threading=multi variant=release --layout=system --stagedir=stage_clang stage 2>&1 | tee stage.log

7.  古いやつがあれば消しておく:

        % rm -rf ~/local/boost*

8.  インストール (`prefix` によっては要 `sudo`):

        % ./b2 -j2 toolset=gcc-5 cxxflags="-std=c++11" link=static runtime-link=shared threading=multi variant=release --layout=system --stagedir=stage_gcc install --prefix=${HOME}/local/boost-gcc 2>&1 | tee install.log
        % ./b2 -j2 toolset=clang cxxflags="-std=c++11 -stdlib=libc++" linkflags="-stdlib=libc++" link=static runtime-link=shared threading=multi variant=release --layout=system --stagedir=stage_clang install --prefix=${HOME}/local/boost-clang 2>&1 | tee install.log

9.  ヘッダーは1か所でもよさそう:

        % rsync -auv boost ~/local/include

------------------------------------------------------------------------

Macなら [Homebrew]({{< relref "mac/homebrew.md" >}}) でもインストールできる:

    % brew install boost --c++11 --without-single

ただし `--layout=tagged` になっているため、
リンクするときは末尾に `-mt` が必要になる。

Ubuntuなら [ppa:boost-latest/ppa](https://launchpad.net/~boost-latest/+archive/ppa)
リポジトリを加えて `libboost*-dev` を適当にインストール:

    % sudo add-apt-repository ppa:boost-latest/ppa
    % sudo apt-get update
    % apt-cache search libboost.*-all-dev
