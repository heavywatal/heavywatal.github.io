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

### 使うとき

インストールした場所とリンクするライブラリをコンパイラに伝える必要がある。
コマンドを直打ちするなら:
```sh
clang++ -I${HOME}/local/include -L${HOME}/local/lib mysource.cpp -lboost_iostreams-mt -o a.out
```

Makefileの変数でいうと:
```make
CPPFLAGS = -I${HOME}/local/include
LDFLAGS = -L${HOME}/local/lib
LDLIBS = -lboost_iostreams-mt
```


## [math](http://www.boost.org/doc/libs/release/libs/math/doc/html/)

### [distribution](http://www.boost.org/doc/libs/release/libs/math/doc/html/dist.html)

確率分布に従った乱数生成はC++11から `<random>` でサポートされるようになったが、
確率密度関数(PDF)や累積密度関数(CDF)はまだ標準入りしてない。

```c++
// #include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/chi_squared.hpp>

namespace bmath = boost::math;

bmath::normal_distribution<> dist(mean, sd);

bmath::mean(dist);
bmath::median(dist);
bmath::standard_deviation(dist);

bmath::pdf(dist, x);
bmath::cdf(dist, x);
bmath::quantile(dist, p);
```

{{%div class="warning"%}}
右側の裾が欲しいときは精度を保つために `complement()` を使う。
(Rでいう `lower.tail=FALSE`)

```c++
// good
bmath::quantile(bmath::complement(dist, p));
// bad
bmath::quantile(dist, 1.0 - p)

// good
bmath::cdf(bmath::complement(dist, x));
// bad
1.0 - bmath::cdf(dist, x);
```
{{%/div%}}


## [iostreams](http://www.boost.org/doc/libs/release/libs/iostreams/doc/)

要ビルド＆リンク `-lboost_iostreams-mt`

### gzip 圧縮と展開

https://github.com/heavywatal/cxxwtils/blob/master/zfstream.hpp

```c++
#include <iostream>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

int main() {
    namespace bios = boost::iostreams;
    {
        bios::filtering_ostream ost;
        ost.push(bios::gzip_compressor());
        ost.push(bios::file_descriptor_sink("hello.txt.gz"));
        ost << "Hello world!";
    }
    {
        bios::filtering_istream ist;
        ist.push(bios::gzip_decompressor());
        ist.push(bios::file_descriptor_source("hello.txt.gz"));
        std::string buffer;
        std::getline(ist, buffer, '\0');
        std::cout << buffer << std::endl;
    }
}
```

- gzipフィルタはコンストラクタで渡してもよい。
- `file_descriptor` はfailビットが立つとすぐ例外を投げて
  "No such file or directory" などを知らせてくれるので便利。
  標準streamのような沈黙を求める場合は代わりに
  `std::ifstream` などを `push()` することも可能。


## [program_options](http://www.boost.org/doc/libs/release/doc/html/program_options.html)

要ビルド＆リンク `-lboost_program_options-mt`


## [coroutine2](http://www.boost.org/doc/libs/release/libs/coroutine2/doc/html/)

要ビルド＆リンク `-lboost_context-mt`

Pythonの`yield`みたいなことをC++でもできるようになる。

[Fibonacci generator on gist](https://gist.github.com/heavywatal/e9c4d705b5617e4fc6ea32452db18860)

{{%div class="warning"%}}
オブジェクトの寿命に注意。

- `yield`返しはmoveなので、
  次の処理で再利用するつもりならコピーコンストラクタ越しに新品を返す。
- generator的なものを返す関数を作ると、
  それを抜ける時に寿命を迎えるオブジェクトがあることに注意。
{{%/div%}}
