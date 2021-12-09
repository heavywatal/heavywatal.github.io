+++
title = 'gcc'
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

https://gcc.gnu.org/

## Usage

直接コマンドを実行するなら:

    g++ -g -Wall -Wextra -O3 main.cpp -o a.out

[make]({{< relref "make.md" >}}) や [CMake]({{< relref "cmake.md" >}})
などのツールを使ってビルドすることが多く、
ターミナルからコンパイラを直接実行することはほぼ無い。

### Options

- https://gcc.gnu.org/onlinedocs/gcc/Invoking-GCC.html
- https://gcc.gnu.org/onlinedocs/gcc/Option-Summary.html

#### [出力オプション](https://gcc.gnu.org/onlinedocs/gcc/Overall-Options.html)

`-c`
:   コンパイルするだけでリンクしない

`-o {file}`
:   出力先のファイル名を指定

#### [C/C++オプション](https://gcc.gnu.org/onlinedocs/gcc/C-Dialect-Options.html)

`-std=c++11`
:   2011年のISO標準でコンパイルする。

`-std=c++14`
:   2014年のISO標準でコンパイルする。
    g++-6 ではこれのGNU方言である `gnu++14` がデフォルトになった。

#### [警告オプション](https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html)

`-Wall`
:   基本的な警告

`-Wextra`
:   Wallより厳しい警告

`-Werror`
:   警告をエラー扱いにする

#### [デバッグオプション](https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html)

`-g`
:   gdb で使えるデバッグ情報を埋め込む

#### [最適化オプション](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html)

`-O1`
:   軽い最適化

`-O2`
:   プログラムサイズを大きくしない最適化をすべて実行

`-O3`
:   インライン展開なども行う

`-Os`
:   コードサイズが最小になるように

#### [プリプロセッサオプション](https://gcc.gnu.org/onlinedocs/gcc/Preprocessor-Options.html)

`-D {name}`, `-D {name=definition}`
:   ソースコードで `#define` するのと同じようにマクロを定義する

`-I {dir}`, `-isystem {dir}`
:   `#include <***>` のインクルードパスの先頭に `{dir}` を追加。

`-iquote {dir}`
:   `#include "***"` のインクルードパスの先頭に `{dir}` を追加。

`-M`, `-MM`
:   ファイルの依存関係を書き出す。
    前者はシステムヘッダーを含み、後者は含まない。

#### [リンクオプション](https://gcc.gnu.org/onlinedocs/gcc/Link-Options.html)

`-l {library}`
:   リンクするライブラリを指定する。
    ファイル名が `libsfmt.a` のときは `-lsfmt`

`-L {dir}`
:   ライブラリのサーチパスを追加する


## Installation

https://gcc.gnu.org/install/

### Mac

Macにある `/usr/bin/gcc` は
`gcc` の顔をした `clang` なので、
本物の `gcc` が欲しいときは別途インストールが必要。
[Homebrew]({{< relref "homebrew.md" >}}) で入れるのが楽チン:

```sh
brew install gcc
```

### Ubuntu

元から入ってるけど最新版をソースからインストールしたい。

1.  gmp, mpfr, mpc をインストール:

        sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev

1.  `--with-system-zlib` でビルドするには zlib が必要:

        sudo apt-get install zlib1g-dev

1.  システムに既に入ってるやつの `configure` オプションを
    `gcc -v` 見るでみて、それを参考に `configure`

    -   `--enable-languages` は使う言語だけに
    -   `--prefix` は `/usr/local`
    -   既存のバージョンとごっちゃになるのを避けるため
        `--program-suffix` を設定
    -   他言語サポートは不要なので `--disable-nls`
    -   `--disable-multilib` を付けないと
        `fatal error: gnu/stubs-32.h: No such file or directory. stub-32.h`
        などと怒られて止まる

    ```
    get -O- http://ftp.gnu.org/gnu/gcc/gcc-4.9.1/gcc-4.9.1.tar.bz2 | tar xj
    kdir build
    d build/
    ./gcc-4.9.1/configure --enable-languages=c,c++,go,fortran --prefix=/usr/local --program-suffix=-4.9 --enable-shared --enable-linker-build-id --without-included-gettext --enable-threads=posix --disable-nls --with-sysroot=/ --enable-clocale=gnu --enable-libstdcxx-debug --enable-libstdcxx-time=yes --enable-gnu-unique-object --disable-libmudflap --enable-plugin --with-system-zlib --disable-browser-plugin --enable-gtk-cairo --with-arch-directory=amd64 --enable-multiarch --disable-werror --with-abi=m64 --disable-multilib --with-tune=core2 --enable-checking=release --build=x86_64-linux-gnu --host=x86_64-linux-gnu --target=x86_64-linux-gnu
    ```

1.  コンパイル。ただしデフォルトではコンパイラオプションが `-g O2`
    となっていて膨大なデバッグシンボルがついた状態で出来上がってしまうらしいので、
    オプションを上書きしてそれを回避する:

        make BOOT_CFLAGS="-O2" CFLAGS="-O2" CXXFLAGS="-O2" bootstrap

    `fatal error: bits/predefs.h: No such file or directory`
    って怒られたら:

        sudo apt-get install libc6-dev-i386

    `/usr/bin/ld: cannot find crti.o: No such file or directory`
    って怒られたら:

        sudo ln -s /usr/lib/x86_64-linux-gnu /usr/lib64

1.  `sudo make install` でインストール。の前に
    `make DESTDIR=${HOME}/tmp/gcc-test install`
    などとしてテストするといいらしい。
    ディスクレスクラスタの場合、子ノードにもインストール:
    `sudo make DESTDIR=/nfsroot install`
