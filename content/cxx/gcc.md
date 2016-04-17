+++
title = 'gcc'
[menu.main]
  parent = "cxx"
+++

[<http://gcc.gnu.org/>](http://gcc.gnu.org/)

## Usage

直接コマンドを実行するなら:

    % g++ -g -Wall -Wextra -O3 main.cpp -o a.out

{{%div class="note"%}}
[/dev/make]({{< relref "dev/make.md" >}})

コンパイルの度にオプションとかファイルを打ち込むのは面倒なので
`Makefile` にルールを記述しといて `make` コマンドでビルドする。
{{%/div%}}

### Options

<http://gcc.gnu.org/onlinedocs/gcc-4.9.0/gcc/Invoking-GCC.html>

**出力オプション** <http://gcc.gnu.org/onlinedocs/gcc-4.9.0/gcc/Overall-Options.html>
:   

    `-c`
    :   コンパイルするだけでリンクしない

    `-o {file}`
    :   出力先のファイル名を指定

**C/C++オプション** <http://gcc.gnu.org/onlinedocs/gcc-4.9.0/gcc/C-Dialect-Options.html>
:   

    `-std=c++11`
    :   2011年のISO標準でコンパイルする。
        早くこれがデフォルトになればいいのに...。
        2014年予定の新標準も `-std=c++1y` で試験的に有効にできる。
        <http://gcc.gnu.org/projects/cxx1y.html>

**警告オプション** <http://gcc.gnu.org/onlinedocs/gcc-4.9.0/gcc/Warning-Options.html>
:   

    `-Wall`
    :   基本的な警告

    `-Wextra`
    :   Wallより厳しい警告

    `-Werror`
    :   警告をエラー扱いにする

**デバッグオプション** <http://gcc.gnu.org/onlinedocs/gcc-4.9.0/gcc/Debugging-Options.html>
:   

    `-g`
    :   gdb で使えるデバッグ情報を埋め込む

**最適化オプション** <http://gcc.gnu.org/onlinedocs/gcc-4.9.0/gcc/Optimize-Options.html>
:   

    `-O1`
    :   軽い最適化

    `-O2`
    :   プログラムサイズを大きくしない最適化をすべて実行

    `-O3`
    :   インライン展開なども行う

    `-Os`
    :   コードサイズが最小になるように

**プリプロセッサオプション** <http://gcc.gnu.org/onlinedocs/gcc-4.9.0/gcc/Preprocessor-Options.html>
:   

    `-D {name}`, `-D {name=definition}`
    :   ソースコードで `#define` するのと同じようにマクロを定義する

    `-I {dir}`, `-isystem {dir}`
    :   `#include <***>` のインクルードパスの先頭に `{dir}` を追加。

    `-iquote {dir}`
    :   `#include "***"` のインクルードパスの先頭に `{dir}` を追加。

    `-M`, `-MM`
    :   ファイルの依存関係を書き出す。
        前者はシステムヘッダーを含み、後者は含まない。

**リンクオプション** <http://gcc.gnu.org/onlinedocs/gcc-4.9.0/gcc/Link-Options.html>
:   

    `-l {library}`
    :   リンクするライブラリを指定する。
        ファイル名が `libsfmt.a` のときは `-lsfmt`

    `-L {dir}`
    :   ライブラリのサーチパスを追加する

## Installation

<http://gcc.gnu.org/install/>

### Mac

Xcode についてくる `gcc-4.2` は
`gcc` の顔をした `clang` なので、
本当の `gcc` が欲しいときは別途インストールが必要。
[Homebrew]({{< relref "mac/homebrew.md" >}})` か `[MacPorts]({{< relref "mac/macports.md" >}}) で入れるのが楽チン:

    % brew install gcc --without-multilib

### Ubuntu

元から入ってるけど最新版をソースからインストールしたい。

1.  gmp, mpfr, mpc をインストール:

        % sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev

2.  `--with-system-zlib` でビルドするには zlib が必要:

        % sudo apt-get install zlib1g-dev

3.  システムに既に入ってる `gcc` の `configure` オプションを見る:

        % gcc -v
        Using built-in specs.
        COLLECT_GCC=gcc
        COLLECT_LTO_WRAPPER=/usr/lib/gcc/x86_64-linux-gnu/4.8/lto-wrapper
        Target: x86_64-linux-gnu
        Configured with: ../src/configure -v --with-pkgversion='Ubuntu 4.8.2-19ubuntu1' --with-bugurl=file:///usr/share/doc/gcc-4.8/README.Bugs --enable-languages=c,c++,java,go,d,fortran,objc,obj-c++ --prefix=/usr --program-suffix=-4.8 --enable-shared --enable-linker-build-id --libexecdir=/usr/lib --without-included-gettext --enable-threads=posix --with-gxx-include-dir=/usr/include/c++/4.8 --libdir=/usr/lib --enable-nls --with-sysroot=/ --enable-clocale=gnu --enable-libstdcxx-debug --enable-libstdcxx-time=yes --enable-gnu-unique-object --disable-libmudflap --enable-plugin --with-system-zlib --disable-browser-plugin --enable-java-awt=gtk --enable-gtk-cairo --with-java-home=/usr/lib/jvm/java-1.5.0-gcj-4.8-amd64/jre --enable-java-home --with-jvm-root-dir=/usr/lib/jvm/java-1.5.0-gcj-4.8-amd64 --with-jvm-jar-dir=/usr/lib/jvm-exports/java-1.5.0-gcj-4.8-amd64 --with-arch-directory=amd64 --with-ecj-jar=/usr/share/java/eclipse-ecj.jar --enable-objc-gc --enable-multiarch --disable-werror --with-arch-32=i686 --with-abi=m64 --with-multilib-list=m32,m64,mx32 --with-tune=generic --enable-checking=release --build=x86_64-linux-gnu --host=x86_64-linux-gnu --target=x86_64-linux-gnu
        Thread model: posix
        gcc version 4.8.2 (Ubuntu 4.8.2-19ubuntu1)

4.  それを参考にして `configure`

    -   `--enable-languages` は使う言語だけに
    -   `--prefix` は `/usr/local`
    -   既存のバージョンとごっちゃになるのを避けるため
        `--program-suffix` を設定
    -   他言語サポートは不要なので `--disable-nls`
    -   `--disable-multilib` を付けないと
        `fatal error: gnu/stubs-32.h: No such file or directory. stub-32.h`
        などと怒られて止まる

    <!-- -->

        % wget -O- http://ftp.gnu.org/gnu/gcc/gcc-4.9.1/gcc-4.9.1.tar.bz2 | tar xj
        % mkdir build
        % cd build/
        % ../gcc-4.9.1/configure --enable-languages=c,c++,go,fortran --prefix=/usr/local --program-suffix=-4.9 --enable-shared --enable-linker-build-id --without-included-gettext --enable-threads=posix --disable-nls --with-sysroot=/ --enable-clocale=gnu --enable-libstdcxx-debug --enable-libstdcxx-time=yes --enable-gnu-unique-object --disable-libmudflap --enable-plugin --with-system-zlib --disable-browser-plugin --enable-gtk-cairo --with-arch-directory=amd64 --enable-multiarch --disable-werror --with-abi=m64 --disable-multilib --with-tune=core2 --enable-checking=release --build=x86_64-linux-gnu --host=x86_64-linux-gnu --target=x86_64-linux-gnu

5.  コンパイル。ただしデフォルトではコンパイラオプションが `-g O2`
    となっていて膨大なデバッグシンボルがついた状態で出来上がってしまうらしいので、
    オプションを上書きしてそれを回避する。:

        % make BOOT_CFLAGS="-O2" CFLAGS="-O2" CXXFLAGS="-O2" bootstrap

    `fatal error: bits/predefs.h: No such file or directory`
    って怒られたら:

        % sudo apt-get install libc6-dev-i386

    `/usr/bin/ld: cannot find crti.o: No such file or directory`
    って怒られたら:

        % sudo ln -s /usr/lib/x86_64-linux-gnu /usr/lib64

6.  インストール。の前に `make DESTDIR=$HOME/tmp/gcc-test install`
    などとしてテストするといいらしい。 :

        % sudo make install

        # ディスクレスクラスタの場合、子ノードにもインストール
        % sudo make DESTDIR=/nfsroot install
