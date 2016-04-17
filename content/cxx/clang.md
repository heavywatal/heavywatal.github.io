+++
title = 'clang / llvm'
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

C language family frontend for LLVM

<http://clang.llvm.org/>

## Installation

<http://clang.llvm.org/get_started.html>

### Mac

Xcode の Command Line Tools にほぼ最新版が含まれている:

    % clang -v
    Apple LLVM version 5.1 (clang-503.0.40) (based on LLVM 3.4svn)
    Target: x86_64-apple-darwin13.1.0
    Thread model: posix

### Ubuntu 12.04 (Binary)

1.  ダウンロードして展開:

        % wget -O- http://llvm.org/releases/3.3/clang+llvm-3.3-amd64-Ubuntu-12.04.2.tar.gz | tar xz

2.  適当なところに置いてパスを通す:

        % sudo mv clang+llvm-3.3-amd64-Ubuntu-12.04.2 /usr/local/llvm-3.3
        % sudo ln -s llvm-3.3 /usr/local/llvm
        % export PATH=/usr/local/llvm/bin:$PATH

3.  [libc++をインストールする]({{< relref "#libc" >}})

### Ubuntu 12.04 (Source)

1.  必要な道具を揃える ([Requirements](http://llvm.org/docs/GettingStarted.html#requirements)):

        % sudo apt-get install build-essential subversion texinfo

    `llvm/test/` のテストを実行したい場合は:

        % sudo apt-get install dejagnu tcl expect

2.  LLVM をダウンロード:

        % svn checkout http://llvm.org/svn/llvm-project/llvm/tags/RELEASE_33/final/ llvm

3.  その中の `tools/` に Clang を、
    `projects/` に Compiler-RT をダウンロード:

        % cd llvm/tools/
        % svn checkout http://llvm.org/svn/llvm-project/cfe/tags/RELEASE_33/final/ clang
        % cd ../projects/
        % svn checkout http://llvm.org/svn/llvm-project/compiler-rt/tags/RELEASE_33/final/ compiler-rt

4.  外に `build` ディレクトリを作って、そこに入る:

        % cd ../..
        % mkdir build
        % cd build

5.  オプションを確認してビルド:

        % ../llvm/configure --help
        % ../llvm/configure --enable-optimized --disable-assertions --enable-targets=x86_64 --enable-shared --enable-libffi
        % make
        % sudo make install

    {{%div class="note"%}}
裸で configure するとデバッグモードでビルドされる
    {{%/div%}}

6.  [libc++をインストールする]({{< relref "#libc" >}})

### `libc++`

<http://libcxx.llvm.org/>

1.  `libsupc++` の位置を確認:

        % echo | g++-4.8 -Wp,-v -x c++ - -fsyntax-only

2.  `libcxx` をダウンロード:

        % svn checkout http://llvm.org/svn/llvm-project/libcxx/tags/RELEASE_33/final/ libcxx

3.  一時ディレクトリに移ってビルド:

        % mkdir build-libcxx
        % cd build-libcxx
        % CC=clang CXX=clang++ cmake -G "Unix Makefiles" -DLIBCXX_CXX_ABI=libsupc++ -DLIBCXX_LIBSUPCXX_INCLUDE_PATHS="/usr/local/include/c++/4.8.1/;/usr/local/include/c++/4.8.1/x86_64-linux-gnu/" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local ../libcxx
        % make
        % sudo make install

4.  `Makefile` で以下のようなオプションと共に使う:

        CPPFLAGS := -isystem /usr/local/include/c++/v1
        CXXFLAGS := -std=c++11 -stdlib=libc++
        LD_FLAGS := -L/usr/local/lib
        LDLIBS := -lsupc++
