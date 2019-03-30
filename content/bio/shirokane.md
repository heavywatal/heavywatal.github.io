+++
date = 2018-06-14T10:00:00+09:00
title = 'SHIROKANE'
subtitle = 'HGCスパコン'
tags = ["job"]
[menu.main]
  parent = "bio"
+++

https://supcom.hgc.jp/

## 利用開始

えくせる、ハンコ、郵送。

無料プランはファイル総数1000個という縛りがキツすぎて実質使えない。


## 環境整備

- [ハードウェア](https://supcom.hgc.jp/japanese/sys_const/system-main.html)
- [ソフトウェア](https://supcom.hgc.jp/internal/cgi/version_up_s3/select.cgi)
    - RHEL 6.9
    - `/usr/local/bin/cmake` 3.8.1
    - `/usr/local/bin/gcc` 4.9.3
    - `/usr/local/bin/git` 1.7.7.1
    - `/usr/local/package/boost/1.64.0/` (ただしC++11不可)
    - `/usr/local/package/gcc/7.3.0/`
    - `/usr/local/package/python/3.6.4/`
    - `/usr/local/package/r/3.5.0/` (ただしC++14不可)


## tumopp on R

http://heavywatal.github.io/rtumopp/

提供されている R-3.5 はバージョンとしては十分に新しいが、
古いコンパイラでビルドされているためC++14のライブラリを使えない。
そのため `/usr/local/package/gcc/7.3.0`
を使ってRをビルドするところから始める必要がある。

Linuxbrewでインストールした `cmake`, `gcc`, `glibc`
などがあると干渉してエラーを起こすので、
アンインストールするかPATHから外すなどしておく。

1.  `~/.bashrc` などで以下のように環境変数を設定する:
    ```
    module unload intel_env
    module unload gcc
    module load /usr/local/package/gcc/modulefiles/gcc/7.3.0
    export PATH=${HOME}/local/bin:${PATH}
    ```

1.  シェルを再起動してこれを反映: `exec $SHELL -l` 。
    `g++ -v` で `gcc version 7.3.0` が表示されることを確認。

1.  Rをソースコードからビルドしてインストール
    (configureオプションは `/usr/local/package/r/3.5.1/lib64/R/etc/Makeconf` を参考に):
    ```
    export GCC_PREFIX=/usr/local/package/gcc/7.3.0
    wget -O- https://cran.r-project.org/src/base/R-3/R-3.5.1.tar.gz | tar xz
    cd R-3.5.1/
    ./configure --prefix=${HOME}/local --enable-R-shlib --enable-shared --enable-memory-profiling --with-tcl-config=/usr/local/lib/tclConfig.sh --with-tk-config=/usr/local/lib/tkConfig.sh CC=${GCC_PREFIX}/bin/gcc CXX=${GCC_PREFIX}/bin/g++ CPPFLAGS="-I${GCC_PREFIX}/include -I/usr/local/include -I/usr/include" LDFLAGS="-L${GCC_PREFIX}/lib64 -L${GCC_PREFIX}/lib -L/usr/local/lib64 -L/usr/local/lib" F77=${GCC_PREFIX}/bin/gfortran FC=${GCC_PREFIX}/bin/gfortran JAVA_HOME=/usr/local/package/java/jdk1.8.0_162_64 PKG_CONFIG_PATH=/usr/local/lib64/pkgconfig:/usr/local/lib/pkgconfig:/usr/local/share/pkgconfig:/usr/lib64/pkgconfig:/usr/lib/pkgconfig:/usr/share/pkgconfig
    make -j4
    make install
    ```
    NEWS.pdfが無い云々と怒られたら適当に `touch doc/NEWS.pdf` などして凌ぐ。

1.  インストールしたRを起動し、パッケージをインストール:
    ```r
    install.packages("remotes")
    remotes::install_github("heavywatal/rtumopp")
    ```
    最新版にアップデートするのもこのコマンド。


### Linuxbrew

https://github.com/Linuxbrew/brew/wiki/CentOS6

まず `mpfr` インストール時の `make check` でコケるが、
`brew edit mpfr` で該当箇所をコメントアウトして対処。

しかし `glibc`/`gcc` のインストールがうまく行かない。
`LIBRARY_PATH` や `LD_LIBRARY_PATH` をいじっても
`libiconv.so.2` を見つけてもらえない:
```
HOMEBREW_NO_AUTO_UPDATE=1 HOMEBREW_BUILD_FROM_SOURCE=1 brew install gcc --without-glibc
==> Installing gcc --without-glibc
==> Downloading https://ftp.gnu.org/gnu/gcc/gcc-5.5.0/gcc-5.5.0.tar.xz
############################################################################################################################################################ 100.0$
==> ../configure --with-isl=/yshare1/home/watal/.linuxbrew/opt/isl@0.18 --with-bugurl=https://github.com/Linuxbrew/homebrew-core/issues --prefix=/yshare1/home/wat$
==> make
Last 15 lines from /home/watal/.cache/Homebrew/Logs/gcc/02.make:
/tmp/gcc-20180614-19441-uv2cp4/gcc-5.5.0/build/./gcc/xgcc -B/tmp/gcc-20180614-19441-uv2cp4/gcc-5.5.0/build/./gcc/ -dumpspecs > tmp-specs
/tmp/gcc-20180614-19441-uv2cp4/gcc-5.5.0/build/./gcc/xgcc: error while loading shared libraries: libiconv.so.2: cannot open shared object file: No such file or di$ectory
make[3]: *** [specs] Error 127
make[3]: *** Waiting for unfinished jobs....
/bin/sh ../../gcc/../move-if-change tmp-attrtab.c    insn-attrtab.c
/bin/sh ../../gcc/../move-if-change tmp-dfatab.c     insn-dfatab.c
/bin/sh ../../gcc/../move-if-change tmp-latencytab.c insn-latencytab.c
echo timestamp > s-attrtab
rm gcc.pod
make[3]: Leaving directory `/tmp/gcc-20180614-19441-uv2cp4/gcc-5.5.0/build/gcc'
make[2]: *** [all-stage1-gcc] Error 2
make[2]: Leaving directory `/tmp/gcc-20180614-19441-uv2cp4/gcc-5.5.0/build'
make[1]: *** [stage1-bubble] Error 2
make[1]: Leaving directory `/tmp/gcc-20180614-19441-uv2cp4/gcc-5.5.0/build'
make: *** [all] Error 2
```

シムリンクを張ることで一応通るようになる:
```sh
cd ~/.linuxbrew/lib/
ln -s /usr/local/lib/libiconv.so.2
cd ~/.linuxbrew/include/
ln -s /usr/local/include/gd.h
ln -s /usr/local/include/gdfontl.h
ln -s /usr/local/include/gdfonts.h
```

が、このglibcを使うと随所で不具合が...

