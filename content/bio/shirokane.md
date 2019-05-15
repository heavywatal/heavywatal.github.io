+++
date = 2018-06-14T10:00:00+09:00
title = 'SHIROKANE'
subtitle = 'HGCスパコン'
tags = ["job"]
[menu.main]
  parent = "bio"
+++

https://supcom.hgc.jp/

## 環境整備

- [ハードウェア](https://supcom.hgc.jp/japanese/sys_const/system-main.html)
- [ソフトウェア](https://supcom.hgc.jp/internal/cgi/version_up_s3/select.cgi)
    - RHEL 7.6
    - `/usr/bin/cmake` 2.8.12
    - `/usr/bin/cmake3` 3.13.4
    - `/usr/bin/gcc` 4.8.5
    - `/usr/bin/python` 2.7.15
    - `/usr/bin/python3` 3.4.9
    - `/usr/local/package/boost/1.67.0/`
    - `/usr/local/package/gcc/7.3.0/`
    - `/usr/local/package/python/3.6.5/`
    - `/usr/local/package/r/3.5.0/`


## tumopp

### on R

http://heavywatal.github.io/rtumopp/

Shrokane5 で提供されている R-3.5 は
`/usr/local/package/gcc/7.3.0`でビルドされているため
C++14のライブラリも問題なく使える。

1.  普通に `cmake` とするとバージョン2.8.12のほうが参照されてしまうので、
    優先的にパスの通ってるところに
    `ln -s /usr/bin/cmake3 ~/local/bin/cmake`
    などとしてバージョン3.8以上が使えるようにする。

1.  Rを起動し、パッケージをインストール:
    ```r
    install.packages("devtools")
    devtools::install_github("heavywatal/rtumopp")
    ```
    最新版にアップデートするのもこのコマンド。

1.  パッケージを読み込んで実行:
    ```r
    library(tumopp)
    result = tumopp("-D3 -Chex -N256")
    ```


### command-line

http://heavywatal.github.io/tumopp/

Shrokane5 (RHEL 7) では問題なくLinuxbrewを使える。

```sh
brew install heavywatal/tap/tumopp
exec $SHELL -l
tumopp -D3 -Chex -N256
tumopp -h
```
