+++
title = 'C++コマンドライン引数'
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

## 理想

-   ヘッダ1つの `#include` で使える
-   オプションの管理を `main()` でなくクラスのソースファイルなどに分散できる
-   オプション定義時に格納先の変数をアドレス渡しで指定できる
-   ファイルや `std::string` からも読み込める
-   ヘルプを自動生成してくれる
-   読み込める形式で全オプションと値を `std::string` あるいはファイルに書き出せる
-   マクロではなく `template` など真っ当な C++ で書かれている

## GNU `getopt`

<http://www.gnu.org/s/libc/manual/html_node/Getopt.html>

-   `gcc` ならインストール不要、`#include <getopt.h>` するだけで使える
-   オプションの管理を各ソースファイルに分散できない
-   ヘルプなど自動生成してくれない
-   使い方がわかりにくい

## `boost::program_options`

<http://www.boost.org/doc/html/program_options.html>

-   格納する変数のアドレスを渡せる
-   ファイルから読み込みをサポート
-   読み込み可能なファイルの出力方法は用意されてないので自分で書く必要がある
-   Boostライブラリのヘッダだけでなくビルドとリンクが必要 cf. [boost]({{< relref "boost.md" >}})
-

    各ソースファイルに分散するには
    :   1.  各クラスの `static` メソッドで `options_description` オブジェクトを生成
        2.  `main()` の `options_description` オブジェクトに `add()` する。

## `gflags`

<https://gflags.github.io/gflags/>

-   テンプレートではなくマクロをふんだんに使って実装されているのでちょっと怖い
-   `main()` 関数まわりをほとんど変更することなく、 各ソースファイルで自由にオプションを定義できる(しかもたった1行で)
-   格納されるのは頭に `FLAGS_` という接頭辞がついた変数
-   入力可能なファイルを出力することも可能
-   ライブラリのビルドとリンクが必要

### Usage

`main()` 関数に書く必要があるのはこれだけ

```c++
##include <gflags/gflags.h>

int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("This is a program to test gflags");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    // do something
    return 0;
}
```

あとは個々のソースファイルでオプションを追加。`namespace` にも入れられる。

```c++
##include <gflags/gflags.h>

namespace tapiola {
    DEFINE_uint64(sibelius, 0, "string that is displayed with --help flag");
}

void func(){
    std::cout << tapiola::FLAGS_sibelius << std::endl;
}
```

## cmdline

<https://github.com/tanakh/cmdline>\
<http://d.hatena.ne.jp/tanakh/20091028>

-   たった1つのヘッダファイルで構成されてるのでインストールもビルドも楽チン
-   直感的でC++らしいデザインなので分かりやすい
-   ファイルの読み取りはサポートしていないが `std::string` は読める
-   変数に直接格納することはできず、パーサのメソッドで値を取得:
    `template <class T> const T &parser::get(const std::string &name)`
-

    各ソースファイルに分散できるっちゃできる？
    :   1.  `main()` あたりで `parser` オブジェクトを定義
        2.  それを各クラスの `static` メソッドに渡し、中で `parser::add()`
        3.  `main()` で `parser::parse(argc, argv)`
        4.  再び各クラスの `static` メソッドを呼んで `parser::get()` から変数に代入

## TCLAP

<http://tclap.sourceforge.net/>

-   "Templatized C++ Command Line Parser Library"の名のとおり
    `template` で書かれておりヘッダだけで構成される
-   が、`configure` と `make install` というインストール手順を踏む
-   ファイルの読み取りはサポートしていないが、`std::string` からパース可能
-   読み込めるファイルの出力は用意されていないが、自分で書くのは簡単そう
-   格納する変数は指定できず、`*Arg` オブジェクトの `getValue()` メソッドで値を取得
-

    各ソースファイルに分散するには
    :   1.  `main()` あたりで `CmdLine` オブジェクトを定義
        2.  それを各クラスの `static` メソッドに渡し、中で `*Arg` オブジェクトを生成
        3.  `main()` で `CmdLine::parse(argc, argv)`
        4.  各 `*Args` オブジェクトに格納されている値を `getValue()` メソッドで取得

## getoptpp

<http://code.google.com/p/getoptpp/>

-   `istream` っぽく `>>operator` を使う
-   格納する変数を指定できる
-   基本的にはライブラリをビルドして使うが、ちょっといじればヘッダの `#include` だけでも使える
-   ファイルの読み込みやヘルプの生成は一切手伝ってくれない
-

    各ソースファイルに分散するには
    :   1.  `main()` あたりで `(argc, argv)` を引数に `GetOpt_pp` オブジェクトを生成
        2.  それを各クラスの `static` メソッドに渡し、中で値を取得:
            `args >> Option('i', "long_option", var)`
