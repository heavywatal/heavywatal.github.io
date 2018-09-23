+++
title = 'C++コマンドライン引数'
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

## 理想

-   ビルド不要でヘッダ1つ
-   標準ライブラリのみに依存していてポータブル
-   ヘルプを自動生成してくれる
-   オプション定義時に格納先の変数を紐付けできる
-   `(argc, argv)` だけでなく `std::string` とかからも読み込める
-   読み込める形式で全ての値を書き出せる
-   マクロではなくtemplateやlambdaなど真っ当なC++
    (できればC++11以降の簡潔なスタイル) で書ける

## GNU `getopt`

<http://www.gnu.org/s/libc/manual/html_node/Getopt.html>

-   UNIX的な環境ならインストール不要だがC/C++標準ではない
-   ヘルプなど自動生成してくれない
-   C++というよりCなので手作業が多い

## `boost::program_options`

<http://www.boost.org/doc/html/program_options.html>

-   格納先の変数を紐付け可能
-   ファイルからも読み込める
-   読み込み可能なファイルの出力方法は用意されてないので自分で書く必要がある
-   ビルドとリンクが必要で大掛かり cf. [boost]({{< relref "boost.md" >}})

## `gflags`

<https://gflags.github.io/gflags/>

-   テンプレートではなくマクロをふんだんに使って実装されているのでちょっと怖い
-   `main()` 関数まわりをほとんど変更することなく、 各ソースファイルで自由にオプションを定義できる(しかもたった1行で)
-   接頭辞 `FLAGS_` のついた変数が自動的に定義されて、そこに値が格納される
-   入力可能なファイルを出力することも可能
-   要ビルド＆リンク

### Usage

`main()` 関数に書く必要があるのはこれだけ

```c++
#include <gflags/gflags.h>

int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("This is a program to test gflags");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    // do something
    return 0;
}
```

あとは個々のソースファイルでオプションを追加。`namespace` にも入れられる。

```c++
#include <gflags/gflags.h>

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

-   ヘッダファイル1つ
-   直感的でC++らしいデザインなので分かりやすい
-   demangle機能のためにポータビリティが犠牲に
-   `std::string` から読める
-   変数に直接格納することはできず、パーサのメソッドで値を取得:
    `template <class T> const T &parser::get(const std::string &name)`

## TCLAP

<http://tclap.sourceforge.net/>

-   "Templatized C++ Command Line Parser Library"の名のとおり
    `template` で書かれておりヘッダだけで構成される
-   が、`configure` と `make install` というインストール手順を踏む
-   `std::string` からパース可能
-   格納する変数は指定できず、`*Arg` オブジェクトの `getValue()` メソッドで値を取得

## getoptpp

<http://code.google.com/p/getoptpp/>

-   `std::istream` っぽく `>>operator` を使う
-   格納する変数を指定できる
-   基本的にはライブラリをビルドして使うが、ちょっといじればヘッダの `#include` だけでも使える
-   ファイルの読み込みやヘルプの生成は一切手伝ってくれない

## clipp

<https://github.com/muellan/clipp>

:   理念がしっかりしていて、かなり柔軟に使える。
:   そのまま使うには少し難しかったり、
    値を一括して取得する機能が欠けたりという問題はある。
:   [nlohmann/json](https://github.com/nlohmann/json)
    を使ってそのへんをうまくやる補助ライブラリ
    [clippson](https://github.com/heavywatal/clippson)
    を作って利用中。


## そのほか

Githubで上位に出てくるこれらもそのうち試したい:

https://github.com/docopt/docopt.cpp
:   ヘルプを自動生成するのではなく、ヘルプからパーサを構築する
:   元々はPython用に作られ、それから多言語に移植されてる実績
:   `bool`か`std::string`でゲットするしかないので、手動でキャストして代入
:   `std::vector<std::string>`から読める
:   要ビルド＆リンク (header-only化しようとしてる雰囲気はある)

https://github.com/jarro2783/cxxopts
:   Lightweight C++ command line option parser
:   `parse(argc, argv)` だけ

https://github.com/Taywee/args
:   名前空間が `args` という大胆さ
:   同じ名前を何回も書かなきゃいけないような、少々やぼったいインターフェイス
:   ヘッダ1つ、ヘルプ自動生成なのは良い

https://github.com/adishavit/argh
:   A minimalist argument handler.
:   ハイフンの数を区別できないし、ヘルプ自動生成も無い。
