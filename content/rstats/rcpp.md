+++
title = 'Rcpp'
subtitle = "RからC++を使う"
tags = ["r", "c++", "package"]
[menu.main]
  parent = "rstats"
  weight = -45
+++

プログラムの書き方によって速度やメモリ効率は大きく変わる。
Rでは大抵、生のforループを避けて、R標準のベクトル演算やちゃんとしたパッケージの関数を使っていれば大丈夫。
でも、どうしても、さらに速度を追い求めたい場合にはRcppが有用となる。

長さnの調和級数を求める例:

```r
r_for = function(n) {
  s = 0; for (i in seq_len(n)) {s = s + 1 / i}; s
}

r_vec = function(n) sum(1 / seq_len(n))

Rcpp::cppFunction("double rcpp(int n) {
  double s = 0; for (int i = 1; i <= n; ++i) {s += 1.0 / i;} return s;
}")  # Compilation takes a few seconds here

n = 1000000L
rbenchmark::benchmark(r_for(n), r_vec(n), rcpp(n))[,1:4]

#       test replications elapsed relative
# 1 r_for(n)          100   3.968   29.835
# 2 r_vec(n)          100   0.473    3.556
# 3  rcpp(n)          100   0.133    1.000
```

## Documentation

- Project Home: http://www.rcpp.org/
- CRAN: https://cran.r-project.org/package=Rcpp
    - [Rcpp-JSS-2011.pdf](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-jss-2011.pdf):
      原典。
    - [Rcpp-introduction.pdf](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-introduction.pdf):
      なぜRcppを使うのか。
    - [Rcpp-attributes.pdf](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-attributes.pdf)
      現在主流となっているRcppの使い方全般。
      それに "Rcpp Attributes" という名前がついていて、
      "inline" という古いパッケージのやり方を置き換えたらしい。
    - [Rcpp-modules.pdf](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-modules.pdf):
      関数やclassをRにexportする。
    - [Rcpp-package.pdf](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-package.pdf):
      自作RパッケージでRcppを使う。
    - [Rcpp-extending.pdf](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-extending.pdf):
      自作classをRcppで扱う。
    - [Rcpp-sugar.pdf](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-sugar.pdf):
      ベクトル化とlazy評価が効くRの記法をC++側で使う。
    - [Rcpp-quickref.pdf](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-quickref.pdf)
    - [Rcpp-FAQ.pdf](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf)
- GitHub: https://github.com/RcppCore/Rcpp
    - 上のPDFは分量が多いわりに意外と網羅的ではない。
      ざっくり読んでなんとなく分かってきたら、
      さらなるドキュメントを求めてネットの海を彷徨うよりソースコードに当たったほうが早い。特に
      [`inst/unitTests`](https://github.com/RcppCore/Rcpp/tree/master/inst/unitTests)
      はかなり参考になる。
- API: http://dirk.eddelbuettel.com/code/rcpp/html/
- Advanced R: [Rewriting R code in C++](https://adv-r.hadley.nz/rcpp.html)
- <span class="fragment" data-fragment-index="1">
  [みんなのRcpp](https://teuder.github.io/rcpp4everyone_ja/) and
  [Rcpp for everyone](https://teuder.github.io/rcpp4everyone_en/)
  by 津駄@teuderさん
  </span>

## Rスクリプトの途中で使う

ファイルあるいは文字列をコンパイルして使う:
```r
Rcpp::sourceCpp("fibonacci.cpp")

Rcpp::sourceCpp(code='
  #include <Rcpp.h>
  // [[Rcpp::plugins(cpp14)]]
  // [[Rcpp::export]]
  int fibonacci(const int x) {
    if (x < 1) return 0;
    if (x == 1) return 1;
    return fibonacci(x - 1) + fibonacci(x - 2);
  }
')
fibonacci(9L)
# [1] 34
```

いろいろな準備を任せて、関数をひとつだけ定義するショートカット:
```r
Rcpp::cppFunction(plugins = c("cpp14"), '
  int fibonacci(const int x) {
    if (x < 1) return 0;
    if (x == 1) return 1;
    return fibonacci(x - 1) + fibonacci(x - 2);
  }
')
fibonacci(9L)
# [1] 34
```


## Rパッケージで使う

- [Rcpp-package.pdf](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-package.pdf) by Dirk Eddelbuettel and Romain François
- [Compiled code (`src/`) - R packages](http://r-pkgs.had.co.nz/src.html) by Hadley Wickham

### 準備手順

-   まずRcppコード以外の部分を作っておく。
    [devtoolsページを参照]({{< relref "devtools.md" >}})。

-   `usethis::use_rcpp()` を実行して設定を整える。
    `DESCRIPTION` や `src/.gitignore` などが書き換えられる。

-   `R/***-package.R` など適当なとこに `@useDynLib` の設定を追加:

    ```r
    #' @useDynLib hello, .registration = TRUE
    #' @importFrom Rcpp sourceCpp
    #' @keywords internal
    "_PACKAGE"
    ```
    `@importFrom Rcpp sourceCpp` を省くと、パッケージ利用時に
    `'enterRNGScope' not provided by package 'Rcpp'`
    のようなエラーが出る場合がある
    (明示的に `library(Rcpp)` するなどして既にRcppロード済みの環境では動く)。

-   同じところに `.onUnload` も定義しておく:
    ```r
    .onUnload = function(libpath) {
      library.dynam.unload("hello", libpath)
    }
    ```
    すると `unloadNamespace("hello")` したときに共有ライブラリもちゃんと外れるようになる。
    ちなみに `devtools::unload()` はこれを省略してもちゃんとリロードしてくれる。

-   外部ライブラリのリンクに関する設定など、
    開発者側で指定すべきビルドオプションは `src/Makevars` に指定:
    ```
    CXX_STD=CXX14
    PKG_CPPFLAGS=-DSTRICT_R_HEADERS -I/usr/local/include
    PKG_LIBS=-L/usr/local/lib -Wl,-rpath,/usr/local/lib -lthankyou
    ```
    `STRICT_R_HEADERS` を定義しておくことで余計なマクロ定義を防げる。
    `configure` や [CMake]({{< relref "cmake.md" >}}) を使って
    `src/Makevars.in` から生成する手もある。

    参考: [Japan.R 2018 LT "Rcppパッケージで外部C++ライブラリを使う"](https://heavywatal.github.io/slides/japanr2018/)

-   どうしてもユーザ側で指定すべきオプションがある場合は
    `~/.R/Makevars` に書いてもらう。
    例えばMPI依存パッケージをmacOSでビルドでしようとすると
    `clang: error: unsupported option '-fopenmp'`
    と怒られるので `brew install llvm`
    で別のコンパイラを入れて下記のように指定する:
    ```
    LLVM_LOC=/usr/local/opt/llvm
    CC=$(LLVM_LOC)/bin/clang
    CXX=$(LLVM_LOC)/bin/clang++
    ```

-   `src/` 以下にソースコードを書く。


### ソースコード `src/*.cpp`

```c++
// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>

//' First example
//' @param args string vector
//' @export
// [[Rcpp::export]]
int len(const std::vector<std::string>& args) {
    return args.size();
}
```

- `std::abort()` や `std::exit()` は呼び出したRセッションまで殺してしまう。
  例外は投げっぱなしで拾わなくても大丈夫で、
  `std::exception`の派生クラスなら`what()`まで表示してもらえる。
- グローバル変数やクラスのstaticメンバは `dyn.unload()` されるまで生き続ける。
  `parallel::mclapply()` とかでフォークした先での変更は子同士にも親にも影響しない。


## 詳細

アタリがついてる場合は
[namespace Rcpp](http://dirk.eddelbuettel.com/code/rcpp/html/namespaceRcpp.html)
とかからブラウザのページ内検索で探すのが早い。

Rcppで楽ができるとはいえ、R本体の内部情報もいずれ知ることになる。。。

- https://github.com/hadley/r-internals
- https://cran.r-project.org/doc/manuals/r-release/R-ints.html
- https://github.com/wch/r-source/blob/trunk/src/include/Rinternals.h


### 型

`SEXP`: S Expression
:   Rのあらゆるオブジェクトを表すC言語上の型。
    これの扱いは Rcpp が肩代わりしてくれるので基本的には直接触らない。

[`Rcpp::RObject`](http://dirk.eddelbuettel.com/code/rcpp/html/classRcpp_1_1RObjectMethods.html)
:   `SEXP` の thin wrapper であり Rcpp から R の変数を扱う際の基本クラス。
    明示的に `PROTECT`/`UNPROTECT` を書かずにRAIIで済ませられる。


[`Rcpp::Vector<T>`](http://dirk.eddelbuettel.com/code/rcpp/html/classRcpp_1_1Vector.html)
:   [`vector/instantiation.h`](http://dirk.eddelbuettel.com/code/rcpp/html/instantiation_8h_source.html)
    抜粋:

    ```c++
    typedef Vector<LGLSXP>  LogicalVector;
    typedef Vector<INTSXP>  IntegerVector;
    typedef Vector<REALSXP> NumericVector; // DoubleVector
    typedef Vector<STRSXP>  StringVector;  // CharacterVector
    typedef Vector<VECSXP>  List;          // GenericVector
    ```

    クラスのメンバとして生の配列ではなくそこへの参照を保持する。
    しかし `std::vector` とは異なり、
    このオブジェクトをコピーしてもメモリ上の中身はコピーされず、
    ふたつとも同じ生配列を参照する。

    C++関数がRから呼ばれるとき
    `Rcpp::Vector<>` 受け取りの場合はうまく参照渡しになるが、
    `const std::vector<>&` 受け取りの場合はコピーが発生する。


[`Rcpp::DataFrame`](http://dirk.eddelbuettel.com/code/rcpp/html/classRcpp_1_1DataFrame__Impl.html)
:   Rの上では強力だけどC++内では扱いにくい。
    出力として使うだけに留めるのが無難。


関数オーバーロードもテンプレートもそのままRにexportすることはできない。
実行時の型情報で振り分ける関数で包んでexportする必要がある。
http://gallery.rcpp.org/articles/rcpp-return-macros/

### タグ

`[[Rcpp::export]]`
:   これがついてる関数は `RcppExport.cpp` を介してライブラリに登録され、
    <code>.Call(`_{PACKAGE}_{FUNCTION}`)</code>
    のような形でRから呼び出せる様になる。
    それを元の名前で行えるような関数も `RcppExport.R` に自動で定義してもらえる。
:   `[[Rcpp::export(".new_name_here")]]`
    のように名前を変更することもできる。
    ドットで始まる名前にしておけば
    `load_all(export_all=TRUE)`
    の状態での名前空間汚染を多少調整できる。
:   Rパッケージの `NAMESPACE` における `export()` とは別物。

`[[Rcpp::plugins(cpp14)]]`
:   `Makevars` に `CXX_STD=CXX14` を書くのとどう違う？
    ほかにどんなプラグインが利用可能？

`[[Rcpp::depends(RcppArmadillo)]]`
:   ほかのパッケージへの依存性を宣言。
    たぶんビルド時のオプションをうまくやってくれる。
    `#include` は自分で。

`[[Rcpp::interfaces(r,cpp)]]`
:   `Rcpp::export` するとき、どの言語向けにいろいろ生成するか。
    何も指定しなければ `r` のみ。
    `cpp` を指定すると、ほかのパッケージから
    `Rcpp::depends` できるようにヘッダーを用意してくれたりするらしい。

`[[Rcpp::init]]`

`[[Rcpp::internal]]`

`[[Rcpp::register]]`


## 自作C/C++クラスをRで使えるようにする

"Rcpp Modules" の機能を使う。

1.  `RcppExports.cpp` に自動的に読み込んでもらえるヘッダー
    (e.g., `src/{packagename}_types.h`)
    で自作クラスの宣言と
    `Rcpp::as<MyClass>()` / `Rcpp::wrap<MyClass>()` の特殊化を行う。

    ```c++
    #include <RcppCommon.h>

    RCPP_EXPOSED_CLASS(MyClass);
    // これで as<MyClass> / wrap<MyClass> の特殊化が定義される
    // 必ず #include <Rcpp.h> より前に来るように

    #include "myclass.hpp"
    // 自作クラスの宣言
    ```


1.  どこかのソースファイルでモジュールを定義

    ```c++
    #include <Rcpp.h>`

    RCPP_MODULE(mymodule) {
      Rcpp::class_<MyClass>("MyClass")
        .constructor<int>()
        .const_method("get_x", &MyClass::get_x)
      ;
    }
    ```

1.  `{packagename}-package.R` でモジュールを読み込む。
    関数やクラスを全てそのまま公開するか、
    `Module` オブジェクト越しにアクセスさせるようにするか。

     ```r
     Rcpp::loadModule("mymodule", TRUE)`
     # obj = MyClass$new(42L)

     modulename = Rcpp::Module("mymodule")
     # obj = mymodule$MyClass$new(42L)
     ```

パッケージを読み込むといくつかのRC/S4クラスが定義される。

`Rcpp_MyClass`
:   `C++Object` を継承した Reference Class (RC)。

`C++Object`
:   R上でC++オブジェクトを扱うための親S4クラス。
    Rコンソール上での表示はこれの `show()` メソッドがデフォルトで利用される。

`C++Class`
:   コンストラクタをR側にexposeするためのクラスで、
    `MyClass$new(...)` のようにして新規オブジェクトを生成する。
    ただしデフォルト引数を扱えないのでファクトリ関数を普通に
    `[[Rcpp::Export]]` したほうが簡単かも。
:   staticメソッドも同様に扱えれば一貫性があったんだけど今のところ無理そう。
    `C++Function` としてならexposeできる。

`C++Function`
:   わざわざModule機能でexposeした関数を扱うS4。
    普通に `[[Rcpp::Export]]` する場合と比べたメリットは？

`Module`
:   `environment` を継承したS4。

RC/S4関連文献

- `?setRefClass` or https://stat.ethz.ch/R-manual/R-devel/library/methods/html/refClass.html
- https://adv-r.hadley.nz/s4.html
- http://adv-r.had.co.nz/OO-essentials.html#rc


### マクロ

http://dirk.eddelbuettel.com/code/rcpp/html/module_8h.html

`RCPP_EXPOSED_AS(MyClass)`
: `as<MyClass>` を定義してくれるマクロ。参照型やポインタ型もやってくれる。

`RCPP_EXPOSED_WRAP(MyClass)`
: `wrap<MyClass>` を定義してくれるマクロ。

`RCPP_EXPOSED_CLASS_NODECL(MyClass)`
: 上の2つを同時にやってくれるショートカット。

`RCPP_EXPOSED_CLASS(MyClass)`
: それらの前にさらに `class MyClass;` の前方宣言もする。




## 関連書籍

<a href="https://www.amazon.co.jp/dp/1461468671/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=ba6c791d1fab3179e6a351c2347bbdc9&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1461468671&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=1461468671" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/B0748CFJL3/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=2a5507d53af7170cbd1a42f540f72351&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=B0748CFJL3&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=B0748CFJL3" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
