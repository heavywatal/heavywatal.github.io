+++
title = 'Rcpp'
subtitle = "RからC++を使う"
tags = ["r", "package"]
[menu.main]
  parent = "rstats"
  weight = -45
+++

- https://cran.r-project.org/package=rcpp
- http://www.rcpp.org/
- <span class="fragment" data-fragment-index="1">
  [みんなのRcpp](https://teuder.github.io/rcpp4everyone_ja/) and
  [Rcpp for everyone](https://teuder.github.io/rcpp4everyone_en/)
  by **津駄@teuder** さん
  </span>


### ベンチマークの一例

長さnの調和級数:

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

多くの場合、生のforループを避け、R標準のベクトル演算やちゃんとしたパッケージの関数を使うことで十分な速度が出る。
それよりもさらに速度や省メモリを追い求めたい場面で、Rcppが有用となる。


## 使い捨てコードで使う

ファイルあるいは文字列をコンパイルして使う:
```r
# Rcpp::sourceCpp("fibonacci.cpp")

Rcpp::sourceCpp(code='
  #include <Rcpp.h>
  // [[Rcpp::plugins(cpp14)]]
  // [[Rcpp::export]]
  int fibonacci(const int x) {
    if (x == 0) return 0;
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
    if (x == 0) return 0;
    if (x == 1) return 1;
    return fibonacci(x - 1) + fibonacci(x - 2);
  }
')
fibonacci(9L)
# [1] 34
```


## 自作パッケージで使う

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

-   [Rcpp](https://cran.r-project.org/web/packages/Rcpp)
    が型変換などをスムーズにしてくれる:
    ```c++
    // [[Rcpp::plugins(cpp14)]]
    #include <Rcpp.h>

    //' First example
    //' @param args string vector
    //' @export
    // [[Rcpp::export]]
    Rcpp::IntegerVector len(const std::vector<std::string>& args) {
        return {static_cast<int>(args.size())};
    }
    ```

その他の注意点

- `std::abort()` や `std::exit()` は呼び出したRセッションまで殺してしまう。
  例外は投げっぱなしで拾わなくても大丈夫で、
  `std::exception`の派生クラスなら`what()`まで表示してもらえる。
- グローバル変数やクラスのstaticメンバは `dyn.unload()` されるまで生き続ける。
  `parallel::mclapply()` とかでフォークした先での変更は子同士にも親にも影響しない。
