+++
title = 'めざせC++高速プログラム'
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

## はじめに

-   速いプログラムで得られるメリットを超えるようなコストを払わないように。
    まずは動くプログラムを書いて目的を達成することが大事。
    自分律速じゃなくてプログラム律速だなと感じた段階でリファクタリングを考える。
-   プログラム本来の意図が読み取れなくなりそうなマニアックな高速化は避ける。
    清く正しくメンテナンスしやすいプログラムを書くほうが結局は生産的。
-   学習目的でない限り、車輪の再発明を避ける。
    やろうとしていることはきっと既に誰かが実現し、
    再利用可能な形で公開してくれているはず。
    まずは標準ライブラリとBoostを探してみる。

## 頑張れコンパイラ

Intelの icc でビルドされたプログラムは速いらしい。
gcc や clang の最適化技術も着々と進歩しており、
新しいコンパイラを使うほうがその恩恵を受けられる、はず。

### 最適化オプション

<http://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html>

```sh
% g++ -O3 main.cpp
```

### コンパイル時定数

定数同士の四則演算などがコード中に転がっていると、
コンパイルの段階でコンパイラが計算して１つの定数に置き換えてくれる。
つまり、実行時の計算を減らすことができる。
したがって、定数と変数が混ざった計算があるときは、定数をまとめておくとよさそう。

四則演算の中でも、割り算は特に遅いらしい。
コンパイル時に逆数のかけ算になるような書き方のほうが高速になる可能性がある。
(自動的にやってくれたりしないのかな?)

```c++
// 毎回割り算
rand() / 4294967296.0;

// カッコ内はコンパイル時に計算されて下のようになるはず
(1.0 / 4294967296.0) * rand();
(2.328306436538696e-10) * rand();
```

C++11では明示的に `constexpr` を指定することで、
確実にコンパイル時定数にできるので積極的に使うべし。
それ以外にも、テンプレートメタプログラミングで
なるべくコンパイル時に計算させることができる。

### インライン展開

関数の呼び出しにはコストがかかり、
コードが短くて何度も呼ばれる関数では特にそのコストの割合が高くて馬鹿にできない。
呼び出される予定のところにコンパイルの段階で関数を展開してやること
(=インライン展開) でそのコストを無くせる。

ファイルをまたいで使う関数の場合は定義宣言を
`ソース.cpp` ではなく `ヘッダ.h` に置き、頭に `inline` を付ける。
メンバ関数の場合はクラス定義の中で定義宣言にするだけでよい。
ただしそのように書いてもコンパイラへのヒントになるだけで、
インライン展開に必要でも十分でもない。

`-O2` では `-finline-small-functions` がオンになり、
`-O3` では `-finline-functions` がオンになる。
実行ファイルのサイズがデカくなるし、コンパイルは遅くなるので、バランスを考える。

```c++
inline double rexp(const double lambda=1.0) {
    return -log(1.0 - runif()) / lambda;
}
```

### 関数オブジェクト

関数オブジェクトとは、`operator()` を定義した `class` または `struct` のこと。
`std::for_each()` や `std::transform()` の最後の引数に関数を渡す場合など、
何度も呼び出される小さい関数については、
普通の関数や関数ポインタで渡すよりも関数オブジェクトを使ったほうが高速。
インライン展開などコンパイラによる最適化がしやすいらしい。
あと、メンバ変数として値を保持できるので、
引数のやり取りやメモリの確保・開放が少なくてすむという場面もありそう。

```c++
class isAboveThreshold {
  public:
    isAboveThreshold(const double t): threshold_(t);
    bool operator()(const double x) const {
        return x > threshold_;
    }
  private:
    const double threshold_;
};

std::transform(v.begin(), v.end(), result.begin(), isAboveThreshold(3.14));
```

C++11からはラムダ式が便利

```c++
const double threshold = 3.14;
std::transform(v.begin(), v.end(), result.begin(),
               [threshold](int x) -> bool {return x > threshold});
```

## 余計な一時オブジェクトを作らない

普通にプログラムを書くと、思わぬところで余計な一時オブジェクトが作られることになる。
メモリもcpu時間ももったいないので、なるべく避けよう。

### `const` 参照渡し or ポインタ渡し

`int` とか `double` のような単純な型の場合は気にしなくてもいいが、
そうでないオブジェクトの値渡しにはコピーコンストラクタが呼ばれるコストが生じる。
したがって、`std::string` とかSTLコンテナとか自作クラスなどは参照渡し
（そいつのアドレスだけを渡す）のほうが遥かに高速。
参照渡しされた仮引数の中身を関数内で変更すると、
呼び出し元の実引数の値も変更されることになる。
変更するつもりがない場合は明示的に `const` をつけて受け取るようにすると安全

```c++
void do_something(const std::vector<int>& huge_array) {
    // コピーせず本体への参照のみ受け取る。
    // constがついているので、変更しようとするとコンパイルエラー
}
```

引数の中身を変更するつもりがある場合はポインタを受け取るようにすると、
呼び出すときに `&` が必要となるので意図が伝わりやすい

```c++
void do_something(std::vector<int>* huge_array) {
    // コピーせずポインタを受け取る。
    // huge_arrayの変更は呼び出し元にも影響する。
}

do_something(&some_array);  // ああ、この配列は変更されるんだな
```

### 複合代入演算子

`int` や `double` とかならまだしも、`std::string` とかだとかなり違うはず。

```c++
// (a+b)の結果を持つ一時オブジェクトが作られ、xに代入される。
T x = a + b;

// xをaで初期化して、bを足す。一時オブジェクトは作られない。
T x(a);
x += b;
```

### ムーブ (C++11)

## コンテナ

### 用途に合わせる

`std::vector`
:   メモリ上に連続した領域を確保する。
    そのため、途中に `insert()` するのは遅い
    （メモリ領域を別のところに再確保して全部コピーしなければならないので）。
    その代わり、要素へのアクセスは高速で、インデックス `[ ]` によるランダムアクセスも可能。

`std::deque`
:   `vector` とほぼ同じだが、`reserve()` ができない。
    代わりに、`push_front()` と `pop_front()` が使える。
    つまり、先頭に対する追加・削除が必要な場合だけ使うコンテナ。

`std::list`
:   上記２つとは違って、飛び飛びの領域にまたがって存在できる。
    これは、各要素の値とともに前後の要素を示すイテレータを格納することで実現されている。
    そのため、途中への `insert()` はイテレータを書き換えるだけなので高速。
    ただし、その分メモリは余分に食うし、ランダムアクセスできないし、イテレータで総なめするのも遅い。

### メモリは一気に確保

`std::vector` の `push_back()` は勝手にメモリ領域を確保してくれるので、
大きさの心配をする必要がなくて便利。
これは、領域が足りなくなる都度「別のところに倍の領域を確保してコピー」
という処理をすることによって行われる。
始め1、次2、4、8、16、、、という具合に。
なので、10000個の要素を格納すると分かっている場合には、始めからそれだけ確保しておくと速い。

```c++
std::vector<int> v;
v.reserve(10000);
```

## ループ

`for` や `while` を回すときの継続条件式は回る度に評価されるので、
そのコストをよく考えなければならない。

### 継続条件

イテレータで回すときの `v.end()` も毎回呼び出されるので、
ループの中身が軽い処理の場合には無視できない差になるかもしれない。

```c++
for (size_t i=0; i<v.size(); ++i) {
    // このv.size()は毎回呼び出されてしまう。
}

for (size_t i=0, n=v.size(); i<n; ++i) {
    // size()の呼び出しは一度だけで済み、nの名前はforの中しか汚さない。
}

for (vector<int>::iterator it=v.begin(), v_end=v.end(); it!=v_end; ++it) {
    // これも、end()の呼び出しは一度きり。
}

for (const auto& x: v) {
    // コンテナ全体を舐めるなら C++11 range-based for
}
```

### 前置インクリメント

`for` でカウンタをインクリメントするだけなら
後置インクリメント `i++` でも
前置インクリメント `++i` でも計算結果は同じになる。
が、`i++` だと前の値を記憶してからプラスする（一時オブジェクトを作るという無駄が発生する）ので、
`++i` のほうがいいらしい。特にイテレータのとき。

```c++
for (std::vector<int>::iterator it=v.begin(); it!=v.end(); ++it) {
    // do something
}
```

### 出来る限り外で処理

`if-else` や `try-catch` のブロックをループの外で大きく取れないか、
一時変数の定義や演算などをループの外で予めやっておけないか、確認すべし。

## 入出力

### まとめて書き出す

ディスクレスクラスタで計算させる場合などは
特に通信やディスク書き込みのオーバーヘッドが大きい。
ファイルに逐次追記していくのではなく、
`std::ostringstream` や `std::vector<std::string>`
などでメモリ上にためておき、最後の最後でまとめて書き出すほうがいい。

### 標準入出力

Cストリーム(`std::printf` とか)とC++ストリーム(`std::cout` とか)
が混在するプログラムでもちゃんと関数が呼ばれた順に入出力を行うため、
デフォルトではこれら２つのストリームが同期するようになっている。
必要がなければ切っておく。

`std::cin` はデフォルトで `std::cout` に結びつけられてて、
`std::cin` される度に `std::flush` されてしまうらしいので、
そうならないように切り離す。

```c++
##include <iostream>

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
}
```

## 参考図書

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4774157155/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51Pwt65tXnL._SX150_.jpg" alt="C++ ポケットリファレンス" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4894714515/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51TFT3FMS1L._SX160_.jpg" alt="Effective C++ 原著第3版 (ADDISON-WESLEY PROFESSIONAL COMPUTING SERIES)" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4894714108/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41W5R878R7L._SX160_.jpg" alt="Effective STL―STLを効果的に使いこなす50の鉄則" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4894712458/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51BC2Q8XFRL._SX160_.jpg" alt="Efficient C++パフォーマンスプログラミングテクニック" /></a>
