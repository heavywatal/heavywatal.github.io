+++
title = '擬似乱数生成器'
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

-   `<cstdlib>` の `std::rand()` は乱数の質も悪く、速度も遅いので非推奨。
    C++11 から標準ライブラリに追加された `<random>` を使う。
-   `<algorithm>` の `std::random_shuffle()` は引数省略で
    `std::rand()` が使われてしまうので非推奨。
    C++11 で追加された `std::shuffle()` に生成器を明示的に渡して使う。
-   非標準の生成器としてはSFMTやdSFMTが高速で高品質。
    徐々にPCGやXorshift系の利用も広がってきている。


## `<random>`

-   <https://en.cppreference.com/w/cpp/numeric/random>
-   <https://cpprefjp.github.io/reference/random.html>
-   <https://cplusplus.com/reference/random/>

C++11 ではまともに使える乱数ライブラリが追加された。
乱数生成エンジンと分布関数オブジェクトを組み合わせて使う。

```c++
#include <iostream>
#include <random>

int main() {
    // seed
    std::random_device rd;
    const auto seed = rd();

    // engine
    std::mt19937 rng(seed);

    // probability density distribution
    const double mean = 0.0;
    const double sd = 1.0;
    std::normal_distribution<double> dist(mean, sd);

    // generate!
    for (size_t i=0; i<8; ++i) {
        std::cout << dist(rng) << std::endl;
    }

    return 0;
}
```

## Mersenne Twister

<http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/mt.html>

松本眞と西村拓士によって開発された高速・高品質な擬似乱数生成器。
標準の `<random>` でも利用可能になっており、
パラメータ定義済みの `std::mt19937` がよく使われる。
[Rのデフォルト](https://stat.ethz.ch/R-manual/R-patched/library/base/html/Random.html)にも採用されている。

### SFMT

<http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/>

Mersenne Twisterを松本眞と斎藤睦夫がさらに改良したもの。
SIMD命令を利用して、速度も品質も向上したらしい。
標準には含まれず、ソースからのビルドが必要。

double型を直接生成する亜種
[dSFMT](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/#dSFMT)
もある。
整数乱数から除算で変換するよりも良質で高速。
v2.1からは整数も出力可能。
かつてはJuliaのデフォルトに採用されていた。

本家の実装はC言語の関数として提供されていて、
C++標準 `<random>` と組み合わせて使う設計にはなっていない。
そこで `std::mt19937` と同じように使えるようにしたラッパークラス
[`wtl::sfmt19937`](https://github.com/heavywatal/sfmt-class/) を書いた。
[CMake]({{< relref "cmake.md" >}})
でSFMTを簡単に導入するためのインストーラをとしても使える。


## PCG

<https://www.pcg-random.org/>

線形合同法(LCG: linear congruential generator)の出力をpermutationする。
周期はMTに比べれば短いけど、省メモリ、高速で、品質テストの結果も優秀。
[`numpy.random.default_rng`](https://numpy.org/doc/stable/reference/random/generator.html)
としても採用されている。
作者O'Neillのブログは読み応えがあっておもしろい。

本家C++実装 [pcg-cpp](https://github.com/imneme/pcg-cpp) は
[CMake]({{< relref "cmake.md" >}}) 非対応なので使いにくい。
[issue](https://github.com/imneme/pcg-cpp/issues/43) にも挙がり
[PR](https://github.com/imneme/pcg-cpp/pull/44) も提出されたが、
巨大な差分を生じる悪手によりマージされず頓挫したらしい。
とりあえず自分のフォーク
[heavywatal/pcg-cpp](https://github.com/heavywatal/pcg-cpp)
で最低限の `CMakeLists.txt` を書いて対処。

また、本家pcg-cppは実験的な亜種クラスを大量に含んでいてソースが複雑化・肥大化している。
そこで `pcg32` と `pcg64` に必要な部分だけを抜き出し、
single-headerで手軽に使える軽量版
[pcglite](https://github.com/heavywatal/pcglite)
を作ってみた。

- [pcg-cpp](https://github.com/heavywatal/pcg-cpp):
  ~3600 lines / 3 headers
- [pcglite](https://github.com/heavywatal/pcglite):
  ~250 lines / 1 header


## Xorshift family

<https://prng.di.unimi.it/>

[Marsaglia (2003)](https://www.jstatsoft.org/article/view/v008i14)
から始まって、別の開発者(主にVigna)によっていくつかの改良版が派生している。
Xoshiro256++が[Juliaのデフォルト](https://docs.julialang.org/en/v1/stdlib/Random/)に採用されている。

PCG作者との論争を見る限り、こちらの作者は口が悪い。
乱数の品質でどちらが優れているのかは分からないけど。


## Seed

### `/dev/urandom`

`/dev/random` は擬似乱数を生成するデバイスで、
環境ノイズから乱数を生成するため自然乱数に近いが、遅い。
`/dev/urandom` は十分なノイズが蓄積していなくても
内部プールの再利用によってすぐに生成してくれるが、
そのために `/dev/random` ほど安全ではない。
長期に渡って使われる暗号鍵の生成以外の目的、
すなわちランダムっぽい乱数シードの生成のような目的では
`/dev/urandom` の利用で十分らしい。
が、今やこれを直接読むコードは必要なく、次のように標準ライブラリを利用する。

### `std::random_device`

C++11から標準の `<random>` で提供される。
基本的には `/dev/urandom` から生成するらしい。

```c++
std::random_device seeder;
const auto seed = seeder();
```

大量のシミュレーションを回すときなど、
エントロピーがもっと必要な場合は `std::seed_seq` を利用する。
