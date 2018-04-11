+++
title = "ビット演算"
subtitle = "Bitwise operation"
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

https://github.com/heavywatal/scribble/blob/master/cxx/bitwise.cpp

## 基本

### 型

きっちりサイズが定義されているものと、
下限だけが定義されてて環境依存のものがある。

https://github.com/heavywatal/scribble/blob/master/cxx/sizeof.cpp

- 1byte == 8bits
- `sizeof(char)` == 1
- `sizeof(bool)` >= 1
- `sizeof(int)` >= 4
- `sizeof(uint32_t)` == 4
  [`<cstdint>`](http://en.cppreference.com/w/cpp/header/cstdint)
- `sizeof(uint128_t)` == 16
  [`<boost/multiprecision/cpp_int.hpp>`](https://boostjp.github.io/tips/multiprec-int.html)

符号あり(signed)型の場合、左端のビットが符号を司る。

```c++
// uint8_t
0b00000000 //    0
0b10000000 //  128
0b11111111 //  255

// int8_t
0b01111111 //  127
0b10000000 // -128
0b10000001 // -127
0b10000010 // -126
//
0b11111101 //   -3
0b11111110 //   -2
0b11111111 //   -1
```

`0b` は
[binary literal](http://en.cppreference.com/w/cpp/language/integer_literal)
(C++14) の接頭辞

### 演算子 operator

```c++
// AND, OR, XOR
x & y
x | y
x ^ y

// NOT, complement of 1
~ x

// negation, complement of 2, (~x | 1)
- x

// shift
x << n
x >> n
```

```c++
// assuming 8-bit
x = 5   // 0b00000101
x & 1   // 0b00000001 =  1
x | 2   // 0b00000111 =  7
x ^ 3   // 0b00000110 =  6
~ x     // 0b11111010 = -6 or 250u
- x     // 0b11111011 = -5 or 251u
x << 1  // 0b00001010 = 10
x >> 1  // 0b00000010 =  2
```

- `<<` は右から0を詰めて左端を捨てる
- `unsigned` に対する `>>` は左から0を詰めて右端を捨てる
- `signed` の負数に対する `>>` は未定義だが、だいたい左端をコピーして右端を捨てる(i.e., 符号を維持しつつ2で割る)
- 8-bit型に対する操作が必ずしも8-bitで返ってくるとは限らないし、
  暗黙の整数型は大概32-bitとかなので、両辺とも明示的に型を指定するほうが安心。
  特に上位ビットも変化する `~` や `-` を使って比較するときは要注意。


## 応用

符号反転 `-x` は、全ビット反転して1加えることに相当する (~x | 1) 。
つまり、符号反転して全ビット反転すると、1加えるのと同じ。

```c++
x = x + 1
x += 1
++x
x = -~x
```

### `std::vector<bool>`

- 特殊化されているため普通のSTL vectorではない
- 省メモリ
- イテレータあり e.g., `std::next_permutation()`
- `data()` メンバが無い。つまりナマのビット列としてアクセスできない？


### `std::bitset<N>`

- http://en.cppreference.com/w/cpp/utility/bitset
- http://www.cplusplus.com/reference/bitset/bitset/
- https://cpprefjp.github.io/reference/bitset.html

ビット数がコンパイル時定数


### `boost::dynamic_bitset<>`

http://www.boost.org/doc/libs/release/libs/dynamic_bitset/dynamic_bitset.html

可変サイズ版bitset

`find_first()`, `find_next()`, `is_subset_of()`
など便利な補助関数もあるが、速度的な最適化はされてないっぽいので、
実行速度がシビアな場面では普通に `operator[]` で自前ループを書いたほうがいいかも。
`operator&` とかも意外と遅い。
