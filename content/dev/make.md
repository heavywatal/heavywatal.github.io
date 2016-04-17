+++
title = 'make'
[menu.main]
  parent = "dev"
+++

-   <http://www.gnu.org/software/make/>
-   <http://www.gnu.org/software/make/manual/make.html>

あらかじめコンパイルの命令を `Makefile` に書いておくことで、
複雑な `gcc` コマンドを何度も打つのを避けられる。

## Makefile

cppソースコードと同じディレクトリに入れるだけでいきなり使える
`Makefile` の例。

```makefile
####### Options

PROGRAM := ./a.out
CXX := g++
CC := ${CXX}
CXXFLAGS := -O3 -std=c++11
CPPFLAGS := -Wall -Wextra -fno-strict-aliasing -iquote/usr/local/include -iquote${HOME}/local/include ${DBGFLAGS}
LDFLAGS := -L/usr/local/lib -L${HOME}/local/lib
LDLIBS := -lboost_program_options
TARGET_ARCH := -m64 -msse -msse2 -msse3 -mfpmath=sse

####### Dependencies

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

-include Dependfile
Dependfile:
        ${CXX} -MM ${CPPFLAGS} ${CXXFLAGS} ${TARGET_ARCH} ${SRCS} > Dependfile

####### Targets

.DEFAULT_GOAL := all
.PHONY: all clean run debug

all: ${PROGRAM}
        @:

clean:
        ${RM} ${OBJS} ${PROGRAM}

run:
        ${PROGRAM}

debug:
        ${MAKE} all DBGFLAGS="-g -DDEBUG"

${PROGRAM}: ${OBJS}
        ${LINK.o} $^ ${LOADLIBES} ${LDLIBS} ${OUTPUT_OPTION}
```

### Automatic Variables

<http://www.gnu.org/software/make/manual/make.html#Automatic-Variables>

`$@`
:   ターゲット

`$^`
:   必須項目のスペース区切り

`$<`
:   必須項目の先頭

### Implicit Variables

<http://www.gnu.org/software/make/manual/make.html#Name-Index>

<http://www.gnu.org/software/make/manual/make.html#Implicit-Variables>

`CC`
:   Cコンパイラ `cc`

`CXX`
:   C++コンパイラ `g++`

`COMPILE.cpp`
:   `$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c`

`LINK.cpp`
:   `$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)`

`LINK.o`
:   `$(CC) $(LDFLAGS) $(TARGET_ARCH)`

`OUTPUT_OPTION`
:   `-o $@`

`RM`
:   `rm -f`

`CPPFLAGS`
:   プリプロセッサ用オプション。 e.g. `-Wall -Wextra -fno-strict-aliasing -DNDEBUG -iquote ${HOME}/include`

`CXXFLAGS`
:   C++コンパイラ用オプション。 e.g. `-O3 -std=c++11`

`LDFLAGS`
:   ライブラリパスを指定する。 e.g. `-L/usr/local/lib -L{HOME}/local/lib`

`LDLIBS`, `LOADLIBES`
:   リンクするライブラリを指定する。 e.g. `-lboost_program_options -lpthread`

`TARGET_ARCH`
:   マシン依存なオプションを指定する。 e.g. `-march=core2 -m64 -msse -msse2 -msse3 -mfpmath=sse`

### Functions

<http://www.gnu.org/software/make/manual/make.html#Functions>

文字列関連 <http://www.gnu.org/software/make/manual/make.html#Text-Functions>
:   `$(subst {from},{to},{text})`

    `$(findstring {find},{in})`

    `$(filter {pattern...},{text})`

ファイル名 <http://www.gnu.org/software/make/manual/make.html#File-Name-Functions>
:   `$(dir {names...})`

    `$(notdir {names...})`

    `$(basename {names...})`

    `$(addprefix {prefix},{names...})`

    `$(wildcard {pattern})`

    `$(abspath {names...})`

条件分岐 <http://www.gnu.org/software/make/manual/make.html#Conditional-Functions>
:   `$(if {condition},{then},{else})`

    関数じゃない `ifeq`, `ifneq`, `ifdef`, `ifndef`, `else`, `endif`
    もある。
    <http://www.gnu.org/software/make/manual/make.html#Conditional-Syntax>

その他
:   

    `$(foreach {var},{list},{text})`
    :   `{list}` の中身をそれぞれ `{var}` に入れて
        `{text}` を実行

    `$(file {op} {filename},{text})`
    :   `{text}` の結果をファイルに書き出す

    `$(call {variable},{params...})`
    :   `$(1) $(2)` などを使って定義しておいた
        `variable` を関数のように呼び出す

    `$(origin {variable})`
    :   変数がどう定義されたかを知れる

    `$(error {text...})`, `$(warning {text...})`, `$(info {text...})`
    :   エラーや警告をプリントする

    `$(shell {command...})`
    :   シェルを呼び出す

### Target

<http://www.gnu.org/software/make/manual/make.html#Standard-Targets>

`all`
:   ディレクトリ内のcppソースをコンパイル

`clean`
:   コンパイル済みオブジェクトを一掃

`run`
:   プログラムを走らせてみる

`debug`
:   デバッグモードでコンパイル

`Dependfile`
:   依存関係を読み取って `Dependfile` に書き出す

`' '`
:   v3.81以降であれば `.DEFAULT_GOAL` が効くので `make all` と同じ

<!-- -->

    % make clean
    % make
    % make run

### Option

<http://www.gnu.org/software/make/manual/make.html#Options-Summary>

`-f file`
:   `Makefile` じゃない名前のファイルを指定したければ

`-j jobs`
:   並列コンパイル。コア数+1くらいがちょうどいいらしい

`-C directory`
:   そのディレクトリに行って `make`

`-p`
:   自動的に作られるものも含めてすべての変数を表示
