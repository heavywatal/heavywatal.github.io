+++
title = 'make'
tags = ["c++"]
[menu.main]
  parent = "dev"
+++

- https://www.gnu.org/software/make/
- https://www.gnu.org/software/make/manual/make.html

あらかじめコンパイルの命令を `Makefile` に書いておくことで、
複雑なコマンドを何度も打つのを避けられる。

## Makefile

C++ソースコードと同じディレクトリに入れるだけでいきなり使える
`Makefile` の例。

```makefile
## Options

PROGRAM := a.out
CXX := clang++
CC := clang
CXXFLAGS := -O3 -std=c++14
CPPFLAGS := -Wall -Wextra -I/usr/local/include -I${HOME}/local/include
LDFLAGS := -L/usr/local/lib -L${HOME}/local/lib
#LDLIBS := -lboost_program_options
TARGET_ARCH := -m64 -msse -msse2 -msse3 -mfpmath=sse

## Dependencies

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

-include Dependfile
Dependfile:
        ${CXX} -MM ${CPPFLAGS} ${CXXFLAGS} ${TARGET_ARCH} ${SRCS} > Dependfile

## Targets

.DEFAULT_GOAL := all
.PHONY: all clean

all: ${PROGRAM}
        @:

${PROGRAM}: ${OBJS}
        ${LINK.cpp} ${OUTPUT_OPTION} $^ ${LDLIBS}

clean:
        ${RM} ${OBJS} ${PROGRAM}
```


### [Rule](https://www.gnu.org/software/make/manual/make.html#Rules)

コロンとタブを使って以下のような形式でルールを書くのが基本。
この`Makefile`があるところでターミナルから `make TARGET` と打つと、
ターゲットよりもソースファイルが新しい場合にコマンドが実行される。

```make
TARGET : SOURCE1 SOURCE2
        COMMAND

a.out : main.cpp sub.cpp
        g++ -O2 main.cpp sub.cpp
```

下記のようなパターンルールが予め定義されている。
(see [Pattern Rule](https://www.gnu.org/software/make/manual/make.html#Pattern-Rules))

```make
%.o : %.c
        $(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
%.o : %.cpp
        $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
```

以下に紹介するように、ほかにも様々な変数や関数が用意されていて、
個別のファイル名などをいちいち入力しなくても済むようになっている。


### [Implicit Variables](https://www.gnu.org/software/make/manual/make.html#Implicit-Variables)

<https://www.gnu.org/software/make/manual/make.html#Name-Index>

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
:   プリプロセッサ用オプション。
    e.g., `-Wall -Wextra -DNDEBUG -I${HOME}/local/include`

`CXXFLAGS`
:   C++コンパイラ用オプション。 e.g., `-O3 -std=c++14`

`LDFLAGS`
:   ライブラリパスを指定する。 e.g., `-L/usr/local/lib -L{HOME}/local/lib`

`LDLIBS`
:   リンクするライブラリを指定する。
    昔は`LOADLIBES`も同じ機能だったが非推奨。
    e.g., `-lboost_program_options -lz`

`TARGET_ARCH`
:   マシン依存なオプションを指定する。
    e.g., `-march=native -m64 -msse -msse2 -msse3 -mfpmath=sse`


### [Automatic Variables](https://www.gnu.org/software/make/manual/make.html#Automatic-Variables)

`$@`
:   ターゲット

`$<`
:   必須項目の先頭

`$^`
:   必須項目のスペース区切り
:   重複してても削らずそのまま欲しい場合は `$+`
:   新しく更新があったファイルだけ欲しい場合は `$?`


### [Functions](https://www.gnu.org/software/make/manual/make.html#Functions)

[文字列関連](https://www.gnu.org/software/make/manual/make.html#Text-Functions)
:   `$(subst FROM,TO,TEXT)`
:   `$(findstring FIND,IN)`
:   `$(filter PATTERN...,TEXT)`

[ファイル名](https://www.gnu.org/software/make/manual/make.html#File-Name-Functions)
:   `$(dir NAMES...)`
:   `$(notdir NAMES...)`
:   `$(basename NAMES...)`
:   `$(addprefix PREFIX,NAMES...)`
:   `$(wildcard PATTERN)`
:   `$(abspath NAMES...)`

[条件分岐](https://www.gnu.org/software/make/manual/make.html#Conditional-Functions)
:   `$(if CONDITION,THEN,ELSE)`

    [関数じゃない条件分岐](https://www.gnu.org/software/make/manual/make.html#Conditional-Syntax)
    (`ifeq`, `ifneq`, `ifdef`, `ifndef`, `else`, `endif`) もある。

その他
: `$(foreach VAR,LIST,TEXT)`:
  `LIST` の中身をそれぞれ `VAR` に入れて `TEXT` を実行
: `$(file op FILENAME,TEXT)`:
  `text` の結果をファイルに書き出す
: `$(call VARIABLE,PARAMS...)`:
  `$(1) $(2)` などを使って定義しておいた `VARIABLE` を関数のように呼び出す
: `$(origin VARIABLE)`:
  変数がどう定義されたかを知れる
: `$(error TEXT...)`, `$(warning TEXT...)`, `$(info TEXT...)`:
  エラーや警告をプリントする
: `$(shell COMMAND...)`:   シェルを呼び出す

### [Targets](http://www.gnu.org/software/make/manual/make.html#Standard-Targets)

`all`
:   ディレクトリ内のcppソースをコンパイル

`clean`
:   コンパイル済みオブジェクトを一掃

`_`
:   v3.81以降であれば `.DEFAULT_GOAL` が効くので `make all` と同じ

```sh
% make clean
% make
```

## [Options](https://www.gnu.org/software/make/manual/make.html#Options-Summary)

`-f file`
:   `Makefile` じゃない名前のファイルを指定したければ

`-j jobs`
:   並列コンパイル。コア数+1くらいがちょうどいいらしい

`-C directory`
:   そのディレクトリに行って `make`

`-p`
:   自動的に作られるものも含めてすべての変数を表示
