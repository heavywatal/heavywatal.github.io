+++
tags = ["package", "c++"]
title = "CMake"
subtitle = "Cross-platform Make"
[menu.main]
  parent = "dev"
+++

[autotools]({{< relref "autotools.md" >}}) の `configure` 的な位置づけで、
環境に合わせた [Makefile]({{< relref "make.md" >}}) を自動生成する。

https://cmake.org/cmake/help/latest/


## `CMakeLists.txt`

ディレクトリ毎に配置して、階層的に制御するのが好まれる。
プロジェクトのトップに置くものは、以下のようなコマンドで始める必要がある。

```cmake
cmake_minimum_required(VERSION 3.1)
project(helloworld CXX)
```

`add_*()` でビルドターゲットを作成し、
`target_*()` でそれらの設定を整えて、
`install()` でインストールする対象や行き先を指定する、というのが基本の流れ。

```cmake
add_executable(a.out hello.cpp)
target_compile_definitions(a.out PRIVATE -DNDEBUG)
target_compile_options(a.out PRIVATE -Wall -Wextra -Wpedantic -O3)
install(TARGETS a.out
  RUNTIME DESTINATION bin
)
```


## [Commands](https://cmake.org/cmake/help/latest/manual/cmake-commands.7.html)

### Scripting commands

- `foreach()`
- `include(GNUInstallDirs)`
- `macro()`
- `message(STATUS Hello world)`
- `option()`
- `set(name value)`
- `file(GLOB srcfiles *.cpp)`:
  グロブにマッチするファイルを列挙して変数に格納


### [Project commands](https://cmake.org/cmake/help/latest/manual/cmake-commands.7.html#id4)

サブディレクトリを利用する:

- `add_subdirectory(source_dir [binary_dir] [EXCLUDE_FROM_ALL])`
- `aux_source_directory(<dir> VAR)`

ターゲットを定義する:

- [`add_executable(<name> [EXCLUDE_FROM_ALL] ...)`]
  (https://cmake.org/cmake/help/latest/command/add_executable.html)
- [`add_library(<name> [STATIC|SHARED|OBJECT] [EXCLUDE_FROM_ALL|IMPORTED] ...)`]
  (https://cmake.org/cmake/help/latest/command/add_library.html)

ターゲットのプロパティを追加する:

- [`target_compile_definitions(<target> <INTERFACE|PRIVATE|PUBLIC> ...)`]
  (https://cmake.org/cmake/help/latest/command/target_compile_definitions.html)
- [`target_compile_options(<target> [BEFORE] <I|P|P> ...)`]
  (https://cmake.org/cmake/help/latest/command/target_compile_options.html)
- [`target_include_directories(<target> [SYSTEM] [BEFORE] <I|P|P> ...)`]
  (https://cmake.org/cmake/help/latest/command/target_include_directories.html)
- [`target_link_libraries(<target> <I|P|P> ...)`]
  (https://cmake.org/cmake/help/latest/command/target_link_libraries.html)
- ターゲットなしの `include_directories()` `link_directories()` `link_libraries()`
  などは全体に影響が及ぶ亜種。

[インストールするものや宛先を指定する。](https://cmake.org/cmake/help/latest/command/install.html)

- `install(TARGETS)`
- `install(<FILES|PROGRAMS>)`
- `install(DIRECTORY)`
- `install(EXPORT)`: 外部のCMakeから使いやすくするためのconfigをインストールする。
  似て非なる `export()` はビルドツリーにあるものを使わせるための謎コマンド。



### オプション

`PRIVATE`
: このターゲットをビルドするときだけ使い、これを利用するときには参照させない。
  例えば、このプロジェクトのライブラリをビルドするにはBoostが必要だけど、
  これを利用するときにそれらのパスを知る必要はない、とか。

`INTERFACE`
: このターゲットでは使わないけど、これを利用するときには参照させる。
  例えば、ヘッダーライブラリを作る場合とか。

`PUBLIC`
: このターゲットにもこれを利用するターゲットにも使う。使う場面あるかな？

`EXCLUDE_FROM_ALL`
: `make [all]` から外れて、明示的なターゲット指定でのみビルドされるようになる。


## [Variables](https://cmake.org/cmake/help/latest/manual/cmake-variables.7.html)

### Variables that Provide Information

- `CMAKE_BINARY_DIR`: ビルドツリーの最上階。
- `CMAKE_COMMAND`: フルパス`cmake`
- `CMAKE_CURRENT_*_DIR`: `add_subdirectory()` の先で使う。
- `CMAKE_PROJECT_NAME`: `project()` で設定したやつ。
- `CMAKE_SKIP_RPATH`
- `CMAKE_SOURCE_DIR`: 最上階の`CMakeLists.txt`があるとこ。
- `CMAKE_VERBOSE_MAKEFILE`: とりあえず `TRUE`

### Variables that Change Behavior

- `CMAKE_BUILD_TYPE`: Debug, Release, RelWithDebInfo, MinSizeRel.
- `CMAKE_INSTALL_PREFIX`: configureでの`--prefix`に相当
- `CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT`
- `CMAKE_PREFIX_PATH`: パッケージやファイル探索 `find_*()` の候補パスを追加する

### Variables that Describe the System

`APPLE`, `UNIX`, `WIN32`

### Variables that Control the Build

- `CMAKE_BUILD_RPATH`
- `CMAKE_BUILD_WITH_INSTALL_NAME_DIR`
- `CMAKE_BUILD_WITH_INSTALL_RPATH`
- `CMAKE_INSTALL_NAME_DIR`
- `CMAKE_INSTALL_RPATH`
- `CMAKE_INSTALL_RPATH_USE_LINK_PATH`
- `CMAKE_MACOSX_RPATH`

### C++

```cmake
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-g)
else()
  add_definitions(-DNDEBUG)
endif()
add_compile_options(-O3 -march=native -Wall -Wextra -Wpedantic)
```

- `CMAKE_CXX_COMPILER`

## Modules

[ExternalProject](https://cmake.org/cmake/help/latest/module/ExternalProject.html)

### [Boost](https://cmake.org/cmake/help/latest/module/FindBoost.html)

```cmake
set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost REQUIRED COMPONENTS program_options iostreams filesystem system)
message("Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")

target_link_libraries(mytarget ${Boost_LIBRARIES})
```

`BOOST_ROOT`
: 探索パス。e.g. `$(brew --prefix)`, `${HOME}/local`


### [CTest](https://cmake.org/cmake/help/latest/module/CTest.html)

```cmake
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
```

`enable_testing()` と書くほうが短いけど
`BUILD_TESTING=TRUE` 固定なのでオプションで切れない。

```cmake
# test/CMakeLists.txt
add_executable(test-gene gene.cpp)
add_test(NAME gene COMMAND $<TARGET_FILE:test-gene>)
```

`ctest -V` で実行。
一部のテストのみ実行したいときは `-R <pattern>` で絞る。


## CLI

### [`cmake`](https://cmake.org/cmake/help/latest/manual/cmake.1.html)

いろんな中間ファイルができる上に `cmake clean` は無いので、
ビルド用の空ディレクトリを外に作って out-of-source で実行するのが基本。
やり直したいときは、そのディレクトリごと消す。

```sh
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/local -DCMAKE_BUILD_TYPE=Debug /path/to/project
```

デフォルトでは `Makefile` が書き出されるので
`make && make install` のように実行してもいいけど、
`cmake` からそれを実行することもできる:

```sh
cmake --build . -- -j2
cmake --build . --target install
```

オプション

`-DCMAKE_XXX=YYY`
: コマンドラインから変数を設定する。

`-G <generator-name>`
: Makefile, Ninja, Xcode, etc.

`-E <subcommand>`
: `chdir <dir> <cmd>`
: `make_directory <dir>`

`-H <dir>` (*undocumented*)
: ソースツリーを指定する。

`-B <dir>` (*undocumented*)
: ビルドツリーを指定する。

`-L`
: キャッシュされている変数をリストアップ。
  `H` をつけると説明文も。
  `A` をつけるとadvancedな変数も。
  見るだけなら `-N` オプションと共に。
