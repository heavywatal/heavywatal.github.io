+++
tags = ["package", "c++"]
title = "CMake"
subtitle = "Cross-platform Make"
[menu.main]
  parent = "dev"
+++

環境に合わせた [Makefile]({{< relref "make.md" >}}) を自動生成する。
似たようなことをする `configure` スクリプトと比べて動作が高速で、
ライブラリの依存関係なども簡潔・柔軟に記述できる。

`configure` ではそれを生成する開発者だけが
[autotools]({{< relref "autotools.md" >}}) を使うのに対して、
CMakeでは開発者と利用者の双方がCMakeをインストールして使う。

https://cmake.org/cmake/help/latest/


## 基本

`CMakeLists.txt` を各ディレクトリに配置して、階層的に管理する。
プロジェクトのトップに置くものは、以下のようなコマンドで始める必要がある。

```cmake
cmake_minimum_required(VERSION 3.15)
project(helloworld
  VERSION 0.1.0
  LANGUAGES CXX)
```

`add_executable()` や `add_library()` でビルドターゲットを作成し、
`target_*()` でそれらの設定を整えて、
`install()` でインストールする対象や行き先を指定する、というのが基本の流れ。

```cmake
add_executable(a.out hello.cpp)
target_compile_options(a.out PRIVATE -Wall -Wextra -pedantic)
install(TARGETS a.out
  RUNTIME DESTINATION bin
)
```

[`cmake` コマンドの使い方は後述](#cli)


## Commands

<https://cmake.org/cmake/help/latest/manual/cmake-commands.7.html>

### Scripting commands

- `configure_file(<input> <output> [COPYONLY] [@ONLY])`:
  ファイルの一部を置換しつつ複製する。
  例えば `@PROJECT_VERSION@` などを含む `config.hpp.in` に値を埋め込んで `config.hpp` を生成するとか。
- `function(<name> [args...])`
- `foreach(var IN LISTS list)`
- `message(STATUS "Hello world!")`
- `option(VARIABLE "Message" ON)`
- `set(VARIABLE value)`
- [`cmake_path()`](https://cmake.org/cmake/help/latest/command/cmake_path.html)
  パス操作全般。
  3.20 より古い環境では `get_filename_component()` 。


### Project commands

<https://cmake.org/cmake/help/latest/manual/cmake-commands.7.html#project-commands>

外部のCMakeプロジェクトを利用する:

- [`find_package(<name> [version] [REQUIRED] ...)`](https://cmake.org/cmake/help/latest/command/find_package.html):
  ビルド・インストール済みのを取り込む。
- [`add_subdirectory(source_dir [binary_dir] [EXCLUDE_FROM_ALL])`](https://cmake.org/cmake/help/latest/command/add_subdirectory.html):
  インストール前のソースツリーを取り込む。

ターゲットを定義する:

- [`add_executable(<name> [EXCLUDE_FROM_ALL] ...)`](https://cmake.org/cmake/help/latest/command/add_executable.html)
- [`add_library(<name> [STATIC|SHARED|OBJECT] [EXCLUDE_FROM_ALL|IMPORTED] ...)`](https://cmake.org/cmake/help/latest/command/add_library.html)
  | option | Mac | Linux | Win | description |
  | ------ | --- | ----- | --- | ----------- |
  | `STATIC` | `.a` | `.a` | `.lib` | オブジェクトをまとめて目次をつけたarchive。静的リンクで複製は生じるが最適化される。 |
  | `SHARED` | `.dylib` | `.so` | `.dll` | オブジェクトをまとめて実行可能ファイル一歩手前まで加工したもの。動的リンクでストレージにもメモリにも1つだけあればいいので効率はいいけどRPATHのお世話が必要。 |
  | `OBJECT` | `.o` | `.o` | `.obj` | コンパイル単位のオブジェクト。ライブラリを経由しないぶん`STATIC`より効率的かと思いきや、リンク時に最適化されず、使わないコードまで保持される。使い道はライブラリを複数作りたい時くらい？ |


ターゲットのプロパティを追加する:

- [`set_target_properties(<target> PROPERTIES p1 v1 p2 v2 ...)`](https://cmake.org/cmake/help/latest/command/set_target_properties.html)
- [`target_compile_definitions(<target> <INTERFACE|PRIVATE|PUBLIC> ...)`](https://cmake.org/cmake/help/latest/command/target_compile_definitions.html)
- [`target_compile_features(<target> <P|P|I> ...)`](https://cmake.org/cmake/help/latest/command/target_compile_features.html)
- [`target_compile_options(<target> [BEFORE] <I|P|P> ...)`](https://cmake.org/cmake/help/latest/command/target_compile_options.html)
- [`target_sources(<target> <I|P|P> ...)`](https://cmake.org/cmake/help/latest/command/target_sources.html)
- [`target_include_directories(<target> [SYSTEM] [BEFORE] <I|P|P> ...)`](https://cmake.org/cmake/help/latest/command/target_include_directories.html):
  次の関数があるおかげでこれを直接使うことは意外と少ない。
- **[`target_link_libraries(<target> <I|P|P> ...)`](https://cmake.org/cmake/help/latest/command/target_link_libraries.html):**
  この関数でターゲット間の依存関係を繋げていくのがCMakeの肝。
  ライブラリ側 `-L -l` だけではなく、インクルード側 `-I` のオプションもお世話してくれる。
- ターゲットなしの `include_directories()` `link_directories()` `link_libraries()`
  などはディレクトリ単位で影響が及ぶ亜種で、非推奨。

[インストールするものや宛先を指定する。](https://cmake.org/cmake/help/latest/command/install.html)

- `install(TARGETS)`
- `install(<FILES|PROGRAMS>)`
- `install(DIRECTORY)`
- `install(EXPORT)`: 外部のCMakeから使いやすくするためのconfigをインストールする。
  似て非なる `export()` はビルドツリーにあるものを使わせるための謎コマンド。



### オプション

`PRIVATE`
: このターゲットをビルドするときだけ使い、これを利用するときには参照させない。
  例えば「このプロジェクトのライブラリをビルドするにはBoostヘッダーが必要だけど、
  これを利用するときにそれらのパスを知る必要はない」とか。

`INTERFACE`
: このターゲットでは使わないけど、これを利用するときには参照させる。
  例えば、ヘッダーライブラリを作る場合とか。

`PUBLIC`
: このターゲットにもこれを利用するターゲットにも使う。使う場面あるかな？

`EXCLUDE_FROM_ALL`
: `make [all]` から外れて、明示的なターゲット指定でのみビルドされるようになる。


## Variables

<https://cmake.org/cmake/help/latest/manual/cmake-variables.7.html>

### Variables that Provide Information

- `PROJECT_SOURCE_DIR`, `PROJECT_BINARY_DIR`:
  直近の `project()` におけるソースツリー、ビルドツリーの最上階。
- `CMAKE_SOURCE_DIR`, `CMAKE_BINARY_DIR`:
  ソースツリー、ビルドツリーの最上階。
  他のプロジェクトから `add_subdirectory()` で使われる場合、
  `PROJECT_*` よりも上になる。
- `CMAKE_CURRENT_SOURCE_DIR`, `CMAKE_CURRENT_BUILD_DIR`:
  ソースツリー、ビルドツリーにおける現在地。
- `PROJECT_NAME`, `PROJECT_VERSION`: `project()` で設定したやつ。
- `CMAKE_SKIP_RPATH`
- `CMAKE_VERBOSE_MAKEFILE`: とりあえず `ON`

### Variables that Change Behavior

- `BUILD_SHARED_LIBS`:
  STATIC/SHAREDが明示されていない `add_library()` でどっちをビルドするか。
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


### Environment variables

<https://cmake.org/cmake/help/latest/manual/cmake-env-variables.7.html>

変数として自動的に利用可能になったりはせず、
[`$ENV{VAR}`](https://cmake.org/cmake/help/latest/variable/ENV.html)
みたいな形で参照するのが基本。

ただし、一部の環境変数はCMakeの変数の初期値として採用される。
e.g., `CMAKE_PREFIX_PATH`, `CXX`, `<PackageName>_ROOT`, etc.

何が渡っているかは `cmake -E environment` で確認できる。

[`CMAKE_EXPORT_COMPILE_COMMANDS`](https://cmake.org/cmake/help/latest/envvar/CMAKE_EXPORT_COMPILE_COMMANDS.html)
: 定義しておくとコンフィグ時に `compile_commands.json` を生成してもらえる。
  これで [clangd](https://clangd.llvm.org/) にコンパイルオプションを伝えられる。
  ソースファイルの親ディレクトリを辿るだけでなく
  `build` という名のサブディレクトリも探してくれるので
  `-B build` の慣習に従っていればコピーやシムリンクも不要。


### C++

```cmake
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
  WINDOWS_EXPORT_ALL_SYMBOLS ON
)
target_compile_options(common PRIVATE
  -Wall -Wextra -pedantic
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>:-march=native>
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>:-march=armv8.3-a+sha3>
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
cmake_print_variables(CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS_DEV "-O2 -g")
```

`CMAKE_CXX_*` のようなグローバル設定を使わず
`target_*()` でターゲットごとに設定するのが今後の主流。
[`CMAKE_CXX_KNOWN_FEATURES`](https://cmake.org/cmake/help/latest/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html)
に `cxx_std_17` などの便利なメタタグが導入されたのは CMake 3.8 から。

Predefined variable              | default
---------------------------------|----
`CMAKE_CXX_FLAGS`                |
`CMAKE_CXX_FLAGS_DEBUG`          | `-g`
`CMAKE_CXX_FLAGS_MINSIZEREL`     | `-Os -DNDEBUG`
`CMAKE_CXX_FLAGS_RELEASE`        | `-O3 -DNDEBUG`
`CMAKE_CXX_FLAGS_RELWITHDEBINFO` | `-O2 -g -DNDEBUG`

`#ifndef NDEBUG` なコードを残しつつ、
そこそこ速くコンパイル＆実行したい、
という組み合わせ `-O2 -g` は用意されていないので自分で定義する。
`CMAKE_CXX_FLAGS_???` を適当に作れば
`-DCMAKE_BUILD_TYPE=???` をcase-insensitiveに解釈してもらえる。


## Generator expressions

<https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html>

ビルド時の状態に応じてに変数を評価する仕組み。
コンフィグ時に評価される `if()` とは使い方が異なる。

例えば、プロジェクト内のビルド時と、外部パッケージとして利用される時とで、
インクルードパスを使い分ける。
```cmake
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
```

コンフィグ時には評価前の文字列でしかないので、
`if()` などの条件分岐にも使えないし、
`cmake_print_variables()` や `message()` しても中身は見えない。
実際にどんな値が入ったかを確かめるには
`file(GENERATE OUTPUT <outfile> CONTENT <content>)`
などでビルド時に書き出すことになる。


## Modules

<https://cmake.org/cmake/help/latest/manual/cmake-modules.7.html>

`include()` や `find_package()` から使う。

### include(CMakePrintHelpers)

変数の名前と中身を表示してくれる。
名前を2回書かずに済む。
次の二つは等価:
```cmake
cmake_print_variables(VAR)
message(STATUS VAR="${VAR}")
```

### GNUInstallDirs

<https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html>

インストール先のディレクトリを指定するときの標準的な値を決めてくれる。

Variable                    | Value
--------------------------- | -----
`CMAKE_INSTALL_BINDIR`      | `bin`
`CMAKE_INSTALL_INCLUDEDIR`  | `include`
`CMAKE_INSTALL_LIBDIR`      | `lib`, `lib64`
`CMAKE_INSTALL_DATADIR`     | `share`
`CMAKE_INSTALL_FULL_<dir>`  | `${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_<dir>}`


### CMakePackageConfigHelpers

<https://cmake.org/cmake/help/latest/module/CMakePackageConfigHelpers.html>

他のプロジェクトから以下のように利用されるライブラリを作りたい。
外部プロジェクトであることを明確にするため
`名前空間::ターゲット` という形でリンクするのが筋:
```cmake
project(OtherProject CXX)
find_package(MyLib)
target_link_libraries(OtherTarget PRIVATE MyLib::MyLib)
```

こうやって使ってもらうためには `*config.cmake` ファイルを
[`find_package()` の探索先](https://cmake.org/cmake/help/latest/command/find_package.html#config-mode-search-procedure)
のどこかにインストールする必要がある。
選択肢があり過ぎて悩ましく、公式ドキュメントにも推奨などは書かれていないが、
[識者のコメント](https://discourse.cmake.org/t/what-should-the-destination-be-for-a-header-only-librarys-cmake-config-file/8473)
によれば次の二択で使い分ける方針が良さそう:

- `share/cmake/${PROJECT_NAME}`:
  header-only, architecture-independent.
  `find_package()` は `${CMAKE_INSTALL_DATADIR}` に追従せず `share` を読みにいく。
- `${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}`:
  共有ライブラリあり、architecture-dependent.

`install(TARGETS)` の中で `EXPORT` のためのターゲット定義し、
`install(EXPORT)` でその設定を行う:
```cmake
project(MyLib
  VERSION 0.1.0
  LANGUAGES CXX)
# ...

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-config
)

set(config_destination ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(EXPORT ${PROJECT_NAME}-config
  DESTINATION ${config_destination}
  NAMESPACE ${PROJECT_NAME}::
)
```

バージョン情報も同じところに送り込む:
```cmake
set(version_file ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(${version_file}
  COMPATIBILITY AnyNewerVersion
)
install(FILES ${version_file}
  DESTINATION ${config_destination}
)
```

`find_package()` から使うにはここまでの設定で十分だが、
`add_subdirectory()` からでも同じ形で使えるようにするため、
`ALIAS` を設定しておいたほうがいい:
```cmake
add_library(${PROJECT_NAME})
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
```

上記のように直接 `EXPORT *-config` するのは簡易版。
ライブラリ利用時に必要な変数をインストール時に設定してあげたいときなどは、
同じ内容を `EXPORT *-targets` のような名前で書き出しておき、
`configure_package_config_file(Config.cmake.in, ...)`
のようなコマンドで `*-config.cmake` を生成し、
その中から `include(*-targets)` する。


### FetchContent

<https://cmake.org/cmake/help/latest/module/FetchContent.html>

外部ライブラリをコンフィグ時に取ってくる。
CMake 3.11 から。

```cmake
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
cmake_print_variables(FETCHCONTENT_SOURCE_DIR_IGRAPH)
FetchContent_Declare(
  igraph
  GIT_REPOSITORY https://github.com/igraph/igraph.git
  GIT_TAG ${PROJECT_VERSION}
  GIT_SHALLOW ON
)
FetchContent_MakeAvailable(igraph)
cmake_print_variables(igraph_SOURCE_DIR, igraph_BINARY_DIR)
```

`FetchContent_Declare()`
: まずこれで依存関係を宣言する。
  複数ある場合、先に全部宣言してからまとめてMakeAvailableを呼ぶのが推奨。
: ソースに関するオプションは
  [ExternalProject](https://cmake.org/cmake/help/latest/module/ExternalProject.html)
  とほぼ同じ。
: `FIND_PACKAGE_ARGS`: (3.24+)\
  `EXCLUDE_FROM_ALL`: (3.28+)

`FetchContent_MakeAvailable()`
: 宣言された依存ライブラリを利用可能な状態にする。(3.14+)
: これひとつを実行することが推奨されているが、各段階を手動で書くこともできる。
  1. `find_package()` を試みて、見つからなければ次に進む。(3.24+)
  1. `FetchContent_GetProperties()` で過去にPopulateしたものがあるか確認。
     `<lowercaseName>_POPULATED` が定義されていなければ次に進む。
  1. `FetchContent_Populate()` でソースコードを取得する。
     `FETCHCONTENT_SOURCE_DIR_<uppercaseName>`
     が定義されている場合はfetchせずそこにあるものを使う。
     成功したら3つの変数をセットする:
     - `<lowercaseName>_POPULATED`
     - `<lowercaseName>_SOURCE_DIR`
     - `<lowercaseName>_BINARY_DIR`
  1. `add_subdirectory()`
     でプロジェクトに取り込む。

似て非なる `ExternalProject` はビルド時に実行されるので
`add_subdirectory()` の対象にできず、
`execute_process()` で git を直接叩くなどして凌いでいた。

### FindThreads

<https://cmake.org/cmake/help/latest/module/FindThreads.html>

`-lpthread` とか自分で書かない。

```cmake
find_package(Threads)
target_link_libraries(MyTarget PRIVATE Threads::Threads)
```

### FindBoost

<https://cmake.org/cmake/help/latest/module/FindBoost.html>

Boostの特別扱いはdeprecatedになった。
普通のCMakeパッケージとして探すには、
`CONFIG` とか `NO_MODULE` オプションを足して明示的にConfigモードを使う。

```cmake
find_package(Boost CONFIG REQUIRED COMPONENTS context)
target_link_libraries(MyTarget PRIVATE Boost::context)
```

ヘッダーだけでいい場合は `Boost::boost` をリンクする。


### CTest

<https://cmake.org/cmake/help/latest/module/CTest.html>

```cmake
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
```

```cmake
# test/CMakeLists.txt
add_executable(test-gene gene.cpp)
add_test(NAME gene COMMAND $<TARGET_FILE:test-gene>)
```

`ctest -V` で実行。
一部のテストのみ実行したいときは `-R <pattern>` で絞る。


[`include(CTest)`](https://github.com/Kitware/CMake/blob/master/Modules/CTest.cmake)
は勝手に[CDash](https://www.cdash.org/)の設定をして
`DartConfiguration.tcl` を生成する。
次のように書き換えればそのへんをスキップできる:

```cmake
option(BUILD_TESTING "Build the testing tree." ON)
enable_testing()
```

## CLI

### `cmake`

<https://cmake.org/cmake/help/latest/manual/cmake.1.html>

ビルド用のディレクトリを別に作って out-of-source で実行するのが基本。
やり直したいときは、そのディレクトリごと消す。
3.0以降 `cmake --target clean` はあるが、
CMakeのバージョンを上げたときなどcleanしたい場面で使えない。

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=${HOME}/local -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j2
cmake --install build
```

`-S <dir>`
: ソースツリーを指定する。
  3.13から。それまではundocumentedで `-H<dir>` という形だった。

`-B <dir>`
: ビルドツリーを指定する。
  3.13から。それまではundocumentedだった。

`-D <var>=<value>`
: コマンドラインから変数を設定する。

`-G <generator-name>`
: Makefile, Ninja, Xcode, etc.
: デフォルトでは `Makefile` が書き出されるので
  `make && make install` と書いてもいいけど、
  generator非依存の `cmake --build` を使ったほうがいい。
: `cmake --install <dir>` が使えるのは3.15以降。
  それまでは `cmake --build build --target install` と明示する必要があった。

`-E <command>`
: シェルの違いを気にせず基本的なコマンドが使えるように。e.g.,
: `chdir <dir> <cmd>`
: `make_directory <dir>`

`-L`
: キャッシュされている変数をリストアップ。
  `H` をつけると説明文も。
  `A` をつけるとadvancedな変数も。
  見るだけなら `-N` オプションと共に。


## Versions

- 3.28: `FetchContent_Declare(... EXCLUDE_FROM_ALL)`,
  [C++20 modules](https://cmake.org/cmake/help/latest/manual/cmake-cxxmodules.7.html),
  [Ubuntu 24.04 noble](https://launchpad.net/ubuntu/noble/+source/cmake)
- 3.24: `FetchContent_Declare(... FIND_PACKAGE_ARGS)`, `--fresh`
- 3.23: `FILE_SET`
- 3.22: [Ubuntu 22.04 jammy](https://launchpad.net/ubuntu/jammy/+source/cmake)
- 3.20: [`cmake_path()`](https://cmake.org/cmake/help/latest/command/cmake_path.html), `cxx_std_23`
- 3.19: [`CMakePresets.json`](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)
- 3.17: [`CMAKE_EXPORT_COMPILE_COMMANDS`](https://cmake.org/cmake/help/latest/envvar/CMAKE_EXPORT_COMPILE_COMMANDS.html)
- 3.16: [Ubuntu 20.04 focal](https://launchpad.net/ubuntu/focal/+source/cmake)
- 3.15: `cmake --install`
- 3.13: `cmake -S . -B build`, `target_link_directories`, `target_link_options`
- 3.12: `cxx_std_20`, linking `OBJECT` libraries
- 3.11: `include(FetchContent)`
- 3.8: `cxx_std_17`