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
cmake_minimum_required(VERSION 3.1)
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


### Project commands

<https://cmake.org/cmake/help/latest/manual/cmake-commands.7.html#project-commands>

サブディレクトリを利用する:

- [`add_subdirectory(source_dir [binary_dir] [EXCLUDE_FROM_ALL])`](https://cmake.org/cmake/help/latest/command/add_subdirectory.html)

ターゲットを定義する:

- [`add_executable(<name> [EXCLUDE_FROM_ALL] ...)`](https://cmake.org/cmake/help/latest/command/add_executable.html)
- [`add_library(<name> [STATIC|SHARED|OBJECT] [EXCLUDE_FROM_ALL|IMPORTED] ...)`](https://cmake.org/cmake/help/latest/command/add_library.html)

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

### C++

```cmake
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)
set_target_properties(${PROJECT_NAME} PROPERTIES
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
  WINDOWS_EXPORT_ALL_SYMBOLS ON
)
target_compile_options(${PROJECT_NAME} PRIVATE
  -march=native -Wall -Wextra -pedantic
)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_DEV "-O2 -g")
```

`CMAKE_CXX_*` のようなグローバル設定を使わず
`target_*()` でターゲットごとに設定するのが今後の主流。
[`CMAKE_CXX_KNOWN_FEATURES`](https://cmake.org/cmake/help/latest/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html)
に `cxx_std_14` などの便利なメタタグが導入されたのは CMake 3.8 から。

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

文脈に応じて変数を評価する仕組み。

プロジェクト内のビルド時と、外部パッケージとして利用される時とで、
インクルードパスを使い分ける。
```cmake
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
```

## Modules

<https://cmake.org/cmake/help/latest/manual/cmake-modules.7.html>

`include()` や `find_package()` から使う。

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

Linuxでは `lib64` が標準的だが、Linuxbrewで使うには `lib` にする必要がある。
```cmake
if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()
```

### CMakePackageConfigHelpers

<https://cmake.org/cmake/help/latest/module/CMakePackageConfigHelpers.html>

他のプロジェクトから以下のように利用されるライブラリを作りたい。
外部プロジェクトであることを明確にするため
`名前空間::ターゲット` という形でリンクするのが筋:
```cmake
project(otherproject CXX)
find_package(mylib)
# add_subdirectory(mylib)
target_link_libraries(othertarget PRIVATE mylib::mylib)
```

`${CMAKE_INSTALL_PREFIX}/share/` らへんにconfigファイルを送り込む必要がある。
`install(TARGETS)` の中で `EXPORT` のためのターゲット定義し、
`install(EXPORT)` でその設定を行う:
```cmake
project(mylib
  VERSION 0.1.0
  LANGUAGES CXX)
# ...
install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-config
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(EXPORT ${PROJECT_NAME}-config
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
  NAMESPACE ${PROJECT_NAME}::
)
```

バージョン情報も同じところに送り込む:
```cmake
set(VERSION_CONFIG ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${VERSION_CONFIG} COMPATIBILITY AnyNewerVersion
)
install(FILES ${VERSION_CONFIG}
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)
```

`find_package()` から使うにはここまでの設定で十分だが、
`add_subdirectory()` からでも同じ形で使えるようにするため、
`ALIAS` を設定しておいたほうがいい:
```cmake
add_library(${PROJECT_NAME} SHARED)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
```

### FetchContent

<https://cmake.org/cmake/help/latest/module/FetchContent.html>

外部ライブラリを取ってきて配置する。
前からあった
[ExternalProject](https://cmake.org/cmake/help/latest/module/ExternalProject.html)
はビルド時に実行されるため
`add_subdirectory()` の対象にできないなどの問題があったが、
こちらはコンフィグ時に実行される。

```cmake
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
message(STATUS "FETCHCONTENT_SOURCE_DIR_IGRAPH: ${FETCHCONTENT_SOURCE_DIR_IGRAPH}")
FetchContent_Declare(
  igraph
  GIT_REPOSITORY https://github.com/igraph/igraph.git
  GIT_TAG ${PROJECT_VERSION}
  GIT_SHALLOW ON
)
FetchContent_MakeAvailable(igraph)
message(STATUS "igraph_SOURCE_DIR: ${igraph_SOURCE_DIR}")
```

CMake 3.11 からの新機能なので、もう少し普及するまでお預け。
当面は `execute_process()` で凌ぐ:
```cmake
find_package(Git)
execute_process(COMMAND
  ${GIT_EXECUTABLE} clone --recursive --depth=1 --branch=master https://github.com/USER/REPO.git ${SUBDIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
add_subdirectory(${SUBDIR} EXCLUDE_FROM_ALL)
```

### FindThreads

<https://cmake.org/cmake/help/latest/module/FindThreads.html>

`-lpthread` とか自分で書かない。

```cmake
find_package(Threads)
target_link_libraries(mytarget PRIVATE Threads::Threads)
```

### FindBoost

<https://cmake.org/cmake/help/latest/module/FindBoost.html>

```cmake
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED COMPONENTS filesystem)
target_link_libraries(mytarget PRIVATE Boost::filesystem)
```

ヘッダーだけでいい場合は `Boost::boost` ターゲットをリンクする。

探索パスを追加するには `BOOST_ROOT` を設定する。


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

ビルド用の空ディレクトリを外に作って out-of-source で実行するのが基本。
やり直したいときは、そのディレクトリごと消す。
3.0以降 `cmake --target clean` はあるが、
CMakeのバージョンを上げたときなどcleanしたい場面で使えない。

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=${HOME}/local -DCMAKE_BUILD_TYPE=Debug
```

デフォルトでは `Makefile` が書き出されるので
`make && make install` のように実行してもいいけど、
`cmake` からそれを実行することもできる:

```sh
cmake --build build -j 2
cmake --build build -j 2 --target install
```

3.15以降は `cmake --install <dir>` が使える。

`-S <dir>`
: ソースツリーを指定する。
  3.13から。それまではundocumentedで `-H<dir>` という形だった。

`-B <dir>`
: ビルドツリーを指定する。
  3.13から。それまではundocumentedだった。

`-DCMAKE_XXX=YYY`
: コマンドラインから変数を設定する。

`-G <generator-name>`
: Makefile, Ninja, Xcode, etc.

`-E <command>`
: シェルの違いを気にせず基本的なコマンドが使えるように。e.g.,
: `chdir <dir> <cmd>`
: `make_directory <dir>`

`-L`
: キャッシュされている変数をリストアップ。
  `H` をつけると説明文も。
  `A` をつけるとadvancedな変数も。
  見るだけなら `-N` オプションと共に。
