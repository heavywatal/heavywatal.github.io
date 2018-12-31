+++
title = 'Qt'
subtitle = "多OS対応GUIフレームワーク"
tags = ["c++"]
[menu.main]
  parent = "dev"
+++

<http://qt-project.org/>

<http://qt-project.org/doc/qt-5/reference-overview.html>

## インストール (Mac OS X)

<http://qt-project.org/doc/qt-5/macosx-building.html>

<http://doc.qt.digia.com/qtcreator-extending/getting-and-building.html>

1.  Qt Library を `developer` オプション付きでインストール。
    例えば Homebrew なら:

        brew install qt5 --developer
        brew linkapps

1.  Qt Creator のソースコードを
    [公式サイト](http://qt-project.org/downloads)
    から落として展開:

        wget -O- http://download.qt-project.org/official_releases/qtcreator/3.1/3.1.0/qt-creator-opensource-src-3.1.0.tar.gz | tar xz

1.  Qt Library でインストールされた qmake を使ってビルド:

        cd qt-creator-opensource-src-3.1.0/
        mkdir build
        cd build/
        qmake -r ..
        make
        mv bin/Qt\ Creator.app /Applications/

1.  Qt Creator を立ち上げて環境設定を開き
    `Build and Run --> Qt Version --> Add...`
    から qmake のパスを指定。例えば:

        ${HOME}/.homebrew/opt/qt5/bin/qmake

1.  隣のタブの `Kits --> Manual --> Desktop (default) --> Qt version`
    を確認する。そのほかの項目も適当に。

{{%div class="note"%}}
古い `QMAKEFEATURES` 設定が残ってたりすると失敗するので
再インストールを試みるときなどは注意

```sh
Failed to process makespec for platform 'macx-clang'
ASSERT: "fileName.isEmpty() || isAbsolutePath(fileName)" in file */qtbase/qmake/library/ioutils.cpp, line 61
```
{{%/div%}}

## qwt --- Qt用のグラフ描画プラグイン

<http://qwt.sourceforge.net/>

### インストール

<http://qwt.sourceforge.net/qwtinstall.html>

<http://qt-project.org/doc/qt-5.1/qtdoc/deployment-plugins.html>

<https://qt-project.org/doc/qtcreator-3.1/adding-plugins.html>

{{%div class="note"%}}
Qt Library のバージョンを変えたらプラグインもビルドし直すべし。
{{%/div%}}

1.  [プロジェクトページ](http://sourceforge.net/projects/qwt/files/qwt/)
    からソースコードをダウンロードして展開:

        wget -O- http://sourceforge.net/projects/qwt/files/latest/download?source=files | tar xj

1.  インストール済みのQtを使ってビルド:

        cd qwt-6.1.0/
        qmake qwt.pro
        make

1.  `/usr/local/qwt-{VERSION}/` にインストールされる:

        sudo make install

1.  Qt Creator から見えるところにプラグインをシムリンク:

        ln -s /usr/local/qwt-6.1.0/plugins/designer/libqwt_designer_plugin.dylib $(brew --prefix)/opt/qt5/plugins/designer

1.  dylibの依存関係を確認し `qwt` の部分を絶対パスに書き換える:

        otool -L /usr/local/qwt-6.1.0/plugins/designer/libqwt_designer_plugin.dylib
        sudo install_name_tool -change qwt.framework/Versions/6/qwt /usr/local/qwt-6.1.0/lib/qwt.framework/qwt /usr/local/qwt-6.1.0/plugins/designer/libqwt_designer_plugin.dylib
        otool -L /usr/local/qwt-6.1.0/plugins/designer/libqwt_designer_plugin.dylib

1.  Qt Creator で適当なプロジェクトのuiファイルを開き、
    Design タブに Qwt Widgets が登場していれば成功。

    無い場合はその画面のメニューバーから
    `Tools --> Form Editor --> About Qt Designer Plugins...`
    を見て原因を探す。

    {{%div class="note"%}}
例えばこんなエラーなら、
プラグインそのものは Qt Creator から見えているが、
`qwt` ライブラリが読み込めていない

```sh
Library not loaded: qwt.framework/Versions/6/qwt
Referenced from: /usr/local/qwt-6.1.0/plugins/designer/libqwt_designer_plugin.dylib
Reason: image not found
```

これすらも無い場合はプラグインがそもそも見えていない。
    {{%/div%}}

### ビルド方法

1.  `.pro` ファイルにこんな感じで追加:

        include ( /usr/local/qwt-6.1.0/features/qwt.prf )

1.  Qt Creator 左下の緑三角でビルド＆ランしてみると怒られる:

        dyld: Library not loaded: qwt.framework/Versions/6/qwt
          Referenced from: /Users/***/build-MyApp-Desktop-Release/MyApp.app/Contents/MacOS/MyApp
          Reason: image not found
        The program has unexpectedly finished.

1.  実行ファイルを見てみると、こいつも相対パスでライブラリ参照している:

        otool -L /path/to/MyApp.app/Contents/MacOS/MyApp
                qwt.framework/Versions/6/qwt (compatibility version 6.1.0, current version 6.1.0)

1.  絶対パスに修正:

        install_name_tool -change qwt.framework/Versions/6/qwt /usr/local/qwt-6.1.0/lib/qwt.framework/Versions/6/qwt build-MyApp-Desktop-Release/MyApp.app/Contents/MacOS/MyApp

1.  Qt Creator 左下の緑三角で再びランしてみる。

    {{%div class="note"%}}
この段階で、開発マシンで起動できるようにはなっているはず。
しかしQtも含めあらゆるライブラリを絶対パスで参照してしまっているので、
このままほかのマシンにコピーしても起動できない。
    {{%/div%}}

1.  ライブラリをアプリケーションバンドル内にコピーし、参照先を変更する:

        macdeployqt build-MyApp-Desktop-Release/MyApp.app
