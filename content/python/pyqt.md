+++
title = 'PyQt'
subtitle = "Python bindings for Qt framework"
tags = ["python"]
[menu.main]
  parent = "python"
+++

<http://www.riverbankcomputing.co.uk/software/pyqt/>
<http://pyqt.sourceforge.net/Docs/PyQt5/index.html>

クロスプラットフォームなGUIアプリケーションを作るためのC++ツールキットである
`Qt` をPythonから利用するためのライブラリ。
ほかのPythonバインディングとしては
`Qt` プロジェクト公式の `PySide` があるけど、
なぜかそれよりもこちらのほうが早く Qt 5 や Python 3 に対応したりしてる。
単独で実行可能なファイルに固めるには `PyInstaller` を使う。

## Installation

Linuxなら `apt-get` を使ってもよい。

1.  `Qt` Library をインストール
2.  SIP をインストール <http://www.riverbankcomputing.co.uk/software/sip/download>:

        % tar xzf sip-4.13.2.tar.gz
        % cd sip-4.13.2/
        % /usr/local/bin/python2.7 configure.py --arch=x86_64
        % make
        % sudo make install

    Macで `--enable-framework` を付けずにビルドして使いたい場合は、
    `sipconfig.py` と `siputils.py` にある以下の行をコメントアウトする:

        #            if "Python.framework" not in dl:
        #                error("SIP requires Python to be built as a framework")

3.  `PyQt` をインストール <http://www.riverbankcomputing.co.uk/software/pyqt/download>:

        % tar xzf PyQt-mac-gpl-4.9.1.tar.gz
        % cd PyQt-mac-gpl-4.9.1/
        % /usr/local/bin/python2.7 configure.py --use-arch=x86_64
        % make
        % sudo make install

## Example

    #!/usr/bin/python
    # -*- coding: utf-8 -*-
    """
    """
    import sys

    from PyQt4 import QtGui, QtCore

    class MainWindow(QtGui.QMainWindow):
        def __init__(self, *args):
            super(MainWindow, self).__init__(*args)
            self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
            self.setWindowTitle("Window Title")
            self.resize(480, 240)

            self.statusBar()

            self.centralwidget = QtGui.QWidget(self)
            hbox = QtGui.QHBoxLayout(self.centralwidget)
            hbox.addWidget(QtGui.QPushButton("Close", self, clicked=self.close))
            hbox.addWidget(QtGui.QPushButton("Quit", self, clicked=app.quit))
            hbox.addWidget(QtGui.QPushButton("Exit", self, clicked=sys.exit))
            self.setCentralWidget(self.centralwidget)

    if __name__ == '__main__':
        app = QtGui.QApplication(sys.argv)
        win = MainWindow()
        win.show()
        win.raise_()
        sys.exit(app.exec_())

## QThread

アプリケーションをマルチスレッド化するには
`QtCore.QThread` クラスを継承して `run` をオーバーライドする。
素直にそうするとMacでは以下のようなエラーが出てくる。:

    Python[31533:53e7] *** __NSAutoreleaseNoPool(): Object 0x11e01f6c0 of class NSConcreteMapTable autoreleased with no pool in place - just leaking

`PyQt` に限らずCocoaアプリのマルチスレッド化には
`NSAutoreleasePool` の呪文が必要らしいので以下のようにする。:

    import QtCore

    if sys.platform.startswith("darwin"):
        import Foundation

    class Thread(QtCore.QThread):
        def __init__(self, parent=None, target=None):
            super(Thread, self).__init__(parent)
            self.target = target

        def run(self):
            if sys.platform.startswith("darwin"):
                pool = Foundation.NSAutoreleasePool.alloc().init()
            self.target()
            return

## i18n

1.  ソースファイル内の翻訳対象文字列を `self.tr()` で囲む
    (`tr` は `QObject` を継承してるクラスなら持ってるはず):

        self.tr("Gene")

2.  プロジェクトファイル(.pro)を作り、ソースファイルと言語ファイルを指定:

        SOURCES = gui.py
        TRANSLATIONS = translations/ja_JP.ts

3.  プロジェクトファイルに書かれたソースから翻訳対象を抜き出して
    `ts` ファイルを生成:

        % pylupdate4-2.7 pyqt4.pro

4.  Qt Linguist で `ts` ファイルを開き、翻訳。
    ただのXMLなのでテキストエディタで直接編集しても良い
5.  `ts` ファイルをコンパイルして `qm` ファイルを作る
    (Linguist のメニューからもできる):

        % lrelease translations/ja_JP.ts

6.  `qrc` ファイルに `qm` ファイルを登録する:

        <!DOCTYPE RCC><RCC version="1.0">
        <qresource>
        <file>translations/ja_JP.qm</file>
        </qresource>
        </RCC>

7.  `qrc` ファイルを `py` モジュールに変換する:

        % pyrcc4-2.7 translations.qrc -o translations.py

8.  ソースファイルで `py` モジュールを `import` し、
    ファイル名をリソースとして記述。
    実際にはシステムあるいは環境変数のロケールを反映するために
    `QtCore.QLocale().system().name()` などを使うと良い。:

        if __name__ == '__main__':
            app = QtGui.QApplication(sys.argv)

            import translations
            translator = QtCore.QTranslator()
            translator.load(":translations/ja_JP.qm")

            app.installTranslator(translator)

            win = MainWindow()
            win.show()
            win.raise_() # for Mac
            sys.exit(app.exec_())

## `PySide`

`Qt` プロジェクト公式のPythonバインディング。
既に広く使われている `PyQt`
とのAPI互換性を維持しつつLGPLで公開されている。
最近は開発が進んでいないようで、Qt5とかPython3への対応はまだ。

-   <http://qt-project.org/wiki/PySide>
-   <http://qt-project.org/wiki/PySideDocumentation>
-   <http://pyside.github.io/docs/pyside/>

### Installation

<http://qt-project.org/wiki/PySideDownloads/>

Macなら公式サイトからインストーラ `pyside-*.pkg`
をダウンロードしてコマンドから実行:

    % cd ~/Downloads/
    % sudo installer -pkg pyside-0.4.1-qt47-py26apple.pkg -target "/"

Ubuntuなら [PPA](https://launchpad.net/~pyside/+archive/ppa)
からインストール:

    % sudo add-apt-repository ppa:pyside
    % sudo apt-get update
    % sudo apt-get install python-pyside

## `PyInstaller`

<http://www.pyinstaller.org/>\
<http://pythonhosted.org/PyInstaller/>

Pythonスクリプトを単独で実行可能なアプリケーションに変換するプログラム。
Linux, Mac, Windowsに対応。

昔はインストールも使用もいろいろ難しい手続きが必要だったけど
今はかなり使いやすくなってる。

[pip]({{< relref "pip.md" >}}) で一撃インストール:

    % pip install pyinstaller

使う時もコマンドひとつ:

    % pyinstaller --windowed myqtapp.py
