+++
title = 'DnaSP'
tags = ["genetics", "windows"]
[menu.main]
  parent = "bio"
+++

DNA Sequence Polymorphism

<http://www.ub.es/dnasp/>

## Installation on Mac

普通のWindowsアプリと同じように WineBottler と Wine
でインストールして使えるが、
インストール後にActiveXの設定をしないと以下のように怒られて
強制終了になってしまう（起動はできる）:

    Error number: 429
    ActiveX component can't create object

### はじめに

1.  [WineBottler]({{< relref "mac/winebottler.md" >}}) および Wine をインストール。
2.  [DnaSP作者のウェブサイト](http://www.ub.es/dnasp/) からインストーラ
    dnasp51001.msi をダウンロード。
3.  微妙に違う2つのやり方があるのでどっちにするか決める。

### DnaSP.appを作る方法

1.  WineBottler を使って `/Applications` に DnaSP.app を作る。
2.  DnaSP.app を一度だけ起動して、すぐ終了する。
    `~/Library/Application Support/Wine/prefixes/DnaSP`
    にPrefixができる。
    WineBottler で“DnaSP”じゃないアプリケーション名を指定して作った人は
    `prefixes/` の下の名前も変わるので注意。
3.  ActiveX の設定をする。
    1.  <http://www.dll-files.com/> らへんから
        `scrrun.dll` を入手し、
        `~/Library/Application\ Support/Wine/prefixes/DnaSP/drive_c/windows/system32/`
        にコピー
    2.  Terminal.app で以下のコマンドを実行:

            cd ~/Library/Application\ Support/Wine/prefixes/DnaSP/drive_c/windows/system32
            cp ../../Program\ Files/DnaSP\ v5/mfc40.dll ./
            WINEPREFIX=${HOME}/Library/Application\ Support/Wine/prefixes/DnaSP /Applications/Wine.app/Contents/MacOS/startwine regsvr32.exe scrrun.dll mfc40.dll threed32.ocx

        wget をインストールしてある人は以下のコマンドで一発のはず:

            cd ~/Library/Application\ Support/Wine/prefixes/DnaSP/drive_c/windows/system32
            wget http://www.dll-download-system.com/dlls/scrrun.zip
            unzip scrrun.zip scrrun.dll
            rm scrrun.zip
            cp ../../Program\ Files/DnaSP\ v5/mfc40.dll ./
            WINEPREFIX=${HOME}/Library/Application\ Support/Wine/prefixes/DnaSP /Applications/Wine.app/Contents/MacOS/startwine regsvr32.exe scrrun.dll mfc40.dll threed32.ocx

4.  インストールされた DnaSP.app を実行

### DnaSP.appを作らず直接dnasp5.exeを実行する方法

1.  Wine でPrefixを `~/.wine` に設定。
    違う場所にした人はこれ以下を適宜読み替えること。
2.  インストーラ dnasp51001.msi を実行し、
    Next 連打で普通にインストールする。
3.  Terminal を起動して以下のコマンドを入力し、
    ActiveX の設定をする。
    コマンドでの作業が嫌いな人はその下の手順でどうぞ:

        cd ~/.wine/drive_c/windows/system32/
        wget http://www.dll-download-system.com/dlls/scrrun.zip
        unzip scrrun.zip scrrun.dll
        rm scrrun.zip
        cp ../../Program\ Files/DnaSP\ v5/mfc40.dll ./
        WINEPREFIX=${HOME}/.wine /Applications/Wine.app/Contents/MacOS/startwine regsvr32.exe scrrun.dll mfc40.dll threed32.ocx

    1.  `~/.wine/driver_c/Program Files/DnaSP v5/` にある
        `mfc40.dll` を
        `~/.wine/drive_c/windows/system32/` にコピー
    2.  <http://www.dll-files.com/> らへんから
        `scrrun.dll` をダウンロードし、
        `~/.wine/drive_c/windows/system32/` にコピー
    3.  右上メニューバーの
        `ワイングラス > DOS Prompt` で出てきた黒い窓に
        以下のコマンドを1行ずつ入力して実行:

            cd windows\system32
            regsvr32 scrrun.dll mfc40.dll threed32.ocx

        DOS窓が出てこなかった人は Terminal から以下のコマンドを実行する:

            cd ~/.wine/drive_c/windows/system32/
            /Applications/Wine.app/Contents/MacOS/startwine regsvr32 scrrun.dll mfc40.dll threed32.ocx

4.  インストールされた dnasp5.exe を実行
