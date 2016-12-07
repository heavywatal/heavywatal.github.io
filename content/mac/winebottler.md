+++
title = 'WineBottler'
tags = ["mac", "windows"]
[menu.main]
  parent = "mac"
+++

## run Windows apps on OS X

-   Wine: [<http://www.winehq.org/>](http://www.winehq.org/)
-   WineBottler: [<http://winebottler.kronenberg.org/>](http://winebottler.kronenberg.org/)

WineはLinuxなどのUNIX系OS上でWindowsアプリを実行できるようにするプログラム群。
WineBottlerはMac OS XでWineを使うための環境構築を手伝ってくれる。

## Note

-   仮想的なC:ドライブなどが含まれるひとつのWindows環境をPrefixと呼ぶ。
-   実行するWindowsアプリがPrefix内にある必要は無いが、
    実行する環境としてひとつのPrefixを指定する必要がある。
-   Windowsアプリ.exeから単独のMacアプリ.appを作る方法と、
    PrefixをひとつだけMacに作ってWindowsアプリをそこで実行する方法の2通りの使い方ができる。
    -   前者の場合、アプリそれぞれに対してPrefixが作られるため
        Mac上に複数のWindows環境が作られることになるが、
        いかにもMacネイティブっぽいアプリケーションができあがるので見栄えが良い。
    -   後者の場合、一度Prefixを作ってしまえばあとは何も意識せず
        その環境で直接Windowsアプリを実行したりインストールしたりできる。
    -   両者を使い分けても良い。
        とりあえずWindowsアプリを実行できる環境としてPrefixをひとつ作っておくが、
        よく使うアプリは\*.appに固める、とか。
-   Wine や WineBottler を使うときはOSの言語設定を英語にする
    (要、一旦ログアウトして再ログイン)。
    基本的にはもう日本語には戻さないつもりで。（と思ったけど、
    `Wineメニュー --> Configuration --> Desktop Integration --> Appearance`
    からフォントを適切に設定すれば日本語のままでも大丈夫かも。知らんけど）

## Installation

1.  [作者のウェブサイト](http://winebottler.kronenberg.org/) から
    `WineBottlerCombo_*.dmg` をダウンロード
2.  それを開いて、Wine.app と WineBottler.app を
    `/Applications/` 以下にコピー

## make Mac Application.app from Windows Application.exe

Windowsアプリ.exe またはその インストーラ.msi から
単独の Macアプリ.app を作る方法

1.  WineBottler.app を起動
2.  左サイドバーから `Create Custom Prefixes` を選択
    1.  Install File [\_\_\_\_] に元となるWindowsアプリまたはそのインストーラを入力。
        (Select File...) から選んでもよい。
    2.  それがインストーラではなく単体で動くアプリケーションの場合は
        その下の [\_] Copy only. にチェック
    3.  prefix template は new prefix のままで
    4.  Winetricks はここでは無視。ランタイム環境が必要な場合はここで指定しとくのかな？
    5.  Self-contained [\_] にチェックすると、
        Wine がインストールされてないMacでも動くアプリケーションを作れる。
        が、200MBくらい大きくなるので、自分のMacで使うだけならチェックしないほうがいい。
    6.  その下の3つはとりあえず無視。
    7.  Runtime Version [\_\_\_\_] に元のWindowsアプリのバージョンを入力。
        分からなきゃそのまま1.0で。
    8.  Install ボタンを押す。

3.  アプリケーションの名前と置く場所を聞かれるので、
    `/Applications` などの適当な場所に `適当な名前.app` を指定してSave
4.  Windowsっぽいインストーラが立ち上がったら基本的に全部Yesで進む。
5.  Select Startfile と言われたら、
    `Program Files/` 以下にインストールされた目的のアプリを指定してOK
6.  できあがったアプリケーションをダブルクリックして実行

この方法で作ったアプリケーションを一度起動すると、
`~/Library/Application\ Support/Wine/prefixes/`
以下にアプリケーションと同名のPrefixが作られ、
以後そのアプリケーションはそのPrefix環境で実行される。
(e.g. Notepad.app を立ち上げると
`~/Library/Application\ Support/Wine/prefixes/Notepad`
というPrefixが作られ、その環境で実行される)

## execute Application.exe on One Prefix

実行環境としてPrefixをひとつ作っておき、`Windowsアプリ.exe` を直接実行する方法。
Prefixは複数作れるが、ひとつだけ `~/.wine` というPrefixを作って
常にそれを使うのが分かりやすい。（LinuxにおけるWineの動作）

### Prefix

1.  Wine.app を起動
2.  画面右上にある `ワイングラスのメニュー --> Change Prefix...`
    -   Add... ホーム以下の適当な場所に適当な名前で作る
        (ホーム直下に `.wine` という名前で作ることを推奨)
    -   間違って作ったPrefixなどはマイナスアイコンで消せる
    -   虫眼鏡アイコンをクリックするとそのPrefixがどこにできたかFinderで確認できる
    -   設定したいPrefixを選択してOK

3.  右上ワイングラスのメニューで current prefix がちゃんと設定されてるか確認

### Execution

1.  上記の手順で目的のPrefixに設定されてることを確認。
2.  FinderなどからWindowsアプリ
    (notepad.exe とか setup.msi とか) を実行
    -   普通はダブルクリックで行ける
    -   ダメなら右クリックから Wine.app を指定

3.  Wineダイアログが出てきたら、
    Run directly in [...] でPrefixを指定してGo
    -   Don’t show this dialog again
        にチェックすれば次回からこのステップを飛ばせる

4.  そのアプリがインストーラなら、手順に従ってインストールする
    -   全部Yesで進むと大概 `C:\\Program Files\\` 以下にインストールすることになるが、
        それらはMac上では `<Prefix>/drive_c/Program Files` に相当する。

### misc

-   OSの言語設定を日本語にした状態でWineを使うと文字化けしたりPrefix内の設定が壊れたりするらしい。
-   X11 (または XQuartz) というウィンドウシステムを使って
    Windowsアプリを表示するので、Wineが動いてるうちはそれらが起動した状態になる。
    （何か変なの起動したとか思わないように）
-   PrefixはLinuxに倣って `~/.wine` とするのがおすすめだけど、
    これだと不可視ファイルになってしまうので混乱する人も出てきそう。
    Finderで開くには2つの方法がある
    -   右上の `ワイングラス --> Change Prefix... --> 右端の虫眼鏡アイコン`
    -   Terminalを開いて `open ~/.wine` を実行
-   インストールされたWindowsアプリを立ち上げる度にPrefixを開いて
    `Program Files` まで行くというのは面倒なので、
    `/Applications/` 以下にエイリアスを作るなどする。
    Finder で `command + alt` + drag-and-drop するか、
    コマンドで以下のようにする。:

        % ln -s ~/.wine/drive_c/Program\ Files /Applications/

-   wine コマンドは
    `/Applications/Wine.app/Contents/Resources/bin/wine` にあるが、
    パスなどをうまく設定するため
    `/Applications/Wine.app/Contents/MacOS/startwine`
    というスクリプトから動かしたほうがいい。
    `.zshrc` などに以下のようなエイリアスを記述するとよいだろう:

        if [ $(uname) = Darwin ]; then
              alias wine=/Applications/Wine.app/Contents/MacOS/startwine
        fi

    これでWindowsアプリをコマンドラインから簡単に起動できる。
    ちなみに `<Prefix>/drive_c/windows/system32/` 以下にある
    実行ファイルにはパスが通っており、拡張子.exeも省略可能:

        % wine notepad

-   シェル変数 `WINEPREFIX` によって
    wine コマンドが実行されるPrefixを指定できる。
    相対パスではなく絶対パスで。デフォルトは `~/.wine`:

        % WINEPREFIX=${HOME}/Wine wine notepad
