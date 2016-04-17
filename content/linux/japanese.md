+++
title = 'Linux日本語環境'
tags = ["linux", "writing"]
[menu.main]
  parent = "linux"
+++

なんとなく、日本語Remix版のUbuntuをインストールするのではなく、
英語版Ubuntuに必要な日本語パッケージを導入するという形にしている。

## Google日本語入力

1.  [Japanese packages for testers](https://launchpad.net/~japanese-testers/+archive/ppa)
    のPPAリポジトリを追加してパッケージリストを更新:

        % sudo add-apt-repository ppa:japanese-testers/ppa
        % sudo apt-get update

2.  Synapticあるいはコマンドラインから\`ibus-mozc\`とその依存パッケージをインストール:

        % sudo apt-get -u install ibus-mozc

    もしくは Language Support を起動して
    `Install / Remove Languages... --> Japanese`
    にチェックして Apply。
    `Keyboard input method system --> ibus`

3.  一旦ログアウトなどしてiBusを再起動
4.  Keyboard Input Method (IBus Preferences) を起動して
    `Input Method --> Select --> Japanese --> Mozc`

ここまででとりあえず使える。次に「ことえり」ライクなキーバインドに設定する。

1.  言語バー右端の星の隣の
    `設定ボタン --> Property --> Keymap style --> Custom keymap --> Customize...`
2.  `Edit --> Import predefined mapping --> Kotoeri`
3.  `Edit --> Export to file...` で適当な名前で保存
4.  そのファイルをテキストエディタで開き、以下の2行（タブ区切り）を追加:

        DirectInput Henkan  IMEOn
        Precomposition  MuHenkan    IMEOff

5.  `Edit > Import from file...` で今書き換えたファイルを指定
6.  変換キーで日本語モード、無変換キーで半角英数モードに切り替えられる。
