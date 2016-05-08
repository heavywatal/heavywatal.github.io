+++
date = "2016-04-25T17:03:29+09:00"
tags = ["editor", "writing"]
title = "Atom"
subtitle = "最強のテキストエディタ"

[menu.main]
  parent = "dev"
+++

https://atom.io

http://flight-manual.atom.io/

Githubが開発したオープンソースのGUIテキストエディタ。
ChromiumとNode.js(を用いたElectronフレームワーク)でできている。
コミュニティの力により現在も急速に成長中。
当然[Git]({{< relref "git.md" >}})によるバージョン管理とも相性が良い。

```sh
brew cask install atom
```

## 環境設定

http://flight-manual.atom.io/using-atom/sections/basic-customization/

https://github.com/heavywatal/dotfiles/tree/master/.atom

いつもの `cmd-,` キーで設定画面を起動。

設定ファイルは `~/.atom/` 以下に置かれる。
設定画面から"Open Config Folder"ボタンを押すとAtom内でそれらを開くことができる。
変更は即時反映される。

`config.cson`
: Core Settings

`keymap.cson`
: Keybindings

`snippets.cson`
: 定型句に名前を付けておいて簡単に呼び出せるようにする。
  http://flight-manual.atom.io/using-atom/sections/snippets/

`styles.less`
: エディタ本体も含めていろんな部分をCSS的にスタイル設定可能。

`cmd-alt-i` でWeb Inspectorを起動させればあらゆる要素を調べることができる。
カーソル位置のスコープを知りたいだけなら `cmd-alt-p` が簡便。


## パッケージ

https://atom.io/packages

https://atom.io/users/heavywatal/stars

環境設定のInstallメニューからインストールし、
Packagesメニューで管理する。

`apm` コマンドを利用してもよい。
```sh
apm install pigments
apm uninstall pigments
apm upgrade
```

## 使い方

とりあえず `cmd-shift-p` でコマンドパレットを呼び出してみる。

矩形(ブロック)選択
: `ctrl-shift-down` / `ctrl-shift-up`
: MacではデフォルトでMission Controlに割り当てられてしまっているので
  システム環境設定からそれを解除しておく。
