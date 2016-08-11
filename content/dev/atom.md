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

いつもの <kbd>cmd ,</kbd> キーで設定画面を起動。

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

<kbd>cmd-alt-i</kbd> でWeb Inspectorを起動させればあらゆる要素を調べることができる。
カーソル位置のスコープを知りたいだけなら <kbd>cmd-alt-p</kbd> が簡便。


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

## Tips

とりあえずコマンドパレット
: <kbd>cmd-shift-p</kbd> で呼び出し、やりたいことを打ち込んでみる

矩形(ブロック)選択
: <kbd>ctrl-shift-down</kbd> / <kbd>ctrl-shift-up</kbd>
: MacではデフォルトでMission Controlに割り当てられてしまっているので
  システム環境設定からそれを解除しておく。

小文字から大文字へ "Editor: Upper Case"
: <kbd>command-k-u</kbd>

大文字から小文字へ "Editor: Lower Case"
: <kbd>command-k-l</kbd>