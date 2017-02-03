+++
date = "2016-04-25T17:03:29+09:00"
tags = ["editor", "writing"]
title = "Atom"
subtitle = "最強のテキストエディタ"

[menu.main]
  parent = "dev"
+++

Githubが開発したオープンソースのGUIテキストエディタ。
ChromiumとNode.js(を用いたElectronフレームワーク)でできており、
どのOSでも同じように動作する。
コミュニティの力により現在も急速に成長中。
当然[Git]({{< relref "git.md" >}})によるバージョン管理とも相性が良い。

- https://atom.io
- https://flight-manual.atom.io/

```sh
% brew cask install atom
```

## Tips

とりあえずコマンドパレット
: <kbd>command-shift-p</kbd> で呼び出し、やりたいことを打ち込んでみる

矩形(ブロック)選択
: <kbd>ctrl-shift-down</kbd> / <kbd>ctrl-shift-up</kbd>
: MacではデフォルトでMission Controlに割り当てられてしまっているので
  システム環境設定からそれを解除しておく。

選択範囲を掴んで移動
: <kbd>ctrl-command-***</kbd>

コメントアウト、解除
: <kbd>command /</kbd> `editor:toggle`

閉じタグを挿入
: <kbd>command-alt .</kbd> `bracket-matcher:close-tag`

小文字から大文字へ
: <kbd>command-k-u</kbd> `editor:upper-case`

大文字から小文字へ
: <kbd>command-k-l</kbd> `editor:lower-case`

プロジェクト内検索
: <kbd>command-shift-f</kbd> `project-find:show`

### Tree view

key  | command
---- | ----
<kbd>ctrl-0</kbd> | tree-view:toggle-focus
<kbd>m</kbd> | tree-view:move
<kbd>d</kbd> | tree-view:duplicate
<kbd>a</kbd> | tree-view:add-file
<kbd>shift-a</kbd> | tree-view:add-folder

矢印キーはそのものでもEmacs/Vim系でも想像通りの挙動

プロジェクト内のファイルを開きたいだけなら
<kbd>command-t</kbd> でインクリメントサーチする癖をつけるほうが早い。


## 環境設定

- http://flight-manual.atom.io/using-atom/sections/basic-customization/
- https://github.com/heavywatal/dotfiles/tree/master/.atom

いつもの <kbd>command ,</kbd> キーで設定画面を起動。

設定ファイルは `~/.atom/` 以下に置かれる。
設定画面から"Open Config Folder"ボタンを押すとAtom内でそれらを開くことができる。
変更は即時反映される。

`config.cson`
: Core Settings

`keymap.cson`
: Keybindings
  デバッグしたいときは <kbd>command .</kbd> でKey Binding Resolverを起動するとよい。

`snippets.cson`
: 定型句に名前を付けておいて簡単に呼び出せるようにする。
  http://flight-manual.atom.io/using-atom/sections/snippets/

`styles.less`
: エディタ本体も含めていろんな部分をCSS的にスタイル設定可能。

<kbd>command-alt-i</kbd> でWeb Inspectorを起動させればあらゆる要素を調べることができる。
カーソル位置のスコープを知りたいだけなら <kbd>command-alt-p</kbd> が簡便。


## パッケージ

- https://atom.io/packages
- https://atom.io/users/heavywatal/stars

環境設定のInstallメニューからインストールし、
Packagesメニューで管理する。

`apm` コマンドを利用してもよい。
```sh
apm install pigments
apm uninstall pigments
apm upgrade
```
