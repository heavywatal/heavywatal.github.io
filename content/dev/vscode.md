+++
title = "VSCode"
subtitle = "Microsoft製テキストエディタ"
date = 2022-06-09T10:52:28+09:00
tags = ["editor", "writing"]
[menu.main]
  parent = "dev"
+++

Microsoftが開発しているGUIテキストエディタ。
[Atom]({{< relref "atom.md" >}}) を継ぐもの。

<https://code.visualstudio.com>


## Keyboard Shortcuts

<https://code.visualstudio.com/docs/getstarted/keybindings>

key  | command | description
---- | ------- | -----------
<kbd>⇧</kbd><kbd>⌘</kbd><kbd>p</kbd> | **Command Palette** | これさえ覚えれば
<kbd>⌥</kbd><kbd>⌘</kbd><kbd>p</kbd> | Open Git Projects | [GPM](https://marketplace.visualstudio.com/items?itemName=felipecaputo.git-project-manager)を使ってディレクトリを開く
<kbd>⌘</kbd><kbd>p</kbd> | Search files by name | プロジェクト内の別ファイルを開く
<kbd>⇧</kbd><kbd>⌘</kbd><kbd>f</kbd> | Find in Files | プロジェクト内をファイル横断で検索
<kbd>⇧</kbd><kbd>⌘</kbd><kbd>e</kbd> | Show Explorer | プロジェクト内のファイル一覧
<kbd>⇧</kbd><kbd>⌘</kbd><kbd>x</kbd> | Show Extensions | 拡張一覧。インストールもここから
<kbd>⌘</kbd><kbd>/</kbd> | Toggle Line Comment | コメントアウトしたり外したり
<kbd>^</kbd><kbd>`</kbd> | Toggle Terminal
<kbd>^</kbd><kbd>⏎</kbd> | Run Selected Text | エディタからターミナルにテキストを送る(要設定)
<kbd>⌥</kbd><kbd>↓</kbd> | Move line down
<kbd>⇧</kbd><kbd>⌥</kbd><kbd>↓</kbd> | Copy line down
<kbd>⌥</kbd><kbd>⌘</kbd><kbd>↓</kbd> | `cursorColumnSelectDown` | 矩形(ブロック)選択
<kbd>⌘</kbd><kbd>k</kbd><kbd>⌘</kbd><kbd>s</kbd> | Keyboard Shortcuts | 一覧
<kbd>⌘</kbd><kbd>k</kbd><kbd>⌘</kbd><kbd>r</kbd> | Keyboard Shortcuts Ref. | 抜粋PDF

設定ファイルは
[`~/.config/Code/User/keybindings.json`](https://github.com/heavywatal/dotfiles/blob/master/.config/Code/User/keybindings.json)

```json
    {
        "key": "ctrl+enter",
        "command": "workbench.action.terminal.runSelectedText",
        "when": "editorTextFocus && resourceExtname != .ipynb"
    },
```


## 環境設定

<https://code.visualstudio.com/docs/getstarted/settings>

いつもの <kbd>⌘</kbd><kbd>,</kbd> キーで設定画面を起動。

設定ファイルは
[`~/.config/Code/User/settings.json`](https://github.com/heavywatal/dotfiles/blob/master/.config/Code/User/settings.json)


## 拡張パッケージ

<https://code.visualstudio.com/docs/editor/extension-marketplace>


### R

<https://code.visualstudio.com/docs/languages/r>

人にはRStudioを薦めるけど個人的にはこっちのほうが使いやすい。

1. [languageserver](https://github.com/REditorSupport/languageserver/)
   をインストールする:
   ```r
   install.packages("languageserver")
   ```
1. [R extension for VSCode](https://marketplace.visualstudio.com/items?itemName=REditorSupport.r)
   をインストール。
1. [R session watcher](https://github.com/REditorSupport/vscode-R/wiki/R-Session-watcher)
   の設定。
   1. VSCodeの設定に `"r.sessionWatcher": true` を追加して再起動。
      "Create R terminal" コマンドで専用コンソールを立ち上げる場合はこれだけ。
   1. ターミナルの[tmux]({{< relref "tmux.md" >}})からRを立ち上げるような場合、
      設定ファイルを明示的に読み込む必要があるので `.Rprofile` に追記:
      ```r
      if (interactive() && Sys.getenv("RSTUDIO") == "") {
        Sys.setenv(TERM_PROGRAM = "vscode")
        source("~/.vscode-R/init.R")
      }
      ```

ここまでやれば `View()` や `help()` などをVSCodeで表示できる。
以下はお好みで。

- [VSCodeの設定](https://github.com/REditorSupport/vscode-R/wiki/Extension-settings)
- [Rの設定](https://github.com/REditorSupport/vscode-R/wiki/R-options)。
  例えば、作図や`View()`の度に Editor Layout "Two Columns" で分割表示されるのが嫌なので:
  ```r
  options(vsc.plot = FALSE)
  options(vsc.view = "Active")
  options(vsc.helpPanel = "Active")
  ```
  ["Two Rows" に固定する方法があればいいんだけど...](https://github.com/REditorSupport/vscode-R/issues/1126)
- VSCode内で図を描画したければ
  [httpgd](https://nx10.github.io/httpgd/) を入れる:
  ```r
  install.packages("httpgd")
  ```
  VSCodeの設定に `"r.plot.useHttpgd": true` を追加。
  <br>ただしこれも前述のレイアウト問題の影響を受ける。
- [IPython]({{< relref "ipython.md" >}})
  のようにコンソール上での補完や色付けなどを強化したければ
  [radian](https://github.com/randy3k/radian)
  を入れる。
  [設定](https://github.com/randy3k/radian#settings)
  はRコードの形で `~/.config/radian/profile` に書く。


### Python

<https://code.visualstudio.com/docs/languages/python>

Microsoft公式パッケージが提供されていて、がっちりサポートされている感じ。

[Python](https://marketplace.visualstudio.com/items?itemName=ms-python.python)
: とりあえずこれを入れれば下記の依存パッケージも自動で入るはず。

[Pylance](https://marketplace.visualstudio.com/items?itemName=ms-python.vscode-pylance)
: 補完や型チェックなどを担う language server.
: 中では [Pyright](https://github.com/microsoft/pyright) が使われている。
: 警告を非表示にするには次のように書く:
  ```py
  # pyright: reportUnusedVariable=false
  ```
  設定項目: <https://github.com/microsoft/pyright/blob/main/docs/configuration.md>

[Jupyter](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter)
: Jupyter notebook support


## Glossary

IntelliSense
: Microsoft製の補完システム。
  各種 language server から情報を受け取る。

CodeLens
: 補助的な情報を表示する機能。
: ファイルの内容と区別しにくい形で挿入されてすごく邪魔なので確実に切る。

[Integrated Terminal](https://code.visualstudio.com/docs/editor/integrated-terminal)
: VSCode内で開けるターミナル。<kbd>^</kbd><kbd>`</kbd>
: シェルを介さず直接 [tmux]({{< relref "tmux.md" >}}) を起動する設定も可能。
