+++
title = "VSCode"
subtitle = "Microsoft製テキストエディタ"
date = 2022-06-09T10:52:28+09:00
aliases = ["/dev/atom.html"]
tags = ["editor", "writing"]
[menu.main]
  parent = "dev"
+++

Microsoftが開発しているGUIテキストエディタ。
[Atom](https://github.com/atom/atom) を継ぐもの。

<https://code.visualstudio.com>


## フォルダ単位で開く

VSCodeはファイルごとにひとつのウィンドウを開くのではなく、
フォルダ単位でウィンドウを開き、複数のファイルをタブで表示するのが基本。
そうは言ってもホームフォルダまるごととかではなく、
次のような機能を効率的に使えるような、ちょうどいい階層のフォルダを開きたい。

- <kbd>⇧</kbd><kbd>⌘</kbd><kbd>e</kbd>: 左側にツリー表示
- <kbd>⇧</kbd><kbd>⌘</kbd><kbd>f</kbd>: ファイル横断検索
- <kbd>⌘</kbd><kbd>p</kbd>: ファイル名を検索して開く
- 後述の [GitHub Copilot](#github-copilot) に文脈を読んで編集してもらう

例えば[Git]({{< relref "git.md" >}})を使っている人はリポジトリ単位で開く。
[Git Project Manager](https://marketplace.visualstudio.com/items/?itemName=felipecaputo.git-project-manager)
という拡張を入れておくと
<kbd>⌥</kbd><kbd>⌘</kbd><kbd>p</kbd>
から選択できる。


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


### GitHub Copilot

- <https://github.com/features/copilot>
- <https://marketplace.visualstudio.com/items?itemName=GitHub.copilot>

AI coding assistant.
プログラミングはもちろんのこと、論文書きやメール書きも含めテキスト仕事なら何でもサポートしてもらえる。
対話形式もできるけど何より補完機能が便利。
[Git]({{< relref "git.md" >}})を触ったことがない人にもおすすめ。
いろいろな利用法があるけどここに書くのはVSCode上の。

Microsoft製のCopilotだけを使うわけではなく(というか選択肢に無い？)、
ChatGPT や Claude のようないくつかの大規模言語モデルから選択できる。


#### 設定

<https://github.com/features/copilot/plans>\
Free でも使えるけどかなり限定的。
GitHub Education に登録すれば Pro 相当の機能を無料で使わせてもらえる。

1. [GitHub](https://github.com) にアカウントを作る。
   このときはGmailでも何でもいいので普段よく使うアドレスを登録。
1. [GitHub Education](https://github.com/education) にも登録する。
   ここでは教育機関のメールアドレスが必要: `@tohoku.ac.jp`, `@hogwarts.edu`, etc.
   学生証や職員証の写真をアップロードする必要があったかも。
1. VSCode を起動して GitHub Copilot 拡張をインストール。
   <kbd>⇧</kbd><kbd>⌘</kbd><kbd>x</kbd> `copilot`
1. 右上の人型アイコン or 上部中央右の
   <img height=16 width=16 src="https://cdn.simpleicons.org/githubcopilot" alt="Copilotボタン"></a>
   から "Sign in with GitHub to Use GitHub Copilot" みたいなメニューを選択。
   ウェブブラウザに飛ばされるので指示に従って認証。
1. 試しにVSCodeでテキストファイルを新規作成 <kbd>⌘</kbd><kbd>n</kbd>
   して適当に書いて改行<kbd>⏎</kbd>:

   > Dear Professor Makino,

   1秒待つと次のような文がグレーで表示される:

   > I hope this message finds you well. I am writing to express my gratitude for your invaluable support and guidance throughout my research journey. Your mentorship has been instrumental in shaping my academic path, and I am truly thankful for the opportunities you have provided me.

   <kbd>tab</kbd> で採用。
   続けて <kbd>⌘</kbd><kbd>i</kbd> で簡易チャットを開き "Translate it to Japanese" と頼んでみる:

   > 私はこのメッセージがあなたに届くことを願っています。私は研究の旅を通じてあなたの貴重なサポートと指導に感謝の意を表したく、この手紙を書いています。あなたの指導は私の学問的な道を形作る上で非常に重要であり、私に提供してくださった機会に心から感謝しています。


設定確認・変更:

- [GitHub上でのCopilotの設定](https://github.com/settings/copilot)
- VSCode上部
  <img height=16 width=16 src="https://cdn.simpleicons.org/githubcopilot" alt="Copilotボタン"></a>
  右のV字をクリック → "Configure Code Completions...":
  - "Status: Ready (Disabled)" になってたらその下の "Enable Completions" で有効化。
  - "Edit Keyboard Shortcuts..."
  - "Edit Settings..."

#### 使い方

大きく3つの方法がある。
いずれにせよ、開いているファイルや選択中のテキストを文脈として読み取ってもらえる。
AIを使う、というと一般にはチャットのイメージが強いだろうけど、
質問や命令を考える必要さえない自動補完こそがズボラな人間の強力な味方。

1. チャット用のサイドバーで質問や編集を行う。
   [3つのモードがある](https://github.blog/ai-and-ml/github-copilot/copilot-ask-edit-and-agent-modes-what-they-do-and-when-to-use-them/):
   - **Ask**: 質問に答えてもらうだけ。\
     e.g., "How do you remove axis ticks in ggplot?", "Explain this code".
   - **Edit**: 編集を提案してもらう。ファイル書き換えユーザーが差分を確認してから。\
     e.g., "Correct grammatical and typographical errors", "Translate it to Japanese".
   - **Agent**: 目的だけを伝えて、ファイルの編集もいちいち確認せず任せる。
2. <kbd>⌘</kbd><kbd>i</kbd> で **Inline Chat** を開き、
   上記 "Ask" や "Edit" に相当することをその場で行う。
3. とにかく自動補完をオンにしてテキストを書き、続きを **Inline Suggestion** してもらう。


key  | command | description
---- | ------- | -----------
<kbd>esc</kbd> | Hide Inline Suggestion | 却下
<kbd>tab</kbd> | Accept Inline Suggestion | 全部採用
<kbd>⌘</kbd><kbd>→</kbd> | Accept Next Word Of ... | 部分採用
<kbd>⌥</kbd><kbd>&bsol;</kbd> | Trigger Inline Suggestion | 自動補完されないときに手動で
<kbd>⌘</kbd><kbd>i</kbd> | Toggle Inline Chat | その場で簡易チャット
<kbd>^</kbd><kbd>⌘</kbd><kbd>i</kbd> | Toggle Chat | 右側にバーを開いてチャット
<kbd>⌘</kbd><kbd>.</kbd> | Set Chat Mode | Ask → Edit → Agent → ...

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
   1. Integrated Terminalから手動でRを立ち上げるような場合、
      設定ファイルを明示的に読み込む必要があるので `.Rprofile` に追記:
      ```r
      if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
        source("~/.vscode-R/init.R")
      }
      ```
      [tmux]({{< relref "tmux.md" >}}) が `TERM_PROGRAM` を上書きすることに注意。

ここまでやれば `View()` や `help()` などをVSCodeで表示できる。
以下はお好みで。

- [VSCodeの設定](https://github.com/REditorSupport/vscode-R/wiki/Extension-settings)
  - `"r.rmarkdown.enableCodeLens": false`: CodeLensは邪魔。
  - `r.session.viewers.viewColumn`:
    [interactive viewers](https://github.com/REditorSupport/vscode-R/wiki/Interactive-viewers)
    をどこで表示するか。
    - `Active`: 既存editor group内の新規タブ。
    - `Two`, `Beside`: 新規editor groupを作って表示。2つの違いは不明。
      左右ではなく上下に分割して表示させるには
      `"workbench.editor.openSideBySideDirection": "down"` 。
      [See vscode-R#1126](https://github.com/REditorSupport/vscode-R/issues/1126).
    - `Disable`: ネイティブに任せる。
    - デフォルト: {plot: Two, browser: Active, viewer: Two, pageViewer: Active, view: Two, helpPanel: Two}
- [Rの設定](https://github.com/REditorSupport/vscode-R/wiki/R-options)。
  ほとんどの項目はVSCode側で設定できる。
  それらをどうしてもRから上書きしたいときに使う。
- 図をどうやって描画するか:
  - デフォルト: `png()` で書き出してVSCode内に表示。
    [ragg](https://github.com/r-lib/ragg/)を使うオプションが欲しいところだけど無い。
    描画領域のサイズ変化に追従する機能なども無い。
  - [httpgd](https://nx10.github.io/httpgd/) を使ってVSCode内に表示。
    `install.packages("httpgd")` した上で
    VSCodeの設定に `"r.plot.useHttpgd": true` を追加。
    SVG形式なので[要素数が多くなるほど急激に重くなり](https://github.com/nx10/httpgd/issues/113)、例えばdiamondsの散布図でも使い物にならない。
    PNGも選べるように[unigd](https://nx10.github.io/unigd/)を開発中らしい。
  - 外部ターミナルでRを起動するのと同じように独立XQuartzなどで表示。
    上記 `viewColumn` を `"plot": "Disable"` に設定する。
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
