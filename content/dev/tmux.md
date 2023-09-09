+++
title = 'tmux'
subtitle = "仮想端末でリモート仕事を安全に"
tags = ["job", "shell"]
[menu.main]
  parent = "dev"
+++

<https://tmux.github.io/>

[GNU screen](https://www.gnu.org/software/screen/)
の後を継ぐ端末多重化ソフト(terminal multiplexer)。

1つの画面の中でウインドウを追加・分割して複数の端末を開く
:   GUIアプリのタブ代わりに。
:   1つのsshセッションで複数の端末を持てる。

ssh切断後も端末丸ごと継続され、後でまた繋ぎ直せる
:   不意のssh切断でも作業が失われない。
:   別の端末から接続しても同じ作業を継続できる。
:   [`nohup`]({{< relref "nohup.md" >}}) とかバックグラウンド化とか考えるより楽チン。

[Homebrew]({{< relref "homebrew.md" >}}) で一発インストール:
`brew install tmux`

## キーバインド

tmux 内で **prefix key** に続けて特定のキーを送信すると、
そのキーに応じたさまざまなコマンドを実行できる
(e.g. <kbd>prefix</kbd><kbd>?</kbd> でキーバインドを列挙)。
prefix keyはデフォルトで <kbd>C-b</kbd>
(<kbd>control</kbd><kbd>b</kbd>の略記、<kbd>^b</kbd>と等価)
だがそれはキャレット左移動に使われるべきなので後述のように変更する。

key          | command | description
------------ | ------- | -----------
<kbd>?</kbd> | `list-keys` |
<kbd>:</kbd> | `command-prompt` |
<kbd>d</kbd> | `detach-client` |
<kbd>c</kbd> | `new-window` |
<kbd>n</kbd> | `next-window` |
<kbd>p</kbd> | `previous-window` |
<kbd>l</kbd> | `last-window` |
<kbd>,</kbd> | `rename-window` |
<kbd>.</kbd> | `move-window` |
<kbd>0 1 2 3</kbd> | `select-window -t :=N` |
<kbd>"</kbd> | `split-window` | 横長・縦並びに分割: 日
<kbd>%</kbd> | `split-window -h` | 縦長・横並びに分割: Φ
<kbd>;</kbd> | `last-pane` | 直前のペイン(往復)
<kbd>o</kbd> | `select-pane -t:.+` | 番号順にペインを巡回
<kbd>↑</kbd><kbd>↓</kbd><kbd>←</kbd><kbd>→</kbd> | `select-pane -U` |
<kbd>C-↑</kbd><kbd>C-↓</kbd> | `resize-pane -U` | ペインサイズ変更
<kbd>C-o</kbd> | `rotate-window` | レイアウトを維持してペインを回す
<kbd>space</kbd> | `next-layout` | レイアウトを変更する
<kbd>!</kbd> | `break-pane` | ペインを独立したウィンドウにする
<kbd>[</kbd> | `copy-mode` |
<kbd>]</kbd> | `paste-buffer` |

### コピーモード

上に戻ってスクロールしたり、その内容をコピーしたいときはコピーモードを使う。
キーボードから <kbd>prefix</kbd><kbd>[</kbd> で入れるほか、
`set -g mouse on` を設定すれば上スクロールで自然に入れる。

コピーモードでのキー操作はデフォルトだとemacs風で、
環境変数 `EDITOR`/`VISUAL` やオプション `mode-keys` からviに変更できる。
シェルを介さずに直接起動する場合も考えると明示的にオプション設定しておくのが無難。

command           | vi    | emacs | description
----------------- | ----- | ----- | -----------
`cancel` | <kbd>q</kbd> | <kbd>esc</kbd> | コピーモード終了
`begin-selection` | <kbd>space</kbd> | <kbd>C-space</kbd> | 選択開始点をマーク
`copy-pipe-and-cancel` | <kbd>enter</kbd> | <kbd>M-w</kbd> | 選択範囲内を `copy-command` に送って終了

コピーと同時に終了せずモードや選択状態を維持したい場合はキーバインドを
`copy-pipe` や `copy-pipe-no-clear` に変更する。

`copy-pipe*` の宛先はデフォルトでtmux内のバッファになっており、
ペーストはtmux内で <kbd>prefix</kbd><kbd>]</kbd> するしかない。
macOSで次のように設定しておけば、
<kbd>⌘command</kbd><kbd>c</kbd>と同じところにコピーして、
アプリを超えて<kbd>⌘command</kbd><kbd>v</kbd>できるようになる:
```
if "command -v pbcopy" "set -s copy-command pbcopy"
```


## 設定

設定ファイル： `~/.tmux.conf`

<https://github.com/heavywatal/dotfiles/blob/master/.tmux.conf>

`prefix <key>`
:   使えるのは `^h` や `^[` のようなASCIIキャレット記法が存在するもの。
    シェルやエディタであまり使わず左手だけで完結できるキーがいい。
    tmux の頭文字で覚えやすい <kbd>C-t</kbd> がよかったけど
    [`fzf`](https://github.com/junegunn/fzf) と衝突。
:   同じキーを `bind <key> send-prefix` に設定しておけば、
    2回押しのうち1回分がtmuxを貫通して伝わる。
    使う頻度の低いキーとの衝突ならこれで乗り切れる。

`aggressive-resize [on | off]`
:   サイズの異なる端末からattachしたときにウィンドウサイズを変更する。

`update-environment <variables>`
:   attachするときに環境変数を親プロセスから持ち込んで既存sessionの値を上書きする。
:   `DISPLAY` や `SSH_AUTH_SOCK` などがデフォルトで含まれているので、
    `-a` オプションで追加するのが無難。
:   `TERM_PROGRAM` は特殊で、指定しても強制的に `tmux` に上書きされる。
    `showenv TERM_PROGRAM` で上書き前の情報がとれるようにはなるので、
    それを使って環境変数を再上書きすることは可能。
    ただしattachの度にそこまで実行できないことには注意。
    [See tmux#3468](https://github.com/tmux/tmux/issues/3468).

デタッチ後しばらくしてシェルを起動すると残存セッションを忘れがちなので、
以下のようなものを `.zshrc` とかに書いておけば表示で気付ける。

```sh
if [ -n "$TMUX" ]; then
  eval $(tmux showenv TERM_PROGRAM)
else
  tmux has-session >/dev/null 2>&1 && tmux list-sessions
fi
```

`open` や `pbcopy` などがうまく働かなくて
[reattach-to-user-namespace](https://github.com/ChrisJohnsen/tmux-MacOSX-pasteboard)
が必要になる問題は既に解消された。


## 利用例

1.  リモートサーバーに ssh ログインし、
    tmux の新しいセッションを開始:
    ```sh
    ssh remote.sample.com
    tmux
    ```
1.  ウィンドウを左右に分割し <kbd>prefix</kbd><kbd>%</kbd>、
    右ペインでPythonインタプリタを起動 `python`:
1.  左ペインにフォーカスを戻し<kbd>prefix</kbd><kbd>o</kbd>、
    ファイルを閲覧したり何だり `less ~/.tmux.conf`
1.  新しいウィンドウを作って <kbd>prefix</kbd><kbd>c</kbd>、
    `root` 仕事をしたり何だり `su -`
1.  ウィンドウを切り替える <kbd>prefix</kbd><kbd>l</kbd>, <kbd>prefix</kbd><kbd>n</kbd>, <kbd>prefix</kbd><kbd>p</kbd>
1.  このセッションをデタッチし <kbd>prefix</kbd><kbd>d</kbd>、
    ログアウトして家に帰る `exit`
1.  家からサーバーに再び ssh ログインして、
    さっきの tmux セッションをアタッチして作業を再開:
    ```sh
    ssh remote.sample.com
    tmux attach -d
    ```
1.  セッション内のすべてのペインで `exit` して終了。
