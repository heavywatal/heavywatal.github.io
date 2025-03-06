+++
title = 'zsh'
tags = ["shell"]
[menu.main]
  parent = "dev"
+++

## The most powerful shell

-   <https://www.zsh.org/>
-   <https://zsh.sourceforge.io/>

## Installation

基本的にはOSに入ってる `/bin/zsh` を使う。
新しいのを入れるなら[Homebrew]({{< relref "homebrew.md" >}})を使うのが楽:

```sh
brew install zsh
brew install zsh-completions
```

## Configuration files

<https://zsh.sourceforge.io/Doc/Release/Files.html>

`$ZDOTDIR` 以下の個人設定ファイルが場合に応じて下記の順で読まれる。
システム全体の設定ファイルとして `/etc/z*` が個人設定ファイルの前に読み込まれる。

`/etc/zshenv`
:   スクリプトの実行時も含めてあらゆる場合に読み込まれ、オプションでも外せない。
    Ubuntuはここで基本的な `PATH` を設定。

`.zshenv`
:   スクリプトの実行時も含めてほぼあらゆる場合に読み込まれる。
    インタラクティブ用の設定などはせず、最低限の記述に留める。
    例えば `ZDOTDIR`, `unsetopt NOMATCH` など。
:   ここで `PATH` を設定したい気もするけど、
    OSによっては次の `/etc/zprofile` で上書きされてしまう。
    `unsetopt GLOBAL_RCS` としてそれを防ぐ手もあるけどやや危険。

`/etc/zprofile`
:   Macでは `/usr/libexec/path_helper` が
    `/usr/bin` などの基本的な `PATH` を設定する。
:   ここで `/etc/profile` を `source` して設定するLinuxもある。

`.zprofile`
:   ログインシェルとして立ち上げるときのみ読まれる。
    `export` する環境変数(`PATH` とか)を設定するのに適している。
    `.bash_profile` との共通設定を `.profile` に書いておいて `source` するとか。
:   例えばローカル環境Mac + リモート環境Linux CUIで開発する場合、
    ターミナルも[tmux]({{< relref "tmux.md" >}})もデフォルトでログインシェルを立ち上げるので、
    `.zshrc` に一本化してしまっても構わない。
:   使い分けるのはLinux GUIを使う場合とか、
    よほど重い初期化をログインシェル1回で済ませたい場合とか。

`.zshrc`
:   ログイン・非ログイン問わず、インタラクティブシェルとして立ち上げるときに読まれる。
    こだわりが無ければだいたいどの設定もこれに書いておけば問題ない。
:   `.zprofile` と使い分けるなら
    `setopt` や `autoload` など、親シェルから引き継がれないものはこちら。
    `alias` などは別ファイルを読み込む形にして `.bashrc` と共有。

`.zlogin`
:   `.zshrc` より後に読まれる以外は `.zprofile` と同じ。使わない。

`.zlogout`
:   ログアウト時にしてほしいことが万が一あれば。

以下は Bash の設定ファイル。
<https://www.gnu.org/software/bash/manual/bash.html#Bash-Startup-Files>

`/etc/profile`
:   ログインシェルの場合ここから開始。
:   OSによって `/etc/bash.bashrc` や `/etc/bashrc` を `source` する。

`~/.bash_profile`
:   ログインシェルの場合の2つめ。
    もし存在しなかったら `~/.bash_login` を読みに行き、
    それも存在しなかったら `~/.profile` を読む。
    どれかが読めたところで自動読み込みは終了。
    つまりこれが存在する場合は `~/.profile` が自動では読み込まれないし、
    いずれにせよ `~/.bashrc` を自動で読みには行かない。
:   ほかのシェルとの共通設定を `~/.profile` と `~/.bashrc` に書き、
    ここではそれらを `source` するだけにしておくと見通しがいい。

`~/.bashrc`
:   非ログインでインタラクティブの場合に自動で読まれる唯一のファイル。
:   ログインシェルの場合は自動では読まれないので
    `~/.bash_profile` から `source` するのが普通。
:   非インタラクティブでも[SSH]({{< relref "ssh.md" >}})越しの
    [`rsync`]({{< relref "rsync.md" >}}) とか `scp`
    とかで読まれて失敗の原因となる。
    環境変数やメッセージを `echo` したくなるかもしれないけど、
    標準出力や標準エラーへの書き出しは厳禁。
    始めに `$-` とか `$PS1` を調べて終了するようにしておくと安心。

`$BASH_ENV`
:   スクリプト実行直前にこの環境変数の指すファイルを読み込む。
:   POSIXモードの `sh` として呼ばれた場合は `$ENV` を読みに行く。
:   スクリプトの中で明示的に `source` するほうが明快で安全なので使わない。


### 起動時間短縮

まずはプロファイリングしてボトルネックを知る:

```sh
# head of .zshenv
zmodload zsh/zprof

# tail of .zshrc
zprof
```

`compinit` とかが遅かったり複数回呼ばれていたりするので順番やオプションを変えてみる。


## Glob

### Qualifiers

ファイルの種類を限定したり順番を変えたりできる。
`setopt EXTENDED_GLOB` で有効になる。

- `*(/)`: directories
- `*(.)`: plain files
- `*(@)`: symbolic links
- `*(*)`: executable plain files
- `*(On)`, `*(^on)`: descending order by name.\
  `*(Om)`, `*(^om)`: ascending order by modification time, oldest first.\
  例えば `file` と `file.backup` を比較したいときにただ
  `diff file*` とするとbackupのほうが後に来てしまうのを解決。
  名前順と時間順でデフォルトの方向が逆なのはいずいけど仕方ない。
