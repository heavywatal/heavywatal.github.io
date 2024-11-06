+++
title = 'プロセス管理'
subtitle = "nohup, disown, kill"
tags = ["job"]
[menu.main]
  parent = "dev"
+++

## 実行中プロセスを知る

`/Applications/Utilities/Activity Monitor.app` を見るのが楽チン。

### `ps`

OSによって挙動がかなり異なるのでややこしいが、
LinuxのほうがいくらかBSDに歩み寄ってくれているらしい。
とりあえず全プロセスを表示するオプション `ax`
とカラムをたくさん表示するオプション `u`
をつけて `less` や `grep` に流すのが基本。

ソートの仕方は共通ではない。
出力形式の調整は共通の `o` オプション。

```sh
ps aux | less

## Sort on Mac
ps auxm | head  # by Memory
ps auxr | head  # by CPU

## Sort on Linux
ps auxk -pmem | head  # by Memory
ps auxk -pcpu | head  # by CPU

## Output format
ps axo user,pid,pcpu,pmem,command
```

### `top`

一度表示して終わりではなく、<kbd>q</kbd> で閉じるまで一定間隔で更新される。
起動してからソート項目を切り替えられる。
プラスで昇順、マイナスで降順。

by     | Linux | Mac
------ | ----- | ----
PID    | <kbd>N</kbd> | <kbd>o -pid</kbd>
CPU    | <kbd>P</kbd> | <kbd>o -cpu</kbd>
Memory | <kbd>M</kbd> | <kbd>o -mem</kbd>
Time   | <kbd>T</kbd> | <kbd>o -time</kbd>


### `jobs`

システム全体のプロセスが見える上記コマンドとは違い、
これで見えるのはそのシェルから実行したジョブだけ。
末尾に`&`をつけてバックグラウンドで走らせたジョブや、
<kbd>ctrl-z</kbd> でsuspendしたジョブを眺めるのに使う（下記）。


## ジョブコントロール

何らかのプログラムを実行:

    top

ここで <kbd>ctrl-z</kbd> を押すと、プロセスはバックグラウンドで一時停止する:

    [1]  + 19310 suspended  top

バックグラウンドのプロセスを確認するには `jobs` コマンド。
左からジョブ番号、カレントジョブか否か、状態、コマンド内容:

    jobs
    [1]  + suspended  top

フォアグラウンドで再開するには `fg` コマンド。
引数としてジョブ番号かプログラム名をパーセントに続けて指定する
(`%1` とか `%a.out` とか)。
[zsh]({{< relref "zsh.md" >}}) なら補完もしてくれる。
引数を省略すると、`jobs` で `+` がついてるカレントジョブが選択される。
:

    fg %1
    [1]  + 19310 continued  top

再び一時停止して、今度はバックグラウンドで再開する。コマンドは `bg`:

    [1]  + 19310 suspended  top
    bg %top

始めからバックグラウンドで走らせるなら末尾にアンド:

    top &

中止させるには `kill` コマンドで `SIGTERM` シグナルを送る。
指定するのはジョブ番号でも生のプロセスIDでもよい。
それでも応答しないやつを強制終了するには
`-9` オプションをつけて `SIGKILL` を送るのが最後の手段。
`killall` は指定した文字列と前方一致するプロセスを
すべて `kill` するショートカット。

    kill %1
    kill -9 %a.out
    kill 19310
    killall a.out


## ログアウト後も継続

バックグラウンドで実行中のプロセスも、
ログアウトするとhangup(`HUP`)シグナルによって終了してしまう。
これを無視して実行し続けるようにプログラムを起動するのが `nohup`。
バックグランド化までは面倒見てくれないので末尾の `&` を忘れずに:

    nohup COMMAND [ARGUMENTS] &

標準出力と標準エラー出力は指定しなければ `nohup.out`
または `~/nohup.out` に追記モードでリダイレクトされる。

- 標準出力の書き出し先を指定するには `>{OUTFILE}`
- 標準エラー出力の書き出し先を指定するには `2>{OUTFILE}`
- 標準エラー出力を標準出力と同じところに流すには `>{OUTFILE} 2>&1`

```sh
nohup COMMAND >out.log 2>err.log &
```

[ssh]({{< relref "ssh.md" >}}) 接続先のサーバーで `nohup` ジョブを走らせるときは
標準入出力をすべて切っておかないと期待通りにsshを抜けられない:

```sh
nohup COMMAND >/dev/null 2>&1 </dev/null &
```

うっかり普通に開始してしまったプロセスを後から `nohup` 状態にするには、
一時停止して、バックグラウンドで再開して、`disown` に渡す:

    ./a.out  # control + z
    [1]  + 19310 suspended  a.out
    bg %1
    [1]  + 19310 continued  a.out
    disown %1

See [tmux]({{< relref "tmux.md" >}}).\
`nohup`, `disown` がプロセス単位で切り離すのに対して、
`tmux` は端末セッション丸ごと切り離し＆復帰することができる。
