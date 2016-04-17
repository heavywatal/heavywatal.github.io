+++
title = 'プロセス管理'
subtitle = "nohup, disown, kill"
[menu.main]
  parent = "dev"
+++

## 実行中プロセスを知る

`top`

`ps`

`jobs`

Activity Monitor

## ジョブコントロール

プログラムを実行:

    % ./a.out

ここで `control + z` を押すと、プロセスはバックグラウンドで一時停止する:

    [1]  + 19310 suspended  a.out

バックグラウンドのプロセスを確認するには `jobs` コマンド。
左からジョブ番号、カレントジョブか否か、状態、コマンド内容:

    % jobs
    [1]  + suspended  a.out

フォアグラウンドで再開するには `fg` コマンド。
引数としてジョブ番号かプログラム名をパーセントに続けて指定する
(`%1` とか `%a.out` とか)。
[zsh]({{< relref "zsh.md" >}}) なら補完もしてくれる。
引数を省略すると、`jobs` で `+` がついてるカレントジョブが選択される。
:

    % fg %1
    [1]  + 19310 continued  a.out

再び一時停止して、今度はバックグラウンドで再開する。コマンドは `bg`:

    [1]  + 19310 suspended  a.out
    % bg %a.out

始めからバックグラウンドで走らせるなら末尾にアンド:

    % ./a.out &

中止させるには `kill [-9] {job}` 。
あるいは `ps` 等で調べて生のプロセスIDを指定してもよい。
`-9` オプションをつけると `SIGKILL` で強制終了。
`killall` は指定した文字列と前方一致するプロセスを
すべて `kill` する。

    % kill %1
    % kill -9 %a.out
    % kill 19310
    % killall a.out

## ログアウト後も継続

バックグラウンドで実行中のプロセスも、
ログアウトするとhangup(`HUP`)シグナルによって終了してしまう。
これを無視して実行し続けるようにプログラムを起動するのが `nohup`。
バックグランド化までは面倒見てくれないので末尾の `&` を忘れずに:

    % nohup COMMAND [ARGUMENTS] &

標準出力と標準エラー出力は指定しなければ `nohup.out`
または `$HOME/nohup.out` に追記モードでリダイレクトされる。

> -   標準出力の書き出し先を指定するには `>{OUTFILE}`
> -   標準エラー出力の書き出し先を指定するには `2>{OUTFILE}`
> -   標準エラー出力を標準出力と同じところに流すには `>{OUTFILE} 2>&1`

    % nohup COMMAND >out.log 2>err.log &

[ssh]({{< relref "ssh.md" >}}) 接続先のサーバーで `nohup` ジョブを走らせるときは
標準入出力をすべて切っておかないと期待通りにsshを抜けられない:

    % nohup COMMAND >/dev/null 2>&1 </dev/null &

うっかり普通に開始してしまったプロセスを後から `nohup` 状態にするには、
一時停止して、バックグラウンドで再開して、`disown` に渡す:

    % ./a.out  # control + z
    [1]  + 19310 suspended  a.out
    % bg %1
    [1]  + 19310 continued  a.out
    % disown %1

{{%div class="note"%}}
[tmux]({{< relref "tmux.md" >}})

`nohup`, `disown` がプロセス単位で切り離すのに対して、
`tmux` は端末セッション丸ごと切り離し＆復帰することができる。
{{%/div%}}
