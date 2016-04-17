+++
title = 'Shell Script'
[menu.main]
  parent = "dev"
+++

<https://www.gnu.org/software/bash/manual/html_node/Shell-Parameter-Expansion.html>

## Parameter Expansion

```sh
string=0123456789abcdef
array=(0 1 2 3 4 5 6 7 8 9 a b c d e f)
echo ${string[1]}  # 0
echo ${array[1]}   # 0
echo ${string[-1]} # f
echo ${array[-1]}  # f
echo ${#string}    # 16
echo ${#array}     # 16
```

### Special Parameters

```sh
$0        # The full path of the script
$1        # First argument
$#        # The number of arguments ($0 is not included)
$*        # 全ての引数。"$*" はひと括り: "$1 $2 $3"
$@        # 全ての引数。"$@" は個別括り: "$1" "$2" "$3"
$@[2,-1]  # 2つ目以降の引数(zshのみ？)
$-        # シェルの起動時のフラグ、setコマンドを使って設定したフラグの一覧
$$        # シェル自身のプロセスID
$!        # シェルが最後に起動したバックグラウンドプロセスのプロセスID
$?        # 最後に実行したコマンドのexit値
```

### Character String

文字列の置換には正規表現を使える `sed` が便利。
でもパスや拡張子のちょっとした操作なら、
わざわざ外部コマンドやパイプを使わずにシェル(`bash`)の機能だけで実現できる。

-   `*` グロブの最短マッチか最長マッチが使える。（正規表現は使えない）
-   `$@` や `$*` に対しては、個々の要素に対して実行されてリストが返される。

```sh
${STR#pattern}        # sed -e "s/^pattern//")     shortest
${STR##pattern}       # sed -e "s/^pattern//")     longest
${STR%pattern}        # sed -e "s/pattern$//")     shortest
${STR%%pattern}       # sed -e "s/pattern$//")     longest
${STR/pattern/repl}   # sed -e "s/pattern/repl/")
${STR//pattern/repl}  # sed -e "s/pattern/repl/g")
${STR/#pattern/repl}  # sed -e "s/^pattern/repl/")
${STR/%pattern/repl}  # sed -e "s/pattern$/repl/")

BASENAME=${PATH##*/}
DIRNAME=${PATH%/*}
EXTENTION=${PATH##*.}
BASENAME=${FILENAME%.jpg}
```

### Other

```sh
## expression  # if VAR is null
${VAR:-WORD}  # return WORD; VAR remains null
${VAR:=WORD}  # return WORD; VAR is assigned WORD
${VAR:?WORD}  # display error and exit
${VAR:+WORD}  # nothing occurs; otherwise return WORD
```

```sh
${!PREFIX*}  # the names of variables begins with PREFIX
${!PREFIX@}
```

## Misc

### Command

コマンドの結果を文字列として受け取る方法は
`` `command ``\` と `$(command)` の2つあるが、
後者のほうが入れ子など柔軟に使える。

### Arithmetic

足し算くらいなら `$((expression))` で

### コア数の取得

Linuxなら `/proc/cpuinfo`、
Macなら `system_profiler` 。でもたぶんIntelとPPCは書き方が違う。

```sh
[[ -r /proc/cpuinfo ]] \
&& CORES=$(grep cpuid /proc/cpuinfo | wc -l)
|| CORES=system_profiler | grep Cores | awk '{print $5}'
```
