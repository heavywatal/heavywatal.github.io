+++
title = 'シェルスクリプト'
tags = ["shell"]
[menu.main]
  parent = "dev"
+++

https://www.gnu.org/software/bash/manual/html_node/


## if-then

- [Conditional Constructs](https://www.gnu.org/software/bash/manual/html_node/Conditional-Constructs.html)
- [Bash Conditional Expressions](https://www.gnu.org/software/bash/manual/html_node/Bash-Conditional-Expressions.html)
- `[` はPOSIXコマンドで `]` は最後の引数。
- `[[ expression ]]` は便利な面もあるけどBash拡張なので注意。

```sh
# basic
if test -e ~; then
  echo 'test -e ~'
fi

# popular
if [ -e ~ ]; then
  echo '[ -e ~ ]'
fi

# bash/zsh extension; not POSIX
if [[ -e ~ ]]; then
  echo '[[ -e ~ ]]'
fi

# shortcut with exit status
test -e ~ && echo 'test -e ~ && echo' || echo 'not printed'
[ -e ~ ] && echo '[ -e ~ ] && echo' || echo 'not printed'
[[ -e ~ ]] && echo '[[ -e ~ ]] && echo' || echo 'not printed'
[[ ! -e ~ ]] && echo 'not printed' || echo '[[ ! -e ~ ]] || echo'
```

AND/OR
```sh
[ -e ~ ] && [ -d ~ ] && echo '[ -e ~ ] && [ -d ~ ]'
[ -e ~ ] || [ -f ~ ] && echo '[ -e ~ ] || [ -f ~ ]'
[[ -e ~ && -d ~ ]] && echo '[[ -e ~ && -d ~ ]]'
[[ -e ~ || -f ~ ]] && echo '[[ -e ~ || -f ~ ]]'
```

文字列
```sh
EMPTY=""
NOTEMPTY="content"
[ -z "$EMPTY" ] && echo '-z "$EMPTY"'
[ -n "$NOTEMPTY" ] && echo '-n "$NOTEMPTY"'
[ "$NOTEMPTY" != "$EMPTY" ] && echo '"$NOTEMPTY" != "$EMPTY"'
[ "$NOTEMPTY" = "content" ] && echo '"$NOTEMPTY" = "content"'
[[ "$NOTEMPTY" =~ "tent$" ]] && echo '"$NOTEMPTY" =~ "tent$"'
```

ファイル、ディレクトリ
```sh
x="${HOME}/.bashrc"
[ -e $x ] && echo "x exists"
[ -d $x ] && echo "x is a directory"
[ -f $x ] && echo "x is a regular file"
[ -s $x ] && echo "x is a regular file with size >0"
[ -L $x ] && echo "x is a symlink"
[ -r $x ] && echo "x is readable"
[ -w $x ] && echo "x is writable"
[ -x $x ] && echo "x is executable"
```

数値比較
```sh
x=1
y=2
[ $x -eq $y ] && echo "x == y"
[ $x -ne $y ] && echo "x != y"
[ $x -lt $y ] && echo "x <  y"
[ $x -le $y ] && echo "x <= y"
[ $x -gt $y ] && echo "x >  y"
[ $x -ge $y ] && echo "x >= y"
```


## Looping Constructs

<https://www.gnu.org/software/bash/manual/html_node/Looping-Constructs.html>

基本形。クオートもカッコも不要:
```sh
for x in 1 2 3; do
  echo $x
done
# 1
# 2
# 3
```

`for x in "1 2 3"` のようにクオートすれば当然1要素扱いになるが、
一旦変数に代入したあとの挙動はシェルによって異なるので要注意:
zshでは1要素扱い、dashやbashではスペース区切りで分割。
```sh
STRING="1 2 3"
for x in $STRING; do
  echo $x
done
```

## Special Parameters

<https://www.gnu.org/software/bash/manual/html_node/Special-Parameters.html>

```sh
$0        # The full path of the script
$1        # First argument
$#        # The number of arguments ($0 is not included)
$*        # 全ての引数。"$*" はひと括り: "$1 $2 $3"
$@        # 全ての引数。"$@" は個別括り: "$1" "$2" "$3"
$-        # シェルの起動時のフラグ、setコマンドを使って設定したフラグの一覧
$$        # シェル自身のプロセスID
$!        # シェルが最後に起動したバックグラウンドプロセスのプロセスID
$?        # 最後に実行したコマンドのexit値
```

## Parameter Expansion

<https://www.gnu.org/software/bash/manual/html_node/Shell-Parameter-Expansion.html>

```sh
string=abcdef
echo ${string}   # abcdef
echo ${#string}  # 6
```

うっかりアンダースコアで変数名をつなげてしまいがちなので、
常に{カッコ}つける癖をつける:

```sh
for CITY in sendai yokosuka; do
  echo hello_$CITY_people     # hello_
  echo hello_${CITY}_people   # hello_sendai_people
done
```

変数が未定義or空だったり、空白や特殊文字を含んでいたりして事故りがちなので、
基本的にダブルクォートを付ける癖をつける:

```sh
EMPTY=""
[ -n $EMPTY ] && echo "Bug! this is printed"
[ -n "$EMPTY" ] && echo "OK, this is not printed."

DIR=${HOME}/directory with space
# with: command not found
echo $DIR
# /Users/watal/directory

DIR="${HOME}/directory with space"
echo $DIR
# /Users/watal/directory with space
cd $DIR
# cd: /Users/watal/directory: No such file or directory
cd "$DIR"
# cd: /Users/watal/directory with space: No such file or directory
```

シングルクォートの中の変数は展開されずそのまま渡される:
```sh
echo "$HOME"  # /Users/watal
echo '$HOME'  # $HOME
```


### 部分列

```sh
# ${parameter:offset}
# ${parameter:offset:length}
echo ${string:0:3}  # abc
echo ${string:1:3}  # bcd
echo ${string:2}    # cdef
echo ${string: -2}  # ef
```

### 置換

文字列の置換には正規表現を使える `sed` が便利。
でもパスや拡張子のちょっとした操作なら、
わざわざ外部コマンドやパイプを使わずにシェル(`bash`)の機能だけで実現できる。

-   `*` グロブの最短マッチか最長マッチが使える。（正規表現は使えない）
-   `$@` や `$*` に対しては、個々の要素に対して実行されてリストが返される。

```sh
# ${STR#pattern}        # sed -e "s/^pattern//")     shortest
# ${STR##pattern}       # sed -e "s/^pattern//")     longest
# ${STR%pattern}        # sed -e "s/pattern$//")     shortest
# ${STR%%pattern}       # sed -e "s/pattern$//")     longest
# ${STR/pattern/repl}   # sed -e "s/pattern/repl/")  longest
# ${STR//pattern/repl}  # sed -e "s/pattern/repl/g") longest
# ${STR/#pattern/repl}  # sed -e "s/^pattern/repl/") longest
# ${STR/%pattern/repl}  # sed -e "s/pattern$/repl/") longest

FILEPATH=/root/dir/file.tar.gz
echo ${FILEPATH##*/}       # file.tar.gz
echo ${FILEPATH##*.}       # gz
echo ${FILEPATH#*.}        # tar.gz
echo ${FILEPATH%/*}        # /root/dir
echo ${FILEPATH%.*}        # /root/dir/file.tar
echo ${FILEPATH%%.*}       # /root/dir/file
echo ${FILEPATH/%.*/.zip}  # /root/dir/file.zip
```

変数が設定されていない場合にどうするか
```sh
              # if VAR is null, then
${VAR:-WORD}  # return WORD; VAR remains null
${VAR:=WORD}  # return WORD; WORD is assigned to VAR
${VAR:?WORD}  # display error and exit
${VAR:+WORD}  # nothing occurs; otherwise return WORD
```


## Misc.

### Command Substitution

<https://www.gnu.org/software/bash/manual/html_node/Command-Substitution.html>

コマンドの結果を文字列として受け取る方法は
`` `command` `` と `$(command)` の2つあるが、
後者のほうが入れ子など柔軟に使えるのでより好ましい。

### Arithmetic Expansion

<https://www.gnu.org/software/bash/manual/html_node/Arithmetic-Expansion.html>

簡単な数値計算は `$((expression))` でできる。

[Shell Arithmetic](https://www.gnu.org/software/bash/manual/html_node/Shell-Arithmetic.html)

### Process Substitution

<https://www.gnu.org/software/bash/manual/html_node/Process-Substitution.html>

1つのコマンドから出力を受け取るだけならパイプ `command1 | command2` で足りるけど、
2つ以上から受け取りたいときはプロセス置換 `command2 <(command1a) <(command1b)` を使う。

標準入力や標準出力を受け付けないコマンドでも、
ファイル名を引数で指定できればプロセス置換で対応できる。
例えばgzip圧縮・展開の一時ファイルを作りたくない場合とか:
```sh
somecommand infile outfile
gunzip -c infile.gz | somecommand /dev/stdin /dev/stdout | gzip -c >outfile.gz
somecommand <(gunzip -c infile.gz) >(gzip -c >outfile.gz)
```

### Arrays

<https://www.gnu.org/software/bash/manual/html_node/Arrays.html>

POSIXでは未定義の拡張機能であり、bashとzshでも挙動が異なる。
```sh
array=(a b c d e f)  # bash        | zsh
echo ${array}        # a           | a b c d e f
echo ${array[0]}     # a           |
echo ${array[1]}     # b           | a
echo ${array[-1]}    # (error)     | f
echo ${array[@]}     # a b c d e f | a b c d e f
echo ${#array}       # 1           | 6
echo ${#array[@]}    # 6           | 6
```

文字列と同様にコロンを使う方法ならzshでもゼロ基準なので安全。
```sh
echo ${array[@]:0:3}  # a b c
echo ${array[@]:1:3}  # b c d
echo ${array[@]:2}    # c d e f
echo ${array[@]: -2}  # e f
```

`array[*]` と`array[@]` はクオートで括られる挙動が異なる。
前者はひと括り、後者はそれぞれ括られる。
```sh
for x in "${array[*]:3}"; do
  echo $x
done
# d e f

for x in "${array[@]:3}"; do
  echo $x
done
# d
# e
# f
```

### オプション

シェル起動時に `bash -u script.sh` とするか、
スクリプト内で `set -u` とする。
途中で `set +u` として戻したりもできる。

後述のように注意は必要だけどとりあえず走りすぎない設定で書き始めたい:
```sh
set -eux -o pipefail
```

`-e`
: エラーが起きたらすぐ終了する。
  これが無ければスクリプトの最後まで走ろうとする。
: 挙動が変わる要因が多すぎて
  (サブシェルか、条件文がついてるか、POSIXモードか、bashバージョンなど)
  逆に難しくなるので使わないほうがいいという見方もある。
  それが気になるほど込み入ってきたらもはやシェルスクリプトの出番ではなく、
  Pythonとかでかっちり書いたほうがよいのでは。
: `source` している中で `exit` するとシェルごと落ちることにも注意。

`-u`
: 定義されていない変数を展開しようとするとエラー扱い。

`-v`
: シェルに入力されるコマンドを実行前にそのまま表示する。
  よくあるverboseオプションなので、
  スクリプトに書いてしまうより必要に応じてコマンドに足すほうが自然か。

`-x`
: 変数を展開してからコマンドを表示する。
  `-v` との併用も可能。

`-o pipefail`
: パイプの左側でエラーになったら終了する。
  読み手側が早めに切り上げて `SIGPIPE` が発生しても止まるので注意。
: POSIXではないがbashでもzshでも利用可能。

`-o posix`
: bash/zsh拡張を切ってPOSIX準拠モードで動く。
  `sh` として呼び出されるとこれになる。


## 関連書籍

<a href="https://www.amazon.co.jp/dp/4797380942?&linkCode=li3&tag=heavywatal-22&linkId=faccc593edc22f1e89bf63f13e924766&language=ja_JP&ref_=as_li_ss_il" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4797380942&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4797380942" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
