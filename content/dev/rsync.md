+++
title = 'rsync'
tags = ["communication"]
+++

ファイルをコピーし、2つのディレクトリを同期する。
更新があったものだけをコピーする、
ひとつのsshセッションで複数のファイルを送受信する、
という使い方が可能なので `cp` や `scp` よりも便利な場合が多い。

- <https://rsync.samba.org/>
- <https://github.com/WayneD/rsync>

## 基本

オプションについては後述するとして、基本は:
```sh
rsync -auvC SRC/ DST/
```

SRC側の末尾のスラッシュの有無によって結果が大きく異なることに注意。
DST側は付けても付けなくても同じ。

```sh
## SRC/ 以下のファイルが DST/ 以下に入る
rsync -auvC SRC/ DST

## ディレクトリ DST/SRC が作られる
rsync -auvC SRC DST

## 結果は同じだが、下のほうがより明示的
rsync -auvC SRC/DIR DST
rsync -auvC SRC/DIR/ DST/DIR/
```

### ネットワーク越しの送受信

SSHを使ってリモートホストとの送受信も可能。
その場合は相手側を `remote_machine:~/dir` のようにコロンから絶対パスで指定する。

```sh
# push
rsync -auvC SRC/ user@example.com:~/DST/

# pull
rsync -auvC user@example.com:~/SRC/ DST/
```

[SSHの設定]({{< relref "ssh.md" >}})をちゃんとやっておくと楽。
- パスワード無しでログインできる。
- ユーザー名を省略したり、ホスト名を短縮した別名をつけたりもできる。
- `RequestTTY yes` を付けてると当然ながら怒られる。
- リモート側の `.bashrc` とかで標準出力に何かを表示するようにしてあるとコケる


## Options

`-a, --archive`
:   バックアップ用途のオプション一式 `-rlptgoD`

`-u, --update`
:   受け手の方が新しいファイルをスキップ

`-v, --verbose`
:   冗長なメッセージ表示

`-n, --dry-run`
:   実際に送受信を行わず試してみる

`-z, --compress`
:   圧縮・展開のCPUコストはかかるけど通信量は減る

`--delete`
:   SRC側に存在しないものがDST側にあったとき削除

`--delete-excluded`
:   除外設定されているファイルが受け手側にあったら削除（危険！）

`--ignore-existing`
:   受信側に存在していたら無視

`--iconv`
:   文字コードを変換する。
    異なるOS/FS間でウムラウトや日本語を含むファイルを送受信するときに使う。
    リモート側だけ指定すれば十分だけどローカルも明示的に指定できる。
    その場合、pushかpullかによらず `={LOCAL},{REMOTE}` の順番で。
    例えばローカルのLinuxマシンとリモートの古いMacでやり取りする場合は
    `--iconv=utf8-mac` or `--iconv=utf8,utf8-mac` 。
    新しいMacのAPFSはLinuxと同じ `utf8` と見なしておけば良さそう...?


### Exclude and include

<https://download.samba.org/pub/rsync/rsync.1#PATTERN_MATCHING_RULES>

- 先に記述したものほど優先。
- 正規表現ではなくglob寄りの独自文法。
- 基本的にはpathの最終コンポーネント(basename)の部分が評価対象。
- 上位ディレクトリから順に評価し、除外されたらそれより下は読みに行かない。
- 末尾が `/` ならディレクトリにのみマッチ。

`--include=<PATTERN>`
:   マッチするファイル・ディレクトリを除外しない

`--exclude=<PATTERN>`
:   マッチするファイル・ディレクトリを除外

`--exclude-from=<FILE>`
:   ファイルに記述した除外パターンを読む

`-C`, `--cvs-exclude`
:   ほぼいつでも無視したいであろうものを除外する。
    基本的にはこれを使う方針で大丈夫そうだが
    `tags`, `core` とかはディレクトリ名として普通に使いそうなので注意。
:   `RCS SCCS CVS CVS.adm RCSLOG cvslog.* tags TAGS .make.state .nse_depinfo *~ #* .#* ,* _$* *$ *.old *.bak *.BAK *.orig *.rej .del-* *.a *.olb *.o *.obj *.so *.exe *.Z *.elc *.ln core .svn/ .git/ .hg/ .bzr/`
:   これに加えて **[`${HOME}/.cvsignore`](https://github.com/heavywatal/dotfiles/blob/master/.cvsignore) も自動で読まれる**というのが便利。
    リモートとのやり取りの場合、ソース側のファイルが読まれる(コマンド実行側とは限らない)ということに注意。


## wget

<https://www.gnu.org/software/wget/manual/wget.html>

`-nv`, `--no-verbose`
: 基本が結構うるさいverboseなので、最低限の表示に抑える。
  `-q`, `--quiet` よりはマイルド。

`-P`, `--directory-prefix`
: 出力先のディレクトリ。デフォルトは `.` カレント。

`-O`
: 単一の出力ファイルを指定。
  シェルの `>` リダイレクトと同じように、既存のファイルはまず消される。
  ハイフン `-` で標準出力。
: 何もつけなければURLの最後の部分がファイル名になる。
  これを繰り返すと `file.1` のように連番が付く。

`-nc`, `--no-clobber`
: 出力先ファイルが存在する場合、サイズや時間によらず連番も取得もしない。
  使い所がわからない。

`-c`, `--continue`
: 出力先ファイルが存在する場合、連番をつけずに部分ファイルとして扱う。
  不足分だけ取得して通信量を節約できるとは限らない。

`-N`, `--timestamping`
: ローカルのファイルが新しい場合は取得しない。
  時間だけでなくサイズの違いも考慮してくれる。

### Recursive retrieval options

これの存在が curl との大きな違い。

`-m`, `--mirror`
: `-r -N -l inf --no-remove-listing`

`-r`, `--recursive`
: 再起的に取得。
  ftpでは普通にディレクトリを潜っていく。
: http(s)ではリンクを辿ることを意味する。
  一見ディレクトリ指定のような
  `https://example.com/dir/` を渡した場合でも、
  `--default-page` (`index.html`) ファイルを取得してそのリンクを辿る。

`-l`, `--level`
: 再帰の最大深度を指定。0で無制限。デフォルトは5。

`-A`, `--accept`, `-R`, `--reject`, `--accept-regex`, `--reject-regex`
: 指定したsuffixあるいはパターンにマッチするものだけを取得、あるいは除外。
  comma-separated list。
: ディレクトリ版は `-I`, `--include-directories`, `X`, `--exclude-directories`

`-H`, `--span-hosts`
: 外部ホストへのリンクも辿る。

`-L`, `--relative`
: 相対リンクのみ辿る。

`-np`, `--no-parent`
: ディレクトリを遡らない。
  `-r` を使うときは大概この挙動を意図してるはず。
: URL末尾のスラッシュ `/` をちゃんとつけないとそのディレクトリが起点にならないので注意。

`-nH`, `--no-host-directories`
: ホスト名のディレクトリを作らない。
  デフォルトでは `example.com/file` のようにホスト名のディレクトリが作られる。
: 逆に確実に作りたい場合は `-x`, `--force-directories` 。

`--cut-dirs`
: 取得元のディレクトリ階層を上から指定した数だけ無視して下層だけ作る。

| option                    | output |
|---------------------------|--------|
| `-r -np -P tmp`           | `tmp/example.com/dist/dir/file` |
| `-r -np`                  | `example.com/dist/dir/file` |
| `-r -np -nH`              | `dist/dir/file` |
| `-r -np -nH --cut-dirs=1` | `dist/file` |
| `-r -nd`                  | `file` |


## curl

<https://curl.se/>

大概どのOSでも最初から入ってるので、
単一ファイルを取得するには手軽で便利。
インストールスクリプトなどで使われることも多い。

`-f`, `--fail`
: 403 Forbidden や 404 Not Found のときにそのページを取得せずエラー終了する。

`-s`, `--silent`
: プログレスバーやエラーの表示を抑制する。
: `-S`, `--show-error` を追加してエラー表示だけ有効にすることが多い。

`-L`, `--location`
: 3xxリダイレクトを辿る。

`-o`, `--output`
: 出力先ファイルを指定。
  デフォルト無指定もしくはハイフン `-` で標準出力。
: URLに `{one,two}` や `[1-3]` のような表現を含めて複数ファイルを取得する場合は
  `file#1.txt` のようなプレースホルダが使える。

`-O`, `--remote-name`
: URLの最後の部分をファイル名として保存。
: サーバー側が指定した名前を採用したい場合は
  `-J`, `--remote-header-name` を使う。

`-R`, `--remote-time`
: リモートファイルのタイムスタンプを取得してローカルファイルに適用。

`--output-dir`
: 出力先ディレクトリを指定。
  デフォルトは `.` カレント。
  存在しないディレクトリを作ってもらうには `--create-dirs` も必要。

`--no-clobber`
: 出力先ファイルが存在する場合、上書きせず連番をつける。

