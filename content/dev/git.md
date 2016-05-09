+++
title = 'Git'
subtitle = '分散型バージョン管理システム'
tags = ["vcs", "writing"]
[menu.main]
  parent = "dev"
+++

https://git-scm.com/

Gitは分散型バージョン管理システムの代表格。
プログラムのソースコードはもちろんのこと、
研究ノートや論文の原稿などあらゆるテキストの管理に使える。

[GitHub](https://github.com)はGitをより便利に使うためのオンラインサービス。
個人的なリポジトリ置き場としてはもちろんのこと、
ほかの人と共有・協力してプロジェクトを進めるプラットフォームとしても使える。

Gitのライバルとして[Mercurial]({{< relref "mercurial.md" >}})もあるが、
[BitBucket](https://bitbucket.org) (GitHubのライバル)
がGit対応した今となってはMercurialを積極的に使う理由は無い気がする。
むしろ[Atom]({{< relref "atom.md" >}})における差分表示など、
Gitでなければ得られない恩恵が大きくなってきている。


## 基本

### 手元の変更を外に伝える

working directory (working tree)
: 手元のファイルの変更はまだリポジトリに登録されていない
: ↓ `add`

staging area (index)
: 次のコミットに含めるファイルをマークする段階
: ↓ `commit`

local repository
: 変更履歴が `.git/` 内に記録されている
: ↓ `push`

remote repository
: GitHubなど別マシンのリポジトリに反映


### 外部の変更を手元に取り込む

remote repository
: ↓ `fetch`

local repository
: 変更が `.git/` に取り込まれたが、見えてるファイルには反映されてない
: ↓ `checkout` or `merge`

working directory
: 手元のファイルが最新版に同期されている


### 用語

blob
: git内部で1つのファイルを指すオブジェクトで、add時に作られる。
  ファイル名などのメタデータは持たず、
  ファイルの内容にのみ依存したハッシュIDを持つ。

tree
: git内部で1つのディレクトリを指すオブジェクトで、commit時に作られる。
  blobやファイル名などのメタデータに依存したハッシュIDを持ち、
  その変化は親に伝播する。

commit
: git内部でroot treeのsnapshotを指すオブジェクト。
  root treeのハッシュID、著者、コメントなどの情報を持つ。
  動詞としては、staging areaの情報をひとつのcommitとしてリポジトリに登録することを指す。

repository
: commitの履歴を保持する拠点。
  `git init` で手元に新規作成するか、`git clone` でリモートから複製する。

`origin`
: remoteリポジトリの典型的なshortname。
  clone時に自動的に追加され、
  push先やfetch元を省略したときにデフォルトで使われる。
  `git remote -v` で確認。

`master`
: デフォルトのブランチの典型的な名前。

`HEAD`
: 現在checkoutしているbranch/commitを指すポインタ。
  基本的にはmasterの最新コミットを指していることが多い。
  ひとつ前は`HEAD^`、ふたつ前は`HEAD^^`。


## よく忘れるコマンド

直前のcommitを修正・上書き (コメント修正やファイルの追加忘れに有用):

    git commit --amend

直前の動作を取り消す:

    # 直前のcommitを取り消す (indexとworkingはそのまま)
    git reset --soft HEAD^

    # 直前のaddを取り消す (workingはそのまま)
    git reset --mixed HEAD

    # working directoryの変更を取り消す, DANGEROUS!
    git reset --hard HEAD
    # それでもuntrackedは残るので、消したければ git clean

    # 直前のresetを取り消す
    git reset --hard ORIG_HEAD

tracking対象から外して忘れさせる(手元のファイルはそのまま):

    git rm --cached <file>

差分を表示:

    # HEAD vs working (staging前のファイルが対象)
    git diff

    # HEAD vs index (staging済みcommit前のファイルが対象)
    git diff --staged

    # HEAD vs working+index (commit前の全ファイルが対象)
    git diff HEAD

    # 最新コミットの変更点
    git diff HEAD^ HEAD


## Submodule

### 既存のリポジトリをsubmoduleとして追加する

```sh
git submodule add https://github.com/mbostock/d3.git

# ブランチを指定する場合:
git submodule add -b gitsubmodule_https https://github.com/heavywatal/x18n.git
```

{{%div class="note"%}}
gh-pagesで公開する場合は参照プロトコルを
`git://` ではなく `https://` にする必要がある。
{{%/div%}}


### submoduleを含むメインリポジトリを使い始めるとき

最初にclone/fetchしてきた時submoduleたちは空なのでまず:

    git submodule update --init

    # 使いたいbranchがmasterではない場合は --remote
    git submodule update --init --remote x18n

    # 歴史があって重いリポジトリはshallowに
    git submodule update --init --depth=5 d3

### submoduleを更新

1.  更新分をまとめて取得:

        git submodule foreach git fetch

2.  好きなとこまでチェックアウト:

        cd d3/
        git checkout v3.5.6

3.  メインリポジトリでその変更をコミット:

        cd ..
        git commit



## GitHub Pages

### ユーザーサイトを作る

1. `USERNAME.github.io` という名前のリポジトリをGitHub上で作成
2. 公開したいウェブサイトをmasterブランチとしてpush
3. `https://USERNAME.github.io` にアクセスしてみる。

例えば本ウェブサイトは
`heavywatal.github.io` というリポジトリの
`source`ブランチでMarkdownテキストを書き、
[Hugo]({{< relref "hugo.md" >}})
で変換・生成したHTMLファイルを`master`ブランチに書き出している。


### プロジェクトサイトを作る

`gh-pages` ブランチの内容が
`https://USERNAME.github.io/PROJECT/` で公開される。

e.g., doxygen生成物を公開

1.  メインリポジトリで `git checkout --orphan gh-pages`
2.  index.html など必要なファイルを用意して最初のコミット:

        git add --all
        git commit -a

3.  `git push origin gh-pages`
4.  メインリポジトリに戻る:

        git fetch
        git checkout master

5.  doxygenの出力先(e.g., html)をサブモジュールにする:

        git submodule add -b gh-pages `git remote -v|grep origin|head -n 1|awk '{print$2}'` html


## プルリクエストを送る

-   masterブランチは更新取得のためだけに使い、開発は別ブランチで行う。
-   大元のリポジトリをupstream、自分のリポジトリをoriginと名付ける。
-   push済みのcommitをrebaseするとIDが変わっちゃうのでダメ。

### 基本の流れ

例えば `user` さんの `project` のドキュメント修正に貢献する場合。

1.  GitHub上で元のリポジトリから自分のリポジトリにフォークする
2.  forkした自分のリポジトリからローカルにcloneする:

        git clone git@github.com:heavywatal/project.git
        cd project/

3.  元のリポジトリにupstreamという名前をつけておく:

        git remote add upstream git://github.com/user/project.git

4.  開発用のブランチを切って移動:

        git checkout -b fix-typo

5.  ソースコードに変更を加える:

        emacs README.md
        git commit -a -m "fix typo in README.md"

6.  この間に起こった元リポジトリの更新をmaster越しに取り込む:

        git checkout master
        git fetch upstream
        git merge upstream/master
        git checkout fix-typo
        git rebase -i master

7.  自分のリポジトリにpush:

        git push origin fix-typo

8.  GitHub上の自分のリポジトリからプルリクエストを送る
9.  修正を求められたらそのブランチで変更して普通にプッシュ


## 問題と対処

### detached HEAD からの復帰

submoduleなどをいじってると意図せずdetached HEAD状態になることがある。
その状態でcommitしてしまった変更をmasterに反映したい。

1. pushしようとして怒られて気付く
   ```
   % git push
   fatal: You are not currently on a branch
   % git status
   HEAD detached from *******
   ```

2. masterに戻ると道筋を示してくれる:
   ```
   % git checkout master
   Warning: you are leaving 2 comits behind, not connected to
   any of your branches
   If you want to keep them by creating a new branch, this may be a good time
   to do so with:

    git branch <new-branch-name> *******
   ```

3. 言われたとおりbranchを作ってmerge
   ```
   % git branch detached *******
   % git merge detached
   ```

4. 不要になったbranchを消す
   ```
   % git branch -d detached
   ```

