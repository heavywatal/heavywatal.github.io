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

### 準備

- ローカルマシンにGitをインストールする。
  MacならCommand Line Toolsに付属のを使うか、[Homebrew]({{< relref "homebrew.md" >}})で新しいのを入れる。
  Linuxなら初めから入ってるのを使うか、Linuxbrewで新しいのを入れる。
- GitHubに個人アカウントを作る。
- [SSH公開鍵を作って]({{< relref "ssh.md" >}})マシンとGitHubに登録する。
  https://help.github.com/articles/connecting-to-github-with-ssh/
- `~/.gitconfig` にユーザ名やアドレスを登録する。
  https://git-scm.com/docs/git-config

    ```sh
    git config --global user.name "Watal M. Iwasaki"
    git config --global user.email "heavy.watalあmail.com"
    less ~/.gitconfig
    ```
- ついでに `pushinsteadof` の設定をしておく。
  httpsで高速にclone/fetch/pullして、
  sshでパスワード無しでpushする、というのが楽ちん。
- 設定例: [`heavywatal/dotfiles/.gitconfig`](https://github.com/heavywatal/dotfiles/blob/master/.gitconfig)


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
: git内部で1つのファイルを指すオブジェクトで、`add`時に作られる。
  ファイル名などのメタデータは持たず、
  ファイルの内容にのみ依存したハッシュIDを持つ。

tree
: git内部で1つのディレクトリを指すオブジェクトで、`commit`した時に作られる。
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
  `clone`時に自動的に追加され、
  `push`先や`fetch`元を省略したときにデフォルトで使われる。
  `git remote -v` で確認。

`master`
: デフォルトのブランチの典型的な名前。

`HEAD`, `@`
: 現在checkoutしているbranch/commitを指すポインタ。
  基本的には`master`の最新commitを指していることが多い。
  1つ前は `HEAD^` か `HEAD~`、
  2つ前は `HEAD^^` か `HEAD~~` か `HEAD~2`。
  (`HEAD^2` は `merge` で複数の親がある場合の2番目)

zshの`EXTENDED_GLOB`が有効になってる場合は
`HEAD^` がパターン扱いされてエラーになるので、
`HEAD\^` のようにエスケープするか `unsetopt NOMATCH` しておいたほうがいい。


## よく使うコマンド

### reset

`git reset <DESTINATION>` は `HEAD` の位置を戻す処理で、
オプションによってindexとworing treeもそこに合わせるように変更される。
`--soft` なら `HEAD` 移動のみ。
`--mixed` なら移動した `HEAD` にindexも合わせる。
`--hard` なら移動した `HEAD` にindexとworiking treeも合わせる。
直前の動作を取り消す用途に絞って使うのが無難:
```sh
# commit直後、それを取り消す (indexとworkingはそのまま)
git reset --soft HEAD^

# add直後、それを取り消す (workingとHEADはそのまま)
git reset --mixed HEAD

# 変更したファイルをHEADの状態に戻す (DANGEROUS!)
git reset --hard HEAD

# reset直後、それを取り消す
git reset --hard ORIG_HEAD

# divergedになってしまった手元のbranchを破棄 (DANGEROUS!)
git reset --hard origin/master
```

直前のcommitをちょっと修正したいだけなら `git commit --amend` が簡単。
それより前のを修正するには `git rebase -i HEAD~3` とかで戻ってrewordやedit。


### diff

差分を表示:
```sh
# HEAD vs working (staging前のファイルが対象)
git diff

# HEAD vs index (staging済みcommit前のファイルが対象)
git diff --staged

# HEAD vs working+index (commit前の全ファイルが対象)
git diff HEAD

# 特定コミットの変更点 (diffじゃない...)
git show [revision]
```

### rm, clean

tracking対象から外して忘れさせる(手元のファイルはそのまま):
```sh
git rm --cached <file>
```

`.gitignore` で無視されてるuntrackedファイルを消す:
```sh
git clean -fdX
```
無視されていないuntrackedファイルも消したい場合は小文字の `-fdx` (危険)。

### tag

特定のコミット(省略すると`HEAD`)にタグ付けする。
lightweightとannotatedの2種類が存在し、後者にはメッセージなどが紐付く。
```sh
git tag v0.1.0 [revision]
git tag -a v0.1.0 -m "Message!"
```

GitHubリポジトリに `git push --tags` するとアーカイブが作られ、
Releasesページに反映される。
annotated tagであれば `git push --follow-tags`
でcommitとtagを同時にpushできる。

https://git-scm.com/book/en/Git-Basics-Tagging


## Submodule

### 既存のリポジトリをsubmoduleとして追加する

```sh
git submodule add https://github.com/mbostock/d3.git

# ブランチを指定する場合:
git submodule add -b gitsubmodule_https https://github.com/heavywatal/x18n.git
```

`gh-pages` で公開する場合は参照プロトコルを
`git://` ではなく `https://` にする必要がある。


### submoduleを含むメインリポジトリを使い始めるとき

最初に`clone`/`fetch`してきた時submoduleたちは空なのでまず:

    git submodule update --init

    # 使いたいbranchがmasterではない場合は --remote
    git submodule update --init --remote x18n

    # 歴史があって重いリポジトリはshallowに
    git submodule update --init --depth=5 d3

### submoduleを更新

1.  更新分をまとめて取得:

        git submodule foreach git fetch

1.  好きなとこまでチェックアウト:

        cd d3/
        git checkout v3.5.6

1.  メインリポジトリでその変更をコミット:

        cd ..
        git commit



## GitHub Pages

https://help.github.com/articles/user-organization-and-project-pages/

### ユーザーサイトを作る

1. `USERNAME.github.io` という名前のリポジトリをGitHub上で作成
1. 公開したいウェブサイトを`master`ブランチとして`push`
1. `https://USERNAME.github.io` にアクセスしてみる。

例えば本ウェブサイトは
`heavywatal.github.io` というリポジトリの
`source`ブランチでMarkdownテキストを書き、
[Hugo]({{< relref "hugo.md" >}})
で変換・生成したHTMLファイルを`master`ブランチに書き出している。

GitHubが勝手にJekyll処理しようとすることがあるので、
`.nojekyll` という空ファイルを作っておく。


### プロジェクトサイトを作る

https://help.github.com/articles/configuring-a-publishing-source-for-github-pages/

リポジトリの内容を `https://USERNAME.github.io/PROJECT/` に公開することができる。
方法は複数あり、リポジトリの設定画面から選択できる。

1. `gh-pages` ブランチの内容を公開 (古い方法)
1. `master` ブランチの内容を公開
1. `master` ブランチの `/docs` ディレクトリのみを公開

前は1番の方法しかなくてブランチの扱いがやや面倒だったが、
今では`master`だけで簡単に済ませられるようになった。

## Pull Request (PR)

-   大元のリポジトリを`upstream`、フォークした自分のリポジトリを`origin`と名付ける。
-   デフォルトブランチ(`master`とか`develop`とか)は更新取得のためだけに使い、変更は新規ブランチで行う。
-   `push`済みのcommitを`rebase`するとIDが変わっちゃうのでダメ。

### 基本の流れ

例えば `USER` さんの `PROJECT` のコード修正に貢献する場合。

1.  `github.com/USER/PROJECT` のForkボタンで自分のGitHubリポジトリに取り込む
1.  forkした自分のリポジトリからローカルに`clone`:

        git clone git@github.com:heavywatal/PROJECT.git
        cd PROJECT/

1.  大元のリポジトリに`upstream`という名前をつけておく:

        git remote add upstream git://github.com/USER/PROJECT.git

1.  PR用のブランチを切って移動:

        git checkout -b fix-typo

1.  コードを変更して`commit`:

        emacs README.md
        git diff
        git commit -a -m "Fix typo in README.md"

1.  この間に`upstream`で更新があったら、それをデフォルトブランチ越しに取り込む:

        git checkout master
        git fetch upstream
        git merge upstream/master
        git push origin/master
        git checkout fix-typo
        git rebase -i master

1.  自分のリポジトリに`push`:

        git push [-f] origin fix-typo

1.  GitHub上に出現する"Compare & pull request"ボタンを押す
1.  差分を確認し、コメント欄を埋めて提出
1.  修正を求められたらそのブランチで変更し、自分のリポジトリに`push`すればPRにも反映される
1.  マージされたらブランチを消す


## 問題と対処

### 用済みブランチの掃除

PRマージ済みのリモートブランチを消したい。
明示的に消すオプションを付けるか、空ブランチをpushするか:

```sh
git push --delete origin issue-42
git push origin :issue-42
```

別のマシンやGitHub PR画面のボタンから消したりすると、
消したはずのリモートブランチの記録が手元のリポジトリに残ることがある。
まず確認:

```sh
git branch -a
git remote show origin
```

こうしたstaleなブランチを刈り取る方法には二通りある:

```sh
git fetch --prune
git remote prune origin
```

それでも消えないローカルブランチは手動で消す:

```sh
git branch -d issue-666
error: The branch 'issue-666' is not fully merged.
If you are sure you want to delete it, run 'git branch -D issue-666'.
```

マージ済みのブランチじゃないと上記のように怒ってくれる。
消してもいい確信があればオプションを大文字 `-D` にして強制削除。


### detached HEAD からの復帰

submoduleなどをいじってると意図せずdetached HEAD状態になることがある。
その状態でcommitしてしまった変更を`master`に反映したい。

1. `push`しようとして怒られて気付く
   ```
   git push
   fatal: You are not currently on a branch
   git status
   HEAD detached from *******
   ```

1. `master`に戻ると道筋を示してくれる:
   ```
   git checkout master
   Warning: you are leaving 2 comits behind, not connected to
   any of your branches
   If you want to keep them by creating a new branch, this may be a good time
   to do so with:

    git branch <new-branch-name> *******
   ```

1. 言われたとおりbranchを作って`merge`
   ```
   git branch detached *******
   git merge detached
   ```

1. 不要になったbranchを消す
   ```
   git branch -d detached
   ```


### 別のリポジトリをサブディレクトリとして取り込む

[Subtree Merging](https://git-scm.com/book/tr/v2/Git-Tools-Advanced-Merging#_subtree_merge)

オプション `-X subtree=${subdir}` を利用してサブディレクトリに入れると、
全体の `git log` では統合されてるように見えるけど、
各ファイルの履歴は途絶えてしまって
`git log --follow ${subdir}/hello.cpp`
などとしても統合前までは辿れない。
予め全ファイルをサブディレクトリに動かすだけのcommitをしておいて、
ルート同士でmergeすると `--follow` が効く状態で取り込める。

```sh
cd /path/to/${subrepo}/
mkdir ${subdir}
git mv $(git ls-tree --name-only master) ${subdir}/
git commit -m ":construction: Move all to ${subdir}/ for integration"

cd /path/to/${mainrepo}/
git remote add ${subrepo} /path/to/${subrepo}
git fetch ${subrepo}
git merge --no-commit --allow-unrelated-histories ${subrepo}/master
git commit
```

- fetchせずにmergeしようとするとブランチ情報が無くて怒られる:
  `not something we can merge`

- 異なる起源をもつリポジトリのmergeは危険なので
  `--allow-unrelated-histories` を明示しないと拒否される:
  `fatal: refusing to merge unrelated histories`
