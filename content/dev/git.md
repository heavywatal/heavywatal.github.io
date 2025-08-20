+++
title = 'Git'
subtitle = '分散型バージョン管理システム'
tags = ["vcs", "writing"]
[menu.main]
  parent = "dev"
+++

[GitHub]: https://github.com

https://git-scm.com/

<iconify-icon inline icon="bi:git"></iconify-icon>
Gitは分散型バージョン管理システムの代表格。
プログラムのソースコードはもちろんのこと、
研究ノートや論文の原稿などあらゆるテキストの管理に使える。

<iconify-icon inline icon="bi:github"></iconify-icon>
[GitHub]はGitをより便利に使うためのオンラインサービス。
個人的なリポジトリ置き場としてはもちろんのこと、
ほかの人と共有・協力してプロジェクトを進めるプラットフォームとしても使える。

Gitのライバルとして[Mercurial](https://www.mercurial-scm.org/)もあるが、
[BitBucket](https://bitbucket.org) (GitHubのライバル)
がGit対応した今となってはMercurialを積極的に使う理由は無い気がする。


## 基本

### 準備

- ローカルマシンにGitをインストールする。
  最初から入ってる場合も多いけど、それが古すぎる場合は
  [Homebrew]({{< relref "homebrew.md" >}}) で新しいのを入れる。
- [GitHub]に個人アカウントを作る。
- [SSH公開鍵を作って]({{< relref "ssh.md" >}})マシンとGitHubに登録する。
  <https://docs.github.com/authentication/connecting-to-github-with-ssh>
- `~/.gitconfig` にユーザ名やアドレスを登録する。
  <https://git-scm.com/docs/git-config>

    ```sh
    git config --global user.name "Watal M. Iwasaki"
    git config --global user.email "heavywatalあmail.com"
    less ~/.gitconfig
    ```
- ついでに `pushinsteadof` の設定をしておく。
  httpsで高速にclone/fetch/pullして、
  sshでパスワード無しでpushする、というのが楽ちん。
- 設定例: [`heavywatal/dotfiles/.gitconfig`](https://github.com/heavywatal/dotfiles/blob/master/.gitconfig)


### 手元の変更を外に伝える

📁 working directory (working tree)
: 手元のファイルの変更はまだリポジトリに登録されていない

&nbsp; ↓ `git add`

<iconify-icon inline icon="bi:git"></iconify-icon> staging area (index)
: 次のコミットに含めるファイルをマークする段階

&nbsp; ↓ `git commit`

<iconify-icon inline icon="bi:git"></iconify-icon> local repository
: 変更履歴が `.git/` 内に記録されている

&nbsp; ↓ `git push`

<iconify-icon inline icon="bi:github"></iconify-icon> remote repository
: GitHubなど別マシンのリポジトリに反映


### 外部の変更を手元に取り込む

<iconify-icon inline icon="bi:github"></iconify-icon> remote repository
: 最新版

&nbsp; ↓ `git fetch`

<iconify-icon inline icon="bi:git"></iconify-icon> local repository
: 変更が `.git/` に取り込まれたが、見えてるファイルには反映されてない

&nbsp; ↓ `git merge` or `git rebase`

📁 working directory
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
: 動詞としては、staging areaの情報をひとつのcommitとしてリポジトリに登録することを指す。

repository
: commitの履歴を保持する拠点。
  `git init` で手元に新規作成するか、`git clone` でリモートから複製する。

`origin`
: remoteリポジトリの典型的なshortname。
  `clone`時に自動的に追加され、
  `push`先や`fetch`元を省略したときにデフォルトで使われる。
  `git remote -v` で確認。
: 他の人のリポジトリをforkして使う場合、
  自分のを `origin`, 元のを `upstream` と名付けるのが一般的。

`master`, `main`
: デフォルトのブランチの典型的な名前。

`HEAD`, `@`
: 現在参照しているbranch/commitを指すポインタ。
  基本的には`main`の最新commitを指していることが多い。
  1つ前は `HEAD^` か `HEAD~`、
  2つ前は `HEAD^^` か `HEAD~~` か `HEAD~2`。
  (`HEAD^2` は `merge` で複数の親がある場合の2番目)

zshの`EXTENDED_GLOB`が有効になってる場合は
`HEAD^` がパターン扱いされてエラーになるので、
`HEAD\^` のようにエスケープするか `unsetopt NOMATCH` しておいたほうがいい。


## よく使うコマンド

### reset

`git reset <DESTINATION>` は `HEAD` の位置を戻す処理で、
オプションによってindexとworking treeもそこに合わせるように変更される。
直前の動作を取り消す用途に絞って使うのが無難。

- `--soft`: `HEAD` 移動のみ。indexとworkingはそのまま。
  ```sh
  # 直前のcommitを取り消し、staged状態まで戻す
  git reset --soft HEAD~
  ```
  直前のcommitをちょっと修正したいだけなら `git commit --amend` が簡単。
  それより前のを修正するには `git rebase -i HEAD~3` とかで戻ってrewordやedit。
- `--mixed`: 移動した `HEAD` にindexも合わせる。
  `git add` を取り消すのに使える:
  ```sh
  # 同義だけど git status で提案されるのは後者
  git reset --mixed HEAD~
  git restore --staged .
  ```
- `--hard`: 移動した `HEAD` にindexとworking treeも合わせる。危険。
  ```sh
  # 変更したファイルをHEADの状態に戻す (DANGEROUS!)
  git reset --hard HEAD

  # reset直後、それを取り消す
  git reset --hard ORIG_HEAD

  # divergedになってしまった手元のbranchを破棄 (DANGEROUS!)
  git reset --hard origin/main
  ```


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

[delta](https://dandavison.github.io/delta/)を使うと見やすくなる。


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
lightweightとannotatedの2種類が存在し、メッセージなどが紐付く後者が推奨。
```sh
git tag v0.1.0 [revision]
git tag -a v0.1.0 -m "Message!"
```

GitHubリポジトリに `git push --tags` するとアーカイブが作られ、
Releasesページに反映される。
annotated tagであれば `git push --follow-tags`
でcommitとtagを同時にpushできる。

https://git-scm.com/book/en/Git-Basics-Tagging


### rebase

ブランチの根本を別のコミットに付け替える。
よく使うのは、開発ブランチを `main` の先頭に追従させるとき。

```sh
git switch some-branch
git rebase main
```

最も近い共通祖先(MRCA)コミットから先を丸ごと移すだけならこのように単純だが、
ブランチのブランチとか、ブランチの一部だけを移したい場合は次のようにする。

```sh
git rebase --onto <newbase> <base> <tip>
```

これで base–tip 間のコミットがnewbaseから伸びる形になる。


### Submodule

- <https://git-scm.com/book/en/v2/Git-Tools-Submodules>
- <https://git-scm.com/docs/git-submodule>

リポジトリの中に別のリポジトリへの参照を保持する仕組み。
個々のファイルの履歴は混ぜずに依存関係を管理したいときに使える。
ソースコードを丸ごとvendorするより冗長性は小さいけど、異物を抱えている感じは否めない。
各分野でパッケージやモジュールを管理する仕組みが発展して不要になることが増えてきた。
e.g., CMake FetchContent, Hugo Modules, etc.

既存のリポジトリを追加する方法は `git clone` と似ている:
```sh
git submodule add https://github.com/rstudio/hex-stickers.git
```

submoduleを含むリポジトリを普通に `git clone` してくると参照しか取れないので、
まず次のようなコマンドで実体を取得する:
```sh
git submodule update --init --recursive
```

submoduleの更新をまとめて取得:
```sh
git submodule foreach git pull
```

submoduleを削除するのは1つのコマンドで完結できない:
```sh
git submodule deinit $SUB_MOD
git rm -r $SUB_MOD
rm -rf .git/modules/$SUB_MOD
```

submoduleに何らかの変更を加えたら親リポジトリでも参照を更新する。


## GitHub Pages

<https://docs.github.com/pages>

### ユーザーサイトを作る

1. `USERNAME.github.io` という名前のリポジトリをGitHub上で作成
1. 公開したいウェブサイトを `main` ブランチとして `push`
1. `https://USERNAME.github.io` にアクセスしてみる。

例えば本ウェブサイトは
[`heavywatal/heavywatal.github.io`](https://github.com/heavywatal/heavywatal.github.io/tree/gh-pages)
というリポジトリの `source` ブランチでMarkdownテキストを書き、
[Hugo]({{< relref "hugo.md" >}})
で変換・生成したHTMLファイルを `gh-pages` ブランチに書き出している。

GitHubが勝手にJekyll処理しようとすることがあるので、
`.nojekyll` という空ファイルを作っておく。


### プロジェクトサイトを作る

<https://docs.github.com/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site>

リポジトリの内容を `https://USERNAME.github.io/PROJECT/` に公開することができる。
方法は複数あり、リポジトリの設定画面から選択できる。

- GitHub Actions を使ってビルドからデプロイまで自動化
- 特定のブランチ、特定のディレクトリを公開。e.g.,
  - `gh-pages` ブランチの `/` ルートを公開
  - `main` ブランチの `/docs` ディレクトリのみを公開

ソースコードやプログラムを更新した結果、
生成されるドキュメントがどう変わるかを一応確認してからデプロイしたいので、
今のところ手でビルド→確認→commit→pushという流れを採用している。

`main` ブランチだけで済ませる方法は楽な面もあるけど、
コミット履歴が汚くなるというデメリットがある。
結局 `gh-pages` ブランチを作る古来の方法がいい塩梅。

その場合 `main` の作業ディレクトリで `git switch --orphan gh-pages` するより、
新しいディレクトリから初期化するほうが簡単かつ安全:
```sh
git init -b gh-pages docs
cd docs/
git remote add origin https://github.com/USERNAME/PROJECT.git
git commit --allow-empty -m ":beer: Create orphan branch"
git push -u origin gh-pages
```


## Pull Request (PR)

-   <https://docs.github.com/en/pull-requests>
-   他人のプロジェクトに変更を提案するための仕組み。
-   リモートリポジトリを2つ参照することになる。
    - 大元のを`upstream`、フォークした自分のを`origin`と呼ぶのが一般的。
    - 大元のを`origin`、フォークを`自分のユーザー名`にすることもある。
      e.g., [Homebrew PR](https://docs.brew.sh/How-To-Open-a-Homebrew-Pull-Request)
-   デフォルトブランチ(`main`とか`develop`とか)は更新取得のためだけに使い、変更は新規ブランチで行う。


### 基本の流れ

例えば `USER` さんの `PROJECT` のコード修正に貢献する場合。

1.  `github.com/USER/PROJECT` の "Fork" ボタンで自分のGitHubリポジトリに取り込む
1.  forkした自分のリポジトリからローカルに`clone`:
    ```sh
    git clone https://github.com/heavywatal/PROJECT.git
    cd PROJECT/
    ```

1.  大元のリポジトリに`upstream`という名前をつけておく:
    ```sh
    git remote add upstream https://github.com/USER/PROJECT.git
    ```

1.  PR用のブランチを切って移動:
    ```sh
    git switch -c fix-typo
    ```

1.  コードを変更して`commit`:
    ```sh
    vim README.md
    git diff
    git commit -a -m "Fix typo in README.md"
    ```

1.  この間に`upstream`で更新があったか確認:
    ```sh
    git fetch upstream
    ```

    必要ならそれを取り込む:
    ```sh
    git switch main
    git merge --ff-only upstream/main
    git switch fix-typo
    git rebase main
    ```

1.  自分のリポジトリに`push`:
    ```sh
    git push -u origin fix-typo
    ```

1.  PR用のURLが表示されるのでそこから飛ぶ。
    もしくはGitHub上に出現する "Compare & pull request" ボタンを押す。
1.  差分を確認し、コメント欄を埋めて提出。
1.  待つ。

    修正を求められたらそのブランチで変更し、自分のリポジトリに`push`すればPRにも反映される。

    この間に `upstream` に更新があった場合、
    `rebase` でPRに取り込むべきか、マージする側でその処理をするか、
    判断は元のプロジェクト管理者に従う。
    `rebase` すると `git push --force` が必要になる。
1.  マージされたらブランチを消す。


### 用済みブランチの掃除

PRマージ済みのリモートブランチを消したい。
一番簡単なのは、GitHub PR画面のdelete branchボタンを押すこと。
手元のgitからやるには、明示的に `--delete` するか、空ブランチをpushするか:

```sh
git push --delete origin issue-42
git push origin :issue-42
```

リモートブランチを消しても手元のリポジトリには残る。
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
```
```
error: The branch 'issue-666' is not fully merged.
If you are sure you want to delete it, run 'git branch -D issue-666'.
```

マージ済みのブランチじゃないと上記のように怒ってくれる。
消してもいい確信があればオプションを大文字 `-D` にして強制削除。

`--single-branch` をつけ忘れて `git clone` したあと
`git branch --remote --delete` でリモートブランチへの参照を消しても、
次の `fetch`/`pull` でまた復活してしまう。
そうならないように消すには:
```sh
git remote set-branches origin main
git fetch --prune
```


## 問題と対処

`git status` で次の一手を提案してくれることが多いので、まずはそれに従う。


### Trailing whitespace 以外の変更だけ add する

大概のIDEには保存時に行末の空白を自動削除するオプションがある。
それによって自分のソースコードは常にきれいに保てるが、
他人の汚いコードや [knitr]({{< relref "knitr.md" >}}) の結果などを編集するときに余計な差分を作ってしまう。
[VSCode]({{< relref "vscode.md" >}}) なら "Save without Formatting"
で設定を変えずに済ませられることは覚えていても、
"Find in Files" で一括編集したときにも空白が削られることは忘れがち。

```sh
git diff --ignore-space-at-eol | git apply --cached
```

上記ワンライナーで大概うまくいくが、
変更箇所が近かったりすると `error: patch failed` と蹴られる。
その場合は次のようにworkaround:

```sh
git diff --ignore-space-at-eol > tmp.diff
git stash
git apply --cached tmp.diff
git stash drop
```


### detached HEAD からの復帰

submoduleなどをいじってると意図せずdetached HEAD状態になることがある。
その状態でcommitしてしまった変更を`main`に反映したい。

1. `push`しようとして怒られて気付く
   ```sh
   git push
   ```
   ```
   fatal: You are not currently on a branch
   ```
   ```sh
   git status
   ```
   ```
   HEAD detached from *******
   ```

1. `main`に戻ると道筋を示してくれる:
   ```sh
   git switch main
   ```
   ```
   Warning: you are leaving 2 comits behind, not connected to
   any of your branches
   If you want to keep them by creating a new branch, this may be a good time
   to do so with:

    git branch <new-branch-name> *******
   ```

1. 言われたとおりbranchを作って`merge`
   ```sh
   git branch detached *******
   git merge detached
   ```

1. 不要になったbranchを消す
   ```sh
   git branch -d detached
   ```

### サブディレクトリを別のブランチやリポジトリに切り分ける

[`git subtree`](https://github.com/git/git/tree/master/contrib/subtree)
は公式サブコマンドとして取り込まれてはいるけどドキュメントには入ってない。
ドキュメントのある
[`git filter-branch`](https://git-scm.com/docs/git-filter-branch)は非推奨になり、
今のところサードパーティの[`git-filter-repo`](https://github.com/newren/git-filter-repo)が推奨。

> [!CAUTION]
>
> [`man git-filter-repo`](https://htmlpreview.github.io/?https://github.com/newren/git-filter-repo/blob/docs/html/git-filter-repo.html)
を読み、歴史改変がいかに危険で面倒なことかを理解する。
>
> 例えば `main` ブランチの `docs/` で公開していた GitHub Pages を
`gh-pages` ブランチに切り分け、 `main` の履歴から消去したいとする。
`main` の履歴の検索性が上がるなどのメリットがあるとして、
散在する全てのクローンやフォークが一旦無効になるリスクとコストに見合うかどうか。
`gh-pages` ブランチに置くのはソースコードではなく生成物なので、
履歴の重要性はそれほど大きくない。
自分ひとりで開発している小規模リポジトリならギリギリ試す価値ありかもしれないけど、
`main` はそのままにして `gh-pages` は新しいorphanブランチを作るほうが無難。

インストール方法はいろいろあるので好きなやつを選ぶ。
Python製なのでuvがセットアップしてあれば早い:
```sh
uv pip install git-filter-repo
```

安全のため手元のリポジトリを使わず、
新たに `git clone` したフレッシュなリポジトリを使う。

`docs/` 以下の履歴だけを保持し、トップレベルに移動する例:
```sh
git filter-repo --path docs/ --path-rename docs/:
ls
git log --graph --all
git remote -v
```

元のリポジトリを誤って上書きしないようにリモートの設定が消去される。
必要に応じて `git remote add origin` で登録し直す。


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
git mv $(git ls-tree --name-only main) ${subdir}/
git commit -m ":construction: Move all to ${subdir}/ for integration"

cd /path/to/${mainrepo}/
git remote add ${subrepo} /path/to/${subrepo}
git fetch ${subrepo}
git merge --no-commit --allow-unrelated-histories ${subrepo}/main
git commit
```

- fetchせずにmergeしようとするとブランチ情報が無くて怒られる:
  `not something we can merge`

- 異なる起源をもつリポジトリのmergeは危険なので
  `--allow-unrelated-histories` を明示しないと拒否される:
  `fatal: refusing to merge unrelated histories`
