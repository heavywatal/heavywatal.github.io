+++
title = "Git入門2018"
date = 2018-10-15T13:30:00+09:00
draft = false
tags = ["vcs", "writing"]
[menu.main]
  parent = "lectures"
+++

2018-10-15 中央水産研究所？


## Gitが無い世界の風景

バージョン違いのファイルがフォルダを埋め尽くす。

```
% ls
analysis.R
analysis2.R
analysis-20180129.R
analysis-20180129改良版.R
analysis-20180210.R
analysis-20180210バグ？.R
analysis佐藤edit.R
analysis佐藤edit田中.R
analysis完全版.R
analysis最終.R
analysis最終改.R
analysis決定版！.R
analysis真・最終.R
plot.R
plot2.R
plot最終.R
plot論文.R
```

## オンラインストレージやバックアップ機能では不十分

- [Dropbox<i class="fab fa-fw fa-dropbox"></i>](https://dropbox.com) とか
  [Google Drive<i class="fab fa-fw fa-google-drive"></i>](https://drive.google.com/)
  では、保存のたびに履歴が残る。
- Time Machine<i class="fas fa-fw fa-clock"></i>では、一定時間間隔で履歴が残る。

でも、バージョン管理や共同作業のためのツールじゃないから...

- いつまでも履歴を保持してもらえるとは限らない。
- いつのバージョンに戻したらいいのか、日時以外の手掛かりが無い。
- 共有ファイルにおける更新の衝突に対処しにくい。


## Git+GitHubを使うのが楽ちん

例: https://github.com/heavywatal/clippson

- 履歴を残すタイミングは任意 = 手動。
- バージョン(リビジョン)ごとにメッセージを残せる。
- 差分を簡単に見られる。
- 共同作業のための機能がちゃんとある。


## Git and GitHub<i class="fab fa-fw fa-github"></i>

[Git](https://git-scm.com/)は分散型バージョン管理システムの代表格。
プログラムのソースコードはもちろんのこと、
研究ノートや論文の原稿などあらゆるテキストの管理に使える。
各自のコンピュータにインストールして使うのはこっち。

[GitHub<i class="fab fa-fw fa-github"></i>](https://github.com)はGitをより便利に使うためのオンラインサービス。
個人的なリポジトリ置き場としてはもちろんのこと、
ほかの人と共有・協力してプロジェクトを進めるプラットフォームとしても使える。


## 類似ツール

- Version Control System (VCS)
    - [Git `git`](https://git-scm.com/)
    - [Mercurial `hg`](https://www.mercurial-scm.org/)
    - その他 svn, cvs, rcs など。
- Hosting Service
    - [GitHub<i class="fab fa-fw fa-github"></i>](https://github.com):
      公開リポジトリは無料。情報も連携も豊富。
    - [Bitbucket<i class="fab fa-fw fa-bitbucket"></i>](https://bitbucket.org/):
      非公開リポジトリも無料。
    - [GitLab<i class="fab fa-fw fa-gitlab"></i>](https://about.gitlab.com/):
      非公開リポジトリも無料。ローカル版もあり。
    - [Gitea](https://gitea.io/en-us/):
      ローカル版のみ。
    - その他 SourceForge, Google Code など。

VCSは基本的にGit一択。<br>
ホスティングサービスは、使い方や予算に応じて選択。

## 基本的な構造と操作

### 手元の変更を外に伝える

<i class="fas fa-fw fa-folder"></i> working directory (working tree)
: 手元のファイルの変更はまだリポジトリに登録されていない
: ↓ `git add`

<i class="fas fa-fw fa-check"></i> staging area (index)
: 次のコミットに含めるファイルをマークする段階
: ↓ `git commit`

<i class="fas fa-fw fa-code-branch"></i> local repository
: 変更履歴が `.git/` 内に記録されている
: ↓ `git push`

<i class="fab fa-fw fa-github"></i> remote repository
: GitHubなど別マシンのリポジトリに反映


### 外部の変更を手元に取り込む

<i class="fab fa-fw fa-github"></i> remote repository
: ↓ `git fetch`

<i class="fas fa-fw fa-code-branch"></i> local repository
: 変更が `.git/` に取り込まれたが、見えてるファイルには反映されてない
: ↓ `git checkout` or `git merge`

<i class="fas fa-fw fa-folder"></i> working directory
: 手元のファイルが最新版に同期されている


## 用語

repository
: commitの履歴を保持する拠点。
  `git init` で手元に新規作成するか、`git clone` でリモートから複製する。

commit
: git内部でroot treeのsnapshotを指すオブジェクト。
  root treeのハッシュID、著者、コメントなどの情報を持つ。
  動詞としては、staging areaの情報をひとつのcommitとしてリポジトリに登録することを指す。

tree
: git内部で1つのディレクトリを指すオブジェクトで、`commit`した時に作られる。
  blobやファイル名などのメタデータに依存したハッシュIDを持ち、
  その変化は親に伝播する。

blob
: git内部で1つのファイルを指すオブジェクトで、`add`時に作られる。
  ファイル名などのメタデータは持たず、
  ファイルの内容にのみ依存したハッシュIDを持つ。

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


## 下準備

- ローカルマシンにGitをインストールする。
  MacならCommand Line Toolsに付属。
  Linuxなら `apt` や `yum` で入れられる。
  [Homebrew](https://brew.sh/)/[Linux](http://linuxbrew.sh/)などで最新版を入れるのも良い。
- [GitHub<i class="fab fa-fw fa-github"></i>](https://github.com)に個人アカウントを作る。
- [SSH公開鍵を作って](/dev/ssh.html)マシンとGitHubに登録する。
  https://help.github.com/articles/connecting-to-github-with-ssh/
- `~/.gitconfig` にユーザ名やアドレスを登録する。
  https://git-scm.com/docs/git-config

    ```sh
    git config --global user.name "Watal M. Iwasaki"
    git config --global user.email "heavy.watalあmail.com"
    cat ~/.gitconfig
    ```
- ついでに `pushinsteadof` の設定をしておく。
  httpsで高速にclone/fetch/pullして、
  sshでパスワード無しでpushする、というのが楽ちん。
- 設定例: [`heavywatal/dotfiles/.gitconfig`](https://github.com/heavywatal/dotfiles/blob/master/.gitconfig)

## Gitの基本

### 新しいリポジトリを作る

1.  GitHubの右上の "+" から "New repository" を選択。
1.  適当に空欄を埋めて "Create repository" を押す。
    いくつかのファイル (`README.md`, `LICENSE`, `.gitignore`)
    をここで作ることもできるけど、今回はとりあえず空っぽのリポジトリを作る。
1.  手元のマシンにローカルリポジトリを作る:

    ```sh
    mkdir YOUR_REPO
    cd YOUR_REPO/
    git init
    git commit --allow-empty -m ":beer: Create repository"
    ```

1.  先程作ったリモートリポジトリを紐付けて、プッシュしてみる:

    ```sh
    git remote add origin https://github.com/YOUR_NAME/YOUR_REPO.git
    git push -u origin master
    ```

### 既存のリポジトリを取ってくる


### 手元の変更をリモートに


### リモートの変更を手元に


## その他よく使うコマンド

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


## チーム作業

- 基本的に、自分のリモートリポジトリにpushできるのは自分だけ。
- コラボレータを設定して、権限を与えることも可能。
  ("Settings > Collaborators")
- でも権限を持つ人が増えすぎると更新の衝突など管理リスクも増える。
- 権限を持たない人はforkからPull Requestを送り、
  権限を持つ少数の人がそれをレビューしてmergeするスタイルが安全。


## Pull Request (PR)

他人のリポジトリに貢献するためのGitHubの機能。

例: https://github.com/Rdatatable/data.table/pull/2807

### 2人1組でPRとmergeを体験

- MANAGER: リポジトリの管理権限を持つ人
- PLAYER: 権限を持たず、PRを送る人

(できれば横に並んで相手の画面も見えるように)

1.  MANAGER: GitHubで新しいリポジトリを作成
1.  MANAGER: 何かtypoを含む `README.md` を作ってpush
1.  PLAYER: 相手のGitHubリポジトリでその `README.md` が見えることを確認
1.  PLAYER: 右上のForkボタンで自分のGitHubリポジトリに取り込む
1.  PLAYER: forkした自分のリポジトリからローカルに`clone`:

        git clone https://github.com/{PLAYER}/PROJECT.git
        cd PROJECT/

1.  PLAYER: 大元のリポジトリに`upstream`という名前をつけておく:

        git remote add upstream https://github.com/{MANAGER}/PROJECT.git
        git remote -v

    ちなみに自分のリポジトリには自動的に `origin` という名前がついている。

1.  PLAYER: PR用のブランチを切って移動:

        git checkout -b fix-typo

1.  PLAYER: コードを変更して `commit`:

        emacs README.md
        git diff
        git commit -a -m "Fix typo in README.md"

1.  PLAYER: この間に`upstream`で更新が無いかどうか確認:

        git fetch upstream

    もしあったら、それをデフォルトブランチ(`master`)越しに取り込む:

        git checkout master
        git merge upstream/master
        git push origin/master
        git checkout fix-typo
        git rebase -i master

1.  PLAYER: 自分のリポジトリに`push`:

        git push origin fix-typo

1.  PLAYER: GitHub上に出現する "Compare & pull request" ボタンを押す。
1.  PLAYER: 差分を確認し、コメント欄を埋めて提出。
1.  MANAGER: 受け取ったPRを確認。必要に応じて修正を要求したり、自分で修正したり。
1.  PLAYER: 修正を求められたらそのブランチに続けてcommitしてまたpush。
1.  MANAGER: 問題が無ければmergeする。
1.  PLAYER: 無事マージされたら作業ブランチを消す。

## GitHubのその他の便利機能

- Issues:
  バグ報告に用いられる。
  チームの課題を列挙するのにも使える。
  例: https://github.com/tidyverse/ggplot2/issues

- Projects:
  プロジェクトのタスク管理のためのツール。
  もちろんissueとも連携可能。
  例: https://github.com/r-lib/pillar/projects/1

- Wiki:
  チーム内のちょっとした情報共有などに。
  でもできればそういう文書もちゃんとGitで管理したほうがいい。
  例: https://github.com/gnab/remark/wiki


## GitHub Pages

https://help.github.com/articles/user-organization-and-project-pages/

### ユーザーサイトを作る

1. `USERNAME.github.io` という名前のリポジトリをGitHub上で作成
2. 公開したいウェブサイトを`master`ブランチとして`push`
3. `https://USERNAME.github.io` にアクセスしてみる。

例えば本ウェブサイトは
`heavywatal.github.io` というリポジトリの
`source`ブランチでMarkdownテキストを書き、
[Hugo](https://gohugo.io/)
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


## 開発環境

最近のテキストエディタはGit+GitHubとの連携機能が付いてたりする。
これらを使えば、ターミナルからコマンドを打たなくても済む。

- [Atom](https://atom.io/): GitHub製
- [VS Code](https://code.visualstudio.com/): Microsoft製
- [RStudio](https://www.rstudio.com/): RStudio製
