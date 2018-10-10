+++
title = "Git入門2018"
date = 2018-10-15T13:30:00+09:00
draft = false
toc = true
tags = ["vcs", "writing"]
[menu.main]
  parent = "lectures"
+++

<style>
.fa-git-square, .fa-share-alt-square {
  transform: rotate(45deg);
}
.fa-git, .fa-git-square, .fa-share-alt-square {
  color: #f03c2e;
}
.fa-github {
  color: #000;
}
.fa-bitbucket {
  color: #2e87fb;
}
.fa-gitlab {
  color: #e04432;
}
.fa-tint {
  color: #949494;
}
.fa-coffee {
  color: #62982f;
}
</style>

2018-10-15 中央水産研究所？

## 事前準備

### OS別

新しめのUNIX系OSが好ましい。

- macOS
    -   10.13 High Sierra以降
    -   Command Line Tools
    -   Gitは標準装備のやつで大丈夫
- Linux
    -   Ubuntu 16.04以降 or CentOS 7以降
    -   build-essential とか "Development Tools" 的なもの
    -   普通にaptやyumで入るGitが古いv1系だったら、どうにかしてv2以降を入れる
- Windows
    -   10 (1803以降)
    -   [Git本体およびGit Bashを公式から入れる](https://git-scm.com/download/)。
        インストール途中でいろいろ選択肢あるけどとりあえずデフォルトでよいのでは。
    -   [RStudio](https://www.rstudio.com/)のターミナルのシェルとしてGit Bashを設定。
        "Tools > Global Options > Terminal > Shell"
    -   RStudioで新しいTerminalを立ち上げて最低限のコマンド操作に慣れる。
        e.g.,  `pwd`, `ls`, `cd`, `mkdir`
    -   Git Bashにおける `HOME` がRStudioにおけるそれと異なることに注意。

        ```sh
        ### Git Bash
        echo $HOME
        # /c/Users/watal
        echo ~
        # /c/Users/watal
        ```
        ```r
        ### R
        normalizePath("~")
        # "C:\\Users\\watal\\Documents"
        ```

### OS共通

-   適当なテキストエディタ(開発環境)を入れておく。
    Git/GitHubとの連携機能が付いていて、
    変更箇所を色付けしてくれたりコマンド入力を肩代わりしてくれたりするのが便利。
    - [Atom](https://atom.io/): GitHub製
    - [VS Code](https://code.visualstudio.com/): Microsoft製
    - [RStudio](https://www.rstudio.com/): RStudio製

-   [GitHub<i class="fab fa-fw fa-github"></i>](https://github.com)に個人アカウントを作る。

-   Git<i class="fas fa-fw fa-share-alt-square"></i>の初期設定をターミナルから行う:

    ```sh
    git --version  # 2.0以上であることを確認
    git config --global user.name "Watal M. Iwasaki"
    git config --global user.email "heavywatalあmail.com"
    git config --global push.default simple
    cat ~/.gitconfig
    ```
    ```ini
    [user]
      name = Watal M. Iwasaki
      email = heavywatalあmail.com
    [push]
      default = simple
    ```


### SSHの設定 (任意)

-   GitHubとの通信に2つの方式がある。
    - HTTPS: 設定不要で高速だが、操作によってパスワード入力が必要
    - SSH: 一旦ちゃんと設定すればパスワードなしで快適
-   ダウンロード操作(clone/fetch/pull)は高速なHTTPSで、<br>
    アップロード操作(push)はパスワード無しのSSHで、というのが楽ちん。
-   [SSH公開鍵を作って]({{< relref "ssh.md" >}})ローカルマシンとGitHubに登録する。
-   設定ファイル `~/.gitconfig` に `pushinsteadof` の設定を追加:

    ```ini
    [url "git@github.com:"]
      pushinsteadof = https://github.com/
    ```

----

以下、本編。

🚧 UNDER CONSTRUCTION 🚧


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

- いつまでも履歴を保持してもらえるとは限らない。<br>
  (ここでDropboxの履歴を例示しようと思ったらエラーで使えなかった)
- オフラインだったりバッテリー駆動だったりすると保存漏れが起きる。
- いつのバージョンに戻したらいいのか、日時以外の手掛かりが無い。
- 共有ファイルにおける更新の衝突に対処しにくい。


## Git<i class="fas fa-fw fa-share-alt-square"></i> and GitHub<i class="fab fa-fw fa-github"></i>

- 履歴を残すタイミングは任意 = 手動。
- オフラインでも作業できる。
- バージョン(リビジョン)ごとにメッセージを残せる。
- 差分を簡単に見られる。
- 共同作業のための機能がちゃんとある。

e.g., https://github.com/tidyverse/stringr/commits/master

### それらはどういう関係？

[Git<i class="fas fa-fw fa-share-alt-square"></i>](https://git-scm.com/)は分散型バージョン管理システムの代表格。
プログラムのソースコードはもちろんのこと、
研究ノートや論文の原稿などあらゆるテキストの管理に使える。
各自のコンピュータにインストールして使うのはこっち。

[GitHub<i class="fab fa-fw fa-github"></i>](https://github.com)はGitをより便利に使うためのオンラインサービス。
個人的なリポジトリ置き場としてはもちろんのこと、
ほかの人と共有・協力してプロジェクトを進めるプラットフォームとしても使える。

### 類似ツール

- Version Control System (VCS)
    - [Git<i class="fas fa-fw fa-share-alt-square"></i> `git`](https://git-scm.com/)
    - [Mercurial<i class="fas fa-fw fa-tint"></i> `hg`](https://www.mercurial-scm.org/)
    - その他 svn, cvs, rcs など。
- Hosting Service
    - [GitHub<i class="fab fa-fw fa-github"></i>](https://github.com):
      公開リポジトリは無料。情報も連携も豊富。
    - [Bitbucket<i class="fab fa-fw fa-bitbucket"></i>](https://bitbucket.org/):
      非公開リポジトリも無料。
    - [GitLab<i class="fab fa-fw fa-gitlab"></i>](https://about.gitlab.com/):
      非公開リポジトリも無料。ローカル版もあり。
    - [Gitea<i class="fas fa-fw fa-coffee"></i>](https://gitea.io/en-us/):
      ローカル版のみ。
    - その他 SourceForge, Google Code など。

VCSは基本的にGit一択。<br>
ホスティングサービスは、使い方や予算に応じて選択。


### GitHub<i class="fab fa-fw fa-github"></i>の使いみち

- 基本: ファイルのバージョン管理<i class="fas fa-fw fa-share-alt-square"></i>
    - プログラムのソースコード:
      e.g., [ggplot2](https://github.com/tidyverse/ggplot2), [rstan](https://github.com/stan-dev/rstan)
    - 論文や本の原稿、サプリ:
      e.g., [R4DS](https://github.com/hadley/r4ds), [Advanced R programming](https://github.com/hadley/adv-r/)

- Issues:
  バグ報告や機能要望に用いられる。
  チームの課題を列挙するのにも使える。<br>
  e.g., https://github.com/tidyverse/ggplot2/issues

- Projects:
  プロジェクトのタスク管理のためのツール。
  もちろんissueとも連携可能。<br>
  e.g., https://github.com/r-lib/pillar/projects/1

- Wiki:
  チーム内のちょっとした情報共有などに。
  でもできればそういう文書もちゃんとGitで管理したほうがいい。<br>
  e.g., https://github.com/gnab/remark/wiki

- GitHub Pages:
  リポジトリの内容をウェブサイトとして公開できる。<br>
  e.g., https://kazutan.github.io/kazutanR/


## 構造と操作の概要

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

https://help.github.com/articles/github-glossary/

repository
: commitの履歴を保持する拠点。
  「ひとつのRパッケージ」とか「1冊の本の原稿」のような単位で作る。
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


## 基本操作を実践

### 既存のリポジトリを取ってくる `clone`

1.  GitHub上の適当なリポジトリをひとつ選ぶ。
    (e.g., <https://github.com/heavywatal/clippson>)
1.  右の方の緑の "Clone or download" ボタンを押す。
1.  SSHではなくHTTPSを選択し、URLをコピー。
1.  ターミナルにコマンドを入力:<br>
    `git clone https://github.com/heavywatal/clippson.git`
1.  中身を眺めてみる:

    ```sh
    cd clippson/
    ls -al
    git log
    git remove -v
    ```

### 新しいリポジトリを作る `init`

1.  GitHubの右上の "+" から "New repository" を選択。
1.  Repository name を例えば `helloworld` として "Create repository" を押す。
    いくつかのファイル (`README.md`, `LICENSE`, `.gitignore`)
    をここで作ることもできるけど、今回はとりあえず空っぽのリポジトリを作る。
1.  手元のマシンにローカルリポジトリを作る:

    ```sh
    mkdir helloworld
    cd helloworld/
    git init
    ls -al
    ```

    リポジトリの本体 `.git/` が作成されたことを確認。

1.  空っぽのコミットを作る:

    ```sh
    git status
    git commit --allow-empty -m ":beer: Create repository"
    git status
    git log
    ```

    事あるごとに `git status` や `git log` を確認すると安心。

1.  先程作ったリモートリポジトリを紐付けて、pushしてみる:

    ```sh
    git remote -v
    git remote add origin https://github.com/YOUR_NAME/helloworld.git
    git remote -v
    git push -u origin master
    git status
    ```

1.  GitHubで履歴を閲覧し、 `git log` と同じになってることを確認。


### 手元の変更をリモートに `push`

1.  上で作ったリポジトリに、適当なファイルを追加:

    ```sh
    echo "# Hello, world!" > README.md
    cat README.md
    git status
    ```

1.  作ったファイルをstaging areaに追加:

    ```sh
    git add README.md
    git status
    git diff --staged
    ```

1.  この変更をcommit:

    ```sh
    git commit -m ":memo: Create README.md"
    git status
    git log
    git show
    ```

1.  リモートにpush:

    ```sh
    git push
    git status
    git log
    ```

### リモートの変更を手元に `fetch`

1.  上のリポジトリでそのまま `git fetch` してみる。
    ローカルとリモートは同じ状態なので当然何も起こらない。

1.  練習のためGitHub上で `LICENSE` ファイルを作成する。
    1.  GitHub上のリポジトリのトップページを開き
        "Create new file" ボタンを押す。
    1.  ファイル名に `LICENSE` と入力。
    1.  右に現れる "Choose a license template" というボタンを押す。
    1.  とりあえず "MIT License" を選択。
    1.  YearとNameを適当に埋めて "Review and submit"。
    1.  "Commit directly to the `master` branch" を選択して "Commit new file"

1.  その変更をローカルリポジトリに取り寄せる:

    ```sh
    git fetch
    git status
    git log --all
    ls -al
    ```
    リポジトリ内部 `.git/` の `origin/master` は更新されたが、
    working directoryにはまだ反映されていない。

1.  `origin/master` の内容を手元のファイルに反映する:

    ```sh
    git merge
    git status
    git log
    git show
    ls -al
    ```

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

push済みのものは改変しちゃダメ！


## チーム作業

- 基本的に、自分のリモートリポジトリにpushできるのは自分だけ。
- コラボレータを設定して、権限を与えることも可能。
  ("Settings > Collaborators")
- でも権限を持つ人が増えすぎると更新の衝突など管理リスクも増える。
- 権限を持たない人はforkからPull Requestを送り、
  権限を持つ少数の人がそれをレビューしてmergeするスタイルが安全。


### Pull Request (PR)

他人のリポジトリに貢献するためのGitHubの機能。

e.g., https://github.com/Rdatatable/data.table/pull/2807

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

1.  PLAYER: `README.md` をテキストエディタで変更して `commit`:

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

## Further reading

- [Pro Git book](https://git-scm.com/book/en/v2)
- [GitHub Learning Lab](https://lab.github.com/)
