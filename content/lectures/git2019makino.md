+++
title = "Git入門2019"
date = 2019-10-30T14:00:00+09:00
draft = false
toc = true
tags = ["vcs", "writing"]
[menu.main]
  parent = "lectures"
+++

2019-10-30 東北大学 生命科学研究科 進化ゲノミクス分野 牧野研

[前半スライド](/slides/makino2019r/5-git.html)

## Environment / 環境

### OS

新しめのUNIX系OSが好ましい。

- macOS
    -   10.14 Mojave以降
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
   -    UNIX派が仕方なくWindowsを使う場合は
        [MSYS2](https://www.msys2.org/)の `pacman` で環境を構築するのが良さそう...?


### Common

-   適当なテキストエディタ(開発環境)を入れておく。
    初期状態でもGit/GitHubとの連携機能が付いていて、
    変更箇所を色付けしてくれたりコマンド入力を肩代わりしてくれたりするのが便利。
    - [Atom](https://atom.io/): GitHub製
    - [VSCode]({{< relref "vscode.md" >}}): Microsoft製
    - [RStudio](https://www.rstudio.com/): RStudio製

-   [GitHub<img height=16 width=16 src="https://cdn.simpleicons.org/github">](https://github.com)に個人アカウントを作る。

-   Git<img height=16 width=16 src="https://cdn.simpleicons.org/git">の初期設定をターミナルから行う:

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

-   Enable macOS Keychain to skip password authentication:

    ```sh
    git config --global credential.helper osxkeychain
    ```


### SSH (任意)

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


## Essential commands / 基本操作

### Fetch existing repositories: `clone`

1.  GitHub上の適当なリポジトリをひとつ選ぶ。
    (e.g., <https://github.com/heavywatal/tumopp>)
1.  右の方の緑の "Clone or download" ボタンを押す。
1.  SSHではなくHTTPSを選択し、URLをコピー。
1.  ターミナルにコマンドを入力:<br>
    `git clone https://github.com/heavywatal/tumopp.git`
1.  中身を眺めてみる:

    ```sh
    cd tumopp/
    ls -al
    git log
    git remote -v
    ```

最新版のスナップショットだけでなく、
履歴もごっそり複製するので、
このあとはオフラインでもいろいろ操作できる。


### Create new repositories: `init`

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


### Export local changes to a remote server: `push`

1.  上で作ったリポジトリに、適当なファイルを追加:

    ```sh
    echo '# Hello, world' > README.md
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

### Import changes from a remote server: `fetch`

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

`git fetch` と `git merge` を一気にやってくれる `git pull` というコマンドもあり、
普段の一人作業ではよく使う。


## Other commands

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

**リモートにpush済みのものは改変しちゃダメ！**


## Collaboration

- 基本的に、自分のリモートリポジトリにpushできるのは自分だけ。
- コラボレータを設定して、権限を与えることも可能。
  ("Settings > Collaborators")
- でも権限を持つ人が増えすぎると競合・衝突など管理リスクも増える。
- 権限を持たない人はforkからPull Requestを送り、
  権限を持つ少数の人がそれをレビューしてmergeするスタイルが安全。


### Pull Request (PR)

他人のリポジトリに貢献するためのGitHubの機能。<br>
e.g., https://github.com/Rdatatable/data.table/pull/2807

1.  貢献したいリポジトリをForkして自分のGitHubアカウントに追加。
1.  Forkした自分のリポジトリを手元にclone。
1.  PR用のブランチを作って、そこでソースコードを編集。
1.  コミットして、ひとまず自分のGitHubアカウントにpush。
1.  GitHub上で大元のリポジトリにPRを送る。
1.  取り込んでもらえたら、用済みのブランチを削除。


### 2人1組でPRとmergeを体験

- 🐸 KING: リポジトリの管理権限を持つ人
- 🐰 PAWN: 権限を持たず、PRを送る人

(できれば横に並んで相手の画面も見えるように)

1.  🐸 GitHubで新しいリポジトリを作成
1.  🐸 何かtypoを含む `README.md` を作ってpush
1.  🐰 相手のGitHubリポジトリでその `README.md` が見えることを確認
1.  🐰 右上のForkボタンで自分のGitHubリポジトリに取り込む
1.  🐰 forkした自分のリポジトリからローカルに`clone`:

        git clone https://github.com/{PAWN}/PROJECT.git
        cd PROJECT/

1.  🐰 大元のリポジトリに`upstream`という名前をつけておく:

        git remote add upstream https://github.com/{KING}/PROJECT.git
        git remote -v

    ちなみに自分のリポジトリには自動的に `origin` という名前がついている。

1.  🐰 PR用のブランチを切って移動:

        git checkout -b fix-typo

1.  🐰 `README.md` をテキストエディタで編集して `commit`:

        git diff
        git commit -a -m ":memo: Fix typo in README.md"

    Git連携機能のあるエディタを使っている場合、
    そこからdiffやcommitをやってみてもよい。
    コードの追加・変更・削除による色分けの便利さも体感しよう。

1.  🐰 この間に`upstream`で更新が無いかどうか確認:

        git fetch upstream

    もしあったら、それをデフォルトブランチ(`master`)越しに取り込む:

        git checkout master
        git merge upstream/master
        git push origin/master
        git checkout fix-typo
        git rebase -i master

1.  🐰 自分のリポジトリに`push`:

        git push origin fix-typo

1.  🐰 GitHub上に出現する "Compare & pull request" ボタンを押す。
1.  🐰 差分を確認し、コメント欄を埋めて提出。
1.  🐸 受け取ったPRを確認。必要に応じて修正を要求したり、自分で修正したり。
1.  🐰 修正を求められたらそのブランチに続けてcommitしてまたpush。
1.  🐸 問題が無ければmergeする。
1.  🐸 自分のローカルリポジトリに pull (fetch+merge) する。
1.  🐰 無事マージされたら作業ブランチを消す。

## Tips

-   習うより慣れる。
    最初はコマンドが多くて難しそう・面倒くさそうに感じるけど、
    だんだん意識しなくても使えるようになる。
-   `git status` やエラー文をちゃんと読む。
    どうすればいいかだいたい書いてくれてるし、
    そのままウェブ検索すればだいたい解決策が見つかる。
-   `--force` とか `-f` のような強制オプションは、
    間違えると取り返しがつかなくなるので基本的に使わない。
-   コミットを簡潔に要約するメッセージを書く。
    [好ましいスタイルについては諸説ある](https://www.google.co.jp/search?q=commit+message+best+practices)けど、
    とりあえず大文字で始まる命令形の一文を書くところから始めたらよいのでは。
    コミットの内容に応じた分類で[先頭に絵文字を入れるスタイル](https://github.com/carloscuesta/gitmoji/)も人気になりつつある。
-   GitHubにpushされたら自動的に
    [Slack](https://slack.com/)や[Twitter](https://twitter.com/)に投稿、
    というような連携が可能。
-   RStudioでもディレクトリを"Project"として扱うことでGitを活用できる。


## Glossary / 用語

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


## Further reading

- [GitHub](https://github.com/): 他人のGit活用事例が見放題。
- [GitHub Learning Lab](https://lab.github.com/):
  公式ボットが手取り足取り教えてくれるらしい。
- [Pro Git book](https://git-scm.com/book/en/v2): Gitの公式？本。
- [Bookdown](https://bookdown.org/): R Markdownで本を書く。
