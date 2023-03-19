+++
title = "Hugo"
subtitle =  "静的サイトを高速生成"
tags = ["writing"]
aliases = ["/dev/hugo.html"]
[menu.main]
  parent = "misc"
+++

Markdown記法のテキストをHTMLに変換する、静的ウェブサイト生成プログラム。
公式ドキュメントが充実しているので基本的にそれに従えば簡単にできる。

- <https://gohugo.io/documentation/>
- <https://github.com/gohugoio/hugo>

高速さとシンプルさに惹かれてSphinxから移行し、
本サイトもHugoでビルドしている。
オフラインの研究ノートとしても有用。

## Quickstart

<https://gohugo.io/getting-started/quick-start/>

-   [インストール方法はいろいろ用意されてる](https://gohugo.io/installation/)。
    例えばMacなら[Homebrew]({{< relref "homebrew.md" >}})で一発:
    `brew install hugo`

    バージョンを細かく指定したりソースコードを改変したりしたい場合はgitから:
    ```sh
    export GOPATH=${HOME}/.go
    export PATH=${PATH}:${GOPATH}/bin
    mkdir ${HOME}/src
    cd ${HOME}/src
    git clone https://github.com/gohugoio/hugo.git
    cd hugo
    go install --tags extended
    ```
    `--tags extended` はSASS/SCSS対応版をコンパイルするオプション。

-   骨組みを作る:
    ```sh
    cd path/to/site
    hugo new site .
    ```

-   ページをMarkdownで書く:
    ```sh
    hugo new about.md
    ```

    ```markdown
    +++
    date = 2016-02-26T19:10:22+09:00
    title = "About"
    +++

    ## Heading

    normal *italic* **bold**
    ```

-   テーマをとりあえず全部インストール:
    ```sh
    git clone --depth 1 --recursive https://github.com/gohugoio/hugoThemes.git themes
    ```

-   適当なテーマでウェブサーバーを走らせる:
    ```sh
    hugo server --theme blank
    ```

-   ブラウザから http://localhost:1313/about にアクセスしてみる。
    `hugo server`, `hugo -w` はファイルを監視していて変更をすぐに反映する。


## 設定

<https://gohugo.io/getting-started/configuration/>

長らく `config.toml` だったが今は `hugo.toml` がデフォルト。
`config/_default/hugo.toml` に置いても同じ。

`config/` 直下のディレクトリ名と `-e/--environment` オプションで切り替え可能。
ただしデフォルトの挙動がわかりにくい罠なので注意。例えば `config/_default/` と `config/production/` を持って `hugo` を実行するとproduction環境になってしまう。
production環境を作らず `config/public/` とかにしておけば明示的に `-e public`
を渡さない限り常にデフォルトのdevelopment環境になるので分かりやすい。

```sh
hugo        # -e production (confusing!)
hugo -w     # -e development
hugo server # -e development
```

## Theme

<https://themes.gohugo.io/>

デフォルトのテーマというものが存在しないのがちょっと厳しい。
ユーザーによっていろいろ投稿されてるけどほとんどがブログ用途。
ということで非ブログ用に簡単なものを自作して使っている:

https://github.com/heavywatal/hugo-theme-nonblog

### Performance

<https://gohugo.io/troubleshooting/build-performance/>

ページによって内容が変わらないテンプレートは `partial` の代わりに
[`partialCached`](https://gohugo.io/functions/partialcached/)
を使う。
ビルドするときに
`--templateMetrics --templateMetricsHints`
オプションを付けるとどのへんを変えたら良いか教えてくれる。


## Content

### Markdown

[CommonMark](https://spec.commonmark.org/)
: "Markdown"の正式な仕様というものが存在せず、
  いくつかの方言(flavor)が乱立していたが、
  現在ではこれが事実上の標準仕様となりつつある。
  [2017年からGFMがこれに準拠することになった](https://githubengineering.com/a-formal-spec-for-github-markdown/)のもよかった。

[GitHub Flavored Markdown (GFM)](https://github.github.com/gfm/)
: CommonMarkに準拠しつついくらかの機能を追加したもの。
  基本的な書き方は[GitHub Helpのページ](https://help.github.com/articles/basic-writing-and-formatting-syntax/)が読みやすい。

[Blackfriday](https://github.com/russross/blackfriday)
: HugoのMarkdownエンジンは長らくこれだった。
  CommonMark準拠じゃないし、
  リストまわりでの不具合が放置されてるし、
  などなど不満が募るうちにGoldmarkに取って代わられた。

[Goldmark](https://github.com/yuin/goldmark/)
: 2019年末からHugoはこっちに移行した。
  基本的にはCommonMark準拠だけど、
  デフォルト設定での生HTMLコードの扱いがちょっと変。


### Front matter

<https://gohugo.io/content-management/front-matter/>

タイトルや日付などのメタデータをファイルの先頭で記述する。
YAMLやJSONでもいいけど、
[TOML](https://github.com/toml-lang/toml)のほうが将来性ありそう。


## 閲覧・公開方法

### Hugo Server

付属の簡易サーバーを起動。
```
hugo server
open http://localhost:1313/
```

### localhost (Mac)

`public/` 以下に生成されるファイルをシムリンクやコピーで
`/Library/WebServer/Documents/` に配置すれば
<http://localhost> で閲覧できる。
ユーザーの `~/Sites/` をドキュメントルートにする方法でもいいが、
単にシムリンクを張るほうが楽ちん。

```sh
ls /Library/WebServer/
sudo mv /Library/WebServer/Documents /Library/WebServer/Documents.orig
sudo ln -s ~/path/to/public /Library/WebServer/Documents
open http://localhost/
```

[5秒待たされる問題](https://stackoverflow.com/questions/70698918/)
を避けるため `KeepAlive Off` を設定しておく。
`/etc/apache2/httpd.conf` を直に書き換えるか、
`/etc/apache2/other/keepalive.conf` のようなファイルを作り、
`sudo apachectl graceful` で変更を反映。


### GitHub Pages

`public/` 以下に生成されるファイルをしかるべきrepository/branchに置くだけ。

See [Git]({{< relref "git.md" >}})
