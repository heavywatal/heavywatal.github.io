+++
title = "Hugo"
subtitle =  "静的サイトを高速生成"
tags = ["writing"]
[menu.main]
  parent = "dev"
+++

Markdown記法のテキストをHTMLに変換する、静的ウェブサイト生成プログラム。
公式ドキュメントが充実しているので基本的にそれに従えば簡単にできる。

- https://gohugo.io/overview/introduction/
- https://discuss.gohugo.io/
- https://github.com/spf13/hugo

高速さとシンプルさに惹かれてSphinxから移行し、
本サイトもHugoでビルドしている。

## Quickstart

http://gohugo.io/overview/quickstart/

* インストール
    ```sh
    % brew install hugo
    ```
  まだまだ開発途上で機能追加や不具合修正も活発なので
  githubをフォローして最新版を入れたほうがいいかも
  https://gohugo.io/overview/installing/
    ```sh
    export GOPATH=${HOME}/.go
    export PATH=${PATH}:${GOPATH}/bin
    go get -v github.com/spf13/hugo
    ```

* 骨組みを作る
    ```sh
    cd path/to/site
    hugo new site .
    ```

* ページをMarkdownで書く
    ```sh
    hugo new about.md
    ```

    ```markdown
    +++
    date = "2016-02-26T19:10:22+09:00"
    title = "About"
    +++

    ## Heading

    normal *italic* **bold**
    ```

* テーマをとりあえず全部インストール
    ```sh
    git clone --depth 1 --recursive https://github.com/spf13/hugoThemes.git themes
    ```

* ウェブサーバーを走らせる
    ```sh
    hugo server -w D -t hyde
    ```

* ブラウザから http://localhost:1313/about にアクセスしてみる。
  監視オプション `-w` を付けておけば、ファイルの変更がすぐに反映される。


## 設定

https://gohugo.io/overview/configuration/

`config.toml`

## Theme

http://themes.gohugo.io/

デフォルトのテーマというものが存在しないのがつらい。
ユーザーによっていろいろ投稿されてるけどほとんどがブログ用途。
ということで非ブログ用に簡単なものを自作して使っている:

https://github.com/heavywatal/hugo-theme-nonblog

## Content

### Markdown

正式な仕様が未だに存在せず、いくつかの方言(flavor)が存在する。

[CommonMark](http://spec.commonmark.org/)
: 標準仕様決定に向けて議論中。

[GitHub Flavored Markdown (GFM)](https://help.github.com/articles/basic-writing-and-formatting-syntax/)
: いま最もよく書かれているのはこれかな？ Atomでも標準サポート。

[Blackfriday](https://github.com/russross/blackfriday)
: HugoのMarkdownエンジンはこれ。
  残念ながら上記2つとも微妙に違う。

### [Front matter](https://gohugo.io/content/front-matter/)

タイトルや日付などのメタデータをファイルの先頭で記述する。
YAMLやJSONでもいいけど、
[TOML](https://github.com/toml-lang/toml)のほうが将来性ありそう。
