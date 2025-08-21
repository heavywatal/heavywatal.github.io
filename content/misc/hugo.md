+++
title = "Hugo"
subtitle =  "静的サイトを高速生成"
tags = ["writing", "web"]
[menu.main]
  parent = "misc"
+++

[Markdown]: {{< relref "markdown.md" >}}

[Markdown]記法のテキストをHTMLに変換する、静的ウェブサイト生成プログラム。
公式ドキュメントが充実しているので基本的にそれに従えば簡単にできる。

- <https://gohugo.io/documentation/>
- <https://github.com/gohugoio/hugo>

高速さとシンプルさに惹かれてSphinxから移行し、
本サイトもHugoでビルドしている。
オフラインの研究ノートとしても有用。

## Quick start

<https://gohugo.io/getting-started/quick-start/>

-   Hugo本体をインストール。
    [方法はいろいろ用意されてる](https://gohugo.io/installation/)。

    - 手動でよければOSに合った[公式prebuilt binary](https://github.com/gohugoio/hugo/releases)をダウンロードしてPATHを通すのが簡単。
    - コマンドで管理するなら[Homebrew]({{< relref "homebrew.md" >}})で一発:
    `brew install hugo`
    - ソースコードを改変したりしたい場合は [git]({{< relref "git.md" >}}) から:
      ```sh
      export GOPATH=${HOME}/.go
      export PATH=${PATH}:${GOPATH}/bin
      mkdir ${HOME}/src
      cd ${HOME}/src
      git clone https://github.com/gohugoio/hugo.git
      cd hugo
      go install -v
      ```
      SCSSのための `--tags extended` オプションは不要になった。
-   SCSSを使う場合は Dart Sass を別途インストールしてPATHを通す。
    方法はいろいろあるけど
    [公式prebuilt binary](https://github.com/sass/dart-sass/releases)
    をダウンロードするか
    [Homebrew]({{< relref "homebrew.md" >}}) を使うのが簡単:
    ```sh
    brew install heavywatal/tap/dart-sass
    # or
    brew install sass/sass/sass
    ```
-   ちゃんとインストールできているか確認: `hugo env`
-   空っぽのディレクトリに移動して骨組みを作る:
    ```sh
    cd path/to/site
    hugo new site .
    ```

-   設定ファイル `hugo.toml` を作ってテーマを指定:
    ```toml
    [module]
    [[module.imports]]
    path = "github.com/theNewDynamic/gohugo-theme-ananke"
    ```

-   適当なページを作る `hugo new about.md`:
    ```markdown
    +++
    title = "About"
    +++

    ## Heading

    normal *emphasis* **strong**
    ```

-   ウェブサーバーを走らせる:
    ```sh
    hugo server
    ```
-   ブラウザから <http://localhost:1313/about> にアクセスしてみる。
    `hugo server`, `hugo -w` はファイルを監視していて変更をすぐに反映する。


## Configuration

<https://gohugo.io/configuration/>

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


## Template

- <https://gohugo.io/templates/>
- <https://gohugo.io/functions/>
- <https://gohugo.io/methods/>

先頭のドット `.` は現在のcontextを表す。
普通のページのテンプレートではページを表す状態から始まり、
`{{ with }}` や ``{{ end }}`` などによって階層的に変化していく。

`.Site.Params.author.name` のように chain してメソッドや変数を参照できる。
TOMLやYAMLではハイフン `-` がキーに入っても問題ないが、
Hugo template では chain できなくなるため非推奨。
[`index`](https://gohugo.io/functions/collections/indexfunction/)
関数を使えばそういうやつでも参照できるけど。

[`{{ .Param "key" }}`](https://gohugo.io/methods/page/param/)
: `.Params.key` が存在すればそれ、しなければ `site.Params.key` を参照する。
: `{{ .Params.key }}`:
  各ページの front matter の `[params]` セクションに書いたものを参照。
: `{{ site.Params.key }}`
  サイト全体の設定 `hugo.toml` の `[params]` セクションに書いたものを参照。


### Performance

<https://gohugo.io/troubleshooting/build-performance/>

ページによって内容が変わらないテンプレートは `partial` の代わりに
[`partialCached`](https://gohugo.io/functions/partialcached/)
を使う。
ビルドするときに
`--templateMetrics --templateMetricsHints`
オプションを付けるとどのへんを変えたら良いか教えてくれる。


## Content

<https://gohugo.io/content-management/>

[Markdown attributes](https://gohugo.io/content-management/markdown-attributes/)
とか
[Shortcodes](https://gohugo.io/content-management/shortcodes/)
とかの機能は便利だけど [Markdown] からの逸脱を最小限に留めたい気もする。


### Front matter

<https://gohugo.io/content-management/front-matter/>

タイトルや日付などのメタデータをファイルの先頭で記述する。
YAMLやJSONでもいいけど、
[TOML](https://toml.io/)のほうが将来性ありそう。


## 閲覧・公開方法

### Hugo Server

付属の簡易サーバーを起動。
```
hugo server
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

`public/` 以下に生成されるファイルを `gh-pages` ブランチにアップロード。

See [Git #github-pages]({{< relref "git.md#github-pages" >}}).
