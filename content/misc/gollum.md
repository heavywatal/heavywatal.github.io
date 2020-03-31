+++
title = "Gollum"
subtitle = "MarkdownとGitで動くWikiエンジン"
date = 2020-03-31T19:21:46+09:00
tags = ["writing"]
[menu.main]
  parent = "misc"
+++

https://github.com/gollum/gollum

研究室内の連絡・情報共有は、フロー型の情報ならSlackで、ストック型の情報ならWikiで。
という方針になったので、学内ネットワークで閲覧編集可能な自前Wikiサーバーを立てる。
開発環境は macOS, 本番環境は Ubuntu 18.04 LTS.

## gollumインストールとWiki新規作成

1.  最新版のRubyを[rbenv](https://github.com/rbenv/ruby-build/wiki)で入れる:

    ```sh
    # sudo apt install autoconf bison build-essential libssl-dev libyaml-dev libreadline-dev zlib1g-dev libncurses5-dev libffi-dev libgdbm6 libgdbm-dev
    brew install rbenv
    rbenv install -l | less
    rbenv install 2.7.0
    rbenv global 2.7.0
    eval "$(rbenv init -)"
    ```

    `eval "$(rbenv init -)"`
    はここでインストールした `ruby` や `bundle` にPATHを通すコマンド。
    新しいシェルを起動するたびに実行する必要があるので
    `.bashrc` 等の設定ファイルに記述したくなるところだが、
    シェルの起動がかなり遅くなってしまうので適宜手動で実行することにする。

1.  Wiki用のリポジトリ(ここでは`labwiki`)を作成して空コミット:

    ```sh
    git init labwiki
    cd labwiki/
    git commit --allow-empty -m ":beer: Create repository"
    ```

1.  `Gemfile` を作成してコミット:

    ```gemfile
    source 'https://rubygems.org'
    gem 'commonmarker'
    gem 'gollum'
    ```

    ---

    gollum本体をいろいろいじくるために自分のフォークを使う場合:

    ```gemfile
    gem 'gollum', :github => 'heavywatal/gollum', :branch => 'custom'
    ```

    gollum が依存する gollum-lib のほうはその `Gemfile` じゃなくて
    `~/.bundle/config` とか `labwiki/.bundle/config` に設定:

    ```sh
    bundle config gollum-lib heavywatal/gollum-lib
    ```

    開発環境ではローカルのクローンを使うように設定:

    ```sh
    SRCDIR=${HOME}/fork
    git clone https://github.com/heavywatal/gollum-lib.git ${SRCDIR}/gollum-lib
    git clone https://github.com/heavywatal/gollum.git ${SRCDIR}/gollum -b custom
    bundle config --local local.gollum-lib ${SRCDIR}/gollum-lib
    bundle config --local local.gollum ${SRCDIR}/gollum
    ```

1.  `bundle install` で諸々インストール。

1.  `bundle exec gollum` でとりあえず走らせてみる。


## 設定

### 基本

`config.rb` を作成して
`bundle exec gollum -c config.rb`
のように指定して読ませる。

```ruby
require 'gollum/app'

wiki_options = {
  page_file_dir: 'source',
  css: true,
  mathjax: false,
  emoji: true
}
Precious::App.set(:wiki_options, wiki_options)
```

例えばこの場合、
ページのソースファイルはリポジトリのルートではなく `source/` から読まれるようになる。

`css: true` により `custom.css` を読み込まれるようになるが、
残念ながらリポジトリルートではなく `source/custom.css` に置かなければならない。
また、ローカルファイルではなくコミット済みのものが読まれることに注意。

なぜかここからは設定できずコマンドからのみ設定可能な項目もある。<br>
e.g., `-b/--base-path`, `--allow-uploads page`.


### ポート番号なしでアクセスする

デフォルトでは <http://example.com:4567> のように4567番ポートのルートで動く。
これを80番ポートの`/wiki`以下で動くように調整して
http://example.com/wiki/ のようにアクセスできるようにする。

1.  Apacheの設定ファイルを新規作成:

    ```sh
    sudo vim /etc/apache2/sites-available/gollum-wiki.conf
    ```

    ```apache
    ProxyRequests Off¬
    <Proxy *>¬
      Order deny,allow¬
      Allow from all¬
    </Proxy>¬
    ProxyPass /wiki http://localhost:4567/wiki¬
    ProxyPassReverse /wiki http://localhost:4567/wiki¬
    ```

1.  その設定ファイルを有効化してApache再起動:

    ```sh
    sudo a2ensite gollum-wiki.conf
    sudo systemctl restart apache2
    ```

1.  `bundle exec gollum -b /wiki` で起動。

    `config.rb` のほうで `Precious::App` に
    `base_path: '/wiki'` を渡してもなぜか効かない。


### BASIC認証でやんわりパスワードをかける

`config.rb` にこんな感じで書くだけ:

```ruby
module Precious
  class App < Sinatra::Base
    passwd = {
      'user1' => '0b14d501a594442a01c6859541bcb3e8164d183d32937b851835442f69d5c94e',
      'user2' => '6cf615d5bcaac778352a8f1f3360d23f02f34ec182e259897fd6ce485d7870d4'
    }
    use Rack::Auth::Basic, 'Private Wiki' do |username, password|
      digested = Digest::SHA256.hexdigest(password)
      if passwd.key?(username) && digested == passwd[username]
        Precious::App.set(:username, username)
      end
    end

    before do
      session['gollum.author'] = {
        name: settings.username,
        email: "#{settings.username}@users.noreply.github.com"
      }
    end
  end
end
```

もっとちゃんとした認証システムにしたほうがいいのかもしれないけど、
大学のファイアウォール外からはアクセス不可能なのでとりあえずこれくらいで...

パスワードのハッシュ値は
`sha256sum <(pbpaste)`
とか
`echo -n 'greatpassword' | sha256sum`
のようなコマンドで計算できる。


### Markdownパーサー/レンダラを設定する

https://github.com/gollum/gollum/wiki/Custom-rendering-gems

Markdownを読んでHTMLに変換するライブラリは
[github-markup](https://github.com/github/markup)
を通して選択できるようになっている。
デフォルトでは
[kramdown](https://kramdown.gettalong.org/)
が利用されるらしいが、
なるべく[CommonMark+GFM](https://github.github.com/gfm/)準拠で高速なのが良いので、
[commonmarker](https://github.com/gjtorikian/commonmarker)
を使うことにする。
例によって `config.rb` に追記:

```ruby
module Gollum
  class Markup
    GitHub::Markup::Markdown::MARKDOWN_GEMS.clear
    GitHub::Markup::Markdown::MARKDOWN_GEMS['commonmarker'] = proc do |content|
      exts = %i[
        table
        tasklist
        strikethrough
        autolink
      ]
      parse_opts = %i[
        UNSAFE
        SMART
      ]
      render_opts = %i[
        UNSAFE
        GITHUB_PRE_LANG
      ]
      doc = CommonMarker.render_doc(content, parse_opts, exts)
      doc.to_html(render_opts)
    end
  end
end
```

ほかにどんなのが利用可能かはこちらを参照:
https://github.com/github/markup/blob/master/lib/github/markup/markdown.rb
