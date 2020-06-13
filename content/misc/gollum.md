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

## ソフトウェア選定

[pukiwiki](https://pukiwiki.osdn.jp/)
:   学生のときの研究室で使ってたので馴染み深い。
    でも独自記法だしphpとか文字コードとか考えたくないので却下。

[crowi](https://github.com/crowi/crowi/)
:   Node.js + MongoDB で動くモダンな Markdown wiki。
:   生のファイルが見えないデータベースっぽいので管理が難しそう。
:   [growi](https://growi.org/) はこれをフォークしたもので、
    機能もドキュメントも強化されてるし、
    docker-compose とかですぐ使えるのも楽ちん。
:   日本語の人しか使わなそう...?

[gitit](https://github.com/jgm/gitit)
:   pandoc + git で動くのでかなり手堅い感じ。
:   Haskell の勉強を兼ねていじくり回す時間があれば...

[gollum](https://github.com/gollum/gollum)
:   Ruby + git で動く Markdown wiki。
:   GitHubやGitLabのWikiにも採用されているのでコミュニティが大きそう。
:   自前サーバーを管理できる人が抜けても内部データを簡単に再利用可能。
:   Rubyはよく知らないけど理解しなくても雰囲気でいじれそう。

[Hugo]({{< relref "hugo.md" >}})
:   静的ウェブサイトを作る用途には最高だけどWiki機能は無い。
:   [Netlify CMS](https://www.netlifycms.org/)
    でWiki-likeなガワを取り付けることは可能だけど、
    編集内容のpush先がGitHubとかになるので、
    それを学内サーバーに即時反映させるのが難しい。


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

    gollum本体をいろいろいじくる場合は自分のフォークを使う:

    ```gemfile
    gem 'gollum-lib', :github => 'heavywatal/gollum-lib', :branch => 'custom'
    gem 'gollum', :github => 'heavywatal/gollum', :branch => 'custom'
    ```

    開発環境ではローカルのクローンを使うように設定:

    ```sh
    SRCDIR=${HOME}/fork
    git clone https://github.com/heavywatal/gollum-lib.git ${SRCDIR}/gollum-lib
    git clone https://github.com/heavywatal/gollum.git ${SRCDIR}/gollum -b custom
    bundle config local.gollum-lib ${SRCDIR}/gollum-lib
    bundle config local.gollum ${SRCDIR}/gollum
    ```

1.  `bundle install` で gollum 及び依存パッケージをまとめてインストール。
    `--local` を付けてこのプロジェクト専用にしてもよい。

1.  `bundle exec gollum` でとりあえず走らせる。
    手元のコンピュータなら <http://localhost:4567> で確認。


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
    ProxyRequests Off
    <Proxy *>
      Order deny,allow
      Allow from all
    </Proxy>
    <Location /wiki>
      ProxyPass http://localhost:4567/wiki
      ProxyPassReverse http://localhost:4567/wiki
    </Location>
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
    use Rack::Auth::Basic, 'Private Wiki' do |username, password|
      users = File.open(File.expand_path('../users.json', __FILE__)) do |file|
        JSON.parse(file.read, symbolize_names: true)
      end
      name = username.to_sym
      digested = Digest::SHA256.hexdigest(password)
      if users.key?(name) && digested == users[name][:password]
        Precious::App.set(:author, users[name])
      end
    end

    before do
      session['gollum.author'] = settings.author
    end
  end
end
```

`session['gollum.author']` にハッシュを渡しておくとコミッターに反映してもらえる。
ユーザー情報は別ファイル(ここでは`users.json`)に分離しといたほうが見通しがいい。

```json
{
  "user1": {
    "name": "First User",
    "email": "user1@example.com",
    "password": "0b14d501a594442a01c6859541bcb3e8164d183d32937b851835442f69d5c94e"
  },
  "user2": {
    "name": "Second User",
    "email": "user2@example.com",
    "password": "6cf615d5bcaac778352a8f1f3360d23f02f34ec182e259897fd6ce485d7870d4"
  }
}
```

もっとちゃんとした認証システムにしたほうがいいのかもしれないけど、
大学のファイアウォール内なのでとりあえずこれくらいで...

パスワードのハッシュ値は
`sha256sum <(pbpaste)`
とか
`echo -n 'greatpassword' | sha256sum`
のようなコマンドで計算できる。


### Markdownパーサー/レンダラを変更する

https://github.com/gollum/gollum/wiki/Custom-rendering-gems

Markdownを読んでHTMLに変換するライブラリは
[github-markup](https://github.com/github/markup)
を通して選択できるようになっている。
デフォルトでは
[kramdown](https://kramdown.gettalong.org/)
が利用されるらしいが、
なるべく[CommonMark/GFM](https://github.github.com/gfm/)準拠で高速なのが良いので、
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

Gollum::Markup.formats.select! { |k, _| k == :markdown }
```

ほかにどんなのが利用可能かはこちらを参照:
https://github.com/github/markup/blob/master/lib/github/markup/markdown.rb


### `systemd` で自動的に開始

```sh
sudo vim /etc/systemd/system/gollum.service
```

```ini
[Unit]
Description=Gollum wiki server
After=network.target

[Service]
Type=simple
User=YOURNAME
WorkingDirectory=/path/to/your/labwiki
ExecStart=bundle exec gollum -c config.rb -b /wiki --allow-uploads dir
Restart=on-abort
StandardOutput=file:/var/log/gollum.log
StandardError=file:/var/log/gollum.log

[Install]
WantedBy=multi-user.target
```

```sh
sudo systemctl start gollum.service
sudo systemctl enable gollum.service
```
