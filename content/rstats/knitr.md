+++
title = 'knitr'
subtitle = 'Markdownにコード実行結果を編み込む'
tags = ['r', 'tidyverse']
[menu.main]
  parent = 'rstats'
  weight = -58
+++

<a href="https://yihui.org/knitr/">
<img src="/_img/hex-stickers/knitr.webp" align="right" width="120" height="139">
</a>

特殊Markdownファイル (.Rmd, .qmd) に含まれるRコードを実行し、
結果を編み込んで汎用Markdownファイル (.md) に変換する。

一般的には [Quarto](https://quarto.org/) や
[R Markdown](https://rmarkdown.rstudio.com/)
といったパイプラインの中で利用されることが多く、
直接触れる必要性はあまりない。
私はMarkdown→HTML処理をPandocではなく[Hugo]({{< relref "hugo.md" >}})に任せたいので
`.Rmd` |> knitr |> Hugo という流れで使う。

````r
## Heading

Paragraph.

```{r, example-code-chunk}
#| fig.height: 5
#| fig.width: 6
answer = 6 * 7
ggplot(mpg) + aes(displ, height) + geom_point()
```

Answer to the ultimate question is `r answer`.
````

- <https://yihui.org/knitr>
- [逆順で理解する R Markdown Presentation](/slides/tokyor102/) 2022-10
  [Tokyo.R #102](https://tokyor.connpass.com/event/262836/)

## Options

<https://yihui.org/knitr/options/>


### Format

chunk header `{r, tag=value, tag=value}` 形式
: knitrではこれが基本という扱い。
: 値はRコードとして評価される。文字列にはquoteが必要。
: 改行は許されない。
: `{r` 直後の `,` はあってもなくてもいい。
  Yihuiは付けているが、RStudioの新規作成RmdやPosit社のcheatsheetでは付いてない。
: 先頭のオプションはtag無しquote無しで `label` 扱いされる。
  明示的に `label` の値を設定してもよい。

chunk内 `#| tag=value, tag=value` 形式
: Rコードも書けるし改行も許される。
: 改行はカンマの代わりにはならない。

chunk内 `#| tag: value` 形式
: Quartoでは[こちらを推奨](https://quarto.org/docs/computations/r.html#chunk-options)。
: 一行一項目で改行する前提。
: 区切り文字としてピリオドよりハイフンが好まれる。
  `^(fig|out)-` は `\1.` に置換された上でheader形式のオプションとmergeされる。
: 値はYAMLの型である必要がある。
  - 論理値は小文字: `true`, `false`
  - 文字列はquote無しでもいいが付けておいたほうがたぶん安全。
    "double" では `\n` が改行扱いされるなどエスケープシーケンスが有効。
    'single' ではそういうのが起こらず文字通り。
  - Rコードを渡すにはひと手間: `#| message: !expr 'NA'`

一括 `knitr::opts_chunk$set(tag = value)` 形式
: 文書内、それより後ろのchunkのデフォルトを変更する。


### Chunk label

区切りにはなぜか**ハイフンが推奨**。ピリオドやアンダースコアではなく。

figureやcacheの出力先に使われるのでuniqueになるように注意する。
文書内での重複はknitrがチェックして怒ってくれるが、
ディレクトリ内の文書間での重複までは見てくれない。

省略すると `unnamed-chunk-%d` のような形で自動的に割り振られる。
入力ファイルによって変わるように
`unnamed.chunk.label` オプションを変更しておいたほうが安全。


### General

- `knitr::opts_chunk$set(error = FALSE)`:
  knitr単独で使う場合、デフォルトではchunk内でエラーが生じても止まらず危険。
  あえてエラーを表示したいchunkでのみ明示的に `error = TRUE` とするように。
- `knitr::opts_chunk$set(comment = "")`:
  出力結果の行頭に何も書かないようにする。
  `NA` でも `NULL` でも同じ。
  デフォルトでは `## ` が入ってしまう。
- `knitr::opts_chunk$set(message = NA)`:
  `message()` からの出力をknitrで捕まえずにconsoleにそのまま流す。
  [knitr 1.42](https://github.com/yihui/knitr/blob/master/NEWS.md#changes-in-knitr-version-142)
  から `message = FALSE` はどこにも出力しないようになって危険。
  `utils::capture.output(..., type = "message")` したいときとかも。
- `knitr::opts_chunk$set(warning = NA)`: 同上。
- `echo: false`: コードを表示しない。
  `echo: -1` のように数値で特定の行だけを除外したりもできる。
- `results: 'markup'`: 標準出力があるたびにコードフェンスに入れて編み込む。(デフォルト)<br>
  `results: 'hide'`: 表示しない。`FALSE` でも。<br>
  `results: 'asis'`: コードフェンスに入れず生の文字列として編み込む。<br>
  `results: 'hold'`: 標準出力をchunk最後まで溜めてまとめて編み込む。<br>
- `eval: false`: 評価しない。
- `include: false`: コードも結果も出力しない。実行だけする。
  空行が1行だけ追加されるのを防ぐのは難しそう。
- `collapse: true`: 実行結果をコードと同じブロックに入れる。
- `ref.label`: ほかのchunkをソースコードとして読み込む。
  順番は関係なく、後ろのコードを前で参照することもできる。
  [複数のchunkに同じlabelをつけて1つだけ`eval=TRUE`にする](https://bookdown.org/yihui/rmarkdown-cookbook/reuse-chunks.html#same-label)
  という手もある。
  [`<<label>>`](https://bookdown.org/yihui/rmarkdown-cookbook/reuse-chunks.html#embed-chunk) で埋め込むのはちょっとお行儀が悪い気がするので使わない。
- `opts.label`: ほかのchunkのオプションを継承する。


### External

- `file: 'setup.R'`: 外部ファイルをchunkのソースコードとして読み込む。
  `code: expr! 'readLines("setup.R")'` も同様。
  `{r}` 以外のengineも指定できる。
  `{embed}` なら拡張子を考慮しつつコードフェンスに入れて埋め込み。
  `{asis}` なら文字通りそのまま。
- `{cat, engine.opts = list(file = "hello.py")}`:
  chunkの内容をファイルに書き出す。
  そのままでは表示されず `echo` も効かない。
  `lang` を追加設定することで言語付きコードフェンスになる。


### Cache

- `knitr::opts_chunk$set(cache = TRUE)`:
  実行結果を保存しておき、chunk内容に変更が無ければ再利用する。
  乱数生成やchunk間の依存関係などによって危険性が跳ね上がる諸刃の剣。
  invalidationにはchunkコードとオプションのMD5 digest値が使われる。
- `knitr::opts_chunk$set(autodep = TRUE)`:
  グローバル変数の依存関係を自動で解決。
- `knitr::opts_chunk$set(cache.path = glue(".cache/{inputname}/"))`:
  名前からの予想に反してpathそのものではなくprefix。
  相対パスの基準はchunk内の `getwd()` ではなく
  `knit()` を呼び出す側の `getwd()`
  (i.e., 推奨に従えばinputではなくoutput側のディレクトリ)。
  デフォルトの `cache/` ではhugoに拾われたりして邪魔なので、
  `.cache/` に変えるとか絶対パスで別の場所を指定するとか。
  複数の文書を同一ディレクトリで扱う場合は入力ファイル名を含めたほうがいい。
  chunk labelの重複を自力で避けたとしても `autodep` の `__global` などが衝突してしまうので。
- `cache.vars`: 保持する変数を明示的に制限する。
- `cache.globals`: そのchunkで作られないグローバル変数を明示的に指定して `autodep` を助ける。
- `cache.rebuild: !expr '!file.exists("some-file")'`:
  強制的にchunkを評価してcacheを作り直す。
  外部ファイルの存在などに依存させたいとき便利。
- `cache.extra`:
  invalidationの依存関係を追加する。
  `tools::md5sum("some-file")` とか `file.mtime("some-file")`
  とかで外部ファイルの内容や変更時間を参照できる。
  そもそもすべてのchunkオプションがinvalidation評価対象なので、
  tagが `cache.extra` である必要はない。
  それゆえにknitrコードにもドキュメントにもちゃんと書かれていなくてわかりにくい。
  [yihui/knitr#405](https://github.com/yihui/knitr/pull/405)
- `dependson`: 依存するchunkをlabelや数値で指定。


### Figure

- `fig.path`: `cache.path` と同様の挙動。
- `fig.width: 10`
- `fig.height: 5`
- `fig.show: hold`
- `fig.show: animate`
  - `animation.hook: gifski` <https://github.com/r-rust/gifski>
  - `interval: 0.25`
- `dpi`: デフォルトは72。
- `fig.retina`: `dpi` に掛け算、 `out.width` に割り算をして
  見かけサイズ据え置きで解像度を変更する。
  `<img>` タグに `width` attribute が追加されることになる。
- `knitr::opts_chunk$set(fig.process = wtl::oxipng)`:
  [PNGを圧縮](https://github.com/shssoichiro/oxipng)したり
  [WebPに変換](https://developers.google.com/speed/webp/docs/cwebp)したりできる。
- `knitr::opts_chunk$set(dev = "ragg_png")`:
  日本語ラベルが文字化けしない、描画が速い、などと言われる
  [ragg](https://github.com/r-lib/ragg/) を使う場合。
  XQuartzのPNGに比べるとなんとなく文字がガチャガチャで美しくない気がするので、
  日本語を使うchunkでだけ設定するほうがいいかも。


### Package Options `knitr::opts_knit`

<https://yihui.org/knitr/options/#package-options>

knitrロード前に `options(knitr.package.verbose = TRUE)` とするか、
`knit()` 実行前に `knitr::opts_knit$set(verbose = TRUE)` とするか。
chunk内からの実行では遅い。
`knitr::opts_chunk` とは別であることに注意。

- `verbose = TRUE`: 帯に短し襷に長し。
  `progress = TRUE` のときはchunkの内容をいちいち表示するので量が多すぎ、
  `progress = FALSE` のときはキャッシュを使った場合にその旨が表示されるだけで情報不足。
- `unnamed.chunk.label = glue("_{inputname}")`:
  複数の文書を同一ディレクトリで扱う場合を考えるとデフォルトの
  `unnamed-chunk` ではやや不安。
- `root.dir`: コード評価のディレクトリ。デフォルト(inputファイルと同じ)を維持。
- `base.dir`: plotの出力先を絶対パスで指定したいことがもしあれば。
  chunk optionではないことに注意。


### Global R Options

<https://yihui.org/knitr/options/#global-r-options>

なぜ `knitr::opts_knit` とまた違うくくりがあるのかは謎。歴史的経緯？

- `knitr.progress.fun`:
  デフォルトではchunk名の確認さえままならないほどコンパクトなので、
  例を参考にして適当に差し替える。


## Functions

[`knitr::knit(input, output = NULL, ...)`](https://github.com/yihui/knitr/blob/master/R/output.R)
: `input`: 呼び出し元とは違うディレクトリにあるファイルを指定しても、
  chunk内の `getwd()` はこのファイルが基準となる。
: `output`: 可能な限り `NULL` のままにしておく。
  `input` とは違うディレクトリに書き出したい場合、
  出力先に予め `setwd()` しておくことが強く推奨されている。
  (個人的には作業ディレクトリをinput側に統一するほうが簡単そうに思えるけど。)
  - Dangerous: `knit("report.Rmd", "outdir/report.md")`
  - Redundant: `setwd("outdir"); knit("../report.Rmd", "report.md")`
  - Good: `setwd("outdir"); knit("../report.Rmd")`

[`knitr::fig_chunk(label, ext, number, ...)`](https://bookdown.org/yihui/rmarkdown-cookbook/fig-chunk.html)
: chunk名から図のパスを取得。
  離れたところで図を使い回せる。

`knitr::current_input()`
: 処理中の入力ファイル名を取得。

`knitr::knit_exit()`
: 文書の途中で終了。


## Hooks

- <https://yihui.org/knitr/hooks/>
- <https://bookdown.org/yihui/rmarkdown-cookbook/chunk-hooks.html>

### Chunk hooks

- chunk前後に評価するコードをsetできる。
- setしたオプションが `NULL` でなければ(`FALSE` でも)評価される。
  cacheによってchunk本体が評価されない場合は評価されない。
- オプションの名前は output hooks と被らなければ何でも。
- cache invalidation より後なのでここで `options` の値を変更しても遅い。
- この関数が返す文字列はasisで埋め込まれる。
- `getwd()` は `knit()` 実行元と同じであり、
  chunk内と同じとは限らないことに注意。

```r
knitr::knit_hooks$set(foo = \(before, options, envir) {
  if (before) {
    # evaluated before chunk code
    options$cache.rebuild = TRUE  # too late!
  } else {
    # evaluated after chunk code
  }
})
```

### Output hooks

- コード評価後、コードや実行結果の文字列をいじる。
- `source`, `output`, `warning`, `message`, `error`,
  `plot`, `inline`, `chunk`, `document`.
- `getwd()` は場合によって異なる謎仕様なので注意。
  Rの場合は `knit()` 実行元と同じ、
  engineを変えた場合はchunkと同じ、、、？

```r
knitr::knit_hooks$set(source = \(x, options) {
  message("hook source ", getwd())
  paste0("src: ", x)
})
```

### Option hooks

- `options` を受け取り、変更し、返す関数をsetしておく。
- setしたオプションの値が `NULL` でさえなければ評価される。
  つまり、元々値を持つ `echo` などにsetすると毎回評価される。
- cache invalidation前に評価される唯一のhookか。
- `getwd()` は `knit()` 実行元と同じであり、
  chunk内と同じとは限らないことに注意。

```r
knitr::opts_hooks$set(fig.square = \(options) {
  stopifnot(is.numeric(fig.square))
  options$fig.width  = options$fig.square
  options$fig.height = options$fig.square
  options
})
```

## Quarto

https://quarto.org/docs/output-formats/hugo.html
