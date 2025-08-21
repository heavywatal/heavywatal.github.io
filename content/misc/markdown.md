+++
title = 'Markdown'
subtitle =  "軽量マークアップ言語"
tags = ["writing", "web"]
[menu.main]
  parent = "misc"
[params]
  toc = true
+++

[Hugo]: {{< relref "hugo.md" >}}

## Introduction

- <https://spec.commonmark.org/current/>
- <https://docs.github.com/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax>

HTMLに変換することを前提に設計された軽量マークアップ言語。
ある程度HTMLを理解しておいたほうが使いやすい。
同類他言語と比べて、変換前のテキストのままでも読みやすいのが特徴。


### Flavors

[CommonMark](https://spec.commonmark.org/)
: Markdownの正式な仕様というものが存在せず、
  いくつかの方言(flavor)が乱立していたが、
  現在ではこれが事実上の標準仕様となりつつある。
  [2017年からGFMがこれに準拠することになった](https://githubengineering.com/a-formal-spec-for-github-markdown/)のもよかった。

[GitHub Flavored Markdown (GFM)](https://github.github.com/gfm/)
: CommonMarkに準拠しつついくつかの機能を追加したもの。
: [`table`](#tables),
  [`tasklist`](#task-lists),
  [`strikethrough`](#emphasis),
  [`autolink`](#links),
  [`tagfilter`](#raw-html)

[Hugo] / [Goldmark](https://github.com/yuin/goldmark/)
: 基本的には CommonMark/GFM 準拠だけど独自の設定や拡張がたくさん。
: [Markdown attributes](https://gohugo.io/content-management/markdown-attributes/)
: [Shortcodes](https://gohugo.io/content-management/shortcodes/)
: [Summaries](https://gohugo.io/content-management/summaries/)
: <https://gohugo.io/configuration/markup/#goldmark>

[Pandoc's Markdown](https://pandoc.org/MANUAL.html#pandocs-markdown)
: 独自の拡張が多くて怖いけど `--from`/`--to` オプションに
  `commonmark` や `gfm` を渡せば大丈夫そう。

[knitr]({{< relref "knitr.md" >}}) / R Markdown / Quarto
: コードブロックのRやPythonを実行し、その結果をMarkdownに編み込む。


### Rules

- 記号は HTML entities に自動変換される。HTML entities を書くとそのまま出力。
  - Markdown `2>&1` → HTML `2&gt;&amp;1` → 2>&1
  - Markdown `2&gt;&amp;1` → HTML `2&gt;&amp;1` → 2&gt;&amp;1
- backslashをASCII記号の前に置くとエスケープ
  - Markdown記法を無効化し、記号そのものとして表示: \*em\* \\
  - 開始する文字だけエスケープすれば足りる: \[text](href) \<a>
  - そのほかの前に置かれたbackslashはそのまま表示: \n\t\1
- [blockのほうがinlineよりも先に処理される](https://spec.commonmark.org/current/#blocks-and-inlines)
  - `one
  - two`


## Leaf blocks

### Horizontal rules

`***`, `---`, `___` だけの行は `<hr>` に変換される。

***
---
___


### Headings

ATX style
: 行頭の `#` の数に応じて `<h1>` から `<h6>` に変換される。
: `### h3` のように半角スペースを挟む必要がある。
: `### h3 #####` のように後ろを修飾してもいい

SETEXT style
: 次の行に `===` または `---` を並べて下線のようにすると、
  前者は `<h1>`, 後者は `<h2>` に変換される。


#### h4 ####

##### h5

###### h6


### Fenced Code blocks

3連backquoteで挟むと `<pre><code> </code></pre>` に変換される。
言語を指定すると
[syntax highlighting](https://gohugo.io/content-management/syntax-highlighting/)
できることが多い。

```python
def hello():
    print("Hello, world!") # コメント
```
```
Hello, world!
```

未定義の言語を指定すると `pre` に `chroma` クラスが付かない:
```unknown-language
print("Hello, world!")
```

`text`, `plain`, `no-highlight` を指定すると `chroma` あり色なし:
```text
print("Hello, world!")
```

3連backquoteを含むコードは4連backquoteで表現できる:
````markdown
```python
def hello():
    print("Hello, world!") # コメント
```
````

基本的には fenced code blocks を書くとして、
行頭にスペース4つ置くと indented code block になることは頭に入れておく。

    インデントのつもりがコードブロックになったりするので。


### HTML blocks

次の条件を満たすところはMarkdownとしては変換されずそのままHTMLに出力される。
自動的に `<p></p>` で囲まれたりもしない。

1. 行頭 `<pre`, `<script`, `<style`, `<textarea` からその閉じタグまで。
   空行では途切れない。
2. `<!--` から `-->` まで
3. `<?` から `?>` まで
4. `<!` とそれに続くASCII文字から `>` まで（Hugoでは無効？）
5. `<![CDATA[` から `]]>` まで
6. 行頭blockタグから空行まで。
   タグは不完全でも閉じタグでも構わない。
   閉じタグでは終わらない。
   - `<div>`, `<p>`, `<blockquote>`
   - `<h1>`, `<h2>`, `<h3>`, `<h4>`, `<h5>`, `<h6>`
   - `<li>`, `<ul>`, `<ol>`, `<menu>`
   - `<dl>`, `<dt>`, `<dd>`
   - `<main>`, `<article>`, `<aside>`, `<footer>`, `<header>`, `<nav>`, `<section>`, `<address>`
   - `<table>`, `<thead>`, `<tbody>`, `<tfoot>`, `<tr>`, `<th>`, `<td>`, `<caption>`
   - `<details>`, `<summary>`
   - `<figure>`, `<figcaption>`
   - `<search>`, `<fieldset>`, `<form>`, `<legend>`
7. inlineでも独自でも、任意のタグだけの行から空行まで。
   前項とは違って、開始タグは一行で完結している必要がある。

type 6 であれば、直前のMarkdownブロックは空行が無くても `</p>` で閉じるし、
開始タグ直後の改行で独立させる必要もない:
<blockquote>*HTML block type 6*</blockquote>
*終了タグ後も空行を入れないとHTML block続行。*

閉じタグ前であっても空行でMarkdownに戻る。開始タグ直後である必要もない:
<blockquote>
*HTML block type 6*

**Markdown block** wrapped in `<p></p>`
</blockquote>

非blockタグで始めるには空行でMarkdownの段落を終わらせたうえ、
開始タグ単独で1行を終える必要がある:

<del>
*HTML block type 7*
</del>

そうじゃないと
<del>
*raw HTML* in **Markdown block** wrapped in `<p></p>`
</del>

<del>*raw HTML* in **Markdown block** wrapped in `<p></p>`</del>

<!-- *HTML comment* -->

<!--
*HTML comment*
-->

<? *Processing Instruction* ?>

<?
*Processing Instruction*
?>

<![CDATA[ *Character Data* ]]>

<![CDATA[
*Character Data*
]]>


### Paragraphs

地の文は `<p></p>` で囲まれて段落になる。

空行をひとつ以上挟むと改段落。
何行空けても出力結果は同じ。


### Tables

[GFM `table` extension](https://github.github.com/gfm/#tables-extension-)

| normal | left  | center | right |
| ------ | :---- | :----: | ----: |
| n      | l     | c      | r     |
| N      | L     | C      | R     |


## Containter blocks

### Block quotes

> Blockquote
>
> second paragraph

#### Alerts

[2023-12-14](https://github.blog/changelog/2023-12-14-new-markdown-extension-alerts-provide-distinctive-styling-for-significant-content/)
GitHub 独自の拡張として
[Alerts](https://docs.github.com/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax#alerts)
が追加された。
GFMに正式に組み込まれたわけではない。

[Hugo] では
[blockquote render hooks](https://gohugo.io/render-hooks/blockquotes/)
を使えば同じことを実現できる。
[community#16925](https://github.com/orgs/community/discussions/16925)
の実装例を見るとHTMLは `<div class="markdown-alert">` 。

> [!NOTE]
> Useful information that users should know, even when skimming content.

> [!TIP]
> Helpful advice for doing things better or more easily.

> [!IMPORTANT]
> Key information users need to know to achieve their goal.

> [!WARNING]
> Urgent info that needs immediate user attention to avoid problems.

> [!CAUTION]
> Advises about risks or negative outcomes of certain actions.

> [!UNKNOWN]
> Unknown alert type shown without color or icon.

> [!type]+ title
> Hogo provides `AlertSign` and `AlertTitle` for the extended syntax like
> [Obsidian callouts](https://help.obsidian.md/callouts).
> This theme may also support it if it gets more popular.


### Lists

行頭に `- `, `+ `, `* ` を置くと `<ul><li>` の unordered list になる:
- Unordered list
- Nested level 1
  - Nested level 2
    - Nested level 3
      - Nested level 4
- Nested level 1

行頭に `1. ` か `1) ` を置くと `<ol><li>` の ordered list になる:
1. Ordered list
1. Nested level 1
   1. Nested level 2
      1. Nested level 3
         1. Nested level 4
1. Nested level 1

ひとつでも段落を持つアイテムがあると、全ての `<li>` の中に `<p>` が入る:
- list with paragraphs

  paragraph 2 in item 1
- item 2 in p\
  second line

リスト開始は段落継続よりも強いので空行を入れる必要はない。
逆に言えば、空行の有無でレイアウトを変えることもできない。


#### Task lists

[GFM `tasklist` extension](https://github.github.com/gfm/#task-list-items-extension-)

- [ ] Unchecked
- [x] Checked
  - nested ul
    - [ ] nested ul


`<ol>` でも行けるっぽい:
1. [ ] Unchecked
2. [x] Checked


#### Definition lists

GFMにも無いけど[Hugo]/Goldmarkではデフォルト有効な拡張。

term
: definition

`<dt>`
: `<dd>`


## Inlines

### Code spans

backtick/backquote 1つか2つで挟むと `` `inline code` ``.

`line endings
are converted
to spaces.`

` <- leading and trailing spaces are removed -> `
ただし両側ともにある場合のみで、削るのは1つずつ。

`space->    <-here`
内側のスペースの数は変わらずHTMLに出力されるけど、
ブラウザで1つに潰されるのがデフォルト。
文字通り表示させるためにCSSで
`code {white-space: pre-wrap;}` としておくのが推奨。


### Emphasis

- *emphasis, typically italic*
- **strong, typically bold**
  - 外側にスペースが無くても有効: left**strong**right
  - 内側のスペースで無効化: ** normal **
  - underscore `_` は外側にスペースが必要: left__normal__right
- ***strong emphasis, typycally bold italic***
- ~strikethrough~, ~~strikethrough~~:
  [GFM `strikethrough` extension](https://github.github.com/gfm/#strikethrough-extension-).
  チルダ1つでも `<del></del>` に変換するのがGFMの仕様かつHugoの初期設定だけど暴発しがち。
  安全のため無効化するか、2つ重ねでのみ変換するように設定する。


### Links

- [Link text](https://spec.commonmark.org/current/#links "optional title"):
  `[text](<href> "optional title")`
- [Link reference definitions]
  ```markdown
  [text]: <href> "optional title"
  ```
- Explicit autolink: `<https://example.com>` <https://example.com>
- [GFM `autolink` extension](https://github.github.com/gfm/#autolinks-extension-):
    `<` と `>` で挟まなくても
  `https://` とか `www.` とかで始まるURLらしきものを認識してリンク生成する。
  https://example.com
- Footnote `[^1]`[^1] `[^name]`[^name]

[Link reference definitions]:
  <https://spec.commonmark.org/current/#link-reference-definitions>
  "Optional title"

[^1]: `[^1]: ` GFMにも無いけどHugo/Goldmarkではデフォルト有効な拡張。
    <https://github.com/yuin/goldmark?tab=readme-ov-file#footnotes-extension>

[^name]:
    `[^name]: ` 番号でない名前もつけられるけど出力は通し番号。

    空行とスペース4つで複数段落も可能。


### Images

`![Alt text](</path/to/file.png> "Optional title")`
![Alt text](/not-found.png "Optional title")

![/favicon.svg](</favicon.svg>)

<img src="/favicon.svg" height=16>

[raw HTML](#raw-html) `<img src="..." alt="..." title="...">`
を書いた方が分かりやすいし他のattributeも使える。
ただし完全に同義というわけではなく、
段落まわりの挙動で違いが出ることがある。
すなわち `![]()` 記法のほうはどう書いても `<p>` の中に入るが、
生タグ独立で [HTML blocks](#html-blocks) になると囲われない。


### Raw HTML

[HTML blocks](#html-blocks) の条件を満たさない使い方をした生HTMLタグは、
Markdownの文脈を切らずにそのままHTMLに出力される。
つまり `<span>` とか `<kbd>` で修飾したりできる。

安全性のためにデフォルトでは無効化されている場合が多い。

- [GFM `tagfilter`](https://github.github.com/gfm/#disallowed-raw-html-extension-):
  `<title>`,
  `<textarea>`,
  `<style>`,
  `<xmp>`,
  `<iframe>`,
  `<noembed>`,
  `<noframes>`,
  `<script>`,
  `<plaintext>`.

- [Hugo `renderer.unsafe = false`](https://gohugo.io/configuration/markup/#rendererunsafe)


### Line breaks

Soft line breaks
: 改行ひとつはHTMLでも改行ひとつのままで、普通はスペースとしてレンダリングされる。
: 改行を改行として見せる設定も可能ではある。
  - CSSでレンダリングを変える:
    [`white-space: pre*`](https://developer.mozilla.org/docs/Web/CSS/white-space)
  - パーサーで `<br>` を入れる:
    [Hugo `renderer.hardWraps`](https://gohugo.io/configuration/markup/#rendererhardwraps)
: スペースを入れたくない言語のためのオプションがある場合もある。
  e.g., [goldmark CJK extension](https://github.com/yuin/goldmark#cjk-extension)

Hard line breaks
: 行末にバックスラッシュ `\` を置くと `<br>` が入って段落内改行。
: 行末にスペース2つ置いても同じ効果を得られることになってるけど、
  trailing whitespace は取り除かれるべきだし見えにくいので使わない。

