+++
title = 'LaTeX'
tags = ["writing"]
aliases = ["/dev/tex.html"]
[menu.main]
  parent = "misc"
+++

## 基本操作

### TeX Live をMacにインストール

とにかく動く環境が欲しい人は、素直にフルのMacTeXをインストールするだけ。
ただし2GB以上の大きなインストーラをダウンロードする必要があるので注意。
BibDesk, LaTeXiT, TeX Live Utility, TeXShop などのGUIアプリが不要で、
必要なパッケージをコマンドラインからインストールできる人は下記の手順で小さくインストールできる。

1.  <http://www.tug.org/mactex/> から
    `BasicTeX.pkg` を入手してインストール。
    あるいは
    `brew cask install basictex`

2.  `/Library/TeX/texbin/` にパスを通す。
    基本的には `/etc/paths.d/TeX` 越しに自動的に設定されるはず。

3.  今後の `tlmgr` 操作で管理者権限を使わなくて済むようにパーミッション設定:
    `sudo chown -R $(whoami):wheel /usr/local/texlive/`

4.  `tlmgr update --self --all` で諸々アップデート。

TeX Liveの新パージョンが年1回リリースされるので、そのときはインストールからやり直す。
(`tlmgr update` では更新できない)


### TeXをPDFにコンパイル

1.  ソースコードを書く

    ```tex
    \documentclass{article}
    \begin{document}
    Hello, World!
    \end{document}
    ```

2.  ターミナルからコンパイル:

    ```sh
    pdflatex hello.tex
    open hello.pdf
    ```

   MacTeXに含まれるTeXShopというアプリを使えば、
   原稿書きからコンパイルまで面倒を見てもらえるので、
   ターミナルにコマンドを打ち込む必要はないらしい。


### tlmgr でパッケージ管理

BasicTeXの場合は最小限のパッケージしか含まれていないので、
大抵いくつか別途インストールする必要がある。
"By default, installing a package ensures that all dependencies of this package are fulfilled"
と言っているが嘘で、実際には依存パッケージも手動でインストールする必要がある。

```sh
tlmgr update --self --all
tlmgr search --global japanese
tlmgr search --global --file zxxrl7z
tlmgr info --list newtx
tlmgr install chktex
```

管理者権限不要の `--usermode` も用意されているが、
なんかうまくいかないので必要に迫られない限り使わないほうがよい。

{{%div class="warning"%}}
パッケージをアンインストールしようと思って
``tlmgr uninstall some-package``
などとするとTeX Live全体が消えてしまうので注意。
正しくは ``tlmgr remove some-package``
{{%/div%}}

## 基本要素

### 数式

とりあえず
[Short Math Guide for LaTeX](ftp://ftp.ams.org/pub/tex/doc/amsmath/short-math-guide.pdf)
(PDF)を読むべし。

-   基本的に [{amsmath}](http://www.ams.org/publications/authors/tex/amslatex) を使う。
    数式環境のデファクトスタンダード。
    アメリカ数学会(AMS)が開発したらしいが、
    ガチ数学じゃなくても数式を書く場合はこれらしい。
    ```tex
    \usepackage{amsmath}
    \usepackage[all, warning]{onlyamsmath}
    ```

-   Inline math:
    ```tex
    if $N_e u \ll 1$, then the population is monomorphic most of the time,
    ```

    `$...$` はTeXの古いやり方で、
    新しいLaTeXでは `\(...\)` を用いるべし、という流れもある。
    けどダラーのほうが書きやすいし読みやすいので、
    しばらくは `chktex -n46` で様子を見る。

-   Display math:
    ```tex
    \begin{equation*}\label{eq:growth}
      N_t = N_0 e^{rt}
    \end{equation*}
    ```

-   改行を含む数式を等号で揃える

    ```tex
    \begin{equation}\label{eq:growth}
    \begin{split}
      N_t &= N_0 e^{rt} \\
          &= N_0 \lambda^t
    \end{split}
    \end{equation}

    \begin{align}
      S &= 4 \pi r^2\label{surface}\\
      V &= \frac {4 \pi r^3} 3\label{volume}\\
    \end{align}
    ```

    `split` 環境では全体でラベルが1つ、
    `align` 環境では行ごとにラベルが生成されるので、
    用途が少し違う。
    古い `eqnarray` 環境でも似たようなことはできるが、
    スペースとかに問題あるらしく非推奨。

-   場合分け
    ```tex
    \begin{equation}\label{eq:heaviside}
      H(x) = \begin{cases}
        & 0 \text{if $x \le 0$} \\
        & 1 \text{if $x > 0$}
      \end{cases}
    \end{equation}
    ```

-   記号:
    <http://www.ctan.org/tex-archive/info/symbols/comprehensive/>
    に網羅されてるけど、
    だいたい Short Math Guide for LaTeX にまとめられてるやつで足りるはず。
    -   "given that" を示す縦棒はパイプ記号 `|` ではなく
        `\mid` を使うのが正しいし適度なスペースが入って読みやすい。
    -   カッコの大きさを変えたいときは `\left(` と `\right)` を使っておけば、
        前後のサイズに応じて自動的にうまいことやってくれる。
        e.g., 分数を挟むとか、カッコの入れ子とか
    -   斜体にしたくない文字を普通にするには `\mathrm dt` 。
        記号じゃないテキストには `\text{otherwise}` 。
        よく使われるやつは定義済み e.g., `\log`, `\exp`

### 図

```tex
\usepackage[final]{graphicx}
%%%
\begin{figure}[htbp]
\includegraphics[width=100mm]{great.pdf}
\caption{some description
}\label{fig:great}
\end{figure}
```

`[hb]` のようにして場所を指定できるが、
結局はいろんな兼ね合いでコンパイラが位置を決める。
デフォルトは `[tbp]`。
複数指定した場合、どの順序で書いても下記の順に優先される。

- `h` here: ソースと同じ位置に
- `t` top: 上部に
- `b` bottom: 下部に
- `p` float page: 図だけを含む独立したページに

ラベル `\label{}` は `\caption{}` 直後が良いらしく、
スペースさえもchktex警告対象なので波カッコ前で改行するという苦肉の策。

キャプションをカスタマイズするには
[{caption}](https://www.ctan.org/pkg/caption)
の `\captionsetup{...}` を用いる。

ひとつの領域に複数の図を貼るには
[{subfig}](https://www.ctan.org/pkg/subfig)
の `\subfloat[caption]{filename}` を用いる。

[{epstopdf}](https://www.ctan.org/pkg/epstopdf)
でEPSを取り込もうとすると(e.g., PLOS)、
ファイルもパッケージも揃ってるはずなのに
`! Package pdftex.def Error: File '*-eps-converted-to.pdf' not found.`
というエラーが出る。
変換プログラム本体である `ghostscript` をHomebrewか何かで入れる必要がある。

GIFアニメをそのまま埋め込むことはできないので、
[{animate}](https://www.ctan.org/pkg/animate)で連番PNGを読み込む。

```tex
% tlmgr install animate media9 ocgx2
\usepackage[autoplay,final,controls=all,type=png]{animate}

\begin{document}
% \animategraphics[options]{frame rate}{file basename}{first}{last}
\animategraphics[scale=0.5]{10}{dir/basename_}{1}{9}
\end{document}
```


### 表

<https://en.wikibooks.org/wiki/LaTeX/Tables>

キャプションやラベルなどをまとめるのが `table` 環境、
表本体が `tabular` 環境

```tex
\begin{table}
\caption{Parameters
}\label{table:parameters}
\begin{tabular}{rl}
  \toprule
  Symbol & Description \\
  \midrule
  $\mu$ & mutation rate per division \\
  $N_0$ & initial population size \\
  \bottomrule
\end{tabular}
\end{table}
```

罫線は標準の `\hline` だと上下スペースが狭苦しいので、
[{booktabs}](https://www.ctan.org/pkg/booktabs)
の `toprule`, `midrule`, `bottomrule` を使う。

ページに合わせて幅をうまいことやるには [{tabulary}](https://www.ctan.org/pkg/tabulary) 。

複数ページにまたがる表を作るには標準 [{longtable}](https://www.ctan.org/pkg/longtable) やその派生。


### 箇条書き

[{enumitem}](https://www.ctan.org/pkg/enumitem)
を使うといろいろなオプションが設定可能になる。

```tex
\usepackage{enumitem}
%%%

\begin{itemize}
  \item Judas Priest
  \item Iron Maiden
\end{itemize}

\begin{enumerate}[nosep,leftmargin=*]
  \item Judas Priest
  \item Iron Maiden
\end{enumerate}

\begin{description}[nextline]
  \item[key1] value1
  \item[key2] value2
\end{description}
```

<http://konoyonohana.blog.fc2.com/blog-entry-58.html>

### 引用

<https://en.wikibooks.org/wiki/LaTeX/More_Bibliographies>

1.  [Bibdesk](http://bibdesk.sourceforge.net/)
    などの文献管理アプリでbibtex形式の文献リストを作る。
    e.g., `mybibdata.bib`
2.  コマンドにcite keyを入れて本文に挿入。
    このとき標準の `\cite` ではなく `{natbib}` のものを使う。
    ```tex
    \usepackage[authoryear,round,sort&compress]{natbib}
    %%%

    \citet{hudson1987g}         # Hudson et al. (1987)
    \citep{hudson1987g}         # (Hudson et al. 1987)
    \citep*{hudson1987g}        # (Hudson, Kreitman and Aguadé, 1987)
    \citep[eq. 5]{hudson1987g}  # (Hudson et al. 1987, eq. 5)
    \citep[see][]{hudson1987g}  # (see Hudson et al. 1987)
    ```

3.  最後の方に文献リストを挿入:
    ```tex
    \bibliographystyle{abbrvnat}
    \bibliography{mybibdata}
    ```

    スタイルを規定する `bst` ファイルはだいたい各Journalで提供してくれる。

4.  元の `tex` をコンパイルして `aux` を生成
5.  `bibtex` に `aux` を渡して `bbl` を生成
6.  再び `tex` をコンパイルすると `bbl` を踏まえて `aux` が更新される。
    (このときPDF出力すると、文献リストはできるけど引用部分はハテナ?になる)
7.  さらにもう1回コンパイルして完成

最初の2回は `pdflatex -draftmode` としてPDF出力を省略すると早い。
適切な `Makefile` を書いて自動化するべし。
[`latexmk`](https://www.ctan.org/pkg/latexmk) を使うともっと楽ちん。


### 文字の修飾

- <https://en.wikibooks.org/wiki/LaTeX/Fonts#Font_styles>
- <https://en.wikibooks.org/wiki/LaTeX/Colors>

```tex
\emph{emphasis}
\textit{italic}
\textbf{bold}
\texttt{monospace}
{\huge huge text}
```

`\tiny`, `\scriptsize`, `\footnotesize`, `\small`,
`\normalsize`,
`\large`, `\Large`, `\LARGE`, `\huge`, `\Huge`

```tex
\usepackage[normalem]{ulem}  % \uline{}, \sout{}
\usepackage{color}           % \textcolor{}
\usepackage{soul}            % \hl{} using {color}
%%%

\uline{underlined text}
\sout{strikethrough}

\textcolor{red}{colored text}
{\text{red} colored text}

\hl{highlighted text}
```

[{ulem}](https://www.ctan.org/pkg/ulem) は
`[normalem]` オプションを付けて読まないと
`\emph` が下線に変更されてしまうので注意。

[{soul}](https://www.ctan.org/pkg/soul) のドキュメントによれば
`\hl{環境}` に `$数式$` を入れられるはずだが
"Extra }, or forgotten $" というエラーで弾かれる。


## Tips

### ダメな使い方を警告してもらう

いろんなパッケージを使ったいろんな書き方がネット上にあふれているが、
中には古すぎたりするため避けたほうがよいものもある。
ファイルの先頭で `nag` を読み込むことで、
コンパイル時にそういうのを警告してもらえる。

```tex
\RequirePackage[l2tabu, orthodox]{nag}

\documentclass{article}  % これよりも前
```

`chktex` コマンドを使えばコンパイルよりも手軽にチェックできる。
[Atom]({{< relref "atom.md" >}}) に
[linter-chktex](https://atom.io/packages/linter-chktex)
を入れれば編集中のファイルの警告箇所を随時ハイライトしてもらえる。


### ligature問題

表示の美しさという点でリガチャは素晴らしいけど、
PDF内の検索やPDFからのコピペ時に問題が発生する。
例えば `fi` が合字になるため `definition` が検索でひっかからない。
`definition` をコピペすると `de nition` になってしまう。

次のコードをプリアンブルの頭の方に記述するといいらしいが、うまく機能しない。。。

```tex
\input{glyphtounicode.tex}
\pdfgentounicode=1
% あるいは
\usepackage{mmap}
```

## 日本語を使う

### LuaLaTeX

- pdfTeXの後継として、今後のスタンダードと目される
- かなり動作が遅い
- [`luatexja`](https://www.ctan.org/pkg/luatexja) が精力的に開発されている
- 最初にフォント関連の問題に遭遇するかも
    - 依存パッケージ:
      `luaotfload`, `adobemapping`, `ctablestack`, `ipaex`
    - `bad argument #1 to 'open'` などと怒られる問題は、
      必要な CMap が LuaTeX から見えていないのが原因

      ```sh
      ## キャッシュを再構築
      % rm -rf ~/Library/texlive/2016basic/texmf-var/luatex-cache/
      % luaotfload-tool -u -v

      ## 見えてるか確認
      % kpsewhich -format=cmap UniJIS2004-UTF32-H
      % kpsewhich -format=cmap UniJIS2004-UTF32-V
      % kpsewhich -format=cmap Adobe-Japan1-UCS2
      ```

### XeLaTeX

-   とにかく日本語入りでコンパイルできればいい、というのであればこれが一番楽
-   BasicTeX に含まれているのでそのまま使える
-   OS上にあるOTFやTTFなどのフォントがそのまま使える
-   日本語に特化したツールは開発されていないので細かい制御がイマイチらしい


### upLaTeX

-   日本語を使えるように LaTeX を改良したもの
-   歴史が長いので日本語組版のための便利な道具が揃っている
-   BasicTeX には含まれていないのでいろいろインストールが必要

    `ptex`, `uptex`
    :   `platex` などを含む

    `jfontmaps`
    :   フォント埋め込みツール `kanji-config-updmap` など

    `jsclasses`
    :   `\documentclass{jsarticle}` などを含む

    `japanese-otf`, `japanese-otf-uptex`
    :   `\usepackage[uplatex,deluxe,jis2004]{otf}` でヒラギノが使える

    `ptex2pdf`
    :   `pdflatex` のように一発でPDFを作るための便利スクリプト。
        例えば `uplatex` を使ってコンパイルするには:

            % ptex2pdf -u -l main.tex

-   フォントマップの設定をする:

        % kanji-config-updmap status
        % kanji-config-updmap auto
        % kanji-config-updmap noEmbed
        % kanji-config-updmap hiragino

## フォント

Computer Modern
:   Knuth先生が作ったデフォルトフォント。

`\usepackage{lmodern}` --- Latin Modern
:   Computer Modern の改良版。

`\usepackage{times}`
:   ローマンとサンセリフにそれぞれ Times と Helvetica を割り当てる。
    数式は Computer Modern のまま。

`\usepackage{txfonts}`
:   `times` の改良版？
    数式も Times にする。
    既にメンテナンスは放棄されている。

[`newtx`](https://www.ctan.org/pkg/newtx)
:   `txfonts` の後継で現役。
    本文と数式を別々に指定できる。
    `\usepackage[libertine]{newtxmath}` とすると Libertine を数式に使える。
    インストールするときは `newtx` だけでなく
    `txfonts` と `boondox` も入れないと
    `Unable to find TFM file` と怒られる。

[`newpx`](https://www.ctan.org/pkg/newpx)
:   `newtx`と同等の機能を美しいPalatinoで。
    `palatino`, `pxfonts`,
    [`tex-gyre-pagella`](https://www.ctan.org/pkg/tex-gyre-pagella),
    [`tex-gyre-math-pagella`](https://www.ctan.org/pkg/tex-gyre-math-pagella) も入れておく。
    **TeX Gyre Pagella** はOpenType志向のPalatinoクローン。

[`libertine`](https://www.ctan.org/pkg/libertine) --- Linux Libertine
:   印刷用途でも通用するよう作られた美しいフリーフォント。
    `fontaxes` もインストールする必要がある。

XeTeX なら OS のフォントをフルネームで指定して使える

```tex
\usepackage{amsthm, amsmath} % must be called ahead of mathspec
\usepackage[all,warning]{onlyamsmath}
\ifxetex
  \usepackage{mathspec} % must be called ahead of fontspec
  \usepackage[math-style=TeX,bold-style=TeX]{unicode-math}
  \usepackage[no-math]{fontspec}
  \setmainfont{TeXGyrePagellaX}
  \setmathfont{texgyrepagella-math.otf}
  % \usepackage{xeCJK}
  % \setCJKmainfont[Scale=0.9,BoldFont=Hiragino Mincho ProN W6]
  %                                   {Hiragino Mincho ProN W3}
  % \setCJKsansfont[Scale=0.9,BoldFont=Hiragino Kaku Gothic ProN W6]
  %                                   {Hiragino Kaku Gothic ProN W3}
  % \setCJKmonofont[Scale=0.9,BoldFont=Hiragino Kaku Gothic ProN W6]
  %                                   {Hiragino Kaku Gothic ProN W3}
\else
  \usepackage[T1]{fontenc}
  \usepackage[utf8]{inputenc}
  \usepackage{newpxtext}
  \usepackage{newpxmath}
  \usepackage{textcomp}
  % \usepackage[uplatex,deluxe,jis2004]{otf}
\fi
```

フォント関連をいじったあと明示的にマップを更新するには `updmap-sys`
