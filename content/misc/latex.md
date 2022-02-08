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

1.  <https://www.tug.org/mactex/> から
    `BasicTeX.pkg` を入手してインストール。
    あるいは
    `brew install --cask basictex`

1.  `/Library/TeX/texbin/` にパスを通す。
    基本的には `/etc/paths.d/TeX` 越しに自動的に設定されるはず。

1.  今後の `tlmgr` 操作で管理者権限を使わなくて済むようにパーミッション設定:
    `sudo chown -R $(whoami):admin /usr/local/texlive/`

1.  `tlmgr update --self --all` で諸々アップデート。

TeX Liveの新パージョンが年1回リリースされるので、そのときはインストールからやり直す。
(`tlmgr update` では更新できない)


### TeXをPDFにコンパイル

1.  ソースコードを書く

    ```latex
    \documentclass[a4paper]{article}
    \begin{document}
    Hello, World!
    \end{document}
    ```

1.  ターミナルからコンパイル:

    ```sh
    pdflatex hello.tex
    open hello.pdf
    ```

   MacTeXに含まれるTeXShopというアプリを使えば、
   原稿書きからコンパイルまで面倒を見てもらえるので、
   ターミナルにコマンドを打ち込む必要はないらしい。

#### `latexmk`

引用文献や図表への参照などを入れたりすると
`pdflatex` を一度実行するだけではPDFが完成せず、
コマンドを何度か繰り返し実行しなければならなくなる。
`latexmk` はそのへんをうまくお世話してくれる。

対象の `.tex` ファイルを明示的に指定してもいいし、
省略してディレクトリ内のすべてのファイルを対象にしてもいい。

```
latexmk -h    # print help
latexmk       # generate document
latexmk -pv   # preview document after generation
latexmk -pvc  # preview document and continuously update
latexmk -C    # clean up
```

設定を `~/.latexmkrc` とかに書いておける:

```perl
$pdf_mode = 1;
$pdflatex = 'pdflatex -file-line-error -halt-on-error -synctex=1 %O %S';
$lualatex = 'lualatex -file-line-error -halt-on-error -synctex=1 %O %S';
```

https://www.ctan.org/pkg/latexmk

#### SyncTeX

ソースコードとPDFの対応箇所を行き来するための仕組み。

- Skim to source: <kbd>shift-cmd-click</kbd>
- Atom to PDF: <kbd>ctrl-alt-s</kbd>
- VS Code to PDF: <kbd>alt-cmd-j</kbd>


### tlmgr でパッケージ管理

BasicTeXの場合は最小限のパッケージしか付いてこないので、必要なものを別途インストールする。
"By default, installing a package ensures that all dependencies of this package are fulfilled"
と言っていて実際tlmgr自体はそういう作りになっているが、
多くのパッケージで依存関係がちゃんと記述されていないので、
実際にはエラーを読みながら依存パッケージを手動でインストールすることになる。

```sh
tlmgr update --self --all
tlmgr search --global japanese
tlmgr search --global --file zxxrl7z
tlmgr info --list newtx
tlmgr install chktex latexmk
```

管理者権限不要の `--usermode` も用意されているが、
なんかうまくいかないので必要に迫られない限り使わないほうがよい。

パッケージをアンインストールしようと思って
`tlmgr uninstall some-package`
などとするとTeX Live全体が消えてしまうので注意。
正しくは `tlmgr remove some-package`


## 基本要素

### 数式

とりあえず
[Short Math Guide for LaTeX](https://mirror.ctan.org/info/short-math-guide)
(PDF)を読むべし。

-   基本的に [{amsmath}](https://www.ams.org/arc/resources/amslatex-about.html) を使う。
    数式環境のデファクトスタンダード。
    アメリカ数学会(AMS)が開発したらしいが、
    ガチ数学じゃなくても数式を書く場合はこれらしい。
    ```latex
    \usepackage{amssymb,amsmath}
    \usepackage[all,warning]{onlyamsmath}
    ```

-   Inline math:
    ```latex
    if $N_e u \ll 1$, then the population is monomorphic most of the time,
    ```

    `$ ... $` はTeXの古いやり方で、
    新しいLaTeXでは `\( ... \)` を用いるべし、という流れもある。
    けどダラーのほうが書きやすいし読みやすいので、
    しばらくは `chktex -n46` で様子を見る。

-   Display math:
    ```latex
    \begin{equation*}\label{eq:growth}
      N_t = N_0 e^{rt}
    \end{equation*}
    ```

    生TeXの `$$ ... $$` を使ってはいけない。
    `\[ ... \]` はOK。

-   改行を含む数式を等号で揃える

    ```latex
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
    ```latex
    \begin{equation}\label{eq:heaviside}
      H(x) = \begin{cases}
        & 0 \text{if $x \le 0$} \\
        & 1 \text{if $x > 0$}
      \end{cases}
    \end{equation}
    ```

-   記号:
    <https://www.ctan.org/tex-archive/info/symbols/comprehensive/>
    に網羅されてるけど、
    だいたい Short Math Guide for LaTeX にまとめられてるやつで足りるはず。
    -   カッコの大きさを変えたいときは `\left(` と `\right)` を使っておけば、
        前後のサイズに応じて自動的にうまいことやってくれる。
        e.g., 分数を挟むとか、カッコの入れ子とか
    -   "given that" を示す縦棒はパイプ記号 `|` ではなく
        `\mid` を使うのが正しいし適度なスペースが入って読みやすい。
        絶対値もパイプではなく `\lvert x \rvert` のようにする。
    -   斜体にしたくない文字を普通にするには `\mathrm dt` 。
        記号じゃないテキストには `\text{otherwise}` 。
        よく使われるやつは定義済み e.g., `\log`, `\exp`

### 図

```latex
\usepackage[final]{graphicx}
%%%
\begin{figure}
\includegraphics[width=\columnwidth]{awesome.pdf}
\caption{some description}%
\label{fig:awesome}
\end{figure}
```

`\begin{figure}[htbp]` のようにして配置場所の優先順位を指定できるが、
結局はいろんな兼ね合いでコンパイラが決める。
デフォルトは `[tbp]`。
複数指定した場合、どの順序で書いても下記の順に優先される。

- **h**ere: ソースと同じ位置に
- **t**op: 上部に
- **b**ottom: 下部に
- float **p**age: 図だけを含む独立したページに

幅の指定には `\textwidth`, `\columnwidth`, `\linewidth` などの変数が使える。
twocolumnの文書内で `\textwidth` まで広げたい場合は
アスタリスク付きの `\begin{figure*}` 環境を用いる。

ラベル `\label{}` は `\caption{}` 直後(または内部)じゃないとおかしくなるらしい。
間に改行やスペースが入るのもダメっぽいので `\caption{}` 行末に `%` を置いておく。
本文からは `\ref{fig:awesome}` のようにラベルで参照しておけば番号に置き換えられる。

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

```latex
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
表本体が `tabular` 環境。

```latex
\begin{table}
\caption{Parameters}%
\label{table:parameters}
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

```latex
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


### 引用

<https://en.wikibooks.org/wiki/LaTeX/More_Bibliographies>

1.  [Bibdesk](https://bibdesk.sourceforge.io/)
    などの文献管理アプリでbibtex形式の文献リストを作る。
    e.g., `mybibdata.bib`
1.  プリアンブルで [{natbib}](https://ctan.org/pkg/natbib) を呼び出し、スタイルを指定。
    ```latex
    \usepackage[authoryear,round]{natbib}
    \bibliographystyle{abbrvnat}
    ```
    スタイルとして `{natbib}` の `abbrvnat` や `unsrtnat` をそのまま使うことはまれで、
    各Journalの提供する、あるいは有志の作る `.bst` ファイルをダウンロードして使う。

    `latex makebst` コマンドから質問に答えて新規作成することも可能。
    後で間違いに気付いた場合、イチからやり直すより中間ファイルの
    `.dbj` ファイルを編集して `latex *.dbj` で `.bst` を生成すると早い。
    `.bst` そのものを読み解いて編集することも不可能ではない。

1.  LaTeX本文にcite keyを挿入。標準の `\cite` は使わない。

    comand | `authoryear,round` | `numbers`
    ---- | ---- | ----
    `\citet{hudson1987g}`        | Hudson et al. (1987)  | Hudson et al. [42]
    `\citep{hudson1987g}`        | (Hudson et al., 1987) | [42]
    `\citealt{hudson1987g}`      | Hudson et al. 1987    | Hudson et al. 42
    `\citealp{hudson1987g}`      | Hudson et al., 1987   | 42
    `\citeauthor{hudson1987g}`   | Hudson et al.         | Hudson et al.
    `\citeyear{hudson1987g}`     | 1987                  | 1987
    `\citenum{hudson1987g}`      | 42                    | 42
    `\citep[eq. 5]{hudson1987g}` | (Hudson et al., 1987, eq. 5) | [42, eq. 5]
    `\citep[see][]{hudson1987g}` | (see Hudson et al., 1987) | [see 42]
    `\citet*{hudson1987g}`       | Hudson, Kreitman, and Aguadé (1987) | Hudson, Kreitman, and Aguadé [42]
    `\citep*{hudson1987g}`       | (Hudson, Kreitman, and Aguadé, 1987) | [42]

    デフォルトは `authoryear,round,semicolon` だが
    `\bibliographystyle{}` 次第で括弧が勝手に `square` になったりする。

    "et al., 1987" の間のカンマを取りたい、といった微調整は
    `\bibpunct{(}{)}{;}{a}{}{,}` のようにする。

1.  最後の方に文献リストを挿入:
    ```latex
    \bibliography{mybibdata}
    ```

1.  元の `.tex` をコンパイルして `.aux` を生成
1.  `bibtex` に `.aux` を渡して `.bbl` を生成
1.  再び `.tex` をコンパイルすると `.bbl` を踏まえて `.aux` が更新される。
    (このときPDF出力すると、文献リストはできるけど引用部分はハテナ?になる)
1.  さらにもう1回コンパイルして完成

最初の2回は `pdflatex -draftmode` としてPDF出力を省略すると早い。
適切な `Makefile` を書いて自動化すると楽で、
[latexmk](#latexmk) を使うともっと楽ちん。


### 文字の修飾

- <https://en.wikibooks.org/wiki/LaTeX/Fonts#Font_styles>
- <https://en.wikibooks.org/wiki/LaTeX/Colors>

```latex
\emph{emphasis}
\textit{italic}
\textbf{bold}
\texttt{monospace}
{\huge huge text}
```

`\tiny`, `\scriptsize`, `\footnotesize`, `\small`,
`\normalsize`,
`\large`, `\Large`, `\LARGE`, `\huge`, `\Huge`

```latex
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


### 行番号

共著者や査読者との議論をスムーズにするため各行に番号を振る。
{amsmath} の `fleqn` オプションを使うと行番号が消えるとか、
`twocolumn` のときに `\pagewiselinenumbers{}` がページワイズにならないとか、
いろいろ不具合はあるものの
[{lineno}](https://www.ctan.org/pkg/lineno)
を使うしかなさそう。

```latex
\usepackage[mathlines,pagewise,switch]{lineno}
\renewcommand\linenumberfont{\normalfont\scriptsize\sffamily\color[gray]{0.5}}%
\setlength\linenumbersep{4truemm}

\linenumbers{}
```

## Tips

### ダメな使い方を警告してもらう

いろんなパッケージを使ったいろんな書き方がネット上にあふれているが、
中には古すぎたりするため避けたほうがよいものもある。
ファイルの先頭で `nag` を読み込むことで、
コンパイル時にそういうのを警告してもらえる。

```latex
\RequirePackage[l2tabu,orthodox]{nag}

\documentclass[a4paper]{article}  % これよりも前
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

```latex
\input{glyphtounicode.tex}
\pdfgentounicode=1
% あるいは
\usepackage{mmap}
```

## 日本語を使う

### LuaLaTeX

-   OS上にあるOTFフォントがそのまま使える
-   pdfTeXの後継として、今後のスタンダードと目される
-   [{luatexja}](https://www.ctan.org/pkg/luatexja) が精力的に開発されている
-   動作が遅い


### XeLaTeX

-   OS上にあるOTFフォントがそのまま使える
-   とにかく日本語入りでコンパイルできればいい、というのであればこれが早い
-   日本語に特化したツールは開発されていないので細かい制御がイマイチらしい


### upLaTeX

-   日本語を使えるように LaTeX を改良したもの
-   歴史が長いので日本語組版のための便利な道具が揃ってるらしいけど未来は無さそう


## フォント

Computer Modern
:   Knuth先生が作ったデフォルトフォント。

{lmodern} --- Latin Modern
:   Computer Modern の改良版。

{times}
:   ローマンとサンセリフにそれぞれ Times と Helvetica を割り当てる。
    数式は Computer Modern のまま。

{txfonts}
:   {times} の改良版？
    数式も Times にする。
    直接は使わない。

[{newtx}](https://www.ctan.org/pkg/newtx)
:   {txfonts} の後継で現役。
    本文と数式を別々に指定できる。
    `\usepackage[libertine]{newtxmath}` とすると Libertine を数式に使える。
    インストールするときは `newtx` だけでなく
    `txfonts` と `boondox` も入れないと
    `Unable to find TFM file` と怒られる。

[{newpx}](https://www.ctan.org/pkg/newpx)
:   {newtx}と同等の機能を美しいPalatinoで。
    {palatino}, {pxfonts}, {newtx},
    [{tex-gyre-pagella}](https://www.ctan.org/pkg/tex-gyre-pagella),
    [{tex-gyre-math-pagella}](https://www.ctan.org/pkg/tex-gyre-math-pagella) も入れておく。
    **TeX Gyre Pagella** はOpenType志向のPalatinoクローン。

[{libertinus}](https://www.ctan.org/pkg/libertinus)
:   美しい[Linux Libertine](https://www.ctan.org/pkg/libertine)の後継プロジェクト。
    type1もOTFも数式もサポートしていて便利だがひと回り小さいことに注意。
    使うときは `\usepackage{libertinus}` でよしなにやってくれるらしいが
    依存パッケージのインストールは例によって手動:
    `libertibnus libertinus-fonts libertinus-type1 libertinust1math libertinus-otf`

LuaTeX/XeTeXならOSのフォントをフルネームで指定して使えるが、
共同執筆とかを考えるとTeX Liveパッケージやプリセットに含まれるものを使うのが安全。

TeX Liveから入れたフォントをOSに認識させるにはシムリンクを張るだけ:
`ln -s /Library/TeX/Root/texmf-dist/fonts/opentype ~/Library/Fonts/texlive-opentype`


```latex
\usepackage{amssymb,amsmath} % must be called ahead of mathspec
\usepackage[all,warning]{onlyamsmath}
\usepackage{iftex}
\iftutex
  \usepackage[math-style=TeX,bold-style=TeX]{unicode-math}
  \usepackage[no-math]{fontspec}
  \setmainfont{TeX Gyre Pagella}
  \setmathfont{TeX Gyre Pagella Math}
  \setsansfont{TeX Gyre Heros}
  % \usepackage{luatexja}
  % \usepackage[hiragino-pron,scale=0.92,deluxe,jis2004,match,nfssonly]{luatexja-preset}
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
