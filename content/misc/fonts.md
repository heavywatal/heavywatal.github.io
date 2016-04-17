+++
title = 'Fonts'
tags = ["writing"]
[menu.main]
  parent = "misc"
+++

使えるフリーフォント

[Noto Fonts](http://www.google.com/get/noto/)
:   -   "No tofu" を目指してGoogleが開発している多言語対応フォント。
    -   Droid fontを継承してるっぽい
    -   [Noto Sans CJK](http://www.google.com/get/noto/help/cjk/)
        の日本語部分はAdobeと共同開発で、
        同じものを違う名前で公開している。
        英数字は Source Sans Pro がリンクされているので、
        どっちかといえば Source Han Sans を使ったほうがいいのかな。
    -   Noto Serif はディスプレイ上での視認性がとてもよい(があまり美しくはない)

[Source Sans Pro](https://github.com/adobe-fonts/source-sans-pro)
:   -   Adobe初のオープンソースフォント
    -   小文字の l が曲がってて良いけど全体的に狭苦しい
    -   Googleと共同開発の
        [Source Han Sans](https://github.com/adobe-fonts/source-han-sans)
        (源ノ角ゴシック)は素晴らしい日本語フォント。
        ウェイトもたくさんある。
    -   Source Han Code JP は普通の固定幅フォントよりも幅広で、
        和文と欧文の幅が3:2になるよう調整されている。
        俺は2:1のUbuntu Monoのほうが好き。

[Ubuntu](http://font.ubuntu.com)
:   -   スクリーン上での視認性を重視してデザインされたフォントで、
        個性的でありながら見やすい
    -   エルの小文字 l が曲がってて良い
    -   Ubuntu Mono は特にプログラミング用のフォントとして最適
    -   普通の欧文 `monospace` よりもひとまわり小さく、
        Osaka-Mono など日本語の半角文字と同じ幅になるのもポイント

[Open Sans](https://www.google.com/fonts/specimen/Open+Sans)
:   -   すっきりゆったりニュートラルで見やすい
    -   ファイルサイズも小さいのでウェブフォントとして使いやすい

[Linux Libertine](http://www.linuxlibertine.org)
:   -   印刷に耐えうるクオリティを目指して開発されたカッコいいセリフ体
    -   近代 LaTeX ではこれがスタンダードになっていくのでは
    -   兄弟の Linux Biolinum も半セリフみたいな感じでカッコいい
    -   ひとまわり小さいのでほかのフォントとバランスとるのが難しい

[Roboto](https://github.com/google/roboto)
:   -   かっちりコンパクト、いかにもシステム向けなフォント
    -   でもちょっと c とかの切れ目が小さすぎる気がする

------------------------------------------------------------------------

[DejaVu](http://dejavu-fonts.org/)
:   -   Bitstream Vera の後を継いで多言語対応が進められている
    -   Linuxにはだいたい入ってる
    -   DejaVu Sans は Verdana のように幅が広く、スクリーン上で読みやすい
    -   DejaVu Sans Mono はプログラミングに使いやすい

[Liberation](https://fedorahosted.org/liberation-fonts/)
:   -   Fedora が Arial, Times New Roman, Courier New からの解放を目指して作ってる
    -   英数字はいいけど、グリフ数は少ない
    -   でも Liberation Mono は Courier New というより Courier な太さ?

[GNU FreeFont](http://www.gnu.org/software/freefont/)
:   -   Helvetica, Times, Courier の置き換えを狙ったフォント
    -   でも FreeMono は Courier というより Courier New な細さ?

[IPA](http://ossipedia.ipa.go.jp/ipafont/)
:   -   十分使える日本語フリーフォントだったが、
        源ノ角ゴシックの登場により役目を終えた感・・・？
    -   公式のIPAフォントは更新サイクルが遅いので
        [Takao](https://launchpad.net/takao-fonts)
        という名前でUbuntuコミュニティが保守している

[MigMix](https://mix-mplus-ipa.sourceforge.jp/migmix/)
:   -   視認性の高いM+フォントに足りない漢字をIPAゴシックから補完した合成フォント。

## serif

<table summary="" border=0 cellpadding=3px cellspacing=3px class="tabel"><tr style="font-family: 'Times', sans-serif" bgcolor="#eeeeee">
  <td rowspan="2">Times</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Times', sans-serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Times New Roman', sans-serif" bgcolor="#eeeeee">
  <td rowspan="2">Times New Roman</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Times New Roman', sans-serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Garamond', sans-serif" bgcolor="#eeeeee">
  <td rowspan="2">Garamond</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Garamond', sans-serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Georgia', sans-serif" bgcolor="#eeeeee">
  <td rowspan="2">Georgia</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Georgia', sans-serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Palatino', sans-serif" bgcolor="#eeeeee">
  <td rowspan="2">Palatino</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Palatino', sans-serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Linux Libertine', sans-serif" bgcolor="#eeeeee">
  <td rowspan="2">Linux Libertine</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Linux Libertine', sans-serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Noto Serif', sans-serif" bgcolor="#eeeeee">
  <td rowspan="2">Noto Serif</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Noto Serif', sans-serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Lora', sans-serif" bgcolor="#eeeeee">
  <td rowspan="2">Lora</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Lora', sans-serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
</table>
## sans-serif

<table summary="" border=0 cellpadding=3px cellspacing=3px class="tabel"><tr style="font-family: 'Helvetica Neue', serif" bgcolor="#eeeeee">
  <td rowspan="2">Helvetica Neue</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Helvetica Neue', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Helvetica', serif" bgcolor="#eeeeee">
  <td rowspan="2">Helvetica</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Helvetica', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Arial', serif" bgcolor="#eeeeee">
  <td rowspan="2">Arial</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Arial', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Lucida Grande', serif" bgcolor="#eeeeee">
  <td rowspan="2">Lucida Grande</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Lucida Grande', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Lucida Sans Unicode', serif" bgcolor="#eeeeee">
  <td rowspan="2">Lucida Sans Unicode</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Lucida Sans Unicode', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Lucida Sans', serif" bgcolor="#eeeeee">
  <td rowspan="2">Lucida Sans</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Lucida Sans', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Tahoma', serif" bgcolor="#eeeeee">
  <td rowspan="2">Tahoma</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Tahoma', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Calibri', serif" bgcolor="#eeeeee">
  <td rowspan="2">Calibri</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Calibri', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Verdana', serif" bgcolor="#eeeeee">
  <td rowspan="2">Verdana</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Verdana', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'DejaVu Sans', serif" bgcolor="#eeeeee">
  <td rowspan="2">DejaVu Sans</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'DejaVu Sans', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Ubuntu', serif" bgcolor="#eeeeee">
  <td rowspan="2">Ubuntu</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Ubuntu', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Open Sans', serif" bgcolor="#eeeeee">
  <td rowspan="2">Open Sans</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Open Sans', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Noto Sans', serif" bgcolor="#eeeeee">
  <td rowspan="2">Noto Sans</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Noto Sans', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Source Sans Pro', serif" bgcolor="#eeeeee">
  <td rowspan="2">Source Sans Pro</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Source Sans Pro', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Roboto', serif" bgcolor="#eeeeee">
  <td rowspan="2">Roboto</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Roboto', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Linux Biolinum', serif" bgcolor="#eeeeee">
  <td rowspan="2">Linux Biolinum</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Linux Biolinum', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
</table>
## monospace

<table summary="" border=0 cellpadding=3px cellspacing=3px class="tabel"><tr style="font-family: 'Courier', serif" bgcolor="#eeeeee">
  <td rowspan="2">Courier</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Courier', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Courier New', serif" bgcolor="#eeeeee">
  <td rowspan="2">Courier New</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Courier New', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Monaco', serif" bgcolor="#eeeeee">
  <td rowspan="2">Monaco</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Monaco', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Menlo', serif" bgcolor="#eeeeee">
  <td rowspan="2">Menlo</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Menlo', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'DejaVu Sans Mono', serif" bgcolor="#eeeeee">
  <td rowspan="2">DejaVu Sans Mono</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'DejaVu Sans Mono', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Source Code Pro', serif" bgcolor="#eeeeee">
  <td rowspan="2">Source Code Pro</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Source Code Pro', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Source Han Code JP', serif" bgcolor="#eeeeee">
  <td rowspan="2">Source Han Code JP</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Source Han Code JP', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Ubuntu Mono', serif" bgcolor="#eeeeee">
  <td rowspan="2">Ubuntu Mono</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Ubuntu Mono', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Osaka-Mono', serif" bgcolor="#eeeeee">
  <td rowspan="2">Osaka-Mono</td>
  <td>ABCDEFGHIJKLMNOPQRSTUVWXYZ</td><td>0123456789</td></tr>
<tr style="font-family: 'Osaka-Mono', serif">
  <td>abcdefghijklmnopqrstuvwxyz</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
</table>
## Japanese

<table summary="" border=0 cellpadding=3px cellspacing=3px class="tabel"><tr style="font-family: 'Hiragino Kaku Gothic ProN', serif" bgcolor="#eeeeee">
  <td rowspan="2">Hiragino Kaku Gothic ProN</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'Hiragino Kaku Gothic ProN', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Hiragino Kaku Gothic Pro', serif" bgcolor="#eeeeee">
  <td rowspan="2">Hiragino Kaku Gothic Pro</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'Hiragino Kaku Gothic Pro', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'YuGothic', serif" bgcolor="#eeeeee">
  <td rowspan="2">YuGothic</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'YuGothic', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Source Han Sans', serif" bgcolor="#eeeeee">
  <td rowspan="2">Source Han Sans</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'Source Han Sans', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Osaka-Mono', serif" bgcolor="#eeeeee">
  <td rowspan="2">Osaka-Mono</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'Osaka-Mono', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'TakaoExGothic', serif" bgcolor="#eeeeee">
  <td rowspan="2">TakaoExGothic</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'TakaoExGothic', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'TakaoPGothic', serif" bgcolor="#eeeeee">
  <td rowspan="2">TakaoPGothic</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'TakaoPGothic', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'TakaoGothic', serif" bgcolor="#eeeeee">
  <td rowspan="2">TakaoGothic</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'TakaoGothic', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'MigMix 1P', serif" bgcolor="#eeeeee">
  <td rowspan="2">MigMix 1P</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'MigMix 1P', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'MigMix 2P', serif" bgcolor="#eeeeee">
  <td rowspan="2">MigMix 2P</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'MigMix 2P', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'MigMix 1M', serif" bgcolor="#eeeeee">
  <td rowspan="2">MigMix 1M</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'MigMix 1M', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'MigMix 2M', serif" bgcolor="#eeeeee">
  <td rowspan="2">MigMix 2M</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'MigMix 2M', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'VL PGothic', serif" bgcolor="#eeeeee">
  <td rowspan="2">VL PGothic</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'VL PGothic', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'VL Gothic', serif" bgcolor="#eeeeee">
  <td rowspan="2">VL Gothic</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'VL Gothic', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'MS PGothic', serif" bgcolor="#eeeeee">
  <td rowspan="2">MS PGothic</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'MS PGothic', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'MS Gothic', serif" bgcolor="#eeeeee">
  <td rowspan="2">MS Gothic</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'MS Gothic', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Hiragino Mincho ProN', serif" bgcolor="#eeeeee">
  <td rowspan="2">Hiragino Mincho ProN</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'Hiragino Mincho ProN', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'Hiragino Mincho Pro', serif" bgcolor="#eeeeee">
  <td rowspan="2">Hiragino Mincho Pro</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'Hiragino Mincho Pro', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'YuMincho', serif" bgcolor="#eeeeee">
  <td rowspan="2">YuMincho</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'YuMincho', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'TakaoExMincho', serif" bgcolor="#eeeeee">
  <td rowspan="2">TakaoExMincho</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'TakaoExMincho', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'TakaoPMincho', serif" bgcolor="#eeeeee">
  <td rowspan="2">TakaoPMincho</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'TakaoPMincho', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'TakaoMincho', serif" bgcolor="#eeeeee">
  <td rowspan="2">TakaoMincho</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'TakaoMincho', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'MS PMincho', serif" bgcolor="#eeeeee">
  <td rowspan="2">MS PMincho</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'MS PMincho', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
<tr style="font-family: 'MS Mincho', serif" bgcolor="#eeeeee">
  <td rowspan="2">MS Mincho</td>
  <td>噌祇逢餅鯖鰯</td>
  <td>１２３４５６７９</td><td>！％＆（）＊＋、ー。・：；＜＝＞？＠「￥」｛｝〜</td>
<tr style="font-family: 'MS Mincho', serif">
  <td>わたるメタル</td>
  <td>0123456789</td><td>!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~</td>
</tr>
</table>
