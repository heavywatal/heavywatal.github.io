+++
title = 'PAML'
[menu.main]
  parent = "bio"
+++

**Phylogenetic Analysis by Maximum Likelihood**

<http://abacus.gene.ucl.ac.uk/software/paml.html>

<https://www.biostars.org/p/5817/>

## Installation

1.  [作者のウェブサイト](http://abacus.gene.ucl.ac.uk/software/paml.html)
    からソースコード `paml*.*.tgz` をダウンロードして展開:

        % wget -O- http://abacus.gene.ucl.ac.uk/software/paml4.8.tgz | tar xz

2.  中に `cd` してWindows用の実行ファイル `*.exe` を削除:

        % cd paml4.8/
        % rm -f bin/*.exe

3.  `src/Makefile` を適宜書き換えてコンパイル:

        % make -C src/

4.  できあがった実行ファイルだけを `bin/` に移動:

        % mv src/*(*) bin/

5.  このフォルダごと `/usr/local/paml` らへんに移動:

        % cd ..
        % sudo mv paml4.8 /usr/local/
        % sudo ln -s paml4.8 /usr/local/paml

6.  パスを通す。例えば `.zshenv` とかに以下のように記述:

        PATH=${PATH}:/usr/local/paml/bin

{{%div class="note"%}}
Macなら [Homebrew]({{< relref "mac/homebrew.md" >}}) で入れることもできる。
{{%/div%}}

## Sub-programs

`baseml`

`basemlg`

`chi2`

`codeml`

`evolver`

`infinitesites`

`mcmctree`

`pamp`

`yn00`
:   配列ペア間の `$d_N$` と `$d_S$` を推定。
    [Yang and Nielsen (2000) MBE](http://www.ncbi.nlm.nih.gov/pubmed/10666704)

## Citation

[Yang (1997) Comput Appl Biosci](http://www.ncbi.nlm.nih.gov/pubmed/9367129)

[Yang (2007) MBE](http://www.ncbi.nlm.nih.gov/pubmed/17483113)
