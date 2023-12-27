+++
title = 'BioPython'
subtitle = "Tools for biological computation"
tags = ["python"]
[menu.main]
  parent = "python"
+++

<https://biopython.org>
<https://biopython.org/DIST/docs/api/Bio-module.html>

## `SeqIO`

-   <https://biopython.org/wiki/SeqIO>
-   <https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec51>
-   <https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html>

`read(file, format, alphabet=None)`
:   1配列しか含まないファイルを読んで `SeqRecord` を返す

`parse(file, format, alphabet=None)`
:   複数の配列を含むファイルを読んで `SeqRecord` のイテレータを返す

`index(filename, format, alphabet=None, key_function=None)`
:   `id` をキーとした辞書likeオブジェクトを返す

`write(sequences, file, format)`

`to_dict(sequences, key_function=None)`

`convert(in_file, in_format, out_file, out_format, alphabet=None)`

    from Bio import SeqIO

    with open('beer.fasta', 'r') as fin:
        for record in SeqIO.parse(fin, 'fasta'):
            print(record.id)
            print(record.seq)

## `SeqRecord`

-   <https://biopython.org/wiki/SeqRecord>
-   <https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec32>
-   <https://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html>

`Seq` とそのほかの情報をひとまとまりにしたクラス。

`id`\
`seq`

`name`\
`description`\
`dbxrefs`\
`features`\
`annotations`\
`letter_annotations`

    >>> for record in SeqIO.parse(fin, 'fasta'):
    ...     print(record)
    ...
    ID: gi|186972394|gb|EU490707.1|
    Name: gi|186972394|gb|EU490707.1|
    Description: gi|186972394|gb|EU490707.1| Selenipedium aequinoctiale maturase K (matK) gene, partial cds; chloroplast
    Number of features: 0
    Seq('ATTTTTTACGAACCTGTGGAAATTTTTGGTTATGACAATAAATCTAGTTTAGTA...GAA', SingleLetterAlphabet())
    ID: gi|186972391|gb|ACC99454.1|
    Name: gi|186972391|gb|ACC99454.1|
    Description: gi|186972391|gb|ACC99454.1| maturase K [Scaphosepalum rapax]
    Number of features: 0
    Seq('IFYEPVEILGYDNKSSLVLVKRLITRMYQQKSLISSLNDSNQNEFWGHKNSFSS...EEE', SingleLetterAlphabet())

## `Seq`

-   <https://biopython.org/wiki/Seq>
-   <https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec17>
-   <https://biopython.org/DIST/docs/api/Bio.Seq.Seq-class.html>

塩基配列・アミノ酸配列のクラス。
標準 `str` とほとんど同じように扱えるほか、
相補鎖とかなんとかを簡単に扱えるメソッドが備わっている。

-   `complement(self)`
-   `reverse_complement(self)`
-   `transcribe(self)`
-   `back_transcribe(self)`
-   `translate(self, table="Standard", stop_symbol="*", to_stop=False, cds=False)`
-   `ungap(self, gap=None)`

## GenomeDiagram

<https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec328>
<https://biopython.org/DIST/docs/GenomeDiagram/userguide.pdf>

## Entrez

-   <https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc90>

### esearch

-   [NCBI ESearch Utility](http://www.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html)
-   <https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93>

    ```py
    >>> from Bio import Entrez
    >>> handle = Entrez.esearch(db="nucleotide", retmax=10, term="Opuntia")
    >>> record = Entrez.read(handle)
    >>> for item in record.items():
    ...     print(item)
    ...
    (u'Count', '390')
    (u'RetMax', '10')
    (u'IdList', ['257359511', '283467266', '246905625', '246905624', '246655205', '246655204', '240253899', '240253897', '240253576', '240253574'])
    (u'TranslationStack', [{u'Count': '200', u'Field': 'Organism', u'Term': '"Opuntia"[Organism]', u'Explode': 'Y'}, {u'Count': '390', u'Field': 'All Fields', u'Term': 'Opuntia[All Fields]', u'Explode': 'Y'}, 'OR', 'GROUP'])
    (u'TranslationSet', [{u'To': '"Opuntia"[Organism] OR Opuntia[All Fields]', u'From': 'Opuntia'}])
    (u'RetStart', '0')
    (u'QueryTranslation', '"Opuntia"[Organism] OR Opuntia[All Fields]')
    ```

### efetch

-   [NCBI EFetch Utility](http://eutils.ncbi.nlm.nih.gov/corehtml/query/static/efetch_help.html)
-   <https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc96>

返り値は結果(XML)へのハンドル。:

    from Bio import Entrez

    handle = Entrez.efetch(db="nucleotide", id="186972394,186972394", rettype="fasta")
    record = SeqIO.parse(handle, "fasta")
    for x in record:
        print(i)

## Installation

[pip]({{< relref "pip.md" >}}) で一発:

    pip install biopython


## 書籍

<a href="https://www.amazon.co.jp/dp/487311845X/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=72a416f5d10a9e84aaab4b3ee9613329&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311845X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=487311845X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873118417/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=6b1a04ec880b6c730bd6e80273e30e9c&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873118417&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873118417" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873117488/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=2181a50362009e68f507d44fc38716b4&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
