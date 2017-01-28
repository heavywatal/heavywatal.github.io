+++
title = 'BioPython'
subtitle = "Tools for biological computation"
tags = ["python"]
[menu.main]
  parent = "python"
+++

<http://biopython.org>
<http://biopython.org/DIST/docs/api/Bio-module.html>

## `SeqIO`

-   <http://biopython.org/wiki/SeqIO>
-   <http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec51>
-   <http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html>

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

-   <http://biopython.org/wiki/SeqRecord>
-   <http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec32>
-   <http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html>

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

-   <http://biopython.org/wiki/Seq>
-   <http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec17>
-   <http://biopython.org/DIST/docs/api/Bio.Seq.Seq-class.html>

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

<http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec328>
<http://biopython.org/DIST/docs/GenomeDiagram/userguide.pdf>

## Entrez

-   <http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc90>

### esearch

-   [NCBI ESerch Utility](http://www.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html)
-   <http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93>

<!-- -->

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

### efetch

-   [NCBI EFetch Utility](http://eutils.ncbi.nlm.nih.gov/corehtml/query/static/efetch_help.html)
-   <http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc96>

返り値は結果(XML)へのハンドル。:

    from Bio import Entrez

    handle = Entrez.efetch(db="nucleotide", id="186972394,186972394", rettype="fasta")
    record = SeqIO.parse(handle, "fasta")
    for x in record:
        print(i)

## Installation

[pip]({{< relref "pip.md" >}}) で一発:

    % pip install biopython


## 書籍

<a href="https://www.amazon.co.jp/dp/479738946X/ref=as_li_ss_il?ie=UTF8&qid=1485612008&sr=8-6&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=5ea5e48ecc83b9439f21406b6f57c062" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=479738946X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=479738946X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/IPython%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9%E3%82%AF%E3%83%83%E3%82%AF%E3%83%96%E3%83%83%E3%82%AF-%E5%AF%BE%E8%A9%B1%E5%9E%8B%E3%82%B3%E3%83%B3%E3%83%94%E3%83%A5%E3%83%BC%E3%83%86%E3%82%A3%E3%83%B3%E3%82%B0%E3%81%A8%E5%8F%AF%E8%A6%96%E5%8C%96%E3%81%AE%E3%81%9F%E3%82%81%E3%81%AE%E3%83%AC%E3%82%B7%E3%83%94%E9%9B%86-Cyrille-Rossant/dp/4873117488/ref=as_li_ss_il?_encoding=UTF8&psc=1&refRID=X16VFSS3W75RMTG7VGCH&linkCode=li3&tag=heavywatal-22&linkId=b79e2290571289b02621392257a4ac1c" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/Python%E3%81%AB%E3%82%88%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E5%88%86%E6%9E%90%E5%85%A5%E9%96%80-_NumPy%E3%80%81pandas%E3%82%92%E4%BD%BF%E3%81%A3%E3%81%9F%E3%83%87%E3%83%BC%E3%82%BF%E5%87%A6%E7%90%86-Wes-McKinney/dp/4873116554/ref=as_li_ss_il?s=books&ie=UTF8&qid=1485612076&sr=1-16&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=f024f5e24f16c38402149f97591c8aab" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873116554&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873116554" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/%E3%82%BC%E3%83%AD%E3%81%8B%E3%82%89%E4%BD%9C%E3%82%8BDeep-Learning-Python%E3%81%A7%E5%AD%A6%E3%81%B6%E3%83%87%E3%82%A3%E3%83%BC%E3%83%97%E3%83%A9%E3%83%BC%E3%83%8B%E3%83%B3%E3%82%B0%E3%81%AE%E7%90%86%E8%AB%96%E3%81%A8%E5%AE%9F%E8%A3%85-%E6%96%8E%E8%97%A4-%E5%BA%B7%E6%AF%85/dp/4873117585/ref=as_li_ss_il?s=books&ie=UTF8&qid=1485612076&sr=1-3&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=cde2e18b945af3ead43dc15e51b00af6" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117585&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117585" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
