# 2장: Biopython으로 무엇을 할 수 있나요?

이 섹션은 Biopython을 빠르게 시작할 수 있도록 설계되었으며, 사용 가능한 기능과 사용 방법에 대한 전반적인 개요를 제공합니다. 이 섹션의 모든 예제는 Python에 대한 기본적인 지식이 있고, 시스템에 Biopython을 성공적으로 설치했다고 가정합니다. Python 지식을 복습해야 한다고 생각되면, Python 공식 웹사이트에서 시작하기에 좋은 무료 문서를 제공합니다 (https://docs.python.org/3/).

컴퓨터를 이용한 생물학 작업의 상당 부분이 인터넷 데이터베이스와의 연결을 포함하기 때문에, 일부 예제를 실행하려면 인터넷 연결이 필요합니다.

이제 준비가 되었으니, Biopython으로 할 수 있는 일들을 살펴보겠습니다.

## 2.1 Biopython이 제공하는 기능의 일반적인 개요

서론에서 언급했듯이, Biopython은 컴퓨터로 작업하는 생물학자들에게 중요한 "것들"을 다룰 수 있는 기능을 제공하는 라이브러리 모음입니다. 일반적으로 이는 최소한 약간의 프로그래밍 경험(물론 Python으로!)이 필요하거나, 적어도 프로그래밍을 배우려는 관심이 있어야 한다는 것을 의미합니다. Biopython의 역할은 재사용 가능한 라이브러리를 제공하여 프로그래머로서 여러분의 작업을 더 쉽게 만드는 것입니다. 이를 통해 특정 파일 형식을 파싱하는 내부 작업에 집중하는 대신 여러분의 관심사인 특정 질문에 답하는 데 집중할 수 있습니다 (물론, 존재하지 않는 파서를 작성하여 Biopython에 기여하고 싶다면 언제든 환영합니다!). 즉, Biopython의 목표는 여러분을 행복하게 만드는 것입니다!

Biopython에 대해 주목할 한 가지는 종종 "같은 일을 하는" 여러 가지 방법을 제공한다는 것입니다. 최근 릴리스에서 개선되었지만, 이는 여전히 Python에서 이상적으로 한 가지 올바른 방법이 있어야 한다는 점에서 좌절스러울 수 있습니다. 그러나 이는 라이브러리에 대한 많은 유연성과 제어를 제공하기 때문에 실제로는 이점이 될 수 있습니다. 이 튜토리얼은 일을 쉽게 처리할 수 있도록 일반적이거나 쉬운 방법을 보여줍니다. 대안적 가능성에 대해 더 알아보려면 Cookbook (22장, 여기에는 몇 가지 멋진 트릭과 팁이 있습니다), 내장 "docstrings" (Python help 명령어를 통해, 또는 API 문서), 또는 궁극적으로 코드 자체를 참조하세요.

## 2.2 시퀀스 다루기

논란의 여지가 있겠지만(물론!), 생물정보학의 중심 객체는 시퀀스입니다. 따라서 우리는 Biopython에서 시퀀스를 다루는 메커니즘인 Seq 객체에 대한 간단한 소개로 시작하겠습니다. 이에 대해서는 3장에서 더 자세히 다룰 것입니다.

대부분의 경우 시퀀스를 생각할 때 우리는 'AGTACACTGGT'와 같은 문자열을 떠올립니다. 다음과 같이 이 시퀀스로 Seq 객체를 만들 수 있습니다 - ">>>"는 Python 프롬프트를 나타내며 그 뒤에 입력할 내용이 옵니다:

```python
>>> from Bio.Seq import Seq
>>> my_seq = Seq("AGTACACTGGT")
>>> my_seq
Seq('AGTACACTGGT')
>>> print(my_seq)
AGTACACTGGT
```

Seq 객체는 지원하는 메서드에서 Python 문자열과 다릅니다. 일반 문자열로는 다음과 같은 작업을 할 수 없습니다:

```python
>>> my_seq
Seq('AGTACACTGGT')
>>> my_seq.complement()
Seq('TCATGTGACCA')
>>> my_seq.reverse_complement()
Seq('ACCAGTGTACT')
```

그 다음으로 중요한 클래스는 SeqRecord 또는 Sequence Record입니다. 이는 시퀀스(Seq 객체로)와 함께 식별자, 이름, 설명 등의 추가 주석을 포함합니다. 시퀀스 파일 형식을 읽고 쓰는 Bio.SeqIO 모듈은 SeqRecord 객체와 함께 작동하며, 이에 대해서는 아래에서 소개하고 5장에서 더 자세히 다룰 것입니다.

이것으로 Biopython 시퀀스 클래스의 기본 기능과 사용법을 다루었습니다. 이제 Biopython 라이브러리와 상호 작용하는 방법에 대한 아이디어를 얻었으니, 생물학적 파일 형식을 다루는 흥미진진한 세계로 뛰어들 시간입니다!

## 2.3 사용 예제

Biopython의 파서와 다른 모든 것들로 바로 뛰어들기 전에, 우리가 하는 모든 것에 동기를 부여하고 삶을 더 흥미롭게 만들기 위한 예제를 설정해 봅시다. 결국, 이 튜토리얼에 생물학이 없다면 왜 읽고 싶어 하겠습니까?

제가 식물을 좋아하기 때문에, 우리는 식물 기반 예제를 사용해야 할 것 같습니다 (다른 생물을 좋아하는 팬 여러분, 죄송합니다!). 최근 지역 온실을 방문한 후, 우리는 갑자기 레이디 슬리퍼 난초에 대한 놀라운 집착을 갖게 되었습니다 (왜인지 궁금하다면, Flickr에서 레이디 슬리퍼 난초 사진을 보거나 Google 이미지 검색을 해보세요).

물론 난초는 보기에 아름다울 뿐만 아니라 진화와 계통학을 연구하는 사람들에게도 매우 흥미롭습니다. 그래서 레이디 슬리퍼 진화에 대한 분자 연구를 위한 연구 제안서를 작성하려고 한다고 가정해 봅시다. 우리는 어떤 연구가 이미 수행되었는지, 그리고 우리가 어떻게 그것에 기여할 수 있는지 알고 싶습니다.

조금 조사해보니 레이디 슬리퍼 난초는 난초과(Orchidaceae)와 컴브레툼아과(Cypripedioideae)에 속하며 5개의 속으로 구성되어 있다는 것을 알게 되었습니다: Cypripedium, Paphiopedilum, Phragmipedium, Selenipedium, Mexipedium.

이제 더 많은 정보를 찾아볼 준비가 되었습니다. Biopython 도구가 어떻게 도움이 될 수 있는지 살펴보겠습니다. 2.4절에서 시퀀스 파싱부터 시작하겠지만, 난초는 나중에도 다시 등장할 것입니다 - 예를 들어 12장에서는 난초에 대한 논문을 PubMed에서 검색하고 GenBank에서 시퀀스 데이터를 추출할 것이며, 13장에서는 특정 난초 단백질에 대한 데이터를 Swiss-Prot에서 추출하고, 6.7.2절에서는 ClustalW를 사용한 난초 단백질의 다중 시퀀스 정렬 작업을 할 것입니다.

## 2.4 시퀀스 파일 형식 파싱하기

생물정보학 작업의 많은 부분은 생물학적 데이터를 저장하기 위해 설계된 다양한 파일 형식을 다루는 것과 관련이 있습니다. 이러한 파일들은 흥미로운 생물학적 데이터로 가득 차 있으며, 이를 어떤 프로그래밍 언어로 조작할 수 있는 형식으로 파싱하는 것이 특별한 과제입니다. 그러나 이러한 파일을 파싱하는 작업은 형식이 꽤 자주 변경될 수 있다는 사실과, 형식에 작은 미묘한 차이가 포함될 수 있어 가장 잘 설계된 파서조차도 깨질 수 있다는 점 때문에 좌절될 수 있습니다.

이제 Bio.SeqIO 모듈을 간단히 소개하겠습니다 - 자세한 내용은 5장에서 확인할 수 있습니다.

우리의 친구인 레이디 슬리퍼 난초에 대한 온라인 검색부터 시작하겠습니다. 이 소개를 간단하게 하기 위해, 우리는 NCBI 웹사이트를 수동으로 사용하고 있습니다. NCBI의 핵산 데이터베이스를 살펴보겠습니다. Entrez 온라인 검색(https://www.ncbi.nlm.nih.gov/nuccore/?term=Cypripedioideae)을 사용하여 Cypripedioideae(레이디 슬리퍼 난초의 아과)라는 텍스트가 언급된 모든 것을 검색합니다.

이 튜토리얼이 처음 작성되었을 때, 이 검색은 단 94개의 결과만을 제공했으며, 우리는 이를 FASTA 형식 텍스트 파일과 GenBank 형식 텍스트 파일로 저장했습니다(ls_orchid.fasta와 ls_orchid.gbk 파일, Biopython 소스 코드의 Doc/examples/ 아래에도 포함되어 있습니다).

오늘 검색을 실행하면 수백 개의 결과를 얻을 것입니다! 튜토리얼을 따라가면서 동일한 유전자 목록을 보고 싶다면, 위의 두 파일을 다운로드하거나 Biopython 소스 코드의 docs/examples/에서 복사하세요. 2.5절에서는 Python 내에서 이와 같은 검색을 수행하는 방법을 살펴볼 것입니다.

### 2.4.1 간단한 FASTA 파싱 예제

좋아하는 텍스트 편집기에서 레이디 슬리퍼 난초 FASTA 파일 ls_orchid.fasta를 열면, 파일이 다음과 같이 시작되는 것을 볼 수 있습니다:

```
>gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAACGATCGAGTG
AATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGTGACCCTGATTTGTTGTTGGG
...
```

이 파일은 94개의 레코드를 포함하고 있으며, 각 레코드는 ">" (큰따옴표) 기호로 시작하는 줄 다음에 한 줄 이상의 시퀀스가 이어집니다. 이제 Python에서 다음을 시도해 보세요:

```python
>>> from Bio import SeqIO
>>> for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
...     print(seq_record.id)
...     print(repr(seq_record.seq))
...     print(len(seq_record))
...
```

화면에 다음과 같은 내용이 출력될 것입니다:

```
gi|2765658|emb|Z78533.1|CIZ78533
Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
740
...
gi|2765564|emb|Z78439.1|PBZ78439
Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
592
```

### 2.4.2 간단한 GenBank 파싱 예제

이번에는 GenBank 파일 ls_orchid.gbk를 로드해 보겠습니다. FASTA 파일에 사용한 코드 조각과 거의 동일하다는 점에 주목하세요. 유일한 차이점은 파일 이름과 형식 문자열을 변경한 것뿐입니다:

```python
>>> from Bio import SeqIO
>>> for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
...     print(seq_record.id)
...     print(repr(seq_record.seq))
...     print(len(seq_record))
...
```

이는 다음과 같은 결과를 출력할 것입니다:

```
Z78533.1
Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
740
...
Z78439.1
Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
592
```

이 경우 seq_record.id로 더 짧은 문자열이 사용되었음을 알 수 있습니다.

### 2.4.3 파싱이 너무 좋아요 - 더 이야기해 주세요!

Biopython에는 많은 파서가 있으며, 각각은 파싱하는 시퀀스 형식과 관련된 고유한 특성을 가지고 있습니다. 5장에서는 Bio.SeqIO에 대해 더 자세히 다루며, 6장에서는 시퀀스 정렬을 위한 Bio.Align를 소개합니다.

가장 인기 있는 파일 형식에 대한 파서는 Bio.SeqIO 및/또는 Bio.AlignIO에 통합되어 있지만, 일부 희귀하고 관심이 적은 파일 형식의 경우 파서가 전혀 없거나 아직 연결되지 않은 오래된 파서가 있을 수 있습니다. 최신 정보는 위키 페이지 http://biopython.org/wiki/SeqIO 와 http://biopython.org/wiki/AlignIO 를 확인하거나 메일링 리스트에 문의하세요. 위키 페이지에는 지원되는 파일 유형의 최신 목록과 일부 추가 예제가 포함되어 있어야 합니다.

특정 파서와 이를 사용하여 멋진 작업을 수행하는 방법에 대한 정보를 찾을 수 있는 다음 장소는 이 튜토리얼의 22장인 Cookbook입니다. 찾고 있는 정보를 발견하지 못했다면, 열심히 일하는 문서 작성자들을 도와 그것에 대한 cookbook 항목을 제출하는 것을 고려해 보세요! (물론 어떻게 하는지 알아낸 후에 말이죠!)

## 2.5 생물학적 데이터베이스에 연결하기

생물정보학에서 매우 일반적으로 해야 하는 일 중 하나는 생물학적 데이터베이스에서 정보를 추출하는 것입니다. 특히 반복적인 작업이 많은 경우 이러한 데이터베이스에 수동으로 접근하는 것은 상당히 지루할 수 있습니다. Biopython은 일부 온라인 데이터베이스를 Python 스크립트에서 사용할 수 있게 함으로써 시간과 노력을 절약하려고 합니다. 현재 Biopython은 다음 데이터베이스에서 정보를 추출하는 코드를 가지고 있습니다:

• NCBI의 Entrez (및 PubMed) - 12장 참조.
• ExPASy - 13장 참조.
• SCOP - Bio.SCOP.search() 함수 참조.

이러한 모듈의 코드는 기본적으로 이러한 페이지의 CGI 스크립트와 상호 작용하는 Python 코드를 쉽게 작성할 수 있게 해주어, 다루기 쉬운 형식으로 결과를 얻을 수 있게 합니다. 일부 경우에는 결과가 Biopython 파서와 긴밀하게 통합되어 정보를 추출하는 것이 더욱 쉬워집니다.

## 2.6 다음으로 할 일

여기까지 오셨다면 Biopython의 기본 사항을 잘 이해하고 있을 것이며 이제 유용한 작업을 시작할 준비가 되었을 것입니다. 지금 할 수 있는 가장 좋은 일은 이 튜토리얼을 끝까지 읽는 것이고, 그 다음에는 원한다면 소스 코드를 살펴보고 자동 생성된 문서를 확인해 보는 것입니다.

무엇을 하고 싶은지, 그리고 Biopython의 어떤 라이브러리가 그 일을 할 수 있는지에 대한 그림이 잡히면 Cookbook (22장)을 살펴보세요. 여러분이 하고 싶은 일과 유사한 작업을 수행하는 예제 코드가 있을 수 있습니다.

무엇을 하고 싶은지는 알지만 어떻게 해야 할지 모르겠다면, 주저하지 말고 Biopython 메인 리스트에 질문을 올려주세요 (http://biopython.org/wiki/Mailing_lists 참조). 이는 우리가 여러분의 질문에 답하는 데 도움이 될 뿐만 아니라, 다음 사람이 여러분이 하고 싶은 일을 할 수 있도록 문서를 개선할 수 있게 해줄 것입니다.

코드를 즐기세요!