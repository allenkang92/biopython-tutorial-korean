# 1장 소개

## 1.1 Biopython이란 무엇인가요?

Biopython은 계산 분자생물학을 위한 [파이썬](https://www.python.org) 모듈 모음집이라고 할 수 있습니다. 파이썬은 과학 컴퓨팅 분야에서 널리 사용되는 인터프리터 기반의 객체지향적이며, 유연한 언어입니다. 파이썬은 배우기 쉽고, 매우 명확한 구문을 가지고 있으며, C, C++ 또는 FORTRAN으로 작성된 모듈을 가지고 쉽게 확장할 수 있습니다. 2000년에 처음 개발된 이후로, Biopython은 전 세계의 많은 봉사자들에 의해 지속적으로 개발되며 유지되고 있습니다.

Biopython 웹사이트([http://www.biopython.org](http://www.biopython.org))에서는 모듈, 스크립트 및 웹 링크를 포함한 온라인 리소스를 제공하며, 이는 파이썬 기반의 생물정보학 소프트웨어를 개발하고 연구하는 개발자들을 위한 것입니다. Biopython은 다양한 생물정보학 파일 형식(BLAST, Clustalw, FASTA, Genbank, ...)의 파서를 제공하며, 온라인 서비스(NCBI, Expasy, ...)에 대한 접근 기능을 포함합니다. 표준 시퀀스 클래스, 시퀀스 정렬 및 모티프 분석 도구, 클러스터링 알고리즘, 구조 생물학을 위한 모듈, 계통학 분석을 위한 모듈도 포함되어 있습니다.

## 1.2 Biopython 패키지에서 찾을 수 있는 기능

Biopython의 주요 릴리스는 다음과 같은 다양한 기능들을 제공합니다:

- 다양한 생물정보학 파일을 파이썬에서 사용할 수 있는 데이터 구조로 파싱할 수 있으며, 다음과 같은 형식을 지원합니다.
  - BLAST 결과 - 독립형 및 웹 기반 BLAST 모두 지원
  - Clustalw
  - FASTA
  - GenBank
  - PubMed 및 Medline
  - ExPASy 파일 (Enzyme, Prosite 등)
  - SCOP ('dom' 및 'lin' 파일 포함)
  - UniGene
  - SwissProt

- 지원되는 파일 형식은 레코드 단위로 반복(iteration)하거나 Dictionary 인터페이스를 통해 인덱싱(indexing)하여 접근할 수 있습니다.
- 다음과 같은 인기 있는 온라인 생물정보학 서비스와 연동할 수 있습니다.
  - **NCBI** – BLAST, Entrez, PubMed 서비스
  - **ExPASy** – Swiss-Prot, Prosite 엔트리 및 Prosite 검색
- 다음과 같은 일반적인 생물정보학 프로그램과의 인터페이스를 제공합니다.
  - **NCBI의 독립형 BLAST**
  - **Clustalw alignment 프로그램**
  - **EMBOSS 커맨드 라인 도구**
- 시퀀스, 시퀀스 ID 및 시퀀스 특징을 처리할 수 있는 표준 시퀀스 클래스를 제공합니다.
- 번역, 전사, 무게 계산과 같은 일반적인 시퀀스 작업을 수행할 수 있는 도구를 제공합니다.
- k-최근접 이웃(k-NN), 나이브 베이즈(Naive Bayes) 또는 서포트 벡터 머신(SVM)을 사용하여 데이터를 분류할 수 있는 코드를 제공합니다.
- 정렬을 처리하는 코드와 치환 행렬(서브스티튜션 매트릭스)을 생성하고 처리하는 표준 방법을 제공합니다.
- 병렬화할 수 있는 작업을 별도의 프로세스로 분할하는 코드를 제공합니다.
- 기본 시퀀스 조작, 번역, BLAST 등을 수행할 수 있는 GUI 기반 프로그램이 포함되어 있습니다.
- 광범위한 문서와 모듈 사용에 관한 도움말을 제공하며, 이 파일, 온라인 위키 문서, 웹사이트 및 메일링 리스트를 통해 도움말을 찾을 수 있습니다.
- BioSQL과 통합되어 있으며, 이는 BioPerl 및 BioJava 프로젝트에서도 지원되는 시퀀스 데이터베이스 스키마입니다.

이러한 다양한 기능을 통해 Biopython을 다운로드하고 사용해보시길 바랍니다!

## 1.3 Biopython 설치하기

Biopython의 모든 설치 정보는 이 문서에서 분리되어 있어 최신 상태를 유지하기 쉽습니다.

간단한 방법은 다음과 같습니다: `pip install biopython` 명령어를 사용하세요. 그 외의 설치 방법은 메인 README 파일[https://github.com/biopython/biopython/blob/master/README.rst]을 참고하세요.

## 1.4 자주 묻는 질문 (FAQ)

1. **Biopython을 과학 출판물(scientific publication)에서 어떻게 인용해야 하나요?**

   Biopython의 주요 인용으로 우리의 응용 프로그램 노트 [5, Cock et al., 2009]를 사용해 주세요. 추가로, 특정 Biopython 모듈에 대한 참고문헌이 필요하다면 아래 목록의 출판물을 참조해 주세요 (자세한 내용은 웹사이트에서 확인할 수 있습니다).
   - 공식 프로젝트 발표: [4, Chapman and Chang, 2000]
   - **Bio.PDB**: [16, Hamelryck and Manderick, 2003]
   - **Bio.Cluster**: [10, De Hoon et al., 2004]
   - **Bio.Graphics.GenomeDiagram**: [35, Pritchard et al., 2006]
   - **Bio.Phylo** 및 **Bio.Phylo.PAML**: [45, Talevich et al., 2012]
   - Biopython, BioPerl, BioRuby, BioJava 및 EMBOSS에서 지원되는 **FASTQ** 파일 형식: [6, Cock et al., 2010]

2. **"Biopython"의 올바른 표기는 무엇인가요? "BioPython"이라고 적어도 되나요?**

   올바른 표기는 "Biopython"입니다. "BioPython"은 맞지 않으니 주의하세요.

3. **Biopython 소프트웨어는 어떤 라이선스로 배포되나요?**

   Biopython은 Biopython License Agreement 하에 배포됩니다. 그러나 Biopython 1.69 이후로는 일부 파일이 Biopython License Agreement 또는 BSD 3-Clause License 중에서 선택하여 사용할 수 있도록 이중 라이선스 체계로 명시되어 있습니다. 이는 향후 모든 Biopython을 이중 라이선스 체계로 제공하기 위한 목적입니다.

4. **Biopython 로고는 무엇이며, 어떤 라이선스로 배포되나요?**

   2017년 7월 Biopython 1.70 릴리스 기준으로 Biopython 로고는 노란색과 파란색 뱀이 이중 나선을 형성하고 그 아래에 소문자로 "biopython"이라고 적혀 있습니다. Patrick Kunzmann이 디자인한 이 로고는 Biopython License Agreement 또는 BSD 3-Clause License 중에서 선택할 수 있는 이중 라이선스 하에 배포됩니다.

   이전에는 두 마리의 노란색 뱀이 이중 나선을 형성하고 그 사이에 "BIOPYTHON"이라고 적혀 있는 로고가 사용되었습니다. 이 로고는 2003년 Henrik Vestergaard와 Thomas Hamelryck가 공개 경쟁의 일환으로 디자인한 것입니다.

5. **릴리스마다 무엇이 새로운지 알려주는 변경 로그가 있나요?**

   소스 코드와 함께 제공되는 `NEWS.rst` 파일을 확인하거나 깃허브의 최신 `NEWS` 파일을 읽어보세요.

6. **내 `print` 명령어에 문제가 생겼어요. 뭐가 잘못된 건가요?**

   Biopython 1.77부터는 Python 3만 지원하기 때문에, 이 튜토리얼도 Python 3 스타일의 `print` 함수를 사용합니다.

7. **내가 설치한 Biopython 버전이 무엇인지 어떻게 알 수 있죠?**

   아래 코드를 사용하세요:
   ```python
   import Bio
   print(Bio.__version__)
   ```

   만약 `import Bio`에서 오류가 발생한다면 Biopython이 설치되어 있지 않은 것입니다. `__version__` 앞뒤로 언더스코어가 두 번씩 들어가야 한다는 점을 유의하세요. 두 번째 라인이 실패한다면, 현재 설치된 버전이 매우 오래된 것입니다.

   버전 문자열 끝에 + 기호가 붙은 경우, 이는 공식 릴리스가 아니라 그 버전이 배포된 후의 개발 코드 스냅샷을 의미합니다. 이 표기법은 Biopython 1.68이 출시된 2016년 6월까지 사용되었습니다.

   버전 문자열 끝에 .dev<number>와 같은 형태로 나타나는 경우, 이는 공식 릴리스가 아니라 그 버전이 배포되기 전의 개발 코드 스냅샷을 의미합니다.

8. **이 문서의 최신 버전은 어디에서 찾을 수 있나요?**

   Biopython 소스 코드 아카이브를 다운로드하면 HTML 및 PDF 형식의 해당 버전을 확인할 수 있습니다. 최신 공식 버전은 다음 링크에서 찾을 수 있습니다:

   - [튜토리얼 (HTML 버전)](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
   - [튜토리얼 (PDF 버전)](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)

9. **서열 비교에 문제가 있어요. 무엇이 잘못된 건가요?**

   Biopython 1.65부터 Seq와 MutableSeq 클래스(및 그 서브클래스)는 문자열 기반 비교를 사용합니다. 이는 `str(seq1) == str(seq2)`와 같이 명시적으로 수행할 수 있습니다.

   이전 버전의 Biopython에서는 Seq 객체에 대한 인스턴스 기반 비교를 사용했으며, 이는 `id(seq1) == id(seq2)`와 같이 명시적으로 수행할 수 있습니다.

   만약 이전 버전의 Biopython을 지원해야 한다면 이러한 명시적 비교를 사용하여 문제를 피할 수 있습니다. 자세한 내용은 섹션 3.10을 참조하세요.

10. **Bio.SeqIO와 Bio.AlignIO에서 어떤 파일 형식을 읽고 쓸 수 있나요?**

    내장된 docstring을 확인하세요. 예를 들어 `from Bio import SeqIO`를 사용한 후 `help(SeqIO)` 명령을 입력해보세요. 또는 [SeqIO 위키 페이지](http://biopython.org/wiki/SeqIO)와 [AlignIO 위키 페이지](http://biopython.org/wiki/AlignIO)에서 최신 목록을 확인할 수 있습니다.

11. **왜 Bio.SeqIO와 Bio.AlignIO의 parse, read 및 write 함수는 파일 이름 대신 핸들을 요구하나요?**

    Biopython 1.54 이후 버전에서는 이러한 함수가 파일 이름이 아닌 핸들을 사용합니다. 이전 버전을 사용 중이시라면 핸들을 명시적으로 사용하세요 (섹션 25.1 참조). 특히 데이터 쓰기가 완료된 후 출력 핸들을 명시적으로 닫아야 합니다.

12. **왜 Bio.SeqIO.write()와 Bio.AlignIO.write() 함수는 단일 레코드나 정렬을 받아들이지 않고 리스트나 반복자(iterator)를 요구하나요?**

    Biopython 1.54 이후 버전에서는 이러한 함수가 단일 레코드 대신 리스트나 반복자를 요구합니다. 단일 레코드를 전달하려면 `[...]`를 사용하여 하나의 요소를 가진 리스트를 생성하세요.

13. **str(...) 함수가 Seq 객체의 전체 서열을 반환하지 않는 이유는 무엇인가요?**

    Biopython 1.45 이상 버전을 사용하세요.

14. **Bio.Blast가 최신의 일반 텍스트 NCBI BLAST 결과와 함께 작동하지 않는 이유는 무엇인가요?**

    NCBI는 BLAST 도구의 일반 텍스트 출력을 지속적으로 조정합니다. 최신 버전의 Biopython을 사용하지 않는다면 업그레이드를 시도해보세요. 하지만 NCBI와 Biopython 양측에서는 XML 출력을 사용할 것을 권장합니다. XML은 컴퓨터 프로그램이 읽도록 설계되어 있습니다.

15. **왜 내 Bio.Entrez.efetch() 스크립트가 더 이상 작동하지 않나요?**

    이는 2012년 2월에 도입된 EFetch 2.0에 따른 NCBI 변경 사항 때문일 수 있습니다. 우선, 기본 반환 모드가 변경되었으므로 `retmode="text"`를 호출에 추가해야 합니다. 둘째, ID 목록을 제공하는 방식에 대해 더 엄격해졌습니다. Biopython 1.59 이후 버전은 리스트를 자동으로 쉼표로 구분된 문자열로 변환합니다.

16. **Bio.Blast.NCBIWWW.qblast()가 NCBI BLAST 웹사이트와 동일한 결과를 제공하지 않는 이유는 무엇인가요?**

    같은 옵션을 사용해야 합니다. NCBI는 웹사이트의 기본 설정을 자주 변경하며, 이 설정이 qblast의 기본값과 일치하지 않을 수 있습니다. 갭(gap) 패널티와 기대값(expectation threshold)과 같은 옵션을 확인하세요.

17. **왜 SeqRecord 객체를 서로 더할 수 없나요?**

    Biopython 1.53 이후 버전을 사용하세요.

18. **Bio.SeqIO.index_db()가 작동하지 않는 이유는 무엇인가요? 모듈은 제대로 임포트되는데 index_db 함수가 없다고 나와요!**

    Biopython 1.57 이상 버전과 SQLite3 지원이 있는 파이썬 버전을 사용하세요.

19. **MultipleSeqAlignment 객체는 어디에 있나요? Bio.Align 모듈은 잘 임포트되지만 이 클래스는 찾을 수 없어요!**

    Biopython 1.54 이상 버전을 사용하세요. 또는 `Bio.Align.Generic.Alignment` 클래스가 일부 기능을 지원하지만, 이 클래스의 사용은 이제 권장되지 않습니다.

20. **왜 애플리케이션 래퍼를 통해 커맨드라인 도구를 바로 실행할 수 없나요?**

    Biopython 1.55 이후 버전에서 지원되지만, 이 기능은 Biopython 1.78부터 사용 중단(deprecated)되었습니다. 대신 파이썬의 `subprocess` 모듈을 직접 사용하세요.

21. **디렉터리에서 코드를 찾으려고 했는데, 해당 작업을 수행하는 코드를 찾을 수 없어요. 코드가 어디에 숨겨져 있나요?**

    한 가지 알아둘 점은 우리는 `__init__.py` 파일에 코드를 넣는 경우가 있습니다. 이 파일에서 코드를 찾는 것에 익숙하지 않은 경우 혼란스러울 수 있습니다. 이렇게 하는 이유는 사용자들이 더 쉽게 임포트할 수 있도록 하기 위함입니다. 예를 들어, `from Bio.GenBank import GenBank`와 같이 "중복된" 임포트 대신 `from Bio import GenBank`와 같이 사용할 수 있습니다.

22. **왜 Bio.Fasta가 작동하지 않나요?**

    Bio.Fasta 모듈은 Biopython 1.51 (2009년 8월)에서 사용 중단(deprecated)되었으며, Biopython 1.55 (2010년 8월)에서 제거되었습니다. 대신 Bio.SeqIO를 사용하도록 기존 코드를 변환하는 간단한 예제가 DEPRECATED.rst 파일에 포함되어 있습니다.

더 일반적인 질문에 대한 답변은 [파이썬 FAQ 페이지](https://docs.python.org/3/faq/index.html)를 참고할 수 있습니다.