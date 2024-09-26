# 4장: 시퀀스 주석 객체

## 4.1 SeqRecord 객체

`SeqRecord`(Sequence Record) 클래스는 `Bio.SeqRecord` 모듈에서 정의됩니다. 이 클래스는 식별자 및 기능과 같은 상위 수준의 기능을 시퀀스에 연결할 수 있도록 하며, `Bio.SeqIO`의 시퀀스 입출력 인터페이스(5장에서 자세히 다루어집니다)를 위한 기본 데이터 유형입니다.

SeqRecord 클래스 자체는 비교적 간단하며, 다음과 같은 속성들을 제공합니다:

- `.seq`: 시퀀스 자체로, 일반적으로 `Seq` 객체로 표현됩니다.
- `.id`: 시퀀스를 식별하는 기본 ID로, 문자열 형태를 가집니다. 대부분의 경우, 이는 액세션 번호와 같은 것일 수 있습니다.
- `.name`: 시퀀스의 일반적인 이름/ID로, 문자열 형태를 가집니다. 어떤 경우에는 액세션 번호와 동일할 수도 있지만, 클론 이름과 같을 수도 있습니다. GenBank 레코드의 LOCUS ID와 유사하다고 볼 수 있습니다.
- `.description`: 시퀀스에 대한 사람이 읽을 수 있는 설명 또는 표현적 이름으로, 문자열 형태를 가집니다.
- `.letter_annotations`: 시퀀스의 각 문자에 대한 추가 정보를 담은 제한된 딕셔너리 형태로, 키는 정보의 이름을 나타내고, 값은 시퀀스와 동일한 길이를 가지는 파이썬 시퀀스(예: 리스트, 튜플 또는 문자열)입니다. 이는 종종 품질 점수나 이차 구조 정보를 위해 사용됩니다.
- `.annotations`: 시퀀스에 대한 추가 정보를 담은 딕셔너리 형태로, 키는 정보의 이름을 나타내고, 값은 해당 정보의 내용입니다. 이를 통해 시퀀스에 더 많은 "구조화되지 않은" 정보를 추가할 수 있습니다.
- `.features`: 시퀀스의 특징에 대한 더 구조화된 정보를 가진 `SeqFeature` 객체의 리스트입니다 (예: 유전체 상의 유전자 위치, 또는 단백질 시퀀스의 도메인). 시퀀스 기능의 구조는 4.3절에서 설명합니다.
- `.dbxrefs`: 문자열 형태의 데이터베이스 상호 참조의 리스트입니다.

---

## 4.2 SeqRecord 객체 생성 및 접근

`SeqRecord` 객체는 `Bio.SeqRecord` 모듈에서 임포트할 수 있으며, 이를 사용하여 시퀀스 정보를 보강할 수 있습니다.

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

simple_seq = Seq("GATC")
simple_seq_r = SeqRecord(simple_seq)
print(simple_seq_r)
위의 코드는 단순한 SeqRecord 객체를 생성하는 예시입니다. 이렇게 생성된 객체는 Seq 객체 외에도 추가적인 주석 정보를 담을 수 있습니다.

python
코드 복사
simple_seq_r.id = "ABC123"
simple_seq_r.description = "간단한 예시 시퀀스"
print(simple_seq_r)
출력 결과는 다음과 같습니다:

perl
코드 복사
ID: ABC123
Name: <unknown name>
Description: 간단한 예시 시퀀스
Number of features: 0
Seq('GATC')
이렇게 SeqRecord 객체에 추가적인 정보를 할당할 수 있습니다. id, name, description 등은 SeqRecord 객체에 흔히 사용되는 속성들입니다.

4.3 시퀀스 기능 (Features) 및 SeqFeature 객체
생물정보학에서 시퀀스 기능(Feature)이란 특정 위치에 존재하는 유전자, CDS, rRNA, 기타 기능적 요소 등을 의미합니다. Biopython에서는 이러한 기능을 SeqFeature 객체로 표현하며, 이는 Bio.SeqFeature 모듈에서 임포트할 수 있습니다.

python
코드 복사
from Bio.SeqFeature import SeqFeature, FeatureLocation
SeqFeature 객체는 location, type, strand와 같은 주요 속성을 가집니다.

location: 해당 기능이 시퀀스에서 차지하는 위치를 나타내는 FeatureLocation 객체입니다.
type: 해당 기능의 종류를 나타내는 문자열로, 예를 들어 'gene', 'CDS', 'repeat_region' 등이 될 수 있습니다.
strand: 해당 기능이 양성 가닥(+1), 음성 가닥(-1), 또는 지정되지 않음(None) 중 어느 것인지 나타내는 정보입니다.
다음은 SeqFeature 객체를 생성하고 SeqRecord 객체에 추가하는 예시입니다:

python
코드 복사
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

simple_seq = Seq("GATCGATCGATCGATC")
simple_seq_r = SeqRecord(simple_seq, id="XYZ789", description="테스트 시퀀스")

# 기능의 위치 지정 (3에서 12까지)
feature_loc = FeatureLocation(start=3, end=12)
simple_feature = SeqFeature(location=feature_loc, type="CDS")

# 기능을 시퀀스 레코드에 추가
simple_seq_r.features.append(simple_feature)
print(simple_seq_r)
출력 결과는 다음과 같습니다:

perl
코드 복사
ID: XYZ789
Name: <unknown name>
Description: 테스트 시퀀스
Number of features: 1
Seq('GATCGATCGATCGATC')
4.4 기능 위치 지정
FeatureLocation 객체는 생물학적 기능이 시퀀스 내에서 어디에 위치하는지를 나타내는 데 사용됩니다. 이는 시퀀스의 서브시퀀스(부분 시퀀스)를 나타낼 때 유용하게 사용됩니다.

위치 정보를 생성할 때 start와 end 값을 이용합니다. 이때 인덱스는 파이썬의 일반적인 규칙을 따르며 0부터 시작하고, 끝 인덱스는 포함되지 않습니다. 예를 들어, FeatureLocation(5, 10)은 시퀀스의 5번째부터 9번째까지의 서브시퀀스를 나타냅니다.

python
코드 복사
from Bio.SeqFeature import FeatureLocation

feature_loc = FeatureLocation(start=5, end=10)
print(feature_loc)
4.5 기능 방향성 (Strand)
SeqFeature의 strand 속성은 해당 기능이 양성 가닥(+1), 음성 가닥(-1), 또는 지정되지 않음(None)인지 나타냅니다.

+1: 양성 가닥 (5' -> 3')
-1: 음성 가닥 (3' -> 5')
None: 방향성 미지정 (해당 기능의 방향성이 중요하지 않은 경우)
4.6 SeqFeature의 중첩
Biopython에서는 한 SeqFeature 안에 다른 SeqFeature가 포함될 수 있습니다. 이는 예를 들어 유전자(gene) 기능 안에 코딩 서열(CDS)과 같은 기능이 포함되는 경우 유용합니다. 이를 위해 SubFeature 리스트를 SeqFeature에 추가할 수 있습니다.

python
코드 복사
from Bio.SeqFeature import SeqFeature, FeatureLocation

# 유전자 기능 생성
gene_feature = SeqFeature(FeatureLocation(0, 20), type="gene")

# CDS 기능 생성 및 gene 기능에 포함시킴
cds_feature = SeqFeature(FeatureLocation(5, 15), type="CDS")
gene_feature.sub_features = [cds_feature]

print(gene_feature)
4.7 다수의 위치를 가지는 기능 (Compound Location)
때때로 시퀀스 기능은 여러 위치에 걸쳐 나타날 수 있습니다. 예를 들어 분절화된 엑손(exon)은 시퀀스의 비연속적인 부분에 위치할 수 있습니다. Biopython에서는 이러한 비연속적 기능을 CompoundLocation 객체로 표현합니다.

python
코드 복사
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation

# 두 개의 위치로 구성된 기능 생성
exon1 = FeatureLocation(0, 50)
exon2 = FeatureLocation(100, 150)
compound_loc = CompoundLocation([exon1, exon2])

# 해당 위치에 대한 기능 생성
multi_exon_feature = SeqFeature(compound_loc, type="exon")

print(multi_exon_feature)
4.8 기능 및 서브시퀀스 추출
SeqFeature는 .extract() 메서드를 사용하여 기능의 위치에 해당하는 서브시퀀스를 추출할 수 있습니다.

python
코드 복사
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

# 시퀀스 생성
my_seq = Seq("AGTACACTGGTACGCGTTACGTAGCT")

# 특정 위치에 대한 기능 생성
feature_loc = FeatureLocation(5, 15)
my_feature = SeqFeature(location=feature_loc, type="domain")

# 기능의 서브시퀀스 추출
sub_seq = my_feature.extract(my_seq)
print(sub_seq)  # 출력: 'ACTGGTACGC'
4.9 SeqRecord의 데이터베이스 교차참조 (dbxrefs)
SeqRecord 객체에는 해당 시퀀스와 관련된 다른 데이터베이스의 레코드에 대한 참조를 저장하는 dbxrefs 속성이 있습니다. 이 속성은 문자열로 구성된 리스트 형태로 저장됩니다.

python
코드 복사
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 데이터베이스 참조 추가
my_seq = Seq("ATGCGTACGTAGC")
my_record = SeqRecord(my_seq, id="ABC123", description="예시 시퀀스")
my_record.dbxrefs = ["GeneID:123456", "UniProtKB/Swiss-Prot:P12345"]

print(my_record)
4.10 특징별 레터 주석 (Letter Annotation)
SeqRecord 객체는 letter_annotations 속성을 통해 시퀀스의 각 염기나 아미노산에 대한 추가 정보를 저장할 수 있습니다. 이 속성은 사전(dictionary) 형태로 저장되며, 키는 정보의 이름이고, 값은 시퀀스와 길이가 동일한 파이썬 시퀀스(예: 리스트, 튜플)입니다.

python
코드 복사
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 시퀀스 생성
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
my_record = SeqRecord(my_seq)

# 각 염기에 대한 품질 점수 추가
my_record.letter_annotations["phred_quality"] = [40, 30, 20, 35, 25, 20, 30, 35, 20, 25,
                                                 35, 20, 30, 20, 25, 20, 30, 35, 30, 40,
                                                 25, 20, 35, 20, 30, 20, 25, 40, 20, 35,
                                                 30, 20, 35, 25, 30]
print(my_record.letter_annotations)
4.11 SeqRecord의 출력
SeqRecord 객체의 내용을 출력할 때는 객체의 주요 정보만 요약하여 표시됩니다. Bio.SeqIO 모듈을 사용하면 시퀀스를 다양한 형식으로 파일에 저장하거나 읽어올 수 있습니다. 이러한 기능은 5장에서 자세히 다룰 예정입니다.

python
코드 복사
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# SeqRecord 생성
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
my_record = SeqRecord(my_seq, id="XYZ789", description="예시 시퀀스")

# FASTA 형식으로 출력
print(my_record.format("fasta"))
출력 결과는 다음과 같습니다:

shell
코드 복사
>XYZ789 예시 시퀀스
GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA
SeqRecord 객체는 .format() 메서드를 통해 다양한 파일 형식으로 변환이 가능합니다.