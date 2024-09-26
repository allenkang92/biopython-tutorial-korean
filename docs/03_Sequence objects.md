# 3장: 시퀀스 객체

생물학적 시퀀스는 생물정보학에서 가장 중요한 객체 중 하나입니다. 이 장에서는 Biopython에서 시퀀스를 다루는 `Seq` 객체에 대해 설명하겠습니다. 4장에서는 `SeqRecord` 객체를 소개할 예정이며, 이는 시퀀스 정보와 주석을 결합하는 데이터 구조입니다. 5장에서는 시퀀스의 입출력에 대해 다룰 것입니다.

시퀀스는 본질적으로 'AGTACACTGGT'와 같은 문자열입니다. 이는 생물학적 파일 형식에서 시퀀스를 가장 일반적으로 나타내는 방식이기 때문입니다.

Biopython의 `Seq` 객체는 표준 파이썬 문자열과 유사하지만, 생물학적 시퀀스 분석에 필요한 추가적인 기능을 제공합니다. 예를 들어, `translate()` 메서드를 통해 시퀀스를 단백질로 번역하거나, `reverse_complement()` 메서드를 사용하여 역상보적 시퀀스를 생성하는 등의 다양한 메서드를 제공합니다.

## 3.1 시퀀스는 문자열처럼 동작한다

대부분의 경우 `Seq` 객체는 파이썬 문자열과 매우 유사하게 동작합니다. 예를 들어, 시퀀스의 길이를 확인하거나 각 염기(letter)를 반복(iterate)할 수 있습니다.

```python
from Bio.Seq import Seq

my_seq = Seq("GATCG")

for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))

print(len(my_seq))  # 시퀀스의 길이 출력
```

출력 결과는 다음과 같습니다:

```
0 G
1 A
2 T
3 C
4 G
5
```

`Seq` 객체의 각 요소에 접근하는 방식도 문자열과 동일합니다 (파이썬의 인덱스는 0부터 시작한다는 점을 유의하세요):

```python
print(my_seq[0])  # 첫 번째 문자: G
print(my_seq[2])  # 세 번째 문자: T
print(my_seq[-1]) # 마지막 문자: G
```

`Seq` 객체의 `.count()` 메서드를 사용하여 특정 염기 또는 아미노산의 개수를 확인할 수 있습니다. 이는 파이썬 문자열의 `.count()` 메서드와 동일한 방식으로 동작합니다.

```python
from Bio.Seq import Seq

print("AAAA".count("AA"))  # 문자열의 경우: 중첩된 부분까지 포함하여 2 반환
print(Seq("AAAA").count("AA"))  # Seq 객체의 경우: 동일하게 2 반환
```

## 3.2 시퀀스 슬라이싱

파이썬 문자열처럼 `Seq` 객체에서도 슬라이싱(slicing)을 사용하여 시퀀스의 부분 문자열을 추출할 수 있습니다.

```python
from Bio.Seq import Seq

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")

print(my_seq[4:12])  # 출력: 'GATGGGCC'
```

이렇게 슬라이싱하면, 시작 인덱스는 포함되고 끝 인덱스는 포함되지 않는 파이썬의 기본 슬라이싱 규칙을 따릅니다.

시퀀스의 끝에서부터 슬라이싱하거나, 특정 부분만 추출하고 싶을 때도 사용 가능합니다:

```python
print(my_seq[-4:])  # 마지막 4개의 염기: 'TCGC'
print(my_seq[:-4])  # 마지막 4개를 제외한 시퀀스: 'GATCGATGGGCCTATATAGGATCGAAAAT'
```

## 3.3 Seq 객체를 문자열로 변환하기

`Seq` 객체를 일반 파이썬 문자열로 변환하려면 `str()` 함수를 사용합니다:

```python
from Bio.Seq import Seq

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")

print(str(my_seq))  # 출력: 'GATCGATGGGCCTATATAGGATCGAAAATCGC'
```

파이썬의 `print()` 함수를 사용하면 `Seq` 객체가 자동으로 문자열로 변환되어 출력됩니다:

```python
print(my_seq)  # 출력: 'GATCGATGGGCCTATATAGGATCGAAAATCGC'
```

## 3.4 시퀀스 연결하기 (Concatenation)

`Seq` 객체는 `+` 연산자를 사용하여 시퀀스를 쉽게 연결할 수 있습니다:

```python
from Bio.Seq import Seq

seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")

print(seq1 + seq2)  # 출력: 'ACGTAACCGG'
```

여러 시퀀스를 한 번에 연결하고 싶다면 `for` 루프를 사용하거나 `.join()` 메서드를 사용하여 간편하게 연결할 수 있습니다:

```python
list_of_seqs = [Seq("ACGT"), Seq("AACC"), Seq("GGTT")]
concatenated = Seq("").join(list_of_seqs)
print(concatenated)  # 출력: 'ACGTAACCGGTT'
```

## 3.5 시퀀스의 대/소문자 변경

`Seq` 객체는 `.upper()`와 `.lower()` 메서드를 제공하며, 이를 통해 시퀀스의 대소문자를 변환할 수 있습니다.

```python
from Bio.Seq import Seq

dna_seq = Seq("acgtACGT")

print(dna_seq.upper())  # 출력: 'ACGTACGT'
print(dna_seq.lower())  # 출력: 'acgtacgt'
```

## 3.6 뉴클레오타이드 시퀀스의 상보적 서열

DNA 시퀀스는 상보적인 염기쌍으로 이루어져 있습니다. Biopython의 `Seq` 객체는 `.complement()`와 `.reverse_complement()` 메서드를 제공하여 이러한 상보적 시퀀스를 생성할 수 있습니다.

```python
from Bio.Seq import Seq

dna_seq = Seq("GATCG")

print(dna_seq.complement())        # 상보적 시퀀스: 'CTAGC'
print(dna_seq.reverse_complement()) # 역상보적 시퀀스: 'CGATC'
```

`complement()` 메서드는 A와 T, G와 C가 상보적인 염기쌍임을 활용하여 상보적 시퀀스를 반환합니다. `reverse_complement()`는 상보적 시퀀스를 뒤집어 역상보적 시퀀스를 반환합니다.

## 3.7 전사 (Transcription)

DNA 시퀀스는 전사 과정을 통해 RNA로 변환됩니다. `Seq` 객체의 `.transcribe()` 메서드를 사용하면 DNA 시퀀스를 RNA로 전사할 수 있습니다.

```python
from Bio.Seq import Seq

coding_dna = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
print(coding_dna.transcribe())  # 출력: 'GAUCGAUGGGCCUAUAUAGGAUCGAAAAUCGC'
```

DNA의 티민(T)은 전사 과정에서 유라실(U)로 변환됩니다.

## 3.8 역전사 (Back-transcription)

RNA 시퀀스는 때로는 역전사 과정을 통해 DNA로 변환될 수 있습니다. `Seq` 객체의 `.back_transcribe()` 메서드를 사용하여 RNA 시퀀스를 DNA로 변환할 수 있습니다.

```python
from Bio.Seq import Seq

m_rna = Seq("GAUCGAUGGGCCUAUAUAGGAUCGAAAAUCGC")
print(m_rna.back_transcribe())  # 출력: 'GATCGATGGGCCTATATAGGATCGAAAATCGC'
```

## 3.9 번역 (Translation)

`Seq` 객체는 `.translate()` 메서드를 통해 DNA 또는 RNA 시퀀스를 단백질로 번역할 수 있습니다.

```python
from Bio.Seq import Seq

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print(coding_dna.translate())  # 출력: 'MAIVMGR*KGAR*'
```

번역 결과의 `*` 표시는 종결 코돈(Stop Codon)을 의미합니다.

`Seq` 객체는 표준적인 생물학적 번역 테이블을 사용하여 시퀀스를 단백질로 변환합니다. 이러한 테이블은 NCBI의 번역 테이블 번호에 따라 사용할 수 있으며, 이를 통해 다양한 코돈 표를 지원합니다.

## 3.10 알 수 없는 문자 (Ambiguous Characters)

생물학적 시퀀스는 때로는 알 수 없는 문자를 포함할 수 있습니다. 예를 들어, 뉴클레오타이드 시퀀스에 "N"이라는 문자는 A, T, G 또는 C를 의미하는 모호한 염기입니다.

Biopython의 `Seq` 객체는 이러한 모호한 문자를 허용하며, `.complement()`, `.reverse_complement()`, `.transcribe()`, `.translate()` 등 대부분의 메서드는 이러한 문자를 처리할 수 있습니다.

```python
from Bio.Seq import Seq

ambiguous_dna = Seq("ATGCGCN")
print(ambiguous_dna.complement())  # 출력: 'TACGCGN'
print(ambiguous_dna.translate())   # 출력: 'MR'
```

## 3.11 서열의 멈춤 코돈 처리

번역할 때 종결 코돈(stop codon)이 나타나면 일반적으로 단백질 시퀀스에 `*`로 표시됩니다. `.translate()` 메서드는 `to_stop=True` 옵션을 사용하여 첫 번째 종결 코돈에서 번역을 멈출 수 있습니다.

```python
from Bio.Seq import Seq

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
protein_seq = coding_dna.translate(to_stop=True)
print(protein_seq)  # 출력: 'MAIVMGR'
```

위의 예제에서 `to_stop=True`는 번역이 첫 번째 종결 코돈(Stop Codon)에서 멈추게 합니다.

## 3.12 전체적인 번역 테이블 (Translation Tables)

Biopython은 다양한 번역 테이블을 지원합니다. NCBI에서 정의한 번역 테이블을 참조하면 코돈 표에 따라 다른 번역 결과를 얻을 수 있습니다.

기본적으로 `.translate()` 메서드는 표준 번역 테이블을 사용하지만, `table` 인자를 통해 다른 번역 테이블을 지정할 수 있습니다.

```python
from Bio.Seq import Seq

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print(coding_dna.translate(table=2))  # 표 2를 사용한 번역
```

## 3.13 Seq 객체와 불변성 (Immutability)

`Seq` 객체는 불변(immutable)입니다. 즉, 생성 후에 시퀀스 데이터를 변경할 수 없습니다. 이는 시퀀스를 안전하게 다루기 위한 설계 원칙입니다.

만약 시퀀스의 일부를 변경하거나 수정하고 싶다면, 새로운 `Seq` 객체를 생성해야 합니다.

```python
from Bio.Seq import Seq

my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
print(my_seq)
# Seq 객체는 불변이므로 일부를 변경하려면 새로 만들어야 함
modified_seq = my_seq[:5] + "T" + my_seq[6:]
print(modified_seq)
```

## 3.14 MutableSeq 객체

Biopython에서는 `Seq` 객체와 다르게 시퀀스를 수정할 수 있는 `MutableSeq` 객체를 제공합니다. `MutableSeq` 객체는 변경 가능하므로, 시퀀스의 특정 부분을 변경할 수 있습니다.

```python
from Bio.Seq import MutableSeq

mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
print(mutable_seq)
mutable_seq[5] = "T"  # 특정 위치의 염기 변경
print(mutable_seq)
```

`MutableSeq` 객체는 여러 편리한 메서드를 제공하며, `Seq` 객체와 마찬가지로 다양한 시퀀스 분석 기능을 사용할 수 있습니다. 다만, `MutableSeq`는 성능 및 안전성 측면에서 다소 제약이 있을 수 있으므로 필요에 따라 신중히 사용해야 합니다.

## 3.15 문자열과 직접 작업하기

이 장을 마무리하며, 시퀀스 객체를 사용하고 싶지 않거나 함수형 프로그래밍 스타일을 선호하는 분들을 위해, Bio.Seq 모듈에는 일반 Python 문자열, Seq 객체 또는 MutableSeq 객체를 받아들이는 모듈 수준의 함수들이 있습니다:

```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate

my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
print(reverse_complement(my_string))
print(transcribe(my_string))
print(back_transcribe(my_string))
print(translate(my_string))
```

그러나 기본적으로는 `Seq` 객체를 사용하여 작업하는 것이 권장됩니다.

이로써 Biopython의 시퀀스 객체에 대한 기본적인 소개를 마칩니다. 이 장에서 다룬 내용을 통해 생물학적 시퀀스를 효과적으로 다룰 수 있는 기초를 갖추게 되었습니다. 다음 장에서는 이러한 시퀀스 객체를 포함하는 더 복잡한 데이터 구조인 `SeqRecord` 객체에 대해 알아보겠습니다.