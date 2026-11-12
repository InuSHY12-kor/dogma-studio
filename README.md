# GeneCanvas

Windows용 플라스미드·DNA 서열 편집 프로그램입니다. 화면 구성과 작업 흐름은 the reference tool과 비슷하게 만들었습니다.

## 다운로드

- **설치 버전**: `GeneCanvas-Setup-<버전>.exe`
- **설치 없이 실행**: `GeneCanvas-<버전>-portable.exe`

두 파일 모두 GitHub **Actions → Build Windows → 최근 실행 → Artifacts**에서 받을 수 있습니다.
`v0.1.0`처럼 `v`로 시작하는 태그를 push하면 **Releases** 페이지에도 자동으로 올라갑니다.

> 코드 서명이 없어서 처음 실행할 때 Windows SmartScreen 경고가 뜰 수 있습니다. **추가 정보 → 실행**을 누르면 됩니다.

## 기능

| 영역 | 내용 |
|---|---|
| 파일 | GenBank(.gb/.gbk), FASTA, .dna 열기 · GenBank/FASTA 저장 · 지도 SVG 내보내기 · 끌어다 놓기로 열기 · 파일 연결(더블클릭으로 열기) · 여러 파일을 탭으로 열기 |
| **NCBI** | **accession 번호(예: L09137, NM_000546)로 서열 바로 받기 · 단어로 검색해서 목록에서 고르기 · 받은 뒤 흔한 feature 자동 인식** |
| Map | 원형/선형 지도, Feature 화살표, 제한효소 절단 위치 라벨, 프라이머, ORF, 지도에서 드래그로 영역 선택, 확대/축소 |
| Sequence | 이중가닥 표시, 10 bp 눈금, Feature 막대, 아미노산 번역, 효소 위치, 직접 타이핑으로 편집, 복사/붙여넣기, DNA·단백질 검색(IUPAC 지원), 위치로 이동 |
| Enzymes | 효소 77종, 절단 횟수·위치·말단 형태(5′/3′/blunt), 원하는 효소를 골라 지도에 표시 |
| Features | 추가·편집·삭제, 종류·색·방향 지정, 흔한 서열 자동 인식(T7, lac operator, loxP, His tag 등) |
| Primers | 선택 영역으로 프라이머 만들기(위/아래 가닥), Tm(nearest-neighbor), GC%, 결합 위치(5′ 꼬리가 있으면 3′ 15 nt로 찾음) |
| **Cloning** | **제한효소 클로닝(벡터·insert를 자르고 ligation, 양쪽 방향 모두 확인) · Golden Gate(Type IIS: BsaI/BsmBI/BbsI/SapI, overhang으로 순서 결정, 내부 사이트 경고) · Gibson/overlap assembly(2~6 조각, 순서·방향 자동 탐색) · Gateway BP/LR(att 자리에서 재조합, attB 생성 여부 확인) · att 자리 표시** |
| **Tools** | **코돈 최적화(8개 organism, 알고리즘 4종, 효소 사이트 추가/제거, 보존 영역) · 서열 비교(여러 플라스미드 정렬) · PCR 프라이머 자동 설계(Tm·GC·hairpin·dimer·3′ 중복결합 검사, 제한효소 꼬리 붙이기) · 시퀀싱 프라이머 자동 배치 · 돌연변이 유발 프라이머 · .ab1 정렬** |
| **Compare** | **여러 서열을 한 화면에서 비교 · 원형은 회전·역상보까지 맞춰서 정렬 · 치환/삽입/결실을 색으로 표시 · 차이 있는 줄만 보기** |
| **Analyze** | **프라이머 구조(헤어핀·self dimer·cross dimer를 ΔG와 그림으로) · RNA 2차 구조(arc 다이어그램, dot-bracket, stem 목록) · 단백질 성질(분자량·pI·흡광계수·hydropathy·전하 곡선·조성) · 외부 구조 예측 툴 연결** |
| **Alignment** | **.ab1 시퀀싱 결과를 서열에 정렬 · 양쪽 가닥 자동 판별 · 낮은 품질 끝부분 자동 제거 · 크로마토그램 표시 · 차이(치환/결실/삽입) 목록 · 여러 read 비교** |
| Actions | 제한효소 절단 + 가상 아가로스 겔, 번역(6 frame), 역상보, 원점 변경, 원형/선형 전환 |
| **인쇄** | **지도·feature 표·프라이머 표·제한효소 목록·전체 서열·정렬 결과를 골라서 인쇄하거나 PDF로 저장** |
| History | 모든 변경 기록, 원하는 단계로 되돌리기, Undo/Redo |

## 기능별 사용법

**NCBI에서 받기** — `Ctrl+D` 또는 툴바의 NCBI 단추. accession 번호를 넣고 Download, 또는 단어를 넣고 Search 후 목록에서 Open.

**제한효소 클로닝** — Cloning → Restriction Cloning. 벡터와 insert를 고르고 각각 자를 효소를 체크하면 생기는 조각이 양 끝 모양(5′-AATT 같은)과 함께 나옵니다. Simulate Ligation을 누르면 붙을 수 있는 조합만 보여주고, Create로 새 탭에 만듭니다. Feature는 그대로 따라갑니다.

**Gibson assembly** — Cloning → Gibson. 조각을 2~6개 고르고 최소 겹침 길이를 정하면, 순서와 방향을 모두 시도해서 원형으로 닫히는 조합을 찾습니다.

**Gateway** — Cloning → Gateway. 두 플라스미드를 고르면 att 자리를 찾아 7 bp core에서 교차시킵니다. 결과 연결 부위가 실제 attB1/attB2가 되었는지 확인해서 알려줍니다.

**프라이머 설계** — Tools → Design PCR Primers. 증폭할 구간, 목표 Tm, 길이 범위를 정하면 후보 5쌍을 Tm·GC·경고와 함께 보여줍니다. 5′ 꼬리(제한효소 자리 등)를 붙일 수 있고, Tm은 꼬리를 뺀 결합 부위로만 계산합니다.

**시퀀싱 결과 정렬** — Tools → Align Sequencing Traces, 또는 .ab1 파일을 창에 끌어다 놓기. 기준 서열에 맞춰 양쪽 가닥을 다 시도하고, 품질 낮은 끝을 자른 뒤 정렬해서 차이를 모두 표시합니다.

**Golden Gate** — Cloning → Golden Gate. 효소(BsaI 등)를 고르면 각 part의 사이트 개수와 생기는 overhang을 보여줍니다. part는 사이트가 **정확히 2개**여야 하고, 더 있으면 "part 안에 사이트가 있다"고 빨간 글씨로 알려줍니다(먼저 코돈 최적화로 없애면 됩니다). 조립이 끝난 산물에는 그 효소 사이트가 남지 않는데, 이것도 결과 표에서 확인할 수 있습니다.

**서열 비교** — Tools → Compare Sequences. 맨 처음 고른 서열이 기준이 되고 나머지가 거기에 맞춰집니다. 원형 플라스미드는 origin이 달라도, 심지어 역상보로 저장돼 있어도 알아서 맞춥니다. 같은 자리는 점으로, 다른 자리는 색으로 표시합니다.

**코돈 최적화** — Tools → Codon Optimisation. 순서대로 정하시면 됩니다.

1. **무엇을 최적화할지**: 열려 있는 서열의 한 구간, 또는 붙여넣은 단백질 서열
2. **ORF 확인**: `Use` 드롭다운에서 CDS feature나 자동으로 찾은 ORF를 고르면 범위가 정확히 들어갑니다. 바로 아래에 `✓ Valid reading frame — 287 codons…`처럼 확인 결과가 나오고, 길이가 3의 배수가 아니거나 중간에 종결코돈이 있으면 무엇이 문제인지 알려줍니다
3. **Organism**: 8종 중 선택. 각 표가 어떤 데이터로 만들어졌는지 함께 표시됩니다
4. **알고리즘**: 4종 중 선택, 각각 언제 쓰는지 설명이 따라붙습니다
5. **효소 사이트**: 없앨 사이트와 **넣을** 사이트를 각각 고릅니다. 지금 그 구간에 들어 있는 사이트는 목록으로 보여주고, 클릭하면 바로 "없앨 목록"에 들어갑니다
6. **보존 영역**: 건드리면 안 되는 구간(선택 영역이나 feature)을 추가하면 그 안의 코돈은 그대로 둡니다
7. 결과로 CAI·GC·GC3·최장 반복·사이트 수가 **전후 비교**로 나오고, 원래 자리에 바꿔 넣거나 새 탭으로 열 수 있습니다

**프라이머 구조 보기** — Analyze → Primer Structures, 또는 Primers 표에서 `Structure`. 헤어핀, self dimer, cross dimer를 사다리 그림과 ΔG로 보여줍니다.

```
5'-GCGCGCGC\
   |||||||| TTTT
3'-CGCGCGCG/
   ΔG -11.92 kcal/mol
```

ΔG는 Tm과 같은 SantaLucia 1998 nearest-neighbor 값에서 계산합니다(37 °C 기준). 대략 −6 kcal/mol보다 약하면 문제되는 일이 드물고, −9보다 강하면 대개 문제가 되며 특히 3′ 끝이 물려 있을 때 그렇습니다. 마지막 5염기 중 몇 개가 짝지었는지도 함께 표시합니다.

**RNA 2차 구조** — Analyze → RNA Secondary Structure. 구간을 정하면 접은 구조를 arc 다이어그램, stem 목록, dot-bracket으로 보여줍니다. T는 U로 읽습니다.

**단백질 성질** — Analyze → Protein Properties. CDS feature나 ORF를 고르거나 단백질 서열을 붙여넣으면 분자량, pI, 280 nm 흡광계수(1 mg/mL A280 포함), 순전하, GRAVY, aliphatic index, Kyte–Doolittle hydropathy 그래프, pH별 전하 곡선, 아미노산 조성을 보여줍니다.

알고리즘 4종:

| 이름 | 내용 |
|---|---|
| Most frequent codon | 항상 그 organism에서 가장 흔한 코돈. CAI가 가장 높지만 같은 코돈이 반복되고 GC가 치우칠 수 있음 |
| Frequency sampling, rare codons dropped | 기준치(기본 10%) 아래 희귀 코돈을 뺀 뒤 자연 빈도대로 섞음. 일반적인 발현용 기본값 |
| Frequency sampling | 자연 빈도 그대로 섞음 |
| Harmonised to the source organism | 원래 organism에서의 희귀/흔함 정도를 목표 organism에서도 비슷하게 맞춤. 번역 속도 리듬이 folding에 중요할 때 |

## 직접 빌드하기

```bash
npm install
npm start            # 개발 모드로 실행
npm run dist:win     # Windows에서 설치 파일과 portable exe 만들기
```

## 구조

```
main.js         Electron 메인 프로세스 (창, 파일 대화상자, NCBI 요청, 인쇄/PDF)
preload.js      화면에서 쓰는 window.gcHost API
app/index.html  화면 틀
app/style.css   스타일
app/dna.js      DNA 로직 (효소, 번역, Tm, ORF, 파일 형식)
app/clone.js    클로닝 (제한효소 ligation, Golden Gate, Gibson, Gateway)
app/trace.js    .ab1 읽기와 서열 정렬
app/primer.js   프라이머 설계
app/msa.js      여러 서열 비교 (앵커 기반 정렬)
app/codon.js    코돈 최적화
app/struct.js   프라이머 구조, RNA 접기, 단백질 성질
app/codon-tables.js  코돈 사용표 (자동 생성 파일)
app/app.js      화면과 상호작용
tools/build-codon-tables.js  NCBI CDS에서 코돈 사용표를 다시 만드는 스크립트
```

## 코돈 사용표에 대하여

`app/codon-tables.js`는 **실제 NCBI 코딩 서열에서 직접 계산한 값**입니다. 손으로 적은 숫자가 아닙니다.

| Organism | 출처 |
|---|---|
| Escherichia coli K-12 | 유전체 NC_000913.3의 전체 CDS (4,269개 유전자, 1,328,211 코돈) |
| Homo sapiens / Mus musculus | RefSeq Select (유전자당 대표 전사체 1개) 598개 |
| S. cerevisiae, CHO, P. pastoris, Drosophila, Arabidopsis | RefSeq mRNA 481~600개 |

시작·종결 코돈은 선택의 여지가 없으므로 집계에서 뺐고, 길이가 3의 배수가 아니거나 중간에 종결코돈이 있는 불완전한 기록도 제외했습니다.
다시 만들려면 `node tools/build-codon-tables.js`를 실행하면 됩니다(인터넷 필요, 몇 분 소요).

## 계산의 근거와 한계

어디까지가 정확한 계산이고 어디부터가 근사인지 구분해 두었습니다.

| 항목 | 상태 |
|---|---|
| Tm, 프라이머 헤어핀·dimer ΔG | **정확**. SantaLucia 1998 nearest-neighbor 값에서 계산하며, 출판된 ΔG°37 값과 0.02 kcal/mol 안에서 일치하는 것을 확인했습니다 |
| 단백질 분자량, 280 nm 흡광계수, GRAVY, aliphatic index | **정확**. 표준 공식 그대로이며 수작업 계산과 대조했습니다 |
| 단백질 pI | pKa 값 묶음에 따라 달라집니다. EMBOSS 묶음을 쓰며 화면에도 표시합니다. 다른 툴과 0.2~0.5 정도 차이 날 수 있는데 이는 정상입니다 |
| 코돈 사용표 | **실측 데이터**. NCBI 실제 CDS에서 계산 (아래 표 참고) |
| **RNA 2차 구조 ΔG** | **근사입니다.** 탐색 알고리즘 자체는 표준 최소자유에너지(Zuker) 방식이지만, 에너지 값은 전체 Turner 파라미터 세트가 아닌 단순화된 모델을 씁니다. 구조의 대략적인 모양을 보는 데는 쓸 만하지만, 논문에 쓸 숫자는 ViennaRNA RNAfold로 확인해 주세요. 앱 안에서도 이 점을 표시하고 링크를 걸어 두었습니다 |
| 단백질·RNA 3차 구조 | 예측하지 않습니다. AlphaFold DB, ESMFold, SWISS-MODEL 등으로 연결만 합니다 |

검증 삼아 효모 tRNA-Phe를 접어 보면 acceptor stem(1–7 / 66–72)을 정확히 찾습니다. 접기 속도는 100 nt 0.14초, 200 nt 0.5초이며, 길이의 세제곱에 비례하므로 앱 안에서는 600 nt로 제한하고 그보다 길면 외부 서버로 안내합니다.

## 알아두실 점

- 예제 플라스미드(pGC-Demo)는 기능을 보여주려고 만든 합성 서열이며, 실제 벡터 서열이 아닙니다.
- NCBI 요청은 공개 E-utilities API를 씁니다. API 키 없이 쓰므로 요청 간격을 자동으로 띄웁니다. 받은 서열 말고는 아무것도 보내지 않습니다.
- 저장은 GenBank 형식입니다. the reference tool `.dna` 형식으로는 저장할 수 없습니다(읽기만 됩니다).
- Tm은 nearest-neighbor 모델(50 mM Na⁺, 250 nM primer) 기준입니다.
- Gateway 시뮬레이션은 att 자리의 7 bp core에서 교차시키는 방식입니다. 표준이 아닌 att 자리를 쓰면 연결 부위가 attB가 아니라고 알려줍니다.
