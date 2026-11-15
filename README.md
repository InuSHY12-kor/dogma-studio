# GeneCanvas

플라스미드·DNA 서열 편집 프로그램입니다. 화면 구성과 작업 흐름은 the reference tool과 비슷하게 만들었고,
**Windows 데스크톱**과 **안드로이드** 두 가지로 나옵니다. 둘은 `app/` 폴더 하나를 같이 쓰기 때문에
기능을 한 번 만들면 양쪽에 함께 들어갑니다.

## 다운로드

**Windows**

- **설치 버전**: `GeneCanvas-Setup-<버전>.exe`
- **설치 없이 실행**: `GeneCanvas-<버전>-portable.exe`

**안드로이드**

- `GeneCanvas-<버전>.apk` — 받아서 설치하면 됩니다(스토어 배포가 아니므로 "출처를 알 수 없는 앱 설치"를 한 번 허용해야 합니다).

세 파일 모두 GitHub **Actions → Build Windows / Build Android → 최근 실행 → Artifacts**에서 받을 수 있습니다.
최신 판은 **Releases** 페이지(`v0.3.0`)에서 바로 받을 수 있습니다. `v`로 시작하는 태그를 push하면 자동으로 올라갑니다.

> 코드 서명이 없어서 Windows에서는 처음 실행할 때 SmartScreen 경고가 뜰 수 있습니다. **추가 정보 → 실행**을 누르면 됩니다.
> 안드로이드 APK도 디버그 서명이라 설치할 때 경고가 뜹니다.

### 폰에서 쓰는 법

마우스가 없는 대신 손가락 동작이 같은 일을 합니다.

| 동작 | 하는 일 |
|---|---|
| 탭 | 그 자리를 고름 — feature, 효소, 커서 위치 |
| 길게 누르기 (feature·효소 위에서) | 그 항목의 메뉴 |
| 길게 누르고 끌기 (빈 곳에서) | 구간 선택 |
| 그냥 끌기 | 스크롤 |
| 두 손가락 오므리기/벌리기 | 맵 확대·축소 |

맵·시퀀스 옆 사이드 메뉴는 엄지가 닿는 오른쪽 아래에 있고, 폰에서는 효소 이름에서 위치 숫자를 빼서
좁은 화면에서도 이름이 겹치지 않게 했습니다(사이드 메뉴 → 줄 번호에서 다시 켤 수 있습니다).
별도 OS 창으로 열리던 Feature·Enzyme 데이터베이스, Cloning, Simulation 창은 폰에서는 같은 화면 안
대화상자로 열립니다.

## 기능

| 영역 | 내용 |
|---|---|
| 파일 | GenBank(.gb/.gbk), FASTA, .dna 열기 · GenBank/FASTA 저장 · 지도 SVG 내보내기 · 끌어다 놓기로 열기 · 파일 연결(더블클릭으로 열기) · 여러 파일을 탭으로 열기 |
| **NCBI** | **accession 번호(예: L09137, NM_000546)로 서열 바로 받기 · 단어로 검색해서 목록에서 고르기 · 받은 뒤 흔한 feature 자동 인식** |
| Map | 원형/선형 지도, Feature 화살표, 제한효소 절단 위치 라벨, 프라이머, ORF, 지도에서 드래그로 영역 선택, 확대/축소 |
| Sequence | 이중가닥 표시, 10 bp 눈금, Feature 막대, 아미노산 번역, 효소 위치, 직접 타이핑으로 편집, 복사/붙여넣기, DNA·단백질 검색(IUPAC 지원), 위치로 이동 |
| Enzymes | 효소 203종, 절단 횟수·위치·말단 형태(5′/3′/blunt), 효소별 표시 토글, Numbers/Lines, noncutters, **단백질을 바꾸지 않고 효소 사이트 넣기·빼기(silent mutagenesis)** |
| Features | 추가·편집·삭제, 종류·색·4방향 지정, **세그먼트와 인트론(점선 표시, 번역이 엑손을 건너뜀)**, cleavage arrow, 길이 0인 point feature, **유전암호 14종**, 일괄 편집, 표시/숨김, 데이터 내보내기, 흔한 서열 자동 인식(T7, lac operator, loxP, His tag 등), **나만의 feature 종류(custom type)** |
| Primers | 선택 영역으로 프라이머 만들기(위/아래 가닥), **효소 사이트·코돈·펩타이드 삽입 도구**, Tm(nearest-neighbor), GC%, MW, **결합 조건 직접 설정(3′ 최소 일치 염기·최소 Tm)**, 표시/숨김, 일괄 편집, 리스트 가져오기, 주문용 내보내기, **프로그램 단위 프라이머 데이터베이스(라벨로 세트 구성)** |
| **Cloning** | **한 창에서 방법만 바꿔 쓰는 Vector · Fragment · Product 창(드롭다운으로 알고리즘 선택, 지도에서 효소·feature·프라이머를 직접 클릭) · 제한효소 클로닝(효소 개수 제한 없음) · Golden Gate(Type IIS: BsaI/BsmBI/BbsI/SapI, overhang으로 순서 결정, 내부 사이트 경고) · Gibson/overlap assembly(2~6 조각, 순서·방향 자동 탐색) · Gateway BP/LR(att 자리에서 재조합, attB 생성 여부 확인) · att 자리 표시** |
| **Tools** | **코돈 최적화(8개 organism, 알고리즘 4종, 효소 사이트 추가/제거, 보존 영역) · 서열 비교(여러 플라스미드 정렬) · PCR 프라이머 자동 설계(Tm·GC·hairpin·dimer·3′ 중복결합 검사, 제한효소 꼬리 붙이기) · 시퀀싱 프라이머 자동 배치 · 돌연변이 유발 프라이머 · .ab1 정렬** |
| **Compare** | **여러 서열을 한 화면에서 비교 · 원형은 회전·역상보까지 맞춰서 정렬 · 치환/삽입/결실을 색으로 표시 · 차이 있는 줄만 보기** |
| **파일 종류** | **DNA · RNA · Protein 파일을 각각 만들 수 있고**, 열려 있는 파일의 종류에 따라 상단 메뉴와 탭이 바뀜 |
| **Import** | 다른 파일에서 Feature/Primer 가져오기 · 최근 파일 30개 + 찾아보기 · **서열을 맞춰서 위치를 찾아 넣음**(좌표가 달라도 됨) |
| **Enzymes** | 효소 **203종** · 한 화면에 다 보이는 선택 창 · **나만의 효소 세트 저장**(파일이 아니라 프로그램에 저장돼서 모든 파일에 적용) · 효소/Feature 데이터베이스 열람과 CSV 내보내기 |
| **Gel** | **아가로스 젤을 별도 탭 문서로** · 0.5~3.0% · TAE/TBE · running 정도 조절 · ladder 15종(NEB, Thermo, Invitrogen, Promega, λ, φX) + **직접 입력한 ladder 저장** · lane 1~20개(개수만큼 젤이 넓어짐) |
| **Analyze** | **프라이머 구조(헤어핀·self dimer·cross dimer를 ΔG와 그림으로) · RNA 2차 구조(arc 다이어그램, dot-bracket, stem 목록) · 단백질 성질(분자량·pI·흡광계수·hydropathy·전하 곡선·조성) · 외부 구조 예측 툴 연결** |
| **Alignment** | **.ab1 시퀀싱 결과를 서열에 정렬 · 양쪽 가닥 자동 판별 · 낮은 품질 끝부분 자동 제거 · 크로마토그램 표시 · 차이(치환/결실/삽입) 목록 · 여러 read 비교** |
| Actions | 제한효소 절단 + 가상 아가로스 겔, 번역(6 frame), 역상보, 원점 변경, 원형/선형 전환 |
| **인쇄** | **지도·feature 표·프라이머 표·제한효소 목록·전체 서열·정렬 결과를 골라서 인쇄하거나 PDF로 저장** |
| History | 모든 변경 기록, 원하는 단계로 되돌리기, Undo/Redo, **제작 이력 트리**, 텍스트로 보기·내보내기·복사, **단계별로 바뀐 염기 색칠**, 이력 다듬기(Trim) |
| **Search** | **`(AG)` `R` `!R` `"R"` 문법 · 정확/유사(미스매치·갭)/단백질/이름 4개 모드 · 일치 항목 색칠 · 미스매치 빨간 박스 · 최근 검색** |
| **보기** | **Split Window(좌우/상하) · Focus on Region(구간만 보기) · GC content 플롯 · 지도 회전 · 글꼴 크기 · Compact/Plain text · Blocks of 3·10 · 고정 줄너비 · 미니맵 · 색 스킴 6종** |
| **Info** | **파일 메타데이터(출처·host·strain·methylation·accession·author·comments) · 참고문헌(PubMed PMID로 가져오기) · 파일 첨부(embedded files)** |

## 기능별 사용법

**DNA / RNA / Protein 파일** — File → New DNA / New RNA / New Protein File. 같은 창에 탭으로 열리고, 보고 있는 파일의 종류에 따라 상단 메뉴가 바뀝니다.

| 파일 | 상단 메뉴 | 탭 |
|---|---|---|
| DNA | File, Edit, View, Features, Primers, Cloning, Tools, Analyze, Actions, Help | Map, Sequence, Enzymes, Features, Primers, Alignment, Compare, History |
| RNA | File, Edit, View, Features, **RNA**, Help | Map, Sequence, Features, **Structure**, History |
| Protein | File, Edit, View, Features, **Protein**, Help | Sequence, Features, **Properties**, History |

RNA 파일은 화면에 U로 보이고 상보가닥을 그리지 않습니다(내부에는 T로 저장해서 ORF·번역·Feature 기능이 그대로 동작합니다). Structure 탭에서 바로 접을 수 있고, Protein 파일은 Properties 탭에서 분자량·pI·흡광계수·hydropathy·전하 곡선·조성을 보여줍니다. Protein 메뉴의 `Back-translate and Codon-optimize…`를 누르면 그 단백질이 코돈 최적화 창에 미리 채워집니다.

**Feature / Primer 가져오기** — Features → Import Features from Another File, Primers → Import Primers. 새 창 위쪽에 **지금 열려 있는 파일**과 **최근 파일 30개**가 나오고 맨 아래 `Browse…`가 있습니다. Feature는 **자기 서열을 대상 파일에서 찾아서** 넣기 때문에 번호가 달라도, origin이 돌아가 있어도, 가닥이 뒤집혀 있어도 제자리에 들어갑니다(원래 좌표 그대로 넣는 선택지도 있습니다).

**Primer 세트 / 효소 세트** — Primers → Primer Sets, Analyze → Enzyme Sets. 이름을 붙여 저장하고, 이름 바꾸기·현재 선택으로 갱신·다른 파일에 적용·삭제가 됩니다. 파일이 아니라 프로그램에 저장돼서 모든 파일에서 쓸 수 있습니다.

**Feature 속성과 나만의 Feature** — Feature를 더블클릭하면 이름, 종류, **Molecule(DNA/RNA/protein)**, 방향(정방향/역방향/양방향), 색, **번역 표시와 reading frame**, **단백질 일치 인식**, note, qualifier(`gene=bla` 형식)를 지정할 수 있습니다. `Save to My Features`로 저장하면 Detect Common Features가 함께 찾아 줍니다. 단백질 일치를 켜 두면 **코돈이 전부 바뀐 서열에서도 같은 단백질이면 인식합니다** — 코돈 최적화한 유전자에서 태그나 CDS를 찾을 때 쓸 수 있습니다.


**효소 세트 저장** — Map/Sequence 탭 왼쪽의 `Choose Enzymes…`, 또는 Analyze → Choose Enzymes. 한 화면에 모든 효소가 보이고 각각 몇 번 자르는지 숫자가 붙습니다. 체크한 뒤 `Save as my set…`으로 이름을 주면 저장되고, **다른 파일을 열어도 그 세트가 그대로 목록에 나옵니다**. 서열에서 자동으로 계산되는 세트(unique cutter, 6+ cutter, Type IIS, blunt 등)도 함께 제공합니다.

**Agarose 시뮬레이션** — 툴바의 **Agarose** 단추. 파일을 연 것처럼 새 탭으로 열립니다. lane은 최대 12개이고, lane마다 서열 하나와 **digest 또는 PCR**을 고릅니다. 창 안 지도에서 효소를 클릭하면 그 효소로 자르고 Shift-클릭으로 double digest, `amplify with`로 바꾸면 프라이머를 클릭하고 마주보는 것을 Shift-클릭해서 증폭 산물을 올립니다. 젤은 그때그때 다시 그려지고, `Open as an Agarose document`로 지금 구성을 아래의 아가로스 젤 문서로 넘길 수 있습니다.

**아가로스 젤(문서)** — Analyze → New Agarose Gel. 새 탭으로 열립니다. 젤 농도, 버퍼(TAE/TBE), running 정도, lane 개수(1~20)를 고르면 lane 수에 맞춰 젤이 넓어집니다. 각 lane에 ladder나 열려 있는 서열을 넣고, 서열 lane은 `Cut with` 단추로 효소를 골라 digest할 수 있습니다. `Add a ladder…`로 원하는 brand의 band 크기를 직접 입력해 저장하면 그 ladder도 목록에 추가됩니다(`3000*`처럼 별표를 붙이면 진한 기준 band).

**지도 옆 사이드 메뉴** — Map·Sequence 탭 오른쪽의 세로 아이콘 줄입니다. 아이콘을 누르면 그 항목을 켜고 끄고, 옆의 작은 펼침 단추를 누르면 그 그룹의 나머지 설정이 **아이콘을 가리지 않게 왼쪽에** 열립니다. 순서는 **효소 · feature · primer · ORF/번역 · 색 · 줄 번호 · alignment · 화면 분할**이고, 효소 세트와 저장된 세트, feature 라벨 위치, 프라이머 hybridisation, ORF 최소 길이와 번역 옵션, 염기 색과 GC 플롯, 눈금자·줄 너비·blocks까지 여기서 끝납니다.

같은 사이드 메뉴가 **클로닝 창과 Agarose 시뮬레이션의 지도 옆에도** 있습니다. 각 패널에서 지금 지도에 그려지는 효소·feature·primer를 체크 목록으로 직접 고를 수 있고(Show all / Hide all 포함), 설정은 탭별·서열별로 따로 기억됩니다.

**제작 이력(History)** — History 탭의 `How it was made`. 이 벡터가 **어떤 기법으로, 어떤 조각이 합쳐져서** 만들어졌는지 트리로 보여줍니다. 각 조각은 미니맵과 크기가 함께 나오고, `Open`을 누르면 그 조각의 **서열이 그대로 새 탭으로 열립니다**. 제한효소 클로닝, Gibson, Golden Gate, Gateway, 코돈 최적화로 만든 결과에 자동으로 기록되며, **GenBank로 저장해도 서열까지 남습니다**(COMMENT에 보관).


**NCBI에서 받기** — `Ctrl+D` 또는 툴바의 NCBI 단추. accession 번호를 넣고 Download, 또는 단어를 넣고 Search 후 목록에서 Open.

**클로닝(한 창에서 전부)** — Cloning 메뉴의 어느 항목을 눌러도 같은 **Vector · Fragment · Product** 창이 열리고, 창 오른쪽 위 드롭다운에서 **제한효소 클로닝 · Golden Gate · Gibson · In-Fusion · Gateway BP/LR**을 고르면 적용되는 알고리즘이 바뀝니다. 방법을 바꿔도 고른 조각은 그대로 두고 접합부만 다시 판정합니다.

- Vector/Fragment 탭에서는 **지도에서 직접** 고릅니다. 효소를 클릭하면 그 효소로 자르고 Shift-클릭으로 효소를 더합니다(개수 제한 없음). feature를 클릭하면 그 구간을 PCR로 증폭하고, 프라이머를 클릭하면 한쪽 끝이 정해집니다.
- Product 탭의 지도도 클릭해서 선택하고 복사할 수 있습니다(선택 구간, 또는 산물 전체).
- Golden Gate는 part의 사이트 개수와 생기는 overhang을 보고 **overhang이 겹치거나 part 안에 사이트가 남으면** 알려줍니다. Gateway는 att core에서 교차시키고 연결 부위가 attB1/attB2가 되었는지, 부산물이 몇 bp인지 적습니다.
- 창은 **완전히 별개의 OS 창**이라 띄워 둔 채 본 창에서 다른 작업을 할 수 있습니다.

**프라이머 설계** — Tools → Design PCR Primers. 증폭할 구간, 목표 Tm, 길이 범위를 정하면 후보 5쌍을 Tm·GC·경고와 함께 보여줍니다. 5′ 꼬리(제한효소 자리 등)를 붙일 수 있고, Tm은 꼬리를 뺀 결합 부위로만 계산합니다.

**시퀀싱 결과 정렬** — Tools → Align Sequencing Traces, 또는 .ab1 파일을 창에 끌어다 놓기. 기준 서열에 맞춰 양쪽 가닥을 다 시도하고, 품질 낮은 끝을 자른 뒤 정렬해서 차이를 모두 표시합니다.

**서열 비교** — Tools → Compare Sequences. 맨 처음 고른 서열이 기준이 되고 나머지가 거기에 맞춰집니다. 원형 플라스미드는 origin이 달라도, 심지어 역상보로 저장돼 있어도 알아서 맞춥니다. 같은 자리는 점으로, 다른 자리는 색으로 표시합니다.

**코돈 최적화** — Tools → Codon Optimization. 순서대로 정하시면 됩니다.

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

ΔG는 Tm과 같은 SantaLucia 1998 nearest-neighbor 값에서 계산합니다. 기본은 37 °C이고, `Calculation parameters…`에서 온도를 바꾸면 그 온도에서 다시 채점합니다(같은 self-dimer가 25 °C에서 −6.5, 60 °C에서 −2.7 kcal/mol). Monomer와 Homodimer를 탭으로 나눠서 봅니다. 대략 −6 kcal/mol보다 약하면 문제되는 일이 드물고, −9보다 강하면 대개 문제가 되며 특히 3′ 끝이 물려 있을 때 그렇습니다. 마지막 5염기 중 몇 개가 짝지었는지도 함께 표시합니다.

**RNA 2차 구조** — Analyze → RNA Secondary Structure. 구간을 정하면 접은 구조를 arc 다이어그램, stem 목록, dot-bracket으로 보여줍니다. T는 U로 읽습니다.

**단백질 성질** — Analyze → Protein Properties. CDS feature나 ORF를 고르거나 단백질 서열을 붙여넣으면 분자량, pI, 280 nm 흡광계수(1 mg/mL A280 포함), 순전하, GRAVY, aliphatic index, Kyte–Doolittle hydropathy 그래프, pH별 전하 곡선, 아미노산 조성을 보여줍니다.

알고리즘 4종:

| 이름 | 내용 |
|---|---|
| Most frequent codon | 항상 그 organism에서 가장 흔한 코돈. CAI가 가장 높지만 같은 코돈이 반복되고 GC가 치우칠 수 있음 |
| Frequency sampling, rare codons dropped | 기준치(기본 10%) 아래 희귀 코돈을 뺀 뒤 자연 빈도대로 섞음. 일반적인 발현용 기본값 |
| Frequency sampling | 자연 빈도 그대로 섞음 |
| Harmonised to the source organism | 원래 organism에서의 희귀/흔함 정도를 목표 organism에서도 비슷하게 맞춤. 번역 속도 리듬이 folding에 중요할 때 |

## the reference tool 사용자 가이드 기준

v0.2.0에서는 the reference tool 사용자 가이드 15개 챕터(131개 문서)를 전부 읽고 빠진 기능을 채웠고,
v0.3.0에서는 그 화면들을 the reference tool 쪽에 더 가깝게 다시 만들었습니다.
챕터별로 무엇이 들어갔는지는 [CHANGELOG.md](CHANGELOG.md)에 적어 두었습니다.

| 챕터 | 다룬 내용 |
|---|---|
| 01 Projects, 02·03·04 Display Options | Split Window, Focus on Region, GC 플롯, 지도 회전, 글꼴 크기, Compact/Plain text, Blocks, 고정 줄너비, 미니맵, 단백질 Region/Site/Bond |
| 05 Info Panel | 파일 메타데이터, 참고문헌(PubMed), 첨부 파일 |
| 06 Enzymes | Enzymes 메뉴 전체, 효소별 표시, noncutters, silent mutagenesis |
| 07·08 Features | 세그먼트·인트론, 번역 옵션, point feature, cleavage arrow, custom feature type |
| 09·10 Primers | 결합 조건, 삽입 도구, 구조 탭, 프라이머 데이터베이스 |
| 11 Translations | 리딩프레임 행, 유전암호 14종, Letter Codes, Genetic Code Tables |
| 12 Colors | 염기·단백질 색 스킴, 직접 칠하기, History colours |
| 13 Searching | 질의 문법, 4개 검색 모드, feature 메타데이터 검색 |
| 14 Focus on Region | 포커스 창과 모든 뷰 연동 |
| 15 History | 텍스트 보기·내보내기·복사·인쇄, 단계별 색칠, Trim |

## 직접 빌드하기

```bash
npm install
npm start               # 데스크톱을 개발 모드로 실행
npm run dist:win        # Windows 설치 파일과 portable exe 만들기
```

안드로이드는 Android SDK와 JDK 21이 있으면 같은 저장소에서 바로 빌드됩니다.

```bash
npm run android:sync    # app/ 폴더를 안드로이드 프로젝트로 복사
npm run android:apk     # 디버그 APK 만들기 (android/app/build/outputs/apk/debug/)
npm run android:open    # Android Studio에서 열기
```

`app/`을 고치면 두 빌드에 모두 반영됩니다. `npm run android:sync`는 그 복사본을 갱신하는 것뿐이라,
웹 부분을 따로 컴파일하는 단계는 없습니다.

## 구조

```
main.js         Electron 메인 프로세스 (창, 파일 대화상자, NCBI 요청, 인쇄/PDF)
preload.js      데스크톱이 화면에 넘겨주는 window.gcHost API
app/mobileboot.js  같은 window.gcHost를 안드로이드에서 채워 넣는 쪽
                   (시스템 파일 선택기, Documents 폴더와 공유 시트, 설정 저장, NCBI 직접 호출)
app/mobile.css  폰 화면 레이아웃 (전부 미디어 쿼리 안이라 데스크톱에는 영향 없음)
app/touch.js    손가락 동작을 앱이 아는 마우스 이벤트로 바꿔 주는 층
capacitor.config.json / android/   안드로이드 껍데기 (Capacitor)
app/index.html  화면 틀
app/style.css   스타일
app/dna.js      DNA 로직 (효소, 번역, Tm, ORF, 파일 형식)
app/clone.js    클로닝 (제한효소 ligation, Golden Gate, Gibson, Gateway)
app/trace.js    .ab1 읽기와 서열 정렬
app/primer.js   프라이머 설계
app/msa.js      여러 서열 비교 (앵커 기반 정렬)
app/codon.js    코돈 최적화
app/silent.js   효소 사이트를 단백질 그대로 넣기·빼기 (silent mutagenesis)
app/search.js   검색 (the reference tool 질의 문법, 정확·유사·단백질 검색)
app/assembly.js 클로닝 엔진 (벡터 준비, overlap 설계, 제품 조립)
app/asmui.js    Vector/Fragment/Product 창 (본 창과 별도 창에서 모두 동작)
app/struct.js   프라이머 구조, RNA 접기, 단백질 성질
app/gel.js      아가로스 젤 이동 모델과 ladder 데이터
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
| Tm, 프라이머 헤어핀·dimer ΔG | **정확**. SantaLucia 1998 nearest-neighbor 값에서 계산하며, 출판된 ΔG°37 값과 0.02 kcal/mol 안에서 일치하는 것을 확인했습니다. 온도를 바꾸면 그 온도의 ΔH·ΔS에서 다시 계산합니다 |
| DNA/RNA 분자량·농도 계산 | **정확**. 뉴클레오타이드 평균 질량에서 결합마다 물 하나를 빼는 방식이며, 20-mer가 표준 oligo 공식과 0.2 Da 안에서 일치하고 ng/µL ↔ nM 변환도 수계산과 대조했습니다 |
| 유전암호 14종 | NCBI translation table 1–6, 9–14, 16, 21을 표준 코드와의 차이로 적어 두었습니다 |
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
