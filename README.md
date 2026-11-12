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
| **Cloning** | **제한효소 클로닝(벡터·insert를 자르고 ligation, 양쪽 방향 모두 확인) · Gibson/overlap assembly(2~6 조각, 순서·방향 자동 탐색) · Gateway BP/LR(att 자리에서 재조합, attB 생성 여부 확인) · att 자리 표시** |
| **Tools** | **PCR 프라이머 자동 설계(Tm·GC·hairpin·dimer·3′ 중복결합 검사, 제한효소 꼬리 붙이기) · 시퀀싱 프라이머 자동 배치 · 돌연변이 유발 프라이머 · .ab1 정렬** |
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
app/clone.js    클로닝 (제한효소 ligation, Gibson, Gateway)
app/trace.js    .ab1 읽기와 서열 정렬
app/primer.js   프라이머 설계
app/app.js      화면과 상호작용
```

## 알아두실 점

- 예제 플라스미드(pGC-Demo)는 기능을 보여주려고 만든 합성 서열이며, 실제 벡터 서열이 아닙니다.
- NCBI 요청은 공개 E-utilities API를 씁니다. API 키 없이 쓰므로 요청 간격을 자동으로 띄웁니다. 받은 서열 말고는 아무것도 보내지 않습니다.
- 저장은 GenBank 형식입니다. the reference tool `.dna` 형식으로는 저장할 수 없습니다(읽기만 됩니다).
- Tm은 nearest-neighbor 모델(50 mM Na⁺, 250 nM primer) 기준입니다.
- Gateway 시뮬레이션은 att 자리의 7 bp core에서 교차시키는 방식입니다. 표준이 아닌 att 자리를 쓰면 연결 부위가 attB가 아니라고 알려줍니다.
