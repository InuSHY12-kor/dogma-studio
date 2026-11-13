# 변경 내역

## v0.2.0

the reference tool 사용자 가이드 15개 챕터(131개 문서)를 전부 읽고 빠진 기능을 채운 판입니다.
작업 단위는 챕터별로 커밋되어 있고, 각 항목은 headless Chromium에서 동작을 확인했습니다.

### Cloning

- **Gibson Assembly / Restriction Cloning을 Vector · Fragment · Product 3탭 창으로 재작성.**
  Vector 탭에서 backbone과 자르는 방식(제한효소 / PCR / 그대로)을 고르고, Fragment 탭의
  **+ / −** 버튼으로 조각을 늘리면 1, 2, 3으로 번호가 붙습니다. Product 탭은 벡터를 검은
  테두리로, 조각을 서로 다른 색으로 그리고 각 접합부의 overlap 길이와 Tm, 조각 아래 방향
  화살표를 함께 보여줍니다. 조건이 맞으면 초록색 **Ready for Gibson Assembly**, 아니면
  노란색으로 이유를 적습니다.
- **Frameshift 경고는 클로닝이 CDS 두 개를 실제로 하나의 ORF로 이을 때만** 뜨고, 어느 CDS와
  어느 CDS가 몇 bp 어긋나는지 tooltip에 적습니다. 서로 다른 ORF는 경고하지 않습니다.
- **Choose Overlapping Primers**로 overlap 길이를 직접 정할 수 있습니다.
- Assembly 창은 **완전히 별개의 OS 창**이라, 띄워 둔 채로 본 창에서 다른 작업을 할 수 있습니다.

### Enzymes (가이드 06)

- Enzymes 전용 메뉴(2단 서브메뉴 지원), 효소별 표시 토글, 사이트 금색 하이라이트
- Enzymes 뷰: Numbers / Lines, 이름·사이트 수·첫 위치·선택영역과의 거리·표시 여부 정렬,
  행을 펼치면 그 효소 단독 절단 조각 크기 표시(Alt-클릭은 전체)
- Noncutters 창
- **Silent Mutations** — 단백질을 바꾸지 않고 효소 사이트를 넣거나 빼기. 넣을 때는 효소 세트·
  기존 사이트 유무·허용 염기 변경 수·인식 길이로 후보를 걸러 목록화하고, 적용 방식 4가지
  (서열 변경 / 변경 + 돌연변이 프라이머 / feature 주석 / 주석 + 프라이머)와 CSV 내보내기.
  뺄 때는 사이트에 걸친 코돈마다 동의 코돈을 사용빈도와 함께 나열하고, 사이트를 남기는
  코돈은 선택 불가로 표시합니다.

### Features (가이드 07, 08)

- **Feature 세그먼트와 인트론** — 점선으로 표시되고 번역이 엑손을 건너뜁니다. 접합부에서
  쪼개진 코돈도 맞게 읽고, 역방향 CDS도 됩니다. GenBank `join(...)`으로 왕복합니다.
- Feature Translation Options — reading frame, **유전암호 14종(NCBI 표)**, 번역 번호 시작값,
  비-ATG 개시코돈을 Met으로 읽기, 세그먼트별 번호 건너뛰기
- 인트론 — 소문자로 쓰기 / 소문자 영역을 인트론으로 주석 / **Splice to Remove Introns**
- Cleavage arrow, **point feature**(길이 0), 4방향 directionality, prioritise display
- Edit Features 일괄 편집, 표시/숨김, Export Feature Data(13개 열), Add to Common Features
- Features 뷰: 표시 토글·체크박스·상세 펼침(세그먼트/인트론/cleavage/산물 크기/qualifier)·
  열 정렬·Find Features(이름·종류·qualifier 검색)
- **Custom Feature Types** — 표준 INSDC 타입 데이터베이스 + 나만의 타입(색·설명·번역 허용·
  qualifier 선택), CSV로 공유, 모르는 타입이 든 파일을 열면 custom 등록/misc_feature 변경 선택

### Translations (가이드 11)

- 리딩프레임 행 — ORFs only / Top 3 / Bottom 3 / All 6 / 개별 프레임,
  정방향 주황·역방향 초록, 행 끝에 `+1 →` 같은 라벨
- 3-letter 아미노산 코드, Translation Options(최소 ORF 길이·개시코돈·가닥·유전암호)
- Tools → Letter Codes, Genetic Code Tables(64코돈 행렬, 표준과 다른 코돈 음영)
- Add Translated Feature, 겹치는 세그먼트로 리보솜 프레임시프트(−1 slippage) 표현

### Colours (가이드 12, 03)

- 염기 색 스킴 — Nucleotide / GC content / Purine·Pyrimidine / Custom / History
- **Set DNA Colour**(가닥별), **History colours**(스냅샷 비교로 그 단계가 바꾼 염기만)
- 단백질 스킴 — Zappo, Clustal X, DNASTAR 분류, Transmembrane tendency + 색상 키

### Searching (가이드 13)

- `app/search.js` 신규. `(AG)` `R` `!R` `!(AG)` `"R"` 문법 전부 지원
- 4개 모드 — DNA 정확 / 유사(미스매치·갭, 미스매치는 빨간 박스) / 단백질 / 이름(효소·feature·primer)
- 전체 일치는 노란색, 현재 항목은 초록색, 상태줄에 가닥·차이 수·좌표
- 최근 검색 목록, Help → Search Tips

### Focus on Region (가이드 14)

- 좌표 입력·미니맵·파란 박스 드래그로 패닝/크기 조절
- 포커스 중에는 Map(선형으로)·Sequence·Enzymes·Features·Primers 전부 그 구간만

### History (가이드 15)

- Show as Text, File → Export History, Edit → Copy History, 인쇄에 History 포함
- History Actions — Erase Selected Operation Plus Ancestors, Trim History Tree(N단계 유지/전체 삭제)
- Operation 라벨을 클릭하면 그 단계의 입력이 제품에서 차지하는 영역을 색칠
- Trim은 Undo로 되돌아갑니다(스냅샷에 provenance 포함)

### Display (가이드 02, 03, 04)

- Sequence Format — Default / Compact / Plain text, Blocks of 3·10, 고정 줄너비, 시퀀스 미니맵
- Show Locations, **Set Origin**(linear은 음수 좌표 가능 — 프로모터 영역 번호 매기기)
- **Split Window**(좌우/상하, 분할선 드래그, 양쪽 선택 연동)
- **GC content 플롯** — 선형맵 상단·원형맵 링, 중점 기준 빨강/파랑, 윈도우·축·중점 설정
- **Rotate Map** 0–359.9°, Font Size, 툴바 표시/스타일
- 단백질 주석을 **Region / Site / Bond**로 분리, 단백질 미니맵, compact format

### Info Panel (가이드 05)

- natural/synthetic, 실험 확인, laboratory host·transformation strain·methylation,
  source organism·sequence class, description·accession·code number·author·comments
- **References + PubMed PMID로 인용정보 가져오기**(NCBI E-utilities)
- **Embedded Files** — 드래그&드롭으로 아무 파일이나 서열 파일 안에 보관, 이름변경·순서·내보내기
- 전부 GenBank 헤더(DEFINITION, ACCESSION, KEYWORDS, SOURCE, ORGANISM, REFERENCE, COMMENT)로 왕복

### Primers (가이드 09, 10)

- **Hybridization Parameters** — 3′ 최소 일치 염기 수와 최소 Tm. 결합 탐색이 3′ 끝에서
  안쪽으로 가장 긴 일치 구간을 찾는 방식으로 바뀌었고, 어닐링된 염기 수와 Tm을 보고합니다.
- Add Primer 재작성 — 선택영역·가닥 선택, **효소 사이트/코돈/펩타이드 삽입 도구**,
  5′ tail 빨간 표시, 프라이머가 템플릿에 어닐링된 모습 표시
- Edit Primers 일괄(색·인산화·대소문자), 표시/숨김
- Primers 뷰 — 상세 펼침(결합 부위마다 어닐링 길이·Tm), 열 정렬·Sort by(Date added 포함), Actions
- 리스트 붙여넣기 import, Export Primer Data(9개 열), IDT OligoAnalyzer 연결
- 프라이머 구조 **Monomer / Homodimer** 탭 + **온도 파라미터**(구조를 실제 온도에서 채점)
- **Primer Database** — 프로그램 단위 저장소, 컬러 라벨로 세트 구성, CSV·TSV·열린 파일·
  붙여넣기 import(무효·중복 경고), JSON 백업/복원, **Annotate From Primer Database**

### 고친 것

- `file://`에서 stylesheet를 읽을 수 없어 내보낸 지도가 검게 나오던 문제
- async 다이얼로그 핸들러가 창을 닫지 않던 버그 5곳
- 두 효소로 벡터를 자를 때 짧은 쪽 arc를 남기던 문제
- Gibson REV 프라이머 꼬리를 반대쪽 조각에서 떼어 오던 문제
- MSA 원형 회전에서 gap이 더 생기던 문제, RNA folding 속도

## v0.1.0

첫 공개 판. Map/Sequence/Enzymes/Features/Primers/History 뷰, GenBank·FASTA·.dna
읽기, NCBI accession 다운로드, 제한효소·Golden Gate·Gibson·Gateway 클로닝, 코돈 최적화,
.ab1 정렬, 아가로스 젤, 프라이머 설계, 인쇄·PDF.
