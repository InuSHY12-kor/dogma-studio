# 변경 내역

## v0.3.0

the reference tool 화면을 기준으로 **지도·클로닝·시뮬레이션 화면을 다시 만든 판**입니다.
요청하신 16가지 항목을 모두 반영했고, 각 항목은 headless Chromium과 Electron에서 동작을 확인했습니다.

### 새로 생긴 것

- **Agarose 시뮬레이션 탭** — 툴바의 `Digest`가 **Agarose**로 바뀌었고, 누르면 파일을 연 것처럼
  **새 탭**으로 열립니다. lane은 최대 12개, 각 lane마다 서열 하나와 **제한효소 digest 또는 PCR**을
  고릅니다. 창 안의 지도에서 **효소를 클릭하면 그 효소로 자르고, Shift-클릭으로 두 번째 효소를
  더해 double digest**가 됩니다. lane을 `amplify with`로 바꾸면 **프라이머를 클릭하고 마주보는
  프라이머를 Shift-클릭**해서 증폭 산물을 올립니다(확인: M13 fwd + M13 rev → 123..225, 103 bp).
  lane이 바뀔 때마다 젤이 다시 그려지고, Enzymes·Primers 표에서도 같은 lane을 조작할 수 있습니다.
  `Open as an Agarose document`로 지금 구성을 기존 아가로스 젤 문서로 넘길 수 있습니다.
- **클로닝 창 하나로 통합** — Cloning 메뉴의 모든 항목이 **같은 Vector · Fragment · Product 창**을
  엽니다. 창 **오른쪽 위 드롭다운**에서 방법을 고르면 적용되는 알고리즘이 바뀝니다:
  제한효소 클로닝 · Golden Gate · Gibson · In-Fusion · Gateway BP/LR. 방법을 바꿔도 고른 조각은
  그대로 두고 접합부만 다시 판정하므로, 같은 설계를 다른 방법으로 바로 시험해 볼 수 있습니다.
  - 엔진에 두 가지 준비 방식이 추가됐습니다: **Type IIS 사이트 사이를 잘라낸 part**와
    **att core 사이 구간**(att 전용 접합 판정 포함).
  - Golden Gate는 **overhang이 겹치거나 part 안에 사이트가 남으면** 알려주고, Gateway는
    **부산물 크기**를 함께 적습니다.
  - 제한효소로 자를 때 **효소 개수 제한이 없습니다**. 한 튜브에 넣은 것처럼 모든 절단을 함께
    계산하고, 이웃한 두 절단 사이에서 가장 긴 조각(또는 `Keep the stretch between the cuts
    instead`를 켜면 가장 짧은 조각)을 씁니다.
- **지도 옆 사이드 메뉴를 클로닝 창과 시뮬레이션 창에도** — 본 창 Map 탭과 **같은 순서·같은
  아이콘·같은 툴팁 위치**(커서에 가리지 않게 왼쪽)로 붙었습니다. 각 패널에서 **지금 지도에
  그려지는 것을 직접 고릅니다**: 효소 세트 + 자르는 효소 전체 체크 목록(세트에 없는 효소를
  체크하면 chosen enzyme으로 유지), feature·primer 체크 목록과 Show all / Hide all, ORF와 최소
  길이, feature 색, 눈금자. 조각을 실제로 자르는 효소는 세트를 아무리 좁혀도 지도에 남습니다.
  설정은 **탭별·서열별로 따로** 기억합니다.

### 바뀐 것

- **선형 지도가 창 가운데에** 놓입니다(확인: 846 px 창에서 위아래 177 px씩). 그림이 창보다 커지면
  위쪽 끝부터 스크롤됩니다.
- **효소·feature 이름이 다른 라벨의 지시선에 가려지지 않습니다.** 지시선을 모두 먼저 그리고 이름을
  나중에 그립니다(원형·선형 모두).
- **왼쪽 패널의 기능이 전부 오른쪽 사이드 메뉴로** 들어갔습니다. 아이콘은 켜고 끄기, 옆의 작은
  **펼침 단추**는 그 그룹의 나머지 설정을 엽니다. 순서는 **효소 · feature · primer · ORF · 색 ·
  줄 번호 · alignment**이고, 글리프는 the reference tool 사이드 툴바를 따랐습니다. 왼쪽 패널은 선택·서열
  수치만 남기고 기본 꺼짐(View → Show Details Panel로 다시 켬).
- **아미노산으로 DNA를 고를 수 있습니다.** 아미노산을 클릭하면 그 코돈 3염기가 선택되고, 드래그
  하면 코돈 단위로 늘어나며, 더블클릭하면 번역된 feature 전체가 선택됩니다. 역방향 CDS에서도
  됩니다. **Ctrl+Shift+C = 역상보 복사, Ctrl+Alt+C = 번역 복사**(Edit 메뉴와 단축키 목록에도 있음).
- **Feature 창을 세그먼트 중심으로 다시 만들었습니다.** 방향은 4개 단추(무방향/정방향/역방향/
  양방향), 위쪽에 세그먼트 색 띠와 전체 길이, 아래에 이름·위치·크기·색 표와 Split / Merge /
  Delete Segment. 표에 적은 내용이 그대로 저장됩니다.
- **Feature 데이터베이스와 제한효소 데이터베이스가 각각 별도의 창**으로 열립니다. 창을 띄워 둔
  채 본 창에서 계속 작업할 수 있습니다. 효소 창은 인식서열 복사 단추, 양쪽 가닥 절단 그림(5′/3′/
  blunt와 Type IIS 모두 올바르게), overhang·palindrome·인식 길이, 제조사별 buffer 활성표, 같은
  인식서열을 갖는 효소와 **붙일 수 있는 말단을 만드는 효소** 목록을 보여줍니다.
- **DNA Calculations 창**을 요청하신 모양으로 다시 만들었습니다: 전체와 선택 구간을 번호 붙인
  블록으로 나누고, MW(지수 표기), 염기 조성, GC, µg↔pmol, 직접 입력하는 양 입력칸 두 개(분자 수와
  nM 표시), 맨 아래 A260 환산.
- **유전암호 표**를 표준 배치(1번째 염기 왼쪽, 2번째 위, 3번째 오른쪽)로 고쳤습니다. 2번째 염기 C
  이후가 어긋나던 문제가 사라졌고, 칸마다 **코돈 · 세 글자 이름 · (한 글자) · 사용빈도 숫자**가
  한 줄로 들어갑니다(가장 흔한 코돈이 1.00).
- **Alignment 창이 빨리 뜹니다.** 정렬을 seed 대각선 주변 띠에서만 계산하고 후보를 강한 순으로
  시도합니다. 900염기 read가 62 ms → 13 ms(동일한 identity·strand·시작 위치·차이 목록).
- 철자: 인터페이스와 문서의 `optimisation` → `optimization`.

### 고친 것

- 지도에서 **효소 이름이 서로 겹쳐 찍히던 문제**(이름 길이를 위치 표기까지 포함해 계산하고,
  줄이 모자라면 그림이 아래로 늘어납니다). 오른쪽 끝 이름은 절단선 왼쪽으로 붙습니다.
- 공용 지도의 **선택 영역이 검게 칠해지던 문제**(fill이 지정돼 있지 않았습니다).
- 클로닝 창의 **`Cut ---- with ----` 첫 칸이 아무 일도 하지 않던 문제** — 첫 칸이 효소 자체가 되고,
  하나 고르면 다음 줄이 열립니다.
- 시뮬레이션 문서에서 선택 정보가 없어 화면이 그려지지 않던 문제, 탭을 닫은 뒤 보기 종류가
  남아 있던 문제.

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
