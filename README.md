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
| Map | 원형/선형 지도, Feature 화살표, 제한효소 절단 위치 라벨, 프라이머, ORF, 지도에서 드래그로 영역 선택, 확대/축소 |
| Sequence | 이중가닥 표시, 10 bp 눈금, Feature 막대, 아미노산 번역, 효소 위치, 직접 타이핑으로 편집, 복사/붙여넣기, DNA·단백질 검색(IUPAC 지원), 위치로 이동 |
| Enzymes | 효소 77종, 절단 횟수·위치·말단 형태(5′/3′/blunt), 원하는 효소를 골라 지도에 표시 |
| Features | 추가·편집·삭제, 종류·색·방향 지정, 흔한 서열 자동 인식(T7, lac operator, loxP, His tag 등) |
| Primers | 선택 영역으로 프라이머 만들기(위/아래 가닥), Tm(nearest-neighbor), GC%, 결합 위치(5′ 꼬리가 있으면 3′ 15 nt로 찾음) |
| Actions | 제한효소 절단 + 가상 아가로스 겔, 번역(6 frame), 역상보, 원점 변경, 원형/선형 전환 |
| History | 모든 변경 기록, 원하는 단계로 되돌리기, Undo/Redo |

## 직접 빌드하기

```bash
npm install
npm start            # 개발 모드로 실행
npm run dist:win     # Windows에서 설치 파일과 portable exe 만들기
```

## 구조

```
main.js         Electron 메인 프로세스 (창, 파일 열기/저장 대화상자)
preload.js      화면에서 쓰는 window.gcHost API
app/index.html  화면 틀
app/style.css   스타일
app/dna.js      DNA 로직 (효소, 번역, Tm, ORF, 파일 형식)
app/app.js      화면과 상호작용
```

예제 플라스미드(pGC-Demo)는 기능을 보여주려고 만든 합성 서열이며, 실제 벡터 서열이 아닙니다.
