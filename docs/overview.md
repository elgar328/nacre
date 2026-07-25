# nacre

순수 Rust CAD 커널. 정확한 b-rep(NURBS/해석적 곡면)을 진실로 보관하면서, append-only 단일 참조 저장과 출처 태그 근사를 결합한 하이브리드 설계. 목표는 정밀 기계 CAD(STEP 입출력 포함), 개인 프로젝트.

**전체 설계는 `docs/design.md`에 있다. 작업 전 반드시 읽을 것.** 이 파일은 그 요약과 불변 규칙이다.

## 이름

nacre = 진주층(자개). 한 겹씩 침착만 되고 이미 만든 겹은 절대 고쳐지지 않으며(append-only), 그 층상 구조 덕분에 원료 광물보다 수천 배 강해진다(robustness). 설계 철학이 이름에 들어 있다.

## 절대 원칙 (설계의 DNA — 위반 금지)

1. **정확 기하가 진실이다.** 평면·원통·구·NURBS는 해석적 형태로 영구 보관. 메시는 파생물(캐시). 절대 "메시를 진실로" 만들지 말 것 — 그것이 Fornjot이 실패하고 우리가 피하려는 바로 그 결정이다.
2. **모든 객체는 append-only 저장소(`Store<T>`)에 한 번만 존재하고 `Handle<T>`로만 참조된다.** 수정·삭제 없음. 동일성은 Handle 비교(정수)이지 좌표 비교(부동소수점)가 아니다. Handle 유효성은 자기 Model 안에서만 성립.
3. **연산 이력을 보존한다.** 모델은 연산 로그의 재생 결과. replay 보장은 "동일 로그·동일 파라미터 → 동일 모델"까지. 파라메트릭 편집(로그 중간 수정)은 v1 비목표.
4. **tolerance는 "발견된" 교차에만 존재한다.** `Origin::Constructed`(구성 시 동일성을 앎 — tolerance 없음) vs `Origin::Discovered { tol }`(교차 계산으로 발견 — 국소 tol 보유)를 타입으로 구분. tol은 상수가 아니라 relaxation의 실측 달성치.
5. **진실/캐시 분리.** 교차곡선=곡면쌍(진실)+스플라인(캐시), 모델=연산로그(진실)+tessellation(캐시), 정점=곡면들(진실)+f64좌표(캐시), 인접=위상store(진실)+Adjacency(캐시). 캐시는 언제든 버리고 재생성 가능.

## 범위 (비목표 — 요청받아도 커널에 넣지 말 것)

순수 기하 커널만: 기하·위상·수치봉합·연산·tessellation·검증까지. **제외:** GD&T/PMI(제작공차), 스타일·색상·레이어·visibility, 제작자·승인·날짜, 제품구조·조립·리비전. 코드-CAD 스크립트가 요구하는 편의(다중 솔리드 값, 다인수 불리언, 패턴·미러, 값 의미론)도 커널 밖 — **편의 레이어**(`nacre-kit`, 워크스페이스 밖 별도 리포)의 몫이다. 이들은 응용 영역 — 커널은 영구 Handle만 제공하고 응용이 사이드카(`HashMap<Handle<_>, 응용데이터>`)로 매단다.

**설탕 vs 커널 판별 기준 (2026-07-26 결정 — 앞으로 모든 편의 기능에 적용).** *편의 레이어는 커널 op을 **조합**만 한다. 새로운 exact 술어 판단(내/외·일치·방향)이 필요해지는 순간 그것은 커널(또는 커널 측) 몫이다.* 근거: exactness 규율(indirect predicates·`Origin`·`validate`)이 전부 커널에 살아서, 앱이 f64로 포함 판정을 시작하면 커널이 없앤 "조용히 틀림"이 한 층 위에서, 소비자마다 다르게 부활한다. 적용 예 — 다인수 fuse/cut/common(fold)·다중 솔리드 값·섬마다 extrude 호출 = **설탕**; 구멍 있는 스케치(`Profile2d` 다중 루프)·링 중첩 판정(`point_in_ring`) = **커널**. nacre-step은 커널 Model을 AP242(Ed2) 형상 서브셋 엔티티로 번역하는 **어댑터** — 직렬화 백엔드는 교체 가능(커널 무지; 개발 중 step-io, 최종 경량 라이터). STEP 가져오기(외부 파일 열기)는 v1 비목표(healing·근사→정확 승격·히스토리 없는 dumb solid를 동반, 앱 레이어 전제).

## 불리언 전략 (핵심)

자체 구현을 커버리지 사다리로 쌓는다. **OCCT에 위임하지 않는다.**
- **M5 `PolyhedralBoolean`**: 평면 솔리드 전용. 평면-평면 교차는 닫힌 형식 직선, 내/외 판정은 exact 술어(orient3d). 여기서 강건 불리언은 연구가 아니라 케이스워크. → "OCCT 없이 직동하는 진짜 커널" 달성 지점.
- **M6 `QuadricBoolean`**: 평면+실린더+구+원뿔. 닫힌 형식 교차(타원·원). 실용 커버리지 대부분.
- **M7 `HybridBoolean`**: 일반 NURBS SSI. 출처태그 메시 → 조합 결정 → 정확 곡면 스냅백. 연구 구간.
- 커버리지 밖 입력은 `BoolError::Unsupported`로 **정직하게 거절**. 조용히 틀리지 말 것.

## OCCT / 외부 코드 규칙

- **OCCT는 오라클(dev 테스트 채점자) 전용.** 제품 경로에 위임·링크 없음. 연동은 out-of-process(STEP 왕복), macOS는 `brew install opencascade` 1순위 / uv+OCP 폴백. 헬퍼는 `tools/occt-helper/`, 프로토콜 `helper <fuse|cut|common> <a> <b> <out>` + exit 0/1/2 + stdout JSON.
- **OCCT 코드는 읽되 베끼지 말 것** — LGPL-2.1이라 번역/차용은 라이선스 오염이고, OCCT 아키텍처가 딸려와 우리 설계와 충돌한다.
- **알고리즘 차용 가능 소스:** Manifold(Apache-2.0), 문헌(Piegl&Tiller *The NURBS Book*, Hoffmann).

## 정밀도 인프라 (실수 최소화 — 코드보다 먼저)

- **판정(부호)** = exact 적응 술어(`geometry-predicates` 크레이트). **반복 구성(좌표)** = double-double(`twofloat`/`qd`). 임의정밀 유리수는 반복에 쓰지 말 것(비트 폭발).
- relaxation은 3단 사다리: f64 Newton → 정체 시 double-double → 실패 시 `RelaxError::Tangential`로 상신(삼키지 말 것).
- 모든 연산 직후 `validate` 실행(디버그). 오일러-푸앵카레는 내부 루프 항 포함형 `V−E+F = 2(S−G) + L_i`.
- proptest로 불변식(멱등성·부피보존·불리언 대수)과 OCCT 오라클 diff.

## 마일스톤 (docs/design.md §8 참조)

M1 뼈대(Store/Handle, Plane/Line, 정육면체, validate, OBJ덤프·STEP출력→기존 뷰어에 위임) → M2 스케치·extrude·replay·STEP → M3 곡선기하(NURBS·Arc·Cylinder, tess 출처태그) → M4 면 위 작업(imprint·pad, 여기까지 전부 Constructed, 오라클 가동) → M5~M7 불리언 사다리(이 즈음 워크스페이스 밖 인터랙티브 디버그 뷰어 — §design 1).

**첫 작업: M1의 `Store<T>`/`Handle<T>`를 proptest와 함께.** 인프라(검증·테스트)가 코드보다 먼저라는 게 이 프로젝트의 방법론이다.

## 작업 스타일

- 큰 덩어리 말고 작은 단위로. 각 단위는 테스트를 남기고 끝낸다.
- 워크스페이스 구조. 크레이트는 마일스톤 따라 `cargo new`로 추가(지금 9개 다 만들지 말 것).
- 커밋은 작게, 의미 단위로.
- 기하 코드 버그는 눈으로 잡는다 — 애매하면 OBJ 덤프해서 확인.
- **셀별 진행 기록은 `docs/dev-log.md`에 append한다** — 설계 규칙·불변은 `design.md`, 오른 사다리·반증된 예측 등 셀 단위 로그는 dev-log로 분리(design.md 재비대 방지).
