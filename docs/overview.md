# nacre

순수 Rust CAD 커널. 정확한 b-rep(NURBS/해석적 곡면)을 진실로 보관하면서, append-only 단일 참조 저장과 출처 태그 근사를 결합한 하이브리드 설계. 목표는 정밀 기계 CAD(STEP 출력 포함), 개인 프로젝트.

이 파일은 불변 규칙이다. 설계의 「왜 이 모양인가」와 타입은 `docs/design.md`, 지금 하는 일과 남은 일은 `docs/todo.md`. **작업 전 셋을 읽는다.**

## 이름

nacre = 진주층(자개). 한 겹씩 침착만 되고 이미 만든 겹은 절대 고쳐지지 않으며(append-only), 그 층상 구조 덕분에 원료 광물보다 수천 배 강해진다(robustness).

## 절대 원칙 (설계의 DNA — 위반 금지)

1. **정확 기하가 진실이다.** 평면·원통·구·NURBS는 해석적 형태로 영구 보관한다. 메시는 파생물(캐시)이다. 메시를 진실로 만들지 않는다.
2. **모든 객체는 append-only 저장소(`Store<T>`)에 한 번만 존재하고 `Handle<T>`로만 참조된다.** 수정·삭제 없음. 동일성은 Handle 비교(정수)이지 좌표 비교(부동소수점)가 아니다. Handle 유효성은 자기 Model 안에서만 성립한다.
3. **연산 이력을 보존한다.** 모델은 연산 로그의 재생 결과다. replay 보장은 「동일 로그·동일 파라미터 → 동일 모델」까지. 파라메트릭 편집(로그 중간 수정)은 v1 비목표.
4. **tolerance는 「발견된」 교차에만 존재한다.** 구성 시 동일성을 아는 점에는 tolerance 개념 자체가 없다. 이 구분은 정점의 실현 캐시가 든다: 좌표에는 증명된 축별 경계가 있거나(`Bounded` — 정의에서 실현됨), 싼 도로가 멈춰 아직 안 된 것이거나(`Ceiling` — 더 물으면 답이 있다), 실현이 뒤에 없다(`Unrealized`). 경계는 좌표가 참값에서 얼마나인지이지 tolerance가 아니다. 발견된 점의 국소 tol은 배열이 측정해 쓰지만(`SeamVertex.tol`) 캐시에 저장하지 않는다 — 잔차는 담체까지의 거리일 뿐 좌표가 참값에서 얼마나인지를 말하지 않는다.
5. **진실/캐시 분리.** 교차곡선 = 곡면쌍(진실) + 스플라인(캐시), 모델 = 연산로그(진실) + tessellation(캐시), 정점 = 곡면들(진실) + f64 좌표(캐시), 인접 = 위상 store(진실) + Adjacency(캐시). 캐시는 언제든 버리고 재생성할 수 있다.

## 범위 (비목표 — 요청받아도 커널에 넣지 않는다)

순수 기하 커널만: 기하·위상·수치봉합·연산·tessellation·검증까지. **제외:** GD&T/PMI(제작공차), 스타일·색상·레이어·visibility, 제작자·승인·날짜, 제품구조·조립·리비전. 코드-CAD 스크립트가 요구하는 편의(다중 솔리드 값, 다인수 불리언, 패턴·미러, 값 의미론)도 커널 밖 — 편의 레이어(`nacre-kit`, 워크스페이스 밖 별도 리포)의 몫이다. 커널은 영구 Handle만 제공하고 응용이 사이드카(`HashMap<Handle<_>, 응용데이터>`)로 매단다.

**설탕 vs 커널 판별 기준.** 편의 레이어는 커널 op을 **조합**만 한다. 새로운 exact 술어 판단(내/외·일치·방향)이 필요해지는 순간 그것은 커널 몫이다. 근거: exactness 규율(간접 술어·정점의 정의·`validate`)이 전부 커널에 살아서, 앱이 f64로 포함 판정을 시작하면 커널이 없앤 「조용히 틀림」이 한 층 위에서 소비자마다 다르게 부활한다. 예 — 다인수 fuse/cut/common(fold)·다중 솔리드 값·섬마다 extrude 호출 = 설탕; 구멍 있는 스케치(`Profile2d` 다중 루프)·링 중첩 판정 = 커널.

`nacre-step`은 커널 Model을 AP242(Ed2) 형상 서브셋 엔티티로 번역하는 어댑터다(직렬화 백엔드는 교체 가능). STEP 가져오기는 v1 비목표 — healing·근사→정확 승격·히스토리 없는 dumb solid를 동반하므로 앱 레이어 전제다.

## 불리언 전략

자체 구현을 커버리지 사다리로 쌓는다. **OCCT에 위임하지 않는다.**

- **M5 평면 솔리드.** 평면-평면 교차는 닫힌 형식 직선, 내/외 판정은 exact 술어. 엔진은 **평면 클래스당 셀 복합체 하나**(`nacre-ops::arrangement`, 문헌의 winding 단일 배열)다. 실패하는 배치마다 경로를 덧대지 않는다 — 새 결함은 그 엔진을 일반화해서 닫는다.
- **M6 평면 + 이차곡면.** 닫힌 형식 교차(타원·원). 실용 커버리지 대부분.
- **M7 일반 NURBS SSI.** 출처태그 메시 → 조합 결정 → 정확 곡면 스냅백. 연구 구간.
- 커버리지 밖 입력은 `BoolError::Rejected`로 **이름 붙여 거절**한다. 조용히 틀리지 않는다.

## OCCT / 외부 코드 규칙

- **OCCT는 오라클(dev 테스트 채점자) 전용.** 제품 경로에 위임·링크 없음. 연동은 out-of-process(STEP 왕복), macOS는 `brew install opencascade` 1순위 / uv+OCP 폴백. 헬퍼는 `tools/occt-helper/`, 프로토콜 `helper <fuse|cut|common> <a> <b> <out>` + exit 0/1/2 + stdout JSON.
- **OCCT 코드는 읽되 베끼지 않는다** — LGPL-2.1이라 번역·차용은 라이선스 오염이고, OCCT 아키텍처가 딸려와 우리 설계와 충돌한다.
- **알고리즘 차용 가능 소스:** Manifold(Apache-2.0), 문헌(Piegl & Tiller *The NURBS Book*, Hoffmann).
- **제품 경로의 의존성은 순수 Rust다 — 성능을 이유로도 예외를 두지 않는다.** 목표 환경이 브라우저(`wasm32-unknown-unknown`)이고, C 의존은 거기서 빌드되지 않는다.

## 정밀도 인프라

- **판정(부호)은 좌표가 아니라 정의로 한다.** 전부 유리수면 exact 적응 술어(`nacre-predicates`)로 끝난다. 모션이 좌표를 무리수로 만들면 **CIP**(`nacre-judge`)가 받는다 — 점을 `정의 + 방향별 tol`로 들고, sound 필터로 거르고, 못 거르면 정밀도를 올린다. 틀린 부호는 결코 내지 않는다.
- **고정밀 층은 `astro-float`(순수 Rust 임의정밀).** 임의정밀 유리수(`Rat`)는 **정의**에 쓰고 반복 계산에는 쓰지 않는다(비트 폭발).
- **정밀도는 상수가 아니라 모델이 정한다.** 필요한 비트를 조건수에서 계산해 **한 번에 점프**한다(배증 아님). 상한 `JUDGE_PREC_CAP`은 정확성이 아니라 비용 한계다.
- **증명하지 못한 판정은 조용한 0이 아니다.** 일치 정밀도보다 가깝다고 *증명되면* 일치로 처리하고 근거를 보고하며(`boolean_with_report`), 못 좁히면 원인에 이름을 붙여 거절한다(`JudgeExhausted`·`DegenerateWitness`·`PrecisionBudget`). 사용자에게 되묻지 않는다.
- **정확한 점은 정의에서 실현하고, 호출자가 고른 정밀도로 올린다**(`nacre_ops::realize_vertex`). f64 좌표를 만드는 통로는 실현 하나다. 실패를 삼키지 않는다.
- `validate` 는 순수 함수이고 테스트·proptest·census·오라클이 연산마다 부른다(연산 자신은 부르지 않는다 — `nacre-ops` 는 `nacre-validate` 에 dev-의존만 한다). 오일러-푸앵카레는 내부 루프 항 포함형 `V−E+F = 2(S−G) + L_i`.
- proptest로 불변식(멱등성·부피보존·불리언 대수)과 OCCT 오라클 diff.

## 마일스톤

M1 뼈대(Store/Handle, 평면·직선, 정육면체, validate, OBJ·STEP 출력) → M2 스케치·extrude·replay·STEP → M3 곡선기하(NURBS·Arc·Cylinder, tess 출처태그) → M4 면 위 작업(pad·pocket, 오라클 가동) → M5~M7 불리언 사다리. 상세는 `design.md` 「마일스톤 사다리」, **현재 위치는 `todo.md` 의 맨 위 항목**이다.

## 작업 스타일

- 큰 덩어리 말고 작은 단위로. 각 단위는 테스트를 남기고 끝낸다. **인프라(검증·테스트)가 코드보다 먼저다.**
- 커밋은 작게, 의미 단위로. **커밋 본문이 그 단위의 기록이다** — 무엇을 재고 무엇이 반증됐는지는 거기 적는다.
- 크레이트는 필요해질 때 추가한다.
- **모듈의 자리.** 여러 단계·관심사를 담은 모듈은 폴더다(`src/arrangement/` — 줄 수는 신호지 규칙이 아니다). `mod.rs` 는 그 모듈의 **어휘와 선언**이다: 둘 이상의 단계가 주고받는 타입, 매크로, `mod` 줄, 재수출. 단계는 파일 하나씩. 타입의 필드는 올리지 않는다 — 다른 파일이 필드를 읽어야 하면 그 타입은 단계 사이의 기록이므로 `mod.rs` 로 간다(private 은 정의 모듈과 그 **후손**에서 보인다). 자기가 한 단계인 타입은 필드를 쥐고 다른 단계가 부르는 메서드만 `pub(super)` 로 연다. 함수의 `pub(super)` 는 「이 단계의 입구」라는 뜻이고, 두 파일이 서로를 부르면 그 둘은 한 단계다. 자식 파일은 `use super::*;` 로 시작한다 — 그래서 `mod.rs` 안의 **정의**는 glob 이 들여온 같은 이름을 가리지만 자식에서 재수출한 이름은 못 가린다(`arrangement` 의 `boolean`·`Wall` 이 `mod.rs` 에 있는 이유). **계측의 자리**: 엔진의 어휘를 이름하는 계측은 그 모듈의 것이고 파일만 `src/tests/probes/` 로 내보낸다(`#[path]` 로 마운트해 모듈은 제자리에 둔다 — 최상위로 올리면 거는 쪽과 이름하는 쪽이 양방향 간선을 만들어 순환이 는다). 아무것도 이름하지 않는 계측만 최상위로 간다(`phase`).
- **파일 경로에 기대는 테스트가 다섯 있다.** 소스를 경로로 읽는 넷 — `par.rs` 의 병렬 스위치 스캔(`src` 전체를 걷는다), `tests/probes/rotation_sweep.rs` 의 두 스캔(모듈을 파일 또는 폴더로 받는다), `tests/instruments/module_graph.rs`(`src` 의 모듈 파일을 읽어 간선을 센다), `tests/instruments/silencers.rs`(`src` 전체를 걷는다) — 과, 가드가 울린 **파일**을 경로로 고정하는 `tests/reject_census.rs`. 파일을 옮기면 스캔은 걷는 범위 밖으로 나간 파일을 조용히 놓칠 수 있고 `reject_census` 는 빨개진다: 옮긴 뒤 위반을 심어 빨개지는지 본다.
- **테스트의 자리.** 단위 테스트는 제품 파일 안이 아니라 `src/tests/<모듈>.rs` 에 산다(제품 모듈의 자식으로 남도록 `#[path]` 로 선언한다 — private 항목을 본다). 2,000줄을 넘으면 `src/tests/<모듈>/` 폴더로 주제별로 나눈다. 200줄 미만이고 proptest 가 없으면 인라인으로 둬도 된다. 테스트를 한 단 깊은 모듈로 옮기면 `super::` 의 뜻이 바뀐다 — 옮긴 뒤 `super::`·`self::` 경로를 훑는다.
- **proptest 시드는 테스트 파일 옆에 자동으로 기록되고 커밋한다**(`<파일>.proptest-regressions`). 단위 테스트의 `proptest!` 블록은 `WithSource` 설정을 달고 `src/tests/` 아래에만 산다 — 기본 설정은 시드를 크레이트 루트의 별도 트리에 써서, 파일을 옮기면 조용히 끊긴다.
- **픽스처는 한 벌** — `crates/nacre-ops/tests/support/`(`fixtures.rs`·`stated.rs`)를 단위·통합 테스트가 함께 읽는다. 본문이 같은 헬퍼만 합친다(같은 이름 ≠ 같은 일). `census.rs`·`reject_census.rs` 는 코퍼스를 동결하려고 자기 픽스처를 든다.
- **프로세스 전역 상태를 읽는 테스트는 자기 바이너리를 갖는다**(`census`·`reject_census`·`wide_datum_cost`) — 같은 프로세스의 어떤 불리언이든 그 카운터를 움직인다. 같은 바이너리 안에서 전역 계측을 읽어야 하면 계측이 내주는 독점 세션을 잡고, 개수가 아니라 존재·전칭만 단언한다.
- 기하 코드 버그는 눈으로 잡는다 — 애매하면 예제로 OBJ를 덤프해서 확인한다(출력은 `target/` 아래).

### 문서 규칙

| 내용 | 집 |
|---|---|
| 불변 규칙·범위·관문 | `docs/overview.md` |
| 설계 — 왜 이 모양인가, 타입, 「가지 말 것」 표 | `docs/design.md` |
| 지금 하는 일·다음·알려진 결함·보류 | `docs/todo.md` |

1. **현재형만 적는다.** 「지금 무엇이 참인가」만. 고친 흔적·날짜·작업 번호·경위를 덧붙이지 않는다 — 바뀌면 그 문장을 **고쳐 쓴다**. 결정을 받치는 측정치는 경위가 아니므로 남긴다.
2. **역사는 git이 든다.** 경위는 커밋 본문에. 재 보고 버린 접근만 `design.md` 의 「가지 말 것」 표에 한 행으로 남긴다(무엇을 · 왜 · 커밋).
3. **`todo.md` 의 끝난 항목은 체크하지 않고 지운다.** 항목은 번호가 아니라 이름을 갖는다.
4. 같은 사실을 두 곳에 적지 않는다. 다른 문서의 내용은 절 **이름**으로 가리킨다.
5. 규칙이 `design.md` 에 살려면 코드에 집이 있어야 한다 — 코드에 없는 이름을 현재형으로 가르치지 않는다. 타입을 바꾸는 변경은 `design.md` 의 그림을 같은 커밋에서 고친다.
6. 소스 주석도 같다: 코드가 **왜 이 모양인지**를 현재형으로. 문서의 절 번호·작업 번호·날짜를 인용하지 않는다.

### 관문 (한 단위를 닫기 전에 도는 것) — 여기가 유일한 원본이다

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --no-fail-fast
cargo test -p nacre-ops --no-default-features
cargo test -p nacre-ops --test census -- --ignored --nocapture | grep '^c '          # 두 프로파일 diff
cargo test -p nacre-ops --release --test census -- --ignored --nocapture | grep '^c '
cargo test -p nacre-ops --test reject_census
cargo test --workspace --no-fail-fast -- --ignored --skip measure_
cargo test -p nacre-ops --release --test perf -- --ignored --nocapture   # 성능은 release로 따로
cargo doc --workspace --no-deps                                          # intra-doc 링크
cargo build --workspace                                                  # 테스트 기능을 끈 제품 빌드 — 경고 0
cargo build --release -p nacre-ops                                       # debug_assertions 를 «끈» 빌드 — 경고 0
```

- **커밋 훅이 못 보는 것 셋.** 훅은 `clippy --all-targets` 와 `test` 를 돈다 — dev-의존이 `test-util` 을 켠 채로고, `cfg(test)` 와 `debug_assertions` 가 **둘 다 켜진** 빌드다. 그래서 (1) 어떤 항목을 `test-util` 뒤로 보내 제품 빌드에서 그 항목만 쓰던 헬퍼가 죽은 코드가 돼도 훅은 초록이고(맨 `cargo build --workspace` 가 경고한다), (2) 그 항목을 가리키던 intra-doc 링크가 깨져도 훅은 모르며(`cargo doc`), (3) **`debug_assertions` 가 꺼져야 드러나는 것**은 훅도 `cargo build --workspace`(디버그)도 못 본다 — `cfg(debug_assertions)` 항목만 담던 glob 재수출이 릴리스에서 비거나, 그런 항목만 쓰던 `use` 가 릴리스에서 죽은 코드가 되는 부류다. **`cargo build --release` 가 그 눈이고, 훅 밖에 있다.** 가시성·`cfg`·게이트를 건드린 변경은 셋 다 돌린다.
- **`cargo doc`**: 이 저장소는 intra-doc 링크를 2,000곳 넘게 쓰고 `fmt`·`clippy`·`test` 중 무엇도 그것을 해석하지 않는다. 경고의 기준선은 **72**(전부 「공개 문서가 비공개 항목을 링크」 부류)이고 규칙은 「기준선보다 늘지 않는다」, unresolved link 는 **0**. 이 계기는 이름이 사라져 깨지는 링크만 잡는다 — 이름이 살아 있는데 가리키는 대상이 바뀐 링크는 개명한 사람이 손으로 훑는다.
- **모듈을 옮겼으면 그래프를 본다**: `cargo test -p nacre-ops --test instruments module_graph`. 두 겹이다 — 불리언 파이프라인이 한 방향으로 흐르는지(`draft` 는 바닥, `assembly` 는 엔진을 안 부른다, 엔진은 앞문을 안 부른다)와 엔진 **밑**의 넷(`bands`·`combinatorics`·`nesting`·`planes`)이 엔진을 안 부르는지와 클래스 표가 그 위의 이름을 안 부르는지(`planes` → `combinatorics`)를 **단언**하고, 나머지 간선은 **얼린 표**로 든다. 표에 없는 간선이 생기면 답은 둘뿐이다: 없애거나, 표에 넣고 **커밋 본문에 왜 그 방향이 맞는지 적는다**. 그래프는 DAG 가 아니다(순환 1쌍 — `construct ⇄ ops`, 같은 층이라 결함이 아니다 — 이 남아 있고 표가 그것을 든다) — 이 관문은 «더 나빠지지 않는다»를 지킨다. 오늘의 표는 `measure_module_graph` 가 찍는다.
- **이름을 개명·은퇴시켰으면 문서도 훑는다**: `python3 tools/deadname-sweep.py` (기본 인자 = 세 문서). 이 계기는 한 물음만 답한다 — 「`crates/` 비주석 사용이 0인가」. 출력은 후보지 작업 목록이 아니다(수식 기호·외부 도구가 섞인다). 값(개수·변종 수)이 바뀐 변경은 이름이 아니라 **옛 숫자**로 문서를 훑는다.
- **이름이 `measure_` 로 시작하는 테스트는 스윕에서 빠진다.** 단언 없이 표를 찍는 계측·스파이크다(스윕에서는 출력이 캡처돼 시간만 쓰고, 전역 타이머를 읽는 것은 병렬이면 숫자가 틀린다). 따로 돌린다: `cargo test -p nacre-ops <이름> -- --ignored --nocapture --test-threads=1`. `#[ignore]` 가 붙었어도 **단언이 있는 테스트에는 이 접두사를 주지 않는다** — 스윕이 그 단언이 도는 유일한 자리다.
- **`perf.rs` 는 스윕에서 빠지고 release 로 따로 돈다** — debug 스윕 시간의 93%가 그 파일의 테스트 둘이고 debug 시간 숫자는 의미가 없다. 다만 그 둘은 큰 회전 fold 를 수십 번 쌓는 커버리지이기도 하므로(release 에서는 `debug_assert!` 가 꺼진다) 가끔 debug 로도 돌린다: `cargo test -p nacre-ops --test perf -- --ignored` (28분).
- **dev·test 프로파일은 `opt-level = 2`, `debug = 1`** 이다. 안쪽 루프가 정확 산술이라 `opt-level = 0` 에서 9배 느리다(스위트 475s → 81s). 증분 컴파일은 느려지지 않고, census 는 비트 동일하며, 패닉의 `file:line` 은 그대로다.
- **카고는 한 번에 하나만** — 겹치면 서로 6배 느려진다. 커밋 훅이 스위트를 돌므로 커밋 중에는 다른 cargo 를 시작하지 않는다. `target/` 은 한 세션에 ~9GB 자라므로 주기적으로 `cargo clean` 한다.
- **성능 비교는 같은 세션 A/B 로만 한다**(`git stash`). 저장된 기준선은 만료된 것으로 본다 — 한 줄도 안 건드린 단계가 +9% 움직인다. 구조적 변화는 trip 수가, 상수 인자는 시간이 말한다.
- **관문 명령은 파이프로 감싸지 말고 파일로 받아 cargo 자신의 `$?` 를 읽는다** — `| tail` 의 exit code 는 `tail` 의 것이다.
- 형제 리포: kit(`../nacre-kit`: fmt·clippy·test), 앱(`../nacre-playground/web`: `npm run wasm:all` + `npx tsc --noEmit` + `npx vitest run`, 그리고 `wasm/` 에 clippy). 커널만 바뀌어도 wasm 재빌드가 필요하다. 거절 이유의 이름을 개명·분할한 커밋은 kit 의 `reason_sentence` 표도 함께 고친다(모르는 이름은 틀리지 않고 덜 풍부해질 뿐이지만, 조용히 그렇게 된다).
