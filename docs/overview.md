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
- **M5 `PolyhedralBoolean`**: 평면 솔리드 전용. 평면-평면 교차는 닫힌 형식 직선, 내/외 판정은 exact 술어(orient3d). → "OCCT 없이 직동하는 진짜 커널" 달성 지점. **★ 여기를 케이스워크로 쌓지 말 것** — 실패하는 배치마다 경로를 덧대는 방식은 2026-07에 폐기했고, 지금은 **평면당 셀 복합체 하나**(`arrangement::boolean`, 문헌의 winding 단일 배열)가 프로덕션이다. 새 결함은 그 엔진을 일반화해서 닫는다.
- **M6 `QuadricBoolean`**: 평면+실린더+구+원뿔. 닫힌 형식 교차(타원·원). 실용 커버리지 대부분.
- **M7 `HybridBoolean`**: 일반 NURBS SSI. 출처태그 메시 → 조합 결정 → 정확 곡면 스냅백. 연구 구간.
- 커버리지 밖 입력은 `BoolError::Rejected`로 **정직하게 거절**. 조용히 틀리지 말 것.

## OCCT / 외부 코드 규칙

- **OCCT는 오라클(dev 테스트 채점자) 전용.** 제품 경로에 위임·링크 없음. 연동은 out-of-process(STEP 왕복), macOS는 `brew install opencascade` 1순위 / uv+OCP 폴백. 헬퍼는 `tools/occt-helper/`, 프로토콜 `helper <fuse|cut|common> <a> <b> <out>` + exit 0/1/2 + stdout JSON.
- **OCCT 코드는 읽되 베끼지 말 것** — LGPL-2.1이라 번역/차용은 라이선스 오염이고, OCCT 아키텍처가 딸려와 우리 설계와 충돌한다.
- **알고리즘 차용 가능 소스:** Manifold(Apache-2.0), 문헌(Piegl&Tiller *The NURBS Book*, Hoffmann).
- **★ 제품 경로의 의존성은 순수 Rust로 간다 — 성능을 이유로도 예외를 두지 않는다.** 고정밀 층에 `astro-float`를 고른 기준(순수 Rust 임의정밀)과 같은 기준이다. **실제로 한 번 시험대에 올랐다(2026-07-29):** 회전 부울 시간의 36.5%가 메모리 할당이고, C로 된 `mimalloc`으로 전역 할당자를 바꾸면 **코드 한 줄에 1.21×**가 나왔다. **채택하지 않는다** — 전역 할당자는 프로그램에서 가장 안전이 중요한 부품이라 남의 unsafe C를 신뢰 기반 한복판에 놓는 일이고, C 컴파일러 의존이 빌드·크로스컴파일에 붙으며, 무엇보다 **`wasm32-unknown-unknown`에서 빌드조차 안 된다**(실측: `wchar.h` 없음) — 목표 환경이 브라우저인데 거기서 못 쓴다. 브라우저용 **순수 Rust** 할당자는 웹앱 착수 때 따로 고르고 재는 별개 항목이다.

## 정밀도 인프라 (실수 최소화 — 코드보다 먼저)

- **판정(부호)는 좌표가 아니라 정의로 한다.** 전부 유리수면 exact 적응 술어(`nacre-predicates`, `geometry-predicates` 기반)로 끝난다. 모션이 좌표를 무리수로 만들면 **CIP**(`nacre-cip`)가 받는다 — 점을 `정의 + 방향별 tol`로 들고, sound 필터로 거르고, 못 거르면 정밀도를 올린다. 틀린 부호는 결코 내지 않는다(no-silent-wrong). 자세히는 design.md §CIP.
- **고정밀 층은 `astro-float`(순수 Rust 임의정밀).** ~~double-double(`twofloat`/`qd`)~~ 은 **실측으로 탈락**했다(H1.5: 영점 근처 cos 오차 ~1.8e-16 — 부호 판정이 일어나는 곳이 곧 영점 근처라 치명적). 임의정밀 유리수(`Rat`)는 **정의**에 쓰고 반복 계산에는 쓰지 말 것(비트 폭발).
- **정밀도는 상수가 아니라 모델이 정한다.** 필요한 비트를 계산해 **한 번에 점프**한다(배증 아님). 상한 `JUDGE_PREC_CAP`은 정확성이 아니라 비용 한계다.
- **증명하지 못한 판정은 조용한 0이 아니다.** 일치 정밀도보다 가깝다고 *증명되면* 일치로 처리하고 **근거를 보고**하며(`boolean_with_report`), 못 좁히면 원인에 이름을 붙여 거절한다(`JudgeExhausted`·`DegenerateWitness`·`PrecisionBudget`). 사용자에게 되묻지 않는다 — 그 정책은 실측이 반박했다(§9 ③).
- 모든 연산 직후 `validate` 실행(디버그). 오일러-푸앵카레는 내부 루프 항 포함형 `V−E+F = 2(S−G) + L_i`.
- proptest로 불변식(멱등성·부피보존·불리언 대수)과 OCCT 오라클 diff.
- *(미구축 — 곡면이 올 때)* relaxation 3단 사다리: f64 Newton → 정체 시 고정밀 → 실패 시 `Tangential`로 상신(삼키지 말 것).
- ★★★★★ **점과 평면의 표현이 어디로 가는지는 `docs/truth-and-cache.md` 가 정한다** — 최종 타입
  구조는 **확정됐다(2026-08-05)**: 그 문서가 도착점(최종 타입·숫자 규칙 7·남은 이행 항목)을
  말하고, 「이행」의 단계표가 진행의 진실이다. **그 방향의 작업(정점 정의·평면 표현·`SurfaceDef`·
  datum 평면·스케치 타입)에 착수하기 전 반드시 읽을 것** — 특히 「반증된 답들」 표(같은 길로 다시
  가지 말 것)와 제약 C1~C8. 지역 코드만 보고 계획을 세우면 문서가 이미 지우기로 한 상태를
  영구화하는 설계를 내놓게 된다(2026-08-04 에 실제로 그랬다).

## 마일스톤 (docs/design.md §8 참조)

M1 뼈대(Store/Handle, Plane/Line, 정육면체, validate, OBJ덤프·STEP출력→기존 뷰어에 위임) → M2 스케치·extrude·replay·STEP → M3 곡선기하(NURBS·Arc·Cylinder, tess 출처태그) → M4 면 위 작업(imprint·pad, 여기까지 전부 Constructed, 오라클 가동) → M5~M7 불리언 사다리(이 즈음 워크스페이스 밖 인터랙티브 디버그 뷰어 — §design 1).

**인프라(검증·테스트)가 코드보다 먼저**라는 게 이 프로젝트의 방법론이다 — M1이 `Store<T>`/`Handle<T>`를 proptest와 함께 시작한 이유이고, 지금도 새 단위마다 같다. **현재 위치는 `docs/dev-log.md`의 마지막 셀이 말한다**(이 파일에 적으면 반드시 낡는다).

## 작업 스타일

- 큰 덩어리 말고 작은 단위로. 각 단위는 테스트를 남기고 끝낸다.
- 워크스페이스 구조. 크레이트는 마일스톤 따라 `cargo new`로 추가(지금 9개 다 만들지 말 것).
- 커밋은 작게, 의미 단위로.
- 기하 코드 버그는 눈으로 잡는다 — 애매하면 OBJ 덤프해서 확인.
- **셀별 진행 기록은 `docs/dev-log.md`에 append한다** — 설계 규칙·불변은 `design.md`, 오른 사다리·반증된 예측 등 셀 단위 로그는 dev-log로 분리(design.md 재비대 방지).

### 관문 (칸을 닫기 전에 도는 것) — **여기가 유일한 원본이다**

★ 이 목록은 지금까지 칸마다 dev-log에 다시 적혀 왔고, 그래서 표류했다. 절차는 여기서 고치고
dev-log에는 **결과**만 적는다.

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --no-fail-fast
cargo test -p nacre-ops --no-default-features
cargo test -p nacre-ops --test census -- --ignored --nocapture | grep '^c '          # 두 프로파일 diff
cargo test -p nacre-ops --release --test census -- --ignored --nocapture | grep '^c '
cargo test -p nacre-ops --test reject_census
cargo test --workspace --no-fail-fast -- --ignored \
    --skip boolean_wall_clock --skip profile_check_wall_clock
cargo test -p nacre-ops --release --test perf -- --ignored --nocapture   # 성능은 release로 따로
```
그리고 kit(`../nacre-kit`: fmt·clippy·test)과 앱(`../nacre-playground/web`: `npm run wasm:all`
+ `npx tsc --noEmit` + `npx vitest run`, 그리고 `wasm/`에 clippy).

★★ **`perf.rs`는 스윕에서 빠지고 release로 따로 돈다** (2026-08-22 실측). 스윕 **1802초 중
1673초(93%)가 그 파일의 테스트 둘**이었고, 그나마 debug 빌드라 시간 숫자는 의미가 없었다.
release로 옮기면 30분 → **6분 20초**(스윕 2분 10초 + perf 4분 8초).
★ **그런데 그 둘은 시간을 단언하지 않을 뿐 «커버리지»다** — `.expect("fuse")`로 큰 회전 fold를
수십 번 쌓으므로 깨지면 패닉한다. release로 옮기면서 잃는 것은 그 fold들의 **`debug_assert!`**
(release에선 꺼진다)이므로, **가끔은 debug로도 한 번 돌린다**:
`cargo test -p nacre-ops --test perf -- --ignored` (28분).

★★★ **스위트가 느린 이유는 «양»이 아니라 «최적화»였다 (2026-09-08).** 워크스페이스 `Cargo.toml`이
`[profile.dev]`·`[profile.test]`에 `opt-level = 2`, `debug = 1`을 든다. 이 커널의 안쪽 루프는
정확 산술이라 `opt-level = 0`에서 **아홉 배** 느리다 — 같은 세션 A/B: ops 라이브러리 테스트
**68.4s → 8.0s**, 전체 스위트 **475s → 80.9s**, census **90s → 1.9s**, 키트 **22.6s → 0.5s**.
**편집-컴파일 루프는 안 느려진다**(한 줄 편집의 증분 컴파일 1.22s → 1.14s: 증분 컴파일이
옵티마이저 비용을 콜드 빌드 한 번으로 몰아준다). census는 **비트 동일** — 최적화가 판정을 바꾸지
않는다. `debug = 1`은 라인 테이블만 남기므로 패닉의 `file:line`(그건 DWARF가 아니라 `Location`)은
그대로고, 사라지는 디버거용 변수 정보가 **산출물 크기의 대부분**이다(바이너리 32MB → 11.7MB).

★★ **카고는 옛 세대를 절대 안 지운다** — `target/debug`가 **한 세션에 9.1GB** 자라는 것이 실측됐고
(`incremental/` 4.2GB + `deps/` 5.9GB, 그중 편집마다 남는 테스트 바이너리 1.3GB), 여러 세션이면
70~90GB가 된다. **주기적 `cargo clean`은 필요하다**(위 `debug = 1`이 주기를 늘려 준다).
★ 그리고 **카고는 한 번에 하나만** — 겹쳐 돌리면 서로 느려진다(실측 68.4s → 407.8s, **6배**).
커밋 훅이 스위트를 돌므로, 커밋 중에는 다른 cargo를 시작하지 않는다.

★ **관문 명령은 파이프로 감싸지 말고 파일로 받아 cargo 자신의 `$?`를 읽는다** — `| tail -N`의
exit code는 `tail`의 것이고, 잘린 출력에는 읽을 것도 안 남는다(2026-08-22에 두 번 당했다).
