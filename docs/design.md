# nacre — 하이브리드 CAD 커널 설계

정확한 기하를 진실로 보관하는 고전 b-rep의 골격 위에, Fornjot에서 검증된 위생 규율(append-only 단일 참조, 위상과 근사의 동시 구축)을 얹은 설계다. 목표는 정밀 기계 CAD(STEP 입출력 포함)이며, 개인 + AI 협업 개발을 전제로 실수 여지를 줄이는 인프라를 처음부터 포함한다.

# 어디서나 닿는 규칙

## 설계 원칙

이 문서의 모든 결정은 다섯 가지 절대 원칙에서 나온다. **원칙과 범위(비목표)의 유일한 원본은 `docs/overview.md` 의 「절대 원칙」·「범위」다** — 여기에 다시 적지 않는다(두 곳에 적으면 한쪽이 상한다). 원칙 4 의 구분을 드는 타입은 정점의 실현 캐시 `PointCache`(`nacre-topo`)이고 변종 셋의 뜻은 그 doc 이다. 파라메트릭 편집이 v1 비목표인 까닭(연산이 원시 Handle 을 참조하면 상류 수정이 하류 번호를 밀어낸다 — topological naming)과 진화 경로 `OpRef` 는 「연산 층」 절.

## 크레이트 구조

의존 방향은 아래에서 위로만 흐른다. 순환 의존 금지. 크레이트마다 맡는 일의 원본은 그 크레이트의 `//!` doc 이다.

```
nacre/                    # 워크스페이스(crates/ 아래). 최상위 `nacre` 크레이트는 파사드(재수출 전용)
├── nacre-store      # typed-index 저장소: Store<T>/Handle<T> (기하·위상 무지, nacre-* 의존 없음)
├── nacre-math       # 자체 선형대수: Point<D>·Vector<D> (nacre-* 의존 없음)
├── nacre-exact      # 유리수·정수로 정확히 답할 수 있는 모든 것: Rat·Angle·Isometry·PlaneName·CylinderDef·QuadVal·QuadRoot·정확 술어 (nacre-* 의존 없음)
├── nacre-predicates # f64 입력의 exact 부호 술어(indirect predicates); geometry-predicates 위·standalone (nacre-* 의존 없음)
├── nacre-judge      # 정확히 답할 수 없는 회전 좌표의 부호 판정(CIP): f64 필터→상승→«증명 또는 미결»
│                    #   kernel(WitnessPoint) + predicate(질의를 이름 경로나 kernel 로 보내는 라우터)
│                    #   + cylinder(두 평면이 원통과 만나는 점 PiercePoint 와 그 질문)
│                    #   ← exact · predicates · math
├── nacre-geom       # f64 기하 «실현»(Surface·Curve)·교차(intersect)·호가 낀 링의 정확 술어(mixed)
│                    #   ← math · exact(유리수 링 술어)
├── nacre-topo       # b-rep 위상: Vertex/Edge/Face/Shell/Solid·half-edge·Model(진실과 그 캐시)   ← store · geom · math · exact
├── nacre-tess       # tessellation: 출처 태그   ← store · topo · geom · math · predicates · props(OBJ 의 면 법선)
│   └── polygon      # 평면 다각형 삼각분할: y-단조 분해 + 단조 삼각분할 + Delaunay 플립
├── nacre-validate   # 불변식 검사: 참조·위상(닫힘·다양체·방향·오일러-푸앵카레)·기하   ← topo · store · math · geom · exact
├── nacre-props      # 질량 특성과 파생 질의: 부피·면적·중심·상자·면(해석적, tess 무관)   ← topo · geom · math · store
├── nacre-ops        # 연산: datum 평면/sketch/extrude/boolean/transform/mirror/copy
│                    #   불리언 = 평면 클래스당 셀 복합체 하나 — 「연산 층」·「원통 — 밴드·게이트·차트·룰링」 절
│                    #   불리언은 한 방향이다: boolean(앞문) → arrangement(엔진) → assembly(조립),
│                    #   셋 다 draft(지어지기 전 면의 어휘)를 읽는다 — module_graph 시험이 단언한다(overview 「관문」)
│                    #   ← topo · geom · math · store · exact · judge (+ 선택적 rayon)
├── nacre-step       # Model→STEP(AP242) 내보내기 어댑터   ← exact · topo · geom · math · store (+ brep-to-step)
├── nacre-oracle     # [dev, publish = false] OCCT 비교 하네스   ← topo · step · math · store
├── nacre            # 파사드   ← exact · geom · math · ops · props · step · store · tess · topo · validate
└── tools/occt-helper/  # 워크스페이스 밖 헬퍼: brew OCCT(1순위) 또는 uv+OCP(폴백) — 「검증·오라클 인프라」 절
                        #   OCCT는 오라클 전용 — 제품 경로에 위임 없음
```

**의존의 규칙.**
- 제품 경로에서 `nacre-ops` 는 `nacre-tess`·`nacre-validate`·`nacre-props`·`nacre-step` 에, `nacre-exact` 는 `nacre-predicates` 에 의존하지 않는다 — 시험 전용 dev-dependency 로만 쓴다(발행되는 라이브러리의 의존 표면이 깨끗하게 남는다). 형제 사이의 dev-dependency 를 path 로만 적는 까닭은 루트 `Cargo.toml` 의 주석.
- 기하는 위상을 모르고(순수 수학), 위상은 기하를 **캐시 값으로** 든다(`SurfaceCache` 의 팔이 geom 의 `Plane`·`Cylinder` 를, `EdgeCache` 가 `Curve` 를). 진실은 topo 자신의 타입이고 geom 은 그 실현이다.
- `Store`/`Handle` 은 최하위 `nacre-store` 에 있어 위의 크레이트가 자유롭게 참조한다. 곡면 아레나가 드는 진실(`Handle<MotionNode>` 를 든 topo 타입)을 topo 아래 크레이트는 이름 지을 수 없으므로 **`nacre-geom` 은 `Handle` 을 갖지 않는다.**
- `nacre-predicates` 는 어떤 `nacre-*` 에도 의존하지 않는다 — 평면 계수·좌표를 평범한 배열로만 받아, nacre 밖으로 떼어 내 단독으로 공개할 수 있게 둔다(「`nacre-predicates` 는 순수 수치층이다」).
- **외부 의존.** exact 유리수는 손수 짓지 않는다 — 손수 지은 유리수의 버그는 exactness 목표를 훼손한다. `nacre-exact` 는 `num-rational`(+`num-traits`·`num-bigint`·`num-integer`)과 `astro-float` 를 쓴다(MIT/Apache·순수 Rust; 그 매니페스트의 주석).
- **라이선스.** 커널은 MIT/Apache-2.0 듀얼이다(Manifold(Apache-2.0) 알고리즘 차용과 호환). indirect predicates 의 참고처는 논문(Attene 2020, arXiv 2105.09772; Shewchuk 1997; Lévy PCK)과 `geometry-predicates` 소스로 한정하고, LGPL 인 Attene 참조 구현의 소스는 **작성 중 열람 금지 / 완성 후 실행 대조만 허용**이다(「마일스톤 사다리」의 indirect predicates 문단). 이차곡면 참고처(QI·CGAL Spherical Kernel)는 소스를 열람조차 하지 않는다 — 「참고 논문과 라이선스 규율」. OCCT 는 overview 「OCCT / 외부 코드 규칙」.

**어느 크레이트에 두나 — 축은 «주제»가 아니라 «정확성»이다.** 「유도된 값 + 그 산술」(유리수·정수로 정확히 답하는 값·이름·술어, 그리고 ℚ 밖이되 정확한 `QuadVal = a + b√c`) → `nacre-exact`, 그 f64 실현 → `nacre-geom`, 「아레나 항목의 진실」 → `nacre-topo`. 주제로 그린 지도(「숫자·기하·위상」)는 경계를 설명하지 못한다 — `nacre-exact` 에는 기하의 이름과 술어(`cylinders_nested`·`SeamOrder`·`PlaneName`·`CylinderDef`)가 든다. 이것을 강제하는 것은 취향이 아니라 의존 그래프다 — **`nacre-judge` 는 `nacre-topo` 에 의존하지 않으므로** 판정이 쓰는 타입은 topo 아래에 있어야 한다(`PlaneName`·`CylinderDef`·`QuadRoot` 가 topo·ops 가 아니라 exact 에 사는 이유 — topo 는 `CylinderDef` 를 재수출한다).
- `nacre-exact` 와 `nacre-predicates` 는 입력의 수로 갈린다: 유리수·정수 위의 값과 술어는 exact, f64 좌표 위의 적응 부호 술어는 predicates 다(exact 의 정수 부호 술어는 predicates 의 f64 쌍둥이와 dev 에서 차등 시험된다). 정확히 답할 수 없는 회전 좌표의 부호는 `nacre-judge` — predicates 의 toleranced 쌍둥이다(「CIP — 허용오차 부호 술어」).
- 견고성이 첨예한 코드는 한 모듈에 있지 않다: f64 닫힌 꼴 교차와 직선 링 술어는 `nacre-geom::intersect`, 호가 낀 링은 `nacre-geom::mixed`, 부호는 위 셋, 셀의 내/외는 `nacre-ops::nesting` 이다. 어디서든 판정(부호)은 좌표가 아니라 정의로 한다(overview 「정밀도 인프라」).
- 크레이트를 나눌지는 「한쪽만 쓰는 소비자가 있는가」로 잰다.

**파사드 `nacre`.** 사용자는 파사드 하나만 의존하고, 약속하는 API 는 파사드뿐이다(`crates/LAYER-README.md`). 스크립트 앱·편의 레이어·디버그 뷰어는 워크스페이스 밖에 둔다.
- 층을 **모듈로** 재수출하고(`nacre::exact`·`geom`·`math`·`ops`·`props`·`step`·`store`·`tess`·`topo`·`validate`), 자주 쓰는 것은 `nacre::prelude` 에 담는다. **평면(flat) 재수출은 하지 않는다** — 층 분리가 이 설계의 뼈대이고, 같은 맨이름이 층마다 다른 것을 뜻한다(`Surface` 의 두 철자 — 「타입이 사는 자리 — 규칙」).
- **공개 API 에 나오는 아래 층의 타입은 누군가 재수출한다.** `Operation::Transform { isometry: Isometry }`·`Mirror { axis: Axis, offset: Rat }` 가 exact 타입을 싣으므로 파사드가 `nacre-exact` 를 재수출한다 — 없으면 소비자가 그 op 를 만들 수조차 없다(`the_whole_pipeline_runs_through_the_facade_alone` 이 컴파일로 잠근다). `nacre-judge` 의 이름은 그것을 싣는 ops 가 재수출하고(`nacre_ops::{Decision, Evidence, Site}`) `nacre-predicates` 의 타입은 공개 API 에 나오지 않으므로, 파사드는 두 크레이트를 모듈로 내보내지 않는다 — 발행은 필요하다(judge 는 ops 의, predicates 는 judge·tess 의 하드 의존).
- **prelude 의 기준은 확인 가능한 성질이다**: «`Operation` 의 모든 변종이 prelude 이름만으로 만들어진다». 내용은 실제 소비자의 import 와, 그 소비자보다 나중에 생긴 표면(스케치 앞문·파생 조회·거절 사유)의 합집합이다 — 낡은 소비자만 보면 최신 API 가 빠진다.
- **문서 예제와 예제는 프로덕션 API 만 쓴다**(`test-util` 이 필요한 예제는 독자가 실행할 수 없다). 예제는 파사드에 둔다 — 층 크레이트의 예제도 crates.io 에 실려 독자가 따라 하는 길이 된다. 검증은 «`nacre::` 경로만으로 끝까지 가기»다: 파사드 doctest, `crates/nacre/tests/facade.rs`, `crates/nacre/examples/plate_with_bore.rs` — 재수출이 빠지면 컴파일이 깨진다.
- **기능.** `parallel`(기본 on)은 `nacre-ops/parallel` 로 전달되고 judge 의 고정밀 캐시를 `Sync` 로 만드는 스위치이기도 하다; `test-util` 은 topo 의 시험 전용 문과 push 계측을 전달한다 — 전달은 파사드 `Cargo.toml` 의 주석, 그 문과 계측(`WIDE_PLANES`·
  `SEEDED_HITS`·`surface_derive_counts` — 제품 빌드는 아무것도 세지 않는다)은 `nacre-topo/Cargo.toml` 의 주석. 소비자가 rayon 을 끌 수 있으려면 파사드는 `nacre-ops` 를 **`default-features = false` 로** 매단다; 이 속성은 매니페스트에 살아 시험이 `Cargo.toml` 을 직접 읽는다(`crates/nacre/tests/facade.rs` 의 `the_ops_dependency_stays_default_free_so_consumers_can_drop_rayon`). **순차 조합은 `-p nacre-ops --no-default-features` 로 확인한다** — 워크스페이스를 통째로 지으면 형제 크레이트(`nacre-oracle`·`-props`·`-step`·`-tess`·`-validate`)의 dev-dependency 와 파사드의 기본이 기본 피처를 도로 켜서 아무것도 안 재게 된다. 그 명령이 순차인 것은 **크레이트가 제 자신에 거는 dev-dependency**(통합 시험이 제 `test-util` 을 쓰려고 거는 한 줄 — `nacre-ops`·`nacre-topo`)**가 `default-features = false` 이고 `default`·`parallel` 을 이름으로도 고르지 않기** 때문이다. 기본을 켜 두면 `--no-default-features` 아래에서도 기본이 도로 켜져 순차 빌드 대신 병렬 빌드가 초록으로 돈다. 자기 dev-dependency 를 새로 거는 크레이트도 이렇게 적는다(`crates/nacre-ops/tests/instruments/manifest.rs` 가 모든 매니페스트에 건다).
- **`Document`(op 로그 + Model + tess 캐시) 번들은 두지 않는다** — 「저장소와 동일성」 절이 이 상위 번들을 지목하지만 의미론을 정하는 새 타입이라, kit 이 무엇을 원하는지 보이기 전에 정하지 않는다(「회전 지원의 확정 결정」 4).

**편의 레이어 `nacre-kit` (워크스페이스 밖, 별도 리포).** 코드-CAD 스크립트와 커널 사이의 층: 다중 솔리드 값(compound), 값 의미론(재사용 시 `Copy` 자동 삽입), 다인수 fuse/cut/common(fold), 프로파일 헬퍼와 섬-분해 호출, 패턴·미러, 에러의 사람용 매핑, 표시 메타데이터(색·투명도 — 커널 비목표라 여기가 제자리). 경계 규칙은 overview 「범위」의 «설탕 vs 커널 판별 기준»이다. **Rust 로 두는 이유:** 헤드리스 `cargo test` 가 되고, 술어 인접 로직이 exactness 도구가 있는 쪽에 남고, 프론트엔드를 바꿔도 살아남고, wasm 경계가 함수 하나로 유지된다. 스크립트 문법은 `nacre-playground` 가 구현하고, 의미론과 그 이유는 `nacre-kit` 의 코드와 주석에 있다(여기에 옮겨 적지 않는다).

**공개 표면 — 밖에서 무엇이 되는가.**
- **`Model` 의 저장소는 비공개이고 문(accessor)으로 읽는다.** 공개 필드가 없고, 전량 순회 문은 두지 않는다 — 아레나는 supersede 된 항목도 들므로 소비자는 live 면을 걷는다.
- **위상 순회와 내부 정보는 밖에서 읽힌다.** 위상 셀의 필드는 `pub` 이고 `Model::reachable` 도 공개다. STEP 에 안 담기는 내부 정보 — 정점의 정의 `Vertex` 와 실현 캐시 `Model::vertex_cache`, `Model::motion`(필드가 아니라 좁은 문), 메시의 출처(`Tessellation` 의 `by_face`·`by_edge`, `TessTriangle.face`, `TessOrigin`) — 도 공개 문으로 읽힌다: 커널 내부를 보는 디버그 뷰어는 커널 크레이트가 아니라 nacre 를 직접 링크하는 워크스페이스 밖 앱의 자리이고(연산 로그를 스텝별로 재생한다 — append-only 라 «N번째까지 replay» 가 공짜다), 그것이 읽을 것은 공개 문에 있어야 한다. 뷰어 없이 하는 시각 확인은 「마일스톤 사다리」의 M1.
- **파생 질의는 값을 돌려주는 읽기 전용이다**(위상 순수성 유지): `nacre-props` 의 `bounds`·`centroid`·`face_props`, `nacre-ops` 의 `face_plane`, `Model::he_start`. 각 질의의 약속 — 곡선을 인지하는 상자, f64 측정이라 구성에 되먹이지 않는 중심, 곡면에서 `None` 인 법선, 생김새로 면을 골라 위상 명명 문제를 비켜 가기 — 은 그 doc 이다.
- **면의 스케치 프레임**(`face_sketch_frame`)은 면이 어디 있는가로 정하지 평면이 어떻게 저장됐는가(모션을 점에 옮겨 적었나, 노드로 기록했나)로 정하지 않는다. 세계 방정식이 진술된 평면(`Model::world_plane_name` — 모션 없음, 또는 접히는 사슬)은 바깥 법선의 세계 arbitrary-axis 프레임을 갖는다: 원점은 세계 원점의 정사영, `+u = ẑ × n`(법선이 수직이면 `ŷ × n`; `plane_frame_default`). 그 밖의 평면(사슬에 프레임 노드나 사분각 밖 회전이 든 것, 넓은 이름, 이름 없음)은 제 `Canonical` 프레임을 사슬로 나른다. 어느 쪽이든 원점은 평면이 정하고 면(꼭짓점·외곽선)을 읽지 않는다. 버린 두 규약(면에서 읽는 원점, `any_perpendicular` 축)은 「가지 말 것」.
- **`face_plane` 은 그 프레임 위의 `Extrude` 가 프로파일을 놓는 바로 그 프레임이다** — `face_sketch_frame` 의 프레임을 실현한 것이고 둘째 유도를 짓지 않는다. 갈라지면 앱이 계산한 자리와 보스가 어긋난다(두 축으로 비대칭인 프로파일로 잡는다 — `face_plane_is_the_frame_pad_places_profiles_in`; 실현이 `face_plane` 과 비트 동일하다는
  계약은 `tests/invariants/sketch_frame_contract.rs`).
- **`flip` 은 호출자가 진술하지 않는다.** 생성자는 `false` 를 두고, `ŵ` 가 어느 쪽을 보는지는 쓰는 연산이 진실에서 `flip_toward` 한 곳으로 정한다(프레임은 `frame_toward` 로 만든다). 프레임 노드는 `push_frame_node` 로만 push 한다 — 그 문을 extrude 와 offset datum 이 함께 쓴다.
- **로그는 평면을 핸들로 싣는다**(`Extrude { frame: SketchFrame, .. }`; 평면의 진술을 싣는 변종은 `DatumPlane` 하나). 로그 속 평면 핸들이 가리키는 것은 씨앗(`Model::new` 이 심는 세계 축 평면 셋, `seeded_plane`)·이미 지은 면의 평면·`DatumPlane` 이 넣은 평면뿐이다 — 로그 밖에서 지은 평면을 가리키면 replay 가 그 핸들을 재현하지 못해 로그가 자기완결이 아니다.

## 제약 — 어떤 답이든 만족해야 하는 것

| | | 이 구조에서 |
|---|---|---|
| C1 | 진실은 하나, 나머지는 캐시 | 생산자는 점(또는 정점 핸들)만 진술, 이름·좌표·계수는 전부 유도 |
| C2 | base case 는 «수» | `PlanePoints::Known` — 모든 사슬이 유리수 세 점에서 끝난다 |
| C3 | 발견된 좌표는 `Rat` 에 안 들어간다 → 가리킨다 | `PlanePoints::Through` |
| C4 | 모션 섞인 판정은 f64 필터 + 상승 | `Judge` 의 f64 필터 → 고정밀 상승 + `Decision` |
| C5 | append-only ⇒ 참조는 과거로만 (DAG) | Through(평면→정점)는 datum 이 정점보다 나중이라 성립; 벽은 값(`Known`)이라 순환 자체가 없다 |
| C6 | 깊이는 사용자 조작당 1단 | 불리언은 평면을 만들지 않는다 — Frame·Through 재귀가 그래서 유계다 |
| C7 | 거절은 정직하게, 이름 붙여 | `Decision` 의 `Exhausted`·`Degenerate` → 이름 붙은 거절 |
| C8 | `Rat` 을 넓히지 않는다 | 넓어진 것은 이름의 **그릇**(동일성, 그리고 `WideFrame` 으로 서는 프레임)과 중간값뿐 |

---

## 숫자 규칙 — 일곱 개

1. **진실은 «유리수»이거나 «핸들»이다.** `Rat`(i128 분수)에 들어가는 값은 적고, 안 들어가는 값(넓은 이름의
   풀이, 계수 곱)은 저장하지 않고 **가리킨다**. 가리키는 사슬의 끝은 언제나 수다(C2·C3 의 동시 만족).
   **폭은 사용자가 쓴 기하와 그것이 놓인 자리가 정하지 연산 횟수가 정하지 않는다** — 불리언은 평면을 새로
   만들지 않고 피연산자의 평면을 재사용하며(`a_boolean_mints_no_surface`), 모션은 진술로 옮겨 적으면 닿은
   자리에 직접 지은 것과 같은 폭의 이름이 되고 기록하면 계수를 건드리지 않으며, 오프셋은 법선을 고정한 채
   d 만 옮겨 분모가 곱이 아니라 lcm 으로 자란다. `nacre-ops/tests/instruments/point_width.rs`
   (`a_stacked_name_is_as_wide_as_the_place_it_lands`)가 잠그고, 발견된 좌표의 폭도 그 파일이 잰다.
2. **화살표는 한 방향뿐이다.** `Rat → f64`(캐시), `Rat → BigFloat`(판정 상승). `f64 → 진실`은 없다 — f64 를
   들어올려 진실로 삼는 순간 반올림이 진실에 구워진다(「가지 말 것」). 진실의 입구는 `Rat::from_decimal`
   하나다(들어올리기가 아니라 사용자가 *쓴* 십진수의 결정적 정규형). `Rat::try_from_f64` 는 `f64` 를 정확히
   읽어 비교하거나 그림을 나르는 문이지 진실의 입구가 아니다(그 doc).
3. **무리수는 값이 아니라 정의로 존재한다.** 회전은 `Angle`(유리수 degree), 기울어진 프레임은
   `Motion::Frame`(법선을 평면 핸들로 가리킨다)이 들고, cos/sin/√ 는 실현 시점에 정확 반올림으로 나온다.
   **모션은 곡면이 든다**(`Model::surface_motion`) — 정점·모서리는 곡면을 가리켜 정의상 따라온다.
4. **실현 통로는 하나다.** 정의를 호출자가 고른 정밀도로 실현하는 길(`realize_vertex(model, v, Precision)` —
   증인점에서는 `WitnessPoint::coord_hp(prec)`)을 f64 캐시(표시·tess, 그리고 STEP — 사다리를 오르며
   `round_to_f64` 로 정확 반올림), 판정 상승이 전부 탄다. 읽어 내는 규칙도 하나다(`Realization::to_f64`):
   최근접 f64, 또는 0 에서 일치 정밀도 `max(1, 점의 최대 |좌표|)·2⁻¹⁸⁰` 안이 증명된 좌표는 `+0.0` 이다(참값 0 은
   반올림으로는 어느 정밀도에서도 정해지지 않는다). 문턱의 크기는 점 자신의 것이라 모델 크기를 넘지 않으므로
   규칙 5 의 모델 크기 문턱보다 엄격한 쪽이다. 반올림보다 먼저 물으므로 어느 단에서 결정해도 답이 같다.
   캐시는 값과 함께 «무엇이 증명됐나»를 든다 — 정점은 `PointCache`(`Realized` 만 축별
   경계), 곡면은 `CacheStanding`; 저장된 tol 은 없다(`PointCache` doc). **f64 좌표는 즉석에서 다시 풀지
   않는다 — 표에서 읽는다**(판정마다 정의를 다시 지으면 회전 불리언 시간의 77% — `plane_def` doc,
   「가지 말 것」). 메시가 곡선 위에 새로 찍는 샘플점은 이 통로 밖이다 — 보기의 수라 libm 으로 짓는다.
   최근접 f64 는 STEP 의 약속이고, STEP 은 원을 곡선으로 내보내 샘플을 쓰지 않는다.
5. **판정은 3단이고, 결과는 3갈래다.** 술어마다 제가 닿는 평면만 보고 정확한 길부터 묻는다. 이름 정수 도로
   (`name_rows` — 그것을 묻는 술어만)는 모든 평면에 이름이 있고 `is_rotated` 가 모두 거짓이거나 전부 한 모션일 때
   열린다(`cmp_coord` 는 앞의 것만): f64 필터 → `nacre_exact::int_*`. 몇몇 술어는 제 정확한 길을 따로 둔다(예:
   `cmp_coord` 의 `cancel_cmp_coord`, `planes_coplanar` 의 세계 이름 비교와 `coplanar_by_composed_rotation`) — 각 술어의
   doc. 그 밖은 커널: 구간 필터 → (직접 `orient3d` 이고
   네 점이 한 모션을 공유하고 바탕이 `f64` 이며 구간 필터가 못 가르면) 모션 전 좌표의 `nacre_predicates::orient3d` →
   astro-float 상승. 커널에 가지 않고, 증명하지 못하면 모른다(`None`)고 답하는 술어도 있다(`planes_parallel`). 결과는:
   - 부호가 **증명되면** → 그 부호(`Sign`).
   - 부호는 못 갈랐지만 **«일치 정밀도보다 가깝다»가 증명되면** → 일치로 처리하되 근거를 실어
     **보고**한다(`Decision::Coincident{within}`, `boolean_with_report`). 일치 정밀도
     (`ProofStandard::coincidence` = 모델 크기 × 2⁻¹⁸⁰ — f64 출력 해상도 2⁻⁵² 보다 두 워드 아래)는
     **모델에서 유도**되는 값이지 사용자 설정이 아니다(유도와 아래로 미는 비대칭은 `coincidence_for` doc).
     전역 tolerance(*"가까우면 붙여라"*)와 방향이 반대(*"가깝다고 **증명되어야** 일치"*)다.
   - 둘 다 증명 못 하면 → **이름 붙은 거절**(`JudgeExhausted`·`DegenerateWitness`·
     `PrecisionBudget`). 조용한 0은 없다(C7).

   정밀도는 상수가 아니라 모델이 정하고, 필요한 비트를 계산해 **한 번에** 점프한다.
6. **동일성은 정준 이름의 `==` 다 — 이름의 그릇은 임의정밀이다.** 이름은 하나만 저장한다:
   `PlaneName = Narrow([Rat;4]) | Wide([BigInt;4])`, i128 에 들어가면 반드시 `Narrow`(정규화 불변식 —
   `plane_name_exact` 가 중간 계산을 임의정밀로 닫아 지킨다). 동일성은 enum 전체의 `==` 다. 산술 지름길은
   `narrow()` 투영을 빌려 읽고 `Wide` 에서 물러서지만, 프레임은 넓은 이름에도 임의정밀 길(`WideFrame`)로 선다.
   좁은 형태(법선·`n·n`)가 없다고 **진실을 거절하지 않는다** — 필요하면 점 셋에서 임의정밀로 실현한다.
   «항상 이름을 갖는다»는 유리수 닫힘인 평면에만 적용된다 — 세 정점이 한 프레임에서 만나지 않는(혼합 프레임)
   `Through` 평면은 이름이 없다(세계 계수가 `z = 0` 처럼 유리수여도). 그런 평면은 **진술 키**
   (`SurfaceKey::Through` — 정렬 삼중항 + 모션)로 intern 되어 «같은 진술 = 한 핸들»은 지켜지고, 같은 평면이
   이름 있는 핸들과 이름 없는 핸들 둘이 될 수 있다. 기하 동일성은 술어가 답한다: 클래스 발견이 핸들이 다른
   평면 쌍마다 `planes_coplanar` 를 물어(정확 0, 또는 일치 정밀도 안이 증명되면 보고와 함께) 한 클래스로 합치고,
   결과의 기하는 한 핸들일 때와 비트 동일하다 — 느릴 뿐이다(`one_plane_two_handles.rs`,
   `a_plane_stated_twice_booleans_as_it_does_once`).
7. **`Rat` 을 넓히지 않는다(C8).** 넓어지는 것은 중간값(BigInt 이름 유도·`Expansion` 술어·
   BigFloat 상승)뿐이다. 근거는 규칙 1 이다 — 저장되는 값의 폭은 사용자가 쓴 기하가 정하고 연산을 쌓아도
   자라지 않으므로, 넓은 그릇이 필요한 자리는 저장이 아니라 계산이다.


### 수의 이름 — 접미사는 결과가 어떤 수이고 얼마나 믿을 수 있는지를 말한다

같은 값의 다른 숫자 판이 있거나 정확한 값에서 수를 만드는 함수는 이름 끝이 그 결과를 말한다. 숫자 판이 하나뿐인
접근자(`SketchPlane::origin()` 등)와 `nacre-math` 는 달지 않고, `nacre-geom` 은 f64 판을 맨이름으로 두고 `Rat` 판에
`_rat` 을 단다(`three_planes` · `point_in_ring_2d_rat`). 판정의 칸은 모두 `_iv` 이상이다 — `_approx` 는 진실에 대한
부호를 정하는 산술에 들어가지 않는다(대략값·계측·보기, 정밀도 정책의 크기 어림 `standard_from` 의 `scale`, 이름 없는
평면 캐시의 향 맞추기 `align_cache_sense` — 캐시의 향이지 진실의 부호가 아니다).

| 철자 | 뜻 | 예 |
|---|---|---|
| `_rat` | 정확한 `Rat`. 무리수이거나 넘칠 수 있으면 `Option` | `cos_sin_rat`·`sqrt_rat`·`chain_point_rat` |
| `_hp` | 호출자의 `prec` 에서의 `HpBounded` — 값과 건전한 반지름 | `cos_sin_hp`·`inv_sqrt_hp`·`quad_hp`·`coord_hp` |
| `_iv` | f64 구간 `Bounded` — 값과 증명된 반지름 | `orient3d_iv`·`det3_iv`·`coord_iv` |
| `_f64` | 참값의 최근접 `f64`(규칙 4), 또는 반올림 없이 정확한 `f64`. `nacre-exact` 의 근사 실현(`cos_sin_f64`·`inv_sqrt_f64`·`sqrt_f64`·`unit_vector_f64`)은 값이 정해질 때까지 사다리(`ladder_f64`)를 오르고, 0 이 아닌 값의 최근접 `f64` 가 정상수가 아니면(정상 범위 밖) `NaN` 이다 — 정확 갈래도 같다(`nearest_normal_f64`). ops 의 `Realization::to_f64` 는 이 규칙 밖이다(정확 갈래는 subnormal 을 통과시키고 근사 갈래는 `+0.0` — 그 doc). 앞의 낱말이 약속을 바꾸면 그 낱말이 말한다 — `upper_f64` 는 그 경계 이상인 `f64` | `cos_sin_f64`·`sqrt_f64`·`unit_vector_f64`·`Rat::to_f64`·`Mag::upper_f64` |
| `_approx` | 약속 없는 `f64` — 여러 번 반올림했거나 캐시·대략값을 `f64` 로 계산했다 | `QuadVal::to_f64_approx`·`Isometry::point_approx`·`coord_approx` |
| `_big` | 좁은 쌍둥이가 있을 때의 넓은 입력판, 숫자 타입 뒤에 둔다. BigFloat 출력은 `_bigfloat` | `sqrt_rat_big`·`inv_sqrt_hp_big`·`rat_to_bigfloat` |
| `_uncached` | memo 문 뒤의 본체 | `cos_sin_hp_uncached`·`coord_hp_uncached` |
| `_exact` | 정확함에 대한 말 — 반올림 없이 결정된 답, 정확한 변종이나 반지름 0 의 구간, 정확한 피연산자를 전제로 한 연산. 반올림하는 함수에는 쓰지 않는다 | `plane_name_exact`·`is_exact`·`HpBounded::exact`·`div_exact` |

정밀도 인자는 exact·judge 에서 `prec`, ops 에서 `bits`, 십진 자릿수는 `places` 다.

**바닥.** 다음 자리는 `f64` 로 읽은 오차 경계를, 0 이 아니면 `f64::MIN_POSITIVE`(가장 작은 정상수) 아래로 내리지
않는다: `Mag::upper_f64`·`bigfloat_upper_f64`(둘이 `pow2_upper_f64` 를 나눠 쓰고 정확한 0 은 0 으로 읽는다),
`Bounded::round_off`(0 도 내지 않는다 — 반올림이 정확했는지 아는 호출자가 0 을 준다), `cyl_gate::round_slack`,
`Realization::to_f64` 의 `half_ulp`, STEP 불확도의 `half_ulp`, 계측 `audit_rim_centre`. 비정규수 구역의 반올림 오차는
값에 비례하지 않고 고정된 반 칸(`2⁻¹⁰⁷⁵`)이라 상대 오차로 쓴 식이 거기서 성립하지 않기 때문이다: `|x|·ε/2` 는 0 이
되고(정확하다는 거짓) `× (1 + kε)` 는 대개 값을 키우지 못한다.

**바닥 없는 식은 정상 범위에 기댄다.** f64 구간 산술 — `Bounded` 의 `add`·`sub`·`mul`(값에 `2ε·|value|`, 반지름의 자기
반올림에 `× (1 + 4ε)` — 연산마다 증명된 한계다; 증인점의 tol 전파도 이것을 걷는다) — 과 `inv_sqrt_error_of` 의 정확한 차는
값이 정상 범위에 머문다는 전제에 기댄다. 정확한 차는 0 이 아닌 `Rat` 이라 크기가 `2⁻¹²⁷` 위에 있어 전제가 타입으로 서고,
구간 산술의 입력 좌표는 i128 분수라 크기가 대략 `2⁻¹²⁷` 위에 머물며 몇 번의 곱으로는 정상 범위를 벗어나지 않아 전제가
선다(논증).

**동사.** `realize_*` 는 모델 안의 정의(정점·담체 곡면)를 값으로 바꾸는 문이고 정밀도를 늘 명시한다(`Precision`·`bits`·
`places` — 기본값 없음). 정밀도 대신 push 의 싼 예산을 묻는 문은 이름이 그렇게 말한다(`realize_vertex_cache`·
`realize_surface_cache` — 계측이 쓰는 문이다). 값(스칼라·quad·계수)과 한 타입의 자기 읽기는 명사 + 접미사로 쓴다
(`quad_hp`, `WitnessPoint::coord_hp`). 깔때기는 `push_*_realized`, 캐시를 더 정확하게 올리는 것은 ops 의 `raise_*`(실현해서
쓴다)와 topo 의 `refine_*_cache`(쓰기 문), 내보내기 전의 문은 `refine_caches(_of)`, 읽고 대조만 하는 계측은 `audit_*` 다.
예산은 `Budget::{Cheap, Paid}` — push 가 내는 싼 예산과 refine 문이 내는 유료 예산이다. `_at` 은 이 커널에서 위치의
낱말(`point_at`·`normal_at`·`surface_handle_at`)이라 정밀도에 쓰지 않는다. 판정은 실현 문을 부르지 않고 자기
어휘(`coord_hp`·`trial_bound`·`precision_for`)를 쓴다.

---

## 타입이 사는 자리 — 규칙

- 두 길이나 두 상태의 기하를 비교하는 회귀 시험은 **위상 정확 일치 + 좌표 ε(모델 크기 상대, `2⁻⁴⁰` — 여유 2000배)** 를
  단언하고, 합격/불합격이 아니라 **최대 편차 숫자를 찍는다**(누적 드리프트 감시 — `plane_anchor.rs`). 비트 동일은 보너스
  신호다.
- **비트 동일 관문은 그 안의 population 에 대해서만 보증한다.** 바꾸려는 것이 닿는 population 을 대장에 먼저 넣고
  (17자리 `fw`·기울어진 `tp` 가족이 그래서 있다 — `census.rs` doc), *"어느 population 인가"* 는 추측하지 말고 계측이
  이름을 대게 한다. 좌표 관문은 «답은 같은데 더 나쁜 길로 갔다»를 원리적으로 못 보므로 **"정확 경로를 탔는가"를 직접
  단언하는 테스트**를 함께 둔다.
- **두 `Surface` 의 철자**: 한 파일에서 맨 `Surface` 는 한쪽만 뜻하고 다른 쪽은 경로로 적는다. 진실을 다루는 파일은
  **맨이름이 진실**(topo, `Handle<Surface>` 가 이름 짓는 것)이고 실현을 `nacre_geom::Surface` 로 적는다; 캐시를 읽는
  `nacre-tess/src/lib.rs` 와 픽스처 `nacre-ops/tests/coverage/features.rs` 는 거꾸로 맨이름이 geom 이고 진실을 `nacre_topo::Surface` 로 적는다.
  한 파일에서 둘 다 맨이름이면 `Handle<Surface>` 는 진실을, 맨 `Surface::Plane(p)` 는 캐시를 뜻해 한 단어가 표지 없이 두
  가지가 된다. 섞으면 **하드 오류**다(geom 은 tuple 변종, 진실은 struct 변종 ⇒ E0532) — 조용히 틀릴 자리가 0.
- 모델에서 판정 표를 짓는 타입(`WorkingPlane`·`WorkingCyl`·`FaceInfo`)은 `nacre-judge` 와 `nacre-topo` 를 둘 다 아는
  `nacre-ops` 에 산다. 어느 크레이트에 두나는 「크레이트 구조」 절.
- **병렬 불변식**: 병렬 구간에서 **곡면 push 금지**(재생 결정성). 곡면의 문은 `Model::push_plane`·
  `push_plane_through`·`push_cylinder`(사설 `push_plane_raw`·`push_cylinder_raw` 로 모인다). 워크스페이스의 진짜 병렬은
  `nacre-ops/src/par.rs` 한 파일(`(0..n).into_par_iter().map(f)`)에만 살고, 그 클로저는 「인덱스 → 값」(`Fn(usize) -> R`)이라
  **`&mut Model` 을 들 수 없다** ⇒ 병렬 구간의 곡면 push 는 구조적으로 불가능하다 — 이것이 잠금이다. `clippy.toml` 은
  rayon 의 병렬 입구를 그 파일에 묶는다(*"병렬이 두 곳에 살면 조용히 깨진다"*;
  아래 「캐시에 쓰는 문은 깔때기만 부른다」). `par.rs` 의 허용이 그 파일에서 끄는 범위는 `disallowed.rs` doc.

---

## 문의 이름 — `Model` 의 접근자 명명 규칙

저장소가 비공개이므로 밖은 **문(=`Model` 의 메서드)** 으로만 읽는다. 규칙:

> **`x(h)` 는 «진실»(아레나 항목 그 자체), `x_cache(h)` 는 «캐시»(실현). 실체 이름을 단 맨 메서드가
> 캐시를 돌려주는 일은 없다.**

문의 가족(서명은 `nacre-topo/src/model/`):

- **진실** — 실체 이름 그대로: `surface`·`vertex`·`edge`·`face`·`shell`·`solid`·`motion`, 진실의 조각
  `surface_motion`(곡면이 진술된 모션), 루트 `live_solids`, 심은 세계 평면 `seeded_plane`.
- **개수·핸들 복원** — `x_count()`·`x_handle_at(index) -> Option<Handle<X>>`.
- **캐시** — `surface_cache(h) -> &SurfaceCache`(종류마다 팔; 조각은 체이닝 — `.plane()`·`.cylinder()`·`.name()`·
  `.standing()`·`.to_geom()`), `vertex_cache(v) -> &PointCache`(`Realized | Ceiling | Unrealized`; `.coord()`·`.bound()`).
  조각의 지름길: `surface_cache_standing`·`surface_name`·`vertex_point`·`edge_curve`·`line_direction_cache`.
- **세계 진술** — 진술된 프레임에서 세계로, 모션 사슬을 접어 정확히: `world_plane_name`·`chain_plane_coeffs`·
  `chain_point_rat`·`chain_dir_rat`.
- **쓰기 문** — `push_*`(곡면은 진실을 진술하며 들어온다 — 그림과 그것이 아는 것을 함께; `push_motion` 은 intern),
  캐시의 둘째 쓰기 문 `refine_vertex_cache`·`refine_surface_cache`(오르기만 한다 — 정점은 `Ceiling` 을 `Realized` 로),
  통째 재생 `rebuild_edge_cache`.

분류·비교·갈래는 진실에 묻고, 표시·tess·경계·측정은 캐시에 묻는다(`Model::surface` doc).

**캐시에 쓰는 문은 깔때기만 부른다.** f64 캐시에 쓰는 topo 의 문(`push_vertex`·`refine_vertex_cache`·`push_plane`·
`push_plane_through`·`push_cylinder`·`refine_surface_cache`·`push_edge`·`rebuild_edge_cache`)과 픽스처용 셋(`*_unregistered`
둘·이름을 심는 `plant_surface_name`), 그리고 깔때기가 다 지은 캐시를 받는 `push_vertex_prerealized` 를 워크스페이스 루트의
`clippy.toml`(`disallowed-methods`)이 막는다 — 크레이트 경계 때문에 `pub` 인 문을 ops 의 깔때기(`realize.rs`) 밖에서 부르면
clippy 가 거절한다. 그래서 «STEP 이 읽는 `f64` 는 실현 길을 지났거나, 지나지 않았다고 표시돼 있다»가 약속이 아니라
잠금이다. 제품에서 부르는 자리(깔때기 모듈, 빈 모델의 세계 평면 씨앗 — 이름이 캐시를 유도한다, seam 표가 실현한 캐시를
넘기는 `reconstruct.rs`)는 `#[expect]` 에 이유를 적고, 테스트 빌드는 크레이트마다 허용한다(픽스처는 문으로 짓는다).
같은 `clippy.toml` 이 rayon 의 병렬 입구를 `par.rs` 에, 링 술어 셋을 제 주인 아래 묶는다 — `point_in_ring` 은
`nesting::cell_inside` 의 재시도, `side_of` 는 `ring_against_plane` 의 걷기, `point_in_mixed_ring` 은 `rational_ring_side` 의
원인 순서. 허용은 하나하나 `tests/instruments/disallowed.rs` 의 표에 이유와 함께 오른다. 그 시험이 거절하는 것(표 밖의
허용, lint 묶음 허용, `clippy.toml` 에서 빠진 경로, 매니페스트·cargo 설정의 lint 수준, 루트 밖의 `clippy.toml`)과 잠금이
닿지 않는 자리(훅 밖의 커밋·사용자 설정·다른 저장소), clippy 가 보는 기능 조합은 그 doc 이다.

**진실과 캐시를 이름으로 가르는 이유.** `Plane` 은 `Surface::Plane`, 즉 실체 이름이라 `Model` 에 `plane(h)` 같은 문을 달면
진실을 줄 것처럼 읽힌다. 캐시를 주는 문은 `_cache`(또는 캐시의 조각 이름 `_point`·`_curve`·`_name`)로 끝나 **「캐시의」라고
말하고 들어온다**. 조각을 돌려받은 타입의 메서드로 꺼내는 것(`vertex_cache(v).bound()`)은 `nacre_geom::Surface` 가
`distance()`·`normal_at()` 을 자기 메서드로 갖고 안에서 가르는 것과 같은 관용구다 — 구·원뿔이 와도 `Model` 의 문은 늘지 않는다.

성능은 이 명명의 고려사항이 아니다. 진실과 캐시의 문은 어느 철자든 «디버그 가드 + 배열 색인 하나 + 판별자 읽기»이고
할당도 복사도 없다(`#[inline]`) — 평면의 이름도 실현과 같은 캐시 칸에 있다.

# 층과 하위 체계

## 저장소와 동일성

Fornjot에서 그대로 가져오는 부분. 수정·삭제 없는 append-only `Vec` 기반 저장소와, 인덱스+타입만 가진 Handle.

```rust
pub struct Store<T> { items: Vec<T> }

impl<T> Store<T> {
    pub fn push(&mut self, item: T) -> Handle<T> { /* index를 Handle로 */ }
    pub fn get(&self, h: Handle<T>) -> &T { /* 항상 유효 — 삭제가 없으므로 */ }
}

/// 번호표일 뿐(u32 인덱스). T는 "어느 store를 가리키는가"의 타입 라벨.
/// derive를 쓰지 않는다 — `#[derive(Hash/Eq/Ord)]`는 T에 불필요한 바운드를
/// 자동 추가하고, 그러면 `HashMap<Handle<Edge>>`(「위상 층」 절의 Adjacency)가 `Edge: Hash`를
/// 요구해 컴파일 실패한다. 바운드는 **T와 무관해야** 하고, T가
/// 무엇이 되든 Handle 이 그 바운드를 강요하면 안 된다. index만 비교/해시하는
/// 수동 impl로 T 바운드를 끊는다. Copy 필수(정수 번호표), PhantomData<fn()->T>로
/// T와 무관하게 Send+Sync·공변성 확보.
/// (디버그 빌드에서는 두 타입 모두 아래의 store-id 필드를 하나 더 든다.)
pub struct Handle<T> { index: u32, _t: PhantomData<fn() -> T> }

impl<T> Clone for Handle<T> { fn clone(&self) -> Self { *self } }
impl<T> Copy for Handle<T> {}
impl<T> PartialEq for Handle<T> { fn eq(&self, o: &Self) -> bool { self.index == o.index } }
impl<T> Eq for Handle<T> {}
// Hash / PartialOrd / Ord / Debug 도 index만 보는 수동 impl.
```

동일성은 `Handle` 비교로 끝난다. "이 두 정점이 같은 점인가?"는 `h1 == h2`이며 부동소수점이 개입하지 않는다. 삭제가 없으므로 Handle은 영구히 유효하고, "지워진 객체를 가리키는 참조" 계열의 버그가 원천 차단된다. 불리언 등으로 객체가 소비되어도 항목을 지우지 않고, 결과 Solid가 새 항목들을 참조할 뿐이다.

주의: Handle의 유효성은 **자기 Model 안에서만** 성립한다. replay가 새 Model을 반환하므로 모델 두 개가 공존하는 순간이 실제로 생기고, A 모델의 Handle을 B에 쓰면 조용히 엉뚱한 객체가 나온다. 디버그 빌드에서 Handle 이 자기를 발급한 store 의 id 를 들고 `Store::get` 이 검사한다(`#[cfg(debug_assertions)]` — 릴리즈에서는 zero-cost로 제거). id 는 프로세스 로컬 `AtomicU64` 카운터로 `Store` 생성 때마다 발급한다(replay마다 새 값). 따라서 **Handle 자체는 직렬화하지 않는다** — 영속화 대상은 연산 로그(「연산 층」 절)이고, replay가 인덱스를 결정적으로 재생성한다.

**로그 속 핸들은 «인덱스 어휘»다.** `Operation` 의 변종은 `Stated` datum 하나를 빼고 전부 핸들을 싣는데(`Extrude` 의 프레임 평면 — 앞선 연산이 만든 면의 평면일 수 있다 · `DatumPlane` 의 정점/프레임 · `Boolean`·`Transform`·`Mirror`·`Copy`), 그 핸들이 가리키는 셀은 로그가 기록된 모델의 것이지 replay 가 짓고 있는 모델의 것이 아니다. 핸들에서 모델을 건너 살아남는 부분은 **인덱스뿐**이므로, replay 는 op 를 적용하기 직전 그 인덱스를 **자기 아레나의 핸들로 재고정(rebind)** 한다. 범위 밖이면 `LogHandleOutOfRange { cell, index }` — **존재만** 답하고, live 여부·합법성은 여전히 op 자신의 몫이다(`SolidNotLive` 등).

재고정은 **op 단위 just-in-time** 이다(op N 의 핸들은 op N−1 이 만든 셀을 가리키므로 사전 일괄 변환은 성립하지 않는다). 덕분에 원자성이 공짜로 따라온다: 재고정 실패는 그 op 가 아직 아무것도 push 하기 전에 일어나고, 실패한 replay 의 지역 모델은 통째로 버려진다.

**`apply` 는 재고정하지 않는다.** `apply` 의 모델은 호출자의 것이고 따라서 그 핸들도 호출자의 것이다. 거기서 재고정하면 외래 핸들을 **조용히 세탁**해 위 store-id 가드의 가치를 파괴한다. 재고정이 정당한 것은 「모델을 자기가 처음부터 짓는」 replay 하나뿐이다.

**유효 범위 — 오독하기 쉬운 자리.** 인덱스 어휘가 성립하는 것은 **로그를 처음부터 통째로 재생할 때**뿐이다(그때 인덱스는 재생이 스스로 만든 것이다). 로그 **중간을 수정한** 재생에는 성립하지 않는다 — 상류 수정이 하류 인덱스를 밀어낸다(v1 비목표; 「연산 층」 절의 계획 `OpRef` 참조). 그리고 이것은 「핸들을 직렬화해도 된다」는 뜻이 **아니다**: 직렬화 대상은 여전히 로그이고, 그 안의 참조가 인덱스일 뿐 핸들이 아니다(재생이 매번 새로 발급한다).

**자기완결성 전제.** 재고정은 「로그가 그 모델의 전체 이력이다」를 **대체하지 않고 요구한다**. 로그 밖에서 셀을 넣은 모델(예: 로그 없이 `apply` 를 부르는 테스트 픽스처, 심기용 `push_*` 문)은 replay 로 재현되지 않으며, 범위 검사는 인덱스가 **범위 안이면서 다른 셀**을 가리키는 경우를 잡지 못한다.

**거절은 live 모델에 원자적이지만 아레나 인덱스에는 아니다.** 거절에는 두 종류가 있다:

- **이른 거절** — 요청만 보고 판정(`ZeroDistance`, 창 밖, `Profile2d::check`, `NonPlanarFace`, `SolidNotLive`…). 아레나 Δ = 0.
- **늦은 거절** — 기하를 지어 봐야 알 수 있는 판정(`Boolean(_)`). 배열이 지은 셀이 **append-only 아레나에 남는다**.

live 모델은 어느 쪽이든 **거절 전 상태로 복원된다** — 커밋 후 거절은 없다. 부울에서 이것은 우연이 아니라 **구조**다: 조립의 피연산자 은퇴는 「결과가 받아들여진 뒤 한 자리」에서만 일어나고(`assemble_fuse_cut` → `reconstruct`), 그와 별도로 `boolean` 이 엔진의 **모든** `Err` 에 live set 스냅샷을 복원한다(`restore_live`; `coverage/rejects.rs` 와 `self_touch.rs` 가 단언한다). 그러나 아레나 길이는 되돌지 않으므로, **거절 뒤에도 기록을 이어가려면 모델을 로그로 다시 지어야 한다**(`model = replay(&log)`). 이 규율을 어긴 세션은 replay 가 재현할 수 없는 인덱스를 적게 되고, 결과는 **이름 붙은 거절이거나 발산**이며 — 패닉은 아니다(예: `LogHandleOutOfRange{Solid, 4}`).

이 절의 주장은 전부 `crates/nacre-ops/tests/invariants/replay.rs` 가 측정한다.

저장소는 도메인별로 나눈다.

**`Model` 의 필드 목록은 여기 그리지 않는다 — 「타입 구조」 부의 「최종 타입 — 진실」·「캐시」 절이 진실이다.** 진실/캐시 어휘를 두 곳에 그리면 **반드시** 한쪽이 상한다. 도메인 구분만 남긴다:

- **진실**(아레나, append-only): 곡면 · 정점 · 간선 · 면 · 셸 · 솔리드 · 모션(interned) · live 솔리드 목록
- **interning 표**(같은 진술이 두 핸들이 되지 않게 한다; 순회하지 않는다): 모션 · 곡면(열쇠는 셋 — 평면의 정준 이름 + 모션 · 이름 없는 평면의 진술 자체 · 원통의 진술 자체)
- **캐시**(핸들 인덱스 병렬, 버리고 재생 가능 — 정점 좌표는 정의에서, 정의가 못 내는 구성 수치는 로그에서): 정점 좌표 · 곡면 실현 · 간선 곡선 · 역방향 인덱스(「위상 층」 절)
- **순수 가속기**(진실도 캐시도 아님 — 비어 있어도 항상 옳다): 모션 사슬 접두의 고정밀 값 표

타입 이름·필드·공개 여부는 그 절들을 본다.

**Model은 "진실"만 담는다 — tess도 ops도 Model의 필드가 아니다.** 모델은 연산 로그의 재생 **결과**이므로 `ops: Vec<Operation>`은 Model을 *만들어내는 입력*이지 Model 안에 든 게 아니고, tessellation은 Model에서 *뽑아낸 파생 캐시*다(진실/캐시 분리). 게다가 레이어링상 Model은 `nacre-topo`에 사는데 `Tessellation`·`Operation`은 topo보다 위 크레이트라, Model에 필드로 넣으면 topo→tess/ops 순환이 된다. 그래서 op 로그와 tess 캐시는 상위 레이어(ops/파사드)에서 Model과 **나란히** 보관하고, "동일 로그·동일 cfg → 동일 모델 + 함께 재생성되는 tess"라는 커플링은 그 상위 번들이 책임진다.

(점 저장소는 두지 않는다 — 두 Vertex가 한 Point를 공유하는 상황은 설계상 존재하면 안 된다. `Vertex` 는 자기 **정의**만 들고, 좌표와 그것에 대해 증명된 것은 인덱스 병렬 캐시(`Model::vertex_cache`, 좌표 조각은 `Model::vertex_point`)에 산다.)

### 편집 연산의 supersede 의미론 — live 도달가능성 (partially persistent)

편집 연산(boolean·transform·mirror)이 기존 위상을 바꿀 때, append-only라 옛 셀을 **지우지 못한다**. 그래서 Model은 **`live_solids: Vec<Handle<Solid>>`(살아있는 솔리드 목록)를 진실로 보유**하고, 편집 연산은 옛 셀을 남긴 채 새 셀을 push한 뒤 **live 목록이 새 결과 Solid만 참조하도록 갱신**한다(소비된 입력 Solid는 목록에서 빠진다). **"살아있는 모델"의 정의 = live_solids에서 하향 참조로 도달 가능한 셀의 폐포(reachable closure).** supersede된 옛 셀은 아레나에 남되 어떤 live solid도 안 가리키므로 자동으로 "안 보인다".

**`Mirror` — 반사는 루프 역순으로 표현한다.** `Operation::Mirror { solid, axis, offset }`는 좌표평면 `axis = offset` 기준 반사다(**3개 방향 × 임의 위치**). 길이를 보존하고 **손잡이를 뒤집는다** — `Isometry`가 proper motion 전용이라 회전·이동의 어떤 조합으로도 못 만드는 유일한 강체 운동이다("−1배 스케일"은 미러와 스케일을 뭉개는 트릭이고 커널에 스케일은 없다). `transform`처럼 입력을 supersede하므로 원본을 남기려면 `copy`와 짝짓는다.

**축 정렬로 한정하는 이유는 회전과 같다** — 커널의 회전축도 X/Y/Z뿐이다. 임의 평면 반사 `x − 2((x−p)·n)n`은 무리수 법선 때문에 좌표가 근사가 되고, 그러면 모션 포리스트에 tol을 가진 반사 노드가 필요해진다. **임의 축 회전과 임의 평면 미러는 같은 항목**이며 그 설계가 들어올 때 함께 열린다(계획).

**구현의 핵심은 방향 대수다.** 반사는 det = −1이라 `R(a) × R(b) = −R(a × b)`이므로, 좌표만 반사하면 루프의 감김이 함의하는 법선이 뒤집혀 **솔리드가 안팎으로 뒤집힌다**. 그래서 모든 루프(outer + inner)를 **역순으로 감고**(`Loop::reversed` — 순서 역전과 `forward` 반전을 함께 하므로 모서리 짝이 유지된다), `Orientation` 플래그는 **건드리지 않는다**: 반사가 내적을 보존하므로 `sign(plane.normal · n_out)`이 그대로다. `nacre-geom`의 `Plane::mirrored`는 단위 법선의 한 성분 부호만 바꾸므로 비트 그대로 단위이고 다시 정규화하지 않는다. **`validate`의 그물은 둘이다**: 면 단위 orientation↔감김 불일치는 `FaceMisoriented`가 잡지만(반사는 루프 역순과 법선 반사가 함께 가서 이 검사에 불변), 플래그와 감김이 **함께** 뒤집힌 전-셸 반전은 면 단위 검사가 원리적으로 못 보므로, 평면 셸이면 셸의 부호 있는 부피(`ShellMisoriented`)가 잡고 곡면이 낀 셸은 불리언과 OCCT 오라클이 그물이다.

**반사도 모션이다 — 사슬에 들어간다.** `Motion::Mirror { axis, offset }`이 사슬 노드다. 기록 여부는 다른 모션과 **같은 규칙**을 따른다(「모션 사슬 — `MoveNode`」의 옮겨 적기 규칙) — `2c − x` 는 모든 유리수 `c` 에서 유리수라 새 솔리드의 반사는 진술로 옮겨 적힌다.

**반사는 비고유(`det = −1`)라 지름길에 패리티 규약이 필요하다.** exact 지름길들은 *"모션이 행렬식을 보존한다"* 를 근거로 **사전-모션 프레임**에서 답한다. 반사는 그것을 뒤집고, 게다가 실패가 보수적이지 않다 — 네 점이 같은 사슬을 지니면 행렬식이 **균일하게** 뒤집혀 확신을 가진 오답이 된다. 자리마다 인자를 유도하는 대신 **base 데이터의 손잡이를 한 번 맞춘다**: 사슬의 반사 수가 홀수면 base 데이터에 정확한 반사를 한 번 더 걸어 건넨다 — 점에는 `shared_base` 가, 평면(이름)에는 `name_stored_ints` 가. 그러면 base와 moved가 **고유 모션으로 관계**되어 어떤 행렬식 질문이든 등식이 그대로 선다. **base 평면의 부호는 그 전역 부호까지 센다** — 외적이 유사벡터라 반사는 삼각형의 회전을 뒤집고 평면의 기록은 안 뒤집으므로, 평면을 따로 반사하면 부호가 하나 남고 `plane_pair_dir_sign`이 바로 그 부호를 읽는다. 그래서 이름의 계수에는 `x` 유지·`y/z/d` 부정을 건다 — 반사(`x` 계수 부정)와 그 전역 부호(`det = −1`)를 함께 센 것이다. 예외는 `cmp_coord`의 지름길 하나뿐이다: 그것은 행렬식이 아니라 **축별 비교**이고 반사는 그 축을 뒤집으므로, 반사가 든 사슬에서는 끈다(보수적 미스).

**미구현(정직 거절):** 곡면(원/원통)의 반사 — 반사가 매개화 손잡이를 뒤집으므로 곡면 기하의 결정 사항이다. `MirrorNotPlanar` 로 거절한다.

**예외 하나 — `copy`는 supersede하지 않는다.** `Operation::Copy { solid }`는 솔리드의 독립 쌍둥이를 만들고 **원본을 live 목록에 남긴다** — 즉 **live 목록에 추가만 하는 유일한 연산**이다(다른 모든 편집 연산은 "옛 것 빼고 새 것 넣기"). 필요한 이유: `transform`·`boolean`이 입력을 live에서 빼므로, copy 없이는 같은 공구를 두 번 쓰거나 원본을 남긴 복사본을 두는 것이 불가능하고, 편의 레이어의 패턴·미러가 성립하지 않는다.

**구현:** `transform::transform_solid`(결정적 깊은 복제 walker: 곡면 → 정점 → 간선 → 면 → 셸 → 솔리드)를 **항등 아이소메트리**로 부르고 원본을 live에서 빼지 않는다 — `transform`과의 차이가 그 한 줄의 부재뿐이다. 위상 셀은 **반드시 복제**한다(두 live solid가 모서리를 공유하면 `edge_uses`가 면 4개 사용으로 읽어 manifold 검사가 깨진다). **곡면은 복제되지 않는다**: 항등 모션은 기록할 노드가 없고 진술이 비트 동일하므로, 같은 진술을 다시 push 하면 interning 이 **원본의 `Handle<Surface>`** 로 돌려준다(같은 진술은 두 핸들이 되지 않는다). 정점의 정의는 복제 후에도 살아남고, 좌표 캐시는 옮겨 받는 것이 아니라 새 정점의 정의를 실현해 채운다. 비-live 입력은 거절한다(소모된 핸들의 부활을 막고, 소비자 버그를 드러낸다). 이름은 **copy**(Handle의 clone은 얕은 참조 복사라는 반대 뜻).

**결과로, 모델을 소비하는 모든 코드는 store 전체가 아니라 live 도달가능 집합만 순회·카운트한다** — `validate`의 Euler/manifold, `Adjacency`(「위상 층」 절), `nacre-props`(부피·면적), `nacre-step` export, 오라클. 편집 중에는 아레나를 정리(compact)하지 않는다 — Handle 인덱스 안정성(=replay 결정성)을 유지하기 위해.

**대안 대비 (왜 도달가능성인가).** (a) *통째 복사*(편집마다 모든 셀을 새 독립 Solid로 복제) — 안 바뀐 셀까지 복사해 **정점 공유(single-reference)를 깨고** 메모리·오염을 부른다. 탈락. (b) *tombstone*(무효 Handle 집합을 손으로 관리) — 표시 누락 시 조용히 오염. 도달가능성은 **"아무도 안 가리키면 죽음"이라 누락 실수가 구조적으로 불가능**. 채택. 이는 partially persistent data structure(과거 버전은 읽기만, 최신만 수정)의 표준 CAD 적용이며, 이 절 위의 "소비돼도 안 지우고 결과 Solid가 새 항목 참조" 원칙의 구체화다.

**상용과 의도적으로 다른 선택 — 근거.** Parasolid·ACIS·OCCT 등 상용 커널은 40년간 "기하 불변 + **위상 가변**(Euler operator로 in-place 편집)"이 정설이다. nacre는 "기하 불변 + **위상도 불변**(append-only supersede)"으로 간다. 가변 위상은 dangling 참조·transient topology numbering·persistent naming 문제를 낳고 그걸 길들이는 데 40년 케이스워크가 필요한데 우리에겐 그 경험이 없다. 불변 위상은 그 버그 부류를 **구조적으로 제거**해 경험 없이도 안전하고, undo·결정론·fuzzing이 공짜이며, Rust 소유권 모델과도 정합적이다(가변 위상 그래프는 borrow checker와 싸워야 해 악명 높게 어렵다). 대가는 메모리인데 2020년대엔 작다. 즉 "상용을 이기는 방식"이 아니라 **"우리 제약(경험 없는 개인+AI 프로젝트·Rust·결정론 중시)에 최적화된 방식"**이다.

**성능·전제 (구현은 프로파일링 후).** 매 순회의 도달가능성 계산이 병목이면 generation/epoch 태그나 reference counting을 **파생 캐시**로 도입할 수 있다(진실=도달가능성, 이건 캐시 — Adjacency가 "진실=위상store, 파생=역인덱스"인 것과 같은 층 분리). ref-counting은 하향 참조가 acyclic이면 죽음 판정에 바로 쓸 수 있다. **전제**: 도달가능성·ref-counting 둘 다 "Solid→Shell→Face→Loop→Edge→Vertex 부모→자식 단방향, 순환 없음"을 요구한다 — 「위상 층」 절의 위상 참조가 이를 만족하고(모든 셀이 하향 Handle만 보유), 역방향 인덱스 `Adjacency`는 진실이 아니라 모델 밖의 파생 인덱스라 무관하다.

**계획(v2) — 세션 중 메모리 관리: compact보다 재구축(rebuild-from-log) 우선.** 저장·세션·undo 형식을 설계할 때의 방향이며, 그 전에는 구현·상세 확정을 하지 않는다("경계는 지금, 내용은 그때"). 로그가 진실인 설계에선 compact가 주력일 필요가 없다. 배경: append-only라 편집·브랜치를 반복하면 도달불가(dead) 셀이 누적된다.
- **정리는 compact보다 재구축이 낫다.** compact(dead 셀 물리 제거 + 인덱스 재배치)는 "현재 브랜치의 과거 단계 셀"도 도달불가라 함께 지워 **undo를 죽인다**. 재구축(현재 live 브랜치의 로그를 0부터 재생해 새 조밀 아레나 생성)은 로그에 모든 과거 단계가 있어 **undo 히스토리째로 되살아난다**. 로그가 진실이라 결정론적으로 동일 결과 보장. "dead를 골라 지우기"보다 "live를 로그로 새로 짓기"가 더 단순·안전하며 undo까지 지킨다.
- **저장 시 compact 불필요.** 스냅샷 직렬화 때 dead를 파일에 안 담으려 아레나를 compact할 이유가 없다 — 직렬화 시점에 **live 도달가능 폐포만 순회해 쓰면 된다**(아레나 불변 → 저장 후 계속 작업). "validate·adjacency·step은 live 폐포만 순회"·"STEP은 live만 export"와 동일 패턴. STEP 출력·로그 저장 모두 아레나를 안 건드리므로 어떤 저장도 compact를 트리거하지 않는다.
- **compact가 유효한 단 하나의 자리 = 메모리 위기 시 비상 회수.** 아주 길고 무거운 세션에서 dead가 실제로 메모리를 압박하면, compact로 dead·과거단계를 회수하고 현재 상태 스냅샷 + 로그를 유지한 채 가볍게 이어간다. undo는 일시 상실되나 로그가 살아 **재구축으로 복구 가능**(영구 손실 아님). 즉 compact = "undo 즉시성을 메모리와 맞바꾸는, 재구축으로 되돌릴 수 있는 스위치". 전제: **로그는 절대 버리지 않는다**(compact는 상태만 버리고 로그는 지킴).
- **역할 정리:** 재구축 = 주력 정리(undo 유지), live 폐포 필터링 = 저장(아레나 불변), compact = 비상 메모리 회수(undo는 재구축으로 복구). 셋은 대립이 아니라 안전망 관계.

**계획(파라메트릭 편집 v2).** 편집 대상 참조를 transient Handle 인덱스로 두면 kernel numbering에 종속돼 상류 수정에 깨진다(「연산 층」 절의 topological naming). 정석은 **저장된 기하 cue(3D 참조점 등)로 tolerance 내 재탐색**해 참조를 의미로 resolve하는 것(HistCAD 계열). 「연산 층」 절의 계획 `OpRef`는 구현하지 않되, 참조를 "순수 인덱스"가 아니라 "기하 cue 포함"으로 확장할 자리를 로그 포맷에 남긴다.

## 해석 기하 층 (`nacre-geom`)

이 크레이트는 **해석적**이지만 **정확하지 않다** — 계수가 f64 다. 곡면의 진실은
`nacre_topo::Surface` 이고(평면은 `PlanePoints` — **유리수 세 점**(`Known`) 또는 **모델 정점 셋을
지난다**(`Through`, datum 평면 — 좌표로는 그 평면을 말할 수 없다) — 원통은 `CylinderDef`, 둘 다 모션
슬롯을 곁에 둔다), 여기 있는 `Surface`·`Curve` 는 그 **실현**이다. 산술이 닫힌 형태라는 것과 값이
정확하다는 것은 다른 말이다. 모든 포함 질의는 호출자의 epsilon 을 받고, 이 층은 tolerance 를 저장하지
않는다.

**이 크레이트는 `nacre-store` 에 의존하지 않는다 — `Handle` 을 들지 않는다.** 핸들은 아레나 항목의
이름이고 아레나가 드는 것은 곡면의 진실(`nacre_topo::Surface`, 모션 핸들을 든다)이라, topo 아래의
타입은 그것을 이름 지을 수 없다.

```rust
pub enum Surface {
    Plane(Plane),
    Cylinder(Cylinder),
    // 계획: `Nurbs(NurbsSurface)`(평가기 `NurbsSurface` 는 이미 있다)·`Sphere` 는 생산자가
    // 붙을 때 변종이 된다. Cone·Torus·회전면·스윕면은 마일스톤 따라. 그래서 이 타입은
    // 처음부터 `Copy` 가 아니다 — 올 `Nurbs` 변종이 힙 제어점을 소유한다.
}

pub enum Curve {
    Line(Line),
    /// 림의 전체 carrier — 평면이 원통과 만나는 곡선, 원통의 θ 로 매개화: 중심(축∩평면)/축/
    /// **ref_dir**/반경 + 기울기 `tilt: Option<[t₀, t₁]>`(축 방향 높이 `t₀ cos θ + t₁ sin θ`;
    /// `None` = 원). `ref_dir` 이 θ=0 앵커이고, 원통의 seam 은 **+ref_dir** 에 고정된다.
    /// 원인지는 정확한 축 관계가 정하고 공식(점·활꼴·모멘트·상자·STEP 의 수 `written`)은 모두
    /// 메서드다 — 원의 갈래는 원의 식을 그대로 돈다. 호 = 이 담체 + 서로 다른 두 끝점 정점,
    /// 온전한 원 = 같은 정점이 양 끝(`[v, v]`, 위상 층 절의 rim 규칙) — Line 과 같은
    /// carrier-vs-trim. brep-to-step `Curve::Circle`·`Curve::Ellipse`(끝점 동일 여부로 온전/호 구분)과 정합.
    Ellipse(Ellipse),
    // 계획: `Nurbs(NurbsCurve)` 도 생산자가 붙을 때 온다(평가기 `NurbsCurve` 는 이미 있다).
}
```

`Curve` 는 간선 기하의 **캐시 쪽**이다 — 간선이 기록하는 담체와 끝점이 진실이고, 곡선은 거기서
유도된다(위상 층 절).

**정의에서 좌표를 실현하는 일은 `nacre_ops::realize_vertex` 가 한다** — 정점의 **정의**를 받아 호출자가 고른
정밀도에서 좌표를 만들고, **끝에서 딱 한 번 반올림한다**. 길은 `Vertex` 변종이 아니라 **값**이 고른다:
세 평면의 만남은 정확한 비율이라 긴 나눗셈 하나로 모든 자릿수가 점 자신의 자릿수이고, 모션이나 근호가
닿는 것은 요구 비트에서 접근해 **구간이 자릿수를 결정할 때만** 찍는다. **실패는 삼키지 않고 위로
올린다** — 특히 접선·퇴화처럼 «판정이 아니라 정책»이 필요한 구역은 이름을 달아 보고한다. 조용히 넘기는
순간 디버깅 불가능한 커널이 된다.

정밀도 분업 원칙: **판정(부호)에는 적응 정밀 술어**(`geometry-predicates` — 쉬운 케이스 f64, 아슬아슬할 때만 확장, 부호는 항상 정확), **반복 구성(좌표)에는 «확장» 부동소수점**. 임의 정밀 유리수(rug/malachite)는 술어에는 완벽하지만 Newton 반복에 넣으면 비트 길이가 반복마다 폭발하므로 구성에는 쓰지 않는다. 이 투자는 평균이 아니라 꼬리를 산다 — 호출 빈도가 낮은 악조건 케이스만 정확히 개선되고, tolerance가 "상수"가 아니라 "실측 보증값"이 된다. 경계는 「반복이면 확장 산술, 아니면 유리수」이지 「불리언이면 f64」가 아니다.

`realize_vertex` 가 **호출자가 고른** 정밀도에서 사다리를 오른다(`Precision` 은 *"always stated, never defaulted"*). 판정 층의 정밀도는 모델이 정하고(정밀도 인프라 절), 실현은 호출자가 정한다 — **둘은 다른 손잡이다.**

**정점의 캐시는 태어날 때부터 실현값이다.** 연산이 정점을 push 하는 자리 다섯(불리언 셋·압출·이동)은 전부 ops 의 깔때기 `push_vertex_realized(model, def, figure)` 를 지난다: 정의를 사다리의 첫 두 단(128비트, 못 정하면 256비트)에서 실현해 `PointCache::Realized { coord, bound }` 로 넣는다. 못 넣었을 때 서는 것은 «구성의 폴백»이 아니라 «이름 붙은 이유»다: 깔때기는 `realize_vertex_cache` 의 `Err` 를 망라적으로 갈라 `Ceiling`(첫 단의 비트가 모자랐거나, 재생할 이력이 **비용 한계** `CHEAP_REPLAY_COST_CAP` 을 넘어 아예 안 걸었거나 — 둘 다 «더 물으면 답이 있다») 과 `Unrealized`(도로가 없거나, 정의가 점을 이름하지 않는다)로 나누고, 좌표만 구성 자리의 값에서 가져온다. 그래서 «실현 통로는 하나»가 구조다 — 표시·tess·STEP·물성·validate 가 읽는 f64 는 realize 의 메모다. **쓰기 문은 둘**이다: 짓는 `push_vertex` 와, `Ceiling` 만 올리는 `Model::refine_vertex_cache` — 그것을 부르는 것은 내보내기 직전의 «비싼 문» `nacre_ops::refine_caches` 하나다(정점·곡면·간선을 차례로 올리고 남은 것을 이유별로 보고한다; 로그 끝의 문이다 — 정점 캐시는 뒤 연산의 진실을 정하지 않지만, 평면 캐시는 뒤 불리언의 seam 표가 `Ceiling` 정점에 쓰는 구성 수치라 그 표의 거절이 달라질 수 있다). 내보내기는 `&Model` 을 받으므로 그 문을 부르지 않는다 — 부르는 것은 호출자다. 계약은 census 가 행마다 단언한다: 답 ⇔ 변종, 값은 비트 동일. 거절 인구는 이름을 든다(`NoMeet` 등).

판정 술어는 다시 둘로 나뉜다: **명시 좌표점(f64 격자 위)은 direct 술어**(orient3d 등, 좌표를 직접 받음)로, **교차로 정의된 점은 indirect 술어**(implicit point = "어느 원시 요소들의 교차인지"라는 정의를 받아 좌표를 만들지 않고 부호를 정확히 계산 — Attene 2020)로 판정한다. 이유: 교차점을 f64 좌표로 만드는 순간 오차가 끼고, "정확한 direct 술어 × 부정확한 입력 = 부정확한 답"이 되기 때문. indirect 술어는 그 구멍을 닫는다(정의째 받으므로 부정확한 중간 좌표가 없음). 둘은 같은 확장 산술 바닥(`geometry-predicates`)을 공유하며, indirect 층은 `nacre-predicates`가 그 위에 쌓는다 — `indirect_cmp_coord`·`indirect_plane_side` 가 그것이다. indirect predicate가 작동하려면 점이 "정의를 보유"해야 하므로, **점이 정의를 갖는다**는 결정(위상 층 절)이 그 전제다 — 판정이 필요하면 indirect 술어, 좌표가 필요하면 **정의에서 실현한다**(위 `realize_vertex`; 둘은 대체가 아니라 역할 분담).

**STEP 연결(계획).** STEP 의 `surface_curve`/`intersection_curve`(3D curve + 곡면별 pcurve + master 지정)에 대응하는 것은 `Curve` 의 변종 하나가 아니다 — 그런 변종은 `Handle<Surface>` 를 들어야 하므로 이 크레이트가 표현할 수 없다. 그 내용은 진실/캐시 선을 따라 갈린다: **«어느 두 곡면인가»는 이미 간선의 진실**(`Edge::surfaces`, 담체 두 면)이고, **근사 곡선과 그 오차는 간선의 캐시**(`EdgeCache`, 지금은 `curve` 하나)다. 계획: 마친(marched) 곡면 교차를 실제로 만들 때 `EdgeCache` 를 `{curve, err}` 로 키우고, pcurve 가 갈 자리도 거기다. STEP 쪽 엔티티 묶음의 보존·복원은 그때 검증 대상이 된다(`nacre-step` 의 커버리지는 **평면 + 원통**이고 그 엔티티는 0건).

## 위상 층 (`nacre-topo`)

b-rep 의 셀이 사는 층이다. **모든 셀은 정확 기하를 `Handle` 로만 참조하고, 위상은 좌표를 보지 않는다.**
셀은 모든 필드가 핸들이거나 플래그일 때만 `Eq`/`Hash` 를 유도한다 — 셀이 f64 를 들 수 있는가에 대한
규칙이고, 좌표를 드는 셀은 없다.

진실 타입의 상세 모양 — 곡면(`Surface`·`PlanePoints`·`CylinderDef`)·정점(`Vertex` 의 변종들)·간선,
그리고 **모션 삼종**(`Motion`·`FramePlacement`·`MotionNode`) — 은 「타입 구조」 부가 그린다. 여기는 그
위에 선 **규칙**과, **규칙의 주어가 되는 셀**을 적는다. 그래서 `Edge` 는 거기서도 여기서도 그려진다:
rim 이 `[v, v]` 라는 것과 담체를 유도할 수 없다는 것은 **필드 값에 대한 문장**이라, 필드를 안 그리면
«무엇이 둘인지»를 확인할 데가 없다. 모션 삼종은 이 층에 살고 아레나가 **진실로** 들지만
(`motions: Store<MotionNode>`) 규칙의 주어가 아니라 이름만 짓고 넘긴다 — 아래 「서피스는 자기
leaf 에서 이어진다」가 그 노드 위에 선다.

**`Model` 은 진실과 진실에서 유도되는 것을 든다**: 정확 기하의 진실(`surfaces`·`motions`), 위상 store 다섯
(`vertices`·`edges`·`faces`·`shells`·`solids`)과 `live_solids`, 그리고 진실에서 유도되는 곁표들
— 인덱스 병렬 캐시 셋(`surface_cache`·`edge_cache`·`vertex_cache`), interning 표,
가속기 `prefix_hp`. tess 와 ops 는 들지 않는다.

**정확 기하는 정점만이 아니라 «면»에도 적용된다.** 비-90° 회전은 평면의 계수를 **무리수로** 만든다 — f64 로 저장된 `Plane{origin, normal, raw}` 은 그 순간 진실이 아니라 **반올림된 상**이고, 그 사실을 적을 곳이 없으면 커널은 기본값으로 "정확하다"고 답하게 된다. 그래서 **곡면의 진실이 자기 `motion` 을 든다**: 유리수 진술은 `motion` 이 이름 짓는 틀에서 말하고, 모션 이력이 그것을 세계로 나른다(`None` = 세계).

**정점은 자기 정의이고 그것뿐이다.** `Vertex` 는 정의의 enum 이다(`ThreePlane`·`OnSeam`·`Pierce` —
변종마다 자기 진실을 말하고 불변도 변종별이다). 좌표와 «그 좌표에 대해 증명된 것»은 정점 store 와
인덱스 병렬인 점 캐시가 든다(`PointCache` — `Realized { coord, bound }`·`Ceiling`·`Unrealized`;
`Model::vertex_cache`/`Model::vertex_point` 로 읽고 `Model::push_vertex` 가 채운다 — 정의가 먼저,
좌표가 둘째). 캐시의 변종은 캐시가 좌표에 대해 **아는 것**을 말하지 정점이 어디서 왔는지를 말하지
않는다. 간선의 캐시와 달리 버리고 다시 쓰는 문이 없다: 정점은 push 될 때 실현되고, 뒤에 일어날 수 있는
것은 `Ceiling` → `Realized` 의 **정련**뿐이다. 다시 쓸 것이 없어서다 — `Realized` 는 이미 정의가 다시 답하는 값이고
(로그가 끝난 뒤 아레나의 모든 정점에 물어도 비트 동일 — `vertex_caches_regenerate_from_their_definitions_at_the_end_of_a_log`;
census 는 결과 정점마다 묻는다), `Ceiling`·`Unrealized` 의 좌표는 구성 자리의 수치라 정의가 아니라 **로그**를
재생해야 다시 나온다(재생 잠금이 좌표 비트·변종·경계를 비교한다).

**간선의 진실은 «담체 둘 + 끝점 둘»이다.** 실현된 곡선은 필드가 아니라 곁의 캐시(`Model::edge_cache`,
`EdgeCache`)다 — 담체와 끝점이 그것을 정한다(`Model::derive_edge_curve`). `Model::push_edge` 가 캐시를
함께 채우고, `Model::rebuild_edge_cache` 가 통째로 버리고 재생한다. 값 가운데 topo 가 진실에서 못 내는 것 —
세계 이름 없는 평면이 낀 직선의 방향, 세계 진술 없는 원통의 림 중심(그 사슬을 재생해야 한다) — 은 유도가 자기 길로 답하지 못할 때만
push 하는 쪽에 묻는다(`EdgeGiven`; ops 깔때기 `push_edge_realized` 가 실현해 답한다). «진실이 여기서
답하나»는 그래서 topo 한 곳에서만 정해진다. **곡선의 종류도 진실이 정한다** —
평면 × 원통 간선이 원(평면이 축을 가로지름)인지 룰링(축을 따름)인지는 평면 이름의 법선과 원통 `def.dir()` 의
정확한 관계(`nacre_exact::AxisRelation`, `Model::plane_cylinder_relation`)이고, 비스듬하면 같은 원을 그 평면에
기울인 타원이다(원통의 θ 그대로 — 평면 캐시의 법선이 기울기를 준다; 수가 기울기를 말하지 못하면 `EdgeDecline::Degenerate`).
캐시가 드는 것은 그 곡선의 값(중심·프레임·기울기)뿐이라 재생이 종류까지 같고, 간선 곡선의
종류를 읽는 판정은 진실을 읽는다. 관계는 두 방향의 물음이고 강체 운동이 보존하므로 두 진술을 세계까지 옮길
필요가 없다: 각 방향을 자기 사슬의 선형부로 정확히 옮길 수 있는 만큼 옮기고(`LineDir` — 이동은 그대로, 거울·
사분각은 치환, 회전은 자기 축을 그대로, 프레임은 `ẑ` 를 그 평면의 법선으로), 남은 운동이 **값으로** 같거나 한쪽이
더 갔으면 그만큼 되돌려 한 프레임에서 묻는다. 원통 인구 게이트(`cyl_gate`)는 같은 규칙(`nacre_exact::axis_relation`)을
부르고, 클래스가 원을 싣는가는 그 `Across` 팔 하나(`class_carries_circle`)로 묻는다.

```rust
pub struct Edge {
    /// **이 간선이 경계 짓는 두 면의 곡면** — 담체.
    pub surfaces: [Handle<Surface>; 2],
    /// 끝점 정점 둘 — 경계. **모든 간선이 둘 다 든다**.
    /// **rim 규칙**: 닫힌 솔리드의 원형 rim 은 seam 정점을 써 `[v, v]`(start == end)다
    /// — 그래야 b-rep 이 유효 CW-복합체로 남아 validate 의 오일러-푸앵카레를 통과한다
    /// (실린더 rim 이 대표). 끝점 없는 독립 전체 원은 **이 타입이 표현하지 않는다**
    /// — 와이어프레임/열린 면 요소는 v1 비목표다.
    pub vertices: [Handle<Vertex>; 2],
}
```

**담체는 «끝점의 면집합»에서 유도할 수 없다.** 한 점에 평면 넷이 모이면 두 끝점의 삼중 교집합이
그 간선이 **타지 않는** 평면을 이름 지을 수 있고, **조용히** 그런다(4평면 동시성). 그래서 **생산자마다
자기가 짓고 있는 면에서 담체를 진술하고, 절대 유도하지 않는다.** 「선을 담기만 하는 평면」도 담체가
아니다 — 세 평면이 한 **선**을 공유할 때 셋째 평면은 그 간선을 담지만 경계 짓지 않는다. 담체는 «인접»의
답, 즉 **그 간선을 쓰는 두 면**이다. 저장은 **핸들 오름차순**(집합이지 순서가 아니다 — `Edge::carrier_pair`
가 그 한 철자다)이고, 두 곡면은 둘이다 — 어느 간선도 한 곡면을 자기 자신과 가르지 않는다(원통 옆면은 림 두
루프로 닫히고 seam 간선이 없다; `push_edge` 는 한 원통 둘의 쌍을 `TwoCylinders` 로 거절한다). 어기면
`validate` 가 `EdgeCarrierMismatch` 로 잡는다(self-pair 는 곡면 종류와 무관하게 위반이다).

**출처가 기록되지 않은 표면은 «타입이» 표현 불가능하게 만든다.** 곡면의 진실은 아레나 항목 그 자체이고
(`Model::surface`), 진실 없는 곡면은 **핸들이 이름 지을 수 없다** — 그래서 「정의 없는 곡면」이라는
위반 부류가 존재하지 않는다. 그렇다고 `validate` 검사가 없는 것은 아니다: `VertexCarrierMismatch`·
`EdgeCarrierMismatch`·`VertexOffDefinition` 은 돈다. 타입으로 승격한 것은 표면의 **정의 유무**뿐이다.

**곡면을 미는 문은 여럿이되 «하나의 비공개 깔때기»로 모인다** — `push_plane`·`push_plane_through`·
`push_cylinder`·`push_plane_unregistered` 가 전부 사설 `push_plane_raw`·`push_cylinder_raw` 를 지나고,
그 자리에서 **진실과 f64 캐시가 한 동작으로 함께 들어간다**(*"the truth and its cache enter together or not at all"*). 둘이
떨어질 수 없게 하는 것이 요점이지 문의 개수가 아니다. `surfaces` store 는 비공개이고 통째 순회자는
일부러 없다(아레나는 supersede 된 곡면을 남긴다).

**정의가 진실이고 좌표는 캐시다 — validate 가 그 둘을 맞댄다.** `validate` 는 **모든 살아 있는 정점**이
자기 정의가 주장하는 모든 곡면 위에 앉아 있는지 검사한다(`VertexOffDefinition`). 출신별로 다른 규칙이
아니다: 물음은 하나이고, tol 도 하나 — 모든 정점에 상수 `EPS_CONSTRUCTED`(1e-9)를 쓴다. 캐시의
`bound` 는 여기에 쓰지 않는다: 그것은 **좌표**가 진실에서 얼마나 먼가를 말하고, 이 검사는 캐시된 점이
**담체**에서 얼마나 먼가를 묻는다. 비용을 적어 둔다: 코퍼스의 담체 거리 최대는 **1.07e-14** 라, 이
상수는 그 인구를 다섯 자릿수 느슨하게 본다.

**`Reachable` 은 surfaces 를 추적하지 않는다** — 그래서 곡면 순회는 **면을 통해** 돈다. (시드된 세계
평면 셋은 일부러 면이 없다.)

**서피스는 솔리드의 leaf가 아니라 자기 leaf에서 이어진다.** 한 솔리드에 회전 이력이 하나라는 전제는 거짓이다 — 서로 다른 각도로 회전한 두 피연산자의 불리언 결과는 벽마다 다른 이력을 갖는다. 그리고 불리언 결과의 정점은 자기 모션을 말하지 않으므로(모션은 **면**의 것이고 면이 스스로 기록한다) 정점에게 "이 솔리드는 어느 회전에 있나"를 물으면 답이 없고, 그렇게 물은 다음 회전은 **뿌리부터 다시 시작**해 회전 이전 witness에 두 번째 회전만 태우게 된다(존재하지 않는 평면). 그래서 `transform`은 변환마다 노드 하나가 아니라 **서로 다른 parent leaf마다 노드 하나**를 만든다.

**모델이 들고 다니는 tol 은 없다.** 구성된 점이든 불리언이 발견한 점이든 정점은 태어날 때 정의에서 실현되어 `Realized` 경계를 들고, 저장된 잔차는 없다. 그렇다고 **tolerance 가 사라지는 것은 아니다** — 검사는 그대로 돌고(위 `EPS_CONSTRUCTED`), 사라지는 것은 **모델이 들고 다니는 tol** 이지 비교가 아니다: Fornjot 의 "무-tolerance 단순함"은 «저장된 tol 이 없다»는 뜻으로만 재현된다.

**정의가 진실 → indirect predicate·봉합 우회.** 점이 좌표가 아니라 정의(implicit point)를 진실로 삼는다는 것은, 부호 판정을 좌표를 만들지 않고 indirect 술어로 정확히 하기 위한 전제다(해석 기하 층 절; Attene 2020). 이 결정이 실제로 해소하는 것은 과장 없이 셋: (1) **조합적 결정의 불일치**(orient 부호가 여기선 +, 저기선 −로 모순돼 위상이 깨지는 것 — 불리언 실패의 주원인), (2) "정확한 술어에 부정확한 입력" 구멍, (3) 판정에서의 tolerance 튜닝 취약성(부호 판정엔 tolerance 값 자체가 안 쓰임). **봉합 문제도 대부분 해소된다 — 단 «한 불리언이 민팅하는 코너»에 대해서다**: 세 평면이 만나는 꼭짓점을 하나의 `Vertex`로 만들고 정의를 "P₁∩P₂∩P₃"로 두면, 세 엣지가 그 하나의 Vertex Handle을 공유하므로 "세 곡선이 한 점에서 만나는가?"가 좌표 비교가 아니라 Handle 비교(`h==h==h`)가 되고 — 봉합의 본체(어긋난 세 후보를 tolerance로 화해)가 사라진다. 배열의 병합도 같은 말을 한다: `merge_coincident` 의 키는 **(벽, 정렬된 끝점 «이름» 쌍)** — 둘 다 정준으로 접은 뒤다 — 이지 **좌표가 아니다**(평면이 동시적이면 두 생산자가 **벽도 끝점 이름도** 다르게 적을 수 있어 **두 접기가 다 있어야** 키가 같아진다). 끝점 이름은 세 평면만이 아니다 — 관통점(`NodeId::Pierce`)도 이름이다. **연산을 가로지르는 동일성은 별개 물음이다** — 정점은 interning 하지 않는다(인접 면 집합이 자라므로 저장 시점에 키를 못 잡는다). 좌표를 f64로 뽑을 때도 하나의 Vertex 좌표만 뽑으므로 세 어긋난 좌표가 안 생긴다. 잔여물은 "그 하나의 f64 좌표가 세 평면 위에 정확히 안 놓인다"는 반 ulp 오차뿐(f64 저장의 바닥, 대부분 무해, 극단적으로 조밀한 형상에서만 문제). **해소 못 하는 것**: 접선·퇴화 교차의 위상 정책(부호 0일 때 "접하냐/스치냐"는 **판정이 아니라 정책**이고, 그 구역은 삼키지 않고 위로 올린다), 좌표 정밀도 자체, 공간 인덱스·SSI·성능 등 다른 축. 일반 곡면으로의 적용도는 미결이다.

half-edge는 교과서 구조를 따르되 Fornjot 신설계처럼 공유 `Edge` + 방향 참조로 둔다:

```rust
/// `forward == true` 는 start→end(`vertices[0]` → `vertices[1]`).
pub struct HalfEdge { pub edge: Handle<Edge>, pub forward: bool }
pub struct Loop { pub half_edges: Vec<HalfEdge> }

/// 2-셀 — 곡면 조각. 외곽 루프 하나 + 구멍 루프 N개.
pub struct Face {
    pub surface: Handle<Surface>,
    pub outer: Loop,
    pub inner: Vec<Loop>,
    pub orientation: Orientation,
}

/// 닫힌 면 집합 — 하나의 경계 곡면을 이룬다.
pub struct Shell { pub faces: Vec<Handle<Face>> }

/// 솔리드 = 바깥 껍질 + 내부 공동(cavity)들.
pub struct Solid {
    pub outer: Handle<Shell>,
    pub cavities: Vec<Handle<Shell>>,
}

/// 역방향 인덱스 — 원본이 아니라 파생물 (재스캔으로 언제든 재구성 가능).
/// 검증("모든 엣지는 정확히 두 면에서 반대 방향으로 사용")과 인접 탐색이 이걸 쓴다.
pub struct Adjacency {
    pub edge_uses: HashMap<Handle<Edge>, Vec<(Handle<Face>, bool)>>,
    pub vertex_edges: HashMap<Handle<Vertex>, Vec<Handle<Edge>>>,
}
```

**닫힌 rim 의 `forward` 는 끝점이 정하지 못한다** — 시작과 끝이 같은 정점이기 때문이다. 거기서
`forward` 는 *곡선 자신의 매개화를 따라*, 즉 원의 법선(= 원통의 축 방향)에 대해 CCW 로 읽는다. 읽을 수
있는 유일한 방식이고 모든 생산자가 그렇게 쓴다(아래 캡의 rim 은 `false` 이고 그 면은 `−axis` 를 바깥으로
말한다). `validate` 의 면 방향 검사가 생산자를 거기에 묶는다.

(`Adjacency` 의 값은 평범한 `Vec` 이다 — 다양체 케이스를 힙 없이 담는 `SmallVec<[_; 2]>`/`[_; 4]` 는 의존성 하나를 아끼려 쓰지 않았고, 인라인 저장 이득이 프로파일링으로 정당화되면 도입한다. 정확성과 무관한 최적화다.)

**`validate` 가 지키는 위상 불변.** 참조 무결성이 먼저 돌고(댕글링 핸들이 있으면 거기서 멈춘다), 그
뒤는 전부 **live 도달가능 셀** 위에서, 스스로 지은 `Adjacency` 로 돈다:

- **루프 닫힘** — 각 루프의 half-edge 가 끝→시작으로 이어진다(`OpenLoop`).
- **다양체** — 모든 간선은 정확히 두 번, 반대 방향으로 쓰인다(`NonManifoldEdge`·`NonOpposedEdge`);
  정점 둘레도 다양체다(`NonManifoldVertex`).
- **담체 일치** — 정점의 정의가 이름 짓는 곡면(`VertexCarrierMismatch`)과 간선이 진술한 담체 쌍
  (`EdgeCarrierMismatch`)이 그것을 실제로 쓰는 면들과 맞는다.
- **방향** — 면의 감김과 진술한 바깥 법선이 맞고(`FaceMisoriented`), 껍질의 부호 있는 부피가 제 역할의 부호다 —
  바깥 껍질은 양, 공동은 음(`ShellMisoriented`).
- **기하 접속** — 캐시된 점이 자기 정의의 곡면(`VertexOffDefinition`)·간선의 곡선(`VertexOffCurve`)·
  면의 곡면(`VertexOffSurface`) 위에 tol 안으로 앉아 있다. 원통은 진실과 f64 캐시가 맞는지도 본다
  (`CylinderTruthCacheMismatch`) — 세계에 진술된 원통의 캐시는 문이 그 진술에서 실현하므로, 진술을 생산자가 따로
  옮긴 값과 대조하는 일은 문(`push_cylinder_raw` 의 `debug_assert`)이 덮어쓰기 전에 한다.
- **오일러-푸앵카레** — `V − E + F = 2(S − G) + L_i`(`L_i` = 면의 내부 루프 수, `S` = 껍질 수,
  `G` = genus). **기준은 `V−E+F=2` 가 아니다**: `validate` 는 **χ = V − E + F − L_i** 를 세어 **짝수인지**
  (`EulerParity`)와 **G = S − χ/2 ≥ 0** 인지(`NegativeGenus`)를 본다. 뚫린 솔리드는 genus 1 이라 χ = 0
  이고, 커널이 실제로 그것을 만든다. rim 이 `[v, v]` 인 것은 이 셈이 맞게 하기 위해서다.

`Adjacency`는 진실이 아니라 파생 인덱스다 — 위상 store들이 진실이고, 인덱스는 쓰는 쪽이 필요할 때 **통째로** 짓는다 — `validate` 가 이 `Adjacency` 를 짓고, 불리언의 결과 검사와 tess 는 제 간선 표를 든다(`edge_uses`/`vertex_edges` 를 증분으로 건드리는 코드 0줄). 모델은 인덱스를 들지 않는다(「가지 말 것」). 증분 갱신은 **필요가 증명될 때** 짓는다.

**Adjacency는 store 전체가 아니라 live 도달가능 셀만 인덱싱한다(supersede 의미론).** `rebuild`는 `model.faces`/`model.edges` 전체가 아니라 `live_solids`에서 도달 가능한 면·엣지만 순회한다 — 그래야 supersede된 옛 면의 엣지가 `edge_uses`를 오염시켜 manifold 검사(엣지 정확히 2회)를 깨뜨리지 않는다. 위상 참조가 하향 단방향·비순환이라 도달가능성 순회는 유한·안전하다. `Adjacency`는 **위상 셀의** 역방향 인덱스이고 캐시이므로 도달가능성 진실을 침해하지 않는다. interning 표 넷(`motion_ids`·`interned` — 키 → `Handle`)은 «값으로 핸들을 찾는» 표라 성격이 다르고, **`Adjacency` 와 달리 도달가능성으로 걸러지지 않는다.** 어느 것도 순회하지 않는다(`HashMap` 의 순서가 결과에 닿으면 안 된다).

**곁표에는 «진실도 캐시도 아닌» 셋째 부류가 하나 있다 — 가속기.** `Model::prefix_hp` 는 모션 체인이 이미 접은 **접두사**의 고정밀 값을 들고, 다음 모션이 base 부터 다시 접는 대신 거기서 이어 접는다(이것이 없으면 깊이 n 인 솔리드를 짓는 데 1+2+…+n 이 든다). **캐시와 다른 점**: 캐시는 진실의 메모라 «있으면 읽는다»가 계약이지만, 이것은 비어 있어도 답이 **비트까지 같다** — 없으면 정의에서 다시 접을 뿐이다. **interning 표와 다른 점**: 핸들 정체를 지지 않으므로 **언제든 비울 수 있다**(시험 빌드의 `clear_prefix_hp` 가 그 문이다). 열쇠가 정의 `(base, leaf, prec)` 라 **낡을 수 없고**(무효화 규칙이 없다), **적중이 곧 축출**이라 용량이 이력 길이가 아니라 **살아 있는 꼭짓점 수**에 묶인다. 값은 「접은 체인 노드 수 + 점」이다 — `Motion::Frame` 하나가 여러 노드로 펼쳐지므로 parent 걸음 수로는 접미사를 자를 수 없다. 쓰기는 «접두사 사슬을 잇는» 자리(`transform`)에서만 일어나고, 불리언의 결과 정점은 표를 건드리지 않는다.

## Tessellation 층 (`nacre-tess`) — 출처 태그 파생물

메시는 진실이 아니라 파생물이고, 모든 요소가 출처를 안다. **길은 하나다**: `tessellate(model, cfg)` →
`Tessellation`. 고른 솔리드만 지을 때는 `tessellate_solids(model, solids, cfg)` — 같은 길이고 걷는 범위만 다르다
(간선 표본과 면 삼각분할은 하나씩이고 두 솔리드는 간선을 나누지 않으며 θ 집합은 제 솔리드의 면으로만 자라므로, 면의 삼각형은 곁에 무엇을 짓든 같다). 간선을 샘플해 공유 폴리라인으로 두고(`finish_polylines`) 모든 면을 자기 곡면의 차트에서 삼각분할한다.
`Adjacency` 처럼 **처음부터 통째로** 짓는다.

```rust
pub enum TessOrigin {
    OnVertex(Handle<Vertex>),
    OnEdge { edge: Handle<Edge>, t: f64 },       // 곡선 파라미터
    OnFace { face: Handle<Face>, uv: [f64; 2] }, // 곡면 파라미터
}

pub struct TessVertex { pub pos: Point3, pub origin: TessOrigin }

pub struct TessTriangle {
    pub vertices: [Handle<TessVertex>; 3],
    pub face: Handle<Face>,   // 출처 태그 — 하이브리드 불리언의 핵심 재료
}

/// tess 컨테이너. 면/엣지별 버킷을 둔다.
pub struct Tessellation {
    pub vertices:  Store<TessVertex>,
    pub triangles: Store<TessTriangle>,
    pub by_edge:   HashMap<Handle<Edge>, Vec<Handle<TessVertex>>>, // 엣지 polyline
    pub by_face:   HashMap<Handle<Face>, Vec<Handle<TessTriangle>>>,
}
```

**닫힌 간선의 폴리라인은 첫 점을 반복하지 않는다 — 닫힘은 암묵이다.** 온전한 rim 은 제 꼭짓점에서
한 바퀴를 돌며 샘플하고 그 꼭짓점 앞에서 멈춘다(핸들을 반복하면 면이 퇴화 삼각형을 얻는다). 폴리라인을 **쌍으로** 걷는 소비자는 `n` 각
링에서 `n − 1` 걸음을 얻으므로 닫는 걸음을 스스로 더해야 하고, 분기할 사실은
`edge.vertices[0] == edge.vertices[1]` 이다.

계획: 증분 갱신. 면/엣지별 버킷은 연산이 면 몇 개만 바꿨을 때 그 버킷만 무효화·재생성하기 위한 모양이고,
무효화 표시(`stale` 집합)와 「연산이 위상에 항목을 추가하면 같은 트랜잭션에서 tess 에도 추가한다」
(Fornjot 의 "함께 쌓기")가 그 도착점이다. tolerance 를 바꾼 재계산은 tess 만 통째로 재생성하고
위상·기하는 불변이다.

이 층의 용도는 셋이다. 뷰어가 이걸 그대로 그린다(별도 tessellation 경로 없음). 메시 내보내기(`Tessellation::to_obj`; 계획: STL/3MF)가 이걸 그대로 쓴다. 그리고 계획: 강건 메시 불리언의 입력이 되어 "조합적 결정 → 출처 태그로 정확 기하 스냅백" 파이프라인의 앞단이 된다. uv/t 파라미터를 들고 있으므로 스냅백 시 Newton 초기값이 공짜로 나온다.

참고(성능): 불리언·교차 전에 면쌍 후보를 AABB/BVH로 컬링하는 공간 인덱스가 필수 과제다 — tess의 면별 버킷에서 AABB가 거의 공짜로 나온다. (Truck은 이게 없어 전수 비교에 가깝고, 저자도 "BSP 등 최적화는 미래 과제"로 명시 — 반면교사.)

**틀린 메시는 없느니만 못하다 — 못 그리는 면은 거절한다.** `TessError` 에는 일부러 폴백 변종이 없다:
`DegenerateRing`(링이 「형제 구멍을 가진 다각형」이 아니다 — 세 정점 미만, 넓이 0, 인덱스로 공유·반복된
정점, 스파이크, 서로 지나가는 두 선분, 한 점을 두 번 지나며 그 점을 가로지르거나 한 부채꼴이 다른 부채꼴 안에
드는 링; 구멍이 외곽 링을 가로지르는 것도 잡는다), `SelfTouchingBoundary`
(아래), `HoleWinding`(구멍이 외곽과 같은 방향으로 감겼다 — 고칠 메시가 아니라 깨진 솔리드), `OverBudget`
(아래). `tessellate`·`tessellate_solids` 는 받은 범위를 걷다 첫 거절에서 멈춘다.

**삼각분할 알고리즘 — 브리징 없는 스윕.** 구멍 있는 면은 **y-단조 분해 + 단조 삼각분할**(de Berg §3)로 자른다. 구멍을 외곽 링에 브리지(폭 0 슬릿)로 꿰매 귀 자르기를 하지 않는 이유는 하나다: 브리지가 정점을 **반복**시키고, 반복 정점이 있는 링은 단순 다각형이 아니며, Meisters 의 two-ears 정리는 단순 다각형에만 성립하므로 그런 링에는 **귀가 아예 없을 수 있다**. 스윕은 **링을 절대 합치지 않으므로** 퇴화 링 자체가 안 생긴다 — 구멍 간선이 보통 간선이고, 구멍의 최상단이 split·최하단이 merge 정점이 되어 대각선이 자동으로 구멍을 잇는다. 판정은 전부 부호이고 **tolerance가 없다**: 전부 exact `orient2d`(`nacre-predicates`)이고, 이것은 **어떤 f64 입력에도 정확**하다. (평면의 차트는 축 드롭이라 `uv`가 원본 f64 그대로지만, 원통의 차트는 `atan2`로 **계산된** 값이다 — 그래도 술어가 정확하므로 알고리즘은 같다. 면 사이에서 `uv`를 비교하는 일이 없으므로 — 면끼리 만나는 자리는 **엣지 폴리라인**이 이미 정한다 — 차트가 면마다 달라도 무방하다.) 경계 정점을 추가하지 않으므로 삼각형 수는 `V + 2H − 2L`(`H` 는 시계 고리, `L` 은 반시계 고리 — 보통 하나, 한 점에서만 붙은 조각과 구멍 안의 섬마다 하나 더)로 불변이고 crack-free 규칙도 그대로다. 스윕은 자기 전제를 **가정하지 않고 검사한다**: b-rep 이 보장하지만 이 층은 보장할 수 없는 것이라, 자기 교차 링을 받은 분해는 그러지 않으면 확신에 찬 틀린 메시를 낸다.

**경계가 자기에게 닿는 면 — 접촉 자리의 공유 간선 다리.** 구멍이 다른 링에 정확히 한 점에서 접하는
평면 면은 유효한 솔리드지만(면이 집혔을 뿐 곡면은 그 점에서 2-다양체다 — 집힌 면의 두 엽이 이웃 곡면을
돌아 이어진다) 두 링으로는 분해할 수 없다: 닿는 정점은 **곡선** 간선의 샘플이고, 그것이 내려앉은 직선
간선은 두 끝만 든다. `tessellate` 의 사전 패스(`bridge_shared_edges`)가 그 **이미 있는** 정점을 직선
간선의 폴리라인에 끼워 넣고, 두 링을 그 점에서 하나로 이으며, 스윕은 생겨난 일치 쌍을 기호적으로
순서 짓는다(`polygon::sos`). `by_edge` 가 공유 구조라 이웃 면도 같은 정점을 보므로 crack-free 가 유지된다
(면이 혼자 경계를 바꾸는 것이 아니라 **간선의 폴리라인**이 바뀐다). 접촉은 어떤 폴리라인도 바꾸기 전에
**모든** live 면에서 먼저 모으고, 한 간선의 삽입은 뒤 위치부터 한다; 패스는 멱등이다. 간선→면 표는
live 집합에서 여기서 짓는다.
**한 링이 한 메시 정점을 두 번 지나는 면**(원판에 판 홈의 끝이 테두리에 닿은 캡, 구멍과 보스 자국이 한 꼭짓점에서
만난 판의 윗면)은 이을 것이 없다 — 두 방문이 이미 한 링의 쌍둥이이고(`polygon::Revisit::Pinch`; 다리가 이어 붙인
고리의 두 방문도 같은 규칙이 이름 짓는다), 꼬집힌 점이 여럿이면 점마다 한 쌍이다. 각 방문이 상대의 반직선을 하나도
품지 않으면(자기 재료 부채꼴만 감싼다) 엄밀 볼록한 부채꼴은 자기 반직선을 따라 미끄러지고 볼록이 아닌 하나는
제자리에 둔다(`polygon::pinch`). 서로의 반직선을 둘 다 품으면 재료가 그 점에서만 붙은 두 엽이다 — 방문의 짝을 바꿔
(각 방문이 상대의 들어오는 변을 받는다) 링을 그 점에서 둘로 나누면 그 점은 위의 집힘이 되고, 스윕은 모든 변의 왼쪽이 재료인
링이면 바깥 링이 여럿이어도, 섬이 구멍 안에 있어도 한 번에 그린다(판 띠의 윗면을 핀의 반원 발이 먼 모서리 한 점까지 깎은 면).
다리가 이은 링의 두 엽은 나누면 이음을 되돌릴 뿐이라 `SelfTouchingBoundary` 다. 그 밖은 그 점을 지나가는 링이거나 겹친
부채꼴이다(`DegenerateRing`). 꼬집힌 점이
고리의 극점이어도 링 방향은 그 점의 모든 반직선 가운데 바깥 틈 다음의 첫 반직선이 정한다(`ring_orientation` —
커널 감김과 같은 극점 판독).
이 길이 덮지 않는 모양은 고치지 않고 **세어서** 거절한다(`Declined`: `CurvedEdge`·`CurvedNeighbour`·
`SharedSegment`·`MultiSegment`·`NoEdge`; 보고는 `bridge_report` → `BridgeReport`) — 다른 링의 **정점**에
닿는 접촉, 닿은 선 위에 정확히 놓인 이웃, 다리가 둘 이상이거나 다리와 꼬집힘이 함께인 면, 세 번 지나는 점도
그렇다. 그것들이 `SelfTouchingBoundary`
로 돌아온다. 이 거절은 솔리드가 틀렸다는 말이 아니다.

**곡면도 같은 길을 탄다 — 면마다 «차트» 하나.** **모든 면이 같은
길**을 탄다: 경계 루프를 **그 곡면의 차트**에 올리고 위의 단조 스윕 + 플립을 돌린다.
분기는 하나이고 **고르는 것은 알고리즘이 아니라 차트**다 — 평면은 축 드롭(Newell 법선 → 축 드롭 →
손잡이 보정, `planar_chart` 한 곳), 원통은 `(z, r·θ)`.
원뿔·구·NURBS는 **차트 하나씩** 더하면 된다(새 곡면이 특수 경로를 곱하지 않는 이유).
- **`θ`가 아니라 `r·θ`**: 원통은 가전개면이라 이 스케일이 **등거리 사상**이고, 그래야 플립이
  「그림 속」이 아니라 **곡면 위에서** 메쉬를 개선한다.
- **손잡이 보정은 교환이 아니라 `u` 부호 뒤집기**: `v`가 스윕 축이므로 교환하면 스윕이 축을 따라
  가고, 대각선이 넓은 호를 가로질러 **곡면을 벗어나는 현**이 된다(교환 보정에서는 현이 32°를 가로질러
  sag가 예산의 254배인데 watertight도 삼각형 수도 **그대로 통과**한다).
  그래서 「삼각형은 곡면을 벗어나지 않는다」가 **별도 잠금**으로 서 있다.
- **θ는 루프를 따라 풀어서** 쓴다(절대각이 아니라). 밴드는 두 림으로 닫혀 있어 펼치면 열린 선 둘이므로
  한 생성선에서 잘라 잇는다(`cut_seamless_bands` — 「seam 엣지 방식」): 절단점이 θ 와 θ+2π 에 한 번씩 서서
  펼친 밴드가 **직사각형**이 된다. 남은 구멍은 온전한 회전만큼 바깥 링의 창으로 옮긴다.
- **삼각형 수는 유도대로**: 밴드의 차트 다각형은 정점 `2n+2` ⇒ `V+2H−2 = 2n`
  (원통 전체 `4n−4` 잠금과 맞는다).

**모든 간선은 예산 안의 현이고, «차트»가 그것을 가능하게 하는 점을 댄다.** `TessConfig`의 두 예산은
경계에서는 구성으로 지켜지지만(아래 «곡선을 얼마나 쪼갤지는 예산 «둘»이 정한다» 글머리), 면의 **내부**는 자기
경계가 우연히 대 준 샘플을 물려받는다. 경계가 곡률 방향을 아크로 덮는 동안은 공짜지만, 덮지
않으면 **조용히 틀린다**: 병합된 옆면(지워진 가짜 이음매가 한 구간의 유일한 아크였던 면)은 내부 점
없이는 반 바퀴를 한 번에 건너는 현으로 그려져 **넓이의 1/6**을 잃는다. 아크가 모든 θ 를 덮어도 룰링이 두
림 사이에서 끊기면(다른 면의 간선이 룰링 중간에 꼭짓점을 낸다) 그 꼭짓점의 높이에는 아크가 없어, 거기서 나가는
내부 간선이 각 예산을 넘는다.

그래서 규칙과 공급을 **두 층으로** 나눈다.

- **규칙**(일반·자기 검증): *모든 내부 간선과 모든 면이 선언된 예산 안이다* — 간선의 선형은
  `Surface::distance(3D 현의 중점) ≤ tol`, 각도는 두 끝의 **곡면 법선 사이 각** `≤ max_angle_deg`; 면은
  면 법선과 무게중심의 곡면 법선 사이 각에, 반올림이 그 면 법선을 얼마나 돌렸을 수 있는지를 더해 `≤ max_angle_deg`
  — 법선을 제 점들이 정하지 못하는 면(넓이 0, 외적이 반올림뿐인 실 삼각형)은 밖이다(`facet_in_budget`). 간선이 다
  예산 안이어도 면은 기울 수 있다 — 원통 축에 비스듬한 곡선 위에 늘어선 세 표본은 그 곡선의 평면 안에 눕는 면이다.
  판정은 메시가 내보내는 점(경계는 공유 폴리라인의 정점, 새 점은 차트를 되돌린 점)에서, 원소마다 모서리를 번호 순으로
  놓은 한 철자로 한다. 면에 거리 항이 없는 까닭: 평면의 면은 평면 안이고, 원통에서 반 바퀴 안의 면이 곡면에서 가장
  먼 점은 그 변의 중점이다(`facet_in_budget` 의 doc). 모두 `Surface`에 전수 대응이라 평면은 반올림이 법선을 정하지 못하는 실 삼각형만
  밖이고(코퍼스·스위트 0) 새 곡면은 컴파일 에러다. 넘으면 `TessError::OverBudget` — 고치지 않고 **정직하게
  거절**한다(틀린 캐시는 없느니만 못하다). 경계 간선은 동어반복이라 보지 않는다: 곡선은 현마다 그 자리에서 두 항의 **최댓값**을 지키게
  나뉘고(`walk` — 면의 집합이 `FINEST_CIRCLE` 열을 넘으면 `OverBudget`), 직선 엣지는 거리 0이고 양 끝 법선이 같다. 곡선 안의 꺾임은 법선이 묻지 않으므로 곡선이 그 자리의
  회전율로 칸을 요구하고(`walk`) — 옆면의 곡선은 그 면의 한 집합 위에 선다(`FaceThetas`) — 시험이 잰다(`kinked_edges`).
- **공급**(곡면을 아는 층): 원통 차트가 **자유 파라미터 없는 격자**를 내놓는다. 곡률 방향은
  그 면의 **θ 집합**(`FaceThetas`) — 그 면의 곡선들이 표본을 두는 바로 그 열이라, 줄은 곡선들이 선 생성선 위에
  서고 줄과 곡선 사이의 띠도 두 곡선 사이처럼 생성선을 품는다 — 그리고 열은 차트의 **풀린** `v` 범위 전체를 덮는다
  (구멍은 온전한 회전만큼 옮겨 다닌다). 축 방향은 **경계 «원 엣지»마다 그 원의 자리**(`station` — 중심의
  축 좌표) 하나(타원 간선의 점은 곡선 자신의 높이 `curve_place` 에 서고 줄을 놓지 않는다 — 올라가는 선이라 격자가
  룰링처럼 대한다), 그리고 **룰링이 끊긴 점마다 그 점을 끊는 평면의 자리**(`ruling_break_place` — 면에 그
  평면의 원 간선이 있으면 그 곡선의 `station`, 없으면 그 평면이 원통에서 자를 원의 것 — 끊는 두 평면 가운데 축을
  가로지르는 것이 있으면 저장 순서와 무관하게 그것이, 없을 때만 비스듬한 평면의 타원 높이가 답한다; 한 평면 위의 두
  끊김과 그 평면의 원은 한 `f64` 다): 곡선에서 읽지 샘플된 점에서 읽지 않는다(한 림의 점들은 그 좌표가
  ulp만큼 다르므로 중복 제거하면 림 하나가 181개가 된다). **그 수는 한 철자다** — 띠를 여는 전처리가 생성선 위 점을 그것으로
  정렬하고, 차트는 원 위의 모든 점에 그것을 `u` 로 주고, 격자는 그 높이에 줄을 놓는다. 그래서 줄 가운데
  아크 위의 점은 아크 현 위에 **비트로** 놓이고 삼각분할의 정확 술어가 버린다. 원 위의 점이나 끊긴 점의
  좌표에서 `u` 를 따로 읽으면 기울어진 축에서 몇 ulp 갈려, 그 점이 재료 쪽으로 들어가 옆면에서 기운 조각
  삼각형이 된다(37°·23° 에서 옆면 하나에 25° 넘게 기운 것 최대 50, 최악 90° — 캡 평면에 눕는다). **룰링과의 겹침은 다른 물음이다**: 사용자의 모서리 각이 격자 칸의 배수에 놓이는
  것(벽 `x = 0.5`, `r = 1` → 30° = 2° 의 열다섯째 칸)은 한 수를 두 번 계산한 것이 아니라 기하가 같게 만든
  두 수라 한 곳으로 모을 수 없다 — 올라가는 경계 선분에서 칸의 `1e-6` 배 안의 점은 대지 않는다(`ON_RULING` — 그 doc 이
  문턱의 세 독자와 잰 틈을 든다; census 메시는 1e-12 ~ 1e-3 에서 같다).
  **차트의 `u` 범위 끝에 선 스테이션에는 줄을 놓지 않는다** — 원만 면을 두를 때 그 둘은 면 자신의 두 림이라
  줄의 점이 모두 림 현 위이거나 면 밖이다(놓으면 보통 밴드가 358점을 받고 전부 버려진다 — 메시는 같고 일만 는다);
  타원은 제 원의 자리 너머까지 닿으므로 범위는 스테이션의 순서가 아니라 차트에서 읽는다. 평면 차트는 아무것도
  대지 않는다.

삽입은 후보마다 **놓을 자리**를 찾는다 — 삼각형의 **엄격한 내부**(셋으로 쪼갬)이거나 **내부 간선
위**(둘로 쪼갬). 두 번째가 예외가 아니라 **주된 길**이다: 후보는 면의 모양이 바뀌는 선 위에 앉고
삼각분할은 이미 그 선을 따라 간선을 갖고 있다. **그러나 두 길 다 프로덕션에서 산다** — 네 벽과 코너에 선
보스의 옆면 다섯(노치가 구멍인 밴드)에서 내부/간선이 `−x` 136/42 · `+x` 163/15 · `−y` 130/48 · `+y` 163/15 ·
코너 156/112 이다. 그 비율은 사용자가 그린 모양이 아니라 밴드를 자른 생성선의 자리가 정하므로 어느 갈래도
도달 불가로 지울 수 없다.
**제약 간선은 절대 안 쪼갠다** — 이웃이 모르는 점은
T-정점이고 그건 크랙이다. 자리가 없는 후보는 **버린다**(강제하지 않는다).
`ChartMap`이 되돌리는 길이다(차트는 전단사다) — 발명된 점은 위치를 읽어 올 모델 정점이 없고, 두 팔
모두 정방향이 건 **손잡이 보정을 되돌려야** 한다.

**계약과 품질을 한 뒤집기 패스가 든다.** 단조 삼각분할은 정확하지만 슬리버가 많다. 자유(비제약) 간선에 **Lawson 플립**을 돌려 제약 Delaunay로 보내면 최소각이 최대화된다(정리 — Lawson 1977, Chew 1989). 새 의존성 0(`incircle`이 이미 있다). 경계는 제약이라 안 뒤집히므로 watertight·개수·면적 전부 불변. **예산이 품질보다 앞선다** — 뒤집기 하나는 원소 셋(대각선과 곁의 두 면)을 바꾸고, 그 가운데 예산 밖인 수가 줄면 뒤집고, 같으면 `incircle` 이 엄격히 양일 때 뒤집는다(`delaunay::refine` 의 `within`, 판정은 `within_budget` 과 한 철자 `edge_in_budget`·`facet_in_budget`). 원소마다 판정이 한 철자라 메시 전체의 «밖» 수는 그 차이만큼 움직여 늘지 않고, 같은 수의 뒤집기는 Delaunay 리프팅 퍼텐셜을 줄이므로 (수, 퍼텐셜) 이 사전식으로 줄어 끝난다. 면이 모두 안이면 이것은 간선의 규칙이다 — 예산 밖 대각선은 예산 안 대각선으로 `incircle` 과 무관하게, 둘 다 안이거나 둘 다 밖이면 `incircle` 대로. 키 큰 띠의 두 림이 수평이 아니면(원 곁의 타원) 그 외접원이 표본 하나를 건너뛰는 대각선을 부르고(비스듬히 자른 봉에서 4°, 예산 2°), 곡선 위 연속한 세 표본의 «귀»는 외접원이 그 곡선의 곡률 반지름이라 Delaunay 가 고른다 — 둘 다 계약이 뒤집는다. 평균 최소각은 2배가 된다(원통 캡 tol 1e-4에서 0.51°→1.05°). **다만 캡의 *최악*값은 안 움직이고, 그건 알고리즘 탓이 아니다** — 정n각형 정점은 **공원**이라 `incircle`이 전부 0을 답하고(모든 삼각분할이 Delaunay), 게다가 내각 178.9°와 two-ears가 0.57° 이하 슬리버를 **강제**한다(측정 0.457° = 하한). 고치려면 점을 **추가**해야 하는데, 그건 *캡의* 얘기다 — 캡의 경계는 이미 곡률을 다 샘플하므로 더 넣을 이유가 **품질뿐**이고, 품질은 이 층의 결정이 아니다. (금지되는 것은 «점 추가» 자체가 아니라 **경계** 정점의 추가다 — 위 «모든 간선은 예산 안의 현이고, «차트»가 그것을 가능하게 하는 점을 댄다» 글머리를 보라.)

**곡선을 얼마나 쪼갤지는 예산 «둘»이 정한다.** 원의 분할 수는
`n = max(⌈360°/Δθ⌉, ⌈π/acos(1 − tol/r)⌉)`이다. 뒤의 항(절대 sagitta)만 쓰면 `n ≈ π√(r/2tol)`이라
**작은 원일수록 적게 쪼갠다** — 절대 오차 예산은 지키지만 반지름 0.2가 10각형이 되고, 확대하면 그대로
드러난다. 앞의 항(각도)은 `s/r = 1 − cos(π/n)`에서 반지름이 사라지는 **상대** 예산이라 크기와 무관하게
같은 품질을 준다. 정n각형에서 **중심각 = 선분의 꺾임각 = 인접 면 법선의 벌어짐**이 모두 같으므로 이 한
숫자가 윤곽과 음영을 동시에 묶는다. 기본값 `Δθ = 2°`(=180분할,
반지름의 0.0152%), `tol = 1e-2` — **어느 반지름에서도 한쪽만 쓸 때보다 나빠지지 않는다**(어느 항이
이기는지는 `per_turn`). **곡선의 모든 현은 그것이 놓인 자리에서 예산을 지킨다** — θ 폭 `h` 인 현 위에서 가장 빠른
회전율 `R`(옆면 법선의 `1` 과 곡선 접선의 큰 쪽)로 `h·R ≤ Δθ`, 가장 큰 처짐 반지름 `ρ` 로 `ρ(1 − cos(h/2)) ≤ tol`
(곡선의 몫은 `Ellipse::pace_over` 가 구간마다 낸다). 현은 가운데 매개변수의 접선과 평행해 이웃 현의 꺾임은 반 칸마다
`R·h/2` 이하다. 현을 더 쪼개도 규칙은 남아 꼭짓점·엇갈림 열이 얹혀도 선다. 열은 원통의 `ref_dir` 에서 걷는 걸음이
짓는다(`walk`) — 원만 있는 면은 고른 `τk/N` 그대로, 타원은 빨리 돌거나 많이 처지는 곳에서만 잘아진다. OCCT
`GCPnts_TangentialDeflection` 도 현마다 국소 기준이다(그 각은 현 하나의 반 회전이라 식은 우리 것).
**옆면의 곡선은 제 칸이 아니라 그 면의 θ 집합 위에 선다**(`FaceThetas`) — 그 면의 곡선들 모두를 읽은 걸음의
열, 그 면 모든 꼭짓점의 θ(원통 틀에서 꼭짓점마다 한 번, 한 θ 에 선 꼭짓점은 `ON_RULING` 안이면 한 열), 엇갈린
현이 청한 열(아래). 열의 수는 하나다 — 표본은 그 수에서 짓고 출처 `t` 에 적는다(`FaceThetas::columns`). 두 곡선이 같은
생성선에 표본을 두어야 그 사이 띠의 면이 생성선을 품는다 — 제 칸으로 따로 나누면 띠가 어긋난 표본의 처짐 / 띠
높이만큼 기울어, 각 예산 / 처짐보다 얇은 띠는 `OverBudget` 이었다(반지름 50 원판의 림을 깎은 호, 기울기 1/20 의
타원과 그 원). 범위가 원통이 아니라 면인 것은
불리언이 원통 핸들을 재사용하기 때문이다 — 핸들 단위면 곁에 지은 솔리드의 꼭짓점이 면의 표본을 바꾼다
(`tessellate_solids`). 가파른 타원의 잔 칸은 같은 면의 원 림·캡·격자 줄도 받는다(기울기 2 의 봉에 정거장 하나:
삼각형 904).

**틈 없음(crack-free) 규칙:** 면의 삼각분할은 반드시 `by_edge`의 공유 polyline 정점들을 자기 경계로 소비해야 한다. 인접한 두 면이 각자 독립적으로 엣지를 샘플링하면 공유 엣지에서 정점이 어긋나 T-junction이 생기고 watertight가 깨진다 — 엣지 polyline이 먼저, 면 삼각분할이 그걸 경계 조건으로.
**금지되는 것은 «경계» 정점의 추가다.** 면이 자기 **내부**에 점을 넣는 것은 이웃이 볼 일이 없으므로 크랙과 무관하고, 위 «모든 간선은 예산 안의 현이고, «차트»가 그것을 가능하게 하는 점을 댄다» 글머리가 바로 그것이다. 폴리라인을 더 잘게 하는 것도 원리상 금지가 아니다 — `by_edge`가 공유 구조이므로 **엣지 샘플러가** 촘촘히 내면 양쪽 면이 함께 본다(OCCT도 `BRepMesh_ModelHealer`에서 그렇게 한다; 위 접촉 다리도 같은 원리다). 금지는 **면이 혼자** 경계를 바꾸는 것이다.

**엇갈린 현은 간선에서 나눈다.** 면의 경계 곡선은 만나지 않지만 원호의 현은 원보다 최대 sagitta 안쪽이라, 그보다
가까운 두 곡선의 현은 서로 지나갈 수 있다 — 테두리에서 sagitta 안까지 온 구멍이나 홈의 끝, 테두리 위의 끝에서
접선에 반 칸(1°)보다 가깝게 떠나는 홈의 변. 현을 나누지 않으면 면은 유효한데도 스윕이 그 링을 거절한다. 그래서 면을 삼각분할하기
전에 스윕 자신의 정확한 교차 검사(`polygon::meets`)로 경계가 스스로를 지나는 곳을 묻고(`crossing_chords`), 엇갈림에 든
원호의 현이 **전부** 그 호의 중점(현이 잇는 두 열의 수 사이 — 한 옆면의 두 림이 한 자리에서 엇갈리면 한 수다)에 열
하나를 청한다 — 폴리라인이 아니라 **그 호가 두르는 옆면의 θ 집합**에. 그래서 나눈 θ 는
그 옆면의 모든 곡선과 격자가 함께 받는다(호 하나에만 넣으면 짝 림과의 얇은 띠가 기운다 — `FaceThetas`). 바퀴마다 모든 간선을 집합에서 새로 표본하고 모든 평면 면의 차트를 다시 짓는다 — 한 캡이 청한 열이 옆면을 타고 다른 캡의
경계를 옮긴다 — 어느 집합도 자라지 않는 바퀴까지(`finish_polylines`). 멈춤은 면의 집합이 담는 가장 많은 열(`FINEST_CIRCLE` = 4096 — 걸음이 거절하는 수와
한 문장, `FaceThetas::full`)이다: 둘레 전체가 `1e-8` 로 붙은 거의 동심 구멍은 원 하나에 2만 점 넘게 필요해 거기 닿고,
이름으로 거절된다.
**순서가 정확성의 일부다** — 띠 자르기 뒤다(원호에 표본을 넣어 현을 곡선 쪽으로 미니 엇갈림을 만들 수 있다).
그리고 **한 차트는 한 번 묻는다** — 검사는 차트의 점과 링만 읽으므로 앞 바퀴와 비트로 같은 차트는 앞의 답을 쓰고, 마지막
바퀴의 답을 다리 사전 패스와 그 면의 삼각분할도 읽는다(다리는 직선 위에 정확히 놓인 점만 넣어 모양을 안 바꾼다; `polygon::Asked`). 한 함수가 그
순서를 든다(`finish_polylines`). 평면 면만 묻는다 — 원통 차트에서 원은 수평선, 룰링은 수직선, 타원은 θ 의 그래프라 샘플링이
엇갈림을 못 만든다.

## 진실과 캐시 — 타입 구조

아레나는 append-only 이므로 **진실만 들어간다.** 캐시는 곁표에 두어 언제든 버리고 재생한다(정점의 구성 수치만은
정의가 아니라 로그에서 — 「위상 층」).
정의가 불변이므로 **캐시는 낡을 수 없다** — 무효화라는 개념이 없고, «버리고 재생»만 있다.


## 진실의 타입

```rust
// ─── nacre-topo ── 진실 (아레나, append-only) ──────────────────────────

/// 아레나(`surfaces: Store<Surface>`)가 이것을 든다. 실현은 색인 평행 곁표
/// `surface_cache: Vec<SurfaceCache>`(쓸 수 있어야 캐시이므로 `Store` 가 아니라 사설 `Vec`).
/// ⇒ **맨이름 `Surface` = 이것**이고, geom 의 f64 실현은 `nacre_geom::Surface` 로 적는다.
pub enum Surface {
    Plane {
        points: PlanePoints,
        motion: Option<Handle<MotionNode>>,   // None = 세계, Some = 점들이 적힌 모션 전 프레임의 이력
        sense: Orientation,                   // 평면의 향: 점들의 세계 방향 `parity · L((p₁−p₀)×(p₂−p₀))` 그대로(Forward)/반대(Reversed)
                                              //   면의 바깥 = 평면의 향 × 면의 Orientation. interning 열쇠에 안 든다(첫 진술의 향)
    },
    /// 성분형이 옳은 이유·`ref_dir` 원시 규칙은 원통 절이 상세히 적는다.
    Cylinder {
        def: CylinderDef,                     // { origin, dir, ref_dir }: Rat + r2: BigRat — 반지름은 «제곱», 폭은 임의정밀
        motion: Option<Handle<MotionNode>>,
    },
}

/// 평면의 세 점 — 진술의 종류가 곧 변종이다. `Known`(Rat 아홉, 288 B)은 인구의 거의 전부라
/// 상자에 넣지 않는다.
pub enum PlanePoints {
    /// 구성 평면(벽·캡·사용자 진술 평면): **모션 전 프레임**의 유리수 세 점(공선 아님).
    /// base case — 여기서 끝난다(C2).
    Known([[Rat; 3]; 3]),
    /// datum 평면: 모델 정점 셋을 지난다. 발견된 좌표는 Rat 에 안 들어갈 수 **있으므로**
    /// (타입이 그것을 허용한다 — 코퍼스는 안 때린다) **가리킨다**(C3). 필요할 때 실현한다.
    /// 생산자는 `DatumDef::ThroughVertices`. 능력의 핵심은 폭이 아니라
    /// **좌표로는 그 평면을 말할 수 없다**는 것 — 발견된 정점에서 읽은 좌표는 반올림된 값이고
    /// 그것으로 지은 평면은 «다른» 평면이다(기울어진 기하에서 220/220, `point_width.rs`).
    /// **핸들은 정렬해 저장한다** — 같은 셋 = 같은 진술. 순서가 다르다고 같은 평면이 두 핸들이
    /// 되는 것을 구성 시점 정규화로 막는다. 법선 방향은 점 순서가 아니라 정준 부호 규약과
    /// `Face::orientation`/`flip` 이 들므로 잃는 정보가 없다.
    /// **모션은 대체가 아니라 합성이다** — 정점들이 자기 프레임에서 평면을 박고 `motion` 이
    /// 그것을 밖으로 나른다. 옮기는 쪽은 핸들을 그대로 두고 노드만 기록한다(핸들도 옮기면 두 번 움직인다).
    /// **`motion` 의 뜻은 이름이 있느냐로 갈린다**: 이름 있는 진술의 이름은 세 정점이 만나는 프레임 F 의
    /// 진술이라 `motion` 이 F 를 **담는다**(F 이거나 F 를 잇는 사슬); 이름 없는(혼합 프레임) 진술은 정점마다
    /// 자기 사슬로 세계에 놓이고 `motion` 은 그 **뒤**의 이동이다. 둘 사이의 진술(이름은 서는데 `motion` 이 F 를
    /// 잇지 않는)은 어느 뜻인지 모호하므로 문(`push_plane_through`)이 release 에서도 단언으로 멈춘다.
    Through([Handle<Vertex>; 3]),
}

/// 정점 = 자기 **정의**, 그리고 정의가 곧 타입이다 — `Surface` 와 같은 모양. 좌표는 진실이 아니라
/// 캐시다(`PointCache`, 색인 평행 `vertex_cache`).
/// 단일형이 아니라 **네 변종**이다: 원통 seam 정점은 «세 평면»으로 적을 수 없다.
/// 불변식은 변종별이고 그것이 이 구조의 근거다.
pub enum Vertex {
    /// 세 평면의 교점 — 이름이 곧 점. D ≠ 0(세 평면이 한 점에서 만난다)은 이름이 한 프레임을 말하는 곳에서
    /// `Model::vertex_named_meet` 한 곳이 정확히 답한다(`NamedMeet::NoPoint`) — `validate` 가 그 정의를
    /// `VertexNoPoint` 로 보고하고 실현은 `RealizeError::NoPoint` 로 거절한다. 이름이 한 프레임을 말하지 않는
    /// 정점(혼합 프레임 `Through` 담체 · 세계 이름 없이 갈리는 사슬 — 세 평면 push 의 census 3.3%, 스위트
    /// 4.0%)은 묻지 않는다.
    ThreePlane([Handle<Surface>; 3]),
    /// 두 곡면의 교차 «곡선» 위의 점 — 원통 seam(테두리의 θ=0). 쌍은 곡선을 박고, 그 위의
    /// 점을 고르는 것은 원통 진실의 `ref_dir`(`CylinderDef`)이다 — 그래서
    /// `OnSeam([cylinder, cap])` 은 「rim ∩ `origin + r·ê` 를 지나는 생성선」으로 **정확히 지정된 한 점**이다.
    /// 정의에서 좌표를 재생하는 기계는 `realize_vertex`(`seam_point`)이고, push 깔때기가 같은 도로로
    /// 실현해 캐시에 넣는다.
    OnSeam([Handle<Surface>; 2]),
    /// 이 점은 도법기하의 **관통점**(piercing point)이다 — 두 평면의 교선이 원통 옆면을 뚫는 자리.
    /// 담체 종류를 **구조로** 말한다(닮은 핸들 셋이 아니라). `rank` 는 `CylinderPair` 와 같은
    /// 순서, 곧 원통 둘레 seam 먼저 각도 순서에서 몇 번째인가다(아래). 한 직선의 두 교점은 각이
    /// 같을 수 없고(같으면 그 선이 룰링), 접하면 점 하나가 `rank: 0` 이다. 순위는 원통 진술만
    /// 읽으므로 두 평면의 순서도 법선의 향도 이름에 닿지 않는다 — 평면을 다시 정렬하는 쪽도
    /// 순위를 다시 말할 일이 없다.
    Pierce {
        planes: [Handle<Surface>; 2],         // 오름차순 핸들 — 이름은 읽지 않는다
        cylinder: Handle<Surface>,
        rank: u8,                             // 원통 둘레 각도 순서 자리; 접하면 0 하나
    },
    /// 평면 하나가 두 원통과 만나는 점 — 프리즘 캡 위에서 서로 다른 원의 두 호가 만나는 모서리.
    /// `rank` 는 점이 `cylinders[0]` 둘레 seam 먼저 각도 순서에서 몇 번째인가(축이 평행할 때
    /// 그것이 무엇을 고르는지는 「원통의 진실」 7). 원통 둘레 순서의 철자는 하나다
    /// (`nacre_exact::quad::seam_ranked`). 거울은 각도의 방향을 뒤집지만 원통이 있는 솔리드의
    /// 거울은 거절된다(`MirrorNotPlanar`) — `Pierce` 도 같다. 변환은 두 핸들을 제자리에서 바꿔
    /// 적으므로 아무도 순위를 다시 읽지 않는다 — 정렬하지 않는 까닭이다.
    CylinderPair {
        plane: Handle<Surface>,
        cylinders: [Handle<Surface>; 2],      // cylinders[0] 이 순위의 기준 — 정렬하지 않는다
        rank: u8,                             // 공유 룰링의 각도 순서 자리; 접하는 두 원은 0 하나
    },
}

/// 곡면 집합은 «담체»(어디 위에 있나)를 정하고, 경계가 «어느 조각»인지를 정한다.
/// **곡면 집합만으로 남김없이 정해지는 것은 정점뿐이다** — 차원 0, 평면 셋이 자유도 셋을 다 먹는다.
/// 간선(1차원)은 담체 위 «어느 조각»인지 끝점 둘이, 면(2차원)은 «어느 영역»인지 루프가 더 정한다.
/// 그래서 정점의 진실은 핸들 셋으로 끝나고, 좌표는 어느 실체의 진실에도 없다(전부 캐시).
/// 간선의 `vertices` 는 «집합»이 아니라 «순서»다 — 원 담체에서 `[A, B]` 는 A→B 를 축에 대해
/// 반시계로 가는 호이고 `[B, A]` 는 그 나머지다(`surfaces` 는 오름차순 정렬한 집합인 것과 대비).
pub struct Edge {
    pub surfaces: [Handle<Surface>; 2],       // 담체 (두 끝점 면집합의 교집합으로 파생 불가)
    pub vertices: [Handle<Vertex>; 2],        // 경계
}
// Store<Curve> 는 없다 — 곡선은 진실이 아니라 EdgeCache 다(`edge_cache`, `rebuild_edge_cache` 가 버리고 재생).

// ─── 모션 (nacre-topo — Handle 이 필요하다) ────────────────────────────

/// 모션은 교환되지 않으므로 이력은 종류별 저장소가 아니라 **순서 있는 사슬 하나**다.
pub enum Motion {
    Rotate    { axis: Axis, pivot: [Rat; 3], angle: Angle },
    Translate { offset: [Rat; 3] },
    Mirror    { axis: Axis, offset: Rat },    // det = −1, 사슬 패리티
    /// 평면 자신의 프레임으로의 기저 변경 — 기울어진 면 위 스케치를 정확하게 만든다.
    /// 법선을 성분이 아니라 **이름으로** 든다 — 재귀는 프레임 없는 평면(세계)에서 끝난다.
    /// `flip` 이 노드를 온전한 프레임으로 만든다: 정준 계수에는 방향이 없으므로(첫 0 아닌 성분
    /// 양수) 한 평면의 두 면이 반대로 향할 수 있고, `flip` 이 계수의 부호를 뒤집어 감각까지 이름한다.
    /// 계수는 배치를 읽기 **전에** 뒤집는다 — `Canonical` 은 `+u` 를 법선에서 유도하므로 ŵ 와 û 가 함께
    /// 돌고 v̂ 는 남으며(벽에서 위는 위), `Named` 는 `+u` 가 진술이므로 ŵ 와 v̂ 가 돌고 û 가 남는다.
    /// 어느 쪽이든 반바퀴라 Proper(det = +1) — 패리티에 기여하지 않는다.
    Frame     { plane: Handle<Surface>, placement: FramePlacement, flip: bool },
}

/// 프레임의 원점과 u 방향 — 사용자가 정하면 값, 아니면 정준 유도. `PlanePoints` 와 같은 이분법.
pub enum FramePlacement {
    /// **기본값.** 평면의 순수 함수로 유도되는 정준 프레임 — 원점 = 세계 원점의 수선의 발
    /// (`p = (−d/n·n)·n`), 축 = **Arbitrary Axis**(DXF/AutoCAD 규약: `u = ẑ×n`, 법선이
    /// 정확히 수직이면 `ŷ×n` — 갈래를 정확히 가른다. `n` 은 프레임이 향하는 법선, 곧 `flip` 을
    /// 쓴 계수다 — 이름 없는 도로는 계수의 부호를 증명하지 못하므로 정준 이름의 부호에 묶인
    /// 규약은 세 도로가 함께 지킬 수 없다). 실현 시점에 필요한 정밀도로 계산되고
    /// 아무것도 저장하지 않는다 ⇒ 정준값이 `Rat` 을 넘치는 평면(분모 n·n 제곱, 코퍼스 1.6%)도
    /// **Through 평면**(거대 계수 위엔 i128 에 드는 유리수 점이 일반적으로 없다)도 같은 규약을
    /// 그대로 받는다 — 원점 위치가 오버플로 여부와 무관하게 사용자 기대대로다.
    /// 한 평면 위의 스케치들이 자동으로 한 노드를 공유한다.
    /// 유도 규약은 스펙으로 **동결**한다 — 바뀌면 기존 스케치가 조용히 돈다.
    /// 구현: 좁으면 유리수 `PlaneFrame`(비트 보존), 넘치면 판정층의
    /// `MoveNode::FrameWide`(BigInt 쌍둥이 — 넘침이 존재하지 않는 실현).
    Canonical,
    /// 호출자가 **명시적으로** 이름 붙인 값 — 정준 규약과 다른 프레임을 원할 때만
    /// (`world_zx` 의 `+u = ẑ` 처럼 유도값과 다른 규약, `through_points`·`with_origin`).
    /// **그 평면의 `points` 가 적힌 좌표계**의 3D 유리수다(모션 없으면 세계. 2D 가 아닌 이유:
    /// 평면 위 2D 좌표는 무리수인 u·v 축을 전제한다 — Frame 이 정의하려는 바로 그것).
    /// origin 은 평면 **위**의 점(정확 검사 = 방정식 대입, 성분 모양이 아니다 — C1),
    /// ref_dir 은 법선과 평행하지만 않으면 되고(단위 길이일 필요도 없다) 실현이 평면으로 사영한다.
    /// 저장된 값이 진실이라 replay 는 저장이 보장한다. 이름 있는 평면 전용.
    /// 검사 실패는 **구성 시점의 이름 붙은 거절**이다(`OriginNotOnPlane`·
    /// `RefDirParallelToNormal`·십진 창 밖 `FrameOutsideDecimalWindow`) — 조용히 Canonical 로
    /// 대체하지 않는다(사용자가 말한 곳과 다른 곳에 앉히는 것이 곧 조용히 틀림이다). Through 평면에
    /// Named 를 주면 유리수 점이 그 위에 정확히 놓일 수 없어 같은 검사에 자연히 걸린다 ⇒ 안내는
    /// «인자 없는 기본(Canonical)을 쓰라». 특정 정점에 원점을 앉히는 요구가 증명되면 값이 아니라
    /// **가리키는** 변종(`origin: Handle<Vertex>`)을 추가한다(규칙 1 의 연장 — 계획, todo 「어휘 밖의 요구가 오면」).
    Named { origin: [Rat; 3], ref_dir: [Rat; 3] },
}

pub struct MotionNode {
    pub motion: Motion,
    pub parent: Option<Handle<MotionNode>>,   // 숲 — 이력의 꼬리를 공유. interned(`motion_ids`).
}
// 모션이 기여하는 tol 은 적용 지점에 달려 있어 노드에 저장하지 않는다 — 판정이 뿌리까지 걸어 계산한다.
```

### 정점 — interning 하지 않는다

인접 면 집합은 **자란다**(나중 연산이 그 점을 지나는 면을 더 만든다) — 어떤 키를 잡아도 같은
점의 키가 바뀌므로 점의 동일성은 인접으로 키잡을 수 없다. 저장은 **정하는 세 면**만 들고(나머지
인접은 위상이 안다 — 적으면 두 번째 설명), 점의 동일성은 기하 질문으로 남는다(`merge_coincident`).
*"같으면 같은 핸들"* 이 구성 시점에 일어나는 것은 **곡면뿐**이다. 모서리·면도 interning 하지
않는다(두 면이 여러 토막에서 만나고, 한 평면 위에 떨어진 면이 둘일 수 있다 — 곡면 집합만으로는
«어느 것»을 못 말한다).

`validate` 의 일: 정의하는 세 면 위는 구성상 자명 — 남는 진짜 질문은 **정의하지 않는 나머지
인접 면 위에도 있나**(4-평면 동시성이 그 자리)다.

**저장은 «정하는 셋»이되, 이름은 인접 전부에서 정준으로 고른다.** 불리언 안에서
점을 부르는 이름은 그 점을 지나는 평면 집합의 함수 하나(`canonical_triple` — 사전식 최소의 독립 삼중)
이고, 피연산자 정점은 인접 면의 클래스 전부(위상이 안다)로 그 함수를 부른다. 저장된 정의 셋은 그중
하나의 철자일 뿐 다음 불리언의 이름을 정하지 않는다 — 정준은 불리언마다의 클래스 번호 공간에서 고르는
값이라 저장할 수 있는 것이 아니다.

---

## 스케치

```rust
// ─── 프로파일 (nacre-ops) — 프레임을 모르는 순수 2D. 구성 시점에 Rat 으로 확정된다. ──

pub struct Profile2d {                        // 한 재료 영역
    outer: Ring2d,
    holes: Vec<Ring2d>,
}
pub struct Ring2d {                           // 정점 + 변. edges[i] = vertices[i] → vertices[i+1 mod n]
    vertices: Vec<[Rat; 2]>,                  //   온전한 원 = 정점 1 + Arc 1 (정점 = 솔기)
    edges: Vec<Edge2d>,                       //   필드는 사설 — 손으로 지은 Arc 가 문을 우회하지 못한다
}
pub enum Edge2d {                             // 조각 «하나» — 순서 있는 고리라 시작=앞 꼭짓점 (nacre_geom::mixed)
    Line,
    Arc { center: [Rat; 2], r2: Rat, ccw: bool },      // r2 = 반지름의 «제곱» — 문에서 |start − center|² 로 한 번 유도해 저장
}
// **커널의 문은 «고리»다, «펜»이 아니다** — `Ring2d::new(vertices, edges)`(검증: 짝 맞음·
//    영길이·호의 r²·다음 정점이 원 위)·`Ring2d::circle`·`Ring2d::polygon_decimal` + 단계 문 `arc_turns(_rat)`
//    (변, 끝점)·`arc_to_rat`, 그리고 `from_paths(Vec<Ring2d>)`. 펜은 kit 의 것이고, 제약 해석기·외부 데이터가
//    와도 같은 문으로 온다 — 커널은 «어떻게 그렸나»를 모른다.
// 링 더미 → 짝수 깊이 = 재료(even-odd) → 섬마다 Profile2d 하나 — from_rings/from_paths.
// 정규형(`Ring2d::normalized`): 공선 중간점은 녹이고, 한 원·한 방향의 인접한 두 호는 한 호다.

// ─── 배치 — 어떤 평면 위 + 배치(기본은 정준 유도, 명시하면 값). ──

pub struct SketchFrame {                      // 필드는 사설, 생성자가 검증한다 — 리터럴이 검사를 돌아가지 못한다
    plane: Handle<Surface>,
    placement: FramePlacement,                // `canonical(plane)` = Canonical(유도) / `named(..)` = Named(값, 검사됨)
    flip: bool,                               // 사용자 몫이 아니다 — 생성자는 false, 연산이 쓰임새가
}                                             //   향할 쪽으로 진실에서 정한다(`frame_toward`)
// `SketchFrame::world(model, axis)` 는 이 struct 를 채우는 설탕이다 — 타입 변종이 아니다.
// (XY·YZ 는 Canonical 과 일치, ZX 만 +u = ẑ 가 유도값 −x̂ 과 달라 Named 를 쓴다 — 그 예외가 사는 곳은 이 함수 하나.)
```

| 규칙 | |
|---|---|
| **호의 진실은 끝점이다, 각도가 아니다** | 3D 경계 표현이 요구하는 것은 **정점과 원**이고 각도는 어디에도 안 쓰인다(`Arc` 에 각도 필드 없음). 끝점 표현이 더 넓다(무리수 도의 호도 끝점이 유리수면 정확 — 3-4-5). 스케치는 임의 각을 받는다: `arc_to_rat(center, start, end, ccw)` 이 끝점을 받고 «원 위인가»만 검사하지 각도를 제한하지 않는다; `arc_turns` 는 90° 배수 끝점을 `(x,y)→(−y,x)` 로 유도해 주는 **설탕**일 뿐이다. **제한은 «스케치»가 아니라 «돌출»에 있다**: `refuse_non_quarter_arcs`(`ops`)가 프리즘 빌더에서 비-사분 호를 `ArcSweepNotQuarterTurn` 으로 거절한다 — 빌더의 정확 감김(`Ring2d::winding_sign`)과 원통 옆면의 정확 표현이 사분 가족에서만 서기 때문이고, 스케치 입력의 한계가 아니다. 계획: 호의 중심·반지름을 정의 기반으로 옮기는 것과 곡선 마일스톤이 그 벽을 연다(todo) |
| **원통 원시체는 없다 — 원은 스케치다** | 원통은 원 프로파일의 Extrude 다. 직선–호 꼭짓점 = `Vertex::Pierce`(두 평면의 저장된 이름 × 원통, 점은 원통 둘레 순위로 고른다), 온전한 원 = `OnSeam` |
| **경계는 f64, 진실은 구성 시점에** | 공개 API 는 f64 그대로다. `Rat::from_decimal` 왕복을 `Profile2d` **생성자**에서 한다 — `check()`(자기교차·중첩·포함)가 진실 위에서 정확 술어로 돌고(`orient2d_rat`: narrow 우선 → BigInt 전역 부호, geom 의 `_rat` 워커 쌍둥이), 십진 창(1e38/1e-22) 밖 치수가 구성 시점의 이름 붙은 에러다(`ProfileOutsideDecimalWindow` / sketch 층 `OutsideDecimalWindow`). f64 부호 ≠ 십진 부호다(0.1·0.2·0.3 공선이 이진에선 굽음) — 그래서 check 가 진실 위에서 도는 것은 정확성의 문제이고, 방향은 항상 "작성자가 쓴 수가 이긴다". 생성자는 `check()` 를 돌리지 않는다(O(n²)) — 프로파일을 소비하는 **모든 연산이 먼저 돌린다**; `from_rings`·`from_paths` 가 검사하는 생성자다 |
| **공선 중간점은 생성자가 지운다** | 공선 정점의 양옆 벽은 한 평면 → 교차가 직선이라 정의 불가, 그리고 비-2-manifold 퇴화다. 제거는 형상 불변(무손실 정규화, 거절 아님 — 정리된 프로파일의 프리즘은 깨끗한 쌍둥이와 비트 동일, 모든 코너가 3-평면 정의 보유). 판정은 Rat 위 정확 orient2d, 제거 조건은 **엄격 내부**(중복점·스파이크는 생존해 각자의 이름 붙은 에러로 보고된다 — 경계 포함 판정을 재사용하면 작성자의 실수를 조용히 지운다) |
| **프레임은 평면에서 뜨지 않는다** | `Named` 의 `origin` 은 참조 평면 **위**의 점(정확 검사 — C1). `SketchFrame::named` 가 구성 시점에 검사한다 — `plane_residual_sign`(exact, **전역**: Narrow 는 Rat 대입 → 넘치면 BigInt, Wide 는 BigInt — fail-open 없음) ≠ 0 이면 `OriginNotOnPlane`, `WideFrame::named_of` None 이면 `RefDirParallelToNormal`(그 함수는 폭에 전역이라 None 은 평행/영벡터뿐), 십진 창 밖은 `FrameOutsideDecimalWindow`, 이름 없는 평면은 `PlaneWithoutExactForm`. 좌표계에 주의: «정의하려는 프레임의 (u,v,w)» 가 아니라 **그 평면의 `points` 가 적힌 좌표계**의 3D 점이다(모션 없으면 세계 — 상자 윗면 z=1 이면 origin 은 `[1,1,1]` 같은 점이지 셋째 성분 0 이 아니다). "평면 위" 는 성분이 아니라 **방정식 대입**(`a·x+b·y+c·z+d = 0`, 정확)으로 검사한다. `Canonical` 은 유도라 검사할 것이 없다. 면 위 스케치의 밑캡이 대상 면의 surface 핸들을 공유하는(flush 접촉 = 핸들 비교) 전제이기도 하다. 평면에서 d 떨어진 스케치가 필요하면 origin 을 띄우는 것이 아니라 **오프셋 평면**을 만든다 — 같은 프레임 안 `(0,0,d),(1,0,d),(0,1,d)` 유리수 세 점의 `Known` 평면(datum 가족). 담체가 실제 평면 핸들로 남아 이름·interning·flush 규칙이 그대로 성립한다 |
| **정확 유리수 프레임은 노드를 만들지 않는다** | 프레임의 세계 기저가 유리수이면(`exact_frame` — 평면 이름에서 `RatFrame::of_plane_frame` 으로 묻고, 평면에 모션이 있으면 그 사슬의 접기로 옮긴다(`RatFrame::carried` — 이동·사분각·축 거울은 접히고, 프레임 노드·사분각 밖 회전은 안 접혀 거절). 세계 축 평면 셋이 전부 여기 든다, `Canonical` 이든 ZX 의 `Named` 든) `Motion::Frame` 노드를 생성하지 않고 세계 유리수 산술로 구성한다(모션 = None) — 빈 사슬 지름길(`shared_base` 상쇄·축별 tol 0·정수 Shewchuk)이 그대로 산다. 거울이 든 사슬로 옮긴 기저는 왼손이라 `RatFrame` 이 손방향을 들고 `normal()` 은 옮긴 ŵ 를 돌려준다(쓸기·오프셋이 ŵ 쪽으로, 감김·호의 회전은 노드 도로가 사슬 parity 로 하듯 그 손방향으로). 변환의 옮겨 적기 규칙(「모션 사슬 — `MoveNode`」)과 같은 구성 시점 정규화다. 지반 잠금: 씨앗 평면의 canonical 프레임이 단위축 위에 **비트 정확**으로 실현된다(`a_seeded_planes_canonical_frame_is_the_world_basis_exactly`) — 생략은 중복 제거지 값 손실이 아니다. 실현한 축을 `from_decimal` 로 되들어올리는 게이트는 정규화가 필요한 축(`0.6000000000000001`)과 짧은 십진이 아닌 원점을 틀리게 읽는다(「가지 말 것」) — `SketchPlane::exact` 는 호출자가 *쓴* 축을 읽는 테스트 오라클로만 남는다. |
| **세계 축 평면 셋은 `Model::new()` 가 심는다** | 핸들이 결정적(0·1·2, 법선 축 순서 Z(XY)·X(YZ)·Y(ZX))이라 replay 가 자명하다. 캐시 방향은 **−축**(밑캡의 감각), `Default` 는 `new()` 로 위임(무씨앗 뒷문 없음). 씨앗은 정의상 영구 orphan — 정의·프레임만 가리키는 평면을 순회·직렬화가 따라가야 한다(todo 「정의만 가리키는 평면의 순회·직렬화」) |
| **스케치는 모델에 저장되지 않는다** | 연산 로그가 스케치의 진실이고, 모델에 남는 것은 벽·캡 평면(`Known` 점 + 프레임 모션)이다. 치수의 흐름: 변 `(x1,y1)→(x2,y2)` + 스윕 → 벽의 세 점 `(x1,y1,0),(x2,y2,0),(x1,y1,d)` — 산술 없이 유리수 그대로. 계수(이름)에만 곱이 있고, 그것이 진실이 점이어야 하는 이유다(계수 `i128` 적합 25.8% vs 점 100%) |
| **기본 프레임은 정준 유도(`Canonical`)다** | 원점 = 세계 원점의 수선의 발(부호·스케일 불변), 축 = Arbitrary Axis(DXF) — 사용자가 아무것도 안 정하면 이 규약이고, 규약은 스펙+테스트로 동결한다. 사용자가 명시한 프레임만 `Named` 로 **값을 저장**한다(유도값과 다른 규약을 유도로 흉내내면 규칙 변경에 조용히 돈다) |

연쇄가 닫힌다: 평면에 스케치 → 솔리드 → 그 면 위에 스케치 → … 매 단계가
`SketchFrame{plane: 이전 단계의 평면}` 으로 **사용자 조작당 1단**만 깊어지고(C6 — 불리언은
평면을 만들지 않는다, `a_boolean_mints_no_surface`), 좌표는 끝까지 유리수이며, 무리수는 프레임
실현(`1/√`)에만 산다.

---

## 이름과 interning

이름은 **하나만 저장**하고(곡면 캐시의 평면 팔 `SurfaceCache::Plane` 의 `name`, 문은 `Model::surface_name`), 좁은 형태는 저장이 아니라 **투영**이다:

```
Surface (진실 — 평면이면 점 셋, 원통이면 def, + 모션)
   │ 유도 (정준화: 분모 털기 → gcd → 부호 규약)  ※ 평면이고 유리수 닫힘일 때만
   ▼
PlaneName = Narrow([Rat;4]) | Wide([BigInt;4])    ← 저장은 이것 하나
   ├─ narrow() → Option<&[Rat;4]>    산술·프레임·지름길 — 사본이 아니라 빌려 읽는다
   └─ + 모션  → SurfaceKey::Name     interning 표의 키 — «합친다»

Surface (진실) ─── 이름이 유도 «안 될» 때 ──→ 진실 그대로가 키 — «안 합친다»
   혼합 프레임 Through 평면 → SurfaceKey::Through
   원통                       → SurfaceKey::Cylinder
```

```rust
/// 정준 계수 — 진실(점 셋)에서 유도한 «이름». 그릇이 임의정밀이라 유리수 닫힘인 평면은
/// **항상** 이름을 갖는다. 역할 둘: 동일성(`==`, enum 전체) · 산술 지름길(`narrow()`).
///
/// 정규화 불변식: **i128 에 들어가는 값은 반드시 `Narrow` 로 저장된다** — 유일한
/// 생성자(`plane_name_exact`, 점에서 유도)가 강제한다. 같은 값 = 같은 변종 = 같은 비트라
/// `Eq`/`Hash` 가 구조적으로 성립한다("정준 여부는 값이 말한다"와 같은 원칙). ~99% 가
/// Narrow(인라인, 힙 없음)이고 Wide 는 ~1%. `Wide` 는 **동일성만** 든다.
pub enum PlaneName {
    Narrow([Rat; 4]),
    Wide([BigInt; 4]),   // Box 없음 — BigInt≈32B 라 [BigInt;4]=128B=[Rat;4], enum 크기가
                         // 같아 Box 는 할당+간접만 더한다
}
impl PlaneName {
    /// `Rat` 계수가 필요한 읽기(프레임 유도·모션 운반·실현)가 읽는다.
    /// `None`(Wide)이면 그 읽기만 못 타고 넓은 그릇(정수)이나 일반 경로로 — 느릴 뿐 틀리지 않는다.
    pub fn narrow(&self) -> Option<&[Rat; 4]> { … }
    /// 어느 폭이든 정수 계수로 — 정수 부호 술어가 읽는 그릇 하나. 원시성은 «네» 계수에 대한
    /// 것이라 법선만 읽어 곱하는 소비자는 자기 gcd 로 먼저 나눠야 한다.
    pub fn coeff_ints(&self) -> [BigInt; 4] { … }
}

// ─── interning 표의 키 (nacre-topo) — 표 하나, 열쇠 셋 ──────────────────────

/// 같은 진술이라도 세계(모션 없음)와 모션 전 프레임은 다른 곡면이므로 **모션이 늘 함께** 든다.
/// 평면의 향(`sense`)은 **어느 열쇠에도 없다** — 한 평면을 `+n`/`−n` 으로 진술하면 한 핸들이고
/// 문이 `flipped` 로 알린다.
pub enum SurfaceKey {
    /// **유도된 정준형** — 다르게 진술해도 같은 평면이면 한 핸들(**기하 동일성**).
    Name(PlaneName, Option<Handle<MotionNode>>),
    /// 이름이 없는 `Through` 평면 — 정렬 삼중항(**문자 동일성**: 같은 진술 = 한 핸들).
    Through([Handle<Vertex>; 3], Option<Handle<MotionNode>>),
    /// 원통 — 진술 전체(`ref_dir` 포함). **일부러 보수적**이다(아래). 드물고 큰 쪽이라 상자에 든다
    /// (안 그러면 표의 모든 이름 열쇠가 그 크기를 차지한다).
    Cylinder(Box<CylinderDef>, Option<Handle<MotionNode>>),
}
```

### 표는 하나, 열쇠는 셋 — «합치는가, 아닌가»

| 열쇠 | 무엇 | 실제로 무엇인가 |
|---|---|---|
| `SurfaceKey::Name` | `(PlaneName, motion)` | **유도된 것** — 합친다 |
| `SurfaceKey::Through` | `([Handle<Vertex>;3], motion)` | `Surface::Plane{points: Through(vs), motion}` **그대로**(향 빼고) |
| `SurfaceKey::Cylinder` | `(CylinderDef, motion)` | `Surface::Cylinder{def, motion}` **그대로** |

뒤의 둘은 «어느 종류냐»로 갈려 있을 뿐 같은 관계다 — 열쇠 필드가 **진실의 필드와 같다**, 평면의 향만 빼고.
구분의 뜻은 «어느 종류»가 아니라 **«합치는가 아닌가»** 다. 향을 빼는 것은 평면의 정체에 방향이 없기
때문이다 — 같은 진술이 반대 향으로 오면 한 핸들이고 문이 `flipped` 로 알린다(「가지 말 것」 «향을 든
진실을 interning 열쇠로»).

원통의 **일부러 약한** 보장이 여기서 나온다: `ref_dir` 이 다르면 seam(θ=0 이음매)이 갈라지므로 —
통째 림의 seam 정점이 그 곡면을 담체로 인용하고 유리수 차트의 배제점이 거기다 — 기하가 같아도 **합치면 안 된다.** 기하 동일성은
술어가 물을 때마다 답한다(규칙 6). **원통에 정준 이름을 주지 않는 것은 «없어서»가 아니라 «결정»이다** —
원통에 4계수형이 없는 것은 사실이지만 지키는 것은 정책이다: 누군가 원통의 정준형을 만들어 합치기
시작하면 기하가 같은 원통이 **조용히 합쳐지고 seam 이 갈라진다.** 새 곡면 종류(구·원뿔)도 «합칠
것인가»를 먼저 정하고 나서야 이름을 얻는다. 이 결정이 막는 행은 없다 — 한 원통면을 두 솔리드가 공유하는
배치는 진술이 같아 한 핸들이어도 차트의 전제(«원통 클래스는 한 솔리드의 곡면»)로 거절되므로, 두 핸들로 갈라 둔
값은 거절의 자리만 바꾼다(인구는 todo 「한 원통면을 두 솔리드가 공유한다」). 방향은 todo 「곡면 둘 이상이 만나는
점과 seam 담체」가 든다.

### 열쇠는 진실이 고른다

열쇠는 **진술된 진실의 함수** 하나가 고른다(`Model::surface_key`): `Known` 점은 정준 이름(공선이면 열쇠 없음),
`Through` 는 세 정점이 한 프레임에서 이름을 주면 이름 아니면 진술 자체, 원통은 진술 자체 — 어느 문으로
들어왔는지가 아니다. 세 문(`push_plane`·`push_plane_through`·`push_cylinder`)은 진실을 진술하고 깔때기
하나(`intern`: 열쇠·이미 발급된 핸들과 `flipped` 회신·종류별 아레나 push(캐시와 이름이 함께 선다)·interning 표·계측)를 부른다. 열쇠가 읽는
것(점·정점 정의와 그 담체의 이름·모션 접힘)은 모두 한 번 쓰이고 열쇠보다 먼저 있으므로, 나중에 다시 물어도
push 때의 답과 같다. `push_plane` 의 `bool` 은 돌려준 곡면이 건넨 진술(점·향)과 **반대**를 향함을 말한다 — 정준형에는
방향이 없으므로(`[0,0,1,−3]` 과 `[0,0,−1,3]` 은 한 평면) 같은 평면이 핸들을 공유하면 방향을 어딘가에서
맞춰야 하고, 그 자리는 호출자다(면을 `Orientation::flipped()` 로 기록). 두 진실(점의 방향 × 향)의 비교이지
캐시의 비교가 아니다. 휴면 경로가 아니다 — 판 위에 앉은 보스는 한 평면을 양쪽에서 진술한다(census 의 수는
`push_plane` 의 doc).

**interning 표는 캐시다** — 아레나를 돌며 진실마다 `surface_key` 를 다시 매기면 통째로 재생된다(test-only
`push_plane_unregistered` 는 표를 **건너뛰므로** — 「한 평면을 두 핸들로」 픽스처가 그것을 필요로 한다 —
test-util 아래에서는 재유도가 표와 다르다).

**«같은 곡면 = 한 핸들»은 구조가 아니라 `validate` 가 지킨다**: 표를 쓰는 곳은 깔때기 하나지만, 종류별 raw
push 는 `pub(super)` 라 `model/` 안 어디서나 불릴 수 있다(test-only 문이 그렇게 표를 건너뛴다). `validate` 가
아레나의 모든 곡면에 `surface_key` 를 다시 매겨 같은 열쇠의 두 핸들을 `Violation::DuplicateSurface` 로 보고한다 —
표가 아니라 진실에서 다시 유도하므로 표를 건너뛴 곡면이 보이고, 그래서 표를 버리고 재생할 잠금이 따로 필요 없다.
**평면의 이름도 같은 열쇠로 잰다** — 저장된 이름이 열쇠의 이름(열쇠가 이름이 아니면 없음)과 다르면
`Violation::SurfaceNameMismatch` 다. `Through` 의 열쇠는 담체의 저장된 이름을 읽지만 담체는 늘 더 낮은 핸들이고
같은 고리가 그것도 재므로, 전부 같다는 것은 이름을 버리고 핸들 순서로 다시 지어도 같다는 것이다 — 재생 잠금이
따로 필요 없다. `push_plane_unregistered` 로 든 평면은 이름이 없어 이것에 걸린다.
아레나 전체인 이유는 참조 무결성 검사와 같다(곡면은 대체되지 않는 공유 정의다). 비용은 아레나 크기에 비례한다 —
사분각 회전 2,000번 모델(곡면 12,022)에서 `validate` 19 ms 중 16 ms, 회전 fin fold 80(곡면 325)에서 4%; 스위트
벽시계는 그대로였고, 앱(kit·playground)은 `validate` 를 부르지 않는다.

- 문을 넘나드는 경우는 잠겨 있다: `a_through_plane_and_a_known_plane_that_are_one_plane_share_a_handle` —
  `Through` 평면과 `Known` 평면이 같은 평면이면 한 핸들이다(변종은 그 평면의 존재가 무엇에 근거하는지를
  적지, 누가 먼저 물었는지를 적지 않는다).
- 유리수 점 셋으로 말한 평면은 계수가 **반드시** 유리수라 항상 이름이 있고(⇒ `SurfaceKey::Name`), 진술 키로는
  이름이 계산 안 되는 진술만 간다 — 그러나 **한 평면이 두 표에 걸칠 수는 있다**: 혼합 프레임 `Through` 가 이름
  있는 평면(예: `z = 0`)을 다시 진술하면 그 진술은 진술 키로 간다(아래 «핸들 `!=`», 규칙 6).
- 진실 안에 `f64` 는 **0개**(전부 유리수·핸들)다. `CylinderDef`·`MotionNode`·`PlaneName` 은 `Eq, Hash` 를
  파생하고, `Surface`·`PlanePoints` 는 `PartialEq` 만 파생한다.
- **정준화**: 분모 털기(lcm) → 내용(gcd) 나누기 → 부호 규약(첫 0 아닌 성분 양수). 정준 여부는
  플래그가 아니라 **값이 말한다**(정수인가·서로소인가·부호 규약인가).
- **생산자는 점만 진술한다**(`push_plane(cache, points, motion)`) — 이름은 커널이 유도하므로
  (`plane_name_exact`) *"한 평면을 두 가지로 진술한다"* 가 표현 불가능하다. 이름은 **사용자가 쓴 치수에서**
  짓지 f64 계수에서 들어올리지 않는다(들어올림은 반올림을 보존해 같은 평면의 두 벡터가 계속 다르다).
- 핸들 `==` 는 «같은 평면»을 증명하지만 **핸들 `!=` 는 아무것도 증명하지 않는다** — 이름 없는 평면
  (혼합 프레임 `Through`)은 진술로만 intern 되고, 모션이 다른 두 핸들이 한 세계 평면일 수 있으므로
  기하 술어가 뒤를 받친다. **세계 이름**(`world_plane_name`)은 둘 다 있으면 `==`·`!=` 둘 다 증명이다 —
  평면 클래스 병합이 그것을 먼저 읽는다. interning 은
  결정적이다(먼저 넣은 쪽이 이긴다 — interned 평면은 첫 pusher 의 세 점과 캐시를 유지한다).
- **interning 의 값은 속도가 아니라 동일성이다**(속도 기여 3.7%) — 같은 평면이 두 핸들로 갈리지 않는 것.
- **큰 정수 계수의 정확 부호**: `2^53` 은 f64 *하나*의 한계일 뿐 — 이름의 정수는 넓어도 정확하다.
  판정의 이름 정수 도로가 정수의 최근접 f64 로 필터를 돌린 뒤 못 가른 것을 정수로 답한다 — 정수가 하나하나 f64
  그 자체인 이름이든 아니든 한 길이다(「판정 = 필터 + 계산된 상승 + 실현 캐시」의 «이름 정수 도로»). 정규화
  실패한 유리수는 술어 안에서 분모를 턴다.

---

## 캐시

```rust
// crates/nacre-topo/src/lib.rs
pub struct Model {
    // 진실 — append-only. 저장소는 모두 비공개다: 밖은 「문의 이름」 절의 문으로만 읽고,
    // push 는 생성자(`push_plane`·`push_cylinder`·`push_vertex`·`push_edge`·…) 경유만 한다.
    // 전량 순회는 없다: 아레나엔 superseded 도 있어 소비자는 live face 를 걷는다.
    surfaces: Store<Surface>,
    vertices: Store<Vertex>,
    edges:    Store<Edge>,
    faces:    Store<Face>,
    shells:   Store<Shell>,
    solids:   Store<Solid>,
    motions:  Store<MotionNode>,              // interned — 쓰기는 `push_motion` 하나

    // 루트 — 살아 있는 솔리드 핸들 목록(도달 집합의 뿌리). 저장소가 아니라 «무엇이 현재 모델인가».
    // 읽기는 `live_solids()`, 쓰기는 `push_solid`·`supersede_live`·`restore_live`.
    live_solids: Vec<Handle<Solid>>,
    world_planes: Vec<Handle<Surface>>,       // `Model::new` 가 심은 세계 평면 셋(길이 늘 3)

    // 캐시 — 핸들 인덱스 병렬, 통째로 버리고 재생 가능. 비공개다: 밖에서 벡터를 만지면
    // index-parallel 불변식을 아무도 못 지킨다. `Store` 가 아니라 `Vec` 인 이유: `Store` 는
    // append-only 로 봉인돼 있어 더 높은 정밀도로 고쳐 쓸 수 없고, 고쳐 쓸 수 없으면 캐시가 아니다.
    vertex_cache:  Vec<PointCache>,
    surface_cache: Vec<SurfaceCache>,         // 실현. 아레나가 진실을 든다
    edge_cache:    Vec<EdgeCache>,            // 평가 가능한 곡선
    given_by_pair: HashMap<[Handle<Surface>; 2], EdgeGiven>, // 담체 쌍 → push 하는 쪽이 실현한 조각
    prefix_hp: HashMap<PrefixKey, PrefixValue>, // 모션 사슬 접두의 고정밀 실현 메모(실현 가속기)
    chain_affines: Vec<Option<AxisAffine>>,      // 모션 핸들 인덱스 평행 — 노드가 태어날 때 사슬 전체를 접은 유리수 사상

    // 동일성 — 「같은 곡면 ⇒ 같은 핸들」이 되게(그래야 동일성이 정수 비교다). 곡면의 표는 하나, 열쇠는 셋
    // (`SurfaceKey::{Name, Through, Cylinder}` — 이름 있는 평면 · 이름 없는 `Through` · 일부러 보수적인 원통).
    interned:            HashMap<SurfaceKey,  Handle<Surface>>,
    motion_ids:          HashMap<MotionNode,  Handle<MotionNode>>,
}
// interning 표는 순회하지 않는다 — `HashMap` 의 순서가 결과에 닿으면 안 된다.

// 수치 층 — 「값 + 그 값이 갇힌 오차」. 같은 것의 다른 정밀도는 `Hp` 접두사 하나로만 다르다.
// crates/nacre-exact/src
pub struct Mag { m: f64, e: i64 }                                  // f64 밖 범위의 보수적 크기(0에서 먼 쪽)

// ── 원자 — 값 하나 + 그 반경. «참값 ∈ value ± error» 가 증명된 것만 담는다. ──
pub struct Bounded   { pub value: f64,      pub error: f64 }
pub struct HpBounded { pub value: BigFloat, pub error: Mag }

// ── 경계 지어진 점 — 원자 셋(축별). `WitnessPoint` 의 `coord_hp`(`HpBounded` 셋)와 `coord_iv`(`Bounded` 셋)가
//   이 모양이다. 값과 오차를 «묶어» 든다(따로 든 배열 둘이 아니라). 그래야 «참값 ∈ value±error» 가 구조로 서고
//   값만 옮기고 오차를 안 옮기는 desync 가 불가능하다.

// ── 정점의 캐시 ──
pub enum PointCache {
    Realized   { coord: Point3, bound: [Mag; 3] },                //   정의에서 실현 — 참값 ∈ coord ± bound
    Ceiling    { coord: Point3 },                                 //   싼 도로가 멈췄다 — 비싼 문이 답할 수 있다
    Unrealized { coord: Point3 },                                 //   실현이 뒤에 없다 — 도로가 없거나, 안 물었거나
}
// 연산이 만든 정점은 push 시점에 실현되므로(`push_vertex_realized`) 거의 전부 `Realized` 다
//   (census 결과 정점 5,082/5,082). 모션 이력이 갈린 담체의 코너는 담체마다의 증인 삼각형으로 실현된다
//   (담체마다 `surface_witness_triangle` → `nacre_judge::plane_hp`, 세 평면의 교점 `meet_hp`, 한 묶음 동안 평면은
//   `PlaneMemo` 가 한 번만 실현한다); 그 도로는 증인 점 아홉을 재생하므로 비용 한계가
//   그 합을 센다.
// **`Ceiling` 의 천장에는 벽이 둘이고, 둘은 같은 말을 한다** — (a) 사다리의 두 단(128·256비트)이 못 정했다 — 0 이 아닌
//   좌표는 비용 한계 안에서 256 이 정하므로(한 회전에 약 1비트, 192 + 53 < 256) 인구는 참값 0 이 사분각 밖
//   회전 수십 번보다 깊이 묻힌 좌표다(일치 정밀도 2⁻¹⁸⁰ 까지 반경이 줄지 않는다 — 회전당 약 1비트에서 예상 76번,
//   37° 왕복 픽스처에서 60번은 결정·80번은 못 함),
//   (b) 재생해야 할 이력이 비용 한계(`CHEAP_REPLAY_COST_CAP` = 192)보다 깊어 싼 도로가 아예 안 걸었다(접기가
//   답하는 사슬은 깊이와 무관하게 읽는다). 읽는 쪽에도
//   되찾는 문(`nacre_ops::refine_caches`)에도 뜻이 같다: «더 물으면 답이 있다» ⇒ 한 변종이 맞다.
//   `Unrealized` 는 그 반대 — 도로가 없다(증인 삼각형이 없는 담체 `NoMeet`, 벽을 원통 틀로 읽을 수 없는 관통점 `NoCurvedPoint`),
//   또는 정의가 점을 이름하지 않는다(세 평면이 한 점에서 만나지 않는 `NoPoint` — validate 가 보고하는 생산자
//   결함). 그 갈림은 사람이 읽을 것을
//   위해 있다: «이 모델이 비싸다»와 «실현할 수 없다»는 다른 보고다. 변종을 고르는 것은 호출자가
//   아니라 깔때기(`push_vertex_realized`)이고, 철저한 `match` 라 새 이유가 옛 이름에 조용히 접히지 않는다.
// **어디에도 잔차를 싣지 않는다.** 잔차는 담체까지의 거리 «하나»이고 좌표가 참값에서 얼마나인지를
//   말하지 않으므로 «증명된 것»이 아니다 — 배열의 seam 표도 캐시(`SeamVertex.cache`)를 들고, 자기접촉 체는
//   `Realized` 의 경계로 상자를 부풀리며 경계가 없는 변종은 무한 상자로 읽는다(후보를 버리지 않고 정확 판정에
//   넘긴다). ⇒ validate 는 모든 정점에 구성 ε 을 쓴다. 대가: 잔차를 들던 인구의 검사가 1.07e-14 에서 1e-9 로
//   다섯 자릿수 느슨하다.
// **«잰 것이 있는지»를 변종으로 말한다** — 경계가 없는 변종엔 오차 자리가 아예 없다(있을 수 없는 것을
//   표현하지 않는다).
// **변종 이름은 «출처»(구성/발견)가 아니라 «캐시가 아는 것»이다.** 출처로 이름하면 기록된 이동을 받은
//   발견 정점을 «구성»이라 철자하게 된다 — 침묵은 거짓말이 아니지만 이름은 거짓을 말한다.
// **`bound` 는 `[Bounded; 3]` 의 오차 자리와 같은 성질이지 잔차가 아니다.** 잔차(담체 거리의
//   max)를 «참값 ∈ 값±오차» 자리에 넣으면 타입이 증명 안 된 포함을 주장한다(거의 퇴화한
//   담체 교차는 모든 담체에 가까우면서 정확한 코너에서 멀 수 있다).
// 고정밀 점 실현은 이름 붙은 타입이 아니다 — 연산 하나짜리 실현이라 «Cache» 가 아니고, `[HpBounded; 3]`
//   (메모는 `(usize, [HpBounded; 3])` = 정밀도와 값)으로 든다.

// ── 경계 점 + 정의 — 판정 전용(nacre-judge). 경계 지어진 점의 상위집합이지 중복이 아니다. ──
//   `WitnessPoint = [Bounded;3] + 정의(base:[Rat;3]·chain) + hp(메모)`. 그 «정의» 가 정확 단계의
//   입력이라 캐시(`PointCache`)와 «같은 것»이 아니다 — 이 차이가 진실/캐시 경계 그 자체다.
//   그래서 이것은 캐시로 접지 않는다(접으면 경계가 지워진다). 정의는 「판정」 절.

// ── 곡면의 캐시 — 종류마다 팔 하나, 그 종류의 캐시가 아는 것만 든다. `PointCache` 와 같은 모양이다. ──
pub enum SurfaceCache {
    Plane    { plane: nacre_geom::Plane, name: Option<PlaneName>, standing: CacheStanding },
    Cylinder { cylinder: nacre_geom::Cylinder, standing: CacheStanding },
}   // .plane() · .cylinder() · .name() · .standing() · .to_geom() -> nacre_geom::Surface
pub enum CacheStanding { Unrealized, Ceiling, Realized }          // 캐시가 값에 대해 아는 것 — 순서가 앎이고, 오르기만 한다
//   팔마다 제 종류의 조각을 든다 — geom 의 `Surface` 를 통째로 들면 평면 팔에 원통이 들어갈 수 있다. 기하의
//   메서드(`distance`·`normal_at`·`translated` …)는 `nacre_geom::Surface` 에 남고, 다형 질의는 `to_geom()` 이
//   그 조각에서 새로 짓는 값에 묻는다(문이 geom 의 값을 돌려주지 않는 것은 그것이 힙을 쓰는 변종의 자리를
//   두려고 `Copy` 가 아니어서다). 표지는 `PointCache` 의 세 변종과 같은 말이다(경계를 들 것이 없어 표지만):
//   `Unrealized` 는 아무도 실현하지 않았다(안 물었거나 도로가 없다 — 읽을 수 없는 사슬), `Ceiling` 은 싼 도로가 멈췄다(비용 한계를 넘는 사슬·두 단에서 미결정 — 값을 치르는 문이 올린다),
//   `Realized` 는 모든 부분이 진실의 최근접 f64 다. 둘째 쓰기 문 `refine_surface_cache` 는 표지를 내리지 못한다.
//   평면의 f64 실현은 문이 진실에서 유도한다(`derive_surface_cache`) — 세계 이름이 있으면 단위 법선은 이름 ×
//   향의 정확 반올림(`nacre_exact::unit_vector_f64`, 참 0 성분은 `+0.0`), 앵커는 진실의 첫 점(`Through` 는
//   꼭짓점들이 만나는 프레임의 첫 만남을 평면의 사슬로 옮긴 점; 진실이 그 점을 유리수로 못 놓는 넘침·넓은 만남은
//   받은 앵커). 문이 유도하지 못한 것은 ops 깔때기(`push_plane_realized`·`push_plane_through_realized`)가 push 한
//   뒤 진실에서 실현해 올린다(`raise_surface_cache` → `Model::refine_surface_cache`) — 문이 먼저 정하므로 «진실이 여기서
//   답하나»는 topo 한 곳에 있다. 앵커는 첫 점을 사슬로 재생해(넓은 만남이면 첫 꼭짓점 자신의 실현), 법선은 모션
//   전 이름 법선을 사슬로 옮겨(이름이 `Rat` 보다 넓으면 진술의 세 점을 재생해 그들이 펼치는 법선으로 — 판정의
//   증인이 프레임 탐침이면 그 순서는 진술의 것이 아니라 쓰지 않는다), 둘 다 같은 읽기 규칙으로. 답하지 않는 반쪽은
//   생산자 값이고 표지가 그것을 말한다: 재생이 비용 한계를 넘는 사슬·256 에서 미결정은 `Ceiling`. 이름이 없는
//   (혼합 프레임) `Through` 는 이름 법선이 없으니 증인 삼각형(`plane_witness` — 판정 표가 감는 그 삼각형)이 펼치는
//   법선을 그 삼각형과 함께 오는 회전의 향(`turn`)으로 세우고, 앵커는 정의 점 0 의 실현을 평면 사슬로 옮긴 점이다;
//   비용 한계는 그 평면이 재생하는 정점 담체의 사슬까지 센다(`Model::plane_replay_depth`). 원통도 같은 두 단이다 — 세계에 진술되면 문이 원점·단위 축·단위 `ref_dir`·반지름을
//   정확 반올림하고(`cyl_unit_frame_f64`), 아니면 깔때기(`push_cylinder_realized`)가 사슬을 재생해 같은 넷을 실현하며,
//   비용 한계를 넘는 사슬은 생산자 값에 `Ceiling` 이다. 사설
//   `push_plane_raw`·`push_cylinder_raw` 가 진실과 캐시를 같은 자리에서 채운다.
//   **평면의 향은 진실이다**(`Surface::Plane.sense`): 면의 바깥 = 평면의 향 × `Face.orientation`, 둘 다
//   진실이다. 생산자는 향을 정확히 진술하고(캐시에서 읽지 않는다 — `f64 → 진실` 이다), `flipped` 는 두
//   진실의 비교다. 캐시의 향은 진실을 따른다: 세계 이름이 있으면 법선 자체가 이름 × 향이고,
//   이름 없는 평면은 증인 삼각형과 그 회전의 향(`turn`)에서 실현된다 — 문이 먼저 점들의 f64 외적으로 생산자
//   수치를 맞추고(`align_cache_sense`, 움직이지 않은 것만), 그 수치는 실현이 멈춘 곳(`Ceiling`)에만 남는다. 향
//   잠금(`audit_plane_senses` — census, 그리고 이름 없는 평면을 드는 스위트의
//   `every_plane_faces_the_way_its_truth_says`)이 결과를 지킨다. 원통은 기준이
//   진실 안에 있다(축에서 바깥으로). `figure` 인자가 드는 것은 실현값(앵커·단위 법선)뿐이다.
//   **향을 묻는 판정은 전부 이 진실을 읽는다**: 평면이 보는 쪽 = `sense` × 점들의 방향(이름과의 관계는
//   `Model::plane_name_sense`·`world_plane_name_sense`), 증인 삼각형의 순서 = `sense × orientation`, 두 면의
//   관계 = `face_facing`. 캐시 법선은 어느 방향 판정에도 들지 않는다. (「가지 말 것」의 «평면의 방향을 점
//   순서에서 뽑기» 와 어긋나지 않는다 — 그 행이 버린 것은 점 순서 **단독**이고, 여기서는 `sense` 가 그 관계를
//   진술한다.)
//   평면의 이름(`PlaneName`)은 `plane_name_exact(세 점)` / `plane_name_through(세 정점)` 로 언제든 다시
//   나오므로 캐시이고, 평면 팔이 든다 — 원통 팔에는 그 자리가 없다. `None` 은 세 정점을 한 프레임이 못 놓는
//   `Through`(와 시험 문)뿐이고, `validate` 가 진실에서 다시 유도해 잰다(`SurfaceNameMismatch`). 칸을 짓는 것은
//   raw push 둘뿐이고 다시 쓰는 자리는 제자리에서 고친다(`SurfaceCache::refigure` — 이름이 인자가 아니므로
//   기하를 고쳐 쓰다 이름을 떨어뜨릴 수 없다).

pub struct EdgeCache     { curve: Curve }                          // 평가 가능한 담체 곡선
//   담체와 끝점이 곡선을 정한다(`derive_edge_curve`). `rebuild_edge_cache` 가 통째로 버리고 재생한다.
//   topo 가 진실에서 못 내는 조각은 유도가 push 하는 쪽에 묻는다(`EdgeGiven` — 세계 이름 없는 평면이 낀
//   직선의 방향, 세계 진술 없는 원통의 림 중심). ops 의 답은 담체 쌍마다 `given_by_pair` 에 남아 같은 쌍의
//   다음 간선이 읽고, 재생은 그것도 버려 다시 묻는다.
// **캐시는 «실현값 (+ 필요하면 그 오차)»** 다: 캐시가 드는 오차는 증명된 경계뿐이고
//   (`PointCache::Realized{coord, bound}`), 오차 필드는 «소비자가 있으면» 붙지 대칭으로 붙지 않는다 —
//   `EdgeCache`·`SurfaceCache` 에는 오차 필드가 없다.
```

**두 «오차»를 섞지 말 것.**

f64 는 둘이고 둘은 만나지 않는다 — 모델의 캐시(표시·tess·내보내기용 실현값)와 판정의 작업
사본(연산 하나 동안만 산다). 표의 두 줄은 한 사다리의 두 단이 아니라 두 수명이다.

| 무엇 | 성질 | 어디 사나 |
|---|---|---|
| `Plane::distance_eps(p)` | **이 거리 계산**의 f64 반올림, `3ε·Σ\|pᵢ−oᵢ\|` — `p` 에 의존(두 연산자). 원통판은 `4ε·(축거리+반지름)` | 호출할 때마다 계산, **저장 불가** |
| `PointCache::Realized.bound` | 실현의 축별 경계 — 좌표가 참값에서 얼마나(결정된 좌표는 반 ulp, 일치 정밀도로 0 이 된 좌표는 그 구간), **잔차가 아니다** | 모델과 같이 |

「곡면당 저장하는 tol 은 틀린 양」은 첫째 줄에만 맞다: 거리 계산의 반올림은 `p` 에 달려 있어 곡면의
성질이 아니다 — 평면 위에 정확히 있는 점도 0이 아닌 잔차를 낼 수 있고, 그 잔차가 주장된 허용오차와
**정확히 같아** 비교가 엄격하다는 것 하나로 통과하던 사례가 있었으므로 `distance_eps` 가 그 항을 따로 든다.
`Plane` 은 **원점 + 단위 법선**만 든다(`from_point_normal` 은 한 번 정규화하고, `from_point_unit_normal` 은 받은 단위 법선을 비트 그대로 싣는다 — 단위 벡터를 다시 정규화하면 비트가 움직인다). 판정은 이것을 읽지 않는다.

좌표가 결과인 STEP 출력은 정의에서 실현한 값을 정확 반올림해 낸다 — 정점의 캐시는 태어날 때부터 그
실현값이다(`push_vertex_realized`: 사다리의 낮은 단이 정하는 좌표는 어느 높은 단이 댈 답과 같아
`realize_vertex(…, NearestF64)` 와 비트 동일하고, `nacre-step` 은 `vertex_point` 를 읽는다). 자릿수를 더
원하는 호출자의 문은 `realize_vertex_decimal` 이다 — 캐시를 더 길게 찍으면 점이 아니라 반올림이 찍힌다.

| | 모델 캐시 (`PointCache`·`SurfaceCache`·`EdgeCache`) | 실현 (`WitnessPoint.coord_iv`·`WorkingCyl.realized`) |
|---|---|---|
| 누가 만드나 | 정점: push 깔때기가 **정의에서** 실현한다(실현이 거절하면 구성 자리가 계산한 f64 가 선다). 평면: 세계 이름이 있으면 문이 진실에서 앵커와 단위 법선(정확 반올림)을 유도하고, 문이 유도하지 못한 반쪽은 ops 깔때기가 push 한 뒤 진실에서 실현해 올린다. 답하지 않는 반쪽엔 생산자가 계산한 f64 가 서고, 표지(`CacheStanding`)가 그것을 말한다. 원통: 세계에 진술되면 문이 진술에서, 사슬이 안 접히면 ops 깔때기가 사슬에서 원점·단위 축·단위 `ref_dir`·반지름을 정확 반올림한다(비용 한계를 넘는 사슬은 생산자 값에 `Ceiling`). 간선: 담체·끝점에서 유도 — 원은 원통 캐시의 프레임을 비트 그대로 싣고, 중심은 진술된 원통이면 진실의 만남(정확 반올림), 아니면 ops 깔때기가 사슬로 재생한 축과 캡 평면의 만남을 실현한 것(그 길이 답하지 않을 때만 캐시끼리의 f64 만남); 직선은 시작 정점을 지나 세계 이름 있는 두 평면이면 이름 법선의 정확한 외적의 최근접 단위, seam·룰링이면 원통 캐시의 축, 평면 하나라도 세계 이름이 없으면 ops 깔때기가 두 평면의 계수(세계 이름, 아니면 증인 삼각형 — `PlaneMemo`)의 외적을 실현한 최근접 단위를 따르고, 그 길이 답하지 않을 때만(비용 한계를 넘는 사슬·증인 없는 평면) 끝점의 차를 따른다 | 판정이 **정의에서** 정밀도를 불러 만든다 |
| 오차 | 정점만 축별 경계 셋(`Realized` 변종일 때). 곡면·간선엔 없다 | 축별 경계 셋 — «참값 ∈ 값 ± 경계» 증명됨 |
| 누가 읽나 | 테셀레이션·STEP·물성·validate · ops 판정 표의 seam 표 — 실현 도로가 없는 점의 구성 수치 | 판정 1단(f64 필터) → 2단 정수 → 3단 상승 |
| 수명 | 모델과 함께 | **연산 하나** |
| 판정 경로에 | `nacre-judge` 는 읽을 수 없다(`nacre-topo` 에 의존하지 않는다). 정확 지름길은 평면의 이름을 읽는다(「판정」의 `WorkingPlane`). ops 의 향 판정도 캐시를 읽지 않는다(`sense`·`WorldName`·`face_facing`) — 판정 표가 평면 캐시를 읽는 것은 seam 표에서 실현 도로가 없는 점의 구성 수치뿐이다 | 그 자체 |


| 캐시 | 키 | **수명** | |
|---|---|---|---|
| 모델 캐시 (f64) | 핸들 인덱스 (밀집) | **모델과 같이** | push 시점에 채운다. 판정의 정확 지름길은 읽지 않는다(평면의 이름을 읽는다) |
| 고정밀 메모 | (정의, 정밀도) (희소) | **연산 하나** | 판정이 실제로 만든 점에만 — 위 «실현». `trial_bound` 는 일부러 이 메모를 우회한다 |
| 각도·√ 표 | `(Angle, prec)` / `(Rat, prec)` | 프로세스(스레드) | 값이 키의 순수 함수 — `cos 37°`·`1/√(n·n)` 은 어디서나 같다. 정확 반올림이라 실현이 유일하고 tol 도 유일 |
| 평면 계수 메모 (`PlaneMemo`) | (곡면 핸들, 정밀도) | **연산 하나**(불리언·변환·프리즘) | 모션 이력이 갈린 코너의 도로(seam 표)와 간선의 직선 방향이 담체 평면의 고정밀 계수를 한 번만 실현한다 — 값이 키의 순수 함수라 답을 바꾸지 못한다. 없으면 회전 핀 80 fold 가 1.1 → 1.8 s |
| 간선의 주어진 조각 (`given_by_pair`) | 담체 쌍 | **모델과 같이** | push 하는 쪽이 실현한 직선 방향·림 중심을 쌍마다 남긴다 — 두 진실의 함수(정확 반올림이라 유일)이고 진실은 append-only 라 낡지 않는다. 불리언은 결과의 간선을 매번 다시 push 하므로(회전 핀 80 fold: 이름 없는 평면이 낀 간선 45,888, 걸음마다 새 쌍은 열 남짓) 없으면 그 fold 가 1.33 → 1.94 s |
| 사슬 접기 | 모션 핸들 (밀집) | **모델과 같이** | 노드가 태어날 때(`push_motion`) 부모의 접기에 자기 하나를 합성 — 이동·사분각 회전·축 거울은 «부호 있는 축 치환 + 유리수 이동»(`AxisAffine`)으로 닫힌다. `Frame`·사분각 밖 회전이 끼면 `None`. 질의마다 사슬을 걸으면 깊이의 제곱이다 |

여섯은 인덱스 공간이나 수명이 달라 합치지 않는다.

---

## 판정 (연산 동안만 산다)

판정 타입은 두 크레이트에 나뉘어 산다. **`nacre-judge`** 는 정의를 든 점(`WitnessPoint`)·모션의
펼침(`MoveNode`)·증명 기준(`ProofStandard`)·결과(`Decision`)와, 평면 표에 던지는 질문(`Judge<W: Witness>`)을
갖고, 두 평면이 원통과 만나는 점(`PiercePoint` — `Judge::pierce` 가 이름에서 푼다)과 그 점에 묻는 질문도 같은 표 위에서 갖는다. **`nacre-ops`** 는 연산 하나 동안 모델에서 그 표를 짓는다(`pub(crate)` `FaceInfo`·`WorkingPlane`·
`WorkingCyl`) — `nacre-judge` 는 `nacre-topo` 를 모르므로 핸들은 이 경계에서 해소된다.

이름 규칙: (1) **`Working*` = 모델 타입의 판정층 쌍둥이이자 표의 뿌리**(수명 = 연산 하나 — `Model` 에
`Working…` 이 담기면 읽는 즉시 이상해 보여야 한다). 접두사는 수명 표식이 아니라 **뿌리 표식**이다 —
부품(`WitnessPoint`·`MoveNode`)은 `Working*` 컨테이너 안에 살며 수명을 상속하므로 접두사를 반복하지
않고, 이름은 역할을 말한다. (2) 평면 정의를 이루는 유리수 점은 **증인 점 `WitnessPoint`** 다.

**면과 평면 클래스는 두 타입이다**(`FaceInfo` / `WorkingPlane`). 이 분리가 하중을 진다: `orient_sign`
(이 면의 저장 법선이 바깥을 향하나) 과 `frame_sign`(그 클래스 뿌리의 저장 법선이 뿌리 면의 바깥과
맞나)은 다른 인덱스 공간의 다른 사실이고, 한 이름으로 합치면 규약을 단언할 수 없는 자리가 생긴다.

```rust
// crates/nacre-ops/src/planes/frames.rs — 판정용 평면 클래스. `Surface::Plane` 의 쌍둥이.
pub(crate) struct WorkingPlane {
    pub(crate) plane: Plane,                          // f64 실현
    pub(crate) base_rat:  Option<[Rat; 4]>,           // 클래스 뿌리의 정확 계수(모션 이전 프레임)
    pub(crate) world: Option<WorldName>,              // 같은 평면의 세계 이름과 그 향 — 판정층이 좁은 투영(`Judge::world_rat` — 원통 도로가 세계 축과 비교하는 유일한 서술)과 향을 곱한 방향(`Judge::stored_world_rat`)으로 읽는다
    pub(crate) surf: Handle<Surface>,                 // 대표 곡면 — 결과 정점의 담체로 기록된다
    pub(crate) tri_pt3: [WitnessPoint; 3],            // 같은 증인의 정확한 정의 — 모든 술어가 빌린다
    pub(crate) rotated: bool,                         // 술어 경로 선택 신호. `tri_pt3` 와 함께 복사돼 어긋날 수 없다
    pub(crate) frame_sign: i8,                        // 라벨 프레임이자 차트의 프레임(n_out = frame_sign · n_P)
    pub(crate) chain_id: u64,                         // 모션 잎의 정체(0 = 없음) — 공유 모션 지름길의 열쇠
    pub(crate) name_ints: Option<NameInts>,           // 저장 방향으로 접힌 이름 정수 + 그 최근접 f64 철자
                                                      //   — 판정자가 이름 있는 평면을 읽는 유일한 정확 서술(이름 정수 도로)
}

// 판정용 원통 클래스 — 결정은 전부 `def`(정확)를 읽고, `realized` 는 «측정»만 읽는다.
pub(crate) struct WorkingCyl {
    pub(crate) surf: Handle<Surface>,
    pub(crate) def: nacre_topo::CylinderDef,
    pub(crate) realized: nacre_geom::Cylinder,        // 모델의 캐시가 아니다 — 연산 하나를 산다
    pub(crate) owner: SolidSide,                      // 쌍 루프는 주인이 다른 쌍만 묻는다
    pub(crate) shared: Vec<SharedRuling>,             // 옆면 위에 놓인 두 클래스의 교선 — 게이트의 기록
    pub(crate) records: Vec<(usize, CylRecord)>,      // 게이트가 배열에 넘기는 클래스별 기록 — 절단(`Cuts`), 비스듬한 단면(`Sections`), 접함(`Touches`), 또는 이름 없는 클래스의 말 없음(쪽, 면이 떨어짐, 단면, 스침)
    pub(crate) boxes: Vec<[[Rat; 3]; 2]>,             // 옆면마다 호를 아는 상자 — 이름 없는 클래스가 있을 때만
}
```

**판정 표의 증인(`Witness`)은 진실만 든다** — 세계 이름(`world_name`)과 정의 점(`tri_pt3`). 면 꼭짓점의 좌표
캐시는 표 행에도 트레이트에도 없어 판정이 읽을 수 없다. 증인 삼각형의 순서는 진실의 비트가 정하고(`sense` ×
면의 `Orientation`), 두 면의 바깥이 같은 쪽인지는 같은 곡면 → 향 비트, 같은 세계 이름 → 두 향, 그 밖엔 두 증인의
정의(`normals_agree_judge`)가 답한다(`face_facing`) — 클래스 병합과 같은 갈래 순서다. 평면 클래스 형태
(`PlaneWitness`)가 더하는 것도 정확한 서술뿐이다 — 이름의 행과 모션 전 프레임.

**판정층은 증인 삼각형 하나로 전체(total)다.** 판정 표의 계약은 «평면 위 세 정확한 점»(`Witness::tri_pt3`,
늘 있다)이고, 이름 없는 `Through` 평면도 그 점을 정의상 갖는다 — 정점이 모두 자기 프레임에서 좁은 만남으로 풀리면
진술의 점을 각자의 사슬과 평면의 모션으로 옮긴 것, 아니면(프레임 둘에 걸치거나 만남이 `Rat` 보다 넓은 정점이 있으면)
판정 프레임의 probe(`plane_witness`). 진실의
`Through` 는 판정층에서 변종이 아니라 **유도된 증인 삼각형**으로 나타난다.

**정확 지름길은 평면을 이름으로만 서술한다.** 평면 캐시(원점과 단위 법선 — 그 음함수의 `d` 는 f64 곱 `n·origin`)와 면
꼭짓점의 좌표 캐시는 둘 다 반올림된 상이라 진실 평면에서 `2⁻⁵⁴` 떨어질 수 있고, 둘을 서로 대조해 맞는 것만
싣는 검사는 «반올림된 평면»을 인증할 뿐이다(`3·0.1 = 0.3` 인 모서리를 벽 밖으로 판정해 유효한 교집합을
거절했다). 그래서 판정자가 이름 있는 평면을 묻는 서술은 이름(정의 점에서 반올림 없이 유도한 정준 정수,
`PlaneWitness::name_ints`) 하나이고 — 이름 정수 도로가 그 정수의 최근접 f64 로 필터를 돌리고 못 가른 것을 정수로 답한다 —
`orient3d` 는 넷째 평면도 그 이름으로 묻는다.
한 질문의 네 평면이 한 어휘로 서술되므로 두 서술이 어긋날 자리가 없다. 모션이 있는 평면도 같다 — 공유 모션
지름길은 이름의 base 프레임 정수를 읽고(`name_stored_ints` 가 손잡이까지 맞춘 것), 반올림된 base 삼각형에서 평면을
유도하지 않는다.

```rust
// crates/nacre-judge/src/kernel/frame3/witness.rs

/// 증인 점 — 유리수 base 를 모션 사슬로 나른다. base + chain 이 정의, `coord_iv` 는 f64 실현
/// (축마다 값+경계 한 원자 — 모델 캐시가 아니라 연산 하나의 것), hp 는 고정밀 메모.
/// 같은 정의는 같은 실현(경로 무관)이다. 동등성은 정의(base·chain)만 본다.
pub struct WitnessPoint {
    pub base: [Rat; 3],
    pub chain: HpRc<[MoveNode]>,              // 사슬을 펼친 것. 평면당 한 번 펼치고 증인 점 셋이
                                              // 같은 Rc 를 복제해 든다 — 두 번 펼치지 않으므로
                                              // 어긋날 수 없다. 점이 드는 이유: 자기완결 판정
                                              // (`shared_base`·실현은 점만 받는다).
    pub coord_iv: [Bounded; 3],
    hp: HpCell,                               // HpRc<HpOnce<(usize, [HpBounded; 3])>> — (정밀도, 값)
}                                             // `parallel` 에선 Arc<OnceLock>, 아니면 Rc<OnceCell>

/// 모션의 판정층 펼침 — 핸들 숲은 판정층이 못 푸므로 경계에서 핸들을 해소한다
/// (`Motion::Frame{plane}` → 유리수 벡터를 든 `MoveNode::Frame`). 모션은 교환되지 않으므로
/// 사슬의 **순서가 정의**다.
pub enum MoveNode {
    Rotate { axis: Axis, angle: Angle, pivot: [Rat; 3] },
    Translate { offset: [Rat; 3] },
    Mirror { axis: Axis, offset: Rat },       // 비고유(det = −1) — 공유 모션을 소거할 땐 손방향을 먼저 맞춘다
    Frame { frame: PlaneFrame },              // 평면 자신의 프레임으로의 기저 변환(유리수 벡터 넷·길이 셋)
    FrameWide(WideFrame),                     // 같은 것, 정준 데이터가 i128 을 넘는 평면 — BigInt 그릇
    FrameThrough(Box<FrameThrough>),          // 이름 없는 `Through` 평면의 프레임 — 증인 점을 노드 안에 든다
}

/// 이 연산이 무엇을 증명으로 인정하나 — 전부 모델에서 유도, 설정 없음(nacre-ops 의 `standard_for`).
/// 문턱은 비트 수가 아니라 **모델 단위의 길이**다: 같은 «256비트»가 한 번 돌린 솔리드엔 1e-76,
/// 삼백 번 돌린 솔리드엔 1e+15 의 해상도라, 비트는 커널이 모델마다 계산하는 구현 세부다.
pub struct ProofStandard {
    pub prec: usize,                          // 상승이 정의를 실현하는 정밀도 — 연산 내내 균일(메모가 따뜻하게)
    pub coincidence: Mag,                     // 이보다 가깝다고 **증명**되면 하나다(2⁻¹⁸⁰ 계열)
    pub scale: Mag,                           // 모델 크기 — 길이 한계를 방향 질문의 각도로 바꾼다
    pub cap: usize,                           // 비용 한계이지 해상도 한계가 아니다
}
// `coincidence` 는 전역 tol 이 아니다: tol 은 «이보다 가까우면 모르고 합친다», 이것은 «일치는 이보다
// 가깝다고 증명돼야 한다». 증명 못 하면 상승하고, 그러고도 못 하면 그렇게 말한다.

pub enum Decision {
    Sign(Orient),                             // 증명된 부호 — Zero 는 그것을 증명할 수 있는 경로에서만
    Coincident { within: Mag },               // 증명된 가까움 — 근거를 실어 보고
    Exhausted { at: usize, within: Option<Mag> }, // 한도에서도 0 을 걸쳤다 — 정밀도가 더 있으면 갈린다
    Degenerate,                               // 계량적 답이 없다 — 정밀도는 돕지 못한다
}                                             // 마지막 둘 → 이름 붙은 거절. `orient()` 는 미결을 Zero 로 접는다
```

- **어느 길로 갈지는 평면의 이름과 모션이 고른다** — 「숫자 규칙 — 일곱 개」의 5 와 CIP 의 «길은 tol 이 아니라 이름과
  모션이 고른다» 글머리(C4).
- 모션이 계수에 하는 일: `Translate` 는 `d' = d − n·t`(유리수면 정확), 축과 나란한 법선은
  `R·n = n`(구조적 검사 가능), 그 밖은 실현에 오차. 진실 쪽에는 오차가 없다 — 오차는 실현할 때
  생긴다.
- `Through` 평면의 판정은 두 갈래다.
  **유리수 닫힘이면 `Wide` 이름이 곧 정확 계수다 — 조건부로.** 세 술어(`orient3d`·`cmp_coord`·
  `dir_sign`)에 이름 정수 도로가 있다: 게이트 = 전 평면 이름 있음 ∧ (전부 무이동 **또는** 전부 한 사슬 —
  `cmp` 는 무이동만). 이름은 구성 시 저장 방향으로 σ-접혀(`name_stored_ints`) 필터를 거쳐 BigInt 부호
  술어(`int_*`, exact)로 정확히 답한다. **조건의 이유**: wide datum 의 이름은 세계 이름(담체 = 발견
  정점)이라, 그 평면을 처음 만드는 혼합-프레임 불리언에서는 일부만 잡히고(`wide_datum_cost.rs` 의 표),
  전부-세계 표(2세대)에서는 상승이 0이다 — wide 든 narrow 든
  (`a_second_generation_boolean_on_world_names_never_climbs`). 그릇은 `Expansion` 이 아니라 BigInt 다.
  **이름 없는(혼합 프레임) datum 만** 실현 상승 전용이다: 정점을 자기 세 평면의 고정밀 실현에서
  동차좌표로 만들고 그 위에서 계수를 구간으로 유도한다(**차수 9**) — 정확성 위험이 아니라 비용
  위험이다:
  - **깊이 1 은 필터가 산다.** 생성 200 사례 전부 결정, 계수 800개 중 미결 0, 최악 상대 반경
    6.3e-10. 차수 9 가 `Bounded` 의 여유를 먹지 않는다.
  - **깊이 2 는 필터가 없다.** 담체가 또 `Through` 면 차수가 **81** 이 되고, 계수가 `f64` 범위를
    8/8 전부 벗어난다(고정밀 쪽은 멀쩡하다). 함수는 그때 `None` 을 돌려 **상승으로 보낸다** — `NaN`
    반경이 우연히 `sign()=None` 을 내는 것에 기대지 않는다. ⇒ 깊이 2 는 «느린 길» 이 아니라 상승
    전용이다. 유계는 C6 이, 비순환은 C5 가 준다; 남은 것은 비용뿐이다.
  - **배율의 부호는 값 안에서 없앤다.** join 은 행에 대해 다중선형이라 결과가 참 평면의
    `D0·D1·D2` 배이고, 음수면 **평면 방향이 뒤집힌다**(`frame_sign`·바깥 법선·라벨 프레임이
    전부 그 위에 있다). 계수 `[Bounded;4]`·`[HpBounded;4]` 에 부호를 실을 자리가 없으므로 — 실을 곳
    없는 값은 아무도 안 쓰는 값이다 — 함수가 스스로 정규화한다(`canonical_plane_coeffs` 와 같은 관용구).
  - `D` 가 0 을 품으면 `None` → 상승 → 안 갈라지면 이름 붙은 거절. 조용한 폴백은 없다.
  - 이 동차 경로의 고정밀 판이 `plane_hp_through` 이고 판정 프레임의 동차 경로가 소비한다. f64
    필터 판은 없다 — 판정 표는 구간-계수 경로를 필요로 하지 않는다(판정 평면의 증인 = 프레임 probe).
    잠금은 Pure-대-Meet 차등이다.

---

## 표현력의 경계 — 어휘 밖의 사용례와 확장 경로

대표 사용례를 타입에 대입해 훑으면 **조용히 틀리는 시나리오는 없다** — 아래는 전부
«현 어휘로 진술 불가 ⇒ 구성 시점 거절»이고, 각각 확장 경로가 있다.

| 사용례 가족 | 왜 막히나 | 확장 경로 (계획) |
|---|---|---|
| **각도-계열 평면**: 각진 datum(모서리 축 θ)·드래프트 돌출·일반 챔퍼 | 진실이 *"이름 붙은 직선 둘레로 유리수 각 θ 돌린 평면"* 인데 그 어휘가 없다 — 기운 벽의 셋째 점은 `dist·tanθ` 라 무리수(`Known` 불가), `Motion::Rotate` 는 축 정렬 전용 | `Motion::Frame` 과 같은 «이름으로 들기»: 모서리(정점 핸들 둘)를 축으로 지목하는 회전 변종. 실현은 기존 기계(각도 캐시 + `1/√`) 그대로. **챔퍼가 이 어휘를 요구한다** |
| 임의 축 회전·임의 평면 미러(솔리드 이동) | `Rotate`/`Mirror` 축 정렬 전용 — 비목표. 새 피처의 임의 방향 배치는 `Frame` 으로 된다 | 위와 같은 계열 — 그 설계가 들어올 때 함께 |
| 혼합 정의 datum (정점 2 + 좌표 1) | `Through` 는 핸들 3 전용, `Known` 은 값 3 전용 | 계획: `Through` 원소를 「좌표 값 \| 정점 핸들」 참조로 일반화 — 평면의 점 데이터 수준이라, 정점 아레나에 좌표 변종을 두는 것과는 다른 자리다 |
| 중간 평면(midplane) | 같은 프레임 평행 면 사이는 유리수 평균이라 `Known` 으로 계산 가능. **프레임이 다른** 두 면 사이는 세계 좌표가 무리수 + 어느 정점도 안 지난다 | 계획: 필요해지면 곡면 핸들 둘을 드는 정의 변종 |
| 스케치 내부의 정확한 각도(정육각형·30° 변) | 꼭짓점에 √3·tan30 — `[Rat; 2]` 로 못 적는다. 앱이 계산한 f64 의 십진수가 진실이 된다(의도된 동작) | **정확한 각도는 좌표가 아니라 모션으로 표현된다**(회전·프레임) — 결함이 아니라 이 숫자 시스템의 정의적 성질 |

막히지 않는 것: datum 이 참조한 솔리드의 이동(이동본 정점으로 재지시), superseded
참조(append-only), 원형 패턴(`360/n` 은 항상 유리수), 깊은 체이닝(비용만 — 회전 3200회),
datum 낀 불리언(비용만), 4+평면 정점(validate 몫), 오프셋의 오프셋(연산당 1단),
STEP 출력, undo/replay.

---

## 연산 층 (`nacre-ops`) — 재계산 가능한 로그

`ops/mod.rs`(`Operation`·`DatumDef`·`OpOutput`·`SketchFrame`), `ops/apply.rs`(`apply`·`replay`), `ops/datum.rs`(`datum_plane`),
`ops/frame.rs`(`frame_toward`·`face_sketch_frame`·`push_frame_node`), `ops/prism.rs`(`extrude_on_frame`·`exact_frame`·`build_prism`),
`ops/profile.rs`, `construct.rs`, `sketch.rs`, `error.rs`, `reject_census.rs`. 변종마다의 약속과 거절은 그 doc 이다.

```rust
pub enum Operation {
    DatumPlane { def: DatumDef },  // Stated(SketchPlane) | Offset { frame, dist } | ThroughVertices([Handle<Vertex>; 3])
    Extrude { frame: SketchFrame, profile: Profile2d, dist: f64 },   // dist 는 부호 있는 두께, 0 은 거절; 스케치는 이 인자다(별도 Sketch op 없음)
    Boolean { kind: BoolKind, a: Handle<Solid>, b: Handle<Solid> },   // BoolKind = Fuse | Cut | Common
    Transform { solid: Handle<Solid>, isometry: Isometry },
    Mirror { solid: Handle<Solid>, axis: Axis, offset: Rat },
    Copy { solid: Handle<Solid> },
}
pub fn replay(ops: &[Operation]) -> Result<Model, OpError>;
```

`replay` 는 `Model` 의 메서드가 아니라 nacre-ops 의 자유 함수다 — `impl Model` 은 inherent impl 이라 `Model` 이 정의된
nacre-topo 에만 놓일 수 있는데(coherence), replay 는 `Operation` 을 다루므로 topo 위 층에 산다. 핸들의 재고정은 `replay`
doc 과 「저장소와 동일성」의 «인덱스 어휘». `replay` 는 진실 `Model` 만 돌려주고, 테셀레이션은 `nacre-tess` 가 모델에서
따로 유도하는 캐시다(`TessConfig`).

**연산은 평면을 «이름 부른다».** 평면을 값으로 싣는 변종은 `DatumPlane` 의 `Stated` 하나뿐이라 원칙 2 가 연산 어휘
전체에서 성립하고, 밑캡은 같은 평면의 둘째 진술이 아니라 공유된 핸들이다(`Operation::Extrude` doc). 프레임의 문은 셋이다
— 세계 평면 `SketchFrame::world`, 기존 면 `face_sketch_frame`, 그 밖의 평면은 `Operation::DatumPlane` 으로 먼저 진술하고
돌려받은 프레임.

**ŵ 는 프레임의 것이고, «어느 쪽»은 `dist` 의 부호다.** `flip` 은 프레임을 만든 쪽이 진실에서 한 곳(`flip_toward`)에서
정하므로 면의 프레임은 바깥을 향하고, «면에 그려 안쪽으로»는 같은 프레임의 음의 돌출이다 — `Operation::Extrude`·
`SketchFrame`·`flip_toward` 의 doc(실현한 ŵ 의 내적으로 정하지 않는 까닭은 「가지 말 것」 «향 관계를 캐시 법선으로
읽기»). 원통 축은 쓸기 부호와 무관하게 프레임의 법선이라 «datum −n 에 +d» 와 «datum +n 에 −d» 는 같은 원통을 두 진술로
적는다 — 원통의 열쇠가 진술 그대로인 것(todo 「곡면 둘 이상이 만나는 점과 seam 담체」)과 같은 부류다. ZX 씨앗의 `+u`
예외는 `SketchFrame::world` doc, 프레임 기저를 실현해서 되묻지 않는 까닭은 `RatFrame::of_plane_frame` doc.

**datum 평면.** 평면을 만드는 연산은 `DatumPlane` 하나이고 불리언은 평면을 만들지 않는다(깊이 불변식,
`a_boolean_mints_no_surface`). 진술의 향 `−(진술된 법선)` 과 함께 뒤집어야 하는 자리들은 `datum_plane` doc, 호출자의
프레임을 함께 돌려주는 까닭은 `OpOutput::DatumPlane` doc, `Offset` 이 한 평면에 두 핸들을 막는 정규화는
`DatumDef::Offset`·`OpError::ZeroOffset` doc 과 그 팔의 주석이다.

- **`ThroughVertices` 는 값 어휘로 말할 수 없는 하나다** — 핸들은 값 타입 `SketchPlane` 에 들어갈 수 없어 어휘가
  `DatumDef` 의 변종으로 자란다. 좌표로 지으면 다른 평면이 되는 측정은 `tests/instruments/point_width.rs`, 정점 순서가
  방향인 것과 원인별 거절은 `DatumDef::ThroughVertices` doc 과 `datum_plane` 의 그 팔, 점과 모션을 한 결정에서 내는 까닭은
  `through_points_by_cause` doc, 이동이 그런 평면에 노드를 남기는 까닭은 `transform.rs` 의 `carry_of` doc.
- **datum 은 솔리드 이동을 따라가지 않는다.** append-only 라 datum 이 가리키는 평면은 남고, 나중 `Transform` 은 새 평면을
  만든다 — datum 은 독립적인 기준이다.
- **쓰이지 않는 datum 은 어디에서도 새지 않는다**: `validate` 는 live 도달 집합만 세고 surface 는 경계 검사만 받으며(씨앗
  셋이 이미 영구 orphan), `nacre-step`·`nacre-tess` 는 surface store 를 돌지 않는다(둘 다 면 경유).

**파라메트릭 편집의 진화 경로 (계획: v2 이후, 기록만).** 연산이 원시 Handle 대신 계보 참조
`OpRef { op: usize, output_slot: usize }`(«연산 N 이 만든 k번째 면»)를 담으면 상류 수정이 하류 Handle 번호를 밀어내도
(topological naming) 참조가 의미로 풀려 편집-재생이 된다. v1 은 짓지 않되 로그 직렬화 포맷에 이 확장이 포맷 파괴 없이 들어갈 자리를 남긴다. 인덱스 어휘는 이를
밀어내지 않는다 — replay 안의 재고정 지점이 곧 `OpRef` 해석이 들어올 자리다(교체이지 경쟁이 아니다). 테셀레이션
tolerance 도 같은 맥락에서 연산별 override(`Operation` 항목의 선택 필드)로 자랄 수 있다.

### 프로파일 — 다중 루프와 계약

`ops/profile.rs`(`Profile2d`·`Ring2d`·`Profile2d::check`), `ops/prism.rs`(`build_prism`·`oriented_ring`·`frame_rings`),
`nacre-geom::intersect::flat_corner`; 시험은 `tests/coverage/sketch.rs`.

`Profile2d` 는 외곽 링 하나와 구멍 링 여럿을 담고, `build_prism` 이 구멍마다 벽을 세우고 두 캡에 내부 루프를 단다 — 구멍
있는 스케치는 불리언 없이 한 번에 선다. 그것이 커널 몫인 까닭은 overview 「범위」의 «설탕 vs 커널 판별 기준»이다.

- **감김은 sweep 을 아는 `build_prism` 한 곳에서 정한다** — `Profile2d`·`build_prism`·`frame_rings` 의 doc. 오일러의 `L_i`
  항이 실제로 필요한 형상은 `a_donut_prism_is_a_valid_genus_one_solid`, 같은 도넛을 잘라 만든 것과의 대조는
  `a_swept_hole_matches_the_same_shape_cut_out`.
- **생성자는 진실을 잡고, 계약은 소비하는 연산이 검사한다.** 십진 리프트·평탄 코너 소멸(중복점·스파이크는 살아남아 제
  이름으로 거절된다)과 계약 넷(단순성·서로소·포함·구멍 속 구멍 없음), 판정이 저장된 유리수 위에서인 까닭은
  `Profile2d`·`Profile2d::check`·`flat_corner` 의 doc. 검사가 없으면 조용히 틀린 답이 나오는 측정은
  `silently_wrong_profiles_are_refused`, 이 검사가 `oriented_ring` 의 전제(넓이 0 없음)를 세우는 것은 그 doc.
- **단순성은 작성된 입력의 계약이지 커널 데이터의 불변식이 아니다** — 불리언이 만드는 윤곽의 안쪽은 배열이 정하고 그
  윤곽은 figure-8 조임이어도 정당하다(`loop_winding`). 그래서 `validate` 에 면 단순성 검사를 넣지 않는다 — 정당한 결과를
  거절하게 된다.

### 스케치 어휘 — 링 중첩·고리·호

`sketch.rs`(`from_rings`·`from_paths`·`Ring2d::new`·`arc_turns`·`arc_to_rat`·`SketchError`), `nacre-geom::{mixed, intersect}`
(`Edge2d`, 링 술어), `nacre_exact::{orient2d_rat, winding_sign_quarter_arcs}`, `ops/prism.rs`(호 벽,
`refuse_non_quarter_arcs`).

사용자는 닫힌 경로만 그리고 무엇이 구멍인지 말하지 않는다 — `from_rings`/`from_paths` 가 포함 깊이의 짝홀로 섬별
프로파일을 돌려주고(채우기 규칙 인자는 없다), 자기교차를 중첩보다 먼저 본다(`from_rings`·`classify`·`SketchError` 의 doc).
커널 `Extrude` 1회는 연결된 덩어리 1개이고 섬마다 부르는 것은 편의 레이어라, `Extrude` 에 다중 바디 출력은 필요 없다.

**커널의 문은 고리다** — 정점과 각 정점을 떠나는 변(`Ring2d::new`)이고, 커널은 «어떻게 그렸나»를 모른다(`sketch` 모듈
doc). ops 에는 정책(깊이 패리티)과 문만 있고, 부호 술어는 `nacre-geom::mixed` 에 산다 — 호까지 ℚ(√c) 의 `QuadVal` 부호로,
광선은 `intersect` 의 반열림 규칙을 옮겨 적고(그 모듈 doc), 다각형 링에서 `intersect` 의 술어와 답이 같음을
`polygon_rings_answer_exactly_as_before` 가 잠근다. 부호 원시는 `nacre_exact::orient2d_rat`(그 doc). 끝점은 정확히 일치해야
하고 스냅은 없다(원칙 4 — `mixed` 모듈 doc, `ArcEndOffCircle`·`ZeroLengthEdge`).

**호와 원.** 조각 타입은 `Edge2d { Line, Arc { center, r2, ccw } }` 하나이고 반지름은 제곱으로만 든다(`Edge2d` doc). 회전각은
커널의 말이 아니다 — 사분각 문(`arc_turns`·`arc_turns_rat`)과 임의 끝점 문(`arc_to_rat`)이 끝점을 `Rat` 으로 준다(`sketch`
모듈 doc). 같은 원·같은 방향 이웃 호의 병합은 `Ring2d` 의 정규형이다. 프리즘은 호마다 원통 조각 벽을 세우고(`wall_surfaces`
— `ref_dir` 은 `Seg3::Arc` doc), 옆면의 향은 `push_walls` 의 주석, 감김은 `winding_sign_quarter_arcs` 라 사분각 배수 아닌
호는 `ArcSweepNotQuarterTurn` 이다. 직선–호 꼭짓점 `Vertex::Pierce` 와 호–호 꼭짓점 `Vertex::CylinderPair` 의 이름은
「원통의 진실」 6·7. 두 호 벽 사이의 옆 간선은 공유 룰링이고, 그런 솔리드를 받은 불리언은 입구에서 `CylinderPairEdge` 로
거절한다 — 차트는 옆면의 룰링이 평면이 자른 것이라고 읽는다.

### 거절 분류

`error.rs`(`BoolError`·`RejectReason`·`RejectClass`·`RejectWhere`·`DeclineKind`·`reject`/`reject_at`), `reject_census.rs`
(`ReasonId`); 시험은 `tests/reject_census.rs`, 파사드의 `tests/facade.rs`.

**거절의 이유는 값에 실려 나가고, 소비자는 `RejectReason::class()` 로 분기한다.** 세 클래스의 뜻, 사유는 자라고 클래스는
안 자라는 까닭(`#[non_exhaustive]` 의 자리), `TraceDeclined` 가 첫 면을 이름 대는 것은 `BoolError`·`RejectClass`·
`RejectReason`·`DeclineKind` 의 doc. 위치(`at: Option<RejectWhere>`)는 사유 옆에 싣는 증인이다 — 범주와 측정을 가르는
까닭, 결정적 최솟값, `Eq` 없음은 `RejectWhere` doc 이고, 위치를 싣는 사유는 `reject_at` 을 부르는 자리다.

**오류의 `Display` 는 로그용 한 줄이다**(std 관례: 소문자, 마침표 없음). `BoolError` 는 문장을 새로 짓지 않고 사유의
식별자를 클래스·위치와 싣고(`boolean rejected, not supported: cylinder_pair_contact`), 작성자에게 보이는 문장은 kit 이
든다. 감싸기만 하는 변종(`OpError::Boolean`·`CacheDecline::Cannot`)은 안쪽의 줄과 `source` 를 그대로 넘기고, 다른 원인을
싸는 변종(`StepError::Io`)은 제 맥락만 말하고 원인을 `source()` 로 준다 — 안쪽을 `Display` 에 싣고 `source` 로도 주면
보고가 같은 말을 두 번 찍는다. 파사드로 닿는 오류 타입이 모두 `std::error::Error + Send + Sync` 인 것은
`every_error_type_is_a_std_error` 가 목록으로 잠근다.

**어떤 사유가 실제로 발화하는지는 상설 census 가 답한다 — `nacre-ops::reject_census`.** `reject`/`reject_at` 이 모든
`Rejected` 를 짓는 유일한 깔때기라 `#[track_caller]` 하나가 모든 호출 지점을 계측하고, 테스트 빌드의 계측이라 제품
빌드에는 없다. 「울린 것」과 「밖으로 나간 것」은 다른 인구이고 그 갭이 남은 삼킴의 지표다 — 그래서 `point_in_component`
의 «이 노드로는 못 정함»은 삼켜질 오류가 아니라 기권(`Said`)이다. 측정과 읽는 길 둘(얼린 코퍼스의 게이트, `reject-trace`
전수 훑기)은 `reject_census` 모듈 doc, 게이트가 집합만 박고 줄 번호·횟수를 걸지 않는 까닭과 못 보는 것은
`tests/reject_census.rs` 의 doc, 키가 `detail` 을 드는 까닭은 `ReasonId` doc.

### 구성 산술은 **사용자가 쓴 십진수**로 한다

`construct.rs`(모듈 doc, `RatFrame`), `ops/prism.rs`(`exact_frame`·`frame_rings`), `nacre_exact::Rat::{from_decimal, to_f64}`;
잠금은 `tests/invariants/pipeline.rs` 의 `a_split_dimension_meets_the_undivided_one`.

f64 의 프리즘 좌표(`origin + x·u + y·v`, `+ normal·dist`)는 결합법칙을 지키지 않아, `7.7` 로 한 번에 올린 캡과 `1.1` 다음
`6.6` 으로 올린 캡이 1 ULP 떨어진다 — 쪼개지지도 틀리지도 않고, 아무도 신고하지 않는 얇은 면을 끌고 다닌다(측정은 그
시험). f64 를 `try_from_f64` 로 올리면 정확 산술이 이진 드리프트를 그대로 재현하므로, 복구하는 것은 사용자가 쓴 십진수다
(`construct` 모듈 doc, `Rat::from_decimal` doc).

**결정.** 공개 API 는 f64 그대로라 호출부·TS/JSON 경계가 바뀌지 않고, 커널이 받는 자리마다 `Rat::from_decimal` 로 왕복하는
최단 십진수를 잡는다 — 프로파일 좌표(`ProfileOutsideDecimalWindow`), extrude·offset 거리(`DistOutsideDecimalWindow`),
`SketchPlane` 의 생성자(`from_axes` 등), `SketchFrame::named`(`FrameOutsideDecimalWindow`). 커널이 실현한 기저는 리프트하지
않는다(`realized_plane` doc). 전제는 `Rat::to_f64` 의 최근접 반올림이다(그 doc). 프레임의 세계 기저는 진실에서 유리수로
묻고(`exact_frame`), 없으면 스케치를 프레임 노드 안에 프레임 좌표로 적는다 — 기울어진 면에서도 쪼갠 치수는 한 평면이다
(`two_caps_reached_by_different_arithmetic_are_one_plane` — 두 캡이 한 평면으로 판정된다). 정확한 꼴이 없는 평면(축이
십진 창 밖이거나 퇴화, 배치 산술의 `i128` 넘침 등 — 그 변종의 doc)은 조용히 f64 로 짓지 않고 `PlaneWithoutExactForm` 으로 거절한다.

**못 고치는 것 (비목표).**

- **스크립트가 계산해 넘긴 값.** 커널이 고치는 것은 자기 누적이고, 스크립트가 f64 로 계산한 값은 이미 다른 수다
  (`Rat::from_decimal` doc).
- **4-평면 동시성.** 유리수 입력에서도 생기는 정확한 0 이라 구성 산술이 고칠 것이 아니다 — `RejectReason::FourPlane` doc;
  판정이 평면의 이름으로 묻는 것은 `every_answer_on_unmoved_named_planes_is_the_names_answer`.

### 불리언 엔진 — 평면당 셀 복합체

불리언 경로는 **하나**다: `nacre-ops::arrangement` 의 배열 엔진(`crate::boolean` 의 공개 `boolean` 이 위임한다) — 평면 클래스마다
독립인 셀 복합체와, 원통 클래스마다 제 차트(「원통 — 밴드·게이트·차트·룰링」). 글로벌 detector 도, regime 별 이중 경로도 없다.
평면 위에 놓인 기하는 퇴화가 아니라 **입력**이다(`arrangement` 모듈 doc). 결과 전체의 판정(자기접촉·닫힌 껍질)은 조립(`assembly`)이 한다.

1. **평면 클래스.** 두 솔리드의 면 평면을 클래스로 묶고, 이후의 일은 **평면 클래스당** 독립이다.
2. **추적(trace).** 면을 모든 클래스 위에 **선분**으로 — Seated·Transversal 두 경우와 면 단위 기록 `declined`(밖으로는
   `TraceDeclined { kind, face }`)는 `arrangement` 모듈 doc. **선 위의 추적은 스윕 하나다**(`sweep_line`): 평면 면의 클래스 선과
   옆면의 룰링이 같은 스윕을 쓴다. 정거장의 사건 순서, 꼬집힌 면의 한 정거장, 두 이름의 동점을 라운드의 표에 기록하는 일
   (`Aliases::record_same_point`)과 그 거절(`CoincidentFeatures`·`CylinderPair`)은 `sweep_line`·`LineEvent`·`record_same_point` 의 doc.
3. **교차 분할.** `split_at_crossings`; 한 클래스의 골격은 `class_skeleton`.
4. **셀 추출·중첩.** `+1` 감김 셀만 host 다(`nest_cells`, `emit_faces` 의 doc).
5. **라벨.** `Label = [A_above, A_below, B_above, B_below]` 를 모서리 마스크(`edge_mask`)로 전파하고 **모든** 모서리에서 검증한다
   (`label_cells`; `UnreachedCell`·`LabelConflict`). «above» 는 클래스의 **저장된** 법선 쪽이다 — `Label` doc(`frame_sign`, 방출의 `flip`).
6. **생존 규칙·방출.** 아래.
7. **청소·조립.** `arrangement::whole_result_faces`(`unify_coplanar_faces` → 옆면 방출 → 결과 전체의 `dissolve_straight_angles`)
   → `assemble_fuse_cut` 이 vertex handle 을 면의 first-appearance 로 배정한다(「병렬 — 판정만 병렬, 순서는 고정」).

**생존 규칙.** 한 챔버(평면 한쪽의 공간)의 `(inA, inB)` 에 대한 술어 하나다(`draft::keep`):

| `kind` | `keep(inA, inB)` |
|---|---|
| `Fuse` | `inA ∨ inB` |
| `Cut` | `inA ∧ ¬inB` |
| `Common` | `inA ∧ inB` |

`+1` 셀은 두 챔버(W 의 위·아래)가 `keep` 아래에서 갈릴 때 결과의 면이고, `flip = keep_above == (frame_sign > 0)` 하나가 어느 쪽이
재료인지를 나른다(`emit_faces` doc). 접촉면의 Requicha 식 표 — 연산 × 상대 법선 — 는 이 한 규칙의 귀결이다: Fuse/opp=`inP⊕inQ`,
Fuse/same=`inP∨inQ`, Cut/opp=`inP`, Cut/same=`inP∧¬inQ`, Common/opp=∅, Common/same=`inP∧inQ`(`inP`·`inQ` 는 각 피연산자의 면이 그 셀을
덮는가). 볼록성 항이 없어 비볼록·다중루프를 half-space 가정 없이 `point_in_ring` parity 로 정확히 처리한다.

- **동일평면 병합은 필수 청소다** — 내부를 지우고 바깥을 꿰매지 않으며, 꼬집힌 점에서는 배열의 셀 걷기와 같은 순환 순서
  (`dir_at`·`angular_order`)로 나가고, 못 하면 그룹을 그대로 내보내 기권한다. 기권이 안전한 까닭(«정돈일 뿐»이 아니라 «패스가
  없던 자리»)과 그 뒤 전역 판정이 모양을 이름 대는 것은 `unify_coplanar_faces`·`merge_component` 의 doc. 평각 정점은 옆면이 붙은
  뒤 결과 전체에서 한 번 녹인다 — `dissolve_straight_angles`·`whole_result_faces` 의 doc.
- **45° 폴드의 거절은 갭이 아니라 «정답»이다** — 한 팔이 제 몸에 공면으로 닿는 두께 0 솔리드라 `SelfTouchingResult`(`Impossible`)
  다. 잠금은 `tests/probes/rotation_sweep.rs` 의 `another_rejection_keeps_its_own_reason` 과 `tests/reject_census.rs` 의 `fold_45`.
- **칼날 접촉**(모서리가 상대 면 평면 위, 두 면이 한쪽으로만 떠남) — 스침 결합은 홀짝이고(`edge_mask` 본문), 소식 없는 모서리는
  골격에서 빠지되 룰링 위 칼날의 두 끝에서 룰링 조각을 먼저 자른다(`carries_news`·`class_skeleton` 의 doc; 「점의 이름」의 «선은
  그 위의 모든 꼭짓점에서 잘린다»). 선·점으로만 닿는 두 몸통은 둘로 나오고, 핸들은 성분이 아니라 **출력 솔리드**(재료 성분 + 그
  공동들) 단위로 배정된다 — `assembly::Grouping`·`NonManifoldResultEdge` 의 doc. 잠금: `tests/probes/contact_separates.rs`,
  `tests/probes/knife_edge.rs`.
- **4평면 동시성의 세 갈래** — 피연산자 꼭짓점의 동시성(`NamedRing::concurrencies` 로 추적 전에 표에 넣고, 런 꼭짓점을 지나는 넷째
  평면은 `third_on_l` 이 본다), 선을 **가로지르는** 두 평면(`split_at_crossings` 의 동점), 선을 **담는** 평면(`Aliases::canon_wall` 의 기록을 그
  선의 분할점마다 되읽는다 — `split_at_crossings` 본문). 놓치면 한 점이 두 이름으로 seam 표에 와 `SeamAlias`(`SuspectedDefect`) —
  유효한 입력을 두고 커널이 자기를 탓한다. 잠금: `tests/probes/concurrent_line.rs`.

**결과 계약: 표면이 자기와 닿는 솔리드는 만들지 않는다.** 명제(모서리가 꼭짓점 아닌 점에서, 또는 꼭짓점이, 같은 솔리드의 자기가
경계하지 않는 면의 닫힌 영역에 닿는가)와 단위(출력 솔리드)는 `assembly::self_touch::self_touch_reject` 의 doc 이고, 이런 몸체가 모든 위상
계수와 `validate` 를 통과하는 까닭은 `SelfTouchingResult` 위의 주석이다. 방향은 Parasolid 를 따르고 OCCT 는 받는다; 대가는 쐐기 포켓의
끝이 마침 반대 벽에 닿는 모델링이 막히는 것 — `tests/probes/self_touch.rs` 모듈 doc. 거절은 «무엇이 잘못됐는지»만 말하고 사람이
읽을 문장은 앱이 `RejectClass` 로 갈라 쓴다(`BoolReport`·`RejectClass`·`SelfTouchingResult` 의 doc).

- **체와 정확 판정** — 면 상자(도달 구간; 호는 그 원 전체 — 체라 호의 방향은 읽지 않는다 — 원을 진술 못 하면 무한)를 격자(`FaceGrid`)에 넣고, 두 끝이 면 평면
  위에 정확히 놓인 것(`side_of`)만 정확 판정한다. 측정(666,653 쌍 → 1,583)과 상자를 부풀리는 까닭은 `self_touch_reject`·`FaceGrid`
  의 doc; 후보 평면을 끝점 이름에서 고르지 않는다(「가지 말 것」).
- **평면 팔** — 끝점이 아니라 **열린 구간**이 면의 «안» 칸과 겹치는가(`combinatorics::segment_meets_face`): 관통 슬롯, 모든 경계를
  한 자루로 읽는 짝-홀, 구멍 원을 읽는 때와 침묵 클래스의 원, 닫힌 영역의 접점(`arc_touch`), 원통이 깎은 링(`ClassLine`·
  `pin_for`), 간선 끝의 고정과 `RingNaming` 은 그 doc. 선 `q ∩ w` 의 `w` 를 간선 자신의 두 면에서 얻는 까닭은 `self_touch_reject` 본문.
- **호 간선은 묻지 않는다** — 평면 면 안의 호는 병합의 몫이고, 꼭짓점 없는 점에서 평면 면에 닿는 호는 접점·접선 판정의 몫이거나
  게이트가 먼저 거절하는 모양이다: 비스듬한 면이 옆면 경계의 원에 한 점으로 접하면 `ObliqueCylinderCut`(게이트가 그 원 모두를
  묻는다 — `oblique_section_meet`), 다른 원통의 호면 `CylinderPairContact`. 이 둘은 이 판정이 게이트에 기대는 자리다
  (`self_touch_reject` 본문의 «The arms» 주석).
- **옆면 팔** — 같은 명제를 옆면에 묻는다: 그 옆면의 노드가 아닌 꼭짓점이 원통 곡면 위이고 옆면의 닫힌 영역에 있는가
  (`combinatorics::named_point_in_face` 의 옆면 갈래), 직선 간선이 끝이 아닌 점에서 그 영역에 닿는가
  (`combinatorics::segment_meets_lateral` — 게이트가 기록한 룰링을 따르는 간선은 `ruling_meets_lateral`, 선은 `SharedRuling::at`;
  그 밖의 선은 원통과의 근). 닿는 점이 접선 판정의 행이 «벽 면이 옆면 위에서 접선과 만난다»고 적은 닫힌 구간
  (`planes::Tangency::contact`) 위이고 그 벽 클래스의 결과 면이 덮으면 선 접촉의 한 점이라 그 판정의 몫이다(`tangent_hand`) —
  넘김의 영역이 판정의 영역이라, 접선에 모서리 한 점으로 닿는 면·옆면 너머에서 접선에 닿는 면은 점 접촉이다. 옆면 경계(림·호) 위의
  점에 그 곡선을 옆면과 함께 쓰는 면의 간선이 닿으면 그 면의 두 경계가 만나는 집힘이다(`shares_a_curve` — 평면 팔이 제 면의
  경계를 묻지 않는 것과 같은 규칙). 체는 옆면 경계의 원들의 상자이고 경계는 체를 지난 첫 후보에서 읽는다; 못 읽은 경계는 그것이
  답을 정하는 후보만 들고, 넘김의 주인을 못 읽으면 든다; 못 읽음은 버린다(수는 `self_touch_reject` 본문).
- **점 팔** — 꼭짓점을, 그것이 경계하지 않는 같은 솔리드의 모든 평면 면에 이름의 종류와 무관하게 묻는다. 경계 위도 세는 대가
  (생산자의 T자 이음)와 버리는 기권은 `self_touch_reject` 본문, 읽기는 `combinatorics::named_point_in_face` 한 함수. 림 위에 놓인
  꼭짓점은 배열이 림의 두 면에 다 용접해 «경계하지 않는 면»이 없고, 집힘으로 결과 검사가 이름 짓는다(`nacre_topo::nonmanifold_vertices`).
- 증인은 이름 순서의 최솟값이고, 세울 수 없는 정확 판정의 오류는 자기접촉이 하나도 없을 때만 낸다 — 증명이 «못 셈»을 이긴다
  (`self_touch_reject` 본문).

### 병렬 — 판정만 병렬, 순서는 고정

`nacre-ops::par` 의 헬퍼(`map_range`·`map_every`·`try_map_range`, `PARALLEL_FLOOR`; `parallel` 피처, rayon)가 평면 클래스당 독립인 두 루프
(`trace_result_faces` 의 pass A 트레이스→병합→분할, pass B 셀→중첩→라벨→방출 — 잰 불리언 시간의 90%)와, 모델의 정밀도를
**점마다** 읽는 `standard_for` 를 나눈다.

- **판정만 병렬, 변이·순서는 단일 스레드 고정.** `Store::push` 순서 = handle 정체성 = 재생 결정론이라 헬퍼는 인덱스 순서로 모으고
  병렬 reduce 하지 않는다 — `par` 모듈 doc. 에러는 인덱스-최초(rayon 의 `collect::<Result<_,_>>()` 를 쓰지 않는다)이고, 거절
  입력에서 뒤 클래스에 잠복한 패닉이 닿게 되는 대가는 `try_map_range` doc. 문턱은 재서 넣은 비용 손잡이다 — `PARALLEL_FLOOR`
  doc 의 측정 표와 그 시험.
- **별칭 표는 라운드 스냅샷으로.** 클래스는 라운드 시작의 표와 제 발견을 보고, 순차와 같은 고정점에 이르는 까닭은 스케줄이 아니라
  `Aliases` 의 성질(최소 대표·단조 누적)이다 — `fixed_point` doc. 별칭 36개가 두 라운드에 서는 4평면 모델을
  한 스레드와 여럿이 같게 푸는 것은 `src/tests/suite/parallel.rs`.
- **정밀도는 max 리덕션이라 안전하다** — `trial_bound`·`standard_for_points` 의 doc.
- **브라우저는 순차 빌드다** — 그 까닭과 거기서 병렬에 드는 것은 `par` 모듈 doc.

### 상대가 닿을 수 없는 평면은 소식을 나르지 않는다

`nacre-ops::reuse` 한 모듈이다: 상대 솔리드와 떨어져 있음이 **증명된** 평면 클래스(`class_plans` → `ClassPlan::{PassThrough, Empty}`)는
셀 패스(pass B)를 돌지 않고 자기 솔리드의 면을 배열의 어휘로 다시 적어 낸다(`pass_through`, `VertexClasses`). `trace_result_faces` 의
출력은 그대로 `Vec<LocalFace>` 라 하류(조립·병합·seam 표·검증기)는 이 지름길을 모른다. 왜 필요한가(80핀 fold 의 마지막 fuse 는 169
클래스를 배열하는데 새 핀이 닿는 것은 27)와 증명이 설계의 전부인 까닭(`WorkingPlane::misses` → `nacre_judge::orient3d_filter`, `None`
은 «배열하라» — 증명 실패는 속도만 잃고 톨러런스가 들어오지 않는다)은 `reuse` 모듈 doc 과 `orient3d_filter` doc.

- 정점만 봐도 되는 까닭(평면 다각형은 제 정점의 껍질 안), 정점 이름을 솔리드에서 읽는 법(넷 이상은 동시성이라 배열한다), 규칙표가
  접히는 모양(양쪽이 면을 가진 클래스는 늘 배열)은 `class_plans`·`VertexClasses` 의 doc.
- **원통 클래스가 하나라도 있으면 reuse 를 통째로 끈다**(`trace_result_faces` 진입의 가드) — 껍질 논증의 평면 전제, 평면 클래스
  단위로만 면을 옮기는 `pass_through`, 밴드가 읽는 원판 라벨이 모두 그 가드에 선다(그 자리 주석). 원통 입력에 reuse 를 되살리는
  일은 성능 작업이 아니라 **밴드의 소속 판정**(「밴드와 라벨」)을 다시 설계하는 일이다.
- **추적 패스(pass A)는 모든 클래스에 돈다** — `PassThrough` 가 정점에 이름을 못 붙이면 배열로 물러나고 그 배열은 그 클래스의
  분할을 읽는다(`settle_rounds` 본문; 건너뛰기는 「가지 말 것」). pass B 를 건너뛰는 결정은 `try_map_range` 의 클로저 안에서
  한다 — 그 자리 주석.
- **차등 검사가 영구히 남는다** — 디버그 빌드는 같은 불리언을 `ClassReuse::{Off, Proved}` 로 두 번 풀어 같은 면을 요구한다(`Off` 는
  폴백이 아니라 기준이다; 정규형 `reuse::canonical`, 참조 실행의 제 `Notes`) — `arrangement::boolean` 본문과 `ClassReuse`·`canonical`
  의 doc. 음성 대조(잰 것): 분리 증명 없이 소유 클래스를 전부 통과시키면 커버리지 167건이 실패한다.
- **잰 값.** 회전 fold 80 순차 9.66 → 7.48 s(1.29×), 병렬 2.38 → 1.97 s(1.21×). 80핀 fold 에서 클래스 방문 8800회 중 7021회가
  지름길을 탄다.

## 원통 — 밴드·게이트·차트·룰링 (`nacre-ops::{bands, arrangement::cyl_chart}`)

원통의 벽은 평면 클래스가 아니라 밴드로 정해진다. 배열 엔진은 **평면 클래스 하나씩** 답하고 원통의 측면은 그 어휘
밖이지만, 별개의 새 엔진 없이 **균일-슬랩 정리**에 선다 — 연속한 두 ⊥ 절단 사이의 열린 원통 슬랩에 상대 솔리드의 경계가
없으면 상대 소속이 슬랩 전체에서 균일하고 증인 하나가 밴드를 정한다(정리와 그것이 서는 인구는 `bands` 모듈 doc).
게이트가 **기록하고 통과시킨** (벽, 원통) 쌍에서는 그 전제가 깨져 밴드가 잘린 원과 룰링으로 끊기고(「룰링」), 옆면은
**자기 차트 위의 셀 복합체**로 답한다(「차트와 셀 방출」).

흐름은 한 방향이다: 게이트가 (평면 클래스, 원통) 쌍마다 답을 기록한다(`planes::cylinder_gate` → `WorkingCyl::records`·
`shared`·`grazed`, `Tangency` 행) → 클래스마다 배열이 그 기록만 믿고 원·룰링·현·단면을 그린다(`arrangement::per_class`) →
옆면은 제 차트에서 셀을 읽는다(`arrangement::cyl_chart`) → 결과 전체의 면을 읽는다(`arrangement::whole_result_faces`,
`draft::HeldRims`) → 조립이 짓고 자기접촉을 판정한다(`assembly::reconstruct`, `assembly::self_touch`).

### 밴드와 라벨

밴드의 소속은 배열의 라벨을 읽고(`cyl_chart::read_cell::read_bits`), 옆면의 존재는 추적의 자국을 읽는다
(`read_cell::face_spans`, `arrangement::ArcLabel`). 밴드의 끝은 방출된 원 경계와 옆면 자신의 rim 이고
(`bands::cyl_rows`), 조립이 rim 표로 캡과 밴드를 한 엣지에 잇는다(`assembly::reconstruct`).

- **소속은 재지 않고 배열에게 묻는다.** 밴드는 원판 셀의 네 비트를 한 철자로 읽고, 차트의 셀은 그 끝과 벽이 모두 일치할
  때만 챔버를 가진다 — `read_cell::read_bits`·`CellRead::chamber` 의 doc. 광선이 아닌 까닭: 평면 엔진은 「이 경계 조각이
  남나」를 **라벨**(`emit_faces`)로, 「이 성분이 재료인가, 어느 성분 안에 드나」를 **광선**(`point_in_component`)으로 답하고,
  밴드는 앞의 물음이다 — 좌표 도로 `point_in_faces_rat` 는 이 경계를 뒤집지 않는다. 라벨이 원통을 양쪽 피연산자에 세우는
  까닭은 `read_bits` doc, 「z-범위 판정」의 반례(L-노치)는 `bands` 모듈 doc.
- **밴드 경계는 평면 배열이 실제로 방출한 원 경계**(+ 측면 자신의 rim)이고 임의의 ⊥ 클래스로 자르지 않는다 — 캡 면의 원
  루프와 밴드의 rim 이 같은 엣지 핸들이어야 닫힌-셸 가드가 사용 횟수 2 를 본다(`reconstruct` 의 rim 표 주석).
- **조립의 rim 기계**(통째 원마다 `OnSeam` 정점과 자기루프 rim, seam 간선 없음, 잘린 원은 노드 사이의 호이고 그 담체는
  `Wall::Arc` 의 `plane`)는 `reconstruct` 의 rim 표 주석과 `Wall::Arc` doc. 성분 결합에는 **둘째 규칙**이 붙는다 — 캡과
  밴드는 Node 를 하나도 공유하지 않으므로 **rim 키**로 잇는다.
- **감김은 규칙으로 유도한다** — `reconstruct` 의 «Which way a rim circle is walked» 주석; 잠금은 부피의 부호(뚫린 상자가
  `8 − πr²h`, `src/tests/bands/slab_bands.rs` 의 winding lock).
- **딱지는 «소속»을 말하지 «존재»를 말하지 않는다.** 구멍 난 옆면이 피연산자로 돌아오면 rim 의 두 섹터가 글자 그대로 같은
  딱지를 들고, 존재는 자국의 **종류**(`Graze`·`Transversal`)가 답한다 — 둘은 한 호를 볼 때 함께 채운다(`ArcLabel` doc).
  차트의 셀이 소속(`read_bits`)과 존재(`face_spans`)를 **나란히** 묻는다.
- **`face_spans` 의 진리표**(`Seated` 건너뛰기 포함)는 그 doc 이다. `Seated` 건너뛰기는 방어가 아니라 작동 중인 경로다 —
  한 호가 옆면 자국과 평면 면의 `Seated` rim 을 함께 들고 오는 인구가 실재한다. 자국 수의 계약 «제 솔리드의 옆면 자국은
  많아야 하나, 없으면 셀도 없다»는 차트 census(`src/tests/probes/chart_census.rs`)가 단언하고, 둘 이상인 호와 자국 없는
  호의 답은 단위 시험(`face_spans_reads_the_trace_not_the_label`)이 손으로 쓴 자국으로 진술한다.

### 인구 게이트 (`planes::cylinder_gate`)

게이트는 (평면 클래스, 원통) 쌍마다 생성기 하나(`cyl_gate::class_meet` → `ClassMeet`)로 답하고, 모든 쌍에 답한 뒤 그 답을
읽어 짓지 못하는 것을 이름으로 거절하며(`ObliqueCylinderCut`·`TangentLineInAnotherPlane`·`CylinderPairContact`, 기하는
`CylinderGateUndecided`, 폭은 `WitnessNotRational`), 남은 것을 기록으로 넘긴다(`CylRecord`·`Silence`, `SharedRuling`,
`GrazeLine`, `Tangency`). 갈래와 거절 이름은 `cylinder_gate`·`class_meet` 의 doc(비스듬한 평면은 `laterals_clear_plane`·
`oblique_section_meet`, 경계 단면 대 평면의 판정은 `boundary_stand` 하나). 조각끼리의 물음은 `pieces_clear` 하나
(`Footprint`, `separating_dirs`, `lateral_reach`, `lateral_box`)이고, 면의 경계 조각은 `planes::corners::corner_of` 가 읽는다.

- **한 클래스가 드는 원통의 닫힌 자취는 축이 평행한 원통들의 것뿐이다** — 원 셀 둘의 중첩을 답하는
  `nesting::disk_in_disk` 가 축이 평행할 때만 참이라서다. 비스듬한 단면을 드는 클래스가 다른 축의 단면을 들거나 다른 원통의 축에 수직이면 그 원통이 클래스에 닿든 말든
  `ObliqueCylinderCut` 으로 거절한다(`carries_unnested`).
- **세계 서술이 없는 클래스(돌린 벽)는 말 없는 기록으로만 통과한다**(`Silence`, 묻는 순서는 `class_meet` doc). 먼저
  «원통 전체가 한쪽인가»를 묻는다(`nameless_stance` — 증인 둘, 돌리기 전 틀의 정확한 물음; 쪽은 `WorkingPlane::misses` 한
  곳에서 만들고, 한 원통의 옆면은 함께 묻는다 — `WorkingCyl::clear_of`). 쪽은 더 강한 증명이라 `Apart` 보다 먼저 묻고
  기록해, 그 원통의 관통점·호의 쪽을 구간 판 대신 그 기록이 답한다. 기록(`WorkingCyl::records` — `is_silent`·`beside`)을
  읽는 자리는 원 싣기(`circle_on_class`), 현(`chord_stations`), 선분과 원·룰링의 교차(`split_circles`·`split_rulings`),
  원판과 꼭짓점의 포함(`nesting`), 광선(`curved_count`·좌표 도로)이다. 상자는 그런 클래스가 있을 때만 짓는다.
  - `Silence::Beside` — 관통점·호를 그 평면에 대해 세우는 자리가 계수 없이 기록된 쪽을 읽는다(`side_of`, `PierceAt`,
    `rays::silent_table`, `arc_crossings_of` 의 doc).
  - `Silence::Apart` — 면이 옆면 상자와 떨어지고(`faces_apart`, `cache_reach`) 원통이 온전한 림만 남기면(`only_whole_rims`)
    그 클래스 위에서 원통을 구멍을 메운 / 보스를 걷어낸 단면으로 읽고, 다른 기록이 있어도 선다. 논증·측정과 다른 기록과의
    교선을 묻지 않는 위험은 `Silence::Apart` doc 과 `cylinder_gate` 본문.
  - `Silence::Sectioned` — 평면이 축을 따라 반지름 안이고 면이 떨어지면 단면을 제 배열에서만 참으로 그린다(룰링과 현
    `chord_stations`, 관통점은 `PiercePoint::Bounded`, 「이차곡면 판정 — 원통 우선」). 메운 읽기보다 먼저다(「가지 말 것」).
    그 룰링은 옆면 차트로 가지 않고, 클래스가 내는 면이 그 룰링에 닿으면 `per_class` 가 거절한다. 다른 기록 클래스와의
    교선은 옆면 면 밖(`line_off_lateral`)이거나 이 클래스의 경첩(`hinge_ruling`, `Model::plane_hinge` — 이름 없는 평면이
    정확히 품는다고 보일 수 있는 세계의 유리수 선은 경첩뿐이다)이어야 하고, 아니면 거절한다. `nesting` 의 원판 지름길이 이 클래스의 현에
    기대는 전제는 `nesting` 그 자리 주석.
  - `Silence::Grazes` — 돌리기 전 틀에서 접한다고 증명된 클래스는 림을 묻지 않고(`Apart` 보다 먼저) 면이 떨어지면
    통과한다. 다른 기록과의 교선은 단면 클래스와 같이 정산하고, 단면 클래스와 공유하는 접선은 `GrazeLine` 으로 기록해
    이름 표가 묶는다(`Aliases::record_graze_line`; 범위 밖의 점은 묶지 않는다 — 「가지 말 것」). 이름 있는 짝은 거절하고,
    그 거절이 막는 결함보다 넓다는 측정은 `cylinder_gate` 본문의 그 자리.
- **게이트의 «기록»이 열고 배열은 기록만 믿는다.** `point_plane_clearance_rat` 한 호출의 삼분법이 생성기의 세 답(무기록 ·
  `CylRecord::Cuts` · `CylRecord::Touches`)이고, 기록 집합은 증명 캐시라 배열의 독자는 전제를 `debug_assert` 로 든다.
  **접선을 절단으로 기록하지 않는다** — 절단의 명제(「평면이 반지름 안쪽을 지난다」)가 거짓이 된다(`class_meet` 본문,
  `Tangency` doc). 면이 닿지 않아도 클래스가 옆면을 자르면 기록하는 까닭은 `class_meet` 본문. 캡이 상대 재료 안에 앉는
  변형은 짓는다(`RejectReason::RulingBoundNotYet` 위의 주석).
- **«두 클래스의 교선이 옆면 면 위에 있다»도 게이트가 한 번 기록한다**(`SharedRuling`, `WorkingCyl::shared`). 짝은 생성기가
  이미 답한 것에서만 나오고, 선은 무한 곡면이 아니라 옆면 **면** 위여야 한다(`line_apart`). 씨앗·룰링 분할의 접기·차트의
  station·접선 판정·결과의 간선 분할이 이 기록을 읽고 스스로 판정하지 않는다 — 짝과 독자, 기록이 드는 점(`at`·`holds`·
  `stated_by`)은 `SharedRuling` doc. `Tangency::line_on_another_face` 는 같은 계산(«그 접선을 품는 클래스들»)에서 출발해 그
  클래스들의 **면**을 읽는다.
- **게이트가 물어야 하는 명제는 «두 클래스가 «면»을 공유하는가»다** — 무한 곡면의 거리가 아니다(「점의 이름」의 첫 원리).
  곡면 거리로 못 비킨 쌍은 옆면이 자신을 대변하고(도달 범위 `lateral_reach` — 수직은 그 일반식의 한 경우이지 분기가
  아니다), 같은 솔리드 안의 쌍은 묻지 않는다. 원통 쌍의 차례(`r₁+r₂`, `same_surface`, 중첩, 면 물음, 남으면
  `CylinderPairContact`)는 `cylinder_gate` doc.
- **«두 조각이 점을 공유하지 않는가»는 물음 하나이고 함수 하나다**(`pieces_clear` — 조각, 증인 셋, 엄격 포함은 그 doc).
  묻는 쪽은 쌍 규칙(면 쌍마다), 배열의 간선 그물(`CircleCrossesRuling`), 그리고 원통 위의 선을 읽는 쪽(룰링 분할의 원통
  위 선분 포함)이다. 선을 읽는 쪽 가운데 «면 위 / 모름»을 갈라야 하는 `SharedRuling` 기록과 접선 행은 같은 증인을
  `line_apart` 로 세 갈래로 읽고 `pieces_clear` 는 그 «떨어짐»만 쓴다 — 증인을 하나 더하면 이 자리들이 함께 연다. 분리
  방향 후보와 근호 비교 `g > √p + √q` 는 `separating_dirs`·`reaches_apart` doc.
- **발자국은 «꼭짓점»이 아니라 «경계 조각의 도달 범위»다.** 꼭짓점(폭 0)·원판(`중심 ± ρ`)·호(`Corner::Round`, 범위는
  `arc_extent` — 원통 옆면과 한 벌, 끝 순서는 담체 축에 대해 CCW)·관통 코너(`X + Y√c`, `rank` 라 재진술이 없다)를
  `corner_of` 한 함수가 읽고, 실패는 `CornerFail` 이 원인별로 이름 짓는다. 못 읽은 면은 «무죄»가 아니라 «모름»이다.
- **피연산자가 굽은 이름을 나를 수 있다 — 연쇄의 규칙.** 결과가 다음 불리언의 피연산자가 되면 게이트의 그물
  `CurvedOperandBoundary` 가 트레이서(`loop_triples`)와 같은 물음을 담체로 묻는다(그 이름 위의 주석). 게이트를 여는 것은 도로를
  짓는 작업의 마지막 걸음이다 — 「그물이 덮는 만큼만 연다」는 등식은 예외를 만난다(`OnSeam` 디스크 면 · subdivision 의
  직진 정점 · 담체 등가의 동일 surface).

### 원·호·조각

스케치의 원·호가 프리즘 빌더에서 원통 옆면과 관통점이 되고(`ops::prism`, kit `build.rs`), 클래스 위의 자취는 축에 수직이면
원, 비스듬하면 단면이다(`planes::cyl_geom` 의 `class_carries_circle`·`trace_band`·`arc_extent`,
`arrangement::cyl_trace::trace_on_face`). 한 클래스가 원과 룰링을 함께 드는 규칙은 `per_class` 의 `ClassEdges::of` 가 진다.

- **`cylinder()` 는 원 스케치 + 돌출의 설탕이다**(kit) — 원시체 연산은 없다. 프리즘 빌더가 원·호를 원통 조각 벽과 관통점
  꼭짓점으로 세우고(호–호 접합은 원통 둘 위의 점), 사분각 밖의 호는 `refuse_non_quarter_arcs` 가 거절한다.
- **호의 진실은 끝점이지 각도가 아니다**: `Edge2d::Arc { center, r2, ccw }` 의 시작·끝은 링의 정점이 말하고(전부 유리수;
  온전한 원은 시작 = 끝), 회전각은 커널이 모르는 말이며 끝점 표현이 더 넓다. 옆면의 향과 시계 방향 호의 저장(`[A, B]`
  반시계 규약)은 `ops::prism` 의 `along`·`push_walls` 주석.
- **구간 돌출은 «스케치가 여기서 시작한다»다** — `[lo, hi]` 를 한 번 쓸고 프레임의 원점을 구간의 시작에 둔다(kit `build.rs`
  의 `Dist::Both`): 두 번 쓸어 fuse 하면 한 곡면이 축의 부호만 다른 두 진술이 되어 거절되고, 쓸고 옮기면 결과에 모션이
  남는다. 계약은 `lo < hi` 만 요구한다. 커널 쪽 「축 부호만 다른 쌍둥이」 정준화는 **짓지 않는다** — 그 인구는 위 규칙이
  없애고, 진술을 합쳐도 「한 곡면을 두 피연산자가 소유」로 옮겨 앉을 뿐이다(차트가 「원통 클래스는 한 솔리드의 곡면」을
  전제로 선다 — 그 전제를 바꾸는 것은 **별개 능력**).
- **원 셀은 «클래스가 축에 수직일 때만» 원이다** — 이름은 `class_carries_circle` 하나이고, 원을 읽는 소비자 목록과 비스듬한
  호를 그 도로들 안에 묶는 게이트의 논증(앉은 원 생산자의 더 강한 근거 포함)은 그 doc. 비스듬한 클래스의 자취는
  **단면**(타원)이고 게이트가 `CylRecord::Sections` 로 기록한 쌍에만 선다. 이 쓰임이 서는 까닭은 **평면 단면이 모두 θ 의
  그래프**라는 것이다(`cyl_chart` 모듈 doc): θ 순서(`circular_order_about_seam`)·호 포함(`theta_between`)은 단면에도 참이고,
  원에만 맞던 «점에서의 반경에 축 성분이 없다»는 `nacre_exact::quad::axial_side`·`ccw_tangent_side` 한 술어로 묻는다
  (극점 평면의 `g*` 는 `combinatorics::winding` 의 그 주석). 호의 상자는 `trace_band`. 추적은 ⊥ 클래스와 한
  덜어내기(`trace_on_face`; 꼭짓점은 `corner_on_section`, 가로지름이 없으면 `Theta::misses_an_angle`)이고, 그리지 못하는
  단면은 기권 하나(`DeclineKind::ObliqueSection`)라 불리언은 게이트와 같은 이름 `ObliqueCylinderCut` 으로 거절한다 — 기권의
  목록은 `trace_on_face`·`cycle_on_class` doc. 게이트는 평면이 옆면 경계의 어느 단면에도 한 점으로 접하지 않을 때 단면을 들이고(면이 옆면에
  닿는지는 묻지 않는다 — `oblique_section_meet`), 단면의 딱지와 절단점(`CutRim`)은 차트로 가 다른 선처럼 칸을 가른다. 결과 면이 단면을
  따라 경계 지어지면 결과에 타원 간선이 남는다(`Model::derive_edge_curve` — 원을 그 평면에 기울인 곡선). 단면을 따라 끝나는 옆면에는
  그 단면을 경계로 든 평면 면이 짝으로 있다(차트 census 가 단언한다). 결과를 읽는 쪽은 단면 위 호를 원과 같은 길로 읽는다 — 덩어리를
  가르는 옆면 패리티는 림 평면의 쪽(`in_faces::rim_up`), 면을 잇는 열쇠는 호의 평면까지(`face_components` 의 `JoinKey::Arc` — 캡의 원호와
  단면의 호가 같은 두 점 사이에 선다). 타원 간선을 든 솔리드는 다시 피연산자가 된다 — **옆면 면은 룰링과 평면 단면으로 경계 지어지고,
  면이 어디 있는지를 읽는 자리는 경계 단면을 두 단면의 한 물음(`section_meet`·`section_height_order`)으로 읽지 정거장(축 위 높이
  하나)으로 읽지 않는다** — 원의 높이 비교는 그 ⊥·⊥ 경우다. 원 셀의 림 증인은 축에 수직인 클래스에만 낸다 — 비스듬한 클래스의 단면에 내면 `Witness::On` 이 거짓이 되어
  거절이 아니라 **조용한 오답**이다(`nesting` 의 그 주석). 「커버리지가 끝난다」(`NotSupported`)와 「닿을 수 없어야 할 백스톱이
  말했다」(`SuspectedDefect`)는 다른 문장이라 거절 이름을 재사용하지 않는다.
- **한 클래스는 원과 룰링을 둘 다 들어도 되고, 둘이 가로지르면 안 된다.** 규칙은 «간선»에 대한 것이다 — 클래스의 간선은
  어떤 면의 경계이므로 서로 다른 원통의 두 간선이 만나면 두 옆면이 한 점을 공유하고, 쌍 규칙(`lateral_faces_clear`)이 그것을
  부정한다. 배열의 그물(`circle_crosses_ruling` → `CircleCrossesRuling`)은 두 간선을 조각으로(원이 아니라 호, 직선이 아니라
  선분) `pieces_clear` 에 먼저 묻는 **backstop** 이고, 출하한다 — 만들어지지 않은 노드는 패닉이 아니라 조용히 틀린
  솔리드다(`ClassEdges::of` 본문, `crossing_measured` doc). release 는 `debug_assert` 를 컴파일하지 않으므로 «되던 것이
  깨졌다»는 조용히 지나가던 것이 거절로 드러난 것일 수 있다. **접선은 «가로지름»이 아니다**: 그물은 경계를 인자로 든 문
  (`cylinder_ruling_reached_extent` 의 `touch_counts`)에 열린 물음을 하고, 이름이 그 명제를 말한다(`CircleCrossesRuling` doc).

### 점의 이름

한 점은 여러 도로에서 여러 이름으로 태어나고, 표가 그것을 하나로 접는다. 짓는 규칙은 `combinatorics::names`
(`canonical_triple`·`names_on_cylinder`·`pin_for`·`restate_pierce`)에, 접는 표는 `arrangement::aliases`(`Aliases`·`Lines`)에
있다. 배열 뒤에는 조립이 결과 정점을 닿는 면들로 다시 이름 짓고(`assembly::naming`) 직선 간선의 열쇠와 담체를
정한다(`assembly::reconstruct`).

- **이름은 곡면이 짓고, «그 면이 여기 있나»는 면이 답한다.** 이름(`NodeId`·`Vertex` 의 정의)은 무한 곡면을
  부른다 — 정확한 좌표와 같음·병합이 거기서 나온다. 꼭짓점은 제 면이 없는 곡면 위에도 우연히 놓이므로 «있음»은
  따로 묻는다 — 배열 전에는 옆면의 발자국(`Footprint`, 「인구 게이트」)이, 배열 뒤에는 결과 면의 링이 답하고, 결과
  정점의 정의는 그 면들의 곡면만 부른다(조립의 끝 검사와 validate 가 같은 이름 `VertexNamesAbsentSurface` 로 묻는다).
- **평면 집합이 지나는 점의 이름은 함수 하나가 짓는다**(`combinatorics::canonical_triple`). 원장은 타입이다 —
  `Canon3` 만 `NodeId::three_planes`·`Def::Three` 가 받으므로 다른 철자는 컴파일 오류다.
- **원통 위의 점도 이름 하나다** — `combinatorics::names_on_cylinder` 가 짓고, 씨앗(피연산자의 코너, 게이트의
  `SharedRuling`·`GrazeLine`)은 모두 `Aliases::record_cylinder_point` 로 간다. 대표는 관통 이름, 그 가운데 정확한
  이름이 먼저다(`Aliases::rep_rank`). 소비자는 `canon_point` 로 비교한다.
- **이름의 같음은 수렴한 표에서만 묻는다.** 라운드 안에서 정확히 증명된 동점은 거절하지 않고 기록하고
  (`Aliases::record_same_point`), 라운드가 거절하는 것은 같음이 아닌 것뿐이다 — 한 쌍의 두 근, 두 원통의 곡면 위의
  점(`Aliases::join_points`), 한 경계가 한 점을 두 번 자름(꼭짓점 닿음 — `cyl_trace` 의 `CylHoleFeature`). 수렴한
  표를 읽는 단계(pass B)는 어느 라운드도 기록하지 않은 동점을 거절한다.
- **핀은 (이름, 선)의 함수다**(`combinatorics::pin_for`; 관통 이름은 제 짝의 두 평면이 다 품는 선에서만
  위치자다 — `names_the_line`). 피연산자의 관통 코너는 방문이 아니라 제 정의로 이름 짓는다
  (`pierce_name_from_def`). 대표 이름이 다른 쌍의 관통 이름인 점(공유 룰링 위)은 세 평면의 교점이라 유리수이고,
  그 근은 유리수로 줄여 읽는다(판별식이 제곱수인 근을 한 걸음의 다른 근호와 섞지 않는다).
- **클래스 원과의 교차는 모든 호 밖으로 증명되면 절단점이 아니다**(`split_circles` 의 원 팔 주석). 증명하지 못한
  교점은 절단점으로 남는다 — 호 밖이면 정점 하나가 늘 뿐이지만, 호 위의 교점을 빼면 정점을 잃는다.
- **룰링이면서 평면 쌍 선인 선은 «옆면이 거기서 끝나는가»가 어휘를 정한다** — 끝나면 평면 어휘
  (`rulings::lateral_half_of_plane_pair_lines`; 접하는 클래스는 `CylRecord::Touches`), 가로지르면 룰링 어휘
  (`split_rulings` 의 접기). **선은 그 위의 모든 꼭짓점에서 잘린다**(`cut_rulings` — pass A 는
  `cut_rulings_along_segments`).
- 그 밖의 이름 자리: «어느 룰링인가»는 점과 벽의 함수(`ruling_side_signed`), 직선 간선의 열쇠는 두 끝이고 담체는
  쓰는 두 면의 곡면(`EdgeKey::Line`·`JoinKey::Line` — `assembly::reconstruct` 의 그 자리 주석; 한 직선의 두 철자는
  `dissolve_straight_angles`), 결과 정점은 닿는 면들의 곡면으로 짓는다 — 관통 정점도(`assembly::naming` 의 결과 정점
  주석).
- **굽은 링의 코너는 이름을 «유도»하지 않고 «재진술»한다** — 정점의 원통 둘레 순위(`rank`)와 class 이름의 교선 위
  순서(`QuadRoot`) 사이의 문은 `pierce_name_from_def`·`Def::vertex`(`seam_ranked`)이고, 담체는 `combinatorics::Wall` 로
  적는다(두 번째 어휘를 만들지 않는다). 다발의 선은 쌍 하나로 부르고(`Lines`, `combinatorics::line_sense`), `ℓ` 의
  방향 반전이 홀수 번이면 root 를 뒤집는다(`QuadRoot::canonical`).
- **법선만 읽는 자는 `PlaneName` 의 네 계수 정준화를 물려받지 못한다** — `primitive_normal`, 그리고 평면 안에서
  직교해야 하는 차트 축(`Chart2dRat` — 버린 세 차트의 측정 포함).

### 결정은 `def`, 측정은 `cache`

원통 행(`WorkingCyl`)과 평면 행은 정확 서술과 f64 실현을 함께 들고, 결정은 서술만 읽는다 — `WorkingCyl::realized`·
`Wall::Ruling`(`side`·`up`)·밴드의 station·축 부호 `plus_t_is_above`·`Judge::stored_world_rat` 이 이 규칙 위에 선다.

- 규칙은 `WorkingCyl::realized` 의 doc 이 든다. `cache` 를 읽어도 되는 것은 측정뿐이다(공차·rim 의 중심).
- **축 부호는 «어디서»가 아니라 «어느 쪽»을 묻는다** — `world_rat` 은 이름이라 방향을 안 나르고, 향
  (`WorldName`·`Model::world_plane_name_sense`)이 나른다. 잰 수는 `assembly::reconstruct` 의 `axis_up` 주석.
- **f64 경로를 정확 경로로 바꿀 때의 규율**: census 비트 동일은 「출력이 안 변했다」이지 「두 경로가 매번 같은 답을
  냈다」가 아니다 ⇒ **차등 탐침을 먼저**(둘 다 계산·어긋나면 패닉), 코퍼스 전량, **탐침 자체를 부정 대조**, 그다음
  f64 제거.

### 분할과 자국

한 클래스의 배열은 세그먼트를 원·룰링과 만나는 점에서 자른다(`arrangement::per_class` → `arrangement::crossings`의
`collect`·`Landing`·`split_segments_at`·`order_distinct`, 어휘 `Split`·`Piece`, 순서 `combinatorics::order_on`). 링을 선에
대어 읽는 걸음은 `combinatorics::ring_walk`(`ring_against_plane`) 하나다.

- 선분 대 곡선은 한 단계를 곡선 종류마다 한 번씩 돈다(원 먼저 — `per_class`); 축을 따르는 선분은
  `planes::line_off_faces`(게이트와 한 술어)가 답한다. 옆면의 구멍은 다른 면과 똑같이 이름 짓는다(`lateral_cycles`·
  `lateral_shape`, 거절 `DeclineKind::CylFaceHole`).
- **선 위 연속 구간은 굽은 간선을 가로질러 잇지 않는다** — 걸음이 구간을 departure 에서 쪼갠다(`EdgeMeet::Departs`).
- **패리티는 σ 와 무관하다(정리).** 원과 평면은 두 점에서 만나므로 호의 내부는 선을 다시 안 만난다 ⇒ 호는 통째로
  한쪽(σ)이다. 교차 수는 `(s_a ≠ σ) + (σ ≠ s_b)` 이고 패리티는 어느 쪽이든 `s_a ≠ s_b` 다. 정리는 개수만 주고 위치는
  안 주므로, 위치가 필요한 자리에서는 σ 를 쓴다 — 오라클이 있기 때문이다(재연산의 split-run 인구; 부호를 뒤집으면
  `LabelConflict`). 떠나는 간선은 σ 쪽 «가상 노드»로 걷기에 들어간다(`EdgeMeet::Departs { side: σ, .. }`, σ = −`ruling_side`).
  선을 따라 떠나는 호(접선 성분 0)의 σ 는 단면의 중심 쪽이고(`combinatorics::arc_leaving_side` — 혼합 링의 접하는 반직선도 읽는다), 닫힌 호면 선을 다시 건너지 않는다
  (`ring_against_plane`) —
  잰 재연산 셋(접점에서 집힌 판 윗면을 그 접선으로 자르는 — `the_pinched_ring_re_cuts_along_the_rims_tangent`)에서는 그 점의 두
  run 이 접촉을 나눠 가질 뿐이라 이 부호를 뒤집어도 답이 같고, 오라클은 기하로 유도한 단위 시험
  `an_arc_leaving_a_line_along_it_lies_on_its_centres_side` 다.
  **오라클 없는 부호는 만들지 않는다** — 부정 대조가 전부 초록인 채로 실제로 틀린 부호가 가능하다.
- **호의 교차 수는 호 자신에게 묻는다**(`ring_against_plane` 의 `crossings`; 평면 면의 도로들은 `ClassLine` 을 거쳐
  `arc_crossings`). «간선이 그 위인가»는 도로마다 다르다 — 평면 면에선 직선, 옆면의 구멍에선 원(`on_meet`).
- `UnorderedEdges` 가 뜨면 케이스를 열지 않고 상류(`merge_coincident`/`Aliases`)를 본다.

### 룰링

게이트가 기록한 (벽, 원통) 쌍(`WorkingCyl::draws_rulings`·`touches`·`SharedRuling`)에서 배열은 밴드를 잘린 원에서
끊고(`split_rims`) 룰링 조각을 낸다(`arrangement::rulings::rulings_on_class`). 차트(`cyl_chart::chart_of`)는 그 조각을
세로선(`ThetaSeg`)으로, 클래스의 자취를 선으로 얹는다(`cyl_chart` 의 모듈 doc).

- **할선 벽들이 지나는 룰링은 정거장 하나다** — `stations_of_lines`·`Station::chords`·`ChordSays`·`Chart::slit_at`.
  공유 룰링은 할선 벽들에 면이 없는 높이에서도 차트의 정거장이다 — 차트는 곡면 전체에 선을 깔고 칸을
  분류하며(Requicha 의 분류, Mäntylä 의 이웃), 그 선에서 옆면이 끝나는지는 기록이 아니라 이 분류가 정한다.
- **옆면 위의 곡선은 룰링이거나 θ 의 그래프다** — 높이 순서 `quad::section_height_order`(래퍼
  `combinatorics::line_order`), θ 순서 `circular_order_about_seam`(두 점을 같은 원 위에서 받지 않으므로 다른 벽의 룰링도,
  위치가 무리수(`√(r²−d²)`)인 오프셋 벽(`0 < d < r`)의 룰링도 같은 산술로 비교한다 — 지름을 가정하지 않는다), 코너의 각
  순서 `departure_key`. 비스듬한 단면 선을 따라 끝나는 칸은 결과에 타원 간선을 남긴다(「원·호·조각」).
- 룰링의 끝은 언제나 선 위이고(`ThetaSeg::on`), 셀은 slab × 틈이다. 룰링 구간은 `rulings_on_class` 가 정하고 차트는
  다시 계산하지 않는다. 이름은 `Chart::ruling_name` 한 번, 동일성은 `ruling_side`(`RulingCarrier::side`)다.
- **접선 룰링은 `side = 0` 이라는 «이름»이다**(`QuadRoot::Double`) — 따라오는 것: `tangent_codirected_turn`,
  `End::Uncovered`, `SegKind::Tangent`.

### 차트와 셀 방출 (`arrangement::cyl_chart`)

옆면은 차트의 «영역»이고 θ-패널이 셀이다. 차트는 셀을 읽고(`Chart::read_cell`) 영역을 걸어(`cyl_chart::regions`)
옆면을 낸다(`emit_lateral`). 결과 전체 위의 정리 한 번(`arrangement::whole_result_faces`)이 옆면과 이웃 면의 점을
함께 지우고, 조립은 정리 뒤의 림(`HeldRims`)만 읽는다.

- **«같은 엔진»은 같은 «단계»이지 같은 «코드»가 아니다**(`regions` doc; DCEL 이식은 「가지 말 것」). **정점은 이웃
  클래스가 가진 곳에만** 선다 — 직진하는 정거장에 정점을 만들면 벽 면의 간선 하나에 옆면 간선 둘이 대응해 용접이
  깨진다. 꺾이지 않는 점은 `whole_result_faces` 가 옆면까지 모든 링에서 함께 지우고, 조립은 그 뒤의 림을
  덩어리마다 읽는다(`HeldRims::cut_per_group`; 나눈 표 `CutRims` 와 타입이 갈린다).
- **다 지워진 림의 철자는 면이 정한다** — 캡은 `Bound::Circle { cyl: k }`, 옆면은 `Bound::Rim { plane: c, ccw }`
  (Parasolid 도 루프의 종류를 그 면의 곡면에서 읽는다; `dissolve_straight_angles` 의 `closed`).
- 칸의 챔버는 양 끝 원의 가로 딱지로 서고, 세로선은 검사자다 — 말하는 변이 모두 합의할 때만 챔버가 선다
  (`Chart::read_cell`). 사건은 `Chart::events`(`combinatorics::theta_groups`), 닫힌 간선 `[v, v]` 는 한 바퀴로
  읽는다(`ring_floor`·`arc_span`·`Model::derive_edge_curve`). 원을 극점에서 단조 조각으로 자르지 않는다(「가지 말
  것」의 seam 엣지 행).
- **호의 원판 쪽 반간선은 `axis_up · frame_sign` 으로 유도하고, 셀의 유리 코너가 감시한다**(`disk_side_probe`) —
  `arrangement::emit` 의 호 라벨 주석.
- 옆면의 패리티는 모든 루프에 대한 `loop_parity` 이고, 광선과 옆면의 교차는 근마다 면에게 묻는다
  (`lateral_face_crossings`; 면의 각을 읽지 않는 갈래는 todo 「광선이 원통 위 선으로 옆면을 스치면 면의 각을 읽지
  않는다」).

### 감김과 방향

링의 방향은 링 자신의 감김에서 읽는다 — 옆면 고리의 구멍 쪽은 `arrangement::cyl_trace` 의 `material_theta_sign`,
옆면 구멍 도로는 `cycle_on_class`, 감김수는 `combinatorics::winding` 의 `loop_winding`, 혼합 링은 `ray_step_crossing`·`Ring::edges`·`point_in_mixed_ring`.

- **「어느 호가 구멍인가」는 링 자신의 감김이 답한다** — «재료는 진행 방향의 왼쪽», 철자는 `material_theta_sign`
  하나. 호의 방향은 생산자가 벽에 싣는다(`Wall::Arc` 의 `ccw`).
- **감김수는 «영역의 극점»에서 규칙 하나로 읽는다**(`loop_winding` — 극점, `angular_order`, `c1`,
  `smooth_extremum_winding`). 국소 모양마다 갈래를 덧대지 않는다(「가지 말 것」의 케이스워크 행).
- 부호 인자의 잠금은 `a_holes_arcs_run_the_way_the_hole_lies`·`a_bores_lateral_with_a_hole_takes_a_second_cylinder`·
  `smooth_extremum_winding` 의 doc.

### 증인과 중점 (`nesting`)

셀이 셀 안인가는 모듈 하나(`nesting.rs`)가 증인 하나·판정 하나·루프 하나로 답하고, 대답의 어휘
(`combinatorics::Said`·`WitnessTally`)를 3차원판(성분의 깊이, `assembly::grouping`)과 나눈다. 증인 공급의 원시
함수는 `combinatorics`(`edge_interior_points`·`conjugate_midpoint`·`circle_centre_rat`·`witnesses`)와 `nesting` 의
`rim_and_centre`·`ring_own_circle` 이다.

- 실패의 이름(`RingHasNoWitness`·`NoClearRay`·`WitnessNotRational`·`CylinderGateUndecided`)과 도로가 대상을 못
  읽는 사건(`Said::Declines`)은 `Said`·`WitnessTally`·`Said::of_ray_error` 의 doc.
- **말 없는 클래스가 낀 물음은 게이트의 증명으로 답하고, 못 하면 게이트의 이름을 든다** — 정확한 도로가 못 읽는
  링 물음은 «링이 든 상자 밖이면 밖» 한 규칙이 2차원·3차원에서 함께 답한다(`combinatorics::ring_box_says`·
  `silent_name_by_box`).
- 링의 닫힘이 암묵인 관례는 「Tessellation 층」 절과 `Tessellation` 의 doc.

### 접선·접촉·집힌 면

접선은 게이트가 기록하고(`tangencies`, `planes::Tangency`) 배열은 무시하며, 조립의 판정 하나
(`assembly::self_touch::tangency_reject`)가 연산을 `keep` 으로 읽는다. 접선이 남긴 집힌 면은 메시가 다룬다
(「Tessellation 층」 절).

- **접선은 «가르지 않는다».** 게이트는 접선 쌍을 거절하지 않고 `tangencies` 로 기록하고, 배열은 그 선에서 어떤
  셀도 나누지 않고 차트의 어떤 섹터에도 정거장을 세우지 않는다. **셋을 가르는 것은 연산이고, 연산은 `keep` 으로만
  들어온다**(`self_touch::pinch`).
- **«갈라진다»가 곧 결함은 아니다** — 판정은 «두 덩이 ∧ 한 출력 솔리드»(`tangency_reject`, `SelfTouchingResult`).
  셋째 클래스의 면이 그 선에 닿으면 게이트가 거절하고(`TangentLineInAnotherPlane`), 그 선이 피연산자 간선이면
  구조에 맡긴다(`Tangency::line_is_an_edge`).
- 집힌 것은 «표면»이 아니라 «면»이고, χ 는 그 증거가 아니다(접선은 정점을 만들지 않는다 — 시험
  `the_touch_is_a_slit` 의 doc). 처방은 b-rep 이 아니라 메시의 기호적 순서다(「가지 말 것」의 집힌 면 분할
  행). tess 의 거절 이름은 메시에 대한 말이고 솔리드 판결은 `validate` 의 것이다.

### 메시로 재기

굽은 면 메시의 오라클을 정한다 — `nacre-tess` 의 띠 시험과 `nacre-ops` 스위트의 넓이 단언이 이 규칙을 따른다.

- **굽은 면의 메시를 보는 눈은 «넓이»다.** `validate`·watertight·정확 부피·면 수는 옆면의 면 구조가 바뀌어
  **메시가 조용히 틀려도 전부 초록일 수 있다**(넓이 7.873 vs 정확 9.425 인 메시가 그 넷을 전부 통과한다). 그걸 잡는
  오라클은 **메시 넓이 ≈ 정확 넓이** 하나다 ⇒ 옆면 표현을 바꿀 때는 **먼저** 원통 코퍼스에 그 오라클이 서 있어야
  한다.

### 강체 운동과의 교환

이 법칙의 오라클은 `nacre-ops` 의 시험 `src/tests/transform.rs`·`src/tests/suite/commuting_oracle.rs` 다.

- **법칙 — 불리언은 강체 운동과 교환한다**: `transform(rigid(R, t)) ≡ transform(T) ∘ transform(R)` — 세계가 같다. 저장은 각 걸음이 옮겨 적을 수 있는 것을 따르므로 진술의 점 셋은 `invariant` 재진술로 갈릴 수 있고, 노드 구조는 넘침에서만 갈린다(`an_overflow_splits_the_representation_not_the_world`).
  **frame −1 의 메커니즘은 회전이 아니라 세계 평면 씨앗이다** — 이동 하나로 ∥ 딱지가 전부 뒤집힌다.

## 원통의 진실과 seam 규약

### seam 엣지 방식

주기 곡면(원통·구·토러스)의 옆면을 b-rep 로 담는 방식은 둘이고, 둘 다 유효한 AP242 이며 유효한 CW-복합체다.

- **B — 림 두 루프(커널의 표현).** seam 엣지 없이 위·아래 원 두 개로만 옆면 경계를 이룬다(옆면이 두 루프 — 한 림이 바깥, 다른 림이 안쪽). 오일러는 `V−E+F−L_i = 2−2+3−1 = 2` 로 통과한다. 커널의 모든 생산자 — 원 프로파일의 돌출(프리즘)과 불리언 조립 — 가 이 모양으로 짓는다. 어느 림이 바깥인지는 생산자가 정한다: 불리언 조립은 차트에서 감김 +1 인 고리(축 방향 아래 림), 프리즘은 스케치 평면의 림(원통의 축이 스케치 법선이라 음의 돌출이면 위 림)이고, 안쪽 루프는 다른 림과 구멍이다. 림은 통째 원이면 닫힌 간선 하나(그 정점이 θ = 0 의 `OnSeam`), 잘렸으면 호와 룰링의 사슬이다. STEP 출력은 내부 표현 그대로다.
- **A — seam 엣지(짓지 않는다).** 위 원 + 아래 원 + **세로 seam 직선 엣지**, 옆면 = 4-엣지 단일 닫힌 루프 `[bottom, seam, top⁻, seam⁻]`(seam 이 같은 면에서 반대 방향으로 두 번 쓰인다 — self-adjacent). 「가지 말 것」에 한 행.
- 정점조차 없는 제3의 형태(끝점 없는 독립 전체 원)는 오일러가 깨지며, **커널은 그 형태를 표현하지 않는다.**

A 를 받치던 다섯 근거가 이 커널에서 서지 않는 까닭:

1. **불리언.** 조합적 b-rep 불리언(면 루프를 순회하며 교차곡선에서 엣지를 쪼갠다)은 주기 면이 실제 엣지로 둘러싸인 유계 영역이기를 요구한다. 이 커널의 옆면은 그렇게 풀리지 않는다 — 자기 차트 위의 셀 복합체로 답하고(`cyl_chart`), 영역의 경계는 차트 고리 그대로이며, 바깥 고리는 고리의 감김(호 단위의 ccw − cw)이 고른다. 각도 비교는 끊는 자리와 무관하다(`SeamOrder::seam_first`). 업계는 갈린다: OCCT 는 seam 이 필수라 seam 없는 면을 읽으면 끼워 넣고(`ShapeFix_Face::FixMissingSeam`), ACIS 는 `periodic_no_seam` 으로 두 루프 면을 허용하고, Parasolid 는 감는 루프 한 쌍이 네이티브다.
2. **tess 정합.** nacre tess 는 공유 엣지 polyline(`by_edge`)을 경계로 소비해 crack-free 를 얻는다(테셀레이션 절). 옆면을 펼치는 절단은 tess 가 스스로 한다 — 감는 두 림을 한 생성선에서 자르고 그 점을 지나는 루프의 공유 polyline 에 끼워(이웃 면도 같은 점을 본다) 직사각형으로 펼친다(`cut_seamless_bands`). 절단선은 가능한 한 적은 루프를 지난다: 두 루프가 모두 닫힌 간선 하나면 각자의 정점에서(새 점 없음), 아니면 각 림을 한 번·구멍을 0 또는 2 번 지나는 생성선 가운데 구멍을 가장 적게 지나고 기존 샘플에서 가장 먼 것; 지나는 구멍은 끊어 절단선을 따라 바깥 고리에 잇는다(seam 없는 모델러의 패싯팅 — Stallings, *Reduced Topology B-Reps* §10). 매개화의 이음은 캐시(메시)의 일이고 진실(위상)에 오르지 않는다.
3. **«가짜 모서리»의 값.** seam 엣지는 기하가 아니라 «각도를 어디서 재기 시작하나»(`ref_dir`)를 위상에 새긴다 — 잘린 림마다 θ = 0 정점과 그 자리의 호 쪼개기, 구멍이 그 자리에 걸리면 바깥 걸음에 잇는 slit, θ = 0 에서만 보는 집힘 그물, self-pair 담체 규칙이 따라붙고, 같은 모양이 seam 의 자리에 따라 다른 위상이 된다(벽 위 보스의 노치가 seam 에 걸리면 바깥 걸음에 이어지고 아니면 구멍이 된다). B 에서 그 자리는 통째 림의 정점 하나에만 닿는다 — 네 벽의 보스가 한 철자로 나온다(`a_boss_on_a_wall_has_one_lateral_face`).
4. **위상 균일.** 옆면은 구멍 있는 평면 면과 같은 틀 — 바깥 루프 하나와 안쪽 루프들 — 이고, validate 의 오일러는 `L_i` 항을 이미 든다. 어떤 간선도 한 곡면을 자기 자신과 가르지 않는다 — 원통의 자기 쌍은 `push_edge` 가 거절하고(`TwoCylinders`), 평면의 자기 쌍 `(P, P)` 는 조립이 못 한 병합으로 거절하며(`CoplanarMerge`), validate 가 둘 다 잡는다. 한 점에서 고리가 집히는 차트 모양은 영역 걸음에서 위치와 무관하게 세어 스위트 32(캡 원 점수판 20 · 안으로 드는 쐐기 12)·census 0 이고, 지어진 결과의 옆면에는 남지 않는다(꼭짓점을 두 번 지나는 옆면 고리 — 스위트·census 0) — 3-D 에서 선으로 닿는 모양이다.
5. **유리수 차트.** 원통 위 점의 각도 매개변수는 단위원의 유리수 반각(Weierstrass) 매개변수 t 이고, 그 차트는 원에서 **정확히 한 점을 못 덮는다**(t=∞). 그 배제점은 그대로 θ = 0 이다 — `ref_dir` 은 모델 기하이고 통째 림의 정점이 거기 있다. 각도 순서는 그 점을 첫째로 읽어(`seam_first`) 앞/뒤 어느 쪽에 끊어도 답이 같다; 사라지는 것은 그 점을 지나던 «간선»뿐이다.

**어댑터 분리.** 내부 표현 ≠ 교환 표현이다. nacre-step 은 내부 표현을 그대로 내보내고, OCCT 는 seam 없는 옆면을 읽으며 필요한 seam 을 스스로 넣고(`ShapeFix_Face::FixMissingSeam`, STEP 읽기의 기본값) checkshape 유효로 읽는다 — 위·아래 림 정점이 한 생성선에 있지 않은 것, 축을 향한 보어, 구멍이 모든 각을 덮는 옆면까지, 면 수와 부피가 nacre 의 것과 같다(`nacre-oracle` 의 `seamless_laterals_read_as_valid_solids`, 헬퍼의 `valid` 는 양쪽으로 보정). 그래서 내보낼 때 seam 을 넣거나 빼는 어댑터는 없다.

### 원통의 진실

1. **진실 형태.** `nacre_topo::Surface::Cylinder { def: CylinderDef, motion }`, `CylinderDef { origin, dir, ref_dir, r2 }` — 전부 유리수다. 반지름은 **제곱** `r2` 로 들고 그 타입은 **넓은** `BigRat` 이다: 축에서 곡면 위 유리수 점까지 거리의 제곱은 유리수지만 반지름 자체는 아닐 수 있고, 어떤 정확 술어도 제곱하지 않은 반지름을 묻지 않는다. 넓은 이유는 진술된 `Rat` 의 제곱이 `i128` 에 안 들어갈 수 있고(16자리 십진 반지름), 진술은 그 제곱의 폭 때문에 거절되지 않기 때문이다. 수로서의 반지름은 실현이다(`radius_f64`, 제곱근이 유리수면 `radius_rat`).
   - `dir`·`ref_dir` 은 **정규화하지 않은 원시**다 — normalize 는 정확형을 파괴한다(`normal_def` 와 같은 규칙). 캐시는 실현이다 — 세계에 진술되면 문이 진술에서 원점·단위 축·단위 `ref_dir`(축에 수직인 부분)·반지름을 각각 정확 반올림하고(`push_cylinder_raw`, `nacre_exact::cyl_unit_frame_f64`), 진술이 없으면(사분각 밖 회전·프레임 노드) ops 깔때기가 사슬을 고정밀로 재생해 같은 넷을 반올림한다(`push_cylinder_realized`; 방향은 재생한 두 점의 차). 원통 위 seam 정점도 같은 재생으로 세계에서 캡과 만나 실현된다(`seam_point_met`).
   - 성분형이 옳은 이유: 사용자 어휘의 직접 리프트이고 곱이 없는 셔플이다. 「유도된 곱(계수)」 모양이 아니라 평면의 *점* 정의에 대응하는 원통판이다.
   - `ref_dir` 은 프리즘이 스케치에서 진술한다(`construct.rs`): 온전한 원은 자기 한 정점 쪽(중심 → 정점, 반지름이 유리수면 그것으로 나눈 단위 길이 — `Ring2d::circle` 의 정점은 프레임 `(r, 0)` 이라 프레임 `+x̂`), 호의 벽은 프레임의 `x̂`. 열쇠가 진술 그대로라(아래 3) 한 원의 호들이 한 곡면으로 intern 되려면 같은 철자여야 한다 — 호가 자기 끝점이 아니라 프레임에서 seam 을 읽는 이유다.
2. **seam 은 모델 기하이고 영구 고정이다**(+`ref_dir` 방향, θ=0). 차트가 자기 배제점을 seam 위에 놓도록 맞춘다 — 차트가 seam 에 적응하지, 그 역은 없다.
3. **interning 은 보수적이다.** 원통의 열쇠(`SurfaceKey::Cylinder`)는 def **문자 동일**(`ref_dir` 포함) + motion 동일만 합친다. 같은 축·반지름·다른 `ref_dir` 을 합치면 seam 이 갈라지므로, 잘못 합칠 위험이 0 인 키다. 기하 동일성은 술어의 몫이다. 평면의 `flipped` 대응물은 필요 없다(문자-동일 키면 캐시 구성도 동일).
4. **`Vertex::OnSeam([원통, 캡])` 의 정의** = 「rim ∩ θ = 0 생성선」(`origin + r·ê` 를 지나는 룰링, `ê` 는 `ref_dir` 의 축에 수직인 단위 부분). `ref_dir` 이 진실에 있으므로 유일점을 정확히 지시하고, 좌표 캐시는 load-bearing 이 아니다. 캡이 축에 수직이든 비스듬하든 한 식이다 — 축과 캡의 만남(유리수 중심)에서 `ê` 를 축 방향으로 캡 위에 올린 걸음(`realize::seam_lift`, 수직 캡에선 `ê` 그대로)만큼.
5. **seam 은 간선이 아니다.** 옆면은 위·아래 림 두 루프로 닫히고(「seam 엣지 방식」 B), `[s, s]` 는 어떤 곡면에도 담체 쌍이 아니다 — `push_edge` 가 한 원통 둘을 거절하고 validate 가 모든 self-pair 를 `EdgeCarrierMismatch` 로 잡는다.
6. **branch 방향.** 평면∩평면∩원통 = 최대 2점이고, `Vertex::Pierce { planes, cylinder, rank }` 변종이 받는다 — 변종이 자기 진실을 말한다(슬롯 재활용 금지). `rank` 는 원통 둘레 순위라(`seam_ranked`) 평면의 순서와 향이 이름에 닿지 않고, 접점은 교점이 하나라 `rank: 0` 이다 — 다시 적을 일이 없으니 접점을 따로 저장하지 않는다. 배열 안의 class 이름은 교선 위 순서 `QuadRoot`(`Lo | Hi | Double`)로 남고, 그 재정렬 규칙은 `QuadRoot` **한 곳**에 산다 — 쌍을 뒤집으면 `canonical`, 같은 선의 다른 쌍이면 `along`(접점 예외는
`flipped(Double) = Double` 에서 저절로 나온다).
7. **원통 위의 꼭짓점은 룰링으로 이름하고 실현한다.** 관통점(`Vertex::Pierce` — 평면∩평면∩원통)과 원통 둘 위의 점(`Vertex::CylinderPair` — 평면∩원통∩원통)은 한 규칙이다: 점은 원통의 룰링 하나 위에 있고, `rank` 는 그 룰링의 원통 둘레 seam 먼저 순서(`nacre_exact::quad::seam_ranked` — 한 철자)이며, 남은 평면이 룰링 위의 점을 고른다. 룰링을 주는 것은 관통점에선 두 평면의 교선을 축 방향으로 쓸어 만든 벽(`realize::wall_rulings` — 교선의 두 만남이 그 벽의 두 룰링 위에 있다; 두 평면 가운데 하나가 축에 평행하면 교선이 그 안에 있어 그 평면이 곧 그 벽이고, 아니면 다발의 평면 `(n_Q·m)·P − (n_P·m)·Q` — `realize::swept_wall`), 원통 쌍에선 근축 평면(`nacre_exact::quad::radical_plane`, `realize::pair_rulings` — 축이 평행해 공유 룰링이 최대 둘)이다. 발은 원통 원점을 지나 축을 가로지르는 평면(`pair_feet_plane`) 위의 점이고 근호 하나다. 순위는 원통 진술만 읽으므로 평면의 순서와 향, 캡의 틀이 이름에 닿지 않는다 — 만드는 쪽(프리즘의 `pierce_def`·`pair_def`)은 꼭짓점을 축 방향으로 그 평면에 내린 발이 몇 번째 룰링의 발인지로 순위를 찾고(`ruling_rank`), 그래서 돌린 몸의 면 위에도 슬롯과 잎이 선다. 순위는 근사로 정하지 않는다 — 한 근이 seam 줄 위에 있으면 근사는 0° 와 359.99…° 를 끝내 가르지 못한다. 세계 진술이 다 있으면 세계에서, 아니면 원통의 진술 틀에서 정확히 정한다. 벽을(다발이면 두 평면을) 그 틀로 당겨 오는 규칙은 `plane_in_cylinder_frame` 하나다: 원통 사슬을 끝에서부터 걸으며 벽도 기록한 운동은 짝짓고, 정확한 운동(사분각 회전·이동·거울)은 되돌리고, 그 밖의 운동은 벽을 제자리에 둬야 한다 — 변환이 평면의 진술을 지키는 것이 바로 그 조건에서다. 실현은 seam 정점처럼 두 길이다(`ruling_point`): 사슬이 접히고 평면이 세계 이름을 가지며 발이 유리수면 발을 축 방향으로 평면까지 옮긴 정확한 점, 아니면 발을 원통 사슬로 재생해 평면과 세계에서 만난다 — 평면이 따로 답하므로 어느 사슬에 앉아도 된다(캡의 법선 둘레로 돈 캡은 진술을 지키고 벽만 노드를 받는다). 원통 쌍 사이의 간선은 공유 룰링, 곧 축 방향 직선이다(`Model::cylinder_pair_relation` 이 정확히 평행을 묻는다). 축이 평행하지 않은 원통 쌍(근호 두 겹)은 아직 만드는 쪽이 없다 — `todo.md` 「원통–원통」.

## 이차곡면 판정 — 원통 우선

이차곡면 교차는 난이도가 갈린다. **평면∩이차곡면**(평면∩원통=타원, 평면∩구=원, 평면∩원뿔=원뿔곡선)은 닫힌 형식이라 쉬운 쪽이고, **이차곡면∩이차곡면**(일반적으로 4차 공간곡선)은 특수 케이스만 닫힌 형식이며 일반은 이미 SSI 에 가깝다. 그래서 이 영역은 「하나의 알고리즘」이 아니라 「어느 곡면쌍까지를 긋고 그 판정을 무엇으로 할지」의 선긋기다. 선은 이렇다: **원통부터**(평면∩원통), 그 뒤 구·원뿔. 일반 이차곡면쌍은 커버리지 밖으로 두고 `Rejected` 로 정직하게 거절한다 — 조용히 틀리는 길이 없으므로 「되는 데까지만 하고 나머지는 거절」이 안전하다.

**판정 방법.** 평면∩원통 교차점은 좌표가 `a+b√c`(a, b, c 유리수) 꼴이다(nacre-exact `QuadVal`).

- 점-대-평면 부호와 같은-판별식 비교는 `a+b√c` 의 부호(`QuadVal::sign`)로 닫힌다.
- 배열에서 관통점은 이름에서 다시 푼 만남(`PiercePoint` — 판정층의 `Judge::pierce` 가 이름에서 풀고, 정확한 만남을 이미 든 자리는 그 쌍으로 짓는다)으로 다니고, 점에 묻는 것(평면의
  쪽, 어느 룰링, 구간 점)은 그 타입의 메서드다. `a + b√c` 비교기가 읽는 정확한 쌍은 `exact` 로 꺼낸다.
- 이름의 평면 하나에 세계 진술이 없으면 만남은 **구간 상자**(`PiercePoint::Bounded`)다 — 같은 규약(방향 `n₁ × n₂`, 근
  `Lo`/`Hi` = `mid ∓ √disc/2A`)의 구간판(`quad::plane_plane_cylinder_hp`)을 그 평면의 저장 향 계수(`Judge::stored_hp` — 이름이
  없으니 정준 향이 곧 저장 향)로 푼다. 구간은 «다르다»만 증명하므로 같음은 이름이 준다: 점은 제 이름의 두 평면 위에 있고,
  원통을 단면하는(`Silence::Sectioned`) 한 이름 없는 평면·그 원통·한 룰링의 두 점은 축 성분 0 인 좌표가 같다(`KeyCmp` —
  비스듬한 평면이면 룰링 쪽은 타원을 둘로 가를 뿐이라 그 규칙은 거짓이다). 접점(`Double`)과 축에 평행한 선은
  구간이 말할 수 없어 그 자리의 이름으로 거절된다. 이 도로는 세계 진술 없는 평면에만 켜진다 — 정확한 쌍이 실패한 이름
  있는 점(넓은 이름, 넘침)은 구간으로 넘어가지 않는다. 규약이 두 도로에서 같다는 것은 debug 빌드가 모든 정확한 관통점을
  구간으로도 풀어 두 상자가 겹치는지로 확인하고, 정확한 점의 쪽과 어느 룰링(`PiercePoint::side`·`ruling_side`)·원 위 세 점의
  구간 회전(`between_on_circle_hp`)·한 선 위 두 관통점의 구간 순서(`Judge::line_order`)·평면 둘의 선이 원통과 만나는
  모양(`Judge::line_meet` 의 정확 가지)도 구간 몸통이 갈리는 곳에서 정확한 답과 맞춰 본다. 판정층 시험은 평면의 향을 모든
  방향으로 돌린 표에서 각 몸통의 부호를 참 기하에 대고 잠근다.
  원통 판정의 구간 팔은 판정층(`nacre_judge::cylinder`)의 문이다 — 점의 상자, 쪽, 어느 룰링, 세 순서(선 위·축 위·세계 축),
  원 위 세 점의 회전, 선이 축을 따라 어느 쪽으로 가는가(`Judge::stored_line_along` — 저장 향의 외적과 축의 곱, 유리수
  방향이 없는 단면 클래스의 선과 룰링이 한 점에서 만날 때), 선이 원통을 지나는가, 축을 따르는 선이 호들 밖인가(방향만
  읽어 선이 옆면 위든 곁이든 같다 — 선분의 한 점으로도, 두 평면의 선으로도 묻고, 뒤의 것은 무한면 밖이어도 밖이다), 게이트의 물음(축 거리 대 반지름); 정준 방향의 향
  («이름 없는 평면의 정준 방향은 저장 방향»)도 한 문(`Judge::canon_sense`)이 말한다. 평면 둘의 선이 원통과 어떻게 만나나는 한
  문(`Judge::line_meet`)이 두 도로를 한 어휘(`LineMeet`)로 답해, 호의 교차·룰링 분할의 교차·룰링 위의 교점·원판 위의 현이 그
  답 하나를 읽는다 — 정확 도로만 말할 수 있는 모양(평행·일치·옆면 위·접점)은 구간 도로에서 «말하지 못함»이다(정확 전용
  `plane_plane_cylinder` 호출은 그 밖에도 남는다). nacre-ops 에 남는 것: 정확 팔(선 위 순서·세계 축·축 좌표의 `a + b√c`
  비교), 점의 계단으로 거절 이름을 고르는 자리(`crossing_on_ruling`), 링 상자 탐침(`ring_box_points_hp`·`plane_sign_hp`·
  `outside_ring_box`, `nesting` 의 말 없는 클래스 탐침), 정확한 점만 받는 구성과 판정(관통점의 유리수 좌표·켤레 중점·섞인 링·
  룰링 접기·별칭 기록·호의 접선 방향). 두 도로가 한 몸통이 될 수 없는 물음은 «한 점의 두 이름»이다 — 원 위 순서
  (`theta_between`)에서 정확 도로만 그것을 말한다(정밀도를 올려도 구간은 0 을 증명하지 못한다 — todo 「원통 판정의 상승 —
  인구가 생기면」). 실현(`realize.rs` 의 고정밀 재생)은 판정이 아니라 이 목록 밖이다.
- 원통 면 위 사건들의 원형 순서는 서로 다른 절단면(서로 다른 판별식 u, v)에서 온 두 점의 비교라 `ℚ(√u,√v)` 의 4항 원소 부호가 필요하다. 이것은 재귀 탑으로 닫힌형이다(`P+√v·Q` 분해 → 상반 부호면 `sign(P)·sign(P²−vQ²)` — nacre-exact `quad.rs` 의 `biquad_sign`).
- 일반 대수적수 기계도, QI 의 pencil 분류도 원통에는 필요 없다.
- **차트 t 의 수치는 실현(캐시)·문서화 전용이다.** t 값 비교는 노름 라디칼(`|d|`, `|e₁|`)까지 끌고 와 3-라디칼이 된다. 원형 순서는 **(seam 반평면 부호 `w·e₂`, 같은 반평면 안 외적 부호 `(w₁×w₂)·d`)** 두 술어로 판정한다(e₁ = `ref_dir` 의 축-수직 성분, e₂ = d×e₁ — 둘 다 유리수라 위 두 부호 원시로 닫힌다). 차트의 배제점 = seam 이고, seam 모선 위의 점(θ=0)은 **이름**(`SeamIncident`)으로 답한다 — 순환 순서의 순위는 `SeamOrder::seam_first` 한 곳이 준다(seam 점이 첫째: 순환 순서는 어디서든 한 번 끊을 수 있다). 호의 범위(`arc_span`)와 노드의 원형 정렬(`circular_order`)이 같은 규칙을 읽는다.

### 참고 논문과 라이선스 규율

일반 이차곡면쌍과 곡면 위 배열(원통 차트의 곡선 어휘)에 갈 때 꺼낼 참고처다.

- **QI 3부작 + 구현 논문**(pencil 분류로 교차 타입을 대수적으로 exact 결정). Dupont·Lazard·Lazard·Petitjean, "Near-optimal parameterization of the intersection of quadrics" — SoCG 2003, JSC 2008 저널 3부작(Part I 생성 알고리즘 43(3):168–191, Part II pencil 분류 43(3):192–215, Part III 특이 교차 43(3):216–232). C++ 구현 논문: Lazard·Peñaranda·Petitjean, "Intersecting Quadrics: An Efficient and Exact Implementation", SoCG 2004 / Comp. Geom. 35(1–2):74–99, 2006.
- **CGAL 3D Spherical Kernel**(2차 대수적수 exact — 다른 계보). de Castro·Cazals·Loriot·Teillaud, Comp. Geom. 42(6–7):536–550, 2009. 간접 술어(Attene, 정의 기반)와 다른 계보라 「확장」이 아니라 「대안 접근」으로만 참고한다.
- **곡면 위 배열.** "Arrangements on Parametric Surfaces" I(Berberich·Fogel·Halperin·Mehlhorn·Wein)·II(Berberich·Fogel·Halperin·Kerber·Setter), Math. Comput. Sci. 4, 2010 — 엔진은 곡선·곡면 종류를 모르고, 곡선 종류마다 술어 묶음(비교·위아래·분할·교점과 중복도)과 곡면 위상(원통에는 무한 면이 없다, seam 은 식별선)만 댄다.
- **이차곡면 교선의 정확 배열.** Berberich·Hemmer·Kettner·Schömer·Wolpert, SoCG 2005 — 교선을 평면에 사상해 배열하고, 퇴화는 경우 코드가 아니라 정확한 0 의 판정이 받는다.
- **원뿔곡선 호의 정확 배열.** Berberich·Eigenwillig·Hemmer·Hert·Mehlhorn·Schömer, "A computational basis for conic arcs and Boolean operations on conic polygons", ESA 2002 — 호는 지지 곡선과 대수적 끝점이고, 끝점은 두 곡선의 근으로 지어 정체가 출처에서 온다(관통 이름과 같은 생각); 원은 같은 기계의 특수형이다. 단면 호를 원통 각 θ 로 순서 매기는 것(축 방향 사영)은 이 출처들 밖의 nacre 유도다.
- **평면 × 원통 단면의 매개화.** Wang·Goldman·Tu, "Enhancing Levin's method for computing quadric-surface intersections", CAGD 20, 2003 — 원통을 `θ` 와 높이로 매개화하면 평면 단면의 높이가 `cos θ, sin θ` 의 일차식이다.

라이선스 규율은 indirect predicates 에 적용한 것과 같고 더 엄격하다.

- **QI 구현**(LORIA/INRIA, gamble.loria.fr/qi): "free for non-commercial use" — 비상업 한정이고 내부 부품은 또 다른 라이선스다. **소스 열람·차용 금지**, 논문만 참고.
- **CGAL Spherical Kernel**(`Circular_kernel_3`): 패키지 오버뷰에 License: **GPL**. 링크하면 nacre 전체가 GPL 에 오염된다. **소스 열람·링크·차용 절대 금지**, 논문만 참고. (CGAL kernel 기반부는 LGPL 이나 Spherical Kernel 은 GPL.)
- **요약:** 논문만 읽고 clean-room 구현, 소스는 열람조차 안 함, 완성 후 실행 대조(dev 전용)만 허용. exact 산술 기반은 `geometry-predicates`(MIT/Apache) 등 자유 라이선스 부품으로.

### 곡면 내/외 판정 후보

「찾은 교차를 어떻게 분류하나」의 후보이며 SSI 해법이 아니다. 1순위 **exact ray casting**(우리 계보와 정합, Cherchi 2022 검증, 메시 경유라 곡면 직접 판정을 우회). 2순위 **GWN**(generalized winding number, Jacobson 2013 — 메시/point cloud in-out 은 성숙 기법이고 libigl·Axom[BSD] 구현이 있다. watertight 무관 강건성이 강점이나 always-closed 에선 덜 필요하고, 느리며 경계 round-off 가 있다. trimmed NURBS 정확 GWN 확장[Spainhour 2024~26]은 검증 진행 중 — 「메시 경유 없이 곡면에서 직접 판정하고 싶어질 때」의 대안). 참고: **graph cuts**(Diazzi/Attene 2021 — 일부 모호한 자기교차 케이스에 우수), **EMBER winding number vector**(Trettner 2022). 채택 전 라이선스를 확인한다 — Cherchi/Attene 계열은 LGPL 이라 위와 같은 규율(논문·MIT/Apache 소스만 참고, LGPL 소스 열람 금지, 실행 대조만 허용).

## `nacre-predicates` 는 순수 수치층이다

간접 술어의 implicit point 는 「어느 원시 요소들의 교차인지」를 정의로 보유한다(점은 자기 정의다). 그 정의가 `Handle<Surface>` 를 담는데, `nacre-predicates` 는 위상 계층보다 **아래** 층이라 그것을 참조하면 의존이 순환한다 — `Store`/`Handle` 을 최하위 `nacre-store` 로 내려 푼 것과 같은 구조의 문제다. 그래서:

- `nacre-predicates` 는 평면 계수·좌표를 **평범한 배열**로만 받는다(커널 타입 무의존, standalone 분리 가능).
- Handle 기반 정의(`Vertex::ThreePlane([Handle<Surface>; 3])`)는 위상 계층(`nacre-topo`)에 살고, 부호 판정 시 위층이 Handle → 계수를 뽑아 술어에 넘긴다.
- `nacre-judge` 도 같은 규율이다: `nacre-exact`·`nacre-predicates`·`nacre-math` 에만 의존하고, 평면 표는 `Witness`/`PlaneWitness` **트레이트**로 받는다(구현은 `nacre-ops` 의 행 타입).

## 면 위 작업(pad/pocket)은 부호 있는 extrude + 불리언

pad·pocket 은 커널 연산이 아니라 편의 레이어(kit)의 조합이다 — 상용 CAD 와 동형인 **tool body + boolean**: 면의 스케치 프레임(`face_sketch_frame`) → `Extrude`(pad 는 `+d`, pocket 은 `−d` — 같은 프레임에서 ŵ 반대로) → `Boolean`(면의 솔리드가 첫 피연산자, pad 는 `Fuse`, pocket 은 `Cut`). 새 정확 술어는 없다; 판정은 전부 불리언의 것이다.

- **프레임은 한 곳에서 정해진다.** `face_sketch_frame` 이 정한 프레임이 곧 `face_plane` 이 호출자에게 약속하는 프레임이고, `Extrude` 는 그 프레임에서 고리를 짓는다(`frame_rings` — 게이트 `exact_frame` 이 세계 기저가 유리수라 하면 세계 도로, 아니면 `push_frame_node`).
- **밑면은 대상 면의 `Surface` 핸들을 공유한다.** `Extrude` 의 밑캡은 늘 프레임의 평면 핸들이라 새 평면을 push 하지 않고 아무것도 진술하지 않는다 — 어느 쪽으로 쓸든 near cap 은 면과 flush 다.
- **공면은 진실로 안다.** 불리언의 평면 클래스 병합(`plane_classes`)은 **먼저 `Surface` 핸들로 묶고**(interning 이 한 평면을 한 핸들로 만들었다 — O(1), 정확, 회전 무관) 서로 다른 핸들의 대표끼리만 `Judge::planes_coplanar` 에 묻는다. 판정은 두 세계 이름이 다 있으면 정준형의 동일성으로 끝나고(같으면 병합, 다르면 증명된 비병합), 없으면(사분각 밖 회전·프레임) 정의 점의 도로 — 합성 회전, 상승 — 로 간다. 캐시(평면 캐시의 계수·면 꼭짓점 좌표)는 어느 갈래에도 없다: 서로 다른 두 평면의 반올림된 상이 비트 단위로 같을 수 있고, 그 병합은 한 몸의 벽을 다른 몸의 평면 위에 놓는다. 증인이 평면을 펴지 못하는 행(이름 없는 `Known` 의 공선 점)은 표(`collect_planes`)가 `DegenerateFace` 로 거절한다. 「참조로 아는 공면은 공짜, 우연한 공면만 계산」 — 상용 커널이 coincident 면을 imprint 로 공유 토폴로지로 승격하는 것의 nacre 판이다.
- **결과 면은 클래스 대표의 핸들을 든다.** 같은 평면이 두 핸들로 진술돼 있으면(기운 면의 프레임 노드 진술) 공구의 먼 캡은 결과에서 몸 쪽 핸들 위에 남는다 — 「내 캡이 어디로 갔나」를 핸들로 묻는 쪽은 그 경우를 놓친다(테스트 픽스처 `fixtures::Feature::cap_face` 는 그래서 평면과 위치로 고른다).
- **프로파일은 면 밖으로 나가도 된다(오버행).** 「프로파일이 면 안에 있다」는 제약이 없다. 면 안에 담긴 footprint 도, 면을 넘는 footprint(boss cantilever, 모서리 슬롯)도 같은 불리언이 답한다. 몸보다 깊은 pocket 은 관통 컷이고, 몸을 가르는 pocket 은 여러 솔리드다.
- **kit 이 결과를 보고 판단하는 전제는 하나다**(깊이 양수·섬 하나는 입력 검사): pad 의 `Fuse` 가 여러 솔리드로 돌아오면(두 단일 셸 솔리드는 닿지 않았을 때만 갈라진다) footprint 가 면을 빗나간 것이고, kit 이 그것을 거절한다 — 솔리드 개수라 정확 술어가 아니다.

**exactness 는 후퇴하지 않는다.** 회전 없는 불리언의 출력 정점은 전부 평면 삼중항으로 이름 붙는다(`Vertex::ThreePlane`; 잠금 = `an_unrotated_boolean_names_every_vertex_by_its_plane_triple`). 상자 모서리는 실제로 세 평면의 교점이므로 그 정의는 참이고, 축정렬 사례에서 그 정점들의 tol 은 0 이라 `EPS_CONSTRUCTED`(1e-9)보다 엄격하다. 그 귀결: **provenance 로 면·정점을 고르는 코드는 성립하지 않는다** — 「컷에서 온 정점」이라는 구별은 출력에 없고, 모든 면에 대해 같은 답이 나온다.

**업계 정합.** 피처 레이어는 pad = tool body + boolean 으로 통합돼 있고(SolidWorks/NX/Creo/Fusion), 커널 레이어는 그래도 coincidence 전담 로직을 갖는다. 상용은 그것을 imprint + tolerance 로 두고, nacre 는 **커널 불리언 안**에 둔다(`Surface` 핸들 공유 + exact 술어) — tolerance 대신 exact 술어로 그 자리를 채우는 소수파다.

## 다중 솔리드와 non-manifold 정책

기준은 OCCT 의 거동이다: 엣지로만 만나는 두 박스의 `Fuse`, 막대를 가르는 `Cut` — OCCT 는 둘 다 **COMPOUND of 2 SOLIDs**(checkshape valid, STEP `2× MANIFOLD_SOLID_BREP`, non-manifold 엔티티 0)로 낸다. 엣지 접촉 Fuse 는 접촉 엣지 + 2 정점을 **공유**하되(V16→14, E24→23) 각 솔리드는 깨끗한 manifold 박스다.

- **불리언은 `Vec<Handle<Solid>>` 를 낸다.** 「Cut 두 동강」「엣지 접촉 Fuse」는 전부 정상 결과다. 양수 성분 전부를 각각 솔리드로 방출하고, `Model` 의 `live_solids` 가 그것들을 든다.
- **manifold 는 솔리드별 판정이다 — 두 솔리드가 경계 요소를 공유해도 각각은 manifold 다.** 공유 엣지가 compound 수준에서 4면에 닿는 것은 각 솔리드가 2면씩이라 위반이 아니다. 각 솔리드 안의 manifold 전제(엣지당 면 2·Euler)는 그대로다.
- **non-manifold 솔리드는 지원하지 않는다 — 불필요하다.** OCCT 도 만들지 않는다(엔티티 이름부터 MANIFOLD).
- **솔리드 간 경계 요소 공유는 허용되고 필요하다.** 솔리드들이 정점·엣지·면 Handle 을 공유할 수 있다. 그 공유는 **불리언이 명시적으로 만든다**(전역 자동병합이 아니다) — 아래 「명시 공유」와 같은 메커니즘이다.

## 회전 지원의 확정 결정

축정렬 전용의 한계를 넘어 회전(비정렬·마름모·각도 스케치)을 지원하는 설계의 결정들이다.

1. **유리수 치수·각도 표현**(`nacre-exact`). 유리수는 입력과 유리수-순수 파생에만 살고, 무리수·복잡 연산·비트 상한이 닿는 순간 캐시는 f64/고정밀 실현으로 내려간다. 정의는 불변이라 필요할 때 재계산으로 정확히 복원된다.
2. **명시 공유 — 전역 자동병합은 없다**(TNP 와 충돌한다). 면 위 스케치 = `Surface` Handle 재사용, 그리고 불리언이 계산한 접촉만 공유한다. 우연히 같은 좌표는 별개로 남고 불리언 시점에 판정된다. 「같다」고 판정되면 Handle 재사용으로 추이성 붕괴를 구조적으로 막는다.
3. **판정 정책 — 묻지 않고, 근거를 붙여 보고한다.** 자세한 것은 CIP 절의 「판정 결과와 보고」.
4. **계획: op-log 소유.** 상위 `Document` 가 `Vec<Operation>`(진실) + 파생 `Model` 을 소유하고, 피처 트리는 편의 레이어의 단계(pad·pocket 같은 조합)를 기록한다. TNP 는 「치수 변경 = 자동 replay, 위상 변화 = 수동 재지정」. (`Document` 타입은 코드에 없다.)
5. **다중 솔리드 출력** — 위 절.

## CIP — 허용오차 부호 술어

`nacre-judge`(Certified Indirect Predicates)는 회전이 들어왔을 때의 부호 판정 층이다. 층의 지도는 `lib.rs`·`kernel/mod.rs`·
`kernel/frame3/mod.rs` 의 doc 이고, 구간 산술(`Bounded`·`HpBounded`)과 각도의 실현은 `nacre-exact`(`bounded.rs`·`angle.rs`·
`realize.rs`)에 있다.

### 왜 필요한가

`nacre-judge` 의 `kernel`(toleranced 부호)·`predicate`(평면 배열 술어)·`cylinder`(관통점) — `lib.rs` doc; `nacre-predicates` 는 이름 정수 도로의 f64 필터와
공유 모션 `orient3d` 를 댄다.

축정렬 판정은 좌표가 유리수라 exact 다. **회전이 들어오면 좌표가 무리수가 된다** — 유리수 도수의 cos·sin 은 `{0, ±½, ±1}`
밖에서 무리수이고(Niven), 임의 각도면 초월수다. 회전 좌표는 근사일 수밖에 없으므로, 모델링 변환·각도 스케치가 그런 좌표를
만들 때 **판정만은 조용히 틀리지 않게** 지키는 것이 CIP 다.

**핵심 명제.** indirect predicates(Attene 2020)는 선형 요소(선·평면)의 교차, 즉 다항식에만 성립하고 회전은 다항식이 아니다
(수학적 한계). 할 수 있는 일은 하나다: **회전 좌표는 근사하고, 그 위의 판정에 입력 근사(tol)까지 반영해 부호를 인증하거나
정직하게 기권한다.** 세 계보의 하이브리드다 — Attene 2020 간접 술어(implicit point 구조) · Shewchuk 1997/CGAL 필터(f64 필터 →
고정밀 상승) · Guibas 1989 epsilon-geometry(값 + tol, sound 판정 아니면 기권). **포기한 것은 초월수까지의 이론적 완비성이지
건전성이 아니다** — 틀린 부호는 결코 내지 않는다.

**대체가 아니라 흡수다.** 점을 「정의 + 누적 tol」로 들고, 어느 도로로 갈지는 tol 이 아니라 이름과 모션이 가른다 — 평면의
이름이 정확히 답하는 것은 이름 정수 도로가, 나머지는 필터 + 고정밀 상승이 답한다(`predicate` 모듈 doc). 술어의 계산 구조(3-평면 정의, 행렬식
부호, 좌표 무독)는 그대로 두고 「평면 계수가 정확하다」만 「계수에 tol 이 있다」로 넓힌다 — 회전된 평면은 계수에 tol 이
붙은 평면일 뿐이다.

### 점 = 정의 + 방향별 tol — `WitnessPoint`

`kernel/frame3/witness.rs`(`WitnessPoint`·`at`·`coord_iv`·`coord_approx`·`tol`·`coord_hp`), `kernel/frame3/mod.rs`(`HpCell`),
`nacre_exact::Bounded::of_rat`; 판정 표가 증인을 받는 문은 `nacre-ops::rotated_vertex::plane_witness`. 타입 그림은
「판정 (연산 동안만 산다)」.

- `base` + `chain` 이 정확한 **정의**이고 잃지 않는다. `coord_iv` 는 축마다 f64 값과 그 증명된 오차를 한 원자(`Bounded`)에
  묶은 캐시다 — 회전을 f64 산술로 나르므로 값은 최근접 f64 가 아니고, 참값까지의 거리는 `tol()` 이 든다. 모션 없는 유리수
  점과 모션을 나른 점이 같은 틀(`base`+`chain`, 축마다의 `tol()`)이라 특수 경로가 없다. 같은 정의는 `coord_hp` 에서 같게 실현되고(**경로 무관** — 건전성 논증의
  뿌리) 동등성은 정의만 본다 — `WitnessPoint` 와 그 `PartialEq` 의 doc. `HpRc`·`HpCell` 의 선택은 `HpCell` doc.
- **tol 은 축마다다.** 스칼라(구) tol 은 정보를 버린다 — z축 둘레 회전은 z 오차가 0 인데 구는 모든 방향에 실어
  「z = 5 평면 위인가?」를 애매하게 만들고 상승시킨다. 행렬식은 방향마다 다른 계수로 증폭하므로(2D orient 에서 `a` 의 x-tol 은 `(by−cy)`
  와, y-tol 은 `(bx−cx)` 와 곱해진다) 하나로 뭉치면 증폭을 따라갈 수 없다. 모션의 단계와 술어의 필터(`det3_iv`)가 같은 구간
  산술을 걷는다(「tol 의 전파」).
- **생성자는 `at(base)` 하나다**(`hp` 가 비공개라 크레이트 밖에 다른 길이 없다). 좌표마다 최근접 f64 와 그 반올림 한계를
  계약에서 든다 — 재지 않는 까닭과 한계의 느슨함이 상승으로만 가는 것은 `at`·`Bounded::of_rat` 의 doc, 고정 정밀도의 측정이
  tol 0 으로 읽는 유리수는 `a_rational_nearer_its_f64_than_a_measurement_resolves_is_still_charged`. 유리수에서 온 점이
  술어에 닿는 길은 이것뿐이다 — 반올림된 f64 를 `Rat` 로 다시 들어올리면 다른 점을 tol 0 으로 이름 붙인다(«f64 로 정확하다»를
  전제로 받는 생성자는 「가지 말 것」).

**판정 표의 계약은 «평면 위 정확한 세 점»이다.** 이름 없는 평면은 정점이 모두 자기 프레임에서 좁은 만남으로 풀리면 그 점들
(각자의 사슬 + 평면의 모션)을, 아니면 판정 프레임의 probe 셋을 증인으로 든다(`plane_witness` doc). 그래서 `Through` 평면
전용 판정 기계는 없다 — 구간-계수 경로 없이 기존 세 술어가 답한다(`plane_hp_through` doc 의 끝 문단).

### 모션 사슬 — `MoveNode`

`kernel/frame3/mod.rs`(`MoveNode`), `frames.rs`(`WideFrame`·`FrameThrough`·`JudgedPoint`·`judged_basis`·`anchor_approx`),
`indirect.rs`(`plane_hp_through`·`cramer_hp`·`meet_hp`), `det3.rs`(`chain_parity`·`shared_base`), `chain.rs`(`fold_suffix`).
topo 의 `Motion` 노드를 펼치는 문은 `nacre-ops::rotated_vertex::motion_chain`, 노드를 남길지는 `nacre-ops::transform`
(`carry_of`·`Carry`·`chain_motion`). 변종마다의 실현과 약속은 그 doc 이다.

모션은 교환되지 않으므로 **사슬의 순서가 곧 정의다.** 사슬은 tol 계산만이 아니라 정의의 합성이다 — 이동·반사도 정의에 들어야
`x = 1` 을 `7/11` 옮긴 벽과 `x = 0` 을 `18/11` 옮긴 벽이 한 평면임을 말할 수 있고(f64 로는 1 ULP 갈려 한 부품이 두 몸통이 된다 —
`tests/coverage/placement.rs`), `R` 다음 `T` 를 적을 수 있다. 교차는 담지 않는다(`Vertex` 의 몫).

- **옮겨 적기 — 진술이 정확히 옮겨지면 노드를 만들지 않는다.** 솔리드의 모든 면의 진술(점 셋·원통 정의)이 모션 전체로
  `i128` 안에서 옮겨지면 진술에 옮겨 적고, 하나라도 못 옮기면(`Through` 평면, 넘침) 솔리드 전체가 기록한다 — 면마다 가르면
  그 사이 모서리가 혼합 프레임 정점이 된다. 이력이 있는 면은 언제나 기록한다. 판단은 진실만 읽는다(캐시로 고르는 길은
  「가지 말 것」). 두 상태(`Carry::{Full, None}`)뿐이라 한 면에서 모션이 진술과 노드로 갈리는 일은 없다 — `carry_of`·`Carry`·
  `chain_motion` 의 doc.
- **`Mirror` 만 improper 다**(`det = −1`; 프레임 노드 셋은 `v̂ = ŵ × û` 라 proper). 공유 모션을 행렬식에서 상쇄하는 지름길은
  모션 전 데이터를 먼저 같은 손방향으로 가져온다 — 홀수 반사 사슬은 모든 입력의 행렬식을 일률적으로 뒤집어, 보정이 없으면
  보수적 miss 가 아니라 확신에 찬 오답을 낸다. 점은 `shared_base` 가 반사하고, 평면 계수는 반사한 점에서 유도한 평면의
  부호를 따라야 한다 — 외적이 pseudovector 라 법선만 반사하면 전역 부호가 갈린다(`PlaneWitness::name_ints`·`name_stored_ints`
  의 doc).
- **프레임 노드 셋**은 평면 자신의 프레임으로의 기저 변환이다. `Frame` 은 입력이 전부 정확한 유리수이고 유일한 무리수
  걸음이 길이로 나누기라 판정할 수 있다. `FrameWide` 는 `Rat` 에 안 드는 평면(`Wide` 이름, 또는 제곱 길이가 `i128` 을 넘는
  narrow)의 같은 것이다. `FrameThrough` 는 **이름이 아예 없는** 평면(혼합 프레임 datum)의 것으로, 정의점 셋
  (`JudgedPoint::{Pure, Meet}` — `Meet` 는 좌표가 어디에도 없는 세 담체의 교점)을 싣고 정준 기저를 구간으로 유도한다.
  `v̂` 를 실현하는 두 길, `vertical` 을 고정 rung 에서 한 번 정하는 까닭(`flip` 은 호출자가 넘긴다), `FrameThrough::of` 의 `None` 이 퇴화의 증명이
  아닌 것은 `MoveNode`·`PlaneFrame`·`FrameThrough`·`FrameThrough::of` 의 doc. 동차점을 평면으로 잇는 사영 join 이 아무것도
  나누지 않고 배율의 부호를 값 안에서 없애는 것은 `plane_hp_through` doc, 판정 프레임의 meet 를 나누는 한 자리는
  `anchor_approx` doc.

**실현은 나눠도 된다. 금지되는 것은 술어의 부호 질문에서의 나눗셈이다** — 거기서 몫은 술어가 믿어야 하는 무리수를 제조한다.
실현의 나눗셈은 자기가 만든 반경을 돌려주고, 제수가 0 에 닿을 수 있으면 `None` 이다(호출자는 상승하거나 이름으로 거절한다):
`HpBounded::div`(구간 제수)·`inv_sqrt`(구간 입력의 `1/√x`)·`div_exact`(정확한 제수)의 doc.

### tol 의 전파

`kernel/frame3/witness.rs`(`rotated`·`mirrored`·`translated`·`framed`·`framed_wide`·`framed_through`, `sum`·`diff`·`prod`·`in_frame`),
`kernel/frame3/mod.rs`(`narrow_hp`), `indirect.rs`(`plane_iv`), `nacre_exact::{Bounded, Angle::cos_sin_error_of, inv_sqrt_error_of}`;
사슬 접기는 `nacre-topo` 의 `Model::chain_point_rat`. 건전성 시험은 `kernel/tests/frame3/witness_points_and_frames.rs`.

- **오차 한계는 최악 상한이다**(삼각부등식 — 실제보다 작아지지 않는다). 확률 전파는 보장이 아니다(「가지 말 것」). 커짐은
  다른 길로 다룬다 — 유리수 구간은 tol 0, 애매하면 필요한 만큼 정밀도를 올린다.
- **단계는 `Bounded` 로 걷는다.** 입력이 전부 구간이고(좌표, 쓴 `cos`/`sin` 쌍에서 잰 오차, 유리수의 `of_rat`, 프레임의 `1/√`
  의 잰 오차나 고정밀 실현을 좁힌 구간) 덧셈·뺄셈·곱셈이 `Bounded` 의 것이라 반지름이 연산마다 증명되고 손으로 적는 항이
  없다. 기존 tol 의 방향 회전(`|R|·기존 tol`)과 좌표 혼합(`u` 가 한 축에선 `cos` 의, 다른 축에선 `sin` 의 오차를 만난다)은
  그 산술에 저절로 든다 — `rotate_about`·`rotated` 의 doc, 회전이 기존 tol 을 나르는 것은 시험
  `seeded_tol_transports_through_rotation`. 접선형 tol 이 건전하지 않은 것은
  `the_tangential_tol_is_unsound_and_the_mixing_one_is_not`(「가지 말 것」). 각도는 정확한 유리수 도수라 각도 자체의
  불확실성 항은 없다 — 유도 치수(구속 솔버·비유리수 각)가 들어오면 그 불확실성은 `cos`/`sin` 의 구간으로 들어가야 한다.
- **정확한 연산은 반지름을 그대로 넘긴다** — 정확한 0 과의 합·차, 정확한 0·2의 거듭제곱과의 곱. 그래서 원점 피벗의 사분각
  회전·원점을 지나는 거울·원점 축정렬 프레임은 tol 0 이고(`quadrantal_origin_chain_is_tol_zero`·
  `an_axis_aligned_frame_through_the_origin_is_tol_zero`) 그런 모델의 커널 필터는 입력이 정확한 만큼 좁다 — 어느 길로
  갈지는 tol 이 아니라 이름과 모션이 고른다(아래 «길은 tol 이 아니라 이름과 모션이 고른다» 글머리). 규칙이 `Bounded` 안이
  아니라 `witness.rs` 의 `sum`·`diff`·`prod` 에 사는 까닭은 그 앞 주석과 「가지 말 것」.
- **회전 기여는 적용점에 의존한다** — 회전은 방향을 정하므로 그 뒤의 이동이 기여를 키운다. 그래서 점은 tol 값만이 아니라
  자기 사슬(조상 모션)을 들고 다닌다.
- **유리수 구간은 접는다.** 기준은 「직선/회전」이 아니라 「유리수/무리수」다. 이동·사분각 회전·축 거울의 사슬은 노드가
  태어날 때 한 유리수 사상으로 접히고(「캐시」의 사슬 접기), 접히는 사슬 아래의 점은 실현·모서리 판독이 정확한 유리수로
  읽는다(`realize.rs`). 증인점의 사슬(`motion_chain`)은 노드를 하나씩 걸으므로 그 도로에서는 연속 이동도 걸음마다 반올림한다.
  같은 축·같은 피벗의 연속 유리수-각도 회전을 각도 합으로 한 번 실현하는 묶기는 짓지 않는다
  (todo 「같은 축 회전의 각도 합성」). 묶지 않은 증분 실현은 걸음마다 행합 `|cos|+|sin| ≥ 1` 을 곱해 한계가 지수로
  자란다(2D 실험에서 30걸음에 실제 오차의 약 20만 배; 실제 오차는 `R` 이 노름을 보존해 평평하다;
  `bundling_is_tighter_than_incremental` 은 방향만 잠근다). 자라도 sound 한 최악 보장이라 상승·거절로 갈 뿐 조용히
  틀리지 않는다.
- **평면 계수의 tol 은 점 tol 의 따름정리다** — 계수는 세 점의 뺄셈·외적이다(`plane_iv` doc, `plane_coefficient_tol_soundness`).
  판정 경로에 sqrt 섭동은 없다 — 판정은 정규화하지 않은 원시 계수를 읽고(`plane_iv`), 단위 법선은 평면 캐시의 것이다
  (`nacre-topo` `surface_cache.rs` 의 `world_unit_normal_f64`).

### 고정밀 층 — `astro-float`

`nacre_exact::{HpBounded, Mag, Angle::cos_sin_hp, Angle::cos_sin_f64, round_to_f64}`, `realize.rs`(`ladder_f64`), `lib.rs` 의
메모(`TRIG`·`F64_TRIG`·`F64_ERR`).

고정밀 층은 **`astro-float`**(순수 Rust 임의정밀)다. 정밀도가 다이얼이라 고정 폭의 천장이 없고,
「더 넓은 고정 폭으로 상승하거나 정직하게 거절하거나」의 딜레마가 「필요한 정밀도로 올린다」 하나가 된다. double-double 을 쓰지 않는 까닭은
「가지 말 것」.

**이 층은 상승 경로 전용이 아니다.** f64 캐시도 이 층의 실현을 정확 반올림한 사본이다 — `Angle::cos_sin_f64` 는 libm 을
부르지 않고 128비트에서 시작하는 사다리(`ladder_f64`)로 최근접 f64 를 정한다. 정확 반올림된 값은 유일하므로 플랫폼·프로파일·
호출처에 무관하게 같은 비트다. `f64::cos` 가 각도의 함수조차 아닌 측정, 90°-계열 분기가 종료 조건인 까닭, 각도당 메모의 값은
`cos_sin_f64`·`F64_TRIG` 의 doc(libm 은 「가지 말 것」).

이 문서의 「고정밀 상승」은 전부 이 층을 가리킨다.

### 판정 = 필터 + 계산된 상승 + 실현 캐시

`nacre-judge` 의 `predicate`(`Judge::orient3d`·`Judge::cmp_coord`·`Judge::plane_pair_dir_sign`·`Judge::planes_coplanar`·
`Judge::name_rows`)와 `kernel::frame3`(`orient3d_judge`·`indirect_orient3d_judge`·`indirect_cmp_coord_judge`·
`dir_sign_judge`·`escalate`·`WitnessPoint::coord_hp`), 이름 정수 도로의 필터 `nacre-predicates`(`indirect_plane_side_filter`·
`indirect_cmp_coord_filter`·`det3_sign_filter`)와 정수 쌍둥이 `nacre-exact`(`int_plane_side`·`int_cmp_coord`·`int_dir_sign`),
표를 짓는 `nacre-ops`(`planes::collect_planes`·`rotated_vertex::plane_witness`·`tolerant`). 타입은 「판정 (연산 동안만 산다)」.

**길은 tol 이 아니라 이름과 모션이 고른다.** 술어마다 자기가 닿는 평면만 보고(`predicate` 모듈 doc): ① 산술 없이 답이
정해지는 질문은 조합 항등식이 0 으로 답하고, ② 몇몇 술어는 제 정확한 길을 따로 둔다 — `cmp_coord` 는 이름 도로보다
먼저 `cancel_cmp_coord`(여섯 평면이 한 축 둘레의 한 모션을 공유하고 그 축의 비교가 이름으로 정해질 때 — 조건은 그 doc),
`planes_coplanar`·`planes_parallel` 은 세계 이름 비교와 합성 회전(`coplanar_by_composed_rotation`·
`parallel_by_composed_rotation`) 등, 각 술어의 doc, ③ `name_rows` 를 묻는 술어는 모든 평면이 이름을 갖고 `is_rotated` 가 모두
거짓이거나 전부 한 사슬이면(`cmp_coord` 는 앞의 것만 — 좌표 비교는 모션 불변이 아니다) **이름 정수 도로**, ④ 그 밖은
커널 — 증인 점의 구간 필터, (직접 `orient3d` 에서만) 공유 모션 상쇄, astro-float 상승. `planes_parallel` 은 커널에 가지 않고,
증명하지 못하면 `None` 이다. tol 0 인 점도 이름 없는 평면이나 섞인 모션에 걸리면 커널로 간다(`Witness::is_rotated` doc).

**이름 정수 도로.** 이름 있는 평면의 질문을 그 이름의 정수로 정확히 답한다 — 폭과 무관하게 한 길이다: 정수의 최근접
f64 로 `nacre-predicates` 의 필터를 돌리고(입력이 반올림이어도 답하면 옳다 — 여유는 그 필터 doc 의 유도 표), 못 가르면
`nacre-exact` 의 BigInt 쌍둥이가 답한다. 같은 좌표의 두 점(어떤 필터도 못 넘는 정확한 0)은 대부분 이 도로에 오지
않는다 — 조합 항등식이 먼저 답한다(lib 스위트(순차)에서 이 도로에 온 `cmp_coord` 질문의 0 은 2.52M 중 10,514, 항등식을 끄면
4.98M 중 2.47M — `Judge::cmp_coord` 본문). 게이트는 `Judge::name_rows`, 이름을 저장 방향(σ × `frame_sign`)과 정준
손잡이(홀수 반사 사슬의 `parity·C`)로 접는 것은 구성 때 한 번 `name_stored_ints`, 한 사슬을 든 질문이 모션 전 이름으로
정확히 답하는 까닭은 `PlaneWitness::chain_id`·`PlaneWitness::name_ints` 의 doc(이름이 없을 때 오르는 측정 포함).

**필터의 모양은 입력이 정한다.** 이름 정수 도로의 입력은 정수의 최근접 f64 라 술어마다 상수 오차 한계(`32·U` 류)를
유도해 둔다. 커널의 입력은 회전 tol 을 든 증인 점이고 그 tol 은 경로마다 편차가 커서, 상수 하나로 잡으면 과대면
무용하고 과소면 위험하다 — 그래서 커널의 필터는 구간 하나다(`Bounded` — `frame3::orient3d_iv`·`det3_iv`·`cramer_iv`,
손으로 유도한 반경 없음). tol 0 인 점은 반경 0 이라 그 구간이 저절로 좁다.

**조합 항등식이 먼저다.** `orient3d` 는 `j ∈ {p, q, r}`, `cmp_coord` 는 두 점이 공유하는 평면이 그 좌표를 고정할 때
(`Judge::shared_planes_hold` — `forced_rank_deficient`, Frobenius–König), `plane_pair_dir_sign` 은 법선이 번호나 알려진
0 만으로 종속일 때(`Judge::normals_forced_dependent`) 0 이다 — 알려진 0 은 `Judge::world_normal_zeros`. 논증·전제·측정과,
0 을 술어에 보내면 무엇이 되는지는 그 함수들과 `Judge::orient3d` 의 doc.

**상승은 계산하고, 캐시는 값이다.** 다음 정밀도는 모자란 비트 `log₂(간격/목표)` 로 한 번에 점프하고(`escalate` doc; 배증은
「가지 말 것」), 상승 중에 구별하는 상태들은 `Gap` doc. 캐시는 판정 결과가 아니라 **점의 값**이다 — 연산의 정밀도에서 한 번
실현해 여러 판정이 빌린다(`WitnessPoint::coord_hp`; 결과 메모는 「가지 말 것」). 평면의 구간 계수와 법선의 알려진 0 은
`Judge` 가 첫 사용 때 채우고, 트레이트 밖에 두는 까닭은 `Judge` 의 `iv`·`zeros` 필드 doc.

**평면 증인은 정의를 소유하고, 「회전됨」은 별개 사실이다.** `Witness::tri_pt3`(늘 있는 캐시)와 `Witness::is_rotated`(경로
신호)는 서로를 유도하지 않는다 — 까닭은 두 doc, 정의를 판정마다 다시 지으면 드는 비용은 `plane_def` doc. 둘은 `nacre-ops` 의
`collect_planes` 가 함께 정해 싣는다(`FaceInfo::rotated` doc; 증인 삼각형의 철자는 `rotated_vertex::plane_witness`).

### 정밀도는 모델이 정한다

`nacre-judge::kernel::frame3`(`ProofStandard`·`trial_bound`·`precision_for`·`coincidence_for`), `nacre-ops::planes::standard`
(`standard_for`·`JUDGE_PREC_CAP`·`CLIMB_HEADROOM`), 거절 `RejectReason::PrecisionBudget`(`arrangement::setup`).

판정의 오차 반경은 `C·2⁻ᵖʳᵉᶜ` 이고 `C` 는 모델의 성질(회전 이력 1회당 약 1비트, 좌표 크기)이라, 연산마다 한 번 읽어
필요한 비트를 계산한다 — 점마다 `trial_bound`, 그 최댓값을 `precision_for` 가 비트로, `standard_for` 가 표 전체에서.
고정 정밀도가 회전 이력의 길이를 조용히 정하는 것은 `precision_for` doc 과 「가지 말 것」 «고정 판정 정밀도».

- 문턱은 모델 단위의 길이이지 비트 수가 아니고, `coincidence` 는 tolerance 가 아니라 «가깝다고 증명되어야 일치»다 —
  `ProofStandard` doc. 값(f64 출력 해상도 `scale·2⁻⁵²` 의 두 워드 아래)의 유도는 `coincidence_for` doc, 규칙은 「숫자 규칙
  — 일곱 개」 5.
- **모델의 깊이와 한 판정의 얇기는 별개 예산이다.** 깊이는 `JUDGE_PREC_CAP`(정확성이 아니라 비용 한계 — 넘으면 판정 전에, 평면 표가 서자마자
  `PrecisionBudget`), 한 판정의 얇기는 `CLIMB_HEADROOM`(모델 정밀도에 대한 상대값)이 든다 — 값의 근거·측정과 공유하면
  안 되는 까닭은 두 상수의 doc, 비용 측정은 `a_rotation_history_past_the_budget_is_rejected_by_name` 의 doc.
- 설정으로 노출할 후보는 일치 정밀도 하나뿐이고 나머지는 거기서 유도된다(`standard_for` doc).

### 판정 결과와 보고

`Decision`·`Decision::orient`(`kernel::frame3`), `Judge::record`·`Notes`·`Evidence`·`Site`(`predicate`),
`arrangement::undecided_reject`, `RejectReason::{JudgeExhausted, DegenerateWitness}`, `boolean_with_report`·`BoolReport`
(`nacre-ops::boolean`). 규칙은 「숫자 규칙 — 일곱 개」 5.

판정은 `Decision` 으로 답한다 — 갈래마다의 뜻(증명된 부호, 증명된 일치, 비트가 모자람과 그때 세운 상한, 비트로 풀리지
않는 퇴화)은 `Decision` doc. 기하는 `Decision::orient` 의 부호만 읽어, 결론 없는 결과는 `Zero` 로 접혀 «같다고 취급하라»가
되고 보고 채널은 기하 결정을 하나도 건드리지 않는다. 증명된 부호가 아닌 판정은 `Judge::record` 가 연산의 `Notes` 에
적고(전역도 thread-local 도 아니고 진단 전용 — `Notes`·`Judge` doc), 연산 끝에서 이렇게 읽힌다:

- `undecided_reject` 가 `Exhausted`·`Degenerate` 를 원인 이름의 거절(`JudgeExhausted`·`DegenerateWitness`)로 바꾼다 —
  하류 증상보다 **먼저** 본다(`arrangement::boolean` 끝의 주석). 두 원인의 차이는 그 변종들의 doc.
- `BoolReport` 가 일치들을 근거와 함께 돌려준다 — 첫 항목은 평면 클래스 병합이다(`BoolReport::merges` doc). 병행 진입점
  `boolean_with_report` 를 두는 까닭은 그 doc.

**묻지 않고, 근거를 붙여 보고한다.** 사용자에게 확인하는 길은 「가지 말 것」 «「애매하면 사용자에게 묻는다」». 사람이
아는 것은 *의도*이고 의도는 이미 스크립트에 있다 — 같은 회전을 공유하면 상쇄로 정확히 답하고, `Surface` 핸들을
공유하면 핸들 비교다. 보고는 사용자가 의도치 않은 일치를 알아차리고 설계 치수를 고치게 하는 진단이지 커널의 결정을
대신하는 프롬프트가 아니다(`BoolReport` doc).

### 근본 한계

어느 한 모듈의 것이 아니라 CIP 전체의 범위다 — overview 「절대 원칙」, 그것을 지키는 장치는 `ProofStandard`·`Decision`(위 절들).

회전이 들어가면 **서로 다른 경로로 같은 위치에 도달해도 실현이 다를 수 있고, exact 술어로도 그 다름을 그대로 반영한다** — 술어는 「주어진 입력에 대해」 정확할 뿐 입력의 근사를 고치지 못한다. **어떤 유한 정밀도로도 보장은 불가**하다(초월수의 상등은 유한 정밀도로 결정 불가) — 정밀도를 올리는 것은 답을 얻는 방법이지 완비성을 사는 방법이 아니다. 다만 **실무적으로는 충분하다**: 상용 CAD 의 uncertainty 가 약 1e-9 인데 f64 는 1000mm 점을 10000번 회전해도 약 1e-10 이고, 커널이 실제로 오르는 폭(회전 3200회에서 3456비트)은 그보다 수백 자리 아래다. **이론적 완벽함은 포기하고 실무적 충분함을 tol 로 보장한다** — 그 tol 을 **측정된 값**으로 들고 다니는 것이 tolerance-fudge 와 다른 점이다. 부호를 확정 못 하는 잔여를 커널은 조용히 0 으로 추측하지 않는다 — 증명되면 일치로 처리해 보고하고, 좁힐 수 없으면 원인을 이름 붙여 거절한다.

### 도메인별 적용

`nacre-judge::cylinder`(`PiercePoint`·`Judge::pierce`), 「이차곡면 판정 — 원통 우선」, 「계획: M7 SSI — 작은 loop·접선 탐지 전략」.

- **평면·이차곡면:** 판정이 다항식이므로 CIP 가 그대로 얹힌다. **CIP 의 가능 여부 = 간접 술어의 가능 여부**다(CIP 는 술어 위의 층이다). 커널의 반경은 구간 산술(`Bounded`·`HpBounded`)이 내므로 술어가 고차가 되어도 손으로 유도할 공식이 없고, 상수 필터를 둘 때만 그 술어의 오차 한계를 유도한다(`nacre-predicates` 필터 doc 의 표).
- **자유곡면(M7):** 메시 조합 판정의 필터 + 뉴턴 스냅백 점의 tol 추적·고정밀 재수렴(폭은 같은 방식으로 계산해 정한다). **SSI(교차를 찾는 것)는 CIP 밖**이다 — 위상 존재 문제이지 tol 문제가 아니다.

**계보.** CGAL Lazy_kernel(구간 근사 + DAG 연산 이력 + 애매하면 exact 재평가 + 캐싱)과 구조가 같다. nacre 는 DAG 대신 **점이 드는 모션 사슬**을 쓰고, 구간은 같다(`Bounded` — 좌표·평면 계수·행렬식이 한 구간 산술을 걷는다). 특화 커널의 이득은 술어의 구조에서 온다 — 술어가 소수·고정이고, Attene 2020 이 선·평면 교차에서 CGAL lazy 보다 빠른 이유와 같다.

## 계획: M7 SSI — 작은 loop·접선 탐지 전략

「내/외 판정」은 *찾은* 교차를 어떻게 분류하나이고, 이 절은 *찾는* 일 — M7 의 진짜 도박 — 이다. 방법은 M7 진입 시 그 시점의 최신으로 비교해 정한다.

**SSI 난제 넷 중 「곡면 결함」과 「연산 대상」을 가른다.** M7 이 도박인 이유는 부호 판정이 아니라 **교차 곡선(불리언 seam)을 처음부터 찾는 SSI** 다.

- **① 작은 loop** — 두 곡면이 좁은 영역에서 작은 닫힌 고리로 교차. 격자·메시 간격보다 작으면 놓친다 → seam 누락 → 불리언 오답. **최대 난제.**
- **② 접선 교차** — 곡면이 가로지르지 않고 스치듯 닿음(구가 평면에 접). 교차 유무가 모호하고 거리 0 판단이 f64 로 어렵다.
- **③ 특이점/cusp** — 교차 곡선이 갈라지거나 꺾임(직교 원통 X자). marching 이 가지를 놓치거나 branch jumping.
- **④ 자기교차·겹침** — 곡선이 아닌 면으로 포개짐, 또는 곡면이 스스로 꼬임.

**처리 정책.** ③④ 중 **곡면 자체 결함**(곡면 자기교차, 곡면 고유 특이점)은 생성 시 검사해 거부한다 — 연산 대상이 아니다. **겹침**(두 정상 곡면이 면으로 포개짐, coplanar 류)과 **교차 cusp**(정상 곡면들의 교차에 생기는 특이)는 곡면 결함이 아니라 정상 연산 대상이라 거부할 수 없다 — 별도 케이스로 처리한다. 실질 난제는 ①② + 겹침/교차 cusp.

**①② 탐지 — 계보와 강도 순(약 → 강).**

- **거리 함수 임계점**(기준선, Patrikalakis 계보). 한 곡면 → 다른 곡면의 방향 거리 함수를 만들고 그 gradient 벡터장의 위상(회전수·Poincaré index)으로 임계점을 탐지한다. 정통이나 (a) 거리 함수 계산이 무겁고 (b) 임계점들이 가까우면 강건성이 떨어진다.
- **법선 평행점**(보장 있음). 정리: 두 비특이 곡면이 닫힌 loop 로 교차하면 **양쪽 법선이 평행한 직선이 반드시 존재**한다. 법선 평행점을 전부 찾아 그 지점에서 세분하면 **loop 를 놓칠 수 없음이 증명된다**(heuristic 이 아닌 보장). Bézier normal vector surface 로 계산한다.
- **Winding number**(유력). winding number 이론 + 세분 결합. 위상 불변량이라 **격자 해상도와 무관하게** loop 존재를 판정한다 → 격자보다 작은 loop 도 방어. 방향 거리장 + gradient 를 격자에 두고 고립 임계점·비고립 임계곡선을 계산, 그로부터 각 분기 시작점을 뽑아 곡선을 추적한다. (ACM TOG 2026, "A Robust and Efficient Intersection Algorithm for NURBS Surfaces: Handling Small Loops and Tangent Intersections")
- **하이브리드가 최신 흐름:** winding number(탐지) + 거리/법선(정밀 위치) + 적응 세분(국소화).

**파이프라인 스케치.**

1. 두 곡면 간 거리를 영역별로 계산한다(**GP 아님** — 곡면을 해석적으로 아니 직접 평가가 정확·빠르다; GP 는 함수 미지 시의 도구라 여기선 근사만 더한다).
2. 거리 **국소 최소**(줄었다 늚) 영역 = 교차 후보. 「나란히 가까이 붙은 넓은 영역」(거리는 작지만 최소가 아님)은 제외해 세분 폭발을 막는다.
3. 후보 영역만 적응 세분, 재귀(간격 ≥ 새 최대 변까지). 세분 리미트를 둔다.
4. 리미트 도달 → 뉴턴법: 수렴 = 교차 존재, 발산 = 미교차. 뉴턴 스냅백 점은 tol 있는 점이므로 **CIP 에 통합**한다(수렴 잔차 = tol, 애매하면 정밀도를 올려 재수렴).
5. 리미트에서도 확신 불가(작은 loop 는 「다 찾았다」의 수학적 보장이 근본적으로 불가) → **`Rejected` 정직 거부.**

**출발점:** Li·Yang·Jia, "Advances and challenges in surface–surface intersection computation — An overview", Computer-Aided Design 193:104039, 2026 — SSI 분야 전체 개관. (관련: Li·Jia·Chen, "Fast Determination and Computation of Self-intersections for NURBS Surfaces", ACM TOG 44(2), 2025 — ④ 자기교차 판정·거부용.)

**CIP 와의 관계.** CIP 는 **부호 판정** 층이다. M7 에서 CIP 가 닿는 곳은 (a) 메시 조합 판정(내/외 분류)의 필터, (b) 뉴턴 스냅백 점의 tol 추적·정밀도 상승 — 둘 다 **정밀화·판정**이다. **SSI** 자체는 위상 존재 문제라 CIP 밖이다. 「찾은 것을 정밀하게」(CIP·스냅백)와 「못 찾은 것을 찾기」(SSI)는 다른 일이고, SSI 가 성공한 뒤라야 스냅백·판정이 의미 있으며 SSI 실패 시 `Rejected` 다.

## 비목표: STEP 가져오기

단순 미구현이 아니라 nacre 정체성과 맞는지부터 물어야 하는 항목이다.

1. **경계가 커널과 앱을 가로지른다.** 파일 파싱·워크플로 결정은 애플리케이션 레이어이고, 파싱된 형상을 유효한 b-rep 으로 재봉합(healing)하는 것만 커널 몫이다. import 는 커널 단독으로 완결되지 않는다 — GD&T·메타데이터를 커널에서 뺀 것과 같은 논리.
2. **export 보다 근본적으로 비싸다.** 넓은 스키마(AP203/214/242 각 에디션) 파싱 + 진짜 본체인 healing(외부 파일은 면들이 각 CAD 의 tolerance 로 느슨히 붙어 있어, 공유 엣지·단일 참조 정점 위상으로 들이려면 어긋난 면을 꿰매야 한다 = tolerant modeling 문제 그 자체)과 근사 → 정확 승격(외부 근사 곡면·교차곡선을 nacre 의 정확 기하 진실로 올리기).
3. **히스토리 없는 dumb solid.** nacre 의 정체성은 「모델 = 연산 로그를 재생한 결과」인데 외부 형상에는 그 로그가 없다. 재생성·히스토리 편집이 불가한 이등 시민이 된다 — 상용 커널도 STEP import 를 dumb solid 로 취급해 보기·측정 수준으로 제한한다. topological naming 을 포함한 비목표들을 한꺼번에 끌어들인다.

**결론:** import 는 healing·근사→정확 승격·naming 을 동반하는 별도 대형 과제이며 앱 레이어를 전제한다. v1 비목표다. 그 뒤에도 「구현 가능한가」가 아니라 「히스토리 없는 dumb solid 를 받아들이는 것이 nacre 핵심 가치와 맞는가」를 먼저 묻는다. (자기 출력을 step-io 리더로 되읽는 라운드트립 테스트는 「방금 내가 쓴 것」을 읽는 것이라 healing 문제가 없다 — 외부 파일 import 와는 별개 사안이고 오라클 검증용으로 유효하다.) 입력이 필요해지면 step-io 를 별도 리더 크레이트로 붙이며, 출력과 입력은 별개 크레이트다.

## 비목표: 와이어프레임/서피스 모델링 — 단 확장 접합면은 지킨다

nacre 는 solid-first 이며 always-closed(모든 결과가 닫힌 솔리드)를 전제한다. 사용자가 점·곡선·열린 면(open shell) 같은 개방 형상을 독립적으로 만들고 편집하는 wireframe/surface 모델링은 v1 범위 밖이다: (a) 개방 요소 + sewing 은 always-closed 와 정면충돌하고(열린 면은 validate 가 즉시 불량 판정), (b) 자유곡면 sewing 은 healing 문제라 난이도가 곡면 불리언급이며, (c) 「정밀 기계 부품 → STEP」 목표에는 솔리드 모델링으로 충분하다. (점·선·면은 커널의 내부 벽돌로서는 전부 존재한다. 비목표인 것은 「사용자가 개방 요소를 독립적으로 생성·편집하는 기능」이다.)

이것은 「ops 에 함수를 추가하면 자연히 되는」 것이 아니다. 개방 형상 지원은 (1) 떠 있는 점·곡선·열린 면을 담는 자료구조(topo 확장), (2) validate 의 유효성 이원화(「닫힌 솔리드 = 엄격 검사」 vs 「개방 요소 = 느슨한 검사」), (3) sewing/healing 연산을 전제한다.

**지키는 규율**(미래의 「전면 재작성」을 「국소 확장」으로 바꾸는 장치):

- **(a)** always-closed 가정을 코드 전반에 흩뿌리지 않는다 — 오직 validate 의 명시적 규칙으로만 둔다. 다른 함수는 「면은 무조건 어떤 솔리드에 속한다」 같은 암묵 가정으로 지름길을 쓰지 말고 Adjacency 등을 통해 조회한다.
- **(b)** validate 는 개별 불변식 검사의 모음이다(하나의 거대 함수 금지) — 유효성 이원화가 「규칙을 켜고 끄는 일」이 되게.
- **(c)** 진실은 연산 로그다 — 형상 종류 확장은 연산(`Operation`) 추가이지 기존 데이터의 마이그레이션이 아니다. append-only 라 옛 연산의 의미는 고정이고 확장이 기존 로그를 건드리지 않는다.
- **(d)** `Model` 최상위가 solids 만이 아니라 향후 개방 요소도 담을 구조적 여지를 갖는다(free-element 컨테이너를 실제로 만들지는 않되, Solid 를 유일 최상위로 못박아 확장을 막지 않는다).

**상세 설계하지 않는 것:** 개방 형상의 정확한 자료구조 표현, sewing 알고리즘의 tolerance 처리, validate 이원화의 구체 규칙 — 곡면 도메인을 실제로 구현하며 얻는 지식(relaxation 이 어디서 깨지는지, tolerance 가 실제로 어떻게 도는지)이 있어야 정확하다. 경계(접합면)는 미리, 내용(구현)은 그때. 이 규율들은 wireframe/surface 와 무관하게도 좋은 설계(깔끔한 검증, 공짜 replay/undo)라 부담이 아니다.

## 미결

- 트리밍 곡면의 pcurve 표현 도입 시점.
- Sketch 제약 솔버의 범위(무제약 프로파일만 있다).
- OpRef 계보 참조의 도입 시점과 직렬화 포맷의 여유분.
- `Store` 스냅샷·직렬화 포맷(자체 vs STEP 재활용). 세션 중 메모리 관리는 **compact 보다 재구축(rebuild-from-log) 우선**이다(재구축 = 주력 정리·undo 유지, live 폐포 필터링 = 저장, compact = 비상 회수).
- OCCT history → 출처 매핑의 실제 충실도(측정 대상).
- 멀티스레딩 경계: `Store` 가 `&mut` 독점인 설계라 연산 단위 병렬은 없다(의도적 단순화).
- tess 의 seam polyline 샘플 밀도와 법선 비분리 규칙의 세부.
- op-log 소유자 `Document`(회전 지원의 확정 결정 4).
- CIP 의 측정 질문: 고정밀 층의 속도, 동적 필터의 실제 성공률(입력 tol 이 있을 때 얼마나 자주 상승하는가), 실무 형상에서 회전 tol 이 판정 경계에 얼마나 근접하는가, 고차(이차곡면) 술어에서의 필터 성공률.
- 구·원뿔, 일반 이차곡면쌍의 범위 선과 판정 방법.
- M7 SSI 탐지 방법의 선택.

## 검증·오라클 인프라 (`nacre-validate`, `nacre-oracle`)

### `validate` — 모델 불변식 검사

```rust
pub fn validate(model: &Model) -> Vec<Violation>;   // 빈 벡터 = 유효
```

`validate` 는 위반을 **전부** 모아 돌려준다(누가 부르는가는 overview 「정밀도 인프라」). 인접 인덱스는 안에서 짓는다. 검사하는 것은 `Violation` 의 변종 그대로다:

| 변종 | 불변식 |
|---|---|
| `DanglingReference` | 모든 handle 이 대상 store 범위 안(`target_index < target_len`). |
| `OpenLoop` | 모든 Loop 가 닫힌다: `end(he[at]) == start(he[at+1])`. |
| `NonManifoldEdge` | 모든 엣지가 half-edge 에 정확히 2회 쓰인다(0 = 고아, 1 = 열린 경계, 3 이상 = 비다양체). |
| `NonOpposedEdge` | 엣지를 쓰는 두 half-edge 의 `forward` 가 반대다(방향 일관성). |
| `NonManifoldVertex` | 정점의 면 링크가 원 하나다 — 모든 엣지가 다양체여도 두 솔리드가 꼭짓점 하나로만 닿는 «핀치»를 잡는다(`nacre_topo::nonmanifold_vertices`). |
| `VertexCarrierMismatch` | 정점 정의의 운반 곡면 종류가 변종의 주장과 맞는다: `ThreePlane` 은 평면 셋, `OnSeam` 은 평면 둘이 아니다(평면 둘의 교선 위 점은 세 평면 교차로 표현되는 진실이다). 불변식은 변종마다 따로다. |
| `EdgeCarrierMismatch` | 엣지가 진술한 `Edge::surfaces` 쌍이 그 엣지를 쓰는 두 면의 곡면 다중집합과 같다 — 두 면이 한 곡면 위에 있으면 쌍이 그 곡면을 든다(다른 하나는 선을 담는 평면). 자기인접 `[s, s]` 는 어느 곡면에도 금지(seam 간선은 없다). 운반 곡면은 유도가 아니라 진술이므로 불일치는 생산자 버그이고, 틀린 운반에서 유도된 곡선 캐시는 조용히 틀린 기하가 된다. |
| `VertexNamesAbsentSurface` | 살아 있는 정점마다 정의가 부르는 곡면이 그 정점을 쓰는 면들의 곡면이다 — 꼭짓점은 제 면들이 만나는 자리이고, 다른 곡면 위에 우연히 놓였다고 그 곡면으로 이름 짓지 않는다(「점의 이름」의 첫 원리). `EdgeCarrierMismatch` 의 정점 짝이되 더 엄격하다: 한 평면 위 두 면 사이의 호 간선은 원통을 담체로 진술하지만, 그 호 끝의 정점은 거기 닿는 면들로 이름 짓는다. |
| `VertexOffCurve` | 엣지에 묶인 정점이 그 엣지 곡선 위에 tol 이내. |
| `VertexOffSurface` | 루프 정점이 그 면의 곡면 위에 tol 이내. |
| `VertexOffDefinition` | **모든** 정점의 캐시 좌표가 자기 정의가 말하는 곡면 하나하나에서 tol 이내. 정의가 진실이고 좌표는 캐시이기 때문이다. tol 은 실측값이 있으면 그것, 없으면 `EPS_CONSTRUCTED`. |
| `EulerParity` · `NegativeGenus` | 오일러–푸앵카레: `V − E + F = 2(S − G) + L_i`. `χ = V − E + F − L_i` 가 홀수면 정수 해가 없고(`EulerParity`), 짝수인데 `G = S − χ/2 < 0` 이면 위상적으로 불가능하다(`NegativeGenus`). |
| `ShellMisoriented` | 셸의 부호 있는 자기 부피가 제 역할(`ShellRole`)의 부호다 — 바깥 셸은 양, 공동(내부 void)은 음. 통째로 뒤집힌 셸은 엣지 대향도 V/E/F/S 도 그대로라 다른 검사가 전부 놓친다. 평면 셸만 검사한다. |
| `FaceMisoriented` | 평면 면의 루프 감김과 진술한 방향의 정합. 루프의 증인(다각형의 뉴웰 면적 벡터, 또는 닫힌 림의 원)과 면이 진술한 바깥 법선(`plane.normal()` × `orientation`)의 단위 내적 `cos` 를 잰다: **외곽 루프는 일치**(`cos ≈ +1`), **내부 루프는 반대**(`cos ≈ −1`). `≈ 0` 은 자기 평면을 펴지 못하는 루프. 불리언의 `collect_planes` 가 `debug_assert` 로 지키는 불변식의 release 쪽 그물이다 — 틀린 플래그·감김은 엣지 대향과 오일러를 모두 통과해 «어느 쪽이 재료인가»로 곧장 흘러든다. |
| `CylinderTruthCacheMismatch` | 원통의 정확한 진실(`CylinderDef`)과 f64 캐시가 같은 원통을 말한다. `field`(`"radius"`·`"origin"`·`"dir"`·`"ref_dir"`)의 성분별로 `\|def − cache\| ≤ CYL_TRUTH_EPS · max(1, \|def\|, \|cache\|)`, `CYL_TRUTH_EPS = 8ε`. 한 곡면의 두 서술은 어긋날 수 있고, 치료는 생산자를 믿는 것이 아니라 소비자 쪽 사후조건이다. 건강한 경로는 절대 1e−17 이내이고, 이 그물이 노리는 결함(seam 타이브레이크 오류로 `ref_dir` 이 ~90° 도는 것, 단위화 누락으로 `dir` 이 길이만큼 커지는 것)은 성분을 ~10¹⁵·ε 움직인다. |
| `VertexNoPoint` | 세 평면 정의의 이름이 한 프레임을 말하는 곳에서 세 평면이 한 점에서 만난다(`Model::vertex_named_meet` 가 `NoPoint` 가 아니다). |
| `DuplicateSurface` | 아레나의 두 핸들이 한 interning 열쇠를 유도하지 않는다 — 열쇠를 진실에서 다시 매긴다(「열쇠는 진실이 고른다」). |
| `SurfaceNameMismatch` | 곡면마다 저장된 이름이 다시 매긴 열쇠의 이름이다(열쇠가 이름이 아니면 없음) — 빠진 이름·남는 이름·틀린 이름. |

오일러 식의 `L_i`(면 내부 루프 수) 항은 빠뜨릴 수 없다. 관통구멍 정육면체로 검산하면 V16 − E24 + F10 = 2 = 0 + L_i(2) 다 — 내부 루프 항이 없는 축약형은 정상 모델을 불량으로 판정한다.

**validate 는 store 전체가 아니라 live 도달가능 셀을 센다**(supersede 의미론). 오일러의 V·E·F·S·L_i 는 `store.len()` 이 아니라 `live_solids` 에서 도달 가능한 정점·엣지·면·셸·내부 루프 수이고, 다양체 조건(엣지 정확히 2회)도 «전역 2회»가 아니라 «도달가능 집합 안에서 2회»다. supersede 된 옛 셀이 아레나에 남으므로, store 길이로 세면 죽은 셀이 오일러를 깨고 죽은 면의 엣지가 다양체 조건을 깬다. store 에 떠 있는 stray 셀은 불량이 아니라 «죽은 아레나 항목»으로 정당하게 무시된다. **참조 무결성 검사만 store 전체를 훑는다** — 살아 있든 죽었든 handle 이 범위 밖이면 버그이기 때문이다. 이 검사는 `.index()`/`.len()` 만 쓰고(`Store::get` 을 부르지 않는다) 가장 먼저 돌며, 위반이 있으면 거기서 단락한다. 그래서 뒤따르는 도달가능성 순회는 항상 범위 안이다.

**테셀레이션 불변식은 `nacre-tess` 가 자기 테스트로 지킨다**(`validate` 의 변종이 아니다): 출처 정합 — `TessOrigin::OnFace { uv }`·`OnEdge { t }` 의 평가값이 `pos` 와 일치, 틈 없음 — 인접 면 삼각분할이 공유 엣지 polyline 정점을 정확히 공유(테셀레이션 절).

### 속성 기반 테스트(proptest)

랜덤 유효 연산열을 생성해 건다: (1) validate 통과, (2) replay 멱등성 — 같은 로그·같은 cfg → 동일 모델, 그리고 더 강한 «`replay(log)` == 세션 모델, 인덱스까지»(`nacre-ops/tests/invariants/replay.rs`), (3) 변환 불변량 — 강체변환 후 부피·면적 보존, (4) 불리언 대수 — `A ∪ A = A`, `A ∩ ∅ = ∅`, `vol(A∪B) + vol(A∩B) = vol(A) + vol(B)`.

### OCCT 오라클(`nacre-oracle`)

nacre 의 결과를 STEP 으로 내보내 OCCT 가 채점한다: 부피·면적·면 개수·바운딩박스. 불리언은 두 입력 솔리드의 STEP 을 OCCT 가 직접 `bfuse`/`bcut`/`bcommon` 한 결과의 값과 nacre 결과의 값을 diff 한다. 정답지를 든 채 개발하는 장치이며, AI 가 생성한 코드의 «그럴듯하지만 틀림»을 잡는 주 방어선이다. **nacre 쪽 부피·면적은 `nacre-props`(정확 기하 발산정리, 해석적)가 계산하고, 오라클은 그 값을 OCCT 가 같은 STEP 에서 독립 계산한 값과 diff 한다** — 완전히 별개인 두 구현의 일치가 양쪽을 교차검증한다. **일치는 형상이 설명한 그것인지 말하지 않는다** — 두 커널이 같은 STEP 을 재므로 픽스처가 어디에 놓이든 맞는다. 그래서 오라클 픽스처는 자기 손 계산 값도 단언하고, 공용 형상(`nacre_ops::fixtures::pocketed_cube`)은 자기 전제를 단언한다.

계획: 허용 편차를 넘은 연산열을 최소화(shrink)해 리포트한다.

**OCCT 는 오라클 전용이며, 연동은 out-of-process 다.** OCCT 를 링크하지 않고, `tools/occt-helper/` 의 헬퍼 프로세스와 STEP 파일로 주고받는다(수송 계층 = `nacre-step`). 이유: C++ 빌드가 Rust 워크스페이스에서 완전히 사라지고, OCCT 가 악조건 입력에서 크래시해도 커널 프로세스가 아니라 헬퍼만 죽는다(크래시 격리 — 수천 회 돌리는 오라클에는 필수). 오라클은 채점만 하므로 history 손실·호출 오버헤드는 무관하다. `nacre-oracle` 은 발행하지 않는 dev 전용 크레이트다.

헬퍼 프로토콜은 고정이다:

```
occt-helper props <in.step>
occt-helper <fuse|cut|common> <a.step> <b.step> [out.step]
```

- `props` 는 솔리드 하나를, 불리언 명령은 **결과**의 성질을 보고한다. `cut` 은 `A − B`. `out` 을 주면 결과를 STEP 으로도 쓴다.
- stdout 은 모든 명령이 같은 key-value 7줄이다: `volume <v>` · `area <a>` · `faces <n>` · `bbox_min <x> <y> <z>` · `bbox_max <x> <y> <z>` · `centroid <x> <y> <z>`(부피 중심) · `valid <1|0>`. `valid` 는 DRAWEXE `checkshape` 를 읽은 형상(`a_1`, 불리언이면 `b_1` 도)과 불리언 결과 하나하나에 돌린 판정이고, 무효인 형상마다 `fault <a|b|r> <BRepCheck_코드>` 줄이 붙는다(`valid 0` 이면 재는 줄은 잰 만큼만 나온다). 판정은 `checkshape` 마다 표지로 감싼 구획에서 읽는다 — 출력 전체를 한 번 grep 하면 «입력 무효 + 결과 유효»가 유효로 읽힌다. 값은 헬퍼 쪽에서 계산한다. 고정 스키마라 파싱 의존성이 0 이고 Rust 쪽 파서는 하나다.
- **무효는 오라클의 문이 거절한다**(`OracleError::Faulty`, 불평마다 `Fault { shape, code }`). 오라클 테스트는 형상마다 방향을 단언하지 않는다 — 면 하나를 뒤집어도 부피는 그대로(정육면체 1)이므로, 문이 거절하지 않으면 부피만 보는 오라클 테스트는 방향 오류를 통과시킨다. OCCT 의 불리언 결과가 무효여도 거절이다: 무효인 정답지로는 채점하지 않는다. 문이 거절함을 `a_wrong_orientation_reads_as_faulty` 가 잠근다(면·간선·경계·원통 옆면·void 를 텍스트에서 뒤집은 파일과 불리언의 두 입력). 오라클 테스트 밖의 인구도 한 번 쟀다 — census 결과 솔리드 349(곡면 168·void 1)를 이 읽기로 재면 전부 유효이고 부피·면 수가 nacre 와 일치하며, 그 파일들의 면·간선·경계·void 를 뒤집으면 잡힌다(11/11). 관문에는 없다.
- **읽기는 방향 수리만 끈다.** OCCT 의 STEP 읽기는 기본으로 형상을 수리하고(`FromSTEP.exec.op : FixShape`), 그 수리가 면·간선·껍질·void 의 방향 오류를 조용히 고친다. 헬퍼는 세션 안 `param` 으로 `FixOrientationMode`·`FixShellOrientationMode`·`FixFaceOrientationMode` 셋만 0 으로 둔다. 나머지 수리는 파일에 없는 PCURVE 를 짓고 seam 없는 옆면에 seam 을 넣으므로(`FixMissingSeam` — 「seam 엣지 방식」) 남긴다 — 수리를 통째로 끄면 곡면 위 면이 늘 무효다(「가지 말 것」). 호출자가 `CSF_STEPDefaults` 로 `exec.op` 를 비우면 FixShape 가 통째로 안 돌고 세 `param` 은 쓰이지 않는다. 어떤 설정으로도 안 보이는 것은 직선 `EDGE_CURVE` 의 `same_sense` 하나이고, brep-to-step 이 직선 방향을 끝점에서 정하므로 nacre 가 틀리게 쓸 수 없다.
- exit code: `0` 성공 · `1` 기하 실패(파일 없음·읽기 불가·빈 형상) · `2` 크래시, 또는 DRAWEXE 가 방향 설정 셋을 다 받지 않음(모르는 이름은 오류를 찍고 exit 0 으로 계속 돌아 수리가 조용히 돌아오므로 헬퍼가 받은 수를 센다).
- 입력은 **단일 솔리드** STEP 이어야 한다(솔리드 하나의 `to_step_solids`). 전송 가능한 루트가 여럿이면 OCCT 가 `x_1, x_2, …` 로 읽고 불리언은 `x_1` 만 쓴다.

구현은 Homebrew OCCT(`brew install opencascade`)의 `DRAWEXE` Tcl 셸 배치 실행이다(`stepread → [bfuse/bcut/bcommon →] vprops/sprops/nbshapes/bounding`, 콘솔 출력 파싱). 콘솔 긁기가 버전 차이로 취약해지면 같은 줄을 내는 얇은 C++ 헬퍼(`STEPControl_Reader` + `BRepGProp` + `BRepBndLib`)나 uv 관리 Python 환경 + OCP 휠로 갈아탄다 — 프로토콜만 지키면 nacre 쪽 코드는 무변경이다.

### indirect predicates 검증

indirect 술어는 틀려도 대부분의 입력에서 맞는 답이 나와 버그가 가장 숨기 쉬운 조각이므로 오라클을 두 겹으로 건다. (1) **direct 술어 대조(property test)**: 좌표를 아는 케이스에서 정의에서 좌표를 실현해(`realize_vertex`) `geometry-predicates` 의 direct 술어(orient3d 등)에 넣은 부호와, 같은 점을 implicit point 로 둔 indirect 술어의 부호가 일치하는지 대조한다. (2) 계획: **Attene C 구현 실행 대조(dev-only)** — 저자 참조 구현을 별개 프로그램으로 실행해 답만 비교한다(아래 라이선스 규율 — `tools/` 격리, 링크 금지, OCCT 오라클과 동일 논리).

### STEP 라운드트립(자기 출력 되읽기)

nacre 가 쓴 STEP 을 step-io **리더로 되읽어** 구조를 비교한다. 외부 파일 import 가 아니라 «방금 쓴 것»이므로 healing 이 필요 없다. 교차 곡선은 일반 SSI 단계(M7)에서 처음 생긴다. 계획: 그 자리는 `EdgeCache` 를 `{curve, err}` 로 키우는 쪽이고, 교차 곡선이 검증 대상이 되는 것도 그때다.

## 마일스톤 사다리

각 단계는 «동작하는 것»을 남기고 끝난다. 3층에 막혀 전체가 멈추는 구조를 피하는 배치다.

**M1 — 뼈대.** `nacre-math`, Store/Handle, 평면·직선만으로 정육면체를 손으로 조립하고 `validate` 를 세운다. 시각 확인은 **기존 뷰어에 위임**한다(전용 뷰어 크레이트를 만들지 않는다): 정상 결과는 STEP 출력 → step-loupe(구조+검증), 중간·깨진 상태는 Tessellation 을 OBJ/STL 로 덤프 → 맥 미리보기(또는 MeshLab·f3d). 커널 내부까지 보는 인터랙티브 디버그 뷰어(면 클릭 → Handle·정의, 법선 화살표, 엣지 polyline·tolerance 공, validate 위반 하이라이트, 로그 스텝별 재생)는 워크스페이스 밖 별도 앱이다 — **발견된 점**과 tolerance 처럼 STEP 에 담기지 않는 내부 정보가 생기는 불리언 단계의 디버깅 생명줄이다.

**M2 — 스케치와 돌출, STEP 내보내기.** 2D 프로파일 → extrude, 연산 로그와 replay, 오일러 연산. STEP 은 어댑터다: 커널 Model → AP242 엔티티 번역만 자체 구현하고 직렬화는 형상 전용 작성기 brep-to-step 에 위임한다(커널은 백엔드 무지 — 헬퍼 프로토콜과 같은 격리). 머리말의 시각은 호출자가 준다 — 커널은 시계를 읽지 않는다(같은 모델과 같은 시각이면 같은 바이트이고, wasm 에는 시계가 없어 시스템 시각을 물으면 패닉한다). 파일의 불확도(`UNCERTAINTY_MEASURE_WITH_UNIT`, «이보다 가까운 두 점은 한 점»)는 내보내는 수 위에서 증명한 맞물림의 위 한계다(`nacre_step::uncertainty` — 꼭짓점과 그 곡선·곡면, 간선 곡선과 그 면들의 곡면 사이 틈을 위로만 반올림한 구간으로 재고, 적힌 좌표의 f64 간격의 반 밑으로는 적지 않는다). 참값과의 거리는 그 값이 아니라 `refine_caches` 의 보고가 말하고, 계산이 끝나지 않으면(내보낼 솔리드가 없거나 거리가 정해지지 않으면) 작성기의 관례 값 1e-7 mm 를 적는다. 리더마다 이 값을 다르게 쓴다 — OCCT 는 1e-7 판과 같은 형상으로 읽되 «Poor result from projection» 경고를 붙일 수 있고, Rhino 는 파일을 열 때 이 값을 문서 허용오차로 삼는다. 어댑터가 내보내는 엔티티를 AP242 커널 형상 집합(cartesian_point·direction·line·circle·b_spline_curve/surface·plane·cylindrical_surface·advanced_face·closed_shell·manifold_solid_brep·surface_curve 등)으로 명시 타입 고정한다. 커버리지는 마일스톤을 따라 넓힌다(평면 → 곡선·곡면 → 트리밍). 목적은 기능만이 아니라 검증이다: 작은 형상에서 좌표계·방향·면 orientation 함정을 먼저 밟는다. 검증 뷰어는 두 겹이고 상보적이다 — **step-loupe**(step-io 리더 기반 웹 뷰어: report 가 드롭·고아·비표준 엔티티를 표시해 어댑터 버그를 구조적으로 잡는다; 자기 출력 되읽기의 GUI 판) + **FreeCAD 등 OCCT 기반 뷰어**(step-loupe 와 리더를 공유하지 않는 교차검증). 두 뷰어 모두 작성기와 코드를 나누지 않는다. STEP 출력은 오라클의 수송 계층이기도 하다. 작성기의 출력이 step-io 작성기의 것과 구조가 같다는 것은 brep-to-step 의 차등 테스트가 든다.

**M3 — 곡선 기하.** Arc, Cylinder, NurbsCurve/Surface 평가(The NURBS Book 기준 구현 + 수치 미분 대조 테스트). tess 출처 태그, tolerance 재계산(같은 모델, tol 3단). proptest.

**M4 — 면 위 작업.** 면의 스케치 프레임 위 돌출 + 불리언(pad·pocket 은 kit 의 조합). 여기까지 모델의 모든 점이 **구성된 점**이다(실측 tol 0). OCCT 오라클이 이 단계부터 돈다.

**M5 — 평면 솔리드 불리언.** 모든 면이 평면인 솔리드 간 fuse/cut/common. 평면–평면 교차는 닫힌 형식의 직선이다(SSI 행진·Newton·캐시 불필요). 꼭짓점은 평면 3장의 교차를 **좌표로 만들어 병합하지 않고 implicit point 로 두고**, 내/외·orientation 부호는 그 정의를 **indirect orient3d**(좌표를 만들지 않는다)에 넣어 정확히 판정한다 — 좌표를 만드는 순간의 오차·불일치를 원천 차단한다(Attene 2020). 세 엣지가 하나의 Vertex Handle 을 공유해 봉합 문제의 본체를 우회하고, 정점은 자기 정의를 든다(`Vertex::ThreePlane`). 이 세계에서 강건 불리언은 연구가 아니라 꼼꼼한 케이스워크(공면, 엣지–엣지 퇴화)다. 선·평면(다항식) 교차점은 indirect predicate 이론이 가장 깔끔하게 도는 영역이다. 커버리지 밖 입력은 `Rejected` 로 정직하게 거절한다. **발견된 점**의 경로·국소 tolerance·정의에서의 실현이 여기서 실전에 들어가고, 이 단계에서 «OCCT 없이 직동하는, 평면 위주 기계 부품을 STEP 으로 내보내는» 커널이 된다. 알고리즘 참고: Manifold(Apache-2.0 — 차용·번역 가능), Hoffmann 등 문헌. Truck 대비 벤치마크·정밀도 비교(전역 1e-6 폴리라인 vs 정점별 실측 tol + 닫힌 형식)는 수치로 보일 수 있는 차별점이다 — 공개 지표 후보.

**indirect predicates 는 자체 구현(clean-room)이다.** `nacre-predicates` 가 implicit point 표현과 indirect 술어를 자체 구현하고, 확장 산술 바닥은 `geometry-predicates`(elrnv, MIT/Apache)를 재사용한다 — 이 크레이트는 finished `orient3d` 뿐 아니라 Shewchuk expansion primitive(`two_product`·`two_sum`·`expansion_sum`·`scale_expansion_zeroelim` 등, `predicates` 모듈에 공개)를 노출하므로 expansion 산술을 자체 구현할 필요가 없다. primitive 는 `[lo, hi]` 순서이고, orient3d 부호 규약은 `det[a−d, b−d, c−d]` 다(golden 테스트로 고정). **라이선스 엄수** — Attene 참조 구현(LGPL)은 «돌려서 답 비교는 자유, 열어서 코드를 보는 건 MIT/Apache 소스만»: (a) **작성 중 소스 열람 금지**(LGPL C++ 를 열어 함수 대응·로직 흐름을 따라가면 2차적 저작물). 참고처는 논문(Attene 2020, arXiv 2105.09772; Cherchi 2020 mesh arrangements·2022 interactive booleans; Lévy 2024; Shewchuk 1997)과 `geometry-predicates` 소스로 한정한다. (b) **완성 후 실행 대조는 허용·권장**(dev-only, `tools/` 격리, 우리 크레이트에 링크 금지 — OCCT 오라클과 동일 논리).

**M6 — 평면+이차곡면 불리언.** 평면∩실린더(타원), 평면∩구(원), 평면∩원뿔 — 여전히 닫힌 형식이라 행진이 필요 없다. 실제 기계 부품 면의 대다수가 평면+실린더+원뿔이므로 여기까지로 실용 커버리지의 대부분을 확보한다. 원통 쌍 가운데 **면이 만나지 않는** 쌍(십자로 관통한 두 스터드처럼 축은 만나도 남은 면은 비킨 것, 큰 원통 속에 가로로 든 작은 원통)은 새 기하가 아니라 게이트의 명제 문제다 — 면이 방향 하나로 비켰거나 한 면이 다른 원통 안에 들었음을 증명하면 배열이 그대로 조립한다. 새 기하가 필요한 것은 면이 **만나는** 쌍(4차 곡선)뿐이다. 직교 같은 특수 배치부터 하나씩 열지 않는다 — 옆면 차트의 곡선 어휘와 교선이 배열에 들어오는 길을 먼저 넓힌다(`todo.md` 「원통–원통」).

**M6/M7 난이도 절벽.** M6(평면+이차곡면)는 닫힌 형식 교차라 SSI 의 지옥(위상 판정·시작점 검출·watertightness)을 **대부분 회피**한다 — 여기까지가 «확실히 되는 실용 커널»이다. M7(일반 NURBS)은 60년 미해결 문제로 **직접 진입**한다. 그래서 «**M6 까지가 실용적 종착점, M7 은 별도의 장기 연구 트랙**»이라는 선을 분명히 긋는다. M7 이 안 돼도 M6 까지로 실제 기계 부품 대부분을 커버한다.

**M7 — 일반 NURBS SSI 불리언(연구 구간, 계획).** `nacre-geom::intersect` 에 SSI 행진을 구현하고 일반 곡면쌍의 하이브리드 파이프라인을 완성한다: 출처태그 tess 에 강건 메시 불리언(exact predicates) → 조합 결정 추출 → 살아남은 면은 정확 곡면 유지, 신규 엣지는 국소 SSI 스냅백. OCCT 오라클과 상시 diff. 조합 결정이 메시 해상도에 의존할 수 있으므로 «위상 결과의 tolerance 불변»을 목표 불변식으로 삼고, 위반 사례는 실패 코퍼스에 쌓는다.

**M7 내/외 분류 = exact ray casting(1순위 후보).** 하이브리드 파이프라인에서 «이 patch 가 최종 솔리드 안인가 밖인가»를 분류하는 단계에 exact ray casting 을 쓴다. 메커니즘: 레이가 삼각형 **내부**를 지나면 삼각형 방향(정점 순서)으로 정확 판정하고, 꼭짓점·엣지·접선(coplanar) 같은 애매한 케이스는 레이를 **수치 섭동**해 항상 «내부 통과»로 되돌린다(Simulation of Simplicity, Edelsbrunner–Mücke 1990). 1순위인 이유: (a) **우리 계보와 정합** — 핵심이 orient 술어(채택) + 섭동(퇴화 구역을 **판정이 아니라 정책**으로 올린다는 방침의 정석)이라 새 수학을 들이지 않는다; (b) **검증됨** — Cherchi 2022(Interactive and Robust Mesh Booleans)가 핵심으로 채택, 수백만 삼각형·수백 입력 variadic 까지 테스트하며 GWN 류를 명시적으로 제쳤다; (c) **파이프라인에 그대로 꽂힌다** — 하이브리드는 이미 메시를 경유하므로 분류를 메시 단계에서 한다. 곡면 직접 판정이 필요 없어 GWN 의 trimmed NURBS 확장(최신·검증 진행 중)을 우회한다.

**오해 방지 셋.** (1) **Truck 방식이 아니다** — 표면적으로는 메시를 쓰나 정반대다: Truck 은 교차를 폴리라인 근사로 표현하고 그 근사가 **최종 결과**다(전역 1e-6). nacre 는 메시를 **분류용 임시 도구로만** 쓰고 exact 술어로 분류한 뒤 **정확 곡면으로 스냅백**한다. 차이는 «메시를 쓰느냐»가 아니라 «메시가 최종이냐 임시냐»다. (2) **indirect predicates 의 대체·확장이 아니다** — indirect predicates 는 M5 의 **점** 부호 판정(평면·국소·교차 계산 중), exact ray casting 은 M7 의 **patch** 내/외 분류(전역·분류 단계). 대상도 마일스톤도 다른 별개 부품이다. (3) **«어려운 케이스에만 켜는 정밀 모드»가 아니다** — 실현 사다리·적응 술어는 케이스 단위로 «어려우면 더 정밀하게»이고, ray casting 은 M7 에서 **상시 도는 분류 단계**다. «선택적»의 단위는 케이스가 아니라 **마일스톤**이다(M5·M6 은 닫힌 형식이라 불필요, M7 에서 켠다).

**한계 — ray casting 은 분류를 풀지 SSI 를 풀지 않는다.** 내/외 분류를 검증된 방법으로 채워도, 진짜 도박인 **SSI(곡면 교차 곡선 계산)**는 미해결로 남는다. ray casting 은 «조각을 어떻게 분류하나»만 정확히 할 뿐 «곡면 교차 곡선을 어떻게 정확히 뽑나»를 풀지 않는다. **«스냅백»과 «SSI»는 다른 일이다.** SSI(교차 곡선을 처음부터 **찾는** 것)가 미해결 문제이고, 스냅백(이미 찾은 교차 근처에서 정확 곡면으로 **되당기는** 것)은 그다음의 별개 작업이다. M7 의 위험은 스냅백이 아니라 «SSI 로 교차를 정확히 찾을 수 있는가»에 있다. SSI 가 실패하면 스냅백까지 갈 것 없이 그 연산은 `Rejected` 로 거절된다.

**SSI 는 검증된 해법이 없는 60년 미해결 문제다.** 일반 곡면 교차는 1960년대부터 연구됐으나 CAD 가 받아들일 만큼 강건·신뢰할 해가 없다 — **Parasolid·SISL·IRIT 등 최고 커널조차 특정 위상 케이스(작은 loop, 접선, cusp, 미세 자기교차)에서 실패**하고, 2025년 논문들도 watertightness 보장이 «여전히 challenging»이라 적는다. M7 의 «도박» 표기는 과장이 아니라 정확한 현실 인식이다. 아래는 **«채택할 검증된 기법»이 아니라 «M7 진입 시 읽을 프론티어 후보»**다(저자 벤치마크 주장일 뿐 독립 검증·프로덕션 채택 전이고, 미리 상세 검토하는 것은 «정보 없이 미리 설계» 함정이다): winding number + subdivision 시작점 검출(작은 loop·접선 branch 놓침 완화, 2026), Dixon matrix tracing(branch jumping/missing 을 root-solving 으로, Chen 2025), interval algebraic topology analysis(SSI 위상을 4D 대수계로 분류, Cheng 2023), lower-dimensional formulation gap control(Wang 2025).

M7 은 열린 연구다. M6 까지가 «확실히 되는» 영역이고, M7 은 이 커널의 존재 이유이자 도박이다. OCCT 의 역할은 전 구간에서 **오라클(dev 전용 시험관)뿐**이다 — 제품 경로에 OCCT 위임은 없다. OCCT 소스는 «상용급이 이 케이스를 어떻게 다루나» 열람용 참고서로만 쓴다: LGPL-2.1 이므로 번역·차용은 라이선스 오염이고, 무엇보다 OCCT 의 위상·tolerance 아키텍처가 딸려 들어와 nacre 설계와 충돌한다. 읽되 베끼지 않는다.

## 가지 말 것 — 재 보고 버린 길

이미 지어 보거나 재 보고 버린 접근이다. 다시 시도하려면 아래 숫자를 뒤집을 새 근거가 있어야 한다.
근거가 된 경위는 git 이력에 있다 — 태그 `pre-docs-restructure` 의 개발 로그.

### 성능

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
| 재생 비용 한계(`CHEAP_REPLAY_COST_CAP`)의 깊이에서 프레임 노드를 한 번만 세기(«그 펼침은 짧아 이력이 아니다») | 프레임 노드는 평면을 가리키고 재생은 그 평면의 사슬 전체를 걷는다 — 7° 와 이동 200 번 상자의 벽에 스케치한 프리즘은 깊이 1 로 세어 싼 예산이 201 단을 재생했다(모서리가 `Realized`). 노드 = 1 + `Model::plane_replay_depth` | `919fd282` |
| f64 필터 없는 정확 술어를 판정의 주 도로에 | `indirect_plane_side` 를 모션 없는 `orient3d` 에 쓰자 축정렬 fold 가 14배 느려졌고, `indirect_cmp_coord`(≈2.7µs, 호출마다 할당 — 필터 있는 인증 도로는 ≈1µs)는 지름길을 타는 평면이 늘자 fold 를 11~17% 늦췄다(한 스레드에선 빨랐다 — 병렬에서 할당자가 줄을 세운다). 둘 다 필터를 달자 기준보다 약 30% 빠르다. 필터는 정확 술어와 한 벌이다 | — |
| 공간분할(모델을 리전으로 잘라 잎에서 불리언) | 전부 초록인 채로 느렸다: 회전 fold 80 7.4→8.6s, 축정렬 1.7→2.4s. 정점을 용접하고 모서리를 공유하는 b-rep 에서는 분할면에 캡이 필수이고, 자르고 다시 꿰매는 일이 순비용이다(다각형 수프를 내는 EMBER 는 이 비용이 없다) | 되살릴 자리 `d3e9984`·`bfd01e6`·`12ec34d`·`9d69786` |
| 모션 사슬을 질의마다 걸어 접기 | 질의가 변환마다·평면마다 나오므로 깊이의 제곱이다: 기록된 사분각 회전 2,000번 0.20초 → 6.7초, 37° 4,200번 1.44초 → 2.14초. 노드가 태어날 때 부모의 접기에 자기 하나를 합성해 둔다(`chain_affines`) — 둘 다 기준선 | 실험, 되돌림 |
| C 로 된 전역 할당자(mimalloc) | 1.21× 이지만 `wasm32-unknown-unknown` 에서 빌드되지 않고(`wchar.h`), 제품 경로는 순수 Rust 다. 할당 36.5% 중 27포인트가 할당자 품질 — 브라우저용 순수 Rust 할당자는 별개의 미측정 항목 | — |
| 판정 결과 메모 캐시 | 1.14× 뿐이고(`WitnessPoint` 의 hp 캐시가 비싼 부분을 이미 흡수) 워커 공유 캐시가 병렬 5× 를 갉는다. 일을 «없애는» 레버는 병렬과 곱해지고 «저장하는» 레버는 나눠진다 | — |
| 서피스별 판정 한계를 transform 시점에 미리 계산해 저장 | O(n²): transform 4,200번에 1.58억 마디, 판정에 닿는 서피스는 마지막 6개뿐. «거절이 일보다 먼저 온다»는 속성을 깬다. 언제 계산하나는 누가 그 답을 읽나로 정한다 | 구현 후 되돌림 |
| 회전 좌표의 f64 산술을 FMA·보정 산술로 조이기 | 상한이 0 이다: 회전의 tol 기여를 정확히 0 으로 둬도 상승 131,510→131,392(−0.09%), 벽시계 무변화. 상승은 좌표 tol 이 만들지 않는다 | 코드 무변경 |
| n-ary fuse 를 큰 레버로 보기 | 누적 순서가 맞으면(큰 몸체 먼저) n-ary 는 1.04× 뿐이다. 순서 자체가 2.3×(OCCT 8.27→3.64s), nacre 도 1.48× | — |
| 「OCCT/FreeCAD 가 빠른 건 병렬 덕」이라는 전제 | OCCT 병렬은 CPU 4배로 1.05× 를 산다. 빠른 것은 단일 스레드 로직이고, 비교는 1스레드·같은 답에서만 한다. nacre 의 병렬은 실제로 4.25× 다(평면 클래스마다 독립) | — |
| 절약을 «개수»나 «식의 모양»으로 추정 | 없앤 연산이 싼 쪽이면 개수는 시간을 말하지 않는다: 방문 94% 헛일 → 1.3×, 반복 질문 92.8% → 1.14×, 쌍 3.3배 감소 → 1.4×, 「확실한 이득」이라던 단계 → 1.9% | — |
| rayon 스레드별 타이머를 합산해 비중 계산 | 합(3,719ms)이 벽시계(2,540ms)보다 컸다. 비중은 순차 빌드(`--no-default-features`)로 불리언 한 번을 끝에서 끝까지 회계한다 | — |
| 루프 불변 캐싱의 이득을 미시 벤치마크로 재기 | 최적화기가 재려던 것을 들어낸다: 41~47% 가 실종. 판단은 fold 의 같은 세션 A/B 로 | — |
| 간선의 진실 조각(이름 없는 평면의 직선 방향)을 push 마다 새로 실현하기 | 불리언은 결과의 간선을 매 걸음 전부 다시 push 한다 — 회전 핀 80 fold 에서 그런 push 45,888 에 새 담체 쌍은 걸음마다 열 남짓, fold 가 1.33 → 1.94 s 였다. 값은 담체 쌍 진실의 함수라 모델이 쌍마다 기억한다(`given_by_pair`) | 첫 시안, 되돌림(`784a35f2` 본문) |
| push 하는 쪽(ops)이 «topo 가 진실에서 답하나»를 push 전에 미리 판정하기 | 같은 이름 계산을 두 번 하고(사분각 회전 2000번 build 122 → 147 ms) 그 판정의 철자가 둘이 된다. topo 의 유도가 자기 길로 못 낼 때만 push 하는 쪽에 묻는다(간선은 `EdgeGiven` 닫힘, 곡면은 push 한 뒤 올리기) | 첫 시안, 되돌림(`784a35f2` 본문) |
| 링 감김(`loop_winding`)에서 링 평면이 고정하는 축을 사전식 비교에서 건너뛰기 | 링의 모든 노드가 그 평면 위라 그 축의 비교는 늘 0 이지만, 판정자가 공유 평면 항등식을 들면 더 주는 것이 없다 — perf.rs 순차 명령어 항등식만 −10.0%, 항등식 + 건너뛰기 −9.6%, 건너뛰기만 −3.9%. 0 은 묻는 쪽이 아니라 판정자가 안다(링 평면 밖의 공유 벽도 같은 0 을 낸다) | 실험, 되돌림 |
| `plane_pair_dir_sign` 에 «같은 사슬의 평행한 두 평면»(같은 핀의 마주 보는 두 벽) 항등식 | 참이고 불린다(기울인 핀 24 fold 의 «일치» 782 → 606, perf 에 남은 «일치»의 27%) — 그러나 질문마다 두면 이름 정수 외적이 perf.rs 순차 명령어를 +5.8%(축정렬 fold 80 1.02 → 1.79 s) 늘리고, 상승 직전에만 두면 기준선끼리의 흔들림(0.6%) 안이다. 남은 «일치»는 수는 많아도 시간을 거의 안 쓴다 | 실험, 되돌림 |
| 상대가 닿을 수 없다고 증명된 클래스의 추적 패스(pass A)까지 건너뛰기 | 순차는 더 빨랐다 — 회전 fold 80 9.66 → 4.78 s(2.02배; pass B 만이면 7.69 s), 허브먼저 ring 80 1.95배, 축정렬 fold 80 1.85배. 그러나 건전하지 않다: pass B 의 `PassThrough` 가 네 평면이 모이는 정점에 이름을 못 붙이면 배열로 물러나는데, 그 배열이 읽는 `splits[wc]` 가 비어 클래스의 면이 사라진다. 상대가 못 닿는 자리에 동시성이 있는 모델이 코퍼스에 없어 차등 검사(`ClassReuse::{Off, Proved}`)도 울릴 수 없었다. 클래스 단위 차등(`flip` 을 뒤집으면 커버리지 4건이 실패했다)은 추적하지 않은 클래스에 비교할 배열이 없어 불리언 전체 차등으로 올라갔다(census 130 행은 두 피처 조합에서 비트 동일). 남은 일이 169 클래스 중 27 에 몰려 병렬은 5% 나빠졌다(pass B 만일 때 1.98 s → 2.09 s, 효율은 pass B 만일 때 약 3.9배에서 2.3배로). pass A 는 모든 클래스에 돈다(`trace_result_faces`) | `059ef1c1`, 되돌림 `496ec741` |
| 판정마다 평면의 증인 정의를 다시 짓기(회전 증인을 복제하고 축정렬 증인을 다시 만들며 그 tol 을 120비트 astro-float 로 재기) | 회전 지느러미 25개를 허브에 fuse 하고 구멍 둘을 뚫는 불리언에서 실행 시간의 77%였다(백만 번 불렸고, 큰 절반이 tol 측정). 정의는 매번 같으므로 증인이 들고 술어는 빌린다(`plane_def`) — 같은 불리언이 2.16초 → 0.37초(native). 첫 걸음은 축정렬 증인의 tol 0 을 재지 않고 이름 붙이는 생성자 `Pt3::exact`(뒤의 `WitnessPoint::exact`)로 섰고, 그 생성자는 뒤에 은퇴했다(아래 「정밀도·판정」 표의 `WitnessPoint::exact` 행). 같은 작업에서 브라우저는 ~25초 → 0.4초였지만 그 절반 넘게는 디버그 wasm 을 내던 빌드 설정이라 이 레버의 수가 아니다 | `b34b3b41`·`5971a844`·`3622b65e`; 생성자 은퇴 `76449d63` |

### 정밀도·판정

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
| 이름이 다른 두 점의 좌표를 구간만으로 비교하기 | 구간은 «다르다»만 증명하고 «같다»는 못 한다 — 벽 `x = 9` 위 링의 모든 노드, 같은 벽 위의 상자 꼭짓점과 구멍의 관통점은 그 좌표가 정확히 같아 영원히 못 가르고, 기록된 벽 위의 37° 상자가 감김에서 거절됐다(피연산자를 바꾼 순서는 모든 다른 조각을 연 뒤에도). 두 점이 다 정확히 진술하는 좌표는 정확히 비교한다(`KeyCmp::both_stated`) | 스파이크 |
| 향 관계를 캐시 법선으로 읽기(내적·«첫 0 아닌 성분»의 부호 비교) | 내적은 크기 논증으로 건전해도 판정이 진실을 안 읽고, 첫 성분 비교는 크기 논증도 아니다 — 참값이 아주 작은 성분(10⁻¹⁷)의 부호를 정당한 반올림 상이 뒤집으면 `outward_fix` 가 반대를 답했다(단위 심기). census·스위트·perf 에서 진실과 불일치 0(첫 성분이 최대가 아닌 호출은 census 에 0). 향은 `sense`·`Orientation`·세계 이름의 향(`WorldName`)과 증인 정의(`normals_agree_judge`)가 답한다 | — |
| 간선 곡선의 종류(원·룰링)를 평면·원통 캐시의 허용오차로 고르기(`(n·d)² ≤ 1e-18·\|n\|²\|d\|²`) | 진실이 비스듬한데 캐시는 원을 답했다 — 무리수 축의 테스트 원통(캡을 실현한 f64 에서 들어올린 문)의 림 689: 진실의 캡은 축에 비스듬해 타원인데 간선 캐시는 원이었고, 판정 경로(`loops.rs`·`corners.rs`·원판 캡)가 그 종류를 읽었다. 축에서 `1e-10` 라디안 기운 평면을 심으면 원, 평행에서 `1e-10` 벗어난 평면은 룰링 선이 됐다. 종류는 정확한 관계(`AxisRelation`)가 정한다 | — |
| 교선 위 점을 비교할 축과 부호를 평면 캐시로 고르기 | 참 0 성분을 반올림이 0 에서 벗어나게 하면 그 축이 골라지고, 교선 위 점은 그 좌표가 모두 같아 서로 다른 두 관통점이 «같다»로 읽힌다 — 기울어진 원통(십진 직교 틀) 심기에서 Fuse/Cut/Common 두 순서 여섯이 전부 `CoincidentFeatures` 로 거절됐고, 유리수 방향(`stored_line_dir`)으로는 여섯 모두 해석값 부피다. census(호출 58,322)에는 그 인구가 없었다 | — |
| 접히는 사슬 아래의 참값 0 좌표를 접기 대신 일치 정밀도 규칙(과 싼 도로의 두 번째 단)으로 풀기 | 풀리기는 하지만 사분각 회전 2000번 빌드가 74 → 352ms(4.8배)였다 — 그 사슬은 노드마다 유리수 사상으로 접혀 있어 정점이 정확한 유리수로 읽히고(`chain_point_rat`), 그 도로는 오히려 빠르다(74 → 64ms). 규칙은 접히지 않는 사슬의 몫이다 | 실험, 되돌림 |
| 사슬이 안 접히는 평면의 캐시를 진실의 세 점을 f64 로 재생(`replay`)해 짓기 | 생산자의 f64 보다 참값에서 멀다 — 스위트 22,759 평면에서 256비트 참조에 생산자 값이 더 가까운 것 11,457, 재생이 더 가까운 것 2,439, 같은 것 8,863(census 출력 34행 추가 이동). 둘 다 반올림이다 | 실험, 되돌림 |
| 모션 노드를 기록할지를 캐시 좌표(꼭짓점·평면 원점·원통 축 원점)가 f64 에서 정확히 옮겨지는지로 고르기 | 두 정확 표현(옮겨 적기·사슬) 가운데 하나를 캐시가 골랐고, 표현을 읽는 규약이 그 선택을 물려받았다 — 같은 자리의 면이 `0.5` 로 옮기면 세계 진술, `1/3` 으로 옮기면 사슬이 되어 거울 뒤 스케치 프레임의 손방향이 갈렸고, 같은 벽이 `(name, motion)` intern 에서 두 핸들이 되어 병합이 고정밀 판정에 기댔다. census 이동 419 중 186, 스위트 약 2,765 번이 그 판정에서 갈렸다(바꾼 뒤 census 출력은 비트 동일). 진술이 `i128` 안에서 옮겨지는지만 묻는다 | — |
| 평면 클래스를 캐시(평면 캐시 계수의 비례 · 면 꼭짓점 좌표의 정확 `orient3d`)로 합치기 | 반올림된 상이 같은 서로 다른 두 평면이 한 클래스가 된다 — 심은 두 배치에서 `Common` 의 경사면이 남의 평면을 들었고, 닿지도 않는 두 몸의 Fuse 에서 한 몸의 벽이 다른 몸의 평면을 들었다(`validate` 무위반). census 코퍼스(쌍 30,009)에는 그 인구가 없었다(불일치 0). 병합은 세계 이름, 없으면 정의를 읽는다 | — |
| 정확 지름길의 자격을 캐시끼리의 대조로 주기(평면 캐시의 계수 ↔ 면 꼭짓점의 좌표 캐시) | 둘 다 반올림된 상이라 «반올림된 평면»을 인증한다 — `3·0.1 = 0.3` 인 모서리를 벽 `y = 3x` 밖으로 판정해(진실과 272번 다른 답) 유효한 교집합을 `ZeroLengthEdge` 로 거절했다. census 코퍼스(판정 186만)는 이 배치를 안 때렸다 — 반올림이 유리수 일치를 깨는 배치를 손으로 심어야 보인다. 지름길은 평면의 이름을 읽는다 | `b1d39469` |
| 잔차(점에서 담체 `f64` 평면까지의 거리)를 좌표의 경계로 읽어 자기접촉 체의 상자를 부풀리기 | 잔차는 거리 «하나»라 좌표가 참값에서 얼마나인지를 말하지 않는다 — 참값의 최근접 f64 까지를 못 덮은 seam 점이 스위트 3,415/9,494(다른 것 중), 그중 155 는 천 배 넘게 모자랐고 census 190. 체는 정확 판정 «앞에서» 후보를 버리므로 모자란 상자는 진짜 접촉을 조용히 버릴 수 있다. 체는 `Realized` 의 증명된 경계를 읽고, 경계가 없으면 버리지 않는다 | — |
| 「애매하면 사용자에게 묻는다」 | 회전 코퍼스에서 상승 2,490 중 2,260 이 판정 불가(연산당 수백 건), 가장 느슨한 한계가 2⁻²²⁹, 그리고 f64 출력은 「떨어져 있다」는 답을 표현하지 못한다. 선택지가 하나뿐인 질문은 질문이 아니다 — 근거를 붙여 보고한다 | — |
| 판정의 정밀도 사다리(256→512→1024 배증) | 256·512·1024·2048 비트에서 결과가 비트 동일 — 애매함은 정밀도 부족이 아니라 구조적 0 이다. 필요한 비트는 조건수 `C` 에서 계산해 한 번에 점프한다 | `cc37d86` · `f7a4fec` |
| 최근접 `f64` 실현 사다리를 256비트에서 멈추고 그 구간의 중점(충실 반올림)을 내기 | 반올림 경계에 256비트로 못 보는 만큼 가까운 값이 틀린 쪽으로 간다 — `sqrt_f64((1+2⁻⁵³)² + 2⁻⁴⁰⁰)` 이 `1.0` 을 냈다(최근접은 `1 + 2⁻⁵²`; 256비트 실현이 정확히 경계로 반올림돼 중점의 짝수 쪽). 정확 갈래가 경계·범위 끝이 될 수 있는 값을 먼저 거르고 0 은 반지름 0 으로 오면, 남는 값은 정해질 때까지 오르면 끝난다 — 얼마나 깊이 오르는지는 √ 계열이 분리 한계(`separation_cap`)로 증명하고, cos·sin 은 한계가 없어 비용 경보(`TRIG_LADDER_CAP`)를 둔다 | `c8ffa902` |
| √ 계열의 사다리에도 고정 상한(2¹⁵)을 경보로 두기 | 폭 제한 없는 입력이 정당하게 닿는다 — `sqrt_f64((1+2⁻⁵³)² + 2⁻⁴⁰⁰⁰⁰)` 의 근은 경계에서 약 `2⁻⁴⁰⁰⁰¹` 라 4만 비트가 필요했고 상한에서 패닉했다. 상한은 입력의 폭에서 증명한다(`separation_cap`) | `34e20e07` · `14bb6637` |
| 고정 판정 정밀도 | `C` 는 회전 1회당 1비트 자란다: 256비트에서 245회 회전부터 거절, 그것도 사용자가 대응할 수 없는 하류 증상의 이름으로 | `f7a4fec` |
| f64 `tol` 을 조건수로 읽어 임의정밀 시행 패스를 생략 | 2000/2000 에서 하회했다 | 코드 무변경 |
| double-double(`twofloat`/`qd`)를 고정밀 층으로 | 영점 근처 cos 오차 ~1.8e-16 — 부호 판정이 일어나는 곳이 영점 근처다. 천장도 약 1e-32 로 고정이다; `astro-float` 는 160비트에서 cos/sin 오차 약 1e-58 이고 정밀도가 다이얼이다 | `79172670` |
| 접선형 회전 tol(`da·\|y\|`, `da·\|x\|`) | 축 근처에서 건전하지 않다(cos/sin 이 독립적으로 반올림돼 반경 방향 오차 `\|x·δc − y·δs\|` 가 생기고, 좌표 하나가 0 근처면 접선형은 그것을 과소평가한다). 2D 실험의 랜덤 10,000 에서 220회 상한이 무너졌고 좌표 혼합형은 0회였다. `frame3` 의 테스트 `the_tangential_tol_is_unsound_and_the_mixing_one_is_not` 가 `!tangential_sound` 를 단언한다 | `79172670` |
| libm 의 sin/cos 를 단일 진실원으로 | 같은 프로세스·같은 각에서 호출부마다 1 ulp 다르다(리터럴 각은 LLVM 이 컴파일 타임에 접는다). 정확도 계약도 없다(잰 오차 `< 5ε`). 고정밀 실현을 정확 반올림하는 쪽이 느리지도 않았다 — 회전 fold 80 에서 libm 731 ms, 정확 반올림 715 ms | `009490c7` |
| SoS(ε-nudge)로 퇴화를 없애기 | 첫 벽(64건)은 사라지지만 exact 삼중 기판의 둘째 벽(69건)으로 옮겨갈 뿐, 완주 0. 퇴화는 이름 붙여 거절하거나 엔진이 구조로 답한다 | 스파이크, 프로덕션 무변경 |
| 도달 가능한 실패를 `debug_assert` 로 막기 | release 에서 꺼진다 — 체인 접촉 절단이 그 단언을 반증했다. 도달 가능하면 이름 있는 거절이다 | — |
| 향이 걸린 값(평면 캐시의 법선)을 판정 증인 삼각형으로 펼치기 | 이름 있는 넓은 만남 평면의 증인은 그 프레임의 탐침이라 회전 방향이 진술 점의 것이 아니다 — `sense` 를 곱하자 캐시가 진실과 반대를 봤다(`a_datum_through_wide_meets_keeps_its_name` 에서 감김 `debug_assert` 가 잡음). 법선은 모션 전 이름 법선을 사슬로 운반하고, 증인은 계수(만남·판정)에만 쓴다 | 첫 시안(`e9a6e3ab` 본문) |
| 오차 한계를 확률 전파(RSS)로 | 실제 오차에 더 가깝지만(최악 상한보다 약 100배 작다는 어림, 측정 아님) 보장이 아니다 — 드물게 부호를 뒤집어 조용히 틀린다. 한계는 최악 상한(삼각부등식)이고, 커짐은 유리수 구간의 tol 0 과 정밀도 상승으로 다룬다 | `5bf0ea36`(논증) |
| 증인점에 «좌표가 f64 로 정확하다»를 전제로 받는 생성자(`WitnessPoint::exact`) | 그 전제는 어떤 호출자도 검사할 수 없고, 반올림된 캐시를 건네면 다른 점을 tol 0 으로 이름 붙인다 — 그렇게 깬 호출처가 있었다. 표현 가능한 점에서는 `at` 이 같은 tol 0 을 말하고 그 밖에서도 정직하다 | `76449d63` |

### 타입·진실

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
| 모델이 인접 인덱스를 캐시로 들기(`Model::adj` + `rebuild_adjacency`) | 읽는 쪽이 topo 자기 시험뿐이었다 — `validate` 는 제 `Adjacency` 를, 조립과 tess 는 제 간선 표를 지었다. 시험과 kit 의 재건 호출 1,500여 곳이 시간만 썼다 | `dcb966bc` |
| 원과 타원을 간선 곡선의 두 변종(`Curve::Circle` + `Curve::Ellipse`)으로 | 짓지 않았다 — 논증이다. 둘은 한 명제(평면과 원통이 만나는 곡선)이고, 소비자(tess·props·validate·STEP·ops)마다 점·활꼴·모멘트·감김·상자를 두 철자로 쓰게 된다. 원은 기울기 없는 타원이라 한 타입(`nacre_geom::Ellipse`)이 공식을 메서드로 들고, 원의 갈래는 메서드 안 한 곳에서 원의 식을 연산 그대로 돈다(원의 갈래를 들인 그 커밋에서 census 메시 343 비트 동일, 잰 것). 종류는 정확한 축 관계가 정한다 | `3ab2aad9`·`746e082f` |
| `Vertex` 를 세 평면 단일형으로 | 원통 seam 정점은 «원통 ∩ 캡» 원 위의 한 점이라 세 평면이 없다. 슬롯을 중복으로 채우면 판정의 `D=0` 이 조용한 0 이 된다 — 변종별 enum | — |
| 향을 든 진실을 interning 열쇠로(`Verbatim(Surface)` — `Surface`·`PlanePoints` 에 `Eq, Hash` 를 붙여 진술 그대로를 열쇠로) | 평면의 진실은 향(`sense`)을 들고 평면의 정체에는 방향이 없다 — 같은 진술이 반대 향으로 오면 한 핸들 + `flipped` 여야 하는데 두 핸들이 된다(스위트에서 이름 없는 `Through` 적중 4 중 1). 열쇠는 향을 뺀 진술이다(`SurfaceKey`) | — |
| 정점에 정확성 태그(`Constructed`/`Discovered`)를 달기 | 그 태그는 정확성을 뜻한 적이 없다 — 정확 경로로 만든 정점의 좌표도 반올림이다. 구분은 실현 캐시가 든다 | `76a07b2` |
| 평면의 진실을 계수로 | 계수는 두 점 차의 곱이라 `i128` 적합 25.8%(점은 100%). 분모를 먼저 날려도(lcm 자체가 104~106비트) 중간 타입만 넓혀도(최종 계수 중앙 159비트) 안 든다. 진실은 점 셋이다 | — |
| f64 계수를 `Rat` 으로 들어올리기 | 반올림이 진실에 구워지고, 같은 벽이 두 경로로 두 평면 클래스가 된다 | — |
| 실현한 프레임 기저(캐시 법선의 축·캐시나 정확 투영의 `to_f64` 원점)를 `from_decimal` 로 되들어올려 도로와 점을 정하기 | 짧은 십진이 아닌 좌표는 그 f64 가 찍는 십진으로 돌아온다 — `1/3` 들어 올린 면의 pad 먼 캡이 `43/30` 대신 `28666666666666667/2·10¹⁶` 위에 섰고(6/6), 같은 면의 오프셋도 그랬다. 정규화가 필요한 축은 `0.6000000000000001` 로 돌아와 유리수 3-4-5 면이 프레임 노드 도로로 빠지고, pad 윗면이 같은 평면의 extrude 캡과 다른 핸들이 됐다. 프레임의 세계 기저는 이름에서 묻는다(`exact_frame`) | — |
| `n·n` 같은 파생값을 저장 | 필요하면 정의에서 다시 실현한다. 캐시가 안 만들어진다고 진실을 거절하던 것이 병이었다 | — |
| 모션을 정점에 달기 | 면이 든다 — 정점은 세 면의 교점이라 따라온다(정점에 달면 3중 중복) | — |
| 모션 노드를 합쳐서 사슬을 접기 | 노드 핸들이 `SurfaceKey` 의 세 열쇠 모두에 든 interning 열쇠라 합치면 동일성 판정이 사라진다. 접는 것은 아레나가 아니라 «읽기»다 | — |
| 코퍼스 최대값을 상한으로 쓰기 | 폭은 타입에서 유도한다(이름 ~2²²⁹¹) — 실측 최대는 표본이다 | — |
| 정준 원점이 넘치면 `points[0]` 으로 폴백 | 스케치의 (0,0) 이 오버플로 여부에 따라 달라진다 | — |
| 좁은 이름을 별도 곁표로 | ~99% 가 중복 저장 — `PlaneName = Narrow \| Wide` 한 enum + `narrow()` 투영 | — |
| wide 계수를 f64 `Expansion` 에 | 조각도 f64 라 지수 상한 ~2¹⁰²³, wide 이름은 ~2²²⁹¹. 그릇은 BigInt 정수 산술이다 | — |
| 53비트를 넘는 좁은 이름을 f64 조각 여럿으로 `Expansion` 에 | 필터를 못 넘는 꼬리가 BigInt 보다 느리다 — release, 52비트 입력에서 1조각 `Expansion` 2,880/1,361 ns(cmp/side) 대 BigInt 891/548 ns, 조각이 늘면 더 느리다. 그 질문의 절반 가까이(스위트 cmp 935k/1.98M)가 정확한 0이고 0은 어떤 필터도 못 넘으므로, 꼬리의 값이 곧 도로의 값이다. 조사한 커널(Geogram·Attene·CGAL)도 넓은 입력을 조각으로 정적 필터에 넣지 않고 구간 → 정확으로 간다 | — |
| BigInt 쌍둥이를 필터 없이 | 모든 평면에 쓰면 상승 도로와 같은 값이다(축정렬 작은 불리언 200개: BigInt 436 ms · 지름길 끔 449 ms · 지름길 219 ms), o3d 는 호출당 870 ns 대 toleranced 61 ns. 필터가 먼저다 | — |
| 판정층에 `Through` 전용 평면 기계를 짓기 | 판정 표의 계약은 «평면 위 정확한 세 점»이고 유도된 증인 삼각형(진술의 점, 또는 판정 프레임의 probe)이 이미 준다 — 새 기계 0 으로 열렸다 | — |
| 관통점을 두 평면의 교선이 원통과 만나는 점으로 실현하되 두 평면을 원통 틀에서 함께 읽을 때만 | 벽(축에 평행한 평면)이 있는 관통점은 벽만 틀에서 읽고 다른 평면은 어느 사슬에 있든 세계에서 만나는데, 그 도달을 잃는다(독립 검토). 교선을 축 방향으로 쓴 벽의 룰링 위 점이라 하면 벽이 곧 그 쓸린 벽이고, 둘 다 벽이 아닐 때만 두 평면을 읽는다 — 잰 것: 벽 없음·못 읽음 인구 0, 벽 없음·읽음 스위트 111·kit 1 이 새로 실현 | `c1830a53` |
| 뒤집힌 `Canonical` 프레임의 기준 방향(`ẑ×n`)을 뒤집기 **전** 정준 계수에서 유도하기 | `û` 가 이름의 정규화 부호(첫 0 아닌 성분 양수 — «같은 평면인가»의 열쇠 관례)를 따라, 좁은 도로만 `û` 를 두고 넓은·이름 없는 도로는 `û` 를 돌렸다 — `z=0` 을 `−z` 로 보는 프레임에서 `[0,1.5]²` 가 좁은 도로로 `[0,1.5]×[−1.5,0]`, 이름 없는 도로로 `[−1.5,0]×[0,1.5]`. 이름 없는 도로는 계수의 부호를 증명하지 못해 그 규약에 철자가 없다. 스위트의 그런 프레임 ~156 이 잠금 없이 움직였다(census 0) | `6b80011d` |
| 평면의 방향을 점 순서에서 뽑기 | 점 순서는 향을 말하지 않는다: 씨앗·datum·원통 밑캡은 관례로 점과 반대인 향을 들고(비씨앗 `Known` 342, 씨앗 전부), `Stated` datum 은 순서를 프레임 배치(`ref_dir`)로 읽고, `Through` 삼중은 정렬된 열쇠다. 향은 진실의 한 비트(`Surface::Plane.sense`)가 든다. 이름은 「어디」에 답하지 「어느 쪽」엔 못 답한다 | — |
| 평면 캐시를 «수선의 발 앵커 + 정준 행»으로 | census 32행이 움직이고 8건이 넘친다(`n·n` 이 계수를 제곱한다). 살아남은 것은 «진실 삼중의 첫 점» 앵커뿐 | `8b53a2c` |
| `rotated` 플래그를 풀어 회전 클래스를 원통 관문에 들이기 | 그 플래그는 정확 지름길이 평면의 이름을 세계 진술로 읽는 허가다 — 모션 있는 평면의 이름은 모션 전 프레임을 말한다. 풀었을 때는 클래스 병합이 면 꼭짓점 캐시의 삼각형 판정으로 가서 공유 벽이 두 몸으로 갈렸다(그 판정은 은퇴했다 — 병합은 세계 이름을 읽는다) | 되돌림(`planes` 주석) |
| 접미사만 다른 형제 함수를 «중복»으로 세어 합치기 | 일곱 갈래를 열어 보니 진짜 중복은 0 이었다. `det3`·`orient3d_iv` 계열(정의 12개)·`three_planes_rat`/`_big` 은 층마다 숫자 타입(`Expansion`·`f64`·`Bounded`·`HpBounded`·`BigInt`)이 다른 **정밀도 사다리**이고, 나머지는 트레이트 기본 구현+impl, 인자가 다른 별개 함수, doc 에 사유가 있는 한 줄 편의 문이었다. 세기 전에 **인자 타입**을 먼저 본다 | — |
| 원통 옆면을 seam 엣지(A — 세로 seam 간선 `[c, c]` 를 두 번 쓰는 한 루프)로 짓기 · 원을 극점에서 단조 조각으로 자르기(CGAL 의 표현) | 각도의 시작점(`ref_dir`)이 위상에 새겨진다: 잘린 림마다 θ = 0 정점과 호 쪼개기(census 154·스위트 1,980), seam 에 걸린 구멍을 잇는 slit 걸음, θ = 0 에서만 보는 집힘 그물, self-pair 담체 규칙이 붙고, 같은 모양이 seam 의 자리에 따라 다른 위상이 된다. 단조 조각도 같은 이유다 — 자르는 극점이 좌표틀에 따라 위상에 새겨진다(이쪽은 잰 것이 아니라 논증이다). 옆면은 림 두 루프로 닫고(「seam 엣지 방식」 B), 닫힌 간선 `[v, v]` 는 한 바퀴로 읽는다 | `e60ab2b9`, `497b54b5`; 단조 조각 `278484eb` |
| 면의 스케치 원점을 면에서 읽기(외곽 꼭짓점의 평균, 면 영역의 면적중심) | 꼭짓점 평균은 직선 도중에 꼭짓점이 하나 늘면 움직인다 — 면의 모양은 그대로인데 그 위에 pad 한 보스가 다른 자리에 놓인다. 면적중심은 이산화에는 무관하지만 면의 `f64` 꼭짓점으로 계산한 원점을 `Rat::from_decimal` 로 올리므로 반올림된 캐시가 진실이 된다 — 같은 footprint 를 두 번 pad 하면 둘째 프로파일이 한 ulp 옆에 놓여 넓이 `2.2e-16` 의 면이 남았고 `validate` 는 알리지 않았다(반올림을 줄인 고침도 새는 길을 닫지 못했다). 둘 다 같은 평면 위 다른 외곽선이면 원점이 달라진다. 원점은 평면만의 성질로 유도한다: 세계 방정식이 진술된 평면(`Model::world_plane_name`)은 세계 원점의 정사영, 그 밖의 평면은 사슬로 나르는 제 `Canonical` 프레임의 원점(`the_sketch_origin_does_not_depend_on_the_outline_at_all`) | `0615e95a`·`00ec61e`·`c652ef4` |
| 면의 스케치 축(`+u`)을 `any_perpendicular` 로(법선과 가장 덜 정렬된 세계 축과 외적) | 고르는 축이 법선을 따라 바뀐다 — 상자 뚜껑이 `u = −ŷ, v = +x̂` 로 `SketchPlane::world_xy` 와 90° 어긋나 한 평면이 두 철자를 가졌다. 불연속도 두 극이 아니라 가장 작은 두 성분이 같은 호 전체에 놓인다. 세계 방정식이 진술된 평면(`Model::world_plane_name`)의 축은 arbitrary-axis 규약(`ẑ × n`, 수직이면 `ŷ × n` — `frame_axes`, 정확히는 `plane_frame_default`)이고, 그 밖의 평면은 제 `Canonical` 프레임을 사슬로 나른다. `any_perpendicular` 는 seam·STEP `ref_dir` 처럼 «아무 수직»이면 되는 자리의 것이다 | `20dae82` |

### 정리

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
| 긴 함수를 길이로 쪼개기 | 쪼갤 근거는 «안의 한 단계를 다른 호출자가 이름으로 부를 만한가»다. 300줄이 넘는 열셋 가운데 그런 호출자를 든 단계는 **하나**였고(감사가 손으로 다시 짓던 방향 분할 → `direction_partition`), 나머지는 블록이 바깥에서 받아야 할 지역변수가 1~32개다 | `abaa95d`·`9ab9ac3` |
| 글자 스캔 `comparators.rs` 를 clippy `disallowed-methods` 로 | 비교를 받는 정렬은 옮겨진다(`slice::sort_by` 경로가 풀리고, clippy 가 찾은 9곳이 허용 표의 정렬 9와 같았다). 그러나 같은 규칙의 다른 반쪽 — 손으로 쓴 `Ord`·`PartialOrd`(`nacre-store` 2) — 은 clippy 에 금지 수단이 없어 한 규칙이 두 장치로 갈린다. 글자 스캔이 두 반쪽을 한 표로 본다 | `e42511b6` |
| 글자 스캔 `silencers.rs` 를 `clippy::allow_attributes_without_reason` 으로 | lint 하나의 범위로 좁혀지지 않는다 — 규칙은 `allow(dead_code)` 6곳 위의 이유인데, 그 lint 는 제품의 `allow` 68개 전부(`too_many_arguments` 33·`large_enum_variant` 11…)에 `reason` 을 요구한다 | `e42511b6` |
| 주석 비율(36~45%)을 줄이기 | 주석의 대부분은 사용자에게 가는 doc 15,182줄이고, 본문 주석의 가장 긴 덩어리 열 개(294줄)에서 되풀이는 **14줄(5%)** 뿐이었다 — 나머지는 불변식·측정치·거절한 대안·위험이다. 판단 규칙은 `overview.md` 문서 규칙 7 | `82f7b6c`·`56e6756` |

### 엔진·연산

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
| 정확 비교의 실패를 플래그로 세우고 `Equal` 로 읽어 `sort_by` 에 태우기 | 그런 비교는 전순서가 아니고, 표준 정렬(Rust 1.81~)은 길이가 20을 넘으면 그것을 알아채 패닉할 수 있다 — 길이 40·원소 열에 하나가 끼는 비교마다 실패에서 2000번 중 1896번, 배열의 분할점 정렬에서는 실패하는 점 하나로. 이름 붙은 거절이 아니라 충돌이다. 실패할 수 있는 비교는 `combinatorics::try_sort_by`(실패에서 멈추고 슬라이스를 순열로 남긴다)로 정렬하고, `tests/instruments/comparators.rs` 가 비교를 받는 정렬과 손으로 쓴 `Ord` 를 허용 목록에 묶는다 | — |
| 자기접촉 판정의 후보 평면을 두 끝점 이름 삼중의 교집합으로 고르기 | 이름은 그 점을 지나는 평면 셋만 담는다. 간선 직선을 평면 셋이 품는 다발(기운 막대가 상자 윗면에 간선으로 얹힘)에서는 이름이 그중 둘만 담고, 어느 둘인지는 클래스 번호 — 피연산자 순서 — 가 정해 한 순서에서 그 면을 묻지 않았다(`Ok(1)`, `validate` 깨끗). 두 끝이 면의 평면 위에 정확히 놓이는지를 묻는다 | — |
| 림 위에 놓인 꼭짓점의 접촉을 옆면 점 팔(꼭짓점 × 경계하지 않는 옆면의 닫힌 영역)로 잡기 | 배열이 그 꼭짓점을 림의 두 면(캡이나 판 윗면, 그리고 옆면)에 다 용접해 «경계하지 않는 면»이 없다 — 상자 모퉁이를 보스 림에 얹은 결과의 옆면 윗림이 그 점에서 갈린 닫힌 호였다. 팔을 지어 본 스파이크에서 스위트·census·무시 스윕·kit 의 새 답 0. 그 점은 두 부채가 한 꼭짓점을 나누는 집힘이고, 링크의 마디를 모서리 끝으로 읽는 집힘 검사(`nacre_topo::nonmanifold_vertices`)가 답한다 | 스파이크 |
| 원통이 깎은 면을 가로지르는 현을 좌표 도로의 쌍둥이(`point_in_mixed_ring` 의 구간 판 — 선 위 사건을 ℚ(√c) 차트 좌표로 모아 교대를 읽기)로 묻기 | 짓지 않았다 — 논증이다. 선이 두 원통의 호를 지나면 두 근호의 사건을 서로 비교해야 하는데 `QuadVal` 산술은 근호가 어긋나면 멈추고(디버그 빌드에서는 패닉), 정렬·run·교대 규칙을 한 벌 더 쓰게 된다. `segment_meets_face` 는 이름 도로라 추적기와 같은 원시(`arc_crossings`, 전체 원통 표)로 읽고, 두 근호의 비교는 `order_on` 이 이미 정확히 한다 | — |
| 실패하는 배치마다 경로를 덧대기(케이스워크) | 배관(detector·생존표)이 다음 엔진으로 전달되지 않고 경우의 수가 폭발한다. 평면당 셀 복합체 하나(winding 단일 배열)를 일반화해서 닫는다 | — |
| 원통 차트의 선을 «축을 따르지 않는 모든 평면 클래스»로 | 칸은 양 끝 선의 자취(딱지)로 읽히는데, 자취 없는 단면이 선이 된다: 게이트가 옆면과 떨어졌다고 통과시킨 비스듬한 단면이 필렛 사분원 밖에서 span 의 띠를 지나고, 그런 두 선 사이의 칸은 span 이 «있음»이라 하는데 말할 변이 없어 `CylinderStagesDisagree`(`a_chord_across_a_face_trimmed_by_a_cylinder_is_refused`). 선은 배열이 그린 자취다 | — |
| 평면 DCEL(`walk_cells`·`nest_cells`·`label_cells`)을 원통 차트에 이식해 셀 기계를 하나로 | 단계는 같아도 도구가 같은 단계가 **0** 이다: 바깥/구멍은 감김수 `#(−1) == 성분 수` 대 «가장 낮은 칸의 밑변» — 차트를 감는 띠에는 평면 감김수가 없다, 딱지는 무한 뿌리에서 전파 대 이웃 클래스에서 읽기 — 환면에는 무한 영역이 없다. 한 코너에서 떠나는 방향의 순서도 출처가 다르다 — 평면은 정점의 각도 술어, 차트는 칸을 세우려고 이미 물은 slab 의 높이 순서(`Chart::slabs`)를 읽고 새 술어를 묻지 않는다 | — |
| 볼록 전용 fast path | 「같은 답을 더 빨리」가 거짓이었다: 일반 경로가 답하는 입력을 거절하고 진리를 f64 로 정한다 | — |
| 구멍을 외곽에 브리지로 꿰맨 뒤 귀 자르기 | 브리지가 정점을 반복시켜 two-ears 정리가 깨진다 — 무작위 구멍 면 4,000 표본의 27.6% 실패, 자기교차 후에도 삼각형을 뱉는 조용한 오답의 여지. 스윕으로 간다 | — |
| 한 점을 두 번 지나며 방문마다 제 재료 부채꼴만 감싸는 고리(꼬집힌 면)를 그 점에서 두 고리로 잘라 그리기 | 원판에 판 홈의 끝이 테두리에 닿으면, 잘라도 «구멍이 외곽에 정점으로 닿는» 같은 모양이 남는다 — 스윕은 여전히 한 점의 반직선 넷을 순서 지어야 하고, earcut 도 결국 대각선을 하나 더 긋는다. 두 방문을 쌍둥이로 두고 볼록한 부채꼴만 자기 간선을 따라 미끄러뜨린다(`polygon::pinch`). 두 방문이 서로의 반직선을 품는 두 엽은 잘린 조각이 서로 떨어진 고리라 자르는 것이 답이다 | — |
| 접선 접촉에서 집힌 면을 b-rep 에서 쪼개기 | 집힌 것은 «면»이지 «표면»이 아니다 — 접점 둘레의 링크가 원 하나라 표면은 2-다양체이고 결과는 유효한 솔리드다(OCCT 도 같은 합집합에서 부피·면적·면 수가 같고 집힌 면을 쪼개지 않는다). 쪼개기는 유효한 솔리드의 b-rep 을 바꾸는 churn 이다. 집힌 면은 메시가 그린다 — 공유 간선 다리와 스윕의 기호적 순서(`polygon::sos`, 「Tessellation 층」) | `3cf20151` |
| 곡면 tess 에서 예산을 넘는 간선을 차트 중점에서 쪼개며 라운드 돌기 | 발산한다(179 → 386 → 856 → 1952 → 4412). 이탈이 Δv 에만 의존하고 영역이 이방적이라 어디서 쪼개도 새 간선이 비슷한 Δv 를 갖는다 — v 방향의 띠가 필요하고 그것은 점을 미리 대야 얻어진다 | — |
| 곡면 tess 에서 지름원 encroachment 로 침범 후보를 버리기 | 직선 경계는 폴리라인이 두 점뿐이라 길이 4 솔기의 지름원이 한 면의 후보 360개 전부를 삼킨다. 자리 없는 후보는 삽입 탐색이 이미 거절한다 | — |
| 원통 차트의 격자 후보를 «경계 어디서든 칸의 `1e-6` 배 안이면 버린다» 한 문턱으로 거르기 | 결과는 같았다(지느러미 0, census `c` 행 그대로). 그러나 두 원인 가운데 높이 쪽은 **한 수(원의 자리)를 두 곳에서 계산**한 것이라 한 곳으로 모아야 고쳐진다 — 문턱 하나가 그것을 덮으면 차트의 두 철자가 남는다. 문턱은 기하의 우연(룰링 각 = 격자 칸)에만 쓴다 | — |
| 평면 캡의 경계를 «원호 조각의 볼록 껍질이 서로 떨어질 때까지» 미리 나누기(Li 외, JCDE 2026) | 짓지 않았다 — 논증이다. 나중에 점을 더 넣어도 안 엇갈린다는 보장은 좋지만, 껍질 검사를 f64 로 하면 정확한 접함과 ulp 겹침을 못 가른다: 구멍이 직선 간선에 정확히 접하는 다리 인구(`bridge_shared_edges`)는 껍질이 닿은 채라 자리마다 상한까지 나뉠 것이고(막으려면 문턱이 든다), 그려지던 메시가 바뀐다. 스윕의 정확한 교차 검사로 엇갈림을 찾아 그 옆면에 열을 청하면 그려지던 메시는 하나도 안 움직인다(잰 것 — `crossing_chords`) | — |
| 엇갈린 현을 한 바퀴에 하나씩 나누기 | 회전한 seam 의 거의 동심 구멍(둘레 전체가 근접)에서 23분 넘게 안 끝났다 — 엇갈림 수백 자리 × 자리마다 여러 번 × 바퀴마다 면 전체의 제곱 검사. 한 바퀴에 전부 나누니 바퀴 1~8 | — |
| 엇갈린 현 나누기를 다리 사전 패스와 별개의 패스로 | 두 패스가 같은 폴리라인에 같은 `meets`(면 전체의 제곱 검사)를 묻는다 — 따로 돌리자 census 메시 시간의 7.7%(630.6 ms 중 48.7 ms), census 벽시계 +2.5~3.8%(같은 세션 A/B 두 번). 마지막 바퀴의 검사가 다리의 접촉 후보를 함께 읽으니 잡음 안 | — |
| 바퀴의 재검사를 «θ 집합이 자란 곡선이 있는 면만»으로 고르기 | 짓지 않았다 — 논증이다. 경계가 앞 바퀴와 같은가를 의존(곡선의 집합, 띠 생성선, 다리 삽입)으로 추론하는 둘째 철자이고, 의존 하나를 빠뜨리면 낡은 답이 조용히 남는다. 검사는 차트의 점과 링만 읽으므로 그 차트를 비트로 비교하는 것이 명제 그 자체다(`polygon::Asked`) — 원판 100 곁 9바퀴 모델에서 바퀴마다 새로 묻기 711 ms, 같은 차트의 답 쓰기 457 ms | — |
| 엇갈린 현의 중점을 두 끝 표본의 좌표(`Ellipse::angle_of`)에서, 또는 출처의 풀린 자리(`θ + τ`)에서 읽기 | 좌표에서 읽으면 구멍의 두 림이 제 중심에서 읽은 θ 가 1~2 ulp 갈려 `1e-16` 간격의 두 열이 길이 0 에 가까운 현을 남긴다 — 반지름 10 원판에 구멍을 뚫거나 기둥을 합친 270 개(축 넷)에서 33 이 메시 거절(`SelfTouchingBoundary` 25·`OverBudget` 8), 모두 영 성분이 있는 3-4-5 축. 풀린 자리로 읽으면 원점에서는 0 이지만 몸이 `(1000, −2000, 3000)` 에 있으면 57 — 시작 각이 0 바로 아래로 읽힌 곡선은 모든 열을 한 바퀴 위에 세운다. 열의 수(`[0, τ)`)에서 읽으면 세 자리 모두 0 | 스파이크 |
| 타원 림이 낀 키 큰 옆면 띠의 예산 밖 대각선을 타원 높이의 극값에 격자 줄을 놓아 고치기 | 비스듬히 자른 봉에서 Delaunay 뒤집기가 표본 하나를 건너뛰는 4° 대각선(예산 2°)을 골랐고, 극값 두 줄을 더해도 같은 `OverBudget` 이었다(잰 것). 띠를 등방으로 채우는 촘촘한 줄은 `insert_interior` 가 점 × 삼각형이라 비싸다. 예산이 계약이고 Delaunay 는 취향이라, 뒤집기는 예산 밖 원소의 수를 늘리지 않고 그 수를 줄이는 뒤집기는 `incircle` 과 무관하게 한다(`delaunay::refine` 의 `within`) — 간선만 묻던 첫 판에서 census 메시 24/343 이 대각선만 바뀌고 census 명령어 −2.5% | `daacf795` |
| 면의 기울기를 뒤집기의 거부권으로(뒤집은 뒤 대각선과 두 면이 다 예산 안일 때만 뒤집기), 또는 밖인 수가 줄 때만 뒤집고 같은 수에서는 둘 다 0 일 때만 Delaunay | 타원 림을 제 회전율로 나눈 봉 끝의 귀는 사라지지만 거부권은 스위트 24 시험, «줄 때만»은 23 시험에서 넓이 0 이거나 2°의 배수(4°–70°)로 기운 면을 남겼다 — 단조 분해가 낸 기운 면을 치우는 길이 밖인 수가 같은 뒤집기를 지나간다. 같은 수에서도 `incircle` 을 따르고 면을 반올림 한계와 함께 판정하면(아래 행) 0 | 스파이크 |
| 면을 반올림 한계 없이 판정하기 | 생성선에서 몇 ulp 떨어진 격자점 둘이 생성선 위 경계 점 둘과 실 삼각형을 이루면 그 외적은 반올림뿐이라 법선이 아무 방향이다 — 판정이 동전 던지기가 되어 수를 세는 뒤집기가 멈추고 12°·30° 기운 면이 남았다(쐐기 시험 둘; 그 판정을 차트의 되돌린 점에서 하면 9). 법선이 돌았을 수 있는 각을 더하고 정하지 못하면 밖(`facet_in_budget`) | 스파이크 |
| 룰링이 끊긴 점의 격자 높이를 그 점의 축 좌표에서 읽기 | 한 평면 위의 두 끊김(상자가 두 룰링을 같은 높이에서 끊는다)이 37°·23° 기운 축에서 `1.0`·`0.9999999999999999` 로 한 ulp 갈려 두 줄이 서고 옆면에 기운 삼각형이 230 생긴다(평판 0). 끊는 평면이 원통에서 자를 원의 자리로 읽으면(`ruling_break_place`) 두 끊김과 그 평면의 원이 한 `f64` 다 | — |
| 옆면의 곡선을 면의 격자 열(`τk/N`)에만 표본하기(꼭짓점 θ 없이) | 호의 끝은 짝 곡선의 격자 열에 서지 않는다 — 반지름 50 원판의 림을 상자 벽에서 깎은 호는 그 끝과 아래 원의 열 사이 띠가 여전히 기울어 `OverBudget` 이었다(타원과 그 원의 봉만 고쳐졌다). 집합은 격자와 그 면 모든 꼭짓점의 θ 다(`FaceThetas`) | 스파이크 |
| 곡선만 면의 θ 집합에, 격자 줄은 격자 열에만 | 줄과 림 사이 띠에 꼭짓점 열이 끼면 줄에는 그 생성선이 없다 — 윗림 0.02 아래 줄과 밑 노치(높이 5~100)의 원판에서 줄과 림 사이 면이 1.16° 까지 기울었고(예산 2°), 줄도 같은 집합을 읽자 0.333°(칸의 1/6) 다. 한 집합을 두 독자가 읽는다 | 스파이크 |
| 옆면의 θ 집합을 원통 원의 고른 격자 `τk/N` 에 두고 예산을 넘는 칸만 국소로 나누기(이분, 또는 같은 조각) | 원에 가까운 타원에서 칸마다 «1 → 2» 가 강제되어 장축 끝의 칸으로 고르게 나누기보다 나쁘다 — 회전율 1.03 의 반지름 2 림 186 → 272, 반지름 20,000·1.03 은 3,189 → 6,094 로 `FINEST_CIRCLE` 을 넘어 새로 거절된다(모형으로 잼). 바탕 격자의 칸이 국소 요구와 정렬되지 않기 때문이다. 원통의 `ref_dir` 에서 그 자리가 허락하는 칸으로 걷는 걸음(`walk`)은 같은 모형에서 183·3,166 | 모형 |
| Sutherland–Hodgman 으로 비볼록 링 클립 | 교차를 링 순서로 짝지어 떨어진 조각을 클립선을 따라 이어 버린다. 선 순서로 짝짓는다 | `d3e9984` |
| 반사를 켤레 회전으로 나르기(`M∘R = (M R M⁻¹)∘M`) | 정확했지만 한 질문에 두 메커니즘이 생기고 사슬이 방금 한 반사를 스스로 말하지 못한다 — 미러로 도달한 벽과 이동으로 도달한 같은 벽이 1 ULP 갈라져 한 부품이 두 몸통이 됐다. 반사는 사슬의 노드다 | — |
| 빈 불리언 결과를 오류로 | 빈 결과는 실패가 아니다 — `Ok(vec![])` | — |
| 이름 없는 단면 클래스를 `crossings` 에 넣기 | 그 집합의 명제는 «유리수 clearance 로 증명된 반지름 안쪽»이고 독자들이 클래스의 유리수 계수를 곧바로 읽는다 — 버리는 스파이크에서 이 길은 30/30 을 빌드했지만, 캡 클래스가 W 의 선분을 필렛 원과 다시 묻고(200회) `chord_stations` 가 W 를 기록된 짝처럼 다룬다. 게이트의 말 없는 갈래 셋째(`Silence::Sectioned`)로 두면 «증명» 독자들이 그대로 W 의 선분을 건너뛴다 | 스파이크 v4 |
| 세계 이름 없는 클래스를 원통에 대해 무조건 «말 없음»으로 통과시키기 | 떨어짐을 증명하지 않으면 틀린 답이 조용히 나간다 — 돌린 상자가 핀을 감싸면 Common 이 50(참 `50 + 3π`), 구멍을 가로지르면 메셔가 받지 않는 몸, 판 위로 뜬 상자의 Fuse 는 광선 자리의 `debug_assert` 패닉(임시 계기). 게이트가 옆면 상자로 비킴이나 떨어짐을 증명한 클래스만 지나고, 못 가르는 독자는 게이트의 이름으로 거절한다 | 스파이크, 되돌림(`cd849c28` 본문) |
| 말 없는 클래스의 «쪽»을 `orient3d_filter(t₀, t₁, t₂, v)` 의 부호 그대로 기록하기 | 판정(`Judge::orient3d`)의 쪽과 부호가 반대다(`v` 가 행렬의 끝에서 처음으로 가는 홀순열) — 관통 꼭짓점과 세 평면 꼭짓점이 섞인 링(필렛 하나인 판의 윗면)이 뒤집힌 쪽을 없는 교차로 읽어 결과가 여러 덩어리로 갈라지고 `StraightAngle`; 링이 관통 꼭짓점뿐인 판(필렛 넷)은 한결같은 뒤집힘이라 못 본다. 쪽은 `WorkingPlane::misses` 한 곳에서 `side_of` 관례로 만든다 | 스파이크, 되돌림(`77a8eb81` 본문) |
| 공유 룰링 정거장을 건너 옆면 셀을 잇기 | 두 벽이 모두 원통으로 들어가는 쐐기의 `A − B` 가 한 몸통으로 나왔다 — 옆면이 선을 통과하는 한 면이고 두 벽 면만 그 선을 간선으로 나눠, 스스로 닿는 솔리드가 `validate` 깨끗·부피 항등식 성립인 채로 나갔다. 옆면의 재료 쪽에서 두 벽 사이의 방이 결과에 없으면 옆면을 그 선에서 끊는다(`slit_at`) | — |
| 옆면의 제 간선 위 정거장을 차트에서만 짝의 마디로 자르기(배열의 룰링 조각은 그대로) | 반원기둥과 접선·할선 벽이 그 간선에서 함께 끝나는 삼각형 96 가운데 60 이 서고 36(Fuse 두 순서·`A − B`)이 `EdgeCurveUnderived`: 옆면은 그 선을 상대 꼭짓점에서 끊어 쓰는데 평면 면의 고리는 한 간선으로 써서 간선 짝이 안 맞는다. 마디는 배열(면 고리)과 차트가 한 벌로 가진다 — 룰링 조각을 배열에서 자른다 | — |
| 룰링 위의 칼날 모서리를 골격에 남겨 접기로 룰링에 넘기기 | 같은 삼각형이 서지만 무소식 거르기의 예외와 접기의 무소식 버림, 두 예외로 돌아간다; 예외를 모든 클래스로 넓히면 세계 틀 24 가 `StraightAngle` — 룰링을 그리지 않는 벽 클래스에 남은 칼날이 면 고리의 안테나가 되고 고리의 극점이 그 끝이다. 규칙은 «선은 그 위의 꼭짓점에서 잘린다»라 칼날은 버리고 끝만 버리기 전에 자른다 | — |
| 정거장에서 선 위에 결과 면이 하나라도 있으면 옆면을 끊기(재료 아닌 쪽의 면까지) | 스스로 닿는 몸(보어를 건너 판까지 드는 쐐기, 원통을 감아 관통하는 C자)은 거절되지만, 선으로 닿는 두 몸(선 접촉·접선 옆면의 합집합)의 옆면에 그 선이 간선으로 생겨 `EdgeCurveUnderived` 로 깨졌다(점수판 두 테스트). 한 몸인지는 껍질의 연결이 정한다 — 결과 계약의 옆면 팔 | — |
| 재료가 원통 바깥쪽인 옆면에서 정거장을 끊지 않기(안쪽 방만 읽고, 바깥쪽이면 «이어짐») | 원통을 감싼 상자에 꼭짓점이 옆면 위인 V 홈을 낸 몸에서 원통을 빼면 홈 꼭짓점 선에서 결과가 스스로 닿는데, 옆면이 그 선을 통과하는 한 면으로 나가 `validate` 깨끗·부피 일치인 채로 지어졌다. 바깥쪽이면 두 벽의 바깥 반직선 사이 방을 읽는다 — 이웃한 두 벽의 바깥 칸 라벨(`RulingExtent::outer`, 정거장에서는 `ChordSays::Piece` 의 `outer`) | — |
| `shared_rulings` 의 짝을 두 클래스의 면이 그 선에 닿는 것으로 좁히기 | 짝은 정거장을 묶는 독자(씨앗·접기·룰링 분할·접선 판정 — 「인구 게이트」의 짝 문단)가 면과 무관하게 읽는다. 면이 옆면 구간에서 선에 닿지 않는 짝(스위트 할선 둘 28 · 접선 39, census 접선 15)을 빼자 스위트 다섯이 거절됐다 — 할선 둘 짝만 빠진 둘(원통 림 위 아치의 Fuse `CoincidentNodes`, 집힌 캡 피연산자 `RulingBoundNotYet`), 접선 짝만 빠진 둘(선분 끝에 놓인 림 위 교차 `CoincidentNodes`, 면 밖 셋째 평면 곁 접선 `RulingBoundNotYet`), 둘 다 빠진 하나(먼 계단 곁 구멍을 지나는 접선 `CoincidentNodes`). census `c` 행은 같았다 | 스파이크 |
| 이름 있는 클래스의 `ClassMeet::Clear` 를 하류에 «말 없음»으로 넘기기 | 여는 것이 없다 — 잰 독자에서 정확 도로가 그 쌍에 늘 같은 답을 냈다: 분할의 교점 후보 스위트 768·census 468 전부 조각 밖, 원판의 이름 꼭짓점 스위트 6,870·census 3,108·무시 스윕 3,240 전부 밖, 좌표 증인의 상자 대체 0. 바뀌기만 한다 — 광선이 기권한 이름 꼭짓점에서 다음 증인 대신 상자 탐침이 답하는 불리언이 스위트 1,385·census 81·무시 스윕 1,272(상한), 그리고 재지 않은 독자의 갈래: 말 없는 클래스의 면 판정(`point_in_faces_rat` — 못 보이면 `CylinderGateUndecided`)·광선에 건네는 원통 표(`silent_table`)·거절 이름(`undecided_reason`) | — |
| 접선 공유의 이름 묶음을 옆면 면 밖의 점(무한 원통 위)까지 넓히기 | 반원기둥 위에 돌린 상자의 모서리(`z = 20`, 옆면은 `z ≤ 10`)가 관통 이름을 대표로 얻어, 이름 없는 평면으로 그 결과 꼭짓점을 짓지 못하고 `CylinderGateUndecided` 로 거절했다(`draft.rs`, 상자 − 반원기둥). 묶음은 정거장이 옆면 면의 축 범위 안인 클래스에만 | `8e983cf0` |
| 이름 없는 클래스를 `Silence::Apart` 로 통과시키되 그 원통에 다른 기록(절단·접선)이 있으면 불리언 전체를 거절하기 | 근거(«메운 단면이 관통점을 없는 것으로 읽는데 기록은 다른 클래스의 선 위에 그 점을 둔다»)는 그 클래스의 면이 없는 자리 — 옆면 상자 안 — 의 일이었다. 막대 ∪ 37° 상자 가족 여덟(구멍을 자르거나 접함, 핀을 자름, 벽 평면이 막대의 룰링을 지남, 소유 뒤집기, 상자 면이 절단 평면을 상자 밖에서 가로지름 — 막대가 구멍에서 멀어도, 비스듬한 벽) × 두 순서 × 세 연산 48 이 거절됐고, 거절을 걷자 47 이 유도 부피(1e-9)·validate 0 으로 지어지고 하나는 막대 단독과 같은 `SelfTouchingResult` 였다. 꺼도 스위트·census 가 초록이라 이 거절을 지키는 시험은 없었다. 잠금은 일곱 가족 × 두 순서 | — |
| 참 읽기(`Silence::Sectioned`)를 온전한 림이 있는 원통에도 열되 메운 읽기(`Silence::Apart`)를 먼저 묻기 | 온전한 림만 남긴 핀은 메우고 그 핀을 꿰는 구멍(밑면 턱이 아래 림을 끊는다)은 참으로 읽어, 축을 지나는 한 벽 클래스가 두 읽기를 안팎에 겹친다 — 두 순서 × 세 연산 여섯 칸이 모두 `LabelConflict` 다. 단면을 먼저 물으면 여섯 모두 유도 부피로 선다 | `1b2154be` |
| 동일평면 면을 이어붙여서 정리 | 비스포크였다. 경계 대수로 일반화한다 — 이어붙이지 말고 내부를 지운다 | — |
| 면 위 직접 구성(`ImprintSketch`·`raise_region`) | 소비자 0 이었고 커널 유일의 동일평면 인접면 생산자였다. pad = extrude + Fuse, pocket = extrude + Cut. Split Face 의 소비자(구역별 재질·FEA 경계조건·파팅라인)가 생기면 두 제거 커밋의 revert 가 출발점이다 | — |
| pad/pocket 을 커널 연산(`PadOnFace`·`PocketOnFace`)으로 | 새 정확 술어 없는 «면 프레임 + 돌출 + 불리언» 래퍼였다(설탕 판별 기준). 커널에 남은 이유는 «면에서 안쪽으로»를 기본 어휘가 말하지 못해서였다 — 돌출이 부호 있는 거리를 받자 사라졌다. 바닥 면을 돌려주려는 캡 회수가 **모든** 불리언에 입력 평면 면마다 클래스·향 표(`ClassOf`, 면마다 `face_facing`)를 계산하게 했고, 깊은 pocket 을 `PocketNotBlind` 로 거절하게 했다. kit 은 캡을 읽지 않았고, 가른 몸의 나머지 조각을 값에서 잃었다 | — |
| kit 의 중심 피벗(`rotate…(…, { pivot: "center" })` — 몸의 `props::centroid` 를 피벗으로) | 중심은 f64 측정이고 피벗은 구성의 입력이다 — kit 이 측정값을 짧은 십진수(최대 17자리)로 올려 정확 구성에 되먹이는 유일한 문이었다. 곡면 몸의 중심은 대개 π 가 든 초월수이고 평면 몸도 `5/6` 이 `0.8333333333333334` 가 되므로, 대칭 부품을 중심 둘레로 180° 돌려 원래 몸과 겹치기를 기대하면 1e-15 어긋난다 — 그 한 옵션만 «정확하지 않다»고 기억해야 하는 예외였다. 피벗은 원점이거나 스크립트가 적은 점이다 | kit `02f0789` |
| 정확한 연산(정확한 0 과의 합, 0·2의 거듭제곱과의 곱)의 규칙을 `Bounded` 의 `add`·`sub`·`mul` 안으로 | 판정 필터도 그 경우를 자주 만나지만(덧셈의 42–51%, 곱의 31–40%) 판정이 하나도 바뀌지 않는다 — 두 모양 모두 상승 수·census `c` 행 동일. 검사는 필터 연산마다 붙어 perf 명령어가 늘었다(바로 반환 +8.2%, 청구만 고르기 +16.9%; census +0.2~0.5%). 규칙은 tol 0 이 계약인 증인점 단계(`witness.rs` 의 `sum`·`diff`·`prod`)에 산다 | — |
| 직선각 정리를 평면 면에 한 번, 옆면을 붙인 뒤 결과 전체에 한 번 더 | 결과는 한 번 정리와 같았다(스위트 6,066·census 390 실행에서 면·림 동일). 그러나 같은 전체 물음을 두 번 묻고, 첫 벌이 옆면 없이 내린 결정을 둘째가 덮으며, 옆면은 첫 벌 뒤의 림 표를 따로 읽어야 한다. 정리는 옆면이 붙은 뒤 한 번(`whole_result_faces`) — 덩어리별 정리(`naming`)는 몸 하나씩 묻는 다른 물음이다 | 스파이크 |
| 옆면 룰링을 정리된 평면 면이 쥔 점에서만 자르기(`HeldRims` 의 룰링판) | «결과의 이 선 위 꼭짓점»이라는 한 명제의 곡선 종류별 둘째 철자다 — 비스듬한 단면의 호가 셋째를 부른다. 정리가 옆면까지 보면 림과 룰링이 한 규칙으로 선다 | — |
| 비스듬한 클래스를 그 면이 모두 옆면 상자에서 떨어질 때만 들이기(`faces_apart`) | 접촉만인 세 가족(구멍 위에 선 거싯, 접선 모서리 위 쐐기, 홈 입구 쐐기)의 18칸을 거절했고, 그 요구가 다른 일 둘을 대신 지고 있었다 — 자기접촉이 옆면 한 점 접촉을 묻지 않는 공백(요구를 걷자 꼭짓점이 보스 옆면에 닿는 몸과 간선이 옆면에 접하는 몸이 `Ok n=1`·validate 0)과, span 끝 림만 읽는 림 접촉 검사(창이 난 옆면의 호에 비스듬한 면이 한 점으로 닿는 몸도 `Ok n=1`). 게이트는 단면을 그릴 수 있는지(어느 원에도 한 점으로 접하지 않는가)만 묻고, 접촉은 결과 계약의 옆면 팔이 맡고, 타원 간선은 결과가 든다 | `8ba299fb` |
| 옆면 팔이 닿는 점을 모두 점 접촉으로 묻기 | 평면 면이 옆면에 닿는 배치(벽·판 윗면에 접하는 스터드·보스·핀, 판 옆면에 앉은 보스의 밑 원이 윗모서리에 접함)에서 간선 근이 «닿음»이 되어, 빌드되던 스위트 12 결과와 census 2 칸을 `SelfTouchingResult` 로 거절했다. 그 근은 두 부류다 — 접선 판정의 행이 말하는 선 접촉의 점(스위트 3 근), 옆면과 이웃 면이 함께 쓰는 림 위에 그 면의 간선이 닿는 그 면의 집힘(스위트 51·census 7·kit 4 근) | `b67a2b8a` 본문 |
| 접평면의 결과 면이 덮는 점을 모두 접선 판정에 넘기기 | 접선 판정은 행 — 벽 면이 옆면 위에서 접선과 길이 있게 만나는 구간 — 만 판단해, 접평면 위 삼각형의 꼭짓점 하나가 보스 옆면에 닿고 기둥으로 판과 이은 몸이 `Ok(1)`·validate 깨끗으로 빌드됐다(그 넘김 앞의 게이트는 `ObliqueCylinderCut`). 넘김의 영역은 행의 구간이다 | `b67a2b8a` 본문 |
| 넘김을 판정의 행 구간만으로 좁히기 | 옆면과 이웃 면이 함께 쓰는 림 위에 그 면의 간선이 닿는 근(판 옆면에 앉은 보스의 밑 원이 판 윗모서리에 접함 — 스위트 16 근, 시험 넷)이 «닿음»이 되어 빌드되던 결과를 거절했다. 그 점은 그 면의 두 경계가 만나는 집힘이다 | `b67a2b8a` 본문 |

### 검증·오라클

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
| OCCT 의 STEP 읽기 수리를 통째로 끄고(`FromSTEP.exec.op` 비우기) 방향을 검사하기 | 곡면 위 면이 늘 무효로 읽힌다 — 파일에 PCURVE 가 없다(원통 r=2.5·h=7: 부피 137.445 → 68.72, valid 0; seam 을 넣은 원통도 같다). seam 없는 옆면은 그 위에 수리의 `FixMissingSeam` 에도 기댄다. 그래서 평면 형상만 검사할 수 있었다. 방향 수리 셋만 끄면 같은 카나리(면·모든 면·간선·경계·void)를 평면과 곡면 모두에서 잡고, 정상 형상은 곡면까지 유효하다 | — |
