# nacre — 하이브리드 CAD 커널 설계

정확한 기하를 진실로 보관하는 고전 b-rep의 골격 위에, Fornjot에서 검증된 위생 규율(append-only 단일 참조, 위상과 근사의 동시 구축)을 얹은 설계다. 목표는 정밀 기계 CAD(STEP 입출력 포함)이며, 개인 + AI 협업 개발을 전제로 실수 여지를 줄이는 인프라를 처음부터 포함한다.

## 설계 원칙

이 문서의 모든 결정은 다섯 가지 절대 원칙에서 나온다. **원칙의 유일한 원본은 `docs/overview.md` 의 「절대 원칙」이다** — 여기에 다시 적지 않는다(두 곳에 적으면 한쪽이 상한다). 범위(비목표)의 원본도 같은 문서의 「범위」다. 이 절은 그 원칙들이 이 설계에서 갖는 정확한 경계만 덧붙인다.

**범위가 커널에 요구하는 것은 하나뿐이다.** GD&T/PMI·스타일·색상·레이어·제작 메타데이터·제품 구조는 커널을 쓰는 응용의 영역이고, 커널은 그들을 위해 단 하나만 제공한다 — 영구 유효한 Handle(append-only의 부산물). 응용은 `HashMap<Handle<Face>, 응용데이터>` 사이드카로 무엇이든 매달 수 있다. nacre-step은 STEP의 형상 서브셋만 커널로 번역하고, 비형상 엔티티는 해석 없이 무손실 패스스루로 보존해 라운드트립을 지킨다.

**「연산 이력을 보존한다」의 보장 범위.** replay는 동일 로그·동일 파라미터에서 동일 모델을 보장하고(undo/redo의 기반), tolerance 변경 재계산은 위상·기하를 불변으로 둔 채 tessellation만 재생성한다(「Tessellation 층」 절). 로그 중간의 파라미터를 수정하는 파라메트릭 편집은 v1 비목표다 — 연산이 원시 Handle을 참조하는 한 상류 수정이 하류 Handle 번호를 밀어내기 때문이다(topological naming 문제). 진화 경로는 「연산 층」 절에 있다.

**「tolerance는 발견된 교차에만 존재한다」를 드는 자리.** 구성 시점에 동일성을 아는 요소에는 tolerance 개념 자체가 없다. 그 구분은 **정점의 실현 캐시**(`PointCache`)가 든다. 캐시에 저장된 tolerance 는 **없다**. 변종 셋은 «캐시가 무엇을 아는가»만 말한다 — 증명된 경계가 있거나(`Bounded`), 싼 도로가 멈췄거나(`Ceiling`: 비트 또는 비용, 둘 다 «더 물으면 답이 있다»), 실현이 뒤에 없거나(`Unrealized`). 발견된 점도 증명된 경계만 든다 — seam 표가 정점이 받을 캐시(`SeamVertex.cache`)를 들고 자기접촉 체가 그 경계를 읽는다.

## 크레이트 구조

의존 방향은 아래에서 위로만 흐른다. 순환 의존 금지.

```
nacre/                    # 워크스페이스(crates/ 아래). 최상위 `nacre` 크레이트는 파사드(재수출 전용)
├── nacre-store      # typed-index 저장소: Store<T>/Handle<T> (기하·위상 무지, nacre-* 의존 없음)
├── nacre-math       # 자체 선형대수: Point<D>·Vector<D>·변환 (nacre-* 의존 없음)
├── nacre-exact      # 유리수·정수로 정확히 답할 수 있는 모든 것: Rat·Angle·Isometry·PlaneName·QuadVal·정확 술어 (nacre-* 의존 없음)
├── nacre-predicates # exact f64 부호 술어(indirect predicates); geometry-predicates 위·standalone (nacre-* 의존 없음)
├── nacre-judge      # 부호 판정: 정확히 답할 수 없는 회전 좌표를 f64 필터→상승→«증명 또는 미결»로 정한다
│                    #   kernel(WitnessPoint) + predicate(질의를 정확 경로나 kernel 로 보내는 라우터)
│                    #   ← exact · predicates · math
├── nacre-geom       # f64 기하 «실현»(Surface·Curve)·교차(intersect 격리)
│                    #   ← math · predicates · exact(Rat 링 술어 쌍둥이)
│                    #   곡면의 진실은 topo 의 Surface — 맨이름은 진실, 실현은 nacre_geom::Surface
├── nacre-topo       # b-rep 위상: Vertex/Edge/Face/Shell/Solid·half-edge·Model   ← store · geom · math · exact
├── nacre-tess       # tessellation: 출처 태그   ← store · topo · geom · math · predicates
│   └── polygon      # 평면 다각형 삼각분할: y-단조 분해 + 단조 삼각분할 + Delaunay 플립
├── nacre-validate   # 불변식 검사: 오일러-푸앵카레·watertight·방향성·참조 무결성   ← topo · store · math · geom
├── nacre-props      # 질량 특성: 부피·면적(해석적, tess 무관)   ← topo · geom · math · store
├── nacre-ops        # 연산: datum 평면/sketch/extrude/boolean/transform/mirror/copy
│                    #   불리언은 한 방향이다: boolean(앞문) → arrangement(엔진) → assembly(조립),
│                    #   셋 다 draft(지어지기 전 면의 어휘)를 읽는다. 관문이 그 방향을 단언한다
│                    #   ← topo · geom · math · store · exact · judge (+ 선택적 rayon)
│                    #   부울 = 면당 평면 arrangement 엔진
│                    #   reuse: 상대가 닿을 수 없다고 증명된 평면 클래스는 배열하지 않는다 (「연산 층」 절)
├── nacre-step       # Model→STEP(AP242) 내보내기 어댑터   ← topo · geom · math · store (+ step-io)
├── nacre-oracle     # [dev] OCCT 비교 하네스   ← topo · step · math · store
├── nacre            # 파사드   ← exact · geom · math · ops · props · step · store · tess · topo · validate
└── tools/occt-helper/  # 워크스페이스 밖 헬퍼: brew OCCT(1순위) 또는 uv+OCP(폴백) — 「검증·오라클 인프라」 절
                        #   OCCT는 오라클 전용 — 제품 경로에 위임 없음
```

제품 경로에서 `nacre-ops` 는 `nacre-tess`·`nacre-validate`·`nacre-props`·`nacre-step` 에 의존하지 않는다(테스트 전용 dev-dependency 로만 쓴다 — 발행되는 라이브러리의 의존 표면이 깨끗하게 남는다). `nacre-exact` 도 `nacre-predicates` 를 dev-dependency 로만 쓴다(정수 부호 술어를 f64 쌍둥이와 차등 시험) — 제품 경로의 방향은 predicates 무의존이다.

`nacre-geom`과 `nacre-topo`의 관계는 한 방향이다. 기하는 위상을 모르고(순수 수학), 위상은 기하를 **캐시 값으로** 든다(`SurfaceCache`·`EdgeCache` 가 geom 의 `Surface`·`Curve` 를 감싼다 — 진실은 topo 자신의 타입이고 geom 은 그 실현이다).

**어느 크레이트에 두나 — 축은 «이름»이 아니라 «정확성»이다.** 주제로 그린 지도(「숫자·기하·위상」)는 경계를 설명하지 못한다 — `nacre-exact` 공개 이름 84개 중 **61개가 기하**(`cylinders_nested`·`SeamOrder`·`PlaneName`…)이고, 그 크레이트가 격리하는 것은 스칼라가 아니라 **«유리수·정수로 정확히 답할 수 있는 모든 것»**(값·술어·이름)이다. `nacre-geom` 은 같은 기하의 **f64 실현**을 맡는다 — 진실/캐시와 같은 축이다. 규칙: 「유도된 값 + 그 산술」→ `nacre-exact`, 「아레나 항목의 진실」→ `nacre-topo`. 이것을 강제하는 것은 취향이 아니라 의존 그래프다 — **`nacre-judge` 는 `nacre-topo` 에 의존하지 않으므로** 판정이 쓰는 타입은 topo 아래에 있어야 한다(`PlaneName` 이 topo 가 아니라 scalar 에 사는 이유). 분리 여부는 「한쪽만 쓰는 소비자가 있는가」로 재고, 그런 소비자는 **0** 이라 분리하지 않는다. robustness가 첨예한 코드(교차·분류)는 전부 `nacre-geom::intersect` 한 모듈에 격리한다. 사용자는 파사드 크레이트 `nacre` 하나만 의존하며, 인터랙티브 스크립트 앱 등은 이 워크스페이스 밖의 별도 프로젝트로 둔다.

**파사드 `nacre`.** 열 개 층을 **모듈로** 재수출하고(`nacre::exact`·`geom`·`math`·`ops`·`props`·`step`·`store`·`tess`·`topo`·`validate`), 자주 쓰는 것은 `nacre::prelude`에 담는다.
- **평면(flat) 재수출은 하지 않는다.** **층 분리가 이 설계의 뼈대**여서 이름공간에 남긴다. 같은 맨이름이 층마다 다른 것을 뜻하기도 한다 — `topo::Surface` 는 진실, `geom::Surface` 는 그 f64 실현이다.
- **`nacre-exact` 재수출은 선택이 아니다.** `Operation::Transform { isometry: Isometry }`·`Mirror { axis: Axis, offset: Rat }`가 scalar 타입을 ops의 공개 API로 새어 보내므로, 없으면 소비자가 그 op을 만들 수조차 없다. 반대로 `nacre-judge`·`nacre-predicates`는 공개 API에 새지 않아 재수출하지 않는다(퍼블리시는 필요 — ops·geom 의 하드 의존).
- **prelude의 기준은 확인 가능한 성질이다**: *"`Operation`의 모든 변이가 prelude 이름만으로 만들어진다."* 내용은 실제 소비자의 import 와, 그 소비자보다 나중에 생긴 표면(스케치 앞문·파생 조회·거절 사유)의 합집합이다 — 낡은 소비자만 보면 최신 API가 빠진다.
- **기능**: `parallel`(기본 on)이 `nacre-ops/parallel`로 전달된다. `nacre-ops`를 **`default-features = false`로** 매달아야 소비자가 끌 수 있고, 그러지 않으면 `--no-default-features`에도 rayon이 들어온다. 이 속성은 매니페스트에 살아 훅이 검사하지 않으므로 **테스트가 `Cargo.toml`을 직접 확인**한다. `parallel`은 **`Sync` 스위치이기도 하다**(cip의 hp 캐시가 `Arc<OnceLock>`↔`Rc<OnceCell>` — 워커들이 평면표를 *공유*해야 하므로 하중은 공유 참조 쪽이다). `test-util`은 topo의 테스트 전용 문(심기용 push·`reversed_shell`·`plane_name_through`)과 push 계측(`WIDE_PLANES`·`SEEDED_HITS`·`surface_derive_counts`)을 전달한다 — 제품 빌드는 세지 않는다. **순차 조합은 반드시 `-p nacre-ops --no-default-features`로 확인한다** — 워크스페이스를 통째로 지으면 형제 크레이트들(`nacre-oracle`·`-props`·`-step`·`-tess`·`-validate`)의 dev-dependency와 파사드의 기본이 기본 피처를 도로 켜서 아무것도 안 재게 된다. 그 명령이 순차인 것은 ops 의 자기 dev-dependency가 기본 피처를 고르지 않기 때문이다(`tests/instruments/manifest.rs` 가 모든 자기 dev-dependency에 건다).
- **검증은 "`nacre::` 경로만으로 끝까지 가기"다** — 크레이트 문서의 예제(doctest)와 `tests/facade.rs`. 재수출이 빠지면 컴파일이 깨진다. 문서 예제는 **프로덕션 API만** 쓴다(`test-util`이 필요한 예제는 독자가 실행할 수 없다).
- **`Document`(op 로그 + Model + tess 캐시) 번들은 두지 않는다** — 「저장소와 동일성」 절이 이 층을 지목하지만 의미론을 정하는 새 타입이므로, kit이 무엇을 원하는지 보이기 전에 정하지 않는다(계획).

**편의 레이어 `nacre-kit` (워크스페이스 밖, 별도 리포).** 코드-CAD 스크립트와 커널 사이의 층: 다중 솔리드 값(compound), 값 의미론(재사용 시 `Copy` 자동 삽입), 다인수 fuse/cut/common(fold), 프로파일 헬퍼와 섬-분해 호출, 패턴·미러, 에러의 사람용 매핑, 표시 메타데이터(색·투명도 — 커널 비목표라 여기가 제자리). **Rust로 두는 이유:** 헤드리스 `cargo test`가 되고, 술어 인접 로직이 exactness 도구가 있는 쪽에 남고, 프론트엔드를 교체해도 살아남고, wasm 경계가 함수 하나로 유지된다. 경계 규칙은 overview.md의 "설탕 vs 커널 판별 기준"이다. 스크립트 문법은 `nacre-playground` 가 구현하고, 의미론과 그 이유는 `nacre-kit` 의 코드와 주석에 있다(여기에 복사하지 않는다 — 두 곳에 같은 내용이 있으면 어긋난다).

**공개 표면 — 밖에서 무엇이 되는가.**
- **`Model` 의 저장소는 비공개이고 문(accessor)으로 읽는다.** 곡면·모션·위상 아레나·캐시·interning 표가 전부 비공개이며, 공개 필드는 `surface_name` 하나다(「문의 이름」 절이 캐시 조각으로 접을 몫). 전량 순회 문은 의도적으로 없다 — 아레나는 supersede 된 항목도 들고 있으므로 소비자는 live 면을 걷는다.
- **위상 순회는 밖에서 된다.** 위상 셀(`Vertex/Edge/Face/Shell/Solid`)의 필드는 `pub` 이고 `Model::reachable()`→`Reachable{vertices,edges,faces,shells}`도 공개다(`shell.faces → face.outer.half_edges → edge.vertices` → 좌표는 `Model::vertex_point(vh)`).
- **내부 정보도 읽힌다.** 피킹용 `Tessellation{by_face,by_edge,…}`·`TessTriangle.face`·`TessOrigin`, 정점의 정의 `Vertex{ThreePlane|OnSeam|Pierce}`·실현 캐시 `Model::vertex_cache`·`Model::motion(h)`(필드가 아니라 **좁은 문**). 디버그 뷰어가 읽어야 할 것은 다 읽힌다. `tessellate`·`Tessellation::to_obj`·`to_step`·`to_step_solid`·`validate`·`mass_props`도 공개.
- **파생 값과 에러 표면.** `nacre-props`에 `bounds`·`centroid`·`face_props`(넓이·중심·법선), `nacre-ops`에 `face_plane`, `nacre-topo`에 `Model::he_start`. **원칙: 값을 돌려주는 읽기 전용 질의**(위상 순수성 유지). `TessConfig`는 `tol` 하나뿐 — 면별 override는 계획.
  - **`bounds`는 곡선을 인지한다.** 원통 옆면은 솔기 정점보다 바깥으로 볼록하므로 꼭짓점 min/max는 **조용히 작은 상자**를 준다. 반지름 `r`·법선 `n̂`인 원은 축 `e` 방향으로 `±r·√(1−(n̂·e)²)`만큼 뻗는다(정확). OCCT `bounding`이 심판하되 **등호로 비교하지 않는다** — DRAWEXE는 상자를 보수적으로 부풀린다(~1e-7).
  - **`centroid`는 새 적분이 아니다.** 솔리드는 기준점에서 각 평면 면으로 뻗은 **원뿔들의 부호합**이고, 원뿔의 중심은 밑면 모양과 무관하게 꼭짓점→밑면중심의 **3/4** 지점이다. 즉 `mass_props`가 이미 계산하는 `(Aᵢ, cᵢ, n̂ᵢ)`만으로 `C = R + Σ Vᵢ·¾(cᵢ−R)/ΣVᵢ`가 나온다. 곡면은 그 논증이 깨지므로 **이름 달고 거절**하고, 그래서 `MassProps`의 필드가 아니라 별도 함수다(곡면 솔리드의 부피·넓이는 계속 살아 있어야 한다).
  - **`face_props.normal`은 `Option`이다.** 원통 면에는 하나의 법선이 없는데, **면 고르기는 모든 면을 훑는 일**이라 실패시키면 필터가 통째로 망가진다. `None`이면 자연스럽게 건너뛴다.
  - **면 고르기가 위상 명명 문제를 우회한다.** 코드-CAD는 면을 번호가 아니라 `filter(법선≈+Z).max_by(중심.z)`처럼 **생김새로** 고르고 매 실행 다시 고른다 — 저장된 참조가 없으니 상류가 바뀌어도 썩지 않는다.

**면의 스케치 좌표계 — `face_plane`.** 면의 `SketchFrame`(`face_sketch_frame`)을 실현한 것이고, **그 프레임 위의 `Extrude`(pad 의 보스, pocket 의 공구)가 실제로 프로파일을 놓는 바로 그 프레임**이다(두 번째 유도가 아니라 한 값의 실현 — 갈라지면 앱이 계산한 위치와 보스가 어긋난다. 테스트가 비대칭 프로파일로 고정한다). 실현은 도로마다 그 도로의 것이다(`frame_basis`): 세계 도로는 유리수 세계 기저를 한 번 반올림하고(그 도로의 꼭짓점이 그렇게 실현된다), 프레임 노드 도로는 사슬을 재생한다.

**공개 스케치 어휘 — `SketchFrame` + `face_sketch_frame`.** `SketchFrame{plane: Handle<Surface>, placement, flip}`은 공개 타입이되 필드는 비공개이고, 생성자가 검증한다: `canonical(plane)`은 유도라 검사 없음, `named(model, plane, origin, ref_dir)`는 구성 시점에 정확 검사해 이름 붙은 거절을 낸다(`FrameOutsideDecimalWindow`·`OriginNotOnPlane` — scalar의 `plane_residual_sign`, Wide 이름은 BigInt 팔 —·`RefDirParallelToNormal`). `face_sketch_frame`은 면의 프레임을 진실에서 고르는 한 문이다(`face_plane`은 그 실현). **면이 어느 프레임을 갖는가는 면이 어디 있는가로 정하지, 평면이 어떻게 저장됐는가로 정하지 않는다** — 이동을 점에 옮겨 적었는지 노드로 기록했는지는 변환이 정확히 진술할 수 있는 것으로 고르므로, 그것을 읽는 규약은 같은 두 면을 두 프레임에 스케치한다. 세계 방정식이 진술되는 평면(`world_plane_name` — 모션 없음, 또는 접히는 사슬)은 **바깥 법선의 세계 arbitrary-axis 프레임**(`plane_frame_default` — 세계 원점의 투영, `+u = ẑ × n`, 수직이면 `ŷ × n`)을 갖고, 그것을 사슬의 역(`chain_point_rat_inverse`/`chain_dir_rat_inverse` — 접히므로 정확)으로 평면의 진술 좌표계에 옮긴 `Named` 로 적는다 — 같은 `flip` 으로 평면 이름이 유도하는 프레임과 같으면 `Canonical`(정규형). 거울이 든 사슬로 옮긴 프레임은 왼손이다(`û` 는 세계의 것, `v̂` 는 반대, `ŵ` 는 바깥). 그 밖의 평면(프레임 노드·사분각 밖 회전·Wide·이름 없음)은 자기 `Canonical` 프레임을 사슬로 운반해 갖는다. 그래서 모든 평면 면이 철자를 가지며, 계약(실현이 `face_plane` 과 비트 동일)은 `tests/invariants/sketch_frame_contract.rs` 가 여덟 배치 × 6면으로 잠근다. flip 은 `frame_toward`, 노드 push는 `push_frame_node` 한 곳으로 통일돼 extrude·face 두 도로가 한 모양이다. `Operation`은 평면을 핸들로 싣는다(`Extrude { frame: SketchFrame, .. }`, `DatumPlane { def: DatumDef }`) — replay 자기완결성: 로그 속 평면 핸들의 합법 표적은 씨앗·기존 면·datum뿐이다. 그리고 `Model::new()`가 세계 축 평면 셋을 심는다(핸들 0·1·2, 캐시 방향 −축, `world_plane(Axis)` 접근자, `Default`는 `new()` 위임) — 세계 평면 위 스케치와 원점 상자의 축 면이 같은 surface 핸들을 공유한다.

**면 프레임의 원점은 세계 원점의 정사영이다 — 면을 전혀 읽지 않는다.** 꼭짓점 평균은 **직선 도중에 꼭짓점이 하나 늘면 움직이고**, 면적중심은 이산화에는 무관하지만 면의 `f64` 꼭짓점을 읽어 반올림된 캐시가 진실이 된다(같은 footprint 를 두 번 pad 하면 넓이 `2.2e-16` 의 면이 남았다) — 그리고 둘 다 같은 평면 위 다른 외곽선이면 원점이 달라진다. 정사영은 평면만의 성질이다(잠금: `the_sketch_origin_does_not_depend_on_the_outline_at_all`). 대가: 원점에서 먼 면에서는 스케치 좌표 `(0, 0)` 이 면 밖이므로, 호출자는 `face_plane` 으로 면의 자리를 읽어 프로파일을 놓는다(kit 의 «`f.center` 근처에 그려라» 안내).

X축은 `any_perpendicular` — **가장 작은 성분의 축과 외적**, Onshape 의 perpendicularVector 와 같은 규칙이다. 가장 작은 두 성분이 **동률일 때 불연속**이고(연속인 선택은 수학적으로 불가능), 그 타이브레이크는 `nacre-math`에서 테스트로 못 박혀 있다.

**남는 한계**: 면의 *모양 자체*가 바뀌면 면적중심도 움직인다(면에 매인 어떤 규칙도 그렇다). 그리고 프레임은 f64다 — 축이 정규화를 거치므로 **유리수 법선에 수직인 단위벡터는 일반적으로 무리수**이고, 면 위 스케치의 정확성은 원점이 아니라 **축**이 벽이다.

**디버그 뷰어는 커널 크레이트가 아니라 워크스페이스 밖 별도 앱이다.** 연산 로그를 입력받아 매 동작을 스텝별로 재생하며(append-only라 "N번째까지 replay"가 공짜), STEP에 안 담기는 nacre **내부 정보**(정점의 정의와 그 실현 캐시·`Handle` 관계·인접 등)까지 시각화하는 인터랙티브 도구. 내부 자료구조에 접근해야 하므로 nacre를 **직접 링크**한다(개발 중 path 의존 → 안정화 후 version 의존, 버전별 디버깅도 자연스러워짐). 뷰어 없이 하는 시각 확인: 정상 결과는 STEP→step-loupe(구조+검증), 중간·깨진 상태는 OBJ 덤프→맥 미리보기.

`Store`/`Handle`은 **최하위 `nacre-store`에 둔다.** typed-index 저장소는 기하·위상을 전혀 모르는 순수 인프라이므로 두 층보다 아래에 격리하고, 위의 크레이트(topo의 아레나들·ops·validate·step·tess·props)가 자유롭게 참조한다. 핸들은 아레나 항목을 이름 짓고, 곡면 아레나는 **진실**(`Handle<MotionNode>`를 든 topo 타입)을 들므로 topo 아래 크레이트는 그 타입을 이름 지을 수 없다 ⇒ geom은 `Handle`을 영구히 갖지 않는다. (라이선스는 MIT/Apache-2.0 듀얼 — Manifold(Apache-2.0) 알고리즘 차용과 호환.)

`nacre-exact`는 **유리수·정수로 정확히 답할 수 있는 모든 것**을 격리한다 — 값(`Rat`·`Angle`)·이름(`PlaneName`·`MeetPoint`)·정확 술어, 그리고 ℚ 를 벗어나되 여전히 정확한 `QuadVal = a + b√c`. 사용자가 입력한 치수·각도를 f64 오차 없이 정확히 보존한다(`1.1`→`11/10`, `1.1×7`=정확히 `7.7` — "얇은 막" 문제의 근본 해결). `Rat`은 `Ratio<i128>` + **checked 산술**이다 — 오버플로는 조용히 감기지 않고 드러나며, 넘친 유도는 임의 정밀도 정수로 다시 한다(`plane_name_exact`; `PlaneName` 은 `Narrow | Wide`). 정점의 정의는 아레나가 불변 보존하고, f64 좌표는 그 실현 캐시다. `Angle`은 유리수 deg를 mod-360 **정확 누적**(한 바퀴가 정확히 0으로 닫힘 → 스케치 닫힘)하고 90°계열은 exact 유리수 cos/sin(회전 tol 0). **`nacre-predicates`와 상보(겹침 아님):** predicates는 기하 행렬식의 **부호를 exact 결정**(exact-부호), nacre-exact는 **입력 값과 유리수-순수 누적을 exact 보존**(exact-값) — 역할이 갈려 이름·층이 분리된다. **외부 의존:** `num-rational`(+`num-traits`·`num-bigint`·`num-integer` — 뒤의 둘은 `num-rational` 기본 피처가 이미 끌어오는 것을 명시한 것)과 `astro-float`. 성숙한 checked 유리수+gcd 약분을 제공하고, exact 유리수를 손수 구현하면 버그가 exactness 목표를 훼손하기 때문이다(MIT/Apache·순수 Rust). 어떤 `nacre-*`에도 의존 않는 **의존 그래프 최하단 순수 토대**. **회전 좌표의 toleranced 부호 판정**(무리수 좌표라 exact 못 하지만 부호는 f64 필터→astro-float 상승→*증명된 일치 또는 정직한 미결*로 sound하게 정함)은 **`nacre-judge::kernel`** 이 맡는다 — topo 아래의 순수 술어층, predicates의 쌍둥이(「판정」 절).

`nacre-predicates`는 **발견된 교차점의 부호 판정(내/외·orientation)을 좌표가 아니라 implicit point(정의)째로 하는 indirect predicates**를 격리한다(「해석 기하 층」 절의 정밀도 분업). 바닥의 적응 정밀 확장 산술은 `geometry-predicates`(MIT/Apache) 재사용, 그 위 implicit point 표현과 indirect 술어만 자체 구현. Rust 최초의 오픈소스 indirect predicates가 되도록 **nacre 밖으로 떼어낼 수 있게**(MIT/Apache 단독 공개 가능) 설계한다 — 그래서 어떤 `nacre-*` 에도 의존하지 않는다. 라이선스 엄수: 구현 참고처는 논문(Attene 2020, arXiv 2105.09772; Shewchuk 1997; Lévy PCK)과 `geometry-predicates` 소스로 한정하고, LGPL인 Attene 참조 구현 소스는 **작성 중 열람 금지 / 완성 후 실행 대조만 허용**.

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
- **캐시**(핸들 인덱스 병렬, 버리고 재생 가능): 정점 좌표 · 곡면 실현 · 간선 곡선 · 역방향 인덱스(「위상 층」 절)
- **순수 가속기**(진실도 캐시도 아님 — 비어 있어도 항상 옳다): 모션 사슬 접두의 고정밀 값 표

타입 이름·필드·공개 여부는 그 절들을 본다.

**Model은 "진실"만 담는다 — tess도 ops도 Model의 필드가 아니다.** 모델은 연산 로그의 재생 **결과**이므로 `ops: Vec<Operation>`은 Model을 *만들어내는 입력*이지 Model 안에 든 게 아니고, tessellation은 Model에서 *뽑아낸 파생 캐시*다(진실/캐시 분리). 게다가 레이어링상 Model은 `nacre-topo`에 사는데 `Tessellation`·`Operation`은 topo보다 위 크레이트라, Model에 필드로 넣으면 topo→tess/ops 순환이 된다. 그래서 op 로그와 tess 캐시는 상위 레이어(ops/파사드)에서 Model과 **나란히** 보관하고, "동일 로그·동일 cfg → 동일 모델 + 함께 재생성되는 tess"라는 커플링은 그 상위 번들이 책임진다.

(점 저장소는 두지 않는다 — 두 Vertex가 한 Point를 공유하는 상황은 설계상 존재하면 안 된다. `Vertex` 는 자기 **정의**만 들고, 좌표와 그것에 대해 증명된 것은 인덱스 병렬 캐시(`Model::vertex_cache`, 좌표 조각은 `Model::vertex_point`)에 산다.)

### 편집 연산의 supersede 의미론 — live 도달가능성 (partially persistent)

편집 연산(boolean·transform·mirror)이 기존 위상을 바꿀 때, append-only라 옛 셀을 **지우지 못한다**. 그래서 Model은 **`live_solids: Vec<Handle<Solid>>`(살아있는 솔리드 목록)를 진실로 보유**하고, 편집 연산은 옛 셀을 남긴 채 새 셀을 push한 뒤 **live 목록이 새 결과 Solid만 참조하도록 갱신**한다(소비된 입력 Solid는 목록에서 빠진다). **"살아있는 모델"의 정의 = live_solids에서 하향 참조로 도달 가능한 셀의 폐포(reachable closure).** supersede된 옛 셀은 아레나에 남되 어떤 live solid도 안 가리키므로 자동으로 "안 보인다".

**`Mirror` — 반사는 루프 역순으로 표현한다.** `Operation::Mirror { solid, axis, offset }`는 좌표평면 `axis = offset` 기준 반사다(**3개 방향 × 임의 위치**). 길이를 보존하고 **손잡이를 뒤집는다** — `Isometry`가 proper motion 전용이라 회전·이동의 어떤 조합으로도 못 만드는 유일한 강체 운동이다("−1배 스케일"은 미러와 스케일을 뭉개는 트릭이고 커널에 스케일은 없다). `transform`처럼 입력을 supersede하므로 원본을 남기려면 `copy`와 짝짓는다.

**축 정렬로 한정하는 이유는 회전과 같다** — 커널의 회전축도 X/Y/Z뿐이다. 임의 평면 반사 `x − 2((x−p)·n)n`은 무리수 법선 때문에 좌표가 근사가 되고, 그러면 모션 포리스트에 tol을 가진 반사 노드가 필요해진다. **임의 축 회전과 임의 평면 미러는 같은 항목**이며 그 설계가 들어올 때 함께 열린다(계획).

**구현의 핵심은 방향 대수다.** 반사는 det = −1이라 `R(a) × R(b) = −R(a × b)`이므로, 좌표만 반사하면 루프의 감김이 함의하는 법선이 뒤집혀 **솔리드가 안팎으로 뒤집힌다**. 그래서 모든 루프(outer + inner)를 **역순으로 감고**(`Loop::reversed` — 순서 역전과 `forward` 반전을 함께 하므로 모서리 짝이 유지된다), `Orientation` 플래그는 **건드리지 않는다**: 반사가 내적을 보존하므로 `sign(plane.normal · n_out)`이 그대로다. `nacre-geom`의 `Plane::mirrored`는 단위 법선의 한 성분 부호만 바꾸므로 비트 그대로 단위이고 다시 정규화하지 않는다. **`validate`의 그물은 절반이다**: 면 단위 orientation↔감김 불일치는 `FaceMisoriented`가 잡지만(반사는 루프 역순과 법선 반사가 함께 가서 이 검사에 불변), 플래그와 감김이 **함께** 뒤집힌 전-셸 반전은 면 단위 검사가 원리적으로 못 보므로 그쪽 그물은 부호 있는 부피와 불리언, 그리고 OCCT 오라클이다.

**반사도 모션이다 — 사슬에 들어간다.** `Motion::Mirror { axis, offset }`이 사슬 노드다. 기록 여부는 다른 모션과 **같은 규칙**을 따른다(「모션 사슬 — `MoveNode`」의 옮겨 적기 규칙) — `2c − x` 는 모든 유리수 `c` 에서 유리수라 새 솔리드의 반사는 진술로 옮겨 적힌다.

**반사는 비고유(`det = −1`)라 지름길에 패리티 규약이 필요하다.** exact 지름길들은 *"모션이 행렬식을 보존한다"* 를 근거로 **사전-모션 프레임**에서 답한다. 반사는 그것을 뒤집고, 게다가 실패가 보수적이지 않다 — 네 점이 같은 사슬을 지니면 행렬식이 **균일하게** 뒤집혀 확신을 가진 오답이 된다. 자리마다 인자를 유도하는 대신 **base 데이터의 손잡이를 한 번 맞춘다**: 사슬의 반사 수가 홀수면 base 점들에 정확한 반사를 한 번 더 걸어 건넨다(`shared_base`, `BaseFrame::of` 두 곳뿐). 그러면 base와 moved가 **고유 모션으로 관계**되어 어떤 행렬식 질문이든 등식이 그대로 선다. **base 평면의 부호는 그 전역 부호까지 센다** — 외적이 유사벡터라 반사는 삼각형의 회전을 뒤집고 평면의 기록은 안 뒤집으므로, 평면을 따로 반사하면 부호가 하나 남고 `plane_pair_dir_sign`이 바로 그 부호를 읽는다. 그래서 f64 로 파생할 때는 보정된 삼각형에서 파생하고, 기록(이름)을 쓸 때는 반사한 기록의 방향을 σ(`witness_name_sense`)의 반대로 잡는다. 예외는 `cmp_coord`의 지름길 하나뿐이다: 그것은 행렬식이 아니라 **축별 비교**이고 반사는 그 축을 뒤집으므로, 반사가 든 사슬에서는 끈다(보수적 미스).

**미구현(정직 거절):** 곡면(원/원통)의 반사 — 반사가 매개화 손잡이를 뒤집으므로 곡면 기하의 결정 사항이다. `MirrorNotPlanar` 로 거절한다.

**예외 하나 — `copy`는 supersede하지 않는다.** `Operation::Copy { solid }`는 솔리드의 독립 쌍둥이를 만들고 **원본을 live 목록에 남긴다** — 즉 **live 목록에 추가만 하는 유일한 연산**이다(다른 모든 편집 연산은 "옛 것 빼고 새 것 넣기"). 필요한 이유: `transform`·`boolean`이 입력을 live에서 빼므로, copy 없이는 같은 공구를 두 번 쓰거나 원본을 남긴 복사본을 두는 것이 불가능하고, 편의 레이어의 패턴·미러가 성립하지 않는다.

**구현:** `transform::transform_solid`(결정적 깊은 복제 walker: 곡면 → 정점 → 간선 → 면 → 셸 → 솔리드)를 **항등 아이소메트리**로 부르고 원본을 live에서 빼지 않는다 — `transform`과의 차이가 그 한 줄의 부재뿐이다. 위상 셀은 **반드시 복제**한다(두 live solid가 모서리를 공유하면 `edge_uses`가 면 4개 사용으로 읽어 manifold 검사가 깨진다). **곡면은 복제되지 않는다**: 항등 모션은 기록할 노드가 없고 진술이 비트 동일하므로, 같은 진술을 다시 push 하면 interning 이 **원본의 `Handle<Surface>`** 로 돌려준다(같은 진술은 두 핸들이 되지 않는다). 정점의 정의는 복제 후에도 살아남고, 좌표 캐시는 옮겨 받는 것이 아니라 새 정점의 정의를 실현해 채운다. 비-live 입력은 거절한다(소모된 핸들의 부활을 막고, 소비자 버그를 드러낸다). 이름은 **copy**(Handle의 clone은 얕은 참조 복사라는 반대 뜻).

**결과로, 모델을 소비하는 모든 코드는 store 전체가 아니라 live 도달가능 집합만 순회·카운트한다** — `validate`의 Euler/manifold, `Adjacency`(「위상 층」 절), `nacre-props`(부피·면적), `nacre-step` export, 오라클. 편집 중에는 아레나를 정리(compact)하지 않는다 — Handle 인덱스 안정성(=replay 결정성)을 유지하기 위해.

**대안 대비 (왜 도달가능성인가).** (a) *통째 복사*(편집마다 모든 셀을 새 독립 Solid로 복제) — 안 바뀐 셀까지 복사해 **정점 공유(single-reference)를 깨고** 메모리·오염을 부른다. 탈락. (b) *tombstone*(무효 Handle 집합을 손으로 관리) — 표시 누락 시 조용히 오염. 도달가능성은 **"아무도 안 가리키면 죽음"이라 누락 실수가 구조적으로 불가능**. 채택. 이는 partially persistent data structure(과거 버전은 읽기만, 최신만 수정)의 표준 CAD 적용이며, 이 절 위의 "소비돼도 안 지우고 결과 Solid가 새 항목 참조" 원칙의 구체화다.

**상용과 의도적으로 다른 선택 — 근거.** Parasolid·ACIS·OCCT 등 상용 커널은 40년간 "기하 불변 + **위상 가변**(Euler operator로 in-place 편집)"이 정설이다. nacre는 "기하 불변 + **위상도 불변**(append-only supersede)"으로 간다. 가변 위상은 dangling 참조·transient topology numbering·persistent naming 문제를 낳고 그걸 길들이는 데 40년 케이스워크가 필요한데 우리에겐 그 경험이 없다. 불변 위상은 그 버그 부류를 **구조적으로 제거**해 경험 없이도 안전하고, undo·결정론·fuzzing이 공짜이며, Rust 소유권 모델과도 정합적이다(가변 위상 그래프는 borrow checker와 싸워야 해 악명 높게 어렵다). 대가는 메모리인데 2020년대엔 작다. 즉 "상용을 이기는 방식"이 아니라 **"우리 제약(경험 없는 개인+AI 프로젝트·Rust·결정론 중시)에 최적화된 방식"**이다.

**성능·전제 (구현은 프로파일링 후).** 매 순회의 도달가능성 계산이 병목이면 generation/epoch 태그나 reference counting을 **파생 캐시**로 도입할 수 있다(진실=도달가능성, 이건 캐시 — Adjacency가 "진실=위상store, 캐시=역인덱스"인 것과 같은 층 분리). ref-counting은 하향 참조가 acyclic이면 죽음 판정에 바로 쓸 수 있다. **전제**: 도달가능성·ref-counting 둘 다 "Solid→Shell→Face→Loop→Edge→Vertex 부모→자식 단방향, 순환 없음"을 요구한다 — 「위상 층」 절의 위상 참조가 이를 만족하고(모든 셀이 하향 Handle만 보유), 유일한 역방향 인덱스 `Adjacency`는 진실이 아니라 캐시라 무관하다.

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
    /// 전체 원 carrier(중심/법선/**ref_dir**/반경 — `ref_dir` 이 θ=0 앵커이고, 원통의
    /// seam 은 **+ref_dir** 에 고정된다). 호 = 이 담체 + 서로 다른 두 끝점 정점,
    /// 온전한 원 = 같은 정점이 양 끝(`[v, v]`, 위상 층 절의 rim 규칙) — Line 과 같은
    /// carrier-vs-trim. step-io CurveInput::Circle(끝점 동일 여부로 원/호 구분)과 정합.
    Circle(Circle),
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

**정점의 캐시는 태어날 때부터 실현값이다.** 연산이 정점을 push 하는 자리 다섯(불리언 셋·압출·이동)은 전부 ops 의 깔때기 `push_vertex_realized(model, def, fallback)` 를 지난다: 정의를 사다리의 첫 두 단(128비트, 못 정하면 256비트)에서 실현해 `PointCache::Bounded { coord, bound }` 로 넣는다. 못 넣었을 때 서는 것은 «구성의 폴백»이 아니라 «이름 붙은 이유»다: 깔때기는 `realize_cache` 의 `Err` 를 망라적으로 갈라 `Ceiling`(첫 단의 비트가 모자랐거나, 재생할 이력이 **비용 한계** `CACHE_REPLAY_COST_CAP` 을 넘어 아예 안 걸었거나 — 둘 다 «더 물으면 답이 있다») 과 `Unrealized`(도로가 없다)로 나누고, 좌표만 구성 자리의 값에서 가져온다. 그래서 «실현 통로는 하나»가 구조다 — 표시·tess·STEP·물성·validate 가 읽는 f64 는 realize 의 메모다. **쓰기 문은 둘**이다: 짓는 `push_vertex` 와, `Ceiling` 만 올리는 `Model::refine_vertex_cache` — 그것을 부르는 것은 내보내기 직전의 «비싼 문» `nacre_ops::refine_caches` 하나다(정점·곡면·간선을 차례로 올리고 남은 것을 이유별로 보고한다; 로그 끝의 문이다 — 정점 캐시는 뒤 연산의 진실을 정하지 않지만, 평면 캐시는 뒤 불리언의 seam 표가 `Ceiling` 정점에 쓰는 구성 수치라 그 표의 거절이 달라질 수 있다). 계약은 census 가 행마다 단언한다: 답 ⇔ 변종, 값은 비트 동일. 거절 인구는 이름을 든다(`NoMeet` 등).

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

**`Model` 은 진실만 든 집합체다**: 정확 기하의 진실(`surfaces`·`motions`), 위상 store 다섯
(`vertices`·`edges`·`faces`·`shells`·`solids`)과 `live_solids`, 그리고 진실에서 유도되는 곁표들
— 인덱스 병렬 캐시 셋(`surface_cache`·`edge_cache`·`vertex_cache`), 역방향 인덱스 `adj`, interning 표,
가속기 `prefix_hp`. tess 와 ops 는 들지 않는다.

**정확 기하는 정점만이 아니라 «면»에도 적용된다.** 비-90° 회전은 평면의 계수를 **무리수로** 만든다 — f64 로 저장된 `Plane{origin, normal, raw}` 은 그 순간 진실이 아니라 **반올림된 상**이고, 그 사실을 적을 곳이 없으면 커널은 기본값으로 "정확하다"고 답하게 된다. 그래서 **곡면의 진실이 자기 `motion` 을 든다**: 유리수 진술은 `motion` 이 이름 짓는 틀에서 말하고, 모션 이력이 그것을 세계로 나른다(`None` = 세계).

**정점은 자기 정의이고 그것뿐이다.** `Vertex` 는 정의의 enum 이다(`ThreePlane`·`OnSeam`·`Pierce` —
변종마다 자기 진실을 말하고 불변도 변종별이다). 좌표와 «그 좌표에 대해 증명된 것»은 정점 store 와
인덱스 병렬인 점 캐시가 든다(`PointCache` — `Bounded { coord, bound }`·`Ceiling`·`Unrealized`;
`Model::vertex_cache`/`Model::vertex_point` 로 읽고 `Model::push_vertex` 가 채운다 — 정의가 먼저,
좌표가 둘째). 캐시의 변종은 캐시가 좌표에 대해 **아는 것**을 말하지 정점이 어디서 왔는지를 말하지
않는다. 간선의 캐시와 달리 버리고-재생이 없다: 정점은 push 될 때 실현되고, 뒤에 일어날 수 있는 것은
`Ceiling` → `Bounded` 의 **정련**뿐이다.

**간선의 진실은 «담체 둘 + 끝점 둘»이다.** 실현된 곡선은 필드가 아니라 곁의 캐시(`Model::edge_cache`,
`EdgeCache`)다 — 담체와 끝점이 그것을 정한다(`Model::derive_edge_curve`). `Model::push_edge` 가 캐시를
함께 채우고, `Model::rebuild_edge_cache` 가 통째로 버리고 재생한다. 값 가운데 topo 가 진실에서 못 내는 것 —
세계 이름 없는 평면이 낀 직선의 방향, 세계 진술 없는 원통의 림 중심(그 사슬을 재생해야 한다) — 은 유도가 자기 길로 답하지 못할 때만
push 하는 쪽에 묻는다(`EdgeGiven`; ops 깔때기 `push_edge_realized` 가 실현해 답한다). «진실이 여기서
답하나»는 그래서 topo 한 곳에서만 정해진다. **곡선의 종류도 진실이 정한다** —
평면 × 원통 간선이 원(평면이 축을 가로지름)인지 룰링(축을 따름)인지는 평면 이름의 법선과 원통 `def.dir()` 의
정확한 관계(`nacre_exact::AxisRelation`, `Model::plane_cylinder_relation`)이고, 비스듬하면(타원) 간선이 없다
(`EdgeDecline::Oblique`). 캐시가 드는 것은 그 곡선의 값(중심·프레임)뿐이라 재생이 종류까지 같고, 간선 곡선의
종류를 읽는 판정은 진실을 읽는다. 관계는 두 방향의 물음이고 강체 운동이 보존하므로 두 진술을 세계까지 옮길
필요가 없다: 각 방향을 자기 사슬의 선형부로 정확히 옮길 수 있는 만큼 옮기고(`LineDir` — 이동은 그대로, 거울·
사분각은 치환, 회전은 자기 축을 그대로, 프레임은 `ẑ` 를 그 평면의 법선으로), 남은 운동이 **값으로** 같거나 한쪽이
더 갔으면 그만큼 되돌려 한 프레임에서 묻는다. 원통 인구 게이트(`cyl_gate`)와 원의 클래스(`class_carries_circle`)가
같은 규칙(`nacre_exact::axis_relation`)을 부른다.

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

**모델이 들고 다니는 tol 은 없다.** 구성된 점이든 불리언이 발견한 점이든 정점은 태어날 때 정의에서 실현되어 `Bounded` 경계를 들고, 저장된 잔차는 없다. 그렇다고 **tolerance 가 사라지는 것은 아니다** — 검사는 그대로 돌고(위 `EPS_CONSTRUCTED`), 사라지는 것은 **모델이 들고 다니는 tol** 이지 비교가 아니다: Fornjot 의 "무-tolerance 단순함"은 «저장된 tol 이 없다»는 뜻으로만 재현된다.

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
뒤는 전부 **live 도달가능 셀** 위에서, `model.adj` 를 믿지 않고 **새로 지은** `Adjacency` 로 돈다:

- **루프 닫힘** — 각 루프의 half-edge 가 끝→시작으로 이어진다(`OpenLoop`).
- **다양체** — 모든 간선은 정확히 두 번, 반대 방향으로 쓰인다(`NonManifoldEdge`·`NonOpposedEdge`);
  정점 둘레도 다양체다(`NonManifoldVertex`).
- **담체 일치** — 정점의 정의가 이름 짓는 곡면(`VertexCarrierMismatch`)과 간선이 진술한 담체 쌍
  (`EdgeCarrierMismatch`)이 그것을 실제로 쓰는 면들과 맞는다.
- **방향** — 면의 감김과 진술한 바깥 법선이 맞고(`FaceMisoriented`), 공동 껍질은 안쪽을 향한다
  (`CavityMisoriented`).
- **기하 접속** — 캐시된 점이 자기 정의의 곡면(`VertexOffDefinition`)·간선의 곡선(`VertexOffCurve`)·
  면의 곡면(`VertexOffSurface`) 위에 tol 안으로 앉아 있다. 원통은 진실과 f64 캐시가 맞는지도 본다
  (`CylinderTruthCacheMismatch`) — 세계에 진술된 원통의 캐시는 문이 그 진술에서 실현하므로, 진술을 생산자가 따로
  옮긴 값과 대조하는 일은 문(`push_cylinder_raw` 의 `debug_assert`)이 덮어쓰기 전에 한다.
- **오일러-푸앵카레** — `V − E + F = 2(S − G) + L_i`(`L_i` = 면의 내부 루프 수, `S` = 껍질 수,
  `G` = genus). **기준은 `V−E+F=2` 가 아니다**: `validate` 는 **χ = V − E + F − L_i** 를 세어 **짝수인지**
  (`EulerParity`)와 **G = S − χ/2 ≥ 0** 인지(`NegativeGenus`)를 본다. 뚫린 솔리드는 genus 1 이라 χ = 0
  이고, 커널이 실제로 그것을 만든다. rim 이 `[v, v]` 인 것은 이 셈이 맞게 하기 위해서다.

`Adjacency`는 진실이 아니라 캐시라는 점이 중요하다 — 위상 store들이 진실이고, 인덱스는 **버리고 재생**한다. `adj` 에 쓰는 자리는 `Model::rebuild_adjacency` **하나**이고 그것은 통째 교체다(`edge_uses`/`vertex_edges` 를 증분으로 건드리는 코드 0줄). 소비자가 **일괄 추가 뒤 손으로 한 번** 부른다(*"the cache is otherwise stale"*). 증분 갱신은 **필요가 증명될 때** 짓는다.

**Adjacency는 store 전체가 아니라 live 도달가능 셀만 인덱싱한다(supersede 의미론).** `rebuild`는 `model.faces`/`model.edges` 전체가 아니라 `live_solids`에서 도달 가능한 면·엣지만 순회한다 — 그래야 supersede된 옛 면의 엣지가 `edge_uses`를 오염시켜 manifold 검사(엣지 정확히 2회)를 깨뜨리지 않는다. 위상 참조가 하향 단방향·비순환이라 도달가능성 순회는 유한·안전하다. `Adjacency`는 **위상 셀의** 역방향 인덱스이고 캐시이므로 도달가능성 진실을 침해하지 않는다. interning 표 넷(`motion_ids`·`interned` — 키 → `Handle`)은 «값으로 핸들을 찾는» 표라 성격이 다르고, **`Adjacency` 와 달리 도달가능성으로 걸러지지 않는다.** 어느 것도 순회하지 않는다(`HashMap` 의 순서가 결과에 닿으면 안 된다).

**곁표에는 «진실도 캐시도 아닌» 셋째 부류가 하나 있다 — 가속기.** `Model::prefix_hp` 는 모션 체인이 이미 접은 **접두사**의 고정밀 값을 들고, 다음 모션이 base 부터 다시 접는 대신 거기서 이어 접는다(이것이 없으면 깊이 n 인 솔리드를 짓는 데 1+2+…+n 이 든다). **캐시와 다른 점**: 캐시는 진실의 메모라 «있으면 읽는다»가 계약이지만, 이것은 비어 있어도 답이 **비트까지 같다** — 없으면 정의에서 다시 접을 뿐이다. **interning 표와 다른 점**: 핸들 정체를 지지 않으므로 **언제든 비울 수 있다**(`clear_prefix_hp`). 열쇠가 정의 `(base, leaf, prec)` 라 **낡을 수 없고**(무효화 규칙이 없다), **적중이 곧 축출**이라 용량이 이력 길이가 아니라 **살아 있는 꼭짓점 수**에 묶인다. 값은 「접은 체인 노드 수 + 점」이다 — `Motion::Frame` 하나가 여러 노드로 펼쳐지므로 parent 걸음 수로는 접미사를 자를 수 없다. 쓰기는 «접두사 사슬을 잇는» 자리(`transform`)에서만 일어나고, 불리언의 결과 정점은 표를 건드리지 않는다.

## Tessellation 층 (`nacre-tess`) — 출처 태그 파생물

메시는 진실이 아니라 파생물이고, 모든 요소가 출처를 안다. **길은 하나다**: `tessellate(model, cfg)` →
`Tessellation`. 간선을 한 번 샘플해 공유 폴리라인으로 두고, 모든 면을 자기 곡면의 차트에서 삼각분할한다.
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

**닫힌 간선의 폴리라인은 첫 점을 반복하지 않는다 — 닫힘은 암묵이다.** 온전한 rim 은 `0 .. (n−1)τ/n` 에서
샘플하고 멈춘다(핸들을 반복하면 면이 퇴화 삼각형을 얻는다). 폴리라인을 **쌍으로** 걷는 소비자는 `n` 각
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
정점, 스파이크, 서로 지나가는 두 선분; 구멍이 외곽 링을 가로지르는 것도 잡는다), `SelfTouchingBoundary`
(아래), `HoleWinding`(구멍이 외곽과 같은 방향으로 감겼다 — 고칠 메시가 아니라 깨진 솔리드), `OverBudget`
(아래). `tessellate` 는 모델 전체를 걷다 첫 거절에서 멈춘다.

**삼각분할 알고리즘 — 브리징 없는 스윕.** 구멍 있는 면은 **y-단조 분해 + 단조 삼각분할**(de Berg §3)로 자른다. 구멍을 외곽 링에 브리지(폭 0 슬릿)로 꿰매 귀 자르기를 하지 않는 이유는 하나다: 브리지가 정점을 **반복**시키고, 반복 정점이 있는 링은 단순 다각형이 아니며, Meisters 의 two-ears 정리는 단순 다각형에만 성립하므로 그런 링에는 **귀가 아예 없을 수 있다**. 스윕은 **링을 절대 합치지 않으므로** 퇴화 링 자체가 안 생긴다 — 구멍 간선이 보통 간선이고, 구멍의 최상단이 split·최하단이 merge 정점이 되어 대각선이 자동으로 구멍을 잇는다. 판정은 전부 부호이고 **tolerance가 없다**: 전부 exact `orient2d`(`nacre-predicates`)이고, 이것은 **어떤 f64 입력에도 정확**하다. (평면의 차트는 축 드롭이라 `uv`가 원본 f64 그대로지만, 원통의 차트는 `atan2`로 **계산된** 값이다 — 그래도 술어가 정확하므로 알고리즘은 같다. 면 사이에서 `uv`를 비교하는 일이 없으므로 — 면끼리 만나는 자리는 **엣지 폴리라인**이 이미 정한다 — 차트가 면마다 달라도 무방하다.) 경계 정점을 추가하지 않으므로 삼각형 수는 `V + 2H − 2`로 불변이고 crack-free 규칙도 그대로다. 스윕은 자기 전제를 **가정하지 않고 검사한다**: b-rep 이 보장하지만 이 층은 보장할 수 없는 것이라, 자기 교차 링을 받은 분해는 그러지 않으면 확신에 찬 틀린 메시를 낸다.

**경계가 자기에게 닿는 면 — 접촉 자리의 공유 간선 다리.** 구멍이 다른 링에 정확히 한 점에서 접하는
평면 면은 유효한 솔리드지만(면이 집혔을 뿐 곡면은 그 점에서 2-다양체다 — 집힌 면의 두 엽이 이웃 곡면을
돌아 이어진다) 두 링으로는 분해할 수 없다: 닿는 정점은 **곡선** 간선의 샘플이고, 그것이 내려앉은 직선
간선은 두 끝만 든다. `tessellate` 의 사전 패스(`bridge_shared_edges`)가 그 **이미 있는** 정점을 직선
간선의 폴리라인에 끼워 넣고, 두 링을 그 점에서 하나로 이으며, 스윕은 생겨난 일치 쌍을 기호적으로
순서 짓는다(`polygon::sos`). `by_edge` 가 공유 구조라 이웃 면도 같은 정점을 보므로 crack-free 가 유지된다
(면이 혼자 경계를 바꾸는 것이 아니라 **간선의 폴리라인**이 바뀐다). 접촉은 어떤 폴리라인도 바꾸기 전에
**모든** live 면에서 먼저 모으고, 한 간선의 삽입은 뒤 위치부터 한다; 패스는 멱등이다. 간선→면 표는
모델의 adjacency 가 아니라 live 집합에서 여기서 짓는다(불리언 안에서는 adjacency 가 낡아 있다).
이 길이 덮지 않는 모양은 고치지 않고 **세어서** 거절한다(`Declined`: `CurvedEdge`·`CurvedNeighbour`·
`SharedSegment`·`MultiSegment`·`NoEdge`; 보고는 `bridge_report` → `BridgeReport`) — 다른 링의 **정점**에
닿는 접촉, 닿은 선 위에 정확히 놓인 이웃, 다리가 둘 이상인 면도 그렇다. 그것들이 `SelfTouchingBoundary`
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
경계에서는 구성으로 지켜지지만(아래 「예산 둘」), 면의 **내부**는 자기
경계가 우연히 대 준 샘플을 물려받는다. 경계가 곡률 방향을 아크로 덮는 동안은 공짜지만, 덮지
않으면 **조용히 틀린다**: 병합된 옆면(지워진 가짜 이음매가 한 구간의 유일한 아크였던 면)은 내부 점
없이는 반 바퀴를 한 번에 건너는 현으로 그려져 **넓이의 1/6**을 잃는다.

그래서 규칙과 공급을 **두 층으로** 나눈다.

- **규칙**(일반·자기 검증): *모든 내부 간선은 선언된 두 예산 안의 현이다* — 선형은
  `Surface::distance(3D 현의 중점) ≤ tol`, 각도는 두 끝의 **곡면 법선 사이 각** `≤ max_angle_deg`.
  둘 다 `Surface`에 전수 대응이라 **평면은 저절로 통과**하고(거리 0·법선 동일) 새 곡면은 컴파일
  에러다. 넘으면 `TessError::OverBudget` — 고치지 않고 **정직하게 거절**한다(틀린 캐시는 없느니만
  못하다). 경계 간선은 동어반복이라 보지 않는다: 아크는 `circle_segments`가 두 항의 **최댓값**으로
  나누고, 직선 엣지는 거리 0이고 양 끝 법선이 같다.
- **공급**(곡면을 아는 층): 원통 차트가 **자유 파라미터 없는 격자**를 내놓는다. 곡률 방향은
  샘플러의 **바로 그 함수**(`circle_segments`)의 칸으로 — 그래야 경계가 이미 깔아 둔 아크와 정렬이
  공짜다 — 그리고 `i`는 차트의 **풀린** `v` 범위 전체를 덮는다(구멍은 온전한 회전만큼 옮겨 다닌다).
  축 방향은 **경계 «원 엣지»마다 그 원의 중심의 축 좌표** 하나: 곡선에서 읽지 샘플된 점에서 읽지
  않는다(한 림의 점들은 그 좌표가 ulp만큼 다르므로 중복 제거하면 림 하나가 181개가 된다).
  **최소·최대 스테이션은 버린다** — 그 둘은 면 자신의 두 림이고 거기 찍은 점은 경계 점의 1-ulp
  복제다(남기면 보통 밴드의 삼각형 수가 움직인다). 평면 차트는 아무것도 대지 않는다.

삽입은 후보마다 **놓을 자리**를 찾는다 — 삼각형의 **엄격한 내부**(셋으로 쪼갬)이거나 **내부 간선
위**(둘로 쪼갬). 두 번째가 예외가 아니라 **주된 길**이다: 후보는 면의 모양이 바뀌는 선 위에 앉고
삼각분할은 이미 그 선을 따라 간선을 갖고 있다. **그러나 두 길 다 프로덕션에서 산다** — 네 벽과 코너에 선
보스의 옆면 다섯(노치가 구멍인 밴드)에서 내부/간선이 `−x` 37/141 · `+x` 162/16 · `−y` 32/148 · `+y` 82/96 ·
코너 147/121 이다. 그 비율은 사용자가 그린 모양이 아니라 밴드를 자른 생성선의 자리가 정하므로 어느 갈래도
도달 불가로 지울 수 없다.
**제약 간선은 절대 안 쪼갠다** — 이웃이 모르는 점은
T-정점이고 그건 크랙이다. 자리가 없는 후보는 **버린다**(강제하지 않는다).
`ChartMap`이 되돌리는 길이다(차트는 전단사다) — 발명된 점은 위치를 읽어 올 모델 정점이 없고, 두 팔
모두 정방향이 건 **손잡이 보정을 되돌려야** 한다.

**품질은 직교하는 별도 패스다.** 단조 삼각분할은 정확하지만 슬리버가 많다. 자유(비제약) 간선에 **Lawson 플립**을 돌려 제약 Delaunay로 보내면 최소각이 최대화된다(정리 — Lawson 1977, Chew 1989). 새 의존성 0(`incircle`이 이미 있다). 경계는 제약이라 안 뒤집히므로 watertight·개수·면적 전부 불변. 평균 최소각은 2배가 된다(원통 캡 tol 1e-4에서 0.51°→1.05°). **다만 캡의 *최악*값은 안 움직이고, 그건 알고리즘 탓이 아니다** — 정n각형 정점은 **공원**이라 `incircle`이 전부 0을 답하고(모든 삼각분할이 Delaunay), 게다가 내각 178.9°와 two-ears가 0.57° 이하 슬리버를 **강제**한다(측정 0.457° = 하한). 고치려면 점을 **추가**해야 하는데, 그건 *캡의* 얘기다 — 캡의 경계는 이미 곡률을 다 샘플하므로 더 넣을 이유가 **품질뿐**이고, 품질은 이 층의 결정이 아니다. (금지되는 것은 «점 추가» 자체가 아니라 **경계** 정점의 추가다 — 위 「차트가 점을 댄다」를 보라.)

**곡선을 얼마나 쪼갤지는 예산 «둘»이 정한다.** 원의 분할 수는
`n = max(⌈360°/Δθ⌉, ⌈π/acos(1 − tol/r)⌉)`이다. 뒤의 항(절대 sagitta)만 쓰면 `n ≈ π√(r/2tol)`이라
**작은 원일수록 적게 쪼갠다** — 절대 오차 예산은 지키지만 반지름 0.2가 10각형이 되고, 확대하면 그대로
드러난다. 앞의 항(각도)은 `s/r = 1 − cos(π/n)`에서 반지름이 사라지는 **상대** 예산이라 크기와 무관하게
같은 품질을 준다. 정n각형에서 **중심각 = 선분의 꺾임각 = 인접 면 법선의 벌어짐**이 모두 같으므로 이 한
숫자가 윤곽과 음영을 동시에 묶는다(OCCT의 angular deflection과 같은 개념). 기본값 `Δθ = 2°`(=180분할,
반지름의 0.0152%), `tol = 1e-2`. 각도 항은 반지름 ≈66까지 이기고 그 위는 sagitta 항이 이어받는다 —
**어느 반지름에서도 한쪽만 쓸 때보다 나빠지지 않는다.**

**틈 없음(crack-free) 규칙:** 면의 삼각분할은 반드시 `by_edge`의 공유 polyline 정점들을 자기 경계로 소비해야 한다. 인접한 두 면이 각자 독립적으로 엣지를 샘플링하면 공유 엣지에서 정점이 어긋나 T-junction이 생기고 watertight가 깨진다 — 엣지 polyline이 먼저, 면 삼각분할이 그걸 경계 조건으로.
**금지되는 것은 «경계» 정점의 추가다.** 면이 자기 **내부**에 점을 넣는 것은 이웃이 볼 일이 없으므로 크랙과 무관하고, 위 「차트가 점을 댄다」가 바로 그것이다. 폴리라인을 더 잘게 하는 것도 원리상 금지가 아니다 — `by_edge`가 공유 구조이므로 **엣지 샘플러가** 촘촘히 내면 양쪽 면이 함께 본다(OCCT도 `BRepMesh_ModelHealer`에서 그렇게 한다; 위 접촉 다리도 같은 원리다). 금지는 **면이 혼자** 경계를 바꾸는 것이다.

## 진실과 캐시 — 타입 구조

아레나는 append-only 이므로 **진실만 들어간다.** 캐시는 곁표에 두어 언제든 버리고 재생한다.
정의가 불변이므로 **캐시는 낡을 수 없다** — 무효화라는 개념이 없고, «버리고 재생»만 있다.


## 숫자 규칙 — 일곱 개

1. **진실은 «유리수»이거나 «핸들»이다.** `Rat`(i128)에 들어가는 값은 적고, 안 들어가는 값
   (발견된 좌표 160–480비트, 계수 곱)은 저장하지 않고 **가리킨다**. 가리키는 사슬의 끝은
   언제나 수다(C2·C3 의 동시 만족).
   **폭은 사용자가 쓴 기하와 그것이 놓인 자리가 정하지 연산 횟수가 정하지 않는다** — 불리언은
   평면을 새로 만들지 않고 피연산자의 평면을 재사용하며(`a_boolean_mints_no_surface`), 모션은
   진술로 옮겨 적으면 닿은 자리의 이름(그 자리에 직접 지은 것과 같은 폭)이 되고 기록하면 계수를
   건드리지 않으며, 오프셋은 법선을 고정한 채 d 만 옮기므로 분모가 곱이 아니라 lcm 으로
   자란다(`1.1 + 6.6 = 7.7`). `nacre-ops/tests/instruments/point_width.rs`
   (`a_stacked_name_is_as_wide_as_the_place_it_lands`)가 이것을 잠근다.
2. **화살표는 한 방향뿐이다.** `Rat → f64`(캐시), `Rat → BigFloat`(판정 상승). `f64 → 진실`은
   없다 — f64 를 들어올려 진실로 삼는 순간 반올림이 진실에 구워진다. 유일한 입구는
   `Rat::from_decimal`(들어올리기가 아니라 사용자가 *쓴* 십진수의 결정적 정규형)이다.
3. **무리수는 값이 아니라 정의로 존재한다.** 회전은 `Angle`(유리수 degree), 기울어진 프레임은
   `Motion::Frame`(법선을 이름으로)이 들고, cos/sin/√ 는 실현 시점에 정확 반올림으로 나온다.
   **모션은 면이 든다** — 정점·모서리는 면을 가리켜 정의상 따라온다.
4. **실현 통로는 하나다.** `realize(정의, bits)` — f64 캐시(표시·tess), STEP 출력
   (`realize(128) → round_to_f64` 정확 반올림), 판정 상승이 전부 같은 통로를 탄다. 읽어 내는 규칙도
   하나다(`Realized::to_f64`): 최근접 f64, 또는 일치 정밀도(`max(1, 점의 최대 |좌표|)·2⁻¹⁸⁰` — 점의 크기는
   모델 크기를 넘지 않으므로 규칙 5 보다 엄격한 쪽) 안이 증명된 좌표는 `+0.0` — 참값 0 은 반올림으로는 어느
   정밀도에서도 결정되지 않는다. 반올림보다 먼저 물으므로 어느 단에서 결정해도 답이 같다. 캐시는
   언제나 **값 + tol 쌍**으로 함께 만들어지고 함께 버려진다. **f64 좌표는 즉석에서 다시 풀지
   않는다 — 표에서 읽는다**(근거: 판정마다 정의를 다시 만들면 그것이 시간의 77%, 25초 대 0.4초).
5. **판정은 3단이고, 결과는 3갈래다.** f64 필터 → (세 평면이 모션을 공유하면) 정수/`Expansion`
   정확 경로 → astro-float 상승. 결과는:
   - 부호가 **증명되면** → 그 부호(`Sign`).
   - 부호는 못 갈랐지만 **«일치 정밀도보다 가깝다»가 증명되면** → 일치로 처리하되 근거를 실어
     **보고**한다(`Decision::Coincident{within}`, `boolean_with_report`). 일치 정밀도
     (`Standard::coincidence` = 모델 크기 × 2⁻¹⁸⁰ — f64 출력 해상도 2⁻⁵² 보다 두 워드 아래)는
     **모델에서 유도**되는 값이지 사용자 설정이 아니다 — 낮게 잡으면 비트만 더 쓰고 높게 잡으면
     진짜로 떨어진 것을 합치는 비대칭이 아래로 밀며, 그 아래의 갈라짐은 어떤 출력에도 살아남지
     못한다. 전역 tolerance(*"가까우면 붙여라"*)와 방향이 반대(*"가깝다고 **증명되어야** 일치"*)다.
   - 둘 다 증명 못 하면 → **이름 붙은 거절**(`JudgeExhausted`·`DegenerateWitness`·
     `PrecisionBudget`). 조용한 0은 없다(C7).

   정밀도는 상수가 아니라 모델이 정하고, 필요한 비트를 계산해 **한 번에** 점프한다.
6. **동일성은 정준 이름의 `==` 다 — 이름의 그릇은 임의정밀이다.** 이름은 하나만 저장한다:
   `PlaneName = Narrow([Rat;4]) | Wide([BigInt;4])`, i128 에 들어가면 반드시 `Narrow`(정규화
   불변식). 동일성은 enum 전체의 `==` 이고, 산술·프레임·지름길은 `narrow()` 투영을 빌려 읽는다
   (`plane_name_exact` 가 중간 계산을 임의정밀로 닫으므로 폭이 이름을 잃게 하지 못한다).
   좁은 형태(법선·`n·n`)가 없다고 **진실을 거절하지 않는다** — 필요하면 점 셋에서 임의정밀로
   실현한다. 이 규칙의 «항상 이름을 갖는다»는 유리수 닫힘인 평면에만 적용된다 —
   세 정점이 한 프레임에서 만나지 않는(혼합 프레임) `Through` 평면은 이름이 없다 — 이름은 한 프레임의 유리수
   풀이에서 나오는데 그런 프레임이 없다(세계 계수는 대개 무리수지만, `z = 0` 처럼 유리수여도 그렇다). 그런 평면의
   저장은 **진술 키**(정렬 삼중항+모션, `SurfaceKey::Through`)로 intern 되어 «같은 진술 = 한 핸들»은 지켜지고,
   같은 평면이 이름 있는 핸들과 이름 없는 핸들 둘이 될 수 있다. 기하 동일성은 술어가 답한다: 클래스 발견이
   핸들이 다른 평면 쌍마다 `planes_coplanar` 를 물어(정확 0, 또는 일치 정밀도 안이 증명되면 보고와 함께) 한
   클래스로 합치고, 결과의 기하는 한 핸들일 때와 비트 동일하다 — 느릴 뿐이다(release 로 잰 세 픽스처에서 3–9배;
   `a_plane_stated_twice_booleans_as_it_does_once`).
7. **`Rat` 을 넓히지 않는다(C8).** 넓어지는 것은 중간값(BigInt 이름 유도·`Expansion` 술어·
   BigFloat 상승)뿐이다. 규칙 1 의 근거가 그대로 이 규칙의 근거다 — 저장되는 값의 폭은 사용자가
   쓴 기하가 정하고 연산을 쌓아도 자라지 않으므로, 넓은 그릇이 필요한 자리는 저장이 아니라 계산이다.

---

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
/// 단일형이 아니라 **세 변종**이다: 원통 seam 정점은 «세 평면»으로 적을 수 없다.
/// 불변식은 변종별이고 그것이 이 구조의 근거다.
pub enum Vertex {
    /// 세 평면의 교점 — 이름이 곧 점. (D != 0 은 좌표 재생이 생기는 자리에서 단언한다.)
    ThreePlane([Handle<Surface>; 3]),
    /// 두 곡면의 교차 «곡선» 위의 점 — 원통 seam(테두리 원의 θ=0). 쌍은 곡선을 박고, 그 위의
    /// 점을 고르는 것은 원통 진실의 `ref_dir`(`CylinderDef`)이다 — 그래서
    /// `OnSeam([cylinder, cap])` 은 「rim ∩ +`ref_dir` 광선」으로 **정확히 지정된 한 점**이다.
    /// 정의에서 좌표를 재생하는 기계는 `realize_vertex`(`seam_point`)다.
    /// 계획: 재생한 값을 캐시에 되쓰는 것 — 그것이 없는 동안 STEP 은 만들 때 나온 f64 를 내보낸다(todo).
    OnSeam([Handle<Surface>; 2]),
    /// 이 점은 도법기하의 **관통점**(piercing point)이다 — 두 평면의 교선이 원통 옆면을 뚫는 자리.
    /// 담체 종류를 **구조로** 말한다(닮은 핸들 셋이 아니라). 교선 방향은 저장 순서(오름차순
    /// 핸들)의 **정준 이름 법선**으로 `ℓ = n₀ × n₁` 이고, `root` 는 ℓ 을 따른 오름차순 매개값,
    /// 접선(중근)은 `QuadRoot::Double` 이다. 두 평면 핸들을 다시 정렬하는 쪽은 `root` 를 함께
    /// 다시 말해야 하고, 그 규칙이 사는 곳은 `QuadRoot::canonical` 하나다.
    Pierce {
        planes: [Handle<Surface>; 2],         // 오름차순 핸들
        cylinder: Handle<Surface>,
        root: QuadRoot,                       // Lo | Hi | Double — 선언 순서(Lo < Hi)가 곧 규약
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
    /// **가리키는** 변종(`origin: Handle<Vertex>`)을 추가한다(규칙 1 의 연장 — 계획, todo).
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
| **원통 원시체는 없다 — 원은 스케치다** | 원통은 원 프로파일의 Extrude 다. 직선–호 꼭짓점 = `Vertex::Pierce`(두 평면의 **저장된 정준 이름** × 원통, 근은 매개 값으로 고른다 — 빌더가 자기 점으로 계수를 다시 만들면 `Lo/Hi` 가 뒤집힐 수 있다), 온전한 원 = `OnSeam` |
| **경계는 f64, 진실은 구성 시점에** | 공개 API 는 f64 그대로다. `Rat::from_decimal` 왕복을 `Profile2d` **생성자**에서 한다 — `check()`(자기교차·중첩·포함)가 진실 위에서 정확 술어로 돌고(`orient2d_rat`: narrow 우선 → BigInt 전역 부호, geom 의 `_rat` 워커 쌍둥이), 십진 창(1e38/1e-22) 밖 치수가 구성 시점의 이름 붙은 에러다(`ProfileOutsideDecimalWindow` / sketch 층 `OutsideDecimalWindow`). f64 부호 ≠ 십진 부호다(0.1·0.2·0.3 공선이 이진에선 굽음) — 그래서 check 가 진실 위에서 도는 것은 정확성의 문제이고, 방향은 항상 "작성자가 쓴 수가 이긴다". 생성자는 `check()` 를 돌리지 않는다(O(n²)) — 프로파일을 소비하는 **모든 연산이 먼저 돌린다**; `from_rings`·`from_paths` 가 검사하는 생성자다 |
| **공선 중간점은 생성자가 지운다** | 공선 정점의 양옆 벽은 한 평면 → 교차가 직선이라 정의 불가, 그리고 비-2-manifold 퇴화다. 제거는 형상 불변(무손실 정규화, 거절 아님 — 정리된 프로파일의 프리즘은 깨끗한 쌍둥이와 비트 동일, 모든 코너가 3-평면 정의 보유). 판정은 Rat 위 정확 orient2d, 제거 조건은 **엄격 내부**(중복점·스파이크는 생존해 각자의 이름 붙은 에러로 보고된다 — 경계 포함 판정을 재사용하면 작성자의 실수를 조용히 지운다) |
| **프레임은 평면에서 뜨지 않는다** | `Named` 의 `origin` 은 참조 평면 **위**의 점(정확 검사 — C1). `SketchFrame::named` 가 구성 시점에 검사한다 — `plane_residual_sign`(scalar, **전역**: Narrow 는 Rat 대입 → 넘치면 BigInt, Wide 는 BigInt — fail-open 없음) ≠ 0 이면 `OriginNotOnPlane`, `WideFrame::named_of` None 이면 `RefDirParallelToNormal`(그 함수는 폭에 전역이라 None 은 평행/영벡터뿐), 십진 창 밖은 `FrameOutsideDecimalWindow`, 이름 없는 평면은 `PlaneWithoutExactForm`. 좌표계에 주의: «정의하려는 프레임의 (u,v,w)» 가 아니라 **그 평면의 `points` 가 적힌 좌표계**의 3D 점이다(모션 없으면 세계 — 상자 윗면 z=1 이면 origin 은 `[1,1,1]` 같은 점이지 셋째 성분 0 이 아니다). "평면 위" 는 성분이 아니라 **방정식 대입**(`a·x+b·y+c·z+d = 0`, 정확)으로 검사한다. `Canonical` 은 유도라 검사할 것이 없다. 면 위 스케치의 밑캡이 대상 면의 surface 핸들을 공유하는(flush 접촉 = 핸들 비교) 전제이기도 하다. 평면에서 d 떨어진 스케치가 필요하면 origin 을 띄우는 것이 아니라 **오프셋 평면**을 만든다 — 같은 프레임 안 `(0,0,d),(1,0,d),(0,1,d)` 유리수 세 점의 `Known` 평면(datum 가족). 담체가 실제 평면 핸들로 남아 이름·interning·flush 규칙이 그대로 성립한다 |
| **정확 유리수 프레임은 노드를 만들지 않는다** | 프레임의 세계 기저가 유리수이면(`exact_frame` — 평면 이름에서 `RatFrame::of_plane_frame` 으로 묻고, 평면에 모션이 있으면 그 사슬의 접기로 옮긴다(`RatFrame::carried` — 이동·사분각·축 거울은 접히고, 프레임 노드·사분각 밖 회전은 안 접혀 거절). 세계 축 평면 셋이 전부 여기 든다, `Canonical` 이든 ZX 의 `Named` 든) `Motion::Frame` 노드를 생성하지 않고 세계 유리수 산술로 구성한다(모션 = None) — 빈 사슬 지름길(`shared_base` 상쇄·축별 tol 0·정수 Shewchuk)이 그대로 산다. 거울이 든 사슬로 옮긴 기저는 왼손이라 `RatFrame` 이 손방향을 들고 `normal()` 은 옮긴 ŵ 를 돌려준다(쓸기·오프셋이 ŵ 쪽으로, 감김·호의 회전은 노드 도로가 사슬 parity 로 하듯 그 손방향으로). 변환의 옮겨 적기 규칙(「모션 사슬 — `MoveNode`」)과 같은 구성 시점 정규화다. 지반 잠금: 씨앗 평면의 canonical 프레임이 단위축 위에 **비트 정확**으로 실현된다(`a_seeded_planes_canonical_frame_is_the_world_basis_exactly`) — 생략은 중복 제거지 값 손실이 아니다. 실현한 축을 `from_decimal` 로 되들어올리는 게이트는 정규화가 필요한 축(`0.6000000000000001`)과 짧은 십진이 아닌 원점을 틀리게 읽는다(「가지 말 것」) — `SketchPlane::exact` 는 호출자가 *쓴* 축을 읽는 테스트 오라클로만 남는다. |
| **세계 축 평면 셋은 `Model::new()` 가 심는다** | 핸들이 결정적(0·1·2, 법선 축 순서 Z(XY)·X(YZ)·Y(ZX))이라 replay 가 자명하다. 캐시 방향은 **−축**(밑캡의 감각), `Default` 는 `new()` 로 위임(무씨앗 뒷문 없음). 씨앗은 정의상 영구 orphan — 정의·프레임만 가리키는 평면을 순회·직렬화가 따라가야 한다(todo) |
| **스케치는 모델에 저장되지 않는다** | 연산 로그가 스케치의 진실이고, 모델에 남는 것은 벽·캡 평면(`Known` 점 + 프레임 모션)이다. 치수의 흐름: 변 `(x1,y1)→(x2,y2)` + 스윕 → 벽의 세 점 `(x1,y1,0),(x2,y2,0),(x1,y1,d)` — 산술 없이 유리수 그대로. 계수(이름)에만 곱이 있고, 그것이 진실이 점이어야 하는 이유다(계수 `i128` 적합 25.8% vs 점 100%) |
| **기본 프레임은 정준 유도(`Canonical`)다** | 원점 = 세계 원점의 수선의 발(부호·스케일 불변), 축 = Arbitrary Axis(DXF) — 사용자가 아무것도 안 정하면 이 규약이고, 규약은 스펙+테스트로 동결한다. 사용자가 명시한 프레임만 `Named` 로 **값을 저장**한다(유도값과 다른 규약을 유도로 흉내내면 규칙 변경에 조용히 돈다) |

연쇄가 닫힌다: 평면에 스케치 → 솔리드 → 그 면 위에 스케치 → … 매 단계가
`SketchFrame{plane: 이전 단계의 평면}` 으로 **사용자 조작당 1단**만 깊어지고(C6 — 불리언은
평면을 만들지 않는다, `a_boolean_mints_no_surface`), 좌표는 끝까지 유리수이며, 무리수는 프레임
실현(`1/√`)에만 산다.

---

## 이름과 interning

이름은 **하나만 저장**하고(`Model::surface_name` 곁표), 좁은 형태는 저장이 아니라 **투영**이다:

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
    /// Shewchuk 정확 술어(`Expansion` 조각)·프레임 유도·지름길이 읽는다.
    /// `None`(Wide)이면 그 지름길만 못 타고 판정이 일반 경로로 — 느릴 뿐 틀리지 않는다.
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
하나(`intern`: 열쇠·이미 발급된 핸들과 `flipped` 회신·종류별 아레나 push·곁표·계측)를 부른다. 열쇠가 읽는
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
- **큰 정수 계수의 정확 부호**: `2^53` 은 f64 *하나*의 한계일 뿐 — `i128` 은 `Expansion` 3조각
  (53비트씩)으로 정확하고, 정규화 실패한 유리수는 술어 안에서 분모를 턴다. 분포: 1조각
  71.9% · 2조각 26.1% · 3조각 1.2% · `i128` 초과 0.9%(임의정밀 이름이 받는다). f64 필터가
  먼저 걸러 expansion 은 애매한 ~1%에서만 돈다.

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
    adj: Adjacency,                           // 요청 시 재구축(`rebuild_adjacency`)
    prefix_hp: HashMap<PrefixKey, PrefixValue>, // 모션 사슬 접두의 고정밀 실현 메모(실현 가속기)
    motion_folds: Vec<Option<AxisAffine>>,      // 모션 핸들 인덱스 평행 — 노드가 태어날 때 사슬 전체를 접은 유리수 사상

    // 평면의 정준 이름 — 유도되는 값(캐시)인데 `SurfaceCache` 안이 아니라 곁표로 산다.
    // `Model` 의 유일한 공개 필드다.
    pub surface_name: HashMap<Handle<Surface>, PlaneName>,

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

// ── 경계 지어진 점 — 원자 셋(축별). 실현(`realize`)의 출력과 `WitnessPoint.realized` 가 이 모양이다.
//   값과 오차를 «묶어» 든다(따로 든 배열 둘이 아니라). 그래야 «참값 ∈ value±error» 가 구조로 서고
//   값만 옮기고 오차를 안 옮기는 desync 가 불가능하다.

// ── 정점의 캐시 ──
pub enum PointCache {
    Bounded    { coord: Point3, bound: [Mag; 3] },                //   정의에서 실현 — 참값 ∈ coord ± bound
    Ceiling    { coord: Point3 },                                 //   싼 도로가 멈췄다 — 비싼 문이 답할 수 있다
    Unrealized { coord: Point3 },                                 //   실현이 뒤에 없다 — 도로가 없거나, 안 물었거나
}
// 연산이 만든 정점은 push 시점에 실현되므로(`push_vertex_realized`) 거의 전부 `Bounded` 다
//   (census 결과 정점 4,972/4,972). 모션 이력이 갈린 담체의 코너는 담체마다의 증인 삼각형으로 실현된다
//   (담체마다 `surface_witness_triangle` → `nacre_judge::plane_hp`, 세 평면의 교점 `meet_hp`, 한 묶음 동안 평면은
//   `PlaneMemo` 가 한 번만 실현한다); 그 도로는 증인 점 아홉을 재생하므로 비용 한계가
//   그 합을 센다.
// **`Ceiling` 의 천장에는 벽이 둘이고, 둘은 같은 말을 한다** — (a) 사다리의 두 단(128·256비트)이 못 정했다 — 0 이 아닌
//   좌표는 비용 한계 안에서 256 이 정하므로(한 회전에 약 1비트, 192 + 53 < 256) 인구는 참값 0 이 사분각 밖
//   회전 수십 번보다 깊이 묻힌 좌표다(일치 정밀도 2⁻¹⁸⁰ 까지 반경이 줄지 않는다 — 회전당 약 1비트에서 예상 76번,
//   37° 왕복 픽스처에서 60번은 결정·80번은 못 함),
//   (b) 재생해야 할 이력이 비용 한계(`CACHE_REPLAY_COST_CAP` = 192)보다 깊어 싼 도로가 아예 안 걸었다(접기가
//   답하는 사슬은 깊이와 무관하게 읽는다). 읽는 쪽에도
//   되찾는 문(`nacre_ops::refine_caches`)에도 뜻이 같다: «더 물으면 답이 있다» ⇒ 한 변종이 맞다.
//   `Unrealized` 는 그 반대 — 도로가 없다(증인 삼각형이 없는 담체 `NoMeet`, 회전된 원통 위 관통점 `NoCurvedPoint`). 그 갈림은 사람이 읽을 것을
//   위해 있다: «이 모델이 비싸다»와 «커널에 구멍이 있다»는 다른 보고다. 변종을 고르는 것은 호출자가
//   아니라 깔때기(`push_vertex_realized`)이고, 철저한 `match` 라 새 이유가 옛 이름에 조용히 접히지 않는다.
// **어디에도 잔차를 싣지 않는다.** 잔차는 담체까지의 거리 «하나»이고 좌표가 참값에서 얼마나인지를
//   말하지 않으므로 «증명된 것»이 아니다 — 배열의 seam 표도 캐시(`SeamVertex.cache`)를 들고, 자기접촉 체는
//   `Bounded` 의 경계로 상자를 부풀리며 경계가 없는 변종은 무한 상자로 읽는다(후보를 버리지 않고 정확 판정에
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

// ── 곡면의 캐시 — geom 의 실현을 감싼 구조체. `EdgeCache` 의 거울이다. ──
pub struct SurfaceCache { realized: nacre_geom::Surface, standing: CacheStanding } // Plane(..) | Cylinder(..) 는 geom 의 enum
pub enum CacheStanding { Unrealized, Ceiling, Realized }          // 캐시가 값에 대해 아는 것 — 순서가 앎이고, 오르기만 한다
//   topo 는 감싸기만 하고 기하의 메서드(`distance`·`normal_at`·`translated` …)는 `nacre_geom::Surface` 에
//   남는다 — 다형 질의는 geom 의 enum 이 안에서 가른다. 감싸는 이유는 캐시가 자랄 자리가 있는 이름 붙은
//   것이 되게 하려는 것이었고, 거기서 `PointCache` 의 세 변종과 같은 말이 자랐다(경계를 들 것이 없어 표지만):
//   `Unrealized` 는 아무도 실현하지 않았다(안 물었거나 도로가 없다 — 혼합 프레임 `Through` 는 증인 삼각형이
//   없다), `Ceiling` 은 싼 도로가 멈췄다(비용 한계를 넘는 사슬·두 단에서 미결정 — 값을 치르는 문이 올린다),
//   `Realized` 는 모든 부분이 진실의 최근접 f64 다. 둘째 쓰기 문 `refine_surface_cache` 는 표지를 내리지 못한다.
//   평면의 f64 실현은 문이 진실에서 유도한다(`derive_surface_cache`) — 세계 이름이 있으면 단위 법선은 이름 ×
//   향의 정확 반올림(`nacre_exact::unit_vector_f64`, 참 0 성분은 `+0.0`), 앵커는 진실의 첫 점(`Through` 는
//   꼭짓점들이 만나는 프레임의 첫 만남을 평면의 사슬로 옮긴 점; 진실이 그 점을 유리수로 못 놓는 넘침·넓은 만남은
//   받은 앵커). 문이 유도하지 못한 것은 ops 깔때기(`push_plane_realized`·`push_plane_through_realized`)가 push 한
//   뒤 진실에서 실현해 올린다(`raise_surface` → `Model::refine_surface_cache`) — 문이 먼저 정하므로 «진실이 여기서
//   답하나»는 topo 한 곳에 있다. 앵커는 첫 점을 사슬로 재생해(넓은 만남이면 첫 꼭짓점 자신의 실현), 법선은 모션
//   전 이름 법선을 사슬로 옮겨(이름이 `Rat` 보다 넓으면 진술의 세 점을 재생해 그들이 펼치는 법선으로 — 판정의
//   증인이 프레임 탐침이면 그 순서는 진술의 것이 아니라 쓰지 않는다), 둘 다 같은 읽기 규칙으로. 답하지 않는 반쪽은
//   생산자 값이고 표지가 그것을 말한다: 재생이 비용 한계를 넘는 사슬·256 에서 미결정은 `Ceiling`, 혼합 프레임
//   `Through` 는 `Unrealized`. 원통도 같은 두 단이다 — 세계에 진술되면 문이 원점·단위 축·단위 `ref_dir`·반지름을
//   정확 반올림하고(`cyl_unit_frame_f64`), 아니면 깔때기(`push_cylinder_realized`)가 사슬을 재생해 같은 넷을 실현하며,
//   비용 한계를 넘는 사슬은 생산자 값에 `Ceiling` 이다. 사설
//   `push_plane_raw`·`push_cylinder_raw` 가 진실과 캐시를 같은 자리에서 채운다.
//   **평면의 향은 진실이다**(`Surface::Plane.sense`): 면의 바깥 = 평면의 향 × `Face.orientation`, 둘 다
//   진실이다. 생산자는 향을 정확히 진술하고(캐시에서 읽지 않는다 — `f64 → 진실` 이다), `flipped` 는 두
//   진실의 비교다. 캐시의 향은 진실을 따른다: 세계 이름이 있으면 법선 자체가 이름 × 향이고,
//   이름 없는 평면은 문이 점들의 f64 외적으로 맞추며(`align_cache_sense`), 회전·프레임이 있는 사슬에서는
//   census 잠금(`audit_plane_senses`)이 지킨다. 원통은 기준이
//   진실 안에 있다(축에서 바깥으로). `fallback` 인자가 드는 것은 실현값(앵커·단위 법선)뿐이다.
//   **향을 묻는 판정은 전부 이 진실을 읽는다**: 평면이 보는 쪽 = `sense` × 점들의 방향(이름과의 관계는
//   `Model::plane_name_sense`·`world_plane_name_sense`), 증인 삼각형의 순서 = `sense × orientation`, 두 면의
//   관계 = `face_facing`. 캐시 법선은 어느 방향 판정에도 들지 않는다. (「가지 말 것」의 «평면의 방향을 점
//   순서에서 뽑기» 와 어긋나지 않는다 — 그 행이 버린 것은 점 순서 **단독**이고, 여기서는 `sense` 가 그 관계를
//   진술한다.)
//   평면의 이름(`PlaneName`)은 `plane_name_exact(세 점)` / `plane_name_through(세 정점)` 로 언제든 다시
//   나오므로 캐시다 — 그런데 곁표 `surface_name` 에 따로 살아, 한 곡면의 캐시가 두 벌이고
//   그릇(`Vec`↔`HashMap`)·이름 규칙(`_cache`↔`_name`)·가시성(비공개↔`pub`)이 셋 다 어긋나 있다.

pub struct EdgeCache     { curve: Curve }                          // 평가 가능한 담체 곡선
//   담체와 끝점이 곡선을 정한다(`derive_edge_curve`). `rebuild_edge_cache` 가 통째로 버리고 재생한다.
//   topo 가 진실에서 못 내는 조각은 유도가 push 하는 쪽에 묻는다(`EdgeGiven` — 세계 이름 없는 평면이 낀
//   직선의 방향, 세계 진술 없는 원통의 림 중심). ops 의 답은 담체 쌍마다 `given_by_pair` 에 남아 같은 쌍의
//   다음 간선이 읽고, 재생은 그것도 버려 다시 묻는다.
// **캐시는 «실현값 (+ 필요하면 그 오차)»** 다: 캐시가 드는 오차는 증명된 경계뿐이고
//   (`PointCache::Bounded{coord, bound}`), 오차 필드는 «소비자가 있으면» 붙지 대칭으로 붙지 않는다 —
//   `EdgeCache`·`SurfaceCache` 에는 오차 필드가 없다.
```

**두 «오차»를 섞지 말 것.**

f64 는 둘이고 둘은 만나지 않는다 — 모델의 캐시(표시·tess·내보내기용 실현값)와 판정의 작업
사본(연산 하나 동안만 산다). 표의 두 줄은 한 사다리의 두 단이 아니라 두 수명이다.

| 무엇 | 성질 | 어디 사나 |
|---|---|---|
| `Plane::distance_eps(p)` | **이 거리 계산**의 f64 반올림, `3ε·Σ|pᵢ−oᵢ|` — `p` 에 의존(두 연산자). 원통판은 `4ε·(축거리+반지름)` | 호출할 때마다 계산, **저장 불가** |
| `PointCache::Bounded.bound` | 실현의 축별 경계 — 좌표가 참값에서 얼마나(결정된 좌표는 반 ulp, 일치 정밀도로 0 이 된 좌표는 그 구간), **잔차가 아니다** | 모델과 같이 |

「곡면당 저장하는 tol 은 틀린 양」은 첫째 줄에만 맞다: 거리 계산의 반올림은 `p` 에 달려 있어 곡면의
성질이 아니다 — 평면 위에 정확히 있는 점도 0이 아닌 잔차를 낼 수 있고, 그 잔차가 주장된 허용오차와
**정확히 같아** 비교가 엄격하다는 것 하나로 통과하던 사례가 있었으므로 `distance_eps` 가 그 항을 따로 든다.
`Plane` 은 **원점 + 단위 법선**만 든다(`from_point_normal` 은 한 번 정규화하고, `from_point_unit_normal` 은 받은 단위 법선을 비트 그대로 싣는다 — 단위 벡터를 다시 정규화하면 비트가 움직인다). 판정은 이것을 읽지 않는다.

좌표가 결과인 STEP 출력은 정의에서 실현한 값을 정확 반올림해 낸다 — 정점의 캐시는 태어날 때부터 그
실현값이다(`push_vertex_realized`: 사다리의 낮은 단이 정하는 좌표는 어느 높은 단이 댈 답과 같아
`realize_vertex(…, NearestF64)` 와 비트 동일하고, `nacre-step` 은 `vertex_point` 를 읽는다). 자릿수를 더
원하는 호출자의 문은 `realize_vertex_decimal` 이다 — 캐시를 더 길게 찍으면 점이 아니라 반올림이 찍힌다.

| | 모델 캐시 (`PointCache`·`SurfaceCache`·`EdgeCache`) | 실현 (`WitnessPoint.realized`·`WorkingCyl.realized`) |
|---|---|---|
| 누가 만드나 | 정점: push 깔때기가 **정의에서** 실현한다(실현이 거절하면 구성 자리가 계산한 f64 가 선다). 평면: 세계 이름이 있으면 문이 진실에서 앵커와 단위 법선(정확 반올림)을 유도하고, 문이 유도하지 못한 반쪽은 ops 깔때기가 push 한 뒤 진실에서 실현해 올린다. 답하지 않는 반쪽엔 생산자가 계산한 f64 가 서고, 표지(`CacheStanding`)가 그것을 말한다. 원통: 세계에 진술되면 문이 진술에서, 사슬이 안 접히면 ops 깔때기가 사슬에서 원점·단위 축·단위 `ref_dir`·반지름을 정확 반올림한다(비용 한계를 넘는 사슬은 생산자 값에 `Ceiling`). 간선: 담체·끝점에서 유도 — 원은 원통 캐시의 프레임을 비트 그대로 싣고, 중심은 진술된 원통이면 진실의 만남(정확 반올림), 아니면 ops 깔때기가 사슬로 재생한 축과 캡 평면의 만남을 실현한 것(그 길이 답하지 않을 때만 캐시끼리의 f64 만남); 직선은 시작 정점을 지나 세계 이름 있는 두 평면이면 이름 법선의 정확한 외적의 최근접 단위, seam·룰링이면 원통 캐시의 축, 평면 하나라도 세계 이름이 없으면 ops 깔때기가 두 평면의 계수(세계 이름, 아니면 증인 삼각형 — `PlaneMemo`)의 외적을 실현한 최근접 단위를 따르고, 그 길이 답하지 않을 때만(비용 한계를 넘는 사슬·증인 없는 평면) 끝점의 차를 따른다 | 판정이 **정의에서** 정밀도를 불러 만든다 |
| 오차 | 정점만 축별 경계 셋(`Bounded` 변종일 때). 곡면·간선엔 없다 | 축별 경계 셋 — «참값 ∈ 값 ± 경계» 증명됨 |
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

## 문의 이름 — `Model` 의 접근자 명명 규칙

저장소가 비공개이므로 밖은 **문(=`Model` 의 메서드)** 으로만 읽는다. 규칙:

> **`x(h)` 는 «진실»(아레나 항목 그 자체), `x_cache(h)` 는 «캐시»(실현). 실체 이름을 단 맨 메서드가
> 캐시를 돌려주는 일은 없다.**

```rust
// ── 진실 — 실체 이름 그대로 ───────────────────────────────────────────
m.surface(h) -> &Surface      m.vertex(v) -> &Vertex      m.edge(e) -> &Edge
m.face(f) -> &Face            m.shell(s)  -> &Shell       m.solid(s) -> &Solid
m.motion(n) -> &MotionNode
m.plane_motion(h) -> Option<Handle<MotionNode>>            // 진실의 조각(곡면이 진술된 모션)
m.live_solids() -> &[Handle<Solid>]                        // 루트
m.world_plane(axis) -> Handle<Surface>                     // 심은 세계 평면

// ── 개수·핸들 복원 — `x_count()` · `x_handle_at(index) -> Option<Handle<X>>` ──
m.surface_count() · m.vertex_count() · m.edge_count() · m.face_count() · m.shell_count() · m.solid_count()
m.surface_handle_at(i) · m.vertex_handle_at(i) · m.edge_handle_at(i) · m.face_handle_at(i) · …

// ── 캐시 ─────────────────────────────────────────────────────────────
m.surface_cache(h) -> &nacre_geom::Surface                 // 포장을 벗긴 실현 — `match` 로 바로 가른다
m.surface_cache_standing(h) -> CacheStanding               // 그 캐시가 값에 대해 아는 것
m.vertex_cache(v)  -> &PointCache                          // 변종(Bounded | Ceiling | Unrealized)
m.vertex_cache(v).coord() -> Point3                        // 세 변종이 모두 드는 좌표
m.vertex_cache(v).bound() -> Option<&[Mag; 3]>             // Bounded 만 — 실현의 축별 경계
m.vertex_point(v)  -> Point3                               // = vertex_cache(v).coord() — 어디서나 답한다
m.edge_curve(e)    -> &Curve                               // 간선 캐시의 조각(곡선)
m.line_direction_cache([p, q]) -> Option<[f64; 3]>        // 두 평면이 만나는 방향 — 이름에서, 아니면 push 하는 쪽이 실현해 남긴 것
m.surface_name                                             // pub 곁표 — 평면의 이름(캐시)

// ── 세계 진술 — 진술된 프레임에서 세계로, 모션 사슬을 접어 정확히 ────────
m.world_plane_name(h) -> Option<PlaneName>                 // 평면의 세계 이름(사슬이 안 접히면 None)
m.chain_plane_coeffs(leaf, c) · m.chain_point_rat(leaf, p) · m.chain_dir_rat(leaf, d)   // 접은 사슬의 세 얼굴

// ── 쓰기 문 ──────────────────────────────────────────────────────────
m.push_plane(fig, standing, ..) · m.push_plane_through(..) · m.push_cylinder(..)   // 곡면은 진실을 진술하며 들어온다 — 그림과 그것이 아는 것을 함께
m.push_vertex(def, cache) · m.push_edge(.., given) · m.push_face(..) · m.push_shell(..) · m.push_solid(..)
m.push_motion(..)                                          // interned
m.refine_vertex_cache(v, coord, bound)                     // 캐시의 둘째 쓰기 문 — `Ceiling` 만 «올린다»
m.refine_surface_cache(h, s, standing)                     // 곡면 캐시의 둘째 쓰기 문 — 표지는 오르기만 한다
m.rebuild_edge_cache(given) · m.rebuild_adjacency()        // 통째 재생 — `given` 은 topo 가 못 내는 조각을 묻는다
```

**진실과 캐시를 이름으로 가르는 이유.** `Plane` 은 `Surface::Plane`, 즉 실체 이름이라 `Model` 에
`plane(h)` 같은 문을 달면 진실을 줄 것처럼 읽힌다. 캐시를 주는 문은 `_cache`(또는 캐시의 조각
이름 `_point`·`_curve`)로 끝나 **「캐시의」라고 말하고 들어온다**. 조각을 돌려받은 타입의 메서드로
꺼내는 것(`vertex_cache(v).bound()`)은 `nacre_geom::Surface` 가 `distance()`·`normal_at()` 을 자기
메서드로 갖고 안에서 가르는 것과 같은 관용구다 — 구·원뿔이 와도 `Model` 의 문은 늘지 않는다.

정제는 `Ceiling` 을 `Bounded` 로 올리는 것뿐이다 — 좌표를 더 정확하게 만들 수만 있다.

성능은 이 명명의 고려사항이 아니다. 어느 철자든 안에서 하는 일은 «디버그 가드 + 배열 색인 하나 +
판별자 읽기»이고 할당도 복사도 없다(`#[inline]`). 다만 평면의 이름과 실현이 둘 다 필요한 자리는 조회가
둘(`Vec` 색인 + `HashMap` 해시)이다.

## 판정 (연산 동안만 산다)

판정 타입은 두 크레이트에 나뉘어 산다. **`nacre-judge`** 는 정의를 든 점(`WitnessPoint`)·모션의
펼침(`MoveNode`)·증명 기준(`Standard`)·결과(`Decision`)와, 평면 표에 던지는 질문(`Judge<W: Witness>`)을
갖는다. **`nacre-ops`** 는 연산 하나 동안 모델에서 그 표를 짓는다(`pub(crate)` `FaceInfo`·`WorkingPlane`·
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
    pub(crate) world: Option<WorldName>,              // 같은 평면의 세계 이름과 그 향 — 좁은 투영 `world_rat()` 이 원통 도로가 세계 축과 비교하는 유일한 서술, 향을 곱한 `stored_world_rat()` 이 방향
    pub(crate) surf: Handle<Surface>,                 // 대표 곡면 — 결과 정점의 담체로 기록된다
    pub(crate) tri_pt3: [WitnessPoint; 3],            // 같은 증인의 정확한 정의 — 모든 술어가 빌린다
    pub(crate) rotated: bool,                         // 술어 경로 선택 신호. `tri_pt3` 와 함께 복사돼 어긋날 수 없다
    pub(crate) frame_sign: i8,                        // 라벨 프레임이자 차트의 프레임(n_out = frame_sign · n_P)
    pub(crate) base: BaseFrame,                       // 회전 이전 쌍둥이
    pub(crate) name_ints: Option<NameInts>,           // 저장 방향으로 접힌 이름 정수 + 53비트 안이면 그 f64 행
                                                      //   — 모션 없는 평면의 정확 지름길이 읽는 유일한 서술
}

// 판정용 원통 클래스 — 결정은 전부 `def`(정확)를 읽고, `realized` 는 «측정»만 읽는다.
pub(crate) struct WorkingCyl {
    pub(crate) surf: Handle<Surface>,
    pub(crate) def: nacre_topo::CylinderDef,
    pub(crate) realized: nacre_geom::Cylinder,        // 모델의 캐시가 아니다 — 연산 하나를 산다
    pub(crate) owner: SolidSide,                      // 쌍 루프는 주인이 다른 쌍만 묻는다
}
```

**판정 표의 증인(`Witness`)은 진실만 든다** — 세계 이름(`world_name`)과 정의 점(`tri_pt3`). 면 꼭짓점의 좌표
캐시는 표 행에도 트레이트에도 없어 판정이 읽을 수 없다. 증인 삼각형의 순서는 진실의 비트가 정하고(`sense` ×
면의 `Orientation`), 두 면의 바깥이 같은 쪽인지는 같은 곡면 → 향 비트, 같은 세계 이름 → 두 향, 그 밖엔 두 증인의
정의(`normals_agree_judge`)가 답한다(`face_facing`) — 클래스 병합과 같은 갈래 순서다. 평면 클래스 형태
(`PlaneWitness`)가 더하는 것도 정확한 서술뿐이다 — 이름의 행과 모션 전 프레임.

**판정층은 증인 삼각형 하나로 전체(total)다.** 판정 표의 계약은 «평면 위 세 정확한 점»(`Witness::tri_pt3`,
늘 있다)이고, 이름 없는 `Through` 평면은 자기 프레임의 probe 로 그 점을 정의상 갖는다 — 진실의 `Through`
는 판정층에서 변종이 아니라 **probe 로 유도된 증인 삼각형**으로 나타난다.

**정확 지름길은 평면을 이름으로만 서술한다.** 평면 캐시(원점과 단위 법선 — 그 음함수의 `d` 는 f64 곱 `n·origin`)와 면
꼭짓점의 좌표 캐시는 둘 다 반올림된 상이라 진실 평면에서 `2⁻⁵⁴` 떨어질 수 있고, 둘을 서로 대조해 맞는 것만
싣는 검사는 «반올림된 평면»을 인증할 뿐이다(`3·0.1 = 0.3` 인 모서리를 벽 밖으로 판정해 유효한 교집합을
거절했다). 그래서 모션 없는 평면의 `exact_coeffs`/`exact_normal` 은 이름(정의 점에서 반올림 없이 유도한 정준
정수)의 f64 행이고 — 53비트 안일 때만 선다 — `orient3d` 는 넷째 평면도 그 행으로 묻는다(`indirect_plane_side`).
한 질문의 네 평면이 한 어휘로 서술되므로 두 서술이 어긋날 자리가 없다. 모션이 있는 평면은 `BaseFrame` 이 같은
원칙(증인 점이 f64 로 정확할 때만, `d` 는 기록에서)을 지킨다.

```rust
// crates/nacre-judge/src/kernel/frame3/witness.rs

/// 증인 점 — 유리수 base 를 모션 사슬로 나른다. base + chain 이 정의, realized 는 f64 실현
/// (축마다 값+경계 한 원자 — 모델 캐시가 아니라 연산 하나의 것), hp 는 고정밀 메모.
/// 같은 정의는 같은 실현(경로 무관)이다. 동등성은 정의(base·chain)만 본다.
pub struct WitnessPoint {
    pub base: [Rat; 3],
    pub chain: HpRc<[MoveNode]>,              // 사슬을 펼친 것. 평면당 한 번 펼치고 증인 점 셋이
                                              // 같은 Rc 를 복제해 든다 — 두 번 펼치지 않으므로
                                              // 어긋날 수 없다. 점이 드는 이유: 자기완결 판정
                                              // (`shared_base`·실현은 점만 받는다).
    pub realized: [Bounded; 3],
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
pub struct Standard {
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

- **정확 경로의 조건 = 세 평면이 같은 모션을 공유할 것**(`None` 포함) — 그 프레임 안 정수
  계수로 Shewchuk. 모션이 섞이면 f64 필터 + 고정밀 상승(C4).
- 모션이 계수에 하는 일: `Translate` 는 `d' = d − n·t`(유리수면 정확), 축과 나란한 법선은
  `R·n = n`(구조적 검사 가능), 그 밖은 실현에 오차. 진실 쪽에는 오차가 없다 — 오차는 실현할 때
  생긴다.
- `Through` 평면의 판정은 두 갈래다.
  **유리수 닫힘이면 `Wide` 이름이 곧 정확 계수다 — 조건부로.** 세 술어(`orient3d`·`cmp_coord`·
  `dir_sign`)에 이름-정수 rescue 가 있다: 게이트 = 전 평면 이름 있음 ∧ 하나 이상 wide ∧ (전부 무이동
  **또는** 전부 한 사슬 — `cmp` 는 무이동만). 이름은 구성 시 저장 방향으로 σ-접혀(`name_stored_ints`)
  BigInt 부호 술어(`int_*`, scalar)로 정확히 답한다. **조건의 이유**: wide datum 의 이름은 세계
  이름(담체 = 발견 정점)이라, 그 평면을 처음 만드는 혼합-프레임 불리언에서는 일부만 잡히고
  (상승 454→330), 전부-세계/통째-이동 표(2세대)에서 완전히 공짜다(139 대 20 — 대조군인 narrow 보다
  싸다). 그릇은 `Expansion` 이 아니라 BigInt 다.
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
| C8 | `Rat` 을 넓히지 않는다 | 넓어진 것은 이름의 **그릇**(동일성 전용)과 중간값뿐 |

---

## 타입이 사는 자리 — 규칙

- 회귀 관문의 기본은 **위상 정확 일치 + 좌표 ε(모델 크기 상대, `2⁻⁴⁰` — 여유 2000배)**, 비트
  동일은 보너스 신호다. 합격/불합격이 아니라 **최대 편차 숫자를 찍는다**(누적 드리프트 감시).
- **비트 동일 관문은 그 안의 population 에 대해서만 보증한다.** 바꾸려는 것이 닿는 population 을
  대장에 먼저 넣고(17자리 `fw`·기울어진 `tp` 가족이 있다), *"어느 population 인가"* 는 추측하지
  말고 계측이 이름을 대게 한다. 좌표 관문은 «답은 같은데 더 나쁜 길로 갔다»를 원리적으로 못
  보므로 **"정확 경로를 탔는가"를 직접 단언하는 테스트**를 함께 둔다.
- **두 `Surface` 의 철자**: 한 파일이 둘 다 필요하면 **맨이름은 진실**(topo, `Handle<Surface>` 가
  이름 짓는 것)이고 실현은 `nacre_geom::Surface` 로 적는다. 반대로 하면 같은 파일에서
  `Handle<Surface>` 는 진실을, 맨 `Surface::Plane(p)` 는 캐시를 뜻해 한 단어가 표지 없이 두 가지가
  된다. 섞으면 **하드 오류**다(geom 은 tuple 변종, 진실은 struct 변종 ⇒ E0532) — 조용히 틀릴 자리가 0.
- **타입을 바꾸는 변경은 「타입 구조」 부의 그림을 같은 커밋에서 고친다.**
- **어느 크레이트에 두나** — 「유도된 값 + 그 산술」은 `nacre-exact`(최하단; `nacre-judge` 가 닿아야
  하므로), 「아레나 항목의 진실」은 `nacre-topo`. 구조적 강제: **`nacre-judge` 는 `nacre-topo` 에
  의존하지 않는다** ⇒ 판정이 쓰는 타입은 전부 topo 아래에 있어야 한다(`PlaneName` 이 scalar 에 사는
  이유). 모델에서 판정 표를 짓는 타입(`WorkingPlane`·`WorkingCyl`·`FaceInfo`)은 둘을 다 아는
  `nacre-ops` 에 산다. 본문은 「크레이트 구조」 절.
- **병렬 불변식**: 병렬 구간에서 **곡면 push 금지**(재생 결정성). 곡면의 문은 `Model::push_plane`·
  `push_plane_through`·`push_cylinder`(사설 `push_plane_raw`·`push_cylinder_raw` 로 모인다).
  워크스페이스의 진짜 병렬은 `nacre-ops/src/par.rs` **한 파일 두 자리**
  (`(0..n).into_par_iter().map(f)`)뿐이고, 그 클로저는 「인덱스 → 값」이라 **`&mut Model` 을 들 수
  없다** ⇒ 병렬 구간의 곡면 push 는 구조적으로 불가능하다. 그 파일에 `push_plane`/`push_cylinder`
  0건, 그리고 *"병렬이 두 곳에 살면 조용히 깨진다"* 는 텍스트 가드가 그 한 파일을 지킨다. 감시는 그
  가드가 한다.

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

```rust
pub enum Operation {
    // 모델이 **호출자가 이름 붙인 평면**을 드는 유일한 연산. 평면을 로그 **밖**에서 만들면
    // 그 모델은 자기완결적이지 않으므로(replay 가 핸들을 재현해야 한다) 연산이어야 한다.
    // 평면은 구성 시점 interning 대상이라 «이미 있으면 그 핸들»이 정답이다(아레나 안 자람).
    DatumPlane { def: DatumDef },  // Stated(SketchPlane) | Offset { frame, dist } | ThroughVertices([Handle<Vertex>; 3])
    // 스케치 평면·프로파일 개념은 Extrude 인자로 흡수된다(별도 Sketch op 없음).
    // 평면은 값으로 싣지 않고 `SketchFrame` 으로 «이름 부른다». `dist` 는 부호 있는 두께(음수는 ŵ 반대쪽), 0 은 거절.
    Extrude { frame: SketchFrame, profile: Profile2d, dist: f64 },
    Boolean { kind: BoolKind, a: Handle<Solid>, b: Handle<Solid> },   // BoolKind = Fuse | Cut | Common
    // 강체 변환: `solid` 를 그 상(像)으로 대체한다. `Isometry` 가 정확한 정의(로그의 진실)이고
    // 기하는 실현된 캐시다.
    Transform { solid: Handle<Solid>, isometry: Isometry },
    // 좌표 평면 `axis = offset` 에 대한 반사(대체). 길이 보존·손지기 반전 — 음수 스케일이 아니다.
    Mirror { solid: Handle<Solid>, axis: Axis, offset: Rat },
    // 제자리 복제, 원본은 live 로 남는다 — `live_solids` 에 빼지 않고 더하기만 하는 유일한 연산.
    Copy { solid: Handle<Solid> },
}

// nacre-ops 의 자유 함수다 — Model 의 메서드가 아니다. `impl Model` 은 inherent impl 이라
// Model 이 정의된 nacre-topo 에만 놓을 수 있는데(coherence), replay 는 Operation 을 다루므로
// topo 보다 위 레이어에 살아야 한다.
///
/// 로그를 처음부터 재생. 보장: 동일 로그 → 동일 모델(인덱스까지 재현).
/// **6변종 전부**에서 성립한다 — replay 가 연산 하나마다 로그의 핸들을
/// 인덱스로 재고정하기 때문이며(「인덱스 어휘」), 잠금은 `tests/invariants/replay.rs` 다.
/// `apply` 는 재고정하지 않는다 — 그 모델은 호출자의 것이고, 거기서 조용히 재고정하면 진짜
/// 남의 핸들을 세탁해 교차-모델 가드를 무력화한다.
/// 로그 중간 파라미터를 수정한 재생은 v1 에서 미지원 — Operation 이 원시 Handle 을
/// 참조하므로 상류 수정이 하류 Handle 번호를 밀어낸다(topological naming 문제).
pub fn replay(ops: &[Operation]) -> Result<Model, OpError>;
```

`replay` 는 진실 `Model` 만 돌려준다. 테셀레이션은 `nacre-tess` 가 모델에서 따로 유도하는 캐시이고(`TessConfig`), `nacre-ops` 는 제품 경로에서 `nacre-tess` 에 의존하지 않는다 — `Model` 이 tess 를 필드로 담지 않는 것과 같은 층 규칙이다.

**연산은 평면을 «이름 부른다».** `Extrude` 가 `SketchFrame`(평면 핸들 + 배치 + 쓰임새가 정한 `flip`)을 받으므로, 평면을 **값**으로 싣는 변종은 `DatumPlane` 의 `Stated` 하나뿐이다. 그래서 원칙 2 가 연산 어휘 전체에서 성립한다 — 평면은 아레나에 있고 연산은 그것을 가리킨다. 밑캡이 같은 평면의 두 번째 진술이 아니라 **공유된 핸들**이 되는 이유다.

| 무엇 | 어디서 |
|---|---|
| 세계 평면 위 스케치 | `SketchFrame::world(&m, Axis)` — 씨앗을 이름 부른다 |
| 기존 면 위 스케치 | `face_sketch_frame(&m, face)` |
| 그 밖의 평면 | `Operation::DatumPlane` 로 **먼저 진술**하고 돌려받은 프레임을 쓴다 |

**ŵ 는 프레임의 것이고, «어느 쪽»은 `dist` 의 부호다.** `flip` 은 `frame_toward` **한 곳에서만** 정하므로 호출자가 뒤집힌 프레임을 진술할 수 없다 — 면의 프레임은 그 면의 바깥을 향한다. 프레임을 뒤집으면 축까지 다시 유도되어(`flip` 은 계수의 부호를 먼저 바꾸고 `+u` 를 그 법선에서 짓는다) 같은 스케치가 거울상으로 떨어지므로, «면에 그려 안쪽으로»는 같은 프레임에서 ŵ 반대로 쓰는 **음의 돌출**이다. 밑캡은 어느 쪽이든 프레임의 평면 핸들이다. 원통 축은 쓸기 부호와 무관하게 프레임의 법선이라, «datum −n 에 +d» 와 «datum +n 에 −d» 는 같은 원통을 두 진술로 적는다 — 원통의 interning 열쇠가 진술 그대로인 것(todo 「곡면 둘 이상이 만나는 점과 seam 담체」의 원통 정체)과 같은 부류다.

- **방향은 `flip` 이 들고, 프레임을 만든 쪽이 진실에서 정한다.** 평면의 정준 이름에는 방향이 없고 평면은 interning 되므로(같은 평면을 `+n`/`−n` 으로 진술하면 **한 핸들**), 프레임이 방향을 담을 수 있는 자리는 `flip` 뿐이다. 답은 두 진실의 곱이다: 평면 자신의 좌표계에서 ŵ 가 평면의 향 쪽인가(`frame_normal_sense` — 이름 있는 평면은 이름 법선, 없는 평면은 판정 점의 순서로 프레임을 지으므로 `plane_name_sense`·`sense` 가 답한다) × 모션의 손방향(`motion_parity`) × 쓰임새가 향할 쪽 — datum 은 **호출자가 진술한 법선**(넣은 진술은 그 반대를 향하고, 문이 `flipped` 로 반대로 앉은 핸들을 알린다), 면은 **바깥**(`orientation`). 실현한 ŵ 와 f64 방향의 내적은 같은 물음을 두 반올림에 묻는다(「가지 말 것」 «향 관계를 캐시 법선으로 읽기»).
- **`world_zx` 는 유도로 만들 수 없다**: arbitrary-axis 규약이 ZX 에 `+u = −x̂` 를 주는데 규약은 `+u = +ẑ` 다(`ŵ` 는 둘 다 `+ŷ`). `SketchFrame::world` 가 그 예외를 **한 곳에** 가둔다 — `canonical(씨앗 ZX)` 로 바꾸면 그 스케치들이 90° 돈다(음성 대조로 잠금).
- **`dist` 는 부호 있는 두께다**(0 은 `ZeroDistance`). ŵ 와 스케치의 축은 프레임이, 몸이 평면의 어느 쪽에 서는지는 부호가 말한다 — `DatumDef::Offset` 의 부호 있는 변위와 같은 자리다. 반대편을 향하는 **스케치**(축까지 뒤집힌 것)가 필요하면 그 방향으로 평면을 진술한다.
- **프레임의 기저는 유리수로 묻는다, 실현해서 되묻지 않는다.** `RatFrame::of_plane_frame` 이 `û = u_raw·inv_sqrt_exact(uu)` 로 답한다 — 실현한 축을 `Rat::from_decimal` 로 다시 들어올리면 정규화가 필요한 축(`(0.6,0.8,0)` → 원시 `(3,4,0)`, `uu=25`)이 `0.6000000000000001` 로 돌아와 직교정규가 깨지고, 그 평면이 **조용히 프레임-노드 도로로 옮겨간다**. `plane_frame_named` 가 `v̂` 에 대해 적어 둔 규칙과 같은 것이다.

**datum 평면의 규칙.** 평면을 만드는 연산은 `DatumPlane` **하나**이고, 불리언은 평면을 만들지 않는다(깊이 불변식, `a_boolean_mints_no_surface`).

**정점을 이름 부르는 datum.** `ThroughVertices([Handle<Vertex>; 3])` 는 값 어휘로는 말할 수 없는 하나다 — 발견 정점의 좌표는 반올림이라, 그 좌표로 평면을 지으면 **다른 평면**이 나온다(기울어진 인구 220/220, 축정렬 음성 대조 552/0; `tests/instruments/point_width.rs`). 진실 쪽은 `PlanePoints::Through` 로 **핸들을 든다**.

- **어휘가 여기서만 자란다** — «새 공개 생성자 없음»은 *값으로 말할 수 있는 것*에 걸리는 원칙이고, 핸들은 순수 값 타입인 `SketchPlane` 에 들어갈 수 없다.
- **정렬은 키에만, 방향은 호출자의 정점 순서.** `dist` 가 양수 전용이라 순서가 방향의 유일한 입구이고, 두 개를 바꾸면 «같은 핸들 + 반대 프레임»이다(`frame_toward` 가 `Stated` 와 같은 기계로 그 쪽을 향하게 한다).
- **이동은 핸들을 그대로 두고 노드를 기록한다.** `transform_solid` 는 정점을 **복제**하므로 가리키던 datum 은 새 복사본을 안 따라간다; 정점이 base 를 정하고 모션이 옮기며 **더해질 뿐 곱해지지 않는다**. 대가로 그런 datum 위의 솔리드는 **정확한 강체 이동에도 노드를 얻는다**.
- **거절은 원인별**(`VerticesInMixedFrames`·`CollinearVertices`·`DuplicateVertex`·`VertexNotThreePlane`)이고, 전부 **무엇이든 push 하기 전에** 일어난다. «담체 셋이 안 만난다»는 거절이 아니라 **단언**이다 — 그 정점이 존재한다는 것이 곧 만났다는 뜻이므로 불변식 위반이지 사용자 오류가 아니다.
- 점과 모션은 `Offset` 처럼 **한 결정에서 함께** 나와야 한다 — «점이 어느 프레임에 적혔나»가 곧 모션이기 때문이다.

datum 일반:

- **캐시 법선은 `−(진술된 법선)`** — 씨앗이 `−축`, extrude 밑캡이 `−plane.normal()` 인 그 규약. 근거 둘: `+축` 씨앗에서는 781 캡의 저장 법선이 뒤집히고, `WorkingPlane::frame_sign` 이 «저장 법선 vs 루트 면의 바깥 법선»인데 밑캡의 바깥이 `−N` 이다.
- **연산은 호출자의 프레임을 함께 돌려준다**(`OpOutput::DatumPlane { plane, frame }`) — 편의가 아니라 **정확성**이다. 호출자가 `SketchFrame::named` 로 다시 말하면 f64 를 거쳐 `Rat::from_decimal` 을 두 번 타고 계산값은 다른 유리수로 떨어질 수 있다. `Stated` 의 배치는 **무조건 `Named`**: ZX 평면의 정준 `+u` 는 `−x̂` 인데 규약은 `+ẑ` 라, 유도로 흉내내면 그 스케치들이 조용히 돈다.
- **`Offset` 은 push 전에 정규화한다 — 셋 다 «한 평면에 두 핸들»을 막는 조항이다.** 신원은 `(평면, 부호 있는 거리)` 의 순수 함수여야 한다(배치의 원점·`+u` 는 평행 평면을 옮기지 않는다): ① `flip` 을 부호로 접고 **평면의 정준 프레임**에서 짓는다 ⇒ 한 평면의 다른 프레임을 든 두 호출자가 같은 핸들을 받는다. ② 정준 기저가 정확히 리프트되면 **세계로 말한다**(노드 생략 규칙을 실현된 기저에 적용) — 프레임 노드 아래 두면 키가 `(이름, Some(노드))` 라 상자 캡의 같은 평면과 **interning 되지 않는다**. ③ **`dist == 0` 은 이름 붙은 거절**(`ZeroOffset`): 기울어진 프레임의 밑 평면은 세계에서 유리수라 ②가 안 걸리므로, 0 만이 남는 충돌이다. 0 이 아닌 오프셋은 계수가 `c ∓ d·√(n·n)` 이라 무리수여서 경쟁 진술 자체가 없다. 세계 되당김이 넘치면 **조용히 노드 도로로 가지 않고** 거절한다 — 그 도로가 곧 중복이 생기는 자리다.
- **datum 은 솔리드 이동을 따라가지 않는다.** append-only 라 datum 이 가리키는 평면은 남고, 나중 `Transform` 은 새 평면을 만든다(의도된 동작 — datum 은 독립적인 기준이다).
- **쓰이지 않는 datum 은 어디에서도 새지 않는다**: `validate` 는 live 도달 집합만 세고 surface 는 경계 검사만 받으며(씨앗 셋이 이미 영구 orphan), `nacre-step`·`nacre-tess` 는 surface store 를 돌지 않는다(둘 다 면 경유).

### 프로파일 — 다중 루프와 계약

**다중 루프 프로파일 — 구멍 있는 스케치를 정확히 세운다.** `Profile2d { outer: Ring2d, holes: Vec<Ring2d> }` 가 외곽 링 하나와 구멍 링 N개(제한 없음)를 담고(`polygon`/`with_holes` 생성자, 필드는 비공개; `inners` 는 필드가 아니라 `with_holes(outer, inners)` 의 **인자** 이름이다), `build_prism` 이 구멍마다 벽을 세우고 두 캡에 내부 루프를 단다. **불리언으로 흉내낼 때와의 차이가 요점이다** — 잘라 만든 도넛은 모든 코너가 `Discovered`(측정된 tol)이지만, 프로파일을 쓸어 만든 도넛은 전부 `Constructed`(tol 없음)다. 원칙 4 가 요구하는 바로 그 차이다. 그래서 **설탕으로 흉내내면 안 된다** — "외곽 extrude → 구멍 프리즘 Cut"은 전부 `Constructed` 였을 모델을 불리언·`Discovered` 경로로 내리므로 원칙 4(tolerance 는 발견된 교차에만)를 스스로 어긴다.

**감김은 한 곳에서만 정한다.** `build_prism` 이 sweep 을 아는 유일한 자리이므로 외곽을 sweep 기준 CCW 로 정규화하고 구멍을 그 반대로 맞춘다. 타입은 링을 그대로 담는다 — 두 곳에서 정하면 음의 돌출처럼 sweep 이 외곽을 뒤집는 경로에서 둘이 갈린다. 이 설계는 그 버그 부류를 구조적으로 없앤다: "정규화된 외곽의 반대"와 "sweep 기준 CW"가 같은 규칙이 되기 때문이다.

검산: 도넛 각기둥은 V16 − E24 + F10 − L_i 2 = 0 = 2(S−G), genus 1 — **오일러의 `L_i` 항이 실제로 필요한 형상**이다.

**프로파일 계약은 커널이 검사한다 — `Profile2d::check`.** 생성자(`polygon`/`with_holes`)는 **진실을 잡는다**: 각 좌표를 `Rat::from_decimal` 로 리프트해 `Ring2d` 의 유리수 정점으로 보관하고(십진 창 밖 = `ProfileOutsideDecimalWindow` 구성 시점 에러), **공선 중간점을 소멸**시킨다(엄격 내부만 — 무손실 정규화; 중복점·스파이크는 생존해 제 이름의 에러로 보고된다. 이로써 모든 프리즘 코너가 3-평면 정의를 보유한다). *계약*(단순성·서로소·포함)은 생성자가 아니라 **프로파일을 소비하는 모든 연산이 먼저 `check()` 를 통과시킨다**(진입점은 `extrude_on_frame` 하나다). 계약은 넷: 모든 링이 **단순 다각형**, 링끼리 서로소, 구멍은 외곽 안, 구멍 속 구멍 없음(그건 섬이고 `from_rings` 가 갈라낸다). 판정은 전부 **저장된 유리수 진실 위**에서다 — f64 이진값의 정확 부호와 십진 진실의 부호는 퇴화 근처에서 실제로 갈리고(0.1·0.2·0.3 공선이 이진에선 굽음), 작성자가 쓴 수가 이긴다.

**왜 "호출자의 약속"으로 둘 수 없는가.** 계약을 어긴 프로파일은 검사가 없으면 전부 `extrude` 성공 · `validate` 위반 0건 · STEP 내보내기 성공이다. 나비넥타이는 부피가 `NaN` 이지만, **구멍이 외곽 밖이면 12.0(정답 16), 구멍 속 구멍이면 20.0(짝-홀 정답 52)** — 눈치챌 단서가 없는 그럴듯한 숫자다. 조용히 틀린 답은 약속으로 둘 수 있는 종류가 아니다.

부수적으로, 이 검사가 `oriented_ring` 의 전제를 세운다: 감김을 부호 있는 면적 벡터로 정하는데 면적이 0 이면 그 판정이 우연에 맡겨진다(대칭 나비넥타이는 두 엽이 정확히 상쇄해 0 이다). 단순 다각형은 면적이 0 일 수 없다.

**단순성은 *작성된 입력*의 계약이지 커널 데이터의 불변식이 아니다.** 그리는 링의 안쪽은 단순성 위에서 짝-홀로 *정의되지만*, 불리언이 만드는 윤곽의 안쪽은 arrangement 가 이미 정했고 그쪽은 **비단순(figure-8 조임)이어도 정당하다**(`loop_winding`). 그래서 `validate` 에 면 단순성 검사를 넣으면 안 된다 — 미루는 게 아니라 정당한 결과를 거절하게 된다.

### 스케치 어휘 — 링 중첩·고리·호

**링 중첩 판정 — `from_rings`.** 사용자는 닫힌 경로만 그리고 무엇이 구멍인지 말하지 않는다. `ops::sketch::from_rings(rings) -> Result<Vec<Profile2d>, SketchError>` 가 포함 **깊이**로 정해 **덩어리(섬)별 프로파일 목록**을 돌려준다: **짝수 = 재료, 홀수 = 구멍**(그래서 구멍 속 링은 다시 재료 = 섬이고 자기 몫의 profile 이 된다), 구멍은 **가장 깊은 포함자**(직계)에 붙는다. **채우기 규칙 인자는 없다** — 짝수-홀수가 유일한 규칙이고, 인자를 두면 그 규칙이 둘이 된다.

**경계:** 커널 `Extrude` 1회 = **연결된 덩어리 1개**(외곽 + 그 구멍들). 섬마다 호출해 결과를 묶는 것은 편의 레이어(overview.md 판별 기준)다. 그래서 `Extrude` 의 다중 바디 출력은 필요 없다. 구멍 있는 스케치(도넛)와 섬이 여러 개인 스케치는 코드-CAD 가 요구하는 것이다.

**분류 술어는 `geom::intersect` 에 있다**(격리 규칙): `point_in_ring_2d`(exact `orient2d` 교차 패리티, `RingSide::{Inside, Outside, OnBoundary}`)와 `rings_cross`(적절 교차 + 접촉). 각 워커에 **Rat 쌍둥이**(`_rat` 접미 + `drop_collinear_midpoints`)가 같은 모듈에 나란히 산다 — 부호 원시는 `nacre_exact::orient2d_rat`(narrow 우선 → BigInt 전역, 분모 청소로 부호 보존)이고, geom 은 scalar 에 의존한다(최하단 토대라 순환 없음). f64 판은 tess·f64 폴백이 쓴다. ops 에는 **정책(깊이 패리티)과 조립**만 남는다. 링이 서로 닿거나 교차하면 `SketchError::RingsMeet` 으로 거절한다 — 접촉도 실패다(엄밀한 안쪽이 없다).

**자기교차는 중첩 분류보다 먼저 본다** — 메시지 품질이 아니라 전제조건이다: `point_in_ring_2d` 는 짝-홀 패리티로 답하고, 그건 단순한 링에서만 "안쪽"을 뜻한다. 술어는 `geom::intersect::ring_self_intersection`(인접하지 않은 변은 접촉만으로 실격, 인접한 변은 공선-겹침일 때만 = 되짚는 스파이크, 인접은 **순환**으로 판정). 스케치 층은 `RingSelfIntersects { ring, at }` 로 **점**(걸린 두 변의 현 중점)을 돌려준다 — 점은 편집기가 표식을 놓을 수 있는 것이고 이 enum 의 다른 거절도 점을 든다.

**커널의 문은 고리다**: `Ring2d::new(vertices, edges)`(정점 + 각 정점을 떠나는 변, 검증: 짝 맞음·영길이·호의 r²·다음 정점이 원 위) · `Ring2d::circle` · `Ring2d::polygon_decimal`, 그리고 `from_paths(Vec<Ring2d>)` 가 중첩을 분류한다. 커널은 «어떻게 그렸나»를 모른다 — 펜·제약 해석기·외부 데이터가 전부 같은 문으로 온다.

**호와 원.** 조각 타입은 **하나** — `Edge2d { Line, Arc { center, r2, ccw } }`(`nacre-geom::mixed`) — 끝점은 고리의 정점이지 조각의 필드가 아니다. **반지름은 저장하지 않는다**: 진실은 `r2 = |start − center|²` 이고 그것은 언제나 유리수다(√2 반지름의 호도 담는다). 정확 술어는 r 을 제곱해서만 쓰고, r 자체는 f64 캐시에서 √ 한다. 회전각은 커널이 모르는 말이다: 끝점 표현이 각도 표현보다 **넓고**(무리수 도의 호도 끝점이 유리수면 정확), 90° 배수 회전은 단계 문이 `Rat` 으로 끝점을 계산해 준다(`arc_turns`·`arc_turns_rat` — (변, 끝점)을 돌려준다; 임의 끝점은 `arc_to_rat`). 링은 **정점 + 변 종류**(`Ring2d { vertices, edges: Vec<Edge2d> }`)라 끝점 중복이 타입에 없고, 정규형이 공선 중간점 제거 옆에서 **같은 원·같은 방향의 이웃 호를 병합**한다(필렛 둘이 만나면 반원, 호 둘로 그린 원은 원 하나; 정점 하나에 호 하나면 온전한 원이고 그 정점이 seam 이다). 분류 술어는 `geom::mixed`(선분–원·호–호·점-in-링을 ℚ(√c) 의 `QuadVal` 부호로, 광선 규칙은 배열의 반열림 철자를 옮겨 적음; 다각형 링에서 f64/Rat 술어와 답이 같음을 proptest 로 잠금). 프리즘 빌더는 호마다 원통 조각 벽을 세우고(축 = 프레임 법선 `w`, `ref_dir` = 프레임 `x̂`, 온전한 원만 시작점 방향), 직선–호 꼭짓점은 **`Vertex::Pierce`**(두 평면의 저장된 정준 이름을 핸들 오름차순으로 `plane_plane_cylinder` 에 — 접하는 벽은 `Double`), 감김은 `2·넓이 = a + b·π` 의 정확 부호(`winding_sign_quarter_arcs`, scalar 의 BigInt). 옆면 방향 `Forward ⟺ 호의 ccw == 스윕이 +w`. 이름 붙여 거절: 90° 배수 아닌 호(`ArcSweepNotQuarterTurn` — 어휘는 그런 호를 정확히 진술하지만 빌더의 감김이 `π/2` 배수를 요구한다), 다른 원의 호–호 접합(`ArcsMeetAtVertex` — 두 원통과 한 평면 위의 점은 정점 정의가 없다).

**끝점은 정확히 일치해야 한다.** 가까우면 붙여주는 스냅은 없다 — tolerance 는 커널이 *발견한* 교차의 것이지 호출자가 *구성한* 것의 몫이 아니다(원칙 4). 고리 문에서 그것은 «다음 정점이 호의 원 위에 있는가»(`ArcEndOffCircle`)·«직선의 두 끝이 같은가»(`ZeroLengthEdge`)로 나타난다.

**파라메트릭 편집의 진화 경로 (계획: v2 이후, 기록만).** 연산이 원시 Handle 대신 계보 참조 `OpRef { op: usize, output_slot: usize }`("연산 N 이 만든 k번째 면")를 담으면, 상류 수정 후에도 참조가 의미로 해석(resolve)되어 편집-재생이 가능해진다. v1 에서는 구현하지 않되, 로그 직렬화 포맷을 설계할 때 이 확장이 포맷 파괴 없이 들어갈 자리를 남긴다. 인덱스 어휘는 이를 **밀어내지 않는다** — 인덱스는 v1 의 참조 기계이고 `OpRef` 는 v2 의 것이며, replay 안의 재고정 지점이 곧 `OpRef` 해석이 들어올 자리다(교체이지 경쟁이 아니다). 테셀레이션 tolerance 도 같은 맥락에서 연산별 override(`Operation` 항목의 선택 필드)로 확장될 수 있다.

OCCT 는 제품 경로에 등장하지 않는다 — 역할은 nacre-oracle 의 채점자뿐이다. 사다리의 각 단은 자기 커버리지 안에서 완전해야 하며, 밖은 조용히 틀리는 대신 에러로 거절한다.

### 거절 분류

**거절의 *이유*는 값에 실려 나간다.** `BoolError::Rejected { reason: RejectReason, at: Option<RejectWhere> }`. **소비자는 `RejectReason::class()` 로 분기한다** — `RejectClass` 는 셋이다: `NotSupported`(커널의 현재 커버리지 밖 — 입력 자체는 판정 안 됨; 「아직」이 아니다. 의도된 거절은 미완성 기능이 아니고, 일부는 영영 더 정확한 거절만 얻는다) / `Impossible`(어떤 마일스톤에서도 유효 솔리드가 없음 — 사실의 진술이지 조언이 아니다) / `SuspectedDefect`(엔진 불변식이 깨짐 — 버그 리포트 감이고 «당신 탓»이라 말하지 않는다). 한 번도 발화가 관측되지 않은 변종의 분류는 잠정이다. 변종 이름은 엔진 어휘라 로그·리포트용 안정 식별자로만 쓴다. 사람이 읽는 문장·현지화는 앱 몫이다.

성장은 사유 어휘에서만 일어나므로 `RejectReason`(과 그 세부인 `DeclineKind`)만 `#[non_exhaustive]` 이고 `BoolError`·`RejectClass` 는 exhaustive 다(소비자가 모든 클래스에 대해 결정하도록 강제한다 — 사유는 자라고 클래스는 안 자란다). 트레이스가 불완전한 경우는 `TraceDeclined { kind: DeclineKind, face: Option<Handle<Face>> }` 가 **무엇을 못 했는지와 어느 피연산자 면에서인지**를 함께 싣는다(한 클래스가 여러 면을 기권할 수 있고, 이것은 첫째를 이름 댄다).

**«어디»는 사유 옆에 실려 나간다.** `at: Option<RejectWhere>`(`Point(Point3)` | `Segment([Point3; 2])`, 월드 f64)는 가드가 발화한 순간 보고 있던 **증인**이다(여럿이면 엔티티 자체 순서의 최솟값 — HashMap 순회로 뽑으면 실행마다 다른 좌표가 나온다). 위치를 `RejectReason` 변종 안에 넣지 않는 이유: 사유는 범주 어휘(census 키·테스트가 이름 대는 것, `Copy+Eq`·const-구성)이고 좌표는 측정값이다 — census 의 "측정값은 키 밖" 규칙의 오류-값 판. 좌표는 진단용 실현(캐시급)이라 `RejectWhere` 는 `Eq` 가 없고, 비교는 근사로만 한다(거절끼리 비교할 때는 `reason` 을 꺼내 비교한다). 표면화 사유 중 `self_touching_result`(위반 모서리)·`non_manifold_result_edge`(위반 모서리)·`non_manifold_vertex`(핀치 정점)·`tangent_line_in_another_plane`(접선 위 한 점)이 싣고, `no_clear_ray`(분산적)·`precision_budget`(전-모델)·`cylinder_face_undecided`는 싣지 않는다(면이 증인인 사유가 필요해지면 `RejectWhere` 에 면 변종을 여는 것이 길이다).

**어떤 사유가 실제로 발화하는지는 상설 census 가 답한다 — `nacre-ops::reject_census`.** `reject()`/`reject_at()` 쌍이 크레이트의 모든 `Rejected` 를 짓는 유일한 깔때기라, 거기 `#[track_caller]` 하나씩이 모든 호출 지점을 한꺼번에 계측한다. 테스트 빌드의 계측이다(`test-util` — overview 「작업 스타일」의 계측의 문): 제품 빌드에는 표도, 기록하는 호출도, `reject` 의 `#[track_caller]` 도 없다.

**두 열은 서로 다른 인구다.** 「울린 것(raise)」과 「밖으로 나간 것(surfaced)」을 따로 적는다. 호출자가 대안 경로로 재시도하며 삼키는 거절이 있으면 한 열 집계는 사용자가 보는 것에 대해 거짓을 말한다(`no_clear_ray` 가 울림의 82% 를 차지하면서 한 번도 표면화되지 않는 구성이 실제로 있었다). 그래서 `point_in_component` 의 «이 노드로는 못 정함»은 `Ok(None)` **기권**이지 삼켜질 오류가 아니다 — 두 열의 갭은 「남은 삼킴의 지표」이고 0 에 가까운 것이 목표 상태다(갭이 벌어지면 새 삼킴의 발견이다). 다만 census 는 **어느 사이트가 표면화했는지는 답하지 않는다** — 거절은 rayon 워커에서 울리고 메인 스레드에서 반환되므로 둘을 잇는 값싼 길이 없고, 「이 가드는 삼켜지는가」는 호출자를 읽으면 되는 **정적** 질문이라 필요도 없다.

읽는 길 둘: **게이트**는 `tests/reject_census.rs` 의 얼린 코퍼스 7종이 `(사유, detail, 파일)` 집합을 박는다(줄 번호는 찍기만 — 거기 걸면 무관한 편집마다 빨개져 아무도 안 읽는다. 횟수도 안 건다: debug 는 트레이스를 두 번 돌리고 parallel 은 실패 입력의 나머지 클래스까지 평가한다). **전수**는 `--features reject-trace` 로 스위트 전체를 훑는다 — 테스트 바이너리 30여 개가 각각 별개 프로세스라 메모리 표가 프로세스와 함께 죽기 때문이고, 그래서 이 기능이 크레이트의 `deny(print_stderr)` 에 대한 유일한 예외다(기본 꺼짐, 제품 빌드엔 출력 코드가 아예 안 들어간다). 게이트의 인구는 그 7개 모양뿐이니 **원통 가드는 여기 안 걸린다** — 그건 전수 훑기가 손으로 메운다.

census 키가 `detail` 을 드는 이유: `TraceDeclined` 의 모든 kind 가 **한 줄**에서 올라와, 이름과 위치만으로 키를 잡으면 한 행으로 뭉개진다 — census 가 잡으려는 결함을 census 가 재현하는 꼴이다. 규칙은 **키가 그 사유 자신의 어휘만큼 잘아야 한다**는 것.

### 구성 산술은 **사용자가 쓴 십진수**로 한다

프리즘의 좌표는 곱 셋과 합 둘로 나온다 — `origin + x·u + y·v`, 그다음 `+ normal·dist`. f64 에서 이건 결합법칙을 안 지키므로 **정확히 일치해야 할 두 치수가 어긋난다**: `7.7` 로 한 번에 올린 블록과 `1.1` 다음 `6.6` 으로 올린 블록의 윗면이 1 ULP 떨어진다.

**쪼개지지도, 틀리지도 않는다 — 더러워진다.** 부피는 정확히 맞고 fuse 는 몸통 하나를 준다. 그 몸통이 넓이 `8.9e-16` 짜리 면과 형상에 없는 면 둘을 달고 있고, 이후 모든 불리언·STEP·메시가 그걸 끌고 다닌다. **거절도 오답도 아니라서 아무도 신고하지 않는다**.

**"API 를 `Rat` 으로 바꾸면 된다"가 답이 아니다.** f64 `1.1` 을 `try_from_f64` 로 *올리면* `11/10` 이 아니라 그 f64 가 담은 이진값이다. 정확한 산술을 해도 이진 드리프트가 그대로 재현될 뿐 일치가 **복구되지 않는다**. 복구하는 것은 사용자가 *쓴 십진수*다 — `11/10 + 66/10 = 77/10`, 그리고 그 실현이 정확히 `7.7`.

**결정:**

- **공개 API 는 f64 그대로.** 호출부·TS/JSON 경계가 안 바뀐다. 커널이 받자마자 각 치수를 `Rat::from_decimal` 로 **왕복하는 최단 십진수**로 잡는다. 이 규칙은 프로파일 좌표, extrude 거리(`DistOutsideDecimalWindow`), 스케치 평면의 축(`from_axes` — 정의는 점 셋 `[o, o+x, o+y]`, `PlaneDef` 는 점 셋 단일 필드로 origin·ref_dir·극성이 구조에서 유도)에 걸린다. 단 **커널이 실현한(반올림된) 기저는 리프트하지 않는다**(내부 `realized_plane` — 이미 정확한 정의를 가진 평면에 두 번째 진실을 만드는 것이 두-정확-기술 결함이다). 의도에 대한 *추측*이 아니라 **결정적 정규형**이다(왕복이 보장되므로 서로 다른 f64 는 서로 다른 유리수, 같은 f64 는 항상 같은 유리수). 십진수를 고르는 이유는 하나 — **사람은 십진수를 타이핑한다**.
- **전제조건: `Rat::to_f64` 가 최근접 반올림이어야 한다.** 분자·분모를 각각 반올림한 뒤 나누면 유효숫자 17자리에서 1 ULP 가 어긋나고, 그러면 십진수 복원이 **무손실이 아니다** — 고치려던 것은 어긋난 누적뿐인데 멀쩡한 입력까지 움직인다. 정수 장제법으로 54비트 몫을 뽑아 **끝에서 한 번만** 반올림한다(`Ratio<i128>` 이라 두 항이 2¹²⁷ 미만인 게 모든 시프트를 닫는다).
- **프레임이 유리수일 때만 이 경로를 탄다.** 축이 `{0, ±1}` 이면(world·축정렬 면) 정확하다. 회전된 면은 축이 무리수라 못 든다 — 성분을 십진수로 잡으면 **단위벡터가 아니게 되어** `normal·dist` 가 틀린 *길이*를 준다(고치려던 것보다 나쁜 오차). 그래서 진입점이 `x·x = y·y = 1`, `x·y = 0` 을 **정확히** 검사한다. 평면(또는 그것을 나르는 프레임 사슬)에 정확한 꼴이 없으면 — 축이 십진 창 밖이거나 퇴화했거나 배치 산술이 `i128` 을 넘치면 — 조용히 f64 로 짓지 않고 **이름 붙여 거절**한다(`PlaneWithoutExactForm`): f64 로 지은 프리즘은 정확한 점을 기록하지 못하고 모션을 견디지 못한다.
- **새 저장 구조가 필요 없다.** 같은 유리수는 항상 같은 f64 로 실현되므로 두 경로의 좌표가 **비트 동일**해지고 평면 계수도 같아진다. provenance 를 새로 다는 일이 없다.

**못 고치는 것 (비목표, 정직하게):**

- **회전된 스케치 프레임.** 축이 무리수라 유리수로 들 수 없다. **f64 삼각함수를 거친 사분 회전도 마찬가지다** — `(90°).to_radians().sin_cos()` 가 `cos = 6.1e-17` 을 주므로 정확성은 이 층이 보기 전에 이미 사라졌다(커널의 정확한 사분 회전은 `Angle`/Niven 에서 온다). 그 경로의 건전성은 CIP 담당이다.
- **스크립트가 계산해 넘긴 값.** `1.1 * 7` 을 f64 로 계산해 넘기면 커널은 `7.700000000000001` 을 받고, 그 최단 십진수는 `7.7` 이 아니다. **올바른 동작이다** — 스크립트가 실제로 다른 수를 계산했다. 고치는 것은 **커널 자신의 누적**이다.
- **4-평면 동시성.** ±1 ULP 로 사라지는 정확한 0 이다. 유리수 입력에서도 생긴다(0.1·0.3 캡과 벽 `y = 3x`) — 구성 산술이 고칠 것이 아니라 판정이 진실에 물을 것이고, 판정은 평면의 이름으로 묻는다(`WorkingPlane`).

### 불리언 엔진 — 평면당 셀 복합체

불리언 경로는 **하나**다: `nacre-ops::arrangement` 의 평면 클래스 배열 엔진(`crate::boolean` 의 공개 `boolean` 이 여기에 위임한다). 글로벌 detector 도, regime 별 이중 경로도 없다. 기하가 평면 위에 놓이는 것은 퇴화가 아니라 **입력**이다.

**단계.**

1. **평면 클래스.** 두 솔리드의 면이 놓인 평면을 클래스로 묶는다. 이후의 일은 면당이 아니라 **평면 클래스당** 독립이다.
2. **추적(trace).** 두 솔리드의 각 면을 모든 평면 클래스 위에 **선분**(닫힌 루프가 아니다)으로 추적한다. 면의 경우는 둘이다. **Seated** — 자기 클래스가 자르는 클래스인 면 — 은 통째로 그 평면에 놓이므로 경계 전체가 추적이다. **Transversal / mixed** — 평면을 가로지르는 면(모서리가 평면에 놓일 수도 있다) — 은 `L = W ∩ fp` 위의 2D 다각형-대-직선 클리핑이고, 3값(−/0/+) 스캔이 선 위의 모서리를 거절하지 않고 접는다. 추적을 못 한 면은 **면 단위 기록**(`declined`)이지 솔리드 전체의 중단이 아니다 — 비어 있지 않은 `declined` 는 "이 추적은 불완전하니 이것으로 결론 내지 말라"는 뜻이고, 밖으로는 `TraceDeclined { kind, face }` 로 나간다.
3. **교차 분할.** 추적된 선분을 서로의 교차점에서 자른다(`split_at_crossings`).
4. **셀 추출·중첩.** 분할된 선분에서 셀을 뽑고, 구멍 셀을 host 셀에 중첩시킨다(`+1` 감김 셀만 host 다 — `−1` 셀은 어떤 host 의 구멍이거나 빈 뿌리다).
5. **라벨.** 셀마다 피연산자별·쪽별 재료 라벨 `Label = [A_above, A_below, B_above, B_below]` 을 단다. 모서리를 건널 때 그 모서리의 마스크가 라벨을 뒤집고, 전파가 끝나면 **모든** 모서리(트리·비트리)에서 뒤집힘 관계를 검증한다(`UnreachedCell`·`LabelConflict`). **"above" 는 클래스의 *저장된* 평면 법선 쪽**이다 — 클래스의 정준 유리수 이름(절반의 클래스에서 반대를 향한다)도, 루트 면의 바깥 법선도 아니다(`frame_sign` 이 둘을 잇고, 방출의 `flip` 이 그것을 되접는다).
6. **생존 규칙·방출.** 아래.
7. **조립·청소.** 엔진이 `assembly` 의 `assemble_fuse_cut` 과 `unify_coplanar_faces` 를 불러 껍질을 짓고 청소한다. 입력은 `NodeId` seam 삼중으로 이름 붙은 면당 출력 `draft::LocalFace` 다. `assemble_fuse_cut` 이 vertex handle 을 면의 first-appearance 로 배정한다.

**생존 규칙.** 규칙은 한 챔버(평면 한쪽의 공간)의 `(inA, inB)` 에 대한 술어 하나다:

| `kind` | `keep(inA, inB)` |
|---|---|
| `Fuse` | `inA ∨ inB` |
| `Cut` | `inA ∧ ¬inB` |
| `Common` | `inA ∧ inB` |

`+1` 셀은 **그 두 챔버(W 의 위·아래)가 `keep` 아래에서 갈릴 때** 결과의 면이다(한쪽은 재료, 한쪽은 빈 공간). 셀의 `−1` 구멍은 내부 링으로 따라간다. `flip = keep_above == (orient_sign(wc) > 0)` 하나가 어느 챔버가 재료인지를 나르고, 조립이 외곽과 내부 링을 함께 뒤집어 결과 법선이 남은 재료의 바깥을 향하게 한다.

접촉면(두 피연산자가 같은 평면에 면을 가진 경우)에 대한 Requicha 식 표 — 연산 × 상대 법선 — 는 이 한 규칙의 귀결이다: Fuse/opp=`inP⊕inQ`, Fuse/same=`inP∨inQ`, Cut/opp=`inP`, Cut/same=`inP∧¬inQ`, Common/opp=∅, Common/same=`inP∧inQ`(`inP`·`inQ` 는 각 피연산자의 면이 그 셀을 덮는가). 볼록성 항이 없어 비볼록·다중루프를 half-space 가정 없이 `point_in_ring` parity 로 정확 처리한다.

**동일평면 병합은 기권하고, 전역 판정이 말한다.** `unify_coplanar_faces` 는 불리언 뒤의 **필수** 청소다: 한 평면 클래스·한 `flip` 의 결과 면들을 **모서리-연결 성분**별로 하나로 합치고 평각 정점을 녹인다. 방법은 «안쪽을 지우고, 바깥을 꿰매지 않는다» — 그룹의 모든 방향 링 모서리에서 마주보는 쌍(`a→b` 와 `b→a`)을 지우고 남은 것을 다시 꿴다. 그래서 한 규칙이 모서리 공유·이웃이 정확히 메운 구멍·여러 이웃이 메운 구멍과 그 연쇄를 다 덮는다. 같은 평면에 있어도 닿지 않는 면은 각자 살아남는다. 병합된 윤곽이 figure-8 이 되거나(다시 꿸 순환 집합이 유일하지 않다) 두 면이 같은 방향 모서리를 주장하면 병합은 **기권**하고 그룹을 그대로 내보낸다. 기권이 안전한 이유는 「병합은 정돈일 뿐」이 아니다 — 병합되지 않은 동일평면 면은 떨어진 벽 위의 코너가 결과에 면이 없는 평면을 이름 대게 만들어 조립이 불리언 전체를 거절한다. 안전한 이유는 더 좁다: 그대로 내보내면 파이프라인이 **이 패스가 돌기 전 자리**에 정확히 서므로, 패스 없이 나왔을 것 말고는 아무것도 나올 수 없다. 그 뒤 **결과 전체에 대한 판정**(`self_touch_reject`, 닫힌 껍질 가드, 모든 정점이 다시 풀리는지 검사)이 자기 문장으로 — 그리고 **증인과 함께** — 그 모양을 이름 댄다. 핀치를 여기서 능력-한계 이름으로 거절하지 않는다: 여기 도달하는 모든 입력에서 그것은 자기접촉 결과이기 때문이다.

**45°/315° 폴드의 거절은 갭이 아니라 «정답»이다.** 그 각도에서 부품의 한 팔이 **자기 자신에** 공면으로 닿아 두께 0 인 솔리드가 나온다 — 존재할 수 없는 물건이고, `SelfTouchingResult`(`Impossible`)가 옳다.

**칼날 접촉**(모서리가 상대 면 평면 위, 두 면이 한쪽으로만 떠남): `edge_mask` 의 스침 결합은 **홀짝**이다(같은 쪽 쌍 = 꼬집힘/노치 = 무소식), all-false 모서리는 골격에서 제외한다(`drop_newsless` — 점 접촉 `touches` 의 1D 판). 모서리-만 접촉 fuse 는 **몸통 둘로 나온다**. 면을 잇는 것은 **다양체 접촉(정확히 두 면이 쓰는 링 모서리)뿐**이라, 선·점으로만 닿는 두 몸통은 서로 다른 성분에 놓이고 각자의 핸들을 받는다. 핸들을 가르는 단위는 성분이 아니라 **출력 솔리드**(재료 성분 + 그 공동들)여야 한다 — 공동이 host 껍질에 닿는 경우 성분별로 가르면 핀치가 셀 수 없어져 두께 0 솔리드가 조용히 통과한다. 자기와 닿는 «한» 몸통은 재료가 접촉을 돌아가므로 성분이 하나로 남아 거절이 그대로 발화한다(잠금: `tests/probes/contact_separates.rs`, `tests/probes/knife_edge.rs`).

**4평면 동시성의 세 갈래.** 한 점에 평면 넷이 모이면 이름이 여러 개 생기고, 배열은 그것을 하나로 접는다(`Aliases`). 발견 경로는 셋이다: ① **입력 면의 꼭짓점**에서의 동시성 — 링을 걸으며 `third_on_l` 이 본다. ② 배열이 만든 점에서 **선을 가로지르는** 두 평면이 겹치는 경우 — `split_at_crossings` 가 「두 handle 이 같은 자리로 정렬된다」로 본다. ③ 배열이 만든 점에서 **그 선을 «담은»** 평면 — 그 평면은 선과 같은 방향 family 라 handle 이 되지 않으므로 ①②가 **구조적으로 못 본다**. 셋째는 따로 찾지 않는다: 벽 둘을 한 선으로 접었다는 기록이 곧 「두 번째 평면이 이 선을 담는다」는 진술이므로, 그 선 위 분할점마다 벽 family 를 되읽어 신고한다. 놓치면 한 점이 두 이름으로 seam 표에 도착해 `SeamAlias`(`SuspectedDefect`)로 거절되는데, **그건 유효한 입력을 두고 커널이 자기를 탓하는 것**이다. 잠금은 `tests/probes/concurrent_line.rs`.

**결과 계약: 표면이 자기와 닿는 솔리드는 만들지 않는다.** 조립이 내놓은 결과 솔리드(바깥 껍질 + 공동)에 대해 *"이 솔리드의 모서리가 같은 솔리드의 어떤 면의 «내부»에 있는가"* 를 묻고, 참이면 `SelfTouchingResult`(`Impossible`)로 거절한다. 이런 몸체는 연결돼 있고 부피가 맞고 **모든 위상 계수를 통과한다** — 접촉선에서 면이 쪼개지지 않아 모든 모서리가 여전히 두 번씩 쓰이므로 `check_result_topology` 도 `validate` 도 보지 못한다. 방향은 Parasolid 기준이다: OCCT 는 이런 몸체를 받지만 따르지 않는다. **대가**는 자연스러운 모델링 동작 하나가 막히는 것이다(쐐기 포켓의 끝이 마침 반대 벽에 닿는 경우). 다른 커널이 *"zero thickness geometry"* 라 부르는 조건과 같다.

**거절은 «무엇이 잘못됐는지»만 말한다.** 어느 일치가 의도치 않은 것이었는지, 답이 다른 치수인지 다른 연산인지 몸체를 둘로 나누는 것인지는 **작성자의 의도**이고 커널은 그걸 볼 수 없다 — 「치수를 조금 옮기라」는 그 의도에 대한 추측이다. `BoolReport` 가 같은 규칙을 적어 둔다: *"It is a diagnosis, not a prompt."* 커널은 사유 이름과 클래스만 내놓고 **사람이 읽을 문장은 앱의 몫**이다(`RejectClass` 로 분기, 변종 이름으로 분기하지 않는다).

자기접촉 판정은 체 셋을 싼 것부터 통과시킨다 — ① 후보 평면 = 두 끝점 이름 삼중의 교집합(정확, 좌표 없음) ② 그 평면 위 면들 사이에서 정점 `tol` 로 **부풀린** 경계상자 선별(한 평면이 면 80개를 일 수 있다; 후보 2,538,703 → 203) ③ 남은 것만 정확 판정. 상자를 부풀리는 것은 선택이 아니다: 상자는 반올림된 실현 좌표에서 나오고 정확 판정 **앞에서 후보를 버리는** 데 쓰이므로, 보수적이지 않으면 진짜 접촉이 조용히 사라진다.

**③은 끝점이 아니라 «열린 구간»을 묻는다.** 끝점만 보면 **관통 슬롯**을 놓친다 — 쐐기를 끝까지 뚫으면 접촉선의 두 끝이 면의 링 «위»에 놓이고(링 위 = 내부 아님), 가운데만 내부를 지나는 **현**이 된다. `combinatorics::segment_meets_face` 가 면의 **모든 링**을 그 선에 대고 읽어 교차점을 모으고, 선을 따라 정렬해 칸이 **밖/안/밖/안…** 으로 번갈아 간다는 사실로 안팎을 읽고, 열린 구간과 겹치는 「안」 칸이 있는지 본다. 갈래가 없으므로 구멍(모든 링을 한 자루 = 짝수-홀수)과 여러 칸을 저절로 처리하고, **끝점 판정을 흡수한다**. 선 `q ∩ w` 의 `w` 는 **모서리 자신의 두 면**에서 얻는다: 링이 기록한 벽은 다발에서 `q` 를 가리킬 수 있고, 이름 삼중은 4평면에서 두 번째 평면을 안 담는다(31쌍에서 확인).

### 병렬 — 판정만 병렬, 순서는 고정

불리언 시간의 90% 는 **면당이 아니라 평면 클래스당** 독립인 두 루프에 있고(트레이스→병합→분할, 셀→중첩→라벨→방출), 나머지 지배항은 모델의 정밀도를 **점마다** 읽는 `standard_for` 다. 셋을 `nacre-ops::par` 의 헬퍼 둘(`map_range`·`try_map_range`)로 병렬화한다(`parallel` 피처, rayon).

**원칙: 판정만 병렬, 변이·순서는 단일 스레드 고정.** `assemble_fuse_cut` 이 vertex handle 을 면의 first-appearance 로 배정하므로 `Store::push` 순서 = handle 정체성 = 재생 결정론이다. 그래서 모든 헬퍼가 **인덱스 순서로 수집**하고, f64/BigFloat 를 병렬 reduce 하지 않는다(항목 내부 계산은 순차). 스케줄에 의존하는 답은 한 번 틀린 답이 아니라 **매번 다른 모델**이다.

- **에러는 인덱스-최초.** 순차 루프는 가장 낮은 클래스에서 반환한다. rayon 의 `collect::<Result<_,_>>()` 는 단축평가라 *어느* 에러가 살아남을지 정해지지 않으므로 쓰지 않는다 — 전부 모은 뒤 순서대로 훑는다. 규칙은 `try_map_range` 안에 **한 번만** 산다.
- **별칭 표는 라운드 스냅샷으로.** 클래스는 4-평면 동시성에서 점의 별칭을 발견한다. 각 클래스는 **라운드 시작 시점의 표 + 자기 발견**을 보고, 라운드 끝에 클래스 순서로 흡수한다. 순차 루프(각 클래스가 그때까지의 발견을 본다)와 같은 고정점에 이르는 이유는 스케줄이 아니라 `Aliases` 의 성질이다 — 병합 대표가 **최소 원소**라 최종 분할이 병합 순서와 무관하고, 발견이 단조 누적이라 늦게 알면 라운드가 하나 더 들 뿐이다(4-평면 모델에서 순차·병렬 모두 2라운드·별칭 36개로 동일).
- **정밀도는 max 리덕션이라 안전하다.** `judge_precision` 은 점당 `trial_bound` 와 `precision_for` 로 쪼개져 있다. 결합이 **최댓값**이라 결합적·정확하다 — 부분 *합*이었다면 재결합이 다른 수를 만들어 이렇게 못 나눈다.
- **대가: 거절 입력에서 일을 더 한다.** 순차는 첫 거절에서 나머지 클래스를 안 돌고 병렬은 전부 돈다. 답은 같고 비용만 오르지만, **뒤 클래스에 잠복한 패닉이 도달 가능**해진다 — 거절 코퍼스가 그래서 관문이다.
- **임계값은 재서 넣는다.** `map_range` 는 항목이 점 하나의 실현이라 작은 모델(면 12장 ≈ 점 36개)에서 분배가 계산보다 비싸다(527→576µs). 측정된 크로스오버로 `PARALLEL_FLOOR = 64` 다(그 위는 fold 중간 불리언이 순차로 떨어져 큰 쪽이 손해). **비용 손잡이지 답 손잡이가 아니다** — 양쪽이 같은 값을 같은 순서로 낸다는 것을 테스트가 못박는다.

**범위 밖: 브라우저.** 플레이그라운드는 `default-features = false` 의 wasm 빌드라 이 이득을 못 받는다. `wasm-bindgen-rayon` + `SharedArrayBuffer` + COOP/COEP 는 배포까지 얽히는 별도 작업이다. 되돌리는 길도 코드 밖에 있다 — 그 피처 한 줄을 끄는 것이 곧 순차 동작이다.

### 상대가 닿을 수 없는 평면은 소식을 나르지 않는다

불리언은 *"이 평면들 전부의 배열이 무엇인가"* 를 풀지만, 답하는 질문은 *"`b` 가 `a` 를 어떻게 바꾸는가"* 다. 증분 fold 에서 둘은 거의 전부만큼 다르다 — 80핀 fold 의 마지막 fuse 는 평면 클래스 **169개를 배열하는데 새 핀이 닿는 것은 27개**이고, 나머지 142개에서는 **누적 솔리드가 이미 가진 면을 처음부터 다시 유도**한다. 그 불리언 시간의 **59%** 다.

그래서 **상대 솔리드와 떨어져 있음이 증명된 평면 클래스**는 배열하지 않고, 자기 솔리드의 면을 **배열 엔진의 어휘(평면 클래스 삼중항)로 다시 적어** 낸다. `nacre-ops::reuse` 한 모듈이고, `trace_result_faces` 의 출력 타입은 그대로 `Vec<LocalFace>` 라 **하류(조립·동일평면 병합·seam 표·검증기)는 이 지름길을 모른다**.

- **건너뛰기는 추정이 아니라 정리(theorem)다.** 분리는 `nacre-judge::orient3d_filter` 로 **증명**한다 — 오차 한계가 0 을 넘길 때만 부호를 답하고, 아니면 `None`. `None` 은 "그냥 배열하라"이고 그건 이 모듈이 없을 때의 동작이다. ⇒ **증명 실패는 속도만 잃고 답은 못 바꾸며, 톨러런스가 들어오지 않는다.**
- **필터는 새 술어가 아니라 기존 술어의 절반이다.** `orient3d_judge` 는 f64 필터 → 공유-모션 정확 경로 → 상승의 3단이고, 앞 두 단이 정확히 여기서 원하는 것이다. 그것이 `orient3d_filter` 이고 **판정이 그것을 호출한다**. 오차 예산을 아는 곳이 둘이 되지 않으므로, 건전성은 검사할 성질이 아니라 **같은 코드라는 사실의 귀결**이다.
- **정점만 봐도 되는 이유, 그리고 그 전제.** 평면까지의 부호거리는 아핀이라 모든 정점이 한쪽이면 **볼록껍질 전체**가 한쪽이고, 평면 다각형은 (비볼록이어도) 자기 정점의 껍질 안에 있다. **이 논증은 모든 면이 평면일 때만 성립한다.** 논증을 떠받치는 것은 **"원통 클래스가 하나라도 있으면 reuse 를 통째로 끈다"**는 가드(`trace_result_faces` 진입)다 — 그래서 reuse 가 도는 동안에는 "모든 면이 평면"이 참이다. 곡면이 있는 입력에서 이 지름길을 쓰려면 껍질이 아니라 **곡면의 경계**를 봐야 한다.
- **그 가드는 삼중 임무다 — reuse 를 만지기 전에 읽을 것.** ① 위 논증(평면 전제)을 살려 두고, ② `pass_through` 가 **평면 클래스 단위**로만 면을 옮기므로(어떤 클래스가 PassThrough 면 원통 측면이 **조용히 사라진다**) 그것을 막고, ③ **원통 밴드가 배열의 라벨을 읽을 수 있게** 한다(원통 밴드 절 — 라벨은 *배열된* 클래스에만 존재하고, PassThrough 클래스는 셀도 라벨도 만들지 않는다). 즉 원통 인구에 reuse 를 되살리는 일은 성능 작업이 아니라 **밴드의 소속 판정을 다시 설계하는 일**이다.
- **정점 이름은 솔리드가 이미 갖고 있다.** 배열은 정점을 정렬된 삼중 `[w, f, third]` — 그 점에서 만나는 세 평면 클래스 — 로 부른다. 솔리드의 정점은 자기 주위 면들의 클래스를 정확히 만나므로, **술어 하나 없이 같은 이름이 나온다**. 넷 이상이면 그건 동시성이고 별칭 표가 이름을 고르는 문제라 재현할 수 없으므로 그 클래스는 배열한다.
- **규칙표에 기하가 남지 않는다.** 상대가 못 닿는 평면은 모든 점이 상대 **바깥**이므로 `keep` 이 접힌다: `Fuse` 는 양쪽 다 통과, `Cut` 은 `a` 쪽만 통과·`b` 쪽은 비움, `Common` 은 둘 다 비움. 양쪽이 다 면을 가진 클래스(동일평면 접촉)는 언제나 배열한다.
- **두 패스를 모두 건너뛴다 — 그리고 그게 이 설계의 유일한 베팅이다.** 셀 패스(pass B)를 건너뛰는 것은 안전하다(면이 맞는지만 문제고, 그건 아래 차등 검사가 답한다). **추적 패스(pass A)까지 건너뛰면 엔진이 *아는 것*이 달라진다** — 추적되지 않은 클래스는 자기가 발견했을 동시성을 보고하지 않고, 커버리지 코퍼스에서 **별칭의 17% 가 분리된 클래스에서만** 발견된다. "그 이름을 아무도 안 쓴다"(건드려진 클래스의 교차점이면 그 클래스가 스스로 찾는다)는 **논증이지 측정이 아니므로**, 믿지 않고 검사한다.
- **차등 검사가 영구히 남는다 — 불리언 전체 수준에서.** 디버그 빌드는 같은 불리언을 `ClassReuse::{Off, Proved}` 두 번 풀어 **같은 면을 요구**한다(`Off` 는 폴백이 아니라 기준이다). 비교는 링 회전·면 순서에 무관한 정규형으로 한다 — 배열의 링은 DCEL 순회가 시작한 곳에서, 솔리드의 링은 저장된 곳에서 시작하므로 그대로 비교하면 차이가 아닌 것을 차이로 읽는다. 클래스 단위로는 비교할 배열이 없기 때문이고, 별칭 질문에 직접 답하는 것도 이 수준뿐이다. 참조 실행은 **자기 `Notes`** 를 받는다 — 증거는 `undecided_reject` 가 읽는 부작용이라 두 번 쌓이면 안 된다.
- **건너뛰기는 클로저 *안에서* 한다.** 범위를 줄이면 안 된다 — 인덱스 공간이 유지되어야 `try_map_range` 의 "최저 인덱스 거절 보고"와 방출 순서(=핸들)의 결정성이 지켜진다.

**두 검사 모두 공허하지 않다**(음성 대조): `flip` 을 뒤집으면 커버리지 **4건**이, 분리 증명 없이 소유 클래스를 전부 통과시키면 **167건**이 즉시 실패한다.

**왜 두 패스를 다 건너뛰는가 — 숫자.** pass A·B 를 모두 건너뛰면 순차 실행이 회전 fold 80 에서 **2.02배**(9.66s → 4.78s), 허브먼저 ring 80 에서 1.95배, 축정렬 fold 80 에서 1.85배 빨라진다; pass B 만 건너뛰면 9.66s → 7.69s 에 그친다. 80핀 fold 에서 클래스 방문 8800회 중 **7021회가 지름길**을 탄다. 결과는 두 피처 조합에서 비트 동일하다.

**병렬은 값을 치른다, 그리고 그게 이 설계의 값이다.** 남은 일이 169개 중 **27개 클래스에 몰려서** 병렬 효율이 약 4.1배에서 2.3배로 떨어진다(회전 fold 80 병렬: 2.38s → 2.09s, pass B 만일 때의 1.98s 보다 5% 나쁘다) — 총 일은 훨씬 더 줄지만 병렬은 그 폭을 못 따라온다. **순차를 택한 이유는 브라우저가 단일 스레드이고, 이 비용을 문제 삼는 자리가 브라우저이기 때문이다.** 병렬 효율 회복은 정확성이 아니라 스케줄링 문제다.

## 원통 — 밴드·게이트·차트·룰링 (`nacre-ops::{bands, arrangement::cyl_chart}`)

원통의 벽은 평면 클래스가 아니라 밴드로 정해진다.

배열 엔진은 **평면 클래스 하나씩** 답한다. 원통의 측면은 평면이 아니므로 그 어휘에 들어오지
않고, 그렇다고 평면 배열과 별개의 새 엔진이 필요하지도 않다. 출발점은 **균일-슬랩 정리**다:

> 연속한 두 ⊥ 절단 사이의 **열린 원통 슬랩**에 상대 솔리드의 경계가 전혀 없으면(⊥면은 절단 자신,
> ∥벽은 인구 게이트가 r 밖으로 밀어냈다) 상대 소속이 그 슬랩 전체에서 **균일**하고, 증인 하나가
> 밴드 전체를 정한다.

게이트가 **기록하고 통과시킨** (벽, 원통) 쌍에서는 그 전제가 깨지고, 그 자리에서 밴드는 잘린 원과
룰링으로 끊긴다(아래 「룰링」). 옆면은 결국 **자기 차트 위의 셀 복합체**로 답한다(아래 「차트와 셀 방출」).

이 절에는 **규칙과 불변만** 적는다.

### 밴드와 라벨

- **소속은 재지 않고 배열에게 묻는다.** 평면 배열은 원 경계마다 **원판 셀**을 만들고 `label_cells`가
  거기에 `[inA 위, inA 아래, inB 위, inB 아래]`를 쓴다 — 「원 안쪽에서 이 평면 위/아래에 각 솔리드의
  재료가 있는가」, 곧 밴드의 챔버다. 밴드는 경계면의 그 네 비트를 읽는다(`cyl_chart::read_cell::read_bits` 한 철자).
  양 끝을 다 읽고 일치를 요구한다 — 다르면 균일-슬랩 정리의 전제가 깨진 것이라 거절한다.
- **왜 광선이 아닌가.** 평면 엔진은 두 종류의 질문을 두 기계로 답한다: 「이 경계 조각이 결과에
  남나」는 **라벨**(`emit_faces`), 「이 성분이 재료인가, 어느 성분 안에 들어있나」는 **광선**
  (`point_in_component` — 한 평면이 답할 수 없는 3D 질문). 밴드는 경계 조각이므로 **첫째**다.
  라벨은 좌표도, 광선 방향도, 기권 재시도도, 폭 천장도 없다 — 그래서 원통이 **양쪽 피연산자**에 설 수
  있다(두-구멍 판). 좌표 도로 `point_in_faces_rat` 는 **둘째 질문**의 기계이고 이 경계를 뒤집지 않는다 —
  밴드는 라벨을 읽는다.
- **「z-범위 판정」은 성립하지 않는다**: L-노치의 안쪽 코너는 모든 벽에서 r보다 멀지만 재료 밖이다.
  그 반례는 좌표 도로와 밴드 패스의 픽스처다.
- **밴드 경계는 평면 배열이 실제로 방출한 원 경계**(+ 측면 자신의 rim)다. 임의의 ⊥ 클래스로
  자르지 않는 이유는 rim 공유다 — 캡 면의 원 루프와 밴드의 rim이 **같은 엣지 핸들**이어야
  닫힌-셸 가드가 사용 횟수 2를 본다.
- **조립**의 rim 기계: 통째 원마다 `OnSeam([측면, 평면])` 정점 하나와 자기루프 rim 엣지
  (`derive_edge_curve`가 원을 파생) — 밴드는 그 림 두 루프로 닫히고 seam 간선은 없다. 잘린 원은 노드 사이의
  호들이고 그 담체는 원통과 벽이 싣는 림의 평면이다(`Wall::Arc` 의 `plane` — 캡이 짓든 옆면이 짓든 한 쌍, 둘째로 걷는 면은
  그 평면을 대조한다). 성분 결합에는 **둘째 규칙**이 붙는다 — 캡과 밴드는
  Node를 하나도 공유하지 않으므로 **rim 키**로 잇는다.
- **감김은 규칙으로 유도한다**: rim 원은 축 방향에 대해 CCW이므로 `sign(면의 바깥 법선 · 축)`이
  외곽의 `forward`를 주고 구멍은 그 반대다. 잠금은 **부피의 부호**(뚫린 상자가 `8 − πr²h`).
- **딱지는 «소속»을 말하지 «존재»를 말하지 않는다.** 딱지는 「어느 쪽이 재료인가」를 말하고
  「이 옆면이 그것을 경계짓는가」는 **다른 질문**이다 — 구멍 난 옆면이 피연산자로 돌아오면 rim 의 두
  섹터가 **글자 그대로 같은 딱지**를 든다. 존재는 자국의 **종류**가 답한다(`Graze` = 면이 여기서
  끝난다 / `Transversal` = 관통한다). `ArcLabels` 가 딱지 옆에 기여 목록을 함께 나르고(한 호를 한 번
  볼 때 같이 채운다), `read_cell::face_spans` 가 `keep` **앞에서** 존재를 먼저 묻는다.
  옆면이 자기 차트에서 셀 복합체를 돌므로 존재와 소속은 **한 자리의 두 물음**이다 — `cyl_chart` 의
  셀이 `chamber`(소속, `read_cell::read_bits`)와 존재 읽기를 **나란히** 묻는다.
- **`Seated` 건너뛰기는 방어가 아니라 작동 중인 경로다** — 한 호가 옆면 자국과 평면 면의 `Seated` rim 을
  함께 들고 오는 인구가 실재한다(그것이 없으면 그 호의 답이 엉뚱한 면에서 온다). 기여가 **없는** 호와
  **둘 이상**인 호는 코퍼스 인구가 0이라 `face_spans` 의 **단위 테스트**(`cyl_chart` 의 테스트)가
  손으로 쓴 자국으로 진술하고,
  잠금은 「1이 아니면 빨강」이다.

### 인구 게이트 (`planes::cylinder_gate`)

- **게이트는 (평면 클래스, 원통) 쌍을 다섯 갈래로 가른다**: ⊥ 절단(통과, 상대 몸통에 flush 하게 앉은
  캡 포함) · r 밖으로 비킨 ∥ 벽(무기록 통과) · 평면이 r 안인 ∥ 벽(기록하고 통과) · 접선 ∥ 벽(접선으로
  기록하고 통과 — 단 접선 판정이 읽을 수 없는 행, 곧 그 선을 셋째 평면이 간선 아닌 채로 품거나 행을
  진술하지 못한 것은 게이트 끝에서 거절한다: `TangentLineInAnotherPlane`, 또는 행이 든 이름) · 비스듬한 평면(옆면의 모든 면이 그 평면을 비킴을 증명하면 통과, 아니면
  `ObliqueCylinderCut` — 타원은 짓지 않았다). 모든 산술은 세계 서술 위의 checked `Rat` 이고, 정확히
  결정하지 못한 것은 무엇이 막았는지로 이름이 갈린다 — 서술을 한 틀의 세계에 놓지 못함(돌린 클래스, 접히지 않는
  사슬 위의 원통)은 기하라 `CylinderGateUndecided`, 세계는 들지만 `Rat` 에 안 담김(`Wide` 이름, 넘침)은 폭이라
  `WitnessNotRational`. 보수적인 정직한 거절이지 추측이 아니다.
- **게이트의 «기록»이 열고 배열은 기록만 믿는다.** `point_plane_clearance_rat` **한 호출**의
  삼분법이 그대로 세 기록이다: `Positive` 무기록 · `Negative` → `crossings`(룰링 둘) · `Zero`(접선) →
  `tangencies`(스치는 선 하나 — 그 선 위에 놓인 옆면 **면**마다, 그 선 위에서 옆면이 있는 구간 안에서 선에 닿는 벽 면만; 구간은 옆면의 모든 루프의 림 원호에서 읽는다 — `lateral_cover_on_ruling`, 창이나 홈이 남긴 틈은 구간 밖이다). 기록 집합은 「평면이 r 안이고 면이 띠를
  비키지 못했다」의 **증명 캐시**라
  배열의 독자들은 전제를 `debug_assert` 로 든다. **접선을 `crossings` 에 넣지 않는 것이 핵심**이다 —
  넣으면 `crossings` 의 명제(「평면이 반지름 **안쪽**을 지난다」)가 거짓이 된다.
- **게이트는 면이 닿지 않아도 «클래스가 옆면을 자르면» 기록한다.** 거절 물음은 면 대 면이지만, 그
  클래스 위의 배열에는 룰링이 필요하다 — 이 솔리드의 면 중 그 클래스를 가로지르는 것은 전부 그 룰링에서
  끝난다(캡의 단면은 rim 호가 평면을 만나는 곳에서 끝난다). 단 **진술할 수 있는 발자국**만 절단을
  증명한다 — 진술 못 하는 발자국(회전된 클래스, span 없음)은 어느 쪽도 증명하지 않고 기록하지 않는다.
- **캡이 상대 재료 안에 앉는 변형은 기록되어 통과한 뒤 하류(챔버)가 거절한다.**
- **«두 ∥ 클래스의 교선이 옆면 면 위에 있다»도 게이트가 한 번 기록한다**(`SharedRuling`, 원통 행의
  `WorkingCyl::shared`). 짝은 게이트가 이미 적은 것에서만 나온다: `crossings` 둘, 또는 `crossings` 하나와
  **원통의 반대 피연산자 면을 가진** 접선 클래스(면 클리어런스와 무관하게 — 캡 위에 선 면도 그 선을
  준다). 같은 솔리드의 접선(필렛)은 그 클래스에 반대 피연산자 면이 없어 빠진다. 술어는 하나
  (`plane_plane_cylinder` 가 `OnRuling`)이고, 선은 **옆면 면의 각 범위 안**(`Footprint::theta_holds_line`,
  끝 포함)이어야 한다 — 무한 곡면 위의 선은 면이 없는 곳이다(필렛 사분원 밖에서 슬래브 모서리가 그렇다).
  **접선 행도 같은 술어로 쓴다**: 반원·슬롯 끝의 곡면이 옆면이 없는 쪽에서 벽에 접하면 아무것도 닿지 않는다.
  씨앗·룰링 분할의 접기·차트의 station·접선 판정·결과의 간선 분할이 이 기록을 읽고 스스로 판정하지 않는다.
  `Tangency::line_in_another_plane` 은 같은 계산(«그 접선을 품는 클래스들»)의 다른 읽기다.
- **게이트가 물어야 하는 명제는 «두 클래스가 «면»을 공유하는가»다** — 무한 곡면의 거리가 아니다.
  곡면 거리로 못 비킨 쌍은 **옆면이 자기 자신을 대변한다**: 한 클래스의 옆면 축 구간 전부가 상대의
  **도달 범위**(`lateral_reach` — 사영 구간 `d·o + s·(d·m)` ± 반경 도달 `r·|d⊥|`, 제곱으로 정확, 근호
  없음)와 서로소면 면을 공유하지 않는다. 수직은 그 일반식에서 `s` 항이 사라지는 경우이지 **분기가
  아니다**. 같은 솔리드 안의 쌍은 묻지 않는다 — 유효 솔리드의 두 면은 교차하지 않는다.
- **옆면의 서술은 차트의 발자국이다** — 축 구간 × 각도 범위, 림 호에서 유리수로. 원통 쌍의 게이트는
  그것으로 「두 클래스가 면을 공유하지 않음」을 증명한다(축이 `r₁+r₂` 보다 멀면 곧장 통과, 평행 쌍은
  「같은 곡면인가」(`same_surface`)를 먼저 묻고, 한쪽이 다른 쪽 안에 엄격히 들면 통과, 그 뒤 같은 문).
  면이 **만나는** 쌍만 진짜 새 기하(quartic 교선)이고 `CylinderPairContact` 로 거절한다.
- **갈라지는 방향은 «축 둘»만이 아니다.** 발자국이 말하는 상자의 경계는 캡 평면(법선 = 축)·
  룰링(방향 = 축)·호뿐이므로 **유리수 분리축 후보는 정확히 셋**: `m_a`, `m_b`, `m_a × m_b`.
  남는 것은 호의 **반경 연속체**이고 그건 공통 단면 차트가 있는 **평행** 쌍에서만 잡힌다 — 그래서
  거절은 참인 문장(「갈라짐을 못 보였다」)이다.
  **규칙은 하나다**: 「한쪽의 span vs 다른 쪽의 reach」를 양방향으로 묻는 대신, 자기 축을 따른 reach 가
  **곧 그 span** 이므로 **방향 목록 하나**로 묻는다(`lateral_faces_clear`). 비교도 한 문이다 — 두 reach 의
  양 끝이 각각 근호를 달므로 물음은 `g > √p + √q` 이고 근호 탑이 답한다.
- **발자국은 «꼭짓점»이 아니라 «경계 조각의 도달 범위»다.** 면을 꼭짓점으로 근사하고
  「면 ⊆ 꼭짓점들의 볼록 껍질」에 기대면 **원판**(외곽 루프가 원 간선 하나인 면)에서 전제가 빈다.
  꼭짓점은 폭 0 인 조각, 원판은 `중심 ± ρ` 이고 판독기는 **한 함수**(`corner_of`)다.
  **폭을 가진 조각은 점이 하나로 답하던 물음을 둘로 답한다** — 띠는 네 답(한 조각이 혼자 양쪽에
  걸친다), 축은 「점이 축의 어느 쪽인가」가 아니라 「`t` 너머 그 쪽으로 닿는가」다.
  **못 읽은 면은 «무죄»가 아니라 «모름»이다** — 읽지 않은 값이 판정에 들어가면 조용히 무죄가 된다.
- **호도 «조각»이다 — 「끝 둘이 이름 짓는 만큼 잘린 같은 원」.** 원판과 호가 한 변종
  (`Corner::Round { centre, rho2, axis, arc: Option<RimArc> }`)이고 세 물음의 답이 전부 **`arc_extent`**
  에서 나온다 — strip 은 `n × m` 방향으로, 축 물음은 `m` 방향으로. 그 규칙은 원통 옆면에 쓴 것과
  **같은 한 벌**이고 평면 면도 그것을 부른다.
  **끝의 순서는 «간선의 저장 순서»**(담체 원통의 축에 대해 CCW)에서 온다 — 면의 법선을 쓰면 **여집합**을
  이름 짓게 되고 그건 상계가 아니라 **다른 집합**이다.
- **게이트는 «원통이 만든 코너»를 근호 하나로 읽는다.** 점은 `line.base() + s·line.dir()` 이고
  상대(반지름·span·평면 계수)는 전부 **유리수**이므로 모든 양이 **`X + Y√c`** 하나 — 부호 탑이
  그대로 답한다. 근사도 평행 가정도 없다.
  - **재진술이 없다**: `QuadRoot` 는 정의가 **자기 두 평면을 자기 순서로** 말한 것에 대한 것이고
    여기서는 그것을 그대로 읽는다(아래 `ℓ` 보정은 class 표 순서로 **건너갈 때**의 값이다).
  - **거절은 원인별로 갈린다**(`CornerFail`): `CurvedOperandBoundary` 는 「뒤의 도로가 못 읽는다」가
    **참인** 자리(이음매 꼭짓점, 무리수 근에 끝나는 호)에만 쓰고, 서술을 세계에 놓지 못함은
    `CylinderGateUndecided`, 숫자가 `Rat` 에 안 담김은 `WitnessNotRational` 이다. 이름이 참을 말해야 다음
    사람이 그 이름을 믿는다.
- **피연산자가 굽은 이름을 나를 수 있다 — 연쇄의 규칙.** 불리언의 **결과가 다음 불리언의 피연산자**가
  되면 트레이서가 평면 어휘로 못 읽는 기하가 들어올 수 있다. 그 자리의 그물은 `CurvedOperandBoundary`
  이고, `loop_triples`가 **클래스로** 읽는 사실을 **담체로** 읽은 것이다(`Edge::surfaces` 의 doc 이
  「두 면의 adjacency 답」이라 말한다). 그물은 **그 도로의 예외(온전한 원)까지 함께** 나른다.
  ⇒ **게이트를 여는 것은 도로를 짓는 작업의 마지막 걸음이다** — 「그물이 덮는 만큼만 연다」는 등식은
  예외를 만난다(`OnSeam` 디스크 면 · subdivision 의 직진 정점 · 담체 등가의 동일 surface).
- **로드맵은 «케이스»가 아니라 «능력»으로 쓴다.** 「무엇이 막혔나」로 줄 세우면 다음 작업이
  그 케이스만 겨누게 되고 거기서 비스포크가 시작된다. 코드가 두 번 그렇게 말한다:
  `QuadRoot::Double` 이 저장 변종이므로 「오프셋 `0<d<r`」과 「접선 `d=r`」은 **한 이차식의 근 2개/
  1개**이고, 엔진이 **만드는** 경계 종은 넷(선분·원·룰링·현)인데 피연산자에서 **읽는** 종은 그
  부분집합일 수 있다 — 축은 「막힌 케이스」가 아니라 **이 넷을 읽는 능력**이다. 능력의 축은 넷이다:
  입력 어휘 = 출력 어휘 · 벽∩옆면은 한 이차식이다 · `validate` 의 비다양체 판정이 곡면을 본다 ·
  원통도 자기 차트에서 arrangement 를 돈다.

### 원·호·조각

- **`cylinder()` 는 원 스케치 + 돌출의 설탕이다** — 원시체 연산은 없다. 프리즘 빌더가 원·호를
  원통 조각 벽과 관통점 꼭짓점으로 세운다. 짓지 않은 것: 임의 각도의 호(점의 이름)·호–호 접합.
- **호의 진실은 끝점이지 각도가 아니다**: 스케치의 호는 `Edge2d::Arc { center, r2, ccw }` 이고 시작·끝은
  링의 정점이 말한다(전부 유리수; 온전한 원은 시작 = 끝). **회전각은 커널이 모르는 말**이고 끝점 표현이
  더 넓다. 옆면 방향은 「호의 ccw == 링의 감김」이 **아니라** **`ccw == 스윕이 +w`** 다(구멍 원의 옆면이
  반대로 뒤집힌다). 시계 방향 호의 캡 간선은 `[A,B]` 반시계 규약에 맞춰 정점을 바꿔 저장하고
  반간선이 되돌린다.
- **구간 돌출은 «스케치가 여기서 시작한다»다.** `[lo, hi]` 는 **한 번** 쓸되 프레임의 원점을
  **구간의 시작**에 둔다 — 두 번 쓸어 fuse 하면 한 곡면이 **두 번 진술**되고(축의 부호만 반대) 커널은
  그 쌍에 이름이 없어 거절한다. 「쓸고 옮기기」가 아닌 이유: 옮기면 결과에 **모션**이 남아 같은
  부피·같은 모양이 **다른 성격**이 된다. 계약은 `lo < hi` 만 요구한다.
  커널 쪽 「축 부호만 다른 쌍둥이」 정준화는 **짓지 않는다** — 그 인구는 위 규칙이 없애고,
  진술을 합쳐도 「한 곡면을 두 피연산자가 소유」로 옮겨 앉을 뿐이다(차트가 「원통 클래스는 한 솔리드의
  곡면」을 전제로 선다 — 그 전제를 바꾸는 것은 **별개 능력**).
- **원 셀은 «클래스가 축에 수직일 때만» 원이다 — 그 전제를 지키는 것이 게이트 하나뿐이므로
  규칙에 이름이 있다**(`class_carries_circle`, 총함수). 깨지면 거절이 아니라 **조용한 오답**이 된다
  (`Witness::On` 이 거짓이 되고 판정이 역방향 없이 답한다).
  앉은 원 생산자엔 더 강한 근거가 있다 — 닫힌 간선 하나에서만 나오므로 경계가 온전한 원이고,
  원은 자기 평면을 정하며 원통의 **원형** 단면은 축에 ⊥ 다(게이트가 바뀌어도 산다).
  **소비자는 각각 다르게 틀린다**: 림은 `û ⊥ axis` 가 `û ⊥ n` 을 함의하지 않아 평면을 벗어나고,
  디스크 포함은 **반지름**이 타원의 폭이 아니라 틀린다. ⇒ 거절 이름을 재사용하지 않는다 —
  「커버리지가 끝난다」(`NotSupported`)와 「닿을 수 없어야 할 백스톱이 말했다」(`SuspectedDefect`)는
  **다른 문장**이다.
- **한 클래스는 원과 룰링을 둘 다 들어도 된다.** 위험은 「둘 다 든다」가 아니라 「둘이 **가로지른다**」
  이다(가로지르면 어느 분할도 그 교점을 안 만들며 점 자체가 `평면 ∩ 원통 ∩ 원통` 이라 이름이 없다).
- **규칙은 «간선»에 대한 것이다**: *클래스의 간선은 어떤 면의 경계다 ⇒ 서로 다른 원통의 두
  간선이 만나면 두 옆면이 한 점을 공유하고, 쌍 규칙(`lateral_faces_clear`)이 그것을 부정한다.* 이 한
  문장이 **아무 split 도 자르지 않는 쌍 셋 전부**(원×원·원×룰링·룰링×룰링)를 덮는다. 그것을 «간선»이
  아니라 «원»에 읽으면 간선이 아닌 「원의 연장선」의 **유령 교차**로 거절한다 ⇒ 물음은 **호**에 대한
  것이고, 호를 실현 못 하면 온전한 원으로 **후퇴**한다(상계이므로 건전).
  그러면 그 검사(`CircleCrossesRuling`)는 **backstop** 이다(쌍 규칙이 먼저 거절한다). **그래도 출하한다** —
  만들어지지 않은 노드는 패닉이 아니라 **조용히 틀린 솔리드**이기 때문이다.
- **접선은 «가로지름»이 아니다.** 「닿음」을 묻는 그물은 설계상 **닫힘**(접촉도 닿음)인데
  배열이 못 만드는 것은 **가로지름**이다 — 접선은 아무것도 나누지 않는다. ⇒ 문이 **경계를 인자로
  들고**(`cylinder_ruling_reached_extent` 의 `touch_counts`) 그물은 **열린** 물음을 한다.
  **이름이 명제를 말하게 한다**: 거절 이름은 `CircleCrossesRuling` 이다 — 「만난다」라고 부르면서
  「가로지른다」를 묻는 표류가 버그의 기제다.
  **release 는 `debug_assert` 를 컴파일하지 않는다** — 「되던 것이 깨졌다」로 보이는 것이 실은
  조용히 지나가던 것이 프로덕션 거절로 드러난 것일 수 있다.

### 점의 이름

- **점의 이름은 «그 점을 지나는 평면 집합의 함수»이고, 그것을 짓는 함수는 하나다.**
  `combinatorics::canonical_triple(jd, s)` — 정렬된 클래스 집합의 사전식 최소 독립 삼중(셋이면 판정
  없이 그대로). 생산자 넷이 그것을 부른다: 피연산자 링(정점이 인접 면의 클래스로 **스스로** 이름
  짓는다) · 별칭 대표 · 선 위의 셋째 평면 코너 · 결과 정점 정의(입사 결과 면의 평면 집합 — 별칭
  대표와 **일부러** 다른 집합이다: 「결과에 없는 면으로 정의」를 피한다).
  **원장은 타입이다**: `Canon3` 뉴타입만 `NodeId::three_planes`·`Def::Three` 가 받으므로 다른 철자는
  컴파일 오류다.
- **원통 위의 점도 이름 하나** — 대표는 **관통점 이름 우선**이다(차트의 θ 기하가 이름에서 원통을
  읽는다). 그 규칙:
  1. **동일성은 «별칭 표»** — 발견은 씨앗 둘(피연산자의 관통 코너 × 클래스, `side_of == 0` —
     그리고 쌍대: 게이트의 `SharedRuling` × 축을 가로지르는 클래스마다, 삼평면 이름 · 할선의 관통 이름 ·
     접선 짝의 `Double` 이름), 병합 층 셋이 모든 이름을 정준화하고 소비자는 `canon_point` 로 비교한다.
     추적 중 발견은 같은 라운드의 다른 클래스가 먼저 거절할 수 있어 **버린다**.
  2. **핀은 (이름, 선)의 함수**(`pin_for`) — 이름이 접히면 **다시 유도**한다. 선 위의 점은 그 핀으로
     위치를 묻는다 — 관통 이름은 자기 짝이 부르는 선에서만 위치자다(`names_the_line`). 대표 이름이
     다른 쌍의 관통 이름인 점(공유 룰링 위)은 세 평면의 교점이라 유리수이고, 그 근은 유리수로 줄여
     읽는다(판별식이 제곱수인 근을 한 걸음의 다른 근호와 섞지 않는다).
  3. **클래스 원과의 교차는 «호 위일 때만» 절단점이다** — θ 순서를 먼저 매기고 덮인 것만 남긴다
     (안 그러면 온전한 원의 근이 진짜 정점과 겹쳐 유령 절단점이 조용히 쌓인다).
  4. **룰링이면서 평면 쌍 선인 선은 «옆면이 거기서 끝나는가»가 어휘를 정한다.** 끝나면(옆면의 제
     간선 — 필렛의 접선 이음) 축을 지나는 **할선 클래스 위에서는 평면 어휘**로 말하고 룰링 어휘로는
     침묵한다(`quiet_nodes`). 옆면 제 접선 벽의 클래스 위에서는 그 끝이 `side = 0` 룰링 정거장이고,
     다른 솔리드의 세그먼트가 그 선에 겹치면 아직 짓지 않는다(`lateral_crossings` 의 `OnRuling` 경비 —
     `RulingBoundNotYet`). 가로지르면
     (다른 솔리드의 모서리) 차트에 정거장이 있어야 하므로 **룰링 어휘**가 말한다 — 그 선 위의 평면
     세그먼트는 끝이 같은 룰링 조각에 기여를 넘기고 사라진다(`split_rulings` 의 접기; 겹치는데 안 맞으면
     이름 거절, 안 겹치면 옆면 밖이라 세그먼트로 남는다).
  5. **«어느 룰링인가»는 점과 벽의 함수**(`ruling_side_signed`)다 — 이름의 **근**에서 읽으면 그 이름이
     그 벽을 짝지을 때만 참이라 틀린다.
  6. **직선 간선은 두 끝이 열쇠다 — 평면 쌍 선이든 룰링이든**(`EdgeKey::Line`, 이음의 `JoinKey::Line`).
     담체는 **그 간선을 쓰는 두 면**에서 읽는다(패널은 원통, 평면 면은 평면) — 두 어휘가 한 식이라
     어느 면이 먼저 만들든 같다: {원통, 평면}이면 접선 이음도 룰링도 그 쌍, {평면, 평면}이면 옆면 없이
     룰링 위에 놓인 결과 간선. 한 평면 둘이면 `(P, P)` 로 적혀 결과 검사가 못 한 병합으로 거절하고
     (`CoplanarMerge`), 한 원통 둘이면 `(C, C)` — 어느 면에도 참이 아닌 철자다(seam 간선이 없다): 서로 다른
     두 패널 사이라면 남은 내부 경계인데 스캔은 곡면만 담아 둘을 가르지 못한다(인구 0, 스위트·census·무시된
     스윕·reject_census). 사용 수가 둘이 아니면 껍질
     가드가 거절할 간선이라 만드는 면의 벽과 제 곡면을 적는다. 끝점 이름의 공유 평면에서 유도하면
     「선을 담기만 하는」 평면을 적는다 — 조립의 자기 검사가 validate 의 규칙(적힌 담체 = 쓰는 두 면의
     곡면)을 그대로 묻고 `EdgeCarrierMismatch` 로 거절한다.
  7. **결과 정점은 닿는 면들의 평면으로 이름 짓는다 — 관통 정점도.** 이름의 두 평면이 다 닿으면
     그대로, 독립 평면 셋이 닿으면 정준 삼중, 아니면 닿는 두 평면의 관통 이름(`restate_pierce`, 근은 원래
     점이 놓인 평면으로 정확히 고른다). 공유 룰링 위의 점은 짝마다 관통 이름이 있고 대표가 결과에 면이
     없는 평면을 부를 수 있다.

  ⇒ **규칙은 하나다: 값은 «그것이 무엇인가»에서 유도하고, 이름은 이름일 뿐이다.**
- **굽은 링의 코너는 이름을 «유도»하지 않고 «재진술»한다.** 정점이 `Vertex::Pierce` 로 이미 들고
  있으므로, handle 공간에서 이 배열의 class 공간으로 옮겨 적고 담체는 결과 쪽이 이미 쓰는
  `boolean::Wall` 로 적는다(두 번째 어휘를 만들지 않는다).
  - **재진술의 규칙 — 두 보정은 «같은 규칙»이다.** `Lo`/`Hi` 는 `ℓ = n₁ × n₂` 방향의 순서이고,
    `plane_plane_cylinder` 는 base 를 `{n₁·x = −d₁, n₂·x = −d₂, ℓ·x = 0}` 로 정한다 — `−ℓ` 이
    **똑같이 만족하는** 조건이다. 그래서 **어느 한 법선의 부호를 뒤집으면** 두 평면도 base 도 그대로고
    `ℓ` 만 뒤집혀 **두 근이 자리를 바꾼다**. 쌍의 교환도 같은 이유다.
    ⇒ **`ℓ` 의 방향 반전이 홀수 번이면 root 를 뒤집는다.** 교환분은 `NodeId::pierce` 가 이미 센다.
    (원통은 보정 없음 — 축 방향의 부호는 곡면을 안 바꾼다.)
  - **대응과 부호를 «한 비교»가 답한다**: def 의 각 평면 handle 의 세계 계수가 어느 후보 클래스의
    계수와 **비례하는지** 보면, 비례하는 쪽이 그 클래스이고 **비례 상수의 부호가 곧 보정**이다.
    따로 구하면 두 곳에서 어긋날 수 있다.
- **정준화는 «네 계수»에 대한 것이고, 법선만 읽는 자는 그것을 물려받지 못한다.**
  `PlaneName` 의 정준화는 분모 털기 → **내용(gcd) 나누기** → 부호 규약이고 그 내용은 **네 계수**
  전체의 것이다 ⇒ 법선(셋)만 읽는 소비자는 `gcd(a,b,c)/gcd(a,b,c,d)` 를 떠안고, 차트가 그것을
  **제곱한다**. ⇒ 법선을 쓰는 자리는 `primitive_normal` 로 **자기 세 성분의 gcd** 를 턴다
  (양수 스칼라배라 방향이 그대로이므로 출력 무변은 측정이 아니라 증명이다).
  남는 인구는 「**원시** 법선 > 2⁶³」 이고 여유가 약 6비트다 — 총함수 답은 끝 질문(패리티=부호)을
  `BigInt` 로 올리는 것이다(계획).
  **차트 축은 «평면 안에서 직교»여야 한다** — 좌표 드롭은 퇴화 패턴을 안 보존하고, `e2 = ê_k` 는
  축을 3-D 평면 방향으로 쓰는 소비자에게 자기 평면 밖의 증인을 만들며, `e2 = ê_j × n` 은 `e1` 과
  직교하지 않아 광선의 등위선을 돌린다.

### 결정은 `def`, 측정은 `cache`

- 규칙은 `WorkingCyl::realized` 의 doc 이 든다(「the f64 twin … while **every decision reads `def`**」).
  축 방향 **결정**은 정확 서술로 내린다: `Wall::Ruling` 의 `up` 은 양 끝을 **끊는** 캡들의 축 매개변수를
  비교하고, 밴드 루프의 station 은 접점이 앉은 **cut circle 의 평면**이 축을 가로지르는 자리를 읽는다.
  `cache` 를 읽어도 되는 것은 **측정**뿐이다(공차·rim 의 중심).
  - **밴드 조립의 `axis_sign` 은 「어디서」가 아니라 「어느 쪽」을 묻는다**: 클래스 **프레임**이 축에 대해
    어느 쪽인가. `world_rat` 는 평면의 **이름**이라 방향을 안 나른다(스위트에서 `world_rat · m` 은 436회
    전부 양수, 답은 203회 반대). 방향은 이름의 **향**이 나른다 — `WorldName`(이름 + `sense`, 진실의 향과
    점들에서 `Model::world_plane_name_sense` 가 읽는다)이고, 축 부호는 `plus_t_is_above` 가 그것으로 정확히
    답한다. 평면 캐시의 법선은 어느 방향 판정에도 들지 않는다.
  - **`world_rat` 은 이름이지 방향이 아니다**: 룰링의 **정체**는 그것으로 읽어도 프레임 무관하지만
    **딱지**는 안 된다 ⇒ 방향이 필요하면 `WorldName` 의 향(`κ`)을 곱한다(`stored_coeffs_rat` 이 그 곱의
    한 철자다).
  - **f64 경로를 정확 경로로 바꿀 때의 규율**: census 비트 동일은 「출력이 안 변했다」이지 「두 경로가
    매번 같은 답을 냈다」가 아니다 ⇒ **차등 탐침을 먼저**(둘 다 계산·어긋나면 패닉), 코퍼스 전량,
    **탐침 자체를 부정 대조**, 그다음 f64 제거.

### 분할과 자국

- **분할 세 패스는 한 어휘를 쓴다.** 분할점은 `Split`(평면 클래스 **또는** 어느 원통의 어느 근),
  순서는 `order_located` **한 규칙**. 호 분할의 구간 판정은 `closed_contains` 이고, 좌표를 읽는 것은
  `segment_meets_cylinder` **필터** 하나뿐이며 양 끝에 좌표가 있을 때만 돈다(필터를 잃으면 풀이 한 번이
  들 뿐 답은 안 변한다). 룰링 분할의 방출부는 `split_segments_at` 으로 호 분할과 **한 함수**다.
- **옆면의 구멍은 «다른 모든 면과 똑같이» 이름 짓는다.** `loop_triples` 에 「면 자신이 원통인」 팔이
  있고 `lateral_cycles` 가 옆면의 루프를 모두 이름 짓는다(통째 림은 `LoopRing::Rim`). 자국이 구멍에 닿으면
  `DeclineKind::CylFaceHole` 로 **면을 지목하며** 거절한다.
- **「선 위 연속 구간」은 굽은 간선을 가로질러 이으면 안 된다.** 직선이면 두 점이 선을 정하니 참이지만
  **호면 거짓**이다(선을 떠났다 돌아온다). 그래서 걸음이 **구간을 departure 에서 쪼갠다** — 소비자
  셋(`trace_transversal_face`·모든 광선·자기접촉 탐지기)이 같은 전제를 공유하므로 규칙은 **걸음 안**에
  있고, 셋의 논리는 그대로다(조각이 곧 그들이 늘 뜻하던 「구간」이다).
- **패리티는 σ와 무관하다(정리).** 원과 평면은 두 점에서 만나므로 호의 **내부는 선을 다시 안
  만난다** ⇒ 호는 통째로 한쪽(σ). 그러면 교차 수가 `(s_a ≠ σ) + (σ ≠ s_b)` 이고 **패리티는 어느
  쪽이든 `s_a ≠ s_b`** 다. 다만 정리는 **개수**만 주고 **위치**는 안 준다. 위치가 필요한 자리에서는
  σ 를 쓴다 — 오라클이 있기 때문이다(재연산의 split-run 인구; 부호를 뒤집으면 `LabelConflict`).
  **떠나는 간선은 σ 쪽 «가상 노드»로 걷기에 들어간다**(`EdgeMeet::Departs(σ)`, σ = −`ruling_side`).
  **오라클 없는 부호는 만들지 않는다** — 부정 대조가 전부 초록인 채로 실제로 틀린 부호가 가능하다.
- **호의 교차 수는 두 끝이 아니라 호 자신에게 묻는다.** 원은 선과 최대 두 점에서 만나므로 두 끝이 반대쪽이면
  정확히 한 번이지만, 같은 쪽이면 0 또는 2(호가 부풀어 선을 두 번 넘는다 — 반원·슬롯의 둥근 끝, 비스듬히 자른
  필렛), 한 끝이 선 위면 0 또는 1 이다. 걸음은 그 두 자리에서만 소비자에게 묻고(`crossings`), 교차로 쪼갠 호의
  조각을 departure 처럼 부호 열에 넣는다 — 선 위 끝점의 옆은 그 조각이다. 평면 면 추적기에서는 함수
  하나(`arc_crossings`)가 수와 이름(호가 엄격히 품는 근, 진행 순서)을 함께 낸다: 평행·빗나감·접선은 0, 선 위
  끝점과 한 점인 근은 그 끝점이다. **호의 교차와 원 구멍의 현 flip 노드는 게이트의 `crossings` 에 실린
  쌍에만** 심는다 — 없는 룰링에 노드만 심으면 스퍼다(룰링 도로와 한 규칙).
- **술어는 통일하면 «버그»다 — 자취가 도로마다 다르다.**

  | 도로 | 클래스의 자취 | 「간선이 그 위인가」 |
  |---|---|---|
  | 평면 면 | **직선** `p ∩ q` | **호가 아닌가** |
  | 옆면의 구멍 | **원** `q ∩ 원통` | **담체가 그 원통인가** — 원 위엔 호가 얹힐 수 있다 |

  ⇒ 공유하는 것은 **규칙**(조각 쪼개기), **테스트는 각자의 것**.
- **`UnorderedEdges` 는 케이스가 아니라 «상류»의 물음이다.** 그 문장은 「한 정점에서 두 간선이
  같은 각으로 떠난다 — 상류의 `merge_coincident`/`Aliases` 가 접었어야 했다」이다 ⇒ 이 이름이 뜨면
  그 케이스를 여는 게 아니라 상류를 본다.

### 룰링

- **기록된 (벽, 원통) 쌍에서 밴드는 잘린 원에서 끊긴다**(`split_rims` 가 셋째 경계원). 양끝이 모두 잘린
  구간은 **θ-패널**로 갈라지고, 섹터별 챔버는 잘린 원의 **호별 라벨**을 읽는다(라벨 원칙 그대로,
  광선 없음). 룰링 엣지의 담체는 [벽면, 원통]이고 곡선은 ∥ 팔이 끝점에서 재생한다.
  소비자: props 는 θ-범위를 적분하고(닫힌 체인 = 정확히 2π), tess 는 열린 rim 체인을 감지 않고
  병합한다.
- **두 할선 벽이 진술하는 룰링은 정거장 하나다**(`chart_of` 의 `fold_shared_stations` — 기록의 짝을
  끝 이름으로 찾는다; 둘로 두면 `circular_order` 가 한 점의 두 이름을 거절한다). 두 벽의 딱지는 같은
  정보가 아니다: 원통 안에서 두 벽은 그 선에서 나가는 현 둘이고, 옆면의 `+θ` 접선에서 안쪽으로 돌 때
  먼저 만나는 현이 `+θ` 쪽 방을, 다른 현이 `−θ` 쪽 방을 닫는다(`StationTwin::bounds_plus` — 현 방향
  외적의 정확한 부호). 그리고 **두 현 사이의 방이 결과에 없으면 옆면은 그 선에서 양쪽 다 끝난다**
  (`Chart::slit_at` — 이웃을 잇지 않고, 한 성분 안에서도 경계로 낸다): 양쪽에 남는 재료는 그 선에서만
  만나므로, 한 솔리드면 네 면이 쓰는 간선으로 조립이 거절하고 따로면 두 솔리드다.
- **직교성은 «측정»이 아니라 «게이트의 귀결»이다.** 게이트의 다섯 갈래 중 옆면에 자국을 남기는 것은
  ⊥ 절단의 **원**(가로)과 ∥ 벽의 **룰링**(세로 — 교차면 둘, 접선이면 하나)뿐이다. 여섯 번째 갈래는
  없다. ⇒ 옆면 배열은 **주기 띠 위의 직교 «선분» 배열**이다(격자가 아니다 — 룰링은 유한 구간).
- **새 수학이 필요 없다** — θ 순서는 `circular_order_about_seam` 이고 그 시그니처가 두 점을 **같은
  원 위에서 받지 않으므로**(`(MeetLine, QuadVal)` 쌍) 서로 다른 벽의 룰링도 그대로 비교된다.
  오프셋 벽(`0 < d < r`)의 룰링 위치가 무리수(`√(r²−d²)`)여도 같은 산술(`QuadVal`·근호 하나)이 답한다 —
  산술은 지름을 가정하지 않는다.
- **룰링의 축-끝은 언제나 축-선이다** — 룰링 노드는 `[통과벽, 다른 평면]` 이고 그 「다른 평면」이
  옆면을 만나려면 ⊥ 일 수밖에 없으며 차트는 **모든** ⊥ 클래스를 모은다. ⇒ 배열이 **축-구간마다
  쪼개지고** 한 구간의 셀은 그냥 **θ-섹터**다(걷기가 필요 없다).
- **이미 정해진 값은 차트가 다시 계산하지 않는다** — 룰링 구간은 `rulings_on_class` 가 정하는
  값이라 차트는 그것을 **얹어** 쓴다. `rulings_on_class` 는 `TraceInput::crossings` 에 실린 쌍에서만
  발화한다.
- **이름은 한 번에 나온다** — `Chart::ruling_name` 한 번이면 되고 **파생 사본은 금지**다.
- **룰링의 identity 는 근이 아니라 `ruling_side`** 다. `RulingCarrier::side` 의 `±1`·`0` 은 **동일성 키**이지
  곱할 부호가 아니다.
- **접선 룰링은 `side = 0` 이라는 «이름»이다** — 캐리어를 만드는 두 자리에서 코너의 근
  `QuadRoot::Double` 로 유도하고 `ruling_side`(부호 술어)는 그대로다. 그 이름에서 따라오는 것들:
  직선–호 접선 이음은 「**반대 방향 = 반바퀴**, 같은 방향 = 곡률이 정하는 무한소 회전」(원호는 제
  중심 쪽으로 휜다: `tangent_codirected_turn` = 원호의 회전 `axis_up · ccw · frame_sign`) ·
  **부분 림의 덮이지 않은 호는 간선이 아니다**(차트 `End::Uncovered`) · 앉은 면이 접선 룰링 조각을
  스스로 낸다(마스크 = 축 쪽 비트) · 셋째 평면 위의 코너는 세 평면 이름.

### 차트와 셀 방출 (`arrangement::cyl_chart`)

- **옆면은 차트의 «영역»이다 — θ-패널이 셀이다.** 평면 쪽이 가는 길(클래스마다 셀 복합체)을 원통도
  자기 차트에서 간다. 단계는 평면 쪽과 같다: 방출 셀 → **연결 성분** → 경계 → 경계의 런을 이웃
  클래스의 조각으로 → 사이클 → 띠+구멍 분류. **«같은 엔진»은 같은 «단계»이지 같은 «코드»가 아니다** —
  차트는 환면·직교이며 딱지를 읽는다(평면 DCEL 을 이식하지 않는다).
  **정점은 이웃 클래스가 가진 곳에만** 만든다 — 직진하는 정거장에 정점을 만들면 벽 면의 간선 하나에
  옆면 간선 둘이 대응해 용접이 깨진다. 림에서 «가진 곳»은 **정리된** 캡 면이 가진 점이다(`HeldRims` —
  `held_rims` 가 정리 뒤의 평면 면에서 읽는다). 배열이 원을 나눈 표(`Curved::split_rims`)는 딱지를 읽는
  표일 뿐이다: 캡에 얹힌 상자를 빼면 정리가 캡을 다시 합쳐 그 점들을 지우는데, 나눈 표로 림을 자르면
  옆면만 거기서 끊겨 한 번 쓰인 간선이 남는다. 둘은 타입이 갈라 조립에는 나눈 표가 들어오지 않는다.
  림의 점이 **전부** 지워지면 그 캡 링은 제 원(`Bound::Circle`)이 된다 — 정리는 점을 정확히 하나만
  남기는 링을 만들지 않고, 원이 될 수 없는 링(띠의 사슬 림)은 점을 지키게 한다(`dissolve_straight_angles`).
- **클래스에는 «면의 구간»이 없다** — 면 단위로 자르는 연산을 클래스에 쓰면 「면 없는 자리에 띠를
  발명」하게 된다. 배열에서 그 틈은 **아무것도 안 남기는 셀**이다.
- **전파도 «셋째 부호»도 필요 없다** — `circle_on_class` 가 면 span 안의 모든 ⊥ 클래스에 원을
  남기므로 면이 있는 셀은 **양끝에 반드시 가로 딱지**가 있고, 챔버는 가로만으로 서며 존재는 가로 끝의
  자국과 행의 구간이 답한다. 병합은 θ 방향도 z 방향도 차트의 영역(연결 성분) 단계가 한다.
- **차트는 자기 «세로선»도 읽는다** — 셀의 챔버를 네 변에서 대칭으로 읽는다. 룰링의 딱지는 「벽
  양쪽의 재료」이므로 림의 딱지가 위/아래를 답하듯 **옆을 답한다**. **말하는 변이 모두 합의할 때만
  챔버가 서고, 어긋나면 이름으로 거절한다** — 어긋남은 배열 딱지의 결함이고, 한쪽을 믿으면 그 결함을
  덮는다. 세로는 가로가 말하는 칸에서는 검사자이고(판독 못 하는 벽은 침묵), 가로가 침묵하는 칸에서만
  답이 된다. 대칭은 챔버의 것이다 — 룰링의 자국은 계측 필드라 세로는 존재를 말하지 않는다.
- **섹터가 정거장을 가로지르면 «걸친 호들을 다 읽고 만장일치일 때만» 답한다.** z-선은 모든 ⊥
  클래스에서 오는데 θ-선은 옆면 스윕의 조각뿐이라 **셀이 실재하는 정거장을 가로지른다** ⇒
  `Chart::arc_around` 가 **가장 가까운 rim 노드 사이의 호들**을 돌려준다(θ 를 모르는 폴백은 과대주장이다).
- **호 라벨은 «셀에 직접 묻는 절대 판정»이다** — 저장 프레임 × 축으로 고르면 안 된다. 셀은 원을
  가로지를 수 없으므로 **유리 코너 하나의 반경 부호가 그 셀의 쪽**이다. 잠금은 결론이 아니라
  **전제**를 감시한다(`disk_side_probe`).
- **옆면도 자기 차트에서 «같은 패리티»로 답한다**: 원통 면의 `outer` 는 **모든 루프**이고
  (환면에는 「안」이 없어 총 패리티다) `loop_parity` 가 +z 광선으로 홀수 번 가르는지 답한다 — rim 은
  위에 있으면, 호는 위에 있고 span 이 θ_X 를 품으면, **룰링은 절대 안 센다.**
- **광선과 옆면의 교차는 근마다 면에게 묻는다**(`lateral_face_crossings`): 각 근이 `loop_parity` 로
  면 위인지를 **반쪽 판정 앞에** 묻는다 — 면 밖의 근은 어느 쪽이든 교차가 아니고, 광선의 원점에서
  면 위인 근은 `Graze` 다. **원통을 비켜간 선은 «0»으로 센다**(`CylinderMeet::Miss` → 교차 없음).
  패널·체인 rim·구멍도 건너뛰지 않고 읽는다.

### 감김과 방향

- **「어느 호가 구멍인가」는 링 자신의 감김이 답한다** — 「담체의 바깥 법선 반대쪽」이 **아니다**.
  `run_body_above` 가 쓰는 「**재료는 진행 방향의 왼쪽**」이고, 그래서 반대 면도 `edge → face` 표도
  인접 캐시도 필요 없으며 계단·U자 구멍까지 같은 걸음이 덮는다. 그 감김 규칙의 철자는
  **`material_theta_sign` 하나**다 — 부호를 뒤집으면 ⊥ 잠금과 ∥ 잠금이 **함께** 빨강이다.
- **호의 방향은 생산자가 벽에 싣는다**(`Wall::Arc` 의 `ccw`) — flank 에서 유도하면 **볼록 구멍에서만** 맞는다.
- **감김수는 «영역의 극점»에서 읽는다.** `loop_winding` 은 사전식 최소 **노드**에서 turn 을
  읽고 「노드 집합의 극점이니 껍질 꼭짓점」이라 논증하는데 **호가 있으면 거짓**이다(호가 모든 노드를
  넘어 부푼다). 영역의 극점은 노드이거나 **호 안의 한 점**이다. 호 안이면 `arc_extremum_winding` 이
  답한다 — 그 점(원이 걸치는 첫 세계축 `ê_a` 의 최소)은 **캡 평면 · 축을 담고 법선이 `m × ê_a` 인 평면 ·
  원통**의 pierce 점이라 노드와 같은 비교기(`cmp_key`)로 순서가 나고, 축이 기울거나 반지름이 무리수여도
  같은 도로다(기운 틀의 270° 부채꼴, `r² = 2` 부채꼴 — 둘 다 사전식 최소 노드가 오목한 중심이다). 그 점을
  못 풀면(폭) 노드의 회전으로 물러나지 않고 이름으로 거절한다. 극점이 노드이고 링이 거기서 **매끄러우면**
  (호–호, 선–호 이음) 회전 대신 **곡률**을 읽는다: 코너가 없어도 hull 꼭짓점은 국소적으로 볼록하고,
  winding = 호 자신의 회전 = **`ccw · axis_up · frame_sign`**(`smooth_extremum_winding` — 호 안의 극점도
  같은 곱이다). `turn` 의 `(Arc,Arc)` 팔은 `0` 그대로다 — **「회전」과 「곡률」은 다른 질문이다.** 극점이
  **첨점**이면(직선과 원호가 접해 링이 되돌아 나간다 — 원통을 뺀 열쇠 구멍의 캡) 영역은 원호의 **볼록한
  쪽**에 있으므로 winding 은 그 반대다.
- **부호가 어디까지 잠겼나**(인자를 하나씩 빼서 확인). 호 방향 둘·`plus_t_is_above`·루트
  재진술·`ccw`·`axis_up` 은 뒤집으면 전부 빨강이다. `frame_sign` 과 옆면의 `orient_sign` 은 **잠기지
  않았다** — 구멍 인구에서 둘 다 `+1` 이고(큐보이드 면은 `Reversed`가 아니고, 보스의 옆면은 바깥을
  본다), **bore 의 옆면에 구멍이 생겨야** 돈다. 잠금 doc 이 그렇게 적는다.
- **모서리 = 반간선의 시작**(2-간선 캡에서 두 정점을 보면 안 된다).
- **광선 위의 코너는 «위로 떠나는 스텝»이 센다**(혼합 링의 반열림 규칙). 직선 변은 다른 끝이
  엄격히 위, 호는 그 끝의 접선이 위 — 이웃 둘이 반대쪽이면 1, 같은 쪽이면 0 또는 2. 결정부는
  `nacre_geom::intersect::ray_step_crossing` **한 철자**이고 평면·혼합 도로 넷이 그것을 읽는다.
- **이름-도로가 혼합 링을 직접 읽는다**: `Ring::edges` 가 walls 에서 캐리어·핀을 짓고, 격리·좌표
  도로는 `ring_is_mixed` 로 갈라 유리 탐침 + `point_in_mixed_ring` **한 철자**로 간다. 병합은
  (방향쌍, 벽) 멀티셋의 **2-패스 지움**이다.

### 증인과 중점 (`nesting`)

- **셀이 셀 안인가는 «증인 하나»가 답한다 — 공급 하나, 판정 하나, 루프 하나.**
  그 물음 전부가 모듈 하나(`nesting.rs`)에 있다. `Witness` 는 **경계 위**냐 **내부**냐를 **타입으로**
  말한다 — 경계 위 증인 하나면 결정되고, 내부 증인은 **역방향을 물어야 한다**(상대의 경계 증인으로
  한 단계만). 공급은 게으른 하나(이름 → 유리수 코너 → 관통 코너 → 현 중점 → 간선 내부 → 내부 증인)다.
  증인 원자는 **모듈 비공개**라 다른 공급은 컴파일러가 막는다.
- **증인의 실패는 세 사건이고 이름이 셋이다**: `RingHasNoWitness`(증인이 애초에 없었다) ·
  `NoClearRay`(가진 증인이 전부 기권했다) · `WitnessNotRational`(값을 못 만들었다). 증인을 **여럿**
  쓰려면 「기권」과 「실패」를 갈라야 한다(「이 점이 링 위」는 다음 증인이 remedy, 「값을 못 만들었다」는
  아니다). 그 사이에 넷째 사건이 있다 — **도로가 대상을 못 읽는다**(`Said::Declines`): 광선 도로는 관통
  꼭짓점의 쪽을 원통 표 없이 물으므로 관통 꼭짓점이 있는 링은 어느 광선 평면에서도 못 읽고, 변이 둘인 섞인
  링(현과 호)은 걷지 못한다. 그 도로의 남은 증인은 같은 물음을 되풀이할 뿐이라 건너뛰고 좌표 증인이 답한다.
  아무도 답하지 못하면 이름은 기권과 같다(`NoClearRay`). 광선 도로의 거절을 읽는 자리는 하나다
  (`Said::of_ray_error`) — 증인의 사정도 도로의 사정도 아닌 것(`RingNaming`, 직선 변만의 `DegenerateRing`)은
  삼키지 않고 올린다.
- **유도가 소비자 둘을 가지면 자리를 옮긴다** — 사본을 만들지도, 「증인 원자는 모듈 비공개」를 열지도
  않는 유일한 길이다(`circle_centre_rat` 은 `nesting` 밖 `combinatorics` 에 있다).
- **「간선이 이름 짓는 점」은 한 규칙이다** — `combinatorics::edge_interior_points` 하나,
  팔 셋(켤레 두 근 · 한 직선 위의 두 관통점 · **양 끝 유리수**)이고 **이터레이터**다(한 간선이 점을
  둘 낼 수 있어 접으면 증인이 사라진다). 성분의 3-D 깊이도 같은 공급을 쓰되 **꼭짓점이 답하면 비용 0**
  인 지연 2단계로 읽는다.
- **링의 간선이 자기 중점을 준다**(`conjugate_midpoint`): `plane_plane_cylinder` 가 두 근을 **한 `mid`
  에서** 만들므로 두 끝이 그 쌍인 간선의 중점은 `base + s.a()·dir` — **온전한 유리수**이고 `disc > 0`
  이라 엄격히 사이다. 켤레성은 값이 아니라 **이름**으로 판정한다(같은 canonical 쌍·같은 `cyl`·뿌리
  `{Lo,Hi}`) — 다른 형상에 잘린 현의 **조각**은 한쪽 끝 노드가 달라 걸러진다.
- **원은 «자기 림»에서 증인을 댄다.** 「원에는 노드가 없다」는 **배열의 노드**에 대해 참이지만
  거기서 「이름 지을 경계 점이 없다」가 따라 나오지 **않는다** — 원의 경계 점은 기하적이고 정확히
  유리수이며, 커널이 이미 그 점을 만든다(온전한 원의 `ref_dir` 가 `(정점 − 중심)/radius` 라
  `centre + radius·ref_dir` 가 **그 링의 정점 그 자체**다). **배열이 버리는 것은 노드지 점이 아니다.**
  ⇒ 디스크 셀은 **림 넷(`On`) → 중심(`In`)** 을 지연 공급한다.
  **넷이고 어느 쌍도 못 뺀다**: 림 점이 중심의 광선 선 위인 것은 `û ∥ e₁` 일 때뿐인데 `û₁ ⊥ û₂`
  가 그 평면을 펼치므로 **최대 한 쌍만** 그 선에 놓인다 ⇒ 항상 결정하는 쌍이 하나 있다.
  프레임은 **투영하지 않고**, `ref_dir ⊥ dir` 이 아니면 **기권한다**(모든 소비자의 건전성이
  `û·dir = 0` 을 «사실»로 요구한다).
- **호로만 닫힌 링은 그 원 자체이고, 증인은 중심이다**(`ring_own_circle`·`circle_centre_rat` —
  축 방향과 무관하게 유리수). 원통에서 잘린 링은 코너가 전부 관통점이라 코너 기반 탐침 목록이 빈다.
- **잘린 캡도 원판이고 자기 내부의 점을 이름 짓는다** — 「이 링이 온전한 원인가」가 아니라
  **「이 면을 어느 원이 두르는가」**를 묻는다. 어느 후보가 면 안인지는 **링에게 묻는다**(유도하려면
  「현의 재료 쪽」 부호가 필요한데 그 도로엔 그 오라클이 없다). 후보는 차트의 두 축 걸음이다: 현이 한
  축과 나란하면 다른 축의 광선이 seam 에서 원과 만나는데, seam 위의 근도 순환 순서의 첫째로 판정되므로
  (`SeamOrder::seam_first`) 대각 걸음이 필요 없다.
- **링의 닫힘은 암묵적이고, 소비자가 그것을 알아야 한다.** 커널의 관례는 「링은 첫 점을 끝에
  반복하지 않는다」이다(메시 정점이라 반복하면 퇴화 삼각형). ⇒ 닫는 구간은 소비자가 만든다
  (`edge.vertices[0] == edge.vertices[1]` 이 그 사실이고, 계약은 `Tessellation` 의 말로 적혀 있어야 한다).
  그 자리의 심은 위반은 **두 겹으로 눈이 멀 수 있다** — 「차수 1 없음」은 **남는 구간**에 구조적으로
  무력하고, 픽스처에 호 간선이 없으면 닫힘 검사 자체가 미검증이다.

### 접선·접촉·집힌 면

- **접선은 «가르지 않는다».** 게이트는 접선 쌍을 거절하지 않고 `tangencies` 로 기록한다. 배열은 거기서
  아무것도 보지 않는다 — 그 선은 이 평면의 어떤 셀도 나누지 않고 차트의 어떤 섹터에도 정거장을 세우지
  않는다. **셋을 가르는 것은 연산이고, 연산은 `keep` 으로만 들어온다** — 접점 근처 재료는 렌즈·쐐기
  둘·건너편 셋이고 각 영역의 `(in_A, in_B)` 는 벽 면의 재료 쪽과 옆면의 `orient_sign` 두 비트가 정한다.
  ⇒ **연산 이름을 한 번도 안 쓰고** 세 답이 나온다.
- **«갈라진다»가 곧 결함은 아니다.** 바깥 접선의 Fuse 는 **선으로 닿는 유효한 몸통 둘**이다.
  결함은 그 둘이 **다른 데서 이어져** 있을 때다 ⇒ 판정은 「두 덩이 ∧ **한 출력 솔리드**에」이고, «한 솔리드»는
  덩이의 모양이 말한다: **쐐기 둘**은 선이 내부를 지나는 한 벽 면이 둘 다 경계하고 접선은 면을 가르지 않으므로 늘 한
  성분이다; **렌즈와 건너편**은 `keep` 상 바깥에서 접한 보스의 합집합에서만 생기고 두 피연산자 각각의 재료이므로,
  피연산자가 한 덩이인 한 결과가 **한 몸통**일 때만 한 솔리드다. 판사는 모양만 읽고 몸통이 어느 클래스의 면을
  드는지는 묻지 않는다 — 행이 참인 접촉을 말하기 때문이다: 게이트는 그 선 위에서 옆면이 있는 구간 안에서 선에 닿는 벽 면에만
  행을 쓰고, `runs_through` 는 선이 그 면의 **내부**(구멍·홈까지 읽은)를 지나는가다(`line_runs_through_face` — 선
  자체를 광선으로 삼은 point-in-polygon, 반열림 규칙은 `ray_straddle` 한 벌). 바깥 루프의 꼭짓점이 선 양쪽에
  있다는 것은 볼록하고 구멍 없는 면에서만 그 명제이고, 클래스 읽기만으로는 먼 곳의 같은 평면 면이 증거가 되어 두
  몸통을 거절했다. 그래서 판사는 게이트가 아니라
  `assembly::self_touch::tangency_reject`(민팅 전, grouping 이 손에 있는 곳)이고 이름은 `SelfTouchingResult` 다.
- **접선의 선이 셋째 평면 안에 놓이면**(`Tangency::line_in_another_plane` — 지름과 같은 폭의 상자,
  열쇠 구멍) 주변은 여섯 영역이고 세 영역 판정은 말할 수 없다. 그 선이 접선 벽 면과 할선 면이 나누는
  **피연산자 간선**이고 할선이 그 선을 룰링으로 진술하면(`Tangency::line_is_an_edge` — 게이트가
  `Edge::surfaces` 로 읽는다) 선은 배열의 간선이다: 옆면·접선 벽 면·할선 면이 모두 거기서 끝나므로
  두 덩어리의 접촉은 네 면이 쓰는 간선이고, 판정은 그 행을 **구조에 맡긴다**. 그 밖의 행은 게이트가 배열 전에
  거절한다(`TangentLineInAnotherPlane`) — 판정이 읽지 못할 행을 들여보내면 배열이 먼저 넘어지고, 어느 백스톱이
  나가는지를 클래스 순서가 정했다.
  옆면의 캡 너머로 뻗은 그 간선은 조립의 분할이 캡의 림 점에서 자른다 — 그 점의 이름은 캡과 한 벽의
  관통 이름이라 «짝 = 간선의 두 평면»이라는 이름 사실로는 안 보이므로, 기록된 선마다 점에게 묻는다.
- **접선 접촉에서 집힌 것은 «표면»이 아니라 «면»이다.** 접점 둘레 작은 구 위에서 경계를
  따라가면 **링크가 원 하나**다 — 집힌 면의 두 로브가 이웃 **곡면**을 돌아 이어진다. ⇒ **표면은
  2-다양체**이고 그 결과는 유효한 솔리드다(OCCT 도 같은 합집합에서 부피·면적·면 수가 일치하고 집힌
  면을 쪼개지 않는다). 잠금은 `bands` 의 `the_touch_is_a_slit`(접점에 위상 정점 0 **∧** 닿는 원이
  온전한 `[v,v]` 림).
- **χ 는 그 물음의 증거가 아니다.** 「핀치 하나는 χ 를 홀수로 만든다」는 **정점인 핀치**의
  문장이다. 접선은 정점을 만들지 않으므로 `euler_counts` 가 그 점을 아예 못 세고, 짝수 χ 는
  「핀치 없음」과 「위상 없는 핀치」 **양쪽과 양립**한다 — 판별력 0. 짝수 χ 를 증거로 읽으면
  거짓 규칙이 된다.
  **처방은 b-rep 이 아니라 메시다.** 면 분할은 **유효한 솔리드의 b-rep 을 바꾸는** churn 이고 OCCT
  도 안 한다 ⇒ 답은 **스윕의 기호적 순서(SoS)**다. 걸린 것이 작지 않다: `tessellate` 는 모델
  **전체**를 돌며 첫 실패에서 멈추므로 **집힌 면 하나가 그 세션의 모든 몸통을 지운다**.
- **집힌 경계를 그리는 데는 «셋이 함께» 필요하다**: ① 접점을 **공유되는 직선 간선의 폴리라인에
  끼워** 양쪽 면이 같이 보게 하고(pre-pass, 정점을 안 민다), ② 두 링을 그 점에서 **인덱스 수준으로
  이어 붙여**(다리 길이 0) 같은 좌표의 쌍둥이 인덱스를 남기고, ③ 스윕의 동점을 **쌍둥이에 한해 인덱스
  인식 기호 순서**로 정한다(`polygon/sos.rs` — 접선 방향 섭동, 1차 항은 `Expansion` 으로 정확,
  **좌표는 안 움직인다**).
- **`TessError::SelfTouchingBoundary` — 어떤 정점도 자기와 무관한 경계 조각 위에 있으면 안 된다.**
  `monotone::link` 가 거절하는 「조합적」 pinch(인덱스 공유)의 **기하학적 쌍둥이**이고 `link` 바로
  다음에 선다(`prev`/`next` 가 인접 제외를 공짜로 준다). 새 술어 0 · tolerance 0.
- **접촉은 자기 자격을 증명한다 — `self_touch` 가 교차를 가리면 안 된다.** 접선 증인 하나:
  접촉 정점의 두 이웃이 접촉 선분의 **같은 쪽**이면 접선, 반대쪽이면 교차(`Interior` 에서만 완전,
  `AtEnd`·공선 이웃은 기권). 교차 루프는 **언제나 돌고** `straddles` 는 엄격하다 — 셋은 **원자적**이다.
- **tess 의 거절 이름은 «메시에 대한 말»로만 쓴다**(교차 = 어떤 삼각분할이든 틀린다 / 접촉 = 이
  분해법에 답이 없다). **솔리드 판결은 `validate` 의 것이다.**

### 메시로 재기

- **메싱은 `tessellate` 로 잰다.** `nacre_tess::to_obj` 는 doc 이 스스로 「bootstrap · 직접 평면
  삼각분할」이라 적은 작성기라 굽은 면을 받지 않는다 — 그것으로 재고 「테셀레이션이 안 된다」고
  적으면 **틀린 도구로 잰** 것이다.
- **굽은 면의 메시를 보는 눈은 «넓이»다.** `validate`·watertight·정확 부피·면 수는 옆면의 면
  구조가 바뀌어 **메시가 조용히 틀려도 전부 초록일 수 있다**(넓이 7.873 vs 정확 9.425 인 메시가 그
  넷을 전부 통과한다). 그걸 잡는 오라클은 **메시 넓이 ≈ 정확 넓이** 하나다 ⇒ 옆면 표현을 바꿀 때는
  **먼저** 원통 코퍼스에 그 오라클이 서 있어야 한다.

### 강체 운동과의 교환

- **법칙 — 불리언은 강체 운동과 교환한다**: `transform(rigid(R, t)) ≡ transform(T) ∘ transform(R)` — 세계가 같다. 저장은 각 걸음이 옮겨 적을 수 있는 것을 따르므로 진술의 점 셋은 `invariant` 재진술로 갈릴 수 있고, 노드 구조는 넘침에서만 갈린다(`an_overflow_splits_the_representation_not_the_world`).
  **frame −1 의 메커니즘은 회전이 아니라 세계 평면 씨앗이다** — 이동 하나로 ∥ 딱지가 전부 뒤집힌다.

## 원통의 진실과 seam 규약

### seam 엣지 방식

주기 곡면(원통·구·토러스)의 옆면을 b-rep 로 담는 방식은 둘이고, 둘 다 유효한 AP242 이며 유효한 CW-복합체다.

- **B — 림 두 루프(커널의 표현).** seam 엣지 없이 위·아래 원 두 개로만 옆면 경계를 이룬다(옆면이 두 루프 `[outer=아래 원, inner=위 원]`). 오일러는 `V−E+F−L_i = 2−2+3−1 = 2` 로 통과한다. 커널의 모든 생산자 — 원 프로파일의 돌출(프리즘)과 불리언 조립 — 가 이 모양으로 짓는다: 바깥 루프는 차트에서 감김 +1 인 고리(아래 림), 안쪽 루프는 위 림과 구멍이다. 림은 통째 원이면 닫힌 간선 하나(그 정점이 θ = 0 의 `OnSeam`), 잘렸으면 호와 룰링의 사슬이다. STEP 출력은 내부 표현 그대로다.
- **A — seam 엣지(짓지 않는다).** 위 원 + 아래 원 + **세로 seam 직선 엣지**, 옆면 = 4-엣지 단일 닫힌 루프 `[bottom, seam, top⁻, seam⁻]`(seam 이 같은 면에서 반대 방향으로 두 번 쓰인다 — self-adjacent). 「가지 말 것」에 한 행.
- 정점조차 없는 제3의 형태(끝점 없는 독립 전체 원)는 오일러가 깨지며, **커널은 그 형태를 표현하지 않는다.**

A 를 받치던 다섯 근거가 이 커널에서 서지 않는 까닭:

1. **불리언.** 조합적 b-rep 불리언(면 루프를 순회하며 교차곡선에서 엣지를 쪼갠다)은 주기 면이 실제 엣지로 둘러싸인 유계 영역이기를 요구한다. 이 커널의 옆면은 그렇게 풀리지 않는다 — 자기 차트 위의 셀 복합체로 답하고(`cyl_chart`), 영역의 경계는 차트 고리 그대로이며, 바깥 고리는 고리의 감김(호 단위의 ccw − cw)이 고른다. 각도 비교는 끊는 자리와 무관하다(`SeamOrder::seam_first`). 업계는 갈린다: OCCT 는 seam 이 필수라 seam 없는 면을 읽으면 끼워 넣고(`ShapeFix_Face::FixMissingSeam`), ACIS 는 `periodic_no_seam` 으로 두 루프 면을 허용하고, Parasolid 는 감는 루프 한 쌍이 네이티브다.
2. **tess 정합.** nacre tess 는 공유 엣지 polyline(`by_edge`)을 경계로 소비해 crack-free 를 얻는다(테셀레이션 절). 옆면을 펼치는 절단은 tess 가 스스로 한다 — 감는 두 림을 한 생성선에서 자르고 그 점을 지나는 루프의 공유 polyline 에 끼워(이웃 면도 같은 점을 본다) 직사각형으로 펼친다(`cut_seamless_bands`). 절단선은 가능한 한 적은 루프를 지난다: 두 루프가 모두 닫힌 간선 하나면 각자의 정점에서(새 점 없음), 아니면 각 림을 한 번·구멍을 0 또는 2 번 지나는 생성선 가운데 구멍을 가장 적게 지나고 기존 샘플에서 가장 먼 것; 지나는 구멍은 끊어 절단선을 따라 바깥 고리에 잇는다(seam 없는 모델러의 패싯팅 — Stallings, *Reduced Topology B-Reps* §10). 매개화의 이음은 캐시(메시)의 일이고 진실(위상)에 오르지 않는다.
3. **«가짜 모서리»의 값.** seam 엣지는 기하가 아니라 «각도를 어디서 재기 시작하나»(`ref_dir`)를 위상에 새긴다 — 잘린 림마다 θ = 0 정점과 그 자리의 호 쪼개기, 구멍이 그 자리에 걸리면 바깥 걸음에 잇는 slit, θ = 0 에서만 보는 집힘 그물, self-pair 담체 규칙이 따라붙고, 같은 모양이 seam 의 자리에 따라 다른 위상이 된다(벽 위 보스의 노치가 seam 에 걸리면 바깥 걸음에 이어지고 아니면 구멍이 된다). B 에서 그 자리는 통째 림의 정점 하나에만 닿는다 — 네 벽의 보스가 한 철자로 나온다(`a_boss_on_a_wall_has_one_lateral_face`).
4. **위상 균일.** 옆면은 구멍 있는 평면 면과 같은 틀 — 바깥 루프 하나와 안쪽 루프들 — 이고, validate 의 오일러는 `L_i` 항을 이미 든다. 어떤 간선도 한 곡면을 자기 자신과 가르지 않는다(`push_edge` 가 거절하고 validate 가 잡는다). 한 점에서 고리가 집히는 차트 모양은 위치와 무관하게 세어 스위트 16(전부 쐐기 점수판)·census 0 이고 전부 껍질 가드가 거절한다(`NonManifoldResultEdge`) — 3-D 에서 선으로 닿는 모양이기 때문이다.
5. **유리수 차트.** 원통 위 점의 각도 매개변수는 단위원의 유리수 반각(Weierstrass) 매개변수 t 이고, 그 차트는 원에서 **정확히 한 점을 못 덮는다**(t=∞). 그 배제점은 그대로 θ = 0 이다 — `ref_dir` 은 모델 기하이고 통째 림의 정점이 거기 있다. 각도 순서는 그 점을 첫째로 읽어(`seam_first`) 앞/뒤 어느 쪽에 끊어도 답이 같다; 사라지는 것은 그 점을 지나던 «간선»뿐이다.

**어댑터 분리.** 내부 표현 ≠ 교환 표현이다. nacre-step 은 내부 표현을 그대로 내보내고, OCCT 는 seam 없는 옆면(위·아래 림 정점이 한 생성선에 있지 않은 것까지)을 checkshape 유효로 읽으며 부피·불리언 결과가 seam 있는 파일과 같다(DRAWEXE 로 잰 것). 그래서 내보낼 때 seam 을 넣거나 빼는 어댑터는 없다.

### 원통의 진실

1. **진실 형태.** `nacre_topo::Surface::Cylinder { def: CylinderDef, motion }`, `CylinderDef { origin, dir, ref_dir, r2 }` — 전부 유리수다. 반지름은 **제곱** `r2` 로 들고 그 타입은 **넓은** `BigRat` 이다: 축에서 곡면 위 유리수 점까지 거리의 제곱은 유리수지만 반지름 자체는 아닐 수 있고, 어떤 정확 술어도 제곱하지 않은 반지름을 묻지 않는다. 넓은 이유는 진술된 `Rat` 의 제곱이 `i128` 에 안 들어갈 수 있고(16자리 십진 반지름), 진술은 그 제곱의 폭 때문에 거절되지 않기 때문이다. 수로서의 반지름은 실현이다(`radius_f64`, 제곱근이 유리수면 `radius_exact`).
   - `dir`·`ref_dir` 은 **정규화하지 않은 원시**다 — normalize 는 정확형을 파괴한다(`normal_def` 와 같은 규칙). 캐시는 실현이다 — 세계에 진술되면 문이 진술에서 원점·단위 축·단위 `ref_dir`(축에 수직인 부분)·반지름을 각각 정확 반올림하고(`push_cylinder_raw`, `nacre_exact::cyl_unit_frame_f64`), 진술이 없으면(사분각 밖 회전·프레임 노드) ops 깔때기가 사슬을 고정밀로 재생해 같은 넷을 반올림한다(`push_cylinder_realized`; 방향은 재생한 두 점의 차). 원통 위 seam 정점도 같은 재생으로 세계에서 캡과 만나 실현된다(`seam_point_met`).
   - 성분형이 옳은 이유: 사용자 어휘의 직접 리프트이고 곱이 없는 셔플이다. 「유도된 곱(계수)」 모양이 아니라 평면의 *점* 정의에 대응하는 원통판이다.
   - `ref_dir` 은 프리즘이 스케치에서 진술한다(`construct.rs`): 온전한 원은 자기 한 정점 쪽(중심 → 정점, 반지름이 유리수면 그것으로 나눈 단위 길이 — `Ring2d::circle` 의 정점은 프레임 `(r, 0)` 이라 프레임 `+x̂`), 호의 벽은 프레임의 `x̂`. 열쇠가 진술 그대로라(아래 3) 한 원의 호들이 한 곡면으로 intern 되려면 같은 철자여야 한다 — 호가 자기 끝점이 아니라 프레임에서 seam 을 읽는 이유다.
2. **seam 은 모델 기하이고 영구 고정이다**(+`ref_dir` 방향, θ=0). 차트가 자기 배제점을 seam 위에 놓도록 맞춘다 — 차트가 seam 에 적응하지, 그 역은 없다.
3. **interning 은 보수적이다.** 원통의 열쇠(`SurfaceKey::Cylinder`)는 def **문자 동일**(`ref_dir` 포함) + motion 동일만 합친다. 같은 축·반지름·다른 `ref_dir` 을 합치면 seam 이 갈라지므로, 잘못 합칠 위험이 0 인 키다. 기하 동일성은 술어의 몫이다. 평면의 `flipped` 대응물은 필요 없다(문자-동일 키면 캐시 구성도 동일).
4. **`Vertex::OnSeam([원통, 캡])` 의 정의** = 「rim ∩ +`ref_dir` 방향 ray」. `ref_dir` 이 진실에 있으므로 유일점을 정확히 지시하고, 좌표 캐시는 load-bearing 이 아니다.
5. **seam 은 간선이 아니다.** 옆면은 위·아래 림 두 루프로 닫히고(「seam 엣지 방식」 B), `[s, s]` 는 어떤 곡면에도 담체 쌍이 아니다 — `push_edge` 가 한 원통 둘을 거절하고 validate 가 모든 self-pair 를 `EdgeCarrierMismatch` 로 잡는다.
6. **branch 방향.** 평면∩평면∩원통 = 최대 2점이고, `Vertex::Pierce { planes, cylinder, root: QuadRoot }` 변종이 받는다 — 변종이 자기 진실을 말한다(슬롯 재활용 금지). `QuadRoot` 는 `Lo | Hi | Double` **셋**이다: `Double` 은 접점(두 근이 일치해 점이 하나)이고, 저장된 값이 자기가 접점인지 말할 수 있어야 하므로 저장 변종이다. 재정렬 규칙은 `QuadRoot::canonical` **한 곳**에 살고, 접점 예외는 `flipped(Double) = Double` 이라는 원시에서 저절로 나온다(스왑은 중근을 자기 자신으로 보낸다).

## 이차곡면 판정 — 원통 우선

이차곡면 교차는 난이도가 갈린다. **평면∩이차곡면**(평면∩원통=타원, 평면∩구=원, 평면∩원뿔=원뿔곡선)은 닫힌 형식이라 쉬운 쪽이고, **이차곡면∩이차곡면**(일반적으로 4차 공간곡선)은 특수 케이스만 닫힌 형식이며 일반은 이미 SSI 에 가깝다. 그래서 이 영역은 「하나의 알고리즘」이 아니라 「어느 곡면쌍까지를 긋고 그 판정을 무엇으로 할지」의 선긋기다. 선은 이렇다: **원통부터**(평면∩원통), 그 뒤 구·원뿔. 일반 이차곡면쌍은 커버리지 밖으로 두고 `Rejected` 로 정직하게 거절한다 — 조용히 틀리는 길이 없으므로 「되는 데까지만 하고 나머지는 거절」이 안전하다.

**판정 방법.** 평면∩원통 교차점은 좌표가 `a+b√c`(a, b, c 유리수) 꼴이다(nacre-exact `QuadVal`).

- 점-대-평면 부호와 같은-판별식 비교는 `a+b√c` 의 부호(`QuadVal::sign`)로 닫힌다.
- 원통 면 위 사건들의 원형 순서는 서로 다른 절단면(서로 다른 판별식 u, v)에서 온 두 점의 비교라 `ℚ(√u,√v)` 의 4항 원소 부호가 필요하다. 이것은 재귀 탑으로 닫힌형이다(`P+√v·Q` 분해 → 상반 부호면 `sign(P)·sign(P²−vQ²)` — nacre-exact `quad.rs` 의 `biquad_sign`).
- 일반 대수적수 기계도, QI 의 pencil 분류도 원통에는 필요 없다.
- **차트 t 의 수치는 실현(캐시)·문서화 전용이다.** t 값 비교는 노름 라디칼(`|d|`, `|e₁|`)까지 끌고 와 3-라디칼이 된다. 원형 순서는 **(seam 반평면 부호 `w·e₂`, 같은 반평면 안 외적 부호 `(w₁×w₂)·d`)** 두 술어로 판정한다(e₁ = `ref_dir` 의 축-수직 성분, e₂ = d×e₁ — 둘 다 유리수라 위 두 부호 원시로 닫힌다). 차트의 배제점 = seam 이고, seam 모선 위의 점(θ=0)은 **이름**(`SeamIncident`)으로 답한다 — 순환 순서의 순위는 `SeamOrder::seam_first` 한 곳이 준다(seam 점이 첫째: 순환 순서는 어디서든 한 번 끊을 수 있다). 호의 범위(`arc_span`)와 노드의 원형 정렬(`circular_order`)이 같은 규칙을 읽는다.

### 참고 논문과 라이선스 규율

일반 이차곡면쌍에 갈 때 꺼낼 참고처다(원통에는 불필요).

- **QI 3부작 + 구현 논문**(pencil 분류로 교차 타입을 대수적으로 exact 결정). Dupont·Lazard·Lazard·Petitjean, "Near-optimal parameterization of the intersection of quadrics" — SoCG 2003, JSC 2008 저널 3부작(Part I 생성 알고리즘 43(3):168–191, Part II pencil 분류 43(3):192–215, Part III 특이 교차 43(3):216–232). C++ 구현 논문: Lazard·Peñaranda·Petitjean, "Intersecting Quadrics: An Efficient and Exact Implementation", SoCG 2004 / Comp. Geom. 35(1–2):74–99, 2006.
- **CGAL 3D Spherical Kernel**(2차 대수적수 exact — 다른 계보). de Castro·Cazals·Loriot·Teillaud, Comp. Geom. 42(6–7):536–550, 2009. 간접 술어(Attene, 정의 기반)와 다른 계보라 「확장」이 아니라 「대안 접근」으로만 참고한다.

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

`nacre-judge`(Certified Indirect Predicates)는 회전이 들어왔을 때의 부호 판정 층이다.

### 왜 필요한가

축정렬 판정은 좌표가 유리수라 exact 다. **회전이 들어오면 좌표가 무리수가 된다** — 유리수 각도라도 cos·sin 은 무리수이고(Niven), 임의 각도는 초월수다. exact 산술로 표현할 수 없으므로 회전 좌표는 근사일 수밖에 없다. 모델링 변환·각도 스케치가 회전 좌표를 만들 때 **판정만은 조용히 틀리지 않게** 지키는 것이 CIP 다.

**핵심 명제.** indirect predicates(Attene 2020)는 선형 요소(선·평면)의 교차, 즉 다항식에만 성립한다. 회전은 다항식이 아니므로 간접 술어를 회전에 적용하는 방법은 없다(수학적 한계). 할 수 있는 일은 하나다: **회전 좌표는 근사하고, 그 위의 판정을 exact 로 유지한다.** CIP 는 그 판정에 **입력 근사(tol)까지 반영**하는 층이다. 세 계보의 하이브리드다 — Attene 2020 간접 술어(implicit point 구조) · Shewchuk 1997/CGAL 필터(f64 필터 → 고정밀 상승) · Guibas 1989 epsilon-geometry(값 + tol, sound 판정 아니면 기권). 유리수부는 완전 exact, 초월(회전)부는 sound 오차 한계로 부호를 **인증하거나 정직하게 기권**한다. **포기한 것은 초월수까지의 이론적 완비성이지 건전성이 아니다** — 틀린 부호는 결코 내지 않는다.

**대체가 아니라 흡수다.** 점을 「정의 + 누적 tol」로 들고, **tol = 0 이면 indirect predicates 그대로(exact)**, **tol > 0 이면 필터 + 고정밀 상승**이다. 술어의 계산 구조(3-평면 정의, 행렬식 부호, 좌표 무독)는 그대로 쓰고, 「평면 계수가 정확하다」는 가정만 「계수에 tol 이 있다」로 넓힌다. 회전된 평면은 계수에 tol 이 붙은 평면일 뿐이다.

### 점 = 정의 + 방향별 tol — `WitnessPoint`

```rust
pub struct WitnessPoint {
    pub base: [Rat; 3],
    pub chain: HpRc<[MoveNode]>,
    pub realized: [Bounded; 3],
    hp: HpCell, // memoized high-precision realization
}
```

- `base` + `chain` 이 정확한 **정의**이고 결코 잃지 않는다. `realized` 는 f64 실현(캐시)이며 축마다 **값과 그 오차 한계를 한 원자(`Bounded`)에** 묶는다 — 둘이 서로 어긋날 수 없다. `coord()` / `tol()` 이 벡터로 읽어 준다.
- 평면 교차점(tol = 0)과 회전 점(tol > 0)이 **같은 틀**이다. 특수 경로가 없다.
- **tol 은 xyz 방향별 벡터다.** 스칼라(구) tol 은 정보를 버린다 — z축 둘레 회전이면 z 방향 오차가 0 인데 스칼라는 모든 방향에 실어 「z=5 평면 위인가?」를 불필요하게 애매하게 만들고 상승을 유발한다. 행렬식에서 각 방향의 오차는 **서로 다른 계수로 증폭**되므로(2D orient 에서 `a` 의 x-tol 은 `(by−cy)` 와, y-tol 은 `(bx−cx)` 와 곱해진다) 하나로 뭉치면 정확히 증폭할 수 없다. CGAL Lazy_kernel 이 좌표를 구간으로 드는 이유와 같다. nacre 는 범용 구간 산술 대신 **방향별 tol + 미리 유도한 오차 한계 공식**을 쓴다 — 술어가 소수·고정이라 유도가 가능하다(CGAL 은 술어가 수백 개라 구간을 택했다).
- `hp_coord` 는 정의에서 임의 정밀도로 사슬을 실현한다. 같은 정의를 가진 두 점은 동일하게 실현된다(**경로 무관** — 건전성 논증의 뿌리).
- **정의 동등성.** `WitnessPoint` 의 `PartialEq` 는 `base` 와 `chain` 만 본다. `realized`/`hp` 는 캐시이고 실현은 정의의 순수 함수이므로, 정의가 같은 두 점의 캐시는 정직하게 다를 수 없다. 이 동등성의 소비자는 `shared_base` 의 노드 통째 사슬 비교다 — 두 모션이 *한* 모션이라고 선언하는 것은 정의에 대한 질문이지 캐시에 대한 질문이 아니다. `FrameThrough` 노드가 점을 안에 싣고 있으므로 이 동등성이 노드 비교까지 닿는다.
- `chain` 이 `HpRc` 인 이유: clone 을 refcount 증가로 만들면서 `parallel` 빌드의 스레드 경계를 지킨다(`parallel` 에서 `Arc`, 아니면 `Rc`). 불리언은 모든 worker 에게 같은 `&[WorkingPlane]` 를 건네므로 `WitnessPoint` 는 `Sync` 여야 한다. 고정밀 실현 캐시 `HpCell` 도 같은 선택이다(`Arc<OnceLock>` / `Rc<OnceCell>`) — 두 worker 가 한 셀을 다투어 채워도 같은 값을 계산하므로 답은 누가 이겼는지에 의존하지 않는다.

**생성자는 경우를 이름으로 구분한다.**

- `at(base)` — 유리수 base 의 f64 반올림을 고정밀(120비트)에서 **측정**해 tol 의 씨앗으로 삼는다(f64 로 표현 가능하면 정확히 0). 이후 사슬이 돌리지 않는 축은 정확히 이 값을 유지한다.
- `at_nearest(base)` — 같은 수를 **계약에서** 얻는다: 좌표가 표현 가능하면 정확히 0, 아니면 `|x|·2⁻⁵³`(`Rat::to_f64` 가 「최근접, ties-to-even」을 문서로 약속하므로 `|r − to_f64(r)| ≤ ½ ulp`). `at` 의 측정은 점당 BigFloat 연산 아홉 번이고 이것은 공짜다. base 가 **정의**인 자리(평면 자신의 점, 담체에서 푼 솔리드 정점)에서 `nacre-ops` 가 쓰는 철자다. 한계의 느슨함은 확정적 필터 답을 「상승」으로 바꿀 수만 있고 반대는 없다. **유리수에서 온 점이 술어에 닿는 길은 이것이다** — 반올림된 f64 를 다시 `Rat` 로 들어올리는 것이 아니다(그것은 다른 점을, tol 0 으로 이름 붙인다).
- `at_with_tol(base, tol)` — 이미 tol 을 지닌 뿌리(불리언 seam 정점)의 씨앗. 사슬이 그 tol 을 싣고 간다(`|R|·기존 tol`).
- 「좌표가 f64 로 정확하다」를 전제로 받는 생성자는 없다 — 그 전제는 어떤 호출자도 검사할 수 없고, 반올림된 캐시를 건네면 다른 점을 이름 붙인다. 표현 가능한 점에서는 `at_nearest` 가 같은 tol 0 을 말하고 그 밖에서도 정직하다. (답이 0 인 것을 120비트로 측정하는 것은 비용이다 — 회전 불리언 시간의 절반이 그것이었던 구성이 있다.)

### 모션 사슬 — `MoveNode`

```rust
pub enum MoveNode {
    Rotate { axis: Axis, angle: Angle, pivot: [Rat; 3] },
    Translate { offset: [Rat; 3] },
    Mirror { axis: Axis, offset: Rat },
    Frame { frame: nacre_exact::PlaneFrame },
    FrameWide(WideFrame),
    FrameThrough(Box<FrameThrough>),
}
```

모션은 교환되지 않으므로 **사슬의 순서가 곧 정의다.**

- **트리는 tol 계산만이 아니라 정의의 합성을 한다.** 이동·반사도 정의의 일부다 — 빠지면 `T(7/11)` 한 벽과 `T(18/11)` 한 벽이 같은 평면임을 말할 방법이 없고(1 ULP 갈라짐 → 한 부품이 두 몸통), `R` 다음 `T` 는 표현조차 못 한다. 교차는 담지 않는다(그것은 `Vertex` 의 몫이다).
- **옮겨 적기 — 진술이 정확히 옮겨지면 노드를 만들지 않는다.** 솔리드의 모든 면의 진술(점 셋·원통 정의)이 모션 전체로 `i128` 안에서 정확히 옮겨지면(`transport_points`/`transport_cylinder`) 모션을 진술에 옮겨 적고 노드를 만들지 않는다 — 두 벽은 그 자리의 세계 진술로 한 이름이 된다. 하나라도 못 옮기면(`Through` 평면은 핸들 진술이라 옮길 수 없다, 넘침) 솔리드 전체가 기록한다 — 면마다 가르면 기록한 면과 옮긴 면이 만나는 모서리가 혼합 프레임 정점이 된다. 이력이 있는 면은 언제나 기록한다(그 진술은 사슬 이전의 평면이라 그대로 남으므로). 판단은 진실만 읽는다 — f64 캐시가 정확히 옮겨지는지는 묻지 않는다(「가지 말 것」). 두 상태뿐이라(`Carry::{Full, None}`) 한 면 안에서 모션의 일부는 진술에, 나머지는 노드에 있는 일은 없다.
- **`Mirror` 만 improper 다**(`det = −1`). 공유 모션을 행렬식에서 상쇄하는 판정은 먼저 모션 전 데이터를 같은 손잡이로 가져와야 한다 — `chain_parity` 와 `shared_base`. 홀수 반사 사슬은 모든 입력의 행렬식을 **일률적으로** 뒤집으므로, 보정이 없으면 지름길은 보수적 miss 가 아니라 확신에 찬 오답을 낸다. 평면은 보정된 점에서 **유도**해야 하고 따로 보정하면 안 된다(외적은 pseudovector 라 전역 부호가 갈린다). 프레임 노드 셋은 proper 다(`v̂ = ŵ × û` 구성).
- **`Frame`** — 평면 자신의 프레임으로의 기저 변환. 모든 입력이 정확한 유리수이고 **유일한 무리수 걸음은 길이로 나누는 것**이다: `u_raw = ẑ × n`(또는 `ŷ × n`)은 이미 평면 안에 있어 사영이 필요 없고, 실현 전체가 `1/√(유리수)` 스칼라 둘이다 — 그래서 판정할 수 있다. `v̂` 도 정확하다: `n ⊥ u_raw` 라 `|v_raw|² = |n|²·|u_raw|²` 이 정확하고 역제곱근 하나로 실현된다(f64 에서 `ŵ × û` 로 실현하면 상쇄되지 않는 반올림 둘이 들어 정규직교가 깨진다).
- **`FrameWide`** — 정확 데이터가 `Rat` 에 안 들어가는 평면(`Wide` 이름, 또는 제곱 길이가 `i128` 을 넘는 narrow)의 `Frame`. 임의정밀 정수라 **넘칠 수 없고** 부분형이 없다. 이 변종 때문에 `MoveNode` 는 `Copy` 가 아니다(사슬은 `HpRc` 로 공유되므로 뜨거운 경로가 노드를 복사하지 않는다).
- **`FrameThrough`** — **이름이 아예 없는** 평면(혼합-프레임 datum: 세 정점을 유리수로 푸는 한 프레임이 없어 어느 정확 그릇에도 안 담긴다)의 `Frame`. 노드가 평면의 **정의점 셋을 싣고**, 정준 기저(원점 = 세계 원점에서 내린 수선의 발, 임의 축 규약)를 **구간으로 유도**한다 — 실현이 도는 정밀도에서, 오차를 싣고.

```rust
pub struct FrameThrough {
    pub points: [JudgedPoint; 3],
    pub vertical: bool,
    pub flip: bool,
}

pub enum JudgedPoint {
    Pure(WitnessPoint),
    Meet(Box<[[WitnessPoint; 3]; 3]>),
}
```

  - `vertical` 은 임의-축 분기(`ẑ×n` 대신 `ŷ×n`)이고 **고정 128비트 실현에서 한 번 판정해 저장한다.** 실현 정밀도마다 판정하면 f64 캐시와 상승 사이에서 기저가 바뀔 수 있고, 정밀도에 따라 기저가 바뀌는 프레임은 프레임이 아니다. 이름 있는 평면은 분기가 `n₀ = n₁ = 0` 을 정확히 읽지만 여기서는 정확한 0 을 증명할 수 없으므로 「증명 가능하게 쓸 수 있는 쪽」이다(`|ẑ×n|²` 이 0 에서 떨어져 있음). 어느 분기도 증명 못 하는 평면은 생산자가 **이름으로** 거절한다 — `FrameThrough::of` 의 `None` 은 퇴화의 증명이 아니라 건강함을 증명하지 못한 것이므로, 호출자는 퇴화라고 주장해서는 안 된다(`CollinearVertices` 는 여기서 거짓말이다).
  - `flip` 은 모든 프레임 노드가 드는 같은 측정된 방향이다: 네 계수를 부정한 뒤 기저를 유도한다. 정준 원점은 부호·배율 불변이고 축은 방향만 읽으므로 전역 부정 하나가 자유도의 전부다. 분기 판정은 제곱 길이만 보므로 `flip` 불변이다.
  - **`JudgedPoint`** 는 정의점이 **어떻게 진술됐는지**다. `Pure` = 자기 프레임에서 정확한 점(유리수 base + 사슬). `Meet` = **세 담체 평면의 교점 — 어디에도 좌표가 없는 점.** 각 담체는 자기 증인 삼각형으로 기술되고, 점은 담체들의 Cramer 계(`cramer_hp`)가 낸 동차 `[Dvec : D]` 로 실현된다. `Pure` 는 `[p : 1]` 이라는 특수 경우다.
  - 동차점 셋을 **사영 join**(`plane_hp_through` — 3×4 의 네 3×3 소행렬식, 부호 교대)이 평면으로 잇는다. meet(세 평면 → 점, Cramer)의 **쌍대**라 새 산술이 없다. **여기서는 아무것도 나누지 않는다** — `Dvec/D` 로 아핀 좌표를 만들면 무리수를 제조해 그 반올림을 하류 전체에 굽는다. 차수는 9 다(아핀 외적 경로는 15 — 구간 폭과 요구 정밀도가 그만큼 준다).
  - **배율의 부호는 값 안에서 없앤다.** join 은 행에 대해 다중선형이라 결과가 참 평면의 `D0·D1·D2` 배다. 양의 배는 무해하지만(평면에 묻는 것은 전부 부호 질문) 음의 배는 방향을 조용히 뒤집고, 방향은 `frame_sign`·외향 법선·딱지 프레임 전체의 토대다. 함수가 자기 결과를 부정한다. 어느 `D_i` 의 부호가 미결이면 `None` — 그 세 평면이 한 점에서 만나는지조차 모르므로 호출자가 상승한다.
  - meet 가 나뉘는 자리는 `anchor_coord` 한 곳이다(캐시가 닻을 내릴 아핀 `[f64; 3]`). 제수는 구성상 안전하다 — `FrameThrough::of` 가 rung 유도를 증명했고 거기에 이 점의 `D` 부호가 포함된다.

**실현은 나눠도 된다. 금지되는 것은 술어의 부호 질문에서의 나눗셈이다.** `HpBounded::div` 는 구간 제수로 나누되 먼저 제수를 0 에서 떼어 놓고(`None` = 제수 구간이 0 에 닿을 수 있음 — 호출자는 상승하거나 이름으로 거절한다), 자기가 만든 반경을 돌려준다. 소비자는 프레임 원점 `(−d·n)/(n·n)` 이다. `HpBounded::inv_sqrt` 는 반경을 지닌 입력의 `1/√x` 다(오차 증폭은 도함수에서: `r / (2·L^(3/2))`; `None` = `x` 가 0 에 닿거나 양이 아님 — 퇴화 법선에는 방향이 없고 호출자가 이름으로 말한다). `HpBounded::div_exact` 는 **정확한** 제수만 받는다(wide 프레임의 원점 `num / den`). 부호 질문에서 몫은 술어가 믿어야 하는 무리수를 제조하므로 거기서는 나누지 않는다.

**판정 표의 계약은 «평면 위 정확한 세 점»이다.** 이름 없는 평면의 증인은 판정 프레임의 probe `(0,0,0)·(1,0,0)·(0,1,0)` 이다 — 프레임 사슬을 실은 정확한 정의이고, 정의상 그 평면 위에 있다. 그래서 `Through` 평면 전용 판정 기계는 없다: 구간-계수 경로 없이 기존 세 술어가 그대로 답한다.

### tol 의 전파

- **tol 의 두 축.** 자체 tol = **이산화 오차**(값을 유한 정밀도로 표현할 때; f64 ≈ 1e-16 × 크기, `p` 비트 ≈ 2⁻ᵖ × 크기) + **연산 오차**(유리수 덧셈 = 0, f64 덧셈·곱셈·회전 = 발생). 누적 tol = 앞 tol + 자체 tol.
- **오차 한계는 최악(선형 합)이다.** 분리 계산은 삼각부등식 `|a+b| ≤ |a|+|b|` 로 정당하다(최악 상한이므로 실제보다 작아지지 않는다). 확률 전파(RSS)는 실제 오차에 가깝고 100배 작지만 **보장이 아니다** — 드물게 부호를 뒤집어 조용히 틀린다. 커짐은 다른 방법으로 관리한다(유리수 구간은 tol 0, 애매하면 필요한 만큼 정밀도를 올린다).
- **회전의 tol 갱신**(`WitnessPoint::rotate_about`; 회전 평면의 두 좌표 `i, j`, pivot 기준 `u, v`, 실현된 `c, s`):

```
새 tol_i = |c|·tol_i + |s|·tol_j            [기존 tol 의 방향 회전 — |R| × 기존 tol]
         + |u|·dc + |v|·ds + ε·(|u|+|v|)    [이 회전의 실현 오차 + 곱·결합의 f64 반올림]
         + piv                               [원점이 아닌 pivot 의 산술]
새 tol_j = |s|·tol_i + |c|·tol_j + |u|·ds + |v|·dc + ε·(|u|+|v|) + piv
```

  - `|R|·(기존 tol)` 은 **필수**다 — x 방향 tol 이 z축 90° 회전 뒤 y 로 옮겨 가는데 반영하지 않으면 상한이 깨진다.
  - 새 오차는 **좌표 혼합형**이지 접선형(`da × 모멘트암`)이 아니다. cos/sin 을 독립으로 반올림하면 순수 회전이 아니라 **방사 성분**이 생기고(`err_x = |x·δc − y·δs|`), 접선형은 좌표 하나가 0 근처일 때 그것을 과소평가한다(근거: 랜덤 10000 중 220회 상한 붕괴; 좌표 혼합형은 0회, 최악 약 3.5배 보수).
  - `dc`/`ds` 는 상수가 아니라 **이 각도의 실현을 측정한 값**이다(`Angle::realization_error_of(c, s)` — 위에서 실제로 쓴 그 쌍을 잰다). **축마다 다르게 청구한다**: `i` 는 `u·c − v·s`, `j` 는 `u·s + v·c` 라 `u` 가 한쪽에선 cos 의 오차를, 다른 쪽에선 sin 의 오차를 만난다(356.65° 회전에서 sin 이 cos 보다 8배 멀고, 한 축이 다른 축의 2.4배를 요구했다).
  - **실현이 정확하면 두 항 모두 0 이다.** `cos`/`sin` 이 `0`/`±1` 이면 곱과 결합이 *정확*하므로, 원점 pivot 의 사분각 회전은 tol 0 을 유지하고(`quadrantal_origin_chain_is_tol_zero`) 축정렬 모델은 exact 술어 경로에 남는다. 무조건 청구되는 산술 항은 이 불변을 조용히 깬다.
  - `piv` 는 pivot 의 `Rat → f64` 실현을 **측정**하고 round-to-nearest 걸음(차 둘, 합 하나)을 **센다**. 이분(dyadic) pivot 에서 실현 항은 정확히 0 이고, 원점 pivot 에서는 항 전체가 0 이다.
  - **각도 불확실성 항**(접선: `da_각도 × 회전축까지의 수직거리`)은 유리수 각에서 0 이다. 유도 치수(구속 솔버·비유리수 각)에서 살아나는 항이므로 설계에서 지우지 않고 0 으로 둔다.
- **회전 기여는 거리에 곱해진다.** 회전은 *방향*을 정하므로, 앞선 회전일수록 뒤 이동이 많아 기여가 크다. 이산화·직선 tol 은 적용점과 무관해 점에 값으로 캐시하지만, 회전 기여는 적용점에 의존한다 — 그래서 점은 자기 사슬(조상 모션)을 들고 다닌다.
- **유리수 구간은 묶는다.** 기준은 「직선/회전」이 아니라 **「유리수/무리수」**다. 연속된 유리수 이동을 먼저 합산하면 그 구간의 연산 오차는 0 이고 이산화는 회전과 만나는 지점에서 한 번뿐이다(`(a+b)+c` 는 이산화 2회, `a+(b+c)` 는 1회). 이동·사분각 회전·축 거울의 사슬은 노드가 태어날 때 한 유리수 사상으로 접힌다(`Model::chain_point_rat` 등, 「캐시」의 사슬 접기). 같은 축·같은 피벗의 연속 유리수-각도 회전도 각도를 유리수로 합산하면 누적 각을 한 번만 실현할 수 있지만(`Angle` 의 덧셈이 mod 360 으로 정확하다) 그 합산은 짓지 않았다 — 인접 회전쌍의 인구가 없다(`todo.md`). **묶기가 중요한 이유**: 묶지 않은 증분 실현은 전파 `|R|` 의 행합 `|cos|+|sin| ≥ 1` 을 매 걸음 곱해 tol 한계가 지수로 폭발한다(30걸음에 실제 오차의 약 20만 배; 실제 오차는 `R` 이 노름을 보존하므로 평평하다). 폭발해도 sound 한 최악 보장이라 조용히 틀리지 않고 상승/거절로만 간다. 회전 결과는 무리수라 회전 경계의 실현 오차는 피할 수 없다 — 피할 수 있는 것만 피한다.
- **평면 계수의 tol 은 점 tol 의 따름정리다.** 평면은 세 점으로 정의되고, 회전되면 그 점들이 방향별 tol 을 가지며, 계수는 그 점들의 뺄셈·외적이라 점 tol 이 계수 tol 로 전파된다. 판정 경로에 sqrt 섭동은 없다 — 평면은 정규화하지 않은 원시 계수를 들고 판정은 그것을 읽으며, 정규화된 법선은 크기 소비자에게만 간다.

### 고정밀 층 — `astro-float`

고정밀 층은 **`astro-float`**(순수 Rust 임의정밀)다. double-double(`twofloat`)은 탈락이다: π 상수는 정확하나 삼각함수가 영점 근처에서 f64 수준(cos 오차 약 1.8e-16)이고, 부호 판정이 일어나는 near-degenerate 가 곧 영점 근처다. `astro-float` 는 160비트에서 cos/sin 오차 약 1e-58. **정밀도가 dial 가능**하므로 double-double 의 약 1e-32 천장이 없다 — 「더 넓은 고정 폭으로 상승 or 정직 거절」의 딜레마가 「정밀도 dial」로 단순해진다.

**이 층은 상승 경로 전용이 아니다.** `Angle::cos_sin_f64` 는 libm 을 부르지 않고 128비트 실현을 f64 로 **정확 반올림**해 돌려준다(`round_to_f64`). `f64::cos` 에는 정확도 계약이 없고 각도의 함수조차 아니다(`sin 27°` 가 debug/release 사이, 그리고 한 release 빌드의 두 호출처 사이에서 1 ulp 달랐다 — LLVM 의 컴파일 타임 상수 평가가 런타임 라이브러리와 다르다). 정확 반올림된 값은 **유일**하다: 플랫폼·프로파일·호출처에 무관하게 같은 비트이고, f64 캐시는 진실의 반올림 사본이 된다. 90°-계열 분기는 최적화가 아니라 **종료 조건**이다(`cos 90°` 는 정확히 0 이라 구간이 0 에 걸쳐 어떤 깊이에서도 양 끝이 같은 f64 로 반올림되지 않는다). 실현은 각도당 메모된다(한 번 약 32µs, libm 은 22ns; 판정의 `trial_bound` 와 같은 `TRIG` 항목을 쓰므로 판정까지 가는 모델은 일이 오히려 준다 — 회전 fold 80: 731ms → 715ms).

이 문서의 「고정밀 상승」은 전부 이 층을 가리킨다.

### 판정 = 필터 + 계산된 상승 + 실현 캐시

1. 점들의 tol 로 이번 판정의 오차 한계를 계산한다(변 길이로 증폭). `|행렬식| > 오차 한계` 면 f64 로 확정한다(대부분).
2. 애매하면 상승한다. 관련 점이 tol = 0 뿐이면 **exact 경로**(indirect predicates 그대로), 회전 점이 끼면 **임의 정밀도**로 간다.
3. **다음 정밀도는 계산한다 — 배증하지 않는다.** 간격은 `C · 2⁻ᵖʳᵉᶜ` 꼴이므로 `log₂(간격/목표)` 가 곧 모자란 비트 수다. 한 번에 점프한다(워드 단위 올림).
4. 고정밀로 계산한 **점의 값**을 캐시한다(판정 결과가 아니라 값 — 여러 판정이 재사용한다). 판정에 필요한 점만 정밀화하고 중간 경유점은 정의로만 남는다.

**필터는 동적 하나다.** tol 없는 입력을 가정한 정적 상수 필터(`ε_D ≈ 5u`·`ε_M ≈ 19u` 류)를 따로 두지 않는다. 회전 tol 은 경로마다 편차가 커서 정적 상수 하나로 잡으면 느슨해 무용하거나(과대) 위험하다(과소). 각도를 쓰는 스케치도 회전이라 무회전 케이스는 소수이고, 동적 필터는 tol = 0 인 점을 자동으로 정적급으로 타이트하게 다루므로 정적 분기는 순수 오버헤드다.

**세 술어.** 평면 표의 인덱스로 묻는다(`Judge::orient3d` · `Judge::cmp_coord` · `Judge::plane_pair_dir_sign`).

- **직접 orient3d** — 명시적 네 점. 오차 한계는 `frame3::det3_bound`: 행렬식은 간선 성분의 부호 있는 삼중곱 여섯 개이고 각 성분 `(a−d)[k]` 는 tol `tol_a[k] + tol_d[k]` 를 든다. 여섯 곱의 구간 반경(`Π(|v|+τ) − Π|v|`)의 합 + 행렬식 자체 산술의 f64 반올림 항. `orient3d_judge` 가 필터 → 상승 → 정규화된 간격 판정으로 소비한다.
- **간접 orient3d** — 3-평면 implicit point 대 평면(`indirect_orient3d_judge`). 세 평면의 구간 계수 → Cramer(`cramer_iv`) → `orient3d_from_cramer`.
- **cmp_coord** — 두 implicit point 의 한 축 좌표 비교(`indirect_cmp_coord_judge`).
- 방향 질문(`plane_pair_dir_sign` → `dir_sign_judge`)은 길이 한계를 모델 크기로 각도로 바꿔 쓴다.

**공유 모션은 상쇄한다.** 모션은 이 술어들이 취하는 모든 행렬식을 *자기 행렬식 배까지* 보존한다. 한 판정의 모든 입력이 같은 사슬을 들면(`Witness::chain_id` — 구조적으로 동일한 사슬에만 같은 값) 답은 모션 전 데이터에서의 답이다 — 정확히, tol 없이. 그 판정은 tol 경로를 통째로 떠난다(`Witness::base_tri`, 손잡이 보정은 위 `Mirror` 항).

**wide 이름의 정수 가지.** 세 술어는 평면 이름의 **`BigInt` 정수**로 정확히 답하는 가지를 갖는다(scalar `int_plane_side`·`int_cmp_coord`·`int_dir_sign` — f64 쌍둥이와 같은 규약; `Expansion` 은 지수 상한이 약 2¹⁰²³ 이라 wide[약 2²²⁹¹]를 못 담는다). 게이트 = 전 평면이 이름을 가짐 ∧ 하나 이상 wide ∧ (전부 무이동 또는 전부 한 사슬; `cmp` 는 무이동만). narrow 전용 질문은 기존 경로 그대로다. 이름은 구성 시 `name_stored_ints` 가 σ(정준 ↔ 저장 방향) × `frame_sign` 으로 접어 **`base_coeffs` 와 같은 방향**을 들고, 홀수 미러 사슬은 계수에서 `parity·C`(x 유지·y/z/d 부정 — 외적이 pseudovector)로 보정한다.

**평면 증인은 정의를 소유하고, 「회전됨」은 별개 사실이다.** `Witness::tri_pt3() -> &[WitnessPoint; 3]`(캐시, 항상 존재) + `Witness::is_rotated() -> bool`(술어 경로 선택). 정의의 *존재*가 「회전됨」까지 답하게 하면 정의를 미리 저장할 수 없고(저장하면 모든 평면이 회전으로 읽힌다) 판정마다 정의를 재구성하게 된다(핀 25개 불리언에서 108만 회, 시간의 77%). **플래그는 유도하지 않는다**: `chain.is_empty()` 도 `tol == 0` 도 그것과 등가가 아니다 — 회전된 솔리드의 면이 사슬 없는 점으로 평면을 증언할 수 있고, 90°-계열 회전은 tol 이 정확히 0 이다. 두 값은 `collect_planes` 한 자리에서만 함께 설정한다. 각 평면의 구간 계수는 `Judge` 가 첫 사용 때 짓고 빌려준다 — `PlaneWitness` 는 소비자가 구현하는 트레이트라, 구간을 *생산*할 의무(반경이 실제로 계수를 가둔다는 건전성 계약)가 이 크레이트 밖으로 나가면 안 된다.

### 정밀도는 모델이 정한다

판정의 오차 반경은 `C·2⁻ᵖʳᵉᶜ` 이고 `C` 는 **모델의 성질**(회전 이력 1회당 약 1비트, 좌표 크기)이지 정밀도의 함수가 아니다. 고정 정밀도는 *모델의 회전 이력이 얼마나 길 수 있는지를 조용히 결정*한다(고정 256비트에서 245회 회전한 솔리드는 빌드에 실패한다). 그래서 `C` 를 연산마다 한 번 읽어 필요 비트를 계산한다(`judge_precision`).

```rust
pub struct Standard {
    pub prec: usize,       // the precision escalation realizes at, chosen per model
    pub coincidence: Mag,  // proved-closer-than-this means one thing
    pub scale: Mag,        // model size
    pub cap: usize,        // most bits an escalation may ask for
}
```

- 문턱은 **모델 단위의 길이**이지 비트 수가 아니다. 「256비트」는 회전 1회 모델에선 1e-76, 300회 모델에선 1e+15 를 뜻한다 — 비트로 노출하면 안 된다.
- **`output_precision` = `scale · 2⁻⁵²`** — f64 좌표가 이 모델에서 볼 수 있는 가장 가는 눈금. **`coincidence` = `output_precision · 2⁻¹²⁸`** — 워드 둘 아래.
- **`coincidence` 는 tolerance 가 아니다.** 전역 tol 은 *「이보다 가까우면 붙여라」*(모르는 채 뭉갬)이고, 이것은 *「이보다 가깝다고 **증명되어야** 일치」* 로 방향이 반대다. 무지 위에서 병합되는 것은 없다.
- **`JUDGE_PREC_CAP = 4096`**(nacre-ops `planes/standard.rs`)은 정확성이 아니라 **비용** 한계이고 값의 근거는 측정된 비용이다(`C` 는 회전당 1비트 자란다; 회전 3200회·3456비트에서도 부피는 정확하고 90초; 4096비트 ≈ 약 4000회). 초과는 `RejectReason::PrecisionBudget { needed, cap }` 로 **평면 표가 서자마자** 거절한다 — 모델에는 아무 잘못이 없다.
- **모델의 깊이와 한 판정의 난이도는 별개 예산이다.** 판정이 오를 수 있는 여유 **`CLIMB_HEADROOM = 128`** 은 모델 정밀도에 대한 **상대값**이다(공유하면 회전이 많은 모델에서 여유가 0 이 되어 같은 얇은 증인이 회전 여부에 따라 판정되거나 포기된다). 물리적 의미는 **증인의 얇기 한계** `log₂(1/여인수)` — 모델 크기보다 `2¹²⁸` 배 퇴화한 증인까지 끝까지 판정한다. 측정값이 아니라 워드 둘인 이유는 오차의 비대칭이다: 너무 작으면 답이 있던 판정을 포기하고, 너무 크면 비트만 쓴다(회전 코퍼스와 100·800회 회전 모델에서 모델 정밀도보다 1비트라도 더 요구한 판정은 없다).
- 설정으로 노출할 후보는 일치 정밀도 하나뿐이고, 정밀도·상한·여유는 전부 거기서 유도된다.

### 판정 결과와 보고

```rust
pub enum Decision {
    Sign(Orient),
    Coincident { within: Mag },
    Exhausted { at: usize, within: Option<Mag> },
    Degenerate,
}
```

`Orient::Zero` 하나로는 0 이 증명된 것인지 가정된 것인지 말할 수 없다 — 그것이 커널이 조용히 추측하게 되는 길이다. 네 갈래가 정직한 분할이다.

- `Sign` — **증명된** 부호. `Zero` 는 증명할 수 있는 경로(exact 술어, 상쇄된 회전)에서만 나온다.
- `Coincident { within }` — 둘이 `within` 이내임이 보였고 그것이 일치 한계 이하다. 하나로 취급하며 `within` 이 그 근거다.
- `Exhausted { at, within }` — 상한 `at` 비트에서도 0 에 걸쳐 있고, 그것이 뜻할 수 있는 간격이 일치 한계보다 **크다** — 일치라 부르면 추측이다. 비트가 더 있으면 가를 수 있다. **자기가 세운 상한을 싣는다**: 참값은 `±within` 안이다. `1e-30` 에서 멈춘 판정과 `1e-3` 에서 멈춘 판정은 같은 변종이되 같은 소식이 아니다. `within` 이 `None` 이면 여인수가 상한에서도 안 풀려 인용할 거리 자체가 없다.
- `Degenerate` — 행렬식을 거리로 바꾸는 여인수를 0 에서 뗄 수 없다(퇴화한 증인 삼각형, 만나는 점이 없는 세 평면). **`Exhausted` 와 원인이 다르고 그 차이가 중요하다**: 비트로 해결되지 않는다(8192비트 더 깊어도 행렬식은 비트-정확히 0).
- 상승 중에는 셋을 구별한다: 간격이 아직 넓다 / 여인수가 0 은 아닌데 자기 반경에서 안 떨어졌다(`Unresolved` — **이것도 오른다**; 아래로 뭉개면 충분히 들여다보지 않은 둘을 조용히 병합한다) / 여인수가 정확히 0 이다(`Vanished` → `Degenerate`).
- `Decision::orient()` 는 결론 없는 결과를 전부 `Zero` 로 접는다. 증명된 일치와 고갈된 판정은 배열 엔진에게 *같은 지시*(「같다고 취급하라」)이고, 다른 것은 나중에 무엇을 말할 수 있느냐뿐이다. 그 차이를 제어 흐름 밖에 두므로 보고 채널이 기하 결정을 하나도 건드리지 않는다.

**묻지 않고, 근거를 붙여 보고한다.** 「애매하면 사용자에게 확인」은 성립하지 않는다: (a) 애매한 판정은 연산당 수백 건이고, (b) 그 애매함의 크기는 2⁻²²⁹(좌표 크기 1 인 모델에서 소수점 69자리)라 사람이 판단할 대상이 아니며 — 사람이 아는 것은 *의도*이고 의도는 이미 스크립트에 있다(같은 회전 공유 → 상쇄로 정확히 답함, `Surface` 핸들 공유 → 핸들 비교) —, (c) 「떨어짐」이라는 답은 만들 수조차 없다(출력 좌표가 f64 라 1e-16 아래 두께는 STEP·메시·화면 어디에도 안 나온다). 선택지가 하나뿐인 질문은 질문이 아니다.

- 판정은 일치 정밀도보다 가깝다고 *증명*되면 일치로 처리하고, 그 근거(무엇이 무엇과 몇 이내인지)를 결과와 함께 보고한다 — `boolean_with_report`(병행 진입점; 보고를 원하는 호출자는 서른에 하나 꼴이라 `boolean` 의 시그니처를 넓히지 않는다).
- **수집 채널은 전역도 thread-local 도 아니다.** 연산 단위 컨텍스트 `Judge { planes, standard, notes, .. }` 가 소유한다. 증인 표는 순수한 *설명*으로 남고 술어는 그 컨텍스트의 메서드다 — **연산의 성질은 연산이 갖는다**(행마다 복제하면 자리표시자·2단계 생성·「찍는 걸 잊었을 때」 가드가 딸려온다). 순수 수치 크레이트에 가변 전역이 생기지 않는다. 수집은 진단 전용이다 — 어떤 부호·병합·좌표도 그것을 읽지 않는다.
- 보고의 **첫 항목은 평면 클래스 병합**이다(정점 하나 계산되기 전에 「어떤 평면이 존재하는가」를 바꾸므로).
- **증명하지 못한 판정은 0 이 아니라 거절이다.** `RejectReason::JudgeExhausted`(한 판정이 비트 부족 — 모델 전체가 깊은 `PrecisionBudget` 의 형제) · `RejectReason::DegenerateWitness`(여인수가 사라짐)로 원인별로 나간다. **원인이 증상보다 앞선다** — 하류의 증상이 먼저 울려도 근거를 먼저 본다(한 이름으로 뭉개면 증상이 정밀도 고갈을 가린다). 이것은 사용자가 의도치 않은 일치를 알아차리고 **설계 치수를 고치게** 하는 진단이지, 커널의 결정을 대신하는 프롬프트가 아니다.

### 근본 한계

회전이 들어가면 **서로 다른 경로로 같은 위치에 도달해도 실현이 다를 수 있고, exact 술어로도 그 다름을 그대로 반영한다** — 술어는 「주어진 입력에 대해」 정확할 뿐 입력의 근사를 고치지 못한다. **어떤 유한 정밀도로도 보장은 불가**하다(초월수의 상등은 유한 정밀도로 결정 불가) — 정밀도를 올리는 것은 답을 얻는 방법이지 완비성을 사는 방법이 아니다. 다만 **실무적으로는 충분하다**: 상용 CAD 의 uncertainty 가 약 1e-9 인데 f64 는 1000mm 점을 10000번 회전해도 약 1e-10 이고, 커널이 실제로 오르는 폭(회전 3200회에서 3456비트)은 그보다 수백 자리 아래다. **이론적 완벽함은 포기하고 실무적 충분함을 tol 로 보장한다** — 그 tol 을 **측정된 값**으로 들고 다니는 것이 tolerance-fudge 와 다른 점이다. 부호를 확정 못 하는 잔여를 커널은 조용히 0 으로 추측하지 않는다 — 증명되면 일치로 처리해 보고하고, 좁힐 수 없으면 원인을 이름 붙여 거절한다.

### 도메인별 적용

- **평면·이차곡면:** 판정이 다항식이므로 CIP 가 그대로 얹힌다. **CIP 의 가능 여부 = 간접 술어의 가능 여부**다(CIP 는 술어 위의 층이다). 술어가 고차가 되면 **오차 한계 공식을 그 술어에 대해 재유도**해야 한다.
- **자유곡면(M7):** 메시 조합 판정의 필터 + 뉴턴 스냅백 점의 tol 추적·고정밀 재수렴(폭은 같은 방식으로 계산해 정한다). **SSI(교차를 찾는 것)는 CIP 밖**이다 — 위상 존재 문제이지 tol 문제가 아니다.

**계보.** CGAL Lazy_kernel(구간 근사 + DAG 연산 이력 + 애매하면 exact 재평가 + 캐싱)과 구조가 같다. nacre 는 DAG 대신 **점이 드는 모션 사슬**을, 구간 대신 **미리 유도한 오차 한계**를 쓴다(특화 커널이라 가능하고 더 빠르다 — Attene 2020 이 선·평면 교차에서 CGAL lazy 보다 빠른 이유와 같다).

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

## 계획: STEP 백엔드 교체

step-io 를 경량 AP242 출력 전용 크레이트로 교체한다 — 크리티컬 패스가 아니다. 커널의 STEP 출력은 좁은 슬라이스(AP242, 커널 형상 엔티티만)라 최종적으로 경량 라이터가 이상적이지만, 처음부터 만드는 것은 난이도가 높다(「동작 먼저 → 최적화 나중」). step-io 의 AP242 Ed2 스키마 지식·코드젠 타입 정의는 재활용하되 리더 로직은 배제한다. 교체 검증은 step-io 리더를 오라클로: 「step-io 출력 vs 경량 출력」을 되읽어 비교한다.

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
- M7 SSI 탐지 방법의 선택, STEP 경량 라이터.

## 검증·오라클 인프라 (`nacre-validate`, `nacre-oracle`)

### `validate` — 모델 불변식 검사

```rust
pub fn validate(model: &Model) -> Vec<Violation>;   // 빈 벡터 = 유효
```

`validate` 는 순수 함수이고 위반을 **전부** 모아 돌려준다. 인접 캐시는 안에서 새로 짓으므로 호출자가 `Model::rebuild_adjacency` 를 먼저 부를 필요가 없다. `nacre-ops` 는 이 크레이트를 dev-의존으로만 갖는다 — 연산이 스스로 부르지 않고, 테스트·속성 테스트·오라클이 연산 직후마다 부른다. 검사하는 것은 `Violation` 의 변종 그대로다:

| 변종 | 불변식 |
|---|---|
| `DanglingReference` | 모든 handle 이 대상 store 범위 안(`target_index < target_len`). |
| `OpenLoop` | 모든 Loop 가 닫힌다: `end(he[at]) == start(he[at+1])`. |
| `NonManifoldEdge` | 모든 엣지가 half-edge 에 정확히 2회 쓰인다(0 = 고아, 1 = 열린 경계, 3 이상 = 비다양체). |
| `NonOpposedEdge` | 엣지를 쓰는 두 half-edge 의 `forward` 가 반대다(방향 일관성). |
| `NonManifoldVertex` | 정점의 면 링크가 원 하나다 — 모든 엣지가 다양체여도 두 솔리드가 꼭짓점 하나로만 닿는 «핀치»를 잡는다(`nacre_topo::nonmanifold_vertices`). |
| `VertexCarrierMismatch` | 정점 정의의 운반 곡면 종류가 변종의 주장과 맞는다: `ThreePlane` 은 평면 셋, `OnSeam` 은 평면 둘이 아니다(평면 둘의 교선 위 점은 세 평면 교차로 표현되는 진실이다). 불변식은 변종마다 따로다. |
| `EdgeCarrierMismatch` | 엣지가 진술한 `Edge::surfaces` 쌍이 그 엣지를 쓰는 두 면의 곡면 다중집합과 같다. 자기인접 `[s, s]` 는 어느 곡면에도 금지(seam 간선은 없다). 운반 곡면은 유도가 아니라 진술이므로 불일치는 생산자 버그이고, 틀린 운반에서 유도된 곡선 캐시는 조용히 틀린 기하가 된다. |
| `VertexOffCurve` | 엣지에 묶인 정점이 그 엣지 곡선 위에 tol 이내. |
| `VertexOffSurface` | 루프 정점이 그 면의 곡면 위에 tol 이내. |
| `VertexOffDefinition` | **모든** 정점의 캐시 좌표가 자기 정의가 말하는 곡면 하나하나에서 tol 이내. 정의가 진실이고 좌표는 캐시이기 때문이다. tol 은 실측값이 있으면 그것, 없으면 `EPS_CONSTRUCTED`. |
| `EulerParity` · `NegativeGenus` | 오일러–푸앵카레: `V − E + F = 2(S − G) + L_i`. `χ = V − E + F − L_i` 가 홀수면 정수 해가 없고(`EulerParity`), 짝수인데 `G = S − χ/2 < 0` 이면 위상적으로 불가능하다(`NegativeGenus`). |
| `CavityMisoriented` | 공동(내부 void) 셸의 면 법선이 void 안쪽을 향한다 — 부호 있는 자기 부피가 음수. 통째로 뒤집힌 셸은 엣지 대향도 V/E/F/S 도 그대로라 다른 검사가 전부 놓친다. 평면 공동만 검사한다. |
| `FaceMisoriented` | 평면 면의 루프 감김과 진술한 방향의 정합. 루프의 증인(다각형의 뉴웰 면적 벡터, 또는 닫힌 림의 원)과 면이 진술한 바깥 법선(`plane.normal()` × `orientation`)의 단위 내적 `cos` 를 잰다: **외곽 루프는 일치**(`cos ≈ +1`), **내부 루프는 반대**(`cos ≈ −1`). `≈ 0` 은 자기 평면을 펴지 못하는 루프. 불리언의 `collect_planes` 가 `debug_assert` 로 지키는 불변식의 release 쪽 그물이다 — 틀린 플래그·감김은 엣지 대향과 오일러를 모두 통과해 «어느 쪽이 재료인가»로 곧장 흘러든다. |
| `CylinderTruthCacheMismatch` | 원통의 정확한 진실(`CylinderDef`)과 f64 캐시가 같은 원통을 말한다. `field`(`"radius"`·`"origin"`·`"dir"`·`"ref_dir"`)의 성분별로 `|def − cache| ≤ CYL_TRUTH_EPS · max(1, |def|, |cache|)`, `CYL_TRUTH_EPS = 8ε`. 한 곡면의 두 서술은 어긋날 수 있고, 치료는 생산자를 믿는 것이 아니라 소비자 쪽 사후조건이다. 건강한 경로는 절대 1e−17 이내이고, 이 그물이 노리는 결함(seam 타이브레이크 오류로 `ref_dir` 이 ~90° 도는 것, 단위화 누락으로 `dir` 이 길이만큼 커지는 것)은 성분을 ~10¹⁵·ε 움직인다. |

오일러 식의 `L_i`(면 내부 루프 수) 항은 빠뜨릴 수 없다. 관통구멍 정육면체로 검산하면 V16 − E24 + F10 = 2 = 0 + L_i(2) 다 — 내부 루프 항이 없는 축약형은 정상 모델을 불량으로 판정한다.

**validate 는 store 전체가 아니라 live 도달가능 셀을 센다**(supersede 의미론). 오일러의 V·E·F·S·L_i 는 `store.len()` 이 아니라 `live_solids` 에서 도달 가능한 정점·엣지·면·셸·내부 루프 수이고, 다양체 조건(엣지 정확히 2회)도 «전역 2회»가 아니라 «도달가능 집합 안에서 2회»다. supersede 된 옛 셀이 아레나에 남으므로, store 길이로 세면 죽은 셀이 오일러를 깨고 죽은 면의 엣지가 다양체 조건을 깬다. store 에 떠 있는 stray 셀은 불량이 아니라 «죽은 아레나 항목»으로 정당하게 무시된다. **참조 무결성 검사만 store 전체를 훑는다** — 살아 있든 죽었든 handle 이 범위 밖이면 버그이기 때문이다. 이 검사는 `.index()`/`.len()` 만 쓰고(`Store::get` 을 부르지 않는다) 가장 먼저 돌며, 위반이 있으면 거기서 단락한다. 그래서 뒤따르는 도달가능성 순회는 항상 범위 안이다.

**테셀레이션 불변식은 `nacre-tess` 가 자기 테스트로 지킨다**(`validate` 의 변종이 아니다): 출처 정합 — `TessOrigin::OnFace { uv }`·`OnEdge { t }` 의 평가값이 `pos` 와 일치, 틈 없음 — 인접 면 삼각분할이 공유 엣지 polyline 정점을 정확히 공유(테셀레이션 절).

### 속성 기반 테스트(proptest)

랜덤 유효 연산열을 생성해 건다: (1) validate 통과, (2) replay 멱등성 — 같은 로그·같은 cfg → 동일 모델, 그리고 더 강한 «`replay(log)` == 세션 모델, 인덱스까지»(`nacre-ops/tests/invariants/replay.rs`), (3) 변환 불변량 — 강체변환 후 부피·면적 보존, (4) 불리언 대수 — `A ∪ A = A`, `A ∩ ∅ = ∅`, `vol(A∪B) + vol(A∩B) = vol(A) + vol(B)`.

### OCCT 오라클(`nacre-oracle`)

nacre 의 결과를 STEP 으로 내보내 OCCT 가 채점한다: 부피·면적·면 개수·바운딩박스. 불리언은 두 입력 솔리드의 STEP 을 OCCT 가 직접 `bfuse`/`bcut`/`bcommon` 한 결과의 값과 nacre 결과의 값을 diff 한다. 정답지를 든 채 개발하는 장치이며, AI 가 생성한 코드의 «그럴듯하지만 틀림»을 잡는 주 방어선이다. **nacre 쪽 부피·면적은 `nacre-props`(정확 기하 발산정리, 해석적)가 계산하고, 오라클은 그 값을 OCCT 가 같은 STEP 에서 독립 계산한 값과 diff 한다** — 완전히 별개인 두 구현의 일치가 양쪽을 교차검증한다. `props` 가 양의 해석적 부피를 돌려준다는 것 자체가 «우리 STEP 이 OCCT 에 유효하고, 닫힌 방향 있는 솔리드 하나가 통째로 읽혔다»는 검사다. **일치는 형상이 설명한 그것인지 말하지 않는다** — 두 커널이 같은 STEP 을 재므로 픽스처가 어디에 놓이든 맞는다. 그래서 오라클 픽스처는 자기 손 계산 값도 단언하고, 공용 형상(`nacre_ops::fixtures::pocketed_cube`)은 자기 전제를 단언한다.

계획: 허용 편차를 넘은 연산열을 최소화(shrink)해 리포트한다.

**OCCT 는 오라클 전용이며, 연동은 out-of-process 다.** OCCT 를 링크하지 않고, `tools/occt-helper/` 의 헬퍼 프로세스와 STEP 파일로 주고받는다(수송 계층 = `nacre-step`). 이유: C++ 빌드가 Rust 워크스페이스에서 완전히 사라지고, OCCT 가 악조건 입력에서 크래시해도 커널 프로세스가 아니라 헬퍼만 죽는다(크래시 격리 — 수천 회 돌리는 오라클에는 필수). 오라클은 채점만 하므로 history 손실·호출 오버헤드는 무관하다. `nacre-oracle` 은 발행하지 않는 dev 전용 크레이트다.

헬퍼 프로토콜은 고정이다:

```
occt-helper props <in.step>
occt-helper <fuse|cut|common> <a.step> <b.step> [out.step]
```

- `props` 는 솔리드 하나를, 불리언 명령은 **결과**의 성질을 보고한다. `cut` 은 `A − B`. `out` 을 주면 결과를 STEP 으로도 쓴다.
- stdout 은 모든 명령이 같은 key-value 5줄이다: `volume <v>` · `area <a>` · `faces <n>` · `bbox_min <x> <y> <z>` · `bbox_max <x> <y> <z>`. 값은 헬퍼 쪽에서 계산한다. 고정 스키마라 파싱 의존성이 0 이고 Rust 쪽 파서는 하나다.
- exit code: `0` 성공 · `1` 기하 실패(파일 없음·읽기 불가·빈 형상) · `2` 크래시.
- 입력은 **단일 솔리드** STEP 이어야 한다(`to_step_solid`). 전송 가능한 루트가 여럿이면 OCCT 가 `x_1, x_2, …` 로 읽고 불리언은 `x_1` 만 쓴다.

구현은 Homebrew OCCT(`brew install opencascade`)의 `DRAWEXE` Tcl 셸 배치 실행이다(`stepread → [bfuse/bcut/bcommon →] vprops/sprops/nbshapes/bounding`, 콘솔 출력 파싱). 콘솔 긁기가 버전 차이로 취약해지면 같은 줄을 내는 얇은 C++ 헬퍼(`STEPControl_Reader` + `BRepGProp` + `BRepBndLib`)나 uv 관리 Python 환경 + OCP 휠로 갈아탄다 — 프로토콜만 지키면 nacre 쪽 코드는 무변경이다.

### indirect predicates 검증

indirect 술어는 틀려도 대부분의 입력에서 맞는 답이 나와 버그가 가장 숨기 쉬운 조각이므로 오라클을 두 겹으로 건다. (1) **direct 술어 대조(property test)**: 좌표를 아는 케이스에서 정의에서 좌표를 실현해(`realize_vertex`) `geometry-predicates` 의 direct 술어(orient3d 등)에 넣은 부호와, 같은 점을 implicit point 로 둔 indirect 술어의 부호가 일치하는지 대조한다. (2) 계획: **Attene C 구현 실행 대조(dev-only)** — 저자 참조 구현을 별개 프로그램으로 실행해 답만 비교한다(아래 라이선스 규율 — `tools/` 격리, 링크 금지, OCCT 오라클과 동일 논리).

### STEP 라운드트립(자기 출력 되읽기)

nacre 가 쓴 STEP 을 step-io **리더로 되읽어** 구조를 비교한다. 외부 파일 import 가 아니라 «방금 쓴 것»이므로 healing 이 필요 없다. 교차 곡선은 일반 SSI 단계(M7)에서 처음 생긴다. 계획: 그 자리는 `EdgeCache` 를 `{curve, err}` 로 키우는 쪽이고, 교차 곡선이 검증 대상이 되는 것도 그때다.

## 마일스톤 사다리

각 단계는 «동작하는 것»을 남기고 끝난다. 3층에 막혀 전체가 멈추는 구조를 피하는 배치다.

**M1 — 뼈대.** `nacre-math`, Store/Handle, 평면·직선만으로 정육면체를 손으로 조립하고 `validate` 를 세운다. 시각 확인은 **기존 뷰어에 위임**한다(전용 뷰어 크레이트를 만들지 않는다): 정상 결과는 STEP 출력 → step-loupe(구조+검증), 중간·깨진 상태는 Tessellation 을 OBJ/STL 로 덤프 → 맥 미리보기(또는 MeshLab·f3d). 커널 내부까지 보는 인터랙티브 디버그 뷰어(면 클릭 → Handle·정의, 법선 화살표, 엣지 polyline·tolerance 공, validate 위반 하이라이트, 로그 스텝별 재생)는 워크스페이스 밖 별도 앱이다 — **발견된 점**과 tolerance 처럼 STEP 에 담기지 않는 내부 정보가 생기는 불리언 단계의 디버깅 생명줄이다.

**M2 — 스케치와 돌출, STEP 내보내기.** 2D 프로파일 → extrude, 연산 로그와 replay, 오일러 연산. STEP 은 어댑터다: 커널 Model → AP242 엔티티 번역만 자체 구현하고 직렬화는 step-io 백엔드에 위임한다(커널은 백엔드 무지 — 헬퍼 프로토콜과 같은 격리). 어댑터가 내보내는 엔티티를 AP242 커널 형상 집합(cartesian_point·direction·line·circle·b_spline_curve/surface·plane·cylindrical_surface·advanced_face·closed_shell·manifold_solid_brep·surface_curve 등)으로 명시 타입 고정해, 나중에 경량 라이터가 감당할 범위를 붙박는다. 커버리지는 마일스톤을 따라 넓힌다(평면 → 곡선·곡면 → 트리밍). 목적은 기능만이 아니라 검증이다: 작은 형상에서 좌표계·방향·면 orientation 함정을 먼저 밟는다. 검증 뷰어는 두 겹이고 상보적이다 — **step-loupe**(step-io 기반 웹 뷰어: report 가 드롭·고아·비표준 엔티티를 표시해 어댑터 버그를 구조적으로 잡는다; 자기 출력 되읽기의 GUI 판) + **FreeCAD 등 독립 OCCT 기반 뷰어**(step-io 를 공유하지 않는 교차검증). STEP 출력은 오라클의 수송 계층이기도 하다. 최종 경량 라이터로 교체할 때는 «step-io 출력 vs 경량 출력»을 step-io 리더로 비교해 교체 안전성을 자동 검증한다.

**M3 — 곡선 기하.** Arc, Cylinder, NurbsCurve/Surface 평가(The NURBS Book 기준 구현 + 수치 미분 대조 테스트). tess 출처 태그, tolerance 재계산(같은 모델, tol 3단). proptest.

**M4 — 면 위 작업.** 면의 스케치 프레임 위 돌출 + 불리언(pad·pocket 은 kit 의 조합). 여기까지 모델의 모든 점이 **구성된 점**이다(실측 tol 0). OCCT 오라클이 이 단계부터 돈다.

**M5 — 평면 솔리드 불리언.** 모든 면이 평면인 솔리드 간 fuse/cut/common. 평면–평면 교차는 닫힌 형식의 직선이다(SSI 행진·Newton·캐시 불필요). 꼭짓점은 평면 3장의 교차를 **좌표로 만들어 병합하지 않고 implicit point 로 두고**, 내/외·orientation 부호는 그 정의를 **indirect orient3d**(좌표를 만들지 않는다)에 넣어 정확히 판정한다 — 좌표를 만드는 순간의 오차·불일치를 원천 차단한다(Attene 2020). 세 엣지가 하나의 Vertex Handle 을 공유해 봉합 문제의 본체를 우회하고, 정점은 자기 정의를 든다(`Vertex::ThreePlane`). 이 세계에서 강건 불리언은 연구가 아니라 꼼꼼한 케이스워크(공면, 엣지–엣지 퇴화)다. 선·평면(다항식) 교차점은 indirect predicate 이론이 가장 깔끔하게 도는 영역이다. 커버리지 밖 입력은 `Rejected` 로 정직하게 거절한다. **발견된 점**의 경로·국소 tolerance·정의에서의 실현이 여기서 실전에 들어가고, 이 단계에서 «OCCT 없이 직동하는, 평면 위주 기계 부품을 STEP 으로 내보내는» 커널이 된다. 알고리즘 참고: Manifold(Apache-2.0 — 차용·번역 가능), Hoffmann 등 문헌. Truck 대비 벤치마크·정밀도 비교(전역 1e-6 폴리라인 vs 정점별 실측 tol + 닫힌 형식)는 수치로 보일 수 있는 차별점이다 — 공개 지표 후보.

**indirect predicates 는 자체 구현(clean-room)이다.** `nacre-predicates` 가 implicit point 표현과 indirect 술어를 자체 구현하고, 확장 산술 바닥은 `geometry-predicates`(elrnv, MIT/Apache)를 재사용한다 — 이 크레이트는 finished `orient3d` 뿐 아니라 Shewchuk expansion primitive(`two_product`·`two_sum`·`expansion_sum`·`scale_expansion_zeroelim` 등, `predicates` 모듈에 공개)를 노출하므로 expansion 산술을 자체 구현할 필요가 없다. primitive 는 `[lo, hi]` 순서이고, orient3d 부호 규약은 `det[a−d, b−d, c−d]` 다(golden 테스트로 고정). **라이선스 엄수** — Attene 참조 구현(LGPL)은 «돌려서 답 비교는 자유, 열어서 코드를 보는 건 MIT/Apache 소스만»: (a) **작성 중 소스 열람 금지**(LGPL C++ 를 열어 함수 대응·로직 흐름을 따라가면 2차적 저작물). 참고처는 논문(Attene 2020, arXiv 2105.09772; Cherchi 2020 mesh arrangements·2022 interactive booleans; Lévy 2024; Shewchuk 1997)과 `geometry-predicates` 소스로 한정한다. (b) **완성 후 실행 대조는 허용·권장**(dev-only, `tools/` 격리, 우리 크레이트에 링크 금지 — OCCT 오라클과 동일 논리).

**M6 — 평면+이차곡면 불리언.** 평면∩실린더(타원), 평면∩구(원), 평면∩원뿔 — 여전히 닫힌 형식이라 행진이 필요 없다. 실제 기계 부품 면의 대다수가 평면+실린더+원뿔이므로 여기까지로 실용 커버리지의 대부분을 확보한다. 원통 쌍 가운데 **면이 만나지 않는** 쌍(십자로 관통한 두 스터드처럼 축은 만나도 남은 면은 비킨 것)은 새 기하가 아니라 게이트의 명제 문제다 — 옆면의 축 구간으로 비켰음을 증명하면 배열이 그대로 조립한다. 새 기하가 필요한 것은 면이 **만나는** 쌍(4차 곡선)뿐이다. 직교 같은 특수 배치부터 하나씩 열지 않는다 — 옆면 차트의 곡선 어휘와 교선이 배열에 들어오는 길을 먼저 넓힌다(`todo.md` 「원통–원통」).

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
| f64 필터 없는 정확 술어를 판정의 주 도로에 | `indirect_plane_side` 를 모션 없는 `orient3d` 에 쓰자 축정렬 fold 가 14배 느려졌고, `indirect_cmp_coord`(≈2.7µs, 호출마다 할당 — 필터 있는 인증 도로는 ≈1µs)는 지름길을 타는 평면이 늘자 fold 를 11~17% 늦췄다(한 스레드에선 빨랐다 — 병렬에서 할당자가 줄을 세운다). 둘 다 필터를 달자 기준보다 약 30% 빠르다. 필터는 정확 술어와 한 벌이다 | — |
| 공간분할(모델을 리전으로 잘라 잎에서 불리언) | 전부 초록인 채로 느렸다: 회전 fold 80 7.4→8.6s, 축정렬 1.7→2.4s. 정점을 용접하고 모서리를 공유하는 b-rep 에서는 분할면에 캡이 필수이고, 자르고 다시 꿰매는 일이 순비용이다(다각형 수프를 내는 EMBER 는 이 비용이 없다) | 되살릴 자리 `d3e9984`·`bfd01e6`·`12ec34d`·`9d69786` |
| 모션 사슬을 질의마다 걸어 접기 | 질의가 변환마다·평면마다 나오므로 깊이의 제곱이다: 기록된 사분각 회전 2,000번 0.20초 → 6.7초, 37° 4,200번 1.44초 → 2.14초. 노드가 태어날 때 부모의 접기에 자기 하나를 합성해 둔다(`motion_folds`) — 둘 다 기준선 | 실험, 되돌림 |
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

### 정밀도·판정

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
| 향 관계를 캐시 법선으로 읽기(내적·«첫 0 아닌 성분»의 부호 비교) | 내적은 크기 논증으로 건전해도 판정이 진실을 안 읽고, 첫 성분 비교는 크기 논증도 아니다 — 참값이 아주 작은 성분(10⁻¹⁷)의 부호를 정당한 반올림 상이 뒤집으면 `outward_fix` 가 반대를 답했다(단위 심기). census·스위트·perf 에서 진실과 불일치 0(첫 성분이 최대가 아닌 호출은 census 에 0). 향은 `sense`·`Orientation`·세계 이름의 향(`WorldName`)과 증인 정의(`normals_agree_judge`)가 답한다 | — |
| 간선 곡선의 종류(원·룰링)를 평면·원통 캐시의 허용오차로 고르기(`(n·d)² ≤ 1e-18·|n|²|d|²`) | 진실이 비스듬한데 캐시는 원을 답했다 — 무리수 축의 테스트 원통(캡을 실현한 f64 에서 들어올린 문)의 림 689: 진실의 캡은 축에 비스듬해 타원인데 간선 캐시는 원이었고, 판정 경로(`loops.rs`·`corners.rs`·원판 캡)가 그 종류를 읽었다. 축에서 `1e-10` 라디안 기운 평면을 심으면 원, 평행에서 `1e-10` 벗어난 평면은 룰링 선이 됐다. 종류는 정확한 관계(`AxisRelation`)가 정한다 | — |
| 교선 위 점을 비교할 축과 부호를 평면 캐시로 고르기 | 참 0 성분을 반올림이 0 에서 벗어나게 하면 그 축이 골라지고, 교선 위 점은 그 좌표가 모두 같아 서로 다른 두 관통점이 «같다»로 읽힌다 — 기울어진 원통(십진 직교 틀) 심기에서 Fuse/Cut/Common 두 순서 여섯이 전부 `CoincidentFeatures` 로 거절됐고, 유리수 방향(`stored_line_dir`)으로는 여섯 모두 해석값 부피다. census(호출 58,322)에는 그 인구가 없었다 | — |
| 접히는 사슬 아래의 참값 0 좌표를 접기 대신 일치 정밀도 규칙(과 캐시 도로의 두 번째 단)으로 풀기 | 풀리기는 하지만 사분각 회전 2000번 빌드가 74 → 352ms(4.8배)였다 — 그 사슬은 노드마다 유리수 사상으로 접혀 있어 정점이 정확한 유리수로 읽히고(`chain_point_rat`), 그 도로는 오히려 빠르다(74 → 64ms). 규칙은 접히지 않는 사슬의 몫이다 | 실험, 되돌림 |
| 사슬이 안 접히는 평면의 캐시를 진실의 세 점을 f64 로 재생(`replay`)해 짓기 | 생산자의 f64 보다 참값에서 멀다 — 스위트 22,759 평면에서 256비트 참조에 생산자 값이 더 가까운 것 11,457, 재생이 더 가까운 것 2,439, 같은 것 8,863(census 출력 34행 추가 이동). 둘 다 반올림이다 | 실험, 되돌림 |
| 모션 노드를 기록할지를 캐시 좌표(꼭짓점·평면 원점·원통 축 원점)가 f64 에서 정확히 옮겨지는지로 고르기 | 두 정확 표현(옮겨 적기·사슬) 가운데 하나를 캐시가 골랐고, 표현을 읽는 규약이 그 선택을 물려받았다 — 같은 자리의 면이 `0.5` 로 옮기면 세계 진술, `1/3` 으로 옮기면 사슬이 되어 거울 뒤 스케치 프레임의 손방향이 갈렸고, 같은 벽이 `(name, motion)` intern 에서 두 핸들이 되어 병합이 고정밀 판정에 기댔다. census 이동 419 중 186, 스위트 약 2,765 번이 그 판정에서 갈렸다(바꾼 뒤 census 출력은 비트 동일). 진술이 `i128` 안에서 옮겨지는지만 묻는다 | — |
| 평면 클래스를 캐시(평면 캐시 계수의 비례 · 면 꼭짓점 좌표의 정확 `orient3d`)로 합치기 | 반올림된 상이 같은 서로 다른 두 평면이 한 클래스가 된다 — 심은 두 배치에서 `Common` 의 경사면이 남의 평면을 들었고, 닿지도 않는 두 몸의 Fuse 에서 한 몸의 벽이 다른 몸의 평면을 들었다(`validate` 무위반). census 코퍼스(쌍 30,009)에는 그 인구가 없었다(불일치 0). 병합은 세계 이름, 없으면 정의를 읽는다 | — |
| 정확 지름길의 자격을 캐시끼리의 대조로 주기(평면 캐시의 계수 ↔ 면 꼭짓점의 좌표 캐시) | 둘 다 반올림된 상이라 «반올림된 평면»을 인증한다 — `3·0.1 = 0.3` 인 모서리를 벽 `y = 3x` 밖으로 판정해(진실과 272번 다른 답) 유효한 교집합을 `ZeroLengthEdge` 로 거절했다. census 코퍼스(판정 186만)는 이 배치를 안 때렸다 — 반올림이 유리수 일치를 깨는 배치를 손으로 심어야 보인다. 지름길은 평면의 이름을 읽는다 | `b1d39469` |
| 잔차(점에서 담체 `f64` 평면까지의 거리)를 좌표의 경계로 읽어 자기접촉 체의 상자를 부풀리기 | 잔차는 거리 «하나»라 좌표가 참값에서 얼마나인지를 말하지 않는다 — 참값의 최근접 f64 까지를 못 덮은 seam 점이 스위트 3,415/9,494(다른 것 중), 그중 155 는 천 배 넘게 모자랐고 census 190. 체는 정확 판정 «앞에서» 후보를 버리므로 모자란 상자는 진짜 접촉을 조용히 버릴 수 있다. 체는 `Bounded` 의 증명된 경계를 읽고, 경계가 없으면 버리지 않는다 | — |
| 「애매하면 사용자에게 묻는다」 | 회전 코퍼스에서 상승 2,490 중 2,260 이 판정 불가(연산당 수백 건), 가장 느슨한 한계가 2⁻²²⁹, 그리고 f64 출력은 「떨어져 있다」는 답을 표현하지 못한다. 선택지가 하나뿐인 질문은 질문이 아니다 — 근거를 붙여 보고한다 | — |
| 정밀도 사다리(256→512→1024 배증) | 256·512·1024·2048 비트에서 결과가 비트 동일 — 애매함은 정밀도 부족이 아니라 구조적 0 이다. 필요한 비트는 조건수 `C` 에서 계산해 한 번에 점프한다 | `cc37d86` · `f7a4fec` |
| 고정 판정 정밀도 | `C` 는 회전 1회당 1비트 자란다: 256비트에서 245회 회전부터 거절, 그것도 사용자가 대응할 수 없는 하류 증상의 이름으로 | `f7a4fec` |
| f64 `tol` 을 조건수로 읽어 임의정밀 시행 패스를 생략 | 2000/2000 에서 하회했다 | 코드 무변경 |
| double-double(`twofloat`/`qd`)를 고정밀 층으로 | 영점 근처 cos 오차 ~1.8e-16 — 부호 판정이 일어나는 곳이 영점 근처다 | — |
| 접선형 회전 tol(`da·|y|`, `da·|x|`) | 축 근처에서 건전하지 않다(cos/sin 이 독립적으로 반올림돼 반경 방향 오차가 생긴다). `frame3` 의 테스트가 `!tangential_sound` 를 단언한다 | — |
| libm 의 sin/cos 를 단일 진실원으로 | 같은 프로세스·같은 각에서 호출부마다 1 ulp 다르다(리터럴 각은 LLVM 이 컴파일 타임에 접는다) | — |
| SoS(ε-nudge)로 퇴화를 없애기 | 첫 벽(64건)은 사라지지만 exact 삼중 기판의 둘째 벽(69건)으로 옮겨갈 뿐, 완주 0. 퇴화는 이름 붙여 거절하거나 엔진이 구조로 답한다 | 스파이크, 프로덕션 무변경 |
| 도달 가능한 실패를 `debug_assert` 로 막기 | release 에서 꺼진다 — 체인 접촉 절단이 그 단언을 반증했다. 도달 가능하면 이름 있는 거절이다 | — |
| 향이 걸린 값(평면 캐시의 법선)을 판정 증인 삼각형으로 펼치기 | 이름 있는 넓은 만남 평면의 증인은 그 프레임의 탐침이라 회전 방향이 진술 점의 것이 아니다 — `sense` 를 곱하자 캐시가 진실과 반대를 봤다(`a_datum_through_wide_meets_keeps_its_name` 에서 감김 `debug_assert` 가 잡음). 법선은 모션 전 이름 법선을 사슬로 운반하고, 증인은 계수(만남·판정)에만 쓴다 | 첫 시안(`e9a6e3ab` 본문) |

### 타입·진실

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
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
| 좁은 이름을 별도 곁표로 | ~99% 가 중복 저장 — `PlaneName = Narrow | Wide` 한 enum + `narrow()` 투영 | — |
| wide 계수를 f64 `Expansion` 에 | 조각도 f64 라 지수 상한 ~2¹⁰²³, wide 이름은 ~2²²⁹¹. 그릇은 BigInt 정수 산술이다 | — |
| 판정층에 `Through` 전용 평면 기계를 짓기 | 판정 표의 계약은 «평면 위 정확한 세 점»이고 프레임의 probe 가 이미 준다 — 새 기계 0 으로 열렸다 | — |
| 뒤집힌 `Canonical` 프레임의 기준 방향(`ẑ×n`)을 뒤집기 **전** 정준 계수에서 유도하기 | `û` 가 이름의 정규화 부호(첫 0 아닌 성분 양수 — «같은 평면인가»의 열쇠 관례)를 따라, 좁은 도로만 `û` 를 두고 넓은·이름 없는 도로는 `û` 를 돌렸다 — `z=0` 을 `−z` 로 보는 프레임에서 `[0,1.5]²` 가 좁은 도로로 `[0,1.5]×[−1.5,0]`, 이름 없는 도로로 `[−1.5,0]×[0,1.5]`. 이름 없는 도로는 계수의 부호를 증명하지 못해 그 규약에 철자가 없다. 스위트의 그런 프레임 ~156 이 잠금 없이 움직였다(census 0) | `6b80011d` |
| 평면의 방향을 점 순서에서 뽑기 | 점 순서는 향을 말하지 않는다: 씨앗·datum·원통 밑캡은 관례로 점과 반대인 향을 들고(비씨앗 `Known` 342, 씨앗 전부), `Stated` datum 은 순서를 프레임 배치(`ref_dir`)로 읽고, `Through` 삼중은 정렬된 열쇠다. 향은 진실의 한 비트(`Surface::Plane.sense`)가 든다. 이름은 「어디」에 답하지 「어느 쪽」엔 못 답한다 | — |
| 평면 캐시를 «수선의 발 앵커 + 정준 행»으로 | census 32행이 움직이고 8건이 넘친다(`n·n` 이 계수를 제곱한다). 살아남은 것은 «진실 삼중의 첫 점» 앵커뿐 | `8b53a2c` |
| `rotated` 플래그를 풀어 회전 클래스를 원통 관문에 들이기 | 그 플래그는 정확 지름길이 평면의 이름을 세계 진술로 읽는 허가다 — 모션 있는 평면의 이름은 모션 전 프레임을 말한다. 풀었을 때는 클래스 병합이 면 꼭짓점 캐시의 삼각형 판정으로 가서 공유 벽이 두 몸으로 갈렸다(그 판정은 은퇴했다 — 병합은 세계 이름을 읽는다) | 되돌림(`planes` 주석) |
| 접미사만 다른 형제 함수를 «중복»으로 세어 합치기 | 일곱 갈래를 열어 보니 진짜 중복은 0 이었다. `det3`(정의 10개)·`three_planes_rat`/`_big` 은 층마다 숫자 타입(`Expansion`·`f64`·`Bounded`·`HpBounded`·`BigInt`)이 다른 **정밀도 사다리**이고, 나머지는 트레이트 기본 구현+impl, 인자가 다른 별개 함수, doc 에 사유가 있는 한 줄 편의 문이었다. 세기 전에 **인자 타입**을 먼저 본다 | — |
| 원통 옆면을 seam 엣지(A — 세로 seam 간선 `[c, c]` 를 두 번 쓰는 한 루프)로 짓기 | 각도의 시작점(`ref_dir`)이 위상에 새겨진다: 잘린 림마다 θ = 0 정점과 호 쪼개기(census 154·스위트 1,980), seam 에 걸린 구멍을 잇는 slit 걸음, θ = 0 에서만 보는 집힘 그물, self-pair 담체 규칙이 붙고, 같은 모양이 seam 의 자리에 따라 다른 위상이 된다. 옆면은 림 두 루프로 닫는다(「seam 엣지 방식」 B) | `e60ab2b9`, `497b54b5` |

### 정리

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
| 긴 함수를 길이로 쪼개기 | 쪼갤 근거는 «안의 한 단계를 다른 호출자가 이름으로 부를 만한가»다. 300줄이 넘는 열셋 가운데 그런 호출자를 든 단계는 **하나**였고(감사가 손으로 다시 짓던 방향 분할 → `direction_partition`), 나머지는 블록이 바깥에서 받아야 할 지역변수가 1~32개다 | `abaa95d`·`9ab9ac3` |
| 주석 비율(36~45%)을 줄이기 | 주석의 대부분은 사용자에게 가는 doc 15,182줄이고, 본문 주석의 가장 긴 덩어리 열 개(294줄)에서 되풀이는 **14줄(5%)** 뿐이었다 — 나머지는 불변식·측정치·거절한 대안·위험이다. 판단 규칙은 `overview.md` 문서 규칙 7 | `82f7b6c`·`56e6756` |

### 엔진·연산

| 하지 말 것 | 이유 | 커밋 |
|---|---|---|
| 실패하는 배치마다 경로를 덧대기(케이스워크) | 배관(detector·생존표)이 다음 엔진으로 전달되지 않고 경우의 수가 폭발한다. 평면당 셀 복합체 하나(winding 단일 배열)를 일반화해서 닫는다 | — |
| 평면 DCEL(`walk_cells`·`nest_cells`·`label_cells`)을 원통 차트에 이식해 셀 기계를 하나로 | 단계는 같아도 도구가 같은 단계가 **0** 이다: 셀은 정점의 각도 순서(정확 술어) 대 격자(모서리 넷), 바깥/구멍은 감김수 `#(−1) == 성분 수` 대 «가장 낮은 칸의 밑변» — 차트를 감는 띠에는 평면 감김수가 없다, 딱지는 무한 뿌리에서 전파 대 이웃 클래스에서 읽기 — 환면에는 무한 영역이 없다. 모양이 닮은 것은 사이클 긋기 수십 줄이고 순서의 출처가 다르다(궤도 next 대 좌회전 규칙) | — |
| 볼록 전용 fast path | 「같은 답을 더 빨리」가 거짓이었다: 일반 경로가 답하는 입력을 거절하고 진리를 f64 로 정한다 | — |
| 구멍을 외곽에 브리지로 꿰맨 뒤 귀 자르기 | 브리지가 정점을 반복시켜 two-ears 정리가 깨진다 — 무작위 구멍 면 4,000 표본의 27.6% 실패, 자기교차 후에도 삼각형을 뱉는 조용한 오답의 여지. 스윕으로 간다 | — |
| 곡면 tess 에서 예산을 넘는 간선을 차트 중점에서 쪼개며 라운드 돌기 | 발산한다(179 → 386 → 856 → 1952 → 4412). 이탈이 Δv 에만 의존하고 영역이 이방적이라 어디서 쪼개도 새 간선이 비슷한 Δv 를 갖는다 — v 방향의 띠가 필요하고 그것은 점을 미리 대야 얻어진다 | — |
| 곡면 tess 에서 지름원 encroachment 로 침범 후보를 버리기 | 직선 경계는 폴리라인이 두 점뿐이라 길이 4 솔기의 지름원이 한 면의 후보 360개 전부를 삼킨다. 자리 없는 후보는 삽입 탐색이 이미 거절한다 | — |
| Sutherland–Hodgman 으로 비볼록 링 클립 | 교차를 링 순서로 짝지어 떨어진 조각을 클립선을 따라 이어 버린다. 선 순서로 짝짓는다 | `d3e9984` |
| 반사를 켤레 회전으로 나르기(`M∘R = (M R M⁻¹)∘M`) | 정확했지만 한 질문에 두 메커니즘이 생기고 사슬이 방금 한 반사를 스스로 말하지 못한다 — 미러로 도달한 벽과 이동으로 도달한 같은 벽이 1 ULP 갈라져 한 부품이 두 몸통이 됐다. 반사는 사슬의 노드다 | — |
| 빈 불리언 결과를 오류로 | 빈 결과는 실패가 아니다 — `Ok(vec![])` | — |
| 공유 룰링 위의 점을 피연산자의 세 평면 모서리 × 원통으로 씨앗하기 | 무한 원통 위의 모서리까지 묶는다 — 필렛 사분원 밖에서 슬래브 모서리가 필렛의 무한 원통 위에 있어, 대표가 된 관통 이름이 결과에 없는 평면을 불러 빌드되던 필렛 슬래브가 `VertexNamesAbsentSurface` 로 깨졌다. 그리고 모서리가 아닌 점(프리즘이 캡 너머로 뻗을 때 캡 위의 점)은 못 본다 — 그 배치는 전부 `CoincidentNodes` 로 남았다. 씨앗은 게이트의 기록(옆면 «면» 위의 선) × 축을 가로지르는 클래스다 | — |
| 공유 룰링 정거장을 건너 옆면 셀을 잇기 | 두 벽이 모두 원통으로 들어가는 쐐기의 `A − B` 가 한 몸통으로 나왔다 — 옆면이 선을 통과하는 한 면이고 두 벽 면만 그 선을 간선으로 나눠, 스스로 닿는 솔리드가 `validate` 깨끗·부피 항등식 성립인 채로 나갔다. 두 현 사이의 방이 결과에 없으면 옆면을 그 선에서 끊는다(`slit_at`) | — |
| 동일평면 면을 이어붙여서 정리 | 비스포크였다. 경계 대수로 일반화한다 — 이어붙이지 말고 내부를 지운다 | — |
| 45° 폴드를 위한 figure-8 재봉합 기계 | 지역-한계 사슬의 깊이는 정확히 2겹이었다. 병합은 기권하고 전역 판정(`SelfTouchingResult`)이 발행 전에 말한다 | — |
| 면 위 직접 구성(`ImprintSketch`·`raise_region`) | 소비자 0 이었고 커널 유일의 동일평면 인접면 생산자였다. pad = extrude + Fuse, pocket = extrude + Cut. Split Face 의 소비자(구역별 재질·FEA 경계조건·파팅라인)가 생기면 두 제거 커밋의 revert 가 출발점이다 | — |
| pad/pocket 을 커널 연산(`PadOnFace`·`PocketOnFace`)으로 | 새 정확 술어 없는 «면 프레임 + 돌출 + 불리언» 래퍼였다(설탕 판별 기준). 커널에 남은 이유는 «면에서 안쪽으로»를 기본 어휘가 말하지 못해서였다 — 돌출이 부호 있는 거리를 받자 사라졌다. 바닥 면을 돌려주려는 캡 회수가 **모든** 불리언에 입력 평면 면마다 클래스·향 표(`ClassOf`, 면마다 `face_facing`)를 계산하게 했고, 깊은 pocket 을 `PocketNotBlind` 로 거절하게 했다. kit 은 캡을 읽지 않았고, 가른 몸의 나머지 조각을 값에서 잃었다 | — |
