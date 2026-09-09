# nacre — 하이브리드 CAD 커널 설계

정확한 기하를 진실로 보관하는 고전 b-rep의 골격 위에, Fornjot에서 검증된 위생 규율(append-only 단일 참조, 위상과 근사의 동시 구축)을 얹은 설계다. 목표는 정밀 기계 CAD(STEP 입출력 포함)이며, 개인 + AI 협업 개발을 전제로 실수 여지를 줄이는 인프라를 1일차부터 포함한다.

## 0. 설계 원칙

**범위 (비목표 명시).** nacre는 순수 기하 커널이다: 기하·위상·수치 봉합(tolerance)·연산·tessellation·검증까지만 다룬다. GD&T/PMI(제작 공차), 스타일·색상·레이어·visibility, 제작자·승인·날짜, 제품 구조·조립·리비전은 커널 범위 밖이며, 커널을 사용하는 응용 프로그램의 영역이다. 커널은 이들을 위해 단 하나만 제공한다 — 영구 유효한 Handle(append-only의 부산물). 응용은 `HashMap<Handle<Face>, 응용데이터>` 사이드카로 무엇이든 매달 수 있다. nacre-step은 STEP의 형상 서브셋만 커널로 번역하고, 비형상 엔티티는 해석 없이 무손실 패스스루로 보존해 라운드트립을 지킨다.

이 문서의 모든 결정은 다섯 가지 원칙에서 나온다. 첫째, **정확 기하가 진실이다** — 평면·원통·NURBS는 해석적 형태로 영구 보관하고, 메시는 파생물이다. 둘째, **모든 객체는 append-only 저장소에 딱 한 번 존재하고 Handle로만 참조된다** — 동일성 질문을 좌표 비교(기하)가 아니라 인덱스 비교(명목)로 바꾼다. 셋째, **연산 이력을 보존한다** — 모델은 연산 로그의 재생 결과다. 단, 보장 범위를 정확히 한다: replay는 동일 로그·동일 파라미터에서 동일 모델을 보장하고(undo/redo의 기반), tolerance 변경 재계산은 위상·기하를 불변으로 둔 채 tessellation만 재생성한다(§5). 로그 중간의 파라미터를 수정하는 파라메트릭 편집은 v1 비목표다 — 연산이 원시 Handle을 참조하는 한 상류 수정이 하류 Handle 번호를 밀어내기 때문(topological naming 문제). 진화 경로는 §6에 기록. 넷째, **tessellation은 일회용이 아니라 출처 태그가 달린 1급 부산물이다** — 모든 삼각형·근사점이 자기가 어느 정확한 면/엣지에서 왔는지 안다. 다섯째, **tolerance는 "발견된" 교차에만 존재한다** — 구성 시점에 동일성을 아는 요소에는 tolerance 개념 자체가 없으며, 이 구분을 타입 시스템에 새긴다.

## 1. 크레이트 구조

의존 방향은 아래에서 위로만 흐른다. 순환 의존 금지.

```
nacre/                    # 워크스페이스. 최상위 `nacre` 크레이트는 파사드(재수출 전용)
├── nacre-store      # typed-index 저장소: Store<T>/Handle<T> (기하·위상 무지)
├── nacre-math       # 자체 선형대수: Point<D>·Vector<D>·변환
├── nacre-scalar     # exact 유리수 값 엔진: Rat·Angle·Axis/Rotation/Isometry·Orient
├── nacre-predicates # exact f64 부호 술어(indirect predicates); geometry-predicates 위·standalone
├── nacre-cip        # toleranced 부호 술어(회전): kernel(WitnessPoint 판정) + predicate(평면 배열 술어). predicates의 쌍둥이
├── nacre-geom       # f64 기하 캐시(Surface·Curve)·교차(intersect 격리; scalar 의존 — Rat 링 술어 쌍둥이). 평면의 진실은 topo 의 SurfaceTruth(S6b)
├── nacre-topo       # b-rep 위상: Vertex/Edge/Face/Shell/Solid·half-edge·Model
├── nacre-tess       # tessellation: 출처 태그·증분 갱신
│   └── polygon      # 평면 다각형 삼각분할: y-단조 분해 + 단조 삼각분할 + Delaunay 플립
├── nacre-validate   # 불변식 검사: 오일러-푸앵카레·watertight·방향성·참조 무결성
├── nacre-props      # 질량 특성: 부피·면적(해석적, tess 무관)
├── nacre-ops        # 연산: sketch/extrude/revolve/pad/pocket/boolean. 부울 = 면당 평면 arrangement 엔진
│                    #   reuse: 상대가 닿을 수 없다고 증명된 평면 클래스는 배열하지 않는다 (§6.3)
├── nacre-step       # Model→STEP(AP242) 내보내기 어댑터
├── nacre-oracle     # [dev] OCCT 비교 하네스
└── tools/occt-helper/  # 워크스페이스 밖 헬퍼: brew OCCT(1순위) 또는 uv+OCP(폴백) — §7
                        #   OCCT는 오라클 전용 — 제품 경로에 위임 없음 (§6, §8)
```

`nacre-geom`과 `nacre-topo`가 서로를 모르게 하는 것이 중요하다. 기하는 위상을 모르고(순수 수학), 위상은 기하를 Handle로만 참조한다. robustness가 첨예한 코드(교차·분류)는 전부 `nacre-geom::intersect` 한 모듈에 격리한다. 사용자는 파사드 크레이트 `nacre` 하나만 의존하며, 인터랙티브 스크립트 앱 등은 이 워크스페이스 밖의 별도 프로젝트로 둔다.

**파사드 `nacre` (2026-07-27 구현).** 열 개 층을 **모듈로** 재수출하고(`nacre::topo`·`nacre::ops`…), 자주 쓰는 것은 `nacre::prelude`에 담는다.
- **평면 재수출은 하지 않았다.** 층 전체 공개 이름 93개 중 충돌은 **`Rotation` 하나뿐**(scalar의 정확한 정의 vs topo의 이력 노드)이라 평면화가 가능했지만, **층 분리가 이 설계의 뼈대**여서 이름공간에 남긴다.
- **`nacre-scalar` 재수출은 선택이 아니다.** `Operation::Transform { isometry: Isometry }`·`Mirror { axis: Axis, offset: Rat }`가 scalar 타입을 ops의 공개 API로 새어 보내는 **유일한 지점**이라, 없으면 소비자가 그 op을 만들 수조차 없다. 반대로 `nacre-cip`·`nacre-predicates`는 공개 API에 새지 않아 재수출하지 않는다(퍼블리시는 필요 — ops의 하드 의존).
- **prelude의 기준은 확인 가능한 성질이다**: *"`Operation`의 모든 변이가 prelude 이름만으로 만들어진다."* 내용은 실제 소비자(playground)의 import를 측정한 뒤, 그 코드보다 나중에 생긴 표면(스케치 앞문·파생 조회·거절 사유)을 합집합으로 얹어 정했다 — 낡은 소비자만 보면 최신 API가 빠진다.
- **기능**: `parallel`(기본 on)이 `nacre-ops/parallel`로 전달된다. `nacre-ops`를 **`default-features = false`로** 매달아야 소비자가 끌 수 있고, 그러지 않으면 `--no-default-features`에도 rayon이 들어온다(실측). 이 속성은 매니페스트에 살아 훅이 검사하지 않으므로 **테스트가 `Cargo.toml`을 직접 확인**한다. `parallel`은 **`Sync` 스위치이기도 하다**(cip의 hp 캐시가 `Arc`↔`Rc` — 워커들이 평면표를 *공유*해야 하므로 하중은 공유 참조 쪽이다). `test-util`은 topo의 테스트 전용 `add_cuboid`를 전달한다. **순차 조합은 반드시 `-p nacre-ops --no-default-features`로 확인한다** — 워크스페이스를 통째로 지으면 `nacre-oracle`·`nacre-props`의 dev-dependency가 기본 피처를 도로 켜서 아무것도 안 재게 된다(§6.2).
- **검증은 "`nacre::` 경로만으로 끝까지 가기"다** — 크레이트 문서의 예제(doctest)와 `tests/facade.rs`. 재수출이 빠지면 컴파일이 깨진다. 문서 예제는 **프로덕션 API만** 쓴다(`test-util`이 필요한 예제는 독자가 실행할 수 없다).
- **`Document`(op 로그 + Model + tess 캐시) 번들은 아직 두지 않았다** — §2/§8이 이 층을 지목하지만 의미론을 정하는 새 타입이므로, kit이 무엇을 원하는지 보이기 전에 정하지 않는다.

**편의 레이어 `nacre-kit` (워크스페이스 밖, 별도 리포 — 2026-07-26 결정, 미착수).** 코드-CAD 스크립트와 커널 사이의 층: 다중 솔리드 값(compound), 값 의미론(재사용 시 `Copy` 자동 삽입), 다인수 fuse/cut/common(fold), 프로파일 헬퍼와 섬-분해 호출, 패턴·미러, 에러의 사람용 매핑, 표시 메타데이터(색·투명도 — 커널 비목표라 여기가 제자리). **Rust로 두는 이유:** 헤드리스 `cargo test`가 되고, 술어 인접 로직이 exactness 도구가 있는 쪽에 남고, 프론트엔드를 교체해도 살아남고, wasm 경계가 함수 하나로 유지된다. 경계 규칙은 overview.md의 "설탕 vs 커널 판별 기준"이고, **문법·의미론과 그 결정 이유는 `nacre-kit` 리포의 `docs/syntax.md`·`docs/decisions.md`에 있다**(여기에 복사하지 않는다 — 두 곳에 같은 내용이 있으면 어긋난다). **전제였던 것 — ✅ 2026-07-27 해소:** 파사드 `nacre`가 채워져 소비자가 **한 줄**로 매단다(아래 §파사드).

**공개 표면 조사 (2026-07-26, 코드 실측 — 다시 조사하지 말 것).** 외부 소비자 관점에서 무엇이 막혀 있는지 훑은 결과.
- **이미 열려 있다(막혀 있다고 오해했던 것들):** `Model`의 모든 필드와 `Vertex/Edge/Face/Shell/Solid`의 모든 필드가 `pub`이고 `Model::reachable()`→`Reachable{vertices,edges,faces,shells}`도 공개라 **위상 순회는 밖에서 된다**(`shell.faces → face.outer.half_edges → edge.vertices → vertex.point`). 피킹용 `Tessellation{by_face,by_edge,…}`·`TessTriangle.face`·`TessOrigin`, 내부 정보 `Vertex::def`(=`VertexDef{ThreePlane|OnSeam}`)·`Model::{vertex_point, vertex_tol}`·`Model::motions`(⇒ 디버그 뷰어가 읽어야 할 것은 이미 다 읽힌다), `tessellate`·`Tessellation::to_obj`·`to_step`·`to_step_solid`·`validate`·`mass_props`도 공개(★ 자유 함수 `to_obj(&Model)`은 **칸 ㉓에서 삭제** — 아래). 플레이그라운드가 bounds를 얻으려 tessellate한 것은 불가능해서가 아니라 번거로워서였다.
- **파생 값과 에러 표면 — ✅ 2026-07-26 공개.** `nacre-props`에 `bounds`·`centroid`·`face_props`(넓이·중심·법선), `nacre-ops`에 `face_plane`, `nacre-topo`에 `Model::he_start`. §6의 거절 이유는 그 앞에 끝났다. **원칙: 값을 돌려주는 읽기 전용 질의**(위상 순수성 유지). `TessConfig`는 `tol` 하나뿐 — 면별 override는 미래.
  - **`bounds`는 곡선을 인지한다.** 원통 옆면은 솔기 정점보다 바깥으로 볼록하므로 꼭짓점 min/max는 **조용히 작은 상자**를 준다. 반지름 `r`·법선 `n̂`인 원은 축 `e` 방향으로 `±r·√(1−(n̂·e)²)`만큼 뻗는다(정확). OCCT `bounding`이 심판하되 **등호로 비교하지 않는다** — DRAWEXE는 상자를 보수적으로 부풀린다(실측 ~1e-7).
  - **`centroid`는 새 적분이 아니다.** 솔리드는 기준점에서 각 평면 면으로 뻗은 **원뿔들의 부호합**이고, 원뿔의 중심은 밑면 모양과 무관하게 꼭짓점→밑면중심의 **3/4** 지점이다. 즉 `mass_props`가 이미 계산하는 `(Aᵢ, cᵢ, n̂ᵢ)`만으로 `C = R + Σ Vᵢ·¾(cᵢ−R)/ΣVᵢ`가 나온다. 곡면은 그 논증이 깨지므로 **이름 달고 거절**하고, 그래서 `MassProps`의 필드가 아니라 별도 함수다(곡면 솔리드의 부피·넓이는 계속 살아 있어야 한다).
  - **`face_props.normal`은 `Option`이다.** 원통 면에는 하나의 법선이 없는데, **면 고르기는 모든 면을 훑는 일**이라 실패시키면 필터가 통째로 망가진다. `None`이면 자연스럽게 건너뛴다.
  - **면 고르기가 위상 명명 문제를 우회한다.** 코드-CAD는 면을 번호가 아니라 `filter(법선≈+Z).max_by(중심.z)`처럼 **생김새로** 고르고 매 실행 다시 고른다 — 저장된 참조가 없으니 상류가 바뀌어도 썩지 않는다.

**면의 스케치 좌표계 — `face_plane` (2026-07-26).** `ops::face_frame`을 공개한 것이고, **`PadOnFace`/`PocketOnFace`가 실제로 프로파일을 놓는 바로 그 프레임**이다(두 번째 유도가 아니라 같은 함수의 사영 — 갈라지면 앱이 계산한 위치와 보스가 어긋난다. 테스트가 비대칭 프로파일로 고정한다).

**공개 스케치 어휘 — `SketchFrame` + `face_sketch_frame` (2026-08-06, S9).** 내부에만 있던 `SketchFrame{plane: Handle<Surface>, placement, flip}`이 공개됐다 — 필드는 비공개, 생성자가 검증한다: `canonical(plane)`은 유도라 검사 없음, `named(model, plane, origin, ref_dir)`는 구성 시점에 정확 검사해 이름 붙은 거절을 낸다(`FrameOutsideDecimalWindow`·`OriginNotOnPlane` — scalar의 전역 `plane_residual_sign` 신설, Wide 이름은 BigInt 팔 —·`RefDirParallelToNormal`). `face_sketch_frame`은 `face_frame`이 내부에서 이미 만들던 값을 버리지 않고 공개한 이음새다(`face_plane`은 같은 프레임의 f64 사영 — 호환 유지). **★ 개정(2026-08-16): world-분기(노드 생략) 면에서는 «받아쓰고 실현으로 검증»한다** — pad 는 그 면들을 canonical 프레임이 아니라 법선 유도 세계축(`frame_axes(n)`)에 스케치하므로, 반환 후보(canonical, 그다음 pad 축의 `Named` 받아쓰기)의 실현이 `face_plane` 과 **비트 동일**할 때만 반환하고, 아무 철자도 검증을 못 통과하면 `FrameNotRepresentable` 로 이름 붙여 거절한다(그 인구 = motion 이 기록된 world-분기 면. ★ 2026-08-17 invariant-plane 재진술이 예고대로 **코드 수정 없이** 주 인구를 소멸시켰다 — 모션이 고정하는 평면은 이제 세계 진술을 유지하고 계약 스윕의 회전 핀이 2→0 실측; 잔여 = 정확-진술-가능하지만 불변 아닌 상·mirror 사슬·2세대 이동). 옛 폴백의 전제(*"canonical 이 같은 세계축으로 실현된다"*)는 flip=false·무모션에서만 참이었고, flip=true 축정렬 면 전부에서 점대칭 프레임을 반환하고 있었다(실측 — 같은 발자국이 pad 와 2.55 떨어진 자리에). 계약은 `tests/sketch_frame_contract.rs` 가 4배치 × 6면으로 잠근다. flip 측정은 `measured_frame`, 노드 push는 `push_frame_node` 한 곳으로 통일돼 extrude·face 두 도로가 한 모양이다. `Operation`이 평면 핸들을 싣는 어휘 교체는 S5(datum op)와 함께다 — replay 자기완결성: 로그 속 핸들의 합법 표적은 씨앗·기존 면·datum뿐이다. 그리고 `Model::new()`가 세계 축 평면 셋을 심는다(핸들 0·1·2, 캐시 방향 −축, `world_plane(Axis)` 접근자, `Default`는 `new()` 위임) — 세계 평면 위 스케치와 원점 상자의 축 면이 같은 surface 핸들을 공유한다. 상세·관문 실측은 `docs/truth-and-cache.md` S9 행.

**원점을 꼭짓점 평균에서 면의 *면적중심*으로 바꿨다.** 옛 규칙은 오목한 면에서 면적중심이 아니었고, 더 나쁘게는 **직선 도중에 꼭짓점이 하나 늘면 움직였다** — 면의 모양은 그대로인데 보스가 다른 자리에 앉는다. 면적중심은 **영역의 성질**이라 이산화에 무관하다. *(월드 원점 정사영(Onshape 방식)도 검토했으나 채택하지 않았다: 안정적이지만 원점에서 먼 면에 `pad`하면 프로파일이 면 밖에 앉아 대개 실패한다. Onshape는 사용자가 스케치를 모서리에 구속으로 붙이지만 스크립트엔 그 단계가 없다.)*

X축은 `any_perpendicular` — **가장 작은 성분의 축과 외적**, Onshape `perpendicularVector()`와 같은 규칙이다. 가장 작은 두 성분이 **동률일 때 불연속**이고(연속인 선택은 수학적으로 불가능), 그 타이브레이크는 `nacre-math`에서 테스트로 못 박혀 있다.

**남는 한계**: 면의 *모양 자체*가 바뀌면 면적중심도 움직인다(면에 매인 어떤 규칙도 그렇다). 그리고 프레임은 f64다 — 축이 정규화를 거치므로 **유리수 법선에 수직인 단위벡터는 일반적으로 무리수**이고, 면 위 스케치의 정확성은 원점이 아니라 **축**이 벽이다(⑦).

**디버그 뷰어는 커널 크레이트가 아니라 워크스페이스 밖 별도 앱이다.** 연산 로그를 입력받아 매 동작을 스텝별로 재생하며(append-only라 "N번째까지 replay"가 공짜), STEP에 안 담기는 nacre **내부 정보**(`Origin`의 `Constructed`/`Discovered`, `Discovered`의 tolerance 실측값, `Handle` 관계·인접 등)까지 시각화하는 인터랙티브 도구. 내부 자료구조에 접근해야 하므로 nacre를 **직접 링크**한다(개발 중 path 의존 → 안정화 후 version 의존, 버전별 디버깅도 자연스러워짐). 만드는 시점은 **`Discovered`/tolerance가 처음 등장하는 M5 즈음** — 그 전(M1~M4는 전부 `Constructed`)의 시각 확인은 정상 결과는 STEP→step-loupe(구조+검증), 중간·깨진 상태는 OBJ 덤프→맥 미리보기로 충분해, 인터랙티브 뷰어는 필요가 증명될 때까지 미룬다.

`Store`/`Handle`은 **최하위 `nacre-store`에 둔다.** geom도 Handle을 쓰기 때문이다 — `Curve::Intersection`(§3)이 `Handle<Surface>`를 담으므로, Handle이 topo에 있으면 geom→topo→geom 순환 의존이 된다. typed-index 저장소는 기하·위상을 전혀 모르는 순수 인프라이므로 두 층보다 아래에 격리하고, 위의 모든 크레이트가 자유롭게 참조한다. (라이선스는 MIT/Apache-2.0 듀얼 — Manifold(Apache-2.0) 알고리즘 차용과 호환.)

`nacre-scalar`는 **회전 오버홀의 근본 표현 — exact 유리수 스칼라**를 격리한다. 사용자가 입력한 치수·각도를 f64 오차 없이 정확히 보존한다(`1.1`→`11/10`, `1.1×7`=정확히 `7.7` — "얇은 막" 문제의 근본 해결). `Rat`은 `Ratio<i128>` + **checked 산술**으로, 오버플로가 §4 강등 **트리거**(값의 캐시를 f64/dd로 내리고 tol을 `Origin`에 기록; 정의는 op-log로 불변 보존). `Angle`은 유리수 deg를 mod-360 **정확 누적**(한 바퀴가 정확히 0으로 닫힘 → 스케치 닫힘)하고 90°계열은 exact 유리수 cos/sin(회전 tol 0). **★ `nacre-predicates`와 상보(겹침 아님):** predicates는 기하 행렬식의 **부호를 exact 결정**(exact-부호), nacre-scalar는 **입력 값과 유리수-순수 누적을 exact 보존**(exact-값) — 역할이 갈려 이름·층이 분리된다. **의존 결정(오버홀 최초 새 외부 dep):** `num-rational`(+num-traits)을 채택 — 성숙한 checked 유리수+gcd 약분을 제공하고, exact 유리수를 손수 구현하면 버그가 exactness 목표를 훼손하기 때문(MIT/Apache·순수 Rust). 어떤 `nacre-*`에도 의존 않는 **의존 그래프 최하단 순수 토대**. **회전 좌표의 toleranced 부호 판정**(무리수 좌표라 exact 못 하지만 부호는 f64 필터→astro-float 상승→*증명된 일치 또는 정직한 미결*로 sound하게 정함, §9 ③)은 이제 **`nacre-cip::kernel`으로 분리**됐다 — nacre-math 독립 순수 술어층, predicates의 쌍둥이(§9). **범위:** exact 값 엔진까지. 통합 값+tol `Scalar`·판정 잔여 정책(§9 ③)·커널 배선은 후속 셀.

`nacre-predicates`는 **발견된 교차점의 부호 판정(내/외·orientation)을 좌표가 아니라 implicit point(정의)째로 하는 indirect predicates**를 격리한다(§3 정밀도 분업, §8 M5). 바닥의 적응 정밀 확장 산술은 `geometry-predicates`(MIT/Apache) 재사용, 그 위 implicit point 표현과 indirect 술어만 자체 구현. Rust 최초의 오픈소스 indirect predicates가 되도록 **nacre 밖으로 떼어낼 수 있게**(MIT/Apache 단독 공개 가능) 설계한다. 라이선스 엄수: 구현 참고처는 논문(Attene 2020, arXiv 2105.09772; Shewchuk 1997; Lévy PCK)과 `geometry-predicates` 소스로 한정하고, LGPL인 Attene 참조 구현 소스는 **작성 중 열람 금지 / 완성 후 실행 대조만 허용**(§8 M5·§7). M1~M4는 `Constructed`만 다뤄 이 크레이트가 불필요하므로 실제 구현은 M5.

## 2. 저장소와 동일성

Fornjot에서 그대로 가져오는 부분. 수정·삭제 없는 append-only `Vec` 기반 저장소와, 인덱스+타입만 가진 Handle.

```rust
pub struct Store<T> { items: Vec<T> }

impl<T> Store<T> {
    pub fn push(&mut self, item: T) -> Handle<T> { /* index를 Handle로 */ }
    pub fn get(&self, h: Handle<T>) -> &T { /* 항상 유효 — 삭제가 없으므로 */ }
}

/// 번호표일 뿐(u32 인덱스). T는 "어느 store를 가리키는가"의 타입 라벨.
/// derive를 쓰지 않는다 — `#[derive(Hash/Eq/Ord)]`는 T에 불필요한 바운드를
/// 자동 추가하고, 그러면 `HashMap<Handle<Edge>>`(§4 Adjacency)가 `Edge: Hash`를
/// 요구해 컴파일 실패한다(Edge는 f64 tol을 품어 Hash/Eq 불가). index만 비교/해시하는
/// 수동 impl로 T 바운드를 끊는다. Copy 필수(정수 번호표), PhantomData<fn()->T>로
/// T와 무관하게 Send+Sync·공변성 확보.
pub struct Handle<T> { index: u32, _t: PhantomData<fn() -> T> }

impl<T> Clone for Handle<T> { fn clone(&self) -> Self { *self } }
impl<T> Copy for Handle<T> {}
impl<T> PartialEq for Handle<T> { fn eq(&self, o: &Self) -> bool { self.index == o.index } }
impl<T> Eq for Handle<T> {}
// Hash / PartialOrd / Ord / Debug 도 index만 보는 수동 impl.
```

동일성은 `Handle` 비교로 끝난다. "이 두 정점이 같은 점인가?"는 `h1 == h2`이며 부동소수점이 개입하지 않는다. 삭제가 없으므로 Handle은 영구히 유효하고, "지워진 객체를 가리키는 참조" 계열의 버그가 원천 차단된다. 불리언 등으로 객체가 소비되어도 항목을 지우지 않고, 결과 Solid가 새 항목들을 참조할 뿐이다(사용되지 않는 항목은 스냅샷 저장 시 mark-and-compact로만 정리).

주의: Handle의 유효성은 **자기 Model 안에서만** 성립한다. replay가 새 Model을 반환하므로 모델 두 개가 공존하는 순간이 실제로 생기고, A 모델의 Handle을 B에 쓰면 조용히 엉뚱한 객체가 나온다. 디버그 빌드에서 Handle에 model-id를 넣어 `get` 시 검사한다(릴리즈에서는 zero-cost로 제거). model-id는 프로세스 로컬 `AtomicU64` 카운터로 발급한다(replay마다 새 값). 따라서 **Handle 자체는 직렬화하지 않는다** — 영속화 대상은 연산 로그(§6)이고, replay가 인덱스를 결정적으로 재생성한다.

**로그 속 핸들은 «인덱스 어휘»다.** `Operation` 7변종 중 6개가 핸들을 싣는데(`PadOnFace`·
`PocketOnFace`·`Boolean`·`Transform`·`Mirror`·`Copy`), 그 핸들이 가리키는 셀은 로그가 기록된
모델의 것이지 replay 가 짓고 있는 모델의 것이 아니다. 핸들에서 모델을 건너 살아남는 부분은
**인덱스뿐**이므로, replay 는 op 를 적용하기 직전 그 인덱스를 **자기 아레나의 핸들로 재고정
(rebind)** 한다. 범위 밖이면 `LogHandleOutOfRange { cell, index }` — **존재만** 답하고,
live 여부·합법성은 여전히 op 자신의 `SolidNotLive`/`FaceNotInLiveSolid` 몫이다.

재고정은 **op 단위 just-in-time** 이다(op N 의 핸들은 op N−1 이 만든 셀을 가리키므로 사전
일괄 변환은 성립하지 않는다). 덕분에 원자성이 공짜로 따라온다: 재고정 실패는 그 op 가 아직
아무것도 push 하기 전에 일어나고, 실패한 replay 의 지역 모델은 통째로 버려진다.

**`apply` 는 재고정하지 않는다.** `apply` 의 모델은 호출자의 것이고 따라서 그 핸들도 호출자의
것이다. 거기서 재고정하면 외래 핸들을 **조용히 세탁**해 위 model-id 가드의 가치를 파괴한다.
재고정이 정당한 것은 「모델을 자기가 처음부터 짓는」 replay 하나뿐이다.

**유효 범위 — 오독하기 쉬운 자리.** 인덱스 어휘가 성립하는 것은 **로그를 처음부터 통째로
재생할 때**뿐이다(그때 인덱스는 재생이 스스로 만든 것이다). 로그 **중간을 수정한** 재생에는
성립하지 않는다 — 상류 수정이 하류 인덱스를 밀어낸다(v1 비목표, §6 `OpRef` 참조). 그리고
이것은 「핸들을 직렬화해도 된다」는 뜻이 **아니다**: 직렬화 대상은 여전히 로그이고, 그 안의
참조가 인덱스일 뿐 핸들이 아니다(재생이 매번 새로 발급한다).

**자기완결성 전제.** 재고정은 「로그가 그 모델의 전체 이력이다」를 **대체하지 않고 요구한다**.
로그 밖에서 셀을 넣은 모델(예: 테스트 전용 `add_cuboid`)은 replay 로 재현되지 않으며,
범위 검사는 인덱스가 **범위 안이면서 다른 셀**을 가리키는 경우를 잡지 못한다.

**거절은 live 모델에 원자적이지만 아레나 인덱스에는 아니다.** 거절에는 두 종류가 있다:

- **이른 거절** — 요청만 보고 판정(`NonPositiveDistance`, 창 밖, `Profile2d::check`,
  `NonPlanarFace`, `SolidNotLive`…). 아레나 Δ = 0.
- **늦은 거절** — 기하를 지어 봐야 알 수 있는 판정(`PadMissesFace`, `PocketNotBlind`,
  `Boolean(_)`). 도구 프리즘의 셀이 **append-only 아레나에 남는다**. 실측(2026-08-07):
  `PadMissesFace` +84 셀(v24/e36/f18/sh3/so3), `PocketNotBlind` +80(v24/e36/f16/sh2/so2),
  `Boolean(NonManifoldVertex)` +63(v19/e30/f12/sh1/so1).

live 모델은 어느 쪽이든 **거절 전 상태로 복원된다** — 커밋 후 거절은 없다. ★ 부울에서는 이것이
2026-08-13 부터 **우연이 아니라 구조**다: 조립의 피연산자 은퇴가 세 분기에 복사돼 있던 것을
「결과가 받아들여진 뒤 한 자리」로 모았고(`assemble_fuse_cut` → `reconstruct`), 그와 별도로
`boolean`/`boolean_with_classes` 가 엔진의 **모든** `Err` 에 live set 스냅샷을 복원한다. 전에는
「조립의 모든 거절이 마침 은퇴 앞에 있다」는 사실에 기대고 있었고, 그걸 단언하는 테스트는 없었다
(지금은 `coverage/rejects.rs` 와 `self_touch.rs` 가 단언한다). 그러나 아레나 길이는
되돌지 않으므로, **거절 뒤에도 기록을 이어가려면 모델을 로그로 다시 지어야 한다**(`model =
replay(&log)`). 이 규율을 어긴 세션은 replay 가 재현할 수 없는 인덱스를 적게 되고, 결과는
**이름 붙은 거절이거나 발산**이며 — 패닉은 아니다(측정: `LogHandleOutOfRange{Solid, 4}`).

이 절의 주장은 전부 `crates/nacre-ops/tests/replay.rs` 가 측정한다.

저장소는 도메인별로 나눈다:

```rust
pub struct Model {
    // 정확 기하 (진실)
    pub surfaces: Store<Surface>,
    pub curves:   Store<Curve>,
    // 위상 (기하를 Handle로 참조)
    pub vertices: Store<Vertex>,
    pub edges:    Store<Edge>,
    pub faces:    Store<Face>,
    pub shells:   Store<Shell>,
    pub solids:   Store<Solid>,
    // 파생 캐시 — 역방향 인덱스 (§4)
    pub adj:  Adjacency,
}
```

**Model은 "진실"만 담는다 — tess도 ops도 Model의 필드가 아니다.** §0에서 "모델은 연산 로그의 재생 **결과**"라 했으니 `ops: Vec<Operation>`은 Model을 *만들어내는 입력*이지 Model 안에 든 게 아니고, tessellation은 Model에서 *뽑아낸 파생 캐시*다(진실/캐시 분리). 게다가 레이어링상 Model은 `nacre-topo`에 사는데 `Tessellation`(§5)·`Operation`(§6)은 topo보다 위 크레이트라, Model에 필드로 넣으면 topo→tess/ops 순환이 된다. 그래서 op 로그와 tess 캐시는 상위 레이어(ops/파사드)에서 Model과 **나란히** 보관하고, "동일 로그·동일 cfg → 동일 모델 + 함께 재생성되는 tess"라는 §5·§6의 커플링은 그 상위 번들이 책임진다.

(점 저장소는 두지 않는다 — 두 Vertex가 한 Point를 공유하는 상황은 설계상 존재하면 안 되고, Point는 작아서 간접화의 실익이 없으므로 `Vertex`에 인라인한다.)

### 편집 연산의 supersede 의미론 — live 도달가능성 (partially persistent)

**결정 (M4에서 확정).** 편집 연산(imprint·pad, 이후 boolean)이 기존 위상을 바꿀 때, append-only라 옛 셀을 **지우지 못한다**. 그래서 Model은 **`live_solids: Vec<Handle<Solid>>`(살아있는 솔리드 목록)를 진실로 보유**하고, 편집 연산은 옛 셀을 남긴 채 새 셀을 push한 뒤 **live 목록이 새 결과 Solid만 참조하도록 갱신**한다(소비된 입력 Solid는 목록에서 빠진다). **"살아있는 모델"의 정의 = live_solids에서 하향 참조로 도달 가능한 셀의 폐포(reachable closure).** supersede된 옛 셀은 아레나에 남되 어떤 live solid도 안 가리키므로 자동으로 "안 보인다".

**★ `Mirror` — 반사는 루프 역순으로 표현한다 (2026-07-26 구현).** `Operation::Mirror { solid, axis, offset }`는 좌표평면 `axis = offset` 기준 반사다(**3개 방향 × 임의 위치**). 길이를 보존하고 **손잡이를 뒤집는다** — `Isometry`가 proper motion 전용이라 회전·이동의 어떤 조합으로도 못 만드는 유일한 강체 운동이다("−1배 스케일"은 미러와 스케일을 뭉개는 트릭이고 커널에 스케일은 없다). `transform`처럼 입력을 supersede하므로 원본을 남기려면 `copy`와 짝짓는다.

**축 정렬로 한정하는 이유는 회전과 같다** — 커널의 회전축도 X/Y/Z뿐이다. 임의 평면 반사 `x − 2((x−p)·n)n`은 무리수 법선 때문에 좌표가 근사가 되고, 그러면 회전 포리스트에 tol을 가진 반사 노드가 필요해진다. **임의 축 회전과 임의 평면 미러는 같은 항목**이며 그 설계가 들어올 때 함께 열린다.

**구현의 핵심은 방향 대수다.** 반사는 det = −1이라 `R(a) × R(b) = −R(a × b)`이므로, 좌표만 반사하면 루프의 감김이 함의하는 법선이 뒤집혀 **솔리드가 안팎으로 뒤집힌다**. 그래서 모든 루프(outer + inner)를 **역순으로 감고**(`Loop::reversed` — 순서 역전과 `forward` 반전을 함께 하므로 모서리 짝이 유지된다), `Orientation` 플래그는 **건드리지 않는다**: 반사가 내적을 보존하므로 `sign(plane.normal · n_out)`이 그대로다. 정확성은 `nacre-geom`의 `Plane::mirrored`가 `raw`(비정규화 법선)를 함께 반사해 지킨다 — ops에서 `from_point_normal`으로 다시 만들면 정확 계수가 반올림된 단위 법선으로 대체된다. **`validate`의 그물은 절반이다**: 면 단위 orientation↔감김 불일치는 `FaceMisoriented`가 잡지만(2026-08-17 신설 — 반사는 루프 역순과 `raw` 반사가 함께 가서 이 검사에 불변), 플래그와 감김이 **함께** 뒤집힌 전-셸 반전은 면 단위 검사가 원리적으로 못 보므로 그쪽 그물은 여전히 부호 있는 부피와 불리언, 그리고 OCCT 오라클이다.

**★ 반사도 모션이다 — 사슬에 들어간다 (2026-07-28 개정).** 원안은 반사를 **켤레**로 날랐다(`M ∘ R = (M R M⁻¹) ∘ M` — 축 정렬 반사와 X/Y/Z 회전축에서 켤레는 같은 축·반사된 중심·부호 반전된 각도의 회전이다). 그 방법은 정확했지만 **틀린 자리에 있었다**: 한 질문에 두 메커니즘이 생기고, 사슬이 방금 수행한 반사를 스스로 말하지 못한다. 그래서 `Constructed` 표면은 켤레할 사슬이 없어 **그 반사상이 무조건 정확하다고 선언**됐고, 미러로 도달한 벽과 이동으로 도달한 같은 벽이 1 ULP 갈라져 한 부품이 두 몸통이 됐다.

이제 `Motion::Mirror { axis, offset }`이 사슬 노드다. 기록 여부는 다른 모션과 **같은 규칙**을 따른다 — `2c − x`가 모든 좌표를 정확한 f64에 남기면(이진 미러면) 노드를 안 만들고, 그 면제는 **이력이 이미 있으면 소멸한다**. 이동에 대해 세운 규칙이 여기서 처음 실전이 된다(비-90° 회전은 가수를 다 쓰므로 이동 쪽에서는 도달 불가였다).

**★ 반사는 비고유(`det = −1`)라 지름길에 패리티 규약이 필요하다.** exact 지름길들은 *"모션이 행렬식을 보존한다"* 를 근거로 **사전-모션 프레임**에서 답한다. 반사는 그것을 뒤집고, 게다가 실패가 보수적이지 않다 — 네 점이 같은 사슬을 지니면 행렬식이 **균일하게** 뒤집혀 확신을 가진 오답이 된다. 자리마다 인자를 유도하는 대신 **base 데이터의 손잡이를 한 번 맞춘다**: 사슬의 반사 수가 홀수면 base 점들에 정확한 반사를 한 번 더 걸어 건넨다(`shared_base`, `BaseFrame::of` 두 곳뿐). 그러면 base와 moved가 **고유 모션으로 관계**되어 어떤 행렬식 질문이든 등식이 그대로 선다. **base 평면은 보정된 삼각형에서 파생한다** — 평면을 따로 보정하면 외적이 유사벡터라 전역 부호가 하나 남고, `plane_pair_dir_sign`이 바로 그 부호를 읽는다(차등 시험이 실제로 잡았다). 예외는 `cmp_coord`의 지름길 하나뿐이다: 그것은 행렬식이 아니라 **축별 비교**이고 반사는 그 축을 뒤집으므로, 반사가 든 사슬에서는 끈다(보수적 미스).

**미구현(정직 거절):** 곡면(원/원통 — 반사가 매개화 손잡이를 뒤집어 M6 결정 사항).

**★ 예외 하나 — `copy`는 supersede하지 않는다 (2026-07-26 구현).** `Operation::Copy { solid }`는 솔리드의 독립 쌍둥이를 만들고 **원본을 live 목록에 남긴다** — 즉 **live 목록에 추가만 하는 유일한 연산**이다(다른 모든 편집 연산은 "옛 것 빼고 새 것 넣기"). 필요한 이유: 그전 커널에는 *이동*만 있고 *복사*가 없어(`transform`·`boolean`이 입력을 live에서 뺀다) 같은 공구를 두 번 쓰거나 원본을 남긴 복사본을 두는 것이 불가능했고, 그래서 편의 레이어의 패턴·미러가 성립하지 않았다.

**구현:** `transform::transform_solid`(결정적 7패스 깊은 복제 walker)를 **항등 아이소메트리**로 부르고 원본을 live에서 빼지 않는다 — `transform`과의 차이가 그 한 줄의 부재뿐이다. 위상 셀은 **반드시 복제**한다(두 live solid가 모서리를 공유하면 `edge_uses`가 면 4개 사용으로 읽어 manifold 검사가 깨진다). **기하도 공유가 아니라 복제한다** — 설계 단계의 "핸들 공유" 결정을 뒤집었다: 순수 평행이동은 평면의 exact `raw`를 보존해 사본 계수가 **비트 단위로 동일**하므로 공유의 정확성 이득이 0이고, 두 live 솔리드가 `Handle<Surface>`를 공유하는 상태는 전례가 없어 위험만 남는다(공유는 나중 메모리 최적화). **★ `Discovered.tol`을 그대로 옮기는 것이 copy에서는 근사가 아니라 정확하다** — 점도 평면도 비트 동일이라 "평면에서 벗어난 거리"라는 측정값이 변할 수 없다(회전·이동에서는 stage 2의 tol 전파가 갚아야 할 빚이지만, copy는 갚을 것이 없다). 비-live 입력은 거절한다(소모된 핸들의 부활을 막고, 소비자 버그를 드러낸다). 이름은 **copy**(Handle의 clone은 얕은 참조 복사라는 반대 뜻).

**결과로, 모델을 소비하는 모든 코드는 store 전체가 아니라 live 도달가능 집합만 순회·카운트한다** — `validate`의 Euler/manifold, `Adjacency`(§4), `nacre-props`(부피·면적), `nacre-step` export, 오라클. (M1~M3엔 supersede가 없어 도달가능 집합 == store 전체였기에 store 길이로 세는 단순화가 맞았고, M4가 그 전제를 처음 깬다.) 정리(compact)는 **스냅샷/직렬화 시점에만** 하며(§2 위), 편집 중엔 안 한다 — Handle 인덱스 안정성(=replay 결정성)을 유지하기 위해.

**대안 대비 (왜 도달가능성인가).** (a) *통째 복사*(편집마다 모든 셀을 새 독립 Solid로 복제) — 안 바뀐 셀까지 복사해 **정점 공유(single-reference, §2 상단)를 깨고** 메모리·오염을 부른다. 탈락. (b) *tombstone*(무효 Handle 집합을 손으로 관리) — 표시 누락 시 조용히 오염. 도달가능성은 **"아무도 안 가리키면 죽음"이라 누락 실수가 구조적으로 불가능**. 채택. 이는 partially persistent data structure(과거 버전은 읽기만, 최신만 수정)의 표준 CAD 적용이며, §2 상단의 "소비돼도 안 지우고 결과 Solid가 새 항목 참조" 원칙의 구체화다.

**상용과 의도적으로 다른 선택 — 근거.** Parasolid·ACIS·OCCT 등 상용 커널은 40년간 "기하 불변 + **위상 가변**(Euler operator로 in-place 편집)"이 정설이다. nacre는 "기하 불변 + **위상도 불변**(append-only supersede)"으로 간다. 가변 위상은 dangling 참조·transient topology numbering·persistent naming 문제를 낳고 그걸 길들이는 데 40년 케이스워크가 필요한데 우리에겐 그 경험이 없다. 불변 위상은 그 버그 부류를 **구조적으로 제거**해 경험 없이도 안전하고, undo·결정론·fuzzing이 공짜이며, Rust 소유권 모델과도 정합적이다(가변 위상 그래프는 borrow checker와 싸워야 해 악명 높게 어렵다). 대가는 메모리인데 2020년대엔 작다. 즉 "상용을 이기는 방식"이 아니라 **"우리 제약(경험 없는 개인+AI 프로젝트·Rust·결정론 중시)에 최적화된 방식"**이다.

**성능·전제 (지금은 기록, 구현은 프로파일링 후).** 매 순회의 도달가능성 계산이 병목이면 generation/epoch 태그나 reference counting을 **파생 캐시**로 도입할 수 있다(진실=도달가능성, 이건 캐시 — Adjacency가 "진실=위상store, 캐시=역인덱스"인 것과 같은 층 분리). ref-counting은 하향 참조가 acyclic이면 죽음 판정에 바로 쓸 수 있다. **전제**: 도달가능성·ref-counting 둘 다 "Solid→Shell→Face→Loop→Edge→Vertex 부모→자식 단방향, 순환 없음"을 요구 — §4 위상 참조가 이를 만족함을 확인했고(모든 셀이 하향 Handle만 보유), 유일한 역방향 인덱스 `Adjacency`는 진실이 아니라 캐시라 무관.

**⟳ 재검토 예정 (v2) — 세션 중 메모리 관리: compact보다 재구축(rebuild-from-log) 우선.** *지금 확정·구현이 아니라, v2에서 저장·세션·undo 형식을 실제 설계할 때 꺼낼 방향. 위 "스냅샷 시 mark-and-compact로만 정리"(§2 상단·이 절)를 그때 이 방향으로 재검토한다. 로그가 진실인 우리 설계에선 compact가 주력일 필요가 없다.* 배경: append-only라 편집·브랜치를 반복하면 도달불가(dead) 셀이 누적된다.
- **정리는 compact보다 재구축이 낫다.** compact(dead 셀 물리 제거 + 인덱스 재배치)는 "현재 브랜치의 과거 단계 셀"도 도달불가라 함께 지워 **undo를 죽인다**. 재구축(현재 live 브랜치의 로그를 0부터 재생해 새 조밀 아레나 생성)은 로그에 모든 과거 단계가 있어 **undo 히스토리째로 되살아난다**. 로그가 진실이라 결정론적으로 동일 결과 보장. "dead를 골라 지우기"보다 "live를 로그로 새로 짓기"가 더 단순·안전하며 undo까지 지킨다.
- **저장 시 compact 불필요.** 스냅샷 직렬화 때 dead를 파일에 안 담으려 아레나를 compact할 이유가 없다 — 직렬화 시점에 **live 도달가능 폐포만 순회해 쓰면 된다**(아레나 불변 → 저장 후 계속 작업). M4의 "validate·adjacency·step은 live 폐포만 순회"·"STEP은 live만 export"와 동일 패턴. STEP 출력·로그 저장 모두 아레나를 안 건드리므로 어떤 저장도 compact를 트리거하지 않는다.
- **compact가 유효한 단 하나의 자리 = 메모리 위기 시 비상 회수.** 아주 길고 무거운 세션에서 dead가 실제로 메모리를 압박하면, compact로 dead·과거단계를 회수하고 현재 상태 스냅샷 + 로그를 유지한 채 가볍게 이어간다. undo는 일시 상실되나 로그가 살아 **재구축으로 복구 가능**(영구 손실 아님). 즉 compact = "undo 즉시성을 메모리와 맞바꾸는, 재구축으로 되돌릴 수 있는 스위치". 전제: **로그는 절대 버리지 않는다**(compact는 상태만 버리고 로그는 지킴).
- **역할 정리(v2 방향):** 재구축 = 주력 정리(undo 유지), live 폐포 필터링 = 저장(아레나 불변), compact = 비상 메모리 회수(undo는 재구축으로 복구). 셋은 대립이 아니라 안전망 관계.
- **지금 할 것: 없음** — M1~M5는 세션·undo·저장형식이 아직 없어 무관. 이 메모는 v2에서 저장·세션·undo·뷰어를 실제 설계할 때 꺼낸다. 지금 구현·상세확정 금지("경계는 지금, 내용은 그때").

**미래(파라메트릭 편집 v2) 기록.** 편집 대상 참조를 transient Handle 인덱스로 두면 kernel numbering에 종속돼 상류 수정에 깨진다(§6 topological naming). 정석은 **저장된 기하 cue(3D 참조점 등)로 tolerance 내 재탐색**해 참조를 의미로 resolve하는 것(HistCAD 계열). §6 `OpRef`를 지금 구현하진 않되, 참조를 "순수 인덱스"가 아니라 "기하 cue 포함"으로 확장할 자리를 로그 포맷에 남긴다.

## 3. 정확 기하 층 (`nacre-geom`)

핵심 결정: **교차 곡선을 1급 variant로, 이중 표현으로 둔다.** 상용 커널과 STEP이 공유하는 그 구조다.

```rust
pub enum Surface {
    Plane(Plane),
    Cylinder(Cylinder),
    Sphere(Sphere),
    Nurbs(NurbsSurface),
    // Cone, Torus, 회전면, 스윕면은 마일스톤 따라 추가
}

pub enum Curve {
    Line(Line),
    /// 전체 원 carrier(중심/법선/반경). 호 = Edge(circle) + 서로 다른 두 끝점,
    /// 닫힌 엣지(bounds:None) = 전체 원 — Line+Edge와 동일한 carrier-vs-trim.
    /// step-io CurveInput::Circle(끝점 동일 여부로 원/호 구분)과 정합.
    Circle(Circle),
    Nurbs(NurbsCurve),
    /// 두 곡면의 교집합 — 절차적 정의(진실) + 근사 캐시(힌트)
    Intersection {
        surfaces: [Handle<Surface>; 2],
        /// 표시·초기값용 근사 스플라인. 진실이 아니라 캐시.
        cache: NurbsCurve,
        /// 캐시가 진짜 교차에서 벗어난 최대 거리 (SSI 행진 시 산출)
        cache_err: f64,
    },
}
```

`Intersection`의 평가는 두 단계다: cache에서 대략적 점을 얻고, 두 곡면 위로 Newton relaxation하여 머신 정밀도로 수렴시킨다. 이 함수 하나가 "언제든 정확한 점을 되찾는" 능력의 전부다.

```rust
/// cache 위 파라미터 t의 점을 실제 교차점으로 정련한다 — 적응 정밀도 사다리.
/// 반환: 정련된 점과 "실제 달성한 정확도" (이 값이 정점 캐시의 측정 tol 이 된다).
pub fn relax_to_intersection(
    p0: Point<3>, s1: &Surface, s2: &Surface,
) -> Result<(Point<3>, f64), RelaxError>;
// 1단: f64 Newton — 횡단(정상) 교차의 99%. 머신 정밀도 근처까지 수렴.
// 2단: 수렴 정체 감지(잔차가 √ε_f64 ≈ 1e-8 근처에서 멈춤 = 접선/스침 의심) 시
//      double-double(twofloat/qd 크레이트)로 재시도. √ε_dd ≈ 1e-16이므로
//      악조건에서도 f64 저장 한계까지 꽉 채운 좌표를 얻는다.
// 3단: 그래도 판정 불가 → RelaxError::Tangential로 분류해 위로 보고.
//      좌표가 아니라 위상 정책(접촉 처리, 기호적 섭동 등)이 필요한 구역.

pub enum RelaxError {
    Tangential { best: Point<3>, residual: f64 }, // 접선/퇴화 — 정책 결정 필요
    Diverged,                                     // 초기값 불량 — cache 재행진 대상
}
```

`RelaxError`는 삼키지 않고 위로 올린다. 3층 실패는 여기서 시작되므로, 실패를 조용히 넘기는 순간 디버깅 불가능한 커널이 된다.

정밀도 분업 원칙: **판정(부호)에는 적응 정밀 술어**(`geometry-predicates` — 쉬운 케이스 f64, 아슬아슬할 때만 확장, 부호는 항상 정확), **반복 구성(좌표)에는 고정폭 확장 부동소수점**(double-double). 임의 정밀 유리수(rug/malachite)는 술어에는 완벽하지만 Newton 반복에 넣으면 비트 길이가 반복마다 폭발하므로 구성에는 쓰지 않는다. 이 투자는 평균이 아니라 꼬리를 산다 — 호출 빈도가 낮은 악조건 케이스만 정확히 개선되고, tolerance가 "상수"가 아니라 "실측 보증값"이 된다. (**유리수 경계 규칙**(회전 오버홀 §9): 유리수는 **입력 표현**일 뿐 계산 매체가 아니다 — 입력·유리수-순수 파생에만 살고, 불리언·반복·혼합이 닿는 순간 f64/double-double로 강등해 이 비트 폭발 경고와 정합한다.)

판정 술어는 다시 둘로 나뉜다: **명시 좌표점(f64 격자 위 — `Constructed`)은 direct 술어**(orient3d 등, 좌표를 직접 받음)로, **교차로 생긴 점(`Discovered`)은 indirect 술어**(implicit point = "어느 원시 요소들의 교차인지"라는 정의를 받아 좌표를 만들지 않고 부호를 정확히 계산 — Attene 2020)로 판정한다. 이유: 교차점을 f64 좌표로 만드는 순간 오차가 끼고, "정확한 direct 술어 × 부정확한 입력 = 부정확한 답"이 되기 때문. indirect 술어는 그 구멍을 닫는다(정의째 받으므로 부정확한 중간 좌표가 없음). 둘은 같은 확장 산술 바닥(`geometry-predicates`)을 공유하며, indirect 층은 `nacre-predicates`가 그 위에 쌓는다(§1). indirect predicate가 작동하려면 점이 "정의를 보유"해야 하므로, `Discovered`가 정의를 갖는다는 결정(§4)이 그 전제다 — 판정이 필요하면 indirect 술어, 좌표가 필요하면 relaxation으로 정의에서 뽑는다(둘은 대체가 아니라 역할 분담).

STEP 연결: 이 variant는 `surface_curve`/`intersection_curve`(3D curve + 곡면별 pcurve + master 지정)와 1:1로 대응한다. pcurve가 필요해지는 시점(트리밍 구현 시)에 `pcurves: [NurbsCurve2d; 2]`를 같은 variant에 추가한다. 기존 stepio 코드젠에서 이 엔티티 묶음의 보존·복원을 우선 검증 대상으로 삼는다.

## 4. 위상 층 (`nacre-topo`) — tolerance를 타입에 새기기

이 설계에서 가장 의견이 들어간 결정. 정점·엣지의 "출신"을 enum으로 구분한다.

```rust
pub enum Origin {
    /// 구성 시점에 정의됨 (스케치 점, 스윕 결과 등).
    /// 동일성은 Handle로 완결 — tolerance 개념이 없다. 여기선 point가 곧 진실.
    Constructed,
    /// 교차 계산으로 발견됨. **정의(implicit point — "어느 원시 요소들의
    /// 교차인지")를 보유하며, 그 정의가 진실이고 Vertex.point 좌표는 캐시다**
    /// (진실/캐시 분리에 점이 뒤늦게 합류 — §3). 부호 판정은 좌표가 아니라
    /// 이 정의를 indirect 술어에 넣어 얻고(§3, §8 M5), 좌표가 필요할 때만
    /// 정의에서 relaxation으로 뽑는다. 국소 tolerance도 가진다: tol은 임의 상수가
    /// 아니라 relax_to_intersection이 반환한 실측 달성 정확도에서 온다.
    /// 정의의 구체 형태 = **위상 계층 enum `VertexDef`**(M5-prep 확정). M5 다면체
    /// 꼭짓점은 평면 3장 교차라 초기 형태 `ThreePlane([Handle<Surface>; 3])`.
    /// 부호 판정 시 geom이 Handle→평면 계수를 뽑아 `nacre-predicates`(순수 수치,
    /// 계수만 받음 — §9)에 넘긴다. 실제 `definition` 필드 추가는 Discovered 점이
    /// 처음 생기는 M5-c(PolyhedralBoolean) — 그전엔 원칙만 고정.
    Discovered { tol: f64 /*, definition: VertexDef (M5-c) */ },
}

pub struct Vertex {
    /// f64 좌표. Constructed면 이게 진실, Discovered면 정의에서 뽑은 캐시 (§2·§3).
    pub point: Point<3>,    // 인라인 — 점 저장소 없음 (§2)
    pub origin: Origin,
}
```

**★ 정확 기하는 정점만이 아니라 면에도 적용된다 (`SurfaceDef`).** 위 원칙은 오래 **정점에만** 지켜졌고 서피스는 아무 말도 하지 않았다. 그런데 비-90° 회전은 평면의 계수를 **무리수로** 만든다 — 저장된 `Plane{origin, normal, raw}`은 그 순간 진실이 아니라 **반올림된 상**이 되는데, 그 사실을 적을 곳이 없으니 커널은 기본값으로 "정확하다"고 답했다. 그것이 거짓이 되는 지점이 정확히 문제가 되는 지점이다: 불리언 **결과**는 회전 이력을 하나도 안 들고 다니므로(결과 정점은 전부 `Discovered`), 그 결과의 모든 면이 **반올림된 삼각형을 정확하다고 선언**하고, 한 벽의 두 사본이 ~1e-16 어긋난 채 *확정적으로* "공면 아님" 판정을 받아 **한 평면이 두 클래스**가 됐다. 그 아래로 가짜 교선·1e15 거리의 정점·반대칭을 잃은 `order_along`이 줄줄이 따라왔다(핀 배열 대장 8/90).

```rust
pub enum SurfaceDef {
    /// 계수가 곧 진실 — 구성된 표면, 또는 정확성을 보존하는 운동(이동·90°계열)의 상.
    Constructed,
    /// 모션 이력의 상. 진실은 `(witness, motion)`이고 계수는 캐시다.
    /// `witness`는 **회전 이전**의 비공선 세 점(그 면에서 그 자리에 붙잡는다 —
    /// 정의가 모델의 다른 부분이 살아남는 데 의존하면 안 된다),
    /// `motion`은 `Model::motions` 포리스트의 leaf — 모션은 면이 든다(정점은 자기 평면을 따라간다).
    Moved { witness: [Point3; 3], motion: Handle<MotionNode> },
    /// **정확히 기술할 수 없다** — 오늘은 `push_surface`를 우회해 출처가 기록되지 않은 표면뿐이다.
    /// (S6b 에서 소멸 — 정확한 형태가 없는 표면은 표현 불가능해졌고,
    /// f64 폴백은 이름 붙은 거절이 됐다.)
}

pub struct Edge {
    pub curve: Handle<Curve>,
    /// 끝점 정점 둘. **닫힌 솔리드의 원형 rim은 seam 정점을 써 `[v, v]`(start == end)로
    /// 둔다** — 그래야 b-rep이 유효 CW-복합체(V−E+F=2)로 남아 validate의 오일러-푸앵카레를
    /// 통과한다(실린더 rim이 대표 사례). ★ 설계 초안의 `bounds: Option<…>`(끝점 없는 독립
    /// 전체 원 = `None`)은 구현되지 않았다 — 와이어프레임/열린 면 요소는 v1 비목표(§9)라
    /// 타입이 그 상태를 표현하지 않는다(2026-09-02 정정: 코드의 필드명은 `vertices`).
    pub vertices: [Handle<Vertex>; 2],
    pub origin: Origin,
}

pub struct Face {
    pub surface: Handle<Surface>,
    pub outer: Loop,            // half-edge 순환
    pub inner: Vec<Loop>,       // 구멍
    pub orientation: Orientation,
}
```

`SurfaceDef`는 곁표(`Model::surface_defs`)에 산다 — `Surface`는 `nacre-geom` 타입이라 `Handle<MotionNode>`를 이름 붙일 수 없다. 입구는 `Model::push_surface` 하나이고, **강제는 privacy가 아니라 불변식으로** 한다: `validate`가 *"live 모델의 모든 면은 정의를 가진 표면 위에 있다"* 를 검사하므로(디버그에서 매 연산 뒤 도는 검사) 기록을 빠뜨리면 즉시 터진다. `Reachable`이 surfaces를 추적하지 않으므로 순회는 **면을 통해** 돈다.

**★ 서피스는 솔리드의 leaf가 아니라 자기 leaf에서 이어진다.** 한 솔리드에 회전 이력이 하나라는 전제가 틀렸다 — 서로 다른 각도로 회전한 두 피연산자의 불리언 결과는 벽마다 다른 이력을 갖는다. 그리고 결과의 정점은 전부 `Discovered`라 정점에게 "이 솔리드는 어느 회전에 있나"를 물으면 `None`이 나오고, 다음 회전이 **뿌리부터 다시 시작**해 회전 이전 witness에 두 번째 회전만 태우게 된다(존재하지 않는 평면). 그래서 `transform`은 변환마다 노드 하나가 아니라 **서로 다른 parent leaf마다 노드 하나**를 만든다.

**따라오는 단순화.** 서피스가 스스로 말하게 되자 그 침묵을 메우려 있던 코드가 통째로 없어졌다 — 면의 평면을 정점 provenance에서 사냥하던 `face_plane_witness`·`plane_pts`, 솔리드 단위 추측 `solid_is_rotated`, 그리고 그 사냥이 실패했을 때의 거절 `RotatedUnderdetermined`. 회전 뒤 이동은 그 뒤 **모션 포리스트**가 이름을 갖게 되어 거절이 아니라 결과가 됐다(아래).

효과가 세 가지다. 첫째, 케이스 A(만나는 자리를 아는 구성 연산)만 쓰는 한 모델 전체가 `Constructed`로만 이루어지고, tolerance 코드 경로가 아예 실행되지 않는다 — Fornjot의 "무-tolerance 단순함"이 그 안전지대 안에서 그대로 재현된다. 둘째, `Discovered`가 처음 등장하는 지점이 곧 3층 코드가 개입한 지점이므로, 디버깅·검증에서 "어디부터 위험한가"를 데이터가 스스로 말해준다. 셋째, 검증 규칙을 출신별로 다르게 걸 수 있다(`Constructed`는 정확 일치 요구, `Discovered`는 tol 이내 요구).

**정의가 진실 → indirect predicate·봉합 우회(무게중심 이동).** `Discovered`가 좌표가 아니라 정의(implicit point)를 진실로 삼는다는 것은, 부호 판정을 좌표를 만들지 않고 indirect 술어로 정확히 하기 위한 전제다(§3, §8 M5; Attene 2020). 이 결정이 실제로 해소하는 것은 과장 없이 셋: (1) **조합적 결정의 불일치**(orient 부호가 여기선 +, 저기선 −로 모순돼 위상이 깨지는 것 — 불리언 실패의 주원인), (2) "정확한 술어에 부정확한 입력" 구멍, (3) 판정에서의 tolerance 튜닝 취약성(부호 판정엔 tolerance 값 자체가 안 쓰임). **봉합 문제도 대부분 해소된다**: 세 평면이 만나는 꼭짓점을 하나의 `Vertex`로 만들고 정의를 "P₁∩P₂∩P₃"로 두면, 세 엣지가 그 하나의 Vertex Handle을 공유하므로 "세 곡선이 한 점에서 만나는가?"가 좌표 비교가 아니라 Handle 비교(`h==h==h`)가 되고 — 봉합의 본체(어긋난 세 후보를 tolerance로 화해)가 사라진다. 좌표를 f64로 뽑을 때도 하나의 Vertex 좌표만 뽑으므로 세 어긋난 좌표가 안 생긴다. 잔여물은 "그 하나의 f64 좌표가 세 평면 위에 정확히 안 놓인다"는 반 ulp 오차뿐(f64 저장의 바닥, 대부분 무해, 극단적으로 조밀한 형상에서만 문제). **해소 못 하는 것**: 접선·퇴화 교차의 위상 정책(부호 0일 때 "접하냐/스치냐"는 판정이 아니라 정책 — `RelaxError::Tangential`로 올리는 구역 그대로), 좌표 정밀도 자체, 공간 인덱스·SSI·성능 등 다른 축. 곡면 적용도(M6~7)는 열어둔다(§9).

half-edge는 교과서 구조를 따르되 Fornjot 신설계처럼 공유 `Edge` + 방향 참조로 둔다:

```rust
pub struct HalfEdge { pub edge: Handle<Edge>, pub forward: bool }
pub struct Loop { pub half_edges: Vec<HalfEdge> }

/// 닫힌 면 집합 — 하나의 경계 곡면을 이룬다.
pub struct Shell { pub faces: Vec<Handle<Face>> }

/// 솔리드 = 바깥 껍질 + 내부 공동(cavity)들.
pub struct Solid {
    pub outer: Handle<Shell>,
    pub cavities: Vec<Handle<Shell>>,
}

/// 역방향 인덱스 — 원본이 아니라 파생물 (replay/재스캔으로 언제든 재구성 가능).
/// 검증("모든 엣지는 정확히 두 면에서 반대 방향으로 사용")과 인접 탐색이 이걸 쓴다.
pub struct Adjacency {
    pub edge_uses: HashMap<Handle<Edge>, SmallVec<[(Handle<Face>, bool); 2]>>,
    pub vertex_edges: HashMap<Handle<Vertex>, SmallVec<[Handle<Edge>; 4]>>,
}
```

(구현은 M1에서 `SmallVec` 대신 평범한 `Vec`으로 시작한다 — 의존성 하나를 아끼고, 인라인 저장 이득은 프로파일링으로 정당화되면 나중에 도입한다. 다양체 케이스를 힙 없이 담는 최적화라 정확성과 무관.)

`Adjacency`는 진실이 아니라 캐시라는 점이 중요하다 — 위상 store들이 진실이고, 인덱스는 연산과 함께 증분 갱신되며 불일치 시 재구성한다(M1은 from-scratch 재구축만, 증분 갱신은 ops 로그가 생기는 M2). (진실/캐시 분리 패턴의 네 번째 반복.)

**Adjacency는 store 전체가 아니라 live 도달가능 셀만 인덱싱한다(§2 supersede 의미론).** `rebuild`는 `model.faces`/`model.edges` 전체가 아니라 `live_solids`에서 도달 가능한 면·엣지만 순회한다 — 그래야 supersede된 옛 면의 엣지가 `edge_uses`를 오염시켜 manifold 검사(엣지 정확히 2회)를 깨뜨리지 않는다. 위상 참조가 하향 단방향·비순환이라(§2 전제) 도달가능성 순회는 유한·안전하다. `Adjacency`가 유일한 역방향 인덱스이지만 캐시이므로 도달가능성 진실을 침해하지 않는다.

## 5. Tessellation 층 (`nacre-tess`) — 출처 태그 파생물

메시는 진실이 아니지만 일회용도 아니다. 연산마다 증분 갱신되며, 모든 요소가 출처를 안다.

```rust
pub enum TessOrigin {
    OnVertex(Handle<Vertex>),
    OnEdge { edge: Handle<Edge>, t: f64 },       // 곡선 파라미터
    OnFace { face: Handle<Face>, uv: [f64; 2] }, // 곡면 파라미터
}

pub struct TessVertex { pub pos: Point<3>, pub origin: TessOrigin }

pub struct TessTriangle {
    pub vertices: [Handle<TessVertex>; 3],
    pub face: Handle<Face>,   // 출처 태그 — 하이브리드 불리언의 핵심 재료
}

/// tess 컨테이너. 면/엣지별 버킷이 있어야 증분 갱신이 가능하다
/// (연산이 면 몇 개만 바꿨을 때 그 버킷만 무효화·재생성).
pub struct Tessellation {
    pub vertices:  Store<TessVertex>,
    pub triangles: Store<TessTriangle>,
    pub by_face:   HashMap<Handle<Face>, Vec<Handle<TessTriangle>>>,
    pub by_edge:   HashMap<Handle<Edge>, Vec<Handle<TessVertex>>>, // 엣지 polyline
    pub stale:     HashSet<Handle<Face>>, // 무효화 표시 — 접근 시 또는 명시적 flush 시 재생성
}
```

이 층의 용도는 셋이다. 뷰어가 이걸 그대로 그린다(별도 tessellation 경로 없음). STL/3MF 내보내기가 이걸 그대로 쓴다. 그리고 마일스톤 후반에서, 강건 메시 불리언의 입력이 되어 "조합적 결정 → 출처 태그로 정확 기하 스냅백" 파이프라인의 앞단이 된다. uv/t 파라미터를 들고 있으므로 스냅백 시 Newton 초기값이 공짜로 나온다.

참고(성능): 불리언·교차 전에 면쌍 후보를 AABB/BVH로 컬링하는 공간 인덱스가 필수 과제다 — tess의 면별 버킷에서 AABB가 거의 공짜로 나온다. (Truck은 이게 없어 전수 비교에 가깝고, 저자도 "BSP 등 최적화는 미래 과제"로 명시 — 반면교사.)

**★ 삼각분할 알고리즘 — 브리징 없는 스윕 (2026-07-28 개정).** 원래는 구멍을 외곽 링에 **브리지**(폭 0 슬릿)로 꿰맨 뒤 귀 자르기(ear clipping)를 했다. 그 브리지가 정점을 **반복**시키고, 반복된 정점이 있는 링은 단순 다각형이 아니며, **Meisters의 two-ears 정리는 단순 다각형에만 성립**한다 — 즉 그렇게 만든 링에는 **귀가 아예 없을 수 있다**. 실측: 무작위 rectilinear 면(구멍 1~4개)의 **27.6%**가 메시되지 않았고(구멍 1개 0% / 2개 15% / 3개 39% / 4개 57% — 구멍이 하나면 브리지도 하나라 핀치가 구조적으로 불가능), 그 위에 *조용한 오답*의 여지까지 있었다(귀 판정이 정점 인덱스로 제외해서 슬릿 간선이 삼각형을 관통해도 못 봄 — 링이 자기교차한 뒤에도 삼각형 2개를 더 뱉는 것을 관측).

⇒ **y-단조 분해 + 단조 삼각분할**(de Berg §3). 스윕이 **링을 절대 합치지 않으므로** 퇴화 링 자체가 안 생긴다 — 구멍 간선이 보통 간선이고, 구멍의 최상단이 split·최하단이 merge 정점이 되어 대각선이 자동으로 구멍을 잇는다. 판정은 전부 부호이고 **tolerance가 없다**: 전부 exact `orient2d`(`nacre-predicates`)이고, 이것은 **어떤 f64 입력에도 정확**하다. (평면의 차트는 축 드롭이라 `uv`가 원본 f64 그대로지만, 원통의 차트는 `atan2`로 **계산된** 값이다 — 그래도 술어가 정확하므로 알고리즘은 같다. 면 사이에서 `uv`를 비교하는 일이 없으므로 — 면끼리 만나는 자리는 **엣지 폴리라인**이 이미 정한다 — 차트가 면마다 달라도 무방하다.) 정점을 추가하지 않으므로 삼각형 수는 `V + 2H − 2`로 불변이고 crack-free 규칙도 그대로다.

**★★★★ 곡면도 같은 길을 탄다 — 면마다 «차트» 하나 (2026-08-25, 칸 ⑥).** 원통 면은 손으로 쓴
특수 경로(두 림을 θ로 병합하는 quad 스트립)를 썼고, 그래서 `face.inner`를 **아예 못 읽었다** —
구멍 있는 곡면(불리언의 곡면 병합이 만드는 밴드+노치)을 그릴 수 없었다. 이제 **모든 면이 같은
길**을 탄다: 경계 루프를 **그 곡면의 차트**에 올리고 위의 단조 스윕 + 플립을 돌린다.
분기는 하나이고 **고르는 것은 알고리즘이 아니라 차트**다 — 평면은 축 드롭, 원통은 `(z, r·θ)`.
원뿔·구·NURBS는 **차트 하나씩** 더하면 된다(M6-3/M7이 특수 경로를 곱하지 않는 이유).
- **`θ`가 아니라 `r·θ`**: 원통은 가전개면이라 이 스케일이 **등거리 사상**이고, 그래야 플립이
  「그림 속」이 아니라 **곡면 위에서** 메쉬를 개선한다.
- **손잡이 보정은 교환이 아니라 `u` 부호 뒤집기**: `v`가 스윕 축이므로 교환하면 스윕이 축을 따라
  가고, 대각선이 넓은 호를 가로질러 **곡면을 벗어나는 현**이 된다. 실측(반증된 보정을 넣어): 현이
  32°를 가로질러 sag가 예산의 254배 — 그런데 watertight도 삼각형 수도 **그대로 통과**한다.
  그래서 「삼각형은 곡면을 벗어나지 않는다」가 **별도 잠금**으로 서 있다.
- **θ는 루프를 따라 풀어서** 쓴다(절대각이 아니라). 밴드의 경계는 솔기를 **두 번** 걷고, 그 두 번이
  차트에서 θ=0과 θ=2π가 되어 펼친 밴드가 **직사각형**이 된다. 구멍은 온전한 회전만큼 바깥 링의
  창으로 옮긴다.
- **삼각형 수는 유도대로**: 밴드의 차트 다각형은 정점 `2n+2` ⇒ `V+2H−2 = 2n`, 옛 walk와 같다
  (기존 `4n−4` 잠금이 그대로 초록).

**★★★★★ 모든 간선은 예산 안의 현이고, «차트»가 그것을 가능하게 하는 점을 댄다 (2026-08-25,
칸 ⑧).** `TessConfig`의 두 예산은 경계에만 집행돼 있었고(§ 아래 「예산 둘」), 면의 **내부**는 자기
경계가 우연히 대 준 샘플을 물려받았다. 경계가 곡률 방향을 아크로 덮는 동안은 공짜지만, 덮지
않으면 **조용히 틀린다**: 병합된 옆면(지운 가짜 이음매가 한 구간의 유일한 아크였던 면)이 반 바퀴를
한 번에 건너는 현으로 그려져 **넓이의 1/6**을 잃었다.

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
  복제다(실측: 남기면 보통 밴드의 삼각형 수가 움직인다). 평면 차트는 아무것도 대지 않는다.

삽입은 후보마다 **놓을 자리**를 찾는다 — 삼각형의 **엄격한 내부**(셋으로 쪼갬)이거나 **내부 간선
위**(둘로 쪼갬). 두 번째가 예외가 아니라 **주된 길**이다: 후보는 면의 모양이 바뀌는 선 위에 앉고
삼각분할은 이미 그 선을 따라 간선을 갖고 있다. ★ **그러나 두 길 다 프로덕션에서 산다.** 다섯 개의
병합된 옆면 실측: `−x`·`+x`·`+y` 벽과 코너는 내부 **0** / 간선 180·180·180·270, 그런데 **`−y` 벽만
내부 102 / 간선 76**이다 — `ref_dir`이 θ=0을 노치 **안**에 두는 유일한 벽이라 노치가 다리로
이어지지 않고 **진짜 구멍**으로 남고, 구멍의 이웃에는 후보가 앉을 긴 대각선이 없다. 어느 갈래도
도달 불가로 지울 수 없다 — 어느 쪽이 도는지는 사용자가 그린 모양이 아니라 **솔기의 우연**이 정한다.
**제약 간선은 절대 안 쪼갠다** — 이웃이 모르는 점은
T-정점이고 그건 크랙이다. 자리가 없는 후보는 **버린다**(강제하지 않는다).
`ChartMap`이 되돌리는 길이다(차트는 전단사다) — 발명된 점은 위치를 읽어 올 모델 정점이 없고, 두 팔
모두 정방향이 건 **손잡이 보정을 되돌려야** 한다.

☑ **반증 기록 둘.** (a) 「예산을 넘는 간선을 차트 중점에서 쪼개고 라운드를 돈다」는 **발산한다**
(179 → 386 → 856 → 1952 → 4412; 플립을 빼도 179 → 526 → 1520 → 4180 → 10148). 이탈이 **오직 Δv**에만
의존하는데 영역이 극단적으로 이방적이라, 어디서 쪼개도 새 간선이 비슷한 Δv를 갖는다 — 필요한 것은
v 방향으로만 잘게 나뉜 «띠»이고 그건 **점을 미리 대야** 얻어진다. (b) Shewchuk의 **지름원
encroachment**를 「침범하는 후보를 버린다」로 뒤집어 쓰는 것: 직선 경계는 폴리라인이 두 점뿐이라
세그먼트가 통째로 길고, **길이 4짜리 솔기의 지름원(반지름 2)이 한 면의 후보 360개를 전부** 삼켰다.
문헌의 규칙은 「세그먼트를 쪼갤 수 있어서 그 길이가 곧 해상도인」 알고리즘의 것이다. 자리 없는
후보는 삽입 탐색 자신이 이미 거절하며, 실측상 **같은 집합**이다.

☑ **칸 ⑥의 「Ruppert/Chew 정제는 필요 없음」은 «품질»에 관한 실측이었다**(범위 정정): 그때 잰 것은
실제 인구에서 삼각형 변의 최대 각폭이 정확히 샘플 간격 2.00°라는 것인데, 그 인구에는 **경계가
곡률을 덮지 않는 면이 없었다**. 위 결함이 그 면이다. 각도 개선(quality refinement)은 여전히 이 층의
빠뜨린 단계가 아니지만, **내부 점 공급**은 빠져 있었다.

**품질은 직교하는 별도 패스다.** 단조 삼각분할은 정확하지만 슬리버가 많다. 자유(비제약) 간선에 **Lawson 플립**을 돌려 제약 Delaunay로 보내면 최소각이 최대화된다(정리 — Lawson 1977, Chew 1989). 새 의존성 0(`incircle`이 이미 있다). 경계는 제약이라 안 뒤집히므로 watertight·개수·면적 전부 불변. **실측: 평균 최소각 2배**(원통 캡 tol 1e-4에서 0.51°→1.05°). **다만 캡의 *최악*값은 안 움직이고, 그건 알고리즘 탓이 아니다** — 정n각형 정점은 **공원**이라 `incircle`이 전부 0을 답하고(모든 삼각분할이 Delaunay), 게다가 내각 178.9°와 two-ears가 0.57° 이하 슬리버를 **강제**한다(측정 0.457° = 하한). 고치려면 점을 **추가**해야 하는데, 그건 *캡의* 얘기다 — 캡의 경계는 이미 곡률을 다 샘플하므로 더 넣을 이유가 **품질뿐**이고, 품질은 이 층의 결정이 아니다. (금지되는 것은 «점 추가» 자체가 아니라 **경계** 정점의 추가다 — 위 「차트가 점을 댄다」를 보라.)

**★ 곡선을 얼마나 쪼갤지는 예산 «둘»이 정한다 (2026-08-19, T1).** 원의 분할 수는
`n = max(⌈360°/Δθ⌉, ⌈π/acos(1 − tol/r)⌉)`이다. 뒤의 항(절대 sagitta)만 쓰면 `n ≈ π√(r/2tol)`이라
**작은 원일수록 적게 쪼갠다** — 절대 오차 예산은 지키지만 반지름 0.2가 10각형이 되고, 확대하면 그대로
드러난다. 앞의 항(각도)은 `s/r = 1 − cos(π/n)`에서 반지름이 사라지는 **상대** 예산이라 크기와 무관하게
같은 품질을 준다. 정n각형에서 **중심각 = 선분의 꺾임각 = 인접 면 법선의 벌어짐**이 모두 같으므로 이 한
숫자가 윤곽과 음영을 동시에 묶는다(OCCT의 angular deflection과 같은 개념). 기본값 `Δθ = 2°`(=180분할,
반지름의 0.0152%), `tol = 1e-2`. 각도 항은 반지름 ≈66까지 이기고 그 위는 sagitta 항이 이어받는다 —
**어느 반지름에서도 한쪽만 쓸 때보다 나빠지지 않는다.**

**틈 없음(crack-free) 규칙:** 면의 삼각분할은 반드시 `by_edge`의 공유 polyline 정점들을 자기 경계로 소비해야 한다. 인접한 두 면이 각자 독립적으로 엣지를 샘플링하면 공유 엣지에서 정점이 어긋나 T-junction이 생기고 watertight가 깨진다 — 엣지 polyline이 먼저, 면 삼각분할이 그걸 경계 조건으로. (validate가 이를 검사한다: §7.)
★ **금지되는 것은 «경계» 정점의 추가다.** 면이 자기 **내부**에 점을 넣는 것은 이웃이 볼 일이 없으므로 크랙과 무관하고, 위 「차트가 점을 댄다」가 바로 그것이다. 폴리라인을 더 잘게 하는 것도 원리상 금지가 아니다 — `by_edge`가 공유 구조이므로 **엣지 샘플러가** 촘촘히 내면 양쪽 면이 함께 본다(OCCT도 `BRepMesh_ModelHealer`에서 그렇게 한다). 금지는 **면이 혼자** 경계를 바꾸는 것이다.

갱신 규칙: 연산이 위상에 항목을 추가하면 같은 트랜잭션에서 tess에도 해당 항목을 추가한다(Fornjot의 "함께 쌓기"). tolerance를 바꾼 재계산은 tess만 통째로 재생성하고 위상·기하는 불변.

## 6. 연산 층 (`nacre-ops`) — 재계산 가능한 로그

```rust
pub enum Operation {
    // M5(S5(i)-a): 모델이 **호출자가 이름 붙인 평면**을 드는 유일한 연산. 로그가 이름 부를 수
    // 있던 평면은 씨앗 셋과 «이미 만든 면의 평면» 둘뿐이었고, 그 밖은 전부 Extrude 안에 값으로
    // 살았다. 평면을 로그 **밖**에서 만들면 그 모델은 자기완결적이지 않으므로 연산이어야 한다.
    // 평면은 구성 시점 interning 대상이라 «이미 있으면 그 핸들»이 정답이다(아레나 안 자람).
    DatumPlane { def: DatumDef },  // Stated(SketchPlane) | Offset { frame, dist } | ThroughVertices([Handle<Vertex>;3])
    // M2: 스케치 평면·프로파일 개념이 Extrude 인자로 흡수된다(별도 Sketch op 없음). Extrude는
    // 프로파일에서 새 솔리드를 만들므로 앞선 op의 면을 참조할 필연이 없다 — op간 Handle 참조는
    // 기존 면 위에 작업하는 M4(PadOnFace/PocketOnFace)에서 비로소 필연적으로 도입된다.
    // M2는 매 Extrude 결과가 닫힌 솔리드라 validate가 빈틈없이 걸린다.
    Extrude { frame: SketchFrame, profile: Profile2d, dist: f64 },   // ← 평면을 «이름 부른다»(S5(i)-b)
    // ★ 칸 ⑨(2026-09-05): 원통 원시체 `Cylinder { frame, center, radius, dist }`는 **은퇴했다** — 원이
    // 스케치의 어휘(`Edge2d::Arc`, `start == end`)가 되어 `Extrude`가 나른다. 원 프로파일의 Extrude 가
    // 원시체와 **위치 정준 비트 동일**(정점 비트·솔기 정의·평면/원통 진술·면·간선·부피 비트·삼각형)임을
    // 잰 뒤 걷었다. 상자가 «네 점 + 돌출»인 것과 같은 구조. topo 의 `add_cylinder_exact`는 test-util
    // 문으로 남아 그 잠금의 대조군이다.
    Revolve { plane: SketchPlane, profile: Profile2d, axis: Axis, angle: f64 },
    PadOnFace { face: Handle<Face>, profile: Profile2d, dist: f64 }, // M4: 기존 면 위 — Handle<Face> 참조
    Boolean { kind: BoolKind, a: Handle<Solid>, b: Handle<Solid> },
    // ...
}

**연산은 평면을 «이름 부른다» (S5(i)-b).** `Extrude` 가 `SketchFrame`(평면 핸들 + 배치 +
측정된 `flip`)을 받으면서, 평면을 **값**으로 싣는 변종은 사라졌다. 그러면 원칙 2 가 연산 어휘
전체에서 성립한다 — 평면은 아레나에 있고 연산은 그것을 가리킨다.

| 무엇 | 어디서 |
|---|---|
| 세계 평면 위 스케치 | `SketchFrame::world(&m, Axis)` — 씨앗을 이름 부른다 |
| 기존 면 위 스케치 | `face_sketch_frame(&m, face)` |
| 그 밖의 평면 | `Operation::DatumPlane` 로 **먼저 진술**하고 돌려받은 프레임을 쓴다 |

**방향은 프레임의 것이고, `flip`은 재는 값이다.** `dist > 0`은 두께이고 어디로 가는지는 프레임의
ŵ이다. `flip`은 소비 연산이 `measured_frame`에서 **한 곳에서만** 재므로(S9) 호출자가 뒤집힌 프레임을
진술할 수 없다 — 면의 프레임은 그 면의 바깥을 향한다. 그래서 «면에 그려 안쪽으로» 파는 것은
`toward`를 안쪽으로 재는 **복합 연산**(`PadOnFace`/`PocketOnFace`가 `extrude_and_boolean`의 부호 있는
sweep으로 하는 일)이고, 원시체 연산의 어휘가 아니다. 원통에 대응하는 복합(면 정박 드릴)은 따로 없다 —
원이 스케치 어휘이므로(칸 ⑨) 면 위 원 스케치의 `PocketOnFace`가 그것이다.

- **방향은 `flip` 이 들고, 프레임을 만든 쪽이 잰다.** 평면의 정준 이름에는 방향이 없고 평면은
  interning 되므로(같은 평면을 `+n`/`−n` 으로 진술하면 **한 핸들**, 실측), 프레임이 방향을
  담을 수 있는 자리는 `flip` 뿐이다. datum 은 **호출자가 진술한 법선**에 대해, 면은 **바깥
  법선**에 대해 잰다. S9 의 «flip 은 진술이 아니라 측정» 규칙 그대로이고, 측정자가 늘었을 뿐이다.
- ★ **`world_zx` 는 유도로 만들 수 없다**: arbitrary-axis 규약이 ZX 에 `+u = −x̂` 를 주는데
  규약은 `+u = +ẑ` 다(`ŵ` 는 둘 다 `+ŷ`). `SketchFrame::world` 가 그 예외를 **한 곳에** 가둔다 —
  `canonical(씨앗 ZX)` 로 바꾸면 그 스케치들이 90° 돈다(음성 대조로 잠금).
- **`dist` 는 두께다**(`NonPositiveDistance` 유지). 방향은 프레임이 말한다 —
  `DatumDef::Offset` 의 `dist` 가 **부호 있는 변위**인 것과 대비되며, 그쪽은 부호만이 «어느 쪽»
  을 말하기 때문이다. 반대편을 향하는 스케치는 **그 방향으로 평면을 진술**한다.
- **프레임의 기저는 유리수로 묻는다, 실현해서 되묻지 않는다.** `RatFrame::of_plane_frame` 이
  `û = u_raw·inv_sqrt_exact(uu)` 로 답한다 — 실현한 축을 `Rat::from_decimal` 로 다시 들어올리면
  정규화가 필요한 축(`(0.6,0.8,0)` → 원시 `(3,4,0)`, `uu=25`)이 `0.6000000000000001` 로 돌아와
  직교정규가 깨지고, 그 평면이 **조용히 프레임-노드 도로로 옮겨간다**(실측). `plane_frame_named`
  가 `v̂` 에 대해 이미 적어 둔 규칙과 같은 것이다.

**datum 평면의 규칙 (S5(i)-a).** 평면을 만드는 연산은 `DatumPlane` **하나**이고, 불리언은
여전히 평면을 만들지 않는다(C6 깊이 불변식 무손상, `a_boolean_mints_no_surface`).

★★★ **정점을 이름 부르는 datum (S5(ii)-1, 2026-08-08).** `ThroughVertices([Handle<Vertex>;3])`
는 값 어휘로는 말할 수 없는 하나다 — 발견 정점의 좌표는 반올림이라, 그 좌표로 평면을 지으면
**다른 평면**이 나온다(기울어진 인구 **220/220**, 축정렬 음성 대조 552/0;
`tests/point_width.rs`). 진실 쪽은 `PlanePoints::Through` 로 **핸들을 든다**.

- **어휘가 여기서만 자란다** — S5(i)-a 의 «새 공개 생성자 없음» 은 *값으로 말할 수 있는 것*에
  걸리는 원칙이고, 핸들은 순수 값 타입인 `SketchPlane` 에 들어갈 수 없다.
- **정렬은 키에만, 방향은 호출자의 정점 순서.** `dist` 가 양수 전용이라 순서가 방향의 유일한
  입구이고, 두 개를 바꾸면 «같은 핸들 + 반대 프레임» 이다(`measured_frame` 이 재는 것은
  `Stated` 와 같은 기계).
- **이동은 핸들을 그대로 두고 노드를 기록한다.** `transform_solid` 는 정점을 **복제**하므로
  가리키던 datum 은 새 복사본을 안 따라간다; 정점이 base 를 정하고 모션이 옮기며 **더해질 뿐
  곱해지지 않는다**. 대가로 그런 datum 위의 솔리드는 **정확한 강체 이동에도 노드를 얻는다**.
- **거절은 원인별**(`VerticesInMixedFrames`·`CollinearVertices`·`DuplicateVertex`·
  `VertexNotThreePlane`)이고, «담체 셋이 안 만난다» 는 거절이 아니라
  **단언**이다 — 그 정점이 존재한다는 것이 곧 만났다는 뜻이므로 불변식 위반이지 사용자 오류가
  아니다. (`VertexPointTooWide` 는 열린 항목 17 로 은퇴 — 폭 무관 이름 유도
  `plane_name_from_meets` 가 그 인구를 named 도로에 태운다; 2026-08-09.)
- ★★★★★ **정정(2026-08-08)**: 출하 직후 검토가 결함 하나를 찾았다 — 이 팔이 세 정점이
  공유하는 **프레임을 계산해 놓고 버리고** 모션을 `None` 으로 저장했다. 기울어진 프레임 위
  프리즘의 먼캡은 그 프레임에서 `w = dist`, 즉 정준 이름이 **세계 평면 `z = dist` 와 같다** —
  둘을 가르는 것이 `SurfaceKey` 의 모션인데 그것을 버렸으므로, 그런 datum 이 **모델에 있던
  상자 윗면 핸들로 돌아왔다**(고치기 전 테스트가 그렇게 죽는다). 점과 모션은 `Offset` 처럼
  **한 결정에서 함께** 나와야 한다 — «점이 어느 프레임에 적혔나» 가 곧 모션이기 때문이다.
- ★★★★ **가장 위험했던 자리는 컴파일러가 못 본 곳이다.** 변종 추가가 낸 비망라 에러는 2개뿐이고
  `transform::points_move` 의 `let`-`else`(원통용 폴백)는 거기 없었다 — 원통은 진실이 기하를
  안 들어 «나를 것 없음» 이 참이지만 `Through` 는 기하를 **참조로** 든다. 그대로 뒀다면 노드
  없는 경로로 가 **캐시만 움직이고 진실은 제자리**에 남았을 것이다.

- **캐시 법선은 `−(진술된 법선)`** — 씨앗이 `−축`, extrude 밑캡이 `−plane.normal()` 인 그 규약.
  근거 둘: S9 가 `+축` 씨앗으로 **781 캡의 저장 법선이 뒤집힘**을 실측했고, `WorkingPlane::
  frame_sign` 이 «저장 법선 vs 루트 면의 바깥 법선» 인데 밑캡의 바깥이 `−N` 이다.
- **연산은 호출자의 프레임을 함께 돌려준다**(`OpOutput::DatumPlane{plane, frame}`) — 편의가
  아니라 **정확성**이다. 호출자가 `SketchFrame::named` 로 다시 말하면 f64 를 거쳐
  `Rat::from_decimal` 을 두 번 타고 계산값은 다른 유리수로 떨어질 수 있다. `Stated` 의 배치는
  **무조건 `Named`**: ZX 평면의 정준 `+u` 는 `−x̂` 인데 규약은 `+ẑ` 라, 유도로 흉내내면 그
  스케치들이 조용히 돈다.
- **`Offset` 은 push 전에 정규화한다 — 셋 다 «한 평면에 두 핸들»을 막는 조항이다.** 신원은
  `(평면, 부호 있는 거리)` 의 순수 함수여야 한다(배치의 원점·`+u` 는 평행 평면을 옮기지 않는다):
  ① `flip` 을 부호로 접고 **평면의 정준 프레임**에서 짓는다 ⇒ 한 평면의 다른 프레임을 든 두
  호출자가 같은 핸들을 받는다. ② 정준 기저가 정확히 리프트되면 **세계로 말한다**(S9 의 노드
  생략 규칙을 실현된 기저에 적용) — 프레임 노드 아래 두면 키가 `(이름, Some(노드))` 라 상자
  캡의 같은 평면과 **interning 되지 않는다**. ③ **`dist == 0` 은 이름 붙은 거절**(`ZeroOffset`):
  기울어진 프레임의 밑 평면은 세계에서 유리수라 ②가 안 걸리므로, 0 만이 남는 충돌이다. 0 이
  아닌 오프셋은 계수가 `c ∓ d·√(n·n)` 이라 무리수여서 경쟁 진술 자체가 없다. 세계 되당김이
  넘치면 **조용히 노드 도로로 가지 않고** 거절한다 — 그 도로가 곧 중복이 생기는 자리다.
- **datum 은 솔리드 이동을 따라가지 않는다.** append-only 라 datum 이 가리키는 평면은 남고,
  나중 `Transform` 은 새 평면을 만든다(의도된 동작 — datum 은 독립적인 기준이다).
- **쓰이지 않는 datum 은 어디에서도 새지 않는다**: `validate` 는 live 도달 집합만 세고 surface
  는 경계 검사만 받으며(씨앗 셋이 이미 영구 orphan), `nacre-step`·`nacre-tess` 는 surface
  store 를 돌지 않는다(둘 다 면 경유).

pub struct TessConfig { pub default_tol: f64 /* , 면별 override 등 */ }

// nacre-ops의 자유 함수다 — Model의 메서드가 아니다. `impl Model`은 inherent impl이라
// Model이 정의된 nacre-topo에만 놓을 수 있는데(coherence), replay는 Operation·Tessellation을
// 다루므로 topo보다 위 레이어(ops는 topo·tess 위)에 살아야 한다. Model 메서드로 두면
// §2에서 tess/ops를 Model 필드에서 뺀 것과 같은 층 위반이 된다.
///
/// 로그를 처음부터 재생. 보장: 동일 로그·동일 cfg → 동일 (모델, tess) (인덱스까지 재현).
/// **7변종 전부**에서 성립한다(S5(i)-b 이후 값-전용 변종은 없다) — replay 가 인덱스를 재고정하기
/// 때문이며(§2 「인덱스 어휘」), 측정은 `tests/replay.rs`(생성 세션 + 6변종 전수)다.
/// 로그 중간 파라미터를 수정한 재생은 v1에서 미지원 — Operation이 원시 Handle을
/// 참조하므로 상류 수정이 하류 Handle 번호를 밀어낸다 (topological naming 문제).
/// 반환 쌍: 진실 Model + 함께 생성되는 tess 캐시(§5 커플링; Model은 tess를 필드로 담지 않으므로 §2).
pub fn replay(ops: &[Operation], cfg: TessConfig) -> Result<(Model, Tessellation), OpError>;
// M2 현재형: tess가 아직 없으므로 `replay(ops: &[Operation]) -> Result<Model, OpError>`(cfg·Tessellation 없음).
// M3에서 tess 도입과 함께 위 (Model, Tessellation)·TessConfig 시그니처로 확장한다.
```

**다중 루프 프로파일 — 구멍 있는 스케치를 정확히 세운다 (2026-07-26 구현).** `Profile2d`가 외곽 링 하나와 구멍 링 N개를 담고(`polygon`/`with_holes` 생성자, 필드는 비공개), `build_prism`이 구멍마다 벽을 세우고 두 캡에 내부 루프를 단다. **불리언으로 흉내낼 때와의 차이가 요점이다** — 잘라 만든 도넛은 모든 코너가 `Discovered`(측정된 tol)이지만, 프로파일을 쓸어 만든 도넛은 전부 `Constructed`(tol 없음)다. 원칙 4가 요구하는 바로 그 차이다.

**감김은 한 곳에서만 정한다.** `build_prism`이 sweep을 아는 유일한 자리이므로 외곽을 sweep 기준 CCW로 정규화하고 구멍을 그 반대로 맞춘다. 타입은 링을 그대로 담고, 프레임 2D 면적으로 한 번 더 정규화하던 `placed_profile_unchecked`의 판단은 **삭제**했다(두 곳에서 정하면 pocket처럼 sweep이 외곽을 뒤집는 경로에서 둘이 갈린다). **측정 결과 이 설계는 그 버그 부류를 구조적으로 없앤다** — "정규화된 외곽의 반대"와 "sweep 기준 CW"가 같은 규칙이 되기 때문이다.

검산: 도넛 각기둥은 V16 − E24 + F10 − L_i 2 = 0 = 2(S−G), genus 1 — **오일러의 `L_i` 항이 처음으로 실제로 필요한 형상**이다.

**프로파일 계약은 커널이 검사한다 — `Profile2d::check` (2026-07-26; S3 개정 2026-08-05).** 생성자(`polygon`/`with_holes`)는 S3부터 **진실을 잡는다**: 각 좌표를 `Rat::from_decimal`로 리프트해 `Ring2d{points: Vec<[Rat;2]>}`에 보관하고(십진 창 밖 = `ProfileOutsideDecimalWindow` 구성 시점 에러 — 과거의 조용한 f64 폴백 폐지), **공선 중간점을 소멸**시킨다(엄격 내부만 — 무손실 정규화; 중복점·스파이크는 생존해 제 이름의 에러로 보고된다. 이로써 모든 프리즘 코너가 3-평면 정의를 보유한다 — truth-and-cache Q2 ②). *계약*(단순성·서로소·포함)은 여전히 생성자가 아니라 **프로파일을 소비하는 모든 연산이 먼저 `check()`를 통과시킨다**(진입점은 `extrude`와 pad·pocket이 공유하는 `extrude_and_boolean` 둘뿐). 계약은 넷: 모든 링이 **단순 다각형**, 링끼리 서로소, 구멍은 외곽 안, 구멍 속 구멍 없음(그건 섬이고 `from_rings`가 갈라낸다). 판정은 전부 **저장된 유리수 진실 위**에서다 — f64 이진값의 정확 부호와 십진 진실의 부호는 퇴화 근처에서 실제로 갈리고(0.1·0.2·0.3 공선이 이진에선 굽음), 작성자가 쓴 수가 이긴다.

**왜 "호출자의 약속"으로 둘 수 없었나 — 실측.** 계약을 어긴 프로파일은 전부 `extrude` 성공 · `validate` 위반 0건 · STEP 내보내기 성공이었다. 나비넥타이는 부피가 `NaN`이었지만, **구멍이 외곽 밖이면 12.0(정답 16), 구멍 속 구멍이면 20.0(짝-홀 정답 52)** — 눈치챌 단서가 없는 그럴듯한 숫자다. 조용히 틀린 답은 약속으로 둘 수 있는 종류가 아니다.

부수적으로, 이 검사가 `oriented_ring`의 **적혀 있지 않던 전제를 복원한다**: 감김을 부호 있는 면적 벡터로 정하는데 면적이 0이면 그 판정이 우연에 맡겨진다(대칭 나비넥타이는 두 엽이 정확히 상쇄해 0이다). 단순 다각형은 면적이 0일 수 없다.

**단순성은 *작성된 입력*의 계약이지 커널 데이터의 불변식이 아니다.** 그리는 링의 안쪽은 단순성 위에서 짝-홀로 *정의되지만*, 불리언이 만드는 윤곽의 안쪽은 arrangement가 이미 정했고 그쪽은 **비단순(figure-8 조임)이어도 정당하다**(`loop_winding`). 그래서 `validate`에 면 단순성 검사를 넣으면 안 된다 — 미루는 게 아니라 정당한 결과를 거절하게 된다.

**링 중첩 판정 — `from_rings` (2026-07-26 구현).** 사용자는 닫힌 경로만 그리고 무엇이 구멍인지 말하지 않는다. `ops::sketch::from_rings`가 포함 **깊이**로 정한다: **짝수 = 재료, 홀수 = 구멍**(그래서 구멍 속 링은 다시 재료 = 섬이고 자기 몫의 profile이 된다), 구멍은 **가장 깊은 포함자**(직계)에 붙는다. 채우기 규칙 파라미터는 두지 않는다 — 짝수-홀수가 유일한 규칙이다.

**분류 술어는 `geom::intersect`에 있다**(§1의 격리 규칙): `point_in_ring_2d`(exact `orient2d` 교차 패리티, `RingSide::{Inside, Outside, OnBoundary}`)와 `rings_cross`(적절 교차 + 접촉). S3부터 각 워커에 **Rat 쌍둥이**(`_rat` 접미 + `drop_collinear_midpoints`)가 같은 모듈에 나란히 산다 — 부호 원시는 `nacre_scalar::orient2d_rat`(narrow 우선 → BigInt 전역, 분모 청소로 부호 보존), geom이 scalar 의존을 얻었다(최하단 토대라 순환 없음). f64 판은 tess·f64 폴백이 계속 쓴다. ops에는 **정책(깊이 패리티)과 조립**만 남는다. 링이 서로 닿거나 교차하면 `SketchError::RingsMeet`으로 거절한다 — 접촉도 실패다(엄밀한 안쪽이 없다).

**자기교차는 중첩 분류보다 먼저 본다** — 메시지 품질이 아니라 전제조건이다: `point_in_ring_2d`는 짝-홀 패리티로 답하고, 그건 단순한 링에서만 "안쪽"을 뜻한다. 술어는 `geom::intersect::ring_self_intersection`(인접하지 않은 변은 접촉만으로 실격, 인접한 변은 공선-겹침일 때만 = 되짚는 스파이크, 인접은 **순환**으로 판정). 스케치 층은 `RingSelfIntersects { ring, at }`로 **점**을 돌려준다 — `from_edges`는 순회 순서로 링을 만들어 변 인덱스가 작성자의 입력과 무관하기 때문이다(`OpenChain`·`BranchingVertex`와 같은 규약).

**선 뭉치 — `from_edges` (2026-07-26 구현).** 생성 코드나 외부 데이터가 내놓는 모양 그대로, **순서도 방향도 무관한** 선들을 받아 끝점으로 이어 링을 만들고 `from_rings`에 넘긴다. ★ **호와 원 (2026-09-05, 칸 ⑨).** `Edge2d::{ Line { from, to }, Arc { center, radius, start, end, ccw } }` — 전부 `Rat`, `start == end`는 온전한 원(시작점 = 솔기). 회전각은 커널이 모르는 말이다: 끝점 표현이 각도 표현보다 **넓고**(무리수 도의 호도 끝점이 유리수면 정확), 90° 배수 회전은 생성자가 `Rat`으로 끝점을 계산해 준다(`arc_turns`·`arc_turns_rat`; 반지름은 `|start−center|²`의 **완전제곱** 검사 — `rat_sqrt_exact`). 링은 **정점 + 세그먼트 종류**(`Ring2d { vertices, segs: Vec<Seg2d> }`)라 끝점 중복이 타입에 없고, 정규형이 공선 중간점 제거 옆에서 **같은 원·같은 방향의 이웃 호를 병합**한다(필렛 둘이 만나면 반원, 호 둘로 그린 원은 원 하나). 분류 술어는 `geom::mixed`(선분–원·호–호·점-in-링을 ℚ(√c)의 `QuadVal` 부호로, 광선 규칙은 배열의 반열림 철자를 옮겨 적음; 다각형 링에서 옛 술어와 답 동일을 proptest 로 잠금). 프리즘 빌더는 호마다 원통 조각 벽을 세우고(축 = 프레임 법선 `w`, `ref_dir` = 프레임 `x̂`, 온전한 원만 시작점 방향), 직선–호 꼭짓점은 **`VertexDef::Branch`**(두 평면의 저장된 정준 이름을 핸들 오름차순으로 `plane_plane_cylinder`에 — 접하는 벽은 `Double`), 감김은 `2·넓이 = a + b·π`의 정확 부호(`winding_sign_quarter_arcs`, scalar 의 BigInt). 옆면 방향 `Forward ⟺ 호의 ccw == 스윕이 +w`. 이름 붙여 거절: 90° 배수 아닌 호(`ArcSweepNotQuarterTurn`), 다른 원의 호–호 접합(`ArcsMeetAtVertex` — 두 원통 위의 점은 정점 정의가 없다), 반지름 무리수(`ArcRadiusNotRational`). **원통 원시체는 이 어휘의 설탕이 되어 은퇴했다**(위 `Operation` 목록).

**끝점은 정확히 일치해야 한다.** 가까우면 붙여주는 스냅은 없다 — tolerance는 커널이 *발견한* 교차의 것이지 호출자가 *구성한* 것의 몫이 아니다(원칙 4). 대신 벌어진 끝점과 **가장 가까운 다른 자유 끝점까지의 거리**를 `OpenChain { at, gap }`으로 돌려준다(고치는 데 필요한 숫자가 그것이다). 그 밖의 거절: 영길이·중복 선, 그리고 세 개 이상이 만나는 `BranchingVertex`. **더 구체적인 결함을 먼저 보고한다** — 분기는 항상 어딘가에 홀수 끝점을 남기므로, 순서를 반대로 하면 늘 모호한 쪽(열림)이 보고된다(`check_result_topology`와 같은 원칙). 구멍 있는 스케치(도넛)와 섬이 여러 개인 스케치를 코드-CAD가 요구한다. **설탕으로 흉내내면 안 된다** — "외곽 extrude → 구멍 프리즘 Cut"은 전부 `Constructed`였을 모델을 불리언·`Discovered` 경로로 내리므로 원칙 4(tolerance는 발견된 교차에만)를 스스로 어긴다. 커널은 이미 대부분 준비돼 있다: `Face { inner: Vec<Loop> }` 존재, `nacre-tess::polygon`이 구멍 여럿을 네이티브로 다루는 삼각분할, `validate` 오일러의 `L_i` 항. 막는 것은 입력 타입 하나(`Profile2d { points: Vec<Point2> }` = 폴리곤 하나)다.
- `{ outer, inners }`(구멍 N개, 제한 없음) + `extrude`가 구멍 벽면과 뚜껑 내부 루프를 함께 생성.
- `Profile2d::from_rings(rings, fill_rule) -> Vec<Profile2d>` — 링 목록의 중첩을 exact `point_in_ring`으로 판정해(술어이므로 커널) **덩어리(섬)별 프로파일 목록**을 돌려준다(짝수-홀수 깊이: 0=재료, 1=구멍, 2=구멍 속 섬…). 채우기 규칙 선택은 호출자.
- **경계:** 커널 `Extrude` 1회 = **연결된 덩어리 1개**(외곽 + 그 구멍들). 섬마다 호출해 결과를 묶는 것은 편의 레이어(overview.md 판별 기준). 그래서 `Extrude`의 다중 바디 출력은 필요 없다.

**파라메트릭 편집의 진화 경로 (v2 이후, 지금은 기록만):** 연산이 원시 Handle 대신 계보 참조 `OpRef { op: usize, output_slot: usize }`("연산 N이 만든 k번째 면")를 담으면, 상류 수정 후에도 참조가 의미로 해석(resolve)되어 편집-재생이 가능해진다. v1에서 이를 구현하지 않되, 로그 직렬화 포맷을 설계할 때 이 확장이 포맷 파괴 없이 들어갈 자리를 남긴다. §2의 인덱스 어휘는 이를 **밀어내지 않는다** — 인덱스는 v1의 참조 기계이고 `OpRef`는 v2의 것이며, replay 안의 재고정 지점이 곧 `OpRef` 해석이 들어올 자리다(교체이지 경쟁이 아니다). TessConfig의 tolerance도 같은 맥락에서 연산별 override(`Operation` 항목의 선택 필드)로 확장될 수 있다.

불리언은 처음부터 trait 뒤에 둔다 — 커버리지 사다리의 코드화:

```rust
pub trait BooleanEngine {
    fn boolean(&self, m: &mut Model, kind: BoolKind,
               a: Handle<Solid>, b: Handle<Solid>) -> Result<Handle<Solid>, BoolError>;
}

pub struct PolyhedralBoolean; // M5: 평면 솔리드 전용 — exact 술어로 진짜 강건.
                              //   커버리지 밖 입력은 BoolError::Rejected로 정직하게 거절
pub struct QuadricBoolean;    // M6: 평면+이차곡면 (닫힌 형식 교차)
pub struct HybridBoolean;     // M7: 일반 곡면 — 출처태그 메시 → 조합 결정 → 스냅백.
                              //   주의: 조합 결정이 메시 해상도에 의존할 수 있으므로,
                              //   "위상 결과의 tolerance 불변"을 목표 불변식으로 삼고
                              //   위반 사례를 실패 코퍼스에 축적한다
```

OCCT는 제품 경로에 등장하지 않는다 — 역할은 nacre-oracle의 채점자(§7)뿐이다. 사다리의 각 단은 자기 커버리지 안에서 완전해야 하며, 밖은 조용히 틀리는 대신 에러로 거절한다.

**거절의 *이유*는 값에 실려 나간다 (2026-07-26 구현; 변종명 `Unsupported`→`Rejected` 개명 2026-08-16 — 세 등급 중 `Impossible`·`SuspectedDefect` 둘에 대해 «미지원»이라는 이름이 거짓이었다).** `BoolError::Rejected { reason: RejectReason }` — 이름 붙인 거절이 크레이트 밖에서 하나의 불투명한 에러로 붕괴하던 것을 값으로 옮겼다(그전엔 태그가 `#[cfg(test)]` thread-local에만 기록됐고, 여러 지점이 거절을 울린 뒤 삼키므로 애초에 건전하지도 않았다). **소비자는 `RejectReason::class()`로 분기한다** — `NotSupported`(커널의 현재 커버리지 밖 — 입력 자체는 판정 안 됨; 변종명 `NotSupportedYet`→`NotSupported` 개명 2026-08-17, 이름의 «Yet»이 지원 예정이라는 예측으로 읽혀서) / `Impossible`(어떤 마일스톤에서도 유효 솔리드가 없음) / `SuspectedDefect`(엔진 불변식이 깨짐). 변형 이름은 엔진 어휘라 로그·리포트용 안정 식별자로만 쓴다. 사람이 읽는 문장·현지화는 앱 몫이다.

성장은 `RejectReason`에서만 일어나므로 그것만 `#[non_exhaustive]`이고 `BoolError`·`RejectClass`는 exhaustive다(소비자가 완전히 처리할 수 있게). 트레이스가 불완전한 경우는 `TraceDeclined { kind, face }`가 **무엇을 못 했는지와 어느 피연산자 면에서인지**를 함께 싣는다 — 예전엔 서로 다른 10가지 사유가 전부 `COPLANAR_PAIR` 하나로 나가 커널이 틀린 말을 했다.

**«어디»는 사유 옆에 실려 나간다 (F, 2026-08-17).** `BoolError::Rejected { reason, at }` — `at:
Option<RejectWhere>`(`Point` | `Segment`, 월드 f64)는 가드가 발화한 순간 보고 있던 **증인**이다
(여럿이면 엔티티 자체 순서의 최솟값 — HashMap 순회로 뽑으면 실행마다 다른 좌표가 나온다). 위치를
`RejectReason` 변종 안에 넣지 않는 이유: 사유는 범주 어휘(census 키·테스트가 이름 대는 것, `Copy+Eq`
·const-구성)이고 좌표는 측정값이다 — census 의 "측정값은 키 밖" 규칙의 오류-값 판. 좌표는 진단용
실현(캐시급)이라 `RejectWhere` 는 `Eq` 가 없고, 비교는 근사로만 한다. 표면화 사유 중
`self_touching_result`(위반 모서리)·`non_manifold_result_edge`(위반 모서리)·`non_manifold_vertex`
(핀치 정점)가 싣고(`coplanar_pinch` 는 변종과 함께 은퇴 — 2026-08-17, 아래 폴드 절), `no_clear_ray`(분산적)·
`precision_budget`(전-모델)·`cylinder_face`(발화 인구 0 — 필요해지면 `RejectWhere::Face` 변종으로
연다)는 싣지 않는다.

**어떤 사유가 실제로 발화하는지는 상설 census가 답한다 — `nacre-ops::reject_census` (2026-08-15).**
`reject()`/`reject_at()` 쌍이 크레이트의 모든 `Rejected`를 짓는 유일한 깔때기라, 거기
`#[track_caller]` 하나씩이 모든 호출 지점을 한꺼번에 계측한다. 무조건부(선례: `nacre-topo::WIDE_PLANES`, `nacre-cip::climb_census`
— `#[cfg(test)]`는 통합 테스트가 비-test 빌드를 링크해서 못 쓴다), 비용은 raise 1회당 **~11ns**(실측).

★ **두 열은 서로 다른 인구다.** 「울린 것」과 「밖으로 나간 것」을 따로 적는다: E 이전(2026-08-15)
워크스페이스 전 스위트 실측으로 **raise 151 / surfaced 26**이었고, 그중 122건이 `combinatorics.rs`
한 자리의 `no_clear_ray`인데 **거기서는 한 번도 표면화되지 않았다**(호출자가 링의 다음 노드로
재시도한다). 한 열로 세면 *"거절의 82%가 no_clear_ray"*라고 보고하게 되고, 그건 사용자가 보는 것에
대해 거짓이다. ★ **E(2026-08-17)가 그 인구를 의도적으로 없앴다**: `point_in_component`의 «이 노드로는
못 정함»은 이제 `Ok(None)` 기권이지 삼켜질 오류가 아니다 — 두 열의 갭은 「남은 삼킴의 지표」가 됐고
0 에 가까운 것이 목표 상태다(갭이 다시 벌어지면 새 삼킴의 발견이다). ★ 다만 census는 **어느 사이트가 표면화했는지는 답하지 않는다** — 거절은 rayon
워커에서 울리고 메인 스레드에서 반환되므로 둘을 잇는 값싼 길이 없고(옛 thread-local 태그가 폐기된
이유), 「이 가드는 삼켜지는가」는 호출자를 읽으면 되는 **정적** 질문이라 필요도 없다.

읽는 길 둘: **게이트**는 `tests/reject_census.rs`의 얼린 코퍼스 7종이 `(사유, detail, 파일)` 집합을
박는다(줄 번호는 찍기만 — 거기 걸면 무관한 편집마다 빨개져 아무도 안 읽는다. 횟수도 안 건다: debug는
트레이스를 두 번 돌리고 parallel은 실패 입력의 나머지 클래스까지 평가한다). **전수**는
`--features reject-trace`로 스위트 전체를 훑는다 — 테스트 바이너리 30여 개가 각각 별개 프로세스라
메모리 표가 프로세스와 함께 죽기 때문이고, 그래서 이 기능이 크레이트의 `deny(print_stderr)`에 대한
유일한 예외다(기본 꺼짐, 제품 빌드엔 출력 코드가 아예 안 들어간다). 게이트의 인구는 그 7개 모양뿐이니
**M6 가드는 여기 안 걸린다** — 그건 전수 훑기가 손으로 메운다.

census 키가 `detail`을 드는 이유: `TraceDeclined`의 11개 kind가 전부 **한 줄**에서 올라와, 이름과
위치만으로 키를 잡으면 한 행으로 뭉개진다 — census가 잡으려는 결함을 census가 재현하는 꼴이다. 규칙은
**키가 그 사유 자신의 어휘만큼 잘아야 한다**는 것.

### 6.0 구성 산술은 **사용자가 쓴 십진수**로 한다 (2026-07-29)

프리즘의 좌표는 곱 셋과 합 둘로 나온다 — `origin + x·u + y·v`, 그다음 `+ normal·dist`. f64에서 이건 결합법칙을 안 지키므로 **정확히 일치해야 할 두 치수가 어긋난다**: `7.7`로 한 번에 올린 블록과 `1.1` 다음 `6.6`으로 올린 블록의 윗면이 1 ULP 떨어진다.

**쪼개지지도, 틀리지도 않는다 — 더러워진다.** 부피는 정확히 맞고 fuse는 몸통 하나를 준다. 그 몸통이 넓이 `8.9e-16`짜리 면과 형상에 없는 면 둘을 달고 있고, 이후 모든 불리언·STEP·메시가 그걸 끌고 다닌다. **거절도 오답도 아니라서 아무도 신고하지 않는다** — 이 셀이 다룬 결함 중 가장 조용한 부류다.

**★ "API를 `Rat`으로 바꾸면 된다"가 답이 아니다.** f64 `1.1`을 `try_from_f64`로 *올리면* `11/10`이 아니라 그 f64가 담은 이진값이다. 정확한 산술을 해도 이진 드리프트가 그대로 재현될 뿐 일치가 **복구되지 않는다**. 복구하는 것은 사용자가 *쓴 십진수*다 — `11/10 + 66/10 = 77/10`, 그리고 그 실현이 정확히 `7.7`.

**결정:**

- **공개 API는 f64 그대로.** 호출부·TS/JSON 경계가 안 바뀐다. 커널이 받자마자 각 치수를 `Rat::from_decimal`로 **왕복하는 최단 십진수**로 잡는다. S3이 프로파일 좌표에, S6a가 스케치 평면의 축(`from_axes` — 정의는 점 셋 `[o, o+x, o+y]`, `PlaneDef`는 점 셋 단일 필드로 origin·ref_dir·극성이 구조에서 유도)에 이 규칙을 구현했다. ★ 단 **커널이 실현한(반올림된) 기저는 리프트하지 않는다**(내부 `realized_plane` — 이미 정확한 정의를 가진 평면에 두 번째 진실을 만드는 것이 두-정확-기술 결함의 재발이다). 의도에 대한 *추측*이 아니라 **결정적 정규형**이다(왕복이 보장되므로 서로 다른 f64는 서로 다른 유리수, 같은 f64는 항상 같은 유리수). 십진수를 고르는 이유는 하나 — **사람은 십진수를 타이핑한다**.
- **전제조건: `Rat::to_f64`가 최근접 반올림이어야 한다.** 분자·분모를 각각 반올림한 뒤 나누면 유효숫자 17자리에서 1 ULP가 어긋나고, 그러면 십진수 복원이 **무손실이 아니다** — 고치려던 것은 어긋난 누적뿐인데 멀쩡한 입력까지 움직인다. 정수 장제법으로 54비트 몫을 뽑아 **끝에서 한 번만** 반올림한다(`Ratio<i128>`이라 두 항이 2¹²⁷ 미만인 게 모든 시프트를 닫는다).
- **프레임이 유리수일 때만 이 경로를 탄다.** 축이 `{0, ±1}`이면(world·축정렬 면) 정확하다. 회전된 면은 축이 무리수라 못 든다 — 성분을 십진수로 잡으면 **단위벡터가 아니게 되어** `normal·dist`가 틀린 *길이*를 준다(고치려던 것보다 나쁜 오차). 그래서 진입점이 `x·x = y·y = 1`, `x·y = 0`을 **정확히** 검사하고 아니면 오늘의 f64 경로로 강등한다. i128 넘침도 같은 답이다.
- **새 저장 구조가 필요 없다.** 같은 유리수는 항상 같은 f64로 실현되므로 두 경로의 좌표가 **비트 동일**해지고 평면 계수도 같아진다 — "계수가 곧 진실"이라는 그때의 서술이 그대로 유효하다. 모션 셀처럼 provenance를 새로 다는 일이 없다.

**못 고치는 것 (비목표, 정직하게):**

- **회전된 스케치 프레임.** 축이 무리수라 유리수로 들 수 없다. **f64 삼각함수를 거친 사분 회전도 마찬가지다** — `(90°).to_radians().sin_cos()`가 `cos = 6.1e-17`을 주므로 정확성은 이 층이 보기 전에 이미 사라졌다(커널의 정확한 사분 회전은 `Angle`/Niven에서 온다). 그 경로의 건전성은 CIP 담당이다.
- **스크립트가 계산해 넘긴 값.** `1.1 * 7`을 f64로 계산해 넘기면 커널은 `7.700000000000001`을 받고, 그 최단 십진수는 `7.7`이 아니다. **올바른 동작이다** — 스크립트가 실제로 다른 수를 계산했다. 고치는 것은 **커널 자신의 누적**이고, 그게 실측된 결함이었다.
- **4-평면 동시성.** ±1 ULP로 사라지는 정확한 0이라 유리수 입력과 무관하다(§9 해당 항목).

**실측(전후 동일 코퍼스):** 쪼갠 십진 치수 600표본에서 슬리버 **11.0% → 0%**(정수 0%·소수1자리 13.0%→0%·소수2자리 20.0%→0%). 정수가 원래 0%인 것이 진단이다 — 정확히 표현되는 값은 합도 정확하므로, 결함은 **전적으로 십진수의 문제**였다. 대장 130행 중 3행만 바뀌었고 셋 다 `ex` 계열의 `(7.7, 1.1, 6.6)` 한 케이스다(`0.1+2.9`·`0.3+0.7`은 f64에서 이미 정확히 떨어져 고칠 게 없었다). `add_cuboid` 121행은 산술이 없으므로 비트 동일.

### 6.1 M5 불리언 — 두 메커니즘을 regime로 라우팅 (합성 아님)

> **★ 방향 전환(2026-07-20, 사용자 결정 — 이 절은 이행기 서술).** 아래 "공면 접촉 유무로 배타 라우팅"(detector + 공면 생존표)은 케이스가 늘수록 **경우의 수가 폭발**한다(same_ground·single-shared·containment로 실증). 그래서 design.md가 **M7**(하이브리드, §426 exact ray casting)에 두던 **arrangement/winding 통합 분류 방식을 M5 평면으로 앞당긴다** — 평면은 메시·SSI 불요라 M7 분류법을 exact plane-triple에 그대로 적용. **M5의 목표 = winding 기반 단일 arrangement 엔진**(면당 세분 → sub-face를 in/out/on 분류 → op별 Requicha keep 균일 적용; 글로벌 detector·이중 경로 소거, Requicha 규칙은 유지). 이 방식은 커토버로 프로덕션 엔진(`nacre-ops::arrangement`)이 됐다 — 면당 평면 cell-복합체(트레이스→분할→셀→중첩→라벨→방출→조립). **잔여 커버리지 갭:** ~~containment(seam 없는 포함)·cavity 접촉(현재 outer shell 순회)~~ — 둘 다
닫혔다(2026-08-15 실측 확인). containment 는 fuse/cut/common 셋 다 정확한 부피를 내고
`coverage/coplanar.rs` 와 `coverage/invariants.rs` 의 proptest 가 상설로 잰다; 공동은 안으로
파고들기·벽에 공면으로 앉히기·정확히 메워 공동을 없애기가 모두 정확하다. 회전된 «두 몸통» 사이의
공면 접촉도 정상 동작한다(0·17·30·45° 실측). ⇒ **M5 커버리지 갭은 남아 있지 않다.**

★★ **45°/315° 폴드의 거절은 갭이 아니라 «정답»이다.** 그 각도에서 부품의 한 팔이 **자기 자신에**
공면으로 닿아 두께 0 인 솔리드가 나온다 — 존재할 수 없는 물건이고, 거절이 옳다.
**이름의 여정(2026-08-16→17)**: 하루는 `CoplanarPinch`(병합 가드가 실제로 보는 사실 — 능력
한계, `NotSupported`)였고, 다음 날 스파이크가 지역-한계 사슬의 깊이를 끝까지 재자(정확히 2겹,
그리고 핀치 코너의 이름-실패는 결함이 아니라 **비다양체 기하의 증상**) 수리가 착지했다: 병합은
핀치 그룹에 **기권**(`Ok(None)` — 정리 단계는 못 정리하면 통과한다)하고, `self_touch_reject` 가
발행 전으로 옮겨져 **`SelfTouchingResult`(`Impossible`) + 접촉선 세그먼트 증인**이 답이다.
`CoplanarPinch` 는 생후 2일에 고아화되어 삭제 — 증상 이름이 죽고 진실 이름이 대신하는, 오류
사다리의 목표 상태. figure-8 «재봉합» 기능은 소비자가 없다(전역 검사 둘을 다 통과하는 핀치
인구가 나타나는 날이 그 첫 소비자다). 잠금: `rotation_sweep.rs::the_other_rejections_are_untouched`
(세그먼트의 캡-스팬·수직성), `reject_census` fold-45. ~~회전 접촉~~ — 실체는 회전이 아니라 **칼날 접촉**(모서리가 상대 면 평면 위, 두 면이 한쪽으로만 떠남)이었고 2026-08-11 에 닫혔다: `edge_mask` 의 스침 결합이 합집합→**홀짝**(같은 쪽 쌍 = 꼬집힘/노치 = 무소식), all-false 모서리는 골격에서 제외(`drop_newsless` — 점 접촉 `touches` 의 1D 판). 모서리-만 접촉 fuse 는 **몸통 둘로 나온다**(2026-08-14) — §824 가 *"엣지 접촉 Fuse가 전부 정상 결과"* 로 이미 정해 뒀던 줄의 이행이다. 면을 잇는 것은 **다양체 접촉(정확히 두 면이 쓰는 링 모서리)뿐**이라, 선·점으로만 닿는 두 몸통은 서로 다른 성분에 놓이고 각자의 핸들을 받는다. 핸들을 가르는 단위는 성분이 아니라 **출력 솔리드**(재료 성분 + 그 공동들)여야 한다 — 공동이 host 껍질에 닿는 경우 성분별로 가르면 핀치가 셀 수 없어져 두께 0 솔리드가 조용히 통과한다(실측). 자기와 닿는 «한» 몸통은 재료가 접촉을 돌아가므로 성분이 하나로 남아 기존 거절이 그대로 발화한다(잠금: `tests/contact_separates.rs`, `tests/knife_edge.rs`).
>
> **★★ 갱신(F2 — detector 붕괴 완료).** 위 "경우의 수 폭발"의 실체가 **라우팅뿐**이었음이 실증됐다: 케이스별 detector 7종이 **전부 같은 `coplanar_result_unified`를 호출** — 결과 생성기는 이미 통합돼 있었다. 그래서 dispatch를 **하나의 exact 질문**(`coplanar_contact_count >= 1`, "진짜 공면 접촉인가")으로 붕괴시키고 detector 12종·지원 타입 ~775줄을 삭제했다. 부수적으로 그 좁은 게이트들이 막던 케이스가 열렸고(비볼록 오버행 footprint, 관통 slot/corner cut), 라우팅이 넓어지며 드러난 "성공하되 열린 셸" 한 건은 **`assemble_fuse_cut`의 닫힘 가드**(모든 모서리 정확히 2회 사용, 위반 시 `NON_MANIFOLD_EDGE` 정직 거절)로 차단했다. ∴ **아래 배타-라우팅 서술은 여전히 유효하되 "detector 다발"이 아니라 "질문 하나"이며**, 두 메커니즘(seam / coplanar)의 공존은 폐기 대상이 아니라 **구조적 필연**이다 — (단계4 SoS Cell 4) 실측이 "seam 하나로 통일"을 반증했다(공면 접촉 모서리는 공유 평면에 통째로 누워 transversal seam 자체가 부재 = 구조적 공면성, 섭동으로 해소 불가). 남은 winding 작업은 엔진 대체가 아니라 **커버리지 확장**(공면 arm의 회전·cavity·containment)이다.

**★ 4평면 동시성의 세 갈래(2026-08-14).** 한 점에 평면 넷이 모이면 이름이 여러 개 생기고, 배열은
그것을 하나로 접는다(`Aliases`). 발견 경로는 셋이고 마지막이 오래 비어 있었다: ① **입력 면의
꼭짓점**에서의 동시성 — 링을 걸으며 `third_on_l` 이 본다. ② 배열이 만든 점에서 **선을 가로지르는**
두 평면이 겹치는 경우 — `split_at_crossings` 가 「두 handle 이 같은 자리로 정렬된다」로 본다.
③ 배열이 만든 점에서 **그 선을 «담은»** 평면 — 그 평면은 선과 같은 방향 family 라 handle 이 되지
않으므로 ①②가 **구조적으로 못 본다**. 셋째는 새로 찾을 필요가 없었다: 벽 둘을 한 선으로 접었다는
기록이 곧 「두 번째 평면이 이 선을 담는다」는 진술이므로, 그 선 위 분할점마다 벽 family 를 되읽어
신고한다. 놓치면 한 점이 두 이름으로 seam 표에 도착해 `SeamAlias`(`SuspectedDefect`)로 거절되는데,
**그건 유효한 입력을 두고 커널이 자기를 탓하는 것**이었다. 실측·잠금은 `docs/dev-log.md` 같은 날 항목과
`tests/concurrent_line.rs`.

**★ 결과 계약(2026-08-13): 표면이 자기와 닿는 솔리드는 만들지 않는다.** 조립이 내놓은 결과
솔리드(바깥 껍질 + 공동)에 대해 *"이 솔리드의 모서리가 같은 솔리드의 어떤 면의 «내부»에 있는가"*
를 묻고, 참이면 `SelfTouchingResult`(`Impossible`)로 거절한다. 이런 몸체는 연결돼 있고 부피가
맞고 **모든 위상 계수를 통과한다** — 접촉선에서 면이 쪼개지지 않아 모든 모서리가 여전히 두 번씩
쓰이므로 `check_result_topology` 도 `validate` 도 보지 못한다. 방향은 Parasolid 기준(사용자 결정,
2026-08-12): OCCT 는 이런 몸체를 받지만 따르지 않는다. **대가**는 자연스러운 모델링 동작 하나가
막히는 것이다(쐐기 포켓의 끝이 마침 반대 벽에 닿는 경우). 다른 커널이 *"zero thickness
geometry"* 라 부르는 조건과 같다.

★ **거절은 «무엇이 잘못됐는지»만 말한다.** 어느 일치가 의도치 않은 것이었는지, 답이 다른 치수인지
다른 연산인지 몸체를 둘로 나누는 것인지는 **작성자의 의도**이고 커널은 그걸 볼 수 없다 — 「치수를
조금 옮기라」는 그 의도에 대한 추측이다. `BoolReport` 가 같은 규칙을 이미 적어 뒀다:
*"It is a diagnosis, not a prompt."* 커널은 사유 이름과 클래스만 내놓고 **사람이 읽을 문장은 앱의
몫**이다(`RejectClass` 로 분기, 변주 이름으로 분기하지 않는다).

판정은 체 셋을 싼 것부터 통과시킨다 — ① 후보 평면 = 두 끝점 이름 삼중의 교집합(정확, 좌표 없음)
② 그 평면 위 면들 사이에서 정점 `tol` 로 **부풀린** 경계상자 선별(한 평면이 면 80개를 일 수 있다;
실측 후보 2,538,703 → 203) ③ 남은 것만 정확 판정. 상자를 부풀리는 것은 선택이 아니다:
상자는 반올림된 실현 좌표에서 나오고 정확 판정 **앞에서 후보를 버리는** 데 쓰이므로, 보수적이지
않으면 진짜 접촉이 조용히 사라진다.

**★ ③은 끝점이 아니라 «열린 구간»을 묻는다(2026-08-14).** 끝점만 보면 **관통 슬롯**을 놓친다 —
쐐기를 끝까지 뚫으면 접촉선의 두 끝이 면의 링 «위»에 놓이고(링 위 = 내부 아님), 가운데만 내부를
지나는 **현**이 된다. `combinatorics::segment_meets_face` 가 면의 **모든 링**을 그 선에 대고 읽어
교차점을 모으고, 선을 따라 정렬해 칸이 **밖/안/밖/안…** 으로 번갈아 간다는 사실로 안팎을 읽고,
열린 구간과 겹치는 「안」 칸이 있는지 본다. 갈래가 없으므로 구멍(모든 링을 한 자루 = 짝수-홀수)과
여러 칸을 저절로 처리하고, **끝점 판정을 흡수한다**(그래서 그 함수는 지웠다 — 실측 1,696쌍 전부
일치). 선 `q ∩ w` 의 `w` 는 **모서리 자신의 두 면**에서 얻는다: 링이 기록한 벽은 다발에서 `q` 를
가리킬 수 있고, 이름 삼중은 4평면에서 두 번째 평면을 안 담는다(31쌍 실측).

관문·실측은 `docs/dev-log.md` 의 같은 날 항목.

M5 `PolyhedralBoolean`은 **능력이 겹치는 두 메커니즘을 "공면 접촉 유무"로 배타 라우팅**한다(한 부울에 하나만 돈다).

- **일반 seam 엔진**(`general_boolean`, ray-cast `classof`): 공면 접촉이 **없는** 순수 transversal 전용. 공면은 door에서 정직 거절(`COPLANAR_PAIR`/`VERTEX_ON_FACE_PLANE`) — 평면 위 정점은 `classof`(3D winding)가 미정의라서다.
- **통합 공면 처리기**(`coplanar_result`): 공면 접촉이 **있는 모든 경우**를 통째로 만든다(벽·바닥·관통까지 self-contained; 두 엔진 합성 아님 — 합성은 on-plane `classof` 블로커로 불가). **일반 규칙(케이스가 늘어도 코드가 안 는다):** 각 면 F를 상대 solid의 **F-평면 단면(`section_of_solid`)** 에 대해 `coplanar_reconstruct`로 클립(`clip_face_to_section`). a-면=`a∖b`(Cut)/`a∩b`(Common)·b-면=`b∩a`(flip은 Cut만)·접촉면=상대 footprint 직접(mouth). 접촉면(π)은 단면이 퇴화하므로 footprint를, 벽/바닥만 단면을 쓴다.

**생존 규칙(접촉면 keep/flip = 연산 × 상대법선):** Fuse/opp=`inP⊕inQ`, Fuse/same=`inP∨inQ`, Cut/opp=`inP`, Cut/same=`inP∧¬inQ`, Common/opp=∅, Common/same=`inP∧inQ`. 볼록성 항이 없어 비볼록·다중루프를 half-space 가정 없이 `point_in_ring` parity로 정확 처리.

**pair 선택:** 평면-coplanar 면쌍 중 **풋프린트가 실제로 겹치는(포함 또는 경계 교차)** 쌍만 genuine 접촉(`footprints_overlap`) — 우연히 같은 평면에 있으나 떨어진 쌍(예: slot top ∥ 먼 포켓 바닥)은 배제. **면-평면 분기:** 면이 상대 면과 coplanar면 단면이 퇴화 → 단면 대신 **면 centroid `point_in_solid`**(disjoint라 strict)로 whole/drop.

**스코프 접촉면-flush(공변 — 접촉면 위 두 모서리가 한 선에 포개짐, cut-overhang의 본질):** ① **R0** — `section_of_solid`의 raw triple을 `plane_classes` canon으로 remap(안 하면 mouth·breach벽·b-벽의 공유 코너가 어긋나 조용히 non-manifold), fold 충돌은 `SECTION_TRIPLE_COLLISION` 거절. ② **R1** — 접촉면 위 attachment graze(section 코너가 벽 모서리 strictly 내부)를 거절 대신 F-edge-split crossing으로. ③ **R2** — 접촉면 chord 위 ∂P 정점을 covered로(삼킨/유지 코너), chord 끝점 4-평면은 `FLUSH_VERTEX_COINCIDENT` 거절. 이 셋으로 벽 네 구성(Middle/slab/Corner/Shorten)이 한 규칙으로, 코너 기둥이 canon-triple로 자동 용접.

**★ glue vs clip 경계(정직):** 위는 전부 **clip**(겹치는 풋프린트를 오림). **coincident**(동일 풋프린트 스택)는 **glue**(맞닿은 면 둘 제거 + 옆벽 splice + 인터페이스 정점 remap) — 구조가 달라 통합 안 함, `coincident_merge`로 별도 유지. **정직 거절(후속):** 회전 공면(toleranced 일치 판정), 다중-loop 단면(슬롯이 cavity 관통), 임의 fan degree-≥3 flush-shared-boundary(`turn_at` 미구현), 다중 genuine 접촉, P⊂Q 대칭. 전부 named 태그 — silent-wrong 0(DNA).

**은퇴:** clip bespoke 5경로의 **detector는 "진짜 공면 접촉인가(blind·not-pierces·convex)" gate로 유지**(transversal을 `general_boolean`으로 걸러 silent-wrong 방지)하되 결과는 `coplanar_result`가 만든다 — builder+헬퍼(proj2/clip/splice 계열)는 삭제.

### 6.2 병렬 — 판정만 병렬, 순서는 고정 (2026-07-29)

부울 시간의 90%가 **면당이 아니라 평면 클래스당** 독립인 두 루프에 있었고(트레이스→병합→분할, 셀→중첩→라벨→방출), 나머지 지배항은 모델의 정밀도를 **점마다** 읽는 `standard_for`였다. 셋을 `nacre-ops::par`의 헬퍼 둘로 병렬화한다. 옛 seam 엔진에도 같은 원칙의 병렬이 있었고 커토버(2026-07-21)가 엔진과 함께 들어냈다 — 되살린 것이지 새로 정한 것이 아니다.

**원칙: 판정만 병렬, 변이·순서는 단일 스레드 고정.** `assemble_fuse_cut`이 vertex handle을 면의 first-appearance로 배정하므로 `Store::push` 순서 = handle 정체성 = 재생 결정론(§2)이다. 그래서 모든 헬퍼가 **인덱스 순서로 수집**하고, f64/BigFloat를 병렬 reduce하지 않는다(항목 내부 계산은 순차). 스케줄에 의존하는 답은 한 번 틀린 답이 아니라 **매번 다른 모델**이다.

- **에러는 인덱스-최초.** 순차 루프는 가장 낮은 클래스에서 반환했다. rayon의 `collect::<Result<_,_>>()`는 단축평가라 *어느* 에러가 살아남을지 정해지지 않으므로 쓰지 않는다 — 전부 모은 뒤 순서대로 훑는다. 규칙은 `try_map_range` 안에 **한 번만** 산다.
- **별칭 표는 라운드 스냅샷으로.** 클래스는 4-평면 동시성에서 점의 별칭을 발견하고, 순차 루프는 각 클래스에 그때까지의 발견을 보여줬다. 이제 **라운드 시작 시점의 표 + 자기 발견**을 보고, 라운드 끝에 클래스 순서로 흡수한다. 같은 고정점에 이르는 이유는 스케줄이 아니라 `Aliases`의 성질이다 — 병합 대표가 **최소 원소**라 최종 분할이 병합 순서와 무관하고, 발견이 단조 누적이라 늦게 알면 라운드가 하나 더 들 뿐이다(실측: 4-평면 모델에서 순차·병렬 모두 **2라운드·별칭 36개로 동일**).
- **정밀도는 max 리덕션이라 안전하다.** `judge_precision`을 점당 `trial_bound`와 `precision_for`로 쪼갰다. 결합이 **최댓값**이라 결합적·정확하다 — 부분 *합*이었다면 재결합이 다른 수를 만들어 이렇게 못 나눈다.
- **대가: 거절 입력에서 일을 더 한다.** 순차는 첫 거절에서 나머지 클래스를 안 돌았고 병렬은 전부 돈다. 답은 같고 비용만 오르지만, **뒤 클래스에 잠복한 패닉이 도달 가능**해진다 — 거절 코퍼스가 그래서 관문이다.
- **임계값은 재서 넣는다.** `map_range`는 항목이 점 하나의 실현이라 작은 모델(면 12장 ≈ 점 36개)에서 분배가 계산보다 비쌌다(527→576µs). 측정된 크로스오버로 64를 골랐다(그 위는 fold 중간 부울이 순차로 떨어져 큰 쪽이 손해). **비용 손잡이지 답 손잡이가 아니다** — 양쪽이 같은 값을 같은 순서로 낸다는 것을 테스트가 못박는다.

**실측(14코어, 같은 커밋의 순차 빌드 대비):** 7 fins 64.9→20.8ms(3.13×) · 13 fins 235.6→59.9ms(3.93×) · 25 fins **904.9→178.9ms(5.06×)** · 축정렬 소형 1.2ms→517µs(2.37×). 대장 130행 비트 동일, OCCT 101/101.

**범위 밖: 브라우저.** 플레이그라운드는 `default-features = false`의 wasm 빌드라 이 이득을 못 받는다. `wasm-bindgen-rayon` + `SharedArrayBuffer` + COOP/COEP는 배포까지 얽히는 별도 작업이다. 되돌리는 길도 코드 밖에 있다 — 그 한 줄이 곧 순차 동작이다.

### 6.3 상대가 닿을 수 없는 평면은 소식을 나르지 않는다 (2026-07-29)

부울은 *"이 평면들 전부의 배열이 무엇인가"* 를 풀지만, 답하는 질문은 *"`b`가 `a`를 어떻게 바꾸는가"* 다. 증분 fold에서 둘은 거의 전부만큼 다르다 — 80핀 fold의 마지막 fuse는 평면 클래스 **169개를 배열하는데 새 핀이 닿는 것은 27개**이고, 나머지 142개에서는 **누적 솔리드가 이미 가진 면을 처음부터 다시 유도**한다. 실측으로 그 불리언 시간의 **59%**다.

그래서 **상대 솔리드와 떨어져 있음이 증명된 평면 클래스**는 배열하지 않고, 자기 솔리드의 면을 **배열 엔진의 어휘(평면 클래스 삼중항)로 다시 적어** 낸다. `nacre-ops::reuse` 한 모듈이고, `trace_result_faces`의 출력 타입은 그대로 `Vec<LocalFace>`라 **하류(조립·동일평면 병합·seam 표·검증기)는 한 줄도 안 바뀐다**.

- **★ 건너뛰기는 추정이 아니라 정리(theorem)다.** 분리는 `nacre-cip::orient3d_filter`로 **증명**한다 — 오차 한계가 0을 넘길 때만 부호를 답하고, 아니면 `None`. `None`은 "그냥 배열하라"이고 그건 이 모듈이 없던 시절의 동작이다. ⇒ **증명 실패는 속도만 잃고 답은 못 바꾸며, 톨러런스가 들어오지 않는다.**
- **★ 필터는 새 술어가 아니라 기존 술어의 절반이다.** `orient3d_judge`가 이미 f64 필터 → 공유-모션 정확 경로 → 상승의 3단이었고, 앞 두 단이 정확히 우리가 원하는 것이다. 그래서 그것을 `orient3d_filter`로 꺼내고 **판정이 그것을 호출하게** 했다. 오차 예산을 아는 곳이 둘이 되지 않으므로, 건전성은 검사할 성질이 아니라 **같은 코드라는 사실의 귀결**이다.
- **정점만 봐도 되는 이유, 그리고 그 전제.** 평면까지의 부호거리는 아핀이라 모든 정점이 한쪽이면 **볼록껍질 전체**가 한쪽이고, 평면 다각형은 (비볼록이어도) 자기 정점의 껍질 안에 있다. **이 논증은 모든 면이 평면일 때만 성립한다.** ★ **그 전제를 떠받치던 문장이 M6-2a에서 죽었다**: 옛 근거는 "`collect_planes`가 비평면 면을 거절하므로"였는데, C1이 원통 면을 표에 앉히면서 그 거절이 사라졌다. 지금 논증을 떠받치는 것은 **"원통 클래스가 하나라도 있으면 reuse를 통째로 끈다"**는 가드(`trace_result_faces` 진입)다 — 그래서 reuse가 도는 동안에는 "모든 면이 평면"이 여전히 참이다. 곡면이 들어온 뒤에도 이 지름길을 쓰려면 껍질이 아니라 **곡면의 경계**를 봐야 한다.
- ★★ **그 가드는 이중 임무다 — reuse를 만지기 전에 읽을 것.** ① 위 논증(평면 전제)을 살려 두고, ② `pass_through`가 **평면 클래스 단위**로만 면을 옮기므로(어떤 클래스가 PassThrough면 원통 측면이 **조용히 사라진다**) 그것을 막고, ③ **M6-2a의 밴드가 배열의 라벨을 읽을 수 있게** 한다(§6.4 — 라벨은 *배열된* 클래스에만 존재하고, PassThrough 클래스는 셀도 라벨도 만들지 않는다). 즉 원통 인구에 reuse를 되살리는 일은 성능 작업이 아니라 **밴드의 소속 판정을 다시 설계하는 일**이다.
- **★ 정점 이름은 솔리드가 이미 갖고 있다.** 배열은 정점을 `sorted3([w, f, third])` — 그 점에서 만나는 세 평면 클래스 — 로 부른다. 솔리드의 정점은 자기 주위 면들의 클래스를 정확히 만나므로, **술어 하나 없이 같은 이름이 나온다**. 넷 이상이면 그건 동시성이고 별칭 표가 이름을 고르는 문제라 재현할 수 없으므로 그 클래스는 배열한다.
- **규칙표에 기하가 남지 않는다.** 상대가 못 닿는 평면은 모든 점이 상대 **바깥**이므로 `keep`이 접힌다: `Fuse`는 양쪽 다 통과, `Cut`은 `a`쪽만 통과·`b`쪽은 비움, `Common`은 둘 다 비움. 양쪽이 다 면을 가진 클래스(동일평면 접촉)는 언제나 배열한다.
- **★ 두 패스를 모두 건너뛴다 — 그리고 그게 이 설계의 유일한 베팅이다.** 셀 패스(pass B)를 건너뛰는 것은 안전하다(면이 맞는지만 문제고, 그건 아래 차등 검사가 답한다). **추적 패스(pass A)까지 건너뛰면 엔진이 *아는 것*이 달라진다** — 추적되지 않은 클래스는 자기가 발견했을 동시성을 보고하지 않고, 커버리지 코퍼스에서 **별칭의 17%가 분리된 클래스에서만** 발견된다. "그 이름을 아무도 안 쓴다"(건드려진 클래스의 교차점이면 그 클래스가 스스로 찾는다)는 **논증이지 측정이 아니므로**, 믿지 않고 검사한다.
- **★ 차등 검사가 영구히 남는다 — 불리언 **전체** 수준에서.** 디버그 빌드는 같은 불리언을 `ClassReuse::{Off, Proved}` 두 번 풀어 **같은 면을 요구**한다(링 회전·면 순서에 무관한 정규형으로 — 배열의 링은 DCEL 순회가 시작한 곳에서, 솔리드의 링은 저장된 곳에서 시작하므로 그대로 비교하면 차이가 아닌 것을 차이로 읽는다). 클래스 단위로는 더 이상 비교할 배열이 없기 때문이고, 별칭 질문에 직접 답하는 것도 이 수준뿐이다. 참조 실행은 **자기 `Notes`** 를 받는다 — 증거는 `undecided_reject`가 읽는 부작용이라 두 번 쌓이면 안 된다.
- **건너뛰기는 클로저 *안에서* 한다.** 범위를 줄이면 안 된다 — 인덱스 공간이 유지되어야 `try_map_range`의 "최저 인덱스 거절 보고"와 방출 순서(=핸들)의 결정성이 지켜진다.

**두 검사 모두 공허하지 않다 — 음성 대조로 확인했다.** `flip`을 뒤집으면 커버리지 **4건**이, 분리 증명 없이 소유 클래스를 전부 통과시키면 **167건**이 즉시 실패한다.

- ★ **(2026-08-09, 열린 항목 0) 정확 도로의 산술 천장이 사라졌다.** `solid_points` 의 공유-모션 팔이 읽는 `three_planes_rat` 는 중간값 `i128` 오버플로로 거절하곤 했고, 십진 프레임 피연산자의 구성 코너가 정확히 그 인구였다(8/8 넘침 — 8/8 드는 점, 첫 실패가 솔리드 전체를 포기시키므로 프레임 피연산자는 이 모듈을 통째로 잃었다). 이제 그 함수가 정수 코어로 폴백해 답하고(`scalar`), 프레임 프리즘이 이 모듈에 복권됐다 — census 는 비트 동일(새로 발화한 계획의 산출이 배열 산출과 같았다). 벽-프레임 도구 자신은 혼합-프레임 꼭짓점이라 여전히 정직하게 declines(열린 항목 13/14 의 인구).

**실측(셀 시작 대비):**

| | 셀 전 | pass B만 | **양쪽** |
|---|---|---|---|
| 회전 fold 80 순차 | 9.66s | 7.69s | **4.78s (2.02×)** |
| 허브먼저 ring 80 순차 | 10.04s | 7.89s | **5.16s (1.95×)** |
| 축정렬 fold 80 순차 | 2.59s | 1.91s | **1.40s (1.85×)** |
| 회전 fold 80 **병렬** | 2.38s | 1.98s | **2.09s (1.14×)** |

대장 130행 **양쪽 피처 조합 비트 동일**, 커버리지 169, OCCT 101/101. 80핀 fold에서 클래스 방문 8800회 중 **7021회가 지름길**을 탄다.

**★ 병렬은 5% 나빠졌고, 그게 이 설계의 값이다.** 남은 일이 169개 중 **27개 클래스에 몰려서** 병렬 효율이 3.9×→2.3×로 떨어진다 — 총 일은 훨씬 더 줄지만 병렬은 그 폭을 못 따라온다. **순차를 택한 이유는 브라우저가 단일 스레드이고, 이 셀을 부른 불만이 브라우저에서 나왔기 때문이다.** 병렬 효율 회복은 정확성이 아니라 스케줄링 문제로 별도 항목이다.

### 6.4 원통의 벽은 평면 클래스가 아니라 밴드로 정해진다 (`nacre-ops::bands`, 2026-08-18, M6-2a C4b·K1)

배열 엔진은 **평면 클래스 하나씩** 답한다. 원통의 측면은 평면이 아니므로 그 어휘에 들어오지
않고, 그렇다고 새 배열이 필요하지도 않다 — M6-2a의 인구 게이트가 ⊥ 절단과 **r보다 멀리 떨어진 ∥
벽**만 통과시키기 때문에 다음이 성립한다(**균일-슬랩 정리**):

> 연속한 두 ⊥ 절단 사이의 **열린 원통 슬랩**에는 상대 솔리드의 경계가 전혀 없다(⊥면은 절단 자신,
> ∥벽은 게이트가 밀어냈다). 그러므로 상대 소속이 그 슬랩 전체에서 **균일**하고, 증인 하나가
> 밴드 전체를 정한다.

- **소속은 재지 않고 배열에게 묻는다.** 평면 배열은 원 경계마다 **원판 셀**을 만들고 `label_cells`가
  거기에 `[inA 위, inA 아래, inB 위, inB 아래]`를 쓴다 — "원 안쪽에서 이 평면 위/아래에 각 솔리드의
  재료가 있는가", 곧 밴드의 챔버다. 밴드는 경계면의 그 네 비트를 읽는다(양 끝을 다 읽고 일치를
  요구한다 — 다르면 균일-슬랩 정리의 전제가 깨진 것이라 정직하게 거절한다).
- ★ **왜 광선이 아닌가 — 이 마일스톤에서 가장 비싼 교훈.** 평면 엔진은 두 종류의 질문을 두 기계로
  답한다: *"이 경계 조각이 결과에 남나"*는 **라벨**(`emit_faces`), *"이 성분이 재료인가, 어느 성분
  안에 들어있나"*는 **광선**(`point_in_component` — 한 평면이 답할 수 없는 3D 질문). C4a는 밴드의
  질문을 **둘째로 분류**해 좌표를 받는 새 유리 광선 도로를 지었고, K1이 그것을 지웠다. 밴드는
  경계 조각이므로 **첫째**이고, 답은 한 층 아래에 이미 적혀 있었다. 라벨은 좌표도, 광선 방향도,
  기권 재시도도, 폭 천장도 없다 — 그리고 그래서 원통이 **양쪽 피연산자**에 설 수 있다(두-구멍 판).
  ★ **그 좌표 도로가 2026-08-21에 돌아왔는데, 그것이 이 경계를 뒤집는 게 아니라 확인한다**:
  돌아온 자리는 **둘째 질문**(`point_in_faces_rat` — 꼭짓점이 하나도 없는 성분이 캡 원판의 중심을
  증인으로 들 때)이고, 밴드는 여전히 라벨을 읽는다. 지운 이유는 도로가 틀려서가 아니라
  **밴드에게 틀린 모양**이어서였다는 것을 삭제 note가 이미 적어 뒀다.
- ★ **"z-범위 판정"은 반증됐다**: L-노치의 안쪽 코너는 모든 벽에서 r보다 멀지만 재료 밖이다.
  그 반례는 도로의 픽스처이자 밴드 패스의 픽스처로 남아 있다.
- **밴드 경계는 평면 배열이 실제로 방출한 원 경계**(+ 측면 자신의 rim)다. 임의의 ⊥ 클래스로
  자르지 않는 이유는 rim 공유다 — 캡 면의 원 루프와 밴드의 rim이 **같은 엣지 핸들**이어야
  닫힌-셸 가드가 사용 횟수 2를 본다.
- **조립**은 원통 b-rep 몸통(`cylinder_solid`)의 rim 기계를 그대로 쓴다: `OnSeam([측면, 평면])` 정점, 자기루프 rim
  엣지(`derive_edge_curve`가 원을 파생), 밴드의 seam은 `[측면, 측면]`. 성분 결합에는 **둘째
  규칙**이 붙는다 — 캡과 밴드는 Node를 하나도 공유하지 않으므로 **rim 키**로 잇는다.
- 감김은 규칙으로 유도한다: rim 원은 축 방향에 대해 CCW이므로 `sign(면의 바깥 법선 · 축)`이
  외곽의 `forward`를 주고 구멍은 그 반대다. 잠금은 **부피의 부호**(뚫린 상자가 `8 − πr²h`).
- ★ **룰링 도로 (M6-2 나머지, 2026-08-24)**: 게이트가 **기록-하고-통과** 팔을 얻었다 — 명확
  증명에 실패한 (벽, 원통) 쌍 중 **축이 그 벽 평면 위에 정확히 있는**(`n·o+d = 0`, 유리수
  부호) 것만 기록되어 통과한다. 기록된 쌍에서는 균일-슬랩 정리의 전제(「∥벽은 게이트가
  밀어냈다」)가 깨지므로, 밴드는 **잘린 원에서 끊기고**(`cut_rims`가 셋째 경계원), 양끝이
  모두 잘린 구간은 **θ-패널**(`Ring{호, 룰링, 호, 룰링}`)로 갈라진다 — 섹터별 챔버는 잘린
  원의 호별 라벨(`ArcLabels`)을 읽는다(라벨 원칙 그대로 — 광선 없음). 룰링 엣지의 담체는
  [벽면, 원통], 곡선은 ∥ 팔이 끝점에서 재생. 소비자: props는 θ-범위를 적분하고(닫힌 체인 =
  정확히 2π — 전-2π 공식으로 환원), tess는 열린 rim 체인을 감지 않고 병합한다. 부호 명단
  (패널 감김·디스크-쪽 선택자·룰링 턴·chord sense 전역)은 관통-보스 부피 오라클이 감시한다
  (항목별 뒤집기 실측 — 전부 빨강). **기록되지 않는 것**: 접선(d=r — 리프트 실측이 부피-정확
  제로-두께 핀치를 조립해 보였다; validate의 핀치 감지기는 원통에서 기권하므로 게이트가 정직한
  저지선)은 게이트에서 `WallMeetsLateral`로 남는다; 오프셋 교차(0<d<r)도 그렇게 남아 있었으나
  **칸 ③(2026-09-03)이 열었다** — 기록의 조건이 «축 위»에서 «r 안»으로. 캡이 상대 재료 안에 앉는 변형은 **기록되어 통과한 뒤** 하류(챔버)가
  `RulingBoundNotYet`으로 거절한다 — 거절 이름이 사다리의 것인 이유다.
- ★★★★★ **피연산자가 굽은 이름을 나를 수 있다 — 연쇄의 벽 (2026-08-26, 칸 ⑨).** 불리언의 **결과가
  다음 불리언의 피연산자**가 되면 트레이서가 평면 어휘로 못 읽는 기하가 들어온다. 실측한 경계:
  구멍→구멍·구멍→위 보스·위 보스→위 보스·위 보스→벽 보스·구멍→벽 보스는 **빌드되고**,
  **벽 보스 뒤에는 어떤 원통 작업도 안 된다**(위에 얹는 보스만이면 48개까지 붙는다).
  ★ **오늘 코퍼스에 원통을 연쇄하는 픽스처가 하나도 없었다** — 그래서 이 경계가 이 칸에서 처음
  측정됐고, 아홉 픽스처로 섰다.
  - **인구**(판 12×4×2 + 벽 보스 r=0.5, 면마다 실측): 판의 위·아래 캡은 **아크 1 + branch 정점 2**,
    벽 y=0은 룰링에 **둘로 쪼개져** 각각 branch 정점 2, 보스의 두 캡은 **온전한 원**(이미 처리 —
    `LoopRing::Circle`), 나머지 벽 셋은 손 안 탐, 합쳐진 옆면은 트레이서가 건너뜀.
  - **벽 둘**: ① 게이트의 `face_clears_footprint`가 **branch 정점**에서 `vertex_meet` 없음으로
    기권한다(첫 진단 「아크가 막는다」는 **반증**됐다). ② 그 뒤에는 **패닉**이 있었다 —
    `loop_triples`가 평면 전용 두 질문(옮겨 실은 wall 클래스 · 정점의 세-평면 이름)을 `ClassIx::plane()`
    으로 물어서. 그 접근자의 패닉은 **다른 40여 호출자에 대해 옳고**, 여기가 필터를 놓을 자리다.
  - **이 칸이 한 것**: 그 자리에 **그물**(`CurvedOperandBoundary`로 정직하게 거절)을 치고, 게이트는
    **답은 그대로 두고 이름만** 그것으로 바꿨다. 조건은 「outer 루프에 원통을 담체로 갖는 간선이
    있고, 루프가 간선 하나짜리(온전한 원)가 아니다」 — `loop_triples`가 클래스로 읽는 **같은 사실**을
    담체로 읽은 것이고(`Edge::surfaces`의 문서가 「두 면의 adjacency 답」이라고 말한다), **그 road의
    예외(온전한 원)까지 함께** 나른다(빼면 크로스 보어의 캡 면이 거짓인 이름을 쓴다 — 실측).
  - ★★ **게이트는 «열지 않았다».** 「그물이 덮는 만큼만 연다」는 등식을 정확히 쓰려는 시도가 **세 번
    다 예외를 만났다**: `OnSeam` 디스크 면은 그물 밖(그리고 그 자리는 「무리수 현이 어떤 면도 닿지
    않는 셀에 조용히 떨어지는」 기록된 상태다) · **subdivision**이 만드는 직진 branch 정점은 원통
    이웃이 없다 · 담체 등가는 두 면이 **같은 surface**일 때(동일평면 원 경계) 예외다.
    ⇒ 게이트 열기는 **도로를 짓는 칸의 마지막 걸음**이다(룰링 사다리가 그렇게 올랐다).
  - **다음 칸이 지을 것**: ① 피연산자 branch 정점의 **class-공간 번역** — 정점은 이름을
    `VertexDef::Branch`로 이미 들고 있고 `QuadRoot::canonical`이 인덱스 공간에 제네릭이지만,
    **한 클래스가 여러 surface를 담아서** handle↔class 대응을 새로 세워야 한다(그 타입의 문서가 이미
    「그 대응은 두 번째로 확립되어야 한다」고 경고한다). ② outer 루프의 **아크 어휘**. ③ 그다음 게이트.
- ★★★★★ **①이 섰다 — 피연산자의 굽은 링이 «서술»된다 (2026-08-26, 칸 ⑩).** 링의 코너는 이름을
  **유도하지 않는다**(정점이 `VertexDef::Branch`로 이미 들고 있다) — **handle 공간에서 이 배열의
  class 공간으로 «재진술»**하고, 담체는 결과 쪽이 이미 쓰는 `boolean::Wall`로 적는다(두 번째 어휘를
  만들지 않는다).
  - **재진술의 규칙 — 두 보정은 «같은 규칙»이다.** `Lo`/`Hi`는 `ℓ = n₁ × n₂` 방향의 순서이고,
    `plane_plane_cylinder`는 base를 `{n₁·x = −d₁, n₂·x = −d₂, **ℓ·x = 0**}`로 정한다 — `−ℓ`이
    **똑같이 만족하는** 조건이다. 그래서 **어느 한 법선의 부호를 뒤집으면** 두 평면도 base도 그대로고
    `ℓ`만 뒤집혀 **두 근이 자리를 바꾼다**. 쌍의 교환도 같은 이유로 그렇다.
    ⇒ **`ℓ`의 방향 반전이 홀수 번이면 root를 뒤집는다.** 교환분은 `NodeId::branch`가 이미 세고,
    이 다리는 **부호 둘**을 더한다. (원통은 보정 없음 — 축 방향의 부호는 곡면을 안 바꾼다.)
  - **대응과 부호를 «한 비교»가 답한다**: def의 각 평면 handle의 세계 계수가 어느 후보 클래스의
    계수와 **비례하는지** 보면, 비례하는 쪽이 그 클래스이고 **비례 상수의 부호가 곧 보정**이다.
    따로 구하면 두 곳에서 어긋날 수 있다.
  - **계측: 두 도로, 한 점.** 낸 이름을 `branch_point`로 **실현**해 그 정점 자신의 좌표와 맞춘다 —
    이름은 **이번** 불리언의 클래스, 좌표는 **이전** 불리언이 자기 클래스로 실현한 것이다. 근을
    틀리면 반대편 룰링의 점이 나온다. 벽 보스는 벽에 대해 **대칭**이라 코너 보스를 함께 쓴다.
  - ★★★★★ **오래된 빚 하나를 갚았다**: `assemble_fuse_cut`의 `QuadRoot::canonical` 호출은
    「미실행 — 클래스가 반대 순서로 도착하는 보스가 오면 픽스처를 준다」고 적혀 있었다. **그 인구가
    왔다**(두 번째 불리언은 클래스를 새로 짓는다): 순서 보정을 끄면 잠금이 **빨갛다**.
    ☑ 반면 **부호 보정은 오늘 인구에서 미실행**이다(끄면 초록) — 불필요가 아니라 미실행으로 기록한다.
  - 굽은 링은 이제 `plane_ring`에서 선다 — 코너면 `Branch`, 담체면 `CurvedWall`(서로 다른 명제).
    **게이트는 여전히 안 열렸으므로** 실제 불리언은 그 코드에 닿지 않고, census는 비트 동일이다.
- ★★★★★ **「결정은 `def`, 측정은 `cache`」 — 규칙은 코드가 이미 썼고, 자리는 훑어서 찾는다
  (2026-08-26, 칸 ⑪).** `WorkingCyl::cache`의 doc이 규칙이다(「the f64 twin … while **every
  decision reads `def`**」). 축 방향 **결정**을 realization으로 내리던 자리는 **셋**이다 — 다리의
  `Wall::Ruling::up`(정점 좌표를 읽었다)과, `cache`를 읽는 넷 중 결정인 둘(`band_loop`의 station과
  `axis_sign`). ★ `cache`를 읽는 나머지 둘(`branch_vertex_tol`의 공차, rim의 `centre`)은 **측정**이라
  규칙대로다 — **그대로 두는 것도 결과물**이다.
  **셋 중 둘**을 정확 서술로 옮겼다: `up`은 양 끝을 **끊는** 캡들의 축 매개변수를 비교하고,
  station은 접점이 앉은 **cut circle의 평면**이 축을 가로지르는 자리를 읽는다. 둘 다 **새 산술 0줄**,
  **같은 부등식이 무엇을 읽는가**만 바꾼다(`>`·`<=` 그대로 — 동률에 거절을 만들지 않는다).
  셋째는 아래대로 **반증**됐다.
  - ★★★★★ **`axis_sign`은 «같은 규칙»이 아니다 — 반증됨.** 그 자리는 「평면이 축을 **어디서**
    가로지르나」(위치, 계수 부호와 무관)가 아니라 「클래스 **프레임**이 축에 대해 **어느 쪽**인가」를
    묻는다. `world_rat`는 평면의 **이름**이라 그 방향을 안 나른다 — 실측 **436회 중 203회 불일치**,
    그리고 `world_rat · m`은 **전 호출에서 양수**(정규화된 상수)다. `WorkingPlane`이 이미 써 뒀다:
    「there is no such thing as *the* plane's outward normal」. 클래스의 유일한 방향 사실은 프레임이고
    `axis_sign`은 `plane.normal()` + `frame_sign`으로 **이미 그걸 읽는다**. ⇒ 코드는 그대로 두고
    **주석과 측정치**를 그 자리에 남겼다.
  - **계측의 규율**: census 비트 동일은 「출력이 안 변했다」이지 「두 경로가 매번 같은 답을 냈다」가
    아니다 ⇒ 매 커밋 **차등 탐침을 먼저**(둘 다 계산·어긋나면 패닉), 코퍼스 전량, **탐침 자체를
    부정 대조**, 그다음 f64 제거. `axis_sign`의 반증이 나온 것이 그 순서 덕이다.
- ★★★★★ **게이트가 «원통이 만든 코너»를 읽는다 — 근호가 하나뿐이라 정확하다 (2026-08-26, 칸 ⑫).**
  연쇄를 막던 것은 게이트의 `face_clears_footprint`가 branch 정점에서 멈추는 것이었다
  (`vertex_meet`은 2차 무리수 좌표에 유리 meet가 없어 거절한다). ⇒ **같은 세 질문을 그 서술로**
  묻는다: 점은 `line.base() + s·line.dir()`이고, 상대(반지름·span·평면 계수)는 전부 **유리수**이므로
  모든 양이 **`X + Y√c`** 하나 — M6-1의 부호 탑이 그대로 답한다. **근사도 평행 가정도 없다.**
  - ★★ **재진술이 없다**: `QuadRoot`는 정의가 **자기 두 평면을 자기 순서로** 말한 것에 대한 것이고
    여기서는 그것을 그대로 읽는다(칸 ⑩의 `ℓ` 보정은 class 표 순서로 **건너갈 때**의 값이다).
  - ★★ **거절이 원인별로 갈렸다**: `CurvedOperandBoundary`는 「뒤의 도로가 못 읽는다」가 **참인**
    자리(호 간선·`OnSeam`)에만 남고, 산술·이름 실패는 `CylinderGateUndecided`다. 이름이 참을
    말해야 다음 사람이 그 이름을 믿는다.
  - **실측**: 네 연쇄 픽스처가 게이트를 통과하고 `wall_faces_clear`가 **`Ok(false)`**(보스가 벽에
    서 있다) ⇒ **`d=0` record-and-pass**가 받아 `crossings`에 기록된다(그게 tracer의 룰링 팔이
    뛰는 조건). 새 거절은 한 층 뒤 **`TraceDeclined { BranchNode }`** — 정찰이 미리 잰 그대로다.
  - ★★★★★ **이 칸의 답을 재는 것이 크레이트에 하나도 없었다**: branch 코너의 답을 통째로 버려도
    스위트가 초록이었다(316). 그래서 게이트를 직접 부르는 잠금을 세웠고, 거기서 **더 정직한 사실**이
    나왔다 — 오늘 결정적인 것은 코너의 **값**이 아니라 **읽힌다는 사실**이다(판정이 다른 분리축
    `along`에서 나므로). **값이 결정적인 인구는 아직 없다**고 그 자리에 기록했다.
- ★★★★★ **M6의 남은 것은 «케이스»가 아니라 «능력»이다 (2026-08-26).** 「무엇이 막혔나」로 줄 세우면
  다음 칸이 그 케이스만 겨누게 되고 거기서 비스포크가 시작된다. 코드가 두 번 그렇게 말한다:
  `QuadRoot::**Double**`이 저장 변종이므로 「오프셋 `0<d<r`」과 「접선 `d=r`」은 **한 이차식의 근 2개/
  1개**이고, 엔진이 **만드는** 경계 종은 넷(`segs`·`circles`·`rulings`·`chords`)인데 피연산자에서
  **읽는** 종은 둘이었다. 그래서 목록을 이렇게 읽는다:

  | 능력 | 오늘 닫혀 있는 것 | 상태 |
  |---|---|---|
  | **A. 입력 어휘 = 출력 어휘** | ✅ 배열의 **«분할» 단계까지 닫혔다**(칸 ⑰). ★★★★★ **그리고 A의 남은 구간은 D다** — 아래 「A의 마지막은 D였다」 | 칸 ⑨~⑰ 완료. 트레이서는 굽은 피연산자를 통과시키고 분할 세 패스가 전부 그것을 받는다. 첫 정지는 `loop_winding`의 **`CurvedStraightRun`**이지만 그건 **표면**이고, 뚫고 보면 `label_cells`의 **`LabelConflict`**이며 그 원인이 D다 |

  ★★★★★ **하류 세 패스의 지도 — 셋 다 닫혔다 (2026-08-27).**

  | 패스 | 상태 |
  |---|---|
  | `split_at_crossings` (옛 「평면 전용」 overlay) | ✅ **칸 ⑯** — 분할점이 `Split`(평면 클래스 **또는** 어느 원통의 어느 근), 순서는 `order_located` 한 규칙 |
  | `split_circles` (호 분할) | ✅ **칸 ⑰** — 구간 판정은 `closed_contains`, 정렬은 `order_located`. 좌표를 읽는 것은 `segment_meets_cylinder` **필터** 하나뿐이고 양 끝에 좌표가 있을 때만 돈다 |
  | `split_rulings` (룰링 분할) | ✅ **칸 ⑰** — 방출부는 `split_segments_at`으로 호 분할과 **한 함수**이고, branch 끝 세그먼트를 건너뛰던 자리가 사라졌다(☑ 세그먼트 200·교차 392가 새로 검사되고 census는 비트 동일 — 옛 논증의 **결론**은 맞았다) |

  ☑ 그 칸이 딛은 사실은 ⑮에서 미리 쟀다: seated 링의 branch 코너 이름과 `circle_crossings`가 같은 점에 붙일 이름이 **같은 `NodeId`다**(32/32, 실현 거리 0). ☑ 칸 ⑰이 그것을 인구로 확인했다 — 끝점 중복 제거 **400회**, 그와 갈라지는 옛 거절(`CoincidentNodes`, 세 평면 이름의 끝점) **1회**.

  ★★★★★ **D의 첫 조각이 들어갔다 (칸 ⑲a, 2026-08-28).** 옆면 행이 **안쪽 루프(구멍)**를 들고 오고,
  자국이 구멍에 닿으면 **`DeclineKind::CylFaceHole`로 «면을 지목하며» 거절**한다. 열린 것은 없고,
  **조용히 틀리던 것이 정직해졌다** — 거절이 세 단계 뒤의 증상(`CurvedStraightRun`)에서 **거짓 문장이
  만들어지는 자리**로 옮겨갔다. ☑ 실측: 옆면 행 **287**개 중 구멍 있는 것 **9**, 서술 불가 **2**.
  census 두 프로파일 비트 동일 — 그리고 그건 재기 전에 **논증**이었다: 구멍 난 옆면은 **벽 보스
  결과물**에만 생기고 그걸 받는 작업은 이미 전부 거절된다.

  ★★★★★ **D가 열렸다 — 자국이 각도마다 답한다 (칸 ⑲b0·b1·b2·d, 2026-08-28).** 넷 다 관문 전량,
  census 두 프로파일 비트 동일(`78b06e50…`).

  - **⑲b0** 거절이 구멍의 ⊥ 간선 「앞의 둘」만 보던 것을 **전부의 min/max**로. ☑ 오늘 안쪽 루프
    **10개가 전부 직사각형**이라 답은 안 바뀐다 — 미래를 위한 보수화.
  - **⑲b1 옆면의 구멍을 «다른 모든 면과 똑같이» 이름 짓는다.** `loop_triples`에 「면 자신이
    원통인」 팔이 생기고, `trace_input`이 옆면의 안쪽 루프를 이름 짓는다(바깥 루프는 솔기 간선이
    self-adjacent이라 여전히 `span`이 말한다). ⑲a가 `planes.rs`에 넣었던 `holes`/`CylLoopEdge`와
    헬퍼 둘은 **소비자가 사라져 걷혔다**. ☑ 옆면 **577** 중 구멍 **5**, 이름 못 붙인 것 **0**
    (옛 도로는 **2**를 서술 못 했다).
  - **⑲b2 `circle_on_class`가 «구간 목록»을 답한다.** 걸음은 `ring_against_plane`을 **부른다**
    (새로 쓰지 않는다). 매듭 「어느 호가 구멍인가」의 답은 계획이 적었던 「담체의 바깥 법선 반대쪽」이
    **아니라 링 자신의 감김**이다 — `run_body_above`가 이미 쓰는 «재료는 진행 방향의 왼쪽». 그래서
    반대 면도, `edge → face` 표도, `model.adj` 캐시 문제도 **전부 사라졌고** 계단·U자 구멍까지 같은
    걸음이 덮는다. 교차점의 이름은 링 노드의 룰링을 이 클래스로 **재진술**한다(`ℓ = n₀ × n₁` 순서라
    ⊥ 둘의 「정렬순 × 축 부호」가 같으면 같은 root, 아니면 한 번 뒤집는다).
  - **⑲d 매끄러운 극점은 «곡률»로 감긴다.** `lo`는 hull 꼭짓점이라 코너가 없어도 국소적으로
    볼록하고, winding = 호 자신의 회전 = **`ccw · axis_up · frame_sign`**. `turn`의 `(Arc,Arc)`
    팔은 `0` 그대로 둔다 — 「회전」과 「곡률」은 다른 질문이다.

  ★★★★★ **그래서 벽이 배열 «밖으로» 나갔다.** 연쇄 원통 피연산자가 winding 걸음·셀 걸음·
  **`label_cells`를 전부 통과**하고 이제 **조립**에서 `OpenResultShell`로 선다. ☑ 매달린 간선을
  이름했다: 구멍의 rim 호 둘(z=0·z=2)과 그 현 둘이 **한 번만** 쓰인다 — 조립이 그 반대쪽을 주장할
  옆면 조각을 안 내놓는다. **그것이 다음 칸이고, 원인은 아직 진단되지 않았다.**

  ★★★★★ **첫 진단이 틀렸고 그 사실을 남긴다.** 「`bands::CylRow`도 옆면을 바깥 span 하나로 서술하니
  띠를 통째로 다시 내보내는 것」이라고 적었다(커밋 `2f89be3`의 메시지에도 그렇게 들어갔다).
  **두 근거가 다 무너진다**: (가) 구멍 난 옆면을 «엮는» 것은 `bands.rs`가 아니라 `boolean.rs`의
  `merge_curved_group`이고 그것은 **안쪽 링을 쓴다** · (나) ☑ 실측하니 바로 이 작업에서 **패널 도로가
  같이 발화한다**(보스 클래스에 패널 12 · 띠 8, 픽스처 넷 합). [[one-probe-cannot-rate-severity]]

  ★★★★★ **칸 ⑳의 조사가 그 자리를 정확히 채웠다 (2026-08-28) — 원인은 «셋»이고 둘이 서로
  상쇄하고 있었다.** 첫 융합과 둘째 작업의 간선을 나란히 재니 귀속이 끝났다:

  | 간선 | 첫 융합(옳음) | 둘째 작업 | 원인 |
  |---|---|---|---|
  | 룰링 x=1.5·2.5 | 벽 + **원통** | 벽 + **유령 벽** | (A) 룰링 자국이 전 높이 `Transversal` |
  | 현 | **간선이 아님** | 유령 벽만 = 1 | (A′) 평면 도로가 호를 직선처럼 훑음 |
  | rim 호 z=0·z=2 | 캡 + **원통** | 캡만 = 1 | (B) 패널의 「살릴까」가 한 비트 |

  ★ z=0의 호는 애초에 **하나뿐**이다(y<0 쪽은 경계가 아니다) — 「호 둘 중 하나만 매달렸다」는
  의문이 여기서 풀렸다.

  ★★★★★ **(A)는 칸 ⑳a가 닫았다** (`a248739`). 옆면이 축을 품은 벽에 남기는 자국이 이제
  **축 구간마다** 말한다 — 관통·**닿음**·관통. ☑ 걸음은 여기서 못 쓴다(**13/13 `AllOn`**: ⊥는
  곡선 하나, ∥는 룰링 «둘»이라 `side_of`가 못 가른다) ⇒ **담체로** 읽고 이유를 코드에 적었다.
  정렬은 **유리수**(코너가 이름하는 ⊥ 클래스의 `t`)라 θ보다 쉽다.
  ★ 감김 규칙(「재료는 진행 방향의 왼쪽」)이 세 번째로 쓰여서 **`material_theta_sign` 하나로
  뽑았다** — 그 부호를 뒤집으면 **⊥ 잠금과 ∥ 잠금이 «함께» 빨강**이다(같은 규칙이라는 증거).
  ★★ `world_rat`은 **이름이지 방향이 아니다**: 룰링의 **정체**는 그것으로 읽어도 프레임 무관하지만
  **딱지**는 안 된다 ⇒ `world_rat_sense`를 `side_of`에서 뽑아 쓴다.

  ★★★★★ **그래서 벽이 «앞으로» 옮겨갔다: `OpenResultShell` → `LabelConflict`.** 뒤로 간 것이
  아니라, (A)와 (A′)가 **서로 상쇄**하고 있어서 딱지가 맞아떨어져 보였을 뿐이다. 이제 배열이
  **자기가 불완전하다고 정직하게 말한다.**
  ★★ 그리고 칸 ㉑이 (A′)를 닫자 벽이 **다시 `OpenResultShell`**이 됐다 — 이번엔 **배열이 온전한 채로**.
  남은 것은 **(B) 하나**이고 그건 조립의 것이다.

  ★★★★★ **(A′)가 칸 ㉑이고, 그건 ⊥ 도로가 이미 배운 규칙이었다** (`4d39407`·`05e00d3`).
  `ring_against_plane`의 「선 위 연속 구간」이 **양 끝이 선 위인 굽은 간선을 가로질러** 이어진다
  — 직선이면 두 점이 선을 정하니 참이고 **호면 거짓**이다(선을 떠났다 돌아온다). 그래서 판의 캡이
  **자기 경계가 아닌 현 위에** `Graze`를 냈다.

  ★★★★★ **고침이 «걸음 안»인 이유: 소비자가 셋이고 셋 다 같은 결함을 갖는다.**
  `trace_transversal_face`(없는 경계에 graze) · `every_ray`(`PointOnRing` 거짓 거절) ·
  자기접촉 탐지기(`events`의 한 칸을 오판). ⇒ 걸음이 **구간을 departure에서 쪼개고**, 셋은
  **논리를 한 글자도 안 바꾼다** — 조각이 곧 그들이 늘 뜻하던 「구간」이기 때문이다.

  ★★★★★ **이 칸을 작게 만든 정리 — 패리티는 σ와 무관하다.** 원과 평면은 두 점에서 만나므로 호의
  **내부는 선을 다시 안 만난다** ⇒ 호는 통째로 한쪽(σ). 그러면 교차 수가
  `(s_a ≠ σ) + (σ ≠ s_b)`이고 **패리티는 어느 쪽이든 `s_a ≠ s_b`** = 오늘의 `flanks_differ`.
  ⇒ **네 번째 부호(σ)를 안 만들었다.**
  ★★ 다만 정리는 **개수**만 주고 **위치**는 안 준다 ⇒ 위치가 필요한 셋(뒤집어야 하는 departure 구간 ·
  양쪽이 departure인 조각 · 노드는 전부 위인데 간선이 떠나는 링)은 **`CurvedDeparture`로 거절**한다.
  ☑ 실측 **0회**. ★ σ를 «지금» 안 만드는 이유는 어려워서가 아니라 **잴 수 없어서**다 — 이 사다리가
  오라클 없는 부호를 이미 두 번 실증했다(⑲b2의 `frame_sign` 미실행, ⑳의 `body_above`는 부정 대조
  넷이 초록인 채로 **실제로 틀려 있었다**).
  ☑ **(E3-c, 2026-08-31) σ가 생겼다 — 잴 수 있게 되어서.** 재연산의 split-run 인구가 오라클이 됐고
  (반전 → `LabelConflict`), 이탈은 `EdgeMeet::Departs(σ)`의 가상 노드로 걷기에 들어가며
  `CurvedDeparture`는 세 자리 모두 삭제됐다(서술 없음은 `Unnameable`).

  ★★★★★ **술어는 통일하면 «버그»다 — 자취가 도로마다 다르다.**
  | 도로 | 클래스의 자취 | 「간선이 그 위인가」 |
  |---|---|---|
  | 평면 면 | **직선** `p ∩ q` | **호가 아닌가** |
  | 옆면(`hole_on_class`) | **원** `q ∩ 원통` | **담체가 `wc`인가** — 원 위엔 호가 얹힐 수 있다 |
  ⇒ 공유하는 것은 **규칙**(조각 쪼개기), **테스트는 각자의 것**. 통일했으면 ⑲b2가 연 것을 도로 닫았다.

  ★★ `Feature::Run::flank`이 **`Option<i8>`**이 됐다 — 호를 지나 돌아온 조각에는 「선 밖 이웃」이
  없고 그 값은 σ다. **타입이 그렇게 말하게** 한다.
  ☑ (E3-c) σ가 생기자 값이 항상 있어 **`i8`로 돌아갔다** — 같은 원리의 다음 걸음.

  ☑ **실측**: 쪼개진 구간 **30**(전부 `len=4`·굽은 간선 1·뒤집기 불필요) · `CurvedDeparture` **0**
  (세 자리 전부, 계측) · ★ 새로 배선한 두 소비자(`every_ray`·자기접촉)의 링에는 **호가 없다** ⇒
  그 배선은 **미실행** ·
  광선 캐스터가 `AllOn`으로 링을 건너뛴 횟수 **0**(잠복 무동작) · census 두 프로파일 비트 동일 ·
  회귀 가드(부피 정확) 초록 · 쪼개기를 끄면 빨강.

  ★★ **부호가 어디까지 잠겼나** (인자를 하나씩 빼서 실측). 호 방향 둘·`plus_t_is_above`·루트
  재진술·`ccw`·`axis_up` = **전부 빨강**. `frame_sign`과 옆면의 `orient_sign` = **초록(미실행)** —
  오늘 구멍 인구에서 둘 다 `+1`이고(큐보이드 면은 `Reversed`가 아니고, 보스의 옆면은 바깥을 본다),
  **bore의 옆면에 구멍이 생겨야** 돈다. 잠금 doc에 그대로 적었다.
  ★ **계획의 예측 하나가 반증됐다**: `ArcDir::axis_up`의 「코퍼스 전 클래스에서 `true`, 미실행」이라는
  **계측 주석이 만료**됐다 — 새 규칙은 그것을 발화의 절반에서 `false`로 읽는다.

  ☑ **⑲c(∥ 룰링 자국)는 미실행 — 픽스처 우연이 아니라 게이트로 닫혀 있다.** `rulings_on_class`는
  `TraceInput::crossings`에 실린 쌍에서만 발화하는데 그 필드는 프로덕션에서 **항상 빈다**(그 doc이
  그렇게 쓴다). 같은 `span` 거짓을 품고 있지만 오늘 아무도 그것을 묻지 않는다.

  ★★★★★ **A의 마지막은 D였다 — 증거 사슬 (2026-08-28, 전부 코드·실측).** 연쇄 원통 피연산자의 첫
  정지는 `loop_winding`의 `CurvedStraightRun`이지만, 그것을 뚫고 보면 진짜 벽은 `label_cells`의
  `LabelConflict`이고 원인은 **옆면 서술에 각도가 없다는 것**이다:

  | # | 사실 | 어디서 |
  |---|---|---|
  | ① | `CylFaceInfo`는 옆면을 **`span: [Rat;2]`(축 방향 구간)로만** 서술한다. **θ 필드가 없다** | `planes.rs`, grep 0 |
  | ② | `circle_on_class`는 「클래스의 z가 그 면의 span 안인가」만 묻고, 통과하면 **원 전체**를 자국으로 낸다 — **각도 범위를 안 본다** | `arrangement.rs` |
  | ③ | 그런데 판 모서리에 걸친 보스는 z∈[판] 구간에서 옆면의 **절반이 판 속에 묻혀** 있다 — 그 절반은 A의 껍질이 아니다 | 픽스처 기하 |
  | ④ | 그래서 자국은 「이 원 **전체**가 A의 껍질이고 A가 이 면을 관통한다」고 말하고, **절반에서 거짓**이다 | ☑ 실측: `Arc cyl=0 merged=[(A, Transversal{mat:1})]` |
  | ⑤ | 호를 건널 때의 뒤집힘이 반쪽만 맞아 한 바퀴 돌면 딱지가 어긋난다 | ☑ `LabelConflict`, 셀 5개의 딱지를 찍어 확인 |

  ★ **이 사슬은 ⑲b·⑲d가 끊었다(위).** ①②가 거짓이 아니게 됐고 ⑤의 `LabelConflict`는 **통과**한다 —
  즉 이 표는 이제 «해결된 진단»의 기록이지 현재 상태가 아니다.

  ⇒ **「θ-패널이 셀이 된다」가 D의 문장이고, 지금 없는 것이 정확히 그것이다.** 필드 하나 더하는
  문제가 아니라 **평면 쪽이 이미 간 길(클래스마다 셀 복합체)을 원통 쪽도 가는 일**이다.
  ★ 그래서 로드맵의 순서가 바뀐다: **A → D**가 아니라 **A의 마지막 구간이 D**다.

  ★★★★★ **연쇄 원통의 벽이 열렸다 (칸 ㉒, 2026-08-28, `8782927`).** 마지막 거짓 문장은
  **딱지가 두 질문에 답한다**는 것이었다. 딱지는 「어느 쪽이 재료인가」(**소속**)를 말하지 「이 옆면이
  그것을 경계짓는가」(**존재**)를 말하지 않는다. 구멍 난 옆면이 피연산자로 돌아오면 rim의 두 섹터가
  **글자 그대로 같은 딱지**를 들고, 패널 도로가 둘 다 살려 구멍을 메웠다. 답은 이미 자국에 있었다 —
  ⑲b2 이후 호마다의 **종류**(`Graze` = 면이 여기서 끝난다 / `Transversal` = 관통한다).
  `ArcLabels`가 딱지 옆에 그 기여 목록을 함께 나르고(`ArcLabel`, **한 호를 한 번 볼 때 같이 채운다**),
  `bands::face_spans`가 `keep` **앞에서** 존재를 먼저 묻는다.
  ☑ **실측**: 벽 보스 위의 둘째 원통 작업 **네 픽스처 전부 빌드** · `validate` **0**(= 매달린 간선
  0, 워터타이트) · 부피가 손으로 유도한 오라클과 **일치**(최대 오차 3e-14) · 존재 때문에 버려진 섹터
  **이 네 픽스처에서 4**(작업당 하나, 묻힌 반쪽 — 바이너리 전체로는 6, 아래) ·
  census 두 프로파일 **비트 동일** · 스윕 128 초록 ·
  perf 257.4s vs 255.7s(잡음) · 부정 대조 둘 다 빨강(게이트를 끄면 **옛 `OpenResultShell`로 정확히
  복귀**).
  ★★ **그런데 D의 본체는 아직 남아 있다.** 이 칸은 「정보가 이미 있고 소비자가 하나라」 놓은 계단이고,
  본체(옆면이 자기 차트에서 셀 복합체를 도는 것)를 하면 존재와 소속이 **한 질문**이 된다. 그때까지
  두 도로의 취급이 갈린다: 패널은 자국에게 묻고 띠(`chamber`)는 딱지만 본다(구멍이 거기 안 닿는
  논증은 `chamber`의 doc에).
  ★★★★★ **처음 적은 「미실행 셋」은 «네 픽스처만» 재고 쓴 것이라 틀렸다 — 전 스위트로 다시 쟀다.**
  ☑ 278번의 rim 읽기 중 **정확히 1회** 한 호가 옆면 자국과 **평면 면의 `Seated` rim을 함께** 들고
  온다(`bands::tests::the_gate_still_refuses_what_the_road_does_not_serve`) ⇒ **`Seated` 건너뛰기는
  방어가 아니라 «작동 중»이다.** 그것이 없으면 그 호의 답이 엉뚱한 면에서 온다.
  ★★ 미실행으로 «실제로» 남은 둘 — 기여가 **없는** 호(0회) · 기여가 **둘 이상**인 호(0회,
  `CylinderFaceUndecided`도 따라서 0회). 278개 전부가 옆면 자국 **정확히 하나**를 든다. ⇒ 그 둘은
  `face_spans`의 **단위 테스트**가 손으로 쓴 자국으로 진술하고, 잠금이 「1이 아니면 빨강」으로
  **들리게** 해 둔다(그날은 결함이 아니라 인구가 도착한 것).
  ★★ **존재 때문에 버려진 섹터도 4가 아니라 바이너리 전체로 6이다**(네 연쇄 픽스처 + 구멍 탐침
  테스트 둘 — 전부 같은 벽-보스 계열). ☑ **여섯 전부 「딱지였으면 살아남았을」 섹터**이므로 게이트가
  실제로 6개의 판정을 바꿨고, 그 여섯은 전부 이 칸 전까지 조립에서 거절되던 작업의 것이다
  (census 비트 동일이 그것을 뒷받침한다).
  ★★★★★ **「테셀레이션이 안 된다」고 적었던 것은 «틀린 도구로 잰» 것이다 (2026-08-28 정정).**
  `nacre_tess::to_obj`는 doc이 스스로 「bootstrap · 직접 평면 삼각분할」이라 적은 옛 작성기이고 굽은
  면을 아예 안 받는다(`NonPlanarFace`). 앱과 프로덕션이 쓰는 것은 `tessellate`이며 **네 결과 전부
  메싱된다** — ☑ 정점 910~1450 · 삼각형 1816~2896 · **비다양체 간선 0**(워터타이트). OBJ로 눈으로도
  확인했다. ⇒ **연쇄 원통 결과를 볼 수 없다는 제약은 없다.**

  ☑ **가 보고 확인한 것** (그 벽을 임시로 뚫은 탐침, 커밋 안 함): 가드를 때리는 링은 **단 둘**이고
  둘 다 **꽉 찬 원**(`n=2`, 한 실린더) · 세 인자 부호 공식(`ccw · axis_up · frame_sign[p]`)은
  `RingOrientation`(윤곽 수 + 오일러)을 **통과**한다 · 그러나 부호를 뒤집어도 **같은 자리**에서 서므로
  그 검사는 여기서 **눈이 먼다**(쌍둥이 반쪽이 같이 뒤집힌다).

  ★★★ **별개의 빚 하나를 세어 두었다**: `loop_winding`은 사전순 최소 **노드**를 hull 정점으로 쓰는데,
  호가 있으면 **노드 집합의 hull ≠ 영역의 hull**이다(호가 모든 노드보다 바깥으로 부풀 수 있다).
  ☑ 전 스위트 **183,306 링 중 274개**가 오늘 **영역의 극점이 아닌 자리에서** 회전을 읽는다(스위트는
  초록이므로 답은 맞고 있으나 **전제가 거짓**이다). 답은 그 가드의 doc이 이미 지정했다 — 「**영역의
  극점**에서 읽어라, 호의 내부일 수 있다」. ☑ 다만 **연쇄 픽스처에서는 0회**라 최전선이 아니다.

  ★★★★★ **셀 단계에 남은 「평면 전용」 자리 — 지도 (2026-08-27 grep 실측).** `CurvedStraightRun`
  이 먼저 막고 있어 대부분 아직 안 닿지만, 능력 A를 셀 단계까지 밀면 이 다섯이 차례로 나온다:
  `combinatorics::point_on_ring`(양 끝의 `.class()`를 요구 → `RingNaming`) ·
  `combinatorics::point_in_component`·`cmp_key`·`Chart::ring`(전부 `node_coords_rat`) ·
  `arrangement::point_in_mixed_ring`·`node_in_circle`(같음). ☑ 분할 단계에는 **하나도 안 남았다** —
  거기서 좌표를 읽는 유일한 자리는 `segment_meets_cylinder` 필터이고 그것도 선택적이다.

  ★★★★★ **칸 ⑰이 만난 새 벽은 「미실행이라고 스스로 적어 둔 가드」다.** `RejectReason::CurvedStraightRun`의 doc이 *"an unfired guard by construction … it fires the day that stops being true"* 라고 쓴다. `loop_winding`은 링이 **직진하는** 노드를 지나쳐 걸어 극점의 회전을 읽는데, **한 원의 두 호는 접선 연속**이라 걷기가 그것도 지나친다. 원이 호로 잘리는 인구가 생긴 것이 그날을 만들었다.
  ★★ **⑲d가 그 자리를 열었지만 답의 모양은 달랐다.** 「영역의 극점에서 회전을 읽는」 일반 규칙이 아니라 — **노드의 극점에서 회전 대신 «곡률»을 읽는** 것이다: `lo`는 hull 꼭짓점이므로 코너가 없어도 국소적으로 볼록하고, 그러면 부호는 호 자신의 회전이다. **위의 274개 빚(노드 hull ≠ 영역 hull)은 그대로 열려 있다** — 그건 `lo`가 «영역의 극점이 아닌» 경우이고, ⑲d는 `lo`가 극점인 채로 **매끄러운** 경우를 답한 것이다.

  ★★★★★ **칸 ⑯이 성능을 처음으로 «축»으로 삼았다.** `split_at_crossings`는 불리언의 **12~24%**
  (회전 fold에서 24%, collect 320만 trip)이고, 그 속도는 hoist 둘(`ImplicitPoint`의 Cramer 캐시 ·
  `dir_sign`)에서 나온다. 그래서 순서 규칙의 hoist를 **두 번째 구현이 아니라 인자**로 만들었다
  (`Located`/`OnLine`/`order_located`). 능력을 연 값은 **불리언의 약 1.4%**로 실측·기록.
  ☑ 계기 규율: **trip 수는 결정적이라 날카롭고, 시간은 세션 간에 표류하므로 A/B는 같은 세션에서**
  (`git stash`). 아침 기준선을 몇 시간 뒤 비교에 쓰면 안 된다 — 안 건드린 단계가 +9% 움직였다.
  | **B. 벽∩측면은 «한 이차식»이다** | 오프셋 `0<d<r` · 접선 `d=r`(`WallMeetsLateral`) | ✅ **오프셋은 칸 ③(2026-09-03)으로 닫혔다** — 게이트의 clearance 한 호출이 유일한 철자, 배열은 기록만 믿는다. 남은 것은 **접선**이고, ★ **그것은 C가 아니라 배열의 사다리다**(칸 ⑤ 실측, 2026-09-04 — 아래 「B가 쪼개지는 곳」의 정정) |
  | **C. validate의 비다양체 판정이 곡면을 본다** | 접선의 두께-0 접촉 | 불리언이 아니라 **계측**의 일반화. ★ **정정(칸 ⑤)**: C가 B의 접선을 막고 있는 것이 **아니다**(B의 벽은 배열이다), 그리고 **칸 ㉓의 물음은 답이 나왔다** — 그 두 결과는 유효한 솔리드다. C의 진짜 인구는 **두께 0**이고 그것은 B가 열려야 생긴다 |
  | **D. 원통도 자기 차트에서 arrangement를 돈다** | ✅ **D5(2026-09-02)로 닫혔다** — 방출기가 차트의 «영역»(방출 셀의 연결 성분과 그 경계 사이클)을 걷고, 띠/패널 디스패치·D4 병합 패스는 삭제됐다. 남은 것은 조립의 슬릿 걷기(`band_loop`)뿐 | 아래 「여덟째 계단」 |
  | ─ 실린더 쌍(`CylinderPairContact`) | ✅ **면이 비킨 쌍은 열렸다(칸 ⑧, 2026-09-05)** — 게이트가 옆면의 축 구간으로 «두 클래스가 면을 공유하지 않음»을 증명하면 통과한다(사용자의 십자 스터드). ★ **넷 중 넷(칸 ⑩, 2026-09-06)**: 같은 솔리드의 쌍은 묻지 않고, 평행 쌍은 «같은 곡면인가» 뒤 같은 문(구간·호), 비스듬 평면 팔도 옆면마다 묻는다; 옆면의 서술은 차트의 발자국(축 구간 × 각도 범위, 림 호에서 유리수로) | 면이 **만나는** 쌍만 **진짜 새 기하**(quartic) = M6b 마일스톤 |
  | ─ 스케치 원·호 어휘 | ✅ **열렸다(칸 ⑨, 2026-09-05)** — `Edge2d::Arc`(끝점·유리수), 4분원 배수의 호와 온전한 원이 프리즘 빌더에서 원통 조각 벽·`Branch` 꼭짓점으로 선다; `cylinder()`는 원 스케치 + 돌출의 설탕이 되어 원시체 연산이 은퇴했다. 남은 것: 임의 각도(점의 이름)·호–호 접합. ~~같은 솔리드의 동축 원통 둘(링·와셔)의 불리언~~ → 칸 ⑩이 열었다(와셔·슬롯·필렛 판이 상자와 셋 다 짓는다) | 기능, 엔진과 직교 |

  ★★★★★ **B가 쪼개지는 곳 (2026-08-28, 코드가 직접 말한다).** `WallMeetsLateral` 하나에 **세 가지**가
  들어 있고, 그 doc이 각각을 이렇게 쓴다:
  - **접선**(`d = r`) — 코드가 *「들어올리면 **부피 정확한 solid가 조립된다** … 막는 것은 C다」*라
    적었고 이 표는 그래서 «접선과 C는 한 항목»이라 결론지었다. ❌ **칸 ⑤(2026-09-04)이 실측으로
    반증했다**: 게이트의 접선 팔을 들어내면 12칸(보스 안/밖·보어 슬랩 × 3연산) 중 **한 칸도
    조립되지 않는다**. 그 문장도 칸 ③의 오프셋 문장과 같은 **D5 이전 측정** 위에 서 있었다.
    막는 것은 배열이고, 사다리를 한 계단씩 채워 가며 잰 실측은 이렇다 —
    ① `chord_on_class`가 `CylinderMeet::Pair`를 요구(접선은 현이 아니라 **한 점**) →
    ② `crossing_on_ruling`도 `Pair`를 요구(접선은 룰링이 **하나**) →
    ③ 걷기가 그 노드에서 간선을 못 세운다(`UnorderedEdges`) · 보어 쪽은 스캔의 노드 정렬이
    겹친 두 branch 노드를 못 센다. ①②는 `Tangent` **팔 한 줄**씩이고, **③이 진짜 어휘 결정**이다
    — `RulingCarrier::side`의 `±1`이 다섯 구조체의 동일성 키인데 접선은 룰링이 하나다. ⇒
    **접선은 B의 마지막 인구이고 C와 별개다.**
  - **오프셋**(`0 < d < r`) — 코드는 「**cell D**의 비유리 룰링선」이라 썼다. ★★ **그 「cell D」는
    옛 칸 레터링**(offset = D · tangent = E)이고 **이 표의 「D. 원통 차트 arrangement」와 다른 D다.**
    ✅ **칸 ③으로 닫혔다** — 그 문장(과 «`OpenResultShell`로 돌아온다»)은 D5 이전의 측정이었고,
    무리수 룰링은 배열이 이미 다루고 있었다(아래 문단).
  - 게이트가 못 치운 면 — 충분조건이 약한 것. 별개.

  ★★★★★ **칸 ㉓ — 삼각분할 도로가 «둘»이었고, 살아남은 도로가 2174 중 2를 못 그린다
  (2026-08-29, `b19aebd`·`fe835a1`·`27baaac`).**

  - **부트스트랩을 지웠다.** `to_obj(&Model)`은 M1 작성기이고 반쪽 간선의 **시작 정점만** 읽어
    **호를 현으로 바꿨다**. 크레이트 헤더가 스스로 *「Two paths coexist **during M3** … kept until
    the existing planar callers migrate」*라 적어 둔 임시 상태였고, 이행이 **M6에 세 마일스톤 늦게**
    왔다. ☑ 지운 것 셋(자유 `to_obj`·`triangulate_polygon`·`TessError::NonPlanarFace`) · 호출자
    **9곳**이 `tessellate(..)?.to_obj()`로 · 평면 차트 규칙(newell → drop_axis → 방향 복구)이
    **한 벌**만 남음(`planar_chart`). census 두 프로파일 **비트 동일**.
    ★ **경고는 dev-log에 두 번 적혀 있었는데도 세 번째로 밟혔다** — 문서로는 못 막았고, 함수를
    지우니 막혔다.
  - **`TessError::SelfTouchingBoundary`.** 어떤 정점도 자기와 무관한 경계 조각 위에 있으면 안
    된다 — `monotone::link`가 이미 거절하는 «조합적» pinch(인덱스 공유)의 **기하학적 쌍둥이**이고,
    `link` 바로 다음에 선다(`prev`/`next`가 인접 제외를 공짜로 준다). **새 술어 0 · tolerance 0.**
  - **계기**: `boolean_with_report`의 유일한 성공 출구에서 결과를 메싱해 본다. ☑ **2174건 중 2건**
    실패, 둘 다 접선 — 보스 발자국이 보어 rim에 접하고(`a_segment_tangent_to_the_rim_builds`),
    보스 rim이 판 윗모서리에 접한다(`a_turned_boss_tangent_to_the_plate_top_builds`). 그 면의
    **내부가 한 점에서 집힌다.**
  - ★★★★ **그 둘이 «유효한 솔리드인가»는 이 칸이 답하지 않는다 — 능력 C의 물음이다.** `validate`는
    초록이고 부피도 1e-9 안에서 정확하지만, 면의 내부가 한 점에서 집힌다는 것은 그 점에서 면의 두
    쪽이 만난다는 뜻이고 **표면의 비다양체 점일 수 있다**. `validate`는 못 본다 — 거기엔 위상 정점이
    없다. ⇒ 이름은 「이 분해가 못 그린다」까지만 말한다.
  - ✅ **답: 유효하다 (칸 ⑤, 2026-09-04).** 접점 둘레 작은 구 위에서 경계를 따라가면 **링크가 원
    하나**다 — 집힌 면의 두 로브가 이웃 **곡면**(보어 옆면·보스 옆면)을 돌아 이어진다. 픽스처 1의
    구간표(`p = (11,10,5)`, 보어 벽 `u ≈ −v²/6`): 보스 면 반원 → 플랩 미소호(v ≈ −ρ) → 보어 벽
    반원 → 플랩 미소호(v ≈ +ρ) → 닫힘. 재료도 국소적으로 한 덩이다(면적을 가진 접촉면으로 이어짐).
    **독립 확인**: OCCT가 같은 피연산자로 같은 합집합을 만들어 **부피·면적이 일치**하고
    (32.785398 / 3906.628331), **면 수도 nacre와 같으며**(8 / 12) 접선을 ε 비킨 쌍둥이와도 **같다**
    — 상용 커널도 몸통으로 받고 집힌 면을 쪼개지 않는다(`nacre-oracle`의 두 `#[ignore]` 측정).
    ⇒ **표면은 2-다양체이고 집힌 것은 면이다.** 잠금: `bands`의 `the_touch_is_a_slit`(접점에 위상
    정점 0 **∧** 닿는 원이 여전히 온전한 `[v,v]` 림).
  - ★★★★★ **χ는 이 물음의 증거가 아니다.** 두 결과의 χ는 0·2로 짝수이고 «핀치 하나는 χ를 홀수로
    만든다»(`check_result_topology`)가 참이지만, 그 문장은 **정점인 핀치**의 것이다. 접선은 정점을
    만들지 않으므로 `euler_counts`가 그 점을 아예 못 세고, 짝수 χ는 «핀치 없음»과 «위상 없는 핀치»
    **양쪽과 양립**한다 — 판별력 0. 짝수 χ를 증거로 읽으면 거짓 규칙이 된다(칸 ⑤이 한 번 그럴
    뻔했다).
  - ★★ **그러면 처방은 무엇인가 — b-rep이 아니라 메시다.** 재봉합(§«figure-8»)은 한 링이 한 정점을
    두 번 지나게 만들어 `monotone::link`의 **인덱스 공유 핀치**로 다시 거절되고(좌표 접촉 → 인덱스
    접촉), 그 다섯 철자는 칸 ㉓이 이미 `DegenerateRing`으로 실측했다. ★ 정확히는 **반복 인덱스** 철자에
    대한 문장이다 — **쌍둥이 인덱스**(같은 좌표, 다른 인덱스) 철자는 `link`를 통과하고 `self_touch`·
    `sweep`의 동점 가드·`emit`의 넓이 0이 막는다(칸 ⑦c가 그 셋을 푼다). 면 분할은 **유효한 솔리드의
    b-rep을 바꾸는** churn이고 OCCT도 안 한다(위 면 수). ⇒ 남은 능력은 아래 「남은 능력」이 이미
    이름 붙인 **스윕의 기호적 순서(SoS)**다. 그리고 그 손실은 작지 않다: `tessellate`는 **모델
    전체**를 돌며 첫 실패에서 멈추므로 **집힌 면 하나가 그 세션의 모든 몸통을 지운다**(소비자는
    플레이그라운드; STEP 내보내기는 tess를 안 타 무사하다).
  - ★★ ~~**남은 능력**: 「집힌 경계를 그린다」~~ → ✅ **열렸다(칸 ⑦c, 2026-09-05).** ☑ 싼 철자
    다섯(정점 공유 · 간선 위 · 쪼갠 간선 · 반복 인덱스 병합 · 쌍둥이 정점 병합)이 `DegenerateRing`이던
    것은 맞지만, 그 이유는 «하나가 부족해서»가 아니라 **셋이 함께여야 해서**였다: ① 접점(원의 샘플)을
    **공유되는 직선 간선의 폴리라인에 끼워** 양쪽 면이 같이 보게 하고(pre-pass, 정점을 안 민다),
    ② 두 링을 그 점에서 **인덱스 수준으로 이어 붙여**(다리 길이 0) 같은 좌표의 쌍둥이 인덱스를 남기고,
    ③ 스윕의 동점을 **쌍둥이에 한해 인덱스 인식 기호 순서**로 정한다(`polygon/sos.rs` — 접선 방향
    섭동, 1차 항은 `Expansion`으로 정확, 좌표는 안 움직임). 아래 칸 ⑦b·⑦c 항목.
  - ★★ ~~**남은 능력**: 「원통 쌍의 게이트가 면을 읽는다」~~ → ✅ **열렸다(칸 ⑧, 2026-09-05).**
    사용자의 `fuse(cuboid, z-스터드, y-스터드)`가 둘째 불리언에서 `CylinderPairContact`로 거절되던 것 —
    게이트의 원통–원통 팔이 **무한 곡면**의 거리(축 거리 vs 반지름 합)를 물었고, 축이 원점에서 만나므로
    거리 0이었다. 배열이 필요로 하는 명제는 «두 클래스가 **면**을 공유하는가»이고, 같은 게이트의
    평면–원통 팔(`face_clears_footprint`·`lateral_spans`)과 평면 배열(`trace_transversal_face`의 다각형
    클리핑)은 이미 면 단위로 답하고 있었다. ⇒ 곡면 거리로 못 비킨 비평행 쌍은 **옆면이 자기 자신을
    대변한다**: 한 클래스의 옆면 축 구간 전부가 상대 옆면의 **도달 범위**(`lateral_reach` — 사영 구간
    `d·o + s·(d·m)` ± 반경 도달 `r·|d⊥|`, 제곱으로 정확, 근호 없음)와 서로소면 두 클래스는 면을 공유하지
    않고 배열은 볼 것이 없다(어느 방향이든 하나면 충분). 수직은 그 일반식에서 `s` 항이 사라지는
    경우이지 분기가 아니다. 아래 칸 ⑧ 항목.
  - ★★ **스케치가 원과 호를 그리고, `cylinder()`는 설탕이다 (칸 ⑨, 2026-09-05).** 상자는 «네 점
    프로파일 + 돌출»인데 원통만 별도 도로(`Operation::Cylinder` → `add_cylinder_exact`)였다 — 스케치에
    원이 없어서. 이 칸이 그 자리를 채웠다: 커널 진실은 `Arc { center, start, end, ccw }`(전부 유리수,
    `start == end` = 온전한 원, 회전각은 커널이 모르는 말 — 끝점 표현이 더 넓다), 앱 문법은 펜
    `arc({ center, sweep })`·`circle({ center, r | d })`·`{ fillet: r }`(직각·축 정렬만, 합 == 변이면
    반원으로 병합)·뭉치 `arc({ center, start, sweep })`, `close()`는 펜이 시작점이면 선을 안 그린다.
    ★ 설탕 주장은 **먼저 잰 뒤** 원시체를 걷었다: 원 프로파일의 Extrude 와 원시체가 위치 정준으로
    비트 동일(정점 비트·`OnSeam` 둘·평면/원통 진술·면·간선·부피 비트·삼각형 집합, 첫 시도 적중).
    ★ 실측이 가른 것: 옆면 방향은 «호의 ccw == 링의 감김»(계획)이 아니라 **`ccw == 스윕이 +w`**
    (구멍 원의 옆면이 계획대로면 뒤집힌다); 감김의 `a + b·π` 는 16각형에서 i128 이 넘쳐 BigInt 로
    (scalar 에 두고 ops 는 BigInt 무의존); 시계 방향 호의 캡 간선은 `[A,B]` 반시계 규약에 맞춰 정점을
    바꿔 저장하고 반간선이 되돌린다. ★ census 가족 `arcprofile`(299행, 옛 287행 무변): D 보스 × 상자
    셋 Ok, 구멍 판 fuse/cut Ok(common 은 상자 면이 구멍에 접해 `SelfTouchingResult` — 픽스처의 접촉),
    **링(annulus) × 상자는 셋 다 `CylinderPairContact`** — 동축 쌍을 평행 팔이 거절한다(칸 ⑧이 남긴
    핸들 인턴, 흔한 모양에서 바로 보임), 슬롯 × 상자는 `TraceDeclined { OuterRing }`(접선 룰링 —
    `ruling_side` 가 0). ★ 사용자의 첫 조립체(2026-09-06)가 이 둘에 **비스듬 팔**까지 세 벽을 한 번에
    만났다 — 필렛 판 × 슬롯 창 × 삼각 거싯: 평행 쌍(축 거리 ≤ r₁+r₂)·접선 룰링(먼 상자와도)·거싯
    빗변 vs 구멍(`ObliqueCylinderCut`, 비켰음을 묻지 않음). 셋 다 닿지 않는 면 사이의 거절이다.
    다음 칸의 입력이다. 실측·대화·함정은 dev-log 「칸 ⑨」.
  - ★★★★★ **게이트는 모든 자리에서 «면»을 읽고, 접선 룰링은 이름을 얻는다 (칸 ⑩, 2026-09-06).**
    칸 ⑨의 세 벽은 셋 다 **닿지 않는 면** 사이의 거절이었다. 뿌리 셋: 원통 쌍 루프가 **같은 솔리드
    안의** 쌍까지 곡면 거리로 물었다(유효 솔리드의 두 면은 교차하지 않으므로 증명이 필요 없는 쌍 —
    이제 클래스 표가 소속을 들고 교차 솔리드 쌍만 묻는다; 평행 팔의 첫 물음은 «같은 곡면인가»이고
    그 뒤는 비평행과 **같은 문**); 옆면의 서술이 축 구간뿐이라 각도 범위를 몰랐다(발자국
    `Footprint { span, theta }` — 림 호의 두 끝을 유리수로 읽고, 반경 도달의 극값은 호가 `d⊥`를
    품을 때만 근호, 평행 쌍의 구간 겹침은 공통 단면의 두 **호**가 점을 공유하는가); 비스듬 평면 팔이
    비켰음을 묻지 않았다(옆면마다 `lateral_reach(d = n)`). 그리고 트레이서의 어휘: 접선 룰링은
    `side = 0`이라는 **이름**이다 — 캐리어를 만드는 두 자리에서 코너의 근 `Double`로 유도하고
    `ruling_side`(부호 술어)는 무변. ★ 그 이름이 서자 **여덟 계단**이 드러났다: 직선–호 접선 이음은
    평면 배열에서 «반대 방향 = 반바퀴, 같은 방향 = 곡률 물음»(`tangent_travel_agrees`) · 부분 림의
    덮이지 않은 호는 간선이 아니다(`split_circles`, 차트 `End::Uncovered`) · 앉은 면이 접선 룰링
    조각을 스스로 낸다(`SegKind::Tangent { axis_above }`, 마스크 = 축 쪽 비트) · 셋째 평면 위의
    코너는 세 평면 이름(`third_on_l`) · 증인 둘(유리수 branch 코너, 두 근 사이의 유리수 점) · 게이트는
    면이 닿지 않아도 클래스가 옆면을 **자르면** 기록. ★ 거절이 가리던 배열의 구멍 셋(S1: 두 디스크
    미중첩의 조용한 `Ok(2)`·원형 바깥이 구멍을 안 읽음·병합이 원형 바깥을 못 말함)은 관대화의
    잘못이 아니라 그 인구가 처음 도달한 자리였다. ☑ 실측: `arcwalls` 30행 중 `p1p2`·`p2p3`·
    `fillet15-far`·`slot-far`·`holes-gusset`·`stacked`·`bushing` Ok, 사용자의 `p1 + p2 + p3`(거싯
    `x = 1.6`) **한 몸**·부피 = 부분 합, 와셔 × 상자 Ok, 슬롯 × 상자 common·cut Ok(fuse는 픽스처의 접촉), 반원 D 보스·슬롯 pad/pocket Ok(칸
    ⑨의 «오늘의 이름» 잠금 둘이 양성으로), 칸 ⑨ 이전의 287행 무변. 남은 정직 거절: ~~거싯 옆면이 필렛 축을
    정확히 지나 접선 룰링을 두 평면과 원통이 **공유**하는 치수(사용자의 `x = 1.5`)~~ · ~~접선 벽과
    동일평면인 상자 면이 접점을 지나 뻗음(`rrect-box`)~~ → 칸 ⑫이 열었다 · ~~4-fold의 `+p4`(`FourPlane`)~~ →
    칸 ⑪이 열었다. dev-log 「칸 ⑩」.
  - ★★★★★ **점의 이름은 한 함수가 짓는다 — 네 평면 정점의 여섯 얼굴을 한 규칙으로 (칸 ⑪, 2026-09-06).**
    사용자의 조립체가 스크립트 순서로는 넷째 fuse에서 거절되고 슬롯 판을 마지막에 붙이면 짓는다는
    것이 입력. 기록을 읽으니 같은 명제 «점의 이름 = 그 점을 지나는 평면 집합의 함수»를 07-27부터
    **여섯 번** 다른 얼굴에서 지역적으로 고쳤고(배열의 run 정점 → 별칭 표 → «면 넷 정점은 대장에
    없다» → 결과 정의의 종속 건너뛰기 → 선을 담은 넷째 평면 → 끝점 삼중의 교집합), 07-27이 «점을
    지나는 평면 전부를 알려면 O(V·P) — 순환»이라 미룬 전제가 피연산자 정점에 거짓이었다(위상이
    인접 면을 안다). ★ 계기가 거절 뒤의 **조용한 오식별**을 찾았다: 피연산자 링이 정점을 면마다
    이름 지어 네 평면 정점이 이름 넷(그중 하나는 세 평면이 선을 공유하는 종속 삼중)으로 들어오고,
    판정기가 그 «점»을 판정 불능 → 0으로 읽어 «모든 클래스 위»로 보고, 거짓 집합 속 평행 평면이
    «선 공유»로 오인돼 벽 가족이 합쳐지고, 11평면 기록이 모든 이름을 벽 끝 모서리로 접었다 —
    `DegenerateWitness`가 막았을 뿐이다. ⇒ **규칙 하나** `combinatorics::canonical_triple(jd, s)`
    (정렬된 클래스 집합의 사전식 최소 독립 삼중; 셋이면 판정 없이 그대로 — 비트 무변)을 생산자 넷이
    부른다: 피연산자 링(정점이 인접 면의 클래스로 **스스로** 이름 짓는다; `EdgeFaces`가 간선 걷기에서
    정점 → 면 표를 함께 짓는다) · 별칭 대표(피연산자 씨앗 `seed_from_operands`) · `third_on_l`(«면의
    평면이 이름에 있나»를 버리고 «`wc`가 이름에 있나» 하나) · 결과 정점 정의(입사 결과 면의 평면 집합 —
    별칭 대표와 **일부러** 다른 집합: 08-12의 «결과에 없는 면으로 정의» 사고). **원장은 타입**: `Canon3`
    뉴타입(생성자 = 집합에서 고르는 `canonical_triple`·구성이 셋을 주는 `Canon3::three`)만
    `NodeId::three_planes`·`Def::Three`가 받아 일곱 번째 철자는 컴파일 오류. 완전성: 평면 넷이면 모든
    기록이 같은 집합이라 성분 하나(다섯 이상은 코퍼스에 없고 중단 조건). ☑ 실측: census 가족
    `fourplane` 여덟 쌍 전부 Ok, **두 fold 순서의 결과가 비트 동일**, 옛 329행 무변; 감사
    (`operand_vertex_audit`)가 정점마다 이름 하나·종속 0·접힘 전부 정점 위, 모션 군 전부에서도(강체
    운동 **오라클** 행은 못 넣었다 — 기울어진 평면이 낀 삼중의 f64 실현은 두 경로가 마지막 비트에서
    갈린다, 실측); 사용자의 4부품 fold가 스크립트 순서로 한 몸. `reuse::triple`의 |S| ≠ 3 기권은 남긴다
    (그 패스엔 판정기가 없고 pass-through 면의 이음이 미측정). dev-log 「칸 ⑪」.
  - ★★★★★ **원통 위의 점도 이름 하나 — 축을 지나는 평면의 룰링이 접선 룰링과 겹친다 (칸 ⑫, 2026-09-06).**
    칸 ⑪이 남긴 벽 둘(사용자 거싯 `1.5`/`−2.5`, `rrect-box`): 거싯 옆면 평면이 필렛 축을 지나 그 클래스의
    룰링 하나가 필렛의 **접선 룰링**이고, 그 끝의 접선 코너가 이름 셋(자기 `Branch … Double`, 세 평면
    `[cap, t, wc]`, 클래스의 룰링 교차 `Lo/Hi`)으로 들어와 옆면 스윕이 `theta_between`의 일치에서
    거절했다. ★ 계획의 «인식 하나로 짓는다»는 틀렸다 — 실측 사다리는 여덟 계단이고, 앞의 여섯은 동일성 → 핀 → 유령 절단점 →
    병합 층의 이름 → 마스크 → 짓는다. **규칙**: (1) **동일성은 별칭 표** — 발견은 씨앗 하나
    (`seed_from_operands`: 피연산자의 branch 코너 × 클래스, `side_of == 0` →
    `Aliases::record_on_cylinder`가 이름 셋을 유니온; 추적 중 발견은 같은 라운드의 다른 클래스가 먼저
    거절할 수 있어 버림), 대표는 **`Branch` 우선**(원통 위의 점은 원통 이름이 대표 — 차트의 θ 기하가
    이름에서 원통을 읽는다), 병합 층 셋(`merge_coincident`·`merge_circles`·`merge_rulings`)이 모든
    이름을 정준화하고 분할·스윕·차트 정거장(`Curved.aliases`)이 `canon_point`로 비교한다(소비자 5 →
    12; 지역 거절 여섯 자리는 «표가 모르는 일치»의 바닥으로 남김). (2) **핀은 (이름, 선)의 함수**
    `combinatorics::pin_for` — 이름이 접히면 다시 유도. (3) **클래스 원과의 교차는 호 위일 때만
    절단점**(`split_circles`가 θ 순서를 먼저 매기고 덮인 것만 남김) — 온전한 원의 근이 진짜 정점(슬래브
    윗면·거싯 발 `y = 0 = −2 + r`)과 겹쳐 `CoincidentNodes`였고, 1.6 대조군엔 op당 4개·`rrect-box`엔
    16개의 유령 절단점이 조용히 있었다(옛 353행 0, 다이제스트 무변). (4) **옆면은 자기 쪽을 진술한다**
    — 접선 선은 평면 쌍 선이라 평면 어휘로(`Seg{wall: t, Graze}`; 룰링 어휘로는 침묵), 바닥면의 Graze와
    `merge_coincident`에서 합쳐져 `edge_mask`가 두 비트(판의 물질이 그 선 양쪽에 있다); 쪽은 룰링 run과
    같은 식(`body_side`), `node_axis_param`은 이름의 어느 평면이든 ⊥인 것을 읽는다. ☑ 실측: census
    가족 `tangentline`(`fillet-slab`·`far`·`user 1.5`)과 `arcwalls p1p3`·`rrect-box` 전부 Ok — p1p3
    `55.5 / 48 / ∅`, rrect-box `261.7 / 197.7 / 32`, 최소 픽스처 `(52 + π) + 3`, 옛 353행 무변; 사용자의
    4부품 fold가 **스크립트 원래 치수**(`1.5`/`−2.5`)로 양 순서 비트 동일; 감사가 접선 코너마다 대표
    하나(`Branch`)·다른 근은 다른 점·거절 0, 모션 군 전부에서. 별건으로 남김: 슬래브 면이 접선 벽과
    동일평면이며 접점 너머로 뻗는 `slab wall 1.5`(`StraightAngle`)·`1.6` cut(`OpenResultShell`).
    ★ **잠금이 계단 둘을 더 찾았다** — 최소 픽스처를 피연산자 순서를 바꿔 돌리자 같은 결함 모양이
    두 번 더 나왔다: (5) **«어느 룰링인가»는 점과 벽의 함수**(`ruling_side_signed`) — 두 자리가 그
    `0`을 이름의 **근**(`Double`)에서 읽고 있었고, 그것은 이름이 그 벽을 짝지을 때만 참이다;
    (6) **룰링 간선의 담체 평면은 그 간선을 경계 짓는 면의 것**(`pair_surfs` 스캔을 룰링 키로 확장) —
    끝점 이름의 공유 평면에서 유도하면 대표가 다른 쌍일 때 «선을 담기만 하는» 평면을 적어
    `EdgeCarrierMismatch`가 잡는다. 규칙은 하나다: **값은 그것이 무엇인가에서 유도하고, 이름은 이름일
    뿐이다.** dev-log 「칸 ⑫」.
  - ★★★★★ **셀이 셀 안인가는 «증인 하나»가 답한다 — 공급 하나, 판정 하나, 루프 하나 (칸 ⑬, 2026-09-07).**
    사용자 스크립트(네 코너가 전부 필렛인 판 + 보어 넷 + 리브 둘)가 **첫 fuse**에서 죽었다:
    `WitnessNotRational`. 상대 솔리드는 무관했다 — 200 떨어진 상자와도 실패한다. 판 자신의 캡 링이
    간선 여덟(직선 넷 + 필렛 호 넷)이고 **여덟 코너가 전부 접선 코너**라, 「이 링이 저 보어 원 안인가」를
    묻는 술어가 증인을 삼평면 이름에서만 찾아 0개를 얻고 정직하게 기권한 것이다. 그 점들의 좌표는
    **여덟 중 여덟이 유리수**였고(실측), 형제 술어는 그 공급을 이미 쓰고 있었다. 즉 계산이 틀린 게 아니라
    **규칙이 여러 벌이고 그중 하나가 잘렸다**: 증인 공급 다섯 · 점 판정 넷 · 루프 셋 · 네 갈래 디스패치 **둘**
    (배열과 동일평면 병합이 각각 한 벌씩; 병합 쪽은 내부-증인의 역방향 절까지 빠져 있었다).
    ⇒ **모듈 하나** `nesting.rs`가 그 물음 전부를 갖는다. `Cell`(링 또는 디스크) · `Witness`
    (**경계 위**냐 **내부**냐 — 경계 위 증인 하나면 결정되고 내부 증인은 역방향을 물어야 한다는 규칙이
    주석이 아니라 **타입**이다; 그 역방향은 상대의 경계 증인으로 **한 단계만** 판다) · 게으른 공급 하나
    (이름 → 유리수 코너 → branch 코너 → 현 중점 → 간선 내부 → 내부 증인) · 판정 하나 · 루프 하나와
    **세 이름이 한 규칙에서**(값을 못 만듦 → `WitnessNotRational`, 빈 목록 → `RingHasNoWitness`, 전부
    기권 → `NoClearRay`). 두 도로가 같은 문을 쓰고, 옛 철자 넷(`ring_in_ring`·`node_in_circle`·
    `circle_center_in_ring`·`ring_in_ring_by_witness`)은 호출자가 없어져 **삭제**됐다. 증인 원자는 모듈
    비공개라 여섯째 공급은 컴파일러가 막는다. 새 기하 판정은 **0개**. ☑ 실측: census 가족 `roundplate`
    아홉 행이 열리고 부피가 손계산과 정확히 일치(판 `51895.2213`, 상자와 `81895.2213`, 리브와
    `76295.2213`), 옛 371행 **비트 동일**, 두 프로파일 동일, release perf 같은 세션 A/B 무변, 코퍼스의
    nesting 물음 **33,781**개 중 삼평면 이름이 없는 링 684개는 전부 상대가 링이라 잘린 팔을 **한 번도**
    밟지 않았다(그래서 거절 수로는 찾을 수 없었다). 사용자 스크립트의 fuse가 지어지고, `cut`은 **다음 벽**
    — 한 클래스가 원과 룰링을 둘 다 드는 경우 — 에 안착한다(실측: 판 자신의 옆면 `x = 45`, 코너 필렛에
    접하면서 가로 원통에 잘린다). 별건으로 얼림: 필렛 외곽이 90° 회전 여섯에서 `NoClearRay`(이 칸 이전과
    같고 보어와 무관 — 옛 트리에서 확인). 병합 도로의 역방향 절은 그 물음을 밟는 픽스처를 **만들지
    못했다**(그 도로는 간선으로 이어진 그룹 안에서만 주인을 찾는다) — 논증으로만 선 계단이라고 적었다.
    dev-log 「칸 ⑬」.

  ★★★★★ **한 클래스는 원과 룰링을 둘 다 들어도 된다 — 칸 ⑭ (2026-09-07).**
  `ClassEdges::of`가 그 조합을 통째로 `RulingBoundNotYet`으로 거절하던 것을 **지웠다**. 거절이
  막고 있던 위험은 «둘 다 든다»가 아니라 «둘이 **만난다**»다 — 만나면 어느 분할도 그 교점을
  안 만들고(`split_circles`는 원을 선분하고만, `split_rulings`는 룰링을 선분하고만 자른다) 점 자체가
  `평면 ∩ 원통 ∩ 원통`이라 이름이 없다. 그런데 그 인구는 **두 겹으로 이미 비어 있다**: 원은 온전할
  때만 간선이 되고(아니면 `PartialCircleUncut`) 룰링은 면의 선분이므로, 둘이 만나면 옆면 둘이 점을
  공유하고 그것은 `lateral_faces_clear`가 부정하는 바로 그 명제라 `CylinderPairContact`가 배열 전에
  거절한다(다른 솔리드 쌍); 같은 솔리드 쌍은 쌍 루프가 건너뛰지만, 축이 나란하지 않은 두 원통면을 한
  몸에 담으려면 불리언을 거쳐야 하고 그 순간 다시 이 게이트를 만난다. **그래서 검사를 출하하지
  않는다** — 소비자가 닿을 수 없는 장치가 되기 때문이다. 대신 유도를 주석으로 적고
  `#[cfg(debug_assertions)]` 술어(`circle_meets_ruling_rat`, 부호 케이스 하나와 제곱 둘로 유리수에서
  닫힌다) 위의 `debug_assert`로 지키며, **의무 한 줄을 쌍 게이트의 거절 자리에 적었다**(그 거절을 여는
  칸이 이 검사를 출하 코드로 진다). ☑ 실측: 코퍼스 전체 디버그 스위트에서 교차 **0건**(사용자 모델의
  쌍 8개는 «못 잼» 0), 옛 census 389행 **비트 동일**, 두 프로파일 동일. 사용자 스크립트의 `cut`은
  **다음 벽**으로 옮겨 앉았다 — `boolean::tangency_reject`의 `line_in_another_plane`(계획이 «게이트의 호
  검사»라 재구성했던 것은 실측에 반증됐다; 그 호 거절은 `tangency_rows`가 의도적으로 삼키는 자리였다).
  ★ 곁가지 둘: 원 중심 유도 `circle_centre_rat`이 **소비자 둘**(nesting의 증인 공급, 이 칸의 그물)을
  갖게 돼 `nesting` 밖 `combinatorics`로 **옮겼다** — 사본을 만들지도, 칸 ⑬이 세운 «증인 원자는
  모듈 비공개» 를 열지도 않는 유일한 길이다. 그리고 「원이 띠 **안**에 통째로 드는」 배치(걷기가
  처음 볼 모양)는 감사가 **0**을 세고 손으로 지어도 피연산자 게이트가 먼저 기권한다 —
  계기의 실패가 아니라 아직 안 들어오는 모양이다.
  dev-log 「칸 ⑭」.

  ★★★★★ **게이트의 면 판독기가 «조각의 도달 범위»를 읽는다 — 칸 ⑮ (2026-09-07). 사용자 스크립트가
  끝까지 지어진다.** 발자국 규칙은 면을 **꼭짓점**으로 근사하고 «면 ⊆ 꼭짓점들의 볼록 껍질»에
  기댔다. 원판(외곽 루프가 원 간선 하나인 면 — 원통의 캡)은 꼭짓점이 없어 그 전제가 비고, 도로는
  «못 읽는다»로 거절했다. 그 장벽은 건전했지만 **비껴가는 원판까지 함께 막았고**, `tangency_rows`가
  그 실패를 «안 비껴갔다»로 접어 **있지도 않은 접촉 행**을 썼다 — 그 행이 `line_in_another_plane`을
  만나 사용자의 `cut`을 멈추고 있었다. 고침은 규칙을 «경계 **조각**의 도달 범위»로 읽는 것이다:
  꼭짓점은 폭 0인 조각, 원판은 `중심 ± ρ`. `Corner::Disk` 한 철자, 판독기 **한 함수**
  (`corner_of` — 두 소비자가 각자 짓고 있었고 둘째가 축약본이라 접선 코너를 못 읽었다), 그리고
  **폭을 가진 조각은 점이 하나로 답하던 물음을 둘로 답한다**: 띠는 네 답(`StripSide::Crosses` 추가 —
  한 조각이 혼자 양쪽에 걸친다), 축은 `axis_side` 대신 «`t` 너머 그 쪽으로 닿는가»(하나의 `Orient`로는
  걸친 원판이 위·아래 둘 다라고 말할 수 없다). 산술은 **하나만 새로**: `cylinder_strip_side`의 유리수
  몸통이 여유 `ρ`를 받고 오늘의 것은 그 `ρ = 0` 문이 됐다(부호 케이스 하나와 제곱 둘로 유리수에서
  닫힌다). 칸 ⑭의 디버그 그물은 그 문의 «한 경계» 물음이라 **삭제**됐다 — 몸통 수가 하나 줄었다.
  ☑ 실측: 코퍼스의 접선 행 21 = 판정 9 + **기권 12**였고 기권 12는 **전부** 못 읽은 원판이었다;
  뒤에는 **기권 0**, 못 읽은 면 **0**, 판정 15. census 398행 중 움직인 것은 사용자 스크립트 셋뿐
  (거절 → 지어짐), 부피 `93783.7175`가 손계산과 일치하고 모션 뒤에도 같다. 곁가지: 그 장벽을 재던
  유일한 픽스처(`a_disk_face_on_a_wall_plane_…`)의 원판은 실은 띠를 **비껴가고** 있었다 — 이제
  지어지고 부피가 `8000 − 210π`로 맞는다. 그리고 «못 읽은 면은 무죄가 아니라 모름»을 함께 고쳤다
  (`tangency_rows`의 doc이 그렇게 적어 두고 코드는 안 그랬다 — 읽지 않은 `straddles`가 판정에
  들어가 조용히 무죄가 됐다). dev-log 「칸 ⑮」.

  ★★★★★ **링의 닫힘은 암묵적이다 — 소비자가 그것을 알아야 한다 — 칸 ㉕ (2026-09-09).**
  사용자가 앱에서 **원통 테두리가 끊어진 것**을 보고했고 진단까지 맞았다. 커널의 관례는
  「링은 첫 점을 끝에 반복하지 않는다」이고(`sample_edge`의 온전한 테두리 팔이 `0 … (n−1)τ/n`에서
  멈추고 `boundary_ring`이 wrap-around 중복을 버린다 — 메시 정점이라 반복하면 퇴화 삼각형),
  앱이 `windows(2)`로 **열린 사슬**을 그려 구간 하나가 빠졌다.
  - **고침은 소비자 쪽 한 줄**(`edge.vertices[0] == edge.vertices[1]`이면 닫는 구간 추가).
    커널은 안 고쳤다 — 그 관례를 커널 소비자 셋이 공유한다. **진짜 구멍은 계약이 `Tessellation`의
    말로 없었다는 것**이고 거기 적었다.
  - ☑ **세고 나서 고쳤다**: 원통 선분 358 → **360**, 차수 `[[1,4],[2,356]]` → `[[2,360]]`,
    상자는 무변(`[[3,8]]`). ⚠ 인구는 「평범한 원통」이 아니라 **잘리지 않은 모든 테두리** —
    보어는 위·아래에 원형 입구 둘을 남기므로 **구멍 있는 거의 모든 부품**이다(예측 정정).
  - ⚠ **심은 위반이 두 겹으로 눈이 멀어 있었다**: 「차수 1 없음」은 **남는 구간**에 구조적으로
    무력하고, 픽스처에 **호 간선이 없어** `closed` 검사 자체가 미검증이었다 ⇒ 히스토그램 고정 +
    **필렛 프리즘** 픽스처. 커널 무변(census 398행 비트 동일). dev-log 「칸 ㉕」.

  ★★★★★ **「간선이 이름 짓는 점」은 한 규칙이고, 성분에게는 그것이 없었다 — 칸 ㉔ (2026-09-09).**
  `diamond-void`(네 꼭짓점이 벽에 닿는 마름모 공동)가 `NoClearRay`로 거절됐다 — 평면 엔진에서
  **표면화되는 유일한** `NotSupported` 벽. 성분의 3-D 깊이는 자기 점 하나를 다른 성분에 물어 정하는데
  `probes_of`가 «이름이 있으면 **꼭짓점만**» 내놓고, 두 폴백(`corner_probes`·`coord_probes`)은
  원통 전용이라 **평면 성분에는 폴백이 없다**(실측 `named=24 corner=0 coord=0`). 마름모의 꼭짓점
  여덟이 전부 벽 위 ⇒ 공급 고갈.
  - ★ **렌즈가 첫 안을 물렸다**: 「간선 중점 공급자를 새로 짓는다」는 그 규칙의 **다섯째 철자**였다.
    저장소가 그 수를 이미 적어 뒀다(`Row` doc: *"the five supplies **today's four spellings** draw
    from"*). 그리고 `chord_midpoint_rat`은 이름과 달리 **호를 거부**하므로 그것도 «직선 간선의 내부
    점»이다 ⇒ 두 생산자는 처음부터 **한 규칙의 두 팔**이었고, **빠진 팔은 가장 평범한 것**
    (양 끝이 유리수인 직선 간선)이었다.
  - **고침**: `combinatorics::edge_interior_points` 하나 — 팔 셋(켤레 두 근 · 한 직선 위의 두
    branch · **양 끝 유리수**), **이터레이터**(팔 a·b가 한 간선에 동시에 맞아 오늘 점 «둘»을 내므로
    접으면 증인이 사라진다; ☑ 실측 **284**건(실행 간 안정), a 단독은 **0**건, c는 **~10⁴**건으로 최대
    공급 — 뒤 수치는 proptest가 닿아 실행마다 흔들리므로 자릿수로만 적는다).
    소비자 셋: `nesting`의 두 철자(코너까지 포함해 `edge_witness_points` 한 줄로 — 첫 고침은 «내부»
    절반만 접어 코너 절반이 두 벌로 남았고 감사에서 마무리했다), 그리고 **새로** 3-D 성분 공급
    (`first_deciding`의 지연 2단계 — 꼭짓점이 답하면 비용 0; 그쪽은 코너를 이미 `Probe::Named`로
    가지므로 «내부»만 읽는다).
  - ☑ **census 398행 두 프로파일 비트 동일** · workspace 1,283/0 · ignored 133/0 · kit 초록 ·
    perf 잡음.
  - **`diamond-void`가 진실한 이름을 얻는다**: `SelfTouchingResult` at 코너 간선. ⚠ **어느** 코너인지는
    기하의 명제가 아니다 — 마름모는 네 벽에 다 닿고, 형제 ⑤b의 삼각형은 코너가 하나뿐이라 그것이
    나온다. all-graze 형상 넷으로 잠갔고 **가설이 확인됐다**: all-graze `Cut`은 본질적으로 퇴화라
    전부 거절이며, 같은 입력의 `Fuse`·`Common`은 **정확한 부피**로 지어진다.
  - ★★ **예측 못 한 이득**: `bores_change_nothing_about_a_rounded_plate_under_rigid_motion`이
    「필렛 외곽이 모션 여섯에서 아직 `NoClearRay`」를 **인구 수로** 얼려 뒀는데, 그 수가 **6 → 0**이다.
    그 테스트 자신의 `Ok` 팔이 한 몸·`validate` 깨끗·**정확한 부피**를 검사하므로 여섯이 옳게 지어짐이
    이미 잠겨 있다. ⇒ `NoClearRay`는 이제 전-스위트에서 **0~1회**만 울리는 인구 없는 백스톱이다.
  - ⚠ **빗나간 예측 둘**: 「폴백을 `named` 앞에 두면 census가 움직인다」 → **비트 동일**이라 그 심은
    위반은 **잡히지 않는다**(뒤에 두는 이유는 비용과 증명이지 관측이 아니다) · 증인 선분을 특정
    코너로 못 박은 것도 과잉이었다. dev-log 「칸 ㉔」.

  ★★★★★ **정준화는 «넷»에 대한 것이고, 법선만 읽는 자는 그것을 물려받지 못한다 — 칸 ㉓ (2026-09-08).**
  「구멍 뚫린 정육면체」가 `WitnessNotRational`로 거절됐다 — 좌표가 **계산된(비종료) 십진수이고**
  대략 1e-4 이하일 때. `1/3`·`1/3·1e-3`·`π·1e-2`·`√2·1e-3`은 지어지고 **`1/3·1e-4`·`1/7·1e-3`**은
  안 됐다. 기하는 그 선을 넘나들며 아무것도 안 변한다 — 변한 것은 **평면 오프셋의 분모**다.
  - **뿌리**: `PlaneName`의 정준화는 *분모 털기 → **내용(gcd) 나누기** → 부호 규약*이고 그 내용은
    **네 계수** 전체의 것이다. ⇒ 법선(셋)만 읽는 소비자는 `gcd(a,b,c)/gcd(a,b,c,d)`를 떠안는다.
    축 정렬 평면 `z = s`이면 계수가 `(0,0,D,−N)`이고 `gcd(D,N)=1`이라 4-벡터는 정준이지만
    **법선은 `(0,0,D)`**다. 그리고 `Chart2dRat::of_normal`의 `e2 = n × e1`이 그것을 **제곱한다**
    (실측 `n=(0,0,5·10²²)` → `2.5·10⁴⁵`).
  - **고침**: `primitive_normal` — 법선을 자기 세 성분의 gcd로 나눈다. **양수 스칼라배**라 두 축의
    방향이 그대로고, 그래서 이 칸은 census가 **비트 동일**함을 *측정*이 아니라 **증명**으로 안다.
  - ★★ **더 싼 차트 셋을 검토하고 전부 물렸다**(그리고 그 이유를 `Chart2dRat`의 doc에 적었다):
    **좌표 드롭**은 census를 398 → **389행**으로 후퇴시킨다(패리티는 아핀 불변이지만 **퇴화
    패턴은 아니다**) · **`e2 = ê_k`**는 평면 위의 좌표로는 완벽한데(`p·e2_old = |n|²p_k + n_k d`,
    양의 아핀) `ring_interior_candidates`가 축을 **3-D 평면 방향**으로 써서 캡 증인을 민들므로
    자기 평면 밖의 `Witness::On`이 된다 — 칸 ㉒가 방금 닫은 그 병 · **`e2 = ê_j × n`**은 1차이고
    평면 위지만 `e1`과 **직교하지 않아** 광선의 등위선을 돌린다. 오늘의 틀은 *"not orthonormal"*이
    아니라 **평면 안에서 직교**이고, `axes()`의 문장이 그것에 기댄다.
  - ☑ **실측**(lib + census + reject_census 전수): census **398행 두 프로파일 비트 동일** · 형상 4 × 크기 3 = **12건이 한 몸·`validate`
    0·부피 상대오차 ~2e-16** · 코퍼스의 유일한 폭 거절(5회) 소멸 · reject-trace **이름 집합 17개
    무변** · perf 잡음(±1.5%) · 심은 위반 넷 전부 빨강.
  - ⚠ **벽은 닫힌 게 아니라 옮겨졌다**: 남는 인구는 «**원시** 법선 > 2⁶³»이고 여유가 약 6비트다.
    총함수 답은 이 도로의 끝 질문(패리티=부호)을 `BigInt`로 올리는 것 — 게이트가 이미 그렇게
    총함수가 됐다. 이 칸은 그 사다리의 한 계단이다.
  - ⚠ **계기 정정**: reject-trace 히스토그램의 **수**는 실행마다 흔들린다(같은 코드 두 실행에서
    `witness_not_rational` **1115 vs 165**) — proptest가 무작위 케이스를 돈다. **이름 집합만**
    계기다. dev-log 「칸 ㉓」.

  ★★★★★ **원 셀은 클래스가 축에 수직일 때만 원이다 — 칸 ㉒ (2026-09-08).**
  칸 ㉑이 림 증인 넷을 주면서 **아무도 강제하지 않는 전제**에 소비자를 하나 더 얹었다. 깨지면
  `Witness::On`이 거짓이 되고 `inside`가 **역방향 없이** 답한다 — 거절이 아니라 **조용한 오답**.
  - **전제를 지키는 것은 게이트 하나였다**: `cylinder_gate`가 (평면 클래스 × 실린더) 전부를 훑어
    비스듬한 쌍을 `ObliqueCylinderCut`(타원, M6-3)으로 거절한다 ⇒ **클래스가 원을 담는다 = 평면이
    옆면과 만난다** ⇒ 비스듬했으면 이미 거절됐다. 앉은 원 생산자엔 더 강한 근거가 있다 —
    `LoopRing::Circle`은 **닫힌 간선 하나**에서만 나므로 경계가 온전한 원이고, 원은 자기 평면을
    정하며 원통의 **원형** 단면은 축에 ⊥다(게이트가 바뀌어도 산다).
  - **소비자는 둘이고 각각 다르게 틀린다**: 림은 `û ⊥ axis`가 `û ⊥ n`을 함의하지 않아 **평면을
    벗어나고**, `disk_in_disk`는 **반지름**이 타원의 폭이 아니라 틀린다(중심 둘은 평면 위라 거리는
    옳다). ★ `ask`의 `radial`은 **기울어져도 옳다** — 축까지 거리를 묻는데 타원 셀이 정확히
    「평면 위 + 축에서 r 이내」다.
  - **고침**: 규칙에 이름을 준다(`class_carries_circle`, 총함수) · 림은 **죽이고 이름은 유지** ·
    `disk_in_disk`은 **새 이름 `ObliqueCircleClass`**(`SuspectedDefect`)로 거절한다.
    ⚠ `ObliqueCylinderCut`을 재사용하지 않았다 — 그것은 `NotSupported`(*"커버리지가 끝난다"*)인데
    이 자리의 문장은 *"닿을 수 없어야 할 백스톱이 말했다"*이다. 선례는 `PartialCircleUncut`.
  - ★ **`disk_in_disk`은 조용히 틀리는 절반이 실재한다**(대수로): `n=(1,0,0)`·`m=(3,4,0)`이면 단면이
    ŷ로 늘어나고, `ra=1`·`rb=2`·**장축으로 `4/3`** 오프셋이면 진실은 inside인데 오늘 `Ok(false)`다.
    ⇒ 가드는 그 절반에서 **정확성을 사고**, 나머지에서 건전했던 답을 판다.
  - ☑ census **398행 무변**(두 프로파일) · 표준 잠금·`reject_census`·`--no-default-features` 초록 ·
    ignored 133/0 · kit 초록 · 전-스위트 reject-trace에서 새 이름 **1회 — 그리고 그 하나가 이 칸의
    잠금 자신이다**(유일한 raise 자리를 그 잠금만 밟는다) ⇒ **코퍼스·프로덕션 경로는 0**이고,
    도화선은 살아 있다.
    심은 위반 둘이 빨강이고, 림 쪽은 debug에서 **사후조건 단언이 먼저** 운다. dev-log 「칸 ㉒」.

  ★★★★★ **원은 자기 림에서 증인을 댄다 — 칸 ㉑ (2026-09-08).**
  사용자의 브래킷(XY와 ZX에 각각 필렛 외곽 + 보어 둘, `fuse`)이 `no_clear_ray`로 막혔다.
  **실패 프레임을 계측으로 찍었다**: `a=disk b=ring may_ask_interior=true asked=true` — 출처 셀이
  `Cell::Disk`이고 그 공급은 **중심 하나**뿐이다.
  - **기제**: 보어의 중심이 **필렛 호의 중심과 같은 점**이고(코너 `(−55,−45)` r 12의 호가
    `(−43,−33)`에 중심), 클래스 차트의 패리티 광선은 `p_x`를 따라 간다 ⇒ 그 유일한 증인의 광선이
    **필렛의 접점을 정통으로** 지나 기권한다. p2 클래스(직교)에서도 같은 모양이다.
  - **뿌리는 참인 전제의 거짓 귀결**: 「원에는 노드가 없다」는 **배열의 노드**에 대해 참이지만,
    거기서 「이름 지을 경계 점이 없다」가 따라 나오지 않는다. 원의 경계 점은 **기하적**이고 정확히
    유리수다 — 그리고 커널은 이미 그 점을 민든다(온전한 원의 `ref_dir`가 `(정점 − 중심)/radius`라
    `centre + radius·ref_dir`가 **그 링의 정점 그 자체**다). 배열이 버리는 것은 **노드**지 점이 아니다.
  - **고침**: `nesting::rim_and_centre` 한 벌이 `Cell::Disk`에 **림 넷(`On`) → 중심(`In`)**을 지연
    공급한다. 프레임은 `nacre_scalar::cyl_unit_frame` — **투영하지 않고**, `ref_dir ⊥ dir`이 아니면
    **기권한다**(모든 소비자의 건전성이 `û·dir = 0`을 «사실»로 요구하기 때문이다).
  - ★★ **넷이고, 어느 쌍도 못 뺀다**: 림 점이 중심의 광선 선 위인 것은 `û ∥ e₁`일 때뿐인데
    `û₁ ⊥ û₂`가 그 평면을 펼치므로 **최대 한 쌍만** 그 선에 놓인다 ⇒ 항상 결정하는 쌍이 하나 있다.
  - ★★ **닫혀 있던 길도 열린다**: `inside`의 역방향 패스는 내부 증인을 **플래그도 없이** 건너뛰므로
    디스크는 **언제나** `RingHasNoWitness`였다(코퍼스 인구 0이라 아무 거절도 그것을 보인 적이 없다).
  - ☑ **실측**: 브래킷 **한 몸 · 부피 `118488`**(π가 상쇄하는 정확한 유리수 오라클) · `validate`
    깨끗 · **면 16**(곡면 8 + 평면 8) · bounds `[−55,−45,0]…[55,0,62]`. census **398행 무변**
    (두 프로파일) · 표준 잠금·`reject_census` 무변 · ignored 133/0 · kit 91/0.
    심은 위반 둘이 정확한 이름으로 빨강(림 제거 → `no_clear_ray`와 `RingHasNoWitness`; 반경 `R/2`
    → 민드는 자리의 사후조건). dev-log 「칸 ㉑」.

  ★★★★★ **접선은 가로지름이 아니다 — 칸 ⑳ (2026-09-08).**
  사용자의 부품(블록 + 교차하는 보스 둘, 슬랩과 보어 둘로 절단)이 `circle_meets_ruling`으로 막혔다.
  슬랩의 벽이 `x = ±12`에 서고, 거기서 `d 30` 원통의 룰링이 `y = ±√(15²−12²) = ±9` — 그 벽에 놓인
  `d 18` 보어 **원의 반지름과 정확히 같다**. 곡선이 곡선에 **접한다**.
  - **그물이 «닿음»을 물었다**: `cylinder_ruling_reached`는 설계상 **닫힘**(접촉도 닿음)인데,
    배열이 못 만드는 것은 **가로지름**이다. 접선은 아무것도 나누지 않는다 — 칸 ⑤가 정점의 접선에
    대해 쓴 바로 그 문장이다. ⇒ 문이 **경계를 이름으로 들고**(`touch_counts`), 그물은 **열린**
    물음을 한다.
  - ★★ **이건 «되던 것이 깨진» 게 아니다**: 칸 ⑰ 트리에서 같은 스크립트를 돌리면 칸 ⑭의
    `debug_assert`에 **그대로 패닉**한다. 앱은 **release** 빌드라 그 단언이 컴파일되지 않아 조용히
    지나갔을 뿐이고, 칸 ⑱이 그 조용함을 프로덕션 거절로 바꾸면서 드러났다. 고친 뒤에는
    **디버그에서도** 지어진다(한 몸, `validate` 깨끗, 부피 `88813.84639076446`).
  - **어휘 하나**: 「조각의 strip 도달 범위」가 `nacre_scalar::StripReach` 한 타입이 되어 strip 문과
    룰링 문이 같은 단어를 쓴다.
  - ★ **이름도 고쳤다**: `CircleMeetsRuling` → **`CircleCrossesRuling`**(`circle_crosses_ruling`).
    이름이 「만난다」인데 묻는 것이 「가로지른다」인 표류가 이 버그의 기제였으므로, 이름이 명제를
    말하게 한다.
  - ☑ census **무변**(두 프로파일) · 칸 ⑲의 잠금 전부 무변 · 심은 위반(접촉을 다시 세기)에서
    사용자 부품이 **빨강**. dev-log 「칸 ⑳」.

  ★★★★★ **규칙은 «간선»에 대한 것이다 — 칸 ⑲ (2026-09-08).**
  칸 ⑱이 넣은 `CircleMeetsRuling`이 **틀린 것을 묻고 있었다**: 클래스가 받는 것은 필렛의 **1/4 호**
  인데 검사는 **온전한 원**에게 물었다. 원(중심 `(5,5)`, r 5)은 `y ∈ [0,10]`이라 드릴의 룰링
  `y = ±√4.25`와 교차하지만, 실제 호는 `y ∈ [5,10]`이라 **닿지도 않는다** ⇒ 간선이 아닌 «원의
  연장선»의 유령 교차로 거절하고 있었다. 검사를 끄면 그 스크립트가 **정확한 부피로** 지어진다
  (실측 `1811.5896781896167`, 손계산 `…172`).
  - **규칙은 하나이고 «간선»에 대한 것이다**: *클래스의 간선은 어떤 면의 경계다 ⇒ 서로 다른 원통의
    두 간선이 만나면 두 옆면이 한 점을 공유하고, 쌍 규칙이 그것을 부정한다.* 이 한 문장이 **아무
    split도 자르지 않는 쌍 셋 전부**(원×원·원×룰링·룰링×룰링)를 덮는다. 칸 ⑭이 옳게 적었고,
    칸 ⑱이 그것을 «간선»이 아니라 «원»에 읽어 틀렸다.
  - **고침은 물음을 옮기는 것**: 기여가 진술하는 각도 범위(=면이 덮는 곳, `RimArc`와 같은 CCW
    관례)를 `node_coords_rat`/`branch_coords_rat`로 실현해 **호로** 묻는다. 실현 못 하거나 기여가
    없으면 **온전한 원으로 후퇴**(상계이므로 건전).
  - **문도 넓혔다**: `cylinder_ruling_reached`가 «구간의 끝 둘»을 받는다(점·원판은 대칭 특수형).
    `W²`를 다시 유도하지 않고 `strip_scale` 위에 얹었고, 판정은 칸 ⑱의 `sqrt_exceeds_root_sum`
    한 줄에 인자를 자리바꿈해 넣는 것뿐이다(그 항등식은 이제 **경계를 인자로** 든다).
    「호 → 방향 d를 따른 두 끝」 변환도 `arc_ends_along` 한 자리로 빠져 strip 문과 룰링 문이 같이 쓴다.
  - ★★ **그러자 검사는 backstop이 된다**: 룰링은 언제나 그 원통 축에서 `r`보다 가까이(`√(r²−h²)`)
    있으므로 호가 룰링에 닿으면 쌍 규칙이 **먼저** 거절한다 — 필렛을 6·7·8로 키워 실측했고 전부
    `CylinderPairContact`였다. 그래도 **출하는 유지**한다(민들어지지 않은 노드는 패닉이 아니라
    조용히 틀린 솔리드이고, 이 논증은 이미 한 번 잘못 읽혔다).
  - ☑ **실측**: 사용자 스크립트 Ok(두 몸, 부피 일치) · 필렛 1·2·6·7·8 **무변** · census **무변**
    (두 프로파일) · 원판 1215 케이스에서 넓힌 문이 옛 문과 **완전 일치**. dev-log 「칸 ⑲」.

  ★★★★★ **호도 조각이다 — 칸 ⑱ (2026-09-08).**
  판보다 굵은 드릴이 하는 일은 **캡 평면을 자기 반지름 안에 들이는 것**뿐이고, 그러면 게이트가
  평면 대신 **면**을 묻는다. 그때 `corner_of`가 원 조각을 **루프 전체가 원 하나(디스크)일 때만**
  읽어서, 직선+호가 섞인 루프 — 필렛 outline, 가장 흔한 것 — 이 **호의 위치와 무관하게**
  안 읽혔다(먼 코너 필렛 하나로도 막힘, 실측). 사각 판은 같은 드릴로 잘 뚫린다.
  - **호는 「끝 둘이 이름 짓는 만큼 잘린 같은 원」이다.** `Corner::Disk`와 `Arc`가 한 변종
    (`Round { centre, rho, axis, arc: Option }`)이 되고, 세 물음의 답이 전부 **`arc_extent`**에서
    나온다 — strip은 `n × m` 방향으로, 축 물음은 `m` 방향으로. 그 규칙은 칸 ⑩이 원통 옆면에
    쓴 것과 **같은 한 벌**이고, 이제 평면 면도 그것을 부른다.
  - ★ **끝의 순서는 간선의 저장 순서**(담체 축에 대해 CCW)에서 온다. 면의 법선을 쓰면 **여집합**을
    이름 짓게 되고 그건 상계가 아니라 다른 집합이다 — 단위 테스트가 그 차이를 얼린다.
  - **문도 넓어졌다**: strip 문이 「대칭 margin」 대신 **구간의 끝 둘**을 받는다(점=폭 0,
    원판=대칭, 호=비대칭 + 유리수 shift). 그 과정에서 `√A > √B + √C`가 스칼라 한 줄
    (`sqrt_exceeds_root_sum`)로 접히고 칸 ⑰의 `exceeds_root_sum`도 그것을 부른다.
    ★ 통화는 `MeetPoint` + `Rat` 그대로라 **Wide 지원이 유지**된다.
  - ★★ **그리고 칸 ⑭의 의무가 만기됐다.** 판독기를 열자 「한 클래스가 원과 룰링을 들고 그 둘이
    만난다」가 도달 가능해진다(실측: 필렛 5). 칸 ⑭의 「도달 불가」 논증은 *"원은 통째로 간선이
    되고 면은 그 일부만 쓴다"*는 구멍이 있었고, 칸 ⑰이 바로 그 쌍을 통과시켰다. 그래서
    `debug_assert`가 **`RejectReason::CircleMeetsRuling`**으로 졸업했다(`None`도 같은 이름 —
    `unmeasured`는 세 코퍼스 전부 0). 그 점은 `평면 ∩ 원통 ∩ 원통`, 이름 없는 중첩 근호다.
  - ☑ **실측**: 필렛 1·2 + 굵은 드릴이 **두 몸으로** 지어지고 부피가 손계산과 일치. census는
    S0이 예측한 **3행만** 이동(`curved_operand_boundary` → `cylinder_pair_contact`, 그 픽스처가
    원래 받으려던 답). 접선 행(`tangency_rows`) 20행 **무변**. dev-log 「칸 ⑱」.

  ★★★★★ **갈라지는 방향은 축 둘만이 아니다 — 칸 ⑰ (2026-09-08).**
  필렛 판을 가로질러 드릴을 뚫는 사용자 스크립트가 `cylinder_pair_contact`로 막혔다. 필렛 축과
  드릴 축의 거리(5)가 반지름 합(5 + 7/2)보다 **작아** 무한 곡면 rung이 못 자르고, 면 판독기는
  **두 축**으로만 물어 둘 다 실패했다 — 그런데 필렛의 **1/4 호**는 `y ∈ [5, 10]`에만 있고 드릴은
  `[−7/2, 7/2]`뿐이라 **공통수선으로는 깨끗이 갈라진다**.
  - **후보 집합을 닫았다**(`planes::separating_dirs`). 발자국이 말하는 상자의 경계는 캡 평면(법선 =
    축)·룰링(방향 = 축)·호뿐이므로 유리수 분리축 후보는 **정확히 셋**: `m_a`, `m_b`, `m_a × m_b`.
    남는 것은 호의 **반경 연속체**이고, 그건 공통 단면 차트가 있는 **평행** 쌍에서만
    (`cross_sections_clear`) 잡힌다 — 그래서 거절은 여전히 참인 문장(「갈라짐을 못 보였다」)이다.
  - **셋째 방향은 위 rung의 면 버전이다.** 전체 호·span 없음이면 `separated`가
    `|W·C| > (r_a+r_b)|C|` — `cylinders_clear`의 skew 팔과 글자 그대로 같다(테스트가 얼린다).
  - **규칙이 하나가 됐다.** 옛 `slab_clears`는 「한쪽의 **span** vs 다른 쪽의 **reach**」를 양방향으로
    묻는 두 벌이었는데, 자기 축을 따른 reach가 **곧 그 span**이라 방향 목록 하나로 접혔다.
    비교도 한 문이다: 두 reach의 양 끝이 각각 근호를 달고 있으므로 물음은 `g > √p + √q`이고,
    `nacre-scalar::exceeds_root_sum`이 근호 탑(`biquad_sign`)으로 답한다 — `reach_clears`는 한쪽이
    유리수 구간인 퇴화형이다.
  - ☑ **인구 실측**: 면 쌍 물음은 기본 스윕 **67** + census 코퍼스 **102** = **169**. 두 축이 다
    실패하는 것이 20(그중 평행 9는 단면 도로), 어긋난 11은 **공통수선으로도 안 갈라진다** ⇒
    셋째 방향이 답했을 쌍 **0**, **census 398행 무변**. 소비자는 사용자 스크립트이고, 커널 자신의
    끝단 픽스처가 «셋째 방향만 답하는» 쌍 넷을 이 인구에 처음 넣는다(기본 스윕 67 → 71).
  - ✗ **이 칸이 열지 않는 것 — 그리고 그건 드릴이 아니라 «호»다(실측).** 판보다 굵은 드릴은 판의
    **캡 평면**을 자기 반지름 안에 들이므로, 게이트가 그 평면 대신 **면**을 묻기 시작한다. 평판은
    답한다: 같은 `r = 4.5`가 **사각 판은 뚫고 두 몸으로 자른다**. 그런데 드릴에서 **가장 먼**
    코너에 필렛을 **하나만** 넣어도 같은 컷이 `curved_operand_boundary`다 — 발자국 판독기
    (`corner_of`)가 원 조각을 **루프 전체가 원 하나(디스크)일 때만** 읽기 때문이고, 호의 위치는
    아무 상관이 없다. 즉 벽은 「평면이 원통 면을 만난다」가 아니라 **「직선+호가 섞인 루프를 못
    읽는다」**이고, 칸 ⑮의 「조각은 자기 도달 범위로 답한다」가 **호 조각**에는 아직 안 왔다.
    ★ 얇은 드릴은 캡 평면이 반지름 밖이라 면을 묻지 않으므로 필렛 판도 지어진다. ★ 평판이 정확히
    반두께(`r = 4`)면 드릴이 캡에 닿아 결과가 자기를 만지므로 `SelfTouchingResult` — 다른, 참인
    거절이다. 반대쪽 경계도 옳게 선다: 반두께 6에 `r = 5.9`는 필렛(`y ≥ 5`)까지 닿으므로
    **쌍 규칙이** 여전히 거절한다. dev-log 「칸 ⑰」.

  ★★★★★ **구간 돌출은 «스케치가 여기서 시작한다»다 — 칸 ⑯ (2026-09-07, 키트).**
  `extrude(s, [lo, hi])`가 **두 번 쓸어 fuse**하고 있었고, 그래서 한 필렛 원통면이 **두 번 진술**됐다
  — `o`·`r`·`ref_dir`이 같고 축의 **부호**만 반대. 커널은 한 곡면의 두 진술에 이름이 없어 그 쌍을
  거절하므로, 모서리를 둥글린 프로파일은 **자기 평면을 가로질러 돌출할 수 없었다**. 고침은
  `hi − lo`만큼 **한 번** 쓸되 프레임의 원점을 **구간의 시작**에 두는 것 — fuse도, 뒤집힌 프레임도,
  미러 분류도 없다. ★ 「쓸고 **옮기기**」가 아닌 이유: 옮기면 결과에 **모션**이 남아(원통은 언제나
  얻는다) 같은 부피·같은 모양이 **다른 성격**이 된다. 계약도 넓어져 `lo < hi`만 요구한다
  (`[5, 15]`·`[-15, -5]` 가능 — 스크립트에서 `translate` 한 번이 준다). 기계는 이미 있었다:
  `flipped_world_frame`의 **원점 인자**가 그것이고, 「이 평면을·이 원점에·이 방향으로」 한 문장
  (`world_frame_at`)으로 일반화했다.
  ☑ **커널 정준화는 재고 «짓지 않기»로 결정했다**: 축 부호만 다른 쌍둥이가 코퍼스에 **0**이고, 오직
  한 스케치를 양쪽으로 쓸 때만 생기는데 **그 인구를 이 키트 고침이 없앤다**. 게다가 진술을 합쳐도
  손으로 쓴 `fuse(extrude(s,5), extrude(s,-5))`는 「한 곡면을 **두 피연산자**가 소유」로 옮겨 앉을
  뿐이라 관측상 무변이다(차트가 「원통 클래스는 한 솔리드의 곡면」을 전제로 서 있다 — 그 전제를
  바꾸는 것은 별개 능력). dev-log 「칸 ⑯」.

  ★★★★★ **능력 D의 첫 계단 — 차트의 «선분 집합»이 섰다 (2026-08-29, `e6846a9`).**
  새 모듈 `nacre-ops::cyl_chart`(**`#[cfg(test)]`** — 계기다). 짓고 **재기만** 했고 커토버는 없다.

  - ★★★★★ **직교성이 «측정»이 아니라 «게이트의 귀결»이다.** 인구 게이트가 (평면, 원통) 쌍을 다섯
    갈래로 가르고 옆면에 자국을 남기는 것은 **둘뿐**이다 — ⊥ 절단이 **원**(가로), 축을 정확히 품는
    벽이 **룰링 둘**(세로). ☑ 실측 인구: ⊥ **1128** · 축 통과 **75** · `0<d≤r` 거절 **7** ·
    멀어서 자국 없음 **1494** · 비스듬 거절 **2**. **여섯 번째 갈래는 없다.**
    ⇒ M6-2에서 옆면 배열은 **주기 띠 위의 직교 «선분» 배열**이다(격자는 아니다 — 룰링은 유한 구간).
  - ★★★★ **「새 수학이 필요 없다」가 확인됐다.** θ 순서는 `circular_order_about_seam`이고 그 시그니처가
    두 점을 **같은 원 위에서 받지 않는다**(`(MeetLine, QuadVal)` 쌍) ⇒ 서로 다른 벽의 룰링도 그대로
    비교된다. 걱정거리는 **분기 노드가 없는 룰링**(z도 θ도 이름이 없다)이었는데, ☑ **826개 조각 전부
    양 끝이 이름을 갖는다**(0건).
  - **`Curved`가 룰링 구간을 내보낸다**(`RulingExtent`) — `rulings_on_class`가 이미 정하는 값이라
    차트가 **다시 계산하지 않는다**. `cut_rims`가 이미 타던 그 길에 얹었다.
  - ★★★ **`bands_of`의 «자르기»는 일부러 안 한다.** 그 함수는 `keyed.retain(t ∈ row.span)`으로 **면**의
    구간에 잘라낸다(행이 면이니까). 클래스에는 그런 span이 없고, 배열에서 그 틈은 **아무것도 안 남기는
    셀**이다 — `CylRow`의 doc이 경고하는 「면 없는 자리에 띠를 발명」이 그렇게 제대로 풀린다.
  - ☑ **실측(전 lib 스위트)**: 차트 **263**개 · 거절 **0** · **세로선이 없는 것 198**(75%, 오늘의
    「쉬운 띠」에 대응) · θ 총계 **404 = `Curved`의 404, 손실 0** · 한 클래스가 행 **둘**인 것 4 ·
    가로선 2~6. **그리고 263개 전부에서 「클래스의 z-선이 그 클래스 모든 행의 띠 경계를 덮는다」**
    (사실이 만들어지는 자리에서 단언).

  ★★★★★ **둘째 계단 — 셀, 그리고 오늘의 면과 대조 (2026-08-29, `14d6ac9`·`7afa955`).**
  여전히 `#[cfg(test)]`이고 커토버는 없다.

  - ★★★★★ **갈림길이 «없어졌다» — 걷기가 필요 없다.** 첫 계단이 남긴 (가) 걷기 / (나) 곱 격자를
    재서 고르려 했는데, 조사에서 전제 하나가 둘 다를 지웠다: **룰링의 축-끝은 언제나 축-선이다**
    (263 차트 전부, 어긋난 것 0). 구조적 이유도 있다 — 룰링 노드는 `[통과벽, 다른 평면]`이고 그
    「다른 평면」이 옆면을 만나려면 ⊥일 수밖에 없으며, `chart_of`는 **모든** ⊥ 클래스를 모은다.
    ⇒ **배열이 축-구간마다 쪼개지고**, 한 구간의 셀은 그냥 **θ-섹터**다. face-walk도 꼭짓점 각도
    순서도 없다. (나)는 셀 **1724 vs 949**에 병합이 딱지를 요구하므로(=D2) 스스로 끝내지도 못한다.
  - **θ 순서는 `arrangement::circular_order`를 «부른다»** — 크레이트에 열었을 뿐 새로 쓰지 않았다.
    ★ **룰링당 노드 하나만** 넘긴다: 룰링은 수직선이라 두 끝의 θ가 같고, 둘 다 넘기면 그 함수가
    「한 점에 두 이름」으로 정직하게 거절한다.
  - ★★★★★ **대조의 1순위는 «세는 주장»이 아니라 «절대» 주장이다.** 「셀마다 면 하나」는 매핑이
    **일관되게** 어긋나도 참이다 — `3c1fa56`이 그 모양으로 실측당했다(전역 `axis_up` 반전을 잠금이
    못 잡음). ⇒ 단언은 **방출된 면이 자기 이름으로 지목하는 셀이 차트에 있는가**이고, 세는 것은
    보조다. ★ 그리고 **양쪽 팔을 세어** 초록 0이 「아무도 그 길을 안 갔다」가 아님을 보인다
    (띠 **284** · 패널 **65**).
  - ☑ **실측(전 lib, D1b-a 픽스처 제외)**: 셀 **949** · 룰링 없는 차트 **198** · θ **404** ·
    `bands_of` 구간 **545** 중 차트가 더 잘게 보는 것 **44** · θ 순서 거절 **0** ·
    역순으로만 맞은 패널 **0**.
  - ★★★★★ **계측이 «찾은» 것 둘.**
    ① **홀수 개 룰링이 지나는 구간이 0** ⇒ 「절단 하나면 셀 하나」 가지는 **검증된 게 아니라
    미실행**이다. ② **방출된 띠 80개(263 차트 중 40개)가 θ-섹터를 둘 이상 덮는다.** `band_faces`는
    **rim 둘이 잘렸는지**만 보고 갈라져 안 잘렸으면 **원 전체** 띠를 내는데, 그 자리 주석은
    *「게이트가 경계면을 옆면에서 떼어 놓는다」*에 기댄다 — **기록-하고-통과 팔(`d=0`)은 그렇지
    않다.** 가짜 절단인지 놓친 경계인지는 **딱지의 물음이고 D2 몫**이다.
  - ★★★★ **이 칸이 만들고 잡은 결함**: census를 `band_faces` 뒤로 옮기며 `?`로 그 패스의 거절을
    전파해 **차트 하나를 통째로 잃었다**(263 → 262). 차트는 **배열의** 성질이지, 거기서 면을 못 낸
    패스의 것이 아니다.
  - ★★★★★ **기울어진 축에서 두드렸다 (`14d6ac9`).** 위 유도는 축 방향에 의존하지 않지만 **수는 전부
    +z 인구의 것**이었다 — kit의 `Builder::cylinder`가 world XY 전용이고, 코퍼스의 기울어진 원통 셋은
    **전부 거절 픽스처**라 차트에 닿은 기울어진 원통이 **0개**였다. ⇒ **피타고라스 프레임**
    (`u=(0.6,0.8,0)`, `v=(−0.48,0.36,0.8)`, 법선 `(0.64,−0.48,0.6)`)으로 프리즘을 세우니 **여섯 면
    전부 narrow**이고 축과의 내적이 **정확히 0**(벽) 또는 평행(캡)이다. 관통 보어와 벽이 축을 품는
    반-보어 **둘 다 서고**, `validate` 초록에 부피가 `4−π/4`·`4−π/8`이다.
    ★★★ **그리고 벽은 커널이 아니라 «생성자»였다**: `Model::add_cylinder`는 그 doc이 스스로
    *「the test entry」*라 부르는 길로, f64 축을 정규화하고 캡을 계산된 float에서 되올린다 ⇒ 기울어지면
    narrow 창을 벗어나 게이트가 `CylinderGateUndecided`로 정직하게 거절한다(실측). 프로덕션 도로
    `add_cylinder_exact`에서는 **기울어진 경우가 무리수 경우가 아니다**.

  ★★★★★ **셋째 계단 — 세로선의 «답», 그리고 82건의 판정 (2026-08-29, `6366a8a`·`450f9cd`).**

  - ★★★★★ **간극의 이름**: `Label = [A_above, A_below, B_above, B_below]`는 **평면의 양쪽을 다
    든다** ⇒ ⊥ 클래스의 디스크 딱지가 차트의 **가로** 교차를 이미 통째로 답한다(그래서 `chamber`가
    둘을 읽고 일치를 요구한다). **세로**를 건너는 것은 **벽**을 건너는 것이고 그 답이 없었다.
    쌍대가 없는 이유도 코드가 말한다 — `RulingTrace`는 *「**옆면**이 ∥ 클래스에 남긴 자국」*이다.
  - **룰링이 이제 «원통 안쪽 셀»의 딱지를 나른다.** `per_class`가 `cells`·`face_of`·`labels`를 한
    자리에 들고 있어 조회 한 번이면 된다. ★ 필드는 `#[cfg(test)]` — 유일한 독자가 계기이므로 릴리스는
    **계산도 안 한다**(`Staged`의 `cells`/`labels`가 이미 그 선례다).
  - ★★★★★ **부호는 «유도»했고 «내용»이 검사한다.** `RulingCarrier::side`(안쪽은 `−side` 방향) ·
    보편 경계 규약(*「재료는 진행 방향의 왼쪽」*) · `world_rat_sense`(이름↔저장 법선) 셋이 합쳐져
    **`side·κ`** 하나가 된다 — `ruling_grazes`가 이미 「frame-free 곱」이라 부르는 그것.
    **새 술어 0.** 검사는 독립 서술(*안쪽 셀만 그 원통 솔리드의 재료를 든다*)이고 사실이 만들어지는
    자리에서 단언한다. ☑ **862 조각 전부 딱지, 838 일치·0 모순·24 눈멂**(제 판을 두른 보스는 양쪽에
    같은 재료라 내용이 못 가른다 — 그때는 유도만 말한다).
    ★★ **그리고 이 부호는 부피 오라클이 감시할 수 없다** — 필드가 계기라 프로덕션이 안 읽는다.
    감시자는 그 **내용 검사**이고, `side·κ`를 뒤집으면 **일치가 전부 모순**이 된다.
  - ★★★★★ **82건이 판정됐다 — 셋으로 갈라서.** D1b는 「가짜 절단 / 놓친 경계」 둘로 적었는데, 셋째가
    칸 ㉒의 물음이다(**그 룰링이 이 면 위에 없다**). ☑ 실측: **면에 없음 0 · 놓친 경계 0 · 무해 84**.
    ⇒ **오늘의 띠 도로가 옳다**: 가로지르는 벽이 옆면에서 아무것도 안 바꾼다.
    ☑ 곁수(**구간별 관측 450**이 분모 — 위 862 «조각»과 다른 분모다): 벽이 재료를 바꾸는 것 **122** ·
    그런 구간 **61** · 자국이 전부 graze인 것 **12**(존재 물음도 비어 있지 않다) · 구간의 걷기가
    **닫힌 것 201, 안 닫힌 것 0**.
  - ★★★★ **그 걷기가 `label_cells`의 마지막 검증이고, 일부러 «부호가 없다».** 「이 섹터가 벽의 어느
    쪽인가」를 물으면 **세 번째 부호**를 쓰게 된다. 룰링의 딱지는 벽 **양쪽**을 다 들므로 각 룰링이
    기여하는 것은 「둘이 다른가」뿐이고, 한 바퀴 돌아 제자리로 와야 한다.
    ★★ **그 검사가 못 보는 것도 쟀다**: 모든 구간의 룰링 수가 **짝수**라(`odd_k = 0`) **룰링마다의
    계통 오차는 XOR이 상쇄해 통과한다**(실측 — 초록). 절대 감각을 붙드는 것은 그 걷기가 아니라 위의
    **내용 검사**다.
  - ★★★ **계기가 두 번 틀렸고 둘 다 잠금을 지켜보다 잡았다**: ① 존재 계수기가 「자국 목록이 비었나」를
    물었는데 `MergedRuling`은 **정의상 안 빈다** ⇒ 0은 인구가 아니라 **아무것도 안 본** 것이었다
    (진짜 신호는 `SegKind::Graze`, 12건). ② 닫힘 검사의 첫 부정 대조가 **통과**했다(위 패리티).
  - ★ `Cell::walls`가 이름 대신 **`theta`의 인덱스**를 든다 — 안 그러면 셀이 자기 룰링의 딱지에 닿지
    못한다. 이름은 `Chart::ruling_name` 한 번이면 나온다(파생 사본 금지, D1b가 `side`를 지운 그 이유).

  - ★★★ **계획의 «표제 수»(커토버 예보)는 못 냈다.** 오늘의 답과 셀별로 견주려면 「섹터가 벽의 어느
    쪽인가」 부호가 필요하고, 그건 이 칸이 일부러 안 만든 **세 번째 부호**다. ☑ 부분 답은 된다 —
    **놓친 경계 0 · 안 닫힌 구간 0** ⇒ θ 방향으로 어긋날 자리가 없다. 나머지는 D2b가 그 부호를
    **읽을 때** 나온다.

  ★★ ~~**다음(D2b)**: 딱지를 **전파**로~~ — 넷째 계단(D2b-0, 2026-08-30)이 반증했다: **전파도 셋째 부호도
    필요 없다.** `circle_on_class`가 면 span 안의 모든 ⊥ 클래스에 원을 남기므로 면이 있는 셀은 양끝에 반드시
    가로 딱지가 있고, `Chart::read_cell`(cyl_chart.rs, `#[cfg(test)]`)이 그것만으로 챔버·존재를 읽어 오늘의
    방출과 **셀별로 일치**했다(전 lib 스위트 1056 셀, 어긋남 0 · 44건은 무해 · 디스크 끝 구간은 언제나 통째
    원). 커토버 결정: θ-병합은 차트, z-병합은 `unify_curved_faces` · 존재는 잘린 끝 = marks / 디스크 끝 = span
    (`cyl_rows` 유지) · 룰링 identity는 `Branch` root가 아니라 `ruling_side`. D4 삭제 목록과 실측은 dev-log
    「능력 D, 넷째 계단」. 그다음 D2b(커토버) · D3 · D4(청소).
  ★★ **다섯째 계단(D2b, 2026-08-30) — 커토버 완료.** 프로덕션은 `cyl_chart::emit_lateral`이 옆면을 방출하고,
    `bands::band_faces`는 test 빌드의 **참조 도로**(census가 두 방출을 면 단위로 zip)로 D3까지 남는다. 규칙 넷:
    `Band` iff 띠형 ∧ 양끝 다 잘리지 않음 · 띠는 `bands_of`의 경계(림 클래스 ∪ span 끝)가 아닌 선만 건넌다 · 런은
    림 노드 룰링에서 갈린다 · 런 순서는 첫 섹터. 열린 것: 반높이 보스(Fuse 32000+375π). 바뀐 이름: (0,0) 코너 보스가
    `OpenResultShell` 대신 `CylinderGateUndecided` — 판 클래스의 원판 셀이 B를 잃는 **arrangement 딱지 결함**(발견,
    다음 인구). 실측·부정 대조·A/B는 dev-log 「능력 D, 다섯째 계단」.
  ★★ **여섯째 계단(D3, 2026-08-30) — 참조 도로 삭제.** `bands::{band_faces, bands_of, chamber, panel_faces, panel_probe}`와
    census의 「참조 vs 차트」 대조(항진명제)를 지웠다(프로덕션 무변). 남은 눈: z-선 ⊇ `boundary_lines` ·
    `emitted_faces == whole_emitted + full_runs − z_merge_bandlike + partial_runs`(기록 자리) · `arcs_*`(옛 MARKS 계약).
    ★ `src2_disagree == 0`은 lib 스위트의 artifact였다 — 코퍼스 가족 12를 lib에 넣으며 명제를 「읽지 못한 present 셀 위로 면을
    내지 않는다」(`emit_unknown == 0 || refused`)로 고쳤다. `bands.rs` 4042 → 3591줄. 실측은 dev-log 「능력 D, 여섯째 계단」.
  ★★ **일곱째 계단(D4, 2026-08-30) — 청소 패스가 모든 영역을 말한다.** `Bound::Band { lo: Rim, hi: Rim }`, `Rim::{Circle, Chain}`:
    꿴 사이클을 **감김수**(Σ `seam_step` 부호 — `CutRim.nodes`의 순서 조회, 좌표 없음; `seam_is_node`는 반열림)로 lo(+1)/hi(−1)/
    구멍(0)으로 가르고, 사슬 림은 접촉(θ=0 위의 정점 — 감김과 **다른** 물음)에서 회전해 슬릿에 잇는다. 열린 것: 반높이 보스와
    거울, 옆면 **1**(두 사슬 모양 — 솔기 룰링을 품는 것·감는 호로 지나는 것). ★ `lateral_axis_span` 일반화는 검토가 뺐다
    (`circle_on_class`의 「바깥 답 = span의 온전한 원」 전제를 깨면 조용한 오답). 실측은 dev-log 「능력 D, 일곱째 계단」.
  ★★ **재연산 사다리 E — 첫 계단 E1 (2026-08-30): 결과 면의 링이 전부 이름을 얻는다.** 결과 42개에 «먼 정육면체 Cut» census를
    세우니 Ok 9 · `DegenerateFace` 14 · `CylSpan` 10 · `OuterRing` 5. 원인 넷: 2-간선 캡의 모서리(`shared_vertex`가 두 정점을 봄 →
    **모서리 = 반간선의 시작**, `OpenLoop`의 정리) · 원-경계 면의 삼각형(평면의 진실 점) · **솔기 이음은 모서리가 아니다**(OnSeam
    정점의 두 다리는 한 걸음) · 옆면 외곽 루프(**E2**: `span` 대신 루프들, ⊥ 클래스마다 원 위 스윕) · `CurvedRingWall`(**E3**).
    E1 뒤: `DegenerateFace`·`OuterRing` 0, `CylSpan` 29, Ok 9(Ok는 E2의 몫). 실측은 dev-log 「재연산 사다리 E, 첫 계단」.
  ★★ **E2 (2026-08-30): 옆면의 외곽 루프가 «경계 사이클들»이 되어 트레이서가 읽는다.** 조립의 슬릿 걷기를 거꾸로 읽어(`lateral_cycles`)
    림 원(`LoopRing::Rim`)·사슬·패널·구멍(끼워진 것 복원)을 이름 짓고, 두 옆면 도로는 한 읽기 `band_span`(림 둘 + 구멍들만; 아니면 `CylSpan` —
    차트가 아직 사슬·패널을 못 담는다)으로 `span` 대신 사이클을 읽는다(`CylFaceInfo::t_range`는 범위일 뿐). ★ `hole_on_class`의 Run 팔이
    호의 방향을 flank에서 유도하던 것(볼록 구멍에서만 맞음)을 생산자의 `arc_ccw`로 — 오늘의 구멍 50/50 동치를 재고 바꿨다. 재연산 Ok 9 → **14**
    (끼워진 구멍의 Fuse 5행). 남은 것: **E2-2** 사슬·패널(부재 바깥 답·룰링 t-스윕·차트 `End::Other`) · **E3** `CurvedDeparture`/`CurvedRingWall`.
  ★★ **E3-a (2026-08-30): 평면 스캔이 룰링 위의 교차를 branch 노드로 이름 짓는다.** 곡선 캐리어를 든 링의 교차점 `fc ∩ wc ∩ cyl`은
    `NodeId::branch(wc, fc, cyl, root)`이고, root는 룰링의 `side`를 두 근에 물어 고른다(`crossing_on_ruling` — 옆면의 구멍 재진술
    `hole_on_class`도 같은 문을 쓴다, 한 철자). 그 뒤의 벽 둘도 열렸다: 셀 중첩이 혼합 링을 읽고(`cell_in_cell`이 `point_in_mixed_ring`으로),
    차트 census의 «잘린 끝 = 자국 하나»·«면 있는 셀 = 딱지 있음»은 구멍 안을 지나는 컷에 반증되어 참 명제로 재진술. 가로지르는 census
    (`crossing_census_…`, 210셀·부품 합 부피 오라클): 벽 보스 Fuse × **벽 슬랩** Ok(1). 남은 것: 결과가 두 솔리드로 갈리면 그룹핑 도로
    (`Ring::edges`의 평면 벽 shim)가 혼합 링을 `BranchVertexUnnamed`로 거절 · `coord_key`에 원통 표(현 셀) · **E3-b** 호 교차(E2-2의 룰링 t-스윕 뒤) ·
    **E3-c** `CurvedDeparture`.
  ★★ **E2-2 (2026-08-30): 옆면이 띠가 아니어도 트레이서가 읽는다.** `band_shape` 게이트 대신 `lateral_shape`(범위·림·사이클들); ⊥ 도로의 바깥 답은
    **Option**(범위 엄격 안 Crosses · 림 station 온전한 Graze · 극단 None — run만 말함)이고 모든 사이클을 `cycle_on_class`(감김만 읽음)로 깎는다; ∥ 도로는
    **`ruling_sweep`**(림·θ를 엄격히 품는 호 = 토글, 룰링 런 = Graze + 양끝 호의 θ-flank가 다르면 토글 — `flanks_differ`를 ν로 읽은 것). 재연산 Ok 14 → **21**,
    `CylSpan` 0; 벽 Cut × 벽 슬랩 정확 부피. 계기 명제 셋(세로선 답의 σ·XOR 닫힘의 짝·shim의 표)을 재진술. 남은 것: E3-b 호 교차·E3-c 이탈·그룹핑 shim.
  ★★ **E3-b/c (2026-08-31): 호가 클래스의 선을 떠나고 가로지른다.** 걷기의 `on_meet`가 `EdgeMeet::{On, Departs(σ)}` — 떠나는 간선은 σ 쪽 가상 노드(σ = −ruling_side를
    `outward_fix`로; `CurvedDeparture` 삭제); 스캔은 호 위의 교차를 `crossing_on_arc`(θ-순서, 정확히 한 근)로 짓고 캡의 현이 overlay에 합류(sense는 canonical 근에서); 현 셀은
    원통을 클래스 표에서, 원 구멍의 현 flip 노드는 **게이트의 `crossings`에 실린 쌍에만**(룰링 도로와 한 규칙 — 없는 룰링에 노드만 심으면 스퍼, d=0 보스가 반증). 재연산 Ok 21→34 ·
    가로지르는 census Ok 69→86 · `rul flush` Cut/Common 빌드, Fuse는 결과 전체 판정의 새 가드(stated 평면-자기 간선 → `CoplanarMerge`)가 거절 — 병합 키 확장은 이름-도로 퇴행으로
    철회(그룹핑 팔 캐리어화의 입구). 실측은 dev-log 「재연산 사다리 E, E3-b/c」.
  ★★ **그룹핑 팔 캐리어화 (2026-08-31): 이름-도로가 혼합 링을 읽는다.** `Ring::edges`가 walls에서 캐리어·핀을 직접 짓고(legacy shim과 구성-시점 `BranchVertexUnnamed` 사망), 격리·좌표
    도로는 `ring_is_mixed`로 갈라 유리 탐침 + `point_in_mixed_ring`(combinatorics로 승격, 한 철자)으로, 병합은 (방향쌍, 벽) 멀티셋의 **2-패스 지움**(정확 반대쌍+`Wall::reversed` →
    평면-전용 1:1; 구성상 공허하던 ≥3 검사는 삭제)으로. 가로지르는 census BVU 18→**0**(mid 15는 라벨 광선의 패널/사슬 옆면 벽 — `NoClearRay`, ☑ miss-first 칸이 닫음) · `rul flush` **[Ok×3]** · census의
    접촉-컷 4행 **빌드**(부피 전부 해석값 정확 — `arc turned`는 개통이(구성-시점 BVU), VNAS 3행은 병합이 열었다: 림 솔기 호 쌍이 지워지자 없는 원통을 이름하던 branch 정점째 소멸). 실측은 dev-log 「그룹핑 팔 캐리어화」.
  ★★ **miss-first (2026-08-31): 라벨 광선이 비켜간 원통을 0으로 센다.** 광선 도로가 원통 면의 `Band` 파괴자를 풀이보다 먼저 요구하던 순서를 뒤집어 — `SpanAsk{Band{span,half}|MissOnly}`
    한 철자, MissOnly의 답은 구조적으로 {0, 기권} — 갈린 결과의 라벨이 보스를 비켜가는 광선으로 결정한다(전제: per-solid dissolve가 `(NodeId, ClassIx)` 키로 옆면을 받음). mid 15셀 →
    **Ok(2) 정확 부피**(NoClearRay 56→41·Ok 86→101), corner Cut×mid 잔류 = hit-기권의 검출기, 프로덕션 census **공집합**. 실측은 dev-log 「miss-first」.
  ★★ **호 라벨이 자기 셀을 싣는다 (2026-08-31).** 잘린 원의 per-arc 라벨을 «클래스의 저장 프레임 × 축»(`axis_up`)으로 고르던 것을 **셀에 직접 묻는 절대 판정**으로 바꿨다 — 셀은 원을
    가로지를 수 없으므로 유리 코너 하나의 반경 부호가 그 셀의 쪽이다(전 스위트 스위프 호 4,486 중 **116이 반대편을 싣고 있었다**; 양쪽 다 침묵하는 폴백 200, «양쪽 witness» 0). 두 픽스처로 세운 규칙을
    세 번째가 반증했고(같은 평면·같은 법선·같은 두 호인데 디스크 쪽이 반대 반쪽-간선), 잠금은 결론이 아니라 **전제**를 감시한다(`disk_side_probe`). 함께 룰링의 θ 스테이션이 림 목록과
    무관하게 놓인다(`end_other` −88 = `end_exact` +88). 가로지르는 census Ok 101→**117**·CGU 46→**22**, 재연산 Ok 34→**38**(CGU 0), 프로덕션 census 공집합. 실측은 dev-log 「호 라벨」.
  ★★ **차트가 자기 세로선도 읽는다 (2026-08-31).** `read_cell`이 셀의 **네 변 중 가로 둘**만 읽던 것을 대칭으로 — 룰링의 딱지가 «벽 양쪽의 재료»이므로 림의 딱지가 위/아래를 답하듯 옆을 답한다
    (`RulingExtent::label`/`ThetaSeg::label`이 `cfg(test)`를 잃었다 — 그 doc이 예약한 날). 부호는 새로 안 만들었다: `ruling_sweep`의 인라인 `−(κ·side)`를 `plus_theta_is_above` 원자로 들어올려
    세 소비자가 공유. **채우되 뒤집지 않는다** — 그림자 대조 일치 4,036 / 불일치 218이고 218은 둘 다 이름이 있다(`corner Common`의 코너 딱지 결함 120 · 한 벽의 두 룰링이 다른 딱지 = **차트에
    빠진 룰링**). ★ 함께 디스크 쪽 규칙의 **빠진 인자**를 복구했다(`axis_up · frame_sign`, 기하 감시자와 3,260:0) — 반증된 것은 규칙이 아니라 인자였다. 실측은 dev-log 「차트가 자기 세로선도」.
  ★★ **섹터가 정거장을 가로지르면 걸친 호들을 다 읽는다 (2026-08-31, `90006fb`).** 차트의 z-선은 **모든** ⊥ 클래스에서 오는데 θ-선은 옆면 스윕의 조각뿐이라, 벽의 두 룰링 중 **면이 안 닿는 쪽 선이 없다** ⇒ 셀이 실재하는
    θ 정거장을 가로지른다. `arc_around`가 `None` 대신 **가장 가까운 rim 노드 사이의 호들**을 돌려주고 `read_cell`이 **만장일치일 때만** 답한다(전체-원 팔의 규칙을 일반화; θ를 모르는 `by_span` 폴백의
    과대주장을 대체). 가로지르는 census `CylinderGateUndecided` **22 → 0**, Ok 117 → 137. 차트의 θ 축을 **벽 주도**로 세우는 구조 수리는 별도 칸(2-A/2-B)이라 적었으나 ★ **D5가 반증했다** — RBNY 15칸에서 θ-선은 이미 있었고 없는 것은 «구간에 묶이지 않은 방출»이었다(아래 여덟째 계단). 실측은 dev-log 「섹터가 정거장을」.
  ★★★★★ **잘린 원의 셀은 온전한 원과 같은 질문을 받는다 (2026-09-01, `d861571`).** 중첩 판정의 증인은 **링의 코너**에서만 왔고(`three_plane_probes`는 branch 노드를 버리고 `node_coords_rat`은 `Branch`면 점을 보지도 않는다), 원통에서 잘린 링은 코너가 전부 branch라 **탐침 목록이 빈 채** 「광선이 막혔다」는 이름으로 거절됐다. 게이트가 모든 실린더 쌍을 반지름 합보다 멀리 떼어 놓으므로
    **호로만 닫힌 링은 그 원 자체**이고(`ring_own_circle`), 증인은 **중심**(축 방향과 무관하게 유리수 — `circle_centre_rat`). `node_in_circle`도 같은 규칙을 받았다. 그리고 사실에 이름을 줬다 — **`RingHasNoWitness`**(「증인이 애초에 없었다」)는 `NoClearRay`(「가진 증인이 다 막혔다」)·`WitnessNotRational`(「값을 못 만들었다」)과 다른 사건이다.
    가로지르는 census Ok 146 → **154**, `ring_in_ring`은 이제 **아무것도 안 낸다**. 실측은 dev-log 「잘린 원의 셀은」.
  ★★★★★ **링의 간선이 자기 중점을 준다 (2026-09-01, `997f0b5`).** 남은 20칸은 벽 **패널** 링이라 코너 넷이 다 branch 이름이었다. 증인은 코너가 아니라 **간선**이다 — `plane_plane_cylinder`가 두 근을 **한 `mid`에서** 만들므로(`lo=(mid,−half,disc)`·`hi=(mid,+half,disc)`) 두 끝이 그 쌍인 간선의 중점은 **`base + s.a()·dir`**, 온전한 유리수이고 `disc > 0`이라 엄격히 사이다(`chord_midpoint_rat`; 켤레성은 값이 아니라 **이름**으로 판정 — 같은 canonical 쌍·같은 `cyl`·뿌리 `{Lo,Hi}`).
    그리고 증인을 **여럿 쓰려면 «기권»과 «실패»를 갈라야 한다** — `circle_center_in_ring`의 몸통을 `rational_point_in_ring`(`Result<Option<bool>, _>`)으로 들어올렸다: 「이 점이 링 위」는 **다음 증인이 remedy**고 「값을 정확히 못 만들었다」는 아니다. 들어올린 문은 **가장 센 도로**로 디스패치한다(유리수 다각형은 y에서 half-open인 차트 도로). 가로지르는 census `RingHasNoWitness` **20 → 0**·Ok 154 → **174**(20칸 전부 빌드), 그 이름은 워크스페이스에서 **아무도 안 낸다**. ★ 직전 항의 예고(「`QuadVal` 넓히기가 20칸을 닫는다」)는 **반증**됐다 — 그 노선은 질문을 기권하는 약한 도로로 내린다. 실측은 dev-log 「패널의 간선이」.
  ★★★★★ **잘린 캡도 원판이고, 자기 내부의 점을 이름 짓는다 (2026-09-02, `8fa54ed`).** 3D 성분 깊이 판정의 증인 공급(`coord_probes`)이 **온전한 원**인 캡만 읽어, 벽 보스의 Common을 슬랩이 가른 성분(반원 캡 ×2 · 패널 · lateral, **정점 0개**)이 탐침을 하나도 못 냈다 — 가로지르는 census `NoClearRay` 16 중 **10**이 그것. `face_circle`이 「이 링이 온전한 원인가」(`ring_own_circle`) 대신 **「이 면을 어느 원이 두르는가」**를 묻고(호들이 한 실린더 · 평면이 축에 수직; 호가 없는 패널은 걸러진다), 증인은 그 원판 평면의 점이되 **어느 것이 면 안인지는 링에게 «묻는다»** — 유도하려면 「현의 재료 쪽」 부호가 필요한데 이 도로엔 그 오라클이 없다.
    후보는 중심 + 차트 축·**대각** 여덟 걸음, `λ = r/(|e|²+1)`이라 `λ²|e|² ≤ r²/4`(유리수 부등식 — 탐색 아님). ★ 대각이 필수다: 현이 한 축과 나란하면 그쪽은 경계에 남고 다른 축의 광선은 **정확히 seam**을 때린다(축만으로 16개 면이 전부 기권). `NoClearRay` **16 → 6**·Ok 174 → **184**. ★ 남은 6(corner 가족)의 뿌리는 **재서 확정**됐다 — 전부 `point_in_mixed_ring`의 **「광선 위의 코너」 동률**이고 옆면의 `MissOnly`와는 무관하다(커밋 직후 적었던 반대 진단을 실측이 반증했다 — 소진 호출이 **마지막으로** 죽은 자리의 측정; 코너가 판정이 된 뒤 ②-b가 시도별로 재니 그 앞의 옆면 기권이 남은 4칸의 뿌리였다). 곁다리로 `point_in_faces_rat`의 `DegenerateRing` 가드가 **혼합 포크 위**에 있어 「호+현」 두 간선 링을 거절하던 것을 아래로 내렸다. 그리고 `point_in_mixed_ring`의 **호 팔은 그날까지 한 번도 안 돌고 있었고**(공급이 없어 광선이 안 나갔다), 그래서 세 팔이 한 규약을 쓰도록 먼저 고쳤다(`64a46c9` — seam 팔 둘이 「근이 호 끝」을 조용히 `false`로). 실측은 dev-log 「잘린 캡이」.
  ★★★★★ **감김수는 «영역의 극점»에서 읽는다 (2026-09-01, `fe00e0f`·`035db7a`).** `loop_winding`은 turn을 **사전식 최소 노드**에서 읽고 그 doc이 「노드 집합의 극점이니 링의 껍질 꼭짓점」이라 논증하는데, **호가 있으면 거짓**이다(호가 모든 노드를 넘어 부푼다). 축이 벽 평면 위인 보스는
    바깥 부채꼴이 중심에서 **270° reflex**라 부호가 뒤집히고, 두 셀이 **맞바뀌므로** `#(−1)==성분 수`도 Euler도 못 잡는다 ⇒ void seed가 유계 셀에 꽂혀 한 솔리드의 재료가 **전역으로 뒤집힌다**. 고침은 `CurvedStraightRun`의 doc이 이미 이름 붙인 것(`arc_extremum_winding`).
    실측: 참 위반 1,465/9,102, 그중 답이 어긋나는 184는 **전부 거절되던 corner-lo** ⇒ 출하된 오답 없음. 함께 절대 명제 둘을 세웠다 — **`disk_side_agrees`**(전-간선 검증이 «한 솔리드 두 비트 전역 뒤집힘»에 구조적으로 눈멀어 있다)와 **`first_volume`**(첫 결과의 부피가 자기 자신을 기준 삼고 있었다). 실측은 dev-log 「감김수는 영역의 극점에서」.
  ★★★★★ **옆면은 차트의 «영역»이다 — 여덟째 계단 D5 (2026-09-02, `2fb7960`·`938fdb5`·`b9333de`·`c897e5b`·`f99095c`·`d458730`).** 가로지르는 census의 `RulingBoundNotYet` 15칸이 **전부 한 자리**(방출기의 런 끝 `node_on`)에서 죽었고, 탐침으로 재니 그 코너는 **어느 클래스에도 없는 노드**였다(z-선이 거기서 관통, 벽 노드 0/15). 뿌리는 그 자리가 아니라 **어휘**: 방출기가 «구간마다 띠 아니면 네 노드 패널»로만 말해서 구간을 넘어 이어지는 영역(예: 벽 보스 Fuse × 축 통과 벽 = 8-간선 링 하나)을 철자할 수 없었다. 고침은 평면 쪽의 단계를 차트에서: 방출 셀 → **연결 성분** → 경계 → 경계의 런을 **이웃 클래스의 조각**(잘린 림의 `CutRim.nodes` 사이 호·벽의 룰링 조각)으로 → 사이클 → `classify_cycles`(띠+구멍) 또는 `Bound::Ring`(최저 구간 아래변을 담은 사이클이 외곽). 띠/패널 디스패치·«양 림 절단»·런 사다리·D4 병합 패스(`unify_curved_faces`)가 사라졌고 잘린 림은 **한 철자**(`Rim::Chain`)다. ★ **«같은 엔진»은 같은 단계이지 같은 코드가 아니다** — 평면 DCEL(`walk_cells`의 각순서·`nest_cells`의 `point_in_ring`·`label_cells`의 비유계 루트) 이식은 검토가 코드로 반증했고, 차트는 환면·직교이며 딱지를 읽는다. ★ **정점은 이웃 클래스가 가진 곳에만** — 경계가 직진하는 정거장에 정점을 만들면 벽 면의 간선 하나에 옆면 간선 둘이 대응해 용접이 깨진다. 실측: RBNY **15 → 0**, Ok 184 → **197**, `NoClearRay` 6 → 8(corner 가족의 Fuse × mid — 가른 조각의 nesting), 프로덕션 census 269행 **비트 동일**, 직렬 원장 거절 방출 19 → 0·면 1114 = 성분 1114. 곁다리 둘: 1a에서 정거장을 «다른 선의 노드»로 대표하던 `arc_around`가 `circular_order`의 «한 점 두 이름» 거절을 부르던 것(Q2, 128)을 이름(`crossing_on_ruling`의 canonical `NodeId`) 등식으로 고쳤고(`end_other` 321 → 245, `other_present` 60 → 0, `src0_present` 14 → 0), 남은 245는 전부 면 span 밖 부재 셀의 옳은 «말 안 함»(`OtherWhy::WholeDisagree`)이라 «present 셀은 `Other` 끝을 갖지 않는다»를 기록 자리 단언으로 승격했다. 실측·검토가 바꾼 것·빗나간 예측은 dev-log 「옆면은 차트의 «영역»이다」.
  ★★★★★ **광선 위의 코너는 «위로 떠나는 스텝»이 센다 — 혼합 링의 반열림 규칙 (2026-09-02, `908912f`·`8580b3e`·`2c74c49`).** D5 뒤 가로지르는 census의 `NoClearRay` 8(corner 가족 전부)은 `point_in_mixed_ring`이 «광선 위의 코너»를 **기권**해 3D 깊이 판정의 탐침 공급이 소진된 것이었다 — 평면 도로 `point_in_ring_2d(_rat)`은 같은 코너를 y-반열림으로 늘 **판정**해 왔다. 규칙은 하나다: «광선 위의 코너는 그 코너를 **위로 떠나는** 스텝이 센다»(직선 변은 다른 끝이 엄격히 위, 호는 그 끝의 접선이 위 — 이웃 둘이 반대쪽이면 1, 같은 쪽이면 0 또는 2). 결정부를 `nacre_geom::intersect::ray_step_crossing`(+`ray_straddle`) **한 철자**로 두고 평면 도로 둘·혼합 직선 팔·혼합 호-끝 팔 넷이 읽는다; 호 끝의 이탈은 `arc_departure_side`의 규약(CCW 접선의 쪽 = −`ruling_side`)이라 `ray_step_crossing(Zero, ty, ty)`로 환원되고 둘째 외적은 없다. 남는 기권은 **탐침이 링 위**(코너·광선을 따라 놓인 스텝·변 위·호 근), 접선 광선, seam 근, 수평 접선뿐. 실측: `NoClearRay` 8 → **4** · Ok 197 → **201**(Fuse·Cut × mid 개통), 판정 호출 214 전부 **첫 탐침**(이전 24건이 재시도), census 벽시계 100 s → 21 s, 프로덕션 census 269행 비트 동일. 잠금: 진리표 27행 · 유리수 링 22개 격자 2,278점의 **두 도로 정확 대응** · 코너 물림 픽스처(`with_cut_rings([40,24,20], seam_off)`, seam을 코너에 둔 빌드와 뗀 빌드 둘)에서 **한 끝**만 지나는 광선(3-4-5 점 (37,24); seam-off에선 (40,15)도) — 두 끝이 한 광선이면 전역 부호 오류가 패리티를 지키므로 부호의 유일한 감시자이고, 두 빌드가 호-끝 분기 셋(seam 팔 둘·`(false, false)`)을 다 돌린다. ★ **예측 «8 → 0»은 빗나갔다**: 남은 4(corner × Common × 슬랩)는 규칙이 아니라 **탐침 공급** — 사분 원통의 증인(축 위의 점 = 판 코너 간선, 벽 간선 중점)이 전부 판의 면 링 **위**라 옳게 기권한다(`ProbeAtCorner` 6/6). 다음 계단은 볼록한 잘린 캡의 **코너 평균**을 섹터 내부 증인으로 — ★ **②-b가 반증**: 시도별로 재니 뿌리는 탐침이 아니라 다른 반쪽의 **옆면**이 교차점 물음에 답하지 못하는 것이었다(아래 문단). ★ 단계 0의 계기가 플랜의 «호-끝 팔은 코퍼스가 못 닿는다»를 반쯤 반증했다 — 물린 판은 못 닿고 corner × Common은 닿으며, 직선 팔의 기권이 그것을 가리고 있었다(`ArcRootAtEnd` 0 = «먼저 말한 팔이 있었다»). 실측·검토·빗나간 예측은 dev-log 「광선 위의 코너는」.
  ★★★★★ **옆면도 자기 차트에서 같은 패리티로 답한다 — `BoundEdges::Lateral`이 경계 루프를 든다 (2026-09-03, `3c7d3ed`·`5a97a79`·`630b813`·`9ef2a3e`).** 칸 ②가 남긴 `NoClearRay` 4의 뿌리는 «탐침 공급»이 **아니었다** — 시도마다 재니 축 위 정점의 가로 광선 12/18이 다른 반쪽의 **옆면**에 부딪혀 «이 교차점이 내 면 위인가»를 모른 채 기권(`SpanAsk::MissOnly`)했고 코너 동률은 마지막 시도 6뿐이었다(스위트: MissOnly-적중 69 전부 기권, 띠 팔 적중 **0** — 게이트가 원통 안의 ∥ 벽을 거절하므로 인구 없음이 아니라 거절된 인구). D5가 철자한 경계(`Bound::Ring`의 정거장·`Wall::Arc`/`Ruling`, `Band{Rim,Rim}`, 구멍)를 광선 도로가 버리고 있었다. 이제 원통 면의 `outer`는 **모든 루프**(`LateralLoop::{Circle(class), Ring(edges)}` — rim·사슬·패널·구멍 함께; 환면에는 «안»이 없어 `material_of`가 아니라 **총 패리티**)이고, `loop_parity`가 X에서 +z 광선으로 루프를 홀수 번 가르는지 답한다: rim은 위에 있으면, 호는 위에 있고 span이 θ_X를 품으면(끝이 정확히 θ_X면 칸 ②의 코너 규칙 — `lo`만 센다), 룰링은 절대 안 세며 X가 그 위면 경계(룰링은 `plane_side(fc) == 0 ∧ ruling_side == side`, 정거장의 이름 그대로). 띠 팔은 두 원 루프의 특수 경우라 삭제됐고(격자 2,180 광선 근 전부 동일), span 판정은 `arc_span` 한 철자를 평면·옆면 도로가 같이 부른다. 실측: census `NoClearRay` 4 → **0** · Ok 201 → **205**, 판정 222건 전부 첫 탐침, `ProbeAtCorner` −24(셋째 시도가 안 돈다), 프로덕션 census 269행 비트 동일. 잠금: 옆면 격자 오라클(관통 빌드 + **계단 빌드** — 아래 캡이 판 안이라 아래 경계가 사슬 rim; 세계 좌표의 진리와 근마다 대조, 정거장 열·seam 열 포함). ★ 부정 대조가 플랜을 반증했다: `AtLo`/`AtHi` 맞바꾸기는 반열림의 다른 관례라 닫힌 루프의 패리티를 못 바꾸고(칸 ②의 뒤집기가 빨갰던 건 직선 스텝과의 **불일치**), «양 끝 다 세기»도 노치(두 호가 같은 끝을 공유)에선 안 보인다 — 계단만 본다. 실측·검토·빗나간 예측은 dev-log 「옆면도 자기 차트에서」.
  ★★★★★ **벽∩옆면은 «한 이차식»이다 — 오프셋 벽은 게이트의 기록이 열고 배열은 기록만 믿는다 (2026-09-03, `004b0e1`·`37cd130`·`0bdd144`).** `WallMeetsLateral`의 오프셋 팔(축과 나란한 벽이 축에서 `0 < d < r`, 룰링 `±√(r²−d²)`)은 «들어 올리면 `OpenResultShell`»이라는 **D5 이전** 측정 위에 서 있었다. «축을 지난다»(`class_through_axis`)가 다섯 자리(게이트의 기록 팔·`chord_on_class`·`rulings_on_class`·`crossing_on_ruling`·스캔의 구멍-원 팔)에 철자돼 있었고, 산술(`plane_plane_cylinder`의 QuadVal 근·`ruling_side`·`circular_order`·D5·②-b)은 지름을 전제한 적이 없다 — 뿌리는 **어휘**. 이제 게이트의 `point_plane_clearance_rat(coeffs, o, r)` **한 호출**이 유일한 철자다(`Positive` 비킴 · `Negative` → 면이 판단, 못 비키면 **crossings**에 기록·통과 · `Zero` 접선 → **tangencies**에 기록·통과, 칸 ⑥); 기록 집합은 «평면이 r 안이고 면이 띠를 비키지 못했다»의 증명 캐시라 배열의 네 자리는 기록만 믿고 전제를 `debug_assert`로 든다(검토가 첫 안의 새 술어를 «한 사실의 다섯째 철자»로 반증). `class_through_axis`는 삭제. 실측: 보스 12/12·보어(슬랩 면이 보어를 가로지름·전폭 면이 갈린 보어를 가로지름)·`rul offset`·얼린 픽스처 전부 조립·validate 0·정확 부피(유리수·무리수 룰링, 축이 판 밖·안); 가로지르는 census +3 가족(45칸 전부 Ok, Ok 205 → 250), 재연산 census +3행, 프로덕션 census는 `rul offset` 3행만 이동, reject-trace의 `wall_meets_lateral`은 접선만. **둘째 물음**: 두께 0.2의 현-조각 기둥을 슬랩이 가르면 그 캡은 코너가 branch 이름뿐이고 «중심 + r/2 걸음» 후보가 조각 밖이라 **탐침 목록이 빈다** — 이름을 `RingHasNoWitness`(«증인 없음»; 링 또는 성분)로 가르고, 잘린 캡이 **현마다 둘**을 더 내놓는다: `o − t·n`에서 `q = (n·o + d)/|n|²`, `T = r²/|n|²`, `t_far = 2qT/(q² + T)`(현 너머·원 안), `t_near = q/2` — 근호 없는 유도, 어느 쪽이 면 안인지는 링이 답한다. 부피 오라클 `disk_in_rect`는 원판∩축 정렬 직사각형의 **닫힌 일반식**(현 길이의 적분). 잔여: 현의 발 선 위에 둘째 현이 놓인 사분-조각(«증인 없음»으로 온다), `Probe::Branch`(branch 코너 자체를 증인으로 — ~120–180줄, 막는 원시 술어 = 교차점이 branch 원점의 앞인가 뒤인가). 실측·검토·빗나간 예측은 dev-log 「벽∩옆면은」.
  ★★★★★ **불리언은 강체 운동과 교환한다 — 대조군 하나, 뒤집힌 부호 셋 (2026-09-03, `7fa6f87`·`bd07f99`·`a422cbe`·`cac9b90`).** 원통 불리언의 코퍼스는 «만든 자리 그대로»만 재고 있었다(census의 `rot`는 평면뿐, `trc`는 이동한 원통뿐). 대조군 하나 — `boolean_commutes`: 판·보스를 같은 등거리로 옮겨 «옮긴 결과 == 옮긴 입력의 결과»(이름·정렬 부피·면/간선/정점/공동 수·validate, 그리고 피벗 0 사분 회전이면 **정점 비트 동일**, 강체·dyadic 이동이면 branch 정점을 뺀 비트 동일) — 를 22가족(BOSS 17·enclosed·planar·bore offset/axis·through mid) × 3 연산 × 운동군 15에 걸었더니 990칸 중 150이 갈렸고 세 자리로 귀속됐다. ① ∥ 룰링 딱지 `ruling_interior_is_even`이 `side·κ`만 읽고 **차트 프레임**(`frame_sign`)을 안 건넜다 — ⊥ 도로가 2026-08-31에 얻은 인자의 쌍둥이(105칸; 이제 `plus_theta_is_above ≠ (frame_sign > 0)`, ⊥의 `axis_up ≠ (frame_sign > 0)`과 글자 그대로 대칭). ② `arc_extremum_winding`이 축이 세계 z일 때만 극점을 유리수로 읽어, 축이 ±x/±y로 돌면 부푼 270° 호를 못 보고 `lo`의 turn을 읽었다 — 3/4 원판이 root가 되어 딱지가 안팎으로 뒤집힘(9칸; 이제 «원이 걸치는 첫 세계축» a에서 `m[a]=0`이면 `c − r·ê_a`, 반쪽은 원 자신의 평면 `m × ê_a`로 — seen_ccw 보정 소멸; 직렬 원장의 undecided 5860 → 16). ③ 수송: `motion_is_exact`가 이동의 정확성을 **회전 전** 좌표로 재고, 부정확이면 `chain_motion`이 정확한 회전은 노드 없이 건너뛰고 이동만 기록 — 진실(회전 전 정의 + [T]) ≠ 캐시. 릴리스에선 오프셋 보스가 판과 «안 닿는다»고 답했다(조용한 오답, census `mot xport` 행이 실측). 이제 **법칙** `transform(rigid(R, t)) ≡ transform(T) ∘ transform(R)`: `carry_of`가 솔리드당 한 번 후보 [전부, 회전만, 없음]에 같은 탐침 «실현 == 정확 상»(정점·평면 원점·원통 축 원점)을 물어 첫 통과를 수송하고 나머지를 기록한다(36칸 — planar 쌍의 조용한 발산 9 포함). ★ **frame −1의 메커니즘은 회전이 아니라 세계 평면 씨앗이다**: `Model::new()`의 x=0·y=0·z=0(캐시 −축) 위에 +축 바깥으로 앉는 면 — `t(−4,−4,−2)` 하나로 max-side 9가족의 ∥ 딱지가 전부 뒤집혔고(회전 없음), min-side 가족은 그 이동으로 씨앗에서 떨어져 frame +1이 된다(두 계급이 상보). 잠금: 상시 부분군(계급마다 운동 하나, 396칸 ≈ 42 s)·ignored 전체군(990칸)·`KNOWN` 원장(비어 있음 — 다음 발산이 적힐 모양)·census `mot` 18행·법칙 픽스처(원점 (4.3, 0.5) + rigid(rz90, (−4,8,0)): 두 길 동일). 남는 것: 극점이 무리수인 축(`m[a] ≠ 0`, M6-3 첫 칸 — x축 둘레 기울기 `(0,s,c)`는 이미 덮인다)·`rotated`의 정의(«옮긴 상이면 참»)는 M6-3에서 재정의 후보·원통의 `Mirror`. 실측·검토·빗나간 예측은 dev-log 「불리언은 강체 운동과 교환한다」.
  ★★★★★ **접선 접촉은 결함이 아니었다 — 집힌 것은 표면이 아니라 면이다 (2026-09-04, 칸 ⑤).** 남은 항목 «접선 `d = r`(능력 C)»를 조사했더니 **세 기록이 반증됐고 그중 하나는 이 칸이 세웠던 계획 자신이다**(칸 ⑱과 같은 모양 — 지으려던 것을 계측이 반증했다). ① design.md의 «접선은 들어올리면 조립된다, 막는 것은 C다»가 **D5 이전 측정**이었다: 게이트의 접선 팔을 들어내면 12칸 중 한 칸도 안 조립되고, 사다리를 채워 가며 재니 `chord_on_class`(현이 아니라 한 점) → `crossing_on_ruling`(룰링이 하나) → **걷기의 `UnorderedEdges`**(= `RulingCarrier::side`의 `±1`이 다섯 구조체의 동일성 키인데 접선은 룰링이 하나다)로 이어진다 — 앞의 둘은 `Tangent` 팔 한 줄씩이고 **셋째가 어휘 결정**이다. ⇒ 접선은 **B의 마지막 인구**이고 C와 별개다. ② 칸 ㉓이 열어 둔 «그 둘이 유효한 솔리드인가»의 답은 **유효하다**: 접점의 **링크가 원 하나**이고(집힌 면의 두 로브가 이웃 곡면을 돌아 이어진다) 재료도 국소적으로 한 덩이이며, **OCCT가 같은 피연산자로 부피·면적·면 수가 일치하는 몸통을 돌려준다**(ε 비킨 쌍둥이와도 면 수가 같다 — 상용도 집힌 면을 안 쪼갠다). 집힌 것은 **면**이다. ③ 그래서 «면의 두 경계가 닿으면 거절»이라는 (이 칸이 세웠던) 계획은 **정당한 결과를 거절**한다 — design.md가 이미 *「면 단순성 검사를 넣으면 … 정당한 결과를 거절하게 된다」*라 적어 뒀고, `SelfTouchingResult`의 참 인구는 **두께 0**(45° 폴드)이지 두께가 있는 이 모양이 아니다. ★★★★★ **χ는 증거가 아니다**: 두 결과의 χ는 짝수지만 «핀치 하나는 χ를 홀수로»는 **정점인 핀치**의 문장이고, 접선은 정점을 안 만들므로 `euler_counts`가 그 점을 못 세 판별력이 0이다 — 짝수 χ를 «핀치 없음의 증명»으로 읽으면 거짓 규칙이 되고, 이 칸이 반증을 받아들이는 도중 한 번 그럴 뻔했다. ★★ **처방은 b-rep이 아니라 메시**: 재봉합은 좌표 접촉을 인덱스 접촉으로 바꿔 `monotone::link`가 다시 거절하고(다섯 철자 전부 칸 ㉓이 `DegenerateRing`으로 실측), 분할은 유효한 솔리드를 바꾸는 churn이다 — 남은 능력은 「남은 능력」이 이미 이름 붙인 **스윕의 SoS**이고, `tessellate`가 모델 전체를 돌므로 **집힌 면 하나가 세션의 모든 몸통을 지운다**(소비자는 플레이그라운드; STEP은 무사). 잠금: `bands`의 `the_touch_is_a_slit`(정점 0 ∧ 원이 온전) · `nacre-oracle`의 `#[ignore]` 측정 둘. 실측·검토 여섯 라운드·빗나간 예측은 dev-log 「접선 접촉은 결함이 아니었다」.

  ★★★★★ **접선은 «가르지 않는다» — 게이트는 기록을 멈추고, 판정은 `keep`이 한다 (2026-09-04, 칸 ⑥).** 능력 B의 마지막 인구가 열렸다. **뿌리는 게이트 한 줄**이었다: 엔진은 «a tangent line touches **without separating**»를 이미 여섯 자리에 적어 두었고(`lateral_crossings`·`circle_crossings`+`split_circles`의 `Double` skip·광선 둘·branch 실현 둘), 그와 다른 말을 하는 자리는 `planes.rs`의 거절뿐이었다. ⇒ `point_plane_clearance_rat` 한 호출의 **삼분법이 그대로 세 기록**이다: `Positive` 무기록 · `Negative` → `crossings` · `Zero` → **`tangencies`**. 접선을 `crossings`에 **안** 넣는 것이 핵심 — 그 집합의 명제(«평면이 반지름 **안**»)와 그것을 드는 **세** `debug_assert`(`crossings.contains`의 세 독자가 각각 하나씩)가 전부 그대로 참이고, **배열의 추적·분할·차트·걷기·라벨·방출이 한 줄도 안 바뀐다**(실측: 21칸 중 **15칸 조립**, validate 0, 부피 정확 — 나머지 6칸은 접선을 품는 셋째 평면이 있는 반례 A·B로, 배열이 `CoincidentNodes`로 먼저 거절한다). 칸 ⑤의 «들어올리면 아무것도 조립 안 된다»는 `crossings.insert`를 **남긴 채** 잰 사다리였다.
  **셋을 가르는 것은 연산이고, 연산은 `keep`으로만 들어온다.** 접점 근처 재료는 정확히 세 영역 — 렌즈(원통 안) · **쐐기 둘**(포물선과 평면 사이) · 건너편 — 이고 인접은 `L–W2`·`F–W2`뿐이다(`v = 0`에서 쐐기 둘 사이는 렌즈 안이다). 그래서 «갈라지는가»는 `(W2 ∧ ¬L ∧ ¬F) ∨ (L ∧ F ∧ ¬W2)`이고, 각 영역의 `(in_A, in_B)`는 **벽 면의 재료 쪽**과 **옆면의 `orient_sign`** 두 비트가 정한다. Fuse는 렌즈가 다리, Common은 렌즈만, Cut은 쐐기 둘만 — **연산 이름을 한 번도 안 쓰고** 세 답이 나온다.
  ★★★★★ **«갈라진다»가 곧 결함은 아니다 — 실측이 규칙을 두 번 고쳤다.** ① 바깥 접선의 Fuse(`L ∧ F ∧ ¬W2`)는 **선으로 닿는 유효한 몸통 둘**이고 커널이 이미 그렇게 낸다(`Ok(2)`, validate 0, 메시 OK) — 거절이 아니다. ② 그런데 그 둘이 **다른 데서 이어져** 있으면 한 솔리드의 자기 접촉이다: 노치 벽에 접하면서 `y`로 블록과 겹치는 보스가 `Ok(1)`·validate clean·mesh ok로 **조용히 나왔다**. ⇒ 판정은 «두 덩이 ∧ 벽 클래스와 원통 클래스가 **한 출력 솔리드**에» — 첫째 항은 그 성분 조건이 항상 참(두 쐐기가 같은 두 면을 경계로 쓴다), 둘째 항이 진짜 물음이다. 그래서 판사는 게이트가 아니라 **`boolean::tangency_reject`**, `self_touch_reject` 바로 옆(민팅 전, grouping이 손에 있는 곳)에 선다. 이유는 새 이름이 아니라 **`SelfTouchingResult`**(그 doc이 이제 증인 둘 — 간선 하나, 선 접선 하나 — 을 든다); `WallMeetsLateral`은 생산자가 없어 **삭제**했다.
  «접촉이 선인가 점인가»는 벽 면이 접선을 가로지르는지로 가른다(칸 ⑤이 «점 접촉은 유효»를 확정했으므로). ★ 스팬 안쪽을 함께 요구했더니 얼린 `bore-slab`(슬랩이 보어와 같은 높이)의 코너가 **열린** 구간의 끝에 전부 놓여 전부 버려졌다 — 축 겹침은 「못 비켰다」가 이미 보장하므로 스트립 양쪽만 본다.
  실측: **막힌 스터드 픽스처**(밑면을 `z = 0`에 둔 것)가 Fuse `1 + 0.06π`·Common `0.02π`로 정확하고 Cut은 이름 붙은 거절, **OCCT가 부피·면적·면 수(8)까지 일치**(관통 쌍둥이 `0.35`는 9면·다른 부피 — 판별력 있는 대조). ★ 이 픽스처를 «사용자의 스크립트»라고 불렀는데 **다른 솔리드다**(칸 ⑦b가 잡음): 앱의 `cylinder({center})`는 **중간 높이** 앵커라 사용자의 스터드는 `z ∈ [−1, 1]`로 큐브를 **관통**한다 — 솔리드 1·validate 0·**면 10**·부피 `1 + 0.04π`, 집히는 캡 **둘**(칸 ⑦c 실측). 스위트에서 움직인 것은 **예측한 잠금 셋뿐**. **남는 것**: 이 칸이 «유효하지만 못 그리는» 인구를 키운다 — 사용자의 Fuse도 윗면이 집혀 `SelfTouchingBoundary`(조립된 15칸 중 7칸). 처방은 칸 ⑤이 확정한 **메시의 SoS**다. 접선을 품는 **셋째 평면**이 있으면 국소 영역이 여섯이라 규칙이 안 서고, 그때는 기권한다(구성한 두 경우는 배열이 `CoincidentNodes`로 먼저 거절한다). 실측·반증 다섯·검토는 dev-log 「접선은 가르지 않는다」.

  ★★★★★ **접촉은 자기 자격을 증명한다 — `self_touch`가 교차를 가리지 않는다 (2026-09-05, 칸 ⑦b, `763dc8c`).** `self_touch`는 첫 접촉에서 즉시 반환하고 교차 루프는 접촉이 없을 때만 돌았다 ⇒ ① 접촉이 다른 자리의 교차를 **가리고**, ② 교차점이 홀의 꼭짓점에 정확히 놓이면 **접촉으로 흡수**된다 — 깨진 링 집합이 `SelfTouchingBoundary`(«이 분해법에 답이 없다»)로 나왔다(둘 다 손 픽스처로 먼저 빨갛게). **접선 증인** 하나: 접촉 정점의 두 이웃이 접촉 선분의 같은 쪽이면 접선, 반대쪽이면 교차. `Interior`(상대 링이 그 점에서 직선)에서만 완전하고 `AtEnd`·공선 이웃은 **기권**. 교차 루프는 언제나 돌고 `straddles`는 엄격(0 넷 금지) — 셋은 **원자적**(무조건 루프 + 관대한 straddle이면 접선 픽스처가 거짓 교차, 실측). 거절 이름은 **메시에 대한 말**로만 쓴다(«교차 = 어떤 삼각분할이든 틀린다, 접촉 = 이 분해법에 답이 없다» — 솔리드 판결은 `validate`의 것: 사용자의 원칙). 접촉표 실측: 면 4, 전부 `Interior:Tangent`. 그리고 **사용자의 두 원칙**이 방향을 정했다 — 「정상 형상이면 둘 다 가능해야 한다」(부분 렌더링 접음), 「실제로 접하는 형상은 시각화에서도 접해야 한다」(재샘플 기각). dev-log 「접촉은 자기 자격을 증명한다」.

  ★★★★★ **다리와 순서 — 집힌 캡이 그려진다 (2026-09-05, 칸 ⑦c, `813f933`·`b1179a9`·`d6986a2`).** 세 커밋: **A** pre-pass — `tessellate`의 1·2단계 사이에서 평면 면의 `Interior:Tangent` 접촉을 **전부 모은 뒤** 접촉된 **직선·양쪽 평면** 간선의 폴리라인에 접점 핸들을 끼운다(간선→면은 `Model::adj`가 아니라 살아 있는 면에서 직접 — 불리언 안에서 `adj`는 낡았다; 멱등; 곡선·원통·다중 접촉은 세고 거절). **B1** `polygon/sos.rs` — 쌍둥이 한정 인덱스 인식 순서: `d_o = uv[prev[o]] − T`, `d_h = uv[next[h]] − T`(**접선 방향** — 법선 성분이 있으면 수평 접촉 변에서 «쌍둥이 아닌 점은 정확 비교»가 실제 섭동과 어긋나 실현 불가 순서가 난다), lex 동점은 T가 소거되어 `lex_less(uv[prev[o]], uv[next[h]])`(정확, 동점 불가), `side`의 구조적 0은 일반 1차 항 `(b−a)×(d_c−d_a) + (d_b−d_a)×(c−a)`를 `Expansion`으로. 꿴 자리는 정렬·`left_of`·`insert`·`top`/`bot`·체인 병합뿐이고 `classify`의 turn·접촉 검출·증인·straddles·`emit`·`refine`은 정확 유지. **B2** 스플라이스 — 링 간 **중복 핸들**이 다리, 분할 링은 «이웃이 공선인 링»(기하로, 순서 아님), `a[..=i] ++ b[j+1..] ++ b[..=j] ++ a[i+1..]`, 감김 게이트는 원래 링에, `constrained`·`boundary`는 병합 링에, `decompose`는 쌍둥이의 `AtEnd` 기록 둘만 면제. ☑ **실측**: 사용자의 관통 스터드 **10/10 면·1452 삼각형·누수 0**, 막힌 8/8; `bands`의 두 접선 픽스처가 공유 오라클(`mesh_covers_faces`: 면별 넓이 vs `face_props`, 부피 vs `mass_props`, 상대 `1e-3`)을 통과; census 면제 제거(`r.is_ok()`). tripwire(넓이 0·1차 항 0·단조 `debug_assert`) 무발화 — **두 쌍둥이는 증인이 보증한 볼록 코너**라 `Split`/`Merge`가 아니고 쌍둥이 대각선이 생기지 않는다(구조적 논증). 기호 판정은 정렬 동점 16건뿐, `side_idx`의 1차 분기는 **프로덕션 인구 0**(일반 보증으로 남김). B1의 잠금: **결정적 테스트 82개/1945 테셀레이션에서 위치 정준 삼각형 집합이 바이트 동일** — census 코퍼스 «전체»는 proptest 넷(1536행)이 실행마다 달라 그 형태로는 불가능함을 실측. dev-log 「다리와 순서」.
  ★★★★★ **원통 쌍의 게이트가 면을 읽는다 — 사용자의 십자 스터드 (2026-09-05, 칸 ⑧).** `fuse(cuboid, z-스터드, y-스터드)`의 둘째 불리언이 `CylinderPairContact`였다: 쌍 루프가 `cylinders_clear`(축 거리 vs 반지름 합 — **무한 곡면**의 문장)만 물었고 축이 원점에서 만나 거리 0. 같은 게이트의 평면–원통 ⊥·∥ 팔과 평면 배열은 이미 **면** 단위였고 비스듬 평면 팔만 곡면 단위로 남아 있었으므로, 이 칸은 넷 중 **셋째** 자리를 같은 명제로 맞춘 것이다(엔진 무변, 넷째는 아래). 결정적 실험이 두 얼굴을 보였다 — 수직 쌍을 그냥 통과시키면 사용자 모델은 정확히 조립되지만(면 14·`1+0.08π`·메시) 진짜 교차하는 두 원통도 `Ok(2)`·validate 깨끗으로 **조용히 틀리게** 나온다 ⇒ 완화가 아니라 **증명**. 술어 `lateral_reach(def, span, d)`: 옆면의 도달 범위 = 사영 구간 `d·o + s·(d·m)` ± 반경 도달 `r·|d⊥|`(`|d⊥|² = d·d − (d·m)²/(m·m)`), 제곱으로 정확·근호 없음·`Rat` checked; 수직은 `d·m = 0`이라 `s` 항이 사라지는 경우이지 분기가 아니다. `lateral_faces_clear`는 A의 모든 옆면 구간 × B의 모든 옆면 도달 범위가 서로소면 «면을 공유하지 않음»(어느 방향이든 하나면 충분, 둘 다 묻는다); 평행 팔은 곡면 판정 그대로(`chart_of`의 핸들 인턴이 동축 거절에 기댄다 — 쌓인 동축 원통은 아직 거절). ★ **죽은 팔**: `d·m ≠ 0`은 오늘 프로덕션 도달 불가 — 평면–원통 게이트의 비스듬 팔이 비켰음을 묻지 않고 `ObliqueCylinderCut`으로 거절하며(**넷째 자리**), 비스듬 원통 쌍은 반드시 비스듬 캡을 데려온다. 단위 테스트가 잠그고 넷째 자리가 열리면 살아난다. ☑ 실측: 두 순서 모두 Fuse 면 14·`1+0.08π`, Common 면 3·`0.04π`, Cut `1.0`; 부정 대조(`fuse(b,c)`) 거절 유지; 스터드 위의 십자 `Ok(2)`; census 287행 두 프로파일 동일·⑦c와 바이트 동일; `reject-trace` 히스토그램은 새 픽스처 행 둘만 +1. 효율(사용자의 물음): 배열엔 낭비 작업이 없으나 후보 열거는 전수 — AABB 컬링은 접어 둔 별개 지렛대(~1.3×). dev-log 「원통 쌍의 게이트가 면을 읽는다」.

  ★★★★★ **넷 중 넷, 그리고 접선의 이름 (2026-09-06, 칸 ⑩).** 칸 ⑧이 넷 중 셋째 자리를 «면»으로
  맞췄고, 사용자의 첫 조립체가 남은 자리들을 한 번에 두드렸다 — 같은 솔리드 안의 평행 쌍, 접선
  룰링, 비스듬 평면. 게이트는 이제 모든 자리에서 «상대의 자취가 내 **발자국**(축 구간 × 각도 범위)과
  만나는가»를 묻고, 못 증명한 쌍은 기록(오늘은 거절 한 자리, M6b·M6-3이 그 자리를 «기록을 배열에
  넘긴다»로 바꾼다). 접선 룰링은 `side = 0`이라는 이름을 얻었고, 그 이름이 서자 평면 배열의 각
  순서·부분 림·앉은 면의 빚·셋째 평면 위의 코너·증인 공급이 차례로 드러나 고쳐졌다. 「남은 능력」
  칸 ⑩ 항목과 dev-log 「칸 ⑩」.

  - ★★★★★ **D3의 전제조건 — 여기서 계기 넷이 눈이 멀었던 적이 있다.** 칸 ⑦이 옆면의 면 구조를 바꾸고
    **메시가 조용히 틀린 채** `validate`·watertight·정확 부피·면 수가 **전부 초록**이었다(넓이 7.873 vs
    9.425). 그걸 잡은 눈은 칸 ⑧의 **메시 넓이 ≈ 정확 넓이**다. ⇒ 갈아타기 **전에** 원통 코퍼스에 그
    오라클을 세운다.

  ★ **「D를 하면 오프셋이 떨어져 나온다」는 그럴듯하지만 미측정이다.** 오프셋의 어려움은 룰링 위치가
  무리수(`√(r²−d²)`)라는 것이고 그 산술(`QuadVal`·근호 하나)은 배열이 **이미** 갖고 있다. D 위에서는
  「특수 룰링 케이스」가 아니라 차트 배열의 한 셀 경계가 되므로 떨어져 나올 **가능성이 높지만**,
  재 본 적이 없다. 추측으로 적어 둔다.

  ☑ `UnorderedEdges`(벽 보스 + 공면 캡)는 **아직 진단 안 했다**. 그 문장은 「한 정점에서 두 간선이
  같은 각으로 떠난다 — 상류의 `merge_coincident`/`Aliases`가 접었어야 했다」이므로, 케이스로 열 것이
  아니라 **그 상류를 봐야** 한다.
  ☑ **(E3-b/c) 진단·해소됐다**: E3의 스캔·현 합류가 열리자 그 가족(`rul flush`)은 조립까지 가고,
  Cut/Common은 빌드, Fuse는 2-곤 기권이 남긴 현(stated 평면-자기 간선)을 결과 전체 판정이
  `CoplanarMerge`로 거절한다 — 이름이 코퍼스에서 사라졌다.
- ☑ **다음 자리 — 「원통도 자기 차트에서 arrangement를 돈다」 (칸 ⑦이 남긴 방향, 2026-08-25).**
  주기 면의 바깥 walk(`band_loop`의 슬릿 다리 + 곡면 병합의 `merge_curved_group`)는 이 글을 쓸 때
  **손으로 쓴 케이스워크**였다(D2b가 방출을 차트 셀로, D4가 병합을 감김수 한 규칙으로 바꿨다). 평면 쪽은 이미 클래스마다 **셀 복합체**를 돌려 그 케이스워크를 없앴는데
  (§6.1 arrangement), 원통은 등거리 차트 `(z, r·θ)`를 갖고 있으므로 **같은 엔진을 그 차트에서**
  돌릴 수 있다 — 솔기는 차트의 절단선일 뿐이고 노치·구멍·θ-패널은 전부 그 배열의 셀이 된다.
  tess가 이미 「모든 면이 같은 차트 길을 탄다」로 갔으므로(§5), ops만 남았다. **아직 안 한
  이유**: M6-2의 룰링 도로가 먼저였고, 손으로 쓴 walk가 지금 인구에서는 초록이다.

## 7. 검증·오라클 인프라 (`nacre-validate`, `nacre-oracle`)

1일차부터 CI에 들어가는 것들:

```rust
// 모든 연산 직후 자동 실행되는 불변식 (디버그 빌드에서 강제)
pub fn validate(m: &Model) -> Vec<Violation>;
// - 오일러-푸앵카레: V − E + F = 2(S − G) + L_i  (L_i = 면 내부 루프 수.
//   관통구멍 정육면체로 검산: V16 − E24 + F10 = 2 = 0 + L_i(2) ✓.
//   내부 루프 항을 빠뜨린 축약형을 쓰면 정상 모델을 불량 판정하게 된다)
// - 모든 Loop 닫힘, half-edge 짝 정합, 방향 일관성
// - 면 orientation ↔ 외곽 루프 감김 정합 (FaceMisoriented): 루프의 뉴웰 벡터와
//   면이 진술한 바깥 법선(plane.normal × orientation)의 정렬 — 불리언
//   collect_planes 의 debug_assert 가 지키던 불변식의 release 그물 (2026-08-17)
// - Handle 참조 무결성 + model-id 일치 (디버그)
// - Constructed 정점: 참조 곡선/곡면 위에 정확히 놓임 (eps_machine)
// - Discovered 정점: tol 이내
// - tess 출처 정합: OnFace{uv} 평가값과 pos 일치
// - tess 틈 없음: 인접 면 삼각분할이 공유 엣지 polyline 정점을 정확히 공유 (§5)
```

**validate는 store 전체가 아니라 live 도달가능 셀을 센다(§2 supersede 의미론).** Euler의 V·E·F·S·L_i는 `store.len()`이 아니라 `live_solids`에서 도달 가능한 정점·엣지·면·셸·내부루프 수이고, manifold(엣지 정확히 2회)도 "전역 2회"가 아니라 "도달가능 집합 내 2회"다. 이유: supersede된 옛 셀이 아레나에 남으므로, store 길이로 세면 죽은 셀이 Euler를 깨고 죽은 면의 엣지가 manifold를 깬다. **참조 무결성 검사(dangling handle)만 store 전체를 훑는다** — 안전 게이트이자, 살아있든 죽었든 handle이 OOB면 버그이기 때문. 이 검사가 먼저 단락하므로 이후 도달가능성 순회는 항상 in-bounds라 안전하다. (부수효과: store에 떠 있는 stray 셀은 이제 불량이 아니라 "죽은 아레나 항목"으로 정당하게 무시된다 — M1~M3의 "stray → Euler 불량" 판정은 이 의미론으로 갱신된다.)

속성 기반 테스트(proptest): 랜덤 유효 연산열을 생성해 (1) validate 통과, (2) replay 멱등성 — 같은 로그·같은 cfg → 동일 모델 **✔ 구현(`nacre-ops/tests/replay.rs`)**, 그리고 같은 자리에서 더 강한 「`replay(log)` == 세션 모델, 인덱스까지」까지 함께 건다, (3) 변환 불변량 — 강체변환 후 부피·면적 보존, (4) 불리언 대수 — `A ∪ A = A`, `A ∩ ∅ = ∅`, `vol(A∪B) + vol(A∩B) = vol(A) + vol(B)`.

OCCT 오라클(`nacre-oracle`): 같은 연산열을 자체 커널과 OCCT에 병렬 실행하고 부피·면적·바운딩박스·(가능하면) 면 개수를 diff. 허용 편차를 넘으면 실패한 연산열을 최소화(shrink)해서 리포트. 정답지를 든 채 개발하는 장치이며, AI가 생성한 코드의 "그럴듯하지만 틀림"을 잡는 주 방어선. **nacre 쪽 부피·면적은 `nacre-props`(정확 기하 발산정리, 해석적)가 계산하고, 오라클이 그 값을 OCCT가 같은 STEP에서 독립 계산한 값과 diff한다** — 완전히 별개인 두 구현의 일치로 양쪽을 교차검증(M4 큐브·실린더부터 가동).

**OCCT는 오라클 전용이며, 연동은 out-of-process다.** OCCT를 링크하지 않고, `tools/occt-helper/`의 헬퍼 프로세스와 STEP 파일로 주고받는다(수송 계층 = nacre-step). 이유: C++ 빌드가 Rust 워크스페이스에서 완전히 사라지고, OCCT가 악조건 입력에서 크래시해도 커널 프로세스가 아니라 헬퍼만 죽는다(크래시 격리 — 오라클처럼 수천 회 돌리는 용도엔 필수). 오라클은 채점(부피·면적·위상 수 비교)만 하므로 history 손실·호출 오버헤드가 무관하다.

**indirect predicates 검증(M5).** indirect 술어는 틀려도 대부분 입력에선 맞는 답이 나와 M5에서 버그가 가장 숨기 쉬운 조각이므로 오라클을 두 겹으로 건다. (1) **direct 술어 대조(property test)**: 좌표를 아는 케이스에서 relaxation으로 f64 좌표를 뽑아 `geometry-predicates`의 direct 술어(orient3d 등)에 넣은 부호와, 같은 점을 implicit point로 둔 indirect 술어의 부호가 일치하는지 대조. (2) **Attene C 구현 실행 대조(dev-only)**: 저자 참조 구현을 별개 프로그램으로 실행해 답만 비교(§8 M5 라이선스 규율 — `tools/` 격리, 링크 금지, OCCT 오라클과 동일 논리).

헬퍼 구현 우선순위(macOS 기준): **1순위 — Homebrew** (`brew install opencascade`): DRAWEXE Tcl 스크립트(readstep → bfuse/bcut/bcommon → writestep, 코드 0줄) 또는 brew 라이브러리에 링크하는 얇은 C++ 헬퍼. **폴백 — uv 관리 Python 환경 + OCP 휠**(사전 컴파일이라 C++ 툴체인 불요, uv가 파이썬 버전 고정까지 해결). 어느 구현이든 프로토콜은 동일하게 고정한다: `helper <props|fuse|cut|common> …`, exit code 0/1/2 = 성공/기하 실패/크래시, stdout으로 진단 값(부피·면적·면 개수·bbox — 헬퍼 쪽에서 계산해 전달). 값은 **간단한 key-value 라인**(`volume <v>` 등)으로 낸다 — 5필드 고정 스키마라 파싱 의존성 0(원안의 "JSON"을 이렇게 단순화; 안정성 동일). 프로토콜만 지키면 구현을 갈아타도 nacre 쪽 코드는 무변경. (v0는 `props <in.step>`만 — 불리언 전. `fuse|cut|common`은 M5.)

STEP 라운드트립(자기 출력 되읽기): nacre가 쓴 STEP을 step-io **리더로 되읽어** 구조 비교(외부 파일 import가 아니라 "방금 쓴 것" — healing 불요). `Intersection` 곡선의 surface_curve 이중 표현 보존을 중점 검증.

## 8. 마일스톤 사다리

각 단계는 "동작하는 것"을 남기고 끝난다. 3층에 막혀 전체가 멈추는 구조를 피하는 배치다.

**M1 — 뼈대.** nacre-math, Store/Handle, Plane/Line만으로 정육면체를 손으로 조립. validate 1차 구현. 시각 확인은 **기존 뷰어에 위임**한다(전용 뷰어 크레이트를 만들지 않는다): 정상 결과는 STEP 출력→step-loupe(구조+검증, §7 라운드트립), 중간·깨진 상태는 Tessellation을 OBJ/STL로 덤프→맥 미리보기(또는 MeshLab·f3d). 커널 내부까지 보는 인터랙티브 디버그 뷰어는 워크스페이스 밖 별도 앱으로 M5 즈음 만든다(§1). — M1에서 렌더링에 시간을 쓰지 않는다. (전부 1층. AI 가속 최대 구간.)

**M2 — 스케치와 케이스 A.** 2D 프로파일(선분만) → extrude. 연산 로그와 replay. 오일러 연산 정리. **STEP 내보내기(어댑터).** 최소 평면 b-rep 출력 — 커널 Model→AP242 엔티티 번역만 자체 구현(어댑터), 직렬화는 step-io 백엔드에 위임(커널은 백엔드 무지 — `BooleanEngine` trait·헬퍼 프로토콜과 같은 격리). 어댑터가 내보내는 엔티티를 AP242 커널 형상 집합(cartesian_point·direction·line·circle·b_spline_curve/surface·plane·cylindrical_surface·advanced_face·closed_shell·manifold_solid_brep·surface_curve 등)으로 명시 타입 고정 → 나중에 경량 라이터가 감당할 범위를 개발 내내 붙박음. 커버리지는 마일스톤별 확장(M2 평면 → M3 곡선·곡면 → 이후 트리밍). 목적은 기능일 뿐 아니라 검증: 정육면체를 STEP으로 내보내 뷰어로 열어 좌표계·방향·면 orientation 함정을 작은 형상에서 미리 밟고, `Curve::Intersection`↔surface_curve 대응 가정을 조기 확인. 검증 뷰어는 두 겹 — **step-loupe**(step-io 기반 웹 뷰어: report가 드롭·고아·비표준 엔티티를 표시해 어댑터 버그를 구조적으로 잡고 모양도 시각 확인; §7 "자기 출력 되읽기" 라운드트립의 GUI판) + **FreeCAD 등 독립 OCCT 기반 뷰어**(step-io를 공유하지 않는 교차검증). 둘은 상보적이다. 오라클 연계: nacre-oracle(§7, M4~)이 STEP을 OCCT 헬퍼로 보내는 수송 계층으로 쓰므로 STEP 출력이 M2에 있어야 M4 오라클이 즉시 돈다. 추가로 step-io 리더로 자기 출력을 되읽어 라운드트립 자체검증이 가능하고, 최종 경량 라이터로 교체할 때 "step-io 출력 vs 경량 출력"을 step-io 리더로 비교해 교체 안전성을 자동 검증할 수 있다.

**M3 — 곡선 기하.** Arc, Cylinder, NurbsCurve/Surface 평가(The NURBS Book 기준 구현 + 수치 미분 대조 테스트). tess 출처 태그 완성, tolerance 재계산 데모(같은 모델, tol 3단). proptest 도입.

**M4 — 면 위 작업.** ImprintSketch, PadOnFace — "만나는 자리를 아는" 연산의 완성. 여기까지 모델 전체가 `Constructed`. nacre-oracle 가동(nacre 쪽 부피·면적은 `nacre-props` 해석적 계산, OCCT와 diff). 시각 확인은 M1과 동일하게 기존 뷰어에 위임(STEP→step-loupe, OBJ→맥 미리보기). 커널 내부를 보는 인터랙티브 디버그 뷰어(면 클릭→Handle·Origin, 법선 화살표, 엣지 polyline·tolerance 공, validate 위반 하이라이트, 로그 스텝별 재생)는 워크스페이스 밖 별도 앱으로 **M5 즈음**(`Discovered`/tolerance가 처음 등장해 STEP에 안 담기는 내부 정보의 시각화가 실제로 필요해질 때) 만든다(§1) — 그것이 M6 불리언 디버깅의 생명줄이 된다.

**M5 — 자체 불리언 1단: 다면체.** `PolyhedralBoolean` — 모든 면이 평면인 솔리드 간 fuse/cut/common을 자체 구현한다. 평면-평면 교차는 닫힌 형식의 직선(SSI 행진·Newton·캐시 불필요). 꼭짓점은 평면 3장 교차를 **좌표로 만들어 병합하지 않고 implicit point로 두고**, 내/외·orientation 부호는 그 정의를 **indirect orient3d**(좌표 안 만듦)에 넣어 정확히 판정한다 — 좌표를 만드는 순간의 오차·불일치를 원천 차단(§3·§4; Attene 2020). 세 엣지가 하나의 Vertex Handle을 공유하게 해 봉합 문제의 본체를 우회한다(§4). 이 세계에서 강건 불리언은 연구가 아니라 꼼꼼한 케이스워크(공면, 엣지-엣지 퇴화)다. 하이브리드 파이프라인(출처태그 메시 → 조합 결정 → 스냅백)을 스냅백이 자명한 평면에서 첫 완성. 선·평면(다항식) 교차점은 indirect predicate 이론이 가장 깔끔하게 도는 영역이라 M5가 논문 대표 예시와 정확히 겹친다 — 딱 필요한 만큼 구현.

**indirect predicates 자체 구현(clean-room).** `nacre-predicates`(§1)에 implicit point 표현과 indirect 술어를 자체 구현하되, 확장 산술 바닥은 `geometry-predicates`(MIT/Apache) 재사용. **라이선스 엄수** — Attene 참조 구현(LGPL)은 "돌려서 답 비교는 자유, 열어서 코드 보는 건 MIT/Apache 소스만": (a) **작성 중 소스 열람 금지**(LGPL C++를 열어 함수 대응·로직 흐름을 따라가면 2차적 저작물), 참고처는 논문(Attene 2020, arXiv 2105.09772; Cherchi 2020 mesh arrangements·2022 interactive booleans; Lévy 2024; Shewchuk 1997)과 `geometry-predicates` 소스로 한정. (b) **완성 후 실행 대조는 허용·권장**(dev-only, `tools/` 격리, 우리 크레이트에 링크 금지 — OCCT 오라클과 동일 논리). indirect predicate는 틀려도 대부분 입력에선 맞는 답이 나와 버그가 숨기 쉬우므로 저자 구현을 정답지로 쓰는 게 강력한 검증(§7).

**exact-arithmetic 바닥 확정(M5-prep, 코드로 실증).** `geometry-predicates`(elrnv, 0.3, MIT/Apache)가 finished `orient3d`뿐 아니라 **Shewchuk expansion primitive**(`two_product`·`two_sum`·`expansion_sum`·`scale_expansion_zeroelim` 등, `predicates` 모듈에 공개)를 노출함을 `nacre-predicates` 뼈대가 실제 호출로 확인 → 그 위에 indirect 술어를 쌓을 수 있고 expansion 산술 자체 구현이 불필요. (primitive는 `[lo, hi]` 순서. orient3d 부호 규약 = `det[a−d, b−d, c−d]`, 뼈대 golden으로 고정.) **M5 서브유닛 사다리:** ① `nacre-predicates`(뼈대→implicit point[3-plane]+indirect orient3d, direct 대조 property test) → ② `nacre-geom::intersect`(평면∩평면=닫힌형식 직선, 3-평면 꼭짓점) → ③ `PolyhedralBoolean`(fuse/cut/common, 세 엣지가 한 Vertex Handle 공유로 봉합 우회, `Origin::Discovered.definition` 필드 도입, 커버리지 밖 `Rejected` 거절) + 오라클 부피·불리언대수 proptest.

알고리즘 참고: Manifold(Apache-2.0 — 차용·번역 가능), Hoffmann 등 문헌. 커버리지 밖 곡면 불리언은 명시적 미지원 에러로 정직하게 거절. `Discovered` 경로·국소 tolerance·relaxation 실전 투입. 이 시점에 "OCCT 없이 직동하는, 실용적 기계 부품(평면 위주)을 STEP으로 내보내는" 진짜 커널이 된다. 참고: Truck 대비 벤치마크·정밀도 비교(전역 1e-6 폴리라인 vs 정점별 실측 tol + 닫힌 형식)는 수치로 보여줄 수 있는 차별점 — 공개 지표 후보.

**M6 — 자체 불리언 2단: 이차곡면.** 평면∩실린더(타원), 평면∩구(원), 평면∩원뿔 — 여전히 닫힌 형식이라 행진 불필요. 실린더∩실린더는 특수 케이스(직교 등)부터. 실제 기계 부품 면의 대다수가 평면+실린더+원뿔이므로, 여기까지로 실용 커버리지의 대부분을 확보한다. ★ **정정(칸 ⑧, 2026-09-05)**: 원통 쌍 가운데 **면이 만나지 않는** 쌍(십자로 관통한 두 스터드처럼 축은 만나도 남은 면은 비킨 것)은 새 기하가 아니라 **게이트의 명제** 문제였다 — 옆면의 축 구간으로 비켰음을 증명하면 배열이 이미 조립한다. 여기 남는 것은 면이 **만나는** 쌍(4차 곡선)뿐이다.

**M6/M7 난이도 절벽.** M6(평면+이차곡면)는 닫힌 형식 교차라 SSI의 지옥(위상 판정·시작점 검출·watertightness)을 **대부분 회피** — 여기까지가 "확실히 되는 실용 커널". M7(일반 NURBS)은 그 60년 미해결 문제(아래)로 **직접 진입**. 그래서 "**M6까지가 실용적 종착점, M7은 별도의 장기 연구 트랙**"이라는 선을 분명히 긋는다. 실용 가치는 M5~M6에 있고, M7이 안 돼도 M6까지로 실제 기계 부품 대부분을 커버한다.

**M7 — 자체 불리언 3단: 일반 SSI (연구 구간).** nacre-geom::intersect에 SSI 행진 구현 → 일반 곡면쌍의 `HybridBoolean` 완성: 출처태그 tess에 강건 메시 불리언(exact predicates) → 조합 결정 추출 → 살아남은 면은 정확 곡면 유지, 신규 엣지는 국소 SSI 스냅백. OCCT 오라클과 상시 diff. 실패 케이스 코퍼스 축적.

**M7 내/외 분류 = exact ray casting (1순위 후보).** 하이브리드 파이프라인("출처태그 메시 → 조합 결정 → 정확 곡면 스냅백")에서 "이 patch가 최종 솔리드 안인가 밖인가"를 분류하는 단계에 exact ray casting을 쓴다. 메커니즘: 레이가 삼각형 **내부**를 지나면 삼각형 방향(정점 순서)으로 정확 판정, 꼭짓점·엣지·접선(coplanar) 같은 애매한 케이스는 레이를 **수치 섭동**해 항상 "내부 통과"로 되돌린다(Simulation of Simplicity, Edelsbrunner–Mücke 1990). 1순위 이유: (a) **우리 계보와 정합** — 핵심이 orient 술어(채택)+섭동(`RelaxError::Tangential`로 올리기로 한 퇴화 구역의 정석)이라 새 수학을 안 들인다; (b) **검증됨** — Cherchi 2022(Interactive and Robust Mesh Booleans)가 핵심으로 채택, 수백만 삼각형·수백 입력 variadic까지 테스트하며 GWN류를 명시적으로 제침; (c) **파이프라인에 그대로 꽂힘** — 우리 하이브리드는 이미 메시를 경유하므로 분류를 메시 단계에서 함(Cherchi 검증 형태 그대로), 곡면 직접 판정 불필요라 GWN의 trimmed NURBS 확장(최신·검증 진행 중)을 우회.

**오해 방지 셋.** (1) **Truck 방식이 아니다** — 표면적으로 메시를 쓰나 정반대다: Truck은 교차를 폴리라인 근사로 표현해 그 근사가 **최종 결과**(전역 1e-6), 우리는 메시를 **분류용 임시 도구로만** 쓰고 exact 술어로 분류 후 **정확 곡면으로 스냅백**. 차이는 "메시를 쓰느냐"가 아니라 "메시가 최종이냐(Truck) vs 임시냐(우리)". (2) **indirect predicates의 대체·확장이 아니다** — indirect predicates=M5 **점** 부호 판정(평면·국소·교차 계산 중), exact ray casting=M7 **patch** 내/외 분류(전역·분류 단계). 대상(점 vs 덩어리)도 마일스톤(M5 vs M7)도 다른 별개 부품. (3) **"어려운 케이스에만 켜는 정밀 모드"가 아니다** — relaxation 사다리·적응 술어는 케이스 단위로 "어려우면 더 정밀하게", ray casting은 M7에서 **상시 도는 분류 단계**. "선택적"의 단위는 케이스가 아니라 **마일스톤**(M5·M6엔 닫힌 형식이라 불필요, M7에서 켜짐).

**한계 — ray casting은 분류를 풀지 SSI를 풀지 않는다.** M7 "내/외 분류" 절반을 검증된 방법으로 채워도, 진짜 도박인 **SSI(곡면 교차 곡선 계산)**는 미해결로 남는다. ray casting은 "조각을 어떻게 분류하나"만 정확히 할 뿐 "곡면 교차 곡선을 어떻게 정확히 뽑나"를 풀지 않는다. M7이 도박·연구 구간인 것은 그대로고, ray casting은 그 도박의 한 부품(분류)만 안정화한다. **덧붙임 — "스냅백"과 "SSI"는 다른 일이다.** SSI(교차 곡선을 처음부터 **찾는** 것)가 도박인 미해결 문제이고, 스냅백(이미 찾은 교차 근처에서 정확 곡면으로 **되당기는** 것)은 그다음의 별개 작업이다. 즉 M7의 위험은 "스냅백"이 아니라 "SSI로 교차를 정확히 찾을 수 있는가"에 있다. 스냅백은 SSI가 성공한 뒤라야 의미가 있으므로, SSI가 실패하면 스냅백까지 갈 것도 없이 그 연산은 `Rejected`로 거절된다.

**SSI는 검증된 해법이 없는 60년 미해결 문제.** 일반 곡면 교차는 1960년대부터 연구됐으나 CAD가 받아들일 만큼 강건·신뢰할 해가 여전히 없다 — **Parasolid·SISL·IRIT 등 최고 커널조차 특정 위상 케이스(작은 loop, 접선, cusp, 미세 자기교차)에서 실패**하고, 2025년 논문들도 watertightness 보장이 "여전히 challenging"이라 적는다. M7 "도박" 표기는 과장이 아니라 정확한 현실 인식이다. 아래는 **"채택할 검증된 기법"이 아니라 "M7 진입 시 읽을 프론티어 후보"**로만 걸어둔다(저자 벤치마크 주장일 뿐 독립 검증·프로덕션 채택 전 — 지금 상세 검토는 "정보 없이 미리 설계" 함정이라 M7 진입 시로 미룸): winding number + subdivision 시작점 검출(작은 loop·접선 branch 놓침 완화, 2026), Dixon matrix tracing(branch jumping/missing을 root-solving으로, Chen 2025), interval algebraic topology analysis(SSI 위상을 4D 대수계로 분류, Cheng 2023), lower-dimensional formulation gap control(Wang 2025).

M7은 열린 연구임을 명시한다. M6까지가 "확실히 되는" 영역, M7은 이 커널의 존재 이유이자 도박이다. OCCT의 역할은 전 구간에서 **오라클(dev 전용 시험관)뿐**이다 — 제품 경로에 OCCT 위임은 없다. OCCT 소스는 "상용급이 이 케이스를 어떻게 다루나" 열람용 참고서로만 쓴다: LGPL-2.1이므로 번역·차용은 라이선스 오염이고, 무엇보다 OCCT의 위상·tolerance 아키텍처가 딸려 들어와 nacre 설계와 충돌한다. 읽되 베끼지 않는다.

## 9. 미결 사항 (다음 논의 대상)

트리밍 곡면의 pcurve 표현 시점(M3에 선행 도입 vs M5까지 지연), 닫힌 엣지의 seam 처리(방식 확정 — 아래 "원통 seam: A vs B" 항목), Sketch 제약 솔버의 범위(초기엔 무제약 프로파일만), OpRef 계보 참조의 도입 시점과 직렬화 포맷 여유분, `Store` 스냅샷·직렬화 포맷(자체 vs STEP 재활용) — 이와 함께 **세션 중 메모리 관리: compact보다 재구축(rebuild-from-log) 우선**(§2 "재검토 예정 (v2)" 참조; 재구축=주력 정리·undo 유지, live 폐포 필터링=저장, compact=비상 회수), OCCT history → 출처 매핑의 실제 충실도(M5에서 실측 필요), 멀티스레딩 경계(Store가 &mut 독점인 설계라 연산 단위 병렬은 미지원 — 의도적 단순화).

**원통 seam: A(seam 엣지) vs B(seamless periodic) — A 채택 (M3.3a에서 확정).** 주기 곡면(원통·구·토러스)의 옆면을 b-rep로 담는 두 방식이 있고, 둘 다 유효 AP242·유효 CW-복합체다(초기 판단에서 "B는 오일러가 깨진다"고 봤으나 **오류** — 깨지는 건 정점조차 없는 제3의 변형 C[`bounds:None`, V0]이고, B는 각 원을 seam 정점 `Some([v,v])`로 두고 옆면을 두 루프[outer=아래원, inner=위원]로 담아 `V−E+F−L_i = 2−2+3−1 = 2`로 통과한다).
- **A(채택):** 위 원 + 아래 원 + **세로 seam 직선 엣지**, 옆면 = 4-엣지 단일 닫힌 루프 `[bottom, seam, top⁻, seam⁻]`(seam이 같은 면에서 2회 반대 — self-adjacent). 원통 생성자(`cylinder_solid`)가 이 방식. STEP 출력도 A(OCCT 계열이 생산·기대하는 형태; step-io 검증 테스트와 동형).
- **B(미채택):** seam 엣지 없이 위·아래 원 두 개로만 옆면 경계(옆면이 두 루프). NIST 샘플 계열.

**결정 근거(정직한 저울질):**
1. **불리언(이 커널의 존재 이유) — 결정적.** 조합적 b-rep 불리언은 면 루프를 순회하며 교차곡선에서 엣지를 쪼갠다. B의 주기 곡면은 교차곡선이 u=0/2π를 가로지를 때 쪼갤 엣지가 없어 특수 처리가 필요하다 — 조합 알고리즘이 "모든 면 = 실제 엣지로 둘러싸인 유계 영역"을 요구하기 때문. ★ **정정(2026-08-17 재조사)**: 원문의 *"OCCT·Parasolid·ACIS가 전부 내부적으로 seam을 넣는다"* 는 과장이었다 — 실측: **OCCT만 seam 필수**(없으면 데이터 모델상 invalid), **ACIS는 `periodic_no_seam` 옵션**(seam 없는 두-루프 원통 면 허용), **Parasolid는 winding loop**(주기 매개변수를 감는 루프)로 seamless가 네이티브다. 즉 업계 관행은 근거가 못 되고, **이 선택의 근거는 아래 2~4의 자기 제약 + 5(유리수 차트)로 자립한다**. B는 dumb solid의 교환 표현으로 흔하고, nacre는 연산하는 커널(M6에서 원통이 불리언에 진입)이라 A.
2. **tess 정합.** nacre tess는 "공유 엣지 polyline(`by_edge`)을 경계로 소비"해 crack-free(§5). A는 seam이 진짜 엣지라 그 polyline을 옆면이 u=0·u=2π 양쪽에서 소비 → **기존 메커니즘 그대로, 주기 특수처리 없음**. B는 "이 면은 periodic이니 u-wrap 봉합" 특수처리를 tess에 새로 요구. 즉 **우리 설계에선 A가 tess도 더 단순**("B=tess 간결"은 tess가 주기 곡면 네이티브일 때만 참).
3. **"가짜 모서리" 비용은 흡수됨.** seam은 각진 모서리가 아니라 파라미터 이음매(양쪽이 같은 곡면이라 C¹ 매끄러움)다. 하지만 매끄러운 엣지는 seam 말고도 존재(접선-연속 fillet 등)하므로 nacre는 "sharp vs smooth 엣지" 판정을 **어차피** 갖는다 — 판정식 = 인접 두 면의 곡면 법선이 엣지를 따라 일치하는가. seam은 양쪽 같은 실린더라 자동으로 "smooth"로 판정된다. **seam 전용 플래그 불필요; 일반 술어에 흡수.** 따라서 A가 만든 "엣지 예외"는 새 범주가 아니고, B가 만드는 "면 예외(주기 경계)"는 불리언에서 훨씬 비싸다.
4. **위상 균일.** 모든 면이 실제 엣지의 닫힌 루프로 둘러싸임 → validate·루프순회·불리언이 단일 규칙. self-adjacent seam은 validate가 이미 무수정 통과(엣지 2회·반대만 보고 면 구별 안 함).

5. **★ 유리수 차트와의 시너지(2026-08-17 재조사 — M6 판정 설계의 발견).** M6a의 판정은 원통 위 점을 단위원의 유리수 반각-매개변수 t로 드는 방향이 유력한데, 그 차트는 원에서 **정확히 한 점을 못 덮는다**(t=∞). seam을 그 배제점에 두면 **면 내부 전체가 유한한 유리수 t로 덮인다** — seam이 인공물이 아니라 유리수 차트의 자연 경계가 된다. B(seamless)였다면 그 특이점이 면 한복판에 남고, 주기 방향에 전순서가 없어 정확-조합 순서 비교도 무너진다.

**어댑터 분리(옵션 보존):** 내부를 A로 두어도 nacre-step은 어댑터라 **필요 시 export 시점에 unseam해 B로 내보낼 수 있다**(내부 표현 ≠ 교환 표현). 지금은 A 출력(OCCT 상호운용 안전)이고, 특정 상호운용 요구가 생기면 B 출력을 어댑터에 추가. **tess 층의 남은 세부**(seam polyline 샘플 밀도·법선 비분리 규칙 구현)는 M3 tess 곡면 샘플링에서 확정.

**★ M6-0 규약 (2026-08-17 확정 — 원통의 진실).** M6a(원통-우선) 진입 재조사에서 확정한 여섯:
1. **진실 형태**: `SurfaceTruth::Cylinder { def: CylinderDef, motion }`, `CylinderDef { origin, dir, ref_dir, radius }` 전부 유리수 — dir·ref_dir은 **정규화하지 않은 원시**(normalize가 정확형을 파괴 — `normal_def` 선례), 캐시는 실현. 성분형이 옳은 이유: 사용자 어휘의 직접 리프트 + 곱-없는 셔플 — 반증표가 금지한 «유도된 곱(계수)» 모양이 아니라 평면 *점*의 원통 대응물. ref_dir 원시는 `any_perpendicular` **자신의 규칙을 유리수로**(최소-|성분| 축, 동률 X→Y→Z). ★ **축 선택은 정규화된 `d` — 캐시가 실제로 읽는 값 — 를 읽는다**(구조적 일치): 첫 철자는 원시 성분을 읽고 "양수 배는 순서 보존"을 논증했으나 그건 실수 산술 논증이었다 — f64 나눗셈 반올림이 강부등호를 동률로 붕괴시켜(실측 축 `[0.34, 0.33999999999999997, 1]`: 원시는 Y, `d`는 X) seam이 캐시와 ~90° 어긋난다. 외적은 여전히 원시 정확 성분으로(어느 값이 기저를 골랐든 `ê_k×원시 ∥ ê_k×d` 양의 평행 — seam 방향 정확 보존).
2. **seam은 모델 기하다 — 여기서 영구 고정된다**(+ref_dir 방향, θ=0). M6-1이 정하는 것은 **차트**뿐이고, 차트가 자기 배제점을 기존 seam 위에 놓도록 맞춘다 — 차트가 seam에 적응하지, 그 역은 절대 아니다.
3. **interning은 보수적으로**: `cylinder_ids`(별도 맵), def **문자 동일**(ref_dir 포함) + motion 만 합침 — 같은 축·반지름·다른 ref_dir을 합치면 seam이 갈라지므로 잘못 합칠 위험 0인 키로 시작, 기하 동일성은 규칙 6대로 술어 몫(M6-1). 평면의 `flipped` 대응물 불요(문자-동일 키면 캐시 구성도 동일).
4. **OnSeam 정점의 정의 완성**: OnSeam([원통, 캡]) = "rim ∩ +ref_dir 방향 ray" — ref_dir이 진실에 앉아 유일점을 정확히 지시한다(좌표 캐시의 load-bearing 해제; 좌표-재생 기계는 별도 유예).
5. **seam 담체 `[s,s]` 확정** — "한 면의 매개화 이음매"의 정직한 철자(잠정 딱지 제거), validate «자기-인접 ⇔ 원통» 규칙이 지킨다.
6. **branch 방향**: 평면∩평면∩원통 = 최대 2점 → `VertexDef` **새 변종**으로 받는다(Q3 정신 — 변종이 자기 진실을 말한다; 슬롯 재활용 금지). 구현은 M6-1. ★ **정정(2026-08-21)**: 그 «최대»가 근의 어휘에도 적용된다 — `QuadRoot` 는 `Lo|Hi` 둘이 아니라 **`Lo|Hi|Double` 셋**이다. 접점을 `Lo` 로 적는 M6-1 의 규약은 `Lo` 가 «둘 중 작은 쪽» 인지 «유일한 쪽» 인지 말하지 못하게 만들었고, 그래서 정의를 든 쪽이 재정렬 때 토글해야 하는지 알 수 없었다. 재정렬 규칙은 `QuadRoot::canonical` **한 곳**에 살고, 접점 예외는 `flipped(Double) = Double` 이라는 원시에서 저절로 나온다(스왑은 중근을 자기 자신으로 보낸다).

**패드/포켓 ↔ 불리언 통합 — feature = tool body + boolean (✅ 완료).** `pad`/`pocket`은 상용 CAD와 동형으로 이미 통합됐다: **`pad` = 프로파일 압출 → `Fuse`, `pocket` = 압출 → `Cut`**(`extrude_and_boolean`), 패드/포켓은 불리언의 얇은 sugar다. `raise_region`(M4 직접 구성)은 폐기됐고, "프로파일이 면 안에 있어야 한다"는 컨테인먼트 제약도 제거돼(`extrude_and_boolean`이 미검사) 오버행 패드/포켓이 불리언의 공면-접촉·오버행 경로(`detect_contained_contact`·`detect_pocket_contact`·`detect_overhang_contact`)로 자동 처리된다. 공면 복잡도는 **불리언 하나로 집약**됐다. **남은 격차(후속)**: 오버행 `Cut`/`Common`의 비볼록 kept-solid는 아직 미검증(재보지 않았다). ★ **2026-07-22 갱신 — 셋 중 둘은 닫혔다:** 비볼록 오버행 footprint(옛 볼록 게이트)와 정확 flush-edge(프로파일 테두리 = 면 테두리 공유)를 **한 형상이 동시에** 통과한다 — `a_non_convex_pad_cantilevers_and_runs_flush`(부피 1.625·면적 9.75, 손계산과 일치). corner-flush `Common`(공면 3면, 옛 `vertex_on_face_plane` 거절)도 열렸고 **OCCT와 부피·면적이 일치**한다(`a_corner_flush_common_keeps_the_non_convex_overlap` + `corner_flush_common_matches_occt`). **잔존 직접-구성은 없다 (2026-07-22 은퇴).** `ImprintSketch`(면을 재료 가감 없이 분할하던 별개 op)가 마지막 직접-구성이었는데, `raise_region` 폐기 이후 **소비자가 하나도 없었다**(`region_face`를 읽는 코드 0). 연산과 `imprint`/`prepare_face_split`/`finish_split`/`placed_profile`(엄격 컨테인먼트)·`OpError::ProfileNotContainedInFace`를 제거했다. **그것이 커널 유일의 동일평면 인접면 생산자였으므로**(불리언 출력은 `unify_coplanar_faces`가 항상 정리한다) 불리언의 이음선 처리도 함께 제거했다 — 되살릴 근거(유도·측정)는 dev-log에 남아 있고, 복원은 그 두 커밋을 revert하는 것이 출발점이다. Split Face가 필요해지는 소비자(구역별 재질·FEA 경계조건·금형 파팅라인)는 아직 로드맵에 없다.

- **★ 구조적 공면 = O(1) 참조 인식(핵심).** 면 위 스케치를 압출한 tool body의 밑면은 대상 면의 **surface Handle을 공유**한다(`build_prism`이 대상 면의 `surface_h`를 밑면에 그대로 쓴다). 통합 불리언은 **먼저 Handle 공유를 검사**(float 계산 0)해 공유 면을 seam으로 즉시 채택하고, **공유하지 않는 독립 솔리드만** `planes_coplanar` 기하 감지로 내려간다. 즉 "참조로 아는 공면은 공짜, 우연한 공면만 계산" — 상용 커널이 "coincident 면을 imprint로 공유 토폴로지로 승격"하는 것의 nacre판.
- **★ 게이트(순서 강제).** (a) 불리언 공면 처리를 하나로 흡수·완성(지금 셀 사다리 — Cut 3갈래 통합, Common이 Cut의 detect 재사용) → (b) Handle-공유 O(1) fast path 추가 → (c) 그제서야 pad/pocket을 extrude+불리언 wrapper로 바꾸고 M4 직접 경로 제거. **역순 금지**: 불리언 공면이 robust해지기 전에 M4를 걷어내면 지금 되던 포켓/보스가 커버리지 구멍에 빠진다(M4 직접 경로는 그때까지 신뢰 가능한 fallback).
  - **(c-1) pocket 완료 (실행됨).** `PocketOnFace`는 `build_prism`(inward, top-flush) + `boolean(Cut)`으로 재구현됐다 — contained-coplanar Cut 경로가 빈 seam을 내므로 결과는 전부 Constructed(tolerance 0), 옛 직접 `raise_region`과 등가. 선결로 `detect_pocket_contact`의 볼록성 게이트를 제거(비볼록 kept `a`·비볼록 프로파일 `b` 둘 다 열림, OCCT로 교차검증)했다. `raise_region`은 `pad`용으로 잔존. through-pocket(`dist ≥ 두께`)은 이제 정직히 거절(M4는 미검사 UB였음 — 개선). **게이트 (b) Handle-공유 fast path는 불필요로 폐기**: contained 경로가 이미 빈 seam·전부 Constructed라 참조-fast-path가 더할 exactness가 없다. 다음: pad(=extrude+Fuse), 그다음 M4 직접 경로 완전 제거.
  - **(c-2) pad 완료 (실행됨).** `PadOnFace`도 `build_prism`(outward, top-flush) + `boolean(Fuse)`로 재구현 — pocket과 대칭. 선결로 `detect_contained_contact`의 볼록성 게이트를 제거(비볼록 kept base·비볼록 프로파일 boss 둘 다 열림, OCCT 2종 교차검증). pad·pocket은 이제 부호(±dist)·`BoolKind`·복원-면-부재 에러만 다른 **공통 헬퍼 `extrude_and_boolean` 위의 얇은 wrapper**로 통일됐다("feature = tool body + boolean"의 코드화). `raise_region`은 유일 사용자였던 pad가 떠나며 **삭제**(그 죽은 `Split` 필드도 함께 정리). 남은 M4 직접 기계는 imprint(`imprint`/`prepare_face_split`)뿐 — 완전 제거는 imprint 정리 후 별도 스텝.
  - **(c-3) 오버행 pad/pocket 획득 (실행됨 — 로드맵 payoff).** `placed_profile`을 `placed_profile_unchecked`(CCW+배치)와 strict-containment wrapper로 쪼개고, `extrude_and_boolean`이 unchecked를 쓰게 했다. 이로써 **면 경계를 넘는 프로파일이 기존 오버행 불리언 사이드카로 자동 라우팅**(Fuse: 단일 엣지·코너·spanning slab; Cut: N-wall blind)된다 — "프로파일이 면 안에 있어야 한다"는 제약이 사라졌다. contained 경로는 완전 불변(같은 base_pts). 커버리지 밖(비볼록 오버행·through·far·비축정렬)은 **정직 거절**(`Boolean(Rejected)`) — silent-wrong 아님. 오버행은 contained(빈 seam)와 달리 footprint crossing에서 **Discovered seam 정점**을 만든다(기존 오버행 셀의 성질, 새 tolerance 아님). imprint는 containment **유지**(유효 inner-loop 필요). **알려진 비대칭**: contained pad/pocket은 비볼록을 받으나 오버행은 볼록 게이트로 볼록만 — 비볼록 오버행(오버행 detect 볼록 게이트 제거)이 다음 후보.
  - **(c-4) 오버행 boss가 비볼록 솔리드 수용 (실행됨 — c-3 비대칭 절반 해소).** `detect_overhang_contact`(Fuse)의 whole-solid `is_convex` 게이트를 **접촉면 footprint 게이트**(두 접촉면 outer loop 볼록 + 홀 없음)로 교체. 근거: Fuse 재구성은 **로컬** — arc-split은 (볼록) 접촉면만, `resplit_overhang`은 wall을 edge-local(볼록 무관)로, 나머지 면은 verbatim 재방출. 그래서 **비볼록 솔리드(포켓 파인 부품·부울 결과)에 접촉면만 볼록이면 boss 캔틸레버**가 붙는다. n0로 OCCT-정확 확인(포켓 큐브 옆면 오버행 = 1.17). silent-wrong 원천인 비볼록 **접촉 footprint**(arc 오분류)는 게이트가 거절. **Cut/Common은 whole-solid 게이트 유지** — 그 `clip_bwall_inside_a`가 breached 반평면 SH 클립이라 볼록 kept에서만 "inside a"와 일치(비볼록은 `OVERHANG_ARCS` 정직 거절, n0 실증). 신규 `loop_is_convex_2d`·`face_outer_is_convex`. **후속**: ② 비볼록 footprint 오버행(multi-piece 재구성), 비볼록 솔리드 오버행 Cut(clip_bwall_inside_a 일반화).
- **보존 불변식 — 의도는 유효, 메커니즘 서술은 arrangement가 대체함(커토버에서 갱신).** 원문은 *"(b)의 참조-fast-path가 순수성을 유지해야 한다 — 공유 면 seam은 Discovered로 승격하지 않고 Constructed로 남긴다"* 였으나, **(b)는 이미 폐기됐고**(위 (c-1)), 커토버 뒤의 arrangement는 결과 정점을 **평면 삼중항으로 다시 이름 붙이므로 회전 없는 부울 출력이 전부 `Discovered{ThreePlane}`** 다(`Node::Orig`는 생성되지 않음; 잠금 = `an_unrotated_boolean_names_every_vertex_by_its_plane_triple`). ⇒ **불변식의 의도(*"통합이 exactness를 후퇴시키면 안 된다"*)는 지켜진다** — 측정한 축정렬 사례에서 그 정점들의 `tol = 0`으로 `EPS_CONSTRUCTED`(1e-9)보다 오히려 엄격하고, 상자 모서리는 실제로 세 평면의 교점이라 정의도 참이다. **바뀐 것은 exactness가 아니라 표식이다:** `Origin`은 더 이상 *"어느 정점이 컷에서 왔나"* 를 구별해 주지 않으므로 **provenance로 면·정점을 고르는 코드는 성립하지 않는다**(이 사실을 모른 픽스처가 조용히 엉뚱한 면을 골랐다 — dev-log 참조).
- **업계 정합.** 피처 레이어는 pad=tool body+boolean으로 통합돼 있고(SolidWorks/NX/Creo/Fusion), 커널 레이어는 그래도 coincidence 전담 로직을 보유한다 — 상용은 그것을 imprint+tolerance로 두지만, nacre는 **커널 불리언 안**에 두었다(surface Handle 공유 + exact 술어). nacre는 tolerance 대신 exact 술어로 그 자리를 채우는 소수파.

**완전한 M5 오버홀 (회전 지원) — 확정 설계.** 축정렬 전용의 근본 한계를 넘어 회전(비정렬·마름모·각도 스케치)을 지원하는 실사용 M5의 오버홀. **확정 결정:** ① **유리수 치수·각도 표현**(§1 `nacre-scalar`; 유리수는 입력·유리수-순수 파생에만 살고, 무리수·복잡 연산·비트 상한이 닿는 순간 캐시를 f64/dd로 강등 — 정의는 불변이라 export 시 재계산으로 exact 복원). ② **명시 공유(전역 자동병합 폐기 — TNP 충돌):** 면에 스케치=`Surface` Handle 재사용·불리언이 계산한 접촉만 공유; 우연히 같은 좌표는 별개 유지(불리언 시점에 판정). ③ **CIP 판정 정책 — 묻지 않고 근거를 붙여 보고한다(2026-07-28 개정).** 원안은 *"애매하면 사용자에게 확인(site별·기본 `Reject`)"* 이었고 **실측이 반박했다**: (a) 애매한 판정이 연산당 수백 건이고, (b) 그 "애매함"의 실제 크기가 **2⁻²²⁹**(좌표 크기 1인 모델에서 소수점 69자리)이라 사람이 판단할 대상이 아니며 — 사람이 아는 것은 *의도*이고 의도는 이미 스크립트에 있다(같은 회전 공유 → 상쇄로 정확히 답함, `Surface` 핸들 공유 → 핸들 비교) —, (c) **"떨어짐"이라는 답은 만들 수조차 없다**(출력 좌표가 f64라 1e-16 아래 두께는 STEP·메시·화면 어디에도 안 나온다). 선택지가 하나뿐인 질문은 질문이 아니다. ⇒ 판정은 **일치 정밀도**(아래 ⑥)보다 가깝다고 *증명*되면 일치로 처리하고, 그 근거(무엇이 무엇과 몇 이내인지)를 결과와 함께 **보고**한다 — `boolean_with_report`(병행 진입점; `boolean`은 얇은 래퍼라 기존 호출부 무변경). **수집 채널은 전역도 thread-local도 아니다**: 연산 단위 컨텍스트 `Judge{planes, standard, notes}`가 소유하므로(증인 표는 순수한 *설명*으로 남고, 술어는 그 컨텍스트의 메서드다) 순수 수치 크레이트에 가변 전역이 생기지 않는다 — **연산의 성질은 연산이 갖는다**(행마다 복제하면 자리표시자·2단계 생성·"찍는 걸 잊었을 때" 가드가 딸려온다)(수집은 진단 전용 — 어떤 부호·병합·좌표도 이걸 읽지 않는다). 보고의 **첫 항목은 평면 클래스 병합**이다(정점 하나 계산되기 전에 "어떤 평면이 존재하는가"를 바꾸므로). 그리고 **증명하지 못한 판정은 0이 아니라 거절이다** — `JudgeExhausted`(비트 부족)·`DegenerateWitness`(여인수가 사라짐)로 원인별로 나가고, **원인이 증상보다 앞선다**(트레이서가 먼저 `LoopOrientMismatch`를 울려도 근거를 먼저 본다) — 사용자가 의도치 않은 일치를 알아차리고 **설계 치수를 고치게** 하는 진단이지, 커널의 결정을 대신하는 프롬프트가 아니다. "같다"면 Handle 재사용으로 추이성 붕괴를 구조적 차단하는 것은 그대로다. ⑥ **정밀도는 모델이 정한다.** 판정의 오차 반경은 `C·2⁻ᵖʳᵉᶜ`이고 `C`는 **모델의 성질**(회전 이력 1회당 약 1비트, 좌표 크기)이지 정밀도의 함수가 아니다 — 그래서 고정 정밀도는 *모델의 회전 이력이 얼마나 길 수 있는지를 조용히 결정*한다(실측: 256비트에서 245회 회전한 솔리드가 빌드 실패). `C`를 연산마다 한 번 읽어 필요 비트를 계산하고 워드 단위로 올린다. 두 값은 길이 단위로 분리된다: **`output_precision`**(출력 좌표의 정밀도, 기본 `모델 크기 × 2⁻⁵²` = f64가 이 모델에서 볼 수 있는 가장 가는 눈금)과 **`coincidence_precision`**(이보다 가깝다고 증명되어야 일치, 기본 그보다 워드 2개 아래). 비트로 노출하면 안 된다 — "256비트"가 회전 1회 모델에선 1e-76, 300회 모델에선 1e+15를 뜻한다. **`tolerance`가 아니다**: 전역 tol은 *"이보다 가까우면 붙여라"*(모르는 채 뭉갬), 이것은 *"이보다 가깝다고 **증명되어야** 일치"* 로 방향이 반대다. **판정 결과는 넷으로 갈린다**(2026-07-28 구현): `Sign`(부호 또는 정확한 0 — **증명됨**) · `Coincident { within }`(일치 정밀도보다 가깝다고 증명됨, 근거를 실어 나름) · `Exhausted { at, within }`(상한 `at` 비트에서도 못 가름 — 비트 부족. **자기가 세운 상한을 싣는다**: 참값은 `±within` 안이고 그 크기가 곧 소식의 무게다 — `1e-30`에서 멈춘 판정과 `1e-3`에서 멈춘 판정은 같은 variant이되 **같은 소식이 아니다**. `within`은 `Option`이다: 여인수가 상한에서도 안 풀리면 인용할 거리 자체가 없고, 그때는 없는 채로 나간다) · `Degenerate`(여인수를 0에서 뗄 수 없음 — 비트로 해결되지 않음). 마지막 둘은 **원인이 다르므로 이름도 달라야 한다**(한 이름으로 뭉개면 `LoopOrientMismatch`가 정밀도 고갈을 가리던 함정의 반복). 부족하면 `log₂(간격/목표)`가 곧 모자란 비트 수이므로 **계산해서 한 번에 점프**한다(배증 아님). 상한 `JUDGE_PREC_CAP`은 정확성이 아니라 **비용** 한계이며 값의 근거는 측정이다(회전 3200회·3456비트에서도 부피는 정확하고 90초; 4096비트 ≈ 3900회 ≈ 2분대) — 초과는 `RejectReason::PrecisionBudget`으로 **평면 표가 서자마자** 거절한다. **모델의 깊이와 한 판정의 난이도는 별개 예산이다**: 판정이 오를 수 있는 여유(`CLIMB_HEADROOM`)는 모델 정밀도에 대한 **상대값**이어야 한다(공유하면 회전을 많이 한 모델에서 여유가 0이 되어, 같은 얇은 증인이 회전 여부에 따라 판정되거나 포기된다). 그 여유의 물리적 의미는 **증인의 얇기 한계** `log₂(1/여인수)`다. **설정으로 노출할 후보는 일치 정밀도 하나뿐이고** 정밀도·상한·여유는 전부 거기서 유도된다. ④ **op-log 소유:** 상위 `Document`가 `Vec<Operation>`(진실)+파생 `Model`을 소유하고 sugar(`PadOnFace` 등)를 그대로 기록(피처 트리)·TNP는 "치수 변경=자동 replay, 위상 변화=수동 재지정". ⑤ **다중 솔리드 출력**(OCCT n0 확정: sever는 `Vec<Solid>`·각 manifold). 핵심 메커니즘은 격리 2D/3D 실험으로 선검증 후 이식했다(stage 1~3 대부분 완료 — dev-log).

**불리언 다중 솔리드 출력 + non-manifold 정책 (OCCT n0로 확정).** 실측: 엣지로만 만나는 두 박스의 `Fuse`, 막대를 가르는 `Cut` — OCCT는 둘 다 **COMPOUND of 2 SOLIDs**(checkshape valid, STEP `2× MANIFOLD_SOLID_BREP`, non-manifold 엔티티 0)로 낸다. 엣지 접촉 Fuse는 접촉 엣지+2정점을 **공유**(V16→14·E24→23)하되 각 솔리드는 깨끗한 manifold 박스다. ⇒ 세 결정:
- **`DISCONNECTED_RESULT` 은퇴 → 다중 솔리드 출력.** 불리언이 `Vec<Solid>`를 반환한다. "Cut 두 동강(dev-log.md 5c)"·"엣지 접촉 Fuse"가 전부 정상 결과. `live_solids: Vec<Handle<Solid>>`(§2) 인프라 이미 존재 — 불리언 쪽 제약(양수 셸 정확히 1개 요구)만 걷어 **전 positive 성분을 각 솔리드로 방출**(cavity 부호 분류 (5d) #5와 연동).
- **non-manifold 솔리드 미지원 — 불필요.** OCCT도 안 만든다(엔티티명부터 MANIFOLD). 각 솔리드 안의 manifold 전제(엣지당 면 2·Euler)는 그대로. 공유 엣지가 compound 레벨에서 4면인 것은 각 솔리드가 2면씩이라 위반이 아니다. validate가 **per-solid**면 안 건드림(전역 엣지-면 카운트면 공유 경계 오탐 → 그때만 per-solid 조정).
- **솔리드 간 경계 요소 공유가 필수.** OCCT가 접촉 엣지를 공유하듯, `Vec<Solid>`의 솔리드들이 정점·엣지·면 Handle을 공유할 수 있어야 한다 — 오버홀의 명시 공유와 같은 메커니즘이며, OCCT는 그 공유를 **불리언이 명시적으로 만든다**(전역 자동병합 아님).
- **★ 교정.** 오버홀 논의 중 "요소 공유가 non-manifold를 부른다"는 우려는 **반증됐다** — manifold는 솔리드별 판정이며, 두 솔리드가 경계를 공유해도 각각은 manifold다. 공유는 문제가 아니라 오히려 불리언이 "여기서 만난다"를 계산 없이 알게 해준다(참조-공면의 확장).

### (5d) exactness sweep — 전수조사 표 (진리를 정하는 자리에서 f64를 읽는 곳. 단일 진실원.)

M5 불리언의 **위상 결정**(어느 것이 안/밖·볼록·공면·outer/cavity인가)이 좌표 f64를 읽는 자리를 전수 나열한다. 진리를 정하는 술어만 대상이다 — 부피·validate·OCCT는 net(사후 검산)이라 제외. **정점의 존재론적 tol(`vertex_tol`)은 제외** — Discovered 정점은 "정의 + tol"이 정상 표현이며(원칙상 합법, DNA), 은퇴 대상이 아니다. 은퇴 대상은 "판정이 tol/f64 캐시를 읽어 조용히 틀릴 수 있는 곳"뿐이다.

| # | 자리 | 무엇을 정하나 | 현재 | 뿌리·비고 | 크기 |
|---|---|---|---|---|---|
| 1 | `is_convex` (ops) | 볼록성 → coincident 병합 분기 | **exact `plane_side` ((5d)-1 완료)** | 옛 tolerance가 정점 f64 비공면을 흡수. n0: 도달 피연산자 전부 Constructed 축정렬이라 이미 exact → tolerance 은퇴는 **무동작 청소**(Case A). Discovered 정점이 닿으면 #3로 넘어감 | 완료(작음) |
| 2 | `coplanar` (geom/ops) | 두 면 공면 여부 → 인터페이스 탐지 | **exact rank-1 `planes_coplanar` ((5d)-2 완료)** | 옛 `1e-9` 절대-길이 tolerance는 scale-비불변(근접-공면 false-merge 위험). exact = 두 평면 계수 2×4의 rank-1(여섯 2×2 minor). n0: 유일 불일치는 float-박스 옆면의 같은-법선 공면쌍뿐 → 반대-법선 인터페이스 필터가 배제 → **동작 보존** | 완료(작음) |
| 3 | `is_convex`의 **triple-sourcing** | Discovered 정점을 f64 캐시가 아니라 제 정의(triple)로 | **exact `three_plane_orient3d` ((5d)-3 완료)** | ★ 술어를 새로 안 지음 — `plane_orient(Q1,Q2,Q3,P)=sign(det4)·sign(det3)`는 이미 `three_plane_orient3d`(det3 인수가 술어에 내재, 순서 무관). 발화 픽스처 = **불리언 결과를 피연산자로 넣는 첫 사례**(볼록 `Common` 결과=Discovered 코너 6개 → 셋째 박스 스택). M5 동작 보존(축정렬 교점 f64-exact라 triple=캐시); 소득은 결정에서 캐시-읽기 제거 + ingestion 경화 + CIP 대비 | 완료(중) |
| 4 | `fan_triangles`의 zero-area 드롭 (ops) | 팬 삼각형 버릴지 (winding 경로) | **exact `orient2d` zero-area ((5d)-4 완료)** | winding 부호는 이미 exact(`ray_triangle_cross`=orient3d). 유일한 tol은 `fan_triangles`의 상대 `1e-12` 공선 드롭 — 대좌표 근접-공선 슬리버를 조용히 버려 교차를 놓칠 수 있었다. exact = 외적 세 성분(=세 좌표투영 `orient2d`)이 모두 0. n0: 전 코퍼스에서 nonzero-면적 발화 0건 → **Case A**. ~~광선-free 재설계~~ **기각**: 비볼록에 더 단순한 exact 광선-free 없음(GWN=f64), 방향 휴리스틱은 honest-reject라 exactness 이득 0. `ray_triangle_cross`의 정점 f64-좌표 읽기는 회전 시 CIP의 몫(이 셀 밖) | 완료(작음) |
| 5 | cavity 분류 부호 (ops, (5c)) | 성분이 outer(+)인가 void(−)인가 | **exact extreme-vertex 부호 ((5d)-5 완료)** | 옛 f64 signed-volume flux(sqrt-area·나눗셈-centroid·혼합부호 합)가 위상 라벨을 정했다. exact = 성분의 lex-최소 정점(볼록 코너)에서 `n_x<0`인 인접면이 있으면 outward. **좌표 산술 없음**(순서비교 + 축정렬 법선부호뿐) → `0.3/0.7` 비표현 좌표에서도 exact, **total**(거절 안 함). ~~containment-parity~~ **기각**(감김 재사용이나 total→partial·O(n²)); ~~exact signed-volume 합~~ **기각**(합의 exact 부호 인프라 없음). n0: 전 코퍼스(ops 217 + OCCT 70) 부호 뒤집힘 0 → **Case A** | 완료(중) |
| 6 | `interface_correspondence` (ops) | coincident-merge 정점 대응 (B링 → A링) | **exact 좌표 동일성 `ap==bp` ((5d)-3 완료)** | (5d)-2 완결성 grep이 발견 — coincident 경로의 두 번째 f64-읽기. 인터페이스 정점은 **Constructed라 triple 정의가 없다** → 정답은 좌표 동일성(대응 코너는 같은 점, bit-동일). 견고했으나(순수성 청소) scale-상대 tol 은퇴. 전 fixture Case A | 완료(작음) |
| — | `vertex_tol` (Discovered 정점) | — | tol 보유 | **은퇴 안 함.** 존재론적 tol이며 판정 fudge가 아니다. CIP가 이 tol을 전파·판정에 반영하는 층(회전 도입 시) | 해당 없음 |

**★ "평면 계수 exact화"만으로 #1이 열리지 않는다.** 축정렬 정수 입력의 Constructed 정점은 정규화 전에도 정확히 공면이었다 — (5d)-1의 계수 exact화는 **토대**(그 위에서 `plane_orient`·`coplanar`가 참 평면에 서게)일 뿐, is_convex의 과잉 거절을 실제로 없앤 것은 tolerance 은퇴다. is_convex의 남은 exact화(Discovered 경로)는 **#3 triple-sourcing**이 열고, 그건 발화 픽스처가 있어야 짓는다.

**진행:** #1·#2·#3·#4·#5·#6 **전부 완료 → (5d) exactness sweep 종료.** n0 완결성 감사(nacre-ops 프로덕션 grep)가 clean: **tolerance 리터럴이 하나도 남지 않았고**(`vertex_tol`만, 그건 존재론적 제외 행), 나머지 f64 부호 비교는 전부 exact 술어 부호(`orient3d`/`plane_side`/`planes_coplanar`)·축정렬 구성-방향·축정렬 법선 dot 부호뿐 — tol/누적 기반 위상 판정 0. dev-log.md (5b-0)/(5c) 기록의 "전수조사 표" 참조는 모두 이 표를 가리킨다.

**M6 부호 판정·내외 분류 방식 — 세 마일스톤 중 유일한 미정 (M6 직전 확정).** M5는 방법 확정(indirect predicates), M7은 방법 미정이어도 무방(SSI 실패 시 `Rejected`로 정직하게 거절하는 게 설계에 내장). 반면 **M6만 "확실히 되는 실용 영역"이라 문서가 약속했는데 정작 어떤 방법으로 판정할지가 비어 있다** — M5의 indirect predicates는 선·평면(다항식)에 특화라 이차곡면에 그대로 안 맞고, M7의 exact ray casting은 메시 경유 하이브리드용이라 M6엔 과하다.

미정인 근본 이유: 이차곡면 교차는 난이도가 갈린다. **평면∩이차곡면**(평면∩실린더=타원, 평면∩구=원, 평면∩원뿔=원뿔곡선)은 **닫힌 형식**이라 쉬운 쪽이고, **이차곡면∩이차곡면**(실린더∩실린더 등, 일반적으로 4차 공간곡선)은 특수 케이스(직교 등)만 닫힌 형식이고 일반은 이미 SSI에 가깝다. 따라서 M6는 "하나의 알고리즘"이 아니라 **"어느 곡면쌍까지를 M6로 긋고, 그 판정을 무엇으로 할지"라는 범위+방법이 얽힌 선긋기**다.

**현재 유력 방향 (확정은 M6 직전 재조사):**
- **평면∩이차곡면 교차점 판정** → indirect predicates를 **부분 확장**하는 것이 유력. 교차 곡선이 닫힌 형식이고 이차곡면도 2차 다항식이라, implicit point를 "이 평면과 이 이차곡면의 교차"로 정의하면 indirect 술어의 다항식이 복잡해질 뿐 **여전히 다항식**이라 M5 indirect가 여기까지 늘어날 여지가 있다.
- **일반 이차곡면쌍** → indirect가 버거워지는 지점. exact ray casting(M7 도구)을 앞당겨 쓰거나, 커버리지 밖으로 두고 `Rejected`로 거절. 즉 어려운 이차곡면쌍은 사실상 M7 쪽으로 미룬다.

**M6 이차곡면 교차 — 참고 논문 + 라이선스 규율 (M6 직전 사용, 지금 구현 아님).** 이차곡면 교차 위상의 exact 분류 참고처를 미리 걸어둔다(위 두 불릿의 판정 근거).

- **QI 3부작 + 구현 논문 (pencil 분류로 교차 타입을 대수적으로 exact 결정).** Dupont·Lazard·Lazard·Petitjean, "Near-optimal parameterization of the intersection of quadrics" — SoCG 2003 발표, JSC 2008 저널 3부작(Part I 생성 알고리즘 43(3):168–191, Part II pencil 분류 43(3):192–215, Part III 특이 교차 43(3):216–232). C++ 구현 논문은 별도: Lazard·Peñaranda·Petitjean, "Intersecting Quadrics: An Efficient and Exact Implementation", SoCG 2004 / Comp. Geom. 35(1–2):74–99, 2006. pencil 분류로 교차 타입(원·타원·점·두 원뿔 등)을 부동소수점 오차 없이 대수적으로 결정 — 위 "평면∩이차곡면 indirect 확장"·"일반 이차곡면쌍"의 위상 판정 참고처.
- **CGAL 3D Spherical Kernel (2차 대수적수 exact — 다른 계보, "대안 접근"으로만).** de Castro·Cazals·Loriot·Teillaud, Comp. Geom. 42(6–7):536–550, 2009. 곡면 교차점을 2차 대수적수로 exact 처리. 단 우리 indirect predicates(Attene, 정의 기반)와 다른 계보(2차 대수적수 타입)라 M5 indirect의 "확장"이 아니라 "대안 접근"으로만 참고.

**라이선스 규율 (indirect predicates와 동일, 오히려 더 엄격 — §7·§8 M5 라이선스 규율과 정합).**

- **QI 구현**(LORIA/INRIA, gamble.loria.fr/qi): "free for non-commercial use" — MIT/Apache 아님, 비상업 한정. 내부 부품(Uspensky 실근 분리 등)은 또 다른 라이선스. → **소스 열람·차용 금지**, 위 논문(3부작 + 구현 논문)만 참고.
- **CGAL Spherical Kernel**(`Circular_kernel_3`): 패키지 오버뷰에 License: **GPL** 명시(CGAL 이중 라이선스 중 상위 알고리즘은 GPL). GPL은 강한 copyleft라 링크 시 nacre 전체가 GPL 오염 → **소스 열람·링크·차용 절대 금지**, 논문만 참고. (CGAL kernel 기반부는 LGPL이나, 관심 대상 Spherical Kernel은 GPL.)
- **규율 요약:** 두 참고처 모두 **논문만 읽고 clean-room 구현**, 소스는 열람조차 안 함, 완성 후 실행 대조(dev 전용)만 허용 — Attene 2020(LGPL)에 적용한 그 규율 그대로(§8 M5·§7). exact-arithmetic 기반은 M5처럼 `geometry-predicates`(MIT/Apache) 등 자유 라이선스 부품으로 clean-room.

**★ 재조사 완료(2026-08-17) — M6a 원통-우선 확정.** 예정대로 이 메모를 꺼내 재조사했고 다음이 확정됐다: **M6는 원통부터**(M6a), 사다리는 **M6-0**(원통의 진실 — 유리수 `CylinderDef`, 불리언 없음) → **M6-1**(스칼라 원시 `sign(a+b√c)` + 평면·평면·원통 술어 + `VertexDef` branch 변종) → **M6-2**(축정렬 평면∩원통 불리언) → **M6-3**(타원 + 회전/CIP), 그 뒤 M6b(구·원뿔). 판정 방법: 평면∩원통 교차점은 좌표가 `a+b√c`(a,b,c 유리수) 꼴이다. ★ **정정(2026-08-18, M6-1 조사)**: 원래 "`sign(a+b√c)` 하나로 exact 판정이 닫힌다"고 적었으나 과장이었다 — **점-대-평면 부호와 같은-판별식 비교는 sign1(`a+b√c`)로 닫히지만, 원통 면 위 사건들의 원형 순서는 서로 다른 절단면(서로 다른 판별식 u,v)에서 온 두 점의 비교라 `ℚ(√u,√v)`의 4항 원소 부호(sign2)가 필요하다.** sign2는 재귀 탑으로 닫힌형(`P+√v·Q` 분해 → 상반 부호면 `sign(P)·sign1(P²−vQ²)` — nacre-scalar `quad.rs`). 일반 대수적수 기계도, QI의 pencil 분류도 M6a엔 여전히 불필요(QI 논문은 일반 이차곡면쌍에 가서야 꺼낸다 — 라이선스 규율은 위 그대로). 원통 위 점의 각도 매개변수는 **유리수 반각(Weierstrass) 차트**이되, ★ **차트 t의 수치는 실현(캐시)·문서화 전용이다** — t값 비교는 노름 라디칼(`|d|`, `|e₁|`)까지 끌고 와 3-라디칼이 되므로, 원형 순서는 **(seam 반평면 부호 `w·e₂`, 같은 반평면 안 외적 부호 `(w₁×w₂)·d`)** 두 술어로 판정한다(e₁ = ref_dir의 축-수직 성분, e₂ = d×e₁ — 둘 다 유리수라 sign1/sign2로 닫힘). 차트의 배제점 = seam(§9 원통 seam 절의 근거 5)이고, seam 모선 위의 점(θ=0)은 순위가 아니라 이름(`SeamIncident`)으로 답한다.

**위험도:** 미정이지만 위험하지 않다. M6도 M7과 같은 안전장치(`Rejected` 거절)가 있어, "정한 방법으로 되는 데까지만 하고 나머지는 정직하게 거절"이 가능하다. 차이는 M7은 "도박이라 열어둠", M6는 "방향은 있고 확정만 M6 직전으로 미룸".

**지금 정하지 않는 이유:** M5에서 indirect predicates를 실제 구현해보면 "이것이 이차곡면으로 얼마나 확장되는지"에 대한 감이 생기고, 그 감이 M6 방법 선택을 정확하게 만든다. 지금 확정하면 "M5 구현 경험 없이 미리 상세 설계"라는 함정. M1~M5 진행 중에는 인지만 해두고, M6 직전에 재조사해 확정한다.

**`nacre-predicates`의 geom 타입 참조 — 순환 의존 확인 (M5 직전).** indirect predicate의 implicit point는 "어느 원시 요소들의 교차인지"를 정의로 보유하는데(§4 `Origin::Discovered`), 그 정의가 `Handle<Surface>`(예: 평면 3장 교차 `[Handle<Surface>; 3]`)를 담으면 문제가 생긴다 — `Surface`는 `nacre-geom`에 있고 `nacre-predicates`는 geom보다 **아래** 층(§1)이라, predicates가 `Handle<Surface>`를 참조하면 geom→predicates→geom 순환이 된다. 이는 §1에서 `Store`/`Handle`을 최하위 `nacre-store`로 내려 푼 것과 **동일한 구조의 문제**다. M5 구현 직전에 정한다: (a) predicates가 `Handle<Surface>`를 직접 담지 않고 **좌표·평면 방정식(계수)만 값으로 받는다**(가장 단순 — predicates가 geom을 전혀 모름), 또는 (b) implicit point 정의를 geom 쪽(또는 store 같은 하위 공용 층)에 두고 predicates는 그 위에서 술어만 제공, 또는 (c) store 패턴처럼 공용 최하위 타입으로 분리. 현재 유력안은 (a) — indirect orient3d는 결국 평면 계수들의 다항식 부호이므로, Handle이 아니라 평면 방정식 계수를 넘기면 predicates가 순수 수치 계층으로 남아 순환이 원천 차단된다. **→ 옵션 (a) 확정(M5-prep).** `nacre-predicates`는 평면 계수·좌표 `[f64;N]`만 받는 순수 수치층(커널 타입 무의존, standalone 분리 가능). implicit point의 Handle 기반 정의(`VertexDef::ThreePlane([Handle<Surface>;3])`)는 위상 계층(§4 `Origin::Discovered`)에 두고, 부호 판정 시 geom이 Handle→계수를 뽑아 predicates에 넘긴다. 순환 원천 차단.

**곡면 내/외 판정 후보 (M6 직전 실측 결정).** 셋 다 M6~M7 "내/외 분류"의 후보이며 SSI 해법이 아니다(§8 M7 한계). 1순위 **exact ray casting**(§8 M7 — 우리 계보 정합, Cherchi 2022 검증, 메시 경유라 곡면 직접 판정 우회). 2순위 **GWN**(generalized winding number, Jacobson 2013 — 메시/point cloud in-out은 10년+ 검증된 성숙 기법, libigl·Axom[BSD] 구현 존재; watertight 무관 강건성이 강점이나 우리 always-closed에선 덜 필요, 느림·경계 round-off, trimmed NURBS 정확 GWN 확장은 최신[Spainhour 2024~26, 검증 진행 중] — "메시 경유 없이 곡면에서 직접 판정하고 싶어질 때"의 대안으로 보류). 참고 **graph cuts**(Diazzi/Attene 2021 — 일부 모호한 자기교차 케이스 우수), **EMBER winding number vector**(Trettner 2022). **라이선스**: 채택 전 확인, Cherchi/Attene 계열 LGPL 주의 — indirect predicates와 동일 규율(논문·MIT/Apache 소스만 참고, LGPL 소스 열람 금지, 실행 대조만 허용).

**CIP — Certified Indirect Predicates (회전 시의 부호 판정 층. 구현됨 — 오버홀 stage 1~3 완료, 수식 ②는 아래 :566에서 소진.)**

**왜 필요한가.** 축정렬 판정은 좌표가 유리수라 exact다. 그러나 **회전이 들어오면 좌표가 무리수가 된다** — 유리수 각도라도 cos·sin은 무리수이고(Niven), 임의 각도는 초월수다. exact 산술로 표현할 수 없으므로 회전 좌표는 f64 근사일 수밖에 없다. **회전은 오버홀(stage 1~3)로 지원된다** — 모델링 변환·각도 스케치가 회전 좌표를 만들고, 그때 **판정만은 조용히 틀리지 않게** 지키는 것이 CIP다.

**핵심 명제.** indirect predicates(Attene 2020)는 **선형 요소(선·평면)의 교차**, 즉 다항식에만 성립한다. 회전은 다항식이 아니므로 **indirect predicates를 회전에 적용하는 방법은 없다**(문헌에도 없다 — 수학적 한계). 실제 연구·구현이 하는 일은 하나뿐이다: **회전 좌표는 f64로 근사하고, 그 위의 판정을 exact로 유지한다.** CIP는 그 "위의 판정"에 **입력 근사(tol)까지 반영**하는 층이다. **세 계보의 하이브리드다** — Attene 2020 간접술어(implicit point 구조)·Shewchuk 1997/CGAL 필터(f64 필터→고정밀 상승)·Guibas 1989 epsilon-geometry(값+tol·sound 판정 or 기권). 유리수부는 완전 exact, 초월(회전)부는 sound 오차 한계로 부호를 **인증하거나 정직 기권**(기권은 근거를 실은 보고, §9 ③). **포기한 것은 초월수까지의 이론적 *완비성*이지 *건전성*이 아니다** — 틀린 부호는 결코 내지 않는다(no-silent-wrong). 새 패러다임이 아니라 세 계보의 특정 실현이다.

**★ 대체가 아니라 흡수다.** CIP는 indirect predicates를 지우지 않는다. 점을 "정의 + 누적 tol"로 들고, **tol = 0이면 지금의 indirect predicates 그대로(exact)**, **tol > 0이면 필터 + 고정밀 폴백**이다. 술어의 계산 구조(3-평면 정의, 행렬식 부호, 좌표 무독)는 **그대로 쓰고**, "평면 계수가 정확하다"는 가정만 "계수에 tol이 있다"로 넓힌다. 회전된 평면은 **계수에 tol이 붙은 평면**일 뿐이다.

**구성 요소 (개념 확정).**

1. **모든 점 = 정의 + 누적 tol.** 평면 교차점(정의 = triple, tol = 0)과 회전 점(정의 = 변환 이력, tol > 0)이 **같은 틀**로 통일된다. 특수 경로 없음.
2. **tol의 두 축.** 자체 tol = **이산화 오차**(값을 유한 정밀도로 표현할 때; f64 ≈ 1e-16×크기, `p`비트 ≈ 2⁻ᵖ×크기) + **연산 오차**(fraction 덧셈 = 0, f64 덧셈 = 크기 차로 비트 손실 발생, 곱셈·회전 = 발생). 누적 tol = 앞 tol + 자체 tol.
3. **★ 회전 tol은 거리에 곱해져 전파된다 — 단순 덧셈이 아니다.** 회전은 *방향*을 정하므로, 각도오차 `da`가 **그 회전 이후 판정점까지의 거리**에 곱해져 위치 오차가 된다. 따라서
   ```
   누적 tol = Σ(이산화·직선 연산 오차)  +  Σ_각회전i ( da_i × 거리(회전 i → 판정점) )
   ```
   앞선 회전일수록 뒤 이동이 많아 기여가 크다. **분리 계산은 삼각부등식 |a+b| ≤ |a|+|b|로 정당하다**(최악 상한이므로 실제보다 작아지지 않는다 — 터지지 않는다). **★ 교정(H4-soundness, 아래 미확정①)**: 위 `da × 거리`(접선형)는 **각도 불확실성** 항만 맞고, v1의 **실현 반올림(이산화) 새-오차는 접선형이 아니라 좌표혼합 `(|x|+|y|)·da`**다(독립 cos/sin 반올림의 방사 성분; 접선형은 좌표 0 근처서 과소평가—실측 반증). 유리수 각의 v1은 `da_각도=0`이라 위 둘째 항이 0이고 새 오차는 좌표혼합 항이 담당한다.
4. **★ 캐시 가능한 부분과 아닌 부분이 갈린다.** 이산화·직선 tol은 **적용점 무관**(점의 고정 속성) → 점에 값으로 캐시. **회전 tol은 적용점 의존**(같은 회전도 판정점이 멀면 tol이 커진다) → **캐시 불가, 판정 시 계산**. 그러려면 점의 **부모·조상**을 알아야 하고 조상 중 회전이 어디에 몇 번인지 알아야 한다. → **모션 이력 트리**(노드 = 모션 + 부모 링크; 여러 점이 공통 조상을 공유하는 forest).
   **★ 개정(실측 2026-07-28): 이동도, 반사도 담는다.** 원안은 *"이동·교차는 담지 않는다 — tol에 기여하지 않으므로"* 였다. 전제는 맞다 — 이동은 tol을 거의 늘리지 않는다. 틀린 것은 **결론**이다: 트리가 하는 일은 tol 계산만이 아니라 **정의의 합성**이고, 이동이 빠지면 `T(7/11)`한 벽과 `T(18/11)`한 벽이 같은 평면임을 말할 방법이 없어진다(측정: 1 ULP 갈라짐 → 한 부품이 두 몸통, OCCT는 하나). 그리고 `R` 다음 `T`는 아예 표현 불가라 거절이었다(15/15). 반사도 같은 이유로 들어왔고(§4 `Mirror`), 그래서 모션의 어휘가 **회전·이동·반사로 완결**된다 — 하나의 사슬, 하나의 규칙. 반사만 비고유라 지름길에 패리티 보정이 붙는다(§4). 교차는 여전히 담지 않는다(그건 `VertexDef`의 몫). 판정 시 노드에서 상위로 순회하며 (각도오차 × 판정점까지 거리)를 합산한다.
5. **★ 방향별 tol(x·y·z)이 필요하다.** 행렬식에서 각 방향의 오차가 **서로 다른 계수로 증폭**된다(2D orient에서 `a`의 x-tol은 `(by−cy)`와, y-tol은 `(bx−cx)`와 곱해진다). 하나로 뭉치면 정확히 증폭할 수 없다. **CGAL Lazy_kernel이 좌표를 구간(interval)으로 드는 이유가 이것**이고, 우리는 구간 산술 대신 **방향별 tol 값 + 미리 유도한 오차 한계 공식**을 쓴다(술어가 소수·고정이므로 유도가 가능하다 — CGAL은 범용이라 술어가 수백 개라 유도가 불가능해 구간을 택했다. 우리는 특화 커널이라 Shewchuk/Attene 계열의 "오차 한계 미리"가 더 빠르다). **→ 확정(미확정①): 스칼라 tol 탈락, 모든 tol을 xyz 벡터로. 회전 각도오차 환산 공식은 위 미확정① 참조.**
6. **직선 구간은 fraction 강체로 묶는다.** *(실측 확인 2026-07-28: 아래 예측이 그대로 관측됐다 — `1.0 + to_f64(7/11)`와 `to_f64(18/11)`이 1 ULP 어긋나 한 벽이 두 평면 클래스가 됐다. 다만 고친 방법은 좌표를 다시 실현하는 것이 **아니라** 이동을 정의에 넣는 것이었다: 두 벽의 정의가 같은 값을 실현하면 f64 캐시가 어긋난 채여도 판정이 일치를 증명한다. 캐시는 고칠 대상이 아니었고, 그래서 좌표는 한 비트도 안 움직였다.)* 연속된 직선 이동을 유리수로 먼저 합산하면 그 구간의 **연산 오차 = 0**이고 이산화는 **회전과 만나는 지점에서 1회**뿐이다. 순차로 f64 덧셈하면 이동마다 이산화가 붙는다((a+b)+c는 이산화 2회, a+(b+c)는 1회). 회전 결과는 무리수라 fraction으로 못 담으므로 **회전 경계의 f64 덧셈 오차는 피할 수 없다** — 피할 수 있는 것만 피한다. **(일반화 — 기준은 "직선/회전"이 아니라 "유리수/무리수"다. 직선이라서 tol 0이 아니라 유리수라서 tol 0이며, 강체 묶기는 직선뿐 아니라 **같은 축의 연속 유리수-각도 회전**에도 적용된다[각도를 유리수로 합산해 누적 각을 1회만 실현 — 2D는 항상 한 축, 3D는 같은 축만]. 유도된 무리수[√ 거리·구속 솔버 해·3D 축 변경 합성]가 섞이면 그 구간은 fraction으로 못 묶는다. tol 공식의 거리·각도 항은 지우지 말고 0으로 둔다. **★ 번들링은 필수(H4-amplification 실측): un-bundled 증분 실현은 전파 `|R|`의 행합 `|cos|+|sin|≥1`을 매 스텝 곱해 tol 바운드가 지수 폭발[실측 30스텝에 실제 오차의 ~20만 배; 실제 오차는 R 노름보존이라 평평]. 누적각 1회 실현이 차단. 폭발해도 sound 최악-보장이라 조용히 안 틀리고 상승/거절로만 간다[항목 8].)**
7. **판정 = 필터 + 정밀도 상승 + lazy 캐싱.** 점들의 tol로 이번 판정의 오차 한계를 계산(변 길이로 증폭) → `|행렬식| > 오차 한계`면 f64로 확정(대부분) → 애매하면 정밀도 상승. 관련 점이 **tol = 0(평면 교차점)뿐이면 exact 폴백**(현행 indirect predicates), **회전 점이 끼면 임의 정밀도로 상승**한다 — *초안은 여기를 "f128"이라 적었고 그건 틀렸다*(⑥): 고정 폭은 모델의 회전 이력이 얼마나 길 수 있는지를 조용히 결정하므로, 폭을 고르지 않고 **모자란 비트를 계산해 한 번에 점프**한다(astro-float, 상한 `JUDGE_PREC_CAP`). 한 번 고정밀 계산한 **값**은 캐시해 재사용하고(판정 결과가 아니라 점의 값 — 여러 판정에서 재사용된다), 재계산 시 **회전 지점만** 다시 계산하고 직선 구간은 fraction을 그 정밀도로 이산화해 잇는다. **판정에 필요한 점만 정밀화한다** — 중간 경유점은 정의로만 남긴다.
8. **오차 한계는 최악(선형 합)이다.** 확률 전파(RSS)는 실제 오차에 가깝고 100배 작지만 **보장이 아니다** — tol은 "이보다 클 수 없다"는 보장이어야 하고, 확률 tol은 드물게 부호를 뒤집어 **조용히 틀린다**. 커짐은 다른 방법으로 관리한다(직선 tol = 0이라 N은 회전 수뿐, 애매하면 필요한 만큼 정밀도를 올린다 — ⑥·항목 7).
9. **★ 필터는 동적으로 통일한다.** 지금 (5b-0)의 필터는 **정적**(오차 한계가 상수 `ε_D ≈ 5u`·`ε_M ≈ 19u`)이고, 입력 tol이 없으므로 그것이 옳다. 회전이 들어오면 **동적으로 통일**한다 — 회전 tol은 경로마다 편차가 커서 하나의 정적 상수로 잡으면 느슨해 무용하거나(과대) 위험하다(과소 → 터진다). 그리고 **각도를 쓰는 스케치도 회전이므로 무회전 케이스는 극소수**이고, **동적 필터는 tol = 0인 점을 자동으로 정적급 타이트하게 처리**하므로 정적 분기를 남기는 것이 순수 오버헤드다.

**★ 근본 한계 (정직하게).** 회전이 들어가면 **서로 다른 경로로 같은 위치에 도달해도 f64가 다를 수 있고, exact 술어로도 그 다름을 그대로 반영한다** — 술어는 "주어진 입력에 대해" 정확할 뿐 입력 자체의 근사를 고치지 못한다. **어떤 유한 정밀도로도 보장은 불가**하다(초월수의 상등은 유한 정밀도로 결정 불가) — 정밀도를 올리는 것은 답을 얻는 방법이지 완비성을 사는 방법이 아니다. 다만 **실무적으로는 충분하다**: 상용 CAD의 uncertainty가 ~1e-9인데 f64는 1000mm 점을 10000번 회전해도 ~1e-10이고, 커널이 실제로 오르는 폭(측정: 회전 3200회에서 3456비트)은 그보다 수백 자리 아래다. **이론적 완벽함은 포기하고 실무적 충분함을 tol로 보장한다** — 그리고 그 tol을 **측정된 값**으로 들고 다니는 것이 tolerance-fudge와 다른 점이다. (부호를 확정 못 하는 **잔여**는 커널이 조용히 0으로 추측하지 않는다 — 오차 한계가 `coincidence_precision`보다 좁다고 **증명되면** 일치로 처리하고 그 근거를 보고하며, 좁힐 수 없으면 원인을 이름 붙여 거절한다. 원안의 ask-user는 실측이 반박했다 — §9 ③.)

**마일스톤별 적용.**
- **M5·M6:** 판정이 다항식이므로 CIP가 그대로 얹힌다. **M6에서 CIP의 가능 여부 = 간접 술어의 가능 여부**다(CIP는 술어 위의 층이므로) — 평면∩이차곡면 간접 술어가 서면 CIP도 선다. 단 술어가 고차가 되면 **오차 한계 공식을 그 술어에 대해 재유도**해야 하고, 고차라 필터 성공률이 M5만큼 좋을지는 **측정 대상**이다.
- **M7:** 메시 조합 판정의 필터 + 뉴턴 스냅백 점의 tol 추적·고정밀 재수렴(폭은 M5와 같은 방식으로 계산해 정한다). **SSI(교차를 찾는 것)는 CIP 밖**이다(위상 존재 문제이지 tol 문제가 아니다) — 아래 "CIP와의 관계 정리" 참조.

**아직 정하지 않은 것 (구체화 대상 — 도입 직전에 정한다).**
1. **~~회전 연산의 자체 tol을 무엇으로 표현할지~~ → 확정(H4-soundness 실측 교정).** **모든 tol(회전·직선·이산화)을 xyz 방향별 벡터로 표현·누적**한다. 스칼라(구) tol은 탈락 — 정보를 버린다(z축 둘레 회전이면 z 방향 오차 0인데 스칼라는 모든 방향에 실어 "z=5 평면 위인가?"를 불필요하게 애매 판정→상승 유발). **★ 회전 tol은 세 항의 합이다** (초기 "접선형 새-오차" 안은 H4-soundness에서 반증됨 — 아래 (a)): (a) **새 오차 = 실현 반올림(좌표 혼합)** + (b) **기존 tol 방향 회전** + (c) **각도 불확실성(접선, v1엔 0)**. 회전 시 tol 갱신 =
```
새 tol = |R| · (기존 tol)            [b: 기존 tol의 방향을 회전 — 성분별 절댓값 행렬 |R| × 기존 tol]
       + (|x| + |y|) · da_이산화      [a: 새 오차 — 실현 반올림, 좌표 혼합]      ★ H4 교정
       + da_각도 × |n × (P − 축점)|   [c: 각도 불확실성, 접선 — v1엔 0]
```
**★ (a) 교정(H4-soundness).** 초기 안은 새 오차를 **접선형 `da × 모멘트암`**으로 뒀으나 **반증됐다** — cos/sin을 독립 반올림하면 순수 회전이 아니라 **방사 성분**이 생겨(`err_x = |x·δc − y·δs|`), 접선형은 좌표 하나가 0 근처일 때 그 성분을 **과소평가**한다(실측: 10000 랜덤 중 220회 상한 붕괴). 성분별 새-오차 tol은 **두 좌표를 합친 `(|x|+|y|)·da_이산화`**여야 sound(0으로 안 꺼짐; 실측 0회 붕괴·최악 ~3.5× 이내 보수). `da_이산화`=실현 반올림 한계(f64 ~16 ulp; 상승 시 ~2^-P). (b) `|R|·(기존 tol)`은 **필수** — x방향 tol이 z축 90° 회전 후 y로 옮겨가는데 안 반영하면 상한이 깨진다. (c) **각도 불확실성**(접선)은 유리수 각의 v1엔 0이나(그리고 (a) 여유가 f64 각도-반올림까지 흡수), 유도 치수(구속 솔버·비유리수 각)에 살아나므로 **항은 지우지 말고 0으로 둔다**. `n`=회전축 단위벡터, `P`=판정점. 직선 이동은 (a)(c)를 안 만들지만 앞선 tol을 싣고 간다. 판정 시 각 방향 tol이 행렬식에서 자기 계수로 증폭돼 오차 한계를 이룬다. (CGAL Lazy_kernel의 interval과 같은 이유.)
2. **~~입력 tol → 행렬식 오차 한계 공식~~ → orient3d는 확정·이식(단계 2a-ii): `frame3::det3_bound`**(6개 signed triple-product 구간 반경 `prod_err` + `16ε·mag` f64 반올림; `orient3d_judge`가 필터→astro-float 상승→정규화된 간격 판정으로 소비). exact3d H-a 검증(위반0·tightness 0.12; 프로덕션 재현 0.065·피벗 포함). `plane_side`=explicit orient3d 쌍둥이(같은 공식). **간접 orient3d**(3평면 implicit point·평면 계수)도 확정·이식(단계 2c-i `indirect_orient3d_judge`; exact3d H-b/H-c 검증·프로덕션 재현 위반0). ①이 확정돼 **입력이 명확**하다 — 점 tol = xyz 방향별 벡터. "방향별 점 tol이 행렬식에서 각자 자기 계수로 증폭돼 오차 한계를 이루는"(항목 5) 공식을 술어별로 유도(현행 정적 상수는 tol = 0 가정). **H4-soundness 진행 상황**: **(a) 새-오차·(b) 전파 항 모두 실측 검증 완료.** (a) 접선형 반증→좌표혼합(10000 랜덤 0회 붕괴). (b) `|R|·기존 tol`: 증분 회전 체인 5000회 0회 붕괴(sound·~3× 타이트); 전파의 **필요성은 비대칭 tol**에서 발현(x-tol을 90° 회전 → 오차가 y로 이동, 전파 없으면 tol_y=0으로 과소예측 — 대수적 확인; 비대칭 tol은 v1-후 √거리·구속 솔버에서). 아래 ③(평면 계수 tol 전파)도 같은 유도의 일부. **남은 H4**: consistency(경로-독립 결정성)·amplification(증폭 경계)·속도.
3. **~~평면 계수의 tol 구조~~ → sqrt 섭동은 (5d)-1이 해결, 나머지는 ②로 흡수.** `Plane::through_points`의 정규화(sqrt) 섭동은 미확정이 아니다 — (5d)-1이 Plane에 비정규화 `raw`를 함께 저장해 판정이 쓰는 `coefficients()`가 exact(정의 정점이 정확히 0)이고 `normal()`(sqrt)은 크기 소비자에게만 간다 → **판정 경로에 sqrt 섭동 없음**. 남은 **"회전된 평면의 계수 tol"**은 별개 미확정이 아니라 **①의 따름정리**다: 평면은 세 점으로 정의되고, 회전되면 그 점들이 ①의 방향별 tol을 가지며, 계수는 그 점들의 뺄셈·외적이라 **점 tol이 계수 tol로 전파**된다. ②의 오차 한계 유도의 일부로 함께 다룬다.
4. **~~3D에서 "거리"의 정확한 정의~~ → ①에서 확정.** 새-오차 (a)는 **좌표 크기 `|x|+|y|`**로 실측 확정(H4). 각도-불확실성 (c)의 접선 크기 = `da_각도 × 회전축까지 수직거리(모멘트 암)`(v1엔 0).
5. **~~방향별 tol의 자료구조~~ → 확정: `[f64;3]` xyz 벡터**(단계 2a). `nacre-cip::kernel::frame3::WitnessPoint { base:[Rat;3], chain: HpRc<[RotNode]>, coord:[f64;3], tol:[f64;3], hp }`. (검증 도구였던 2D `Pt2` 프레임은 소비자가 없어 2026-07-28에 삭제 — 불변식은 frame3의 랜덤 체인 tol 테스트와 접선형 반증 테스트가 이어받았다.) **생성자가 세 경우를 이름으로 구분한다**(2026-07-27): `exact`(좌표가 f64로 정확 → tol이 **구조상 0**) · `at`(유리수 base의 반올림을 **측정**) · `at_with_tol`(Discovered seam의 seed). `at`을 exact 경우에 쓰면 답이 0인 것을 120비트로 계산하게 되고, 그게 회전 부울 시간의 절반이었다. `chain`이 `HpRc`인 이유는 clone을 refcount 증가로 만들면서 `parallel` 빌드의 `Send`를 지키는 것(`Rc`는 `Send`가 아니다). (**회전 시 tol 변환**은 ①의 `|R|·기존 tol`로 흡수.) 임의 유리수 피벗까지 exact3d H-f로 검증(피벗 산술 tol 항 추가).
6. **~~고정밀 층의 선택~~ → 확정: `astro-float`(순수 Rust 임의정밀). twofloat(double-double)는 H1.5에서 탈락.** H1.5 실측: twofloat의 π 상수는 정확하나 삼각함수가 "preliminary"라 **영점 근처 cos 오차 ~1.8e-16(f64 수준)**으로 정확도 게이트(~1e-30) 실패 — 부호 판정이 일어나는 near-degenerate가 곧 영점 근처라 치명적. `astro-float`로 교체(160비트에서 cos/sin 오차 **~1e-58**로 통과). **정밀도가 dial 가능**이라 double-double의 ~1e-32 천장이 사라진다(H4-amplification이 tol을 키워도 정밀도를 올리면 됨 → "quad-double 상승 or 정직 거절" 딜레마가 "정밀도 dial"로 단순화). 대가는 속도 **~120µs/call**(상승 경로라 드물고, 각도당 실현을 캐시해 상각; 정밀도를 낮추면 빨라짐) — H4가 속도·상승빈도 실측. **★ 그리고 이 층은 이제 상승 경로 전용이 아니다**: `Angle::cos_sin_f64`가 libm을 부르지 않고 128비트 실현을 f64로 **정확 반올림**해 돌려주므로(`round_to_f64`), 모델을 **구성할 때부터** 임의정밀이 관여한다. 그 대신 f64 캐시가 진실의 반올림 사본이 되고, 실현이 플랫폼·빌드·호출부에 무관해진다 — 실측: `128`비트 실현 1회 ≈ 32µs(libm은 22ns)이지만 각도당 한 번이고 `trial_bound`와 **같은 `TRIG` 항목**을 쓰므로 판정까지 가는 모델은 오히려 일이 준다(회전 fold 80: 731ms → 715ms). **★ 이하 이 문서의 "double-double"·"f128"·"고정밀 상승"은 이 상승 층을 가리키는 일반명이며, 실제 구현은 astro-float다.**
6.5. **평면 witness는 정의를 소유하고, "회전됨"은 별개 필드다** (2026-07-27). `Witness::tri_pt3() -> &[WitnessPoint;3]`(캐시, 항상 존재) + `is_rotated() -> bool`(경로 선택). 전에는 `tri_pt3: Option`의 **존재가 두 질문에 동시에 답**해서 정의를 미리 저장할 수 없었다(저장하면 모든 평면이 회전으로 오인). 그래서 `plane_def`가 판정마다 정의를 재구성했고 — 핀 25개 부울에서 **108만 회** — 그것이 시간의 77%였다. 지금은 순수 접근자다. **플래그는 유도하지 않는다**: `chain.is_empty()`도 `tol == 0`도 `solid_is_rotated`와 등가가 아니며(각각 chain 없는 증언 / 90°계열 회전에서 깨진다) 서로 다른 케이스에서 틀린다. 두 필드는 `collect_planes` 한 자리에서만 함께 설정한다.

7. **회전 이력 트리 → 확정(단계 1c): `Model.rotations: Store<Rotation>` + `parent` 링크 forest.** 남은 **캐시 정책**(hp 실현값 수명·축출)은 단계 2a에서 미실행(WitnessPoint 매번 fresh 계산) → **단계 2b 이후로 defer**(성능 최적화, soundness 무관).

**★ 정리 — ②(오차 한계)가 완전 소진됐다: orient3d/plane_side 직접(2a-ii `det3_bound`)·간접(2c-i `indirect_orient3d_judge`)·cmp_coord(cmp-i `indirect_cmp_coord_judge`) 모두 확정·이식(프로덕션 `nacre-cip::kernel::frame3`).** 남은 §CIP는 ⑦-캐시(성능·defer)뿐. ①·④·⑤·⑥ 확정, ⑦-트리는 1c 완료, ③은 sqrt 부분 (5d)-1이 해결·나머지 ②로 흡수. §CIP의 마지막 수학 ②는 실험 H-a(직접 orient3d)/H-b(평면 계수)/H-c(간접 orient3d)가 검증했고, cmp_coord는 exact3d에 없던 새 수학이라 frame3에서 직접 **H-g**(두 코퍼스·wrong-sign 0)로 검증했다(Cramer 기계는 2c-i 이식분 재사용·최종 부호 결합만 신규). 세 술어 모두 프로덕션에 들어왔다 — **판정층 완성**.

★★★ **평면의 구간을 증인 없이 얻는 길 (S5(ii)-2a, 2026-08-08).** `frame3::{plane_iv_through,
plane_hp_through}` — 동차점 `[Dvec : D]` 셋의 **사영 join**(3×4 의 네 3×3 소행렬식, 부호 교대).
`cramer_iv`(세 평면 → 점, meet)의 **쌍대**라 `det3_iv`/`det3_big` 재사용뿐이고 새 산술이 없다.
소비자는 무리수 모션이 낀 datum 평면(`PlanePoints::Through`)인데, ★★★★★ **그 소비자는 배열
엔진이 아니라 «프레임 실현» 이다**(정정 2026-08-08). 벽은 둘이고 프레임이 먼저다
(`docs/truth-and-cache.md` 열린 항목 16).

★★★★ **16-1 (2026-08-09) — 프레임 벽이 순수-혼합 population 에 열렸다.** 새 부품:
`HpApprox::{div, inv_sqrt}`(구간 나눗셈·역제곱근 — 실현은 나눠도 된다, 금지는 술어의 부호
질문), **`MoveNode::FrameThrough`**(정의점 셋을 싣고 정준 기저를 구간으로 유도 — 분기는 고정
128비트 실현에서 한 번 판정해 저장), topo 의 **진술 interning**(`surface_through_ids`), ops 의
원인 3분기 + `frame_chain` 제3 도로 + `collect_planes` 이종-사슬 가지(`plane_iv`/`plane_hp` 는
원래 세 점을 각자 실현하므로 판정층 새 기계 0). 끝-대-끝 실측: 이름 없는 datum 위 프리즘이
불리언(Common)을 통과, 부피 = 스케치 진술 그대로, 비용 416 climbs(wide 454·narrow 110, 고갈 0).
`WitnessPoint` 는 **정의 동등성**(base+chain)을 얻었다 — `shared_base` 의 전-노드 비교가 그
소비자다.

★★★★ **16-2 (2026-08-09) — 걸친 정점의 datum.** `JudgedPoint { Pure, Meet }`: 정의점이 세 담체
평면의 교점일 수 있다 — 담체들의 `cramer_hp` 가 낸 동차점 셋을 **2a 의 사영 join
(`plane_hp_through`)** 이 평면으로 잇는다(그 고아의 첫 호출자). 잠금은 Pure-대-Meet 차등(같은
세 정점의 두 표기가 한 기저 — 정지·37° 회전). producer 는 straddle 을 `Nameless` 로 합류시키고
(`VerticesInMixedFrames` 는 «담체가 이름 없는 datum» 만 남음), datum 어휘의 분류-수용률이 전
population 100% 가 됐다. ★★★★★ **16-3 (2026-08-09) — 열린 항목 16 완결.** 판정 표는 `WorkingPlaneDef::Through` 가
필요 없었다: 표의 계약(«평면 위 세 정확한 점, n_out 감김»)을 판정 프레임의 probe
`(0,0,0)·(1,0,0)·(0,1,0)` 이 정의상 만족한다 — `collect_planes` 가 그 삼각형을 지으면서
걸친-datum 위의 **불리언이 열렸다**(Common 통과·부피 정확·434 climbs — Wide 와 같은 자릿수).
클래스 병합은 probe 의 동일 사슬 → `shared_base` 증명된 0. 판정층 새 기계 0.
`plane_iv_through` 는 소비자 없이 은퇴(잠금은 Pure-대-Meet 차등으로 이사).
- **나누지 않는다.** `Dvec/D` 로 아핀 좌표를 만드는 것은 무리수를 제조하는 일이고, `Approx` 에
  나눗셈이 없고 `HpApprox::div_exact` 가 반경 0 을 요구하는 것이 그 규율의 집행이다.
- **차수 9**(아핀 외적 경로는 15) — 구간 폭과 요구 정밀도가 그만큼 준다.
- **배율 부호는 값 안에서 없앤다** — 결과가 참 평면의 `D0·D1·D2` 배라 음수면 방향이 뒤집히고,
  `Judge::plane_iv` 시그니처에 부호를 실을 자리가 없다.
- **깊이 1 은 필터가 살고**(미결 0/800, 최악 상대 반경 6.3e-10) **깊이 2 는 차수 81 이라 계수가
  `f64` 밖**(8/8) ⇒ 상승 전용. 그때 `None` 을 **말한다**(`NaN` 이 우연히 옳게 구는 것에 기대지
  않는다).

★★★ **열린 항목 15 (2026-08-09) — wide 이름의 정수 rescue.** 세 판정 술어(`orient3d`·
`cmp_coord`·`plane_pair_dir_sign`)가 이름의 **BigInt 정수**로 정확히 답하는 가지를 얻었다
(scalar `int_plane_side`·`int_cmp_coord`·`int_dir_sign` — f64 쌍둥이와 같은 규약; `Expansion`
은 지수 상한 ~2¹⁰²³ 이라 wide(~2²²⁹¹)를 못 담는다). 게이트 = 전 평면 이름 ∧ 하나 이상 wide ∧
(전부 무이동 또는 전부 한 사슬; `cmp` 는 무이동만) — narrow 전용 질문은 기존 경로 그대로라
census 무접촉. 이름은 구성 시 `name_stored_ints` 가 σ(정준↔저장 방향, interning flipped 실측
1,916) × `frame_sign` 으로 접어 **`base_coeffs` 와 같은 방향**을 들고, 홀수 미러 사슬은
계수에서 `parity·C`(x 유지·y/z/d 부정 — 외적이 pseudovector)로 보정한다. ★ 구현 중 반증:
wide datum 의 이름은 **세계 이름**(담체 = 사슬 없는 발견 정점)이라 1세대 혼합-프레임
불리언에서는 부분만 잡히고(454→330), 전부-세계/통째-이동 2세대 표에서 공짜다(139 대 20 —
narrow 대조군보다 싸다). 상세·잔여는 truth-and-cache.md 항목 15.

**측정해야 아는 것:** 고정밀 층의 속도, 동적 필터의 실제 성공률(입력 tol이 있을 때 정적 대비 얼마나 자주 exact로 떨어지는가), 실무 형상에서 회전 tol이 실제로 얼마나 커지는가(판정 경계에 근접하는가).

**전제:** **(5d) 완료.** 판정이 좌표가 아니라 평면 계수·정의만의 함수가 되어야 그 위에 "계수의 tol"을 얹는 것이 깨끗하다. (5d) 전에 CIP를 시작하면 f64 캐시를 읽는 자리에 tol을 덧칠하는 꼴이 된다.

**계보:** CGAL Lazy_kernel(구간 근사 + DAG 연산이력 + 애매하면 exact 재평가 + 캐싱)과 구조가 같고, **우리는 DAG 대신 이미 있는 연산 로그(원칙 3)를 재사용**하며 구간 대신 **미리 유도한 오차 한계**를 쓴다(특화 커널이라 가능하고 더 빠르다 — Attene 2020이 선·평면 교차에서 CGAL lazy보다 빠른 이유와 같다). **지금 할 것: 없음.** 회전 도입 시 이 항목을 꺼낸다.

**M7 SSI: 작은 loop·접선 탐지 전략 (프론티어 후보, M7 진입 시 확정).** 위 "내/외 판정"은 *찾은* 교차를 어떻게 분류하나이고, 이 절은 *찾는* 일 — M7의 진짜 도박 — 이다(§8 M7 한계·구분 그대로).

**분류 먼저 — SSI 난제 4종 중 "곡면 결함"과 "연산 대상"을 가른다.** M7이 도박인 이유는 부호 판정(내/외 분류; exact ray casting으로 해결됨)이 아니라 **교차 곡선(불리언 seam)을 처음부터 찾는 SSI**다. SSI가 깨지는 대표 상황 넷:
- **① 작은 loop** — 두 곡면이 좁은 영역에서 작은 닫힌 고리로 교차. 격자·메시 간격보다 작으면 놓침 → seam 누락 → 불리언 오답. **최대 난제.**
- **② 접선 교차** — 곡면이 가로지르지 않고 스치듯 닿음(구가 평면에 접). 교차 유무 모호, 거리 0 판단이 f64로 어려움.
- **③ 특이점/cusp** — 교차 곡선이 갈라지거나 꺾임(직교 원통 X자). marching이 가지를 놓치거나 branch jumping.
- **④ 자기교차·겹침** — 곡선 아닌 면으로 포개짐, 또는 곡면이 스스로 꼬임.

**처리 정책:** ③④ 중 **곡면 자체 결함**(곡면 자기교차, 곡면 고유 특이점)은 생성 시 검사해 거부 — 연산 대상 아님. 단 **겹침**(두 정상 곡면이 면으로 포개짐, coplanar 류)과 **교차 cusp**(정상 곡면들의 교차에 생기는 특이, 예: 직교 원통)는 곡면 결함이 아니라 정상 연산 대상이므로 거부 못 함 — 별도 케이스로 처리. 실질 난제는 ①② + 겹침/교차 cusp.

**①② 탐지 — 계보와 강도 순(약→강). M7 진입 시 최신으로 비교 후 확정.**
- **거리 함수 임계점 (기준선, Patrikalakis 계보).** 한 곡면→다른 곡면의 방향 거리 함수(oriented distance field)를 만들고 그 gradient 벡터장의 위상(회전수·Poincaré index)으로 임계점을 탐지해 내부 loop·특이점을 찾음. 정통이나 (a) 거리 함수 계산이 무겁고 (b) 임계점들이 가까우면 강건성↓.
- **법선 평행점 (normal-based, 보장 있음).** 정리: 두 비특이 곡면이 닫힌 loop로 교차하면 **양쪽 법선이 평행한 직선이 반드시 존재**. 따라서 법선 평행점을 전부 찾아 그 지점에서 세분하면 **loop를 놓칠 수 없음이 증명됨**(heuristic 아닌 보장). Bézier normal vector surface로 정확·효율 계산. "간격 < 최대 변" 류 heuristic보다 강함.
- **Winding number (2026 최신, 유력).** winding number 이론 + 세분 결합. winding number는 위상 불변량이라 **격자 해상도와 무관하게** loop 존재를 판정 → 격자보다 작은 loop도 방어. 방향 거리장+gradient를 격자에 두고 고립 임계점·비고립 임계곡선을 계산, 그로부터 각 분기 시작점을 뽑아 곡선 추적. 작은 loop·접선·특이점에 강건. (ACM TOG 2026, "A Robust and Efficient Intersection Algorithm for NURBS Surfaces: Handling Small Loops and Tangent Intersections")
- **하이브리드가 최신 흐름:** winding number(탐지) + 거리/법선(정밀 위치) + 적응 세분(국소화) 결합.

**실무 파이프라인 스케치 (M7 진입 시 위 방법 중 택1로 채움).**
1. 두 곡면 간 거리를 영역별 계산(**GP 아님** — 곡면을 해석적으로 아니 직접 평가가 정확·빠름; GP는 함수 미지 시 도구라 여기선 근사만 더함).
2. 거리 **국소 최소**(줄었다 늚) 영역 = 교차 후보로 표적화. "나란히 가까이 붙은 넓은 영역"(거리 작지만 최소 아님)은 제외해 세분 폭발 방지. ← 단순 "간격<최대변" heuristic의 약점 보완.
3. 후보 영역만 적응 세분, 재귀(간격 ≥ 새 최대 변까지). 세분 리미트 설정.
4. 리미트 도달 → 뉴턴법으로 넘김: 해 수렴=교차 존재, 발산=미교차. 뉴턴 스냅백 점은 tol 있는 점이므로 **CIP에 통합**(수렴 잔차=tol, 애매하면 정밀도를 올려 재수렴).
5. 리미트에서도 확신 불가(작은 loop는 "다 찾았다"의 수학적 보장이 근본적으로 불가) → **`Rejected` 정직 거부.** 조용히 틀리기보다 거부(§8·M7 철학).

**출발점 — M7 진입 시 이것부터 읽는다:** Li·Yang·Jia, "Advances and challenges in surface–surface intersection computation — An overview", Computer-Aided Design 193:104039, 2026. SSI 분야 전체 최신 개관이라 개별 논문 여러 개보다 이 리뷰가 최적 출발점. 그 시점의 최신을 반영해 위 세 계보(거리/법선/winding number)를 재비교 후 채택. (관련: Li·Jia·Chen, "Fast Determination and Computation of Self-intersections for NURBS Surfaces", ACM TOG 44(2), 2025 — ④ 자기교차 판정·거부용.)

**CIP와의 관계 정리 (혼동 방지).** CIP(Certified Indirect Predicates)는 "다항식 판정 + 회전 tol 필터"라 **부호 판정** 층이다. M7에서 CIP가 닿는 곳은 (a) 메시 조합 판정(내/외 분류)의 필터, (b) 뉴턴 스냅백 점의 tol 추적·정밀도 상승 — 둘 다 **정밀화·판정**이다. M7의 도박인 **SSI(교차를 찾는 것)** 자체는 tol 문제가 아니라 위상 존재 문제이므로 CIP 밖이다. "찾은 것을 정밀하게"(CIP·스냅백)와 "못 찾은 것을 찾기"(SSI)는 다른 일. SSI 성공 후라야 스냅백·판정이 의미 있고, SSI 실패 시 `Rejected`.

**STEP 백엔드 교체 (비크리티컬).** 개발 중 백엔드인 step-io를 경량 AP242 출력 전용 크레이트로 교체 — 시점은 M4 이후 임의이며 크리티컬 패스가 아니다. 커널의 STEP 출력은 좁은 슬라이스(AP242, 커널 형상 엔티티만)라 최종적으로 경량 라이터가 이상적이지만, 그걸 처음부터 만드는 건 난이도가 높아 크리티컬 패스에 두지 않는다("동작 먼저 → 최적화 나중", OCCT·불리언과 같은 전략). step-io의 AP242 Ed2 스키마 지식·코드젠 타입 정의는 재활용하되 리더 로직은 배제. 교체 검증은 step-io 리더를 오라클로: "step-io 출력 vs 경량 출력"을 되읽어 비교.

**STEP 가져오기 (v1 비목표 — 핵심 가치 정합성 재검토 대상).** 단순 미구현이 아니라, nacre 정체성과 맞는지 자체를 먼저 물어야 하는 항목이다. (1) **경계가 커널과 앱을 가로지른다** — 파일 파싱·워크플로 결정은 애플리케이션 레이어 기능이고, 파싱된 형상을 유효한 b-rep으로 재봉합(healing)하는 것만 커널 몫. 즉 import는 순수 커널 단독으로 완결되지 않고 앱 레이어를 전제하므로, GD&T·메타데이터를 커널에서 뺀 것과 같은 논리로 v1 범위 밖. (2) **export보다 근본적으로 비싸다** — 넓은 스키마(AP203/214/242 각 에디션) 파싱 + 진짜 본체는 파싱 다음의 healing(외부 파일은 면들이 각 CAD의 tolerance로 느슨히 붙어 있어, nacre의 공유 엣지·단일 참조 정점 위상으로 들이려면 어긋난 면을 꿰매 유효 b-rep으로 복원 = tolerant modeling 문제 그 자체)과 근사→정확 승격(외부 근사 곡면·교차곡선을 nacre 정확 기하 진실로 올리기). (3) **히스토리 없는 dumb solid** — nacre 정체성은 "모델 = 연산 로그를 재생한 결과"인데 외부 형상엔 그 로그가 없다. 재생성·히스토리 편집 불가한 이등 시민이 되며, 이는 CATIA 등 수십 년 상용 커널조차 STEP import를 dumb solid로 취급해 보기·측정 수준으로 제한하는 잘 알려진 현상과 동일. topological naming을 포함한 v1 비목표들을 한꺼번에 끌어들인다. **결론:** import는 healing·근사→정확 승격·naming을 동반하는 별도 대형 과제이며 앱 레이어를 전제한다. v1 비목표. v2에서도 "구현 가능한가"가 아니라 "히스토리 없는 dumb solid를 받아들이는 것이 nacre 핵심 가치와 맞는가"를 먼저 물어야 하는 항목. (참고: 자기 출력을 step-io 리더로 되읽는 라운드트립 테스트는 "방금 내가 쓴 것"을 읽는 것이라 healing 문제가 없으며, 외부 파일 import와는 난이도가 전혀 다른 별개 사안 — 오라클 검증용으로 유효.) 입력이 필요해지면 step-io를 별도 리더 크레이트로 붙이며, 출력과 입력은 별개 크레이트로 둔다.

**와이어프레임/서피스 모델링 (v1 비목표 — 단 확장 접합면은 지금 명시).** nacre는 solid-first이며 always-closed(모든 결과가 닫힌 솔리드)를 전제로 설계한다. 사용자가 점·곡선·열린 면(open shell) 같은 개방 형상을 독립적으로 만들고 편집하는 wireframe/surface 모델링(예: CATIA의 Wireframe & Surface 워크벤치)은 v1 범위 밖이다. 이유: (a) 개방 요소 + sewing은 우리 always-closed 철학과 정면충돌하고(열린 면은 validate가 즉시 불량 판정), (b) 자유곡면 sewing은 healing 문제라 난이도가 곡면 불리언(3층)급이며, (c) "정밀 기계 부품→STEP" 목표엔 솔리드 모델링으로 충분하다. (주의: 점·선·면은 커널의 내부 벽돌로서는 처음부터 전부 존재한다 — 없으면 솔리드 자체가 안 만들어진다. 여기서 비목표인 것은 "사용자가 개방 요소를 독립적으로 생성·편집하는 기능"이다.)

중요 — 이것은 "ops에 함수를 추가하면 자연히 되는" 것이 아니다. 개방 형상 지원은 (1) 떠 있는 점·곡선·열린 면을 담는 자료구조(topo 확장), (2) validate의 유효성 이원화("닫힌 솔리드=엄격 검사" vs "개방 요소=느슨한 검사"), (3) sewing/healing 연산(열린 면을 tolerance 안에서 꿰매 닫는 것)을 전제한다. ops 함수는 이 셋이 먼저 깔린 뒤에야 의미가 있다.

그래서 지금 지킬 규율 (미래의 "전면 재작성"을 "국소 확장"으로 바꾸는 장치):
- **(a)** always-closed 가정을 코드 전반에 흩뿌리지 말 것 — 오직 validate의 명시적 규칙으로만 둔다. 다른 함수는 "면은 무조건 어떤 솔리드에 속한다" 같은 암묵 가정으로 지름길을 쓰지 말고 Adjacency 등을 통해 조회한다.
- **(b)** validate는 개별 불변식 검사의 모음으로 구성한다(하나의 거대 함수 금지) — 나중에 유효성 이원화가 "규칙을 켜고 끄는 일"이 되게.
- **(c)** 진실은 연산 로그다 — 형상 종류 확장은 연산(Operation) 추가이지 기존 데이터의 마이그레이션이 아니다. append-only라 옛 연산의 의미는 고정되어 있으므로, 확장이 기존 로그를 건드리지 않는다.
- **(d)** Model 최상위가 solids만이 아니라 향후 개방 요소도 담을 구조적 여지를 갖는다는 점을 인지하고 설계한다(지금 free-element 컨테이너를 실제로 만들지는 않되, Solid를 유일 최상위로 못박아 확장을 막지는 말 것).

미루는 것 (지금 상세 설계하지 말 것): 개방 형상의 정확한 자료구조 표현, sewing 알고리즘의 tolerance 처리, validate 이원화의 구체 규칙 — 이것들은 곡면 도메인(M6~M7)을 실제로 구현하며 얻는 지식(relaxation이 어디서 깨지는지, tolerance가 실제로 어떻게 도는지)이 있어야 정확하다. 정보가 없는 지금 상세를 확정하면 자산이 아니라 부채가 된다. 경계(접합면)는 지금, 내용(구현)은 그때 — 이것이 원칙이다.

**결론:** 나중에 wireframe/surface를 지원할 때 감수할 것은 "sewing·validate 이원화를 새로 구현하는 국소 작업"이지 "구조 전체를 뒤집어엎는 것"이 아니다. 위 규율 (a)~(d)가 그 차이를 만들며, 이 규율들은 wireframe/surface와 무관하게 지금도 좋은 설계(깔끔한 검증, 공짜 replay/undo)이므로 부담이 아니다.


## 10. 셀 기록 → `docs/dev-log.md`

파일 비대 방지로 셀별 진행 기록·반증된 예측을 **`docs/dev-log.md`**(개발 로그)로 분리했다. 새 셀 기록도 그 파일에 append한다(이 design.md는 설계 규칙·불변만 유지).
