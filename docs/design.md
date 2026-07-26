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
├── nacre-cip        # toleranced 부호 술어(회전): kernel(Pt3 판정) + predicate(평면 배열 술어). predicates의 쌍둥이
├── nacre-geom       # 정확 기하: Surface·Curve·교차(intersect 격리)
├── nacre-topo       # b-rep 위상: Vertex/Edge/Face/Shell/Solid·half-edge·Model
├── nacre-tess       # tessellation: 출처 태그·증분 갱신
│   └── polygon      # 평면 다각형 삼각분할: Newell + 구멍 브리징 + ear clipping
├── nacre-validate   # 불변식 검사: 오일러-푸앵카레·watertight·방향성·참조 무결성
├── nacre-props      # 질량 특성: 부피·면적(해석적, tess 무관)
├── nacre-ops        # 연산: sketch/extrude/revolve/pad/pocket/boolean. 부울 = 면당 평면 arrangement 엔진
├── nacre-step       # Model→STEP(AP242) 내보내기 어댑터
├── nacre-oracle     # [dev] OCCT 비교 하네스
└── tools/occt-helper/  # 워크스페이스 밖 헬퍼: brew OCCT(1순위) 또는 uv+OCP(폴백) — §7
                        #   OCCT는 오라클 전용 — 제품 경로에 위임 없음 (§6, §8)
```

`nacre-geom`과 `nacre-topo`가 서로를 모르게 하는 것이 중요하다. 기하는 위상을 모르고(순수 수학), 위상은 기하를 Handle로만 참조한다. robustness가 첨예한 코드(교차·분류)는 전부 `nacre-geom::intersect` 한 모듈에 격리한다. 사용자는 파사드 크레이트 `nacre` 하나만 의존하며, 인터랙티브 스크립트 앱 등은 이 워크스페이스 밖의 별도 프로젝트로 둔다.

**편의 레이어 `nacre-kit` (워크스페이스 밖, 별도 리포 — 2026-07-26 결정, 미착수).** 코드-CAD 스크립트와 커널 사이의 층: 다중 솔리드 값(compound), 값 의미론(재사용 시 `Copy` 자동 삽입), 다인수 fuse/cut/common(fold), 프로파일 헬퍼와 섬-분해 호출, 패턴·미러, 에러의 사람용 매핑, 표시 메타데이터(색·투명도 — 커널 비목표라 여기가 제자리). **Rust로 두는 이유:** 헤드리스 `cargo test`가 되고, 술어 인접 로직이 exactness 도구가 있는 쪽에 남고, 프론트엔드를 교체해도 살아남고, wasm 경계가 함수 하나로 유지된다. 경계 규칙은 overview.md의 "설탕 vs 커널 판별 기준"이고, **문법·의미론과 그 결정 이유는 `nacre-kit` 리포의 `docs/syntax.md`·`docs/decisions.md`에 있다**(여기에 복사하지 않는다 — 두 곳에 같은 내용이 있으면 어긋난다). **전제:** 남이 `Cargo.toml`에 추가해 쓰려면 파사드 `nacre`가 재수출로 채워져야 한다(현재는 이름만 예약한 placeholder — 소비자가 개별 크레이트를 path로 매다는 상태).

**공개 표면 조사 (2026-07-26, 코드 실측 — 다시 조사하지 말 것).** 외부 소비자 관점에서 무엇이 막혀 있는지 훑은 결과.
- **이미 열려 있다(막혀 있다고 오해했던 것들):** `Model`의 모든 필드와 `Vertex/Edge/Face/Shell/Solid`의 모든 필드가 `pub`이고 `Model::reachable()`→`Reachable{vertices,edges,faces,shells}`도 공개라 **위상 순회는 밖에서 된다**(`shell.faces → face.outer.half_edges → edge.bounds → vertex.point`). 피킹용 `Tessellation{by_face,by_edge,…}`·`TessTriangle.face`·`TessOrigin`, 내부 정보 `Origin::{Constructed, Discovered{tol,definition}, Rotated{base,rotation}}`·`VertexDef`·`Model::rotations`(⇒ 디버그 뷰어가 읽어야 할 것은 이미 다 읽힌다), `tessellate`·`to_obj`·`to_step`·`to_step_solid`·`validate`·`mass_props`도 공개. 플레이그라운드가 bounds를 얻으려 tessellate한 것은 불가능해서가 아니라 번거로워서였다.
- **진짜 빈 것 = 파생 값과 에러 표면:** 솔리드 AABB, centroid(`MassProps`는 volume·area만), 면 법선·면적, 면→`SketchPlane`(`ops::face_frame`이 private → 면 위 스케치 차단), 순회 편의(`he_start`가 `pub(crate)`), 그리고 §6의 거절 이유. 없으면 소비자마다 재구현하며 `Orientation`·평면 `raw` 해석에서 틀릴 수 있다. **원칙: 값을 돌려주는 읽기 전용 질의**(위상 순수성 유지). `TessConfig`는 `tol` 하나뿐 — 면별 override는 미래.

**디버그 뷰어는 커널 크레이트가 아니라 워크스페이스 밖 별도 앱이다.** 연산 로그를 입력받아 매 동작을 스텝별로 재생하며(append-only라 "N번째까지 replay"가 공짜), STEP에 안 담기는 nacre **내부 정보**(`Origin`의 `Constructed`/`Discovered`, `Discovered`의 tolerance 실측값, `Handle` 관계·인접 등)까지 시각화하는 인터랙티브 도구. 내부 자료구조에 접근해야 하므로 nacre를 **직접 링크**한다(개발 중 path 의존 → 안정화 후 version 의존, 버전별 디버깅도 자연스러워짐). 만드는 시점은 **`Discovered`/tolerance가 처음 등장하는 M5 즈음** — 그 전(M1~M4는 전부 `Constructed`)의 시각 확인은 정상 결과는 STEP→step-loupe(구조+검증), 중간·깨진 상태는 OBJ 덤프→맥 미리보기로 충분해, 인터랙티브 뷰어는 필요가 증명될 때까지 미룬다.

`Store`/`Handle`은 **최하위 `nacre-store`에 둔다.** geom도 Handle을 쓰기 때문이다 — `Curve::Intersection`(§3)이 `Handle<Surface>`를 담으므로, Handle이 topo에 있으면 geom→topo→geom 순환 의존이 된다. typed-index 저장소는 기하·위상을 전혀 모르는 순수 인프라이므로 두 층보다 아래에 격리하고, 위의 모든 크레이트가 자유롭게 참조한다. (라이선스는 MIT/Apache-2.0 듀얼 — Manifold(Apache-2.0) 알고리즘 차용과 호환.)

`nacre-scalar`는 **회전 오버홀의 근본 표현 — exact 유리수 스칼라**를 격리한다. 사용자가 입력한 치수·각도를 f64 오차 없이 정확히 보존한다(`1.1`→`11/10`, `1.1×7`=정확히 `7.7` — "얇은 막" 문제의 근본 해결). `Rat`은 `Ratio<i128>` + **checked 산술**으로, 오버플로가 §4 강등 **트리거**(값의 캐시를 f64/dd로 내리고 tol을 `Origin`에 기록; 정의는 op-log로 불변 보존). `Angle`은 유리수 deg를 mod-360 **정확 누적**(한 바퀴가 정확히 0으로 닫힘 → 스케치 닫힘)하고 90°계열은 exact 유리수 cos/sin(회전 tol 0). **★ `nacre-predicates`와 상보(겹침 아님):** predicates는 기하 행렬식의 **부호를 exact 결정**(exact-부호), nacre-scalar는 **입력 값과 유리수-순수 누적을 exact 보존**(exact-값) — 역할이 갈려 이름·층이 분리된다. **의존 결정(오버홀 최초 새 외부 dep):** `num-rational`(+num-traits)을 채택 — 성숙한 checked 유리수+gcd 약분을 제공하고, exact 유리수를 손수 구현하면 버그가 exactness 목표를 훼손하기 때문(MIT/Apache·순수 Rust). 어떤 `nacre-*`에도 의존 않는 **의존 그래프 최하단 순수 토대**. **회전 좌표의 toleranced 부호 판정**(무리수 좌표라 exact 못 하지만 부호는 f64 필터→astro-float 상승→declare-0 사다리로 sound하게 정함)은 이제 **`nacre-cip::kernel`으로 분리**됐다 — nacre-math 독립 순수 술어층, predicates의 쌍둥이(§9). **범위:** exact 값 엔진까지. 통합 값+tol `Scalar`·declare-0 ask-user 정책·커널 배선은 후속 셀.

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
/// 반환: 정련된 점과 "실제 달성한 정확도" (이 값이 Origin::Discovered{tol}의 근거가 된다).
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

pub struct Edge {
    pub curve: Handle<Curve>,
    /// 끝점 정점, 또는 진짜 닫힌 엣지(끝점 없음)면 None.
    /// **닫힌 솔리드의 원형 rim은 None이 아니다** — seam 정점을 써 Some([v, v])
    /// (start == end)로 둔다. 그래야 b-rep이 유효 CW-복합체(V−E+F=2)로 남아
    /// validate의 오일러-푸앵카레를 통과한다(실린더 rim이 대표 사례, add_cylinder).
    /// None은 seam 없는 독립 전체 원(끝점 없는 곡선) — 와이어프레임/열린 면 요소라
    /// v1 비목표(§9)이며 현재 미사용.
    pub bounds: Option<[Handle<Vertex>; 2]>,
    pub origin: Origin,
}

pub struct Face {
    pub surface: Handle<Surface>,
    pub outer: Loop,            // half-edge 순환
    pub inner: Vec<Loop>,       // 구멍
    pub orientation: Orientation,
}
```

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

**틈 없음(crack-free) 규칙:** 면의 삼각분할은 반드시 `by_edge`의 공유 polyline 정점들을 자기 경계로 소비해야 한다. 인접한 두 면이 각자 독립적으로 엣지를 샘플링하면 공유 엣지에서 정점이 어긋나 T-junction이 생기고 watertight가 깨진다 — 엣지 polyline이 먼저, 면 삼각분할이 그걸 경계 조건으로. (validate가 이를 검사한다: §7.)

갱신 규칙: 연산이 위상에 항목을 추가하면 같은 트랜잭션에서 tess에도 해당 항목을 추가한다(Fornjot의 "함께 쌓기"). tolerance를 바꾼 재계산은 tess만 통째로 재생성하고 위상·기하는 불변.

## 6. 연산 층 (`nacre-ops`) — 재계산 가능한 로그

```rust
pub enum Operation {
    // M2: 스케치 평면·프로파일 개념이 Extrude 인자로 흡수된다(별도 Sketch op 없음). Extrude는
    // 프로파일에서 새 솔리드를 만들므로 앞선 op의 면을 참조할 필연이 없다 — op간 Handle 참조는
    // 기존 면 위에 작업하는 M4(PadOnFace/PocketOnFace)에서 비로소 필연적으로 도입된다.
    // M2는 매 Extrude 결과가 닫힌 솔리드라 validate가 빈틈없이 걸린다.
    Extrude { plane: SketchPlane, profile: Profile2d, dist: f64 },
    Revolve { plane: SketchPlane, profile: Profile2d, axis: Axis, angle: f64 },
    PadOnFace { face: Handle<Face>, profile: Profile2d, dist: f64 }, // M4: 기존 면 위 — Handle<Face> 참조
    Boolean { kind: BoolKind, a: Handle<Solid>, b: Handle<Solid> },
    // ...
}

pub struct TessConfig { pub default_tol: f64 /* , 면별 override 등 */ }

// nacre-ops의 자유 함수다 — Model의 메서드가 아니다. `impl Model`은 inherent impl이라
// Model이 정의된 nacre-topo에만 놓을 수 있는데(coherence), replay는 Operation·Tessellation을
// 다루므로 topo보다 위 레이어(ops는 topo·tess 위)에 살아야 한다. Model 메서드로 두면
// §2에서 tess/ops를 Model 필드에서 뺀 것과 같은 층 위반이 된다.
///
/// 로그를 처음부터 재생. 보장: 동일 로그·동일 cfg → 동일 (모델, tess) (인덱스까지 재현).
/// 로그 중간 파라미터를 수정한 재생은 v1에서 미지원 — Operation이 원시 Handle을
/// 참조하므로 상류 수정이 하류 Handle 번호를 밀어낸다 (topological naming 문제).
/// 반환 쌍: 진실 Model + 함께 생성되는 tess 캐시(§5 커플링; Model은 tess를 필드로 담지 않으므로 §2).
pub fn replay(ops: &[Operation], cfg: TessConfig) -> Result<(Model, Tessellation), OpError>;
// M2 현재형: tess가 아직 없으므로 `replay(ops: &[Operation]) -> Result<Model, OpError>`(cfg·Tessellation 없음).
// M3에서 tess 도입과 함께 위 (Model, Tessellation)·TessConfig 시그니처로 확장한다.
```

**★ 다중 루프 프로파일 — `Profile2d`를 `{ outer, inners }`로 일반화 (2026-07-26 결정, 미구현).** 구멍 있는 스케치(도넛)와 섬이 여러 개인 스케치를 코드-CAD가 요구한다. **설탕으로 흉내내면 안 된다** — "외곽 extrude → 구멍 프리즘 Cut"은 전부 `Constructed`였을 모델을 불리언·`Discovered` 경로로 내리므로 원칙 4(tolerance는 발견된 교차에만)를 스스로 어긴다. 커널은 이미 대부분 준비돼 있다: `Face { inner: Vec<Loop> }` 존재, `nacre-tess::polygon`이 구멍 여럿을 브리징하는 삼각분할, `validate` 오일러의 `L_i` 항. 막는 것은 입력 타입 하나(`Profile2d { points: Vec<Point2> }` = 폴리곤 하나)다.
- `{ outer, inners }`(구멍 N개, 제한 없음) + `extrude`가 구멍 벽면과 뚜껑 내부 루프를 함께 생성.
- `Profile2d::from_rings(rings, fill_rule) -> Vec<Profile2d>` — 링 목록의 중첩을 exact `point_in_ring`으로 판정해(술어이므로 커널) **덩어리(섬)별 프로파일 목록**을 돌려준다(짝수-홀수 깊이: 0=재료, 1=구멍, 2=구멍 속 섬…). 채우기 규칙 선택은 호출자.
- **경계:** 커널 `Extrude` 1회 = **연결된 덩어리 1개**(외곽 + 그 구멍들). 섬마다 호출해 결과를 묶는 것은 편의 레이어(overview.md 판별 기준). 그래서 `Extrude`의 다중 바디 출력은 필요 없다.

**파라메트릭 편집의 진화 경로 (v2 이후, 지금은 기록만):** 연산이 원시 Handle 대신 계보 참조 `OpRef { op: usize, output_slot: usize }`("연산 N이 만든 k번째 면")를 담으면, 상류 수정 후에도 참조가 의미로 해석(resolve)되어 편집-재생이 가능해진다. v1에서 이를 구현하지 않되, 로그 직렬화 포맷을 설계할 때 이 확장이 포맷 파괴 없이 들어갈 자리를 남긴다. TessConfig의 tolerance도 같은 맥락에서 연산별 override(`Operation` 항목의 선택 필드)로 확장될 수 있다.

불리언은 처음부터 trait 뒤에 둔다 — 커버리지 사다리의 코드화:

```rust
pub trait BooleanEngine {
    fn boolean(&self, m: &mut Model, kind: BoolKind,
               a: Handle<Solid>, b: Handle<Solid>) -> Result<Handle<Solid>, BoolError>;
}

pub struct PolyhedralBoolean; // M5: 평면 솔리드 전용 — exact 술어로 진짜 강건.
                              //   커버리지 밖 입력은 BoolError::Unsupported로 정직하게 거절
pub struct QuadricBoolean;    // M6: 평면+이차곡면 (닫힌 형식 교차)
pub struct HybridBoolean;     // M7: 일반 곡면 — 출처태그 메시 → 조합 결정 → 스냅백.
                              //   주의: 조합 결정이 메시 해상도에 의존할 수 있으므로,
                              //   "위상 결과의 tolerance 불변"을 목표 불변식으로 삼고
                              //   위반 사례를 실패 코퍼스에 축적한다
```

OCCT는 제품 경로에 등장하지 않는다 — 역할은 nacre-oracle의 채점자(§7)뿐이다. 사다리의 각 단은 자기 커버리지 안에서 완전해야 하며, 밖은 조용히 틀리는 대신 에러로 거절한다.

**거절의 *이유*는 값에 실려 나간다 (2026-07-26 구현).** `BoolError::Unsupported { reason: RejectReason }` — 이름 붙인 거절이 크레이트 밖에서 하나의 불투명한 에러로 붕괴하던 것을 값으로 옮겼다(그전엔 태그가 `#[cfg(test)]` thread-local에만 기록됐고, 여러 지점이 거절을 울린 뒤 삼키므로 애초에 건전하지도 않았다). **소비자는 `RejectReason::class()`로 분기한다** — `NotSupportedYet`(다음 마일스톤이면 됨) / `Impossible`(어떤 마일스톤에서도 유효 솔리드가 없음) / `SuspectedDefect`(엔진 불변식이 깨짐, 리포트 대상). 변형 이름은 엔진 어휘라 로그·리포트용 안정 식별자로만 쓴다. 사람이 읽는 문장·현지화는 앱 몫이다.

성장은 `RejectReason`에서만 일어나므로 그것만 `#[non_exhaustive]`이고 `BoolError`·`RejectClass`는 exhaustive다(소비자가 완전히 처리할 수 있게). 트레이스가 불완전한 경우는 `TraceDeclined { kind, face }`가 **무엇을 못 했는지와 어느 피연산자 면에서인지**를 함께 싣는다 — 예전엔 서로 다른 10가지 사유가 전부 `COPLANAR_PAIR` 하나로 나가 커널이 틀린 말을 했다. 어떤 사유가 실제로 발화하는지는 dev-log의 사유 대장(census) 참조.

### 6.1 M5 불리언 — 두 메커니즘을 regime로 라우팅 (합성 아님)

> **★ 방향 전환(2026-07-20, 사용자 결정 — 이 절은 이행기 서술).** 아래 "공면 접촉 유무로 배타 라우팅"(detector + 공면 생존표)은 케이스가 늘수록 **경우의 수가 폭발**한다(same_ground·single-shared·containment로 실증). 그래서 design.md가 **M7**(하이브리드, §426 exact ray casting)에 두던 **arrangement/winding 통합 분류 방식을 M5 평면으로 앞당긴다** — 평면은 메시·SSI 불요라 M7 분류법을 exact plane-triple에 그대로 적용. **M5의 목표 = winding 기반 단일 arrangement 엔진**(면당 세분 → sub-face를 in/out/on 분류 → op별 Requicha keep 균일 적용; 글로벌 detector·이중 경로 소거, Requicha 규칙은 유지). 이 방식은 커토버로 프로덕션 엔진(`nacre-ops::arrangement`)이 됐다 — 면당 평면 cell-복합체(트레이스→분할→셀→중첩→라벨→방출→조립). **잔여 커버리지 갭:** containment(seam 없는 포함)·cavity 접촉(현재 outer shell 순회)·회전 접촉.
>
> **★★ 갱신(F2 — detector 붕괴 완료).** 위 "경우의 수 폭발"의 실체가 **라우팅뿐**이었음이 실증됐다: 케이스별 detector 7종이 **전부 같은 `coplanar_result_unified`를 호출** — 결과 생성기는 이미 통합돼 있었다. 그래서 dispatch를 **하나의 exact 질문**(`coplanar_contact_count >= 1`, "진짜 공면 접촉인가")으로 붕괴시키고 detector 12종·지원 타입 ~775줄을 삭제했다. 부수적으로 그 좁은 게이트들이 막던 케이스가 열렸고(비볼록 오버행 footprint, 관통 slot/corner cut), 라우팅이 넓어지며 드러난 "성공하되 열린 셸" 한 건은 **`assemble_fuse_cut`의 닫힘 가드**(모든 모서리 정확히 2회 사용, 위반 시 `NON_MANIFOLD_EDGE` 정직 거절)로 차단했다. ∴ **아래 배타-라우팅 서술은 여전히 유효하되 "detector 다발"이 아니라 "질문 하나"이며**, 두 메커니즘(seam / coplanar)의 공존은 폐기 대상이 아니라 **구조적 필연**이다 — (단계4 SoS Cell 4) 실측이 "seam 하나로 통일"을 반증했다(공면 접촉 모서리는 공유 평면에 통째로 누워 transversal seam 자체가 부재 = 구조적 공면성, 섭동으로 해소 불가). 남은 winding 작업은 엔진 대체가 아니라 **커버리지 확장**(공면 arm의 회전·cavity·containment)이다.

M5 `PolyhedralBoolean`은 **능력이 겹치는 두 메커니즘을 "공면 접촉 유무"로 배타 라우팅**한다(한 부울에 하나만 돈다).

- **일반 seam 엔진**(`general_boolean`, ray-cast `classof`): 공면 접촉이 **없는** 순수 transversal 전용. 공면은 door에서 정직 거절(`COPLANAR_PAIR`/`VERTEX_ON_FACE_PLANE`) — 평면 위 정점은 `classof`(3D winding)가 미정의라서다.
- **통합 공면 처리기**(`coplanar_result`): 공면 접촉이 **있는 모든 경우**를 통째로 만든다(벽·바닥·관통까지 self-contained; 두 엔진 합성 아님 — 합성은 on-plane `classof` 블로커로 불가). **일반 규칙(케이스가 늘어도 코드가 안 는다):** 각 면 F를 상대 solid의 **F-평면 단면(`section_of_solid`)** 에 대해 `coplanar_reconstruct`로 클립(`clip_face_to_section`). a-면=`a∖b`(Cut)/`a∩b`(Common)·b-면=`b∩a`(flip은 Cut만)·접촉면=상대 footprint 직접(mouth). 접촉면(π)은 단면이 퇴화하므로 footprint를, 벽/바닥만 단면을 쓴다.

**생존 규칙(접촉면 keep/flip = 연산 × 상대법선):** Fuse/opp=`inP⊕inQ`, Fuse/same=`inP∨inQ`, Cut/opp=`inP`, Cut/same=`inP∧¬inQ`, Common/opp=∅, Common/same=`inP∧inQ`. 볼록성 항이 없어 비볼록·다중루프를 half-space 가정 없이 `point_in_ring` parity로 정확 처리.

**pair 선택:** 평면-coplanar 면쌍 중 **풋프린트가 실제로 겹치는(포함 또는 경계 교차)** 쌍만 genuine 접촉(`footprints_overlap`) — 우연히 같은 평면에 있으나 떨어진 쌍(예: slot top ∥ 먼 포켓 바닥)은 배제. **면-평면 분기:** 면이 상대 면과 coplanar면 단면이 퇴화 → 단면 대신 **면 centroid `point_in_solid`**(disjoint라 strict)로 whole/drop.

**스코프 접촉면-flush(공변 — 접촉면 위 두 모서리가 한 선에 포개짐, cut-overhang의 본질):** ① **R0** — `section_of_solid`의 raw triple을 `plane_classes` canon으로 remap(안 하면 mouth·breach벽·b-벽의 공유 코너가 어긋나 조용히 non-manifold), fold 충돌은 `SECTION_TRIPLE_COLLISION` 거절. ② **R1** — 접촉면 위 attachment graze(section 코너가 벽 모서리 strictly 내부)를 거절 대신 F-edge-split crossing으로. ③ **R2** — 접촉면 chord 위 ∂P 정점을 covered로(삼킨/유지 코너), chord 끝점 4-평면은 `FLUSH_VERTEX_COINCIDENT` 거절. 이 셋으로 벽 네 구성(Middle/slab/Corner/Shorten)이 한 규칙으로, 코너 기둥이 canon-triple로 자동 용접.

**★ glue vs clip 경계(정직):** 위는 전부 **clip**(겹치는 풋프린트를 오림). **coincident**(동일 풋프린트 스택)는 **glue**(맞닿은 면 둘 제거 + 옆벽 splice + 인터페이스 정점 remap) — 구조가 달라 통합 안 함, `coincident_merge`로 별도 유지. **정직 거절(후속):** 회전 공면(toleranced declare-0), 다중-loop 단면(슬롯이 cavity 관통), 임의 fan degree-≥3 flush-shared-boundary(`turn_at` 미구현), 다중 genuine 접촉, P⊂Q 대칭. 전부 named 태그 — silent-wrong 0(DNA).

**은퇴:** clip bespoke 5경로의 **detector는 "진짜 공면 접촉인가(blind·not-pierces·convex)" gate로 유지**(transversal을 `general_boolean`으로 걸러 silent-wrong 방지)하되 결과는 `coplanar_result`가 만든다 — builder+헬퍼(proj2/clip/splice 계열)는 삭제.

## 7. 검증·오라클 인프라 (`nacre-validate`, `nacre-oracle`)

1일차부터 CI에 들어가는 것들:

```rust
// 모든 연산 직후 자동 실행되는 불변식 (디버그 빌드에서 강제)
pub fn validate(m: &Model) -> Vec<Violation>;
// - 오일러-푸앵카레: V − E + F = 2(S − G) + L_i  (L_i = 면 내부 루프 수.
//   관통구멍 정육면체로 검산: V16 − E24 + F10 = 2 = 0 + L_i(2) ✓.
//   내부 루프 항을 빠뜨린 축약형을 쓰면 정상 모델을 불량 판정하게 된다)
// - 모든 Loop 닫힘, half-edge 짝 정합, 방향 일관성
// - Handle 참조 무결성 + model-id 일치 (디버그)
// - Constructed 정점: 참조 곡선/곡면 위에 정확히 놓임 (eps_machine)
// - Discovered 정점: tol 이내
// - tess 출처 정합: OnFace{uv} 평가값과 pos 일치
// - tess 틈 없음: 인접 면 삼각분할이 공유 엣지 polyline 정점을 정확히 공유 (§5)
```

**validate는 store 전체가 아니라 live 도달가능 셀을 센다(§2 supersede 의미론).** Euler의 V·E·F·S·L_i는 `store.len()`이 아니라 `live_solids`에서 도달 가능한 정점·엣지·면·셸·내부루프 수이고, manifold(엣지 정확히 2회)도 "전역 2회"가 아니라 "도달가능 집합 내 2회"다. 이유: supersede된 옛 셀이 아레나에 남으므로, store 길이로 세면 죽은 셀이 Euler를 깨고 죽은 면의 엣지가 manifold를 깬다. **참조 무결성 검사(dangling handle)만 store 전체를 훑는다** — 안전 게이트이자, 살아있든 죽었든 handle이 OOB면 버그이기 때문. 이 검사가 먼저 단락하므로 이후 도달가능성 순회는 항상 in-bounds라 안전하다. (부수효과: store에 떠 있는 stray 셀은 이제 불량이 아니라 "죽은 아레나 항목"으로 정당하게 무시된다 — M1~M3의 "stray → Euler 불량" 판정은 이 의미론으로 갱신된다.)

속성 기반 테스트(proptest): 랜덤 유효 연산열을 생성해 (1) validate 통과, (2) replay 멱등성 — 같은 로그·같은 cfg → 동일 모델, (3) 변환 불변량 — 강체변환 후 부피·면적 보존, (4) 불리언 대수 — `A ∪ A = A`, `A ∩ ∅ = ∅`, `vol(A∪B) + vol(A∩B) = vol(A) + vol(B)`.

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

**exact-arithmetic 바닥 확정(M5-prep, 코드로 실증).** `geometry-predicates`(elrnv, 0.3, MIT/Apache)가 finished `orient3d`뿐 아니라 **Shewchuk expansion primitive**(`two_product`·`two_sum`·`expansion_sum`·`scale_expansion_zeroelim` 등, `predicates` 모듈에 공개)를 노출함을 `nacre-predicates` 뼈대가 실제 호출로 확인 → 그 위에 indirect 술어를 쌓을 수 있고 expansion 산술 자체 구현이 불필요. (primitive는 `[lo, hi]` 순서. orient3d 부호 규약 = `det[a−d, b−d, c−d]`, 뼈대 golden으로 고정.) **M5 서브유닛 사다리:** ① `nacre-predicates`(뼈대→implicit point[3-plane]+indirect orient3d, direct 대조 property test) → ② `nacre-geom::intersect`(평면∩평면=닫힌형식 직선, 3-평면 꼭짓점) → ③ `PolyhedralBoolean`(fuse/cut/common, 세 엣지가 한 Vertex Handle 공유로 봉합 우회, `Origin::Discovered.definition` 필드 도입, 커버리지 밖 `Unsupported` 거절) + 오라클 부피·불리언대수 proptest.

알고리즘 참고: Manifold(Apache-2.0 — 차용·번역 가능), Hoffmann 등 문헌. 커버리지 밖 곡면 불리언은 명시적 미지원 에러로 정직하게 거절. `Discovered` 경로·국소 tolerance·relaxation 실전 투입. 이 시점에 "OCCT 없이 직동하는, 실용적 기계 부품(평면 위주)을 STEP으로 내보내는" 진짜 커널이 된다. 참고: Truck 대비 벤치마크·정밀도 비교(전역 1e-6 폴리라인 vs 정점별 실측 tol + 닫힌 형식)는 수치로 보여줄 수 있는 차별점 — 공개 지표 후보.

**M6 — 자체 불리언 2단: 이차곡면.** 평면∩실린더(타원), 평면∩구(원), 평면∩원뿔 — 여전히 닫힌 형식이라 행진 불필요. 실린더∩실린더는 특수 케이스(직교 등)부터. 실제 기계 부품 면의 대다수가 평면+실린더+원뿔이므로, 여기까지로 실용 커버리지의 대부분을 확보한다.

**M6/M7 난이도 절벽.** M6(평면+이차곡면)는 닫힌 형식 교차라 SSI의 지옥(위상 판정·시작점 검출·watertightness)을 **대부분 회피** — 여기까지가 "확실히 되는 실용 커널". M7(일반 NURBS)은 그 60년 미해결 문제(아래)로 **직접 진입**. 그래서 "**M6까지가 실용적 종착점, M7은 별도의 장기 연구 트랙**"이라는 선을 분명히 긋는다. 실용 가치는 M5~M6에 있고, M7이 안 돼도 M6까지로 실제 기계 부품 대부분을 커버한다.

**M7 — 자체 불리언 3단: 일반 SSI (연구 구간).** nacre-geom::intersect에 SSI 행진 구현 → 일반 곡면쌍의 `HybridBoolean` 완성: 출처태그 tess에 강건 메시 불리언(exact predicates) → 조합 결정 추출 → 살아남은 면은 정확 곡면 유지, 신규 엣지는 국소 SSI 스냅백. OCCT 오라클과 상시 diff. 실패 케이스 코퍼스 축적.

**M7 내/외 분류 = exact ray casting (1순위 후보).** 하이브리드 파이프라인("출처태그 메시 → 조합 결정 → 정확 곡면 스냅백")에서 "이 patch가 최종 솔리드 안인가 밖인가"를 분류하는 단계에 exact ray casting을 쓴다. 메커니즘: 레이가 삼각형 **내부**를 지나면 삼각형 방향(정점 순서)으로 정확 판정, 꼭짓점·엣지·접선(coplanar) 같은 애매한 케이스는 레이를 **수치 섭동**해 항상 "내부 통과"로 되돌린다(Simulation of Simplicity, Edelsbrunner–Mücke 1990). 1순위 이유: (a) **우리 계보와 정합** — 핵심이 orient 술어(채택)+섭동(`RelaxError::Tangential`로 올리기로 한 퇴화 구역의 정석)이라 새 수학을 안 들인다; (b) **검증됨** — Cherchi 2022(Interactive and Robust Mesh Booleans)가 핵심으로 채택, 수백만 삼각형·수백 입력 variadic까지 테스트하며 GWN류를 명시적으로 제침; (c) **파이프라인에 그대로 꽂힘** — 우리 하이브리드는 이미 메시를 경유하므로 분류를 메시 단계에서 함(Cherchi 검증 형태 그대로), 곡면 직접 판정 불필요라 GWN의 trimmed NURBS 확장(최신·검증 진행 중)을 우회.

**오해 방지 셋.** (1) **Truck 방식이 아니다** — 표면적으로 메시를 쓰나 정반대다: Truck은 교차를 폴리라인 근사로 표현해 그 근사가 **최종 결과**(전역 1e-6), 우리는 메시를 **분류용 임시 도구로만** 쓰고 exact 술어로 분류 후 **정확 곡면으로 스냅백**. 차이는 "메시를 쓰느냐"가 아니라 "메시가 최종이냐(Truck) vs 임시냐(우리)". (2) **indirect predicates의 대체·확장이 아니다** — indirect predicates=M5 **점** 부호 판정(평면·국소·교차 계산 중), exact ray casting=M7 **patch** 내/외 분류(전역·분류 단계). 대상(점 vs 덩어리)도 마일스톤(M5 vs M7)도 다른 별개 부품. (3) **"어려운 케이스에만 켜는 정밀 모드"가 아니다** — relaxation 사다리·적응 술어는 케이스 단위로 "어려우면 더 정밀하게", ray casting은 M7에서 **상시 도는 분류 단계**. "선택적"의 단위는 케이스가 아니라 **마일스톤**(M5·M6엔 닫힌 형식이라 불필요, M7에서 켜짐).

**한계 — ray casting은 분류를 풀지 SSI를 풀지 않는다.** M7 "내/외 분류" 절반을 검증된 방법으로 채워도, 진짜 도박인 **SSI(곡면 교차 곡선 계산)**는 미해결로 남는다. ray casting은 "조각을 어떻게 분류하나"만 정확히 할 뿐 "곡면 교차 곡선을 어떻게 정확히 뽑나"를 풀지 않는다. M7이 도박·연구 구간인 것은 그대로고, ray casting은 그 도박의 한 부품(분류)만 안정화한다. **덧붙임 — "스냅백"과 "SSI"는 다른 일이다.** SSI(교차 곡선을 처음부터 **찾는** 것)가 도박인 미해결 문제이고, 스냅백(이미 찾은 교차 근처에서 정확 곡면으로 **되당기는** 것)은 그다음의 별개 작업이다. 즉 M7의 위험은 "스냅백"이 아니라 "SSI로 교차를 정확히 찾을 수 있는가"에 있다. 스냅백은 SSI가 성공한 뒤라야 의미가 있으므로, SSI가 실패하면 스냅백까지 갈 것도 없이 그 연산은 `Unsupported`로 거절된다.

**SSI는 검증된 해법이 없는 60년 미해결 문제.** 일반 곡면 교차는 1960년대부터 연구됐으나 CAD가 받아들일 만큼 강건·신뢰할 해가 여전히 없다 — **Parasolid·SISL·IRIT 등 최고 커널조차 특정 위상 케이스(작은 loop, 접선, cusp, 미세 자기교차)에서 실패**하고, 2025년 논문들도 watertightness 보장이 "여전히 challenging"이라 적는다. M7 "도박" 표기는 과장이 아니라 정확한 현실 인식이다. 아래는 **"채택할 검증된 기법"이 아니라 "M7 진입 시 읽을 프론티어 후보"**로만 걸어둔다(저자 벤치마크 주장일 뿐 독립 검증·프로덕션 채택 전 — 지금 상세 검토는 "정보 없이 미리 설계" 함정이라 M7 진입 시로 미룸): winding number + subdivision 시작점 검출(작은 loop·접선 branch 놓침 완화, 2026), Dixon matrix tracing(branch jumping/missing을 root-solving으로, Chen 2025), interval algebraic topology analysis(SSI 위상을 4D 대수계로 분류, Cheng 2023), lower-dimensional formulation gap control(Wang 2025).

M7은 열린 연구임을 명시한다. M6까지가 "확실히 되는" 영역, M7은 이 커널의 존재 이유이자 도박이다. OCCT의 역할은 전 구간에서 **오라클(dev 전용 시험관)뿐**이다 — 제품 경로에 OCCT 위임은 없다. OCCT 소스는 "상용급이 이 케이스를 어떻게 다루나" 열람용 참고서로만 쓴다: LGPL-2.1이므로 번역·차용은 라이선스 오염이고, 무엇보다 OCCT의 위상·tolerance 아키텍처가 딸려 들어와 nacre 설계와 충돌한다. 읽되 베끼지 않는다.

## 9. 미결 사항 (다음 논의 대상)

트리밍 곡면의 pcurve 표현 시점(M3에 선행 도입 vs M5까지 지연), 닫힌 엣지의 seam 처리(방식 확정 — 아래 "원통 seam: A vs B" 항목), Sketch 제약 솔버의 범위(초기엔 무제약 프로파일만), OpRef 계보 참조의 도입 시점과 직렬화 포맷 여유분, `Store` 스냅샷·직렬화 포맷(자체 vs STEP 재활용) — 이와 함께 **세션 중 메모리 관리: compact보다 재구축(rebuild-from-log) 우선**(§2 "재검토 예정 (v2)" 참조; 재구축=주력 정리·undo 유지, live 폐포 필터링=저장, compact=비상 회수), OCCT history → 출처 매핑의 실제 충실도(M5에서 실측 필요), 멀티스레딩 경계(Store가 &mut 독점인 설계라 연산 단위 병렬은 미지원 — 의도적 단순화).

**원통 seam: A(seam 엣지) vs B(seamless periodic) — A 채택 (M3.3a에서 확정).** 주기 곡면(원통·구·토러스)의 옆면을 b-rep로 담는 두 방식이 있고, 둘 다 유효 AP242·유효 CW-복합체다(초기 판단에서 "B는 오일러가 깨진다"고 봤으나 **오류** — 깨지는 건 정점조차 없는 제3의 변형 C[`bounds:None`, V0]이고, B는 각 원을 seam 정점 `Some([v,v])`로 두고 옆면을 두 루프[outer=아래원, inner=위원]로 담아 `V−E+F−L_i = 2−2+3−1 = 2`로 통과한다).
- **A(채택):** 위 원 + 아래 원 + **세로 seam 직선 엣지**, 옆면 = 4-엣지 단일 닫힌 루프 `[bottom, seam, top⁻, seam⁻]`(seam이 같은 면에서 2회 반대 — self-adjacent). `add_cylinder`가 이 방식. STEP 출력도 A(OCCT 계열이 생산·기대하는 형태; step-io 검증 테스트와 동형).
- **B(미채택):** seam 엣지 없이 위·아래 원 두 개로만 옆면 경계(옆면이 두 루프). NIST 샘플 계열.

**결정 근거(정직한 저울질):**
1. **불리언(이 커널의 존재 이유) — 결정적.** 조합적 b-rep 불리언은 면 루프를 순회하며 교차곡선에서 엣지를 쪼갠다. B의 주기 곡면은 교차곡선이 u=0/2π를 가로지를 때 쪼갤 엣지가 없어 특수 처리가 필요하고, **그래서 OCCT·Parasolid·ACIS가 전부 내부적으로 seam을 넣는다** — 조합 알고리즘이 "모든 면 = 실제 엣지로 둘러싸인 유계 영역"을 요구하기 때문. B는 dumb solid의 교환 표현, A는 연산하는 커널의 내부 표현. nacre는 연산하는 커널(M6에서 원통이 불리언에 진입)이라 A.
2. **tess 정합.** nacre tess는 "공유 엣지 polyline(`by_edge`)을 경계로 소비"해 crack-free(§5). A는 seam이 진짜 엣지라 그 polyline을 옆면이 u=0·u=2π 양쪽에서 소비 → **기존 메커니즘 그대로, 주기 특수처리 없음**. B는 "이 면은 periodic이니 u-wrap 봉합" 특수처리를 tess에 새로 요구. 즉 **우리 설계에선 A가 tess도 더 단순**("B=tess 간결"은 tess가 주기 곡면 네이티브일 때만 참).
3. **"가짜 모서리" 비용은 흡수됨.** seam은 각진 모서리가 아니라 파라미터 이음매(양쪽이 같은 곡면이라 C¹ 매끄러움)다. 하지만 매끄러운 엣지는 seam 말고도 존재(접선-연속 fillet 등)하므로 nacre는 "sharp vs smooth 엣지" 판정을 **어차피** 갖는다 — 판정식 = 인접 두 면의 곡면 법선이 엣지를 따라 일치하는가. seam은 양쪽 같은 실린더라 자동으로 "smooth"로 판정된다. **seam 전용 플래그 불필요; 일반 술어에 흡수.** 따라서 A가 만든 "엣지 예외"는 새 범주가 아니고, B가 만드는 "면 예외(주기 경계)"는 불리언에서 훨씬 비싸다.
4. **위상 균일.** 모든 면이 실제 엣지의 닫힌 루프로 둘러싸임 → validate·루프순회·불리언이 단일 규칙. self-adjacent seam은 validate가 이미 무수정 통과(엣지 2회·반대만 보고 면 구별 안 함).

**어댑터 분리(옵션 보존):** 내부를 A로 두어도 nacre-step은 어댑터라 **필요 시 export 시점에 unseam해 B로 내보낼 수 있다**(내부 표현 ≠ 교환 표현). 지금은 A 출력(OCCT 상호운용 안전)이고, 특정 상호운용 요구가 생기면 B 출력을 어댑터에 추가. **tess 층의 남은 세부**(seam polyline 샘플 밀도·법선 비분리 규칙 구현)는 M3 tess 곡면 샘플링에서 확정.

**패드/포켓 ↔ 불리언 통합 — feature = tool body + boolean (✅ 완료).** `pad`/`pocket`은 상용 CAD와 동형으로 이미 통합됐다: **`pad` = 프로파일 압출 → `Fuse`, `pocket` = 압출 → `Cut`**(`extrude_and_boolean`), 패드/포켓은 불리언의 얇은 sugar다. `raise_region`(M4 직접 구성)은 폐기됐고, "프로파일이 면 안에 있어야 한다"는 컨테인먼트 제약도 제거돼(`extrude_and_boolean`이 미검사) 오버행 패드/포켓이 불리언의 공면-접촉·오버행 경로(`detect_contained_contact`·`detect_pocket_contact`·`detect_overhang_contact`)로 자동 처리된다. 공면 복잡도는 **불리언 하나로 집약**됐다. **남은 격차(후속)**: 오버행 `Cut`/`Common`의 비볼록 kept-solid는 아직 미검증(재보지 않았다). ★ **2026-07-22 갱신 — 셋 중 둘은 닫혔다:** 비볼록 오버행 footprint(옛 볼록 게이트)와 정확 flush-edge(프로파일 테두리 = 면 테두리 공유)를 **한 형상이 동시에** 통과한다 — `a_non_convex_pad_cantilevers_and_runs_flush`(부피 1.625·면적 9.75, 손계산과 일치). corner-flush `Common`(공면 3면, 옛 `vertex_on_face_plane` 거절)도 열렸고 **OCCT와 부피·면적이 일치**한다(`a_corner_flush_common_keeps_the_non_convex_overlap` + `corner_flush_common_matches_occt`). **잔존 직접-구성은 없다 (2026-07-22 은퇴).** `ImprintSketch`(면을 재료 가감 없이 분할하던 별개 op)가 마지막 직접-구성이었는데, `raise_region` 폐기 이후 **소비자가 하나도 없었다**(`region_face`를 읽는 코드 0). 연산과 `imprint`/`prepare_face_split`/`finish_split`/`placed_profile`(엄격 컨테인먼트)·`OpError::ProfileNotContainedInFace`를 제거했다. **그것이 커널 유일의 동일평면 인접면 생산자였으므로**(불리언 출력은 `unify_coplanar_faces`가 항상 정리한다) 불리언의 이음선 처리도 함께 제거했다 — 되살릴 근거(유도·측정)는 dev-log에 남아 있고, 복원은 그 두 커밋을 revert하는 것이 출발점이다. Split Face가 필요해지는 소비자(구역별 재질·FEA 경계조건·금형 파팅라인)는 아직 로드맵에 없다.

- **★ 구조적 공면 = O(1) 참조 인식(핵심).** 면 위 스케치를 압출한 tool body의 밑면은 대상 면의 **surface Handle을 공유**한다(`build_prism`이 대상 면의 `surface_h`를 밑면에 그대로 쓴다). 통합 불리언은 **먼저 Handle 공유를 검사**(float 계산 0)해 공유 면을 seam으로 즉시 채택하고, **공유하지 않는 독립 솔리드만** `planes_coplanar` 기하 감지로 내려간다. 즉 "참조로 아는 공면은 공짜, 우연한 공면만 계산" — 상용 커널이 "coincident 면을 imprint로 공유 토폴로지로 승격"하는 것의 nacre판.
- **★ 게이트(순서 강제).** (a) 불리언 공면 처리를 하나로 흡수·완성(지금 셀 사다리 — Cut 3갈래 통합, Common이 Cut의 detect 재사용) → (b) Handle-공유 O(1) fast path 추가 → (c) 그제서야 pad/pocket을 extrude+불리언 wrapper로 바꾸고 M4 직접 경로 제거. **역순 금지**: 불리언 공면이 robust해지기 전에 M4를 걷어내면 지금 되던 포켓/보스가 커버리지 구멍에 빠진다(M4 직접 경로는 그때까지 신뢰 가능한 fallback).
  - **(c-1) pocket 완료 (실행됨).** `PocketOnFace`는 `build_prism`(inward, top-flush) + `boolean(Cut)`으로 재구현됐다 — contained-coplanar Cut 경로가 빈 seam을 내므로 결과는 전부 Constructed(tolerance 0), 옛 직접 `raise_region`과 등가. 선결로 `detect_pocket_contact`의 볼록성 게이트를 제거(비볼록 kept `a`·비볼록 프로파일 `b` 둘 다 열림, OCCT로 교차검증)했다. `raise_region`은 `pad`용으로 잔존. through-pocket(`dist ≥ 두께`)은 이제 정직히 거절(M4는 미검사 UB였음 — 개선). **게이트 (b) Handle-공유 fast path는 불필요로 폐기**: contained 경로가 이미 빈 seam·전부 Constructed라 참조-fast-path가 더할 exactness가 없다. 다음: pad(=extrude+Fuse), 그다음 M4 직접 경로 완전 제거.
  - **(c-2) pad 완료 (실행됨).** `PadOnFace`도 `build_prism`(outward, top-flush) + `boolean(Fuse)`로 재구현 — pocket과 대칭. 선결로 `detect_contained_contact`의 볼록성 게이트를 제거(비볼록 kept base·비볼록 프로파일 boss 둘 다 열림, OCCT 2종 교차검증). pad·pocket은 이제 부호(±dist)·`BoolKind`·복원-면-부재 에러만 다른 **공통 헬퍼 `extrude_and_boolean` 위의 얇은 wrapper**로 통일됐다("feature = tool body + boolean"의 코드화). `raise_region`은 유일 사용자였던 pad가 떠나며 **삭제**(그 죽은 `Split` 필드도 함께 정리). 남은 M4 직접 기계는 imprint(`imprint`/`prepare_face_split`)뿐 — 완전 제거는 imprint 정리 후 별도 스텝.
  - **(c-3) 오버행 pad/pocket 획득 (실행됨 — 로드맵 payoff).** `placed_profile`을 `placed_profile_unchecked`(CCW+배치)와 strict-containment wrapper로 쪼개고, `extrude_and_boolean`이 unchecked를 쓰게 했다. 이로써 **면 경계를 넘는 프로파일이 기존 오버행 불리언 사이드카로 자동 라우팅**(Fuse: 단일 엣지·코너·spanning slab; Cut: N-wall blind)된다 — "프로파일이 면 안에 있어야 한다"는 제약이 사라졌다. contained 경로는 완전 불변(같은 base_pts). 커버리지 밖(비볼록 오버행·through·far·비축정렬)은 **정직 거절**(`Boolean(Unsupported)`) — silent-wrong 아님. 오버행은 contained(빈 seam)와 달리 footprint crossing에서 **Discovered seam 정점**을 만든다(기존 오버행 셀의 성질, 새 tolerance 아님). imprint는 containment **유지**(유효 inner-loop 필요). **알려진 비대칭**: contained pad/pocket은 비볼록을 받으나 오버행은 볼록 게이트로 볼록만 — 비볼록 오버행(오버행 detect 볼록 게이트 제거)이 다음 후보.
  - **(c-4) 오버행 boss가 비볼록 솔리드 수용 (실행됨 — c-3 비대칭 절반 해소).** `detect_overhang_contact`(Fuse)의 whole-solid `is_convex` 게이트를 **접촉면 footprint 게이트**(두 접촉면 outer loop 볼록 + 홀 없음)로 교체. 근거: Fuse 재구성은 **로컬** — arc-split은 (볼록) 접촉면만, `resplit_overhang`은 wall을 edge-local(볼록 무관)로, 나머지 면은 verbatim 재방출. 그래서 **비볼록 솔리드(포켓 파인 부품·부울 결과)에 접촉면만 볼록이면 boss 캔틸레버**가 붙는다. n0로 OCCT-정확 확인(포켓 큐브 옆면 오버행 = 1.17). silent-wrong 원천인 비볼록 **접촉 footprint**(arc 오분류)는 게이트가 거절. **Cut/Common은 whole-solid 게이트 유지** — 그 `clip_bwall_inside_a`가 breached 반평면 SH 클립이라 볼록 kept에서만 "inside a"와 일치(비볼록은 `OVERHANG_ARCS` 정직 거절, n0 실증). 신규 `loop_is_convex_2d`·`face_outer_is_convex`. **후속**: ② 비볼록 footprint 오버행(multi-piece 재구성), 비볼록 솔리드 오버행 Cut(clip_bwall_inside_a 일반화).
- **보존 불변식 — 의도는 유효, 메커니즘 서술은 arrangement가 대체함(커토버에서 갱신).** 원문은 *"(b)의 참조-fast-path가 순수성을 유지해야 한다 — 공유 면 seam은 Discovered로 승격하지 않고 Constructed로 남긴다"* 였으나, **(b)는 이미 폐기됐고**(위 (c-1)), 커토버 뒤의 arrangement는 결과 정점을 **평면 삼중항으로 다시 이름 붙이므로 회전 없는 부울 출력이 전부 `Discovered{ThreePlane}`** 다(`Node::Orig`는 생성되지 않음; 잠금 = `an_unrotated_boolean_names_every_vertex_by_its_plane_triple`). ⇒ **불변식의 의도(*"통합이 exactness를 후퇴시키면 안 된다"*)는 지켜진다** — 측정한 축정렬 사례에서 그 정점들의 `tol = 0`으로 `EPS_CONSTRUCTED`(1e-9)보다 오히려 엄격하고, 상자 모서리는 실제로 세 평면의 교점이라 정의도 참이다. **바뀐 것은 exactness가 아니라 표식이다:** `Origin`은 더 이상 *"어느 정점이 컷에서 왔나"* 를 구별해 주지 않으므로 **provenance로 면·정점을 고르는 코드는 성립하지 않는다**(이 사실을 모른 픽스처가 조용히 엉뚱한 면을 골랐다 — dev-log 참조).
- **업계 정합.** 피처 레이어는 pad=tool body+boolean으로 통합돼 있고(SolidWorks/NX/Creo/Fusion), 커널 레이어는 그래도 coincidence 전담 로직을 보유한다 — 상용은 그것을 imprint+tolerance로 두지만, nacre는 **커널 불리언 안**에 두었다(surface Handle 공유 + exact 술어). nacre는 tolerance 대신 exact 술어로 그 자리를 채우는 소수파.

**완전한 M5 오버홀 (회전 지원) — 확정 설계.** 축정렬 전용의 근본 한계를 넘어 회전(비정렬·마름모·각도 스케치)을 지원하는 실사용 M5의 오버홀. **확정 결정:** ① **유리수 치수·각도 표현**(§1 `nacre-scalar`; 유리수는 입력·유리수-순수 파생에만 살고, 무리수·복잡 연산·비트 상한이 닿는 순간 캐시를 f64/dd로 강등 — 정의는 불변이라 export 시 재계산으로 exact 복원). ② **명시 공유(전역 자동병합 폐기 — TNP 충돌):** 면에 스케치=`Surface` Handle 재사용·불리언이 계산한 접촉만 공유; 우연히 같은 좌표는 별개 유지(불리언 시점에 판정). ③ **CIP 판정 정책(declare-0):** f128로도 애매하면 커널이 "같다"고 추측하지 않고 **사용자에게 확인**(site별·기본 `Reject`), "같다"면 Handle 재사용으로 추이성 붕괴를 구조적 차단. ④ **op-log 소유:** 상위 `Document`가 `Vec<Operation>`(진실)+파생 `Model`을 소유하고 sugar(`PadOnFace` 등)를 그대로 기록(피처 트리)·TNP는 "치수 변경=자동 replay, 위상 변화=수동 재지정". ⑤ **다중 솔리드 출력**(OCCT n0 확정: sever는 `Vec<Solid>`·각 manifold). 핵심 메커니즘은 격리 2D/3D 실험으로 선검증 후 이식했다(stage 1~3 대부분 완료 — dev-log).

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

**M6 부호 판정·내외 분류 방식 — 세 마일스톤 중 유일한 미정 (M6 직전 확정).** M5는 방법 확정(indirect predicates), M7은 방법 미정이어도 무방(SSI 실패 시 `Unsupported`로 정직하게 거절하는 게 설계에 내장). 반면 **M6만 "확실히 되는 실용 영역"이라 문서가 약속했는데 정작 어떤 방법으로 판정할지가 비어 있다** — M5의 indirect predicates는 선·평면(다항식)에 특화라 이차곡면에 그대로 안 맞고, M7의 exact ray casting은 메시 경유 하이브리드용이라 M6엔 과하다.

미정인 근본 이유: 이차곡면 교차는 난이도가 갈린다. **평면∩이차곡면**(평면∩실린더=타원, 평면∩구=원, 평면∩원뿔=원뿔곡선)은 **닫힌 형식**이라 쉬운 쪽이고, **이차곡면∩이차곡면**(실린더∩실린더 등, 일반적으로 4차 공간곡선)은 특수 케이스(직교 등)만 닫힌 형식이고 일반은 이미 SSI에 가깝다. 따라서 M6는 "하나의 알고리즘"이 아니라 **"어느 곡면쌍까지를 M6로 긋고, 그 판정을 무엇으로 할지"라는 범위+방법이 얽힌 선긋기**다.

**현재 유력 방향 (확정은 M6 직전 재조사):**
- **평면∩이차곡면 교차점 판정** → indirect predicates를 **부분 확장**하는 것이 유력. 교차 곡선이 닫힌 형식이고 이차곡면도 2차 다항식이라, implicit point를 "이 평면과 이 이차곡면의 교차"로 정의하면 indirect 술어의 다항식이 복잡해질 뿐 **여전히 다항식**이라 M5 indirect가 여기까지 늘어날 여지가 있다.
- **일반 이차곡면쌍** → indirect가 버거워지는 지점. exact ray casting(M7 도구)을 앞당겨 쓰거나, 커버리지 밖으로 두고 `Unsupported`로 거절. 즉 어려운 이차곡면쌍은 사실상 M7 쪽으로 미룬다.

**M6 이차곡면 교차 — 참고 논문 + 라이선스 규율 (M6 직전 사용, 지금 구현 아님).** 이차곡면 교차 위상의 exact 분류 참고처를 미리 걸어둔다(위 두 불릿의 판정 근거).

- **QI 3부작 + 구현 논문 (pencil 분류로 교차 타입을 대수적으로 exact 결정).** Dupont·Lazard·Lazard·Petitjean, "Near-optimal parameterization of the intersection of quadrics" — SoCG 2003 발표, JSC 2008 저널 3부작(Part I 생성 알고리즘 43(3):168–191, Part II pencil 분류 43(3):192–215, Part III 특이 교차 43(3):216–232). C++ 구현 논문은 별도: Lazard·Peñaranda·Petitjean, "Intersecting Quadrics: An Efficient and Exact Implementation", SoCG 2004 / Comp. Geom. 35(1–2):74–99, 2006. pencil 분류로 교차 타입(원·타원·점·두 원뿔 등)을 부동소수점 오차 없이 대수적으로 결정 — 위 "평면∩이차곡면 indirect 확장"·"일반 이차곡면쌍"의 위상 판정 참고처.
- **CGAL 3D Spherical Kernel (2차 대수적수 exact — 다른 계보, "대안 접근"으로만).** de Castro·Cazals·Loriot·Teillaud, Comp. Geom. 42(6–7):536–550, 2009. 곡면 교차점을 2차 대수적수로 exact 처리. 단 우리 indirect predicates(Attene, 정의 기반)와 다른 계보(2차 대수적수 타입)라 M5 indirect의 "확장"이 아니라 "대안 접근"으로만 참고.

**라이선스 규율 (indirect predicates와 동일, 오히려 더 엄격 — §7·§8 M5 라이선스 규율과 정합).**

- **QI 구현**(LORIA/INRIA, gamble.loria.fr/qi): "free for non-commercial use" — MIT/Apache 아님, 비상업 한정. 내부 부품(Uspensky 실근 분리 등)은 또 다른 라이선스. → **소스 열람·차용 금지**, 위 논문(3부작 + 구현 논문)만 참고.
- **CGAL Spherical Kernel**(`Circular_kernel_3`): 패키지 오버뷰에 License: **GPL** 명시(CGAL 이중 라이선스 중 상위 알고리즘은 GPL). GPL은 강한 copyleft라 링크 시 nacre 전체가 GPL 오염 → **소스 열람·링크·차용 절대 금지**, 논문만 참고. (CGAL kernel 기반부는 LGPL이나, 관심 대상 Spherical Kernel은 GPL.)
- **규율 요약:** 두 참고처 모두 **논문만 읽고 clean-room 구현**, 소스는 열람조차 안 함, 완성 후 실행 대조(dev 전용)만 허용 — Attene 2020(LGPL)에 적용한 그 규율 그대로(§8 M5·§7). exact-arithmetic 기반은 M5처럼 `geometry-predicates`(MIT/Apache) 등 자유 라이선스 부품으로 clean-room.

**지금 할 것: 없음**(M6 직전까지 인지만). M6 진입 시 이 논문들을 읽고 라이선스 오염 없이 clean-room으로 구현한다. 아래 "M6 직전 재조사" 시점에 이 메모를 꺼낸다.

**위험도:** 미정이지만 위험하지 않다. M6도 M7과 같은 안전장치(`Unsupported` 거절)가 있어, "정한 방법으로 되는 데까지만 하고 나머지는 정직하게 거절"이 가능하다. 차이는 M7은 "도박이라 열어둠", M6는 "방향은 있고 확정만 M6 직전으로 미룸".

**지금 정하지 않는 이유:** M5에서 indirect predicates를 실제 구현해보면 "이것이 이차곡면으로 얼마나 확장되는지"에 대한 감이 생기고, 그 감이 M6 방법 선택을 정확하게 만든다. 지금 확정하면 "M5 구현 경험 없이 미리 상세 설계"라는 함정. M1~M5 진행 중에는 인지만 해두고, M6 직전에 재조사해 확정한다.

**`nacre-predicates`의 geom 타입 참조 — 순환 의존 확인 (M5 직전).** indirect predicate의 implicit point는 "어느 원시 요소들의 교차인지"를 정의로 보유하는데(§4 `Origin::Discovered`), 그 정의가 `Handle<Surface>`(예: 평면 3장 교차 `[Handle<Surface>; 3]`)를 담으면 문제가 생긴다 — `Surface`는 `nacre-geom`에 있고 `nacre-predicates`는 geom보다 **아래** 층(§1)이라, predicates가 `Handle<Surface>`를 참조하면 geom→predicates→geom 순환이 된다. 이는 §1에서 `Store`/`Handle`을 최하위 `nacre-store`로 내려 푼 것과 **동일한 구조의 문제**다. M5 구현 직전에 정한다: (a) predicates가 `Handle<Surface>`를 직접 담지 않고 **좌표·평면 방정식(계수)만 값으로 받는다**(가장 단순 — predicates가 geom을 전혀 모름), 또는 (b) implicit point 정의를 geom 쪽(또는 store 같은 하위 공용 층)에 두고 predicates는 그 위에서 술어만 제공, 또는 (c) store 패턴처럼 공용 최하위 타입으로 분리. 현재 유력안은 (a) — indirect orient3d는 결국 평면 계수들의 다항식 부호이므로, Handle이 아니라 평면 방정식 계수를 넘기면 predicates가 순수 수치 계층으로 남아 순환이 원천 차단된다. **→ 옵션 (a) 확정(M5-prep).** `nacre-predicates`는 평면 계수·좌표 `[f64;N]`만 받는 순수 수치층(커널 타입 무의존, standalone 분리 가능). implicit point의 Handle 기반 정의(`VertexDef::ThreePlane([Handle<Surface>;3])`)는 위상 계층(§4 `Origin::Discovered`)에 두고, 부호 판정 시 geom이 Handle→계수를 뽑아 predicates에 넘긴다. 순환 원천 차단.

**곡면 내/외 판정 후보 (M6 직전 실측 결정).** 셋 다 M6~M7 "내/외 분류"의 후보이며 SSI 해법이 아니다(§8 M7 한계). 1순위 **exact ray casting**(§8 M7 — 우리 계보 정합, Cherchi 2022 검증, 메시 경유라 곡면 직접 판정 우회). 2순위 **GWN**(generalized winding number, Jacobson 2013 — 메시/point cloud in-out은 10년+ 검증된 성숙 기법, libigl·Axom[BSD] 구현 존재; watertight 무관 강건성이 강점이나 우리 always-closed에선 덜 필요, 느림·경계 round-off, trimmed NURBS 정확 GWN 확장은 최신[Spainhour 2024~26, 검증 진행 중] — "메시 경유 없이 곡면에서 직접 판정하고 싶어질 때"의 대안으로 보류). 참고 **graph cuts**(Diazzi/Attene 2021 — 일부 모호한 자기교차 케이스 우수), **EMBER winding number vector**(Trettner 2022). **라이선스**: 채택 전 확인, Cherchi/Attene 계열 LGPL 주의 — indirect predicates와 동일 규율(논문·MIT/Apache 소스만 참고, LGPL 소스 열람 금지, 실행 대조만 허용).

**CIP — Certified Indirect Predicates (회전 시의 부호 판정 층. 구현됨 — 오버홀 stage 1~3 완료, 수식 ②는 아래 :566에서 소진.)**

**왜 필요한가.** 축정렬 판정은 좌표가 유리수라 exact다. 그러나 **회전이 들어오면 좌표가 무리수가 된다** — 유리수 각도라도 cos·sin은 무리수이고(Niven), 임의 각도는 초월수다. exact 산술로 표현할 수 없으므로 회전 좌표는 f64 근사일 수밖에 없다. **회전은 오버홀(stage 1~3)로 지원된다** — 모델링 변환·각도 스케치가 회전 좌표를 만들고, 그때 **판정만은 조용히 틀리지 않게** 지키는 것이 CIP다.

**핵심 명제.** indirect predicates(Attene 2020)는 **선형 요소(선·평면)의 교차**, 즉 다항식에만 성립한다. 회전은 다항식이 아니므로 **indirect predicates를 회전에 적용하는 방법은 없다**(문헌에도 없다 — 수학적 한계). 실제 연구·구현이 하는 일은 하나뿐이다: **회전 좌표는 f64로 근사하고, 그 위의 판정을 exact로 유지한다.** CIP는 그 "위의 판정"에 **입력 근사(tol)까지 반영**하는 층이다. **세 계보의 하이브리드다** — Attene 2020 간접술어(implicit point 구조)·Shewchuk 1997/CGAL 필터(f64 필터→고정밀 상승)·Guibas 1989 epsilon-geometry(값+tol·sound 판정 or 기권). 유리수부는 완전 exact, 초월(회전)부는 sound 오차 한계로 부호를 **인증하거나 정직 기권**(declare-0→ask). **포기한 것은 초월수까지의 이론적 *완비성*이지 *건전성*이 아니다** — 틀린 부호는 결코 내지 않는다(no-silent-wrong). 새 패러다임이 아니라 세 계보의 특정 실현이다.

**★ 대체가 아니라 흡수다.** CIP는 indirect predicates를 지우지 않는다. 점을 "정의 + 누적 tol"로 들고, **tol = 0이면 지금의 indirect predicates 그대로(exact)**, **tol > 0이면 필터 + 고정밀 폴백**이다. 술어의 계산 구조(3-평면 정의, 행렬식 부호, 좌표 무독)는 **그대로 쓰고**, "평면 계수가 정확하다"는 가정만 "계수에 tol이 있다"로 넓힌다. 회전된 평면은 **계수에 tol이 붙은 평면**일 뿐이다.

**구성 요소 (개념 확정).**

1. **모든 점 = 정의 + 누적 tol.** 평면 교차점(정의 = triple, tol = 0)과 회전 점(정의 = 변환 이력, tol > 0)이 **같은 틀**로 통일된다. 특수 경로 없음.
2. **tol의 두 축.** 자체 tol = **이산화 오차**(값을 f64/f128로 표현할 때; f64 ≈ 1e-16×크기, f128 ≈ 1e-34×크기) + **연산 오차**(fraction 덧셈 = 0, f64 덧셈 = 크기 차로 비트 손실 발생, 곱셈·회전 = 발생). 누적 tol = 앞 tol + 자체 tol.
3. **★ 회전 tol은 거리에 곱해져 전파된다 — 단순 덧셈이 아니다.** 회전은 *방향*을 정하므로, 각도오차 `da`가 **그 회전 이후 판정점까지의 거리**에 곱해져 위치 오차가 된다. 따라서
   ```
   누적 tol = Σ(이산화·직선 연산 오차)  +  Σ_각회전i ( da_i × 거리(회전 i → 판정점) )
   ```
   앞선 회전일수록 뒤 이동이 많아 기여가 크다. **분리 계산은 삼각부등식 |a+b| ≤ |a|+|b|로 정당하다**(최악 상한이므로 실제보다 작아지지 않는다 — 터지지 않는다). **★ 교정(H4-soundness, 아래 미확정①)**: 위 `da × 거리`(접선형)는 **각도 불확실성** 항만 맞고, v1의 **실현 반올림(이산화) 새-오차는 접선형이 아니라 좌표혼합 `(|x|+|y|)·da`**다(독립 cos/sin 반올림의 방사 성분; 접선형은 좌표 0 근처서 과소평가—실측 반증). 유리수 각의 v1은 `da_각도=0`이라 위 둘째 항이 0이고 새 오차는 좌표혼합 항이 담당한다.
4. **★ 캐시 가능한 부분과 아닌 부분이 갈린다.** 이산화·직선 tol은 **적용점 무관**(점의 고정 속성) → 점에 값으로 캐시. **회전 tol은 적용점 의존**(같은 회전도 판정점이 멀면 tol이 커진다) → **캐시 불가, 판정 시 계산**. 그러려면 점의 **부모·조상**을 알아야 하고 조상 중 회전이 어디에 몇 번인지 알아야 한다. → **회전 이력 트리**(노드 = 회전의 각도오차·축·위치 + 부모 링크; 이동·교차는 담지 않는다 — tol에 기여하지 않으므로. 여러 점이 공통 조상을 공유하는 forest). 판정 시 노드에서 상위로 순회하며 (각도오차 × 판정점까지 거리)를 합산한다.
5. **★ 방향별 tol(x·y·z)이 필요하다.** 행렬식에서 각 방향의 오차가 **서로 다른 계수로 증폭**된다(2D orient에서 `a`의 x-tol은 `(by−cy)`와, y-tol은 `(bx−cx)`와 곱해진다). 하나로 뭉치면 정확히 증폭할 수 없다. **CGAL Lazy_kernel이 좌표를 구간(interval)으로 드는 이유가 이것**이고, 우리는 구간 산술 대신 **방향별 tol 값 + 미리 유도한 오차 한계 공식**을 쓴다(술어가 소수·고정이므로 유도가 가능하다 — CGAL은 범용이라 술어가 수백 개라 유도가 불가능해 구간을 택했다. 우리는 특화 커널이라 Shewchuk/Attene 계열의 "오차 한계 미리"가 더 빠르다). **→ 확정(미확정①): 스칼라 tol 탈락, 모든 tol을 xyz 벡터로. 회전 각도오차 환산 공식은 위 미확정① 참조.**
6. **직선 구간은 fraction 강체로 묶는다.** 연속된 직선 이동을 유리수로 먼저 합산하면 그 구간의 **연산 오차 = 0**이고 이산화는 **회전과 만나는 지점에서 1회**뿐이다. 순차로 f64 덧셈하면 이동마다 이산화가 붙는다((a+b)+c는 이산화 2회, a+(b+c)는 1회). 회전 결과는 무리수라 fraction으로 못 담으므로 **회전 경계의 f64 덧셈 오차는 피할 수 없다** — 피할 수 있는 것만 피한다. **(일반화 — 기준은 "직선/회전"이 아니라 "유리수/무리수"다. 직선이라서 tol 0이 아니라 유리수라서 tol 0이며, 강체 묶기는 직선뿐 아니라 **같은 축의 연속 유리수-각도 회전**에도 적용된다[각도를 유리수로 합산해 누적 각을 1회만 실현 — 2D는 항상 한 축, 3D는 같은 축만]. 유도된 무리수[√ 거리·구속 솔버 해·3D 축 변경 합성]가 섞이면 그 구간은 fraction으로 못 묶는다. tol 공식의 거리·각도 항은 지우지 말고 0으로 둔다. **★ 번들링은 필수(H4-amplification 실측): un-bundled 증분 실현은 전파 `|R|`의 행합 `|cos|+|sin|≥1`을 매 스텝 곱해 tol 바운드가 지수 폭발[실측 30스텝에 실제 오차의 ~20만 배; 실제 오차는 R 노름보존이라 평평]. 누적각 1회 실현이 차단. 폭발해도 sound 최악-보장이라 조용히 안 틀리고 상승/거절로만 간다[항목 8].)**
7. **판정 = 필터 + 정밀도 상승 + lazy 캐싱.** 점들의 tol로 이번 판정의 오차 한계를 계산(변 길이로 증폭) → `|행렬식| > 오차 한계`면 f64로 확정(대부분) → 애매하면 정밀도 상승. 관련 점이 **tol = 0(평면 교차점)뿐이면 exact 폴백**(현행 indirect predicates), **회전 점이 끼면 f128**. 한 번 고정밀 계산한 **값**은 캐시해 재사용하고(판정 결과가 아니라 점의 값 — 여러 판정에서 재사용된다), 재계산 시 **회전 지점만** 다시 계산하고 직선 구간은 fraction을 그 정밀도로 이산화해 잇는다. **판정에 필요한 점만 정밀화한다** — 중간 경유점은 정의로만 남긴다.
8. **오차 한계는 최악(선형 합)이다.** 확률 전파(RSS)는 실제 오차에 가깝고 100배 작지만 **보장이 아니다** — tol은 "이보다 클 수 없다"는 보장이어야 하고, 확률 tol은 드물게 부호를 뒤집어 **조용히 틀린다**. 커짐은 다른 방법으로 관리한다(직선 tol = 0이라 N은 회전 수뿐, 애매하면 f128).
9. **★ 필터는 동적으로 통일한다.** 지금 (5b-0)의 필터는 **정적**(오차 한계가 상수 `ε_D ≈ 5u`·`ε_M ≈ 19u`)이고, 입력 tol이 없으므로 그것이 옳다. 회전이 들어오면 **동적으로 통일**한다 — 회전 tol은 경로마다 편차가 커서 하나의 정적 상수로 잡으면 느슨해 무용하거나(과대) 위험하다(과소 → 터진다). 그리고 **각도를 쓰는 스케치도 회전이므로 무회전 케이스는 극소수**이고, **동적 필터는 tol = 0인 점을 자동으로 정적급 타이트하게 처리**하므로 정적 분기를 남기는 것이 순수 오버헤드다.

**★ 근본 한계 (정직하게).** 회전이 들어가면 **서로 다른 경로로 같은 위치에 도달해도 f64가 다를 수 있고, exact 술어로도 그 다름을 그대로 반영한다** — 술어는 "주어진 입력에 대해" 정확할 뿐 입력 자체의 근사를 고치지 못한다. f128로도 **보장은 불가**하다(초월수의 상등은 유한 정밀도로 결정 불가). 다만 **실무적으로는 충분하다**: 상용 CAD의 uncertainty가 ~1e-9인데 f64는 1000mm 점을 10000번 회전해도 ~1e-10, f128이면 ~1e-28이다. **이론적 완벽함은 포기하고 실무적 충분함을 tol로 보장한다** — 그리고 그 tol을 **측정된 값**으로 들고 다니는 것이 tolerance-fudge와 다른 점이다. (f128로도 애매한 **잔여 케이스**는 커널이 0으로 추측하지 않고 **사용자에게 확인받는다**(§9 회전 오버홀의 CIP declare-0 정책). 즉 여기 "tol 보장"은 f128이 부호를 확정하는 대다수 경로를 말하고, 확정 못 하는 잔여는 ask-user로 넘긴다.)

**마일스톤별 적용.**
- **M5·M6:** 판정이 다항식이므로 CIP가 그대로 얹힌다. **M6에서 CIP의 가능 여부 = 간접 술어의 가능 여부**다(CIP는 술어 위의 층이므로) — 평면∩이차곡면 간접 술어가 서면 CIP도 선다. 단 술어가 고차가 되면 **오차 한계 공식을 그 술어에 대해 재유도**해야 하고, 고차라 필터 성공률이 M5만큼 좋을지는 **측정 대상**이다.
- **M7:** 메시 조합 판정의 필터 + 뉴턴 스냅백 점의 tol 추적·f128 재수렴. **SSI(교차를 찾는 것)는 CIP 밖**이다(위상 존재 문제이지 tol 문제가 아니다) — 아래 "CIP와의 관계 정리" 참조.

**아직 정하지 않은 것 (구체화 대상 — 도입 직전에 정한다).**
1. **~~회전 연산의 자체 tol을 무엇으로 표현할지~~ → 확정(H4-soundness 실측 교정).** **모든 tol(회전·직선·이산화)을 xyz 방향별 벡터로 표현·누적**한다. 스칼라(구) tol은 탈락 — 정보를 버린다(z축 둘레 회전이면 z 방향 오차 0인데 스칼라는 모든 방향에 실어 "z=5 평면 위인가?"를 불필요하게 애매 판정→상승 유발). **★ 회전 tol은 세 항의 합이다** (초기 "접선형 새-오차" 안은 H4-soundness에서 반증됨 — 아래 (a)): (a) **새 오차 = 실현 반올림(좌표 혼합)** + (b) **기존 tol 방향 회전** + (c) **각도 불확실성(접선, v1엔 0)**. 회전 시 tol 갱신 =
```
새 tol = |R| · (기존 tol)            [b: 기존 tol의 방향을 회전 — 성분별 절댓값 행렬 |R| × 기존 tol]
       + (|x| + |y|) · da_이산화      [a: 새 오차 — 실현 반올림, 좌표 혼합]      ★ H4 교정
       + da_각도 × |n × (P − 축점)|   [c: 각도 불확실성, 접선 — v1엔 0]
```
**★ (a) 교정(H4-soundness).** 초기 안은 새 오차를 **접선형 `da × 모멘트암`**으로 뒀으나 **반증됐다** — cos/sin을 독립 반올림하면 순수 회전이 아니라 **방사 성분**이 생겨(`err_x = |x·δc − y·δs|`), 접선형은 좌표 하나가 0 근처일 때 그 성분을 **과소평가**한다(실측: 10000 랜덤 중 220회 상한 붕괴). 성분별 새-오차 tol은 **두 좌표를 합친 `(|x|+|y|)·da_이산화`**여야 sound(0으로 안 꺼짐; 실측 0회 붕괴·최악 ~3.5× 이내 보수). `da_이산화`=실현 반올림 한계(f64 ~16 ulp; 상승 시 ~2^-P). (b) `|R|·(기존 tol)`은 **필수** — x방향 tol이 z축 90° 회전 후 y로 옮겨가는데 안 반영하면 상한이 깨진다. (c) **각도 불확실성**(접선)은 유리수 각의 v1엔 0이나(그리고 (a) 여유가 f64 각도-반올림까지 흡수), 유도 치수(구속 솔버·비유리수 각)에 살아나므로 **항은 지우지 말고 0으로 둔다**. `n`=회전축 단위벡터, `P`=판정점. 직선 이동은 (a)(c)를 안 만들지만 앞선 tol을 싣고 간다. 판정 시 각 방향 tol이 행렬식에서 자기 계수로 증폭돼 오차 한계를 이룬다. (CGAL Lazy_kernel의 interval과 같은 이유.)
2. **~~입력 tol → 행렬식 오차 한계 공식~~ → orient3d는 확정·이식(단계 2a-ii): `frame3::det3_bound`**(6개 signed triple-product 구간 반경 `prod_err` + `16ε·mag` f64 반올림; `orient3d_judge`가 필터→astro-float 상승→declare-0로 소비). exact3d H-a 검증(위반0·tightness 0.12; 프로덕션 재현 0.065·피벗 포함). `plane_side`=explicit orient3d 쌍둥이(같은 공식). **간접 orient3d**(3평면 implicit point·평면 계수)도 확정·이식(단계 2c-i `indirect_orient3d_judge`; exact3d H-b/H-c 검증·프로덕션 재현 위반0). ①이 확정돼 **입력이 명확**하다 — 점 tol = xyz 방향별 벡터. "방향별 점 tol이 행렬식에서 각자 자기 계수로 증폭돼 오차 한계를 이루는"(항목 5) 공식을 술어별로 유도(현행 정적 상수는 tol = 0 가정). **H4-soundness 진행 상황**: **(a) 새-오차·(b) 전파 항 모두 실측 검증 완료.** (a) 접선형 반증→좌표혼합(10000 랜덤 0회 붕괴). (b) `|R|·기존 tol`: 증분 회전 체인 5000회 0회 붕괴(sound·~3× 타이트); 전파의 **필요성은 비대칭 tol**에서 발현(x-tol을 90° 회전 → 오차가 y로 이동, 전파 없으면 tol_y=0으로 과소예측 — 대수적 확인; 비대칭 tol은 v1-후 √거리·구속 솔버에서). 아래 ③(평면 계수 tol 전파)도 같은 유도의 일부. **남은 H4**: consistency(경로-독립 결정성)·amplification(증폭 경계)·속도.
3. **~~평면 계수의 tol 구조~~ → sqrt 섭동은 (5d)-1이 해결, 나머지는 ②로 흡수.** `Plane::through_points`의 정규화(sqrt) 섭동은 미확정이 아니다 — (5d)-1이 Plane에 비정규화 `raw`를 함께 저장해 판정이 쓰는 `coefficients()`가 exact(정의 정점이 정확히 0)이고 `normal()`(sqrt)은 크기 소비자에게만 간다 → **판정 경로에 sqrt 섭동 없음**. 남은 **"회전된 평면의 계수 tol"**은 별개 미확정이 아니라 **①의 따름정리**다: 평면은 세 점으로 정의되고, 회전되면 그 점들이 ①의 방향별 tol을 가지며, 계수는 그 점들의 뺄셈·외적이라 **점 tol이 계수 tol로 전파**된다. ②의 오차 한계 유도의 일부로 함께 다룬다.
4. **~~3D에서 "거리"의 정확한 정의~~ → ①에서 확정.** 새-오차 (a)는 **좌표 크기 `|x|+|y|`**로 실측 확정(H4). 각도-불확실성 (c)의 접선 크기 = `da_각도 × 회전축까지 수직거리(모멘트 암)`(v1엔 0).
5. **~~방향별 tol의 자료구조~~ → 확정: `[f64;3]` xyz 벡터**(단계 2a). `nacre-cip::kernel::frame3::Pt3 { base:[Rat;3], chain, coord:[f64;3], tol:[f64;3] }` — 2D `frame::Pt2`의 3D 아날로그. (**회전 시 tol 변환**은 ①의 `|R|·기존 tol`로 흡수.) 임의 유리수 피벗까지 exact3d H-f로 검증(피벗 산술 tol 항 추가).
6. **~~고정밀 층의 선택~~ → 확정: `astro-float`(순수 Rust 임의정밀). twofloat(double-double)는 H1.5에서 탈락.** H1.5 실측: twofloat의 π 상수는 정확하나 삼각함수가 "preliminary"라 **영점 근처 cos 오차 ~1.8e-16(f64 수준)**으로 정확도 게이트(~1e-30) 실패 — 부호 판정이 일어나는 near-degenerate가 곧 영점 근처라 치명적. `astro-float`로 교체(160비트에서 cos/sin 오차 **~1e-58**로 통과). **정밀도가 dial 가능**이라 double-double의 ~1e-32 천장이 사라진다(H4-amplification이 tol을 키워도 정밀도를 올리면 됨 → "quad-double 상승 or 정직 거절" 딜레마가 "정밀도 dial"로 단순화). 대가는 속도 **~120µs/call**(상승 경로라 드물고, 각도당 실현을 캐시해 상각; 정밀도를 낮추면 빨라짐) — H4가 속도·상승빈도 실측. **★ 이하 이 문서의 "double-double"·"f128"·"고정밀 상승"은 이 상승 층을 가리키는 일반명이며, 실제 구현은 astro-float다.**
7. **회전 이력 트리 → 확정(단계 1c): `Model.rotations: Store<Rotation>` + `parent` 링크 forest.** 남은 **캐시 정책**(hp 실현값 수명·축출)은 단계 2a에서 미실행(Pt3 매번 fresh 계산) → **단계 2b 이후로 defer**(성능 최적화, soundness 무관).

**★ 정리 — ②(오차 한계)가 완전 소진됐다: orient3d/plane_side 직접(2a-ii `det3_bound`)·간접(2c-i `indirect_orient3d_judge`)·cmp_coord(cmp-i `indirect_cmp_coord_judge`) 모두 확정·이식(프로덕션 `nacre-cip::kernel::frame3`).** 남은 §CIP는 ⑦-캐시(성능·defer)뿐. ①·④·⑤·⑥ 확정, ⑦-트리는 1c 완료, ③은 sqrt 부분 (5d)-1이 해결·나머지 ②로 흡수. §CIP의 마지막 수학 ②는 실험 H-a(직접 orient3d)/H-b(평면 계수)/H-c(간접 orient3d)가 검증했고, cmp_coord는 exact3d에 없던 새 수학이라 frame3에서 직접 **H-g**(두 코퍼스·wrong-sign 0)로 검증했다(Cramer 기계는 2c-i 이식분 재사용·최종 부호 결합만 신규). 세 술어 모두 프로덕션에 들어왔다 — **판정층 완성**.

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
4. 리미트 도달 → 뉴턴법으로 넘김: 해 수렴=교차 존재, 발산=미교차. 뉴턴 스냅백 점은 tol 있는 점이므로 **CIP에 통합**(수렴 잔차=tol, 애매하면 f128 재수렴).
5. 리미트에서도 확신 불가(작은 loop는 "다 찾았다"의 수학적 보장이 근본적으로 불가) → **`Unsupported` 정직 거부.** 조용히 틀리기보다 거부(§8·M7 철학).

**출발점 — M7 진입 시 이것부터 읽는다:** Li·Yang·Jia, "Advances and challenges in surface–surface intersection computation — An overview", Computer-Aided Design 193:104039, 2026. SSI 분야 전체 최신 개관이라 개별 논문 여러 개보다 이 리뷰가 최적 출발점. 그 시점의 최신을 반영해 위 세 계보(거리/법선/winding number)를 재비교 후 채택. (관련: Li·Jia·Chen, "Fast Determination and Computation of Self-intersections for NURBS Surfaces", ACM TOG 44(2), 2025 — ④ 자기교차 판정·거부용.)

**CIP와의 관계 정리 (혼동 방지).** CIP(Certified Indirect Predicates)는 "다항식 판정 + 회전 tol 필터"라 **부호 판정** 층이다. M7에서 CIP가 닿는 곳은 (a) 메시 조합 판정(내/외 분류)의 필터, (b) 뉴턴 스냅백 점의 tol 추적·f128 상승 — 둘 다 **정밀화·판정**이다. M7의 도박인 **SSI(교차를 찾는 것)** 자체는 tol 문제가 아니라 위상 존재 문제이므로 CIP 밖이다. "찾은 것을 정밀하게"(CIP·스냅백)와 "못 찾은 것을 찾기"(SSI)는 다른 일. SSI 성공 후라야 스냅백·판정이 의미 있고, SSI 실패 시 `Unsupported`.

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
