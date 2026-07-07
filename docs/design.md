# nacre — 하이브리드 CAD 커널 1차 설계 스케치

정확한 기하를 진실로 보관하는 고전 b-rep의 골격 위에, Fornjot에서 검증된 위생 규율(append-only 단일 참조, 위상과 근사의 동시 구축)을 얹은 설계다. 목표는 정밀 기계 CAD(STEP 입출력 포함)이며, 개인 + AI 협업 개발을 전제로 실수 여지를 줄이는 인프라를 1일차부터 포함한다.

## 0. 설계 원칙

**범위 (비목표 명시).** nacre는 순수 기하 커널이다: 기하·위상·수치 봉합(tolerance)·연산·tessellation·검증까지만 다룬다. GD&T/PMI(제작 공차), 스타일·색상·레이어·visibility, 제작자·승인·날짜, 제품 구조·조립·리비전은 커널 범위 밖이며, 커널을 사용하는 응용 프로그램의 영역이다. 커널은 이들을 위해 단 하나만 제공한다 — 영구 유효한 Handle(append-only의 부산물). 응용은 `HashMap<Handle<Face>, 응용데이터>` 사이드카로 무엇이든 매달 수 있다. nacre-step은 STEP의 형상 서브셋만 커널로 번역하고, 비형상 엔티티는 해석 없이 무손실 패스스루로 보존해 라운드트립을 지킨다.

이 문서의 모든 결정은 다섯 가지 원칙에서 나온다. 첫째, **정확 기하가 진실이다** — 평면·원통·NURBS는 해석적 형태로 영구 보관하고, 메시는 파생물이다. 둘째, **모든 객체는 append-only 저장소에 딱 한 번 존재하고 Handle로만 참조된다** — 동일성 질문을 좌표 비교(기하)가 아니라 인덱스 비교(명목)로 바꾼다. 셋째, **연산 이력을 보존한다** — 모델은 연산 로그의 재생 결과다. 단, 보장 범위를 정확히 한다: replay는 동일 로그·동일 파라미터에서 동일 모델을 보장하고(undo/redo의 기반), tolerance 변경 재계산은 위상·기하를 불변으로 둔 채 tessellation만 재생성한다(§5). 로그 중간의 파라미터를 수정하는 파라메트릭 편집은 v1 비목표다 — 연산이 원시 Handle을 참조하는 한 상류 수정이 하류 Handle 번호를 밀어내기 때문(topological naming 문제). 진화 경로는 §6에 기록. 넷째, **tessellation은 일회용이 아니라 출처 태그가 달린 1급 부산물이다** — 모든 삼각형·근사점이 자기가 어느 정확한 면/엣지에서 왔는지 안다. 다섯째, **tolerance는 "발견된" 교차에만 존재한다** — 구성 시점에 동일성을 아는 요소에는 tolerance 개념 자체가 없으며, 이 구분을 타입 시스템에 새긴다.

## 1. 크레이트 구조

의존 방향은 아래에서 위로만 흐른다. 순환 의존 금지.

```
nacre/                    # 워크스페이스. 최상위 `nacre` 크레이트는 파사드(재수출 전용)
├── nacre-store      # typed-index 인프라: Store<T>/Handle<T> (기하·위상 무지의 순수 저장소)
├── nacre-math       # 벡터·행렬·변환. nalgebra 래핑 or 자체 (Point<D>, Vector<D>, Transform)
├── nacre-geom       # 정확 기하: Surface, Curve, 평가·미분·국소 교차(SSI relaxation)
├── nacre-topo       # Vertex/Edge/Face/Shell/Solid, half-edge, Model 집계
├── nacre-tess       # 출처 태그 tessellation: TessVertex, TessTriangle, 증분 갱신
├── nacre-ops        # 연산: sketch, extrude, revolve, imprint, boolean(자체 — 커버리지 사다리)
├── nacre-validate   # 불변식 검사: 오일러-푸앵카레, watertight, 방향성, 참조 무결성
├── nacre-step       # STEP 입출력 (기존 코드젠 라이브러리 연결 지점)
├── nacre-viewer     # wgpu 뷰어 (nacre-tess 출력을 그대로 소비)
├── nacre-oracle     # [dev] OCCT 비교 하네스 (out-of-process 헬퍼 경유), proptest 전략
└── tools/occt-helper/  # 워크스페이스 밖 헬퍼: brew OCCT(1순위) 또는 uv+OCP(폴백) — §7
                        #   OCCT는 오라클 전용 — 제품 경로에 위임 없음 (§6, §8)
```

`nacre-geom`과 `nacre-topo`가 서로를 모르게 하는 것이 중요하다. 기하는 위상을 모르고(순수 수학), 위상은 기하를 Handle로만 참조한다. robustness가 첨예한 코드(교차·분류)는 전부 `nacre-geom::intersect` 한 모듈에 격리한다. 사용자는 파사드 크레이트 `nacre` 하나만 의존하며, 인터랙티브 스크립트 앱 등은 이 워크스페이스 밖의 별도 프로젝트로 둔다.

`Store`/`Handle`은 **최하위 `nacre-store`에 둔다.** geom도 Handle을 쓰기 때문이다 — `Curve::Intersection`(§3)이 `Handle<Surface>`를 담으므로, Handle이 topo에 있으면 geom→topo→geom 순환 의존이 된다. typed-index 저장소는 기하·위상을 전혀 모르는 순수 인프라이므로 두 층보다 아래에 격리하고, 위의 모든 크레이트가 자유롭게 참조한다. (라이선스는 MIT/Apache-2.0 듀얼 — Manifold(Apache-2.0) 알고리즘 차용과 호환.)

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
    Arc(Arc),
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

정밀도 분업 원칙: **판정(부호)에는 적응 정밀 술어**(`geometry-predicates` — 쉬운 케이스 f64, 아슬아슬할 때만 확장, 부호는 항상 정확), **반복 구성(좌표)에는 고정폭 확장 부동소수점**(double-double). 임의 정밀 유리수(rug/malachite)는 술어에는 완벽하지만 Newton 반복에 넣으면 비트 길이가 반복마다 폭발하므로 구성에는 쓰지 않는다. 이 투자는 평균이 아니라 꼬리를 산다 — 호출 빈도가 낮은 악조건 케이스만 정확히 개선되고, tolerance가 "상수"가 아니라 "실측 보증값"이 된다.

STEP 연결: 이 variant는 `surface_curve`/`intersection_curve`(3D curve + 곡면별 pcurve + master 지정)와 1:1로 대응한다. pcurve가 필요해지는 시점(트리밍 구현 시)에 `pcurves: [NurbsCurve2d; 2]`를 같은 variant에 추가한다. 기존 stepio 코드젠에서 이 엔티티 묶음의 보존·복원을 우선 검증 대상으로 삼는다.

## 4. 위상 층 (`nacre-topo`) — tolerance를 타입에 새기기

이 설계에서 가장 의견이 들어간 결정. 정점·엣지의 "출신"을 enum으로 구분한다.

```rust
pub enum Origin {
    /// 구성 시점에 정의됨 (스케치 점, 스윕 결과 등).
    /// 동일성은 Handle로 완결 — tolerance 개념이 없다.
    Constructed,
    /// 교차 계산으로 발견됨. 국소 tolerance를 가진다.
    /// tol은 임의 상수가 아니라 relax_to_intersection이 반환한 실측 달성 정확도
    /// (+ 정점 봉합 시 곡선 간 불일치 반경)에서 온다.
    Discovered { tol: f64 },
}

pub struct Vertex {
    pub point: Point<3>,    // 인라인 — 점 저장소 없음 (§2)
    pub origin: Origin,
}

pub struct Edge {
    pub curve: Handle<Curve>,
    /// None = 닫힌 엣지 (원 전체 등 — 실린더 윗면 경계가 대표 사례, M3에서 필요).
    /// 닫힌 엣지는 끝점이 없으며, seam이 필요한 곡면 파라미터화는 tess 층에서 처리한다.
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
    Sketch { plane: SketchPlane, profile: Profile2d },
    Extrude { face: Handle<Face>, dir: Vector<3>, dist: f64 },
    Revolve { face: Handle<Face>, axis: Axis, angle: f64 },
    ImprintSketch { face: Handle<Face>, profile: Profile2d },
    PadOnFace { face: Handle<Face>, profile: Profile2d, dist: f64 },
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
```

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

속성 기반 테스트(proptest): 랜덤 유효 연산열을 생성해 (1) validate 통과, (2) replay 멱등성 — 같은 로그·같은 cfg → 동일 모델, (3) 변환 불변량 — 강체변환 후 부피·면적 보존, (4) 불리언 대수 — `A ∪ A = A`, `A ∩ ∅ = ∅`, `vol(A∪B) + vol(A∩B) = vol(A) + vol(B)`.

OCCT 오라클(`nacre-oracle`): 같은 연산열을 자체 커널과 OCCT에 병렬 실행하고 부피·면적·바운딩박스·(가능하면) 면 개수를 diff. 허용 편차를 넘으면 실패한 연산열을 최소화(shrink)해서 리포트. 정답지를 든 채 개발하는 장치이며, AI가 생성한 코드의 "그럴듯하지만 틀림"을 잡는 주 방어선.

**OCCT는 오라클 전용이며, 연동은 out-of-process다.** OCCT를 링크하지 않고, `tools/occt-helper/`의 헬퍼 프로세스와 STEP 파일로 주고받는다(수송 계층 = nacre-step). 이유: C++ 빌드가 Rust 워크스페이스에서 완전히 사라지고, OCCT가 악조건 입력에서 크래시해도 커널 프로세스가 아니라 헬퍼만 죽는다(크래시 격리 — 오라클처럼 수천 회 돌리는 용도엔 필수). 오라클은 채점(부피·면적·위상 수 비교)만 하므로 history 손실·호출 오버헤드가 무관하다.

헬퍼 구현 우선순위(macOS 기준): **1순위 — Homebrew** (`brew install opencascade`): DRAWEXE Tcl 스크립트(readstep → bfuse/bcut/bcommon → writestep, 코드 0줄) 또는 brew 라이브러리에 링크하는 얇은 C++ 헬퍼. **폴백 — uv 관리 Python 환경 + OCP 휠**(사전 컴파일이라 C++ 툴체인 불요, uv가 파이썬 버전 고정까지 해결). 어느 구현이든 프로토콜은 동일하게 고정한다: `helper <fuse|cut|common> <a.step> <b.step> <out.step>`, exit code 0/1/2 = 성공/기하 실패/크래시, stdout으로 진단 JSON(부피·면적·면 개수 — diff 비교값을 헬퍼 쪽에서 계산해 전달). 프로토콜만 지키면 구현을 갈아타도 nacre 쪽 코드는 무변경.

STEP 라운드트립: 골든 파일 셋에 대해 import → export → import 후 구조 비교. `Intersection` 곡선의 이중 표현 보존 여부를 중점 검증.

## 8. 마일스톤 사다리

각 단계는 "동작하는 것"을 남기고 끝난다. 3층에 막혀 전체가 멈추는 구조를 피하는 배치다.

**M1 — 뼈대.** nacre-math, Store/Handle, Plane/Line만으로 정육면체를 손으로 조립. validate 1차 구현. 시각 확인은 2단계로: 먼저 Tessellation을 OBJ/STL로 덤프하는 함수(기존 뷰어 — MeshLab, f3d 등 — 로 확인, 반나절짜리)로 시작하고, 그다음 최소 wgpu 뷰어(창 하나 + 삼각형 렌더 + 궤도 카메라 + 와이어프레임 토글, 이 이상 금지). 뷰어는 커널을 비추는 거울이지 제품이 아니다 — M1에서 렌더링 품질에 시간을 쓰지 않는다. (전부 1층. AI 가속 최대 구간.)

**M2 — 스케치와 케이스 A.** 2D 프로파일(선분만) → extrude. 연산 로그와 replay. 평면 솔리드에 대한 STEP 내보내기 — 기존 stepio 연결. 오일러 연산 정리.

**M3 — 곡선 기하.** Arc, Cylinder, NurbsCurve/Surface 평가(The NURBS Book 기준 구현 + 수치 미분 대조 테스트). tess 출처 태그 완성, tolerance 재계산 데모(같은 모델, tol 3단). proptest 도입.

**M4 — 면 위 작업.** ImprintSketch, PadOnFace — "만나는 자리를 아는" 연산의 완성. 여기까지 모델 전체가 `Constructed`. nacre-oracle 가동(OCCT와 부피 diff). M3~M4에 걸쳐 뷰어를 디버그 도구로 성장시킨다: 면 클릭 → Handle·Origin 표시, 법선 화살표, 엣지 polyline·tolerance 공 시각화, validate 위반 하이라이트, 연산 로그 스텝별 재생(append-only라 "N번째까지 replay"가 공짜). 예쁜 렌더가 아니라 커널 내부가 보이는 기능만 — 이것들이 M6 불리언 디버깅의 생명줄이 된다.

**M5 — 자체 불리언 1단: 다면체.** `PolyhedralBoolean` — 모든 면이 평면인 솔리드 간 fuse/cut/common을 자체 구현한다. 평면-평면 교차는 닫힌 형식의 직선(SSI 행진·Newton·캐시 불필요), 꼭짓점은 평면 3장 연립(봉합 반경이 머신 정밀도 수준), 내/외 판정은 exact 술어(orient3d)로 문자 그대로 정확 — 이 세계에서 강건 불리언은 연구가 아니라 꼼꼼한 케이스워크(공면, 엣지-엣지 퇴화)다. 하이브리드 파이프라인(출처태그 메시 → 조합 결정 → 스냅백)을 스냅백이 자명한 평면에서 첫 완성. 알고리즘 참고: Manifold(Apache-2.0 — 차용·번역 가능), Hoffmann 등 문헌. 커버리지 밖 곡면 불리언은 명시적 미지원 에러로 정직하게 거절. `Discovered` 경로·국소 tolerance·relaxation 실전 투입. 이 시점에 "OCCT 없이 직동하는, 실용적 기계 부품(평면 위주)을 STEP으로 내보내는" 진짜 커널이 된다. 참고: Truck 대비 벤치마크·정밀도 비교(전역 1e-6 폴리라인 vs 정점별 실측 tol + 닫힌 형식)는 수치로 보여줄 수 있는 차별점 — 공개 지표 후보.

**M6 — 자체 불리언 2단: 이차곡면.** 평면∩실린더(타원), 평면∩구(원), 평면∩원뿔 — 여전히 닫힌 형식이라 행진 불필요. 실린더∩실린더는 특수 케이스(직교 등)부터. 실제 기계 부품 면의 대다수가 평면+실린더+원뿔이므로, 여기까지로 실용 커버리지의 대부분을 확보한다.

**M7 — 자체 불리언 3단: 일반 SSI (연구 구간).** nacre-geom::intersect에 SSI 행진 구현 → 일반 곡면쌍의 `HybridBoolean` 완성: 출처태그 tess에 강건 메시 불리언(exact predicates) → 조합 결정 추출 → 살아남은 면은 정확 곡면 유지, 신규 엣지는 국소 SSI 스냅백. OCCT 오라클과 상시 diff. 실패 케이스 코퍼스 축적.

M7은 열린 연구임을 명시한다. M6까지가 "확실히 되는" 영역, M7은 이 커널의 존재 이유이자 도박이다. OCCT의 역할은 전 구간에서 **오라클(dev 전용 시험관)뿐**이다 — 제품 경로에 OCCT 위임은 없다. OCCT 소스는 "상용급이 이 케이스를 어떻게 다루나" 열람용 참고서로만 쓴다: LGPL-2.1이므로 번역·차용은 라이선스 오염이고, 무엇보다 OCCT의 위상·tolerance 아키텍처가 딸려 들어와 nacre 설계와 충돌한다. 읽되 베끼지 않는다.

## 9. 미결 사항 (다음 논의 대상)

트리밍 곡면의 pcurve 표현 시점(M3에 선행 도입 vs M5까지 지연), 닫힌 엣지의 seam 처리 세부(파라미터화 경계를 tess 층 어디서 끊을지), Sketch 제약 솔버의 범위(초기엔 무제약 프로파일만), OpRef 계보 참조의 도입 시점과 직렬화 포맷 여유분, `Store` 스냅샷·직렬화 포맷(자체 vs STEP 재활용), OCCT history → 출처 매핑의 실제 충실도(M5에서 실측 필요), 멀티스레딩 경계(Store가 &mut 독점인 설계라 연산 단위 병렬은 미지원 — 의도적 단순화).
