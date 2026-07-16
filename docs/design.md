# nacre — 하이브리드 CAD 커널 설계

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
├── nacre-scalar     # [overhaul] exact 유리수 스칼라 + 2D-프레임 판정: Rat(Ratio<i128>, checked=강등 트리거)·Angle(유리수 deg·정확 누적·90°계열 exact·astro-float HP 삼각). `frame` 모듈 = orient2d_judge(f64 필터→astro-float 상승→declare-0). num-rational·astro-float 의존, 커널 타입 무의존(순수 토대)
├── nacre-predicates # [M5] indirect predicates(implicit point 부호 판정) 자체 구현. 바닥 = geometry-predicates(elrnv, MIT/Apache: orient3d + expansion primitive 노출). 순수 수치층(평면 계수·좌표 `[f64;N]`만, 커널 타입 무의존) → geom 무순환 + standalone 분리 가능
├── nacre-geom       # 정확 기하: Surface, Curve, 평가·미분·국소 교차(SSI relaxation)
├── nacre-topo       # Vertex/Edge/Face/Shell/Solid, half-edge, Model 집계
├── nacre-tess       # 출처 태그 tessellation: TessVertex, TessTriangle, 증분 갱신
│   └── polygon      # 평면 다각형 삼각분할: Newell 투영 + 구멍 브리징 + ear clipping
├── nacre-ops        # 연산: sketch, extrude, revolve, imprint, boolean(자체 — 커버리지 사다리). private `arrange` 모듈 = 면당 평면 arrangement(seam 세그먼트 수집 → 셀 추출). private `tolerant` 모듈 = 회전 boolean용 toleranced 술어 래퍼(회전 시 평면을 면 정점 Pt3로 재구성→frame3 judge·overhaul stage 3). `Model`과 geom을 둘 다 쓰므로 geom에 둘 수 없다(§1 의존 방향). nacre-tip 의존(회전 정점 Pt3·읽기전용·tip의 ops는 dev-only라 lib 그래프 무순환)
├── nacre-validate   # 불변식 검사: 오일러-푸앵카레, watertight, 방향성, 참조 무결성
├── nacre-props      # mass properties: 정확 기하 발산정리로 부피·면적(해석적, tess 무관). 소비: 사용자 질의·M5 부피보존 불변식·오라클 diff
├── nacre-tip        # [overhaul §TIP] 회전 정점의 정의+방향별 tol을 forest 순회로 조립(→ scalar::frame3::Pt3)해 orient3d/cmp_coord 판정: tol-0=predicates(Shewchuk exact)·tol>0=filter→astro-float. 단일 Discovered seam=indirect orient3d(2c-ii)·두 Discovered seam=cmp_coord(cmp-ii·ThreePlane→회전 정점 3개로 평면 재구성→frame3 간접 판정)·혼합cmp/이동-후-회전=defer. topo/scalar/predicates/geom 위 read-only 분석
├── nacre-step       # Model→AP242(Ed2) 엔티티 번역 어댑터. 직렬화 백엔드 교체 가능(커널 무지); 개발=step-io, 최종=경량 라이터
├── nacre-oracle     # [dev] OCCT 비교 하네스 (out-of-process 헬퍼 경유), proptest 전략
└── tools/occt-helper/  # 워크스페이스 밖 헬퍼: brew OCCT(1순위) 또는 uv+OCP(폴백) — §7
                        #   OCCT는 오라클 전용 — 제품 경로에 위임 없음 (§6, §8)
```

`nacre-geom`과 `nacre-topo`가 서로를 모르게 하는 것이 중요하다. 기하는 위상을 모르고(순수 수학), 위상은 기하를 Handle로만 참조한다. robustness가 첨예한 코드(교차·분류)는 전부 `nacre-geom::intersect` 한 모듈에 격리한다. 사용자는 파사드 크레이트 `nacre` 하나만 의존하며, 인터랙티브 스크립트 앱 등은 이 워크스페이스 밖의 별도 프로젝트로 둔다.

**디버그 뷰어는 커널 크레이트가 아니라 워크스페이스 밖 별도 앱이다.** 연산 로그를 입력받아 매 동작을 스텝별로 재생하며(append-only라 "N번째까지 replay"가 공짜), STEP에 안 담기는 nacre **내부 정보**(`Origin`의 `Constructed`/`Discovered`, `Discovered`의 tolerance 실측값, `Handle` 관계·인접 등)까지 시각화하는 인터랙티브 도구. 내부 자료구조에 접근해야 하므로 nacre를 **직접 링크**한다(개발 중 path 의존 → 안정화 후 version 의존, 버전별 디버깅도 자연스러워짐). 만드는 시점은 **`Discovered`/tolerance가 처음 등장하는 M5 즈음** — 그 전(M1~M4는 전부 `Constructed`)의 시각 확인은 정상 결과는 STEP→step-loupe(구조+검증), 중간·깨진 상태는 OBJ 덤프→맥 미리보기로 충분해, 인터랙티브 뷰어는 필요가 증명될 때까지 미룬다.

`Store`/`Handle`은 **최하위 `nacre-store`에 둔다.** geom도 Handle을 쓰기 때문이다 — `Curve::Intersection`(§3)이 `Handle<Surface>`를 담으므로, Handle이 topo에 있으면 geom→topo→geom 순환 의존이 된다. typed-index 저장소는 기하·위상을 전혀 모르는 순수 인프라이므로 두 층보다 아래에 격리하고, 위의 모든 크레이트가 자유롭게 참조한다. (라이선스는 MIT/Apache-2.0 듀얼 — Manifold(Apache-2.0) 알고리즘 차용과 호환.)

`nacre-scalar`는 **오버홀(§overhaul)의 §4 근본 표현 — exact 유리수 스칼라**를 격리한다. 사용자가 입력한 치수·각도를 f64 오차 없이 정확히 보존한다(`1.1`→`11/10`, `1.1×7`=정확히 `7.7` — "얇은 막" 문제의 근본 해결). `Rat`은 `Ratio<i128>` + **checked 산술**으로, 오버플로가 §4 강등 **트리거**(값의 캐시를 f64/dd로 내리고 tol을 `Origin`에 기록; 정의는 op-log로 불변 보존). `Angle`은 유리수 deg를 mod-360 **정확 누적**(한 바퀴가 정확히 0으로 닫힘 → 스케치 닫힘)하고 90°계열은 exact 유리수 cos/sin(회전 tol 0). **★ `nacre-predicates`와 상보(겹침 아님):** predicates는 기하 행렬식의 **부호를 exact 결정**(exact-부호), nacre-scalar는 **입력 값과 유리수-순수 누적을 exact 보존**(exact-값) — 역할이 갈려 이름·층이 분리된다. **의존 결정(오버홀 최초 새 외부 dep):** `num-rational`(+num-traits)을 채택 — 성숙한 checked 유리수+gcd 약분을 제공하고, exact 유리수를 손수 구현하면 버그가 exactness 목표를 훼손하기 때문(MIT/Apache·순수 Rust). 어떤 `nacre-*`에도 의존 않는 **의존 그래프 최하단 순수 토대**. **2D-프레임 판정(`frame` 모듈, 이식 #2):** 회전 좌표는 무리수라 exact 못 하지만 **부호**는 사다리로 sound하게 정한다 — `orient2d_judge`가 f64 필터(좌표혼합 tol `(|bx|+|by|)·da`, H4 교정)로 쉬운 케이스를 처리하고, 애매하면 **astro-float**(H1.5, 200비트)로 상승, precision floor 아래면 **declare-0**(`Orient::Zero` 값 반환; 병합/별개 ask-user 정책은 §6·후속). predicates(M5 geometry-predicates, 축정렬 격자)와 여전히 상보 — TIP 상승 판정이라 nacre-scalar에 둔다(부피 커지면 `nacre-tip`으로 분리 후보). **범위(이식 #1·#2):** exact 값 엔진 + 2D-프레임 판정까지. 통합 값+tol `Scalar`·PtChain 전파(3D 체이닝)·declare-0 정책·커널 배선은 후속 셀.

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

정밀도 분업 원칙: **판정(부호)에는 적응 정밀 술어**(`geometry-predicates` — 쉬운 케이스 f64, 아슬아슬할 때만 확장, 부호는 항상 정확), **반복 구성(좌표)에는 고정폭 확장 부동소수점**(double-double). 임의 정밀 유리수(rug/malachite)는 술어에는 완벽하지만 Newton 반복에 넣으면 비트 길이가 반복마다 폭발하므로 구성에는 쓰지 않는다. 이 투자는 평균이 아니라 꼬리를 산다 — 호출 빈도가 낮은 악조건 케이스만 정확히 개선되고, tolerance가 "상수"가 아니라 "실측 보증값"이 된다. (**오버홀의 유리수 경계 규칙은 `docs/overhaul.md §4`**: 유리수는 **입력 표현**일 뿐 계산 매체가 아니다 — 입력·유리수-순수 파생에만 살고, 불리언·반복·혼합이 닿는 순간 f64/double-double로 강등해 이 비트 폭발 경고와 정합한다.)

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
    // 기존 면 위에 작업하는 M4(ImprintSketch/PadOnFace)에서 비로소 필연적으로 도입된다.
    // M2는 매 Extrude 결과가 닫힌 솔리드라 validate가 빈틈없이 걸린다.
    Extrude { plane: SketchPlane, profile: Profile2d, dist: f64 },
    Revolve { plane: SketchPlane, profile: Profile2d, axis: Axis, angle: f64 },
    ImprintSketch { face: Handle<Face>, profile: Profile2d },        // M4: 기존 면 위 — Handle<Face> 참조
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

**패드/포켓 ↔ 불리언 통합 — feature = tool body + boolean (목표, 미착수).** M4 `PadOnFace`/`PocketOnFace`(`imprint`+`raise_region`, 직접 구성)와 M5 불리언 공면 경로(`detect_pocket_contact`·`detect_contained_contact`·overhang 사이드카)는 **부분 중복**이다: top-flush 포함 포켓/보스는 `PocketOnFace`로도 `Cut(솔리드, 프리즘)`으로도 나온다(`cut_a_blind_pocket_into_a_face`). 둘 다 있는 건 **마일스톤 순서의 산물** — M4가 M5보다 먼저라, 불리언이 없던 시절엔 재료 가감이 직접 구성뿐이었다. 목표 구조는 상용 CAD와 동형: **`pad` = 프로파일 압출 → `Fuse`, `pocket` = 압출 → `Cut`**, 패드/포켓은 불리언의 얇은 sugar. "프로파일이 면 안에 있어야 한다"는 제약이 사라지고(오버행 패드 자동 획득), M4 직접 기계(`imprint`/`raise_region`/`prepare_face_split`)를 제거해 코드 표면이 준다. 공면 복잡도는 없어지는 게 아니라 **불리언 하나로 집약**된다(독립-솔리드 불리언 때문에 어차피 필요).

- **★ 구조적 공면 = O(1) 참조 인식(핵심).** 면 위 스케치를 압출한 tool body의 밑면은 대상 면의 **surface Handle을 공유**한다(`imprint`가 이미 원본 `surface_h`를 재사용). 통합 불리언은 **먼저 Handle 공유를 검사**(float 계산 0)해 공유 면을 seam으로 즉시 채택하고, **공유하지 않는 독립 솔리드만** `planes_coplanar` 기하 감지로 내려간다. 즉 "참조로 아는 공면은 공짜, 우연한 공면만 계산" — 상용 커널이 "coincident 면을 imprint로 공유 토폴로지로 승격"하는 것의 nacre판.
- **★ 게이트(순서 강제).** (a) 불리언 공면 처리를 하나로 흡수·완성(지금 셀 사다리 — Cut 3갈래 통합, Common이 Cut의 detect 재사용) → (b) Handle-공유 O(1) fast path 추가 → (c) 그제서야 pad/pocket을 extrude+불리언 wrapper로 바꾸고 M4 직접 경로 제거. **역순 금지**: 불리언 공면이 robust해지기 전에 M4를 걷어내면 지금 되던 포켓/보스가 커버리지 구멍에 빠진다(M4 직접 경로는 그때까지 신뢰 가능한 fallback).
  - **(c-1) pocket 완료 (실행됨).** `PocketOnFace`는 `build_prism`(inward, top-flush) + `boolean(Cut)`으로 재구현됐다 — contained-coplanar Cut 경로가 빈 seam을 내므로 결과는 전부 Constructed(tolerance 0), 옛 직접 `raise_region`과 등가. 선결로 `detect_pocket_contact`의 볼록성 게이트를 제거(비볼록 kept `a`·비볼록 프로파일 `b` 둘 다 열림, OCCT로 교차검증)했다. `raise_region`은 `pad`용으로 잔존. through-pocket(`dist ≥ 두께`)은 이제 정직히 거절(M4는 미검사 UB였음 — 개선). **게이트 (b) Handle-공유 fast path는 불필요로 폐기**: contained 경로가 이미 빈 seam·전부 Constructed라 참조-fast-path가 더할 exactness가 없다. 다음: pad(=extrude+Fuse), 그다음 M4 직접 경로 완전 제거.
  - **(c-2) pad 완료 (실행됨).** `PadOnFace`도 `build_prism`(outward, top-flush) + `boolean(Fuse)`로 재구현 — pocket과 대칭. 선결로 `detect_contained_contact`의 볼록성 게이트를 제거(비볼록 kept base·비볼록 프로파일 boss 둘 다 열림, OCCT 2종 교차검증). pad·pocket은 이제 부호(±dist)·`BoolKind`·복원-면-부재 에러만 다른 **공통 헬퍼 `extrude_and_boolean` 위의 얇은 wrapper**로 통일됐다("feature = tool body + boolean"의 코드화). `raise_region`은 유일 사용자였던 pad가 떠나며 **삭제**(그 죽은 `Split` 필드도 함께 정리). 남은 M4 직접 기계는 imprint(`imprint`/`prepare_face_split`)뿐 — 완전 제거는 imprint 정리 후 별도 스텝.
  - **(c-3) 오버행 pad/pocket 획득 (실행됨 — 로드맵 payoff).** `placed_profile`을 `placed_profile_unchecked`(CCW+배치)와 strict-containment wrapper로 쪼개고, `extrude_and_boolean`이 unchecked를 쓰게 했다. 이로써 **면 경계를 넘는 프로파일이 기존 오버행 불리언 사이드카로 자동 라우팅**(Fuse: 단일 엣지·코너·spanning slab; Cut: N-wall blind)된다 — "프로파일이 면 안에 있어야 한다"는 제약이 사라졌다. contained 경로는 완전 불변(같은 base_pts). 커버리지 밖(비볼록 오버행·through·far·비축정렬)은 **정직 거절**(`Boolean(Unsupported)`) — silent-wrong 아님. 오버행은 contained(빈 seam)와 달리 footprint crossing에서 **Discovered seam 정점**을 만든다(기존 오버행 셀의 성질, 새 tolerance 아님). imprint는 containment **유지**(유효 inner-loop 필요). **알려진 비대칭**: contained pad/pocket은 비볼록을 받으나 오버행은 볼록 게이트로 볼록만 — 비볼록 오버행(오버행 detect 볼록 게이트 제거)이 다음 후보.
  - **(c-4) 오버행 boss가 비볼록 솔리드 수용 (실행됨 — c-3 비대칭 절반 해소).** `detect_overhang_contact`(Fuse)의 whole-solid `is_convex` 게이트를 **접촉면 footprint 게이트**(두 접촉면 outer loop 볼록 + 홀 없음)로 교체. 근거: Fuse 재구성은 **로컬** — arc-split은 (볼록) 접촉면만, `resplit_overhang`은 wall을 edge-local(볼록 무관)로, 나머지 면은 verbatim 재방출. 그래서 **비볼록 솔리드(포켓 파인 부품·부울 결과)에 접촉면만 볼록이면 boss 캔틸레버**가 붙는다. n0로 OCCT-정확 확인(포켓 큐브 옆면 오버행 = 1.17). silent-wrong 원천인 비볼록 **접촉 footprint**(arc 오분류)는 게이트가 거절. **Cut/Common은 whole-solid 게이트 유지** — 그 `clip_bwall_inside_a`가 breached 반평면 SH 클립이라 볼록 kept에서만 "inside a"와 일치(비볼록은 `OVERHANG_ARCS` 정직 거절, n0 실증). 신규 `loop_is_convex_2d`·`face_outer_is_convex`. **후속**: ② 비볼록 footprint 오버행(multi-piece 재구성), 비볼록 솔리드 오버행 Cut(clip_bwall_inside_a 일반화).
- **보존 불변식.** M4 직접 경로는 교차를 아예 계산하지 않아 **전부 Constructed**(tolerance 0)다. (b)의 참조-fast-path는 그 순수성을 유지해야 한다 — 공유 면 seam은 Discovered로 승격하지 않고 Constructed로 남기고, **우연한 공면만** Discovered seam(+tol). 통합이 exactness를 후퇴시키면 안 된다.
- **업계 정합.** 피처 레이어는 pad=tool body+boolean으로 통합돼 있고(SolidWorks/NX/Creo/Fusion), 커널 레이어는 그래도 coincidence 전담 로직(imprint+tolerance)을 보유한다 — 사이드카는 사라지는 게 아니라 **커널 불리언 안으로 들어간다**. nacre는 tolerance 대신 exact 술어로 그 자리를 채우는 소수파.

**완전한 M5 오버홀 (회전 지원) — 전체 설계는 `docs/overhaul.md`.** 축정렬 전용의 근본 한계를 넘어 회전(비정렬·마름모·각도 스케치)을 지원하는 실사용 M5로 가는 오버홀을 별도 문서에 둔다: 유리수 치수·각도 표현, 명시 공유(전역 자동병합 아님 — TNP 충돌), TIP 판정 정책(declare-0), op-log 소유·sugar, 참조-공면. **핵심 novel 메커니즘은 별도 `experiment/exact-2d` 브랜치의 격리 2D 실험으로 선검증**한 뒤 이식한다.

**불리언 다중 솔리드 출력 + non-manifold 정책 (OCCT n0로 확정).** 실측: 엣지로만 만나는 두 박스의 `Fuse`, 막대를 가르는 `Cut` — OCCT는 둘 다 **COMPOUND of 2 SOLIDs**(checkshape valid, STEP `2× MANIFOLD_SOLID_BREP`, non-manifold 엔티티 0)로 낸다. 엣지 접촉 Fuse는 접촉 엣지+2정점을 **공유**(V16→14·E24→23)하되 각 솔리드는 깨끗한 manifold 박스다. ⇒ 세 결정:
- **`DISCONNECTED_RESULT` 은퇴 → 다중 솔리드 출력.** 불리언이 `Vec<Solid>`를 반환한다. "Cut 두 동강(§10 5c)"·"엣지 접촉 Fuse"가 전부 정상 결과. `live_solids: Vec<Handle<Solid>>`(§2) 인프라 이미 존재 — 불리언 쪽 제약(양수 셸 정확히 1개 요구)만 걷어 **전 positive 성분을 각 솔리드로 방출**(cavity 부호 분류 (5d) #5와 연동).
- **non-manifold 솔리드 미지원 — 불필요.** OCCT도 안 만든다(엔티티명부터 MANIFOLD). 각 솔리드 안의 manifold 전제(엣지당 면 2·Euler)는 그대로. 공유 엣지가 compound 레벨에서 4면인 것은 각 솔리드가 2면씩이라 위반이 아니다. validate가 **per-solid**면 안 건드림(전역 엣지-면 카운트면 공유 경계 오탐 → 그때만 per-solid 조정).
- **솔리드 간 경계 요소 공유가 필수.** OCCT가 접촉 엣지를 공유하듯, `Vec<Solid>`의 솔리드들이 정점·엣지·면 Handle을 공유할 수 있어야 한다 — 오버홀의 명시 공유와 같은 메커니즘이며, OCCT는 그 공유를 **불리언이 명시적으로 만든다**(전역 자동병합 아님).
- **★ 교정.** 오버홀 논의 중 "요소 공유가 non-manifold를 부른다"는 우려는 **반증됐다** — manifold는 솔리드별 판정이며, 두 솔리드가 경계를 공유해도 각각은 manifold다. 공유는 문제가 아니라 오히려 불리언이 "여기서 만난다"를 계산 없이 알게 해준다(참조-공면의 확장).

### (5d) exactness sweep — 전수조사 표 (진리를 정하는 자리에서 f64를 읽는 곳. 단일 진실원.)

M5 불리언의 **위상 결정**(어느 것이 안/밖·볼록·공면·outer/cavity인가)이 좌표 f64를 읽는 자리를 전수 나열한다. 진리를 정하는 술어만 대상이다 — 부피·validate·OCCT는 net(사후 검산)이라 제외. **정점의 존재론적 tol(`vertex_tol`)은 제외** — Discovered 정점은 "정의 + tol"이 정상 표현이며(원칙상 합법, DNA), 은퇴 대상이 아니다. 은퇴 대상은 "판정이 tol/f64 캐시를 읽어 조용히 틀릴 수 있는 곳"뿐이다.

| # | 자리 | 무엇을 정하나 | 현재 | 뿌리·비고 | 크기 |
|---|---|---|---|---|---|
| 1 | `is_convex` (ops) | 볼록성 → coincident 병합 분기 | **exact `plane_side` ((5d)-1 완료)** | 옛 tolerance가 정점 f64 비공면을 흡수. n0: 도달 피연산자 전부 Constructed 축정렬이라 이미 exact → tolerance 은퇴는 **무동작 청소**(Case A). Discovered 정점이 닿으면 #3로 넘어감 | 완료(작음) |
| 2 | `coplanar` (geom/ops) | 두 면 공면 여부 → 인터페이스 탐지 | **exact rank-1 `planes_coplanar` ((5d)-2 완료)** | 옛 `1e-9` 절대-길이 tolerance는 scale-비불변(근접-공면 false-merge 위험). exact = 두 평면 계수 2×4의 rank-1(여섯 2×2 minor). n0: 유일 불일치는 float-박스 옆면의 같은-법선 공면쌍뿐 → 반대-법선 인터페이스 필터가 배제 → **동작 보존** | 완료(작음) |
| 3 | `is_convex`의 **triple-sourcing** | Discovered 정점을 f64 캐시가 아니라 제 정의(triple)로 | **exact `three_plane_orient3d` ((5d)-3 완료)** | ★ 술어를 새로 안 지음 — `plane_orient(Q1,Q2,Q3,P)=sign(det4)·sign(det3)`는 이미 `three_plane_orient3d`(det3 인수가 술어에 내재, 순서 무관). 발화 픽스처 = **불리언 결과를 피연산자로 넣는 첫 사례**(볼록 `Common` 결과=Discovered 코너 6개 → 셋째 박스 스택). M5 동작 보존(축정렬 교점 f64-exact라 triple=캐시); 소득은 결정에서 캐시-읽기 제거 + ingestion 경화 + TIP 대비 | 완료(중) |
| 4 | `fan_triangles`의 zero-area 드롭 (ops) | 팬 삼각형 버릴지 (winding 경로) | **exact `orient2d` zero-area ((5d)-4 완료)** | winding 부호는 이미 exact(`ray_triangle_cross`=orient3d). 유일한 tol은 `fan_triangles`의 상대 `1e-12` 공선 드롭 — 대좌표 근접-공선 슬리버를 조용히 버려 교차를 놓칠 수 있었다. exact = 외적 세 성분(=세 좌표투영 `orient2d`)이 모두 0. n0: 전 코퍼스에서 nonzero-면적 발화 0건 → **Case A**. ~~광선-free 재설계~~ **기각**: 비볼록에 더 단순한 exact 광선-free 없음(GWN=f64), 방향 휴리스틱은 honest-reject라 exactness 이득 0. `ray_triangle_cross`의 정점 f64-좌표 읽기는 회전 시 TIP의 몫(이 셀 밖) | 완료(작음) |
| 5 | cavity 분류 부호 (ops, (5c)) | 성분이 outer(+)인가 void(−)인가 | **exact extreme-vertex 부호 ((5d)-5 완료)** | 옛 f64 signed-volume flux(sqrt-area·나눗셈-centroid·혼합부호 합)가 위상 라벨을 정했다. exact = 성분의 lex-최소 정점(볼록 코너)에서 `n_x<0`인 인접면이 있으면 outward. **좌표 산술 없음**(순서비교 + 축정렬 법선부호뿐) → `0.3/0.7` 비표현 좌표에서도 exact, **total**(거절 안 함). ~~containment-parity~~ **기각**(감김 재사용이나 total→partial·O(n²)); ~~exact signed-volume 합~~ **기각**(합의 exact 부호 인프라 없음). n0: 전 코퍼스(ops 217 + OCCT 70) 부호 뒤집힘 0 → **Case A** | 완료(중) |
| 6 | `interface_correspondence` (ops) | coincident-merge 정점 대응 (B링 → A링) | **exact 좌표 동일성 `ap==bp` ((5d)-3 완료)** | (5d)-2 완결성 grep이 발견 — coincident 경로의 두 번째 f64-읽기. 인터페이스 정점은 **Constructed라 triple 정의가 없다** → 정답은 좌표 동일성(대응 코너는 같은 점, bit-동일). 견고했으나(순수성 청소) scale-상대 tol 은퇴. 전 fixture Case A | 완료(작음) |
| — | `vertex_tol` (Discovered 정점) | — | tol 보유 | **은퇴 안 함.** 존재론적 tol이며 판정 fudge가 아니다. TIP가 이 tol을 전파·판정에 반영하는 층(회전 도입 시) | 해당 없음 |

**★ "평면 계수 exact화"만으로 #1이 열리지 않는다.** 축정렬 정수 입력의 Constructed 정점은 정규화 전에도 정확히 공면이었다 — (5d)-1의 계수 exact화는 **토대**(그 위에서 `plane_orient`·`coplanar`가 참 평면에 서게)일 뿐, is_convex의 과잉 거절을 실제로 없앤 것은 tolerance 은퇴다. is_convex의 남은 exact화(Discovered 경로)는 **#3 triple-sourcing**이 열고, 그건 발화 픽스처가 있어야 짓는다.

**진행:** #1·#2·#3·#4·#5·#6 **전부 완료 → (5d) exactness sweep 종료.** n0 완결성 감사(nacre-ops 프로덕션 grep)가 clean: **tolerance 리터럴이 하나도 남지 않았고**(`vertex_tol`만, 그건 존재론적 제외 행), 나머지 f64 부호 비교는 전부 exact 술어 부호(`orient3d`/`plane_side`/`planes_coplanar`)·축정렬 구성-방향·축정렬 법선 dot 부호뿐 — tol/누적 기반 위상 판정 0. §10 (5b-0)/(5c) 기록의 "전수조사 표" 참조는 모두 이 표를 가리킨다.

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

**TIP — Toleranced Indirect Predicates (회전 도입 시의 부호 판정 층. 개념 확정·수식 미확정. 지금 구현 금지.)**

**왜 필요한가.** M5~현재의 판정은 좌표가 유리수라 exact다. 그러나 **회전이 들어오면 좌표가 무리수가 된다** — 유리수 각도라도 cos·sin은 무리수이고(Niven), 임의 각도는 초월수다. exact 산술로 표현할 수 없으므로 회전 좌표는 f64 근사일 수밖에 없다. **회전은 M5의 관심사가 아니다** — 프로덕션에 솔리드 회전이 없고 픽스처에도 없다(`rotated_l_prism`은 기하 회전이 아니라 **정점 순서 순환 이동**이다 — 이름이 오해를 부른다). 문제는 회전이 도입되는 시점(모델링 변환, 각도를 쓰는 스케치 포함)에 발생하며, 그때 **판정만은 조용히 틀리지 않게** 지키는 설계가 TIP다.

**핵심 명제.** indirect predicates(Attene 2020)는 **선형 요소(선·평면)의 교차**, 즉 다항식에만 성립한다. 회전은 다항식이 아니므로 **indirect predicates를 회전에 적용하는 방법은 없다**(문헌에도 없다 — 수학적 한계). 실제 연구·구현이 하는 일은 하나뿐이다: **회전 좌표는 f64로 근사하고, 그 위의 판정을 exact로 유지한다.** TIP는 그 "위의 판정"에 **입력 근사(tol)까지 반영**하는 층이다.

**★ 대체가 아니라 흡수다.** TIP는 indirect predicates를 지우지 않는다. 점을 "정의 + 누적 tol"로 들고, **tol = 0이면 지금의 indirect predicates 그대로(exact)**, **tol > 0이면 필터 + 고정밀 폴백**이다. 술어의 계산 구조(3-평면 정의, 행렬식 부호, 좌표 무독)는 **그대로 쓰고**, "평면 계수가 정확하다"는 가정만 "계수에 tol이 있다"로 넓힌다. 회전된 평면은 **계수에 tol이 붙은 평면**일 뿐이다.

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
6. **직선 구간은 fraction 강체로 묶는다.** 연속된 직선 이동을 유리수로 먼저 합산하면 그 구간의 **연산 오차 = 0**이고 이산화는 **회전과 만나는 지점에서 1회**뿐이다. 순차로 f64 덧셈하면 이동마다 이산화가 붙는다((a+b)+c는 이산화 2회, a+(b+c)는 1회). 회전 결과는 무리수라 fraction으로 못 담으므로 **회전 경계의 f64 덧셈 오차는 피할 수 없다** — 피할 수 있는 것만 피한다. **(일반화 — `docs/overhaul.md §4`: 기준은 "직선/회전"이 아니라 "유리수/무리수"다. 직선이라서 tol 0이 아니라 유리수라서 tol 0이며, 강체 묶기는 직선뿐 아니라 **같은 축의 연속 유리수-각도 회전**에도 적용된다[각도를 유리수로 합산해 누적 각을 1회만 실현 — 2D는 항상 한 축, 3D는 같은 축만]. 유도된 무리수[√ 거리·구속 솔버 해·3D 축 변경 합성]가 섞이면 그 구간은 fraction으로 못 묶는다. tol 공식의 거리·각도 항은 지우지 말고 0으로 둔다. **★ 번들링은 필수(H4-amplification 실측): un-bundled 증분 실현은 전파 `|R|`의 행합 `|cos|+|sin|≥1`을 매 스텝 곱해 tol 바운드가 지수 폭발[실측 30스텝에 실제 오차의 ~20만 배; 실제 오차는 R 노름보존이라 평평]. 누적각 1회 실현이 차단. 폭발해도 sound 최악-보장이라 조용히 안 틀리고 상승/거절로만 간다[항목 8].)**
7. **판정 = 필터 + 정밀도 상승 + lazy 캐싱.** 점들의 tol로 이번 판정의 오차 한계를 계산(변 길이로 증폭) → `|행렬식| > 오차 한계`면 f64로 확정(대부분) → 애매하면 정밀도 상승. 관련 점이 **tol = 0(평면 교차점)뿐이면 exact 폴백**(현행 indirect predicates), **회전 점이 끼면 f128**. 한 번 고정밀 계산한 **값**은 캐시해 재사용하고(판정 결과가 아니라 점의 값 — 여러 판정에서 재사용된다), 재계산 시 **회전 지점만** 다시 계산하고 직선 구간은 fraction을 그 정밀도로 이산화해 잇는다. **판정에 필요한 점만 정밀화한다** — 중간 경유점은 정의로만 남긴다.
8. **오차 한계는 최악(선형 합)이다.** 확률 전파(RSS)는 실제 오차에 가깝고 100배 작지만 **보장이 아니다** — tol은 "이보다 클 수 없다"는 보장이어야 하고, 확률 tol은 드물게 부호를 뒤집어 **조용히 틀린다**. 커짐은 다른 방법으로 관리한다(직선 tol = 0이라 N은 회전 수뿐, 애매하면 f128).
9. **★ 필터는 동적으로 통일한다.** 지금 (5b-0)의 필터는 **정적**(오차 한계가 상수 `ε_D ≈ 5u`·`ε_M ≈ 19u`)이고, 입력 tol이 없으므로 그것이 옳다. 회전이 들어오면 **동적으로 통일**한다 — 회전 tol은 경로마다 편차가 커서 하나의 정적 상수로 잡으면 느슨해 무용하거나(과대) 위험하다(과소 → 터진다). 그리고 **각도를 쓰는 스케치도 회전이므로 무회전 케이스는 극소수**이고, **동적 필터는 tol = 0인 점을 자동으로 정적급 타이트하게 처리**하므로 정적 분기를 남기는 것이 순수 오버헤드다.

**★ 근본 한계 (정직하게).** 회전이 들어가면 **서로 다른 경로로 같은 위치에 도달해도 f64가 다를 수 있고, exact 술어로도 그 다름을 그대로 반영한다** — 술어는 "주어진 입력에 대해" 정확할 뿐 입력 자체의 근사를 고치지 못한다. f128로도 **보장은 불가**하다(초월수의 상등은 유한 정밀도로 결정 불가). 다만 **실무적으로는 충분하다**: 상용 CAD의 uncertainty가 ~1e-9인데 f64는 1000mm 점을 10000번 회전해도 ~1e-10, f128이면 ~1e-28이다. **이론적 완벽함은 포기하고 실무적 충분함을 tol로 보장한다** — 그리고 그 tol을 **측정된 값**으로 들고 다니는 것이 tolerance-fudge와 다른 점이다. (f128로도 애매한 **잔여 케이스**는 커널이 0으로 추측하지 않고 **사용자에게 확인받는다** — `docs/overhaul.md §6`. 즉 여기 "tol 보장"은 f128이 부호를 확정하는 대다수 경로를 말하고, 확정 못 하는 잔여는 ask-user로 넘긴다.)

**마일스톤별 적용.**
- **M5·M6:** 판정이 다항식이므로 TIP가 그대로 얹힌다. **M6에서 TIP의 가능 여부 = 간접 술어의 가능 여부**다(TIP는 술어 위의 층이므로) — 평면∩이차곡면 간접 술어가 서면 TIP도 선다. 단 술어가 고차가 되면 **오차 한계 공식을 그 술어에 대해 재유도**해야 하고, 고차라 필터 성공률이 M5만큼 좋을지는 **측정 대상**이다.
- **M7:** 메시 조합 판정의 필터 + 뉴턴 스냅백 점의 tol 추적·f128 재수렴. **SSI(교차를 찾는 것)는 TIP 밖**이다(위상 존재 문제이지 tol 문제가 아니다) — 아래 "TIP와의 관계 정리" 참조.

**아직 정하지 않은 것 (구체화 대상 — 도입 직전에 정한다).**
1. **~~회전 연산의 자체 tol을 무엇으로 표현할지~~ → 확정(H4-soundness 실측 교정).** **모든 tol(회전·직선·이산화)을 xyz 방향별 벡터로 표현·누적**한다. 스칼라(구) tol은 탈락 — 정보를 버린다(z축 둘레 회전이면 z 방향 오차 0인데 스칼라는 모든 방향에 실어 "z=5 평면 위인가?"를 불필요하게 애매 판정→상승 유발). **★ 회전 tol은 세 항의 합이다** (초기 "접선형 새-오차" 안은 H4-soundness에서 반증됨 — 아래 (a)): (a) **새 오차 = 실현 반올림(좌표 혼합)** + (b) **기존 tol 방향 회전** + (c) **각도 불확실성(접선, v1엔 0)**. 회전 시 tol 갱신 =
```
새 tol = |R| · (기존 tol)            [b: 기존 tol의 방향을 회전 — 성분별 절댓값 행렬 |R| × 기존 tol]
       + (|x| + |y|) · da_이산화      [a: 새 오차 — 실현 반올림, 좌표 혼합]      ★ H4 교정
       + da_각도 × |n × (P − 축점)|   [c: 각도 불확실성, 접선 — v1엔 0]
```
**★ (a) 교정(H4-soundness).** 초기 안은 새 오차를 **접선형 `da × 모멘트암`**으로 뒀으나 **반증됐다** — cos/sin을 독립 반올림하면 순수 회전이 아니라 **방사 성분**이 생겨(`err_x = |x·δc − y·δs|`), 접선형은 좌표 하나가 0 근처일 때 그 성분을 **과소평가**한다(실측: 10000 랜덤 중 220회 상한 붕괴). 성분별 새-오차 tol은 **두 좌표를 합친 `(|x|+|y|)·da_이산화`**여야 sound(0으로 안 꺼짐; 실측 0회 붕괴·최악 ~3.5× 이내 보수). `da_이산화`=실현 반올림 한계(f64 ~16 ulp; 상승 시 ~2^-P). (b) `|R|·(기존 tol)`은 **필수** — x방향 tol이 z축 90° 회전 후 y로 옮겨가는데 안 반영하면 상한이 깨진다. (c) **각도 불확실성**(접선)은 유리수 각의 v1엔 0이나(그리고 (a) 여유가 f64 각도-반올림까지 흡수), 유도 치수(구속 솔버·비유리수 각)에 살아나므로 **항은 지우지 말고 0으로 둔다**. `n`=회전축 단위벡터, `P`=판정점. 직선 이동은 (a)(c)를 안 만들지만 앞선 tol을 싣고 간다. 판정 시 각 방향 tol이 행렬식에서 자기 계수로 증폭돼 오차 한계를 이룬다. (**overhaul.md §6**; CGAL Lazy_kernel의 interval과 같은 이유.)
2. **~~입력 tol → 행렬식 오차 한계 공식~~ → orient3d는 확정·이식(단계 2a-ii): `frame3::det3_bound`**(6개 signed triple-product 구간 반경 `prod_err` + `16ε·mag` f64 반올림; `orient3d_judge`가 필터→astro-float 상승→declare-0로 소비). exact3d H-a 검증(위반0·tightness 0.12; 프로덕션 재현 0.065·피벗 포함). `plane_side`=explicit orient3d 쌍둥이(같은 공식). **간접 orient3d**(3평면 implicit point·평면 계수)도 확정·이식(단계 2c-i `indirect_orient3d_judge`; exact3d H-b/H-c 검증·프로덕션 재현 위반0). ①이 확정돼 **입력이 명확**하다 — 점 tol = xyz 방향별 벡터. "방향별 점 tol이 행렬식에서 각자 자기 계수로 증폭돼 오차 한계를 이루는"(항목 5) 공식을 술어별로 유도(현행 정적 상수는 tol = 0 가정). **H4-soundness 진행 상황**: **(a) 새-오차·(b) 전파 항 모두 실측 검증 완료.** (a) 접선형 반증→좌표혼합(10000 랜덤 0회 붕괴). (b) `|R|·기존 tol`: 증분 회전 체인 5000회 0회 붕괴(sound·~3× 타이트); 전파의 **필요성은 비대칭 tol**에서 발현(x-tol을 90° 회전 → 오차가 y로 이동, 전파 없으면 tol_y=0으로 과소예측 — 대수적 확인; 비대칭 tol은 v1-후 √거리·구속 솔버에서). 아래 ③(평면 계수 tol 전파)도 같은 유도의 일부. **남은 H4**: consistency(경로-독립 결정성)·amplification(증폭 경계)·속도.
3. **~~평면 계수의 tol 구조~~ → sqrt 섭동은 (5d)-1이 해결, 나머지는 ②로 흡수.** `Plane::through_points`의 정규화(sqrt) 섭동은 미확정이 아니다 — (5d)-1이 Plane에 비정규화 `raw`를 함께 저장해 판정이 쓰는 `coefficients()`가 exact(정의 정점이 정확히 0)이고 `normal()`(sqrt)은 크기 소비자에게만 간다 → **판정 경로에 sqrt 섭동 없음**. 남은 **"회전된 평면의 계수 tol"**은 별개 미확정이 아니라 **①의 따름정리**다: 평면은 세 점으로 정의되고, 회전되면 그 점들이 ①의 방향별 tol을 가지며, 계수는 그 점들의 뺄셈·외적이라 **점 tol이 계수 tol로 전파**된다. ②의 오차 한계 유도의 일부로 함께 다룬다.
4. **~~3D에서 "거리"의 정확한 정의~~ → ①에서 확정.** 새-오차 (a)는 **좌표 크기 `|x|+|y|`**로 실측 확정(H4). 각도-불확실성 (c)의 접선 크기 = `da_각도 × 회전축까지 수직거리(모멘트 암)`(v1엔 0).
5. **~~방향별 tol의 자료구조~~ → 확정: `[f64;3]` xyz 벡터**(단계 2a). `nacre-scalar::frame3::Pt3 { base:[Rat;3], chain, coord:[f64;3], tol:[f64;3] }` — 2D `frame::Pt2`의 3D 아날로그. (**회전 시 tol 변환**은 ①의 `|R|·기존 tol`로 흡수.) 임의 유리수 피벗까지 exact3d H-f로 검증(피벗 산술 tol 항 추가).
6. **~~고정밀 층의 선택~~ → 확정: `astro-float`(순수 Rust 임의정밀). twofloat(double-double)는 H1.5에서 탈락.** H1.5 실측: twofloat의 π 상수는 정확하나 삼각함수가 "preliminary"라 **영점 근처 cos 오차 ~1.8e-16(f64 수준)**으로 정확도 게이트(~1e-30) 실패 — 부호 판정이 일어나는 near-degenerate가 곧 영점 근처라 치명적. `astro-float`로 교체(160비트에서 cos/sin 오차 **~1e-58**로 통과). **정밀도가 dial 가능**이라 double-double의 ~1e-32 천장이 사라진다(H4-amplification이 tol을 키워도 정밀도를 올리면 됨 → "quad-double 상승 or 정직 거절" 딜레마가 "정밀도 dial"로 단순화). 대가는 속도 **~120µs/call**(상승 경로라 드물고, 각도당 실현을 캐시해 상각; 정밀도를 낮추면 빨라짐) — H4가 속도·상승빈도 실측. **★ 이하 이 문서·overhaul.md의 "double-double"·"f128"·"고정밀 상승"은 이 상승 층을 가리키는 일반명이며, 실제 구현은 astro-float다.**
7. **회전 이력 트리 → 확정(단계 1c): `Model.rotations: Store<Rotation>` + `parent` 링크 forest.** 남은 **캐시 정책**(hp 실현값 수명·축출)은 단계 2a에서 미실행(Pt3 매번 fresh 계산) → **단계 2b 이후로 defer**(성능 최적화, soundness 무관).

**★ 정리 — ②(오차 한계)가 완전 소진됐다: orient3d/plane_side 직접(2a-ii `det3_bound`)·간접(2c-i `indirect_orient3d_judge`)·cmp_coord(cmp-i `indirect_cmp_coord_judge`) 모두 확정·이식(프로덕션 `nacre-scalar::frame3`).** 남은 §TIP는 ⑦-캐시(성능·defer)뿐. ①·④·⑤·⑥ 확정, ⑦-트리는 1c 완료, ③은 sqrt 부분 (5d)-1이 해결·나머지 ②로 흡수. §TIP의 마지막 수학 ②는 실험 H-a(직접 orient3d)/H-b(평면 계수)/H-c(간접 orient3d)가 검증했고, cmp_coord는 exact3d에 없던 새 수학이라 frame3에서 직접 **H-g**(두 코퍼스·wrong-sign 0)로 검증했다(Cramer 기계는 2c-i 이식분 재사용·최종 부호 결합만 신규). 세 술어 모두 프로덕션에 들어왔다 — **판정층 완성**.

**측정해야 아는 것:** 고정밀 층의 속도, 동적 필터의 실제 성공률(입력 tol이 있을 때 정적 대비 얼마나 자주 exact로 떨어지는가), 실무 형상에서 회전 tol이 실제로 얼마나 커지는가(판정 경계에 근접하는가).

**전제:** **(5d) 완료.** 판정이 좌표가 아니라 평면 계수·정의만의 함수가 되어야 그 위에 "계수의 tol"을 얹는 것이 깨끗하다. (5d) 전에 TIP를 시작하면 f64 캐시를 읽는 자리에 tol을 덧칠하는 꼴이 된다.

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
4. 리미트 도달 → 뉴턴법으로 넘김: 해 수렴=교차 존재, 발산=미교차. 뉴턴 스냅백 점은 tol 있는 점이므로 **TIP에 통합**(수렴 잔차=tol, 애매하면 f128 재수렴).
5. 리미트에서도 확신 불가(작은 loop는 "다 찾았다"의 수학적 보장이 근본적으로 불가) → **`Unsupported` 정직 거부.** 조용히 틀리기보다 거부(§8·M7 철학).

**출발점 — M7 진입 시 이것부터 읽는다:** Li·Yang·Jia, "Advances and challenges in surface–surface intersection computation — An overview", Computer-Aided Design 193:104039, 2026. SSI 분야 전체 최신 개관이라 개별 논문 여러 개보다 이 리뷰가 최적 출발점. 그 시점의 최신을 반영해 위 세 계보(거리/법선/winding number)를 재비교 후 채택. (관련: Li·Jia·Chen, "Fast Determination and Computation of Self-intersections for NURBS Surfaces", ACM TOG 44(2), 2025 — ④ 자기교차 판정·거부용.)

**TIP와의 관계 정리 (혼동 방지).** TIP(Toleranced Indirect Predicates)는 "다항식 판정 + 회전 tol 필터"라 **부호 판정** 층이다. M7에서 TIP가 닿는 곳은 (a) 메시 조합 판정(내/외 분류)의 필터, (b) 뉴턴 스냅백 점의 tol 추적·f128 상승 — 둘 다 **정밀화·판정**이다. M7의 도박인 **SSI(교차를 찾는 것)** 자체는 tol 문제가 아니라 위상 존재 문제이므로 TIP 밖이다. "찾은 것을 정밀하게"(TIP·스냅백)와 "못 찾은 것을 찾기"(SSI)는 다른 일. SSI 성공 후라야 스냅백·판정이 의미 있고, SSI 실패 시 `Unsupported`.

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


## 10. 셀 기록 (오른 사다리 — 반증된 예측 보존)

M5 다면체 불리언 사다리를 셀 단위로 오른 기록이다. **★는 반증된 예측**이며, 이 문서의 최대 자산이므로 지우지 않는다 — 다음 사람이 같은 길을 다시 걷지 않게 한다. 셀 순서는 대체로 시간순이나 의존 관계로 엇갈릴 수 있고, 각 항목의 `(셀 3f-4)`류 꼬리표가 진짜 좌표다. 상호참조는 셀 ID·항목 제목으로 가리킨다(행 번호는 편집에 취약해 걷어냈다).

**볼록 경로 vs ray casting 경로 통합 여부 (M5-d 사다리 완료 후 프로파일링으로 결정).** M5-d1에서 비볼록 ray casting(`point_in_solid`)이 도입됐고, 원리적으로 볼록도 커버 가능하다(볼록 = 비볼록의 특수 경우). 현재는 `boolean()` 셀렉터가 **볼록 → 옛 half-space 경로**, **비볼록 → ray casting 경로**로 분기한다(볼록 회귀 0 보장). **결정 시점: M5-d 서브유닛 2~5 완료 후** — (a) 서브유닛 4에서 `classify_vertex` 은퇴 예정이라 볼록 경로 잔여분이 그때 드러나고, (b) 비볼록 경로가 완성돼야 볼록 완전 대체 가능성을 판단할 수 있으며, (c) 볼록 입력에서 half-space vs ray casting 속도 프로파일링이 필요하다. **판단 기준:** 속도 차이 유의미 + 볼록 흔함 → 두 경로 유지(fast path); 미미 or 유지보수 부담 → ray casting으로 통합. "측정 후 최적화" 원칙(§ SmallVec·reachability 캐시 결정과 동일).

**★ (5b)가 이 질문에 답했다 — 속도가 아니라 커버리지·정확성이 답했다.** fast path의 전제 "같은 답을 더 빨리"가 깨진다: 볼록 경로는 일반 경로가 답하는 입력을 거절하고(`cube_and_notch`·드릴 큐브를 `POKE_THROUGH`로), 진리를 f64로 정한다(`classify_vertex`의 tolerance, `enter_face`의 `t` argmax). 그래서 통합이 아니라 **볼록 `Fuse`/`Cut` 경로 삭제**다. (a)의 예고는 한 셀 늦었다 — `classify_vertex`는 서브유닛 4가 아니라 (5b)에서 죽는다. 속도 프로파일링은 (5b-0)이 했고(위), 필터 뒤 씨임 경로는 볼록의 ~21–27배이되 그 잔차는 exact 산술이 아니라 광선 순회라 (5d)의 몫이다.

**hollow(cavity 보유) 피연산자 — 임시 거절 적용, 완전 처리는 서브유닛 5.** 볼록 경로(`common`/`fuse_cut`)는 "솔리드 = 자기 outer 면 half-space의 교집합"이라 가정하는데 cavity가 이 가정을 깬다. `is_convex`는 틀리지 않았다 — 그것은 "outer shell이 볼록한가"라는 참인 명제를 정확히 답하고, 볼록 경로가 그 답을 "솔리드가 볼록하다"로 잘못 읽는다. **실측:** `[0,3]³`에서 `[1,2]³`를 파낸 hollow 박스(부피 26)에 코너 박스를 `Cut` → `boolean`이 `Ok`를 반환하고 부피 `26.875`(정답 `25.875`), `cavities=0`, **`validate`도 클린**. 거절 없는 오답 — 이 프로젝트가 가장 경계하는 "plausible but wrong". → **`boolean()` 진입부에서 `tag::HOLLOW_OPERAND`로 정직하게 거절한다(입력만; 출력이 hollow인 것은 그대로).** 고칠 곳은 `is_convex`가 아니라 **cavity의 거절 또는 지원**이고, 그 결정은 서브유닛 5(퇴화 경화)에 속한다. 부작용: 이 가드가 `tunnel`을 공개 `boolean` 경로에서 도달 불가로 만든다 — 가드는 defense-in-depth로 남기고, 그 테스트는 `overlap_fuse_cut`을 직접 호출해 커버리지를 유지한다.

**inner loop를 가진 면 — 분류기는 고쳤고, 재구성 경로는 임시 거절 (M5-d3 셀 3a).** `face_points`가 `face.outer`만 읽어 `point_in_solid`·`segment_crosses_face`가 구멍을 메운 것처럼 취급했다. `prepare_face_split`이 `ImprintSketch`/`PocketOnFace`에서 이미 그런 면을 만든다. **실측(pocket된 단위 큐브, 공동 `[0.3,0.7]²×[0.5,1]`):** 공동 안의 점이 `Inside`로 오분류되고, 더 나쁘게는 **분류가 점마다 엇갈린다**(공동 안 박스의 8개 코너가 5 Inside / 3 Outside) — 옆으로 빠져나가는 광선은 뚜껑을 아예 안 만나기 때문. → `face_loops`가 outer 링과 각 구멍 링을 모두 돌려주고 전부 부채꼴 분할한다. 구멍 링은 외향 법선 기준 CW이며(관례가 아니라 매니폴드 조건 — 각 엣지가 인접 벽면 outer 링에서 반대 방향으로 쓰이고 `validate`가 `NonOpposedEdge`로 검사한다), 그래서 부호 붙은 교차가 저절로 상쇄된다. 새 tolerance도 ear-clipping도 없다.

재구성 경로는 구멍을 못 다뤘다(**셀 3f-5에서 해소** — 아래). **실측된 세 가지 파손:** `reconstruct_face`는 구멍을 잃고(`Ok`, 부피 0.96996 vs 0.916625, `validate` 위반 5건), `edge_incidence`는 rim 엣지의 인접 면을 하나만 봐서 `inc[1]`에서 **패닉**하며, `coincident_merge`는 `solid_local_faces`로 구멍을 버린다(`Ok`, 부피 2.0533 vs 2.0). → `overlap_fuse_cut`·`fuse_cut`·`coincident_merge` 진입부에서 `tag::INNER_LOOP_OPERAND`로 거절. 셀 3f-1은 결과에 구멍을 **내지만** 이 가드를 은퇴시키지 않았다 — 구멍 있는 면을 **입력으로 받는** 것은 별개의 일이고, 은퇴는 셀 3f-5다. 주의: `detect_coincident_interface`는 **교차-솔리드 반대법선** 공면 쌍만 세므로, 인터페이스가 아닌 면에 imprint한 볼록 솔리드는 스택으로 보여 merge 경로로 들어온다 — `common`은 `is_convex`+`has_coplanar_pair`가 이미 막아 가드가 불필요하다.

**`indirect_orient3d`는 implicit point를 하나만 받는다 — 그런데 seam 수집에는 그걸로 충분했다 (M5-d3 셀 3b에서 실증).** `nacre_predicates::indirect_orient3d(p: &ThreePlane, q, r, s)`는 explicit 점 셋을 받는다. 애초 예상은 "면당 arrangement의 조합 판정이 seam 정점 둘·셋을 얽으므로 Attene 2020의 LPI/TPI 계열 2·3-implicit 변형이 필요하고, 그때까지 arrangement 코어는 exact orient2d(투영 좌표) + 투영 불확실성 필터로 서야 한다"였다. **seam 세그먼트 수집 단계에 관해서는 틀렸다.** 실제로는 (a) 끝점이 `three_planes(P,Q,R)` — 이미 쓰는 3-평면 implicit point고, (b) "엣지 안 + 상대 면 안" 포함 판정은 `segment_crosses_face`가 exact orient3d로 이미 답하며, (c) 직선 위 두 교점의 전후는 `three_plane_orient3d(P,Q,R_i, R_j.tri)`(1-implicit)와 `det3_sign([n_P,n_Q,n_{R_j}])`(exact)의 곱으로 정해진다. **두 implicit 점이 한 술어에 동시에 들어가는 자리가 없다.** 그래서 2D 투영도, orient2d도, `ARRANGEMENT_DEGENERATE`의 "투영 불확실성 필터"도 필요 없었다(geom에 `plane_pair_dir_sign` 얇은 exact 래퍼 하나만 추가). 다중 implicit 술어가 정말 필요해지는 곳은 **셀 추출**(두 implicit 점의 orient)이며, 그 셀에 이르러 다시 판단한다 — 쓰지도 않을 술어를 미리 만들지 않는다. **그 셀은 3f-3이고, 사야 할 것은 정확히 하나다: 고리의 감김.** 고리들의 감김이 **모두 같으면 서로 포함하지 않는다** — 직접 포함은 깊이를 1 바꾸고 깊이 패리티가 재료의 안팎을 정하므로 감김이 반대가 된다. 그러므로 감김 하나가 중첩 판정과 hole/island 판정을 동시에 준다(그리고 `kept[0]`과 서로를 검산한다). 다중 chord(3e-2)는 조합만으로 풀렸으니 좌표는 여기서 처음 필요하다. **예약은 셀 3h에서 집행됐고, 예상의 절반은 틀렸다 — 아래.**

**`u_and_slab`은 다중 chord와 inner loop를 **둘 다** 요구한다 (셀 3e-1에서 실측).** 손계산은 닫힌 고리 0을 예측했으나 실제는 **2**다: U의 두 기둥이 슬래브의 `y=1.5` 면 **내부에** 사각형 둘을 뚫는다. `multichord`가 U의 base cap에서 먼저 발화할 뿐이고, 그 픽스처는 3e-2(다중 chord)만으로는 열리지 않는다. **셀 3e-2 이후 같은 입력이 `pokehole`로 거절된다** — 태그가 움직인 것이 호 쪽이 끝났다는 측정이다. **셀 3f-3이 나머지 절반을 갚아 세 방향(`4.0` / `12.0` / `6.7`)으로 열었다.** 픽스처마다 `(한 면당 최대 열린 호, 전체 닫힌 고리)`를 `arrangement_agrees_with_todays_seam_bookkeeping`이 못박는다 — 사다리의 어느 칸이 무엇을 여는지가 테스트에서 읽힌다.

**`ARRANGEMENT_DEGENERATE`는 미검증 백스톱이다 (`fourplane`·`seam_branch`·`seam_count_mismatch`·`loop_orient_mismatch`와 같은 줄).** `arrange::seam_segments_on`이 홀수 교차 카운트에서 거절한다(패리티상 짝수여야 한다 — 교점은 enter/exit로 교대하고 직선의 양 끝은 바깥). 그런 퇴화(정점 스침, 엣지 공선)는 그 전에 `CONTACT_DEGENERATE`가 잡을 공산이 크므로, 발화 테스트가 없다. **셀 3c에서 다시 쟀다: 여전히 미발화다.** 조립 골든 넷과 144-인스턴스 proptest(면 859건 통과·5건 거절) 어디서도 뜨지 않았다. **셀 3f-1에서 또 쟀다: 여전히 미발화다** — 불리언 결과가 구멍을 갖는 첫 형상까지 포함해서. (3f-1·3f-2·**3f-3 어느 것도** "대표점 분류"를 쓰지 않았다. 패리티 논증이, 그리고 3f-3에서는 **고리의 감김**이 대신했다. **여기 적었던 "대표점이 필요한 것은 3f-3"은 틀렸다** — 포함 판정이 실제로 필요한 곳은 **호 옆의 고리**를 배치하는 3f-4다. 3f-3에서 다시 쟀고 `ARRANGEMENT_DEGENERATE`는 여전히 미발화다.)

**seam 경로 조립 — 노드 동일성은 증명되고, `edge_seam`의 1:1은 가드에 얹혀 있다 (M5-d3 셀 3c).** `arrange::seam_paths_on`이 seam 세그먼트를 **끝점 triple의 동일성만으로** 이어 붙여 열린 호와 닫힌 고리를 낸다. `∂Y ∩ P ∩ f`는 1-manifold이므로 세그먼트는 끝점에서만 만나고, 그 자리는 Y의 엣지가 `f`를 뚫는 지점이며, 그 엣지를 공유하는 두 Y-면이 양쪽에서 같은 정렬 triple을 낸다. **좌표를 한 번도 읽지 않는다.**

- **노드 키 충돌은 불가능하다(증명).** 경계 노드의 triple은 X-평면 둘, 내부 노드는 Y-평면 둘 — 두 종이 겹칠 수 없다. 서로 다른 `g`는 서로 다른 `Q`를 낸다. 같은 `g` 안에서 세 번째 평면 `R`이 중복되는 유일한 길(비볼록 면의 **공선 엣지 둘**이 같은 이웃 평면 위에 있는 경우)은 `order(R,R) == 0`이라 `fourplane`이 이미 거절한다. 따라서 차수는 1 또는 2뿐이다.
- **`SEAM_BRANCH`(차수 ≥3)는 오늘 도달 불가이지만 죽은 코드가 아니다.** 위 증명의 세 번째 갈래가 `fourplane`에 기대므로, 셀 3e/3h가 그 거절을 완화하면 한 노드에 세 세그먼트가 붙을 수 있다. `tunnel`과 같은 종류의 살아 있는 백스톱.
- **`STRICTARC`의 전제가 사라졌다 — 셀 3e-1에서 은퇴함.** 옛 `reconstruct_face`는 내부 bend를 직선 chord에 투영해 정렬했고, 그 정렬이 접힌 호에서 틀릴 수 있어 `strict`가 거절했다. 인접 순회는 **볼록성과 무관하게 참 호 순서**를 낸다. 실증: `l_and_popup_box`의 박스 바닥면에서 조립된 호가 `(2,0.5)→(2,1)→(1,1)→(1,1.5)`(계단형, 반대로 꺾이는 bend 둘).
- **`edge_seam`의 "엣지당 seam triple 하나"는 타입이 아니라 가드가 지킨다.** 한 엣지가 두 번 교차되면 양 끝점이 같은 쪽이고, 그 분기도 `pierced_face`를 부르므로 `pierced_multi`가 뜬다. straddle 엣지가 세 번 교차돼도 마찬가지. 볼록 경로는 볼록성이 교차를 ≤2로 묶고, 양끝-바깥은 `segment_enters` → `poke_through`. **실측:** `cube_and_notch`(`A=[0,10]³`, `Y=[3,7]×[−1,1.4]×[−1,1.2]`)는 `poke_through`로 거절되고, 같은 입력에서 `seam_paths_on`은 두 경계 노드가 **A의 같은 엣지 위에 있는** 호 하나를 옳게 낸다. → **셀 3d는 splice를 `Handle<Edge>`가 아니라 seam triple로 키잉해야 하고, 한 엣지 위 두 교차의 전후는 `order_along(P, R, ·, ·)`로 정한다**(엣지의 두 평면을 쌍으로 넘기면 되므로 새 술어가 필요 없다). 셀 3e가 그 가드를 걷어내기 **전에** 3d가 끝나 있어야 한다. 순서를 뒤집으면 교차 하나가 조용히 사라지고 `validate`가 통과하는 틀린 면이 나온다. `an_edge_crossed_twice_is_rejected`가 그 가드를 태그로 고정해 먼저 운다. (**셀 3e-3이 비볼록 경로에서 그것을 집행했다** — 예고한 `order_along(P, R, ·, ·)`이 정확히 쓰였고, 순회 방향은 `edge_sign`이 이미 갖고 있었다. 볼록 경로의 `edge_seam`은 남는다.)
- **`insert` 충돌 `debug_assert`는 넣지 않았다 — 공허하기 때문.** `edge_incidence`가 엣지를 중복 제거해 각 엣지를 정확히 한 번 돌려주고 두 솔리드의 엣지 집합은 서로소이므로, `edge_seam.insert`의 이전 값은 입력과 무관하게 항상 `None`이다.
- **경계 끝점의 엣지는 `{P,R}`로 역추적하지 않는다.** 비볼록 면은 같은 이웃 평면 `R` 위에 공선 엣지 둘을 가질 수 있어 엉뚱한 엣지에 붙는다. 생산 시점에 `he.edge`를 `SeamSegment.on_edge`에 기록한다. 내부 노드가 `None`인 것은 정보 손실이 아니다 — 그 점은 상대 솔리드의 면에서 **경계 노드**로 나타나며 거기서 그쪽 엣지를 기록한다.
- **결정성은 자료구조가 아니라 순회 규율에서 나온다.** 인접 맵은 조회에만 쓰고 절대 순회하지 않는다(`HashMap` 순서는 인스턴스마다 다르다). 경로는 가장 낮은 인덱스의 미사용 세그먼트에서 시작한다. 이게 깨지면 골든이 산발적으로 깨지는 데 그치지 않고 셀 3d에서 **연산 로그 재생**(§2 절대 원칙)이 무너진다. `seam_paths_are_deterministic`이 고정한다.
- **검증하지 못한 전제 둘 (정직하게 남긴다).** (a) **triple 동일성 ≠ 기하 동일성**: 서로 다른 triple이 같은 점에 놓이면 경로가 갈라진다. `fourplane`과 `segment_crosses_face`의 graze 검사가 막는다고 보지만 확인하지 않았다. 좌표 근접 검사로 `debug_assert`를 만들지 않는다 — tolerance가 필요하고, tolerance는 `Origin::Discovered` 밖에 설 자리가 없다. (b) **Y의 엣지가 `f`의 평면 `P` 안에 통째로 누우면** 접선 접촉이다. `has_coplanar_pair`는 *면*의 공면만 막는다. `CONTACT_DEGENERATE`가 거절한다고 보지만 확인하지 않았다. (**셀 (5a)가 확인했고 이름을 줬다** — `VERTEX_ON_FACE_PLANE`. 옛 부채꼴도 같은 입력을 거절하고 있었다: `segment_triangle_cross`가 평면 검사를 먼저 하므로 끝점이 평면 위면 모든 삼각형이 `Degenerate`다. 커버리지는 변하지 않았고 태그만 정직해졌다. 전제 (a)는 여전히 미검증이다.)
- **검증의 세 겹 중 두 겹만 걸 수 있다.** 이 셀은 솔리드를 만들지 않으므로 OCCT diff도 부피 항등식도 없다. 대신 (1) 손계산 골든 넷, (2) 독립 오라클 — 차수-1 노드 수 == `reconstruct_face`의 kept/dropped 전이 수(`point_in_solid` 광선 캐스팅으로 테스트가 직접 센다; 전이 수는 `keep` 선택과 무관하므로 `BoolKind`를 고르지 않는다), (3) 구조 불변식 proptest. 오라클의 전제는 "한 엣지가 두 번 교차되지 않음"이고, `cube_and_notch`가 정확히 그 반례라 **테스트로 남겨 도메인 경계를 문서화한다**. 세 겹은 셀 3d(동등 치환)에서 복원된다 — 그때 기존 OCCT diff 17건이 그대로 걸린다.

**arrangement 배선 — 비볼록 경로만, `edge_seam` 소멸 (M5-d3 셀 3d).** `overlap_fuse_cut`이 면을 `arrange::seam_paths_on`으로 재구성한다. 출력과 거절 태그가 전부 불변이고, 3b·3c에서 걸 수 없던 **OCCT diff 17건이 복귀**했다(직접 게이트는 `nonconvex_overlap_{cut,fuse}_matches_occt`).

- **볼록 경로는 배선하지 않았다 — 측정된 차단 요인.** `two_boxes`(`A=[0,1]³`, `B=[0.5,1.5]³`)에서 B의 `y=0.5` 면은 A의 `x=1` 면을 `(1,0.5,0.5)`, 즉 **그 정사각형의 정중앙**에서 가로지른다. 두 대각선이 거기서 만나므로 `segment_crosses_face`가 네 apex 전부에서 스쳐 `contact_degenerate`로 정직하게 거절한다(위 항목). 배선하면 `cut_of_two_cubes`와 그 OCCT 골든이 깨진다. **후퇴가 아니다:** 볼록 경로의 "엣지당 교차 하나"는 가드가 아니라 **볼록성이라는 기하학적 사실**이다(선분은 볼록 경계를 최대 2회 만나고 straddle이면 정확히 1회). 지뢰는 비볼록 경로에만 있었고, 거기서 `edge_seam`은 **재키잉이 아니라 삭제**됐다 — 경계 교차가 생산 시점의 `SeamSegment::on_edge`를 들고 다니므로 `Handle<Edge>`로 splice를 키잉할 이유가 없다. 두 경로 통합은 `segment_crosses_face`의 부채꼴 퇴화 해소(서브유닛 5)에 걸려 있으며, 그건 위 "볼록 경로 vs ray casting 경로 통합"과 **같은 게이트**다. (**셀 (5a)가 그 게이트를 열었다.** `the_centre_of_a_square_face_is_pierced_not_grazed`가 `two_boxes`를 비볼록 경로로 직접 통과시킨다 — `Cut` 부피 `0.875`, `validate` 클린. 배선 자체는 셀 (5b).)
- **`transitions`는 버리지 않고 남겼다.** `kd`/`dk`는 여전히 거기서 나온다. 호의 양 끝에서 뽑으면 두 끝이 **같은 엣지**에 앉을 때(3c의 `cube_and_notch`) 무너진다. 경로로 옮긴 것은 **분류뿐**이고, 3e가 확장할 자리가 바로 거기다.
- ~~**가드 우선순위는 `MULTICHORD` → `POKEHOLE` → splice.**~~ **셀 3e-2에서 거짓이 됐다** — 비볼록 경로의 `MULTICHORD`가 사라져 우선순위 자체가 없다. 호 옆의 고리는 이제 `POKEHOLE` 하나가 잡는다.
- **`SEAM_COUNT_MISMATCH` — 3c의 내부 오라클을 프로덕션 불변식으로 승격.** 처음엔 `debug_assert`로 두려 했으나 **그럴 수 없다**: `n_open`으로 분류한 뒤 `transitions[0]`을 인덱싱하므로 릴리스에서 인덱스 패닉이 된다. 그리고 이 검사는 방어가 아니라 **서로 독립인 두 기계의 합의**다 — `classof`의 광선 캐스팅과 arrangement의 정확한 교차, 그리고 `pierced_face`가 만든 `seam` 목록과 arrangement의 경로. **개수가 아니라 집합 상등**으로 검사한다(개수는 한 chord의 두 끝이 같은 엣지에 앉아도 통과한다). 오늘 미발화 백스톱이며(`fourplane`·`arrangement_degenerate`·`seam_branch`와 같은 줄), 죽은 코드가 아니다: 3e가 `pierced_multi`를 완화하면 가장 먼저 울고, 그 자리는 거절이 아니라 **다중 chord의 진짜 처리**가 되어야 한다. (**셀 3e-3이 그 자리를 채웠다.** 집합 상등이 죽고 run의 **교대 vs `classof`** 가 섰다. 정점 없는 run은 교대만이 정한다 — 예고가 정확했다.)
- **seam에 닿지 않는 면에는 arrangement를 부르지 않는다.** 오늘의 `transitions == 0` 가지 그대로다. 새 코드는 `segment_crosses_face`를 부르는데(옛 코드는 부르지 않았다), 그러면 seam과 무관한 먼 면의 스치는 접촉이 `contact_degenerate`를 내서 **거절 태그를 갈아치울 수 있다**. 비용도 seam에 닿는 면으로 한정된다(`O(F_y·E)`, 나머지는 `O(1)`).
- **arrangement는 기하의 출처가 아니다.** `SeamEnd::point`는 캐시다. 모델에 들어가는 정점은 `seam[..]`의 `point`와 `tol`이다. 셀 3e-1이 `strict`를 지운 뒤로 **`reconstruct_face_paths`는 좌표를 하나도 읽지 않는다** — triple과 핸들이 들어가고 triple과 핸들이 나온다.
- **오늘의 조용한 오답 하나가 닫혔다.** 옛 `bends` 필터는 **평면 위 모든 seam 정점**을 담았다. 그래서 한 면에 chord 하나(`transitions == 2`)와 내부 고리가 동시에 있으면 고리의 노드가 bends에 섞여 **매니폴드이지만 틀린 면**이 나왔다(`POKEHOLE`은 `transitions == 0`일 때만 봤다). 새 코드는 `POKEHOLE`로 거절한다. 단일 볼록 피연산자로는 만들 수 없어 스위트에 없었다.
- **배선 전에 측정했다.** `arrangement_agrees_with_todays_seam_bookkeeping`이 `reconstruct_face`에 도달하는 다섯 픽스처(받아들여지는 셋 + **거절하는 셋**) × 양쪽 솔리드 × 모든 면 = 72면에서 (1) 어떤 면도 거절하지 않고, (2) 경계 노드가 transition 엣지와 **집합으로** 같고, (3) 받아들여지는 입력엔 고리가 없고, (4) 단일 chord의 **bend들**이 chord 방향으로 **강한 단조**임을 확인한다. (4)가 치환의 출력 불변을 증명한다 — 오늘의 `sort_by`는 안정 정렬이라 동률이면 `seam` 삽입 순서로 떨어지고, 그건 arrangement가 모른다. 거절 픽스처를 포함한 이유: 새 코드가 **다른 면에서 먼저** 거절하면 태그가 바뀌는데 받아들여지는 입력만 재서는 그걸 못 본다.
- **측정이 가르친 것: 반사 bend는 chord 구간 밖으로 투영된다.** `l_and_reflex_box`에서 bend들은 끝점 `0.0`·`0.72`에 대해 `−0.24`와 `0.96`에 투영된다. "bend는 두 끝점 사이에 놓인다"는 거짓이고, 참일 필요도 없다 — 오늘 정렬되는 것은 끝점이 아니라 bend뿐이다.

**`segment_crosses_face`는 관통점이 면의 중심에 놓이면 모든 apex에서 스친다.** 사각형 면의 두 대각선이 중심에서 만나므로 fan을 어느 꼭짓점에서 시작해도 관통점이 대각선 위에 놓인다 ⇒ `CONTACT_DEGENERATE`. 정직한 거절이지만, **정육면체 두 개를 축정렬로 겹치면 실제로 걸린다**(M5-d2의 `l_and_corner_box` 주석이 이미 비대칭 좌표를 고른 이유). 픽스처는 변을 서로 다르게 잡는다. ~~근본 해소는 fan 대신 실제 삼각분할(또는 대각선을 피하는 apex 선택)이며 서브유닛 5(퇴화 경화)의 항목이다.~~ **이 예측은 틀렸다 — 셀 (5a)를 볼 것.** 어떤 삼각분할이든 대각선을 만들고, 정사각형의 정중앙은 **모든** apex의 대각선 위에 있다. 부채꼴을 고치는 것이 아니라 없애는 것이 답이었다.

**평면 삼각분할 — 결함이 기록보다 컸고, 3f 앞에 독립 셀로 수선했다 (tess 셀).** §9는 원래 "`triangulate_planar`가 inner loop를 무시해 구멍을 메운다 … 셀 3f에 묶어 고친다"고만 적었다. 3f-1 계획을 검토하다 **세 가지**가 드러나, 약속을 사실에 맞게 고치고 별개 셀로 처리했다.

- **부채꼴은 구멍만이 아니라 비볼록 외곽 링에서도 틀렸다.** `ring[0]`에서 부채꼴로 자르는 것은 링이 **그 정점에서 별 모양일 때만** 옳다. L-프리즘 캡은 `(0,0)`에서, 3e-1이 만든 계단 면은 `(0.5,0.5)`에서 별 모양이라 **우연히** 맞고 있었다. U-프리즘 캡은 아니다 — 부채꼴 변이 노치를 지나 삼각형이 **뒤집혀** 나온다(실측: 캡마다 넓이 `2.6`이 부호 없이 더해져 `+5.2`).
- **`tessellate`와 `to_obj`가 supersede된 죽은 면을 뱉었다.** `Store`는 append-only이고 불리언·pocket·pad는 피연산자를 지우지 않는데, 둘 다 면 저장소를 통째로 순회했다. 실측: stacked-fuse에서 죽은 큐브 둘의 표면적 `+12.0`. 이제 `live_solids`에서 도달 가능한 셀만 메시로 만든다(`Reachable`은 `HashSet`이므로 **저장소 순서로 순회하며 멤버십만 묻는다** — 재생 가능성).
- **★ 사용자가 보는 OBJ는 `tessellate`를 거치지 않았다.** 외부 호출부는 부트스트랩 `to_obj(&Model)`뿐이고, 그건 `face.inner`를 아예 무시하는 **두 번째 부채꼴**이었다. `triangulate_planar`만 고쳤다면 생명줄은 하나도 안 고쳐졌을 것이다. 이제 둘이 `tess::polygon::triangulate_polygon`을 공유한다. (파이프라인 통합 자체는 남았다 — `to_obj`의 "OBJ 인덱스 == 정점 핸들 인덱스" 계약을 바꾸므로 별개 셀이다.)

**게이트: `Σ|삼각형 넓이| == props.area`.** 부채꼴 버그와 구멍 버그를 **동시에** 잡는 등식이다. **부호 있는 합으로는 못 잡는다** — 다각형 밖으로 삐져나간 삼각형은 뒤집혀 있어 상쇄되고 shoelace가 자기 자신과 맞아버린다. 그리고 `watertight`(무방향 변 다중도 2)는 **구멍 버그는 잡지만 부채꼴 버그는 못 잡는다**(부채꼴도 변은 짝을 이룬다). 검사마다 담당이 다르다.

**정직한 실패.** `tessellate`·`to_obj`가 `Result<_, TessError>`를 낸다. 옳은 메시가 없으면(자기교차 링, 감김이 뒤집힌 구멍, 브리징 불가) 에러다 — **부채꼴로 되돌아가지 않는다.** 그 정직함이 곧바로 값을 했다: 새로 쓴 pocket 픽스처가 구멍을 뚜껑 밖에 반쯤 걸쳐 놓았는데(`PocketOnFace`의 프로파일은 면에서 유도된 좌표계에 산다), `validate`는 위상만 보아 통과시켰고 `props`는 링을 그냥 적분했으며 옛 부채꼴은 `face.inner`를 무시했다. **`NoEar`만이 찾아냈다.**

**Newell 법선이 `Orientation` 가정을 없앴다.** `tessellate`의 doc은 "`Orientation::Forward`를 가정한다"고 적혀 있었으나 `assemble_fuse_cut`의 `flip`이 이미 `Reversed` 면을 만든다. 링에서 법선을 뽑으므로 표면 법선을 볼 필요가 없다. 부수적으로 알게 된 것: **링은 언제나 자기 Newell 법선에 대해 CCW다.** 뒤집으면 법선이 뒤집힐 뿐이므로 "시계방향 외곽 링"은 탐지할 수 없고, 탐지하려 들지도 않는다 — 어느 쪽이 바깥인지는 호출자가 안다.

**메시는 캐시이므로 `f64`를 쓴다.** 정확한 기하가 진실이라는 원칙과 충돌하지 않는다.

**곡면 tess는 적응적(곡률 기반 해상도)으로 — M6 곡면 도입 시 검토.** 현재 평면 tess는 귀 자르기로 이미 최소 삼각형(`N + 2h − 2`)을 내므로 적응적 이득이 없다. 그러나 M6에서 곡면(구·원통·원뿔)이 들어오면, 곡면을 삼각형으로 근사할 때 "평평에 가까운 곳은 큰 삼각형, 곡률이 심한 곳만 잘게" 나누는 적응적 삼각화가 **시각 품질을 유지하면서 tess 메시 용량을 줄인다.** FEM용 정삼각형이나 프린터용 방수 같은 특수 용도가 아니라 범용 개선(용량↓ + 매끄러움 유지)이라, 앱이 아니라 **커널 기본 tess로 넣는 것을 고려할 만하다.** 비용도 작다 — 바이너리 용량 증가는 미미(코드 수백 줄), 연산 시간은 삼각형 수 감소와 곡률 판단 비용이 상쇄되어 비슷, 대신 tess 결과 용량이 작아져 커널 내 형상 확인 시에도 메모리가 절약된다.

단 이건 **"측정 후 최적화" 원칙**에 따라 M6 곡면 tess를 실제 설계할 때 결정한다 — 그때 균일 vs 적응적을 곡면 실제 모델로 측정해서, 용량·속도 이득이 복잡도를 정당화하는지 확인하고 채택. 지금(평면 M5)은 균일 귀 자르기로 충분하니 **구현·상세설계 금지.**

*그때의 출발점으로, 오늘 이미 참인 것:* `circle_segments`는 이미 sagitta(현 편차) 기반이라 `TessConfig::tol`과 반지름으로 분할 수를 정한다. 다만 원의 곡률은 **일정**하므로 그건 적응적이 아니라 파라미터 균일이고, 그것으로 충분하다. 원통도 축 방향 곡률이 0이라 `triangulate_cylinder`가 내부 샘플을 아예 두지 않는다. **적응적 이득이 실제로 나오는 것은 곡률이 변하는 곡면이다** — 원뿔(꼭짓점 근처), 구(곡률은 일정하나 uv 격자가 극에서 과다 샘플링), 토러스·NURBS. 그러므로 M6에서 잴 대상은 "곡면 일반"이 아니라 **그 목록**이고, 원통·원은 이미 최적이므로 비교 기준선으로 쓴다.

**불리언 결과가 구멍을 가질 수 있다 — `POKEHOLE` 축소 (M5-d3 셀 3f-1).** `l_and_dimple`(L-프리즘 + 그 윗면에 서 있는 스텁 `[0.3,0.7]²×[0.5,1.5]`)의 seam은 면 **내부의 닫힌 고리**다. 이제 그 고리가 `Face.inner[0]`으로 나간다. `Cut` = 눈먼 포켓(2.92), `Fuse` = 보스(3.08), 둘 다 `validate` 클린이고 같은 면에 같은 구멍을 낸다.

- **여는 모양은 정확히 하나이고, 그 증명은 짧다.** seam을 건널 때마다 분류가 뒤집히므로, 경계가 **전부 유지**(`kept[0]`)인 면 위의 **고리 하나**는 반드시 **버려지는 내부**를 감싼다 — 구멍이지 섬이 아니다. 넓이도 2D 포함 판정도 필요 없다. 나머지 셋은 `POKEHOLE`이 계속 잡는다: 고리 둘(**중첩**일 수 있다 — 구멍 속의 섬), `kept[0] == false`(고리 내부가 유지 영역인 **섬 면**), 열린 호와 고리의 공존. **가릴 수 없는 것은 열지 않는다.**
- **★ 고리의 방향은 부호 세 개의 곱이다 — 좌표를 읽지 않는다.** b-rep의 모든 고리는 **재료를 왼쪽에** 두고 돈다. 평면에서 그 조건은 `n_out_f × t ∝ ±n_out_g`이고, `n_out = s·n_plane`(`s = orient_sign`, `Reversed` 면에서 `−1`)을 넣고 전개하면

  ```
  t  ∝  −keep · s_f · s_g · d ,      d = n_P × n_Q       (keep = +1 ⇔ Outside)
  ```

  `d`는 `seam_segments_on`이 교차를 **이미 정렬해 둔 그 방향**이고, 이웃 두 노드가 `d`를 따라 앞뒤인지는 `order_along(planes, p, q, r_i, r_j)`가 정확히 답한다. `arrange::orient_seam_loop`은 이 규칙으로 **정렬된 고리를 돌려준다** — `-> Result<bool>`로 두면 "그 값을 보고 뒤집어라"는 의무가 호출자에게 남고, 잊어도 아무도 안 잡는다. **방향을 아는 쪽이 방향을 적용한다.**
- **`LOOP_ORIENT_MISMATCH`는 `debug_assert`가 아니라 진짜 거절이다.** 고리의 모든 변이 같은 답을 내야 한다(되감기 변 `n−1 → 0`을 포함해서 — 빠뜨리기 쉽다). 이건 방어가 아니라 **서로 독립인 두 기계의 검산**이다: 정확한 술어 `order_along`과 방향 장부 `PlaneInfo::{n_out, orient}`. `SEAM_COUNT_MISMATCH`와 같은 종(種)이고, 같은 이유로 릴리스에서 살아남는다. 공유 평면 추출 실패와 `order_along == 0`도 같은 태그로 묶인다 — `expect`로 뽑으면 릴리스 패닉이다. **오늘 미발화 백스톱**이며, 3c의 노드 동일성 논증이나 `FOURPLANE`이 완화되면 운다.
  - `seam_paths_on`의 순회 기록을 재사용해 `order_along` 재계산을 피할 수도 있었다. **그러면 안 된다** — 한 기계의 기억을 믿으면 검산할 것이 없다. **재계산이 곧 두 번째 기계다.**
- **전역 부호 뒤집힘은 커널이 못 잡는다.** 모든 변이 일관되게 틀리기 때문이다. 그건 `validate`(위상: `NonOpposedEdge`)와 `tessellate`(기하: `HoleWinding`)의 몫인데 **둘 다 사용자가 부를 수도, 안 부를 수도 있는 검사다.** 그래서 (1) 배선 **전에** `orient_seam_loop`의 순수 골든이 전역 부호를 못박고(`material_outside`를 뒤집으면 결과가 정확히 역순이어야 한다 — 규칙의 유일한 자유도가 그 불리언이므로), (2) `a_flipped_hole_loop_is_caught`가 **이 형상에서** 두 검출기가 실제로 우는지 확인한다. `Store`에 `get_mut`이 없으므로(append-only) 뒤집힌 고리를 가진 새 `Face`를 push하고 새 `Shell`·`Solid`로 갈아끼운 뒤 `live_solids`를 옮긴다 — 옛 면은 `reachable()` 밖으로 떨어진다.
- **부피와 OCCT는 방향에 무감각하다** (`props`가 `|넓이|`를 합한다 — 더 정확히는 **법선을 `Face.orientation`에서 뽑고 링을 읽지 않는다.** 셀 3f-2에서 소스로 확인했고, 그래서 눈멂이 구멍에 국한되지 않는다). 대신 `blind_dimple_{cut,fuse}_matches_occt`는 **넓이도** 비교한다 — 우리 게이트(`props`·메시 넓이 합·손계산 14.8)는 모두 **같은 링을 적분**하므로 구멍의 크기에 대해서는 서로를 검산하지 못한다. 독립 커널만이 그 주장을 세운다. OCCT diff 20 → 22.
- **`flip`과 `inner`의 상호작용은 쓰였으나 미검증이었다 — 셀 3f-5에서 검증됐다.** `Cut`의 B-조각은 모든 고리를 뒤집는데, 당시 스텁 쪽 면에는 고리가 없었고, B의 면이 고리를 가지려면 구멍 있는 피연산자가 필요한데 `INNER_LOOP_OPERAND`가 문 앞에서 거절했다. **가드가 문에 있으면 그 뒤의 코드는 실행되지 않는다.** `cut_slab_by_pocket`이 이제 그것을 밟는다.
- **`INNER_LOOP_OPERAND`는 그대로 남는다.** 결과에 구멍을 **내는** 것과 구멍 있는 면을 **입력으로 받는** 것은 다른 일이다 — 3f-1의 결과를 다시 불리언에 넣으면 정직하게 거절된다. 은퇴는 **3f-5**.

**메시 게이트에 부호 있는 부피를 더했다 — 방향을 보는 두 번째 기계 ((gate) 셀).** 3f-2를 검토하다 그물의 구멍이 드러났고, 추측 대신 **소스를 읽어** 확정했다.

- `props::face_contribution`은 면 법선을 **`Face.orientation`에서** 뽑고 `polygon_area_centroid`는 `|A_vec|`(부호 없는 크기)을 돌려준다. ⇒ `props`는 **어떤 고리 뒤집힘도 읽지 않는다.** 구멍이든 외곽이든. 부피도, 그것으로 채점하는 모든 OCCT diff도 눈멀었다.
- `tess::triangulate_planar`는 `face.orientation`도 표면 법선도 **읽지 않는다.** 삼각형 감김이 링을 따라간다. `HoleWinding`은 inner loop만 본다.

그래서 면의 링을 보는 기계가 `validate` **하나뿐**이었다. 놓치는 결함이 있어서가 아니라 — 뒤집힌 고리는 공유 엣지를 같은 방향으로 두 번 쓰므로 `NonOpposedEdge`가 반드시 뜬다 — **기계가 하나**여서 문제다. `SEAM_COUNT_MISMATCH`와 `order_along` 재계산에서 지켜온 기준에 미달한다.

**`mesh_volume == props.volume`이 두 번째 기계다. 순환이 아니다.** `props`는 `orientation`을 읽고 링을 안 읽는다. `tess`는 링을 읽고 `orientation`을 안 읽는다. 그러므로 두 부피의 일치는 정확히 **"모든 면에서 `Face.orientation`이 제 링의 감김과 일치한다"**는 진술이다. `validate`가 위상(엣지 대립)으로 묻는 것을 기하(표면 발산)로 묻는다.

| 검사 | 잡는 것 | 못 잡는 것 |
|---|---|---|
| watertight | 잃은 면, 공유되지 않은 rim | 다각형 밖으로 나간 부채꼴, 뒤집힌 면 |
| `Σ\|넓이\| == props.area` | 부채꼴, 메워진 구멍 | 뒤집힌 면 (부호가 없다) |
| `mesh_volume == props.volume` | 뒤집힌 면·고리 | 넓이만 틀린 것 |

여섯 픽스처 전부에서 **먼저 쟀다**: 최대 불일치 `1.3e-15`. `only_the_signed_volume_sees_a_reversed_face`가 이빨을 보인다 — 정육면체 한 면의 outer loop를 뒤집으면 삼각형 수·watertight·부호 없는 넓이·`props.volume`이 **전부 침묵**하고 부호 있는 부피만 `⅔`(= `x=1` 면 flux의 두 배) 어긋난다. 정육면체는 `[1,2]³`다 — `[0,1]³`이면 원점을 지나는 세 면의 기여가 `0`이라 그중 하나를 뒤집어도 **공허하게 통과**한다.

*미래 후보(구현 금지):* 이 검사를 `nacre-validate`의 진짜 `Violation`으로 승격할 수 있다. 오늘은 테스트 게이트뿐이라 라이브러리 사용자는 여전히 부를 수도 안 부를 수도 있다. 반쪽엣지 하나짜리 원 고리에서 Newell 법선이 정의되지 않는 문제를 먼저 풀어야 한다.

**면이 섬(island)일 수 있다 — `POKEHOLE` 두 번째 축소 (M5-d3 셀 3f-2).** `Cut(stub, L)`이 열린다: `l_and_dimple`의 피연산자를 뒤집으면 같은 seam 고리가 **경계가 통째로 버려진** 면에 앉는다. 그 면의 유지 영역은 고리의 내부뿐이고 `∂f` 조각이 하나도 없다 — outer loop가 곧 seam ring, 정점 넷 전부 `Discovered`다. 답은 `z=1` 위의 `0.4×0.4×0.5` 박스(0.08), `validate` 클린.

- **증명이 짧아졌다.** 호 없이 고리 하나면 면이 원반 + annulus로 갈리고 `∂f`는 annulus에 있다. 그러니 `kept[0]`이 annulus를 분류하고 원반은 반대 분류를 갖는다. annulus가 살면 **구멍**, 원반이 살면 **섬**. 제3의 답은 없고 넓이는 필요 없다. 빠져나갈 구멍("정점은 다 같은 분류인데 엣지가 두 번 교차")은 그 엣지가 `bnd`에 들자마자 빈 `transitions`와 **집합으로** 달라져 `SEAM_COUNT_MISMATCH`가 먼저 운다 — 3d가 개수 대신 집합을 택한 값.
- **새 부호도, 새 술어도, 새 태그도 없다.** 부호 규칙은 국소적이다: seam 위의 한 점에서 유지 재료가 **어느 쪽인지**만 묻고, 그 고리가 무엇을 **감싸는지**는 묻지 않는다. 그래서 구멍을 CW로 감던 바로 그 `orient_seam_loop` 호출이 섬을 CCW로 감는다. 3f-1의 골든이 `inside == rev(outside)`로 **두 방향을 다 못박아 둔** 것이 한 셀 뒤에 값을 했다. `assemble_fuse_cut`도 손대지 않았다 — `ring` 클로저는 노드 종류를 묻지 않는다.
- **이름이 거짓이 되어 개명했다.** `orient_hole_loop` → `orient_seam_loop`, `HOLE_ORIENT_MISMATCH` → `LOOP_ORIENT_MISMATCH`. doc은 이미 일반적이었다("outer or inner"). ~~`POKEHOLE`은 개명하지 않는다 — 3f-3이 비볼록 발화 지점을 은퇴시키면 볼록 경로의 문자 그대로의 poke-through만 남아 이름이 **다시 정확해진다.**~~ **(5b)에서 뒤집혔다:** 볼록 경로가 통째로 죽어 `POKEHOLE`은 남지 못하고 은퇴했다(아래 셀 (5b)).
- **`flip`이 그런 고리를 처음 만난다.** `Cut`의 B-조각이므로 링이 뒤집혀 면 법선이 `−z`가 되고, 그 면이 박스의 **바닥**이 된다. rim 엣지 넷은 `edge_for`의 정점쌍 dedup 덕에 스텁 옆면 조각이 쓰는 바로 그 엣지이고, `flip` 덕에 두 사용이 대립한다. (`flip` + **inner**는 여전히 미검증.)
- **검출기 표를 실측했다 — 구멍과 정반대다.**

| 오류 | props / OCCT | `validate` | `tessellate` | 부호 있는 메시 부피 |
|---|---|---|---|---|
| **구멍** 고리 뒤집힘 | 눈멂 | `NonOpposedEdge` | `HoleWinding` | (도달 못 함) |
| **섬** 고리 뒤집힘 | 눈멂 (0.08 그대로) | `NonOpposedEdge` | `Ok` — 침묵 | `2·0.16/3` 어긋남 |

  `tessellate`가 섬을 못 보는 이유는 링에서 법선을 **뽑기** 때문이다 — 삼각형이 조용히 안쪽을 향할 뿐이다. (gate) 셀이 없었다면 섬의 방향을 보는 기계가 `validate` 하나였다. `a_flipped_island_loop_is_caught`가 네 검출기에 물어 나온 대로 적었다.

- **`POKEHOLE`이 지키는 것.** 고리 둘은 **중첩**일 수 있고, 그러면 안쪽 고리의 내부가 다시 유지되어 **별개의 면**이 된다 — `Option<LocalFace>` 하나로는 못 돌려준다. 호 옆의 고리도 배치할 수 없다. **"오늘 픽스처에 없으니 도달 불가"라고 쓸 수 없다:** 다면체 **토러스**를 구멍 평면으로 자르면 중첩 고리 둘이 나오고 그 면들은 전부 단순하므로 `INNER_LOOP_OPERAND`도 막지 못한다. 셀 3f-3.
- **발화 테스트는 좁히기 전에 확보했다.** `Cut(slab, u)` — `cut_multi_chord_slab_is_unsupported`의 피연산자를 뒤집으면, 면이 A→B 순으로 도므로(`overlap_fuse_cut`) 슬래브의 `y=1.5` 면(고리 둘)이 U의 base cap(`multichord`)보다 먼저 운다. **실측했고**, 틀린 태그를 먹여 이빨을 확인했다. 태그가 면 열거 순서에 의존한다는 것을 주석에 적었다 — 순서가 바뀌면 소리내어 깨지는 것이 옳다.

**한 면이 여러 chord를 받는다 — `MULTICHORD` 은퇴 (M5-d3 셀 3e-2).** `reconstruct_face_paths`가 `Vec<LocalFace>`를 돌려주고, 한 입력 면이 여러 출력 면으로 갈라질 수 있다. `l_and_notch_bar`가 열린다: `Cut` `2.96`, `Fuse` `3.65`, 둘 다 `validate` 클린. OCCT 25.

- **호들은 순열이다.** 각 호의 두 끝은 하나가 `kd`(kept→dropped 전이), 하나가 `dk`다 — 호가 `f`를 둘로 가르고 두 `∂f` 조각의 분류가 다르므로. `dk[a]`에서 유지 정점을 따라 앞으로 걸으면 어떤 `kd[b]`에 닿고, 그 `b`가 `a`의 후계자다. 그 후계 사상의 **사이클이 곧 유지 영역**이고 사이클 하나가 `LocalFace` 하나다. 호가 하나면 사이클은 `[0]`이고 오늘의 splice와 **글자 그대로 같다**(골든이 못박는다).
- **`arrange::stitch_cycles`는 정수만 주고받는다** — `kept`, `kd`, `dk`. 그래서 골든이 종이 위 계산이 되고, 좌표 무독(無讀) 성질도 유지된다.
- **오늘 가정하던 것을 검사로 바꿨다.** 단일 호에서는 `if kept[t0]`로 kd/dk를 **가정**해도 무해했다. 여럿이면 조용히 틀릴 수 있고, 그 자리는 `classof`의 광선 캐스팅과 arrangement의 정확한 교차가 어긋나는 지점이다 — `SEAM_COUNT_MISMATCH`의 정의 그대로. 새 태그를 만들지 않고 세 검사를 얹었다: `kd`들이 호와 전단사, 중복 `kd` 없음, 사이클이 모든 호를 정확히 한 번. 순회에는 상한을 걸어 릴리스 무한 루프 대신 정직한 거절을 낸다.
- **`MULTICHORD`는 볼록 경로에만 남고 거기서 도달 불가다.** `∂f`는 볼록 폐곡선이고 상대는 볼록 집합이므로 `∂f ∩ B`는 호 하나 — 교차는 0 또는 2. `fourplane`·`tunnel`·`seam_branch`·`seam_count_mismatch`·`loop_orient_mismatch`와 같은 줄의 **미발화 백스톱**이 됐다. `unreachable!()`로 바꾸지 않는다.
- **볼록 경로 `reconstruct_face`는 `Option`으로 남는다.** 위 논증이 곧 이유다. `Vec`으로 넓히면 코드가 없는 자유도를 가진 척한다.

**★ §9가 없다고 단정한 픽스처가 있었다 — 모서리를 물면 된다.** 앞서 다중 chord만 요구하는 형상을 "좁히면 하나, 넓히면 다른 하나"라며 사실상 포기했다. 도망의 원인은 하나다: **한 호의 두 끝이 같은 엣지에 앉으면 그 엣지가 두 번 뚫려** `pierced_multi`/`poke_through`가 먼저 삼킨다. 그러니 절단자가 면의 **모서리**를 물게 하면 각 호의 끝이 서로 다른 두 엣지에 앉는다.

`l_and_notch_bar` = L-프리즘 + 노치에 누운 L자 막대(`z∈[0.5,1.5]`), 두 팔 끝이 캡의 볼록 모서리 `(2,1)`·`(1,2)`를 문다. 캡은 `kept=[T,T,F,T,F,T]`, 전이 넷이 **서로 다른 네 엣지**에 앉는다. `expect = (2, 0)` — 픽스처 전체에 닫힌 고리가 하나도 없다. **한 번에 한 변수.**

- **`Cut`과 `Fuse`가 서로 다른 재구성을 밟는다.** `Cut`의 B-조각에서 막대의 바닥 캡은 `kept=[T,F,F,F,T,F]`로 **두 사이클**을 내어 두 면(두 물린 조각의 바닥)이 되고, `Fuse`에서는 한 고리가 두 호를 쓴다. 그래서 OCCT diff 둘이 다 필요하다.
- **넓이가 이 형상에 눈멀다.** 모서리 물기는 없앤 세 면을 그대로 돌려주므로 `Cut` 결과의 넓이는 L의 `14.0` 그대로다. 부피와 `validate`가 채점한다. 검사는 교환 가능하지 않다.

**★ `n_out`이 한 모서리의 회전이었고, reflex 모서리가 그것을 뒤집었다 (셀 3e-2 픽스처를 만들다 발견).** `outer_tri`는 면의 외곽 링에서 **처음 세 연속 정점**을 잡고 그 회전을 링의 감김으로 믿었다. 그 모서리가 **볼록할 때만** 참이다. 스위트의 모든 픽스처가 우연히 안전했다 — 어느 것도 캡 링을 reflex 모서리 바로 앞에서 시작하지 않았다. L-프리즘의 프로파일을 `(0,0)` 대신 `(2,1)`에서 시작하면(**같은 솔리드**) 캡의 `n_out`이 `−z`로 나온다.

- **어디까지 샜나.** `order_along`은 무사하다 — `tri`를 exact orient3d에 넣고 `dir_sign`(→ `n_out`)을 곱하는데, 둘이 같은 삼각형에서 나오므로 두 번 뒤집혀 상쇄된다. `is_convex`·`classify_vertex`·`segment_enters`는 `n_out`을 볼록 경로에서만 읽고 거기엔 reflex가 없다. **위험했던 것은 `orient_seam_loop`** — `orient_sign` 둘을 곱하는데 상쇄해 줄 짝이 없다. 안쪽을 향한 `n_out`이 구멍이나 섬을 **조용히 뒤집는다**. 디버그에선 `debug_assert`가, 릴리스에선 `validate`와 부호 있는 메시 부피만이 안다.
- **고침:** 링 전체의 Newell 합이 부호를 정하고 삼각형을 그에 맞춰 뒤바꾼다. 한 모서리에 속을 것이 없다. `PlaneInfo::tri`는 이제 "링을 따른다"고 주장하지 않고 소비자가 실제로 요구하는 것 — **RH 법선이 바깥을 향한다** — 만 주장한다.
- **두 겹으로 못박았다.** `outward_normals_agree_with_their_orientation`은 서로 독립인 두 출처(링의 감김 vs b-rep의 `Surface`+`Orientation`)를 네 솔리드의 모든 면에서 맞대본다 — `debug_assert`가 면 하나씩 하던 일이다. `a_rotated_profile_is_the_same_solid_to_the_boolean`은 회전한 L에서 dimple을 잘라 2.92와 클린 `validate`를 요구한다. **둘 다 옛 `outer_tri`에서 실패한다.**
- **교훈:** "우연히 안전"은 `Vec` 하나만큼도 안전하지 않다. 프로파일을 회전해도 솔리드는 같아야 한다는 것이 이 버그를 잡은 최소 진술이었다.

**★ 예약을 집행했다 — 그런데 3-implicit `orient2d`는 끝내 필요 없었다 (셀 3h).** 앞서 "다중 implicit 술어가 정말 필요해지는 곳은 셀 추출(두 implicit 점의 orient)"이라 적었다. **절반만 맞았다.** 필요한 것은 두 implicit 점의 **좌표 비교** 하나였고, orient는 하나도 필요 없었다.

**이유가 우연이 아니라 구조적이다.** seam 고리의 모든 변은 `P ∩ Q_j` 선 위에 있고 방향이 `d_j = n_P × n_{Q_j}`다. 그러면 노드에서의 회전이 점을 만들지 않고 떨어진다:

```
(n_P × n_A) × (n_P × n_B) = n_P · det[n_P, n_A, n_B]        (a×b)×(a×c) = a·det(a,b,c)
  ⇒  turn = s_a · s_b · plane_pair_dir_sign(P, A, B) · orient_sign(P)
```

`s_j`는 `order_along`이 정확히 준다(순회의 기억이 아니라 **재계산**). `det`는 절대 `0`이 아니다 — 노드가 세 평면 위에 있고, 그 점이 존재한다는 것이 곧 법선들의 독립이다. `arrange::turn_at`.

**남는 것은 볼록 껍질 꼭짓점 하나를 고르는 일이고, 거기서만 두 implicit 점이 한 술어에 들어간다.** 사전식 최소 노드는 평면 위 점집합의 극점이므로 껍질 꼭짓점이다. `nacre_predicates::indirect_cmp_coord` — Cramer의 `N/D`를 교차 곱해 `sign(N_a·D_b − N_b·D_a)·sign(D_a)·sign(D_b)`. **Cramer는 공개하지 않는다**(numerator를 넘겨주는 것은 implicit 점을 재료화해 넘기는 것과 같다). `arrange::loop_winding`.

**★ 지름길을 찾았다가 반증했다 — 다음 사람이 같은 길을 다시 걷지 않도록 남긴다.** 고리의 각 변은 `P ∩ Q_j` 위에 있으므로, "다른 모든 노드가 `Q_j`의 한쪽에 있는가"는 **1-implicit** `three_plane_orient3d` 하나로 답한다. 그런 **지지 변**이 있으면 그 끝점이 껍질 꼭짓점이고 새 술어가 필요 없다. **그러나 단순 다각형은 껍질 위에 놓인 변을 하나도 갖지 않을 수 있다** — 정오각형의 각 변을 안쪽으로 얕게 접으면 껍질 꼭짓점은 전부 다각형 꼭짓점이지만 어떤 변도 껍질 위에 없다. 껍질 **꼭짓점**은 언제나 있고 **변**은 없을 수 있다. 지름길은 죽는다.

**감김은 `kept[0]`의 검산이다 — "두 기계"라 쓰면 부정확하다. 세 출처다.** `orient_seam_loop`이 낸 링의 방향은 `keep`에 의존하고 `kept[0]`도 그렇다. `keep`은 양쪽에서 **상쇄된다.** 검사가 실제로 묶는 것은 (1) 국소 부호 규칙(`n_out` 장부 + `order_along`), (2) 고리의 전역 감김(껍질 꼭짓점의 `det3`), (3) `classof`의 광선 캐스팅. 셋 중 어느 것도 옳다고 가정하지 않는다.

- 구멍 ⇔ `kept[0] == true` ⇔ 감김 `−1`(CW, 재료가 고리 밖). 섬 ⇔ `false` ⇔ `+1`.
- 불일치는 `LOOP_CLASS_MISMATCH`로 **릴리스에서 거절**한다. 오늘 도달 불가 — **동작 변화 0이 곧 셋의 합의다.**
- 3f-3에서 같은 등식이 **중첩 탐지기**로 승격된다: 감김이 각 고리의 분류를 정하고, `kept[0]`은 깊이 0에서만 검산한다.

**볼록 고리에서는 껍질 탐색이 무의미하다 — 그래서 비볼록 고리를 만들었다.** 오늘까지 닫힌 고리는 전부 사각형이었고, 볼록 다각형은 모든 꼭짓점의 회전이 같다. `l_and_ell_stub`(L의 캡 안에 선 L자 스텁)의 구멍은 reflex 노드를 갖는다. `Cut` = `3 − 0.1725·0.5`, 넓이 `15.3`, `validate` 클린, OCCT 26. `triangulate_polygon`의 브리징도 **비볼록 구멍을 처음 만나 그대로 통과했다**(리스크로 적어 두었으나 버그는 없었다).

- **★ 반증에 이빨을 붙이는 데 두 번 걸렸다.** `turn_at(reflex) == +1`인데 `loop_winding == −1`임을 단언하는 것은 **참이지만 무력하다** — 이 고리의 `ring[0]`이 우연히 볼록이라, `turn_at(ring[0])`을 쓰는 순진한 구현도 통과한다(실측). 고리는 순환이므로 **감김은 시작점에 불변**이다. reflex 노드로 회전시킨 링을 먹이면 순진한 구현이 `+1`을 낸다. **`outer_tri` 버그와 같은 함정을, 이번엔 코드가 아니라 테스트가 먼저 밟았다.**
- **3c의 미검증 전제("triple 동일성 ≠ 기하 동일성")가 처음 검사된다.** "triple 동일성 ≠ 기하 동일성 … 좌표 근접 검사로 `debug_assert`를 만들지 않는다 — tolerance가 필요하고, tolerance는 `Origin::Discovered` 밖에 설 자리가 없다." 사전식 탐색에서 세 좌표가 **정확히** 같은 서로 다른 triple이 나오면 `LOOP_ORIENT_MISMATCH`로 거절한다. 근접이 아니라 상등이므로 tolerance가 없다 — exact 술어가 그때의 거부 이유를 지웠다. (최소값과 겹치는 경우만 보므로 부분 탐지다.)
- **`SeamEnd::point`는 여전히 읽지 않는다.** 비교는 평면 계수만 본다. "`reconstruct_face_paths`는 좌표를 하나도 읽지 않는다"가 유지된다.

**한 면이 고리를 여럿 가질 수 있다 — 감김 하나가 중첩을 가린다 (M5-d3 셀 3f-3).** `u_and_slab`이 세 방향으로 열린다: `Cut(u,slab)` = `4.0`(섬 둘, `flip`), `Fuse` = `12.0`(같은 면의 구멍 둘), `Cut(slab,u)` = `6.7`(A의 구멍 둘 + U 캡의 사이클 분할). 넓이도 손계산과 맞는다(`18` / `42` / `33.2`). OCCT 29.

- **★ 포함 판정도 대표점도 필요 없다.** 고리 `L`의 감김은 깊이 `d(L)`의 **패리티**가 정한다 — 재료는 언제나 왼쪽이고, 고리를 건널 때마다 분류가 뒤집히므로 `L` 안쪽 영역의 분류는 `kept[0] XOR (d+1 홀수)`다. 직접 포함은 깊이를 1 바꾸므로 감김이 **반대**가 된다. 따라서

  > **모든 고리의 감김이 같다 ⟺ 어떤 고리도 다른 고리를 포함하지 않는다.**

  그리고 그때 전부 깊이 0이므로 감김은 `kept[0]`이 정한다. `classify_loops`가 부호 하나로 **중첩 탐지와 hole/island 판정을 동시에** 한다.

- **★ 그 논증에는 사전조건이 있고, 그것이 3f-4를 정의했다 (셀 3f-4에서 이행).** "깊이 0 영역"이 정의되려면 `∂f`가 **균일하게 분류**되어야 한다 — 즉 면에 **호가 없어야** 한다. 호가 있으면 반례가 있다: **유지 영역 안의 구멍(CW)과 버림 영역 안의 섬(CCW)은 서로 포함하지 않으면서 감김이 반대다.** `classify_loops`를 거기서 부르면 **거짓 `NESTED_LOOPS`**를 낸다. 그래서 "사이클이 하나면 구멍을 거기 붙이면 된다"는 지름길도 틀린다. 호 옆의 고리는 **진짜 포함 판정**을 요구한다.
  사전조건은 공짜로 강제된다: 3d의 집합 상등 덕에 `opens.is_empty() ⟺ transitions.is_empty()`이므로, `!opens.is_empty()`에서 거절하는 것이 곧 사전조건이다.

- **`NESTED_LOOPS`는 도달 가능하고 픽스처가 없을 뿐이다.** 다면체 토러스를 구멍 평면으로 자르면 중첩 고리 둘이 나온다. `fourplane`·`tunnel`의 "미발화 백스톱" 줄에 **넣지 않는다** — 순수 함수 골든이 그것을 직접 발화시킨다. (**토러스는 필요 없었다:** 셀 3f-5가 pocket된 큐브를 뚜껑과 바닥 사이에서 잘라 실물로 발화시킨다. 토러스에 대한 위 문장은 참이다 — 그 면들은 정말 전부 단순하다 — 다만 **훨씬 흔한 길이 `INNER_LOOP_OPERAND` 뒤에 가려져 있었다.** 도달 가능성을 논증할 때 가장 이국적인 예를 들면, 문 앞의 가드가 가리고 있는 평범한 예를 못 본다.)

- **`POKEHOLE`은 은퇴하지 않는다.** 남는 모양은 **호 옆의 고리** 하나뿐이고, `an_arc_beside_a_loop_is_unsupported`가 그것을 고정한다.

**★ 각기둥 보조정리 — 왜 arc+loop 픽스처가 번번이 도망갔는가 (셀 3f-3).** 절단자 `Y`가 각기둥이고 면 `f`가 **압출축에 수직**이면 `Y ∩ P`는 프로파일, 즉 연결집합이다. `프로파일 ∩ f`의 두 성분을 잇는 경로는 `f`를 벗어나고, 경로가 `f` 안에 머문 부분은 시작 성분에 속하므로 그 성분의 폐포가 `∂f`에 닿는다. **따라서 모든 성분이 `∂f`에 닿는다** — 성분이 둘 이상이면 고리가 없다. ∎

그래서 `l_and_staple`은 **축에 평행한** 면을 노린다: XZ 평면의 П자를 `−y`로 압출하면 L의 캡(`z=1`)이 축과 평행해 단면이 둘로 끊긴다. 가까운 다리는 캡 안에 고리를, 먼 다리는 반사 모서리 `(1,1)`을 감싸며 **서로 다른 두 엣지**에 끝을 둔 호를 남긴다(그래서 `pierced_multi`가 조용하다).

- **문 앞의 가드가 손계산을 잡았다.** `has_coplanar_pair`는 **합친** 평면 목록을 훑으므로 **한 피연산자 안의** 두 공면 면도 거절한다. 스테이플의 두 다리 바닥이 처음엔 둘 다 `z=0.5`였다 — `u_prism`이 기둥 높이를 어긋나게 둔 이유가 바로 이것이다. `0.45`로 내렸다.
- **부채꼴이 좌표 하나를 강제했다.** `y∈[0.7,1.3]`이면 관통점 `(1.4, 0.7)`이 캡의 부채꼴 대각선 `y = x/2` 위에 정확히 놓여 `contact_degenerate`다. `0.65`가 그것을 피한다. 서브유닛 5가 부채꼴을 없애면 사라질 제약이다. (**셀 (5a)가 없앴다.** `0.65`는 이제 불필요하지만 되돌리지 않는다 — 한 번에 한 변수.)

**호 옆의 고리 — `POKEHOLE` 은퇴, 그리고 세 번째로 아낀 술어 (M5-d3 셀 3f-4).** `l_and_staple`이 세 방향으로 열린다: `Cut(L,st)` `2.689`(구멍), `Fuse` `3.4495`, `Cut(st,L)` `0.4495`(**섬**). 셋이 `V_∩ = 0.311`로 포함–배제를 정확히 닫는다. OCCT 32.

- **★ 3f-3의 반례가 실물이 됐다.** **같은 면의 같은 고리**가 `Cut(L,st)`에서는 유지 영역 안이라 **구멍**이고, 피연산자를 뒤집으면 유지 영역이 모서리 물린 조각 하나뿐이라 그 **밖**이 되어 **섬**이다. 감김은 둘을 구별하지 못한다(구멍은 CW, 섬은 CCW, 그러나 서로 포함하지 않는다). **위치가 구별한다.**

- **★ 또 새 술어가 필요 없었다 — 광선을 `P ∩ Q_a` 위에 쏜다.** 링의 모든 변은 `P ∩ R` 위에 있고, 고리 노드 `v`는 제 평면 `Q_a`에 대해 `P ∩ Q_a` 위에 있다. 두 선의 교점 `X = {P, Q_a, R}`은 **그 자체로 3-평면 점**이므로, "`X`가 변 안쪽인가"도 "`X`가 `v`보다 앞인가"도 3d가 만든 `order_along`이다. 점을 만들지 않고 좌표를 읽지 않는다.

  | 물어야 할 것 | 술어 |
  |---|---|
  | 두 선이 만나는가 | `plane_pair_dir_sign(P, Q_a, R) ≠ 0` |
  | `X`가 변 안쪽인가 | `order_along(P, R, Q_a, S_i) · order_along(P, R, Q_a, S_j) < 0` |
  | `X`가 `v`보다 앞인가 | `order_along(P, Q_a, R, Q_b)`의 부호 = 광선 방향 |

  3b("1-implicit이면 충분했다")·3h("2-implicit 하나면 충분했다")에 이어 **세 번째로 예상한 술어가 이미 있던 술어로 밝혀졌다.**

- **★ 좋은 광선을 먼저 고르면 특수 경우가 통째로 사라진다.** 링 노드가 광선 **선** 위에 있으면 패리티가 애매하다. `three_plane_orient3d(T, Q_a.tri) == 0`이 그것을 정확히 판정한다. 그런 노드가 없으면 **변의 선이 광선 선일 수 없으므로**(끝점이 선 위에 있을 테니) `det == 0`은 언제나 "평행하고 다른 선" = 교차 없음이다 — 공선 처리가 필요 없고, 모든 교차가 횡단적이라 패리티가 곧 포함이다.

- **후보의 절반이 죽는다 — 실측.** 스테이플의 두 다리는 같은 `y = 0.65`·`y = 1.3` 캡에서 압출됐으므로 고리와 호가 그 평면들을 **공유한다.** 고리의 `y` 선을 따라 쏜 광선은 호의 노드 `(0.8,0.65)`·`(1.4,0.65)`를 정확히 지난다. `x` 선 둘만 살아남는다. **노드의 두 평면을 모두 후보에 넣어야 하는 이유다.** 하나도 없으면 `NO_CLEAR_RAY`(미발화).

- **답은 광선에 무관해야 한다 — 공짜 두 번째 기계.** 링이 단순하므로 모든 깨끗한 광선이 같은 패리티를 내야 하고, 골든이 그것을 단언한다.

- **`SeamPath::Closed`의 미검사 주장이 검사된다.** "닫힌 seam 고리는 `∂f`에 절대 닿지 않는다"를 3f-1 이래 아무도 확인하지 않았다. 교차점 `X == v`가 변 **안쪽**이면 `v ∈ M`이다 → `POINT_ON_RING`(미발화). 이제 모든 면의 모든 고리에서 검사된다.

- **원본 정점도 3-평면 점이다** — `f`의 평면과 인접 두 엣지의 이웃 평면. 그래서 사이클의 링(원본 + seam 혼합)이 seam 고리와 **같은 타입**이고 `point_in_ring`은 차이를 모른다. 전역 triple 유일성은 필요 없다(교차 판정이 변마다 국소적이다). 다만 한 이웃 평면 위에 **공선 엣지 둘**을 가진 면에서는 `ring_edge`가 거절한다 — `SeamSegment::on_edge`의 doc이 경고하던 그 면이고, 오늘 픽스처엔 없다.

- **세 출처의 배역이 바뀐다.** (1) 국소 부호 규칙이 링의 **방향**을, (2) **포함 판정**이 hole/island과 **중첩**을, (3) 감김이 고리마다의 **검산**을 맡는다. `classify_loops`(3f-3의 감김 기반 중첩 추론)는 삭제됐다 — 호가 있으면 틀린다. 잃은 것은 없다: 호가 없으면 `check_loop_class`가 모든 감김을 같게 강제하고, 그것이 3f-3의 정리로 곧 "중첩 없음"이다.

- **`POKEHOLE`은 볼록 경로에만 남는다.** 거기서는 볼록 B가 볼록 면의 내부를 뚫으면 `poke_through`가 먼저 난다고 **본다** — `MULTICHORD`와 달리 그 논증은 볼록 집합의 성질이 아니라 **다른 가드가 먼저 발화한다는 것**에 기댄다. **믿지만 검증하지 않았다.**

- **`NESTED_LOOPS`는 당시 실물 피연산자가 없었다** — 다면체 토러스가 필요한 줄 알았다. **틀렸다: 셀 3f-5가 pocket된 큐브를 뚜껑과 바닥 사이에서 자르면 난다.** 그 피연산자를 막고 있던 것이 `INNER_LOOP_OPERAND`였다.

**구멍 있는 피연산자 — `INNER_LOOP_OPERAND` 은퇴 (M5-d3 셀 3f-5).** 3f-1 이래 불리언은 구멍 있는 면을 **만들 수는** 있어도 **받을 수는** 없었다. `Cut(Cut(a,b),c)`가 안 됐고, pocket된 솔리드도 안 됐다. 계획 단계에서 가설 셋을 세우고 셋 다 코드를 읽어 반증했다.

- **`segment_crosses_face`는 이미 구멍을 옳게 처리한다 — 셀 3a가 고쳤다.** `face_loops`가 outer + `inner`를 모두 돌려주고 구멍 링은 CW이므로 부호 있는 부채꼴 합이 저절로 상쇄된다. `point_in_solid`·`boundaries_intersect`·`pierced_face`도 같은 함수를 쓴다. 부채꼴의 남은 문제는 **퇴화**(대각선 스침)뿐이고 구멍과 무관하다 — 사다리 항목 (5)를 앞당길 이유가 없었다. (셀 (5a)에서 `segment_crosses_face`는 삭제됐다. `point_in_solid`의 **광선** 부채꼴만 남는다.)
- **`nonconvex_seamfree`에는 가드가 없었다.** 구멍 있는 피연산자의 **포함/분리** 불리언은 이미 옳게 통과하고 있었다(`contained_result`가 셸을 통째로 재사용한다). 테스트만 없었다 — `containment_boolean_already_keeps_a_pocket`이 못 박았다(0.919 / 0.92).
- **"볼록 경로는 구멍 있는 면을 못 본다"는 틀렸다.** `imprint`가 region face에 `Surface`를 **재사용**하므로 imprint된 큐브는 **볼록한데 구멍 있는 면을 갖는다**. 옳은 명제는 **구멍 있는 면 ⇒ 오목하거나 공면 쌍이 있다**(rim 엣지의 반대쪽 면은 안쪽 벽이거나 같은 평면의 region face다; 닫힌 셸에서 셋째는 없다). 그래서 `fuse_cut`의 가드는 죽은 코드가 아니라 `has_coplanar_pair`와 **중복**이었고, 지우니 더 정확한 `COPLANAR_PAIR`가 나온다. `fuse_cut`의 옛 주석이 이 논증을 적어 두었으나 **측정된 적이 없었다.** 이제 잰다(`a_pocketed_cube_is_not_convex`, `an_imprinted_cube_is_convex`).

**★ 가드를 문에서 쓰이는 자리로 내리자 새 기계가 필요 없어졌다.** 씨임이 구멍의 **rim**을 건드리지 않으면 `∂f`는 여전히 링 하나다 — `stitch_cycles`의 순열 모형도, `transitions`도, `bnd`도 outer 링 위에서 그대로 성립한다. 그러면 `f`의 구멍은 arrangement가 보는 **또 하나의 닫힌 고리**이고, 3f-4의 `place_loops`에 **이어 붙이기만** 하면 담는 region·seam 고리와의 얽힘(`NESTED_LOOPS`)·모순이 한 호출로 나온다. 세 번의 셀(3b·3h·3f-4)에 이어 **네 번째로 술어를 아꼈다.**

seam 고리와 다른 점은 정확히 둘이다.
1. **구멍은 섬이 될 수 없다.** 재료가 아니므로 `owner == None`이면 그냥 버려진다(그것이 앉아 있던 `∂f` 영역과 함께).
2. **감김 검산이 `owner`와 무관하다.** `f.inner`는 언제나 CW이므로 `check_loop_class(true, w)`다. `check_loop_class(owner.is_some(), w)`를 그대로 쓰면 **버려지는 구멍마다 감김 `+1`을 요구해 실패한다** — 계획 검토가 배선 전에 잡은 버그다.

그리고 **`Node::Orig`로 방출해야 한다.** 구멍 링의 노드는 `f`의 원본 정점이다. seam 고리처럼 `Node::Seam(triple)`로 내면 `assemble_fuse_cut`이 같은 자리에 새 `Discovered` 정점을 만든다 — "arrangement는 **조합의 출처이지 기하의 출처가 아니다**"라는 원칙이 여기서 정점 표현으로 나타난다.

**가드는 이미 코드에 있었고 태그만 틀렸다.** `edge_ix`는 경계 노드의 엣지를 `f.outer.half_edges`에서 찾는다. 씨임이 rim을 가로지르면 못 찾고 **`SEAM_COUNT_MISMATCH`로 거짓말했다.** 이제 inner 링을 뒤져 `SEAM_ACROSS_HOLE_RIM`을 낸다. **완전하다:** `bnd`가 끝점만이 아니라 **모든 경로의 모든 노드**의 `on_edge`를 `edge_ix`에 통과시키고, `Closed` 경로는 `∂f`에 닿지 않는다. 그래서 그 뒤에서는 `∂f`가 링 하나임이 보장되고 `stitch_cycles`가 outer 링만 보아도 옳다.

**배선 전에 재야 했던 것 하나: 구멍 링의 감김.** `Face.orientation`(Forward/Reversed)과 `PlaneInfo.n_out`의 관계 때문에 `loop_winding(rim) == −1`이 보장되지 않는다. 순수 함수 골든(`a_pocket_rim_is_a_clockwise_ring_of_three_plane_points`)이 배선 **전에** 그것을 확정했다. `+1`이었다면 모든 구멍 면이 조용히 뒤집혔을 것이다. **순수 함수 + 손 골든 먼저, 배선은 다음.**

**★ 문 앞의 가드는 그 뒤의 버그를 숨긴다.** `arrange::seam_segments_on`의 두 스윕이 `outer.half_edges`만 돌고 있었다. 상대 면 `g`의 rim 엣지가 `f`의 평면을 `f` 안에서 뚫는 교차를 놓친다. **오답은 아니다** — 그런 엣지는 straddle하므로 `g`가 어차피 거절된다. 피해는 **태그 오염**이다: 그 면에서 홀수 카운트가 나 `ARRANGEMENT_DEGENERATE`가 정직한 가드보다 먼저 운다. **실측:** 스윕 수정을 되돌리면 `overlap_across_a_hole_rim_*` 두 테스트가 **양쪽 피연산자 순서 모두** `arrangement_degenerate`로 떨어진다. 검출기는 서로 대체 가능하지 않다. 같은 자리에서 `edge_incidence`의 반환형을 `Vec<usize>` → `[usize; 2]`로 좁혀 `inc[1]` 맹목 인덱싱(옛 패닉의 근원)을 **타입으로** 없앴다.

**두 좌표가 움직였고 둘 다 측정 결과다.**
- 코너 박스를 `[0.85,1.15]³`로 잡으면 그 수직 엣지가 뚜껑을 `(0.85, 0.85)`에서 뚫는데, 그 점이 뚜껑과 rim이 허용하는 **모든 apex의 부채꼴 대각선 위에** 있다 → `contact_degenerate`. 비대칭 `[0.85,0.8,0.75] → [1.2,1.15,1.1]`이 피한다. **서브유닛 5가 갚을 빚이다.** (**셀 (5a)가 갚았다.** 대칭 박스가 돌아왔고 부피는 `0.9125`가 아니라 아래 "결과"가 처음부터 적어 둔 `0.916625`다.)
- rim 테스트의 박스는 rim 코너 위에 걸치되 큐브 옆벽에 닿지 않는다. 닿는 박스(`[0.55,1.15]² × [0.85,1.15]`)는 pocket 공동을 관통해 벽으로 빠져나가므로 엣지 하나가 면 둘을 뚫고 **`pierced_multi`가 먼저 운다** — 정직하지만 다른 것에 대해서. (**3e-3이 그 가드를 은퇴시킨 뒤 다시 재니 그 아래에 `contact_degenerate`가 있었다** — 부채꼴이지 커버리지가 아니다. 가드 하나를 걷으면 그 아래 가드가 드러난다.)

**결과.** `cut_a_pocket_at_a_corner`(구멍이 `place_loops`로 배치)와 `cut_a_pocket_at_a_bottom_corner`(구멍이 `whole()`로 실려 나감)가 **서로 다른 두 경로로 같은 `0.9125`**를 낸다(셀 (5a)가 대칭 박스를 돌려준 뒤로는 `0.916625`). `coincident_merge_keeps_an_imprinted_hole`은 2.0(옛 오답 2.0533). `cut_slab_by_pocket`은 부피 1.99에 구멍 있는 면 **둘** — `z = 1`의 뒤집힌 뚜껑(`flip` + `inner`, 3f-1 이래 미검증이던 상호작용)과 `z = 0.3`의 발견된 구멍. **실려 들어온 구멍과 발견된 구멍이 한 솔리드에서 만난다.** OCCT diff 여섯(부피 + **넓이**), 32 → 38.

**★ `NESTED_LOOPS`가 처음으로 실물에서 발화한다.** 3f-3 이래 "도달 가능하나 픽스처가 없다"고 적어 왔고 다면체 토러스를 상상했다. 필요 없다 — **pocket된 큐브를 뚜껑과 바닥 사이 높이에서 슬래브로 자르면** 슬래브 밑면 위에 큐브 단면이 구멍으로, 그 안에 pocket 단면이 섬으로 앉는다. 두 kind × 두 피연산자 순서 모두 발화한다. 그것을 막고 있던 것이 `INNER_LOOP_OPERAND`였다. (토러스는 애초에 그 가드에 걸리지 않았다 — 그 면들은 전부 단순하다. 도달 가능성을 논증하며 **가장 이국적인 예**를 들었기에, 문 앞의 가드가 가리고 있던 평범한 예를 못 봤다.)

**남긴 정직한 빚.**
- **`NESTED_LOOPS`가 이제 필요보다 넓게 거절한다.** `f`의 구멍이 seam **구멍** 고리 안이면 옳은 답은 "버린다", seam **섬** 고리 안이면 "그 섬의 `inner`"인데, `place_loops`의 쌍별 검사가 둘 다 거절한다. 잘라 내는 쪽의 단면이 구멍을 **에워쌀 때** 걸린다. 3f-7(3f-6은 rim을 **가로지르는** 별개 경우였고 먼저 열렸다). **(셀 3f-7이 청산 — 깊이 패리티 forest. `Cut` 양쪽이 열렸고, `Fuse`는 pocket을 밀봉해 cavity가 되어 `disconnected_result`로 옮겨갔다((5c) 몫).)**
- **`NON_MANIFOLD_EDGE`·`HOLE_CLASS_SPLIT`는 미발화 백스톱이다.** 전자: `validate`의 `NonOpposedEdge`가 보장하지만 **`boolean`은 `validate`를 부르지 않는다.** 후자: 씨임 없는 면의 rim 정점 분류가 outer와 다르면 상대 경계가 `f`를 가르므로 씨임이 있어야 한다.
- **`solid_local_faces`의 remap 경로는 미검증이다** — merge 픽스처에서 구멍은 A쪽에 있고 remap은 B쪽에만 걸린다.
- **`whole()`의 두 번째 호출부는 도달 불가로 보인다** — `paths.is_empty()`가 집합 등식으로 `touches_seam == false ∧ transitions == ∅`를 강제하고, 그 쌍은 이미 앞에서 반환했다. 방어 가지로 남긴다.
- **`overlap_fuse_cut`·`fuse_cut`로 오는 구멍은 pocket과 불리언 결과뿐이다.** imprint된 솔리드는 `has_coplanar_pair`가 막는다. 두 가지가 그것에 걸려 있다: (a) rim 엣지의 두 인접 평면이 같아져 `order_along`이 0을 내는 퇴화가 도달하지 않는다. (b) `surf_ix`는 `Surface → usize`인데 imprint는 두 면이 **같은 `Surface` 핸들을 공유**하므로 인덱스가 뭉개진다. `coincident_merge`만 imprint를 받고, 거기서는 `solid_local_faces`가 `plane_offset + pos`로 면마다 다른 인덱스를 준다.

**상대를 관통하는 엣지 — `PIERCED_MULTI` 은퇴, 구멍 뚫기 (M5-d3 셀 3e-3).** `Cut(L, rod)`가 **genus 1** 솔리드를 낸다. 부피 `2.94`, 넓이 `14.88`, `validate` 클린, OCCT가 둘 다 확인. 이 커널이 만든 첫 토러스다.

**막고 있던 것은 표현이 아니었다.** `poke_through_hole_is_unsupported`의 주석은 "through-hole은 아직 표현 불가"라고 말했는데, 3f-1이 면에 inner loop를 얹은 이후로 **거짓이었고** 3f-5가 그것을 다시 입력으로 받은 뒤로는 두 겹으로 거짓이었다. 못 따라간 것은 (a) 양 끝점이 같은 쪽인 엣지에 씨임 정점을 하나도 만들지 않는 씨임 목록과 (b) `∂f`를 **엣지로** 인덱싱하는 경계 모형이다. **가드의 주석은 코드보다 오래 산다.**

**★ 인덱스 공간을 엣지에서 교차점으로 옮기면 `stitch_cycles`가 손대지 않고 산다.** `∂f`를 교차점이 `2m`개의 **run**으로 나눈다. 교차점은 언제나 부류를 뒤집으므로 **run의 부류가 교대하고**, 그래서 `succ`의 "kept run을 앞으로 걷기"가 정확히 한 걸음이 된다. `kept[t] && !kept[t+1]`은 "교차점 `c_{t+1}`이 kd"와 동치이므로 호의 슬롯은 **kd 교차점 − 1**이다. 정수만 주고받는 함수는 인덱스가 무엇을 뜻하는지 모르므로 의미를 갈아 끼울 수 있었다 — 그 함수의 **doc은 갈아 끼워야 했다**(정점 → run, 걸음이 하나).

run은 정점을 **0개** 가질 수 있다. 그것이 "한 엣지가 두 번 교차됨"의 정체이고, 로드의 벽면은 그런 run을 둘 갖는다(그 면의 `loop_nodes`는 seam 노드뿐이다).

**한 엣지 위 두 교차의 전후는 `order_along(P, R, ·, ·)`, 순회 방향은 `edge_sign` — 둘 다 이미 있었다.** 3c 기록이 앞의 것을 예고했고, 뒤의 것은 계획이 "엣지의 두 끝점도 그 선 위의 3-평면 점이니…"라고 새로 유도했는데 **유도된 식이 3h가 쓴 `edge_sign`의 정의 그 자체였다**(`R_{i−1} == R_{i+1}` 퇴화 처리까지). 3b·3h·3f-4·3f-5에 이어 **다섯 번째로 아낀 것**이고, 이번에 아낀 것은 술어가 아니라 **함수**다.

**`SEAM_COUNT_MISMATCH`가 3d에서 예고한 자리를 받았다.** 집합 상등 `bnd == transitions`는 죽고 그 자리에 **교대(기계 A) vs `classof`(기계 B)** 가 선다. 정점을 가진 run은 그 정점의 분류를 갖고, 씨앗 하나에서 교대로 전파한 뒤 **모든 정점을 검사**한다. 정점 없는 run은 교대만이 정한다 — 다른 진리 출처가 없다. 덤으로 `n == 2·|호|`가 **호 내부 노드가 `∂f`에 닿는 경우**(오늘 `debug_assert`만 있던 자리)를 릴리스에서 잡고, 호의 두 끝은 부류가 반대여야 한다.

**★ `TUNNEL`이 제 이름을 얻었다 — 패리티 가드.** 선분이 닫힌 곡면을 홀수 번 지나는 것과 양 끝점이 반대쪽인 것은 동치다. `point_in_solid`의 감김과 `pierced_faces`의 정확한 교차가 그 둘을 말하고, 어긋나면 거절한다. 앞서 "`tunnel`의 유일한 진짜 발화 경로는 cavity 비대칭"이라 적어 둔 것이 이제 **remark가 아니라 가드의 형태**다. 셀의 심장은 그보다 작다: `s0 == s1` 가지가 **`continue`를 멈춘다.**

**★ 문 앞의 가드는 뒤에 있는 버그뿐 아니라 뒤에 있는 요구사항도 숨긴다.** `Cut(rod, L)`은 로드를 **둘로 자른다.** 가드를 걷고 재니 `Ok`, 부피 `0.06`, `validate`가 `NegativeGenus { v:16, e:24, f:12, genus:-1 }` — 한 셸 안의 떨어진 두 상자. **`boolean`은 `validate`를 부르지 않으므로 아무도 말하지 않았을 것이다.** `Solid`는 outer shell 하나를 갖고 `boolean`은 핸들 하나를 돌려주므로, 두 성분은 두 솔리드다. 새 가드 `DISCONNECTED_RESULT`(면의 `Node` 공유로 union-find)가 **발화 테스트와 함께** 태어난다. 오늘까지 도달 불가였던 이유가 정확히 `PIERCED_MULTI`다 — A를 가르려면 A의 엣지가 B를 관통해야 하므로.

**잃은 것과 드러난 것, 정직하게.**
- **cavity를 통째로 가로지르는 엣지의 그물이 사라졌다.** 그런 엣지는 outer shell을 **두 번** 뚫으므로 패리티가 맞고(짝수 == 짝수) 통과한다. 옛 `pierced_face`는 hits ≥ 2를 그 자리에서 삼켰다. 당시엔 문 앞의 `HOLLOW_OPERAND`만이 막았다. **((5c-in)이 문을 열어 seam front-end를 all-shells로 바꿨고, (5c-in-2)가 재구성 세그먼트 경로까지 all-shells로 해 void에 *들어가는* cut도 재구성한다 — 열린 void는 outer와 병합. 패리티 TUNNEL은 비대칭이 없어져 은퇴, 잠시 섰던 `SEAM_ENTERS_CAVITY`도 클래스 전체가 정답-이거나-loud라 은퇴.)**
- **`NO_ENTRY_FACE`가 비볼록 경로에서 미발화 백스톱이 됐다** — 기계가 옳으면 `s0 != s1`에 hits가 빌 수 없다. 패리티보다 먼저 검사해 더 구체적인 태그를 남긴다.
- **3f-5의 rim 픽스처 아래에는 `contact_degenerate`가 있었다.** `pierced_multi`를 피하려 고른 좌표였는데, 그 가드를 걷고 재니 부채꼴이 나온다. 서브유닛 5의 빚 다섯째 줄. (**셀 (5a) 뒤에 다시 재니 자연스러운 `[0.55,1.15]² × [0.85,1.15]`가 양쪽 피연산자 순서 모두 `seam_across_hole_rim`까지 도달한다.** 가드 셋을 걷어 낸 자리에 마침내 그 도형에 관한 가드가 섰다.)
- **볼록 경로는 여전히 `POKE_THROUGH`로 거절한다.** `reconstruct_face`가 `edge_seam`(엣지당 triple 하나)을 읽고 `Option<LocalFace>`를 낸다. `cube_and_notch`가 그 자리에 남는다 — **볼록성의 사실이 아니라 그 자료구조의 한계다.** 사다리 (5).
- **`multichord`·`pierced_multi`를 피하려 비틀었던 픽스처들**(`l_and_notch_bar`의 모서리 물기, `l_and_staple`의 축평행 면, `l_and_popup_box`의 여유 높이)은 **그대로 둔다.** 한 번에 한 변수.

**M5-d3의 가드는 모두 은퇴했다.** `SEAM_ACROSS_HOLE_RIM`이 마지막이었고 3f-6이 은퇴시켰다(3f-5가 `INNER_LOOP_OPERAND`를 그것으로 좁혔고, 3e-3이 `PIERCED_MULTI`를 은퇴시켰다). `NESTED_LOOPS`(닫힌 고리 nesting)도 3f-7이 은퇴시켰다. **출력 cavity는 (5c)가 열었다**(씨임 경로가 밀봉 void를 cavity 셸로; `disconnected_result`는 진짜 분리에만 남음). **입력 hollow는 (5c-in)·(5c-in-2)가 열었다**(비껴가는 cut·void에 들어가는 cut 모두 재구성; 셸-병합으로 void가 열린 pocket이 됨). 남은 확장은 (5d) sweep #4·#5.

**선분 부채꼴이 죽는다 — `CONTACT_DEGENERATE` 은퇴, 좌표가 풀린다 (M5 셀 (5a)).** 커버리지 경계가 기하가 아니라 **좌표의 대칭성**에 달려 있었다. 이 커널의 DNA에 어긋난다.

**★ 앞서 "fan 대신 실제 삼각분할(또는 대각선을 피하는 apex 선택)"이라던 예측이 틀렸다.** 둘 다 답이 아니다 — 어떤 삼각분할이든 대각선을 만들고, 정사각형의 정중앙은 **모든** apex의 대각선 위에 있다(이미 실측돼 있던 사실인데 그 예측이 읽지 못했다). 옳은 관찰은 하나뿐이다: **엣지가 면을 뚫는 점은 3-평면 점이다.** 엣지는 평면 `E0`·`E1` 위에, 면은 `Q` 위에 있으니 교점은 `{E0, E1, Q}`이고, **셀 3f-4가 다른 이유로 지은 `point_in_ring`이 바로 그 점의 포함 판정**이다. 3b·3h·3f-4·3f-5·3e-3에 이어 **여섯 번째로 아낀 술어.**

`arrange::edge_crosses_face`는 exact `plane_side`로 straddle을 먼저 묻고(대부분의 (엣지, 면) 쌍이 여기서 걸러진다 — 옛 부채꼴은 **모든** 쌍에서 면 전체를 잘랐다), straddle이면 `det[n_E0, n_E1, n_Q] ≠ 0`이 보장되므로 `X = {E0,E1,Q}`가 존재한다. `point_on_ring`이 `X ∈ ∂g`를, `point_in_ring`이 `X ∈ g`를 답한다. **좌표를 읽는 자리는 `plane_side`의 부호 둘뿐이고 그것도 exact `orient3d`다.**

- **`point_in_ring`만으로는 태그가 거짓말한다 — 검토가 배선 전에 찾았다.** `X`가 면의 **꼭짓점**에 정확히 놓이면 `X`의 두 평면이 둘 다 그 링 노드를 지나 `every_ray`의 광선 후보가 전부 불명확해지고 **`NO_CLEAR_RAY`**가 뜬다. 정직한 답은 `POINT_ON_RING`이다. 그래서 `point_on_ring` 선검사를 세웠다 — 새 술어 없이 `side_of` + `order_along`으로. `an_edge_touching_a_face_is_a_contact_not_a_graze`가 두 자리(변 내부, 꼭짓점)를 모두 못박고 후자에 그 주석을 단다.
- **`plane_side`는 geom의 한 줄이다.** `three_plane_orient3d`는 implicit 점을 받고, 이쪽은 **explicit** 점이 평면의 어느 쪽인지 묻는다. 부호 규약이 뒤집힐 수 있어 **배선 전에 골든과 proptest로 못박았다**(3f-5의 `loop_winding(rim) == −1`이 그랬듯이). geom 67.
- **세 호출부, 하나의 퇴화 정책.** `seam_segments_on`·`pierced_faces`·`boundaries_intersect`가 모두 `edge_crosses_face` 위로 옮겼고, `boundaries_intersect`는 계약이 같으므로 `pierced_faces` 위에 **다시 썼다.**
- **광선 부채꼴은 산다.** `point_in_solid`의 광선 원점은 임의의 explicit 점이라 교점이 3-평면 점이 아니다. `RAY_DEGENERATE`는 남고, 이 셀이 죽인 것은 **선분** 부채꼴뿐이다.

**★ 부채꼴의 죄는 스침을 접촉과 구별하지 못한 것이다.** 이제 `VERTEX_ON_FACE_PLANE`(끝점이 `Q` 위 — 둘 다면 엣지가 `Q` 안에 통째로 누웠다)과 `POINT_ON_RING`(관통점이 `∂g` 위)이 **진짜 접촉**을 이름으로 부르고, 스침은 존재하지 않는다. `common_rejects_non_convex_input`의 박스는 `[0,0,0]`–`[1.5,1.5,0.5]`, 즉 **L의 바닥면과 공면**이다. L의 바닥 엣지들이 통째로 그 평면 안에 눕는다. "모든 apex에서 대각선을 스쳤다"고 부르던 것이 처음부터 접선 접촉이었고, 이제 태그가 그렇게 말한다. 3c의 미검증 전제 (b)가 여기서 검증됐다.

**좌표가 풀렸다.**
- `cut_a_pocket_at_a_corner`가 대칭 `[0.85,1.15]³`로 돌아가 **`0.916625`** — 앞서 "서브유닛 5가 갚을 빚"이라 적고 3f-5가 `0.9125`로 우회했던 그 수. 넓이 `6.8`(직각 코너의 사각 물기는 사분면 셋을 사분면 셋과 맞바꾼다). 두 코너 컷 모두 mesh 게이트 네 검사와 OCCT diff를 통과한다. OCCT 40 유지.
- `overlap_across_a_hole_rim_*`가 자연스러운 `[0.55,1.15]² × [0.85,1.15]`로 돌아가고 **여전히 `SEAM_ACROSS_HOLE_RIM`**이다. 그 도형은 `pierced_multi` → `contact_degenerate` → `seam_across_hole_rim`, **가드 셋 아래**에 있었다.
- 정중앙 관통 게이트가 열렸다 — `overlap_fuse_cut(two_boxes)`가 정중앙을 뚫고 `Cut` 부피 `0.875`, `validate` 클린.
- `l_and_staple`의 `y = 0.65`(3f-4), `cube_and_notch`·`seam_segments_on_two_overlapping_boxes`의 비대칭 extent는 **되돌리지 않았다.** 한 번에 한 변수. 되돌릴 수 있다는 사실만 여기 적는다.

**대가 하나 — imprint된 비볼록 피연산자.** 정확한 포함 판정은 면의 링을 **triple로** 요구하는데, imprint된 면은 rim 엣지의 두 이웃이 **공면**(구멍 있는 면 + 거기서 잘라낸 region 면)이라 triple을 낼 수 없다. `overlap_fuse_cut`·`fuse_cut`·`common`은 이미 `has_coplanar_pair`를 문 앞에 두고 있었고, **`nonconvex_seamfree`만 없었다** — `contained_result`가 셸을 통째로 재사용하며 링을 보지 않았기 때문이다. 이제 그 가드도 선다(피연산자 **각각** 안의 공면 쌍만 본다; 교차-솔리드 공면은 무해하다 — 면의 링은 제 솔리드의 평면만 쓰고, `X`의 `E0`·`E1`은 straddle이 `det ≠ 0`을 보장하므로 `Q`와 평행할 수 없다). **정직한 회귀다:** imprint된 비볼록 피연산자의 포함/분리 불리언이 `COPLANAR_PAIR`로 거절된다. 발화 테스트 `an_imprinted_nonconvex_operand_rejects_as_coplanar_pair`(pocket 바닥에 imprint한 pocketed cube)가 그것을 고정하고, `containment_boolean_already_keeps_a_pocket`이 **pocket은 그대로 통과함**을 고정한다 — 가드가 imprint만 무는지, 모든 구멍을 무는지의 차이다.

**남긴 정직한 빚.**
- **좌표를 읽는 자리가 하나 남았다.** `plane_side`에 넘기는 `p0`·`p1`은 저장된 `Point3`다. 원본 정점이면 그것이 진리이지만, 불리언 **결과를 다시 넣으면**(3f-5 이래 가능) `Discovered` 정점의 점은 반올림된 캐시이고 진리는 triple이다. `face_vertex_triples`로 정점의 triple을 되살려 `three_plane_orient3d`로 부호를 내면 좌표를 아예 안 읽는다. **이 셀은 그러지 않았다** — `point_in_solid`도 같은 좌표를 읽고, 옛 부채꼴도 같은 좌표로 같은 판정을 했으므로 커버리지도 정확도도 변하지 않는다. 다음 기회.
- **ops 스위트가 `1.16s` → `3.85s`(약 3.3×).** `seam_segments_on`이 `(f, g)` 쌍마다 `face_rings`를 다시 만든다. 측정 후 최적화 — 링 triple을 상위에서 캐시하면 된다.
- **"먼 면의 접촉이 거절 태그를 갈아치운다"는 경고는 그대로다.** 씨임과 무관한 먼 면의 접촉이 거절 태그를 갈아치울 수 있다. 이름이 정확해졌을 뿐 문제는 남는다.
- **`face_vertex_triples`의 거절 표면이 전역이 됐다.** 이제 교차 검사에 관여하는 **모든** 면에서 부른다. `LOOP_ORIENT_MISMATCH`(180° 각)와 `ring_edge`의 "한 이웃 평면 위 공선 엣지 둘" 제약이 전역이다. `every_edge_against_every_face_on_every_fixture`가 오늘의 픽스처 전부 × 양쪽 방향에서 교차 수를 못박아 그것을 지킨다.

**다음.** (5b) 볼록 경로 통합(`fuse_cut`·`edge_seam`·`reconstruct_face`·`MULTICHORD`·`POKE_THROUGH` 은퇴 — 3e-3이 poke-through를 열어 뒀고, (5a)가 정중앙을 열었다), 이어서 3f-6.

**부호 층에 부동소수점 필터를 달았다 — `indirect_orient3d`·`det3_sign`, 동작 변화 0 (M5 셀 (5b-0)).** (5b)를 세우며 "두 경로를 나란히 재는 마지막 기회"라 잰 n0가 (5b) 계획의 비용 문장 셋을 반증했다: `overlap_fuse_cut(two_boxes)`가 `fuse_cut`의 **~100배**(5477 µs vs 56 µs, release)였고, 그 시간의 전부가 `three_plane_orient3d` 한 번의 **1.46 µs**였다. Shewchuk `orient3d`는 ~50 ns다.

- **원인은 ops가 아니라 `nacre-predicates`였다.** `indirect_orient3d`가 매 호출 `Expansion`을 만들어 exact 산술을 끝까지 갔다 — 크레이트 doc이 *"no adaptive fast path yet … a later optimization"*이라 예고한 바로 그 부재. Attene 2020의 implicit predicates가 실용적인 이유가 그 필터다.
- **`indirect_orient3d_filter`:** 같은 다항식을 f64로 풀고 반올림 오차 한계(`|fl(x)| > εₓ·x̃`, Higham §3.1의 cancellation-free 척도)로 부호가 확실할 때만 답한다. 유도된 오차는 `ε_D ≈ 5u`·`ε_M ≈ 19u`, 쓰는 상수는 `32u`·`64u`(6×·3.4× 여유). **키우는 것이 안전한 방향** — exact로 더 자주 내려갈 뿐 틀린 부호는 없다.
- **★ "동작 변화 0"은 구조적이다.** 필터는 `|fl(x)| > 0`일 때만 답하므로 **`0`을 절대 주장하지 않는다.** 따라서 모든 영 판정이 여전히 exact를 지나고, 이 커널의 퇴화 가드는 전부 영 판정이다(`FOURPLANE`·`NO_CLEAR_RAY`·`POINT_ON_RING`·`THREE_PLANES`). 거절 태그는 경험이 아니라 구조로 움직일 수 없다. 실측 게이트는 i128 독립 오라클과 **OCCT diff 40**.
- **★ 일곱 번째로 아낀 술어.** 남은 나머지는 `det3_sign`이었다 — `every_ray`에서 `three_plane_orient3d`보다 호출이 많다(`order_along`이 `three_plane_orient3d × dir_sign`, `dir_sign`이 `det3_sign`). 그런데 `det[a,b,c] = orient3d(a,b,c,0)`이고 `geometry_predicates::orient3d`는 **이미 적응형**이라, `det3_sign`을 그 위임 한 줄로 바꾸면 필터가 공짜로 딸려 온다. 그 항등식을 못박은 `prop_det3_sign_matches_orient3d`가 술어보다 먼저 있었다. 3b·3h·3f-4·3f-5·3e-3·(5a)에 이은 일곱 번째.
  - 단, 그 테스트는 `f64::signum(0.0) == 1`이라 **특이 행렬에서 공허**했고(랜덤 연속 입력이 그 자리를 안 밟아 통과), 위임 뒤엔 **동어반복**이 됐다. m2a가 `sign_f64`로 고치고 특이 행렬 proptest를 더했으며(Shewchuk가 정확한 영을 낸다는 이 트리 최초의 확인), m2b가 내용을 `det3(m).sign() == det3_sign(m)`(expansion vs 적응형, 두 기계)으로 되찾았다.
- **실측 (release, `two_boxes` Cut, m1+m2b 후):** `three_plane_orient3d` `1.46 → 0.011 µs`, `overlap_fuse_cut` `5477 → 262 µs`, ops 스위트(debug) `3.76 → 0.63 s` **6×**. **(5a)가 기록한 "ops 3.3× 느려짐"의 원인이 이것이었다** — exact 산술이 아니라 필터의 부재.
- **★ 볼록 vs ray-cast 통합 물음의 답을 이 숫자 위에 다시 쓴다.** 필터 뒤에도 씨임 경로는 볼록 경로의 **~21–27배**다. 그러나 그 잔차는 exact 산술이 아니다 — 술어는 이제 `0.01 µs`로 무시할 수준이고, `point_in_ring`의 `1.26 µs`는 그 안에서 도는 **~60번의 값싼 술어 호출**, 즉 광선 캐스팅의 `O(r)` 순회다. **그것을 더 줄이는 것은 (5d)의 `plane_orient`(광선 없는 `O(F)` 분류기)이지 필터가 아니다.** "measure then decide"의 실측 답: 필터로 끝, 나머지는 알고리즘이라 (5d).
- **필터 발화율도 쟀다** — 축정렬 박스 평면(불리언이 만드는 계수)의 3-평면 정점 × 면 질의 768건에서 **88.3%**. 축정렬은 공면 영이 가장 많은 최악이고, miss 12%는 어차피 exact로 가야 하는 진짜 영이다.
- **크레이트 doc을 고쳤다.** *"no adaptive fast path yet"*은 **`Expansion` 값에 대해서는 여전히 참**이고 `indirect_orient3d`의 **부호에 대해서는 거짓**이 됐다 — 필터가 값이 아니라 부호에 붙는다.

**볼록 Fuse/Cut 경로가 죽었다 — `fuse_cut`·`classify_vertex`·`enter_face`·`POKE_THROUGH` 은퇴 (M5 셀 (5b)).** 애초 "속도 프로파일링 후 결정"이라 미뤄 둔 볼록 vs ray-cast 통합 질문의 답이다(위 "★ (5b)가 이 질문에 답했다" 참조). **속도가 아니라 커버리지와 정확성이 답했다.** fast path의 전제 "같은 답을 더 빨리"가 깨진다.

- **커버리지.** 볼록 경로는 일반 경로가 답하는 입력을 거절했다: `cube_and_notch`(엣지 두 번 교차)와 드릴 큐브(막대 관통)를 `POKE_THROUGH`로. 넷 다 비볼록 쌍둥이가 이미 그린이었고(`cut_notch_bar`·`drill_through_the_l`·`fuse_the_l_and_the_rod`·`cut_rod_by_l_disconnects`), 막던 것은 볼록성이 아니라 `edge_seam`의 엣지당 triple 하나(3d에서 소멸)다. **드릴 큐브는 이 커널의 두 번째 토러스이고, 첫 번째와 달리 두 피연산자가 모두 볼록이다.** 이 셀은 새 기계를 만들지 않았다 — 디스패처 한 줄이 막던 것을 열었다.
- **정확성.** `classify_vertex`의 `1e-9` tolerance와 `enter_face`의 `t = d0/(d0−d1)` argmax가 `Fuse`/`Cut`의 진리를 f64로 정했다. 씨임 경로는 exact `point_in_solid`와 exact `edge_crosses_face`를 쓴다. **`Fuse`/`Cut` 경로의 마지막 부동소수점 술어가 사라졌다** — 커널 전체에서가 아니다(남은 자리는 §9 **"(5d) exactness sweep — 전수조사 표"**: `is_convex`·`coplanar`·cavity 분류 부호 등. `order_ccw`는 3g에서 완료).
- **★ 볼록 vs ray-cast 물음의 (a) 예고가 한 셀 늦었다.** "서브유닛 4에서 `classify_vertex` 은퇴 예정"이라 했으나 (5b)까지 살았다.
- **`common`은 남고 이유가 있다.** 열거(`three_planes`+`three_plane_orient3d`)는 exact라 논거 (b)가 안 걸린다. `is_convex`도 남는다 — `detect_coincident_interface`가 문지기로 쓴다((5b-0) 표대로 (5d) 몫). 비볼록 `Common`(3g)이 `common`을 지운다.
- **태그 넷이 발화 테스트 없이 죽었다** — `MULTICHORD`·`ON_BOUNDARY`·`POKEHOLE`·`OUTSIDE_OR_FOURPLANE`. 디스패치가 `is_convex`를 먼저 검사하므로 볼록 전용 백스톱이었다. `NONCONVEX_OPERAND`도 같은 이유로 미발화지만 `common`이 쓰므로 남는다.
- **★ 3f-2의 예측이 뒤집혔다.** "`POKEHOLE`은 개명하지 않는다 … 볼록 경로의 문자 그대로 poke-through만 남아 이름이 다시 정확해진다"고 적었으나, (5b)가 볼록 경로를 지워 `POKEHOLE`은 살아남지 못하고 죽었다.
- **★ 태그 하나가 움직였다 — 예측보다 하나 더.** `same_ground_overlap`**과 `coincident_merge_rejects_offset_footprint`**가 `COPLANAR_PAIR` → `VERTEX_ON_FACE_PLANE`으로 바뀌었다(n0 인구조사는 고정 픽스처라 하나만 예측했다 — 측정이 둘째를 찾았다). `fuse_cut`의 **결합** 평면 공면 검사가 사라지고 씨임 경로의 **피연산자별** 검사가 통과한 뒤 `plane_side`가 공유 평면을 만난다. **커버리지는 보존**(여전히 거절)되고 태그만 한 평면 덜 구체적이다 — "거절 태그 갈아치움" 경고의 두 번째 실례. 결합 검사를 앞당기면 컨테이너 벽에 닿는 포함 박스를 과잉 거절하므로 **열린 문제로 남긴다.**
- **`fuse_common_inclusion_exclusion`이 두 기계를 잇는 항등식이 됐다** — `V_A+V_B = V_fuse+V_common`에서 `fuse`는 씨임 경로, `common`은 열거 경로. 3g가 `common`을 지울 때까지 교차 검증이다.
- **검증.** 손 골든 `993.28`/`1014.4`/`24`/`29`, 드릴 큐브 genus-1은 `holed_faces == 2`(뚜껑 둘)로(validate 클린이 genus-1을 증명하지 않으므로), 노치 컷의 정점 없는 run을 가진 면은 mesh 게이트 넷으로, `PROPTEST_CASES=4096`이 씨임 경로로 통과(`Cut`의 `.unwrap()`도 `prop_assume` 거절률도 그물). OCCT diff 40 → 44(노치·드릴 Cut/Fuse). ops 스위트는 (5b-0) 덕에 그대로 빠르다.
- **부채 둘.** (1) `boundaries_intersect`는 `overlap_fuse_cut`의 `seam.is_empty()`와 같은 질문을 두 번 묻는다 — 삭제하면 스캔이 한 번이 된다(답·태그 불변). (5b-0)이 속도를 이미 닫아 성능 동기는 없고, 중복 제거가 근거이므로 미뤘다. (2) `edge_crosses_face`의 `VERTEX_ON_FACE_PLANE`은 접촉이 아니라 **평면 위**(공면)면 발화한다 — (5a)가 남긴 과잉 거절이 이제 모든 불리언에 노출됐다. 좁힘은 `point_on_ring`/`point_in_ring`로 가능하다.
- **다음.** 3g(비볼록 `Common`), 3f-6(`SEAM_ACROSS_HOLE_RIM`), (5c) cavity 결정, (5d) exactness sweep.

**비볼록 `Common`을 열고 볼록 `common`을 지웠다 — `order_ccw`의 `atan2` 은퇴 (M5 셀 3g).** (5b)의 테마가 이어졌다: 특수 경로를 지우고 하나의 일반 exact 경로만 남긴다. 이제 `Fuse`·`Cut`·`Common` 셋 다 씨임 경로를 지난다.

- **★ 핵심 새 로직은 keep/flip 한 줄이었다 — De Morgan 쌍대.** `overlap_fuse_cut`의 keep/flip 표에 `Common => (Side::Inside, Side::Inside, false)`. A∩B는 A의 B-안쪽 재료 + B의 A-안쪽 재료, 둘 다 원래 바깥 법선(Cut과 달리 어느 셸도 cavity 벽이 안 되므로 `flip=false`). Fuse `(Outside,Outside,false)`의 쌍대다. 씨임 재구성 기계는 이미 이 픽스처들에서 검증돼 있었다 — `l_and_corner_box`의 `Cut`(2.776)/`Fuse`가 열린 호를, `cut_the_stub_by_the_l`이 island 면을 `keep=Outside`로 돌린다. Common은 **같은 기계를 `keep=Inside`로** 밟는다.
- **★ `orient_seam_loop`의 `Inside` 부호가 처음 밟혔다.** `material_outside = (keep == Side::Outside)` 매핑이 hole(Cut, 재료가 고리 밖 → true)과 island(Common, 재료가 고리 안 → false)를 정확히 가른다 — 공식은 대칭일 뿐 아니라 의미적으로 옳았고 값만 미실행이었다. 발화: `a_common_can_leave_a_closed_seam_loop`(막대가 큐브를 관통, `∩ = [1,2]²×[0,3]` = `3.0`), cube의 캡에서 유지 영역 `[1,2]²`가 씨임으로만 둘러싸인 닫힌 고리. validate 클린.
- **삭제전 증명((5b)의 규율).** `common`이 살아 있는 동안 `the_seam_path_answers_common_overlap`이 `overlap_fuse_cut(Common, two_boxes)`을 직접 불러 `[0.5,1]³ = 0.125`·정점 8·엣지 12·면 6을 못박았다. n2가 `common`을 지운 뒤 `common_of_two_cubes`가 디스패처로 같은 답을 낸다 — **씨임 재구성이 열거와 같은 위상을 낸다는 실측**(카운트가 안 바뀌었다).
- **삭제.** `common`·`enumerate_vertices`·`build_edges`·`build_faces`·`order_ccw`·`assemble` + `ResultEdge`/`ResultVertex`(~240줄). **`order_ccw`의 `atan2`가 커널의 마지막 열거 술어였다.** 태그 여섯이 발화 테스트 없이 죽었다(볼록 전용 백스톱): `NONCONVEX_OPERAND`·`TANGENT_EDGE`·`NONMANIFOLD_FACE`·`DEGENERATE_CENTROID`·`DEGENERATE_RADIUS`·`COLLINEAR_FACE`. `COMMON_OVERLAP`도 두 자리를 열며 죽었다. `FOURPLANE`은 산다(`overlap_fuse_cut`); `is_convex`·`solid_vertices`는 `detect_coincident_interface`가 붙잡아 (5d)로.
- **OCCT diff** 비볼록 Common(`l_and_corner_box`, 부피+넓이), 44 → 45. **동작 보존**은 살아남은 Common 골든(`common_of_two_cubes` 0.125, 포함, proptest 둘)이 씨임 경로로 그대로 냄으로 확인.
- "M5-d3의 남은 가드" 판단은 그대로 유효하다: M5-d3의 남은 가드는 `SEAM_ACROSS_HOLE_RIM`(3f-6) 하나 — 3g는 그것과 독립이다.

**씨임이 구멍 rim을 가로지른다 — `SEAM_ACROSS_HOLE_RIM` 은퇴, M5-d3의 마지막 가드 (M5 셀 3f-6).** 발화 픽스처는 pocket된 큐브에 코너 박스 `[0.55,1.15]²×[0.85,1.15]`를 `Cut` — 박스 발자국이 pocket 구멍과 `[0.55,0.7]²`에서 겹쳐 벽 `x=0.55`·`y=0.55`가 뚜껑의 구멍 rim을 가로지른다. 그러면 `∂(뚜껑)`이 링 하나가 아니라 여럿이 된다.

- **★ 여덟 번째로 술어가 아니라 이미 있던 기계였다 — 이번엔 인덱스 위상.** 3b·3h·3f-4(술어 셋)와 3f-5(가드를 문 자리로 내림)에 이어. `stitch_cycles`·`boundary_runs`·`run_classes`는 3e-3 이래 "엣지가 아니라 **run**"으로 돌고 정수만 주고받았다. 남은 단일 링 가정은 암묵적 `(i+1)%n` 순환 **여섯 자리**뿐이었다(재검으로 넷→여섯 정정 — succ 안 walk 전진 둘을 놓칠 뻔했다). 그것을 명시적 per-ring 후행 맵 `next[]`/`prev[]`로 바꾸니 새 술어도 좌표 읽기도 없이 다중 링이 떨어졌다. `stitch_cycles`의 `succ`가 곧 평면 그래프의 face-walk(교차=노드, 호·run=엣지)라, **호가 서로 다른 링을 자동으로 잇는다** — 호가 outer의 kd에서 rim의 dk로 건너고 `next`가 각 링 안에서만 걷는다.
- **★ 측정이 가르쳤다: 패리티는 링당 짝수다.** "outer→hole 단일 호가 링을 홀수로 만든다"는 애초 우려는 **틀렸다** — 구멍이 그 호를 둘로 쪼개 짝수를 복원한다. n0 실측(뚜껑): 호 둘(각각 outer↔rim), **outer 2교차·rim 2교차**, 둘 다 짝수·비영. 교차 0인 링(일부 링만 교차, 어떤 링은 통째 kept/dropped)은 per-ring `boundary_runs`가 `runs.last_mut()==None`으로 `SEAM_COUNT_MISMATCH` **정직 거절** — 발화 픽스처는 모든 링이 ≥1 교차라 안 닿는다. 위상 논증이 "짝수"라 했지만 판정은 측정이 했다.
- **교차된 구멍은 outer 경계로 흡수된다.** rim이 씨임에 걸리면 그 구멍은 outer 링에 이어져 genus가 떨어지고 inner loop로 남지 않는다. `crossed`(face.inner 인덱스, **단일 진실원**)가 어느 구멍이 흡수될지 정하고, 그 한 집합이 (연결 verts/kept/bnd/by_edge, per-ring `next`/`prev`, `place_loops` 제외, hole 재매핑)을 모두 구동한다. 흡수 구멍은 `place_loops`에서 빠진다(이중 계산 방지 — outer 링이자 placed hole로 두 번 셀 뻔했다); 미교차 구멍은 오늘처럼 placed. 슬라이스 `boundary_runs`가 0-based run 정점 인덱스를 내므로 연결 `verts`/`bnd`를 타려면 `+base` 오프셋이 필수였다(놓치면 조용한 wrong-vertex, 부피·genus·OCCT 삼중 게이트가 잡는다).
- **결과.** `Cut(pc,box)=0.893`(구멍이 바깥으로 열린 코너 노치, genus 0), `Cut(box,pc)=0.081`. 양쪽 validate 클린·OCCT diff(부피+넓이) 확인 — 넓이가 흡수된 rim 크기를 독립 커널로 못박는다. `SEAM_ACROSS_HOLE_RIM` 태그·`edge_ix`의 rim-거절 가지 삭제, `edge_ix`는 이제 crossed rim 엣지를 연결 슬롯으로 해결한다. OCCT 45 → 47.
- **`NESTED_LOOPS`(닫힌 고리 nesting)는 독립이라 다음 셀(3f-7).** 씨임이 구멍을 **가로지르는** 경우(rim 교차 → 흡수, 이 셀)와 씨임 **닫힌 고리가 구멍을 에워싸** "제 구멍을 가진 섬"을 만드는 경우(3f-5의 "NESTED_LOOPS 과잉 거절" 항목)는 코드가 다르다 — 교차 구멍은 `place_loops` 전에 흡수되므로. 한 번에 한 변수: 후자는 3f-7의 nesting-forest `place_loops`.

**닫힌 고리 nesting — `NESTED_LOOPS` 과잉 거절 은퇴 (M5 셀 3f-7).** 발화 픽스처는 pocket된 큐브를 뚜껑(z=1)과 바닥(z=0.5) 사이 z=0.7에서 슬래브로 자른다(`pocket_and_slab(0.7)`). 슬래브 밑면에서 B 단면이 **동심 두 씨임 고리** — 바깥 큐브 벽(구멍)과 그 안 pocket 벽(섬). 고리 안의 고리를 `place_loops`의 쌍별 검사가 무조건 거절했다.

- **★ 3f-4가 설계한 세 기계가 중첩까지 그대로 확장된다.** 3f-4는 "(1) 국소 부호 규칙이 링 방향을, (2) 포함 판정이 hole/island과 중첩을, (3) 감김이 검산을" 맡긴다. 당시 지운 `classify_loops`는 **감김으로 깊이를 추론**해 호가 있으면 틀렸다. 해법은 **깊이는 containment로, base는 region으로** 나누는 것: `is_island = base_kept XOR (depth 짝수)`, `depth`는 `l[0]`을 담는 다른 고리 수. 그러면 고리가 **중첩 깊이를 따라 hole/island를 교대**한다(region 속 구멍, 그 속 섬, 그 속 구멍). 3f-4 반례(유지 영역 구멍 vs 버림 영역 섬, 둘 다 depth 0)는 **per-loop `base`가 가른다** — 감김을 안 쓰므로 안 헷갈린다.
- **`orient_seam_loop`은 중첩을 이미 옳게 방향짓는다.** 방향 공식이 각 엣지의 이웃 평면 `Q`(상대 솔리드 면 법선 `n_out_g`)를 per-edge로 읽으므로, pocket 벽 위 섬은 재료가 안쪽임을, 큐브 벽 위 구멍은 바깥임을 정확히 반영한다. `Cut(pc,slab)`의 B-면 `keep=Inside`에서도 C를 CCW(섬)·P를 CW(구멍)로 방향지어 깊이 패리티와 일치 — 손으로 검증. 감김 검산(`check_loop_class`)이 셋째 독립 기계다.
- **★ `f`의 구멍은 패리티에 안 넣는다 — 언제나 CW.** 저장 방향이 고정 CW라 재분류하지 않는다(`i >= n_rings`). 패리티에 넣으면 미소유 구멍(base거짓·depth0)이 섬으로 오분류돼 오늘 안 나던 `LOOP_CLASS_MISMATCH`를 낸다. forest는 `f` 구멍에게 **섬 owner만** 준다(§9 665의 "섬 고리 안 `f` 구멍" 경로).
- **분해: `nest_loops`(순수 기하 forest) + `classify_nesting`(순수 조합).** 전자는 `point_in_ring`으로 각 고리의 containers·region을, 후자는 좌표 없이 깊이 패리티·owner를 낸다. 조합이 순수라 depth-2 최심 섬 선택·`f` 구멍 강제·감김 불일치를 **합성 forest 골든**으로 배선 전 못박았다(실물 3-링 픽스처 불필요). 섬은 이제 제 구멍을 담는 `LocalFace`가 된다(옛 `inner: vec![]`).
- **결과.** `Cut` 양쪽이 열린다 — `Cut(slab,pc)=1.488`(region 속 구멍 C 곁에 섬 P), `Cut(pc,slab)=0.668`(슬래브 밑면 dropped·regions 비어 **섬 C가 구멍 P를 가짐**, `flip=true`). validate 클린·OCCT diff(부피+넓이) 양쪽. `place_loops`·`NESTED_LOOPS` 삭제. OCCT 47 → 49.
- **★ `Fuse`는 여전히 거절 — 이유가 바뀌었다.** 합집합이 pocket을 밀봉해(`[0.3,0.7]²×[0.5,0.7]`, z=0.7 슬래브가 뚜껑) **닫힌 cavity**를 만든다 — 씨임 경로가 못 만드는 둘째 셸. `disconnected_result`로 거절(옛 `NESTED_LOOPS` 자리). **중첩 자체는 처리됐고**(면 재구성 성공), 거절이 그 너머 cavity로 옮겨갔을 뿐. cavity는 (5c)의 몫. 계획은 Fuse도 열릴 줄 알았으나 **측정이 cavity를 드러냈다.**

**씨임 경로가 갇힌 cavity를 낸다 — `DISCONNECTED_RESULT`를 부호로 가른다 (M5 셀 5c).** 3f-7이 드러낸 Fuse 밀봉 cavity를 연다. `assemble_fuse_cut`이 `cavities: vec![]`을 하드코딩해 못 내던 것을, 재구성 면을 성분 분해해 outer 셸 + N cavity 셸로 조립한다.

- **★ 부호가 곧 "재료 vs void"다.** `DISCONNECTED_RESULT`가 두 의미를 혼동했다 — 진짜 분리(로드가 둘로, 두 솔리드)와 갇힌 void(한 솔리드 + cavity 셸). 재구성 면은 material이 옳은 쪽으로 방향지어졌으므로, 성분의 **부호 있는 부피**가 가른다: 재료 덩어리는 material이 안 → 법선 바깥 → **양수**, void는 material이 밖 → 법선 안쪽 → **음수**. **양수 성분 정확히 1개면 outer+voids(수용), 아니면 거절.** 이 부호는 `validate::check_cavity_orientation`이 재검하는 그 발산식이라 수용=validate 통과. containment 판정도 새 술어도 없이 부호 하나가 판별·검산을 한다(3f-3 "감김 하나가 중첩을 가린다"의 결). **((5d)-5가 이 결정을 flux에서 exact extreme-vertex 부호로 분리 — 이제 결정과 validate net이 서로 다른 방법. 아래.)**
- **★ 측정: void는 이미 안쪽, 진짜 분리는 양수 둘.** n0 실측 — `Fuse(slab,pc)`: 성분 둘, outer `+7.32`(=2.44)·void `−0.096`(=−0.032). **void가 음수라 `reversed_shell` 불필요** — Fuse가 flip 안 하고 pc 벽·바닥이 이미 pocket(=void) 안쪽을 향한다. `Cut(rod,L)`·`convex_cut_can_sever`: 성분 둘 다 양수 → `positives.len()!=1` → 여전히 `DISCONNECTED_RESULT`. 부호를 **materialize된 Face에서 재므로** assemble의 orient/flip 재현 중복이 없다.
- **하류는 이미 cavity-ready — 변경 0.** `Solid { outer, cavities }`·`reversed_shell`·`point_in_solid`·`mass_props`(void 빼기)·`check_cavity_orientation`·`tessellate`/`reachable`·STEP `solid_with_voids`가 이미 `once(&outer).chain(cavities)`를 돈다. 변경은 nacre-ops 조립부뿐(`faces_connected` → `face_components` + `shell_signed_flux`).
- **결과.** `Fuse(slab,pc) = 2.408`(재료), cavity 1(부피 `0.032`), 셸 2, validate 클린, watertight 메시, OCCT `BREP_WITH_VOIDS` 부피+넓이 확인. 양 순서. OCCT 49 → 50.
- **정직한 빚.** (a) **입력 hollow(`HOLLOW_OPERAND`) — (5c-in)·(5c-in-2)가 완전 해소** (seam front-end + 재구성을 all-shells로; 비껴가는 cut·void 진입 cut 모두 재구성). (b) **부호 분류는 f64 — (5d)-5가 해소** (extreme-vertex exact 부호로 은퇴, 아래). 당시엔 거시적 부피라 견고했으나 위상 결정이라 §9 표 #5로 등재했고, 진리는 여전히 씨임의 exact 술어가 정하며 이 부호는 파생 캐시 위 위상 라벨이었다.

**평면 계수 exact + is_convex tolerance 은퇴 (M5 셀 5d-1).** (5d) exactness sweep 시작. `is_convex`의 `1e-9` tolerance가 무엇을 가리는지 조사하다 전제 둘을 정정했다.

- **★ 회전은 M5 관심사가 아니다.** is_convex tolerance는 "constructed 솔리드 정점의 f64 비공면"을 흡수한다는데, 그 뿌리로 의심된 회전 픽스처가 **없다** — `rotated_l_prism`은 정점 순서 순환 이동(정수 축정렬, 기하 회전 아님)이고 프로덕션도 솔리드 회전이 없다. 기울어진 형상이 필요하면 정수/유리수 좌표로 직접 구성해 면을 정확히 공면으로 둔다. **M5 픽스처는 exact 좌표만.** 회전(모델링 변환)은 도입 시점에 그 f64 오차를 tol 추적으로 별도 설계한다. (예외: `common_tilted_prism` proptest가 f64-기울어진 SketchPlane을 쓰나, 그 결과가 coincident 인터페이스가 없어 is_convex 값이 결과에 무관 — 오염 없음.)
- **★ 측정: Discovered 정점은 is_convex에 한 번도 안 닿는다.** 전체 스위트 계측 결과 is_convex 도달 피연산자가 **전부 Constructed**(`Discovered=0`). tolerance의 명목(Discovered 캐시 흡수)은 **미발화**였다. is_convex는 `detect_coincident_interface`(imprint/coincident)에서만 호출되고 불리언 결과는 대개 비볼록이라, "볼록+Discovered+coincident"가 빈 교집합인 것이다.
- **은퇴.** is_convex의 `sd≤scale(1e-9)`를 **exact `plane_side(pi.tri, v) ≤ 0`**로. `pi.tri`가 RH 법선=`n_out` 정렬이라 부호 직접, tolerance·`n_out` 읽기 삭제. 시그니처 불변(Constructed만 닿으므로 handle 불필요). exact가 tolerance와 **모든 피연산자에서 일치**함을 계측이 확인(과잉 거절 0). Discovered가 닿는 날은 (5d)-2의 `plane_orient`(det4: `sign(det4[Q1;Q2;Q3;P])·sign(det3)`)가 triple로 판정 — 그러나 발화 픽스처가 없어 **미리 안 짓는다**(미발화 machinery 회피).
- **토대: 평면 계수 exact.** `Plane`이 `through_points`에서 정규화(sqrt)해 계수가 반올림이던 것을, **비정규화 `raw`를 함께 저장**해 `coefficients()`가 exact(정의 정점에서 `a·v+d=0` 정확). `normal()`은 저장 unit 바이트 동일이라 좌표(`plane_plane`·`three_planes`)·게이트(`PARALLEL_EPS`)·크기(flux·STEP)·`project` 전부 불변. 술어는 scale-invariant라 큰 계수도 부호 불변(`prop_scaling_a_plane_is_invariant`), 회귀 0. 씨임 경로의 Discovered 교점이 이제 참 평면 위에 정의된다 — coplanar·plane_side triple-sourcing·plane_orient가 딛는 토대.
- **전수조사 표 갱신:** `is_convex`는 은퇴(exact `plane_side`, 표 #1 완료). 남은 (5d): `coplanar`(#2)·plane_side triple-sourcing(#3)·`plane_orient` 광선-free 분류기(#4)·cavity 분류 부호(#5). 상세는 §9 **"(5d) exactness sweep — 전수조사 표"**.

**`coplanar` exact화 — rank-1 계수 술어, `1e-9` 은퇴 (M5 셀 5d-2).** 표 #2. 옛 `coplanar`(ops)는 `const EPS=1e-9`를 반올림된 unit 법선 위 `cross`(무차원, 이미 OK)와 **`distance×2`(절대-길이, scale-비불변)**에 걸었다 — 근접-공면을 false-merge할 수 있는 진짜 f64 진리-읽기. 두 평면이 같은 평면 ⟺ 계수 2×4 행렬 rank-1 ⟺ **여섯 2×2 minor 전부 0**(error-free `two_product`, 방향·scale 무관)으로 대체. `nacre-predicates::planes_coplanar` + geom 래퍼(§3 결정층) + ops 세 호출부 배선.

- **★ 계획은 "축정렬뿐 → 무동작(Case A)"을 예상했으나, n0가 광범위한 불일치를 드러냈다 — 다만 전부 무해였다.** 실측: `old=true new=false` **1024개 구별 쌍**. 뿌리는 회전이 아니라 **float 좌표 박스의 옆면**이다 — 두 스택 박스의 같은 `x=x0` 면이 서로 다른 z-삼각형의 `through_points`로 지어져 f64 raw가 비-bit-비례(기하적으론 같은 평면)라, exact가 tolerance와 갈린다. **그러나 1024개 전부 같은-법선**(`n_out.dot>0`)이고 반대-법선은 **0개**. 즉 갈리는 쌍은 정확히 `:2542`의 반대-법선 인터페이스 필터가 배제하는 쌍이고, 실제 coincident 캡 인터페이스(반대 법선, exact 축 계수)는 일치한다. → **인터페이스 검출·거절 판정 불변.** 전체 스위트 + 50 OCCT diff 그대로.
- **★ `coplanar`는 정점이 아니라 평면을 읽는다.** (5d)-1의 `is_convex`(정점→Discovered/Constructed가 축)와 달리 이 판정은 Surface 계수만 본다. 불리언 결과 피연산자(Discovered 정점)도 seam 면이 원래 exact Surface를 유지하므로 Case A다 — #3를 막던 "연쇄 불리언" 문제가 여기엔 없다. 우려하던 tilted-캡 Case B는 발화하지 않았다.
- **정직한 프레이밍.** 현재 버그 수정이 아니라 **잠재적 'plausible but wrong' 오병합 함정 제거 + scale-불변 + 순수성**. 헤드라인 골든은 z=0 vs z=1e-9(정확히 1e-9 떨어짐) — 옛 절대 tolerance는 `distance≤EPS`로 false-merge, exact는 `minor(2,3)=−1e9≠0`으로 분리.
- **완결성 grep이 둘째 tolerance를 찾았다 — `interface_correspondence`(ops:2574).** coincident-merge가 B링→A링 정점을 `1e-9×크기` scale-상대 tol로 매칭한다. 공면(#2)이 아니라 좌표 대응이라 이 셀 밖 — 표 #6으로 기록(이미 scale-상대라 `coplanar`보다 나음; exact 대안은 정점 triple 동일성).
- **범위 밖 f64(의도적).** `:2542`의 `n_out.dot<0`은 방향 라벨(=rank-1 비례상수 λ의 부호, coarse·견고 — 후속 `sign(raw·raw)`로 exact화 가능). `three_planes`의 `COPLANAR_DET_EPS`·`plane_plane`의 `PARALLEL_EPS`는 **구성 게이트**(좌표 캐시 conditioning, 진리-읽기 아님)라 불변.

**is_convex triple-sourcing + interface_correspondence exact (M5 셀 5d-3).** 표 #3·#6. is_convex가 볼록 판정에서 정점의 f64 좌표 캐시(`v.point`)를 읽던 것을 **정의-기반**으로 바꿨다: Constructed 정점은 좌표가 진리라 `plane_side` 그대로, Discovered 정점은 세-평면 교점이므로 그 implicit point를 `three_plane_orient3d`로 면에 대해 직접 판정(캐시 무독).

- **★ 술어를 새로 안 지었다.** 계획했던 `plane_orient(Q1,Q2,Q3,P)=sign(det4)·sign(det3)`가 **이미 있었다** — `three_plane_orient3d`(=`indirect_orient3d`)가 그것이다. (5d)-1 계획이 "det3 인수 필수(det4만 쓰면 D 부호로 뒤집힘)"라 걱정한 그 인수는 **술어 안에 내재**하고, 교점 V가 순서 불변이라 세 평면 순서도 무관. 새 수치 코드 0, 배관뿐(시그니처 `model+handles`, `solid_vertex_handles` 재사용).
- **★ 발화 = 불리언 결과를 피연산자로 넣는 첫 사례.** is_convex는 coincident 게이트에서만 호출되고 Discovered 정점은 불리언 결과에만 있는데, 스위트 어디도 결과를 후속 피연산자로 안 썼다(연쇄 둘은 hollow라 문 앞 거절). n0 픽스처: 축정렬 정수 박스 둘의 `Common`=[1,2]³ → 8코너 중 **6개가 혼합-평면 교점이라 Discovered**(2개는 원래 코너 Constructed). 이를 셋째 박스에 스택 → coincident 병합이 Discovered-코너 솔리드에 is_convex 발화. **ingestion 클린**(부피 2.0, cavities 0, validate 통과) — 잠재 버그 없음.
- **정직한 프레이밍 — "더 정확"이 아니라 "결정에서 캐시-읽기 제거".** M5에선 triple=캐시(축정렬 교점 f64-exact). 진짜 불일치는 exact-계수 평면이 f64-비표현 점에서 만나는 TIP-류라 M5 밖. 소득: **(a) sweep 원칙 집행(결정이 f64 캐시 안 읽음), (b) Discovered-operand 경로 첫 측정·경화, (c) TIP/회전 대비 정의-기반화.** #1·#2처럼 동작 보존(전 스위트 + 50 OCCT 불변).
- **곁들이 #6 (별도 커밋).** `interface_correspondence`의 `1e-9×크기` scale-상대 정점 대응을 **exact 좌표 동일성 `ap==bp`**로. 계획의 "triple 동일성" 제안은 **부적용**(인터페이스 정점은 Constructed라 정의가 없다) — 대응 코너는 같은 점(bit-동일)이라 좌표 동일성이 정답. 견고했던 tol이라 순수성 청소, 전 fixture Case A. #3과 원인 분리 위해 별도 커밋(먼저 착지).

**`fan_triangles` zero-area 드롭 exact화 (M5 셀 5d-4).** 표 #4. `point_in_solid`의 winding 부호는 **이미 exact**(`ray_triangle_cross`=orient3d 부호 정수합)이고, 이 경로에 남은 유일한 tol은 `fan_triangles`가 근사-영면적 팬 삼각형을 버리던 **상대 `1e-12`**뿐이었다. 이를 exact zero-area로 은퇴.

- **★ 판정 = 세 좌표투영 `orient2d`가 모두 0.** 외적 `(t1−t0)×(t2−t0)`의 세 성분이 곧 세 좌표평면 투영의 `orient2d`라, 면적 0 ⟺ 셋 모두 exact 0. **축-독립**(nonzero 면적이면 외적이 nonzero → 최소 한 투영이 비퇴화 → 오탐 없음)이라 법선을 안 넣어 시그니처 불변. `orient2d`는 이미 in-scope, 새 술어 0.
- **★ 정직한 스코프 — winding의 정점 f64-읽기는 안 건드림.** `ray_triangle_cross`는 정점 f64 좌표를 읽는다(Constructed·축정렬 Discovered면 exact, **회전 Discovered면 캐시≠진실**로 silent-wrong 가능 — #3과 같은 범주). 이건 회전 시 **TIP(오버홀)**의 몫이라 이 셀 밖. 즉 winding 경로를 통째로 닫은 게 아니라 **fan tol 하나**만 닫았다.
- **★ 광선-free 재설계는 기각.** 표 #4의 옛 프레이밍("광선-free 분류기")은 과잉이다: 비볼록 평면 솔리드에 exact 광선 캐스팅보다 단순한 광선-free 방법이 없고(GWN은 초월수 solid angle 합=f64), 방향 휴리스틱은 grazing을 exact 검출(`orient3d==0`)해 소진 시 **정직히 `RAY_DEGENERATE` 거절**(조용히 틀림 아님)이라 바꿔도 exactness 이득 0·복잡도만↑. "고려·기각"으로 기록.
- **정직한 동기 — 현 버그 아님, 회전 대비.** 상대 tol이 nonzero 삼각형을 드롭하려면 `|e1||e2|≥1e12`(엣지 ~1e6+)라 축정렬 소좌표 M5엔 안 생긴다 — n0 계측이 전 코퍼스(nacre-ops 215 + OCCT 오라클 70)에서 **nonzero-면적 발화 0건** 확인 → **Case A**(bit-동일). 진짜 슬리버는 회전/오버홀 기하에서 나오므로 이 셀은 **exactness 위생 + 회전 이식 대비**. 회귀 가드는 대좌표 슬리버 단위 테스트(옛 tol이면 드롭·새로는 keep). 곁들여 방향-독립 팬 삼각형을 `'dirs` 루프 밖으로 hoist(1회 계산·재사용).

**cavity 분류 부호 exact화 — extreme-vertex 부호로 flux 은퇴 (M5 셀 5d-5). (5d) sweep 종료.** 표 #5, sweep의 마지막. `assemble_fuse_cut`이 재구성 성분을 outer(+)/void(−)로 가르던 **f64 signed-volume flux**(`shell_signed_flux`: sqrt-area·나눗셈-centroid·혼합부호 합)를 은퇴. 부호가 곧 위상 라벨이라 결정이 f64 누적을 타던 마지막 자리였다.

- **★ 판정 = lex-최소 정점의 바깥향 부호.** 성분의 lex-최소 정점 `v*`는 볼록 껍질 위라 볼록 코너이고, 재료는 좌표 증가 방향에 있다 → **outward ⟺ `v*`에 인접한 면 중 materialized 바깥 법선 `n_x<0`인 면이 하나라도 있다**. 반평행 x-수직 면이 한 정점에 공존 못하므로 이 ∃-형은 "max `|n_x|` 면 선택"과 동치이면서 부호 비교뿐. 주변 `positives`/`[outer_c] else 거절`/cavities 로직은 불변(라벨 소스만 교체).
- **★ "정수라서"가 아니라 "좌표 산술을 안 해서" exact.** 술어는 좌표로 산술을 안 하고 (a) lex 순서 비교 (b) 축정렬 평면 계수 부호(정확히 ±)만 읽는다 — `Fuse(slab,pc)` void의 `0.3/0.7`처럼 **f64 비표현 좌표에서도 exact**(flux는 좌표를 곱·sqrt해 이 강점이 없다). 홀(inner loop)도 무관(극점은 outer 루프에만·판정은 법선 부호뿐). flux가 절대 거절 안 하던 **total**성도 유지(∃가 항상 결정).
- **★ 결정과 net이 이제 서로 다른 방법 — defense-in-depth.** (5c)는 결정 부호와 `validate::check_cavity_orientation`이 **같은 발산식**이었다. #5 후 결정은 extreme-vertex, validate net은 여전히 flux — **독립 교차검증**이 됐다(우연 일치 아님). n0가 전 코퍼스에서 둘 일치 확인.
- **고려·기각.** ~~containment-parity~~(대표점 감김 홀짝): 감김 재사용은 매력이나 `point_in_solid`가 `RAY_DEGENERATE`로 **최종 조립에 유효 입력 거절을 주입**(total→partial)·O(n²)이라 기각. ~~exact signed-volume 합~~: 다항 합의 exact 부호 인프라가 없어 f64 누적 부활, 기각(오버홀 몫).
- **측정·회귀.** n0: `is_shell_outward` vs `shell_signed_flux>0`가 전 코퍼스(nacre-ops 217 + OCCT 70)에서 **부호 뒤집힘 0** → **Case A**(같은 outer·cavity 집합·순서). `Fuse(slab,pc)` cavity 1·`Cut(rod,L)` `DISCONNECTED_RESULT` 불변, OCCT `BREP_WITH_VOIDS` 그대로.
- **★ sweep 종료 — 완결성 감사가 근거.** (5d) 이력이 초기 전수조사의 누락을 반복 증명했으므로(#6은 (5d)-2 grep이, (5c-in) 다섯째 자리는 n0가) 종료 전 nacre-ops 프로덕션을 감사: **tolerance 리터럴 0**(`vertex_tol`만, 존재론적 제외), 나머지 f64 부호 비교는 전부 exact 술어 부호(`orient3d`/`plane_side`/`planes_coplanar`)·축정렬 구성-방향·축정렬 법선 dot 부호 — **tol/누적 기반 위상 판정 0**. #1–#6 전부 완료, sweep 정당하게 종료.

**hollow 입력 operand — seam이 void를 비껴가면 지원 (M5 셀 5c-in). `HOLLOW_OPERAND` 은퇴.** 커버리지 매트릭스의 빨강 하나. cavity는 이미 1급(결과-조립 5c·하류 validate/props/STEP이 `once(&outer).chain(cavities)`로 순회)이라, 갭은 **seam 프론트엔드가 outer 셸만 도는 것**뿐. `solid_shell_handles`(outer+cavities)로 다섯 자리를 all-shells로: `collect_planes`·`solid_vertex_handles`·`edge_incidence`·`solid_face_rings`·**면-재구성 루프**. void가 손 안 대고 결과로 실려 (5c) 성분 분할이 보존.

- **★ 놓친 다섯째 자리 — 재구성 루프.** 4함수만 all-shells로 바꾸고 측정하니 여전히 **silent-drop**(코너 Cut 부피 `26.875`·cavities=0·validate 클린). cavity 면이 `planes`엔 들었으나 `overlap_fuse_cut`의 재구성 루프(`solids.get(solid).outer`만 순회)가 결과로 안 실었다. 그 루프까지 all-shells로 하니 **정답 `25.875`·cavities=1**. 탐색 목록(4함수)이 놓친 자리를 n0 실측이 잡았다.
- **★ silent-wrong이 loud해졌나 — 그게 문을 여는 조건.** 옛 `HOLLOW_OPERAND`가 막던 건 crash가 아니라 조용한 오답이고, 그 조용함을 잡던 게 outer/cavity **비대칭 TUNNEL**이었다. 둘 다 all-shells면 비대칭이 사라져 TUNNEL이 재구성 버그를 못 잡는다. 그래서 **cut이 void에 *들어가면*(seam 정점이 cavity 벽에 앉으면) 재구성 전에 명시적으로 `SEAM_ENTERS_CAVITY`로 거절** — 조용히 틀릴 자리를 원천 차단. `PlaneInfo.is_cavity`로 판정.
- **하위-케이스 사다리.** **(A) far-cut**(seam이 void 비껴감, pocket 부품의 다른 곳 bore) = 지원, 정답. **(B) cut-enters-cavity**(void 진입/관통) = `SEAM_ENTERS_CAVITY` 거절, 재구성은 (5c-in-2). n0가 (A) 통과·(B) 거절(silent 아님)을 배선 전 확정.
- **소득·회귀.** 첫 cavity-생성 연산 뒤에도 **연쇄 불리언**이 된다(hollow 결과를 다시 피연산자로 — (5d)-3 ingestion 위에 얹힘). 비-hollow는 `chain(빈)`이라 불변, 전 스위트 + 50 OCCT 그대로.

**cut-enters-cavity 재구성 — hollow 이야기 완결 (M5 셀 5c-in-2). `SEAM_ENTERS_CAVITY` 은퇴.** (5c-in)이 far-cut을 열고 void에 *들어가는* cut을 좁게 거절한 것을, 이 셀이 마저 연다.

- **★ 배관 갭은 하나였다 — 재구성의 seam-*세그먼트* 경로.** (5c-in)이 seam-*정점* 경로(pierced_faces+solid_face_rings)를 all-shells로 했지만, 재구성의 세그먼트 경로 `reconstruct_face_paths→seam_paths_on→arrange::seam_segments_on`이 `other.outer`만 읽었다(단일 잔여 outer-only, 콜체인 검증). 정점 경로는 cavity seam을 찾는데 세그먼트 경로는 못 만들어 count가 어긋난 게 `seam_count_mismatch`의 뿌리였다. :123을 all-shells로 하니 두 경로가 일치.
- **★ 셸-병합이 그냥 됐다.** 밖에서 void로 뚫으면 outer가 cavity에 연결돼 void가 밀폐를 잃는다(genus 변화). 걱정한 `shell_signed_flux` 병합-부호 버그는 없었다 — `face_components`가 채널을 통해 outer와 cavity 면을 한 성분으로 묶고 병합-flux가 거시적 양수라 outer 라벨·cavities=0. 3f-6의 hole-rim 방향 bookkeeping이 **Reversed cavity 면**의 새 구멍에도 옳게 걸렸다(측정: validate 클린).
- **★ n0가 클래스 전체 안전을 확정.** 가드를 걷으면 모든 하위-케이스가 노출되므로 배선 전 전수 측정: **blind hole 25.96·cavities=0·shells=1**, **through-tunnel 25.92(genus-1)**, l_and_inner stub 2.484 — 전부 정답·validate 클린. **slab이 박스를 가르면 `DISCONNECTED_RESULT`**(loud). **silent-wrong 0, seam_count_mismatch조차 없음.** 클래스가 전부 정답-이거나-loud라 `SEAM_ENTERS_CAVITY`·`PlaneInfo.is_cavity` 완전 은퇴.
- **재사용.** (5c) 성분 분할이 병합을 표현하고, 3e-3 genus-1 drilling(BoundaryRuns·run-class·엣지당 다중 교차)이 관통을 재구성 — 새 코드는 seam_segments_on 한 자리뿐. 회귀 0(비-hollow·far-cut 불변).

**split-face 슬롯 지원 — 공면 가드를 공유-엣지 쌍으로 좁힘 (M5 셀 coplanar-narrow).** `has_coplanar_pair`(문 앞 가드)가 인접성 무시하고 **떨어진 공면 면**까지 거절해, 전폭 홈을 판 바(top이 두 스트립으로 갈라져 한 `Surface`를 공유)가 두 번째 불리언을 못 태웠다. 각인(niche)이 아니라 정상 다중-피처 부품이 진짜 갭이었다(각인은 명시적 `ImprintSketch`만 만든다). 세 개 맞물린 수정으로 연다.

- **★ A′ 정적 감사가 셀의 진짜 깊이를 배선 전에 쟀다.** "공면 면을 distinct 인덱스로" 하면 `order_along`/`side_of`(둘 다 `three_plane_orient3d`)가 **공면 쌍 `(i,j)`에 0을 반환**한다(`V_i=P∩Q∩R_i`가 `R_j`에도 놓임). 그 0을 "축퇴/일치/평면 위"로 읽는 소비자가 **~12곳**(3 FOURPLANE 사이트 + `LOOP_ORIENT_MISMATCH`·`POINT_ON_RING`·`SEAM_COUNT_MISMATCH` 하류) — "한 평면당 한 면" 불변식이 seam 기계에 pervasive하다. **그러나 기하가 대부분을 막는다:** `edge_crosses_face` 봉쇄로 엣지는 공유 평면을 한 점에서만 통과하고 그 점은 최대 한 스트립 안이라, 두 disjoint 스트립은 같은 `third`/ring에 co-appear하지 않는다. **기하 무관하게 발화하는 건 전수-스캔 FOURPLANE(lib:1178) 하나뿐** — n3 gate가 이 예측을 실증(슬롯이 완전 재구성, 남은 11곳 미발화).
- **A. per-face 인덱싱.** `surf_ix`를 `Surface→usize`에서 **`Face→usize`**로 (PlaneInfo에 `face` 필드 추가). 두 코플래너 스트립이 한 `Surface`를 공유해 인덱스가 **덮어써져 붕괴**하던 것을 면마다 distinct 인덱스로. 현 통과 경로의 모든 면은 distinct Surface(전단사)라 **동작 보존**(전 스위트 + 50 OCCT 불변).
- **B. FOURPLANE 좁힘.** seam 정점 삼중 `(e0,e1,entry)`에 대해, 셋 중 하나와 `planes_coplanar`인 평면 `m`은 진짜 4번째가 아니라 같은 평면의 다른 면이므로 skip. 가드가 아직 닫혀 no-op(동작 보존).
- **C. 가드 per-edge 좁힘.** `has_coplanar_pair`(임의 두 공면)를 **`solid_has_coplanar_neighbour_edge`**(엣지를 공유하는 두 공면만)로. 각인 rim·Fuse flat-edge(엣지 공유)는 계속 거절, 슬롯의 두 disjoint 스트립(엣지 미공유)은 통과. `overlap_fuse_cut`은 cross-operand 공면 검사(`cross_coplanar`, 접선 접촉) 유지. `has_coplanar_pair`는 이제 테스트 전용 특성화 술어.
- **정직한 회귀 0.** 각인 거절 3종 불변(rim이 엣지 공유 → 새 가드가 계속 거절; `an_imprinted_cube_is_convex`에 새 가드 단언 추가). 발화 테스트 `a_slotted_bar_chains_through_a_cut`(슬롯 바에 blind pocket) + OCCT diff `slotted_bar_pocket_cut_matches_occt`. **Fuse 결과 인접 공면 병합은 다음 셀(fuse-coplanar-merge)** — 여전히 `COPLANAR_PAIR`.

**coincident Fuse 옆면 병합 + 조합적 dissolve — 붙인 부품이 연쇄된다 (M5 셀 fuse-coplanar-merge).** coplanar-narrow의 새 가드가 세운 다음 갭. 면-대-면으로 붙인 등-footprint 볼록 solid의 Fuse(`coincident_merge`)는 두 interface 면만 버리고 **옆면을 분리 유지**해, 결과에 **가짜 엣지**(공면-인접 면 쌍)가 남았다(`fuse_stacked_cubes` = 10면). 그 결과를 두 번째 seam 불리언에 넣으면 `solid_has_coplanar_neighbour_edge`가 거절 → 붙인 부품을 더 못 가공. (overlap/seam 경로는 `cross_coplanar`가 애초에 막으므로 가짜 엣지 생산자는 `coincident_merge`뿐.)

- **병합 = 순수 Node-링 splice, 조건부.** interface 엣지마다 A-옆면·B-옆면의 loop_nodes를 공유 엣지 자리에서 이어붙인다. `assemble_fuse_cut`의 `edge_for`가 병합 엣지를 Line으로 생성·dedup하므로 새 Vertex/Edge 수작업 0. **`planes_coplanar`일 때만 병합** — 비-공면(기울어진 볼록 solid)은 진짜 이면각이라 분리 유지(무조건 병합은 형상 파괴). 축정렬만 구성 가능한 현재는 항상 공면이나, 조건부라야 모든 볼록 입력에 정확.
- **★ face-only 병합은 연쇄 불가 — gate가 실증, dissolve 필수.** 병합만 하면 6면이 되나 **interface 코너가 straight angle**(쪼개진 수직 엣지 양쪽이 같은 이웃 면)로 남아, `loop_triples`(arrange.rs)가 `a == b`로 거절한다. `face_rings`가 **모든 면**의 triple을 만들므로 cut 위치와 무관하게 `LOOP_ORIENT_MISMATCH` — 상단 pocket cut조차 실패. face-only는 canonical하지도 연쇄 가능하지도 않은 중간 형태였다(커밋 a15c3af가 기록으로 남음).
- **★ dissolve는 조합적 — 부동소수점 collinearity 없이.** interface 정점 `v_i`(인접 엣지 `e_{i-1}`,`e_i`)는 **두 엣지가 다 병합(공면)됐을 때만** collinear다: 그때 A의 수직 엣지 = `Π_i ∩ Π_{i-1}`, B의 = 같은 `Π_i ∩ Π_{i-1}`라 **동일 직선**(증명, exact). 그런 `v_i`를 병합 loop_nodes에서 필터하면(splice 구조상 far corner에 flanked돼 있어) 쪼개진 수직 엣지가 하나로 융합 → **깨끗한 6면/8정점/12엣지 박스**. 새 collinearity 술어 신설 불요, 별도 dissolve 셀 폐지.
- **회귀 0·net.** stacked·side-by-side 다 깨끗한 박스, 연쇄 정상. `fuse_stacked_cubes` 6/8/12·공면-인접 0 단언, `a_fused_stack_chains_through_a_cut`(interface 가로지르는 Cut), OCCT diff `stacked_fuse_then_cut_matches_occt`(fuse+연쇄 부피). imprint hole은 same-normal이라 interface로 안 세어 불변(`coincident_merge_keeps_an_imprinted_hole`). **면 수는 OCCT와 비교 안 함** — OCCT 기본 Fuse는 `UnifySameDomain` 없이 defeature 안 해 10면을 내므로(실측), 6/8/12는 우리 가드가 요구하는 커널 canonical 목표이지 OCCT 파생이 아니다.

**면-국소 스케치 포함 가드 — silent-wrong을 문 앞에서 닫음 (M5 셀 imprint-containment).** M5 진행 중 발견한 M4 갭. `prepare_face_split`(imprint/pad/pocket 공유 관문)이 프로파일을 면 centroid에 놓고 포함 검사 없이 inner-loop 구멍으로 만들었다. 프로파일이 면 region을 벗어나면(경계 가로지름·완전 외부·오버행·경계 접촉) **조용히 무효한 면**을 냈다 — DNA("조용히 틀리는 것을 가장 경계") 정면 위반.

- **★ n0가 창의 실재·크기를 못박았다.** 네 경우 × imprint/pocket 전부: **`apply` Ok · `validate` 통과(위상만) · `props` garbage · `tessellate`의 `NoEar`만 발화**. props는 링을 그냥 적분해 오답을 낸다 — pocket/오버행은 **부피 −3.8(음수)**, pocket/외부 0.952, imprint는 부피-중립이라 **1.0(그럴듯한 오답, 가장 위험)**. props만 읽는 경로는 조용히 틀린다. tess §10 기록(§7 "정직한 실패", "NoEar만이 찾아냈다")의 그 창이 이제 닫힌다. **경계 접촉(flush)도 `NoEar` 발화 → 엄격 포함이 정답**임을 n0가 확정.
- **가드 = exact 엄격 포함.** `prepare_face_split` 입구(모델 변형 전)에서 프로파일이 면 region(outer 안·모든 hole 밖·경계 비접촉)에 엄격 포함인지 판정, 아니면 `OpError::ProfileNotContainedInFace`로 거절. **판정은 exact여야 한다** — false-accept(경계 바로 밖 통과)가 곧 silent-wrong이므로. 면 법선의 최대-성분 축을 **drop한 exact 투영**(좌표를 버릴 뿐 재계산 없음) 후 `orient2d`(Shewchuk, nacre-predicates 신설·geom 재-export)로 point-in-ring·segment-meet. false-reject(정상 거절)는 안전, false-accept만 창을 재개방.
- **계약·회귀 0.** imprint/pad/pocket 계약 = "프로파일 ⊂ 면 region". contained 프로파일 전부 불변(기존 스위트 net). 발화 테스트 `a_profile_reaching_past_the_face_is_rejected`(네 경우 × 세 연산). **★ 오버행/일반 pad은 "미지원"이 아니라 별개 방향** — 실 CAD(FreeCAD 포함)처럼 스케치를 프리즘으로 `Extrude` → `Fuse`/`Cut`(boolean-based pad, 오버행 공짜). 그 관문은 불리언의 **공면 부분 접촉**(프리즘-온-페이스)이라 coplanar-contact 셀들의 다음 스레드(coplanar-narrow→fuse-coplanar-merge→coplanar-contact). 프로파일 **자기교차**(Profile2d "미검사")는 별개 silent-wrong(containment 아닌 프로파일 유효성)이라 별개 셀.

**boss를 면에 Fuse — contained 공면 접촉 (M5 셀 coplanar-contact-boss).** boolean-based pad(오버행) 스레드의 첫 칸. 박스가 다른 박스의 면 **안쪽**에 얹히면(공면·반대법선, footprint가 면에 포함) `detect_coincident_interface`가 **정확 링 일치**를 요구해 지나침 → seam-free `contained_result`에서 `point_in_solid`가 base top 평면 위 boss 바닥 정점에 grazing → **`vertex_on_face_plane` 정직 거절**(n0: silent-wrong 아님).

- **★ coincident_merge를 contained로 일반화 — 직전 두 셀 합성.** 정확-일치 merge(두 interface 면 버리고 stitch)와 달리, 큰 면은 **경계를 유지하고 작은 footprint를 hole로** 얻는다: containing 면 → annular(기존 outer + contained footprint를 inner hole), contained 면 drop, footprint 엣지가 hole과 boss walls를 stitch(`assemble_fuse_cut`의 `edge_for` dedup), 나머지 유지. footprint 영역은 양쪽 solid라 내부(면 없음이 정답). **감지 = `detect_coincident_interface`의 공면·반대법선 쌍 + `profile_strictly_in_region`**(imprint-containment 셀 신설)로 strict 포함 판정. **조립 = `solid_local_faces`·`assemble_fuse_cut`**(coincident_merge와 동일). 두 직전 셀이 building block.
- **★ exact — Profile2d 반올림 회피.** hole = contact 면의 **자기 정점** 그대로(프레임 투영·재구성 안 거침)라 stitch가 bit-exact, remap 불요. **hole 방향 = contact 면 outer loop as-is**: 그 면 법선(-n) 기준 CCW = containing 면 법선(+n) 기준 CW = hole 규약(보는 법선 반전 = CCW↔CW). 뒤집으면 매니폴드 깨짐.
- **범위·net.** **Fuse만**(Cut/Common은 fall-through). contained footprint·볼록. 정확-일치 merge와 strict-포함으로 명확히 갈림(경계 접촉·일치는 이 경로 밖). 발화 `fuse_a_boss_onto_a_face`(부피 1.25·validate·watertight·공면-인접 0 연쇄), OCCT diff `boss_fuse_then_cut_matches_occt`(fuse+연쇄 cut). **겹침·오버행**(footprint가 면 경계 가로지름 — partial split)·**boolean-based pad 통일**은 후속 하위-셀.

**blind pocket을 면에 Cut — boss의 Cut 짝 (M5 셀 coplanar-contact-cut).** boss(반대-법선 Fuse)와 대칭: 프리즘 B가 base A **안**에 top-flush(**같은-법선**, footprint 포함)로 있고 `Cut(A,B)` → blind pocket. **법선 방향이 side를 인코딩** — 반대=B 외부(boss), 같음=B 내부(pocket). 현재 boss와 동일하게 `vertex_on_face_plane` 거절(n0).

- **★ boss 조립을 `cut` 플래그로 일반화.** `contained_contact_result(cc, cut)`: containing 면 → annular(footprint hole), contact 면 drop은 공통. **Cut이면 B 면들 `flip:true`**(제거 영역 경계 — walls 안쪽, far 면 = pocket floor 위쪽; `assemble_fuse_cut`의 seam-Cut flip 재사용) **+ hole 뒤집기**(same-normal이라 B contact outer가 containing 법선 기준 CCW → CW로). boss(반대-법선)는 as-is·flip 없음. 두 반전이 `cut`으로 함께 켜짐.
- **★ 감지 = same-normal + B blind 포함.** `detect_pocket_contact`: 단일 same-normal 공면 쌍, A(Cut의 유지 피연산자) 면이 B footprint 포함, **B가 A에 blind**(접촉-평면 밖 B 정점이 전부 A strict 내부). too-tall B(A를 뚫음)는 seam cut이라 이 경로 밖 — `boundaries_intersect`는 flush 접촉을 오탐하므로 못 씀, 정점-내부 검사로. 발화 `a_pocket_that_punches_through_is_not_a_pocket_contact`.
- **net.** `cut_a_blind_pocket_into_a_face`(부피 0.875·validate·watertight·공면-인접 0), OCCT diff `pocket_cut_then_cut_matches_occt`(pocket + 비볼록 결과에 transversal 연쇄). boss·정확-일치·seam 불변. **겹침·오버행·무-op 조합(반대+Cut, 같음+Fuse)·Common**은 후속.

**오버행 boss를 면에 Fuse — 단일-엣지 공면 부분 겹침 (M5 셀 coplanar-contact-overhang).** boss/pocket은 footprint가 면에 **통째 포함**돼 "면 하나가 통째 hole"로 끝났다. 오버행은 footprint가 면 경계를 **가로질러** 공유 평면이 **부분 겹침**(2D 폴리곤 boolean)이다 — imprint-pad가 못 내는 진짜 새 커버리지(캔틸레버). 첫 하위-셀은 가장 tractable한 **단일 base-엣지 오버행 boss `Fuse`**만. n0: 대표 픽스처(base `[0,1]³`, boss `[0.5,1.5]×[0.25,0.75]×[1,2]`)가 boss/pocket과 동일하게 `vertex_on_face_plane` 정직 거절(silent-wrong 아님).

- **★ single-edge 제약이 2D-boolean을 반평면 1개 클립 + 루프 splice로 축소.** 볼록 footprint Q가 면 P의 **엣지 e 하나만** 가로지르면(proper 교차 정확히 2개, P 꼭짓점이 Q 내부에 없음), 세 영역이 e의 반평면 하나로만 갈린다: overlap=`Q∩e내부`(내부, 면 없음), **cantilever**=`Q∩e외부`(볼록, small 접촉면 −n), **notch**=`P−overlap`(P 루프의 e 구간을 Q 안쪽 호로 치환한 비볼록 단일 폴리곤, +n). 일반 폴리곤 엔진 불요. **P/Q 식별**=containment(오버행엔 부적용) 아닌 "두 교차점이 자기 단일 엣지 위에 놓이는 면=P". NEITHER/BOTH(코너·다중-엣지)면 out-of-scope `None`.
- **★ 교차 정점 = seam 메커니즘 재사용(exact).** 교차점은 P접촉·P벽·Q벽 **3평면 meet**(`three_planes` 좌표 캐시 + `vertex_tol` residual) → `Origin::Discovered{ThreePlane}`을 `Node::Seam(triple)`로 참조, `assemble_fuse_cut`의 기존 Discovered 경로가 vertex/edge dedup·stitch. 발견된 교차라 Constructed 위장 안 함(DNA #4). 벽 조회는 `boolean`이 adjacency 미구축이라 `model.adjacency` 무의존 — 접촉 엣지의 `Handle<Edge>` 공유 면을 shell 직접 순회(`edge_incidence` 방식).
- **★ 벽 re-split은 3D 평면 가드 필수 — T-junction 방지.** boss 벽이 접촉 영역보다 넓어 교차점이 벽 top/bottom 엣지 **내부**에 떨어지므로 그 벽을 쪼개야(안 그러면 non-manifold). 판정은 벽 법선 drop 후 `orient2d==0 && in_bbox`인데, **투영이 축 하나를 버려** 다른 depth의 교차점이 오탐된다(첫 실측: boss top z=2 면이 z=1 교차점을 (x,y) collinear로 흡수 → `VertexOffSurface` residual 1.0). **먼저 `plane_side==0`(벽 평면 위)로 3D 걸러야** 정확. collinear 경계 정점 허용(imprint/slot과 동일).
- **net.** `fuse_an_overhanging_boss`(부피 1.5·validate·watertight·공면-인접 0 — notch/cantilever는 점만 공유해 연쇄 가능), OCCT diff `overhang_fuse_then_cut_matches_occt`(fuse 1.5 + 캔틸레버 관통 드릴 연쇄). 게이트 발화 `a_boss_that_pierces_the_base_is_not_an_overhang`·`a_corner_overhang_is_not_a_single_edge_overhang`. boss/pocket/coincident/seam 불변. **후속(비목표)**: 코너·다중-엣지 오버행, Cut/Common 오버행, 공면 옆벽 일치(공유 엣지), non-axis-aligned·4-plane 동시성, self-intersecting profile — 전부 감지 `None`으로 정직히 흘림.

**오버행 Cut으로 모서리 슬롯 — 벽을 뚫고 나오는 edge-slot (M5 셀 coplanar-contact-overhang-cut).** 오버행 Fuse의 Cut 대칭 짝. 프리즘 B가 base A 면에 top-flush(**같은-법선** 공면)이나 footprint가 A 면의 한 엣지를 가로질러, `Cut(A,B)`가 A의 **모서리에 슬롯**을 판다. 스레드 목표 pad/pocket의 pocket 쪽. n0: 대표 픽스처(A `[0,1]³`, B `[0.5,1.5]×[0.25,0.75]×[0.5,1.0]`)가 **`vertex_on_face_plane` 거절** — **예측(`COPLANAR_PAIR`) 반증**: B가 x=1 벽을 관통해 seam 경로로 가지만, 그 전에 `point_in_solid`가 공면 top(z=1) 위 B 정점에 grazing해 먼저 거절(Fuse 오버행과 같은 태그, 다른 원인).

- **★ R=B∩A는 단일-반평면 클립.** 제거 영역은 볼록 프리즘 B를 **가로지른 A-벽 W의 반평면(inside)**으로 자른 것(Fuse 캔틸레버가 단일-직선 클립이었듯). R의 면이 두 종류로 갈린다: **A 경계에 얹힌 면 → notch**(R-top=A top의 mouth, R-side=A 벽 W의 옆 개구 — 둘 다 R이 A 엣지에 닿는 곳에서 열림), **B 내부 벽 → 슬롯 표면**(x≤1로 클립, `flip:true`). 그래서 Fuse보다: **notch 2개**(mouth+side, 둘 다 `splice_notch` 공유·직전 셀에서 팩터링), seam 정점 **4개**(2개 z=1 접촉엣지, 2개 z=0.5 슬롯 바닥), B 벽은 `clip_face_by_plane`(3D 반평면)으로 재구성. **resplit 불요** — 클립이 정점을 코너에 인라인 넣어 T-junction 없음.
- **★ coplanar-twin 크랙 회피 — 정규 교차 집합.** 접촉 평면 z=1에는 A-top(P)·B-top(Q) 두 coplanar twin이 있어, top 교차점을 mouth는 `{W,B-y,P}`로·B벽 클립은 `{W,B-y,Q}`로 만들면 **같은 점에 정점 2개 → 균열**. `try_overhang`이 top 교차점을 **P 기준으로 canonical하게**(three_planes 첫 인자=P) 만들고 클립은 **점 좌표로 이 집합만 조회**(새 triple 생성 안 함)해 하나의 `Node::Seam`으로 봉합. Fuse엔 접촉 평면이 하나라 없던 문제.
- **★ single-wall breakthrough 게이트.** B가 A를 오직 W로만 뚫음을 **off-plane B 정점 v마다 `point_in_solid(v,A)==Inside` ⟺ `plane_side(W,v)<0`(W inside)** 등가로 exact 판정. B가 바닥/타 벽을 뚫으면 등가가 깨져 거절(`an_edge_slot_through_the_bottom_is_out_of_scope`).
- **net.** `cut_an_edge_slot`(부피 0.875·validate·watertight·공면-인접 0), OCCT diff `edge_slot_cut_then_cut_matches_occt`(슬롯 0.875 + 비볼록 결과에 transversal 드릴 연쇄). boss/pocket/overhang-Fuse/coincident/seam 불변. **후속(비목표)**: through-slot(바닥까지 관통), 코너(2-wall) 오버행 Cut, Common, non-axis-aligned·4-plane 동시성.

**코너 오버행 Fuse — 2-crossing 겹침으로 일반화 (M5 셀 coplanar-contact-overhang-corner).** single-edge Fuse를 **∂P·∂Q가 정확히 2점에서 만나는 볼록 겹침(단일 lens)**으로 일반화. boss footprint가 base 면의 한 변이 아니라 **꼭짓점을 물고 두 변**을 가로지르면(코너 삼킴) L자 캔틸레버. single-edge는 그 특수케이스(삼킨 P 정점 0). n0: 이 픽스처는 M5-d2 기록(§10 위 "태그 하나가 움직였다")의 `coincident_merge_rejects_offset_footprint`와 **동일** — 그때 `vertex_on_face_plane`로 거절되던 것을 이 셀이 **열어 성공으로** 바꾼다(예측 반증; 테스트는 `detect_coincident_interface` None 단위검사로 개명, 성공은 `fuse_a_corner_overhanging_boss`).

- **★ 대칭 2-arc가 핵심.** 두 교차점이 P·Q 경계를 각각 두 호로 가른다: **notch = `P_outside 호` + `Q_inside 호`**, **cantilever = `Q_outside 호` + `P_inside 호`**(교차점이 이음매). single-edge면 `P_inside`가 빈 호라 cantilever=`Q_outside`(옛 동작 보존). corner면 `P_inside`=삼켜진 base 코너(캔틸레버 L의 reflex 꼭짓점, **기존 Node::Orig**라 base 벽과 자동 공유). `split_loop_at_crossings`(각 면을 outside/inside 호로 — **호 선택은 has_outside**; single-edge P_inside가 빈 호라 has_inside론 못 고름) + `stitch_arcs`(두 호를 교차 노드로 이음, splice_notch 일반화).
- **★ 회귀 보호 두 가지.** ① **same-edge param 정렬**: single-edge는 두 교차점이 한 P 엣지라, aug 삽입을 param 순으로 안 하면 `Pv[seg]→C_far`가 `C_near`를 관통해 비매니폴드 — `split_loop_at_crossings`가 같은-엣지 교차점을 엣지 방향으로 정렬. ② **Cut 격리**: `try_overhang` 완화("정확히 2 교차"만, per-crossing p_seg)가 `detect_overhang_cut_contact`로 새지 않게 Cut이 `cs[0].p_seg == cs[1].p_seg` 재검사(single-edge=single-wall 유지). P/Q 배치는 **대칭**(swap해도 같은 두 면).
- **net.** `fuse_a_corner_overhanging_boss`(부피 2.0·validate·watertight·공면-인접 0), OCCT diff `corner_overhang_fuse_then_cut_matches_occt`(2.0 + L 캔틸레버 코너 관통 연쇄). single-edge Fuse·Cut·boss/pocket 불변. 당시 4-교차 슬래브는 out-of-scope로 거절(다음 셀 -multi가 엶). **2-crossing 프레이밍이라 한 변 걸쳐 2코너 삼킴도 자동 커버.** **후속(비목표)**: multi-edge(≥4 교차·가로지르는 슬래브), 코너 Cut, Common, non-axis-aligned·4-plane.

**슬래브 오버행 Fuse — 일반 N-crossing 겹침 (M5 셀 coplanar-contact-overhang-multi).** overhang Fuse를 **∂P·∂Q가 2k점에서 만나는 일반 convex-convex 겹침**으로 일반화. boss가 base 면을 **가로지르면**(슬래브) notch/cantilever가 **여러 조각으로 분리**된다(예: base top이 두 스트립, boss 밑면이 두 처마). single-edge·corner는 1조각 특수케이스. n0: 직전 셀의 `a_spanning_slab_is_not_a_2_crossing_overhang`(4-교차 거절)를 이 셀이 **열어 성공으로** 바꾼다(예측 반증; 성공 `fuse_a_spanning_slab_boss`).

- **★ 볼록 교차점 같은-순환-순서 → 연속-쌍 매칭.** 볼록 P·Q의 교차점은 양 경계에서 같은 순환 순서로 나타난다. ∂P·∂Q를 교차점에서 쪼개 각 호를 inside/outside 라벨하면, **notch 조각 = P-outside 호 + 그 두 끝 교차점을 공유하는 Q-inside 호**, **cantilever 조각 = Q-outside 호 + P-inside 호**. (같은 방향 "연속쌍 P-호+Q-호"는 overlap/external이라 면 아님 — 매칭은 **ends-집합 동일 + inside**. 2-crossing은 두 호가 같은 ends라 inside 필터로 disambiguate → 현 corner 동작 보존.) 조각 수 = 2k. `split_loop_all_arcs`(모든 호+ends+outside) + `stitch_arcs`(직전 셀 **그대로**).
- **★ 회귀·격리.** `try_overhang`을 `[Crossing;2]`→**`Vec`(2k)** 완화. edge-slot Cut은 `try_overhang`을 공유하므로 **`cs.len()==2 && same p_seg`** 재검사로 single-wall 유지(격리). P 홀은 **포함하는 notch 조각에 배정**(다조각에서 모호 회피, 1조각이면 자동). 매칭 실패(비볼록이 is_convex를 빠져나감)는 `OVERHANG_ARCS`로 정직 거절.
- **net.** `fuse_a_spanning_slab_boss`(부피 1.4·validate·watertight·공면-인접 0·연쇄), OCCT diff `spanning_slab_fuse_then_cut_matches_occt`(1.4 + 처마 관통 연쇄). single-edge·corner·edge-slot·boss/pocket/coincident 불변. **후속(비목표)**: P 홀 base 다조각, 코너·multi Cut, Common, non-axis-aligned·4-plane.

**슬래브 채널 Cut — 다벽 오버행 파내기 (M5 셀 coplanar-contact-overhang-slab-cut).** 파내기(Cut) 일반화 첫 칸. edge-slot(1벽)을 **다벽**으로 — 프리즘이 base를 가로질러 **반대 두 벽**으로 나가는 채널 슬롯. edge-slot은 1벽 특수케이스로 통합. n0: 슬래브 Cut 픽스처(A `[0,1]³`, B `[-0.5,1.5]×[0.4,0.6]×[0.5,1]`)가 edge-slot과 동일하게 `vertex_on_face_plane` 거절(공면 top grazing 먼저).

- **★ 세 부분 재구성.** ① **mouth(A top)**: `A_top − overlap` = Fuse의 `split_loop_all_arcs`+notch-piece stitch **그대로 재사용**(cantilever 없이). 슬래브는 2 스트립. ② **각 뚫린 벽 opening**: edge-slot 옆-notch `splice_notch([c,d,d,c])`를 벽마다. ③ **B 벽 클립**: `clip_face_by_plane`를 뚫린 벽 평면들로 fold. 감지는 top 크로싱을 접촉 엣지별로 묶어 **뚫린 벽=엣지**, 1벽 또는 2-평행벽(`n_out.dot<0`)만 허용.
- **★ 실측으로 확정된 세 함정(계획 검토가 배선 전 차단).** ⓐ **mouth split엔 top 크로싱만** — crossings에 top(z=1)·floor(z=0.5) d가 섞여, 전체를 `split_loop_all_arcs(P)`에 주면 d가 P 루프에 오배치. `plane_side(contact)==0` 필터. ⓑ **clip 크로싱은 클립 평면 위로** — 슬래브 B-벽 바닥 엣지가 x=−0.5→1.5로 두 벽 floor 크로싱을 다 지나, x=0 클립에 x=1 크로싱 오삽입. `plane_side(keep_tri)==0` 필터. ⓒ **blind 게이트=모든 off-plane B 정점이 비-뚫린 A 면 전부 strict inside** — through-bottom·비뚫린 벽·코너를 한 조건으로 잡음(이전 "뚫린 벽 하나의 밖"은 through-bottom 놓침).
- **net.** `cut_a_slab_channel`(부피 0.9·validate·watertight·공면-인접 0·연쇄), OCCT diff `slab_channel_cut_then_cut_matches_occt`(0.9 + 채널 밖 관통 연쇄). edge-slot·Fuse 전 계열·boss/pocket 불변. 당시 인접-2벽(코너)은 out-of-scope(다음 셀 -corner-cut이 엶). **후속(비목표)**: 코너 Cut(인접벽·모서리 엣지), through-slot(바닥 관통), >2벽, Common, non-axis-aligned·4-plane.

**코너 슬롯 Cut — 인접 두 벽 파내기 (M5 셀 coplanar-contact-overhang-corner-cut).** 파내기 다음 칸. 프리즘 B가 base A의 **꼭짓점을 물고 인접 두 벽**(x=1·y=1)으로 나가는 코너 슬롯. 슬래브(반대벽)와 달리 두 개구가 **모서리 기둥에서 만나고 벽 옆-엣지까지 접함**. 구조: 각 뚫린 엣지에 top 교차 **1개**(슬래브는 2) + **base 코너 (1,1) 삼킴**(Fuse-코너처럼). n0: `vertex_on_face_plane` 거절(공면 top grazing; 직전 셀의 `a_corner_slab_cut_is_out_of_scope`가 detect None 고정하던 픽스처 → 이 셀이 성공으로 엶, 성공 `cut_a_corner_slot`).

- **★ mouth 재사용 + corner-splice.** mouth(A top L자)는 `mouth_notch_pieces`(slab에서 factor; `split_loop_all_arcs`+`stitch_arcs`) **그대로**. 각 뚫린 벽은 **corner-splice**(`splice_corner`) — 벽 loop에서 **삼킨 base 코너 정점을 `[c(top), d(floor interior), cc(모서리기둥)]` 3노드로 치환**. 두 벽이 `cc={W1,W2,B-bottom}` 공유(코너 아래 z<0.5 기둥 엣지는 남아 두 벽 봉합). 삼킨 코너는 세 면(mouth·두 벽)에서 모두 제거돼 미참조.
- **★ B 벽 클립을 기하 SH로 재작성(핵심 함정).** 슬래브의 canonical-lookup fold는 코너에서 깨진다 — B-bottom `[0.5,1.5]²`을 x=1로 먼저 클립하면 중간 정점 `(1,1.5,0.5)`(A 밖, 비-canonical)가 생기고, 이를 y=1로 다시 클립해야 `cc`가 나오므로 중간점이 필수. `clip_bwall_inside_a`: 뚫린 벽 반평면들로 순차 클립하며 **교차점을 좌표로 계산**(`segment_plane_point`, 중간점 포함), **최종 정점만** Node로 매핑(B 정점=`Node::Orig`, canonical 교차=`Node::Seam`, exact 좌표 매칭). 미매칭은 정직 거절. 슬래브는 최종이 전부 canonical이라 이 방식으로도 동작.
- **★ 감지(slab와 동일 exact blind 게이트).** 인접 두 벽(`n_out.dot==0`; 반대면 slab), 삼킨 코너가 footprint 내부(`point_in_ring2`), off-plane B 정점이 비-뚫린 A 면 전부 strict inside(through-bottom·타 벽·>2벽 배제).
- **net.** `cut_a_corner_slot`(부피 0.875·validate·watertight·공면-인접 0·연쇄), OCCT diff `corner_slot_cut_then_cut_matches_occt`(0.875 + 반대 코너 관통 연쇄). slab-cut·edge-slot·Fuse 전 계열·boss/pocket 불변. 거절 발화 `a_corner_cut_through_the_bottom_is_out_of_scope`. **후속(비목표)**: through-slot(바닥 관통), 뚫린 벽 3개 이상, Common, non-axis-aligned·4-plane.

**파내기 통합 — 하나의 일반 N-wall Cut (M5 셀 coplanar-contact-overhang-cut-general).** 파내기 스레드의 정점. edge-slot·slab·corner 세 갈래(`detect_overhang_cut_contact`+`detect_overhang_corner_cut`, 두 result)를 **하나**(`detect_overhang_cut_general`+`overhang_cut_general_result`)로 흡수하고, 덤으로 **3벽(L·U자)**를 연다. 붙이기가 single→corner→multi를 하나로 합친 것의 파내기 판. **M6에서 곡선 축퇴가 붙을 때 가지-폭발을 막는 토대** — 한 곳에 곡선 analog만 얹으면 됨. n0: 3벽 L 픽스처(A `[0,1]³`, B `[-0.5,1.5]×[0.5,1.5]×[0.5,1]`)가 `vertex_on_face_plane` 거절(공면 top grazing)·두 기존 detect None 확인, 목표 0.75.

- **★ 결과 = 각 A-면을 그 면 위 cap만큼 notch + B 벽(클립·flip).** 각 뚫린 A-벽의 notch를 **top 엣지에 삼킨 base 코너 수(0/1/2)**로 분기: **0**=middle-splice(`splice_notch [c,d,d,c]`; slab), **1**=corner-splice(`splice_corner [c,d,cc]`; corner), **2**=벽이 코너기둥까지 짧아진 사각형(신규 `substitute_corners` — 삼킨 두 코너 정점을 각자 `cc`로 치환, top 엣지가 자동으로 `cc_i→cc_j` B-bottom trace가 됨). mouth·`splice_notch`·`splice_corner`·`clip_bwall_inside_a`(기하 SH)·`assemble_fuse_cut` 전부 재사용.
- **★ 뚫린 벽 = A-top 엣지가 B footprint와 겹치는 엣지(크로싱만으로 판정 금지).** B가 엣지를 통째로 덮으면(양끝 삼킴) 크로싱 0인데도 뚫린다(L의 y=1 벽) — 이걸 놓치면 그 벽이 원본대로 남아 결과가 틀린다. 각 엣지의 삼킨 끝점 수 `n`으로 분류하고 **불변식 게이트 `(엣지 위 크로싱 수) == 2 − n`**(n=0→2, n=1→1, n=2→0)로 볼록-박스 정합성 확인, 안 맞으면 `None`. 코너기둥 `cc = three_planes(W_i, W_j, B-bottom)`는 삼킨 코너의 인접 두 뚫린 벽이 공유(단순 loop이라 코너당 정확히 두 벽).
- **★ 안전 순서로 통합(회귀 0).** 먼저 일반 경로를 dispatch에 **먼저** 걸어 기존 `cut_an_edge_slot`·`cut_a_slab_channel`·`cut_a_corner_slot`이 전부 일반 경로로 그린임을 확인한 뒤, 그때 두 detect + 두 result + 죽은 헬퍼(`clip_face_by_plane`·`crossing_on_segment`·`shared_vertex`)를 제거. edge-slot=1벽 0코너, slab=2 반대벽 0코너, corner=2 인접벽 1코너 each의 특수케이스.
- **net.** `cut_an_l_step`(부피 0.75·validate·watertight·공면-인접 0·연쇄), OCCT diff `l_step_cut_then_cut_matches_occt`(0.75 + 전면 스트립 관통 연쇄). edge-slot·slab·corner이 일반 경로로 재현(회귀 가드)·Fuse 전 계열·boss/pocket 불변. 거절 발화 `an_edge_slot_through_the_bottom_is_out_of_scope`·`a_corner_cut_through_the_bottom_is_out_of_scope`(일반 detect None). **후속(비목표)**: through-slot(바닥 관통), full-face(base⊂커터·교차 0), 한 코너기둥 3면 이상 공유, Common, non-axis-aligned·4-plane 동시성.

**교집합 오버행 — 연산 3종 코너 완성 (M5 셀 coplanar-contact-overhang-common).** Fuse boss·Cut overhang은 이미 통합됐지만 **Common만** coplanar-contact 계열에서 미개방(7개 overhang 셀 후속 목록 전부에 등장)이었다. 파내기와 **동일 픽스처**에서 `Cut`은 A\R을, `Common`은 **볼록 겹침 R=A∩B**를 낸다(볼록∩볼록=볼록이라 결과는 항상 단일 볼록 상자). 감지는 kind-무관 기하라 `detect_overhang_cut_general`을 **통째 재사용**하고, `detect_overhang_common`이 **콘택트 평면 top 크로싱 == 2** 게이트만 얹음(⟺ inside-P arc 1개 ⟺ 1회 stitch로 조립; **뚫린 벽 수 무관** — edge 1벽·corner 2벽·L 3벽 전부 커버, 4+ slab은 `vertex_on_face_plane` 표준 fall-through). n0: edge·corner 픽스처가 그 태그로 거절 확인.

- **★ 결과 = 겹침-top + A 뚫린 벽∩B + B 벽∩A, 전부 flip 없음(R 볼록).** 겹침-top은 mouth의 형제 `overlap_top_piece` — mouth가 **outside**-P arc를 골랐다면 겹침은 **inside**-P arc를 골라 inside-Q arc와 `stitch_arcs`(P∩Q). contained 케이스는 상위에서 걸러지므로 `!outside` arc가 정확히 1개 보장. B 벽∩A는 `clip_bwall_inside_a`(기하 SH, A 뚫린 벽 평면 fold) — Cut의 B-벽 클립과 동일, `flip:false`만 차이.
- **★ A 벽∩B는 clip이 아니라 크로싱에서 정확 구성(구현 중 발견한 함정).** 계획은 A 벽도 `clip_bwall_inside_a`로 B 평면들에 fold하려 했으나, **B의 평행 두 평면(x=0.3·x=0.7)을 연쇄 clip하면 두 번째 교차 좌표가 중간점에서 재계산돼 반올림**된다(`(0.3,1,0.5)`→`(0.30000000000000004,…)` → 저장 크로싱 매칭 실패 → 거절). R이 볼록이라 A 벽∩B의 모든 꼭짓점은 **이미 저장된 크로싱이거나 공유 원본**이므로, `common_wall_face`가 벽의 `WallKind`(Middle `[c,c,d,d]`/Corner `[c,corner,cc,d]`/Shorten `[v,v,cc,cc]`)에서 노드를 **직접** 뽑고 `orient_to_wall`(원본 벽 loop과 signed-area 부호 맞춤)로 winding만 고친다. Cut이 벽을 splice로 정확히 다루는 것과 같은 정신 — clip의 float 손실을 피함. (B 벽은 A의 **수직/단일** 평면으로만 clip돼 exact, Cut과 동일하게 유지.)
- **net.** `common_an_edge_overhang`(0.1)·`common_a_corner_overhang`(0.125, cc 1개)·`common_an_l_step`(0.25, cc 2개 공유; 2-크로싱 게이트가 뚫린 벽 수 무관 증명) — 각 부피·validate·watertight·공면-인접 0·**`is_convex`**(면-방향 오류를 부피 너머로 잡음)·cavity 없음. OCCT diff `edge_overhang_common_matches_occt`·`corner_overhang_common_matches_occt`. Fuse·Cut 전 계열·기존 Common(횡단)·boss/pocket 불변. 거절 `common_a_slab_overhang_is_out_of_scope`(게이트 None + 표준 거절). **후속(비목표)**: spanning-slab Common(k≥2 lens 겹침-top 체이닝), contained Common(B⊂A), through, non-axis-aligned·4-plane, 홀 있는 피연산자.

**M5-d2 거절 가드 — 실측으로 확정된 세 가지 (거절 태그 훅 도입 후).** `assert_eq!(boolean(..), Err(Unsupported))`는 의도와 다른 가드가 발화해도 통과하므로, ops의 모든 `Unsupported` 생성 지점에 `reject(tag)`를 붙이고 테스트가 `assert_rejects(.., tag::…)`로 **어느 가드가 발화했는지**까지 단언하게 했다(`#[cfg(test)]` thread_local, 릴리스 영향 0). 이 계측이 드러낸 것:
- **`strict` seam-arc 가드는 과잉 거절이었다 — 셀 3e-1에서 해소.** 자기교차가 아닌 **단순 계단형 호**(reflex turn 1회 이상)까지 접힌 호와 구별 못 해 거절했다. 그것은 **정렬 인공물**을 막는 가드였고, 셀 3d가 정렬을 없앴다.

  **대체가 아니라 삭제다.** 호는 `∂Y ∩ P ∩ f`의 한 성분, 즉 1-manifold이므로 단순 곡선이고, 양 끝만 `∂f`에 닿는다(내부 노드가 `∂f`에 닿으면 `segment_crosses_face`의 graze 검사가 `CONTACT_DEGENERATE`로 거절한다). 따라서 `유지된 ∂f 구간 + 호`는 회전 방향과 무관하게 **단순 다각형**이다. 그걸 다시 검사하려면 축 투영 + `orient2d`, 즉 이 서브유닛이 3b부터 피해온 2D 기계가 필요하다 — 증명된 것을 검사하려고 들이지 않는다.

  실측: `l_and_popup_box`의 박스 바닥면에서 재구성된 고리는 `(0.5,0.5)→(0.5,1.5)→(1,1.5)→(1,1)→(2,1)→(2,0.5)`, 겹침 발자국(넓이 1.0), `(1,1)`이 반사 꼭짓점. **`strict`는 옳은 면을 거절하고 있었다.** `validate`는 자기교차 면을 통과시키므로(매니폴드, 오일러 성립) 부피와 OCCT diff(`folded_arc_{cut,fuse}_matches_occt`)만이 게이트다.

  **볼록 경로(`reconstruct_face`)는 여전히 투영 정렬을 쓰지만 `strict`가 필요 없다.** 볼록 영역 둘의 경계 교선은 볼록이므로 chord 방향으로 단조다 — 가드가 아니라 기하학적 사실이고, `edge_seam`의 1:1과 같은 종류의 진술이다.
- **`tunnel` 가드의 진짜 역할** (**3e-3에서 가드의 형태가 됐다 — 아래**). 주석이 주장하던 out→in→out 터널은 두 면을 뚫으므로 당시 `pierced_multi`가 먼저 삼켰다. `tunnel`의 **유일한** 발화 경로는 cavity 비대칭이다: `point_in_solid`는 outer+cavity를, `pierced_face`·`collect_planes`·`solid_vertices`는 **outer shell만** 센다. cavity 안에서 끝나는 엣지는 양 끝이 `Outside`로 분류되면서 outer 경계를 1회 관통한다. 살아있는 백스톱이며, 위 hollow 오답을 비볼록 경로에서 막고 있는 것이 바로 이 가드다.
- **`fourplane` 가드는 미검증.** 유리수 좌표로 4-평면 동시성을 인위적으로 맞추는 적대적 픽스처만 존재하고 자연스러운 형상이 없다. 서브유닛 3에서 정상 경로로 확인 가능. (`pokehole`·`tunnel`은 발화 테스트를 갖췄다 — `pokehole`은 3f-1·3f-2에서 두 번 **축소**됐고 이제 `two_loops_on_one_face_are_unsupported`가 남은 범위(고리 ≥2, 호 옆의 고리)를 고정한다; `multichord`는 셀 3e-2에서 비볼록 경로가 은퇴시켜 **미발화 백스톱**이 됐고, `strictarc`는 셀 3e-1에서 은퇴했다. 한 설계 검토는 `multichord`가 "구조적으로 도달 불가능"하다고 결론냈으나 반증됐다 — 한 면의 4-transition은 **서로 다른 두 straddle 엣지**가 각각 1회씩 관통해 만든다.)

**명시 공유 — pad/pocket cap이 면 Surface handle 공유 + Handle 기반 공면-접촉 인식 (오버홀 이식 #3).** 오버홀 §5의 surface 공유를 커널에 켠다. `build_prism`이 밑면 cap을 fresh push하던 것을, `extrude_and_boolean`이 넘긴 **접한 면의 `Surface` handle을 재사용**(월드-plane `extrude`·테스트 caller는 `None`). 그러면 flush 접촉이 **shared-handle 공면쌍**이 되고, 세 감지기(`detect_coincident_interface`·`detect_contained_contact`·`detect_pocket_contact`)가 `shares_or_coplanar`(= `pa.surf == pb.surf || planes_coplanar`)로 **O(1) Handle 동일성**으로 접촉을 인식한다.

- **★ orientation 재조정 — surf와 orient는 함께 이동.** cap의 요구 바깥법선은 `−sweep`인데 공유 surface의 stored normal은 다를 수 있어, `orient = Forward if n_h·(−normal)>0 else Reversed`로 맞춰 **materialized = −sweep 유지**(base_loop 불변). 재구성부(`assemble_fuse_cut`)가 결과 면 orientation을 `planes[idx].orient`(raw)로 정하되 **같은 idx의 surf와 쌍으로** 복사·flip하므로, 이 쌍-이동이 materialized를 flip/non-flip 양쪽에서 보존한다 → **기하 Case A**. (raw orient를 surf와 분리해 쓰는 코드는 없음 — grep 확인.)
- **★ M5에선 handle 분기가 redundant(정직).** 축정렬에서 `surf==` ⟹ `planes_coplanar`(같은 handle=같은 평면)라 OR의 handle 쪽은 결과를 못 바꾼다(순수 redundant). 기하 fallback을 빼면 공유 안 된 공면-접촉(독립 솔리드)이 깨지므로 **fallback 필수**. 즉 handle 검출이 기하와 **관찰상 divergence하는 건 회전에서만** — 이 셀은 **Case A + 메커니즘 실증**(구조: cap이 handle 공유 + `shares_or_coplanar` handle 분기가 비공면 쌍에도 참). 행동 페이오프는 회전 배선(#4) 후.
- **정직한 빚·후속.** 이 셀은 **인식만**, 6개 기하 감지기 **은퇴 0**(§9: 공유가 커버함을 실증한 뒤에만). overhang 감지기 3곳은 기하 게이트 유지. 엣지·정점 공유(§5 2·3층)와 회전 배선은 후속. n0 실측: pad/pocket/boss·coplanar-contact 전 픽스처 부피·validate·OCCT diff 불변, 공유-면(F·cap이 한 Surface) validate/STEP/tess 견딤, 퇴화 ThreePlane 없음(공면-접촉은 특수경로).

**다중 솔리드 불리언 — sever가 `Vec<Solid>`를 낸다 (오버홀 단계 0.4).** `DISCONNECTED_RESULT`(몸체를 두 동강 내는 `Cut(rod,L)`을 거절하던 가드)를 은퇴시키고, **바깥향(재료-포함) 성분이 여럿이면 각각을 솔리드로 반환**한다(OCCT n0 §8이 정한 정책의 구현). `boolean`/`general_boolean`/`assemble_fuse_cut` + 6개 접촉 사이드카가 `Result<Vec<Handle<Solid>>>`를 반환. seam-free 경로(`contained_result`)는 최대 1개라 `vec![…]`로 감쌈. 회전과 무관·독립이라 오버홀 arc보다 앞서 main-병합 가능한 초기 셀.
- **★ 재구성은 이미 성분을 안다 — 가드만 걷었다.** `assemble_fuse_cut`은 진작 결과를 `face_components`(공유 `Node` union-find)로 나누고 `is_shell_outward`((5d)#5의 exact 극단-정점 부호)로 outward/inward를 판별했다. 옛 로직은 "outward 정확히 1개면 수용, 아니면 `DISCONNECTED_RESULT`". 새 꼬리는 `match positives.len()`: **1** ⇒ outer + 나머지 cavity(기존 그대로), **≥2 & inward 없음** ⇒ 각 성분이 제 cavity-free 솔리드(sever), **≥2 & inward 있음** ⇒ `SEVERED_WITH_CAVITY`, **0** ⇒ `NO_OUTWARD_SHELL`.
- **★ 반환 Vec 순서 = 기하-기반 canonical(replay-안정).** 다운스트림 op가 Vec 인덱스/핸들로 솔리드를 고르므로 순서가 사용자 가시적이다. 각 severed 솔리드를 `comp_key`(성분 정점 좌표를 사전순 정렬)로 정렬 — **최종 기하에만 의존**하므로 replay가 같은 기하를 재생성해 구성상 안정(재구성 내부 순서라는 미증명 불변식에 의존 안 함). lex-min 정점 단독은 tie-free가 아니다(DNA 원칙 2: 서로 다른 정점이 같은 좌표 가능) — 전체-좌표-리스트로 동점을 깬다(disjoint 성분은 좌표집합이 반드시 달라 total order). n0: 원시 재구성 순서도 이미 replay 2회 동일이나, 정렬이 이를 기하로 못박음.
- **★ `SEVERED_WITH_CAVITY`는 발화 테스트와 함께 태어난다.** ≥2 outward + cavity 생존은 "어느 outer가 어느 cavity를 소유하는가"에 shell-scoped point-in-shell 판정(=`point_in_solid` winding 코어를 단일 셸로 re-scope)이 필요하나 아직 없어 정직 거절·defer. **도달 가능**(오프셋 void hollow box를 void를 피해 가르는 slab; n0 실측: outward 2 + inward 1) → firing 테스트 `severed_with_cavity_is_rejected`. `NO_OUTWARD_SHELL`(positives==0)은 재료 경계 없는 결과라 **기하적으로 도달 불가·방어 백스톱**(발화 테스트 없음, `fourplane`류).
- **★ pad/pocket는 다중 솔리드를 수용한다(원자성).** `extrude_and_boolean`의 `boolean`은 inspect 전 이미 커밋(입력 은퇴 후 결과 push)하므로, append-only에 롤백이 없어 "커밋 후 거절"은 반쯤 변형된 모델을 Err로 남겨 `model=replay`(DNA #3)를 깬다. 그래서 거절 대신 수용: 반환 솔리드 전부에서 cap을 찾아 cap 보유 primary를 반환하고 나머지 조각은 live로 둔다. pad=Fuse는 결코 sever 안 하고, blind pocket sever는 병리적·드묾 — 예전 `DISCONNECTED_RESULT` 거절이 이제 유효 다중 솔리드 성공으로 바뀌나 비원자적 거절보다 우월.
- **net.** 옛 거절 3테스트가 성공으로 뒤집힘: `cut_rod_by_l_severs_into_two`(2솔리드·validate 클린·부피 합 0.06), `a_convex_cut_severs_its_operand`(2×1.0), `a_slab_splits_a_hollow_box_into_two`(void 열림·cavity 0·합 24.4). OCCT diff `sever_cut_matches_occt_compound`(nacre 2솔리드 vs OCCT COMPOUND-of-2, 부피·면적 합 일치). churn: ~199 단일-솔리드 테스트 사이트를 테스트-전용 `boolean_one` shim으로(먼저 순수 rename n1a, 다음 타입 flip n1b — Rust 컴파일 원자성). **후속(정직한 빚)**: seam-free 다중(disjoint Fuse + edge-touch Fuse, 현 `EmptyResult` 경로)·OCCT shared-edge 토폴로지 재현(§8: §5 명시 공유 + per-solid validate)·multi-outer-with-cavity containment·`Cut(A,A)` empty-Vec.

**`Transform` 연산 — 강체 이동 (오버홀 단계 1a).** 회전 arc의 첫 실제 커널 작업. 강체 회전·이동(isometry)은 곡면을 만드는 revolve가 아니라 솔리드를 통째로 옮기는 것 — 평면 그대로(M6 아님). 1a는 그 **기계적 뼈대**를 exact 이동으로 먼저 세운다(회전은 1b). **이 셀이 커널의 nacre-scalar(오버홀 foundation #1) 최초 소비.**
- **★ 표현 계층은 이미 회전 준비 완료(조사 확정).** `Plane`·topo·validate·tess·STEP 전부 방향 무관(프레임 기반). 축정렬 가정은 `nacre-ops` boolean 단축(is_shell_outward 등, "TIP의 일")에만. ⇒ 이동은 순수 additive(boolean 무영향).
- **★ `Isometry`(nacre-scalar) — 정의·실현 분리.** op가 `Operation::Transform { solid, isometry }`를 싣고, `Isometry`는 유리수 데이터(`translate:[Rat;3]`) + `offset_f64()` 실현만 — **nacre-math 미참조**(frame.rs `Pt2`가 `Point2` 아닌 튜플을 쓰는 선례; nacre-scalar·nacre-math는 형제). `Point3`/`Plane` 적용은 nacre-ops. 1b가 rotation을 struct에 추가해 enum·match churn 0.
- **★ `transform_solid` 다중 패스 복제(재사용 핵심).** shell(outer+cavity)→face→loop(outer+inner 구멍)→edge→vertex를 **결정적 순회**로 복제, old→new dedup 맵(surface·curve·vertex·edge). 패스 순서 = surface→curve→vertex→edge→face→shell→solid(정점 `ThreePlane` remap이 완성된 surface 맵을 읽도록). replay 결정성(DNA 3)·#3 surface 공유·구멍·cavity 보존.
- **★ Discovered 정점 정의 보존.** boolean seam 정점(`Discovered{tol, ThreePlane}`)은 **Origin 종류 유지** + `ThreePlane` surface handle을 이동 surface로 remap(tol 불변). Constructed로 강등하면 정의를 잃고 exactness를 거짓 주장 — 금지. surface가 face surface에 없으면 불변식 위반이라 fail-loud. 이동은 §⑦ 회전 트리에 안 담김(op-log가 정의)·tol에 안 기여; 단계 1은 tol을 기록만·판정 안 함이라 이동 실현의 ~ulp는 무시(엄밀 bump는 단계 2 tol 관리).
- **★ 평면/곡선 = 한 번 변환 + dedup.** `Plane::translated`(origin+offset, **법선·raw 불변 → raw exactness 보존**)·`Line`/`Circle`/`Cylinder` translate. 재구성 아님 → 공유 평면 두 면이 동일 평면(f64 divergence 0). face orientation 불변(이동은 법선 안 돌림). 1b가 회전 변환으로 확장.
- **net.** `transform_translate_cuboid`(부피·면적 불변·bbox=offset 이동·validate/tess/STEP 클린·supersede), `transform_translate_preserves_discovered_definition`(Cut 결과 이동, Discovered 6개 보존·validate 클린·부피 불변), `transform_is_deterministic`(같은 ops→같은 handle), `transform_op_applies`(apply 경로), OCCT diff `translated_solid_cut_matches_occt`(이동 큐브 Cut이 OCCT 일치). boolean 거절 없음(이동은 유효). **후속**: 1b 회전(`Store<Rotation>`·`Origin` 회전 변이·일반각 무리수·90°계열 exact·boolean 거절 가드·리스크 #2/#3), 단계 2 TIP 이식(exact3d), 단계 3 회전 boolean.

**`Transform` 회전 — 축정렬 강체 회전 (오버홀 단계 1b).** 1a 이동 뼈대 위에 **회전**을 얹는다. 솔리드를 좌표축(X/Y/Z) 둘레·유리수 축점을 지나 각도만큼 통째로 기울임 — revolve 아님·곡면 없음(M6 아님). 회전 arc의 하드 셀(리스크 #2/#3): 회전 정의 기록·90°계열 exact·불리언 명시 거절(TIP 전이라 아직 sound 아님). **단계 1 = 표현만·판정 없음**(tol은 기록만·읽지 않음 → 단계 2가 판정). **축 = 축정렬 X/Y/Z 확정** — exact3d가 tol soundness를 축정렬로만 검증(임의 유리수 축 Rodrigues는 "later extension"·미검증)이라 단계 2 판정과 정합. 비-90° 각이라야 기울어진 기하·Rotated origin·거절 경로가 의미 있다(90°계열 축정렬은 여전히 축정렬·exact이라 이동과 동급).
- **★ 회전 정의 = 정의만 저장·tol은 저장 안 함(§⑦).** 회전 tol은 적용점 의존 → 판정 시 계산(단계 2). ⇒ `Origin::Rotated { base:Handle<Vertex>, rotation:Handle<Rotation> }`는 tol 슬롯 없음(현 `Discovered`와 다름). base = 회전 전 정점(arena 잔존; 그 Origin이 회전 전 정의). Discovered→Rotated 강등 아님 — 별도 변이라 거절 가드가 boolean-Discovered와 명확히 구분.
- **★ `Isometry`에 `rotate: Option<Rotation>` 추가(1a 예측 실현).** `Rotation { axis:Axis{X,Y,Z}, point:[Rat;3], angle:Angle }` + `is_exact()`(회전 없음 또는 90°계열 = `try_exact_cos_sin` Some)·`apply_point`/`apply_dir`(f64, math 독립 — nacre-scalar가 nacre-math 미참조 유지). 적용 순서 = 축점 둘레 회전 후 이동. enum churn 0(1a가 예측한 대로 struct 확장뿐).
- **★ 정점만 Rotated 마킹(엣지/곡선은 파생 캐시).** 술어가 판정하는 건 정점(ThreePlane→3평면); 엣지 `Line`은 두 회전 정점을 지나는 파생 캐시라 정점이 회전 정의를 지니면 따라온다. ⇒ 비-exact 회전 시 정점만 `Origin::Rotated`, 엣지 Origin은 1a처럼 보존. 거절 가드는 정점만 스캔해 충분.
- **★ `Store<Rotation>` 공유 노드(§⑦ forest).** 비-exact 회전 1개당 노드 1개를 push, 그 솔리드의 전 회전 정점이 공유(per-vertex Angle 복사 금지). 단일 회전은 `parent:None`; 체인 번들링은 후속 셀(단일 회전은 체인 없어 트리거 없음). **첫 topo→nacre-scalar 의존.**
- **★ 기하는 성분 회전 + 기존 생성자 재구성(geom 무변경).** 회전 있으면 nacre-ops가 평면/곡선 성분(origin/normal/direction)을 `apply_point`/`apply_dir`로 회전시켜 `Plane::from_point_normal` 등으로 재구성(raw exactness는 회전으로 어차피 상실·예상됨). 회전 없으면 1a `translated`(raw 보존). ⇒ geom에 rotated 메서드 불요. `is_exact()`가 Origin 결정 — exact면 1a처럼 보존(Constructed→Constructed·Discovered→ThreePlane remap), 비-exact면 Rotated. 다중패스·dedup·구멍·cavity·결정성은 1a 그대로.
- **★ 불리언 거절 가드(O(1)).** `boolean`(live 체크 뒤) `if solid_is_rotated(a)||solid_is_rotated(b) → reject(ROTATED_UNSUPPORTED)`. Transform은 솔리드 전체를 균일 회전(all-or-nothing)·boolean이 회전 입력 거절 → 혼합 회전 솔리드 없음 ⇒ 정점 하나만 확인해 충분(M5 hot path 무부하). 90°계열은 Rotated 안 만드니(exact·Constructed) boolean 허용(축정렬 유지).
- **★ validate `Origin::Rotated`.** incidence 검사(정점 f64가 회전 곡선/평면 위, 잔차 미미)만 하고 정의-검사는 skip(EPS_CONSTRUCTED 머신-scale). Discovered ThreePlane 정의-검사는 auto-skip. n0 실측: 30° 회전 cuboid·회전 Cut 결과 모두 false-fail 0.
- **net.** `transform_rotate_cuboid_tilts_and_blocks_boolean`(30°Z 유리수 축점: 부피·면적 불변·validate/tess/STEP·모서리 정확 회전 이미지·Rotated origin·boolean 양방향 거절), `transform_rotate_90_is_exact_and_allows_boolean`(90°: Constructed 유지·부피 24 정확·후속 Cut 허용·검증), `transform_rotate_boolean_result_keeps_discovered_base`(Cut 결과 회전 → Rotated{base=Discovered}·validate·부피 불변), `transform_rotate_is_deterministic`·`transform_rotate_op_applies`, OCCT diff `rotated_solid_props_match_occt`(회전 cuboid+회전 Cut 단일-솔리드 export → OCCT volume·area 일치). **defer(후속 셀)**: 번들링(연속 동일축 = 누적각 1회 실현·H4 증폭 방지)·강등(유리수 오버플로)·임의 유리수 축(Rodrigues) — 전부 단계 2 판정 전까지만 필요. 단계 2 TIP 이식(exact3d `Pt3`/`orient3d_judge`), 단계 3 회전 boolean(TIP 소비·거절 은퇴·main 병합).

**재회전 forest — 체인·base=루트 (오버홀 단계 1c).** 1b는 **단일 회전**만 다뤘다. 회전된 솔리드를 다시 회전하면 1b는 조용히 **반쯤-깨진 forest**를 냈다: `rot_node.parent` 항상 None(평평)·`Origin::Rotated.base`가 루트 아닌 **직전 회전 정점**·exact 회전을 회전된 솔리드에 주면 `is_exact→remap`이 노드를 **빠뜨림**(좌표엔 있는데 forest엔 없음). 이 forest는 **단계 2(TIP 판정)가 소비**(회전이력 순회·tol 누적)하므로 지금 **완전하고 sound하게** 세운다. **★ 이 셀은 번들링(최적화)이 아니라 forest 정확성만 고친다 — always-chain**(번들은 stage 2 직전 별도 셀; 그 이득 tol 타이트닝은 stage 1에서 측정 불가·인접-가드 soundness도 그때 검증 가능). always-chain은 exact3d가 sound 확인(un-bundled ~30.8× 느슨할 뿐)·**가드 불요**(이동 개입 여부와 무관하게 각 회전이 자기 노드·자기 da 항 → 항상 sound).
- **★ 관찰 가능한 동작 변화 없음(순수 stage-2 인프라).** 재회전·exact-after-inexact 전·후로 Rotated·boolean 거절·같은 부피/좌표/validate — **내부 forest 구조만** 달라짐(끊긴→완전·base=루트). 방법론("인프라 먼저")·exact3d 검증 모양(`Pt3` base+chain)에 정합. ⇒ 테스트는 필연적으로 **white-box**(`forest_probe`: `Model.rotations` 노드 수·parent 사슬·base 루트 도달성 직접 assert).
- **★ A/B0/B1 결정 트리(`transform_solid`).** **A**(순수 이동)=1a/1b remap. **B0**(fresh 회전, 입력 미회전)=1b 보존: exact→remap(Constructed/Discovered 유지)·inexact→노드 parent None. **B1**(재회전, 입력 이미 회전)=**항상 체인**: 노드 `{axis,point,angle, parent: input_leaf}`·**exact여도 기록**(조상 inexact→합성 inexact·forest 완전; exact 노드는 stage 2 tol 0이나 조상 tol을 `|R|·old`로 transport하므로 필수). `solid_rotation` 헬퍼가 입력 공유 leaf를 peek(대표 첫 정점)·**균일 회전 불변식 `debug_assert`**(fail-loud).
- **★ base = 항상 루트(한 홉).** `Origin::Rotated.base`는 비-Rotated exact 루트(Constructed/Discovered)를 가리키고 `rotation` 부모 사슬이 그 루트로부터의 완전한 회전 이력(exact3d `Pt3.base+chain` 미러). 입력이 Rotated면 새 base=그 base(한 홉), 아니면 vh. **fresh frame 없음**(base를 Rotated 정점에 놓으면 stage 2 재귀가 회전 이중적용/이동손실). 이동은 forest 제외(§⑦)라 회전→이동→회전도 체인 후 base+forest가 이동 미재현이나 **sound**(tol은 노드 축 + 현재 캐시 좌표 모멘트암, base 위치 무관).
- **net.** `rerotate_same_axis_chains`(30Z→20Z: 노드 2개 [Z,Z]·base=루트·불변·거절), `rerotate_different_axis_chains`([Z,X]), `rerotate_exact_after_inexact_records_node`(30Z→90X: exact 90X도 노드 기록·여전히 Rotated), `rerotate_across_translation_chains`(회전→이동→회전: 가드 없이 sound), `rerotate_deep_chain`([Z,X,Y]), `fresh_rotation_of_constructed_unchanged`(B0 회귀 가드), `rerotate_is_deterministic`, OCCT diff `rerotated_solid_props_match_occt`. **★ pre-existing 1b 버그(1c n0 발견 → 셀 1d에서 해소)**: fresh exact 회전을 **Discovered(tol=0)** 솔리드에 주면 `remap`이 tol 0 유지·f64 회전 반올림(~8e-17)이 초과 → `VertexOffSurface`(축정렬 boolean이 tol 0 정점을 냄). B1은 전부 Rotated→EPS_CONSTRUCTED라 클린이라 이 1c 셀 밖이었고, **1d(90°계열 정확 실현)가 근본 해소**. **defer**: 번들링(같은-축 인접 누적각·인접-가드)·강등·Rodrigues.

**90°계열 회전 정확 실현 (오버홀 단계 1d).** 1c n0가 드러낸 pre-existing 1b 버그를 근본 해소. `Isometry::apply_point`/`apply_dir`이 90°계열(0/90/180/270°)에서도 f64 `cos()`/`sin()`을 써서 `cos(90°)≈6.12e-17`(정확히 0 아님) spurious cross-term을 남겼다. 축정렬 boolean이 낸 **`Discovered`(tol=0)** seam 정점을 exact 90° 회전하면 그 ~8e-17 잔차가 tol 0을 초과 → validate `VertexOffSurface`(무효 기하). 이는 §⑦ 원칙("90°계열은 tol 0·축정렬-급 exact")과 코드의 불일치 — exactness 개념은 `is_exact()`·rot_node 결정·boolean 가드에 이미 있는데 **좌표 실현만** f64 spurious를 넣었다.
- **★ `Angle::cos_sin_f64` — f64 실현의 진실 단일화.** `cos()`/`sin()`/`try_exact_cos_sin()` 옆에 신규: quadrantal이면 `try_exact_cos_sin`의 `Rat` 0/±1을 `to_f64()`(정확 0.0/±1.0), 아니면 f64 `cos`/`sin`. `apply_point`/`apply_dir`이 이걸 호출(중복 제거 + 미래 실현 경로가 버그 재도입 못 하도록 한 곳에). exact3d `Pt3::rotate`가 검증한 미러(그쪽은 tol도 다뤄 인라인). `frame.rs::rotated_about_origin`(sign-judgment vehicle)은 무변경.
- **★ 잔차 0 보존(rigid 포함).** exact quadrantal은 순수 swap+sign이라 정점·평면이 동일하게 perturb → vertex−plane dot-product 잔차가 f64 덧셈 **교환법칙**으로 비트-보존(0→0). rigid(회전+이동)는 offset이 vertex−plane_origin에서 상쇄되어 잔차 미도입. ⇒ tol-0 Discovered 정점의 exact 90° 회전이 validate 클린.
- **net.** `nacre-scalar`: `cos_sin_f64_is_exact_for_quadrantal`(0/90/180/270° 정확·비-90°는 cos/sin 동일), `apply_point_exact_for_quadrantal`((3,5)→90°Z→(-5,3) 비트-정확·non-origin 피벗·비-90° 불변). `nacre-ops`: `boolean_result_rotated_90_validates`(Cut 결과 90°·rigid(90°+이동) 둘 다 validate 클린·부피 불변 — 버그 해소·offset-상쇄 가드), `rotate_90_lands_vertices_exactly`(정점 정확 축정렬 격자). 관찰 개선: boolean→90° 회전이 유효(전엔 18 errors)·90° 정점 정확 격자. **n2(OCCT) 없음**: 변화가 좌표 6e-17 정밀도라 OCCT(~6 유효숫자)가 못 봄·회전-불변 props는 1b n2 커버(in-crate 비트-정확 검증이 더 강함). stage 2의 "90°계열 tol 0 정확 실현" 기저경로 확립.

**Pt3 방향별 tol 계산 이식 (오버홀 단계 2a).** 스테이지 2(TIP 코어)의 첫 조각 — 스테이지 1이 세운 tol>0 원천(회전 forest)을 **sound하게 계산**하는 primitive. exact3d가 검증한 `Pt3`(정의+누적 tol)를 프로덕션 `nacre-scalar::frame3`로 이식(2D `frame::Pt2`의 3D 아날로그). **"층만"** — 판정 술어(`orient3d_judge`=2a-ii)·Model forest 순회(nacre-tip=2b)·boolean 무배선. 순수 수치·커널 무의존.
- **★ `Pt3 { base:[Rat;3], chain, coord:[f64;3], tol:[f64;3] }`** — §TIP⑤ 방향별 xyz tol 벡터(자료구조 확정=`[f64;3]`). `at`(base→f64 반올림 tol)·**`at_with_tol`**(Discovered seam tol seed — 2b가 주입, API churn 방지)·`rotate_about`(tol `|R|·old+mix`)·`hp_coord`(astro-float ground truth). `Axis::plane`·`cos_sin_at`·`try_exact_cos_sin` 및 frame.rs 공유(`DA_F64`·`rat_to_big`·`bf_mag`, `pub(crate)` 승격) 재사용 — 복붙 아님.
- **★ 임의 피벗(exact3d에서 유일한 비-verbatim 확장·먼저 검증).** exact3d는 원점 축만 검증했으나 커널 `Rotation`은 유리수 피벗을 가짐 → `rotate_about(axis, angle, point)`가 피벗-상대 실현 + `piv` tol 항(`(|ci|+|cj|+|px|+|py|)·DA`, 원점이면 정확히 0이라 H-d/H-e 비트-불변). **de-risk 규율대로 exact3d를 먼저 확장·검증**(H-f: 5000 혼합-피벗 체인 위반0·tightness 0.500 — 원점과 동급), 그다음 이식. 부차 발견: tol-0 exact 혼합 코퍼스는 512bit GT 실현 노이즈를 노출 → 1e-100 floor로 분리.
- **net.** `nacre-scalar::frame3`: `tol_bounds_error_over_random_chains`(랜덤 피벗·exact/inexact 체인 tol이 실제 f64 오차 상한·GT=512·proptest형), `quadrantal_origin_chain_is_tol_zero`(90°계열 origin tol 0·좌표 정확), `at_seeds_base_rounding_tol`, `seeded_tol_transports_through_rotation`(seed tol이 `|R|`로 전파), `bundling_is_tighter_than_incremental`(H-e). exact3d 무회귀(원점 피벗이 기존 비트-재현). **n2(OCCT) 없음**(순수 수치·boolean 미배선). §TIP⑤ 확정·⑦-캐시 defer. 후속: 2a-ii `orient3d_judge`(②·H-a)·2b nacre-tip forest 브리지·2c 간접 술어.

**orient3d_judge 이식 — 직접 3D 방향 판정 (오버홀 단계 2a-ii).** 2a의 `Pt3.tol`을 **소비하는 판정 술어**. 네 점 f64 행렬식 **필터**(`|det| > det3_bound`면 부호 확정·대부분·~300ns), 애매하면 **astro-float 상승**(정의에서 고정밀 재실현·~1.2ms), 정밀도 floor(`16·scale³·2^-JUDGE_PREC`) 이하 → **declare-0**(`Orient::Zero`·§6 사용자 질문). frame.rs 2D `orient2d_judge`의 3D 아날로그·같은 사다리·`is_positive`(astro-float#44 우회).
- **★ `det3_bound` = §TIP② 유일 남은 수학(직접분)**: 행렬식=6개 signed triple-product(cofactor). 각 edge 성분 tol=`tol_a+tol_d`(apex tol이 뺄셈으로 전파). bound=Σ`prod_err`(구간 곱 반경 `Π(|v|+τ)−Π|v|`)+`16ε·mag`. exact3d H-a 검증(위반0). frame3(nacre-scalar)는 nacre-predicates 무의존.
- **★ tol>0 경로·판정 라우팅**: 설계 §TIP "tol=0이면 지금 predicates 그대로(exact)·tol>0이면 필터+고정밀". 이 judge는 모든 상승을 astro-float로(격리 셀엔 무해) — **tol-0(Constructed)은 상위 층이 `nacre-predicates::orient3d`(Shewchuk exact·더 빠름·declare-0 없음)로 라우팅**(2b/stage 3). frame3가 nacre-predicates 무의존이라 judge가 Shewchuk 못 부름 → 라우팅은 상위 층. `JUDGE_PREC=200`은 튜닝 노브(stage 3 declare-0율 실측 후 상향).
- **★ 공유 dedup·`Orient` 루트화**: `frame::JUDGE_PREC`→`pub(crate)`·`Orient`는 `lib.rs`에 `pub use frame::Orient`(2D/3D/하류 공통 판정 타입 = `nacre_scalar::Orient`)·`bf_mag`/`HP_RM`/`rat_to_big` frame.rs 공유. exact3d 복붙 아님.
- **net.** `frame3`: **full H-a**(`h_a_det3_bound_soundness`·무리수각·**피벗 포함**·5000·위반0·tightness 0.065·**`#[ignore]`**)·`det3_bound_port_check`(빠른 40개·default)·`orient3d_sanity_and_rotation_invariance`(사면체 Positive·거울 Negative·공유 회전 불변)·`orient3d_coplanar_is_zero`(정확 공면 tol0 → Zero). **성능 정책**: 무거운 astro-float soundness(2a `tol_bounds_error...`·2a-ii H-a)는 `#[ignore]`(oracle 패턴·`--ignored` on-demand), default 스위트엔 빠른 것만. **n2(OCCT) 없음**·boolean 미배선("층만"). §TIP② 직접분 소진(남은 수학=간접 2c). 후속: 2b nacre-tip forest 브리지·2c 간접 orient3d.

**nacre-tip: forest→Pt3 브리지 + 커널 정점 orient3d (오버홀 단계 2b).** 수치 판정 층(2a/2a-ii)을 **커널에 연결** — 1c forest의 **첫 소비자**. 회전 정점을 `Pt3`로 조립(base=루트 좌표·chain=forest 순회)해 실제 tol을 계산·판정(end-to-end 관찰 가능). 신규 크레이트 `nacre-tip`(topo/scalar/predicates 위 read-only 분석·nacre-props 패턴).
- **★ `Rat::try_from_f64`(nacre-scalar)**: 커널 f64 좌표를 Pt3.base `[Rat;3]`로 **정확** 변환(비트 분해·i128 오버플로→None=§4). round-trip 검증. f64 캐시→exact 정의 브리지.
- **★ `vertex_pt3(model, vh) -> Result<Pt3, TipError>`**: `Rotated{base,rotation}`→base를 루트(비-Rotated·1c 불변식)까지 따라가 `Pt3::at`, `rotation.parent` 사슬을 `RotNode`로 조립해 `rotate_about` 순차 적용. **★ 순수-회전 가드**: forest는 회전만 담아(§⑦) **회전-후-이동**이면 replay≠커널 coord·`hp_coord`가 이동 잃음 → replay가 커널 `vertex.point`와 **비트-일치**(같은 apply_point 연산) 확인, 불일치 `TranslateInterleaved`(honest defer). **Discovered**(seam)→`IndirectRequired`(2c 간접·직접 술어는 캐시-as-base라 경계서 unsound). soundness는 M5 "Constructed 좌표=진실" 규약 기준.
- **★ `orient3d(model, 4정점) -> Result<Orient, TipError>`**: 전부 tol-0(Constructed 미회전·90°계열)이면 **`nacre_predicates::orient3d`**(Shewchuk exact), 하나라도 tol>0(회전)이면 **`frame3::orient3d_judge`** — §TIP "tol=0=predicates·tol>0=필터+고정밀" 라우팅. 사유별 `TipError`(IndirectRequired·TranslateInterleaved·Downgrade) 뭉치지 않음.
- **net.** `nacre-scalar`: `try_from_f64_is_exact_and_round_trips`. `nacre-tip`: `rotated_vertex_pt3_matches_and_has_tol`(회전 큐브 정점 coord 비트-일치·tol>0)·`unrotated_vertex_is_tol_zero`·`orient3d_unrotated_matches_shewchuk`(tol-0 Shewchuk 경로)·`orient3d_rotated_agrees_on_definite`(judge 경로·definite서 Shewchuk 부호 일치)·`discovered_is_indirect_required`·`translate_after_rotate_is_deferred`. **n2(OCCT) 없음**·boolean 미배선("층만"). 후속: 2c 간접(Discovered seam 진짜 판정·`None` 경로가 여기로)·stage 3 boolean 배선.

**간접 orient3d — implicit point (오버홀 단계 2c-i).** 2b가 `IndirectRequired`로 defer한 **Discovered seam**을 판정하는 술어를 `frame3`에 이식(§3 "Constructed=direct·Discovered=indirect"). seam 정점 `V=∩(3평면)`의 `orient3d(V,q,r,s)` 부호를 **V를 f64로 재료화하지 않고** Cramer `sign(D)·sign(M)`(나눗셈 없음)로 — 회전 boolean이 가장 많이 소비할 술어(우연 공면·횡단 seam).
- **★ `Iv { mid, rad }` 동적 필터(§⑨)**: 대칭 오차 반경 인터벌(`sub`/`add`/`mul` worst-case 반경 + `2ε` f64-반올림 팽창)·`sign()->Option<bool>`(0을 걸치면 `None`=상승). **술어별 오차 한계 공식을 손유도하지 않음** — 인터벌 산술이 구성상 sound. `det3_iv`·`pt_iv`.
- **★ `plane_iv(&Pt3×3)->[Iv;4]`(§② 간접분·§539 따름정리)**: 회전 평면 계수 `n=(p1−p0)×(p2−p0)`·`d=−n·p0`, 점 tol이 뺄셈·외적·내적으로 **계수 tol로 전파**. astro-float 쌍둥이 `plane_hp`(점 `hp_coord` 소비·자체 trig 없음). `indirect_filter`→`indirect_hp` 상승.
- **★ declare-0 floor는 간접 전용 정책**: 항 크기 합 `mag` 기반(`FLOOR_K·mag·2^-prec`·`FLOOR_K=1e6`), explicit orient3d(2a-ii)의 `16·scale³`과 **별개**(섞지 않음). 부호 추출은 2a-ii 인라인 관용구(`is_zero`/`bf_mag`/`is_positive`·astro-float#44 우회) **재사용**(별도 `bf_sign` 없음·frame3에 부호 방식 하나).
- **net.** `frame3`: **full H-b**(`h_b_plane_coefficient_tol_soundness`·평면 계수 tol이 실제 오차 상한·피벗 포함·10000·위반0·tightness 0.410·**`#[ignore]`**)·**full H-c**(`h_c_indirect_orient3d_soundness`·이질 provenance + near-coplanar 공유-피벗 2코퍼스·wrong-sign 0/3341·escalated·**filter_resolved 둘 다 >0**로 두 경로 실증·**`#[ignore]`**)·`iv_arithmetic_is_sound_and_sign_decides`·`indirect_sanity_and_rotation_invariance`(좌표평면 교점 V·q/r swap 부호 반전·공면 s→Zero·공유 회전 불변)·`indirect_port_check`(빠른 16개·default). **n2(OCCT) 없음**·exact3d 무변경(검증=frame3 재현)·boolean 미배선("층만"). §TIP② orient3d **직접+간접 모두 소진**. **정직 구분**: 회전-Pt3 cmp_coord는 오차 한계만 H-b가 커버(함수 구현은 exact3d에도 없음·stage-3 `loop_winding` 후속). 후속: 2c-ii nacre-tip 브리지(Discovered seam→평면 Pt3 조립·`IndirectRequired` 경로가 여기로)·stage 3 boolean 배선.

**nacre-tip 간접 브리지 — Discovered seam 판정 (오버홀 단계 2c-ii).** 2c-i 술어(회전 평면=세 회전 점)를 커널에 연결 — 2b가 `IndirectRequired`로 defer한 Discovered seam이 여기로. `orient3d(a,b,c,d)`에 **단일-Discovered 디스패치**를 접어 stage 3가 direct/indirect를 모른 채 한 진입점만 호출.
- **★ 회전 평면을 f64 계수가 아니라 회전 정점으로 재구성**: 비-90° 회전 시 `Plane::coefficients()`가 무리수라 "정확 술어×부정확 입력" 구멍. Surface는 회전 이력 무보유(origin/normal/raw만) → 정확 provenance는 회전 정점(`Origin::Rotated`)만. `plane_pts`가 `VertexDef::ThreePlane`의 각 surface를 쓰는 면의 **direct·비공선 정점 3개**(`vertex_pt3`→Pt3)로 평면 재구성. **정점 winding 무관**(간접 술어가 평면 법선 방향에 불변). **store 전체 스캔 sound**(append-only·superseded 면 정점도 그 평면 위).
- **★ 부호 치환**: Discovered를 V-슬롯(0번)으로 swap(치환 1회 = `Orient` 1회 반전·`slot!=0`이면 flip). 술어 `judge=orient3d(V,q,r,s)`. **≥2 Discovered**(다중 implicit)·정점 부족은 사유별 `TipError`(`IndirectRequired`·`PlaneUnderdetermined`) honest defer.
- **net.** `nacre-tip`: `indirect_matches_direct_on_rotated_corner`(회전 코너를 간접(3평면)·직접(명시점) 두 방식 판정 등가·**슬롯 0/1 모두**로 flip·no-flip 실증)·`indirect_unrotated_matches_exact_plane`(unrotated seam이 독립 exact-plane `nacre_geom::three_plane_orient3d`와 일치·교차검증)·`underdetermined_plane_is_deferred`·`two_implicit_points_are_deferred`. **n2(OCCT) 없음**·boolean 미배선("층만"·boolean이 아직 비-90° 회전 거절→회전 seam은 합성 테스트). `nacre-geom` dep 추가(`Handle<Surface>` 타입명·테스트 오라클). 후속: stage 3 배선(가드 은퇴·unrotated→exact-plane 성능 라우팅·seam 내외 판정)·다중 implicit `cmp_coord`.

**회전 cmp_coord 술어 — 두 implicit point 축좌표 비교 (오버홀 단계 cmp-i).** §TIP②의 마지막 프리미티브. 두 seam point의 한 축 좌표를 재료화 없이 비교(`sign(a[axis]−b[axis]) = sign(Dvec_a[axis]·D_b − Dvec_b[axis]·D_a)·sign(D_a)·sign(D_b)`·나눗셈 없음·nacre-predicates `indirect_cmp_coord` 규약 일치). exact3d엔 없던 새 수학(orient3d만 검증)이라 frame3에서 직접 검증. stage-3 `loop_winding`·seam 정렬이 소비.
- **★ 사전 chore**: `frame.rs`(2D)→`frame2.rs`로 리네임(3D `frame3`과 대칭·nacre-scalar 내부 전용).
- **★ Cramer 헬퍼 추출(DRY)**: `indirect_filter`/`indirect_hp`의 `D`·`Dvec` 계산부를 `cramer_iv`/`cramer_hp`(+ `mag_dvec`)로 뽑아 orient3d·cmp 공유. 두 함수 시그니처 유지(H-c 테스트가 직접 호출)·본문만 교체 → **bit-identical**(H-c 수치 불변 확인).
- **★ `indirect_cmp_coord_judge`**(pub·반환 `Orient`): `Iv` 필터→astro-float 상승·세 부호곱(`cmp_combine`·짝수 음수→Positive)·`sign_with_floor`/`FLOOR_K` 재사용. **`Zero` 의미**: exact 판의 0("정확히 같음")과 달리 "같음 **또는** floor 아래 분리불가"(stage 3가 다른 축 fallback으로 다뤄야). **winding 불변**(트리플 법선 뒤집으면 `D`·`Dvec` 동시 반전 → cmp-ii가 아무 세 점으로 평면 정의해도 sound·대수 확인).
- **★ `FLOOR_K` 재사용 실측 확인**: cmp의 M은 orient3d보다 고차수(≈13차)라 자동 타당 아님 → **H-g가 위반0으로 확인**(현 1e6 적정).
- **net.** `frame3`: **full H-g**(`h_g_indirect_cmp_coord_soundness`·이질 provenance + near-tie[같은 축 회전으로 tie 보존] 2코퍼스·wrong-sign **0/3177**·escalated 1360·filter_resolved 2640 두 경로 실증·**`#[ignore]`**)·`cmp_sanity_axis_aligned`(축정렬 exact 점 known 순서·swap 반전·같음→Zero)·`cmp_port_check`(빠른 16개)·`cmp_combine_counts_negatives`. **orient3d 무회귀**(리팩터 bit-identical·H-b/H-c·sanity green). **n2(OCCT) 없음**·exact3d 무변경·boolean 미배선("층만"). §TIP② **완전 소진**(orient3d 직접+간접·cmp_coord). 후속: cmp-ii nacre-tip 브리지(두 Discovered→`plane_pts`×2→judge)·stage 3 배선.

**nacre-tip cmp_coord 브리지 — 두 Discovered seam 축좌표 비교 (오버홀 단계 cmp-ii).** cmp-i 술어를 커널에 연결(2c-ii가 orient3d에 한 것과 대칭). `cmp_coord(model, v1, v2, axis)`: 두 정점 모두 Discovered seam이면 각 `ThreePlane`을 `plane_pts`(2c-ii 재사용)로 세 회전 정점 트리플로 재구성→`indirect_cmp_coord_judge`. 반환 `Orient`(Positive=`v1[axis]>v2[axis]`). **winding 무관**(cmp-i 대수 확인)이라 아무 세 점으로 평면 정의해도 sound. 회전 평면 계수가 무리수라 정점 재구성이 핵심(2c-ii와 동일).
- **honest defer**: 비-Discovered 정점→`NotSeam`·surf 정점 부족→`PlaneUnderdetermined`.
- **net.** `nacre-tip`: `cmp_matches_coord_on_rotated_corners`(**독립 변환된 두 큐브** 코너쌍=heterogeneous provenance·여러 축·definite에서 f64 좌표 부호와 wrong-sign 0·`resolved>0`[항상-Zero 버그 방지]·swap antisymmetric)·`cmp_unrotated_matches_exact_plane`(독립 exact-plane `three_plane_cmp_coord`와 교차검증)·`cmp_non_seam_is_deferred`·`cmp_underdetermined_plane_is_deferred`. **n2(OCCT) 없음**·boolean 미배선("층만"·회전 seam 합성 테스트). **★ TIP 판정층 완성**: stage-3 탐색으로 확정 — boolean은 모든 점(seam·원본 정점)을 평면-트리플로 통일해 정렬이 전부 single-implicit orient3d(`order_along`→2c-ii)·two-implicit cmp(`loop_winding`→cmp-ii)뿐이라 **혼합 cmp는 없음**(cmp-iii 불필요). 술어(orient3d 직접+간접·cmp_coord·D부호) 모두 frame3/tip에 있음. 후속: stage 3(배선). **(정정 — 3c-i: 이 "완성"은 arrangement/seam 정렬에 한정됐다.** `point_in_solid`의 **ray casting**은 이 탐색 밖이었고, `orient3d(p, p+d, ·, ·)`의 방향 항이 새 술어 `dir_orient3d`를 요구했다 — §TIP② 새 수학이 아니라 한 열 exact인 기존 Iv 프레임의 적용이라 §TIP② 소진 상태는 유지. 아래 3c-i.)

**stage 3 착수 — toleranced t_orient3d + PlaneInfo Pt3 정의 (오버홀 3a-i).** 회전 boolean 배선의 첫 웨지. **핵심 통찰**(탐색): arrangement(arrange.rs)가 이미 평면-트리플 위 순수 sign-predicate라, 술어에 먹이는 `PlaneInfo.tri`(회전 시 반올림 f64)만 정확 정의로 바꾸면 통째로 rotation-sound. **분해**: 3a(라우팅 층·미배선)→3b(arrange 교체·가드 유지)→3c(point_in_solid toleranced·shell·validate)→3d(가드 은퇴+OCCT·main 병합). 소프트니스 가드는 **부분 은퇴 불가**라 3d가 최종 게이트.
- **★ `PlaneInfo.tri_pt3: Option<[Pt3;3]>`**: 회전 시 `collect_planes`가 면 tri 정점(`outer_tri`가 핸들도 반환→`tip::vertex_pt3`)을 세 Pt3로 **한 번 캐시**(O(N²) order_along가 매 호출 forest 순회 안 하도록). unrotated면 `None`(무비용·exact 경로). `PlaneInfo`는 derive 없는 plain 구조체라 Pt3(Clone) 필드가 Copy 안 깸.
- **★ `tolerant::t_orient3d(planes, p,q,r,j, rotated) -> i8`**: `order_along` 시그니처와 동일(3b drop-in). `!rotated`→기존 `three_plane_orient3d`(핫패스 무변경). `rotated`→`tri_pt3` 세 Pt3로 `frame3::indirect_orient3d_judge`. 라우팅=boolean당 `solid_is_rotated` flag(매 호출 tol 계산 없음). winding 불변.
- **net.** `nacre-ops`: `t_orient3d_rotation_invariant`(**회전 불변 오라클**: 회전 큐브의 frame3 경로 == 같은 큐브 unrotated geom 경로·definite 전부·조립+라우팅 검증·술어 soundness는 frame3 H-c 계승·astro-float GT 불요)·`plane_def_from_face`(tri_pt3 조립·좌표 일치·tol>0)·`t_orient3d_unrotated_forwards_geom`(geom bit-일치). **미배선**(order_along 라이브 무변경·가드 유지)·`t_orient3d`/`tri_pt3`는 3b까지 `#[allow(dead_code)]`. **n2(OCCT) 없음**(3d)·기존 boolean 무회귀. ops→nacre-tip dep(lib 그래프 무순환·tip의 ops는 dev-only). 후속: 3a-ii(cmp/dir_sign/plane_side 래퍼)·3b(배선).

**toleranced t_cmp_coord + t_plane_side (오버홀 3a-ii).** 라우팅 층의 나머지 술어 2종(기존 frame3 judge 재사용).
- **★ `t_cmp_coord`**(loop_winding·두 implicit point 축좌표): `!rotated`→`three_plane_cmp_coord`·`rotated`→각 트리플
  세 Pt3로 `indirect_cmp_coord_judge`. cmp는 좌표라 **회전 불변 아님** → 테스트는 좌표 오라클(cmp-ii 동형).
- **★ `t_plane_side`**(straddle·정점 vs 면): `!rotated`→`plane_side`·`rotated`→정점 Pt3+면 세 Pt3로 `orient3d_judge`
  (직접 judge). orient3d 부호라 **회전 불변** → 테스트는 회전 불변 오라클(3a-i 동형).
- **★ 혼합 회전 버그 수정(래퍼-측)**: 3a-i가 `tri_pt3`를 솔리드별로 채워, 혼합 부울(A 회전·B 축정렬)의 B-plane은
  `None`이라 rotated 경로가 panic할 수 있었음. **`plane_def`가 Some→clone(값싼 chain 복사·forest 순회 없음)·None→
  tri 좌표로 즉석 exact 빌드**(축정렬 좌표는 exact f64·tol 0). collect_planes 시그니처 churn 없이 수정(호출부 ~40 무변경).
- **net.** `nacre-ops`: `t_plane_side_rotation_invariant`·`t_cmp_coord_matches_coord`(좌표 오라클·resolved>0·antisym)·
  둘 다 `_unrotated_forwards_geom`·`mixed_rotation_handled`(혼합서 None plane 즉석 빌드·panic 없음). **미배선**·
  `#[allow(dead_code)]`(3b까지). 후속: 3a-iii(t_dir_sign·frame3 D부호 노출)·3b(배선).

**toleranced dir_sign (frame3 dir_sign_judge + ops 래퍼) — 3a 라우팅 층 완성 (오버홀 3a-iii).** 4번째 술어
`plane_pair_dir_sign`(법선 det·`sign((n_a×n_b)·n_c)`·`order_along`의 `dir_sign`/`turn_at`이 씀) 라우팅.
- **★ frame3 `dir_sign_judge`**(pub): 세 평면 법선 det 부호=cramer `D`(`cramer_iv`/`cramer_hp` 재사용·D만·M 불요).
  **winding 의존적**(법선=`(p1−p0)×(p2−p0)`이라 점 순서가 방향 정함·호출자가 일관 순서).
- **★ ops `t_plane_pair_dir_sign`**: `!rotated`→geom. `rotated`→각 면 **outward tri**의 frame3 `D`를
  `plane_pair_dir_sign`의 **저장-법선 규약**으로 다리 — `det(stored)=orient_sign(p)·orient_sign(a)·orient_sign(b)·
  det(outward)`(`orient_sign` pub(crate)화·평행 단위벡터 dot이라 회전에 robust). 회전 불변 테스트가 이 규약 검증(통과).
- **★ H-i(신규·`#[ignore]`)**: dir_sign soundness를 **법선 near-coplanar 코퍼스**(D≈0·H-a/H-c/H-g 미겨냥 regime)에서
  512-bit stable GT와 wrong-sign **0/1325**(escalated·filter_resolved 양경로). **코퍼스는 퇴화 평면(tiny normal) skip**
  — near-coplanar-*normals* regime이 아니고 boolean 입력도 아님(발견: 그렇지 않으면 false-positive wrong). 판정 자체는
  단일 정밀도(다른 judge와 동일)로 sound.
- **net.** `frame3`: `dir_sign_judge_sanity`·`dir_sign_rotation_invariant`·H-i(`#[ignore]`). `nacre-ops`:
  `t_dir_sign_rotation_invariant`(orient_sign 규약 검증)·`t_dir_sign_unrotated_forwards_geom`. **★ 3a(라우팅 층) 완성**
  — 4술어(orient3d·cmp_coord·plane_side·dir_sign) 전부 toleranced 래퍼. **미배선**·`#[allow(dead_code)]`(3b까지). 후속:
  3b(arrange 배선·`rotated` 스레딩·가드 유지).

**arrange 술어 배선 — index-based 3종·첫 라이브 배선 (오버홀 3b-i).** 3a 래퍼 중 평면-인덱스 기반 3종을 라이브
arrange.rs 호출로 교체 — 회전 boolean으로 가는 첫 라이브 코드 변경. `t_plane_side`(straddle)는 정점 핸들+`model`
스레딩이 필요해 3b-ii로 분리·가드 유지(3d 은퇴).
- **★ per-predicate 파생 라우팅 — `rotated` 스레딩 폐기**: 3a 계획의 "`rotated` flag 스레딩" 대신, 래퍼가 관련 평면의
  `tri_pt3.is_some()`(단일 헬퍼 `any_rotated(planes, &[idx])`)로 **호출마다 자동 분기**하게 리팩터(`rotated: bool`
  파라미터 제거). 그래서 arrange는 시그니처 변경·flag 배관 없이 **호출만 교체**. per-predicate라 혼합 부울(A 회전·B
  축정렬)에서 축정렬-only 술어는 geom(빠름) 유지 — per-boolean flag보다 정밀. **왜 `all`이 아니라 `any`**: 관련 평면
  중 하나라도 무리수 좌표면 그 f64 근사가 orient3d/cmp를 뒤집을 수 있어 frame3 필요. "무리수 여부"는 소비자가 아니라
  **분류기 소유**(`solid_is_rotated`←`Origin::Rotated`←`Isometry::is_exact`) — 90°-계열은 exact `Constructed`라
  `tri_pt3=None`→geom 자동. 향후 exact 케이스 확대는 분류기 층만의 최적화이고 이 헬퍼는 플래그만 읽어 무변경.
- **★ arrange 7 호출부 교체**: `order_along`·`side_of`(`t_orient3d`)·`loop_winding`(`t_cmp_coord`·죽은 로컬 `tri`
  클로저 제거·`ring`이 이미 `[usize;3]`)·`dir_sign`·`turn_at`·`point_on_ring`·`every_ray`(`t_plane_pair_dir_sign`).
  미사용된 geom import 3종 제거·인트라-닥 링크는 전체 경로화(import 제거 시 `unresolved link` vs 유지 시 `unused_imports`
  catch-22 회피). `PlaneInfo.tri_pt3` 필드 allow도 제거(배선으로 라이브 read).
- **★ 지연 의무 — 순서 판정 자리의 declare-0(3d)**: frame3는 `Orient::Zero`(declare-0)를 낼 수 있는데, 이를 조용히
  비-0으로 읽으면 silent-wrong. 7자리 감사: **(A)** `order_along`(전순서 가정)·**(B)** `loop_winding` lex 비교(부분
  declare-0→잘못된 lex-min `lo`→hull 아닌 정점 winding)는 3d 가드 은퇴 시 정직한 거절로 매핑해야 함. **(C)** `side_of`·
  `point_on_ring`·`every_ray`는 0이 설계된 출력. **(D)** `turn_at`은 이미 `det==0` 자기 거절. 이 셀은 가드가 회전을
  막아 라이브 무영향·의무만 기록.
- **net.** `nacre-ops`: 4 래퍼 `rotated` 파라미터 제거·`any_rotated`·기존 tolerant 테스트 10종 인자만 갱신(모든 `true`
  호출이 회전 평면 집합·`false`가 all-None이라 파생이 bit-identical 보존)·`t_plane_side`는 3b-ii까지 `#[allow(dead_code)]`
  유지. **검증**: unrotated 부울 전 스위트 무회귀(전 평면 `tri_pt3=None`→geom 파생→bit-identical)·`cargo doc` 무경고.
  **n2(OCCT) 없음**(3d)·가드 유지(frame3 라이브 미실행). 후속: 3b-ii(`t_plane_side`·`edge_crosses_face` 핸들+model
  스레딩)·seam-gen 술어(4-plane guard·is_convex)·3c.

**t_plane_side 배선 — arrangement 4번째·마지막 술어 (오버홀 3b-ii).** 남은 arrangement 술어 `t_plane_side`(straddle)를
`edge_crosses_face`에 배선. 이로써 arrangement 4술어 전부 배선 완료. `edge_crosses_face`의 링 포함(`point_on_ring`/
`point_in_ring`→`every_ray`)은 이미 3b-i로 toleranced였고, 유일한 f64 잔여였던 straddle(`plane_side(planes[q].tri,
p0/p1)`)만 남았었다.
- **★ 정점 핸들+model 스레딩**: `t_plane_side`가 정점 좌표(geom)와 정의(frame3·`vertex_pt3`) 둘 다 필요해
  `edge_crosses_face`(arrange.rs)·`pierced_faces`(lib.rs) 시그니처를 `p0/p1: Point3`→`v0/v1: Handle<Vertex>`+`model`로.
  **모든 호출부**(seam_segments_on·overlap_fuse_cut·boundaries_intersect·테스트 5)가 이미 `bounds`+model 보유·p0/p1을
  거기서 계산 중이라 핸들 전달로 그 계산이 사라짐(호출부 단순화).
- **★ 왜 인덱스 재구성(`t_orient3d`, 스레딩 불요)이 아닌가**: 끝점을 세-평면 meet `{e0,e1,r}`로 보면 straddle이
  인덱스 전용이 되지만, **Constructed 원본 정점은 세-평면 정의가 없다**(진실=좌표·§4)·세 번째 평면 `r`이 인접 면
  순회·정점당 3+면이라 비유일. 정점 truth 존중하는 `t_plane_side`(Constructed=좌표·Rotated=`vertex_pt3`·
  Discovered=[3c]indirect)가 올바른 접근·스레딩 정당.
- **★ declare-0 자기 보호**: straddle `s0==0||s1==0→VERTEX_ON_FACE_PLANE`. frame3 spurious declare-0도 정직한
  거절로 흡수(turn_at 동형)·별도 3d 정책 불요.
- **★ 정직한 한계(3c)**: `t_plane_side` rotated 경로 `vertex_pt3`는 `Discovered⇒IndirectRequired`라 **Discovered
  끝점이면 패닉**(boolean 결과를 회전 피연산자로 넣으면 그런 코너·design.md 위). 3b-ii는 가드가 회전 차단→geom
  경로라 안전·3c가 Discovered 끝점을 indirect orient3d로 라우팅해야 3d 가드 은퇴 시 패닉 면함(is_convex/
  point_in_solid과 같은 정점-origin 계층). 4-plane 동시성 가드도 orient3d+`planes_coplanar` 함께라 3c.
- **net.** `arrange.rs`: `edge_crosses_face`(model+핸들·straddle→`t_plane_side`·doc 산문/링크 갱신)·`plane_side`
  import 제거(미사용·`plane_side` 링크 전체경로화)·`three_planes` 유지. `lib.rs`: `pierced_faces`(model+핸들)·3 호출부·
  테스트(`pierce` 헬퍼가 `bounds` 반환). `tolerant.rs`: `t_plane_side` allow 제거·모듈 doc(4종 배선). **검증**: unrotated
  무회귀(228+통합 bit-identical)·`cargo doc` 무경고. **n2(OCCT) 없음**(3d)·가드 유지. 후속: 3c(point_in_solid·
  is_convex·t_plane_side Discovered·4-plane 가드·공면 감지기 — 정점-origin 계층·규모 큼).

**frame3 `dir_orient3d_judge` — 방향 orient3d (오버홀 3c-i).** 3c는 정점-origin/분류 계층. 탐색으로 **횡단 회전
Fuse/Cut의 라이브 임계 사이트 셋** 확정(point_in_solid·4-plane 가드·is_shell_outward). headline `point_in_solid`
(내외 분류)을 회전 sound화하려면 ray-삼각형 판정을 반올림 좌표가 아니라 정점 exact 정의·frame3 위에서 해야 하는데,
그 5개 orient3d가 **단 하나의 새 술어**로 환원된다.
- **★ `dir_orient3d_judge(d, base, x, y)`**(frame3 신규 pub): `sign(det[d, x−base, y−base]) = sign(d·((x−base)×
  (y−base)))`. **유도**: ray-삼각형은 `orient3d(p, p+d, ·, ·)`를 묻는데 `p+d`(회전점+rational 방향)는 base+chain
  Pt3로 표현 불가(`R⁻¹d` 무리수) — 그러나 각 판정이 방향 행렬식으로 환원되고 `d`가 **한 열 exact(rad-0)**로 들어간다.
- **★ Iv-직접-열 구현(구성상 sound)**: `det3_iv([d행, x−base, y−base])` 필터(`Iv::sub`가 차분 반올림 추적)→straddle
  시 `det3_big`(hp)+`sign_with_floor`. **`p+d`를 안 만들어** f64 덧셈 반올림 함정 회피. `det3_hp`에서 **임의-벡터
  `det3_big`를 추출**(bit-identical)해 공유. mag는 f64 행 6-곱-절대합(`det3_mag`).
- **★ 검증 — 오라클은 고정밀 GT(회전-불변 아님)**: `d`가 고정이라 점만 돌리면 부호가 바뀌므로 회전-불변 불가 →
  GT-stable 512 vs 640-bit det 부호. 코퍼스 두 영역(랜덤 + **near-grazing**: 회전 후 in-plane edge의 큰 유리수 근사로
  `d`가 회전 normal에 거의 ⊥ → det≈0 → escalation). **wrong-sign 0/1499·escalated 603·filter_resolved 897**(양 경로)·
  sanity(축정렬 known 부호·edge swap 반전·in-plane→Zero). **rotate-together sanity 생략**: `d`가 `[Rat;3]`라 일반
  회전 `R·d`가 무리수(90°만 유리수)라 어색 — GT 코퍼스가 실제 사용(고정 유리수 d + 회전 점)을 직접 검증해 대체.
- **net.** `frame3`: `dir_orient3d_judge`(pub)·`det3_big`(추출)·`det3_mag`·`dir_orient3d_judge_sanity`·
  `h_dir_orient3d_soundness`(`#[ignore]`). **exact3d/기존 frame3 무회귀**(det3_big 추출 bit-identical·H-a/b/c/g/i
  ignored 7건 green). **미배선("층만")**·boolean·point_in_solid 무변경. §TIP② 소진 유지(한 열 exact인 Iv 프레임
  적용). 후속: 3c-ii(ray_triangle_cross 조립)·3c-iii(point_in_solid 배선).

**toleranced ray_triangle_cross 조립 (오버홀 3c-ii).** 3c-i 방향 술어로 `point_in_solid`이 매 면-삼각형마다 부르는
**광선-삼각형 교차 판정**의 toleranced 버전을 ops에 조립·검증(unwired). 전체 배선(핸들 스레딩·퇴화 드롭·라우팅)은 3c-iii.
- **★ 부호 환원**(위험 지점): `ray_triangle_cross`(predicates)의 5 orient3d를 3 `orient3d_ray`(광선-선 edge:
  `orient3d(p,p+d,·,·)`)·1 `orient3d_judge`(s0)·1 `dir_orient3d_judge`(sd)로. declare-0→`Degenerate`(광선 재시도 흡수).
- **★ frame3 `orient3d_ray(base,dir,x,y) = sign(orient3d(base, base+dir, x, y))`**(신규 pub): `p+d`가 base+chain Pt3로
  표현 불가한 것을 열-축소로 **`dir_orient3d_judge(dir, y, x, base)`**(순수 인자 순열·부정 불요·마지막 두 인자 swap이
  부호 흡수)로 환원. 호출부가 f64 원본과 **인자 대 인자 일치**해 실수 최소.
- **★ 건전성 상속·조립만 검증**: 술어(dir_orient3d H-corpus·orient3d_judge H-a/c)가 sound, 조립 로직은 route-무관 →
  **unrotated 정수좌표 bit-identical**(부호 환원 검증)이면 회전에서도 sound(충분성 논거·3b/3c-i 동형).
- **net.** `frame3`: `orient3d_ray`(pub)·`orient3d_ray_matches_materialized`(unrotated서 materialized orient3d와 일치).
  `nacre-ops`: `ray_triangle_cross_tol`(private·`#[allow(dead_code)]`·3c-iii 배선 시 제거)·`ray_triangle_cross_tol_
  matches_f64_unrotated`(정수좌표 bit-identical·6 RAY_DIRECTIONS)·`ray_triangle_cross_tol_is_rotation_equivariant`
  (90° 회전 전체를 돌려 불변·Pt3 chain exercise). **미배선("층만")**·`point_in_solid`·boolean 무변경·가드 유지·기존
  전 테스트 무회귀. 후속: 3c-iii(분류기)·3c-iv(배선).

**point_in_solid_tol — 회전 내외 분류기 (오버홀 3c-iii).** `ray_triangle_cross_tol`로 회전 하 정점 내외 분류기를
ops에 짓고 회전-불변 검증(unwired·배선은 3c-iv). `point_in_solid`(2213·Point3·f64) 미러 — 면 삼각분할해 forward-ray
winding.
- **★ `point_in_solid_tol(model, p: &Pt3, solid)`**: 좌표→정점 exact 정의(Pt3), `ray_face_cross`→`ray_triangle_cross_tol`,
  `RAY_DIRECTIONS`→[Rat;3]. 면 정점 핸들은 신규 `face_loop_verts`(`face_loops` 미러·`he_start`)→`vertex_pt3`(`Discovered`
  면 정점이면 `ROTATED_UNSUPPORTED` 정직 거절·신선 회전 한정). winding≠0은 광선 방향 무관(닫힌 곡면)이라 고정 방향으로도
  회전 프레임서 정확. **두 피연산자 모두(서로 다른 R) 회전까지 지원**(각 Pt3가 참 world 위치·orient3d 프레임 무관).
- **★ 퇴화 in-loop 처리(사전 드롭 없이)**: `Degenerate` 뜰 때만 **`pt3_base_collinear`**로 판단 — 진짜 collinear면 skip
  (넓이 0·기여 0), 아니면 graze라 다음 광선 재시도. f64판의 `fan_triangles`+`triangle_is_degenerate` 사전 드롭보다 별도
  pass·헬퍼 없이 같은 결과.
- **★ `pt3_base_collinear`가 soundness-critical**: 회전 f64 넓이는 near-collinear를 0으로 반올림해 진짜 삼각형을 거짓
  skip→교차 누락→오분류(silent-wrong)((5d)#4 동류). 회전이 collinearity 보존·same-solid=same-chain이라 **`Pt3.base`
  (회전 전 유리수) 3-투영 rational orient2d가 전부 0일 때만** collinear. i128 overflow→not-collinear(최악은 정직한
  `RAY_DEGENERATE`·거짓 skip 없음). unrotated면 base=coord라 현행과 동치.
- **★ 왜 제네릭화 아닌 병렬 함수**: winding 루프 공유하려면 라이브 `point_in_solid`(f64)를 건드려야 해 bit-identical
  무회귀가 위태 → 병렬 함수로 두어 f64 경로 정의상 무변경(안전>DRY).
- **net.** `nacre-ops`: `point_in_solid_tol`·`face_loop_verts`·`pt3_base_collinear`(private·`#[allow(dead_code)]`·3c-iv
  배선 시 제거)·`point_in_solid_tol_is_rotation_invariant`(회전 큐브+오목 L-prism[reflex notch 점 포함] 분류=unrotated
  `point_in_solid` 정답)·`pt3_base_collinear_exact`(회전 collinear→skip·sliver→retry). **미배선("층만")**·`point_in_solid`·
  boolean 무변경·가드 유지·기존 전 테스트 무회귀. 후속: 3c-iv(라우터 배선).

**vertex_in_solid 라우터 — point_in_solid_tol 배선 (오버홀 3c-iv).** 회전 인식 라우터로 일반 부울의 내외 분류를
라이브에 연결(내외 분류의 첫 라이브 회전 경로).
- **★ `vertex_in_solid(model, vh, solid)`**(신규 private): `solid_is_rotated(solid) || 정점 vh Origin::Rotated`면 tol
  (`vertex_pt3(vh)`[Discovered→`ROTATED_UNSUPPORTED` 정직 거절]→`point_in_solid_tol`), 아니면 **현행 `point_in_solid`
  (정점 coord)** — bit-identical. 두 항 다 필요: 대상 솔리드 회전(면 무리수) 또는 질의 정점 회전(혼합 부울·좌표 무리수)
  이면 tol.
- **★ 일반 부울 2 호출부만 배선**: 1602(seam-free contained 분류)·1706(overlap_fuse_cut 분류)→`vertex_in_solid`. **공면
  감지기 point_in_solid(3617/4330)은 배선 안 함** — 그 감지기는 `is_convex`·`n_out.dot`(f64)도 써서 point_in_solid만
  바꾸면 반쪽 배선 → 감지기 전체와 함께 공면-접촉 밀레스톤으로(횡단 부울은 coplanar-pair 게이트서 early-None이라 비임계).
- **★ dead_code allow 제거**: `point_in_solid_tol`·`face_loop_verts`가 라우터→라이브 호출부로 정적 reachable(회전 분기
  미실행이어도)이라 allow 불요(3b-i `tri_pt3` 선례).
- **net.** `nacre-ops`: `vertex_in_solid`·1602/1706 교체·두 allow 제거·`vertex_in_solid_routes_by_rotation`(축정렬→
  `point_in_solid`·회전→`point_in_solid_tol` dispatch 확인). **가드 유지**(라이브는 항상 unrotated 분기→bit-identical·
  tol 분기 미실행까지 3d)·기존 boolean 무회귀·**n2(OCCT) 없음**(3d). 후속: 3c-v(4-plane 가드·is_shell_outward).

**seam 4-plane 동시성 가드 배선 (오버홀 3c-v).** seam 정점 생성(`general_boolean`)의 **4-plane 동시성 가드**를 회전
toleranced로 배선. 가드는 seam 정점(평면 e0·e1·entry의 meet)이 *다른* 4번째 평면 위에 있는 퇴화(4-plane meet)를 잡아
`FOURPLANE` 거절 — 회전 시 평면 계수·`tri`가 반올림 무리수라 f64 `three_plane_orient3d`가 오판→놓친 거절→틀린 seam
(silent-wrong) 가능.
- **★ 한 표현식 교체**: `three_plane_orient3d(&planes[e0/e1/entry].plane, pm.tri[0..2]) == 0` →
  `crate::tolerant::t_orient3d(&planes, e0, e1, entry, m) == 0`(3b-i 인덱스 래퍼 재사용·스레딩 불요·`pm`=`planes[m]`이라
  `pm.tri`=`planes[m].tri`, unrotated 경로 bit-identical). 신규 코드·테스트 없음.
- **★ `planes_coplanar` 트윈-스킵(가드 앞 gate)은 무변경·이 셀서 sound**: 진짜 트윈(같은 Surface 공유 두 면)은
  `t_orient3d`에 도달조차 못 함 → 가드는 **구별되는 평면만** 판정. 트윈 exact의 근거는 **계수 출처** — `collect_planes`가
  평면을 `*model.surfaces.get(face.surface)`에서 읽지(면별 `tri` 재계산 아님) → 같은 Surface 공유 면은 회전 무관
  **byte-identical 계수** → `planes_coplanar`(rank≤1 비례성)에 정확히 비례(tolerance 무관). fuzzy는 *다른*-Surface
  near-coplanar 회전 쌍뿐 → 공면-접촉 밀레스톤(`planes_coplanar` toleranced화)으로 이연(신선 primitive 횡단엔 부재).
- **★ 새 테스트 없음(정당)**: 가드의 `t_orient3d`는 **매 seam 정점 × 다른 모든 평면마다 실행**(`for (m, pm) in
  planes.iter().enumerate()`)되므로 seam 경로를 밟는 **기존 boolean 테스트 전부가 가드 코드를 exercise** — bit-identical
  이라 결과 불변→무회귀가 wiring 검증. 회전 건전성은 `t_orient3d` 자체가 이미 tested(3a-iii H-c·3b-i forwards-geom).
  FOURPLANE **거절 자체는 미발화 백스톱**(§637)이나 가드 코드 경로는 실행됨.
- **net.** `nacre-ops`: 가드 한 줄 교체(import 없음·fully-qualify). `three_plane_orient3d`는 다른 호출부(2068)서 유지.
  `comp_key`(f64 coord 정렬 tie-break·회전에도 결정적)·`is_shell_outward`는 이 셀 아님. 후속: 3c-vi(`is_shell_outward`
  회전-aware·구조적 — lex-min 정점 exact·outward 판정 회전-불변).

**is_shell_outward 회전-exact — 마지막 횡단 f64 자리 종료 (오버홀 3c-vi).** `assemble_fuse_cut`이 결과의 각 연결
컴포넌트를 outer(material·+)/cavity(void·−)로 가르는 `is_shell_outward`((5d)#5의 extreme-vertex 부호)는 **마지막 횡단
임계 f64 자리**였다. 그 답이 곧 위상 라벨이라 회전 시 틀리면 silent-wrong(sever 조각을 cavity로 오판). 두 수치 단계가
회전 하 깨진다: **lex-min v\* 탐색**(f64 좌표 배열 비교)·**outward 판정**(`plane.coefficients()[0]*orient_sign` — 축정렬
법선 가정). **알고리즘은 회전에도 정확**(v\*=min-x 극점·`∃ 인접면 n_x<0`이 outward/void를 가름)이라 두 수치 단계만
exact화. **★ 신규 술어 0개** — 둘 다 이미 배선·검증된 프리미티브로 환원.
- **★ 조립(3c-vi-a)**: `component_is_outward_tol(planes, comp: &[&LocalFace])`(회전-exact 트윈·`model` 무접촉).
  **outward** = `(flip?−1:1)·dir_orient3d_judge([1,0,0], tri_pt3…)`(3c-i·결과 면 outward 법선의 x-성분 부호 = RH-법선
  x-성분×flip; `plane_def`로 혼합-회전 안전). **lex-min** = 각 노드를 `loop_triples`(arrange) 동형 트리플[own+두 loop-edge
  타평면·meet=그 정점·원본/seam 통일]로 표현→`t_cmp_coord`(3b-i) 축별 argmin(`loop_winding` 형태). 비단순/직선각/
  non-manifold→정직 거절. **★ 조립-전 `LocalFace`에서 작동** — 조립이 결과 정점의 회전 provenance를 버리는(Orig→
  Constructed·Seam→Discovered{반올림 계수}) 문제를 구조적 우회(정확 정보는 `planes`+`comp`에). 4테스트(축정렬 f64 등가·
  회전 불변·혼합-회전 무패닉·비볼록 L-prism)·unwired.
- **★ 배선(3c-vi-b)**: `assemble_fuse_cut` positives 계산을 **컴포넌트별 `any_rotated` 인라인 라우팅**(회전→
  `component_is_outward_tol`·축정렬→`is_shell_outward` 무변경)으로. 라우터 추출 안 함(단일-site·dispatch 테스트가
  inverted-routing 못 잡음[tol이 축정렬도 정확]→3c-v식 무테스트). `any_rotated` pub(crate). 비회전 bit-identical
  (라이브 planes 전부 tri_pt3=None→f64 분기·positives 오름차순 c 순서 불변)·전 코퍼스 무회귀. **★ 이로써 횡단 임계 f64
  자리 전부 종료** → 3d(가드 은퇴·회전 부울 end-to-end·OCCT diff·main 병합)가 회전 부울을 처음 라이브 구동. 기각:
  provenance 보존(chained 라우팅 파급)·exact signed-volume(합 부호 인프라 없음, (5d)#5 근거). 후속: 3d.

**회전 부울 첫 라이브 — ROTATED_UNSUPPORTED 진입 가드 은퇴 (오버홀 3d-i).** 3a~3c-vi로 회전 횡단 부울의 모든 임계 f64
술어가 TIP-exact가 됐으므로, `boolean()` 진입부의 blanket 가드(`solid_is_rotated(a)||solid_is_rotated(b)`)를 제거해
완성된 기계를 회전 입력으로 처음 라이브 구동. 유일한 프로덕션 변경은 가드 3줄 삭제. `ROTATED_UNSUPPORTED` 태그·
`solid_is_rotated`는 잔존(`Discovered` 회전 정점[결과 재투입] 거절에서 계속 사용).
- **★ n0 실측 — 전부 DNA-safe(silent-wrong 0)**: (1) **횡단**(같은-iso 회전) 4종 전부 정확·회전-불변 — corner-bite
  Cut 2.776/Fuse 3.700·**sever 2솔리드 Σ0.06**(3c-vi `component_is_outward_tol`이 회전에서 처음 라이브·양 조각 outward
  정확 판정)·containment cavity 1개 2.488, 모두 부피 1e-9 일치·validate clean. (2) **mixed**(한쪽만 회전)도 열림·
  `plane_def`가 None/Some 처리 — 회전 c2 안의 축정렬 d Cut → cavity 1개·vol 23 정확. (3) **회전 공면-접촉**(공면
  밀레스톤 스코프)은 **solved-or-honest-reject**: Z축 회전은 z-접촉면 법선이 정확히 ±Z 유지→감지기가 exact로 잡아 정상
  처리; X축 회전(접촉면 기울임)은 감지기가 놓쳐 general_boolean 도달→boss Fuse는 정확히 풀리고(vol 1.25) pocket Cut은
  정직 거절(Unsupported). **어느 경우도 그럴듯한-오답 없음** — 반올림 계수가 놓치는 공면을 정점-on-평면 TIP 술어
  (`edge_crosses_face`의 `t_plane_side` declare-0)가 base 프레임 exact로 잡기 때문.
- **★ 테스트**: 회전-불변 4종(corner Cut/Fuse·sever·containment)·공면-회전 no-silent-wrong 1종(boss/pocket X축:
  `Err` 또는 정확-clean만 허용). 기존 `transform_rotate_cuboid_tilts_and_blocks_boolean`·`rerotate_same_axis_chains`의
  `ROTATED_UNSUPPORTED` 거절 단언 3곳은 무효화 → **mixed Cut 성공+validate**로 전환(이름 "blocks_boolean" 제거·mixed
  커버리지 보존). 비회전 전 코퍼스 무회귀(가드는 회전에만 발화했으므로 bit-identical).
- **스코프·후속**: fresh-primitive 회전만. near-degenerate/다른-iso 우연공면·회전 결과 재투입(`Discovered`)은 이후.
  **3d-ii** OCCT diff(회전 두 피연산자 `occt_boolean_of`·`#[ignore]`)·**3d-iii** 적대적 회전 코퍼스·**3d-iv** main 병합.

**회전 부울 OCCT 외부 확증 (오버홀 3d-ii).** 3d-i의 회전-불변은 "회전 전후 자기 일관성"이라 같은-부피-다른-위상 버그를
놓칠 여지가 있다. 이 셀은 회전 두 피연산자의 부울을 **OpenCASCADE(OCCT)와 diff**해 독립 커널로 부피·면적을 검산. nacre가
회전 STEP을 export→OCCT가 `bcut/bfuse/bcommon`(rigid라 회전 STEP서 정답)→`approx`(≤1e-6 abs/1e-4 rel) 대조. 프로덕션
변경 0(테스트 전용·`#[ignore]`·기본 스위트 무영향, `-- --ignored`로 실행).
- **★ 5 테스트 전부 통과**(DRAWEXE 로컬 실행): overlap Cut(**Z→X 연쇄 full-tilt**·모든 면 무리수)·Fuse·Common(Z30°),
  **sever Cut(full-tilt)** → 2 솔리드 Σvol/area = OCCT COMPOUND(3c-vi 회전 `is_shell_outward` outward 판정을 외부 truth로
  확증·같은-부피-다른-위상 버그가 숨을 곳), containment Cut → cavity 1개 = OCCT BREP_WITH_VOIDS. `rotated_boolean_matches_occt`
  헬퍼(iso 연쇄 양쪽 적용·OCCT baseline·boolean·Σvol/area approx·결과 반환)로 축약.
- **의의**: 회전 부울이 성숙 커널과 부피·면적 일치 → 3d-i 라이브가 "회전 불변"만이 아니라 **실제 옳은 결과**임을 외부
  검증. 후속: 3d-iii(적대적 스트레스)·3d-iv(main 병합).
