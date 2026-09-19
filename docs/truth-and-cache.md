# 진실과 캐시 — 점·평면·스케치의 최종 타입 구조


아레나는 append-only 이므로 **진실만 들어간다.** 캐시는 곁표에 두어 언제든 버리고 재생한다.
정의가 불변이므로 **캐시는 낡을 수 없다** — 무효화라는 개념이 없고, «버리고 재생»만 있다.

---

## 숫자 규칙 — 일곱 개

1. **진실은 «유리수»이거나 «핸들»이다.** `Rat`(i128)에 들어가는 값은 적고, 안 들어가는 값
   (발견된 좌표 160–480비트, 계수 곱)은 저장하지 않고 **가리킨다**. 가리키는 사슬의 끝은
   언제나 수다(C2·C3 의 동시 만족).
2. **화살표는 한 방향뿐이다.** `Rat → f64`(캐시), `Rat → BigFloat`(판정 상승). `f64 → 진실`은
   없다 — f64 를 들어올려 진실로 삼는 순간 반올림이 진실에 구워진다. 유일한 입구는
   `Rat::from_decimal`(들어올리기가 아니라 사용자가 *쓴* 십진수의 결정적 정규형)이다.
3. **무리수는 값이 아니라 정의로 존재한다.** 회전은 `Angle`(유리수 degree), 기울어진 프레임은
   `Motion::Frame`(법선을 이름으로)이 들고, cos/sin/√ 는 실현 시점에 정확 반올림으로 나온다.
   **모션은 면이 든다** — 정점·모서리는 면을 가리켜 정의상 따라온다.
4. **실현 통로는 하나다.** `realize(정의, bits)` — f64 캐시(표시·tess), STEP 출력
   (`realize(128) → round_to_f64` 정확 반올림), 판정 상승이 전부 같은 통로를 탄다. 캐시는
   언제나 **값 + tol 쌍**으로 함께 만들어지고 함께 버려진다. ★ **f64 좌표는 즉석에서 다시 풀지
   않는다 — 표에서 읽는다**(실측: 판정마다 정의를 다시 만들던 것이 시간의 77%, 25초→0.4초).
5. **판정은 3단이고, 결과는 3갈래다.** f64 필터 → (세 평면이 모션을 공유하면) 정수/`Expansion`
   정확 경로 → astro-float 상승. 결과는:
   - 부호가 **증명되면** → 그 부호(`Sign`).
   - 부호는 못 갈랐지만 **«일치 정밀도보다 가깝다»가 증명되면** → 일치로 처리하되 근거를 실어
     **보고**한다(`Coincident{within}`, `boolean_with_report`). 일치 정밀도
     `same_within = 모델 크기 × 2⁻¹⁸⁰`(f64 출력 해상도 2⁻⁵² 보다 두 워드 아래)은 **모델에서
     유도**되는 값이지 사용자 설정이 아니다 — 낮게 잡으면 비트만 더 쓰고 높게 잡으면 진짜로
     떨어진 것을 합치는 비대칭이 아래로 밀며, 그 아래의 갈라짐은 어떤 출력에도 살아남지 못한다.
     전역 tolerance(*"가까우면 붙여라"*)와 방향이 반대(*"가깝다고 **증명되어야** 일치"*)다.
   - 둘 다 증명 못 하면 → **이름 붙은 거절**(`JudgeExhausted`·`DegenerateWitness`·
     `PrecisionBudget`). 조용한 0은 없다(C7).

   정밀도는 상수가 아니라 모델이 정하고, 필요한 비트를 계산해 **한 번에** 점프한다.
6. **동일성은 정준 이름의 `==` 다 — 이름의 그릇은 임의정밀이다.** 이름은 하나만 저장한다:
   `PlaneName = Narrow([Rat;4]) | Wide(BigInt)`, 들어가면 반드시 `Narrow`(정규화 불변식).
   동일성은 enum 전체의 `==` 이고, 산술·프레임·지름길은 `narrow()` 투영을 빌려 읽는다
   (`plane_name_exact`가 중간 계산을 이미 임의정밀로 닫았으므로 남은 것은 그릇뿐이었다).
   좁은 형태(법선·`n·n`)가 없다고 **진실을 거절하지 않는다** — 필요하면 점 셋에서 임의정밀로
   실현한다. ★ 이 규칙의 «항상 이름을 갖는다»는 `Known` 평면(유리수 닫힘)에만 적용된다 —
   무리수 모션이 낀 `Through` 평면은 정확 계수가 무리수라 이름이 없고, interning 없이 술어가
   매번 답한다(느릴 뿐, 틀리지 않는다).
   저장은 **진술 키**(정렬 삼중항+모션, `surface_through_ids`)로 intern 되어
   «같은 진술 = 한 핸들» 은 지켜지고, 기하 동일성은 예고대로 술어의 몫이다.
7. **`Rat` 을 넓히지 않는다(C8).** 넓어지는 것은 중간값(BigInt 이름 유도·`Expansion` 술어·
   BigFloat 상승)뿐이다.

---

## 최종 타입 — 진실

```rust
// ─── nacre-geom / nacre-topo ── 진실 (아레나, append-only) ─────────────

/// 이 enum 이 `Surface` 이고 아레나
/// (`surfaces: Store<Surface>`)가 그것을 든다. 실현은 `surface_cache: Vec<SurfaceCache>`.
/// ⇒ **맨이름 `Surface` = 이것**이고, geom 의 f64 실현은 `nacre_geom::Surface` 로 적는다
/// (열린 항목 8 의 철자 규칙).
pub enum Surface {
    Plane {
        points: PlanePoints,
        motion: Option<Handle<MotionNode>>,
    },
    /// 성분형이 옳은 이유·`ref_dir` 원시 규칙은 design.md
    /// 원통 절이 상세히 적는다.
    Cylinder {
        def: CylinderDef,                     // { origin, dir, ref_dir } Rat + r2: BigRat — 반지름은 «제곱», 폭은 BigInt(2026-09-15)
        motion: Option<Handle<MotionNode>>,
    },
}

/// 평면의 세 점 — 진술의 종류가 곧 변종이다.
pub enum PlanePoints {
    /// 구성 평면(벽·캡·사용자 진술 평면): **모션 전 프레임**의 유리수 세 점.
    /// base case — 여기서 끝난다(C2). 오늘의 `surface_points` 값 그대로.
    Known([[Rat; 3]; 3]),
    /// datum 평면: 모델 정점 셋을 지난다. 발견된 좌표는 Rat 에 안 들어갈 수 **있으므로**
    /// (타입이 그것을 허용한다 — 실측 코퍼스는 안 때린다) **가리킨다**(C3). 필요할 때 실현한다.
    /// ★ 생산자는 `DatumDef::ThroughVertices`. 능력의 핵심은 폭이 아니라
    /// **좌표로는 그 평면을 말할 수 없다**는 것(실측 220/220).
    /// ★ **핸들은 정렬해 저장한다** — 같은 셋 = 같은 진술. 이름 없는(무리수 모션) datum 은
    /// interning 이 못 받쳐 주므로, 순서가 다르다고 같은 평면이 두 핸들이 되는 것을 구성
    /// 시점 정규화로 막는다. 법선 방향은 점 순서가 아니라 정준 부호 규약과
    /// `Face::orientation`/`flip` 이 들므로 잃는 정보가 없다.
    Through([Handle<Vertex>; 3]),
}

/// 정점 = 자기 **정의**, 그리고 정의가 곧 타입이다 — `Surface` 와 같은 모양. 좌표는 진실이 아니라
/// 캐시다(PointCache).
/// ★ 단일형이 아니라 **세 변종**이다: M5 에 원통 seam 정점이
/// 실재하고 그 점은 «세 평면»으로 적을 수 없다. 변종별 불변식(Q5)이 이 구조의 근거다.
pub enum Vertex {
    /// 세 평면의 교점 — 이름이 곧 점. (D != 0 은 좌표 재생이 생기는 자리에서 단언한다.)
    ThreePlane([Handle<Surface>; 3]),
    /// 두 곡면의 교차 «곡선» 위의 점 — M3 원통 seam(테두리 원의 θ=0). 점을 못 박는 매개
    /// 정보(원통의 `ref_dir`)는 **M6 의 원통 진실과 함께 왔다**(`CylinderDef`) — 그래서
    /// `OnSeam([cylinder, cap])` 은 「rim ∩ +`ref_dir` 광선」으로 **정확히 지정된 한 점**이다.
    /// ✔ 정의는 완성됐고 **재생하는 기계도 섰다**(칸 ㊵ `realize_vertex`, `seam_point`). ⏳ 없는 것은
    /// 그 값을 캐시에 **되쓰는 것**(열린 항목 28) — 그래서 STEP 은 오늘도 만들 때 나온 f64 를 내보낸다.
    OnSeam([Handle<Surface>; 2]),
    /// 이 점은 도법기하의 **관통점**(piercing point)이다.
    /// 담체 종류를 **구조로** 말한다(닮은 핸들 셋이 아니라). `root` 는 정준 선 방향
    /// (`n₀ × n₁`, 저장 순서)을 따른 오름차순이고 접선은 `QuadRoot::Double` 이다.
    Pierce {
        planes: [Handle<Surface>; 2],         // 오름차순 핸들
        cylinder: Handle<Surface>,
        root: QuadRoot,                       // Lo | Hi | Double
    },
    // 원뿔 꼭짓점 `Apex(Handle)` 등은 아직 예고다.
}

/// 곡면 집합은 «담체»(어디 위에 있나)를 정하고, 경계가 «어느 조각»인지를 정한다.
/// **곡면 집합만으로 남김없이 정해지는 것은 정점뿐이다** — 차원 0, 평면 셋이 자유도 셋을 다 먹는다.
/// 간선(1차원)은 담체 위 «어느 조각»인지 끝점 둘이, 면(2차원)은 «어느 영역»인지 루프가 더 정한다.
/// 그래서 정점의 진실은 핸들 셋으로 끝나고, 좌표는 어느 실체의 진실에도 없다(전부 캐시).
/// ★ 간선의 `vertices` 는 «집합»이 아니라 «순서»다 — 원 담체에서 `[A, B]` 는 A→B 를 축에 대해
/// 반시계로 가는 호이고 `[B, A]` 는 그 나머지다(M6-2b; `surfaces` 는 오름차순 정렬한 집합인 것과 대비).
pub struct Edge {
    pub surfaces: [Handle<Surface>; 2],       // 담체 (두 끝점 면집합의 교집합으로 파생 불가)
    pub vertices: [Handle<Vertex>; 2],        // 경계
}
// Store<Curve> 는 없다 — 곡선은 진실이 아니라 EdgeCache 다.

// ─── 모션 (nacre-topo — Handle 이 필요하다) ────────────────────────────

pub enum Motion {
    Rotate    { axis: Axis, pivot: [Rat; 3], angle: Angle },
    Translate { offset: [Rat; 3] },
    Mirror    { axis: Axis, offset: Rat },    // det = −1, 사슬 패리티
    /// 평면 자신의 프레임으로의 기저 변경 — 기울어진 면 위 스케치를 정확하게 만든다.
    /// ★ 법선을 성분이 아니라 **이름으로** 든다 — 재귀는 프레임 없는 평면(세계)에서 끝난다.
    Frame     { plane: Handle<Surface>, placement: FramePlacement, flip: bool },
}

/// 프레임의 원점과 u 방향 — 사용자가 정하면 값, 아니면 정준 유도. `PlanePoints` 와 같은 이분법.
pub enum FramePlacement {
    /// **기본값.** 평면의 순수 함수로 유도되는 정준 프레임 — 원점 = 세계 원점의 수선의 발
    /// (`p = (−d/n·n)·n`), 축 = **Arbitrary Axis**(DXF/AutoCAD 규약: `u = ẑ×n`, 법선이
    /// 정확히 수직이면 `ŷ×n` — 갈래를 정확히 가른다). 실현 시점에 필요한 정밀도로 계산되고
    /// 아무것도 저장하지 않는다 ⇒ 정준값이 `Rat` 을 넘치는 평면(분모 n·n 제곱, 코퍼스 1.6%)도
    /// **Through 평면**(거대 계수 위엔 i128 에 드는 유리수 점이 일반적으로 없다)도 같은 규약을
    /// 그대로 받는다 — 원점 위치가 오버플로 여부와 무관하게 사용자 기대대로다.
    /// 한 평면 위의 스케치들이 자동으로 한 노드를 공유한다.
    /// ★★ 유도 규약은 스펙으로 **동결**한다 — 바뀌면 기존 스케치가 조용히 돈다.
    /// ★ 구현: 좁으면 유리수 `PlaneFrame`(비트 보존), 넘치면 판정층의
    /// `MoveNode::FrameWide`(BigInt 쌍둥이 — 넘침이 존재하지 않는 실현).
    Canonical,
    /// 호출자가 **명시적으로** 이름 붙인 값 — 정준 규약과 다른 프레임을 원할 때만
    /// (`world_zx` 의 `+u = ẑ` 처럼 유도값과 다른 규약, `through_points`·`with_origin`).
    /// **그 평면의 `points` 가 적힌 좌표계**의 3D 유리수다(모션 없으면 세계. 2D 가 아닌 이유:
    /// 평면 위 2D 좌표는 무리수인 u·v 축을 전제한다 — Frame 이 정의하려는 바로 그것).
    /// ★ origin 은 평면 **위**의 점(정확 검사 = 방정식 대입, 성분 모양이 아니다 — C1),
    /// ref_dir 은 법선과 평행하지만 않으면 되고 실현이 평면으로 사영한다.
    /// 저장된 값이 진실이라 replay 는 저장이 보장한다. Known 전용.
    /// ★★ 검사 실패는 **구성 시점의 이름 붙은 거절**이다(`OriginNotOnPlane`·
    /// `RefDirParallelToNormal`·십진 창 밖) — 조용히 Canonical 로 대체하지 않는다(사용자가
    /// 말한 곳과 다른 곳에 앉히는 것이 곧 조용히 틀림이다). Through 평면에 Named 를 주면
    /// 유리수 점이 그 위에 정확히 놓일 수 없어 같은 검사에 자연히 걸린다 ⇒ 안내는 «인자 없는
    /// 기본(Canonical)을 쓰라». 특정 정점에 원점을 앉히는 요구가 증명되면 값이 아니라
    /// **가리키는** 변종(`origin: Handle<Vertex>`)을 그때 추가한다(규칙 1 의 연장).
    Named { origin: [Rat; 3], ref_dir: [Rat; 3] },
}

pub struct MotionNode {
    pub motion: Motion,
    pub parent: Option<Handle<MotionNode>>,   // 숲 — 이력의 꼬리를 공유. interned.
}
```

### 정점 — interning 하지 않는다

인접 면 집합은 **자란다**(나중 연산이 그 점을 지나는 면을 더 만든다) — 어떤 키를 잡아도 같은
점의 키가 바뀌므로 점의 동일성은 인접으로 키잡을 수 없다. 저장은 **정하는 세 면**만 들고(나머지
인접은 위상이 안다 — 적으면 두 번째 설명), 점의 동일성은 기하 질문으로 남는다(`merge_coincident`).
*"같으면 같은 핸들"* 이 구성 시점에 일어나는 것은 **평면뿐**이다. 모서리·면도 interning 하지
않는다(두 면이 여러 토막에서 만나고, 한 평면 위에 떨어진 면이 둘일 수 있다 — 곡면 집합만으로는
«어느 것»을 못 말한다).

`validate` 의 일: 정의하는 세 면 위는 구성상 자명 — 남는 진짜 질문은 **정의하지 않는 나머지
인접 면 위에도 있나**(4-평면 동시성이 그 자리)다.

★ **저장은 «정하는 셋»이되, 이름은 인접 전부에서 정준으로 고른다** (2026-09-06, 칸 ⑪). 불리언 안에서
점을 부르는 이름은 그 점을 지나는 평면 집합의 함수 하나(`canonical_triple` — 사전식 최소의 독립 삼중)
이고, 피연산자 정점은 인접 면의 클래스 전부(위상이 안다)로 그 함수를 부른다. 저장된 정의 셋은 그중
하나의 철자일 뿐 다음 불리언의 이름을 정하지 않는다 — 정준은 불리언마다의 클래스 번호 공간에서 고르는
값이라 저장할 수 있는 것이 아니다.

---

## 스케치

```rust
// ─── 프로파일 — 프레임을 모르는 순수 2D. 구성 시점에 Rat 으로 확정된다. ──

pub struct Profile2d {                        // 한 재료 영역 — 현행 의미 그대로
    outer: Ring2d,
    holes: Vec<Ring2d>,
}
pub struct Ring2d {                           // 정점 + 변. edges[i] = vertices[i] → vertices[i+1 mod n]
    vertices: Vec<[Rat; 2]>,                  //   온전한 원 = 정점 1 + Arc 1 (정점 = 솔기)
    edges: Vec<Edge2d>,
}
pub enum Edge2d {                             // 조각 «하나» — 순서 있는 고리라 시작=앞 꼭짓점 (geom::mixed)
    Line,
    Arc { center: [Rat; 2], r2: Rat, ccw: bool },      // ✔ r² (2026-09-15, 항목 25 · 1단계); 중심의 정의화는 2·3단계
}
// ★ **커널의 문은 «고리»다, «펜»이 아니다** — `Ring2d::new(vertices, edges)`(검증: 짝 맞음·
//    영길이·호의 r²·다음 정점이 원 위)·`Ring2d::circle`·`Ring2d::polygon_decimal` + 단계 문 `arc_turns(_rat)`
//    (변, 끝점)·`arc_to_rat`, 그리고 `from_paths(Vec<Ring2d>)`. 펜은 kit 의 것이고, 제약 해석기·외부 데이터가
//    와도 같은 문으로 온다 — 커널은 «어떻게 그렸나»를 모른다.
// 링 더미 → 짝수 깊이 = 재료(even-odd) → 섬마다 Profile2d 하나 — from_rings/from_paths, 현행 유지.

// ─── 배치 — 어떤 평면 위 + 배치(기본은 정준 유도, 명시하면 값). ──

pub struct SketchFrame {
    pub plane: Handle<Surface>,
    pub placement: FramePlacement,            // 기본 Canonical(유도) / 명시 Named(값, Known 전용)
    pub flip: bool,                           // ★ 사용자 몫이 아니다 — 연산이 실현된 ŵ 를
}                                             //   면의 바깥 법선과 내적해 **측정**으로 정한다
// world_xy()/yz()/zx() 는 이 struct 를 채우는 설탕이다 — 타입 변종이 아니다.
// (xy·yz 는 Canonical 과 일치, zx 만 +u = ẑ 가 유도값 −x̂ 과 달라 Named 를 쓴다.)
```

| 규칙 | |
|---|---|
| **호의 진실은 끝점이다, 각도가 아니다** (칸 ⑨) | 3D 경계 표현이 요구하는 것은 **정점과 원**이고 각도는 어디에도 안 쓰인다(`Arc` 에 각도 필드 없음 — 실측). 끝점 표현이 더 넓다(무리수 도의 호도 끝점이 유리수면 정확 — 3-4-5). ✔ **임의 각은 이미 열렸다**: `arc_rat(center, start, end, ccw)` 이 끝점을 받고 «원 위인가»만 검사하지 각도를 제한하지 않는다; `arc_turns` 는 90° 배수 끝점을 `(x,y)→(−y,x)` 로 유도해 주는 **설탕**일 뿐이다. ⚠ **남은 제한은 «스케치»가 아니라 «돌출»에 있다**: `refuse_non_quarter_arcs`(`ops.rs`)가 프리즘 빌더에서 비-사분 호를 `ArcSweepNotQuarterTurn` 으로 거절한다 — 원통 옆면의 정확 표현이 사분 가족에서만 서기 때문이고(M6), 스케치 입력의 한계가 아니다. ⏳ 그 벽은 열린 항목 25(중심·반지름의 정의 기반화)와 곡선 마일스톤이 함께 연다 |
| **원통 원시체는 없다 — 원은 스케치다** (칸 ⑨) | `Operation::Cylinder` 는 원 프로파일 Extrude 와 위치 정준 비트 동일이 실측된 뒤 은퇴. 직선–호 꼭짓점 = `Vertex::Pierce`(두 평면의 **저장된 정준 이름** × 원통, 근은 매개 값으로 고른다 — 빌더가 자기 점으로 계수를 다시 만들면 `Lo/Hi` 가 뒤집힐 수 있다), 온전한 원 = `OnSeam` |
| **경계는 f64, 진실은 구성 시점에** ✔S3 | 공개 API 는 f64 그대로(§design 6.0). `Rat::from_decimal` 왕복을 `Profile2d` **생성자**에서 한다 — `check()`(자기교차·중첩·포함)가 진실 위에서 정확 술어로 돌고(`orient2d_rat`: narrow 우선 → BigInt 전역 부호, geom 의 `_rat` 워커 쌍둥이), 십진 창(1e38/1e-22) 밖 치수가 구성 시점의 이름 붙은 에러다(`ProfileOutsideDecimalWindow` / sketch 층 `OutsideDecimalWindow`). ★ f64 부호 ≠ 십진 부호가 실측 사실이라(0.1·0.2·0.3 공선이 이진에선 굽음) check 를 진실 위로 옮긴 것이 정확성 변경이고, 그 방향은 항상 "작성자가 쓴 수가 이긴다" |
| **공선 중간점은 생성자가 지운다** ✔S3 | 공선 정점의 양옆 벽은 한 평면 → 교차가 직선이라 정의 불가, 그리고 비-2-manifold 퇴화다. 제거는 형상 불변(무손실 정규화, 거절 아님 — 정리된 프로파일의 프리즘은 깨끗한 쌍둥이와 비트 동일, 모든 코너가 3-평면 정의 보유 = Q2 ② 닫힘). 판정은 Rat 위 정확 orient2d, 제거 조건은 **엄격 내부**(중복점·스파이크는 생존해 각자의 이름 붙은 에러로 보고된다 — 경계 포함 판정을 재사용하면 작성자의 실수를 조용히 지운다) |
| **프레임은 평면에서 뜨지 않는다** ✔S9 | `Named` 의 `origin` 은 참조 평면 **위**의 점(정확 검사 — C1). ★ 구현: `SketchFrame::named` 가 구성 시점에 검사한다 — `plane_residual_sign`(scalar, **전역**: Narrow 는 Rat 대입 → 넘치면 BigInt, Wide 는 BigInt — fail-open 없음) ≠ 0 이면 `OriginNotOnPlane`, `WideFrame::named_of` None 이면 `RefDirParallelToNormal`(그 함수는 폭에 전역이라 None 은 평행/영벡터뿐), 십진 창 밖은 `FrameOutsideDecimalWindow`. ★ 좌표계에 주의: «정의하려는 프레임의 (u,v,w)» 가 아니라 **그 평면의 `points` 가 적힌 좌표계**의 3D 점이다(모션 없으면 세계 — 상자 윗면 z=1 이면 origin 은 `[1,1,1]` 같은 점이지 셋째 성분 0 이 아니다). "평면 위" 는 성분이 아니라 **방정식 대입**(`a·x+b·y+c·z+d = 0`, 정확)으로 검사한다. `Canonical` 은 유도라 검사할 것이 없다. 면 위 스케치의 밑캡이 대상 면의 surface 핸들을 공유하는(flush 접촉 = 핸들 비교) 전제이기도 하다. 평면에서 d 떨어진 스케치가 필요하면 origin 을 띄우는 것이 아니라 **오프셋 평면**을 만든다 — 같은 프레임 안 `(0,0,d),(1,0,d),(0,1,d)` 유리수 세 점의 `Known` 평면(datum 가족, S5). 담체가 실제 평면 핸들로 남아 이름·interning·flush 규칙이 그대로 성립한다 |
| **정확 유리수 프레임은 노드를 만들지 않는다** ✔S9 잠금 | 평면에 모션이 없고 실현된 원점·축이 **정확 유리수 직교**로 증명되면(오늘의 `RatFrame` 게이트 — 세계 축 평면 셋이 전부 여기 든다, `Canonical` 이든 `world_zx` 의 `Named` 든) `Motion::Frame` 노드를 생성하지 않고 세계 유리수 산술로 구성한다(모션 = None) — 빈 사슬 지름길(`shared_base` 상쇄·축별 tol 0·정수 Shewchuk)이 그대로 산다. 미러의 *"정확한 f64 에 남으면 노드를 안 만든다"* 와 같은 구성 시점 정규화다. ★ 지반 잠금(S9): 씨앗 평면의 canonical 프레임이 단위축 위에 **비트 정확**으로 실현됨(`a_seeded_planes_canonical_frame_is_the_world_basis_exactly`) — 생략은 중복 제거지 값 손실이 아니다. 오늘의 게이트 표현식은 `exact()`(★★ 주석 동결 — 바꾸면 노드 인구가 움직인다); `PlaneFrame`+`inv_sqrt_exact` 기반의 더 강한 게이트는 §열린 항목 11 |
| **세계 축 평면 셋은 `Model::new()` 가 심는다** ✔S9 | 핸들이 결정적(0·1·2)이라 replay 가 자명하다. 캐시 방향은 **−축**(밑캡의 감각 — S9 행), `Default` 는 `new()` 로 위임(무씨앗 뒷문 없음). 씨앗은 정의상 영구 orphan — 정의·프레임만 가리키는 평면을 순회·직렬화가 따라가야 한다(§열린 항목 4) |
| **스케치는 모델에 저장되지 않는다** | 연산 로그가 스케치의 진실이고, 모델에 남는 것은 벽·캡 평면(`Known` 점 + 프레임 모션)이다. 치수의 흐름: 변 `(x1,y1)→(x2,y2)` + 스윕 → 벽의 세 점 `(x1,y1,0),(x2,y2,0),(x1,y1,d)` — 산술 없이 유리수 그대로. 계수(이름)에만 곱이 있고, 그것이 진실이 점이어야 하는 이유다(계수 `i128` 적합 25.8% vs 점 100%) |
| **기본 프레임은 정준 유도(`Canonical`)다** | 원점 = 세계 원점의 수선의 발(부호·스케일 불변), 축 = Arbitrary Axis(DXF) — 사용자가 아무것도 안 정하면 이 규약이고, 규약은 스펙+테스트로 동결한다. 사용자가 명시한 프레임만 `Named` 로 **값을 저장**한다(유도값과 다른 규약을 유도로 흉내내면 규칙 변경에 조용히 돈다) |

연쇄가 닫힌다: 평면에 스케치 → 솔리드 → 그 면 위에 스케치 → … 매 단계가
`SketchFrame{plane: 이전 단계의 평면}` 으로 **사용자 조작당 1단**만 깊어지고(C6 — 불리언은
평면을 만들지 않는다, `a_boolean_mints_no_surface`), 좌표는 끝까지 유리수이며, 무리수는 프레임
실현(`1/√`)에만 산다.

---

## 이름과 interning

이름은 **하나만 저장**하고, 좁은 형태는 저장이 아니라 **투영**이다:

```
Surface (진실 — 평면이면 점 셋, 원통이면 def, + 모션)
   │ 유도 (정준화: 분모 털기 → gcd → 부호 규약)  ※ 평면이고 유리수 닫힘일 때만
   ▼
PlaneName = Narrow([Rat;4]) | Wide([BigInt;4])    ← 저장은 이것 하나
   ├─ narrow() → Option<&[Rat;4]>    산술·프레임·지름길 — 사본이 아니라 빌려 읽는다
   └─ + 모션  → SurfaceKey::Name     interning 표(`surface_ids`)의 키 — «합친다»

Surface (진실) ─── 이름이 유도 «안 될» 때 ──→ SurfaceKey::Verbatim(그 진실 그대로)
   무리수 모션의 Through 평면 · 4계수형이 없는 모든 곡면 — «안 합친다»
```

```rust
/// 정준 계수 — 진실(점 셋)에서 유도한 «이름». 그릇이 임의정밀이라 유리수 닫힘인 평면은
/// **항상** 이름을 갖는다. 역할 둘: 동일성(`==`, enum 전체) · 산술 지름길(`narrow()`).
///
/// ★ 정규화 불변식: **i128 에 들어가는 값은 반드시 `Narrow` 로 저장된다** — 유일한
/// 생성자(점에서 유도)가 강제한다. 같은 값 = 같은 변종 = 같은 비트라 `Eq`/`Hash` 가
/// 구조적으로 성립한다("정준 여부는 값이 말한다"와 같은 원칙). ~99% 가 Narrow(인라인,
/// 힙 없음)이고 Wide 는 ~1%.
pub enum PlaneName {
    Narrow([Rat; 4]),
    Wide([BigInt; 4]),   // Box 없음 — BigInt≈32B 라 [BigInt;4]=128B=[Rat;4], enum 크기가
                         // 같아 Box 는 할당+간접만 더한다 (S2 구현에서 실측)
}
impl PlaneName {
    /// Shewchuk 정확 술어(`Expansion` 조각)·프레임 유도·지름길이 읽는다.
    /// `None`(Wide)이면 그 지름길만 못 타고 판정이 일반 경로로 — 느릴 뿐 틀리지 않는다.
    pub fn narrow(&self) -> Option<&[Rat; 4]> { … }
}

/// interning 표의 **키 — 팔이 둘이고 영원히 둘이다** (2026-09-12 확정; 오늘 코드는 아래 ⏳).
///
/// 같은 계수라도 Constructed(세계)와 Moved(모션 전 프레임)는 다른 평면이므로 **모션이 늘 함께**
/// 열쇠에 든다.
pub enum SurfaceKey {
    /// **유도된 정준형** — 다르게 진술해도 같은 평면이면 한 핸들(**기하 동일성**).
    Name(PlaneName, Option<Handle<MotionNode>>),
    /// **진실 그대로** — 글자가 같아야 한 핸들(**문자 동일성**). 이름이 «없는» 모든 경우가
    /// 여기로 온다: 무리수 모션의 `Through` 평면, 그리고 4계수 음함수형이 없는 모든 곡면.
    Verbatim(Surface),
}
```

### ★★★★ 열쇠가 «둘»인 이유 — 그리고 구·원뿔이 와도 안 늘어나는 이유

오늘 코드는 표가 **셋**이다(`surface_ids`·`surface_through_ids`·`cylinder_ids`). 그런데 재 보면
뒤의 둘이 **같은 관계**다 — 열쇠 필드가 **진실의 필드와 같다**:

| 오늘의 표 | 열쇠 | 실제로 무엇인가 |
|---|---|---|
| `surface_ids` | `(PlaneName, motion)` | **유도된 것** — 합친다 |
| `surface_through_ids` | `([Handle<Vertex>;3], motion)` | `Surface::Plane{points: Through(vs), motion}` **그대로** |
| `cylinder_ids` | `(CylinderDef, motion)` | `Surface::Cylinder{def, motion}` **그대로** |

⇒ 둘은 «어느 종류냐»로 갈려 있었을 뿐 관계는 하나다. **진실 자체를 열쇠로** 쓰면 `Verbatim` 한 팔이
둘을 삼키고, **구·원뿔·토러스·NURBS 가 와도 새 팔이 필요 없다** — `Surface::Sphere{..}` 가 생기면
`Verbatim` 이 그냥 받는다. 팔 이름이 «어느 종류»가 아니라 **«합치는가 아닌가»** 를 말하기 때문이다.
☑ 그래서 원통의 **일부러 약한** 보장이 이름에 드러난다: `ref_dir` 이 다르면 seam(θ=0 이음매)이
갈라지므로 기하가 같아도 **합치면 안 된다** — 기하 동일성은 술어가 물을 때마다 답한다(규칙 6).

### 열쇠는 «생산자»가 아니라 «진실»이 고른다

```rust
surface_ids: HashMap<SurfaceKey, Handle<Surface>>     // 표 «하나»

fn surface_key(&self, truth: &Surface) -> SurfaceKey {
    match self.derive_name(truth) {
        Some(n) => SurfaceKey::Name(n, truth.motion()),
        None    => SurfaceKey::Verbatim(truth.clone()),
    }
}
```

오늘은 `push_plane` → 표 ①, `push_plane_through` → 이름 있으면 ① 없으면 ②, `push_cylinder` → ③
으로 **어느 문으로 들어왔는지가** 열쇠를 고른다. 위 모양은 **사설 깔때기**에서 진실을 보고
고르므로, 「이름이 있으면 이름으로, 없으면 그대로」가 **모든 곡면 종류에 자동으로** 적용된다.

★★ **`derive_name` 이 원통에 `None` 을 주는 것은 «없어서»가 아니라 «결정»이다.** 원통에 4계수형이
없는 것은 사실이지만, 그 함수가 지키는 것은 정책이다 — 누군가 원통의 정준형을 만들어 `Some` 을 주기
시작하면 기하가 같은 원통이 **조용히 합쳐지고 seam 이 갈라진다.** 그래서 그 팔은 `None` 을 **명시적으로**
돌려주고 이유를 그 자리에 적는다; 새 종류(구·원뿔)도 «합칠 것인가»를 먼저 정하고 나서야 이름을 얻는다.

★★ **`surface_ids` 는 캐시다** — 아레나를 돌며 진실마다 `surface_key()` 를 다시 매기면 통째로 재생된다.
그래서 `rebuild_surface_names()`(열린 항목 36)의 형제로 «표도 버리고 재생 → 동일» 잠금이 서고,
그 잠금이 `SurfaceCache::Plane.name == 표의 열쇠` 일치까지 증명한다. ⚠ **프로덕션에서만** 그렇다: test-only
`push_plane_unregistered` 는 표를 **건너뛰므로**(실측: `insert` 0건 — 「한 평면을 두 핸들로」 픽스처가
그것을 필요로 한다) test-util 아래에서는 재유도가 오늘 표와 다르다. 그 잠금은 프로덕션 생산자만 거친
모델에서 돈다.

⚠★★★★ **오늘 이 불변식을 «검사하는» 코드가 없다**(실측 2026-09-12): `validate` 에 중복 곡면 검사
**0건**이고, 지키는 것은 생산자 넷(`intern_plane`·`push_plane_through`·`push_cylinder`·
test-only `push_plane_unregistered`)이 **각자 표를 기억하는 규율**뿐이다. ⇒ 다섯째 생산자가 잊으면
컴파일도 되고 테스트도 초록인데 **조용히 중복 핸들**이 생긴다. 문을 하나로 모으면 그 실수가
**구조적으로 불가능**해진다 — 이것이 이 개편이 사는 것이다.
☑ 문을 넘나드는 경우는 이미 잠겨 있다:
`a_through_plane_and_a_known_plane_that_are_one_plane_share_a_handle`.
☑ 인구가 배타적이라는 것도 실측: 유리수 점 셋으로 말한 평면은 계수가 **반드시** 유리수라 항상 이름이
있고(⇒ `Name`), `Verbatim` 으로는 이름이 계산 안 되는 것만 간다 ⇒ **한 평면이 두 팔로 갈라지지 않는다.**

⏳ **실현 조건**(재기 완료): 진실 안에 `f64` 가 **0개**(전부 유리수·핸들)라 `Eq`/`Hash` 파생이
가능하다. ⚠ 다만 오늘 `Surface`·`PlanePoints` 는 `PartialEq` 만 파생하므로 **`Eq, Hash` 를 더해야**
한다(`CylinderDef`·`MotionNode` 는 이미 있다).
⚠ **대가 하나 — 측정 대상**: `Verbatim(Surface)` 는 진실의 사본을 열쇠로 들고, `Surface` 의 크기는
가장 큰 팔(`PlanePoints::Known`, 유리수 아홉 ≈ 288 B)이 정한다. 그런데 `Known` 은 **항상 이름이
있어서 `Verbatim` 으로 갈 일이 없다** ⇒ 쓰지 않는 팔 때문에 열쇠가 커진다. 오늘도 `cylinder_ids` 가
`CylinderDef` 사본을 들고 `PlaneName::Wide` 는 `[BigInt;4]` 라 비슷한 자릿수이지만, **재 보고**
문제면 `Verbatim(Box<Surface>)` 가 탈출구다.

- **정준화**: 분모 털기(lcm) → 내용(gcd) 나누기 → 부호 규약(첫 0 아닌 성분 양수). 정준 여부는
  플래그가 아니라 **값이 말한다**(정수인가·서로소인가·부호 규약인가).
- **생산자는 점만 진술한다**(`push_surface_with_points`) — 이름은 커널이 유도하므로 *"한 평면을
  두 가지로 진술한다"* 가 표현 불가능하다.
- `==` 는 «같은 평면»을 증명하지만 **`!=` 는 아무것도 증명하지 않는다** — 이름 없는 평면
  (무리수 모션의 `Through`)은 interning 되지 않으므로 기하 술어가 뒤를 받친다. interning 은
  결정적이다(먼저 넣은 쪽이 이긴다).
- **interning 의 값은 속도가 아니라 동일성이다**(실측 3.7%) — 같은 평면이 두 핸들로 갈리지 않는 것.
- **큰 정수 계수의 정확 부호**: `2^53` 은 f64 *하나*의 한계일 뿐 — `i128` 은 `Expansion` 3조각
  (53비트씩)으로 정확하고, 정규화 실패한 유리수는 술어 안에서 분모를 턴다. 실측 분포: 1조각
  71.9% · 2조각 26.1% · 3조각 1.2% · `i128` 초과 0.9%(임의정밀 이름이 받는다). f64 필터가
  먼저 걸러 expansion 은 애매한 ~1%에서만 돈다.

---

## 캐시

```rust
pub struct Model {
    // 진실 — append-only. ★★ **모두 비공개**(최종형): 밖은 §문의 이름의 문으로만 읽는다.
    //   ✔ 칸 56 이 닫았다 — 다섯 저장소 + adj + live_solids 가 비공개가 됐고, 남은 공개 필드는
    //   `surface_name` 하나다(§문의 이름이 `surface_cache(h).name()` 으로 접을 몫).
    surfaces: Store<Surface>,                 // 읽기는 좁은 문(surface/surface_count/…), push 는
    vertices: Store<Vertex>,                  //   생성자 경유만. 전량 순회 없음: 아레나엔
    edges:    Store<Edge>,                    //   superseded 도 있어 소비자는 live face 를 걷는다.
    faces:    Store<Face>,
    shells:   Store<Shell>,
    solids:   Store<Solid>,
    motions:  Store<MotionNode>,              // interned

    // 루트 — 지금 살아있는 솔리드 핸들 목록(도달 집합의 뿌리). 저장소가 아니라 «무엇이 현재 모델인가».
    pub live_solids: Vec<Handle<Solid>>,

    // 캐시 — 핸들 인덱스 병렬, 통째로 버리고 재생 가능.
    // ★ **비공개다**(실측: 오늘의 `vertex_cache`·`edge_cache` 에 `pub` 이 없다). 읽기는 좁은 문
    //   (`vertex_point`/`vertex_tol`/`surface`/`edge_curve`)이 이미 열려 있고, 정제도 `Model` 의
    //   좁은 문으로 들어온다 — 밖에서 벡터를 만지면 index-parallel 불변식을 아무도 못 지킨다.
    vertex_cache:  Vec<PointCache>,
    surface_cache: Vec<SurfaceCache>,         // 실현. 아레나가 진실을 든다
    edge_cache:    Vec<EdgeCache>,            // 평가 가능한 곡선

    // 동일성 — 「같은 곡면 ⇒ 같은 핸들」이 되게(그래야 동일성이 정수 비교다).
    // ⏳ 오늘 코드는 표가 «셋»이다(`surface_ids`·`surface_through_ids`·`cylinder_ids`).
    //    최종형은 표 «하나» — 열쇠가 `SurfaceKey{Name|Verbatim}` 이고 **진실에서 유도**된다
    //    (§이름과 interning). 모션도 같은 패턴이다(`motion_ids`, 이미 한 표).
    surface_ids: HashMap<SurfaceKey, Handle<Surface>>,
    motion_ids:  HashMap<MotionNode, Handle<MotionNode>>,
}
// ⏳ **`surface_name` 은 사라진다** — 이름은 유도되므로 캐시이고, `SurfaceCache::Plane.name` 으로
//    들어간다(위 §수치·캐시 타입). 오늘은 `pub HashMap` 곁표라 곡면의 캐시가 두 벌이다.

// 수치 층 — 「값 + 그 값이 갇힌 오차」가 **세 계단**으로 산다(2026-09-13 정리; ⏳ 열린 항목 23 이 집행).
//   같은 것의 다른 정밀도는 `Hp` 접두사 하나로만 다르다.
pub struct Mag { m: f64, e: i64 }                                  // f64 밖 범위의 보수적 크기(0에서 먼 쪽)

// ── 계단 1: «원자» — 값 하나 + 그 반경. ────────────────────────────────
pub struct Bounded   { value: f64,      error: f64 }
pub struct HpBounded { value: BigFloat, error: Mag }

// ── 계단 2: «경계 지어진 점» — 원자 셋(축별). realize 가 이미 `[HpBounded;3]` 로 든다. ──
//   ★ 값과 오차를 «묶어» 든다(따로 든 배열 둘이 아니라). 그래야 «참값 ∈ value±error» 가 구조로 서고
//     transform 이 값만 옮기고 오차를 안 옮기는 desync 가 **불가능**해진다.
pub enum PointCache {
    Bounded { coord: Point3, bound: [Mag; 3] },                   //   정의에서 실현 — 참값 ∈ coord ± bound
    Ceiling { coord: Point3 },                                    //   싼 도로가 멈췄다 — 비싼 문이 답할 수 있다
    Unrealized { coord: Point3 },                                 //   실현이 뒤에 없다 — 도로가 없거나, 안 물었거나
}
// ★ 연산이 만든 정점은 push 시점에 실현되므로 거의 전부 `Bounded` 다(코퍼스 4,578/4,924).
// ★★★★ **`Ceiling` 의 천장에는 벽이 둘이고, 둘은 같은 말을 한다** — ⓐ 사다리 첫 단(128비트)이 못 정했다,
//    ⓑ 이력이 **비용 한계**(`CACHE_REPLAY_COST_CAP` = 192)보다 깊어 싼 도로가 **아예 안 걸었다**. 읽는 쪽에도
//    되찾는 문에도 뜻이 같다: «더 물으면 답이 있다» ⇒ 한 변종이 맞다. `Unrealized` 는 그 반대 — 도로가 없다
//    (`NoMeet` 346 등). 그 갈림은 **사람이 읽을 것을 위해** 있다: «이 모델이 비싸다»와 «커널에 구멍이 있다».
// ⚠ **잔차는 캐시에서 죽었다**(`Measured` 변종과 `PointCache::residual`·`Model::vertex_tol`·`moved_to` 함께).
//    배열은 여전히 재고 셀프터치 체가 읽지만(`SeamVertex.tol`), 캐시는 저장하지 않는다 — 잔차는 담체까지의
//    거리 «하나»이고 좌표가 참값에서 얼마나인지를 말하지 않으므로, 캐시가 드는 «증명된 것»이 아니다.
//    ⇒ validate 는 **모든 정점**에 구성 ε 을 쓴다(`tol_of` 소멸). 잃은 것: 잔차를 들던 인구의 검사가
//    실측 1.07e-14 에서 1e-9 로 다섯 자릿수 느슨해졌다.
//   ★★ **«잰 것이 있는지»를 변종으로** — `tol: Option<f64>` 의 `None`/`Some` 이 나르던 것을 타입이 구조로
//     말한다. 미측정 변종엔 오차 자리가 아예 없다(있을 수 없는 것을 표현 안 함 — `SurfaceCache` 와 같은 원칙).
//   ★★★★ **이름은 «출처»(구성/발견)가 아니라 «캐시가 아는 것»이다** — 쓰는 쪽이 정했다. `transform` 은
//     기록된 이동을 받은 «발견» 정점을 `None` 으로 강등해 왔으므로(*"a recorded move used to demote to
//     `Moved`"*) `None` 은 «구성»이 아니라 «미측정»이었고, 2026-09-13 그림의 `Constructed` 로 이름했다면
//     그 자리가 이동된 발견 정점을 «구성»이라 철자했을 것이다 — `None` 은 침묵했지만 이름은 거짓을 말한다.
//   ★★★★ **`[Bounded; 3]` 가 아니다 — 잔차는 축별 경계가 아니다.** 2026-09-13 그림은 `Discovered([Bounded;3])`
//     라 적었는데, 그것은 아래 「세 «오차»를 섞지 말 것」이 경고한 바로 그 혼동이었다: `PointCache` 의 오차는
//     «저장된 점의 잔차(면에서 얼마나)» **하나**(`pierce_vertex_tol` = 담체 거리의 max)이고, `Bounded.error` 는
//     «참값 ∈ 값±오차» 를 **증명한** 축별 반경이다. 잔차를 거기 넣으면 타입이 증명 안 된 포함을 주장한다
//     (거의 퇴화한 담체 교차는 모든 담체에 가까우면서 정확한 코너에서 멀 수 있다). ⇒ **`PointCache` 는 계단
//     2 가 아니다** — 계단 2(`[Bounded;3]`)는 `realize` 출력과 `WitnessPoint` 에만 산다. «묶음»의 이득(transform
//     이 값만 옮기고 오차를 안 옮기는 desync 불가)은 변종이 그대로 준다.
//   ⏸ 고정밀판은 아직 «타입»이 아니다: 연산 하나짜리 실현이라 이름이 «Cache» 가 아니다(항목 23) —
//     `[HpBounded; 3]` 이거나 아래 `WitnessPoint` 자신이다(튜플 `(usize, [HpBounded;3])` 가 오늘 것).

// ── 계단 3: «경계 점 + 정의» — 판정 전용. 계단 2 의 상위집합이지 중복이 아니다. ──
//   ★★★★ `WitnessPoint = [Bounded;3](계단 2) + 정의(base:[Rat;3]·chain) + hp(메모)`. 그 «정의» 가
//     정확 단계의 입력이라 캐시(`PointCache`)와 «같은 것»이 아니다 — 이 차이가 진실/캐시 경계 그 자체다.
//     그래서 계단 1·2 는 합치되 이것은 캐시로 접지 않는다(접으면 경계가 지워진다).
pub struct WitnessPoint  { base: [Rat; 3], chain: /*모션*/_, realized: [Bounded; 3], hp: /*메모*/_ }

// ⏳ **한 곡면의 캐시 = 변종 «하나»** — 진실과 짝을 이룬다(`Surface::Plane` ↔ `SurfaceCache::Plane`).
//    2026-09-12 확정, 아직 안 지음(오늘 코드는 `{ realized: geom::Surface }` + 곁표 `surface_name`).
pub enum SurfaceCache {
    Plane    { realized: geom::Plane,    name: Option<PlaneName>, tol: ⏸ },
    Cylinder { realized: geom::Cylinder,                          tol: ⏸ },
}
// ★★★★ **왜 변종인가 — 이름이 평면에만 있기 때문이다.** 원통에는 4계수 음함수형이 없다(2차 곡면).
//    단일 구조체로 두면 `name: Option<PlaneName>` 이 **원통에겐 영원히 `None`** 인 필드가 되고,
//    그것은 타입이 「있을 수 없는 것」을 표현 가능하게 두는 것이다. 변종이면 원통에 이름 자리가
//    **아예 없다** — 그리고 구·원뿔이 와도 같은 방식으로 자란다(각자 자기 필드만).
// ★★★★ **왜 이름이 여기 사는가 — 이름은 캐시다.** `plane_name_exact(세 점)` 로 언제든 다시 나온다
//    (`Through` 는 정점을 풀어서 `plane_name_through`). 그런데 오늘은 **`pub surface_name` 곁표**에
//    따로 살아서, 같은 곡면의 캐시가 **두 벌**이고 그릇(`Vec`↔`HashMap`)·이름 규칙(`_cache`↔`_name`)·
//    가시성(비공개↔**`pub`**)이 셋 다 어긋나 있다. 캐시로 접으면 `PointCache`·`EdgeCache` 와 같은
//    「한 실체 = 한 캐시 구조체」가 되고, 재생 문 하나의 **비트 동일 잠금이 「이름은 캐시다」를
//    증명**한다(간선이 만든 선례).
//    재생을 막는 것은 «유도의 부재»가 아니라 **행(`raw`)** 이다 — 유도가 그것을 캐시에서 복사하므로
//    재생이 지우려는 값을 먼저 읽는다. ⇒ «**앵커만** 버리고 재생해도 비트 동일»은 오늘도 서고,
//    전체 재생은 행이 진실에서 나온 뒤에 선다(항목 36 의 같은 정정을 볼 것).
//    ☑ 실측(2026-09-12): `surface_name` 사용 44곳이 전부 «핸들로 조회»(`get` 24·`contains_key` 18·
//    `len`·`iter`·`insert` 각 1) ⇒ `HashMap` 이어야 할 이유가 없다.
// ★★ 접근자 이름은 §문의 이름이 정한다: 진실은 `surface(h)`, 캐시는 `surface_cache(h)`, 조각은
//    그 위에서 체이닝(`.plane()`·`.cylinder()`·`.name()`·`.tol()`). 다형 질의(`distance`·
//    `normal_at`)도 `SurfaceCache` 의 메서드다 — 오늘 geom 의 enum 이 하던 일 그대로.
// ☑ **대가 실측(2026-09-12)**: `m.surface(h)` 호출 **106**곳 중 **58 은 `match`/`matches!`**(변종을
//    바로 가르므로 오히려 나아진다) · **30 은 `let`-`else` 로 한 변종만**(형태만 바뀐다) ·
//    **손봐야 하는 건 18**(enum 통째로 11 + enum 메서드 7), 그중 대부분은 `distance`/`normal_at` 을
//    `SurfaceCache` 로 옮기면 사라진다. 통째로 변환하는 진짜 자리는 `transform.rs` **한 곳**.
pub struct EdgeCache     { curve: Curve }                          // 평가 가능한 담체 곡선
// ★★★ **캐시는 «실현값 (+ 필요하면 그 오차)»**:
//    `PointCache::Bounded{coord, bound}` ✔(칸 54: 캐시가 드는 오차는 **증명된 경계**뿐이다) · `EdgeCache{curve, err?}`(M7 후보) ·
//    `SurfaceCache{realized, tol?}` **동기 미확인**(위 §최종 타입 — 곡면 진실은 정확해 잔차가 없다).
//    ⇒ 오차 필드는 «소비자가 있으면» 붙지 대칭으로 붙지 않는다.
//    그 빠진 필드는 이미 이름이 있다 — design §3 의 `Intersection.cache_err`(*"캐시가 진짜 교차에서
//    벗어난 최대 거리"*). ⏸ 소비자 0 ⇒ 짓지 않고 적어 둔다(열린 항목 19).

// ★★★★ **세 «오차»를 섞지 말 것**
//
// | 무엇 | 성질 | 상태 |
// |---|---|---|
// | `Plane::distance_eps(p)` | **이 거리 계산**의 f64 반올림, `3ε·Σ|pᵢ−oᵢ|` — **`p` 에 의존**(두 연산자) | ✔ 있음, **저장 불가** |
// | `SurfaceCache.tol`       | **이 평면 자신**이 참 평면에서 얼마나 떨어졌나 — 평면의 성질(한 연산자) | ✗ 없음 **+ 동기 미확인**(아래) |
// | `SeamVertex.tol`         | 같은 것의 정점판(*"measured residual"*) — 배열의 작업 데이터, 셀프터치 체가 읽는다 | ✔ 있음 (칸 54: **캐시에서 이사**) |
// | `PointCache::Bounded.bound`     | 실현의 축별 경계 — 좌표가 참값에서 얼마나(반올림·사다리 반경), **잔차가 아니다** | ✔ 칸 52 |
//
// ⚠ 「곡면당 저장하는 tol 은 틀린 양」은 **첫째에만** 맞다. 그 근거는 코드가 적어 뒀다 —
//   *"the tolerance it has to hand is the vertex's … a point that is exactly on the plane can
//   still produce a nonzero residual"*, 그리고 그 사고가 실재했다(*"residual **exactly equal** to
//   the claimed tolerance, saved only by the comparison being strict"*).
// ☑ 반면 **둘째는 `PointCache::Measured.residual` 의 곡면판**이므로 저장하는 게 옳고, **없는 것은 「대체됐다」가
//   아니라 「아직 안 지었다」**다. 만드는 법은 이 절이 이미 적었다(아래: 세 정점의 실현에서 유도).
// ☑ `inv_norm` 도 «없는 게 아니라 다른 배치»다: 오늘 `Plane` 은 **단위 법선**(구성 시 한 번 정규화)
//   + `raw`(정확 계수용)를 들어 같은 정보를 갖는다. 계수 우선 캐시를 고를 때만 필요한 필드다.
```

★★★★★ **f64 가 «둘»이고, 둘은 만나지 않는다**. 아래 표의 세 줄은 **한 사다리의 세 단이 아니라 세 «수명»**이다:

| | 모델 캐시 (`PointCache`·`SurfaceCache`) | 실현 (`WitnessPoint.realized`·`WorkingCyl.realized`) |
|---|---|---|
| 누가 만드나 | 구성·불리언 경로가 f64 로 계산해 **넘겨준다** | 판정이 **정의에서** 정밀도를 불러 만든다 |
| 오차 | `tol` 하나 — 저장된 점의 **잔차**(면에서 얼마나) | 축별 경계 셋 — «참값 ∈ 값 ± 경계» **증명됨** |
| 누가 읽나 | 테셀레이션·STEP·물성·validate | **판정 1단(f64 필터)** → 2단 정수 → 3단 상승 |
| 수명 | 모델과 함께 | **연산 하나** |
| 판정 경로에 | **없다** — cip 이 읽는 곳 0건(실측) | 그 자체 |


| 캐시 | 키 | **수명** | |
|---|---|---|---|
| 모델 캐시 (f64) | 핸들 인덱스 (밀집) | **모델과 같이** | live 도달분만 lazy — `trial_bound` 는 일부러 우회. **판정 경로 밖** |
| 고정밀 메모 | (정의, 정밀도) (희소) | **연산 하나** | 판정이 실제로 만든 점에만 — 위 «실현» |
| 각도·√ 표 | `(Angle, prec)` / `(Rat, prec)` | 프로세스(스레드) | 값이 키의 순수 함수 — `cos 37°`·`1/√(n·n)` 은 어디서나 같다. 정확 반올림이라 실현이 유일하고 tol 도 유일 |

셋은 인덱스 공간이 달라 합치지 않는다. `Through` 평면의 `SurfaceCache` 는 세 정점의 실현에서
유도된다 — 법선 = 점 차의 외적이라 점의 오차 기계를 재사용하고, 이름 실현의 tol 은 지름길에만
쓰이므로 틀려도 조용한 오답이 아니라 느려질 뿐이다.

---

## 문의 이름 — `Model` 은 실체당 «둘», 조각은 체이닝 (2026-09-12 확정)

저장소가 비공개면 밖은 **문(=`Model` 의 메서드)** 으로만 읽는다. 규칙은 **한 줄**이다:

> **`Model` 은 실체당 문이 «둘» — `x(h)` 는 «진실»(아레나 항목 그 자체), `x_cache(h)` 는 «캐시».
> 그 아래 조각은 `Model` 에 문을 더 내지 않고 돌려받은 타입의 «메서드»로 꺼낸다.**

```rust
// ── Model 의 문 — 실체당 둘, 그게 전부다 ──────────────────────────────
m.surface(h)       -> &Surface        m.surface_cache(h) -> &SurfaceCache
m.vertex(v)        -> &Vertex         m.vertex_cache(v)  -> &PointCache
m.edge(e)          -> &Edge           m.edge_cache(e)    -> &EdgeCache
m.face(f) · m.shell(s) · m.solid(s) · m.motion(n)          // 캐시가 없는 실체는 하나뿐
m.surface_count() · m.vertex_count() · …                   // 개수는 예외적으로 Model 에

// ── 조각은 체이닝 ────────────────────────────────────────────────────
m.surface(h).motion()            -> Option<Handle<MotionNode>>  // 진실 조각
m.surface_cache(h).plane()       -> Option<&geom::Plane>        // 캐시 조각
m.surface_cache(h).cylinder()    -> Option<&geom::Cylinder>
m.surface_cache(h).name()        -> Option<&PlaneName>
m.surface_cache(h).tol()         -> ⏸ 「8·캐시 실형」
m.surface_cache(h).distance(p)   -> f64                         // 다형 질의는 캐시가 가른다
m.vertex_cache(v)                -> &PointCache                 // 변종(Bounded | Ceiling | Unrealized) — §캐시
m.vertex_cache(v).bound()        -> Option<&[Mag; 3]>           // Bounded 만 — 실현의 축별 경계
m.refine_vertex_cache(v, c, b)                                 // `Ceiling` 만 «올린다»(둘째 쓰기 문)
//   조각은 그 위에서 `match`. 세 변종이 모두 `coord` 를 들므로 `vertex_point()` 는 어디서나 답한다.
m.edge_cache(e).curve()          -> &Curve                     // (오차 err 은 M7 후보 — 위 §캐시)
```

★★★★ **왜 체이닝인가 — 이름이 «어디의» 조각인지를 스스로 말한다.** 초안은 `plane(h)`·`name(h)` 를
`Model` 에 달려 했는데, `Plane` 은 `Surface::Plane`, 즉 **실체 이름**이라 규칙 ①(실체 이름 = 진실)과
정면으로 부딪친다 — `plane(h)` 가 진실을 줄 것처럼 읽힌다. `m.surface_cache(h).plane()` 은 **이미
「캐시의」라고 말하고 들어왔으므로** 그 오해가 **구조적으로 불가능**하다.
☑ 그리고 새 방식이 아니다 — `geom::Surface` 가 `distance()`·`normal_at()` 을 자기 메서드로 갖고
안에서 가르는 것과 **같은 관용구**를 캐시에 적용하는 것이다.
☑ **`Model` 이 안 자란다**: 구·원뿔·NURBS 가 와도 `SurfaceCache` 에 `sphere()` 가 붙을 뿐이고,
「8·캐시 실형」이 `tol` 을, 수렴 주석이 `EdgeCache::err` 을 더해도 문 개수는 그대로다.

⏳ 오늘의 예외 둘은 **아직 규칙 밖이다**: `plane_motion(h)` → `m.surface(h).motion()`,
`pub surface_name` 곁표 → `m.surface_cache(h).name()`. ★ 칸 56 이 저장소 다섯을 닫으면서
`surface_name` 이 **`Model` 의 마지막 공개 필드**가 됐고, 칸 57 이 그 접기를 **열린 항목 36** 으로
냈다(캐시가 못 보는 조용한 실패 모드가 있어 같은 칸에 넣지 않았다) — 이 절이 정한 **철자**를 그
항목이 **집행**한다. `plane_motion(h)` 은 아직 주인이 없다.
☑ `m.motion(n)` 은 **이미** 이 규칙이다 — 새 규칙을 만드는 게 아니라 넓히는 것이다.


☑ **나머지는 기계적이다** (실측 2026-09-12): `vertex_point` **111** · `edge_curve` **29** ·
`plane_motion` **18** · `surface_name` **44**.

☑ **성능은 이 결정의 고려사항이 아니다.** 어느 철자든 안에서 하는 일은 «디버그 가드 + 배열 색인
하나 + 판별자 읽기»이고 할당도 복사도 없다(`#[inline]`). ★ 오히려 **줄어든다**: 오늘 이름과 실현이
따로 살아 둘 다 필요하면 조회가 둘(`Vec` 색인 + **`HashMap` 해시**)인데, 체이닝이면 `Vec` 색인
**하나**에서 둘 다 꺼낸다.
⚠ 대가는 `PointCache`·`EdgeCache`·`SurfaceCache` 의 필드가 비공개이므로 **조각마다 메서드를 하나씩
지어야 한다**는 것 — 작지만 0은 아니다.

## 판정 (연산 동안만 산다 — nacre-cip)

이름 규칙 셋: ① **`Working*` = 모델 타입의 판정층 쌍둥이이자 표의 뿌리**(수명 = 연산 하나 —
`Model` 에 `Working…` 이 담기면 읽는 즉시 이상해 보여야 한다). ★ 접두사는 **수명 표식이
아니라 뿌리 표식**이다 — 부품(`WitnessPoint`·`MoveNode`·캐시들)은 `Working*` 컨테이너 안에
살며 수명을 상속하므로 접두사를 반복하지 않고, 이름은 역할을 말한다(캐시 타입은 모델 쪽과
일부러 공유된다 — 반올림 사본임이 타입에 보이도록). ② **거울의 변종 이름은 진실과 같다**
(`Known`/`Through`). ③ 평면 정의를 이루는 유리수 점은 **증인 점 `WitnessPoint`** 다.
오늘 코드와의 대응 — **판정층 개명 집행 완료(2026-08-08)**, 항목별 처분:
- ★★★ **`FaceInfo` 는 따라가지 않았다 — 여기 적혀 있던 «둘 다 `WorkingPlane`» 은 정정한다.**
  면/평면 분리가 오늘 하중을 진다: `orient_sign`(이 면의) vs `frame_sign`(그 클래스의)은
  한 함수가 두 종류 인덱스로 불리던 시절 — *규약을 단언할 수 없던 유일한 자리* — 를 가른
  수선이었고, 한 이름으로 합치면 고친 결함의 이름이 되살아난다. 합치기는 개명이 아니라 설계
  작업이며 16 이후 재검토.
- `WorkingPoint` → `WorkingVertex`: 코드에 아직 없음 — 16이 만들 타입이 이 이름으로 태어난다.
- `Standard` → `ProofStandard`: 같은 계열(모양이 다르다 — `same_within` 유도로의 재구성) — 대응
  명시가 없어 유예. **기준: 이 목록에 명시된 것만 기계적 개명이다.**
(점 실현 묶음 `HpPointCache` 는 아직 타입으로 없음 — 오늘은 `(usize, [HpBounded; 3])` 튜플.)

```rust
/// 판정용 평면 — `Surface::Plane` 의 쌍둥이. 정의를 펼쳐 들고 두 실현을 메모한다.
pub struct WorkingPlane {
    pub def: [WitnessPoint; 3],               // 정의 — 증인 삼각형 하나로 전체(16-3 정정)
    pub name: Option<PlaneName>,              // 이름 — Narrow 는 Shewchuk, Wide 는 BigInt(항목 15)
    pub chain: Rc<[MoveNode]>,                // 사슬을 펼친 것. MoveNode = Motion 의 판정층
                                              // 펼침 — 핸들 숲은 판정층이 못 푸므로 경계에서
                                              // 핸들을 해소한다(Frame{plane} → PlaneFrame).
                                              // ★ 평면당 한 번 펼치고, Known 의 증인 점 셋은
                                              // 이 **같은 Rc 를 복제**해 든다(할당 하나,
                                              // 포인터 넷 — 두 번 펼치지 않으므로 어긋날 수
                                              // 없다). 점이 드는 이유: 자기완결 판정
                                              // (shared_base·realize 는 점만 받는다). 평면이
                                              // 드는 이유: Through 는 증인 점이 없고, 평면
                                              // 수준 질문(합성 회전 증명·구조 검사·계수 tol)
                                              // 이 있다.
    pub cache: SurfaceCache,                  // 실현 1 — f64 필터, 항상
    hp: Rc<OnceCell<(usize, HpSurfaceCache)>>,// 실현 2 — 고정밀, 필요할 때만
}

/// ★★★★★ **정정(2026-08-09, 16-3): `Through` 변종은 지어지지 않는다 — 판정층은 증인 삼각형
/// 하나로 전체(total)다.** 판정 표의 계약은 «평면 위 세 정확한 점» 이고, 판정된(이름 없는)
/// 평면은 자기 프레임의 probe 로 그 점을 정의상 갖는다(반증표 참조). 거울 이름 규칙 ② 의
/// 한정: 진실의 `Through` 는 거울에서 **probe 로 유도된 증인 삼각형**으로 나타난다 — 변종이
/// 아니라 유도다.

/// 증인 점 — 유리수 base 를 모션 사슬로 나른다. base + chain 이 정의, realized 는 f64 실현
/// (축마다 값+경계 한 원자 — 모델 캐시가 아니라 연산 하나의 것), hp 는 고정밀 메모.
/// 같은 정의는 같은 실현(경로 무관)이다.
pub struct WitnessPoint {
    pub base: [Rat; 3],
    pub chain: Rc<[MoveNode]>,
    pub realized: [Bounded; 3],
    hp: Rc<OnceCell<(usize, [HpBounded; 3])>>,  // `HpPointCache` 라는 타입은 없다(열린 항목 23)
}

/// 판정용 정점 — `Vertex { surfaces: [3] }` 의 쌍둥이. 세 평면을 가리키기만 하고
/// 좌표를 만들지 않는다. 캐시하는 것은 나누지 않은 동차좌표 [X:Y:Z:W] — 나누면 무리수가
/// 되고 오차가 낀다. 부호 질문은 sign(X − …·W)·sign(W) 처럼 곱셈·뺄셈만으로 답한다
/// (indirect predicates).
///
/// ★ 판정층의 점이 둘(`WitnessPoint`·`WorkingVertex`)인 것은 질문이 둘이기 때문이다 —
/// «값을 실현해 달라»(정의가 값) vs «좌표 없이 부호만 답해 달라»(정의가 교차).
/// 메모 규칙도 여기서 갈린다: **메모는 정의에만 붙는다** — `WitnessPoint`·`WorkingPlane` 은
/// 원천이라 HP 메모를 갖고, `WorkingVertex` 는 파생이라 f64 층만 두고 고정밀은 평면
/// 메모에서 재계산한다(두 번째 원천을 만들면 어긋날 수 있다).
pub struct WorkingVertex<'a> {
    pub planes: [&'a WorkingPlane; 3],
    homog: OnceCell<[Bounded; 4]>,            // f64 층만 — 고정밀은 평면 메모에서 재계산
}

/// 이 연산이 무엇을 증명으로 인정하나 — 전부 모델에서 유도, 설정 없음.
pub struct ProofStandard {
    pub model_size: Mag,                      // 순회 순서 무관 (재생 결정성)
    pub start_bits: usize,
    pub max_bits: usize,                      // 비용 한계이지 해상도 한계가 아니다
}
impl ProofStandard {
    /// model_size · 2⁻¹⁸⁰ — 이보다 가깝다고 증명되면 일치다 (규칙 5).
    pub fn same_within(&self) -> Mag { … }
}

pub enum Decision {
    Sign(Orient),                             // 증명된 부호
    Coincident { within: Mag },               // 증명된 가까움 — 근거를 실어 보고
    Exhausted { at: usize, within: Option<Mag> },
    Degenerate,                               // 마지막 둘 → 이름 붙은 거절
}
```

- **정확 경로의 조건 = 세 평면이 같은 모션을 공유할 것**(`None` 포함) — 그 프레임 안 정수
  계수로 Shewchuk. 모션이 섞이면 f64 필터 + 고정밀 상승(C4 — 기존 규칙이지 새 예외가 아니다).
- 모션이 계수에 하는 일: `Translate` 는 `d' = d − n·t`(유리수면 정확), 축과 나란한 법선은
  `R·n = n`(구조적 검사 가능), 그 밖은 실현에 tol. 진실 쪽에는 오차가 없다 — 오차는 실현할 때
  생기고, 그것이 `tol` 이다.
- `Through` 평면의 판정은 두 갈래다.
  ★★★★ **«유리수 닫힘이면 `Wide` 이름이 곧 정확 계수라 공짜» — 이제 조건부로 참이다**
  (2026-08-09, 열린 항목 15 닫힘). 세 술어(`orient3d`·`cmp_coord`·`dir_sign`)에 이름-정수
  rescue 가 있다: 게이트 = 전 평면 이름 있음 ∧ 하나 이상 wide ∧ (전부 무이동 **또는** 전부
  한 사슬 — `cmp` 는 무이동만). 이름은 구성 시 저장 방향으로 σ-접혀(`name_stored_ints`)
  BigInt 부호 술어(`int_*`, scalar)로 정확히 답한다. **조건의 이유**: wide datum 의 이름은
  세계 이름(담체 = 발견 정점)이라, 그 평면을 처음 만드는 혼합-프레임 불리언에서는 일부만
  잡히고(454→330), 전부-세계/통째-이동 표(2세대)에서 완전히 공짜다(139 대 20 — 대조군인
  narrow 보다 싸다). 그릇은 `Expansion` 이 아니라 BigInt(반증표).
  **무리수 모션이 낀 datum 만** 실현 상승 전용이다: 정점을 자기 세 평면의 고정밀 실현에서
  동차좌표로 만들고 그 위에서 계수를 구간으로 유도한다(**차수 9** — 위 참조) — 정확성 위험이
  아니라 비용 위험이었고, ★ **이제 기계가 있어 쟀다**(S5(ii)-2a, `frame3.rs` 단위 테스트):
  - **깊이 1 은 필터가 산다.** 생성 200 사례 전부 결정, 계수 800개 중 미결 **0**, 최악 상대
    반경 **6.3e-10**. 차수 9 가 `Bounded` 의 여유를 먹지 않는다.
  - ★★★ **깊이 2 는 필터가 없다.** 담체가 또 `Through` 면 차수가 **81** 이 되고, 계수가
    `f64` 범위를 **8/8 전부** 벗어난다(고정밀 쪽은 멀쩡하다). 함수는 그때 `None` 을 돌려
    **상승으로 보낸다** — `NaN` 반경이 우연히 `sign()=None` 을 내는 것에 기대지 않는다.
    ⇒ 깊이 2 는 «느린 길» 이 아니라 **상승 전용**이다. 깊이 제한을 두어야 하는지는 **열린
    항목 16 의 질문**이다(C6 이 유계를, C5 가 비순환을 이미 준다;
    남은 것은 비용뿐이다).
  - ★★ **배율의 부호는 값 안에서 없앤다.** join 은 행에 대해 다중선형이라 결과가 참 평면의
    `D0·D1·D2` 배이고, 음수면 **평면 방향이 뒤집힌다**(`frame_sign`·바깥 법선·라벨 프레임이
    전부 그 위에 있다). `Judge::plane_iv(k) -> [Bounded;4]` 에 부호를 실을 자리가 없으므로 —
    실을 곳 없는 값은 아무도 안 쓰는 값이다 — 함수가 스스로 정규화한다.
  - `D` 가 0 을 품으면 `None` → 상승 → 안 갈라지면 이름 붙은 거절. 조용한 폴백은 없다.
  - ★★★★★ **동차 상승의 종착**(갱신 2026-08-09, 16-3): `plane_hp_through`(고정밀 판)는 판정
    프레임의 동차 경로가 소비하고, `plane_iv_through`(f64 필터 판)는 **소비자 없이 은퇴했다**
    — 판정 표가 계수 경로를 아예 필요로 하지 않았기 때문이다(판정 평면의 증인 = 프레임 probe,
    §열린 항목 16). 잠금은 Pure-대-Meet 차등으로 이사.

---

## 제약 — 어떤 답이든 만족해야 했고, 확정안이 만족하는 것

| | | 확정안에서 |
|---|---|---|
| C1 | 진실은 하나, 나머지는 캐시 | 생산자는 점(또는 정점 핸들)만 진술, 이름·좌표·계수는 전부 유도 |
| C2 | base case 는 «수» | `PlanePoints::Known` — 모든 사슬이 유리수 세 점에서 끝난다 |
| C3 | 발견된 좌표는 `Rat` 에 안 들어간다 → 가리킨다 | `PlanePoints::Through` |
| C4 | 모션 섞인 판정은 f64 필터 + 상승 | `WorkingVertex` 동차좌표 + `Decision` — 기존 규칙 그대로 |
| C5 | append-only ⇒ 참조는 과거로만 (DAG) | Through(평면→정점)는 datum 이 정점보다 나중이라 성립; 벽은 값(`Known`)이라 순환 자체가 없다 |
| C6 | 깊이는 사용자 조작당 1단 | 불리언은 평면을 만들지 않는다 — Frame·Through 재귀가 그래서 유계다 |
| C7 | 거절은 정직하게, 이름 붙여 | 규칙 5 의 3갈래 |
| C8 | `Rat` 을 넓히지 않는다 | 넓어진 것은 이름의 **그릇**(동일성 전용)과 중간값뿐 |

---

## ★★★★★ 반증된 답들 — 같은 길로 다시 가지 말 것

2026-08-04~05 에 실제로 제안됐다가 측정·문서·사용자에게 반박당한 것들. 근거 서사는 git 이력.

| 제안 | 실제 |
|---|---|
| *"`Vertex` 는 세 평면 단일형(`{surfaces:[3]}`)"* (Q3 원안) | ✗ **S7 구현이 반박**: 원통 seam 정점은 «원통 ∩ 캡» = 테두리 **원 전체** 위의 한 점이라 세 평면이 없다. 슬롯을 중복으로 채우는 «잠정 표기»는 타입이 참이 아닌 말을 하게 만들고(세-평면 교점 아님) 판정의 `D=0` 을 조용한 0으로 흘린다 — 사용자 판단으로 **두 변종 enum**(`ThreePlane \| OnSeam`) 채택. 각 변종이 자기 진실을 말하고 불변식은 변종별(Q5), M6 확장(`Branch`·`Apex`)은 변종 추가로 받는다. 단일형의 원래 근거("정점은 술어가 가장 많이 소비")도 이미 소멸 — S6b 이후 판정층은 정점을 읽지 않는다 |
| *"`Origin::Constructed` 에 가드를 달자"* | ✗ 그 태그는 정확성을 뜻하지 않는다 — 정확 경로로 만든 정점도 `Constructed` 이고 그 `point` 는 반올림이다 |
| *"불리언 능력을 포기하자"* | ✗ 설계 위반 — 이 문서는 `Inexact` 를 지우기로 했다 |
| *"평면의 진실을 계수로"* | ✗ 계수는 두 점 차의 **곱**이라 `i128` 적합 25.8%(점은 100%) — 담을 그릇이 없다 |
| *"f64 계수를 `Rat` 으로 들어올리자"* | ✗ 반올림이 진실에 구워지고, 같은 벽이 두 경로로 두 평면 클래스가 된다 |
| *"넓은 이름이 프레임을 연다"* | ✗ 막던 것은 `n·n` **캐시**의 `i128` — 캐시가 안 만들어진다고 진실을 거절하던 것이 병이다 |
| *"`n·n` 을 저장해야 한다"* | ✗ 파생물이다 — 필요하면 정의에서 다시 실현한다 |
| *"모션은 정점이 든다"* | ✗ 면이 든다 — 정점은 세 면의 교점이라 저절로 따라온다. 정점에 달면 3중 중복 |
| *"평면이 정점을 가리키면 무한 재귀"* | ✗ 불리언이 평면을 안 만들어 안 자란다(C6) — 깊이는 사용자 조작당 1단 |
| *"datum 평면은 M6+"* | ✗ M5 다 — 곡면이 없다 (사용자 확정 2026-08-05) |
| *"코퍼스 최대값이 상한"* | ✗ 폭은 타입에서 유도하라(이름 ~2²²⁹¹) — 실측 최대는 표본이다 |
| *"정준 원점이 넘치면 `points[0]` 으로 폴백"* | ✗ 스케치의 (0,0) 위치가 **오버플로 여부에 따라 달라진다** — 유도(`Canonical`)로 가면 실현이 임의정밀로 감당하고 규약이 유지된다 |
| *"좁은 이름을 별도 곁표로"* | ✗ ~99% 가 중복 저장이다 — `PlaneName = Narrow \| Wide` 한 enum + `narrow()` 투영이면 저장은 하나다 |
| *"판정층에 `WorkingPlaneDef::Through` 가 필요하다"* (2a~16-2 내내 내가 가정) | ✗ **판정 표의 계약은 «평면 위 세 정확한 점» 이었고, 판정 프레임의 probe 가 그것을 이미 준다** (원점 = 수선의 발, û·v̂ = 평면 안 — 전부 정의상 평면 위). 16-3 이 그 가정을 실행으로 반증했다 — 판정층 새 기계 0으로 불리언이 열렸다. ★ 교훈: 소비자를 위해 지어 둔 기계(`plane_iv_through`)는 그 소비자와 함께 은퇴했다 |
| *"S5(ii)-2 는 판정층 배선이다"* (2026-08-08 내가 이행표에 적음) | ✗ **벽은 프레임이었고 판정층은 애초에 막고 있지 않았다.** 이름 없는 평면은 `SketchFrame` 을 못 얻어 **판정 표에 도달조차 못 한다** — 그래서 `WorkingPlaneDef::Through` 는 생산자가 «없는» 게 아니라 «있을 수 없다». §열린 항목 16 이 두 벽을 순서대로 적는다. ★ 교훈: «다음 단계가 이 기계를 쓴다» 를 적기 전에 **그 단계가 도달 가능한지** 먼저 확인할 것 |
| *"wide 계수는 `Expansion` 에 넘기면 된다"* (항목 15 원문) | ✗ **f64-expansion 은 wide 를 못 담는다** — 조각도 f64 라 지수 상한이 ~2¹⁰²³ 인데 wide 이름은 ~2²²⁹¹ 실측. `2^53 은 f64 하나의 한계일 뿐` 은 **i128 까지만** 참인 문장이었다. 옳은 그릇은 BigInt 정수 산술 그 자체(항목 15 닫힘) |
| *"wide 상승 887/454 의 인구는 «움직인» wide 평면(담체의 공유 프레임)"* (항목 15 계획, 검토 2회 통과) | ✗ **wide datum 의 담체는 발견 정점이라 사슬이 없다 — 이름은 세계를 말한다**(`the_wide_plane_names_the_world` 가 잠금). 그 평면을 처음 만드는 불리언은 혼합 프레임이라 공유-모션 게이트가 영원히 닫혀 있었을 것이고, 픽스처의 패리티 단언 하나가 구현 중에 이것을 잡았다. ★ 교훈: «어느 프레임이 말하는가» 는 평면마다 담체에게 물을 것 — 조사·검토가 아니라 **표를 찍은 probe** 가 반증했다 |
| *"rat 거절 인구는 reuse 에 도달하지 않아 프로덕션 이득 아마 0"* (항목 0 의 2026-08-08 정정) | ✗ **정정 자체가 인접 명제였다** — 발견 정점(2/3 거절을 잰 인구)이 도달하지 않는 것은 맞지만, 도달하는 인구는 따로 있다: 움직인 솔리드의 **구성** 코너. probe 실측 wf 프레임 프리즘 8/8 넘침·8/8 드는 점 ⇒ 십진 프레임 피연산자는 클래스 reuse 를 통째로 잃고 있었다. ★ 교훈: «인구 X 는 도달하지 않는다» 를 적을 때 **도달하는 인구가 무엇인지**까지 물을 것 — 부정 명제는 반쪽 측정으로도 참처럼 보인다 |

---

## 이행

### 관문 규칙

- 기본 관문은 **위상 정확 일치 + 좌표 ε(모델 크기 상대, `2⁻⁴⁰` — 실측 여유 2000배)**, 비트
  동일은 보너스 신호. 합격/불합격이 아니라 **최대 편차 숫자를 찍는다**(누적 드리프트 감시).
- ★★ **8b — 비트 동일 관문은 그 안의 population 에 대해서만 보증한다.** 바꾸려는 것이 닿는
  population 을 대장에 먼저 넣고(17자리 `fw`·기울어진 `tp` 가족은 이미 있다), *"어느
  population 인가"* 는 추측하지 말고 계측이 이름을 대게 한다. 좌표 관문은 «답은 같은데 더 나쁜
  길로 갔다»를 원리적으로 못 보므로 **"정확 경로를 탔는가"를 직접 단언하는 테스트**를 함께 둔다.
- ★★★ **두 `Surface` 의 철자**(2026-09-11, 영구 규칙 — 열린 항목 8): 한 파일이 둘 다 필요하면
  **맨이름은 진실**(topo, `Handle<Surface>` 가 이름 짓는 것)이고 실현은 `nacre_geom::Surface` 로
  적는다. 반대로 하면 같은 파일에서 `Handle<Surface>` 는 진실을, 맨 `Surface::Plane(p)` 는 캐시를
  뜻해 한 단어가 표지 없이 두 가지가 된다. ☑ 섞으면 **하드 오류**다(geom 은 tuple 변종, 진실은
  struct 변종 ⇒ E0532) — 조용히 틀릴 자리가 0.
- ★★★ **타입 구조를 바꾼 이행 행은 §최종 타입·§캐시의 «그림»도 같이 고친다** (2026-09-12). 이 문서는
  기록 층(이행표·열린 항목)이 자라고 그림 층은 얼어 있는 구조라, 행만 쓰면 그림이 조용히 뒤처진다.
- **어느 크레이트에 두나** — 「유도된 값 + 그 산술」은 `nacre-scalar`(최하단; `nacre-cip` 이 닿아야
  하므로), 「아레나 항목의 진실」은 `nacre-topo`. 구조적 강제: **`nacre-cip` 은 `nacre-topo` 에 의존하지
  않는다** ⇒ 판정이 쓰는 타입은 전부 topo 아래에 있어야 한다(`PlaneName` 이 scalar 에 사는 이유).
  본문은 design.md 「크레이트 구조」.
- **병렬 불변식**: 병렬 구간에서 **곡면 push 금지**(재생 결정성). 오늘의 문은 `Model::push_plane`·`push_cylinder`(사설 `push_plane_raw`·`push_cylinder_raw` 로
  모인다). ☑ **M6 뒤 재측정(2026-09-11)**: 워크스페이스의 진짜 병렬은 `par.rs` **한 파일 두 자리**
  (`(0..n).into_par_iter().map(f)`)뿐이고, 그 클로저는 「인덱스 → 값」이라 **`&mut Model` 을 들 수 없다**
  ⇒ 병렬 구간의 곡면 push 는 **구조적으로 불가능**하다. 그 파일에 `push_plane`/`push_cylinder` **0건**, 그리고
  *"병렬이 두 곳에 살면 조용히 깨진다"* 는 텍스트 가드가 이미 그 한 파일을 지킨다. 규칙은 유지, 감시는
  그 가드가 한다.

---

## 표현력의 경계 — 어휘 밖의 사용례와 확장 경로 (2026-08-05 점검)

대표 사용례를 타입에 대입해 훑은 결과. **조용히 틀리는 시나리오는 없다** — 아래는 전부
«현 어휘로 진술 불가 ⇒ 구성 시점 거절»이고, 각각 확장 경로가 있다.

| 사용례 가족 | 왜 막히나 | 확장 경로 |
|---|---|---|
| ★★ **각도-계열 평면**: 각진 datum(모서리 축 θ)·드래프트 돌출·일반 챔퍼 | 진실이 *"이름 붙은 직선 둘레로 유리수 각 θ 돌린 평면"* 인데 그 어휘가 없다 — 기운 벽의 셋째 점은 `dist·tanθ` 라 무리수(`Known` 불가), `Rotate` 는 축 정렬 전용 | `Frame` 과 같은 «이름으로 들기»: 모서리(정점 핸들 둘)를 축으로 지목하는 회전 변종. 실현은 기존 기계(각도 캐시 + `1/√`) 그대로. ★ **챔퍼가 M6 이므로 M6 가 이 어휘를 요구한다** — 정점 가지 번호 다음의 둘째 항목 |
| 임의 축 회전·임의 평면 미러(솔리드 이동) | `Rotate`/`Mirror` 축 정렬 전용 — 기존 비목표(design.md), 이 설계가 바꾸지 않았다. 새 피처의 임의 방향 배치는 `Frame` 으로 이미 된다 | 위와 같은 계열 — 그 설계가 들어올 때 함께 |
| 혼합 정의 datum (정점 2 + 좌표 1) | `Through` 는 핸들 3 전용, `Known` 은 값 3 전용 | `Through` 원소를 `PlanePointRef(At \| Vertex)` 로 일반화 — 평면의 점 데이터 수준이라, 기각된 `Vertex::At`(정점 아레나)과 다른 자리다 |
| 중간 평면(midplane) | 같은 프레임 평행 면 사이는 유리수 평균이라 `Known` 으로 계산 가능 ✓. **프레임이 다른** 두 면 사이는 세계 좌표가 무리수 + 어느 정점도 안 지난다 | 필요해지면 정의 변종 `Mid([Handle<Surface>; 2])` |
| 스케치 내부의 정확한 각도(정육각형·30° 변) | 꼭짓점에 √3·tan30 — `[Rat; 2]` 로 못 적는다. 앱이 계산한 f64 의 십진수가 진실이 된다(§design 6.0 — 의도된 동작) | **정확한 각도는 좌표가 아니라 모션으로 표현된다**(회전·프레임) — 결함이 아니라 이 숫자 시스템의 정의적 성질. 여기 명시해 둔다 |

검사했고 막히지 않는 것: datum 이 참조한 솔리드의 이동(이동본 정점으로 재지시), superseded
참조(append-only), 원형 패턴(`360/n` 은 항상 유리수), 깊은 체이닝(비용만 — 실측 회전 3200회),
datum 낀 불리언(비용만 — 아래 1), 4+평면 정점(validate 몫), 오프셋의 오프셋(연산당 1단),
STEP 출력, undo/replay.

---

## 열린 항목 — 측정할 것과 알려진 절벽

1. **`Through` 판정 비용** — 두 갈래이고 **오늘 잴 수 있는 것은 하나**다(실측 2026-08-08,
   `nacre-ops/tests/wide_datum_cost.rs`).
   - **(a) 유리수 닫힘** — 문항이 틀렸다: 조각 수가 아니라 **없는 정확 경로의 대가**다(§판정의
     정정 참조). 실측: 이름이 `Narrow` 인 datum 으로 자르면 상승 **278**, `Wide` 면 **887**
     (mean 256비트, 고갈 0). ★ 두 팔은 서로 다른 정점 삼중항을 써야 하므로(그래야 한쪽이
     wide 다) 도구 위치가 다르다 — **3.2배는 계수가 아니라 자릿수**로 읽는다.
   - **(b) 무리수 모션** — **잰다고 적었고, 기계를 지으면서 쟀다**(S5(ii)-2a, 2026-08-08).
     차수는 ~12 가 아니라 **9**(사영 join). **깊이 1: 필터가 산다**(200/200 결정, 미결
     0/800, 최악 상대 반경 6.3e-10). **깊이 2: 필터가 없다**(차수 81, 계수가 `f64` 범위
     밖 8/8) ⇒ 상승 전용. ⇒ 문서의 «비용을 먼저 재고 S5(ii)-2 로 간다» 는 순서가 (b) 에는
     성립하지 않았다 — 잴 기계가 그 단계에서 만들어지므로, 짓되 재면서 지었다.
   필터가 무력하면 `SurfaceCache` 선실현이 다음 수다 — **깊이 2 가 정확히 그 자리**다.
2. **무리수 모션이 낀 datum 은 이름이 없다** — interning 불가, 동일성은 술어가 매번(규칙 6 의
   한정). 느릴 뿐 틀리지 않는다. ⚠ **미측정 주장**(2026-09-12 표지): 이웃 항목(1·2b)은 실측치를
   다는데 이 줄만 근거가 없다. 불리언은 **평면 클래스별로** 셀 복합체를 만들므로, 같은 평면이
   두 핸들이면 «두 클래스»가 되고 — 그것이 정말 비용뿐인지(결과 동일)는 **같은 평면을 두 진술로 넣은
   불리언을 하나 만들어 census 로** 재야 안다. `Verbatim` 열쇠(§이름과 interning)가 이 인구를 낳는다.
2b. **wide 프레임 실현 비용** — `FrameWide` 의 축 실현은 캐시 없이 점마다 돈다(인구가 작아
   수용, S4). 그리고 S4 구현 중 실측 하나: **스케치→돌출 벽의 정준 이름은 ~115비트에 캡**
   (십진 창이 곱을 묶는다) — `Wide` **이름**의 면은 오늘의 구성 경로에서 안 나오고, 첫 생산자는
   datum(S5)이다. `n·n` 넘침(1.6%)은 그 경로에서 실재하며 S4 가 열었다(census `wf` 가족).
4. **정의만 가리키는 평면의 순회·직렬화** — `Model::reachable` 이 면을 통해서만 돈다. 세계 축
   평면·datum 이 가리키는 평면·이동본을 정의 경유로도 따라가야 한다. ★ S9 부로 실재 인구가
   생겼다: 씨앗 셋은 심긴 직후엔 어느 면도 참조하지 않는 orphan 이다(`Reachable` doc 정정
   완료 — validate 는 도달 집합만 검사해 무위반). ★★ **S5(i)-a 는 이것을 강제하지 않았다**
   (확인, 2026-08-07): 스냅샷 포맷이 **존재하지 않고**, `nacre-step`·`nacre-tess` 는 surface
   store 를 **아예 돌지 않으며**(둘 다 면 경유), validate 는 도달 집합만 센다 ⇒ 쓰이지 않는
   datum 은 어디로도 새지 않는다. 일반화는 스냅샷 포맷이 생길 때다.
5. **M6 절벽** — 이차곡면 셋은 최대 8점에서 만나 `[Handle; 3]` 이 «어느 점»을 못 말한다. 가지
   번호(`branch: u8`, 결정적이어야 함) 또는 재명명 — M6 에서 실제 형상을 만나 결정. 통일안을
   채택했으므로 이 문제는 모든 점에 걸린다 — M6 의 가장 큰 항목.
   ★ 모서리에도 같은 자리가 있다: 원통 **seam** 은 두 면의 교차가 아니라 **한 면의 매개화
   이음매**라(양쪽이 같은 원통) 두-면 교차 담체로 적히지 않는다 — S8 은 잠정 표기
   `[h_cyl, h_cyl]` 자기-인접으로 적었고 validate 가 «담체 동일 ⇔ 원통» 을 지킨다
   (`EdgeCarrierMismatch`); 가지 번호·매개 표현은 여전히 M6. S8 은 평면 모서리를
   대상으로 하고, seam 의 담체 표현은 M6 에서 원통의 진실(`ref_dir`)과 함께 결정한다.
6. **폭 주장은 타입 경계에서만 보장이다** — 이 문서의 백분율(0.39%·1.6%·25.8%·71.9%…)은 전부
   **코퍼스 수치**다. 상한이 필요한 자리에는 타입에서 유도한 값을 쓴다.
7. **`check()` 의 진실-위 비용** — 실측(S3, 볼록 링·17자리 좌표): 34ms@100점·3.2s@1000점·
   78s@5000점, 호출당 ~1.7µs(narrow 경로의 gcd 약분이 지배; 이런 좌표는 전-narrow).
   손 스케치(수십 점)는 ms 미만이라 수용, 수천 점 생성기가 실재해지면 그때의 수: (a) Rat
   비교 기반 정확 bbox 선별(교차쌍 대부분 기각), (b) 실현 f64 + 건전 오차 한계 필터 → Rat
   상승(CIP 필터 철학의 2D 판). **둘 다 그 인구가 생기기 전엔 짓지 않는다.**
11. **노드 생략의 더 강한 게이트** — 오늘의 게이트는 `exact()`(축이 정확 유리수 직교로
   리프트되는가)다. `PlaneFrame`+`inv_sqrt_exact` 로 «실현이 정확 f64 에 떨어지는가»를 직접
   묻는 더 강한 게이트가 가능하지만, **표현식을 바꾸면 노드 인구가 움직인다** — 두 ★★ 주석의
   경고 그대로, 교체는 census 관문 동반 필수(S9 에서 기록만).
12. ⏳★ **정점 캐시엔 «버리고 재생» 보증이 없다(S7)** — `rebuild_edge_cache` 의 정점판을 짓지
   않았다. 정점 `D≠0` 의 완전한 유리수 단언(문서 :85 의 «공짜 단언»)도 같은 이유로 유예 —
   오늘은 `push_vertex` 의 핸들 상이성 debug_assert 까지.

   ⚠★★★ **정정(2026-09-11) — 여기 적혀 있던 근거 셋이 전부 지나갔다.** 원문은 *"3b(좌표 재생)가
   ⏸ 이고, 발견 좌표는 배열이 공들여 만든 값(1992 중 238 이 순진 Cramer 와 다름)이며 seam 좌표는
   M6 까지 load-bearing 이다 … 좌표 재생이 생기는 자리(M6/판정 통합)에서 셋을 함께 연다"* 였다.
   - **3b 는 절반 섰다** — §이행 3b 행: 묻는 문(`nacre_ops::realize_vertex{,_decimal}`)은 섰고
     덮어쓰기가 남았다. ⇒ 「⏸ 이므로」는 더 이상 근거가 아니다.
   - **238/1,992 는 «순진 Cramer» 와의 차이였다** — 항목 15 가 그 근거를 반증했고, 칸 ㊵ 가
     정확 반올림 실현으로 실측했다: 축정렬·불리언 인구는 **48/48 · 96/96 일치**, 기울어진 프레임은
     **12/12 불일치(최대 4 ulp)** ⇒ 그쪽에서는 **캐시가 틀린 쪽**이다.
   - **M6 는 도착했다** — `OnSeam` 의 좌표는 `CylinderDef` 로 정확히 지정된다(§최종 타입).
   ⇒ 남은 것은 「보증이 없다」는 **사실**이지 그 근거들이 아니다. ✔ **정점판은 곡면보다 먼저 열렸다(칸 52)** —
   «곡면 실현과 함께»라 미뤄 둔 근거(정점만 하면 자기모순)가 실측 1.07e-14 로 반증됐다. 계약은 한 문장이 됐다:
   *"연산이 만든 모든 정점: `vertex_point(v)` 는 `realize_cache` 의 답과 비트 동일이거나, 그 길이 이름으로 거절한
   것"* — census 가 매 행에서 단언한다. (여기 적혔던 «재생(=`refine_caches`)» 은 지어지지 않았다 — 덮어쓰기가 아니라
   태어날 때 실현.)
   ⚠ **정정(2026-09-17, 칸 55 감사): «덮어쓰기가 아니라» 는 이제 반쪽이다.** 칸 54 가 `refine_vertex_cache` 를
   세웠다 — 계획됐던 전량 `refine_caches` 는 여전히 안 지어졌지만, `Ceiling` 만 `Bounded` 로 올리는 **제한된 둘째
   쓰기 문**은 섰고 그것은 덮어쓰기다(내보내기 직전에 부르는 «비싼 문»; 뒤이은 연산은 캐시를 읽어 진실이 되는 것을
   정하므로 그 뒤로는 연산하지 않는다). ⇒ **이 항목에 남은 것은 «버리고 재생» 보증**이지 «덮어쓰기가 전혀 없다» 가
   아니다. 그래서 ⏳ 로 둔다 — 근거는 다 지나갔고 사실 하나가 남았다.
14. **혼합-프레임 정점은 reuse 가 답하지 못한다(S7 기록)** — 호출자가 세계 좌표로 명시한
   밑캡 위의 프레임-스케치 코너는 세 담체의 모션이 갈려(둘은 프레임, 하나는 세계) 유리수
   pullback 이 없다. 옛 `Origin` 길은 base 정점으로 답했지만 그 정점이 S7 에서 소멸했다.
   실측 인구는 `the_def_road_answers_for_the_populations_it_can_name` 의 ④ 가 핀한다.
   ✔ 칸 52 의 실현 거절 인구(같은 부류): 불리언 결과 정점 4,924 중 `NoMeet` **346**(회전 피연산자가 섞인 결과 —
   두 모션 이력·Wide 이름), 피연산자 6,988 중 `NoMeet` 192 · `NoCurvedPoint` 12(회전된 원통의 seam: 세상에 진술
   못 하는 담체) · 반사된 원통(`world_cylinder_def` 는 순수 이동만 안다) — 전부 구성 폴백으로 서고 census `r` 행이
   센다. 실현 문이 넓어지면 이 수가 준다.
10. **감김(`oriented_ring`)은 아직 f64 다** — 배치된 3D 점의 면적벡터·법선 내적(ops).
   `check()` 가 단순성(≠0 면적)을 진실 위에서 보증하므로 지금은 건전하지만, f64 폴백 소멸
   (S6 이후)과 함께 재검할 것.

22. **`tess` 는 «수치 정밀도» 손잡이가 없다 — 내보내기가 그것을 필요로 한다** (2026-09-10, 논의만;
    ⚠ **번호 정정 2026-09-12 — 이 항목은 «옛 이름»이 있다.** `15.` 로 붙어 있었는데 그 번호는 이미
    쓰이고 있었다(「술어가 이름의 정수를 읽는다」, 2026-08-09 닫힘). 22로 옮긴다.
    ⇒ **「열린 항목 15」라고 적힌 기록은 두 항목으로 갈린다**(실측: 인용 21곳):
    `design.md`·소스 주석 셋(`planes.rs`·`tolerant_tests.rs`·`predicate.rs`)·dev-log 7456·7829 는
    **닫힌 15**를, **dev-log 19734·19752·19806·19811·19835(칸 ㊵)는 «이 항목»**을 가리킨다.
    dev-log 는 그때 참이었던 것을 적는 기록이므로 고치지 않고, 여기에 별칭을 적어 둔다 —
    **칸 ㊵ 기록의 「열린 항목 15」 = 이 항목(22)**.;
    ★ 2026-09-11 칸 ㊵ 로 **정점 쪽 절반이 섰다** — 아래 「칸 ㊵ 가 한 것」).
    ⚠ **두 손잡이를 가르는 것이 이 항목의 전부다.**

    | | 무엇 | 오늘 |
    |---|---|---|
    | **A. 세밀도** | 삼각형 개수 — `TessConfig { tol, max_angle_deg }`. `tol` 은 현 편차(sagitta), `max_angle_deg` 는 「2° 면 800px 원이 0.06px 안」 — 순전히 **보기**의 수 | ☑ **있다** |
    | **B. 수치 정밀도** | 각 좌표를 몇 비트로 **실현**하나 | ✗ **없다** |

    ★★★ **칸 58 의 실측 — 평면 캐시의 «자리»는 B 를 기다리지 않고 닫혔다** (2026-09-18).
    캐시의 앵커는 이제 **진실에서** 나온다(진실 삼중의 첫 점). 인구와 값:

    | 잰 것 | 수 |
    |---|---|
    | 유도가 성립 / 거절 | **5,317 / 532**(거절 = 모션 456 · 이동된 원통 76; **넘침 0 · 이름 없음 0 · wide 0**) |
    | 앵커가 **이미** 진실 첫 점 | **4,678** |
    | 진실의 *다른* 점에 앉아 있던 것 | **45** |
    | 셋 중 **어느 것도 아닌** 점에 앉아 있던 것 | **34** |
    | 실제로 캐시 비트가 바뀐 것 | **61** |
    | census `c` 행: 결과 이동 / 피연산자 이동 | **0 / 10**(debug = release) |

    ⚠ **B(자릿수)가 여는 것은 «방향»이다** — 평면의 단위 법선, 원통의 축과 `ref_dir` 은 전부
    `normalize()`(`self / n2.sqrt()`)를 타므로 유리수가 아니고, STEP 이 내보내는 값이 바로 그것이다.
    ⚠ 그리고 «행(`raw`)» 은 아직 생산자의 것이다: 저장된 `raw` 가 «실현된 진실 점들의 외적»과 같은 것은
    **2,747**, 다른 것이 **2,010**(회전 사슬 456은 답할 수 없음) ⇒ 행을 진실에서 가져오는 칸은
    앵커의 61 과 달리 **약 2,010개를 움직이는 큰 변경**이다.

    실측: `nacre-tess/src/lib.rs:480` 이 `model.vertex_point(v)` 로 **f64 캐시**를 읽는다 —
    그 캐시는 자기 오차를 `PointCache::Measured { residual }` 로 들고 있다(즉 **최근접 f64 가 아닐 수 있다**).
    ⇒ 화면에는 충분하지만, 내보내기(OBJ·STEP)는 **고정밀 실현 → 한 번만 반올림**을 원한다.

    ☑ **기계는 이미 있다**: `nacre-scalar` 의 `HpBounded`(`add`/`sub`/`mul`/`div`/`inv_sqrt` 전부
    `prec: usize` 를 받는다) · `WitnessPoint::hp_coord(prec)`(**`pub(crate)`**) ·
    `judge_precision`/`trial_bound`(정밀도 «고르기»는 `pub`). 그리고 «정확 반올림이라 실현이
    유일»하므로 **같은 정밀도면 답이 하나**다 — 결정적이고 오라클을 갖는다.
    ✗ **없는 것**: `HpPointCache`/`HpSurfaceCache` 라는 «층»(이 문서 §캐시의 도착점) 과, `tess`
    가 좌표를 캐시가 아니라 **실현에서** 받는 길.

    ★★ **가장 값 하는 곳은 정점이 아니라 «곡선 위 샘플점»이다** — 결론은 서지만 ⚠ **여기 적었던
    근거는 2026-09-11 칸 ㊵ 의 감사에 반증됐다.** 「폭 7비트라 f64 에 정확히 든다」가 아니다:
    좌표별 실측으로 `boolean_corner` 의 48개 중 **32개만 dyadic** 인데(3.3·7.7 로 자르니 `33/10`
    은 이진 분수가 아니다) 캐시는 **48개 전부** 최근접이다. 표현 가능성은 기제가 아니었다.

    실측으로 다시 세운 근거: **축정렬·불리언 인구의 정점은 안 움직인다**(48/48 · 96/96 — 두 길이
    같은 유리수를 같게 반올림하니 어긋날 수가 없다). **기울어진 프레임 인구는 움직인다**
    (36 좌표 중 29 — 그쪽 캐시는 더 긴 f64 유도라 아무것의 정확 반올림도 아니다).
    그래도 원·원통을 잘게 나눈 점은 `cos/sin` 평가라 정밀도가 그대로 드러나고 **OBJ 는 거의 전부 그
    샘플점**이므로, 「값은 샘플점에 있다」는 그대로다.

    ⚠ **순서**: B 가 없으면 「STEP 텍스트에 f64 보다 많은 자리」는 **거짓말**이다(f64 의 17자리
    뒤는 이진 반올림의 부산물). ⇒ **B → 그 다음 자릿수 지정.** 거꾸로 하면 없는 정밀도를 출하한다.

    ☑ 그리고 **잠금이 공짜로 딸려 온다**: 실현이 유일하므로 ⑴ 같은 모델은 **바이트 동일**한 파일을
    내고 ⑵ 실현 정밀도를 **두 배로 올려도 f64 결과가 안 바뀌어야** 한다(이미 최근접이므로).
    그 둘이 「이건 정밀도 «기능» 이 아니라 정확성 «고침» 이다」의 증거다.

    #### 칸 ㊵ 가 한 것 (2026-09-11) — **정점 쪽 절반**

    `nacre_ops::realize_vertex{,_decimal}` + 앱의 `v.digits(n)`. 정의에서 실현하고 한 번만
    반올림한다. **캐시는 한 비트도 안 바꿨다**(덮어쓰기는 3b 의 남은 절반).

    ★★★ **그리고 이 문서가 「원리로만」 적어 둔 것이 정점에서 실측됐다.**

    | 인구 | 정점 폭 | 대조 | **캐시가 최근접** |
    |---|---|---|---|
    | 축정렬 상자 · 두 번 불린 | 7비트 | 48 | 48/48 |
    | 기울어진 프레임 위 프리즘 | **59비트** | 12 | **0/12** — 전부, 최대 **4 ulp** |

    ⚠★★★ **가르는 것은 «표현 가능성»이 아니다 — 첫 정리에서 내가 틀렸고 실측이 반증했다.**
    좌표별로 재니 `boolean_corner` 의 48개 중 **32개만 dyadic** 이고(픽스처가 3.3·7.7 로 자르는데
    `33/10` 은 이진 분수가 아니다) **캐시는 48개 전부 일치**한다. ⇒ f64 로 표현 불가한 것이
    어긋남의 원인이 아니다.
    실제로 가르는 것: 좁은 인구는 **두 길이 같은 유리수를 같게 반올림**하니 어긋날 수가 없다.
    기울어진 인구의 캐시는 더 긴 f64 유도에서 나오고 그건 무엇의 정확 반올림도 아니다 —
    **59비트라는 폭은 «산술이 얼마나 있었나»의 대리 지표이지 기제가 아니다.**

    ☑ 새로 선 것: `round_to_digits`(`round_to_f64` 의 십진 쌍둥이 — 못 정하면 «모른다») ·
    `decimals_of_ratio`(유리수는 긴 나눗셈, 찍는 자리가 좌표 자신의 자리) · `nearest_f64_big` ·
    곡선 정점용 유계 산술(`sqrt_bounded`·`realize_quad`·`realize_seam_point`) ·
    `WitnessPoint::realize(prec)`(cip 의 유일한 새 공개 표면, 새 타입 0개).

    ⚠ **고친 벽 하나**: `rad_upper_big` 이 `2^e` 를 f64 로 만들며 `|e| > 1000` 을 거절해서,
    **1024비트 실현이 「모른다」로** 나오고 있었다 — 사다리를 오를수록 답이 나빠졌다. BigFloat 에서
    조립하도록 고쳤다.

    ⚠ **남은 것은 이 항목의 본체**: `tess` 의 손잡이 B 와 **곡선 위 샘플점**. 위 표가 말하듯 정점은
    대부분 안 움직이고, OBJ 는 거의 전부 샘플점이다.

    #### 설계 (2026-09-10 이어진 논의로 확정된 모양) — ✔ 칸 52 (2026-09-15) 가 지었다, 모양은 셋이 달라졌다

    ✔ **지은 모양 — 정제는 없다, push 가 실현한다.** 아래가 그린 명시적 `refine_caches` 는 짓지 않았다: ops 의 push
    깔때기 `push_vertex_realized(model, def, fallback)` 가 정의를 **첫 단(128비트)** 에서 실현해 `Bounded` 로 넣고,
    첫 단이 f64 를 못 정하거나 실현이 이름으로 거절하거나 **모션 깊이가 64 를 넘으면** 구성의 폴백을 든다
    (`realize_cache` — 계기가 같은 물음을 다시 물어 캐시를 잡는다). 아래의 전제 셋이 실측에 뒤집혔다:
    ① *"128비트 실현은 훨씬 느리다 ⇒ 명시적 연산"* — 정점당 3~4 µs(코퍼스 4,924개에 17 ms). 다만 **깊은 모션
    이력은 예외**다: 4,200회 회전 픽스처가 push 마다 사다리를 오르자 5.6초 테스트가 분 단위가 됐다 — 그래서 첫 단
    + 깊이 64 (오차가 회전당 ~1비트라 그 너머는 첫 단이 어차피 못 정한다). ✔ 그 64 는 매직 넘버가 아니라 **두 일을 하나로 묶은 것**이었다 — 칸 54 가 갈랐다(항목 30). 불리언 벤치는 **+~10%**(fold 80 1.00 →
    1.11 s) — 열린 항목 30. ② *"정점만 하면 파일이 자기모순"* — 실현값↔담체 캐시 거리 최대 **1.07e-14**, validate
    ε 의 다섯 자릿수 아래. ③ *"238건이 어디로 가는가"* — 칸 ㊵(축정렬 48/48 동일)와 칸 52 표(527 이동 ≤ 27 ulp).
    ★ 그리고 사용자가 잡은 것: «연산 끝에 덮어쓴다»(둘째 초안)도 필요 없다 — 정의에서 실현할 수 있으면 태어날 때
    넣으면 되고, 캐시는 append-only 로 남는다. 간선 캐시는 끝점 뒤에 서므로 재유도도 없다.

    ⚠ **첫 정리에서 내가 두 번 틀렸고 사용자가 둘 다 잡았다.**
    - *「캐시를 고치면 기하가 망가진다」* → **틀림.** 판정은 캐시를 **1단 f64 필터**로만 쓰고
      (규칙 5) **증명은 진실에서** 나온다 ⇒ **더 정확한 좌표 + 정직한 `tol` 이면 필터가 좁아질
      뿐 답이 나빠질 수 없다**(상승도 준다). 3b 가 보류인 근거는 *"**naive** re-solve 와
      238/1,992 가 다르다"* 인데, **정확 반올림 실현은 naive 재풀이가 아니다**. `OnSeam` 의 doc
      자신이 *"a unique point, **exactly designated** … 정의는 완성됐고 없는 건 재생 기계"* 라
      적었다 ⇒ **불건전해서가 아니라 아무도 짓고 검증하지 않아서** 보류다.
    - *「출력 전용 실현이 유일한 안전한 길」* → **불필요.** 캐시는 아레나가 아니라 덮어써도 되고,
      덮어쓰는 편이 **덤이 크다**(화면·tess·STEP·OBJ 가 한 값을 공유 ⇒ 서로 어긋날 자리가 없음).

    ★ **모양**: 정제를 «내보내기 안»이 아니라 **자기 연산**으로 둔다.
    ```rust
    model.refine_caches(bits);                     // 명시적 — 부르는 쪽이 정한다
    let step = to_step(&model);                    // 서명 무변 (오늘 둘 다 &Model — 실측)
    let obj  = tessellate(&model, cfg)?.to_obj();  // 무변
    ```
    ⚠ 내보내기가 `&mut Model` 을 받으면 **이름이 약속하지 않은 일**을 한다(그 뒤 불리언의 판단이
    달라질 수 있다). ⚠ 그리고 128비트 실현은 f64 캐시보다 **훨씬 느리다** — 화면 갱신마다 부르면
    안 된다. 그 비대칭이 «명시적 연산» 이어야 하는 둘째 이유다.
    ⚠ **세 캐시 전부**: 정점 · 간선(`derive_edge_curve` 가 «carriers and endpoints» 에서 유도 ⇒
    `rebuild_edge_cache` 가 이미 있다) · **곡면**(`SurfaceCache.coeffs` — STEP 의 평면 계수가
    여기서 나온다). 정점만 하면 **파일이 자기모순**이 되고 OCCT 오라클이 *"vertex not on face"* 로 잡는다.

    ★★ **공개 문은 하나, 정밀도는 «이름 있는 인자»로.**
    ```rust
    pub enum Precision { NearestF64, Bits(usize) }   // 나중에 Digits(n) 도
    impl Model {
        pub fn realize_vertex(&self, v: Handle<Vertex>, p: Precision) -> Realized;
        pub fn realize_surface(&self, s: Handle<Surface>, p: Precision) -> RealizedPlane;
    }
    pub struct Realized { coord: [BigFloat; 3], tol: [Mag; 3] }  // 값+tol 한 덩이 (규칙 4)
    ```
    ⇒ 셋이 **같은 문의 소비자**가 된다: `refine_caches` 는 `NearestF64` 로 부르고, 고정밀 STEP 은
    `Bits(n)` 로 불러 **f64 를 안 거치고** 자릿수를 찍고, «점 하나만 아주 정밀히» 도 그 자리에서 열린다.
    ★ 그러면 캐시가 «별개의 진실» 이 아니라 **실현의 «메모»** 가 되고, `refine_caches` 의 계약이
    한 문장이 된다 — **«모든 `vertex_point(v)` 가 `realize_vertex(v, NearestF64)` 와 같아지게 한다»**
    — 그 문장이 곧 잠금이고 멱등도 거기서 따라 나온다.

    ⚠ **이름 둘을 재서 골랐다.** `_at` 은 안 된다 — 이 커널에서 `at` 은 **위치의 낱말**이다
    (`point_at` ×12 · `normal_at` ×6 · `surface_handle_at` ×4 …). 반대로 `realize_*` 는 이미
    «정의를 정밀도로 값이 되게 한다» 는 뜻으로 쓰이고 **전부 `prec` 를 받는다**
    (`realize_cos_sin(prec)` · `realize_inv_sqrt(v, prec)`). ☑ 판정과의 혼동은 **암묵적 기본값**이
    원인이었지 낱말이 아니다 — `Precision` 이 항상 명시되면 기본값이 없다. 판정은 이 문을 안 부르고
    자기 어휘(`hp_coord` · `judge_precision` · `trial_bound`)를 쓴다.

    ★★★ **먼저 잴 것 — 짓기 전에.** *"238 of 1,992 differ from a naive re-solve"* 의 **238건이
    정확 반올림 실현으로는 어디로 가는가**. 같으면 정제는 거의 무동작이고 얻는 것은 곡선 샘플뿐이며,
    다르면 **어느 쪽이 옳은지가 드러난다**(진실에서 나온 쪽이 옳다) ⇒ 그 자체가 버그를 찾는 계기다.
    그 수가 이 일의 크기와 값을 정한다.

25. ⏳★★★★(1단계 ✔ 2026-09-15) **원/원호의 진실은 «실현값»이 아니라 «정의»여야 한다 — 평면의 normal-vs-coefficients 와 같은 갈래**
    (2026-09-13 진단, 곡선 마일스톤 입력).

    오늘 원호는 «실현된 값»으로 저장된다: `Edge2d::Arc{center: [Rat;2], radius: Rat, ccw}`(2026-09-15 까지 `Seg2d`) ·
    `CylinderDef{.., radius: Rat}`. 그래서 **무리수가 되는 두 양**에서 막힌다.

    **① 반지름 — 저장하는 양이 틀렸다.** `radius_of` 는 `r² = |start−center|²`(분수끼리라 **항상 유리수**)를
    구한 뒤 `rat_sqrt_exact(r²)` 로 √ 이 유리수인지 보고, 아니면 `ArcRadiusNotRational` 로 거절한다. 그런데
    정확 술어는 r 을 **제곱해서** 쓴다(원통 게이트 `(n·o+d)² > r²|n|²`, 부피 `πr²h`) — **맨 r 을 유리수로
    요구하는 곳은 `radius_of` 자신뿐**이고, geom `Circle` 은 r 을 **f64**(캐시)로 든다. ⇒ 진실은 `r²: Rat`
    이어야 한다(항상 유리수 = 모든 원을 담는다), r 은 f64 캐시로 √ 해서 쓴다. 이것은 평면이 「유리수 법선을
    저장하려다 무리수·오버플로에 걸려 계수(정의)를 저장하고 법선을 유도」한 것(`normal_def` 재구성, S6b)과
    **같은 실수**다. 막는 인구: 손으로 그린 «비피타고라스» 원호(중심 (0,0), (1,1)→(−1,1), r²=2).

    **② 중심 — r² 트릭으로도 안 풀리는 둘째 축.** 기울어진 코너의 필렛은 **중심 자체가 무리수**다(모서리에서
    이등분선 방향으로 r — 단위 이등분선이 무리수). kit 이 *"축정렬 직각 코너가 아니면 거절"* 하는 진짜 이유.
    이것도 평면과 같다 — 기울어진 평면이 무리수 법선이라 «점 셋(정의)»을 저장하고 실현하듯, 기울어진 원호도
    «만든 방법(두 간선 + 반지름)»을 저장하고 중심을 실현해야 한다.

    ⇒ **뿌리는 하나**: 원/원호가 평면이 S6/S7 에서 지나온 «정의 vs 실현» 분리를 아직 안 지났다. 「다 만들 수
    있어야 한다」가 옳고, 길은 평면이 간 길이다. ⚠ 크기는 doc 한 줄이 아니라 **곡선 마일스톤**(진실 필드
    `radius`→`r²` 를 `Edge2d::Arc`·`CylinderDef` **둘**에서 함께(항목 26 이 조각 타입을 하나로 만들어 셋이 둘이 됐다) + 중심의 정의 기반 표현 신설) — 지금 고치는 게
    아니라 그 계획의 입력이다. ★ `arc_turns`/`arc_rat` 이 갈린 것은 이것과 무관한 편의 설탕(끝점을 90°k 로
    유도 vs 받음)이고, 정의 기반이 되면 둘 다 «정의를 진술하는 한 방법»으로 자연히 정리된다.

    ★★★★ **2026-09-15 설계 방향 (사용자 검토, 칸 ㊿ 뒤 논의) — 원의 타입은 하나, 넓히는 것은 «수의 타입».**
    r² 는 첫 계단일 뿐이다: 중심이 유리수이고 r² 가 유리수인 원만 담기고, 임의 각 코너의 필렛(중심 유도),
    두 호에 접하는 필렛(중심·접점·때로 반지름까지 유도)은 담기지 않는다. 그 답이 «필렛-두-직선»·«한쪽만
    접함»·«두 점 + r²» 같은 **정의 변종의 열거는 아니다** — 구속이 하나 늘 때마다 커널 변종이 하나 느는 벽이다.

    - **다른 CAD.** Parasolid·ACIS·OCCT 는 원을 «중심 f64 + 반지름 f64» 한 타입으로 들고 정확 산술이 없다.
      유도는 커널 밖에서 한다 — 스케치 구속 해석기(SolveSpace·PlaneGCS·Onshape 자체 솔버: 접선·일치·치수를
      방정식으로 세워 뉴턴 계열 수치 해)나 작도 명령(AutoCAD `CIRCLE TTR`·`FILLET`: 오프셋 원 두 개의 교점을 f64
      로). «접한다»는 1e-9 안에서만 참이라 b-rep 이 정점·간선마다 허용오차를 든다(tolerant edge) — 우리가 DNA
      로 거부한 «구성 시점의 tolerance»가 그 대가다. 정확 진영(CGAL 원형 커널, LEDA real, CORE Expr)은 «식을
      저장하고 필요한 정밀도로 구간 평가, 분리 한계로 0 판정» — 이 커널이 정점에서 이미 하는 것과 같은 진영.
    - **핵심 사실: 직선·원 작도는 제곱근으로 닫혀 있다.** 접선·필렛·아폴로니우스는 전부 «오프셋한 직선/원의
      교점»이라 직선끼리는 사칙연산, 원이 끼면 제곱근 **하나**. 유도된 호 위의 필렛은 근호가 중첩된다(차수 2·4·
      8…). 그래도 «자와 컴퍼스로 작도 가능한 수» 한 부류다. 타원·스플라인은 밖(M7 이 그은 선).
    - **저장 타입.** 원은 끝까지 `{ center: [N; 2], r2: N }` 이고, `N` 이 `Rat` 에서 «작도된 수의 정의»로 넓어진다:
      `NumDef { Rat, Add, Sub, Mul, Div, Sqrt }` — 값이 아니라 **정의(식)**, 실현은 `HpBounded` 구간 캐시(정점의
      `Vertex` + `PointCache` 와 같은 분리). `a + b√c` 는 이 식의 깊이 1 특수형이고, 필렛 위의 필렛은 깊이 2 —
      **2단계와 3단계의 저장 타입은 같다.** ⚠ «쪽»(두 교점 중 어느 것)은 정의의 일부다 — 두 중심을 잇는 방향의
      좌/우, 유리수 부호 하나.
    - **경우의 수는 라이브러리에 있다.** 두 직선 접 = `meet(offset(l₁,r), offset(l₂,r))` · 직선 접 + 점 통과 =
      `meet(offset(l,r), circle(p,r²), 쪽)` · 두 원 접 = `meet(circle(c₁,(r₁±r)²), circle(c₂,(r₂±r)²), 쪽)` · 두 점 + r =
      `meet(circle(p,r²), circle(q,r²), 쪽)` · 세 점 = 수직이등분선 둘의 교점(근호 없음) · 접점 = 수선의 발 또는
      `c₁ + r₁/(r₁+r)·(x−c₁)`(근호가 새로 안 생긴다). 「필렛」「접선 호」「TTR 원」은 이 두 연산(오프셋·교점)을 부르는
      **함수**이지 커널 변종이 아니다 — «커널의 문은 고리, 펜은 kit» 과 같은 선.
    - **2단계 ↔ 3단계의 실제 차이 = 부호 판정 엔진.** 커널이 수에게 묻는 건 «이 식의 부호» 하나다.
      (2) 근호 하나: `X + Y√c` 는 제곱해 유리수끼리 비교 — **M6-1 부호 탑, `Pierce`/`QuadRoot` 가 이미 쓴다**,
      새 엔진 없음. (3) 중첩: 제곱해 없애면 깊이마다 차수가 두 배라 폭발 ⇒ 구간 평가 + 정밀도 상승(cip 의
      에스컬레이션 그대로) + **멈추는 규칙 = 분리 한계**(식의 차수·계수 크기에서 «0 이 아니면 |x| ≥ 이만큼»,
      BFMSS) — 오늘 cip 가 `Undecidable` 로 돌려주는 자리가 «0 임이 증명됨»이 된다. 그 앞에 **구조적 0**: 필렛의
      접점은 «수선의 발»로 정의됐으니 직선 위에 있음은 증명이지 측정이 아니다(접하는 두 원의 판별식 0 도) —
      «구조적 0 은 증명되어 온다»가 이미 있어 분리 한계까지 가는 일은 드물다.
    - ⚠ **경계**: 담기는 것은 «순차 작도»(앞서 만든 것에서 다음이 결정됨)로 풀리는 구속이다. 여러 구속을 **연립**
      으로만 풀 수 있는 경우는 3차 이상이 나올 수 있고 그건 제곱근 밖 — 이름 붙여 거절하거나, 훨씬 뒤에 일반
      대수적 수(다항식 + 고립 구간)로. 구속 해석기(«치수를 주면 좌표가 계산되도록»)는 이 위의 프런트엔드로,
      구속을 «작도 순서»로 풀어 `NumDef` 식을 내놓는다 — 닫힘 덕에 «해석기의 답을 커널이 정확히 담는다»가 보장.
    - ⇒ **계단 셋**: (1) r² 진실화(닫힌 형태 유지, 비피타고라스 호 즉시 해방) → (2) `NumDef` 를 세우되 깊이 1 만
      허용, 기존 부호 탑으로 판정(유리수 직선·원에 접하는 필렛 전부) → (3) 같은 타입에서 깊이 제한을 풀고 구간
      정제 + 분리 한계. 저장 타입은 (2) 에 한 번 정해지고 그 뒤 바뀌지 않는다. ⚠ 열린 비용: 같은 원을 다른
      정의로 두 번 만들었을 때의 동일성(정의 해시로는 다르고, 대수적 동일성은 `sign(a−b)=0` 판정 — 인턴은 정의
      기준으로 두고 기하 동일성은 불리언의 일치 판정에 맡기는 것이 평면과 같은 길).

    ✔ **1단계 닫힘 (2026-09-15, 칸 51).** `Edge2d::Arc { center, r2: Rat, ccw }` · `CylinderDef { origin, dir, ref_dir,
    r2: BigRat }` · `ArcSpec.r2`·`QuarterArc.r2`·`Corner::Round.rho2`. 스칼라 문 12 가 r² 를 받고(`&BigRat`), 그중
    반지름의 합/차를 쓰던 셋(`cylinders_clear` 평행 팔·`cylinders_nested`·`skew_axes_clear`)은 이미 있던
    `sqrt_root_sum_cmp`(√a > √b + √c ⟺ a−b−c > 0 ∧ (a−b−c)² > 4bc) 로 — 근호를 만들지 않는다. `realize_seam_point` 와
    신설 `sqrt_f64`(`inv_sqrt_f64` 의 쌍둥이)는 «완전제곱이면 옛 길(비트 동일), 아니면 √ 실현». 스케치 문 `r2_of` 는
    dist² 를 그대로 — `ArcRadiusNotRational`·`radius_of` 소멸. `radius()` 삭제로 컴파일러가 **58 자리**를 짚었고(제품 49 + 소스 내 테스트 9 — 플랜의 «38» 은 표를 더한 어림이었다;
    캐시의 f64 `radius()` 둘은 남는다), 같은 `Rat` 타입이라 못 짚은 자리 **하나**(nesting 의 림 증인 단언이 정확 반지름을 r² 문에 넘김)는 관문이 잡았다 —
    디버그 census 행 수 398 → 160 이 첫 신호.
    ★★ **플랜의 «폭» 주장은 틀렸었다 — 분자만 셌다.** 16자리 십진 반지름은 ~1e-4 아래에서 r² 의 **분모**가 i128 을 넘는다
    (5.000000000000001e-8 → 분모 10²³, 제곱 10⁴⁶). 그 인구는 이미 잠겨 있었다(`a_bored_cube_builds_at_any_size…` 의
    17자리 가수 일곱 자릿수, `reject_census` 의 `cylinder_wide_axis`) — 「폭은 원인이 아니다」의 잠금들. 사용자 결정:
    창을 좁히지 않고 **넓은 타입**으로 — `nacre_scalar::BigRat`(BigInt 유리수, `MeetPoint::Wide` 의 선례;
    `square_of(Rat)` 은 절대 넘치지 않는다). 스케치의 `Edge2d::Arc.r2` 는 `Rat` 그대로(dist² 가 오늘도 `Rat` 이어야
    문을 지난다). ops 에서 `Rat` 산술을 하던 자리 넷(`arc_extent`·격자 탐침의 λ·`t_cap`·`ArcSpec`)은 `narrow()?` 로 —
    옛 `r·r` 이 넘치던 **바로 그 자리**에서 같은 거절. 「유리수 점을 r 로 짓는 자리 여섯」은 `radius_exact()`
    (`rat_sqrt_exact_big`) 뒤에 그대로 서 있다 — 2·3단계의 방문 목록.
    ☑ 잠금: census 398행 HEAD≡debug≡release 비트 동일 · `a_non_pythagorean_arc_is_stated_and_extruded`(활꼴 부피
    π/2−1, 옆면 캐시 √2 비트 동일) · **상자 `cut` 탐침이 초록이라 잠금으로 승격**(`a_non_pythagorean_prism_is_cut_by_a_box`
    — 불리언 길이 r²=2 를 끝까지 지난다) · `sqrt_f64` 완전제곱 250 + 비완전제곱 7 + 넓은 제곱 1 · 비완전제곱 반지름
    쌍(r²=2, 3)의 `cylinders_clear/nested` · 넓은 제곱의 솔기 실현 · kit 의 «√2 반지름 호» 거절 테스트가 수용 잠금으로.
    죽은-이름 스윕: 예측 «은퇴 서술 둘» = 실제.

29. ⏳★★ **pad/pocket 은 kit 의 설탕으로 — 커널은 extrude + boolean 만** (2026-09-15, 칸 52 의 메모).
    design.md §6 「패드/포켓 ↔ 불리언 통합 ✅」대로 커널의 `PadOnFace`/`PocketOnFace` 는 이미 «면 프레임 + 압출 +
    Fuse/Cut» 래퍼이고 새 정확 술어가 없다 — overview 의 설탕 판별 기준(«편의 레이어는 커널 op 을 조합만 한다»)에
    그대로 걸린다. kit 이 두 변종을 부르는 곳은 `build.rs` 한 자리. 제거는 `Operation` 어휘(저장된 로그의 재생)와
    kit 재작성을 함께 재야 하므로 **별도 플랜**.

31. ⏳★★ **픽스처는 제품 도로로 — 실현 통로를 안 지나는 테스트 전용 구성 문을 은퇴시킨다** (2026-09-16 의논으로 확정,
    칸 53 뒤).
    실측: topo 의 `test-util` 문 `add_cuboid` **804곳(51 파일)** · `add_cylinder` **182곳(23 파일)**. topo 는 ops 를 모르니
    둘 다 실현 없이 `Unmeasured` 로 push 하고, census 의 불리언 코퍼스 전부가 이 상자·원통으로 서 있다 — 존재 이유는
    «topo 층 테스트가 ops 없이 솔리드를 원했다»와 «스케치 도로보다 코퍼스가 먼저 있었다». (`Model::new()` 의 씨앗은
    세계 평면뿐, 정점은 없다.)
    ✔ 모양: 테스트 지원 모듈에 **같은 서명의 `cuboid(m, min, max)`·`cylinder(m, …)`** 를 두고 안에서 `apply(Extrude)`
    (사각형·온전한 원 프로파일)를 부른다. 모든 픽스처가 push 깔때기를 지나 `Bounded` 가 되고, «깔때기를 안 지난 push»
    라는 `Declined` 사유가 제품에서 사라지며, census 가 «제품 도로의 census»가 된다. topo 의 `add_cuboid`/`add_cylinder`
    는 은퇴(스윕).
    ★ **남기는 것 하나 — 심기용 날것 문.** validate 의 자기 테스트는 **틀린 모델**(매달린 핸들·뒤집힌 면·정의와 어긋난
    캐시)을 심어야 하고 제품 도로로는 틀린 모델을 만들 수 없다. 그래서 `push_vertex`(날것 캐시)·`push_plane_unregistered`
    같은 문은 `test-util` 아래 topo 에 남는다 — 제품에 없는 것을 테스트용으로 두는 게 정당한 유일한 경우.
    ⚠ 위험은 census 가 움직이는가: 압출 상자의 평면 이름은 `add_cuboid` 와 같은 정준형이라 **비트 동일이 예측**이고,
    원통은 `add_cylinder` 의 seam 모형과 압출 원의 seam 이 같은지 그 칸의 A/B 가 답한다. 986곳 치환은 같은 서명이라
    기계적.

33. ⏳★ **가라앉은 좌표는 경계도 0 으로 실린다** (2026-09-16, 칸 54 에서 찾고 미룸).
    `nearest_f64_big_exact` 가 «가라앉은 값을 정확하다고 말하던» 것은 칸 54 가 고쳤다(`exact = false`). 그런데
    `Realized::to_f64` 의 오차 팔이 `Mag::of(val).times(Mag::pow2(-53))` 이라, `val == 0.0` 이면 **경계도 0** 이다
    ⇒ 0 이 아닌 값에 «참값 ∈ 0 ± 0» 이 실린다. 깃발이 거짓말을 그만둬도 **경계는 여전히 거짓말한다.**
    ✔ **인구 0**: 그러려면 좌표가 2⁻¹⁰⁷⁵ 아래여야 하는데 `Rat` 이름은 i128 비율이고 가장 넓게 잰 `Wide` 이름도
    168비트다. ⇒ 고치려면 `to_f64` 의 오차 모형을 손대야 하고(가라앉은 값의 경계는 최소 비정규수의 반),
    그것은 칸 54 의 범위 밖이라 적어 둔다. 닿으면 `Unrepresentable` 로 거절하는 쪽이 맞을 수도 있다 — 그때 잰다.

35. ⏳★★★ **캐시의 철자 — `enum SurfaceCache` 와 `.plane()`/`.cylinder()` 체이닝** (2026-09-17, 칸 57 이
    범위 밖으로 냈다). 「문의 이름」이 **2026-09-12 에 확정**했고 아직 안 지었다 — 그 미이행 자체가
    항목 18·27 이 벌하는 «문서가 죽은 구조를 현재형으로 가르친다» 부류다.
    ⚠ 칸 57 이 근거로 삼으려던 «재발 방지»는 **인구 0으로 반증**됐다(항목 34 ②) ⇒ 이 항목의 근거는
    «문서가 확정한 구조의 이행»이지 «막는다»가 아니다. 인구: 값 21 + 테스트 44, 메서드 이전 8,
    통째 출구는 `transform.rs` 한 곳. 테스트 30곳의 종류 판별도 이 칸이 한꺼번에 지나간다.
    ⚠★★★★ **가장 위험한 한 줄은 `flipped_against`** — 캐시 normal 로 interning 의 `flipped` 비트를
    정하고 그 비트가 **스위트에서 1,916번** 발화한다. 잘못 재철자하면 면 방향이 조용히 뒤집힌다.

36. ⏳★★★ **곁표 접기 — `pub surface_name` → `SurfaceCache::Plane.name`** (2026-09-17, 칸 57 이 냈다).
    `Model` 의 **마지막 공개 필드**다(칸 56 이 일곱을 닫고 남겼다). 인구 제품 17 · 테스트 53.
    ⚠ **정정(2026-09-18, 칸 58)**: 쓰는 자리는 이제 `intern_plane` 이 아니라 **`push_plane_raw`**
    (`topo:1224`)다 — 칸 58 이 삽입을 그 문으로 옮겼다(«한 자리»인 것은 그대로고 이름이 다르다).
    이름은 유도된다 —
    *"The name is derived, so it cannot disagree with the thing it names."*
    ⚠★★★★ **유일하게 조용한 실패 모드를 가진다**: 이름이 있어야 할 자리에 없거나 없어야 할 자리에
    생겨도 **census 가 못 본다**(그 인구를 안 든다). 그래서 칸 57 이 합치지 않았다.
    ⚠⚠ **정정(2026-09-18, 칸 58): 아래 문장의 근거가 바뀌었다.** *"진실→geom 유도가 트리에 없다"* 는
    **이제 거짓**이다 — `derive_surface_cache` 가 있고 `push_plane_raw` 가 그것을 적용한다. 재생을
    막는 것은 **향이 아니라 행**이다: 유도가 `raw` 를 `surface_cache(h)` 에서 복사하므로 재생이
    «지우려는 값»을 먼저 읽는다. ⇒ 전체 재생은 여전히 못 서지만 **«앵커만 버리고 재생해도 비트 동일»은
    오늘 형태로도 선다** — 다음 칸의 출발점이다.
    ☑ 향을 진실에서 뽑는 길도 **실측으로 닫혔다**: 비씨앗 `Known` 평면 **342개**가 자기 점 순서와
    반대 방향의 캐시를 단다 ⇒ 점 순서는 방향을 名指하지 않는다.
    ⚠ 근거로 삼으려던 `rebuild_surface_cache()` 는 **`realized` 를 재생할 수 없다** — 진실→geom 유도가
    트리에 없다(`push_plane_raw` 의 doc: *"handed in by the producer today"*). ⇒ 이름만 재생하는
    **`rebuild_surface_names()`** + `surface_name_discard_and_regenerate_bit_identical`, 제약 둘:
    (a) **제품 도로 전용**(`push_plane_unregistered` 는 이름을 건너뛰므로 rebuild 가 이름을 **추가**한다)
    (b) **오름차순 핸들 순서**(`Through` 이름이 담체 이름에 의존하는 자기참조, append-only 라 의존이
    항상 낮은 인덱스). ☑ 선례는 `rebuild_edge_cache` + `edge_cache_discard_and_regenerate_bit_identical`.

37. ⏳★★ **STEP 이 곡면을 캐시에서 통째로 읽는다 — `realize_surface` 가 없다** (2026-09-17, 칸 57 이
    발견). 항목 28 은 «정점»의 정확 반올림 통로만 닫았고, `nacre-step` 은 면의 곡면을
    `surface_cache` 에서 그대로 꺼내 `origin`·`normal`·`radius` 를 내보낸다. 항목 8 은 ✔ 이므로 이
    비대칭은 **주인이 없었다**. 문서 §최종형이 `realize_surface(s, Precision)` 를 약속하지만 코드엔
    없고, 유도에는 `PlaneName{Narrow|Wide}` 경유가 필요하다(`Known` 의 세 `Rat` 에서 법선을 유리수로
    유도하면 넘친다). ⇒ 정점은 태어날 때 실현되는데 곡면은 producer 가 넘긴 f64 그대로다.

    ⚠★★★ **정정(2026-09-18, 칸 58).** ① *"곡면은 producer 가 넘긴 f64 그대로다"* 는 **앵커에 대해서는
    거짓이 됐다** — STEP 이 내보내는 `AXIS2_PLACEMENT_3D` 의 location 은 이제 진실에서 나온다(다만
    `normal`·`ref_dir`·`radius` 는 그대로다). ② *"법선을 유리수로 유도하면 넘친다"* 는 **자리를 잘못
    짚었다**: 이름 유도는 `plane_name_exact` → `plane_name_big` 으로 이미 폭을 넘고, 실제로 넘치는 것은
    **수선의 발**(`plane_origin_projection` 이 `n·n` 을 만들어 계수를 제곱한다 — 이 코퍼스에서 8건).
    칸 58 은 그 길을 **안 쓰기로** 했으므로 이 항목이 기다리는 것은 여전히 **방향의 정확 반올림**이고,
    그것은 항목 22 의 «B → 자릿수» 다.

38. ⏳★★★★ **모션 사슬을 «읽을 때» 접는다 — 무엇이 접히는지 원인별로 쟀다** (2026-09-18 실측, 칸 58
    뒤의 조사; 구현은 아직 없다).
    `Model::chain_translation`(`topo:1717`)은 **첫 비이동 노드에서 사퇴**하고, 그 사퇴가
    `world_plane_name`(`topo:1746`)을 통해 **읽는 곳 8자리(호출 11회)**로 퍼진다 — 캐시 유도의 게이트
    (`topo:1326`) · `vertex_meet` 의 `world_road()`(`topo:2060`) · `world_plane_coeffs`
    (`planes.rs:651`) · `planes.rs:2566`·`:2862` · `transform.rs:623` · `combinatorics.rs:5683`·`:5768`.

    **거절 532 를 사슬의 구성으로 쪼갠 실측**(거절이 나는 `decline_reason` 바로 그 자리에서 분류 —
    모델을 도는 탐침은 push 인구가 아니라 «살아남은 곡면»을 세어 명제 옆을 잰다):

    | 부류 | 평면 456 | 원통 76 |
    |---|---|---|
    | `Frame` 노드를 든다 | **169** | 0 |
    | 임의 각 회전을 든다 | **197** | **7** |
    | 이동 + **사분각** 회전뿐 | **24** | **69** |
    | 위 + **거울**까지 | **66** | 0 |

    169+197+24+66 = 456(정확히 분할된다). ⇒ **오늘 있는 기계로 열리는 인구 = 90 평면 + 69 원통 =
    532 중 159(30%)**.

    ★★★★ **접는 것은 아레나가 아니라 «읽기»다.** `self.motions` 를 만지는 자리는 저장소 전체에
    `push` 와 `get` **둘뿐**(수정·삭제 없음)이고, 노드 핸들은 `SurfaceKey`·`ThroughKey`·`CylinderKey`
    의 **interning 열쇠**라 노드를 합치면 사라지는 것이 자식이 아니라 **동일성 판정 전체**다.
    `chain_translation` 의 doc 이 이미 그 규칙이다 — *"it does not license dropping the history."*
    ⇒ 이 항목은 **노드를 하나도 만들지 않고 하나도 지우지 않는다.**

    ★★★ **산술은 이미 다 있다 — 없는 것은 «노드를 차례로 적용하는 걷기»뿐이다.**
    `Isometry::plane_coeffs`(회전은 `try_exact_cos_sin` 경유, 피벗·이동 포함)·`point_rat`·`dir_rat` ·
    `mirror_plane_coeffs`·`mirror_point_rat`. 모양은 기존 가족을 따른다(`chain_fixes` 라는 private
    걷기 → 공개 얼굴 둘): **`chain_plane_coeffs` / `chain_point` / `chain_dir`** 셋이 한 걷기를
    공유하고, 이는 `Isometry` 자신의 삼총사와 같은 절단면이다.
    ⚠★★★ **반환형은 `Isometry` 가 될 수 없다** — 거울은 `det = −1` 이라 `Isometry` 에 자리가 없는데
    그 인구가 **66**이다. ⇒ `SurfaceDeriveCounts` 의 doc 이 적어 둔 처방 *"a motion wants a chain
    folded to an `Isometry`"* 는 **그 66 에 대해 틀렸다**(칸 58 이 남긴 낡은 처방).
    ☑ 거울의 패리티는 **이름이 흡수한다**: `canonical_plane_coeffs` ③ 이 첫 비영 성분을 양으로
    강제하므로 정준 이름은 방향을 말하지 않는다 — 평행이동에서 오늘 그러는 것과 같다.
    ☑ `chain_dir` 은 **거울을 만날 수 없다**: 방향을 옮기는 것은 원통뿐인데(평면은 계수로 움직인다)
    `OpError::MirrorNotPlanar`(`ops.rs:676`)가 미러링을 먼저 거절한다 —
    `edge_carriers.rs:258` 이 *"A mirrored cylinder has no population"* 으로 이미 적어 뒀다.
    ⇒ **새 스칼라 함수는 필요 없다.**

    ⚠★★★★★ **사용자가 낸 규칙 — 같은 축 회전의 각도 합성. 옳지만 오늘 인구가 0 이다.**
    「사슬은 앞부분만이 아니라 **전체**를 봐야 하고, 같은 축·같은 피벗의 회전은 각도를 더해 접을 수
    있으며, 30°+60°=90° 면 **유리수 표현이 되살아나** 주변 이동까지 함께 접힌다」 — `Angle(Rat)` 이
    도 단위 유리수이고 `checked_add` 가 정확히 mod 360 축약하므로 **타입 수준에서 성립한다**.
    그러나 실측: **`pair_same = 0`**(코퍼스 전체에 같은 축·같은 피벗의 **인접 회전쌍이 없다**),
    **`depth_max = 2`**(`chain_motion` 이 한 호출당 최대 둘(회전→이동)을 남기고, 픽스처가 한 솔리드를
    두 번 회전시키지 않는다) ⇒ **접을 인접쌍 자체가 존재하지 않는다.**
    ⚠ 이것은 «픽스처가 안 한다»이지 «사용자가 못 한다»가 아니다 — 같은 축으로 두 번 돌리면 즉시
    생긴다. ⇒ **관측된 결함의 수선이 아니라 «역량»으로 적어 두고**, 인구가 생기면 짓는다. 위 설계가
    그것을 막지 않는다(노드별 적용 위에 «인접 정규화» 한 겹을 얹는 모양이라 되돌릴 것이 없다).
    ⚠ 같은 축이라도 **피벗이 다르면** 각도만 더해지지 않는다 — 합성 피벗에 cos/sin 이 든다.

    **영구히 밖**: `Frame` **169** — 기저가 `1/√유리수`(`chain_fixes` 의 doc 이 못 박아 둔 경계).
    **열린 항목 2 의 인구**: 임의 각 **197** — Niven.

    ⚠★★★★ **칸 58 과 달리 이 일은 «결과를 움직일 수 있다».** `derive_surface_cache` 가
    `world_plane_name` 으로 게이트하므로 **캐시 확장과 이름 확장은 쪼갤 수 없고**, `vertex_meet` 의
    `world_road()` 가 **자동으로** 넓어져 코너가 새 datum 도로를 얻는다(그 doc 의 *"인구를 열 수는
    있어도 움직일 수는 없다"* 는 계약이 바로 이 자리다).
    ☑★★★ **음성 대조군이 이미 트리에 있다** — `topo:3285` 가 **90°** 사슬에 대해
    *"a rotated carrier has no world name to meet with"* 를 단언하고, 바로 위 주석이 스스로를
    반박한다(*"A quarter turn carries x = 0 to y = 0, so the cache is exact … what declines here is
    the rotation"*). ⇒ **그 줄이 빨개지는 것이 «접기가 작동한다»의 증거**다. 다시 겨눌 때 **37° 짝을
    같은 파일에 남긴다**(무엇이 여전히 사퇴하는지 보여 주는 대조가 없으면 「운 좋은 통과」와 구분되지
    않는다).
    ⚠★★★ **«이름이 생긴다» ≠ «길이 열린다» — 둘째 관문이 있다.** `combinatorics.rs:4451` 의 원통
    관문은 *"refuses any class that is **rotated** or has no narrow rational name"* 이라 회전된
    클래스는 이름을 얻고도 계속 거절된다. `rotated` 는 «세계 진술이 없다»와 «행의 `tri` 가 진실이
    아니라 실현본이다»를 **겸하므로 건드리지 않는다**(`planes.rs:644-650` 이 그 실험을 하고 되돌린
    기록 — 공유 벽이 두 몸으로 갈렸다). ⇒ 인구를 `declined_motion` 의 감소로만 세면 **과대평가**다:
    «이름을 얻은 곡면»과 «실제로 답이 달라진 census 행»을 **따로** 센다.
    ⚠ **폴백 순서는 계약이다**(위 §1296-1299 실측): `world_road` 는 **넓히되 승격하지 않는다** —
    먼저 태우면 오늘 답이 있는 코너의 점 철자가 바뀐다.
