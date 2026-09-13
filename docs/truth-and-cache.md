# 진실과 캐시 — 점·평면·스케치의 최종 타입 구조

> ★★★★★ **타입 구조 확정 (2026-08-05).** 이 문서는 리팩토링의 **도착점**이다 — 최종 타입,
> 숫자 규칙, 남은 이행 항목. 여기 오기까지의 유도·측정·반증의 전체 서사는 git 이력(`a464cd0`
> 이전의 이 파일)과 `docs/dev-log.md`에 있다. **진행 상태는 §이행 의 단계표가 진실이다.**


아레나는 append-only 이므로 **진실만 들어간다.** 캐시는 곁표에 두어 언제든 버리고 재생한다.
정의가 불변이므로 **캐시는 낡을 수 없다** — 무효화라는 개념이 없고, «버리고 재생»만 있다.

---

## 숫자 규칙 — 일곱 개

1. **진실은 «유리수»이거나 «핸들»이다.** `Rat`(i128)에 들어가는 값은 적고, 안 들어가는 값
   (발견된 좌표 160–480비트, 계수 곱)은 저장하지 않고 **가리킨다**. 가리키는 사슬의 끝은
   언제나 수다(C2·C3 의 동시 만족).
   > ★★★ **「160–480비트」는 실측되지 않은 추정이었고, 실측은 훨씬 좁다**(2026-08-07,
   > `nacre-ops/tests/point_width.rs`): 발견된 정점(`vertex_tol().is_some()`)의 담체를 정확히
   > 풀어 **좌표마다 기약분수로 줄인 뒤** 잰 폭이 축정렬 불리언 코너 **7비트**, 회전 **4비트**,
   > 기울어진 십진 프레임(census `wf`) **59비트** — **127비트 초과 0건**이고, 이미 껄끄러운
   > 인구에 기능을 하나 더 얹어도 **59→59로 안 움직였다**(기능마다 평면이 십진수로 새로
   > 진술되므로 풀이가 누적되지 않는다). 코퍼스 수치이지 상한이 아니다(모델 5·정점 108).
   >
   > ★★ **규칙의 근거는 «폭» 이 아니라 «가리키기» 이고, 그것은 이미 서 있다.** (처음 여기에
   > 적었던 근거 —`Pt3` 가 base+chain 을 든다 — 는 **틀렸다**: `Pt3` 는 위 이름 대응표의
   > `WitnessPoint`, 즉 **판정층**(연산 하나 동안만 사는) 타입이지 진실이 아니다.) 진실 쪽은 더
   > 강하다 — 정점은 `Vertex::ThreePlane` 으로 **좌표를 아예 안 들고**(세 핸들뿐), 모션은
   > 규칙 3 대로 **면이 든다**.
   >
   > ★★★ **그래서 «연산이 누적되면 결국 넘치지 않나» 의 답은 «아니오» 이고, 논증이 아니라
   > 실측이다**(`tests/point_width.rs`): 40회를 쌓아도 이름 폭이 **1·4·2비트로 상수**이고, 그
   > 사이 surface 는 246·484 개로 늘지만 **distinct 이름은 6(10)개 그대로** — 움직인 면은 전부
   > **같은 이름 + 모션 노드**다(★ 2026-08-17 한정: 모션이 **고정하는** 평면은 노드 대신 세계
   > 진술을 유지하고 원 핸들로 intern-back 한다 — invariant-plane 재진술; distinct-이름 주장은
   > 그대로다). 셋이 맞물린다: 불리언은 `push_plane` 을 **한 번도 부르지
   > 않고**(결과 면은 피연산자 평면을 재사용), 모션은 계수를 **안 건드리며**, 오프셋은 `n` 을
   > 고정한 채 `d` 만 옮기는데 소수 덧셈의 분모는 곱이 아니라 **lcm** 이다(`1.1 + 6.6 = 7.7`).
   > ⇒ 폭은 **사용자가 쓴 기하**가 정하지 **연산 횟수가 정하지 않는다** — 한 번의 도약이지
   > 누적이 아니고, 그 도약은 `from_decimal` 창이 막는다.
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
   매번 답한다(느릴 뿐, 틀리지 않는다). ★ **이 population 은 16-1 로 실재하게 됐다**
   (2026-08-09): 저장은 **진술 키**(정렬 삼중항+모션, `surface_through_ids`)로 intern 되어
   «같은 진술 = 한 핸들» 은 지켜지고, 기하 동일성은 예고대로 술어의 몫이다. 실측 비용은
   상승 416(narrow 110·wide 454) — 절벽이 아니다.
7. **`Rat` 을 넓히지 않는다(C8).** 넓어지는 것은 중간값(BigInt 이름 유도·`Expansion` 술어·
   BigFloat 상승)뿐이다.

---

## 최종 타입 — 진실

```rust
// ─── nacre-geom / nacre-topo ── 진실 (아레나, append-only) ─────────────

/// ✔ **이름과 자리를 얻었다** (2026-09-11, `78d770f`): 이 enum 이 `Surface` 이고 아레나
/// (`surfaces: Store<Surface>`)가 그것을 든다. 실현은 `surface_cache: Vec<SurfaceCache>`.
/// ⇒ **맨이름 `Surface` = 이것**이고, geom 의 f64 실현은 `nacre_geom::Surface` 로 적는다
/// (열린 항목 8 의 철자 규칙).
///
/// ★ **「`SurfaceDef` 를 흡수한다」의 뜻**(원 주석이 짧아 오해를 낳았다): 그 시절 `SurfaceDef` 는
/// **사이드테이블 enum**(`Constructed` / `Rotated{witness, rotation}` / `Inexact`, `a5379a2` 에서
/// 사망)이었고, 「흡수」는 그 **역할**을 진실이 삼킨다는 뜻이다 — `Constructed`/`Rotated` 가
/// `motion: None`/`Some` 으로 표현된다. **struct 냐 enum 이냐를 정한 문장이 아니다.**
/// ⚠ 그러니 `Vertex { def }` 와의 비대칭은 **설계가 아니라 잔재였다**: `Vertex` 는 최초 커밋
/// (`f86e7e3`)에 `point`+`origin`+`def` 를 든 구조체였고 `76a07b2` 가 앞의 둘을 캐시로 보내
/// 껍데기가 남았다. ⏳ **접는다 (2026-09-13 확정)** — `Surface` 와 같은 모양(진실 = enum 그 자체),
/// `m.vertex(v)` 가 바로 `match` 되게. 실측 2026-09-13: `VertexDef` 173자리 · `.def` 218자리,
/// 전부 기계적·census 무영향. 「문의 이름」 칸이 정점 자리를 어차피 전부 지나므로 거기 얹는다(열린 항목 21).
pub enum Surface {
    Plane {
        points: PlanePoints,
        motion: Option<Handle<MotionNode>>,
    },
    /// M6 ✔ **도착했다** — 예고가 아니다. 성분형이 옳은 이유·`ref_dir` 원시 규칙은 design.md
    /// 원통 절이 상세히 적는다.
    Cylinder {
        def: CylinderDef,                     // { origin, dir, ref_dir, radius } 전부 Rat
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
    /// ★ S5(ii)-1 ✔ — 생산자는 `DatumDef::ThroughVertices`. 능력의 핵심은 폭이 아니라
    /// **좌표로는 그 평면을 말할 수 없다**는 것(실측 220/220).
    /// ★ **핸들은 정렬해 저장한다** — 같은 셋 = 같은 진술. 이름 없는(무리수 모션) datum 은
    /// interning 이 못 받쳐 주므로, 순서가 다르다고 같은 평면이 두 핸들이 되는 것을 구성
    /// 시점 정규화로 막는다. 법선 방향은 점 순서가 아니라 정준 부호 규약과
    /// `Face::orientation`/`flip` 이 들므로 잃는 정보가 없다.
    Through([Handle<Vertex>; 3]),
}

/// 정점 = 자기 **정의**, 그리고 정의가 곧 타입이다 — `Surface` 와 같은 모양. 좌표는 진실이 아니라
/// 캐시다(PointCache). `Origin`·`point` 는 소멸. ⏳ 오늘 코드는 `struct Vertex { def: VertexDef }`
/// 껍데기(위 주석)이고, 접는 것은 열린 항목 21.
/// ★ 단일형이 아니라 **세 변종**이다(S7 에서 Q3 수정 — 아래 반증표): M5 에 원통 seam 정점이
/// 실재하고 그 점은 «세 평면»으로 적을 수 없다. 변종별 불변식(Q5)이 이 구조의 근거다.
pub enum Vertex {
    /// 세 평면의 교점 — 이름이 곧 점. (D != 0 은 좌표 재생이 생기는 자리에서 단언한다.)
    ThreePlane([Handle<Surface>; 3]),
    /// 두 곡면의 교차 «곡선» 위의 점 — M3 원통 seam(테두리 원의 θ=0). 점을 못 박는 매개
    /// 정보(원통의 `ref_dir`)는 **M6 의 원통 진실과 함께 왔다**(`CylinderDef`) — 그래서
    /// `OnSeam([cylinder, cap])` 은 「rim ∩ +`ref_dir` 광선」으로 **정확히 지정된 한 점**이다.
    /// ✔ 정의는 완성됐고 **재생하는 기계도 섰다**(칸 ㊵ `realize_vertex`, `seam_point`). ⏳ 없는 것은
    /// 그 값을 캐시에 **되쓰는 것**(칸 ㊸) — 그래서 STEP 은 오늘도 만들 때 나온 f64 를 내보낸다.
    OnSeam([Handle<Surface>; 2]),
    /// ✔ **`Pierce`** (2026-09-13, 열린 항목 24 완료): «branch» 는 «어느 근이냐»(수학의 가지)를
    /// 말하지 «무슨 점이냐»를 말하지 않아 개명 — 이 점은 도법기하의 **관통점**(piercing point)이다.
    /// ★ **M6 ✔ 도착 — 예고했던 것과 모양이 다르다.** 예고는 `Branch{surfaces, branch}` 였고
    /// 실제는 담체 종류를 **구조로** 말한다(닮은 핸들 셋이 아니라). `root` 는 정준 선 방향
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
    Rotate    { axis: Axis, pivot: [Rat; 3], angle: Angle },   // ⏳ 코드는 `point` — 열린 항목 21
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
    /// ★★ 유도 규약은 스펙으로 **동결**한다 — 바뀌면 기존 스케치가 조용히 돈다
    /// (수선의 발·Arbitrary Axis 는 이미 구현·테스트로 고정돼 있다: `c652ef4`·`20dae82`).
    /// ★ 구현(S4 ✔): 좁으면 유리수 `PlaneFrame`(비트 보존), 넘치면 판정층의
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

### 결정 요약 — Q1~Q5 의 답 (2026-08-05)

| | 답 |
|---|---|
| **Q1** 평면의 세 점 | **변종을 정점이 아니라 평면에 둔다.** 구성 평면 = `Known`(값 — 오늘 동작하는 `surface_points` 그대로), datum 평면만 `Through`(핸들). 핸들 쪽을 밀던 힘(발견 좌표는 가리켜야 한다)은 Through 가 받고, 값 쪽을 지키던 힘(좁고 종료 자명)은 Known 이 지킨다 |
| **Q2** `VertexDef::At` | **불필요.** «교점이 아닌 정점» 두 자리가 둘 다 사라지기 때문이다: **스케치 프레임의 base 정점**은 `Origin::Moved` 와 함께 소멸하는 이행기 산물 — 프레임 안 코너는 프레임 모션을 공유하는 세 평면의 교점이고, «공유 프레임에서 풀고 모션 재생 = 비트 동일(8/8)»이 이미 실측됐다. **프로파일의 공선 정점**은 표현할 것이 아니라 **만들지 말 것** — 프로파일 생성자의 공선 중간점 제거(무손실 정규화)로 닫는다 |
| **Q3** `Vertex` 최종형 | `{ def: VertexDef }` — **두 변종**(`ThreePlane` \| `OnSeam`). `Discovered{tol}` → `PointCache.tol`, `Moved` → 면의 모션, `Constructed` 태그는 정확성을 뜻한 적이 없어 잃는 정보가 없다. ★ 원안의 «단일형»은 S7 구현 중 반박됐다(아래 반증표) — 원통 seam 정점이 M5 에 실재하고 세 평면으로 적히지 않는다 |
| **Q4** `SurfaceDef` | `motion` 필드가 변종 안으로 — **«기록 없는 surface»가 표현 불가능**해져 극성 함정이 타입에서 소멸. 선행: `Store<Surface>` 봉인(§이행 S1) + `Inexact` 소멸은 임의정밀 이름(S2)이 선행 |
| **Q5** 원통(M6) | 버틴다 — 불변식은 variant 별(*"모든 **평면**이 점을 갖는다"*), 원통의 진실은 자기 variant 안. `Through` 는 정점의 정의 방식과 무관하게 성립. 알려진 절벽은 정점 쪽(이차곡면 셋 = Bézout 최대 8점 → 가지 번호 `branch: u8` 후보)이고 M6 에서 결정한다 |

**변종이 평면에 있는 이유** (모든 평면 점을 정점 핸들로 통일하고 `VertexDef::At` 을 두는 안의 기각 근거):
- `At` 의 «어느 프레임의 좌표인가»가 애초에 생기지 않는다 — `Known` 의 프레임은 그 평면의
  모션 전 프레임이고, 그 답은 `motion` 필드 바로 옆에 있다. 정점에 프레임을 달면 *"모션은
  정점이 든다"*(반증됨)로 돌아간다.
- append-only 에서 벽 평면이 자기 위상 정점을 가리킬 수 없다(정점이 평면보다 나중) — `At`
  정점을 미리 만들면 같은 점이 아레나에 둘(진술본+위상본), 위상이 재사용하면 «정점 = 세 면의
  교점» 통일이 깨진다.
- 정점은 술어·validate·배열이 가장 많이 소비하는 타입 — 단일형의 가치가 평면 단일형보다 크다.
  평면의 두 변종은 판정층 거울의 `match` 하나로 흡수된다. ⚠ **정정(2026-09-11)**: 그 거울은
  `WorkingPlaneDef` 라는 별도 enum 이 아니다 — 16-3 정정이 변종을 하나로 줄여 `WorkingPlane` 이
  증인 삼각형(`[WitnessPoint; 3]`)을 **직접** 든다.
- 저장 증가·재사용률 질문이 통째로 사라지고, `Known` 은 오늘 동작하는 코드 그대로다.

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
pub struct Ring2d {                           // 정점 + 조각. segs[i] = vertices[i] → vertices[i+1 mod n]
    vertices: Vec<[Rat; 2]>,                  //   온전한 원 = 정점 1 + Arc 1 (정점 = 솔기)
    edges: Vec<Edge2d>,                       //   ⏳ 오늘 코드는 `segs: Vec<Seg2d>` — 열린 항목 26
}
pub enum Edge2d {                             // 조각 «하나» — 순서 있는 고리라 시작=앞 꼭짓점
    Line,
    Arc { center: [Rat; 2], radius: Rat, ccw: bool },  // ⏳25 radius→r²·중심 정의화는 미정 마일스톤
}
// ⏳★★ **오늘 코드는 조각 타입이 둘**이다: 입력 `Edge2d{Line{from,to} | Arc{..start,end..}}`(양끝을 자기가
//    들어 순서 없이 흩어져 들어옴) + 저장 `Seg2d{Line | Arc{center,radius,ccw}}`. 열린 항목 26 이 kit 펜의
//    순서를 살려 둘을 합치고 `start`/`end` 를 없앤다(살아남는 이름 = `Edge2d`). 위 그림이 그 최종형이다.
// 링 더미 → 짝수 깊이 = 재료(even-odd) → 섬마다 Profile2d 하나 — from_rings, 현행 유지.

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
| **원통 원시체는 없다 — 원은 스케치다** (칸 ⑨) | `Operation::Cylinder` 는 원 프로파일 Extrude 와 위치 정준 비트 동일이 실측된 뒤 은퇴. 직선–호 꼭짓점 = `VertexDef::Branch`(두 평면의 **저장된 정준 이름** × 원통, 근은 매개 값으로 고른다 — 빌더가 자기 점으로 계수를 다시 만들면 `Lo/Hi` 가 뒤집힐 수 있다), 온전한 원 = `OnSeam` |
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
으로 **어느 문으로 들어왔는지가** 열쇠를 고른다. 위 모양은 **`push_raw` 한 자리**에서 진실을 보고
고르므로, 「이름이 있으면 이름으로, 없으면 그대로」가 **모든 곡면 종류에 자동으로** 적용된다.

★★ **`derive_name` 이 원통에 `None` 을 주는 것은 «없어서»가 아니라 «결정»이다.** 원통에 4계수형이
없는 것은 사실이지만, 그 함수가 지키는 것은 정책이다 — 누군가 원통의 정준형을 만들어 `Some` 을 주기
시작하면 기하가 같은 원통이 **조용히 합쳐지고 seam 이 갈라진다.** 그래서 그 팔은 `None` 을 **명시적으로**
돌려주고 이유를 그 자리에 적는다; 새 종류(구·원뿔)도 «합칠 것인가»를 먼저 정하고 나서야 이름을 얻는다.

★★ **`surface_ids` 는 캐시다** — 아레나를 돌며 진실마다 `surface_key()` 를 다시 매기면 통째로 재생된다.
그래서 `rebuild_surface_cache()` 의 형제로 «표도 버리고 재생 → 동일» 잠금이 서고, 그 잠금이
`SurfaceCache::Plane.name == 표의 열쇠` 일치까지 증명한다. ⚠ **프로덕션에서만** 그렇다: test-only
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
    //   ⏳ 오늘 코드는 surfaces·motions 만 비공개, vertices·edges·faces·shells·solids 는 `pub` —
    //   열린 항목 21 이 닫는다(면·셸의 push 문 신설 + 전량 순회 60곳 검증이 선행).
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
    //   (이 스케치는 `pub` 으로 그려 뒀었다. 코드가 옳다.)
    vertex_cache:  Vec<PointCache>,
    surface_cache: Vec<SurfaceCache>,         // ✔ 2026-09-11 — 실현. 아레나가 진실을 든다
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
pub struct Bounded   { value: f64,      error: f64 }              // ⏳ 오늘 `Approx`(cip, 50곳)
pub struct HpBounded { value: BigFloat, error: Mag }              // ⏳ 오늘 `HpApprox`(cip, 115곳)
//   ⏳ scalar 의 튜플 `Bounded = (BigFloat, Mag)` 는 `HpBounded` 가 흡수. 산술 20개도 scalar 로 이사.

// ── 계단 2: «경계 지어진 점» — 원자 셋(축별). realize 가 이미 `[Bounded;3]` 로 든다. ──
//   ★ 값과 오차를 «묶어» 든다(따로 든 배열 둘이 아니라). 그래야 «참값 ∈ value±error» 가 구조로 서고
//     transform 이 값만 옮기고 오차를 안 옮기는 desync 가 **불가능**해진다.
pub enum PointCache {                                             // ⏳ 오늘 `{coord: Point3, tol: Option<f64>}`
    Constructed([f64; 3]),                                        //   좌표만 — 측정 없음(checker 가 자기 ε)
    Discovered([Bounded; 3]),                                     //   값+오차 «묶음»(계단 2) — 발견 정점
}
//   ★★ **출처를 변종으로** — 오늘 `tol: Option<f64>` 의 `None`/`Some` 이 나르던 «구성/발견»을 타입이
//     구조로 말한다(reuse 가 그 유무로 «유리수 base 없음»을 읽는다 — `matches!(Discovered)`). 값과 오차는
//     발견 변종 안에서 [Bounded;3] 로 «묶여» desync 불가; 구성 변종엔 오차 자리가 아예 없다(있을 수 없는
//     것을 표현 안 함 — `SurfaceCache` 변종과 같은 원칙). 좌표 중복 없음(발견의 coord = Bounded.value).
//   ⏳ 오늘은 단일 struct 라 구성 정점도 `tol: None` 으로 «오차 자리»를 든다 — 칸 ㊸ 가 변종화하며
//     스칼라 tol 을 축별 `[Bounded;3]` 로 바꾼다.
//   ⏸ 고정밀판은 아직 «타입»이 아니다: 연산 하나짜리 실현이라 이름이 «Cache» 가 아니다(항목 23) —
//     `[HpBounded; 3]` 이거나 아래 `WitnessPoint` 자신이다(튜플 `(usize, [HpBounded;3])` 가 오늘 것).

// ── 계단 3: «경계 점 + 정의» — 판정 전용. 계단 2 의 상위집합이지 중복이 아니다. ──
//   ★★★★ `WitnessPoint = [Bounded;3](계단 2) + 정의(base:[Rat;3]·chain) + hp(메모)`. 그 «정의» 가
//     정확 단계의 입력이라 캐시(`PointCache`)와 «같은 것»이 아니다 — 이 차이가 진실/캐시 경계 그 자체다.
//     그래서 계단 1·2 는 합치되 이것은 캐시로 접지 않는다(접으면 경계가 지워진다).
pub struct WitnessPoint  { base: [Rat; 3], chain: /*모션*/_, coord_tol: [Bounded; 3], hp: /*메모*/_ }

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
//    「한 실체 = 한 캐시 구조체」가 되고, `rebuild_surface_cache()` **하나**가 실현과 이름을 같이
//    재생하며 그 **비트 동일 잠금이 「이름은 캐시다」를 증명**한다(간선이 만든 선례).
//    ☑ 실측(2026-09-12): `surface_name` 사용 44곳이 전부 «핸들로 조회»(`get` 24·`contains_key` 18·
//    `len`·`iter`·`insert` 각 1) ⇒ `HashMap` 이어야 할 이유가 없다.
// ★★ 접근자 이름은 §문의 이름이 정한다: 진실은 `surface(h)`, 캐시는 `surface_cache(h)`, 조각은
//    그 위에서 체이닝(`.plane()`·`.cylinder()`·`.name()`·`.tol()`). 다형 질의(`distance`·
//    `normal_at`)도 `SurfaceCache` 의 메서드다 — 오늘 geom 의 enum 이 하던 일 그대로.
// ☑ **대가 실측(2026-09-12)**: `m.surface(h)` 호출 **106**곳 중 **58 은 `match`/`matches!`**(변종을
//    바로 가르므로 오히려 나아진다) · **30 은 `let`-`else` 로 한 변종만**(형태만 바뀐다) ·
//    **손봐야 하는 건 18**(enum 통째로 11 + enum 메서드 7), 그중 대부분은 `distance`/`normal_at` 을
//    `SurfaceCache` 로 옮기면 사라진다. 통째로 변환하는 진짜 자리는 `transform.rs` **한 곳**.
// ⚠ 버린 스케치(2026-08-05, M6 전): `{ coeffs: [f64;4], tol: [f64;4], inv_norm }` — 평면 전용이라
//    원통으로 확장되지 않는다. `HpSurfaceCache{ coeffs, tol }` 도 같은 이유로 죽었다.
pub struct EdgeCache     { curve: Curve }                          // 평가 가능한 담체 곡선
// ★★★ **캐시는 «실현값 (+ 필요하면 그 오차)»** (2026-09-13 정정 — «셋 다 오차로 수렴»은 과했다):
//    `PointCache{coord, tol}` ✔(정점은 발견 잔차가 실물) · `EdgeCache{curve, err?}`(M7 후보) ·
//    `SurfaceCache{realized, tol?}` **동기 미확인**(위 §최종 타입 — 곡면 진실은 정확해 잔차가 없다).
//    ⇒ 오차 필드는 «소비자가 있으면» 붙지 대칭으로 붙지 않는다.
//    그 빠진 필드는 이미 이름이 있다 — design §3 의 `Intersection.cache_err`(*"캐시가 진짜 교차에서
//    벗어난 최대 거리"*). ⏸ 소비자 0 ⇒ 짓지 않고 적어 둔다(열린 항목 19).

// ★★★★ **세 «오차»를 섞지 말 것** — 이 절을 읽고 한 번 섞였으므로 갈라 적는다(2026-09-11).
//
// | 무엇 | 성질 | 상태 |
// |---|---|---|
// | `Plane::distance_eps(p)` | **이 거리 계산**의 f64 반올림, `3ε·Σ|pᵢ−oᵢ|` — **`p` 에 의존**(두 연산자) | ✔ 있음, **저장 불가** |
// | `SurfaceCache.tol`       | **이 평면 자신**이 참 평면에서 얼마나 떨어졌나 — 평면의 성질(한 연산자) | ✗ 없음 **+ 동기 미확인**(아래) |
// | `PointCache.tol`         | 같은 것의 정점판(*"measured residual"*)                                  | ✔ 있음 |
//
// ⚠ 「곡면당 저장하는 tol 은 틀린 양」은 **첫째에만** 맞다. 그 근거는 코드가 적어 뒀다 —
//   *"the tolerance it has to hand is the vertex's … a point that is exactly on the plane can
//   still produce a nonzero residual"*, 그리고 그 사고가 실재했다(*"residual **exactly equal** to
//   the claimed tolerance, saved only by the comparison being strict"*).
// ☑ 반면 **둘째는 `PointCache.tol` 의 곡면판**이므로 저장하는 게 옳고, **없는 것은 「대체됐다」가
//   아니라 「아직 안 지었다」**다. 만드는 법은 이 절이 이미 적었다(아래: 세 정점의 실현에서 유도).
// ☑ `inv_norm` 도 «없는 게 아니라 다른 배치»다: 오늘 `Plane` 은 **단위 법선**(구성 시 한 번 정규화)
//   + `raw`(정확 계수용)를 들어 같은 정보를 갖는다. 계수 우선 캐시를 고를 때만 필요한 필드다.
```

★★★★★ **f64 가 «둘»이고, 둘은 만나지 않는다** (2026-09-12 — 이 절이 한 사다리처럼 읽혀 생긴 오해를
갈라 적는다). 아래 표의 세 줄은 **한 사다리의 세 단이 아니라 세 «수명»**이다:

| | 모델 캐시 (`PointCache`·`SurfaceCache`) | 실현 (`WitnessPoint.coord`·`WorkingPlane.cache`) |
|---|---|---|
| 누가 만드나 | 구성·불리언 경로가 f64 로 계산해 **넘겨준다** | 판정이 **정의에서** 정밀도를 불러 만든다 |
| 오차 | `tol` 하나 — 저장된 점의 **잔차**(면에서 얼마나) | 축별 경계 셋 — «참값 ∈ 값 ± 경계» **증명됨** |
| 누가 읽나 | 테셀레이션·STEP·물성·validate | **판정 1단(f64 필터)** → 2단 정수 → 3단 상승 |
| 수명 | 모델과 함께 | **연산 하나** |
| 판정 경로에 | **없다** — cip 이 읽는 곳 0건(실측) | 그 자체 |

☑ **이것은 처음부터의 의도다** — 이 문서의 최초판(2026-08-01, `4c51a5e`)이 *"**캐시 이야기라 판정에는
영향이 없고**, 좌표가 결과인 STEP 출력은 `realize(p, 128) → round_to_f64` 로 정확 반올림해 낸다 ⇒ 남는
것은 표시·테셀레이션이 그 f64 로 충분한가뿐"* 이라 적었다. 구현은 그대로다. ⚠ **벗어난 것은 판정이
아니라 내보내기다**: 그 문장이 약속한 «정확 반올림 통로»는 칸 ㊵가 문(`realize_vertex`)만 세웠고 캐시에
**쓰는 것**은 아직이라, STEP 은 오늘도 날것 f64 캐시를 읽는다(기울면 최대 4 ulp) — 칸 ㊸이 닫는다.
⚠ **오해의 뿌리는 어휘다**: 원문 496행이 *"f64 필터가 가장 뜨겁게 읽는 `cache`"* 라 적은 그 `cache` 는
**`WorkingPlane.cache`**, 연산이 끝나면 사라지는 판정용 사본이다 — 모델 캐시가 아니다. 같은 단어가 두
수명을 가리킨다. ⇒ 판정 쪽 «cache» 어휘는 열린 항목 23(`Bounded`)이 «실현»으로 갈아 없앤다.

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
m.vertex_cache(v)                -> &PointCache                 // 변종(Constructed | Discovered) — §캐시
//   조각은 그 위에서 `match` — 발견이면 `[Bounded;3]`, 구성이면 좌표만. (오늘은 coord()/tol() 둘로 갈렸다.)
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

★★★★ **`_truth` 접미사는 사라진다.** `surface_truth` 는 이 커널에서 **유일한** `_truth` 였다
(`vertex_truth`·`edge_truth` 는 없다 — 정점·간선의 진실은 공개 필드 `.get()` 으로 읽혀 왔다).
접미사가 필요했던 것은 **이름을 캐시가 쥐고 있어서**였고, 칸 ㊷이 그 원인을 없앴다.
★ 오늘의 예외 둘도 규칙 안으로 들어온다: `plane_motion(h)` → `m.surface(h).motion()`,
`pub surface_name` 곁표 → `m.surface_cache(h).name()`.
☑ `m.motion(n)` 은 **이미** 이 규칙이다 — 새 규칙을 만드는 게 아니라 넓히는 것이다.

⚠ **`surface(h)` 는 개명이 아니라 «뜻 뒤집기»다** — 오늘은 캐시를 주고 **106곳**이 그것을 쓴다.
☑ 타입이 달라(`&geom::Surface` ↔ `&topo::Surface`) **전부 컴파일 오류로 드러나고 조용히 틀릴 자리가
0**이지만, 106곳을 한 건씩 「진실을 원했나 캐시를 원했나」 판정해야 한다(분포는 §캐시의
`SurfaceCache` 주석: `match` 58 · `let`-`else` 30 · 통째로 18).

☑ **나머지는 기계적이다** (실측 2026-09-12): `vertex_point` **111** · `edge_curve` **29** ·
`plane_motion` **18** · `vertex_tol` **17** · `surface_name` **44**. 저장소가 비공개가 되면
`.get()` **683**곳이 `m.vertex(v)` 꼴로 바뀐다(열린 항목 21).

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
- ✔ `Pt3` → `WitnessPoint` (+위성 `Pt3Error` → `WitnessPointError`)
- ✔ `HpIv` → `HpApprox`, 그리고 짝이 강제한 ✔ `Iv` → `Approx`(수치층 규칙 — 같은 것의 다른
  정밀도는 `Hp` 접두사 하나로만 다르다) + 필드 `mid`/`rad` → `value`/`error`, ✔ `Bound` → `Mag`
  (캐시 절의 철자 그대로 — `HpApprox { value: BigFloat, error: Mag }`)
- ✔ `PlaneGeom` → `WorkingPlane`(모양은 오늘 것 그대로 — 최종 모양은 열린 항목 16이 만든다)
- ★★★ **`FaceInfo` 는 따라가지 않았다 — 여기 적혀 있던 «둘 다 `WorkingPlane`» 은 정정한다.**
  면/평면 분리가 오늘 하중을 진다: `orient_sign`(이 면의) vs `frame_sign`(그 클래스의)은
  한 함수가 두 종류 인덱스로 불리던 시절 — *규약을 단언할 수 없던 유일한 자리* — 를 가른
  수선이었고, 한 이름으로 합치면 고친 결함의 이름이 되살아난다. 합치기는 개명이 아니라 설계
  작업이며 16 이후 재검토.
- `WorkingPoint` → `WorkingVertex`: 코드에 아직 없음 — 16이 만들 타입이 이 이름으로 태어난다.
- `tri_pt3` → `def`: 개명이다. ⚠ **정정(2026-09-11)**: 「타입 변경(`[WitnessPoint;3]` →
  `WorkingPlaneDef`)」이라 적어 뒀는데, 16-3 정정이 `Through` 변종을 없애 **타입은 그대로**
  `[WitnessPoint; 3]` 다(코드 실측: `WorkingPlaneDef` 없음).
- `Standard` → `ProofStandard`: 같은 계열(모양이 다르다 — `same_within` 유도로의 재구성) — 대응
  명시가 없어 유예. **기준: 이 목록에 명시된 것만 기계적 개명이다.**
(점 실현 묶음 `HpPointCache` 는 아직 타입으로 없음 — 오늘은 `(usize, [HpApprox; 3])` 튜플.)

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
// ⚠ **`WorkingPlaneDef` 는 지운다**(2026-09-11) — 바로 위 정정이 변종을 **하나로** 줄였으므로
//    별도 enum 이 존재할 이유가 없다. 코드에도 없다(실측). 판정 평면은 증인 삼각형을 직접 든다:
//        def: [WitnessPoint; 3]      // 유리수 base + 사슬. 오늘의 `tri_pt3` 그대로.

/// 증인 점 — 유리수 base 를 모션 사슬로 나른다. base + chain 이 정의, coord/tol 은
/// f64 캐시(= `PointCache` 모양), hp 는 고정밀 메모(= 평면의 `HpSurfaceCache` 와 대칭).
/// 같은 정의는 같은 실현(경로 무관)이다.
pub struct WitnessPoint {
    pub base: [Rat; 3],
    pub chain: Rc<[MoveNode]>,
    pub coord: [f64; 3],
    pub tol: [f64; 3],
    hp: Rc<OnceCell<(usize, HpPointCache)>>,
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
    homog: OnceCell<[Approx; 4]>,             // f64 층만 — 고정밀은 평면 메모에서 재계산
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
    반경 **6.3e-10**. 차수 9 가 `Approx` 의 여유를 먹지 않는다 — 재기 전에는 몰랐던 것이고,
    이것이 이 단계의 진짜 관문이었다.
  - ★★★ **깊이 2 는 필터가 없다.** 담체가 또 `Through` 면 차수가 **81** 이 되고, 계수가
    `f64` 범위를 **8/8 전부** 벗어난다(고정밀 쪽은 멀쩡하다). 함수는 그때 `None` 을 돌려
    **상승으로 보낸다** — `NaN` 반경이 우연히 `sign()=None` 을 내는 것에 기대지 않는다.
    ⇒ 깊이 2 는 «느린 길» 이 아니라 **상승 전용**이다. 깊이 제한을 두어야 하는지는 **열린
    항목 16 의 질문**이다(«2b» 는 이행표를 떠났다 — C6 이 유계를, C5 가 비순환을 이미 준다;
    남은 것은 비용뿐이다).
  - ★★ **배율의 부호는 값 안에서 없앤다.** join 은 행에 대해 다중선형이라 결과가 참 평면의
    `D0·D1·D2` 배이고, 음수면 **평면 방향이 뒤집힌다**(`frame_sign`·바깥 법선·라벨 프레임이
    전부 그 위에 있다). `Judge::plane_iv(k) -> [Approx;4]` 에 부호를 실을 자리가 없으므로 —
    실을 곳 없는 값은 아무도 안 쓰는 값이다 — 함수가 스스로 정규화한다. 생성 200 중 **172**가
    음수 `D` 를 지나므로 그 이빨은 실제로 물렸고, 정규화를 지우면 대조 테스트가 깨진다(확인).
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

### 완료된 단계 (근거·커밋은 git 이력과 dev-log)

| 단계 | 무엇 | |
|---|---|---|
| 0 | 측정 — 4018 정점 전부 인접 면 정확히 3, 후보 삼중 유일, 병든 삼중 없음 | ✔ 2026-08-01 |
| 1 | 계수 정준화 + interning | ✔ 2026-08-01 |
| 3a | `Vertex.definition` 을 `point` 옆에 추가 — 세 origin 100% 커버 실측 | ✔ 2026-08-01 |
| 2 | `Motion::Frame` — 기울어진 면 위 스케치가 그 평면 자신의 프레임에서 | ✔ 2026-08-03 |
| 1″ | 평면의 진실을 점 셋으로 (계수는 유도된 이름으로 강등) | ✔ 2026-08-04 |
| 3b | 좌표 재생 — **묻는 문은 섰고 덮어쓰기는 아직**(2026-09-11, 칸 ㊵). `nacre_ops::realize_vertex{,_decimal}` 이 정의에서 좌표를 실현하고 한 번만 반올림한다(변종 셋 전부: `ThreePlane`·`OnSeam`·`Branch`). ★ **그리고 재보고가 났다**: 축정렬 인구는 캐시가 이미 최근접 f64(48/48 비트 동일, 7비트 좌표), **기울어진 프레임 인구는 대조 가능한 12개 «전부» 최근접이 아니고 최대 4 ulp**. ⚠ 원인은 «맨티사를 넘음»이 **아니다**(좁은 인구에도 표현 불가 좌표가 16개 있고 전부 일치한다) — 그쪽 캐시가 더 긴 f64 유도에서 나와 아무것의 정확 반올림도 아니기 때문이다. 남은 절반은 캐시 덮어쓰기(`refine_caches`) | ⏸ |
| S1 | `Store<Surface>`·`Store<MotionNode>` 봉인 — 좁은 접근자(`surface`/`surface_count`/`motion`) + `compile_fail,E0616` 잠금. 기록 없는 surface 는 크레이트 밖에서 표현 불가(테스트 전용 입구 `push_surface_unrecorded` 만 예외, `test-util` 게이트) | ✔ 2026-08-05 |
| S2 | **임의정밀 이름** — `PlaneName{Narrow\|Wide}` 그릇(`plane_name_big` 꼬리의 `to_i128` 분기 하나), `surface_coeffs`→`surface_name` 개명, Wide 도 intern. **동일성 팔 절단**: 점 가진 평면의 이름 실패가 공선뿐이 됨 — ★ f64 폴백 연쇄의 **프레임 팔은 S4 몫**(Wide 는 `narrow()=None` 이라 프레임·지름길을 안 연다, 반증표 그대로). 잠금: 스칼라 2(`a_plane_too_wide…`·narrow/wide 합의) + topo 1(`a_wide_plane_interns_but_opens_no_shortcut`) | ✔ 2026-08-05 |
| S4 | **프레임 팔 절단** — `Motion::Frame{plane, placement, flip}` + `FramePlacement{Canonical\|Named}`. `Canonical`(기본)은 사슬 펼침 시 유도: 좁으면 오늘의 `plane_frame_default` 그대로(비트 보존 — 유도의 이사), **넘치면 wide 도로**(`MoveNode::FrameWide` = `PlaneFrame` 의 BigInt 쌍둥이, 실현은 `HpIv` 전 구간 + f64 캐시는 128비트 실현의 좁힘). Wide 이름·`n·n` 넘침(1.6%) 인구의 프레임이 열려 §12 연쇄의 프레임 팔이 닫힘. 잠금: cip 2(on-plane·tol-bounds wide 쌍둥이) + ops 3(`a_wide_plane_hosts_a_canonical_frame`·nn-넘침 개통·**종단** `a_pad_on_a_wall_with_overflowing_squares_takes_the_exact_road`) + census `wf` 가족. ★ 계획의 "심기 + 공개 SketchFrame 통일"은 **S9 로 분리**(아래) | ✔ 2026-08-05 |
| S3 | **`Profile2d` 유리수화 + 공선 중간점 정리** — `Profile2d{outer: Ring2d, holes}`·`Ring2d{points: Vec<[Rat;2]>}`, 생성자가 리프트(창 밖 = `ProfileOutsideDecimalWindow` 구성 시점 에러) + 공선 중간점 소멸(엄격 내부만 — 중복점·스파이크는 생존해 제 이름으로 보고). `check()`·`from_rings` 분류가 진실 위 정확 술어로(`orient2d_rat` scalar + geom `_rat` 쌍둥이 — geom 이 scalar 의존 획득, §design 1 격리 규칙 준수). 잠금: 쌍둥이 프리즘 비트 동일+전 코너 3-평면 정의(Q2 ② 닫힘), 창 에러 2, 십진-이진 부호 분기(십진이 이긴다), 비인접 공선 벽 interning 재핀. census 150줄 비트 동일. ★ 실측: check 는 34ms@100점(호출당 ~1.7µs, gcd 지배 — §열린 항목 7) | ✔ 2026-08-05 |
| S6a | **점 없는 평면의 소멸** — `Inexact` 소멸(S6b 타입 교체)의 전제. `PlaneDef` 를 점 셋 단일 필드로(origin=points[0]·ref_dir=points[1]−points[0]·극성=점 순서 — 불변식이 구조, 좁은-계수 def 실패 계급 사망, `named_plane_points` 은퇴), **`from_axes` 가 축의 십진 진실을 정의로**(#28 «축만 든 호출자» 인구 개통 — 45°급 프레임의 프리즘이 `WideFrame::named_of` 로 정확 경로; 내부의 실현-기저 호출 3곳은 의도적 무-def `realized_plane` 분리 — 두-정확-기술 재발 방지), `add_cylinder` 캡 점 기록(add_cuboid 선례), 넘침-이동은 노드 기록(`motion_is_exact`(칸 ④에서 `carry_of`로 — 법칙 행 S6c) 프로브에 점 수송 포함). 잠금: Named×Wide 프레임 단위 + 축-전용 기울어진 프리즘 종단 + 원통 캡 interning + 넘침-이동 + **전수 관문**(`points_coverage` — 생산 경로별 모델의 live 평면 face 전수가 점 보유). census 150줄 비트 동일 ×3회 | ✔ 2026-08-05 |
| S6b | **타입 교체** — 진실 스토어(`SurfaceTruth{Plane{points: PlanePoints::Known, motion} \| Cylinder{motion}}`)가 캐시 store 와 인덱스-평행으로 탄생, `SurfaceDef`·`surface_defs`·`surface_points`·`push_surface(_with_points/_unrecorded)`·`Violation::UndefinedSurface`·`RejectReason::{InexactSurface, CoordinateOutOfRange}` **사망**. push 는 `push_plane`(interning, flipped 는 f64 캐시 내적 그대로 — 같은 평면이라 부호 정확)·`push_cylinder(motion)`·test-util `push_plane_unregistered`/`set_plane_points_for_test`. f64 프리즘 폴백 → **이름 붙은 거절**(`PlaneWithoutExactForm`·`DistOutsideDecimalWindow` — `Swept::along` 삭제, build_prism 정확-전용). 부수 개선: 이동된 원통이 `Inexact` 강등 대신 모션 기록. ★ 구현 중 반박 1건: normal_def 의 `v = n×u` 곱이 작은-지수 전폭 법선(분모 10²¹→10⁴²)에서 넘침 — proptest 가 폴백 소멸 당일 발견, 원시 방향조차 137비트라 **기저-교차 셔플**(`w = x̂×n` + 대수 부호 `det[ẑ,x̂,n]=n₁`)로 재구성(곱 0개, 전역). 아레나 반전(캐시가 Store·진실이 Vec — `Handle<T>` 타입 매개변수가 강제)은 최종 개명 시 제자리로(§열린 항목) — ✔ **2026-09-11 `78d770f` 이 되돌렸고, 개명도 같은 커밋이다**(당시의 «강제» 는 갈라 커밋할 때만 참이었다: 한 커밋이면 핸들 철자가 안 바뀐다). census 150줄 비트 동일 ×4 | ✔ 2026-08-06 |
| S6c | **수송 법칙**(칸 ④, 2026-09-03) — `transform(rigid(R, t)) ≡ transform(T) ∘ transform(R)`: `carry_of`(옛 `motion_is_exact`)가 솔리드당 한 번, 후보 [전부·회전만·없음]에 같은 탐침 «실현(`Xform::point`) == 정확 상(`point_rat`/`mirror_point_rat`)»을 정점·평면 원점·**원통 축 원점**에 물어 첫 통과를 수송(`transport_points`/`transport_cylinder`)하고 나머지를 기록(`chain_motion(carry)` — 표 없음, parent 있으면 전부). 기록된 노드는 **기록된 부분** 앞의 진술을 말한다(회전 정확·이동 부정확 → 회전된 진술 + [T], `chain_translation` 접힘). 반증한 것: 회전 전 좌표 `p+t`로 재던 탐침(회전 뒤 반올림을 «정확»이라 함)과 «정확한 회전은 노드 없음 ∧ 부정확 이동은 노드»의 반쪽 체인(오프셋 보스가 릴리스에서 «안 닿는다»). | ✔ 2026-09-03 |
| S9 | **공개 스케치 API 통일 + world 평면 사전 심기** — ① `Model::new()` 가 세계 축 평면 셋을 심는다(핸들 0·1·2 = XY·YZ·ZX, points 는 `axis_plane` 삼중 `[0,u,v]`, **캐시 방향은 −축** — extrude 밑캡의 감각과 일치, +축이면 실측 781 캡 flip 재도입; `#[derive(Default)]` 제거 = 무씨앗 뒷문 폐쇄, `world_plane(Axis)` 접근자, `stat seeded_hits` 반증성 다리 신설 = 실측 455). census ε-재기준 1회: 평면 digest 이동 127/150줄, **결과는 143/150 비트 동일 + 나머지 7줄도 부피·면적·centroid 전부 비트 동일**(정점 해시만 이동 — Cramer 가 사실상 스케일-불변으로 반올림, 스칼라 최대 편차 정확히 0), ERR/EMPTY·피연산자 정점 해시 문자 동일. ② 공개 `SketchFrame{plane, placement, flip}`(필드 private + 검증 생성자 — 리터럴 우회 봉쇄): `named()` 가 구성 시점 거절 `FrameOutsideDecimalWindow`·`OriginNotOnPlane`(신규 scalar `plane_residual_sign` — orient2d_rat 급 **전역**, Wide 는 BigInt 팔)·`RefDirParallelToNormal`(판정은 `WideFrame::named_of` 재사용 — 폭에 전역이라 None = 평행뿐), 이름 없는 평면 = `PlaneWithoutExactForm` 재사용. `face_sketch_frame` 신설(이음새 — face_frame 이 만들던 값을 버리지 않고 공개). ③ 내부 통일: flip 측정은 `measured_frame` 한 곳, 노드 push 는 `push_frame_node` 한 곳(extrude·face 두 도로가 한 모양, 게이트 표현식 문자 유지, census 비트 동일). ★ **`Operation` 의 평면-핸들 어휘 교체는 S5 로 유예** — replay 자기완결성: 로그 속 핸들의 합법 표적은 씨앗·기존 면·datum 뿐인데 datum op 가 S5 에야 생긴다. ★ 잠금서 확정 둘: 씨앗 intern 직접 증거(원점 상자 바닥/왼쪽/앞 + z=0 밑캡 = 씨앗 핸들, 아레나 6 유지), ZX 의 canonical 프레임은 `−x̂`(스크립트 삼중과 다름 — Named 로 말할 사례임을 잠금이 명문화) | ✔ 2026-08-06 |
| S8 | **Edge 최종형** — `Edge{surfaces: [Handle<Surface>;2], vertices: [Handle<Vertex>;2]}`: 담체 두 면(오름차순 정렬 쌍) + 경계 두 점, `curve`·`bounds: Option`·`origin` 사망. `Store<Curve>` → `edge_cache: Vec<EdgeCache>`(인덱스-평행 캐시): 유일 입구 `push_edge`(eager 파생, 퇴화 검사는 **팔별** — rim `[v,v]` 는 합법) + `rebuild_edge_cache`(«버리고 재생» 잠금이 비트 동일 증명) + `edge_curve` 접근자·`derive_edge_curve`(직선 = 끝점 through_points, rim 원 = 담체에서 — 신설 geom `line_plane`, seam = 자기-인접 `[cyl,cyl]` 잠정 표기). transform pass 2(곡선 이동) 통째 소멸. validate: 신설 `EdgeCarrierMismatch`(담체 ≠ 인접 관측, `[plane,plane]` 자기쌍 검출) + `UnboundedEdgeInLoop`·`RefKind::EdgeCurve`·`StepError::UnboundedEdge` 순삭. ★ 구현 중 발견 2건: ① **담체는 wall 로 추측하면 틀린다** — 세 평면이 한 직선을 공유하는 인구(해결된 4-평면 동시성)에서 각 면의 arrangement 는 제3의 평면을 wall 로 (옳게) 지목 — 담체는 **전 링 선-주사한 인접성**에서 읽는다(실측: debug_assert 발화가 잡음). ② «전 생산 직선 비트 동일» 주장이 이동 경로에서 반박 — pass 2 는 방향을 직접 회전, 파생은 끝점 차 재정규화라 방향 ~1 ulp(실측 2.8e-16, 직선 83/84 비트 동일, 원 최대 2.2e-16 — 직선 기하는 비관측이라 무해). ③ VertexOffCurve 의 직선 갈래는 **타입상 항진**이 됐다(끝점이 자기 직선 위) — 검사는 원(rim)으로 이빨 유지, `.max(tol_of(edge.origin))` 은 상수 `EPS_CONSTRUCTED` 로 재철자(**무-행동이 아니었다** — `Discovered{tol:0}` 정점의 하한을 edge 항이 받치고 있었음, 실측). census 전 커밋 비트 동일 | ✔ 2026-08-06 |
| S7 | **`Origin` 소멸 — 정점은 자기 정의를 들고, 좌표는 캐시가 된다** — `Vertex{def: VertexDef{ThreePlane([3]) \| OnSeam([2])}}`(Q3 수정: seam 정점이 단일형을 반박 — 반증표), `point`·`Origin`(3변종) 사망, `vertex_cache: Vec<PointCache{coord, tol: Option<f64>}>` 인덱스-평행 + 유일 입구 `push_vertex` + 접근자 `vertex_point`/`vertex_tol`. **`rebuild_vertex_cache` 는 없다**(3b ⏸ — 발견 좌표는 배열이 공들인 값 1992 중 238 이 순진 Cramer 와 다르고, seam 좌표는 load-bearing. ⚠ **그 근거는 2026-09-11 에 지나갔다 — 열린 항목 12 의 정정을 볼 것**; 이 행은 2026-08-07 당시의 기록으로 남긴다): S8 이 모서리에서 얻은 «버리고 재생» 보증은 정점엔 아직 없음을 정직 기록. 소멸한 기계: 스케치 프레임 base 정점(Q2 — 프레임 공유 세 평면의 유리수 Cramer + 사슬 재생이 저장 좌표를 **비트 동일**로 재현, 8/8 실측을 영구 잠금으로 승격)·`remap_origin`·`solid_motion`+정점용 `move_node`(면이 자기 leaf 를 든다 — 규칙 3)·한-홉 base 불변식(타입이 흡수: 중복 적용이 표현 불가)·`exact.rs::base_f64/top_f64`. reuse `solid_points` 는 def 경로로(구성=`Pt3::exact` 문자 동일, 발견=포기 문자 동일, 이동=세 이름의 checked-i128 Cramer→replay; **혼합 프레임은 정직한 decline** = 기록된 유일한 차이). 게이트 `origins_are_remappable`→`defs_are_remappable` 전 정점 확장(발화 0 + **양성 대조**), validate: `tol_of` 1식화·`VertexOffDefinition` **전 정점 확장**(+양성 대조)·신설 `VertexDefCarrierMismatch`(변종 ⇔ 담체 종류). 신설 `nacre_scalar::three_planes_rat`. census **전 커밋 비트 동일**(좌표 verbatim 이사 — 재기준 없음) | ✔ 2026-08-07 |

| M6-0 | **원통의 진실** — `SurfaceTruth::Cylinder { def: CylinderDef, motion }`, `CylinderDef { origin, dir, ref_dir, radius }` 전부 유리수(dir·ref_dir 은 **비정규화 원시** — normalize 가 정확형을 파괴하는 `normal_def` 선례; 성분형은 «유도된 곱» 이 아니라 사용자 어휘의 리프트+셔플이라 반증표 무저촉). ref_dir 원시는 `any_perpendicular` **자신의 규칙을 유리수로**(최소-\|성분\| 축, 동률 X→Y→Z; ★ 축 선택은 **정규화된 `d`** — 캐시의 실제 입력 — 를 읽어 구조적 일치. 첫 철자의 «양수 배는 순서 보존» 은 실수-산술 논증이라 반증됨: f64 나눗셈 반올림이 강부등호를 동률로 붕괴 — 실측 축 `[0.34, 0.33999999999999997, 1]` 에서 원시=Y·`d`=X 로 seam ~90° 어긋남, 회귀 픽스처로 박제. 외적은 원시 정확 성분 그대로라 seam 방향 보존). interning 은 **보수적**(`cylinder_ids` 별도 맵, def 문자 동일 + motion — 다른 ref_dir 을 합치면 seam 이 갈라지므로 위험 0 키로 시작, 기하 동일성은 M6-1 술어 몫). **seam 은 모델 기하라 여기서 영구 고정**(+ref_dir, θ=0) — M6-1 의 유리수 반각 차트가 자기 배제점을 seam 에 맞춘다(역방향 절대 금지). OnSeam = "rim ∩ +ref_dir ray" 로 정의 완성(좌표 캐시 load-bearing 해제 선언 — 재생 기계는 3b 와 함께 유예), `[s,s]` 자기-인접 확정. 규약 전문 design.md §9. 관문: 사전-측정 비트 리터럴 4픽스처(축정렬·피타고라스 비트 동일, 기울면 ≤1ulp), census 원통 가족(+`CylinderFace` 거절 줄), validate `CylinderTruthCacheMismatch`(혼합 절대/상대 — ulp 계량은 전폭 축의 Gram–Schmidt 0-자리 스미어에 반증됨), step-io 독립 왕복 클린 | ✔ 2026-08-17 |

| M6-1 | **판정의 스칼라 원시 + branch 정점** — `nacre-scalar/quad.rs`: `QuadVal`(a+b√c, c≥0, `PartialEq` 비파생 — √8=2√2 정준화는 소인수분해라 값 동등은 sign 탑 몫), sign 탑은 BigInt 코어로 **total**(sign1: c=0 최우선, 상반은 `sign(a)·sign(a²−b²c)` 곱; `biquad_sign`: ℚ(√u,√v) 4항 — 원형 순서가 여기로 환원), `plane_plane_cylinder`(퇴화 사다리 전 변종 명명), **점 = MeetLine + s**(같은 라디칼 공유가 구조적), `circular_order_about_seam`(4계급, seam 모선은 순위 아닌 `SeamIncident` 이름; API 모양은 M6-2 후보). `VertexDef::Branch { planes, cylinder, root: QuadRoot }` — remap은 재정렬+**스왑 시 root 토글**(red 실측 잠금). ★ **정정(2026-08-21)**: 그 토글이 **무조건**이라 접점에서 틀렸다 — 중근은 스왑해도 같은 점(`disc=0` ⇒ `s=mid`, `mid↦−mid` 와 `ℓ↦−ℓ` 가 상쇄)인데 토글하면 한 점이 두 이름을 갖는다. 그리고 저장된 `Lo` 가 「둘 중 작은 쪽」인지 「유일한 쪽」인지 말하지 못해 정의를 든 쪽이 알 방법이 없었다 ⇒ `QuadRoot` 에 **`Double`** 을 더하고 규칙을 `QuadRoot::canonical`(인덱스 공간에 제네릭 — 핸들과 클래스 둘 다 답한다) 한 곳으로 모았다, 흡수 let-else 5곳 명시 match화, `carriers()` 한 철자. 관문: Pell 수렴자 자기-검증 사다리(~1e-73)·512-bit 차등 오라클·1e-18 분리에서 f64 초과 실증·census 비트 동일. ★ 발견: "sign1 하나로 닫힌다"는 과장 — 원형 순서는 sign2(설계 문서 정정) | ✔ 2026-08-18 |

| M6-2b 준비 | ★ **섞인 정점의 비교가 풀렸다**(2026-08-20) — 삼중항 점(`MeetPoint`, 유리수)과 분기점(`a+b√c`)의 좌표 비교가 `nacre_scalar::quad::cmp_coord_meet_branch`/`cmp_coord_branch`로 닫혔다. 둘 다 **BigInt로 들어** 부호 탑의 정수 코어(`sign1_int`·`biquad_sign_int`)에 넘기므로 **세 짝 모두 총체**다. declining을 만들 수 있는 유일한 길은 들어오는 길에 `Rat`로 좁히는 것. ★ 전제: **유리수 인구**(게이트가 `rotated`와 `base_rat` 부재를 둘 다 막는다) — 회전이 열리면 CIP 사다리의 일이다. 남은 것은 산술이 아니라 **이름**(배열은 아직 정점을 평면 삼중항으로만 부른다) | ✔ 2026-08-20 |
| M6-2b 준비 | **정점의 «방향»이 넓어질 자리** (2026-08-20) — 배열의 정점 처리는 좌표를 읽지 않는다(방향 = «어느 벽 평면 + 부호», 회전 = 평면 쌍의 부호). **호에는 탈 평면이 없으므로** M6-2b는 자료구조문제가 아니라 «방향» 개념이 넓어지는 문제다. 그 방향을 쓰는 두 자리가 같은 값을 **다른 철자**로 만들고있어(한쪽은 결과 반전, 한쪽은 인자 교환) `edge_dir`·`turn`으로 각각 한 군데에 모았다 — `turn`의 `0`은소비자마다 뜻이 달라(StraightAngle 거절 vs π 극 버킷) **그대로 돌려주고**, `edge_dir`의 `0`은 둘이같은 답을 요구하므로 원자가 가진다. ★ 넓어질 자리는 **셋**이다: ① 회전 원자(접선은 분기점에 **1차**라`sign1`로 닫힌다 — 접선끼리의 2차는 `CylinderPairContact`로 이미 거절) · ② `angular_order`의 π 극판정(«같은 벽, 반대 부호» → «같은 원, 반대 접선») · ③ `loop_winding`의 공선성 되감기. ④ 그리고 사전순최소 비교는 `nacre_predicates` 확장 산술 vs `QuadVal`로 **표현이 달라 아직 만나지 않는다**. ★★ 실측:회전 원자를 negate하면 83개가 빨개지지만 **단위 골든 둘은 초록**이다 — 읽을 규칙으로 픽스처를 만들어부호가 상쇄된다(규약은 엔드투엔드가 잡는다). `frame_sign = −1`에 닿는 단위 픽스처는 없다 | ✔ 2026-08-20 |
| M6-2a | **원통이 불리언에 들어온다 — ⊥ 절단 인구** — 게이트(평면 class × 원통 쌍마다 정확 판정: ⊥ 통과·∥벽은 `(n·o+d)² > r²|n|²`, **그것이 실패하면 그 클래스의 «면»들이 직접 답한다**(2026-08-19: 정리의 전제가 원래 «상대 몸통의 **경계**»였는데 게이트는 무한 평면을 재고 있었다 — 멀리 선 보스의 벽 평면이 구멍을 지나는 것만으로 거절됐다. 면 판정은 «**발자국 직사각형을 비켜가는가**» — ★ 2026-08-20: 벽이 원통에서 차지하는 자리는 그 평면 위 **직사각형**(가로=띠 `cylinder_strip_side`, 세로=그 lateral **면**의 span `point_axis_side`)이고, 직사각형은 두 띠의 교집합이라 **어느 한 축을 비켜가면 끝**이다 — 검사 둘이 아니라 **한 판정의 분리축 둘**. 면이 여럿이면 직사각형도 여럿이고 **전부**를 비켜가야 한다(간극이 통과되는 근거는 `bands_of`가 절단을 행의 span으로 잘라 **간극엔 밴드를 안 만든다**는 것 — per-span은 엔진이 세우는 슬랩과 같은 범위다). span은 **열린** 구간이라 캡 평면에 닿는 면은 통과하고, 그 뒤의 진짜 장애물은 `CircleMeetsSegment`가 이름을 준다. ★ 빈 span 집합은 «모두 비켰다»가 아니라 «축을 못 쓴다»로 읽는다 — 가드를 빼고 재니 진짜 가로지르는 벽과 접하는 벽이 둘 다 통과했다. ★ 두 축은 **축 정렬 직사각형 면**(=압출이 만든 벽, 챔퍼 포함)에만 완전하고, 앞선 불리언이 베어 문 면엔 보수적이다 — 완전판은 **면의 모서리 법선을 축에 더하는 것**, 같은 판정의 확장이다), 제곱근 없이 정확; outer loop만·직선 모서리만 — 볼록包 논증이 직선에 걸려 있고, 원판 캡 면은 꼭짓점 하나로 «전부 한쪽»을 만족시켜 버린다)·그 외 이름 붙은 거절. ★ **좌석 규칙은 없다**(2026-08-19 삭제): 좌석 원이 어려워지는 것은 그 **경계**가 상대의 경계를 만날 때뿐인데, 그 경계는 평면 위 모서리(→ ∥벽/기운 규칙)이거나 다른 원통의 rim(→ 쌍 규칙)이라 이미 전부 이름이 붙어 있었다 — 규칙은 자기가 말할 수 있는 것보다 넓었고, 실측하니 **울타리 둘까지 가려서**(모서리에 걸친 보스·나란한 두 원통이 `SeatedCylinderCap`으로 발화) 틀린 문장을 주고 있었다). 평면 배열의 **원 요소 2종**(transversal/seated)은 세그먼트가 아니라 **닫힌 셀**(★ 2026-08-20: «원이 세그먼트를 안 만난다»는 그 전제를 이제 `circles_meet_no_segment`가 **셀을 세우기 직전에 직접 확인**한다 — 옛날엔 벽 규칙이 곁다리로 지켜 줬고, 그래서 «원에 대한 약속»이 «벽에 대한 규칙»에 얹혀 있었다. 교차하면 `CircleMeetsSegment` — ★ 2026-08-20: 그 거절이 이제 **어디인지도 말한다**. 위치는 `plane_plane_cylinder(W, 세그먼트의 벽)`의 근이고 그건 곧 `VertexDef::Branch`가 이름 붙이는 `평면∩평면∩원통` 점, 즉 **다음 칸이 원을 호로 쪼갤 분할점**이다(그때 이 거절은 사라지고 위치기가 남는다 — 증인은 오늘의 소비자일 뿐). ★ 탐지기(`segment_meets_cylinder`)는 한 글자도 안 바꿨다: «실체 원통을 만나는가»와 «원을 어디서 가로지르는가»는 다른 질문이고, **원판 안에 통째로 든 세그먼트**가 전자엔 참·후자엔 교차 0이기 때문이다. 이름은 하나로 두고 **증인의 모양**으로 원인을 가른다(가로지름=`Point`, 원판 안 모서리=`Segment`, 길이 0=그 점 자체 — 붕괴 모서리에 위치기를 돌리면 엉뚱한 자리를 가리킨다). 도달 가능한 `CylinderMeet` 변종은 **셋뿐이고 증명된다**(원을 남기는 평면은 축에 ⊥ ⇒ 그 meet 선도 ⊥ ⇒ `OnRuling`·`AxisParallelMiss`·평행/일치 불가능). 증인 선택은 **좌표가 아니라 이름의 최소**다 — 좌표로 고르면 f64 반올림이 증인을 고른다)(유사 half-edge `2n+2i`, 라벨 XOR 무수정). 측면은 `bands.rs`의 **z-밴드**(★ 2026-08-19: 행은 **클래스가 아니라 면**마다다 — 한 lateral 표면이 한 솔리드 안에서 면을 둘 가질 수 있고(보어의 중간을 잘라내면), 옛 «클래스당 첫 면이 말한다»는 둘째 띠의 span을 잘라내 결과 껍질을 열어 놓았다 → `OpenResultShell`. span을 min..max로 합치는 것도 답이 아니다 — 면이 없는 구간에 없는 띠를 지어낸다): 균일-슬랩 정리 + **배열의 원판 셀 라벨**을 읽어 소속을 정한다(C4a의 유리 광선 도로는 K1이 삭제 — 밴드는 "경계 조각이 남나"이지 "점이 솔리드 안인가"가 아니었다). 조립은 rim 테이블 `(group, cyl, plane) → (OnSeam 정점, rim 엣지)` + 성분 결합 **둘째 규칙**(rim 키). 감김은 `sign(면의 바깥 법선 · 축)`에서 유도하고 **부피의 부호**가 그 잠금. ★ 벽의 **own 솔리드 소속도 라벨에서 읽는다** — 드릴은 자기 원통을 채우지만 물려받은 보어 벽은 그 안이 비어 있다(가정으로 박아 두면 두-구멍 판에서 틀린다). 관문: 관통/막힌 구멍(★ 2026-08-19까지 막힌 구멍은 **만드는 것만** 참이었다 — 그 몸통을 다음 불리언의 피연산자로 쓰면 `CylinderGateUndecided`였다. 원인: 원통 lateral 면이 ⊥ 클래스에 **rim으로 닿을 때** `Graze`를 안 냈고, `edge_mask`의 `Graze > Seated`는 정확히 «막힌 보어의 천장» 같은 reflex dihedral을 위한 것이다)·Fuse·Common·**N개 구멍**·드릴 판 flush 절단·**좌석 캡 5종**(양캡 flush 관통·바닥 flush 막힌 구멍·보스 Fuse·전체 Common·포켓 바닥 뚫기), 부피 해석값 일치·genus(1·2)·**OCCT**·STEP 왕복·validate 클린, census 평면 코퍼스 비트 동일. 남은 것: 회전 모델(M6-3)·기운 절단·∥벽이 측면을 가름(M6-2b)·원통쌍 접촉(M6b)·값 경로 i128 천장(`lateral_axis_span`의 t·chart 투영). ★ **곡면 성분의 깊이 판정은 닫혔다**(2026-08-21, `CurvedComponentDepth` 은퇴 — 발화 0곳): 광선이 담체 종류에 대해 일반화됐고(`cylinder_face_crossings` 가 `CylinderMeet` 일곱 변종을 **전부** 답한다 — `circle_crossings` 의 «원을 남기는 평면은 축에 ⊥» 증명이 탐침 광선에는 없다), **꼭짓점이 하나도 없는 성분**은 캡 원판의 **중심**을 증인으로 든다(축 파라미터에서 나오고 면 **위**에 있다 — 축의 내부 점은 다른 종류의 점이라 경계가 맞닿을 때 갈린다). 좌표 도로(`point_in_faces_rat`)는 되살아났지만 **곡면 팔은 이름 도로와 한 벌을 공유한다** — 두 도로가 graze 에 대해 다르게 답하게 될 자리가 없어야 한다. 남은 빚: 좌표 도로에서 곡면 팔에 `k ≠ 0` 을 싣는 픽스처(**기울어진 축**이 필요 — `cylinders_clear` 가 두 원통을 `r₁+r₂` 밖으로 강제하므로 나란한 축에서는 0회 아니면 2회다) | ✔ 2026-08-18 |
| M6-2b | **분기점이 이름을 얻는다** (2026-08-21, `cce9e2e`·`b8de145`·`3c1fa89`) — 「남은 것은 산술이 아니라 이름」의 그 이름. `NodeId::Branch{planes[2](오름차순 **클래스**), cyl, root}` + 유일 생성자 `NodeId::branch`, 문 `three_plane_name -> Option`, 그리고 오늘 분기점을 계산하는 유일한 자리(`circle_crossings`)가 위치용 `u8` 대신 **그 이름을 댄다**. 증인 좌표는 이름에서 다시 푼다(`combinatorics::branch_point`) — 진실은 정의, 좌표는 캐시. ★★★ **`QuadRoot` 가 `Lo|Hi` 둘이 아니라 `Lo|Hi|Double` 셋이 됐다**: M6-1 이 접점을 `Lo` 로 적기로 한 규약이 «둘 중 작은 쪽»과 «유일한 쪽»을 구별 불가능하게 만들었고, 그래서 `transform` 의 **무조건 토글**이 접점에서 틀린 채 아무도 못 고칠 상태였다(스왑해도 같은 점인데 토글하면 한 점이 두 이름). 재정렬 규칙은 `QuadRoot::canonical<T: Ord>` 한 곳(핸들 인덱스와 클래스 인덱스 **둘 다** 답한다)이고 접점 예외는 `flipped(Double)=Double` 에서 저절로 나온다. ★★ **문 뒤가 두 종류**임을 확정: 탐침 목록은 원소를 잃어도 되고(`three_plane_probes` 가 그 licence 를 이름으로 나른다) 링은 안 된다(`Option<Vec<_>>` 로 모아 문법이 막는다). 새 어휘 둘 — `DeclineKind::BranchNode`(트레이서의 기존 와일드카드가 `TraceDeclined{kind,face}` 로 승격하므로 **면 핸들까지** 붙는다)과 `RejectReason::BranchVertexUnnamed`(트레이서 밖 네 자리). `self_touch_reject` 는 **기권**한다(거절 가드가 사퇴하면 «못 봤다»가 «잘못됐다»가 된다 — 자기 doc 의 정책 그대로). `reuse::canonical` 의 키는 `CanonNode` 로 넓어지며 `usize::MAX` 센티널 둘이 사라졌다. ★★★ **실측 넷**: 게이트가 보는 2060쌍 중 463(22.5%)이 내림차순인데 **로케이터에 도달하는 9쌍은 전부 오름차순**이고, 그 이유는 면 push 순서(`add_cuboid` `[−Z,+Z,−Y,+Y,−X,+X]` · `build_prism` 캡→벽)이며, 코퍼스의 원통 축은 **전부 `+Z`** 였다. ⇒ 정준화는 프로덕션에서 한 번도 안 돌았다. 그래서 **축이 `+X` 인 보스**(이 저장소 첫 비-Z축 원통) 픽스처 둘을 지었고, 각각 근이 **정확히 하나**여서 `||` 없는 단언이 된다(둘 다 통과하면 이름 오류와 실현 오류가 증인에서 **상쇄**된다 — 유도 후 그것을 피해 배치). red 실측: 토글을 빼면 교차 픽스처만, `Double` 을 `Lo` 로 흘리면 접점 픽스처만 빨개지고 **기존 울타리 넷은 둘 다에 눈이 멀어 있다**. 관문: 울타리 넷 + census 픽스처 **비트 동일**(도달 9쌍이 오름차순이라 정준화가 항등), 비트 census 두 프로파일 무변화, reject census 는 스왑 인구 행 하나 증가(13→14). **남은 것**: 호 분할 자체, 모서리 «방향»의 arity(호는 두 끝의 접선이 다르다), 그리고 씸 표가 `VertexDef::Branch` 를 주조할 때 **클래스 순서 ≠ 핸들 순서**라 `QuadRoot::canonical` 을 두 번째로 물려야 한다는 것 | ✔ 2026-08-21 |
| M6-2b | **호가 셀을 둘린다** (2026-08-21, `1612335`·`3b54717`·이 칸) — 「넓어질 자리 셋」이 전부 답을 얻었다. ① 방향은 `EdgeDir::{Line, Arc}` 이고 **`(모서리, 노드)` 에서만** 태어난다(`dir_at`, 유일 생성자, `node ∈ {e.node, e.to}` 를 `debug_assert`); ② 회전 원자는 새 원시 없이 닫힌다 — 원 경계의 평면이 축에 ⊥ 이므로 `n_P ∥ m` 이고 BAC-CAB 이 `(d × (m × r))·m = (d·r)(m·m)` 로 무너져 **`sign(d·(x−c))` 한 번**(`quad::plane_side`); ③ lex 최소는 `CoordKey::{Three, Branch}` 로 갈려 `cmp_coord_meet_branch`/`cmp_coord_branch` 를 부른다(cip 무변경). 담체도 사수 `RingEdge.carrier: Carrier::{Plane{wall, sense}, Arc(cyl, def, ccw)}` 로 넓어지고, walk 가 **세 범위**(세그먼트 · 호 · 안 잘린 원의 유사 반모서리)를 돈다. ★★★★ **정준 이름은 저장 법선과 절반의 클래스에서 반대다** — `class_coeffs_rat` 로 `d` 를 만들면 정확히 거기서 거꾸로 읽는다(`plus_t_is_above` 의 「36개가 한꺼번에 빨개졌다」와 같은 함정). `stored_coeffs_rat` 이 그 회전을 한 곳에 모으고 `NameInts::coeff_sign`(cip 에 연 문 하나)과 debug 교차검사한다; **실측 발화**(도달 세 클래스 중 하나). ★★★ **walk-back 이 들고 나오는 값의 두 철자 중 하나만 옳았다**: 직선에서 같고 지름 현에서 정확히 반대 — `break earlier`(공유 노드에서 읽은 것)이고, 굽은 구간을 걸어야 하면 `CurvedStraightRun` 으로 거절한다(미발화). ★★ **red 프로브가 전역 부호를 못 본다**: 같은 원자가 각 순서와 감김에 둘 다 들어가고 walk 가 양쪽 handedness 를 시도하므로 흡수된다 — 잠긴 것은 `ccw` 와 정준→저장 회전뿐이고 `axis_up`·`frame_sign` 은 미행사(둘 다 코퍼스에서 상수). 실측: 두 클래스가 각 **셀 4개**를 첫 시도에 내고(`-1` 개수 = 성분 1), walk-back 걸음 **0**, 코너 없는 링 **0**. **남은 것**: 조립 — `circle_of`·`mask_of`·`bound_of`·원 셀 붙이기가 전부 `2·segs.len()` 로 종류를 가르므로 정지판(`ArcBoundNotYet`)이 `extract_cells` **반환 직전**에 선다 | ✔ 2026-08-21 |
| M6-2b | **정지판이 조립으로 한 걸음** (2026-08-22, `7d6919e`·`1157861`·이 칸) — 반모서리의 «종류» 가 타입이 됐다. `ClassEdges{segs,arcs,circles: Cow}` + `HalfEdgeKind::{Seg,Arc,Circle}`; 다섯 자리가 `he >= 2·segs.len()` 로 두 갈래만 묻던 것이 `kind(he)` 한 곳이 되고, **`edge_at` 이 `RingEdge` 의 유일한 생성자**가 된다. ★★★★ **그 다섯은 「틀린다」가 아니라 «패닉» 했다** — 정지판을 한 단계 밀자 `nest_cells::ring_of` 가 `segs[he/2]` 로 범위 밖 인덱스를 냈다(계획의 목록에 없던 다섯째이자, 바로 그 자리). ★★★ **`RingEdge` 가 두 번 쓰여 있었고 두 철자가 달랐다**: walk 는 분할의 `sense` 를 넘기고 `ring_of` 는 `sense: None` — 오늘 답이 같은 이유는 구조적(분할 산출물이 정지판 밖으로 안 나간다)이고, 이 칸이 그 구조를 바꾼다. `Cow` 로 소유 가능하게 해 전처리의 **두 벌**을 잠금이 아니라 **구조로** 없앴다. ★★ **커밋을 가른 자리가 두 번 움직였다**: 「타입 vs 배관」은 성립 안 함(`of` 가 분할을 도므로 타입이 배관을 요구한다) ⇒ 「배관 전부 vs 정지판 한 줄」 ⇒ A 가 커지자 **위험 hunk**(`ring_of` 접기)에서 다시 갈라 A1/A2. 셋 다 census 로 검산. ★★★★ **잠금 둘**(초안의 「잠금 없다」가 틀렸다): ① 감사(`frame_audit`, 자기 파이프라인을 인라인으로 다시 돈다)와 불리언이 같은 `failed_at` 을 낸다 — `&[]` 로 «컴파일만 되게» 고치면 감사가 끝까지 성공해 `None` 이 되고, 기존 합의 테스트의 픽스처는 일부러 클래스 밖에서 실패하는 입력이라 이 자리를 못 본다; ② `ClassAudit.produced = {cells 4, arcs 2, roots 1, holes 0}` — **전부 유도**(rim 이 모서리를 두 번 가름 ⇒ 호 2·셀 4; `−1` 하나 ⇒ root 1; 바깥이 세 `+1` 과 노드 공유 ⇒ hole 0; 보스가 윗면 위에 서므로 잘리는 클래스 하나). ★★ **못 보는 것을 「미정」이 아니라 유도로 적고 확인**: 분할의 `sense` 뒤집기는 호를 다른 셀에 붙이면서도 네 숫자를 안 건드린다 ⇒ 잠금 통과 (실측 초록). 첫 독자는 `label_cells` 의 keep 결정. **남은 것**: `label_cells`→`emit_faces` 두 걸음, 그다음 위상 조립 — 첫 항목은 `edge_for` 의 키(A·B 사이에 현·안쪽 호·바깥쪽 호 **셋**이 순서 없는 정점 쌍으로 접힌다; `Curve::Circle` 은 이미 있으니 **직선은 순서 없이, 원은 순서 있게**), 그리고 목록에 없던 여섯째 **밴드 패스**(잘린 원은 디스크 라벨이 아예 없다) | ✔ 2026-08-22 |
| M6-2b | **정지판이 라벨을 지난다** (2026-08-22, `b305c22`·`ae55818`; 곁가지 `63100aa`) — 단계 열(`walk→nest→label`)이 `per_class` **한 벌**이 되고 정지판이 그 뒤에 선다. ★★★★ **모양이 중복 제거만큼 중요하다**: 호출부는 정지판을 돌린 뒤에 풀어야 하는데, 단계를 따로 든 `Result` 로 쥐면 그 순서를 **단계마다** 되풀이해야 하고 틀리면 **조용하다**(스위트 초록, 거절의 **이름만** 바뀜 — A2 가 낸 것). `Result` 하나면 `?` 자리가 하나뿐이고 그것이 정지판 뒤다. ★★ 부수 효과 둘: 정지판 이동이 **저절로** 되고(단계를 넣는 순간 그 뒤가 된다), `face_of` 가 `Staged` 에서 **빠진다**(유일한 소비자가 안으로 들어갔다 — 단계가 늘고 인터페이스가 주는 것이 경계를 옳게 그었다는 표시). ★★★★★ **지난 계획의 예상이 반증됐다**: 「`sense` 프로브의 첫 독자는 `label_cells`」— 아니다. 라벨 *순서* 는 움직이지만 순서에 무관한 요약은 **전부 동일**(링 크기·호 인덱스·라벨 삼중항까지 대조). 두 3-셀이 라벨을 맞바꾸는 대칭이라 어떤 다중집합도 못 본다. **첫 독자는 `emit_faces`** — 링이 **역순**으로 나온다 ⇒ `sense` 잠금은 `Ring` 담체 칸의 몫. 잠금: `Produced.pos_labels`(정렬된 `+1` 라벨 셋, **유도값**) + red 둘 (호가 비트를 안 뒤집으면 빨개짐; 라벨 단계를 실패시켜도 이름 **유지** = 가로채기 살아 있음). ★★★ **계획을 열한 번 검토해 잘라낸 것**: 동어반복 단언(`root_label` — `seed` 를 되받음) · 함축된 단언(`kept` — 라벨이 결정) · 근거 잃은 추출 (`is_result_face` — 소비자 하나) · **컴파일도 안 되는** red 프로브 · **정지판을 함수 안으로 + `note` 콜백**(8·9회차가 짓고 10회차가 되돌림) · 자기 모순 셋(요약·위험·파일목록이 본문과 반대). ★ 10회차의 원칙: **「철자 불가능」이 언제나 옳은 것은 아니다 — 그 자리를 이미 보는 가드가 있으면 기계는 값이 아니라 비용**(정지판을 건너뛰는 실패는 `bound_of` 의 assert 가 이미 본다; 안 잡히는 것은 «에러가 먼저 새는» 쪽뿐). ★ 실측 추가: `emit_faces` 도 이미 돌고 면 둘을 낸다(평판 윗면 = 세그먼트 5 + 호 1, 걸터앉은 밑면 = 현 + 호 `flip=true`); **잘린 원의 디스크 라벨 = 0**(목록에 없던 여섯째). **남은 것**: `emit_faces` 한 걸음 — `Ring.walls` 를 담체로(소비자 여덟), 그다음 `edge_for` 의 키와 밴드 패스 | ✔ 2026-08-22 |
| M6-2b | **정지판이 면을 지난다 — `sense` 가 잠긴다** (2026-08-22, `926f31e`·`183a7c9`·이 칸) — `emit_faces` 가 `per_class` 에 들어가 정지판이 저절로 그 뒤가 되고, `Staged` 가 `{cells,nesting,labels,faces,disk_labels}` 로 **넓어진다**. ★★★★ **「단계가 들어가면 인터페이스가 준다」는 법칙이 아니었다** — 지난 칸은 그 필드(`face_of`)의 소비자가 하나였을 뿐이고, 기준은 「줄었나」가 아니라 **「각 필드에 소비자가 있나」**다(앞 셋은 `frame_audit` 이 읽으므로 `#[cfg(test)]` 로 어느 소비자인지까지 타입이 말한다). ★★★★★ **조사가 지난 칸의 「다음」을 둘 다 반증했다**: 첫 벽은 `Ring` 이 아니라 **씸 표**(`BranchVertexUnnamed`)고 소비자는 여덟이 아니라 **여섯**(`:577·615` 는 `NamedRing` 이었다) ⇒ **담체 확장 없이** 한 걸음 더 갔다. ★★★★ **가드는 문장이 아니라 «자리» 를 옮겼다**: 「호가 여기 닿지 않는다」가 거짓이 됐으므로 참인 명제 「sentinel 이 **읽히지** 않는다」를 **읽는 쪽**(`boolean.rs` 의 `planes[r.walls[t]]`)에 세웠다 — red 로 발화 확인. ★★★★★ **잠금은 `ClassAudit::outer_rings`** — 면의 바깥 링을 **좌표 열**로 들고 **회전만 정규화, 역순은 안 접는다**(노드 신원은 클래스 인덱스라 유도 불가 ⇒ 순환; `Produced` 안에는 못 넣는다 — `[f64;3]` 은 `Eq` 가 아니고 무리수 `2−√3/4` 가 `near()` 를 요구한다). red: 분할의 `sense` 뒤집기가 **정확한 역순**을 내며 빨개진다 — **돌린 보스에서만**. ★★★★★ **회전 정규화는 `n ≤ 2` 에서 역순에 눈이 멀다**(역순이 곧 회전) ⇒ 걸터앉은 보스의 2-링은 어떤 프로브로도 무방비이고, 이것이 `sense`-눈멀음과 **독립된** 둘째 이유다. ★★★ **유도의 «전제» 도 단언한다**: 링이 CCW 로 도는 기준 `n_out(wc)` 은 **어느 면이 클래스 루트냐**에 달린 평면 표의 사실이지 픽스처의 숫자가 아니므로, 링 단언 앞에 그것을 박았다(계측: 둘 다 `orient_sign=+1`, `+z`/`+x`). ★★ **계획 8회차가 3-링의 방향이 «역순으로 적혀 있던 것»을 잡았다** — 7회차의 「방향은 유도한다」가 5-링에만 적용돼 있었다; 내부-왼쪽 규칙으로 검산했고 5-링이 음성 대조군이다. ★★★★ **자체 점검 둘이 「고친 것」이 아니라 「고친 것의 설명」에서 나왔다**: ① `three_plane_name` 의 doc 이 *"must be empty"* 라고 적어 둔 게이트를 좌표 도로가 깼다(스위트 전부 초록인 채로) ⇒ `combinatorics::node_point_f64` 로 옮기고, 그 doc 이 **지금 비어 있지 않다**는 것과 남은 셋의 이유(문이 한쪽뿐 — `branch_name` 쌍둥이가 없다)를 적었다; ② **울타리 자신의 doc 이 정반대를 말하고 있었다** — *"this fence is green with that flip"* 인데 이제 빨개진다. ★★ `arc_side` 의 red 표는 **다시 재서** 여전히 참임을 확인했다(전역 부호 뒤집기, ops 582 초록) — 기록된 계측은 새 잠금이 생길 때 다시 재야 한다. **남은 것**: 씸 표의 `VertexDef::Branch` 주조(첫 벽) → `Ring.walls` 담체 + `edge_for` 의 키 → `has_curved_bound`·`merge_component`(같은 키 충돌) → 밴드 패스 | ✔ 2026-08-22 |
| M6-2b | **씸이 분기 정점을 realize 하고 정지판이 배열을 떠난다** (2026-08-22, `3179337`) — 씸 표의 분기 팔(좌표 `branch_point` 재해석 · tol `branch_vertex_tol` = 평면 둘 + 원통면 + **meet line** — 축-평행 평면 ∩ 원통은 자오선 둘이라 닫힌 형식 쌍별은 그 하나뿐, 나머지는 기록하고 뺌), 그리고 정지판이 클래스에서 «만들어지고» 씸 스트레치(unify→밴드→씸) 뒤에서 «터진다». ★★★★★ **가로채기가 두 층에서 같은 모양** — 단계들이 돌고 deferred 가 그들의 답 위에서 이긴다(`Ok` 든 `Err` 든); 조사가 잰 스트레치의 실패 지점 다섯(`unify`·`cyl_rows`·`band_faces`·씸 채움·`SeamAlias` 스캔) 전부가 덮인다. 이름·증인 무변화(기존 울타리 둘 + reject census 두 행 그대로). ★★★★★ **구현이 계획의 잠금 하나를 죽은 것으로 판정**: 「팔이 없으면 `BranchVertexUnnamed` 가 먼저 나가 빨개진다」는 공짜 잠금이 3회차의 가로채기에 **삼켜진다** — 팔을 통째로 꺼도 불리언-층 울타리 전부 초록(실측). ⇒ 씸 채움을 **`seam_table`** 로 추출하고 새 울타리 `the_seam_realizes_a_branch_vertex_and_measures_it` 가 프로덕션이 먹는 바로 그 면 목록으로 직접 부른다(deferred 이름 · 분기 엔트리 둘 · 각 `tol < 1e-12`; 좌표는 `outer_rings` 가 같은 road 로 이미 잠금). ★★★★ **문의 나머지 반쪽이 생겼다**: 초안이 `NodeId::Branch` 를 모듈 밖에서 철자해 문 게이트를 위반 ⇒ doc 이 예고한 쌍둥이 `branch_name` 을 만들고 **기존 위반 셋도 이행**, 게이트 0 히트. red 넷: tol 원통 항 제거(그 단언만 빨강) · 팔 끄기(새 울타리만 빨강 — 분업) · 클래스 층 실패 · 씸 층 실패(둘 다 이름 유지). **남은 것**: 정점 주조 — `edge_faces` 가 밴드 면을 못 봐(`poly_rings()` 만) `def_triple` 이 비고 `StraightAngle` (거짓 문장)로 나온다; `edge_for` 키 충돌 · `Bound::Band` 의 `[v,v]` rim 과 같이 본다 | ✔ 2026-08-22 |
| M6-2b | **모든 결과 정점이 이름을 얻는다 — 주조 직전까지** (2026-08-22, `2b5d952`·`e32f6a5`) — `reconstruct` 의 이름 절반이 `name_result_vertices` 로 추출되고(①, census 비트 동일), def 가 `Def{Three,Branch}` 로 넓어진다(②): **분기 def 는 선언**(이름의 payload — 유도 불필요), **물린 코너는 carried-walls fallback**(far_plane 완전 굶음 한정 + 모든 유도 면의 일치 요구 — 갈리면 결손). 정지판 raise 는 `reconstruct` 안 「이름 뒤·rim 주조 앞」 — 앞이 전부 **시그니처로** 모델 무변이(`group_faces`·`self_touch_reject` 가 `&mut Model` 을 안 받음). 가로채기 셋째 층, 불리언 층은 유지. ★★★★ **reject census 가 계획 밖 사실을 잡았다**: grouping 의 성분 기계(`ring_edges_with_walls`)가 호 링의 분기 노드에서 `branch_vertex_unnamed` 를 울리고 held→deferred 로 삼켜진다 — 이름 불변, census 는 「삼켜진 가드도 울린 가드」. 두 행 blessed. ⇒ **self_touch 는 호 입력에서 건너뛰어진다**(grouping held-Err) — 성분 기계의 분기 팔은 모서리 칸 이후. red 셋: 분기 팔 끄기(완전성 빨강 — fallback 이 MAX-wall 가드로 못 구함이 곧 증명), fallback 끄기(돌린 보스만 빨강), 이름 실패(이름 유지). 새 울타리 `every_result_vertex_of_the_arc_population_is_named`(결손 0 · 선언 둘 · 코너 def 가 `(4,0,2)` 로 realize). **남은 것**: 모서리 칸 — T-정션(쪼개진 쌍둥이: 이웃 클래스가 분기점을 알 수 없음이 실측 ⇒ 주조에서 가른다) · 호 키 접힘(현+호+호 한 키) · `Ring.walls` 담체 · rim `[v,v]`; 그 뒤 `VertexDef::Branch` 주조 + 정준화 둘째 답 | ✔ 2026-08-22 |
| M6-2b | **쪼개진 쌍둥이가 갈라지고 분기 정점이 주조된다** (2026-08-23, `c4b34f0`·`f729fb9`) — ① T-정션 분할: 분기 노드의 평면쌍 == 모서리의 {own,wall} 인 이웃 모서리를 그 노드에서 가른다(이웃 클래스는 그 점을 알 수 없음이 실측 — boolean-side 가 유일한 자리). betweenness = 두 끝이 다른 축의 좌표 betweenness(선 위 단조), 비교자는 `cmp_key` 어휘 그대로(`branch_between`) — **새 정확 기계 0**. 잠금 = edge-key census(양-끝-분기 아닌 세그먼트 키 == 정확히 두 면; 예외에 「분할 자신이 만드는 가운데 조각」 포함). fallback 은 무인구 그물로 은퇴(doc 재작성, 삭제 안 함). ② `node_handle` 분기 팔이 주조가 됐다: 클래스→핸들 + **`QuadRoot::canonical` 둘째 답**(빚 청산; ★ mint-flip 은 오늘 코퍼스 미행사 — canonical off 도 초록, 실측·기록). 정지판 raise 가 정점 실체화 루프 뒤로 (4층 가로채기). ★★ **의도된 거래**: 호 거절이 주조 정점을 store 에 남긴다(live-set 밖 garbage — 아레나 잔여와 같은 부류, replay 규율이 덮음; raise 주석·reject doc 에 명시). 그 거래가 곧 관측: 새 울타리가 거절 후 store 전체 필터로 분기 정점 둘(핸들 오름차순·tol 상한·유도 교차점)을 단언. red 넷: 분할 off(키 복귀) · 팔 백스톱(주조 울타리만 빨강) · canonical off(초록=미행사) · 루프 실패(이름 유지). **남은 것**: 담체 칸 — `Ring.walls`→담체(소비자 여섯) · `edge_for` 키(직선 무순서·원 순서) · sentinel 은퇴; 그 뒤 rim 잘린 판·Band 팔·grouping 분기 팔·밴드 패스·위상 봉합 | ✔ 2026-08-23 |
| M6-2b | **링이 담체를 들고, 호가 모서리로 주조된다** (2026-08-23, `a14d834`·이 칸) — ① `Wall{Plane,Arc{cyl,ccw}}` 가 `usize::MAX` sentinel 을 은퇴시킨다(`bound_of` 가 `edge_at` 의 규약 «짝수 반모서리 = CCW» 를 **인용해 나름**, 새 유도 0; `Ring::edges` 는 legacy shim — names-road 담체화는 grouping-팔 칸; census 비트 동일). 공짜 수확: dissolve 의 `MAX==MAX` 위험 **구조 사멸**(남은 동등쌍 = 같은 원·같은 방향 — 그 경우는 동등이 곧 정답, doc 명시). ② ★★★★★ **새 규약**: 원-담체 모서리의 **정점 순서 = 어느 호** — `[A,B]` = A→B CCW(축 기준), 상보 호 = 역순 쌍(`derive_edge_curve` 의 (Plane,Cylinder) 팔에 문서화 + tess `sample_edge` 0..τ·validate 미인지 빚 명시). `EdgeKey{Line(무순서 핸들쌍),Arc{cyl,from,to}}` — 호 팔은 담체 직접 진술(cap+측면 — 스캔 불참), `pair_surfs` 는 `Wall::Plane` 만(호를 밀면 현의 선-키 오염). rim 표는 cut 집합(분기 def 의 `(g,cyl,plane)` 쌍)을 skip, Circle/Band 팔은 rim 조회 **앞** cut 검사로 `ArcBoundNotYet`(`MissingSeam`=SuspectedDefect 오진 방지), raise = 면 루프 뒤(**5층 가로채기** — 실체화·rim 실패도 `deferred.unwrap_or` 양보). garbage 가 모서리·면까지(같은 부류·규율). ★★★ 프로브가 store 를 직접 쟀다: straddling = 호 둘 역순쌍 + 현 `[PP]` **하나**(현+가운데조각 용접 — 분할 칸의 약속 실현) + `[v,v]` 는 안 잘린 위 rim 뿐; turned = 호 둘 + 현 0(분기쌍이 판 모서리 경유) — 새 울타리 `an_arc_and_its_complement_are_minted_as_two_ordered_edges` 가 이 표를 단언. ★★ **쌍둥이 울타리가 조여졌다**: 「양-끝-분기 제외」 없이 초록(호가 자기 키를 가지므로) — 제외 삭제, 실측 확정. reject census: 잘린-rim 팔의 셋째 울림 blessed(삼켜진 가드도 기록). red 셋: 키 순서 지움(호가 한 모서리로 접힘) · rim skip off(`[v,v]` 복귀) · 호 팔 스텁(이름 유지). raise 이동의 주석 다섯은 계획이 미리 이름 붙여 grep-0. **남은 것**: 밴드 칸 — `band_loop` 조립(호들+씸) · grouping 분기 팔(shim 은퇴) · 잘린 원 디스크 라벨 · 셸/봉합 ⇒ 「걸친 보스가 빌드」; tess 호 팔·validate 혼합 루프 | ✔ 2026-08-23 |
| M6-2b | **밴드가 호들로 조립된다 — 거절 후의 garbage가 이미 닫힌 셸** (2026-08-23, `19d0375`·이 칸) — ① `CutRim{nodes(θ순서),seam_is_node}` 를 분할에서 **나른다**(`SeamIncident` 이름 판정 그대로 — 재유도·f64 비교 0; disk_labels 와 `Curved` 로 묶어 trace 재넓힘 종결, `work[k]` 번역 동일 자리). **`cut` 집합 은퇴**(Def::Branch 초과근사 → `cut_rims` 존재 = 출처 하나). rim 표 cut 팔이 정점을 세운다: 씸≡분기 **별칭** / S 주조(OnSeam·유리점·실측 tol). 면 루프가 품은 호를 S 에서 가른다 — ★★ wrap 매치는 **방향 있는** 판정(carried `ccw` 로 CCW-쌍: 2-노드 원의 상보 호가 같은 무순서 쌍이라 무순서 매치는 둘 다 갈랐다). ② `band_chains` **선주조**(면 순서 결정적; `edge_for` 이중 캡처가 band_loop 안 주조를 금지) + `band_loop` = `[lo 순회, seam↑, hi 역순, seam↓]`(닫힌 모서리/사슬 한 규칙). ★★★ **hi-잘림은 도달 가능** — 플랜 검토가 「미행사」를 반증: **매달린 보스**(`[4,2,-1]`·+Z — 잘린 원이 hi 끝) 게이트 통과 실측, 픽스처 추가(hi+씸분할 조합만 무픽스처 — 사슬 공유). 양끝-잘림 = chamber 뒤 미도달, 한 줄 정직 거절. **raise = 닫힌-셸 가드 뒤(6층)** — 가드가 호 입력에서 돌고 **조용히 통과**, 밴드 팔 울림이 사라져 census 두 행 raised 가 arrangement+combinatorics 로 복귀. 새 울타리 `the_refusals_garbage_is_already_a_closed_shell`(사용 횟수 전부 2 · 루프 정점-연속 · 반모서리 5/6/5 · 씸만 두 번 반대 sense); F1 은 혼합-쌍 필터로 정밀화(씸 `[lat,lat]` 가 걸렸다). red 다섯(분할 off·별칭 off·사슬 뒤집기·씸 생략·스텁→이름 유지). **남은 것 = 모양 뒤**: validate 루프 규약 · tess `sample_edge`(첫 하류, 자기 doc) · 정지판 제거 ⇒ 「걸친 보스가 빌드」(초록 칸); 양끝-잘림 라벨 · grouping shim 은 인구가 올 때 | ✔ 2026-08-23 |
| M6-2b | **거절이 solid를 다 짓고 나서야 말한다** (2026-08-23, `98e9024`·이 칸) — 초록 칸 조사(정지판 리프트 프로브)가 벽 넷을 실측: W1 grouping(진짜 벽 — held-Err) · W2 validate(★ **커널 무죄** — 초승달 flip·감김 손계산 검산; `loop_winding` 현-Newell이 호-지배 루프를 반대로 읽고 **그 doc이 「M6-2b 오면 여기」라고 예약**; straddling·hung 0 이슈 = sense 대향 초록 판정) · W3 props `UnsupportedBoundary` · W4 tess debug_assert. 이 칸 = W1: ★★★ **`JoinKey { Line(무순서) \| Arc{cyl,from,to}(CCW) }`** — 「직선 무순서·원 순서」 **세 번째 등장**(수리 후보 둘 반증: rim-키 = 캡2+밴드=3 사용자; 무순서 쌍 = 2-노드 원에서 현+호1+호2 한 키 접힘, 사용자 6 — 담체 칸 키-접힘의 노드-공간 재현). 캡 걸음은 자기 `ccw`로, 밴드는 CutRim 순환쌍으로 — 조각마다 캡+밴드=2, **현 키 2로 캡↔캡 결합도 처음 열림**(디곤 고립이 n=2의 정체). ★★ 예측 반증: bvu는 self_touch로 이동하지 않는다 — self_touch 첫 가동이 **자기 doc 규칙(분기 끝점 skip)대로 무울림 기권** ⇒ census 두 행이 **한 사이트**(stopper 자신)로 축소. raise 7층 = `Ok(out)` 직전(`grouping?` 양보; `defs_are_remappable`가 garbage에서 성립 — 실측 후 assert 뒤 배치): garbage = 정점·모서리·면·셸·**솔리드**(live-set 밖, retire는 Ok에서만). 울타리 둘: `the_grouping_joins_across_a_cut_rim`(n==1 ×3) · `the_refusal_leaves_a_complete_solid`(솔리드 1 · outer 셸 == 주조 면 · `mass_props`=`UnsupportedBoundary` 기록). red 둘(등록 off → n==2 · grouping 외래 스텁 → 이름 유지 — F3의 구조 단언은 울타리 자신이지 프로브가 아니다). **남은 것 = 소비자 셋 + 정지판**: validate 혼합-루프 팔 · tess θ 부분구간 · props 원호 세그먼트 적분 → 제거 ⇒ 「걸친 보스가 빌드」(초록 칸). 양끝-잘림 라벨 · grouping shim은 인구 올 때 | ✔ 2026-08-23 |
| M6-2b | ★★★★★ **완료 — 걸친 보스가 빌드된다** (2026-08-23, `37b8056`·이 칸) — 소비자 셋이 `[A,B]`-CCW 규약을 배우고 정지판이 은퇴했다. ① `Circle::angle_of`(규약→숫자의 한 철자) + props `planar_face` 혼합 팔(현-부채꼴 + 원호 세그먼트, `Circle::segment_area/centroid`로 승격해 validate와 공유) — **garbage 솔리드 위에서 부피 32+π/4 정확 잠금**(정지판 유지 채; 한-숫자 사각지대는 straddling의 두 z=2 면 16−π/8·π/8이 부호 양방향으로 가름). ② tess `sample_edge` 부분호 + `triangulate_cylinder` **θ-병합 걷기**(닫힌×닫힌 = 동률-hi-전진 퇴화, 기존 메시 무변화 실측) · validate `loop_winding`의 예약된 팔(현-Newell 2×면적 + 세그먼트 벡터 ×2; digon은 <3점 관문 통과) · **정지판 제거**(`arc_stopper`+witness 가족 5항목·`ClassEdges.split` 삭제; **deferred 배관 = 정지판 소켓** — 7층 사다리는 M6-3 타원이 꽂는다; `ArcBoundNotYet`은 백스톱 둘로 존속) · 울타리 전환(마일스톤 `a_boss_overhanging_the_plates_edge_builds` = validate []+watertight; turned = 부피·오일러 + winding red 스위치; F2/F3/mint 개명·구조 단언 유지·피연산자 은퇴 확인; reject census 두 행 삭제) · census 코퍼스 트리오×KINDS(기존 163행 비트 동일 + 새 9행 = t3). ★ 새 실측: **접촉-절단은 정직 거절**(CoplanarMerge/BranchVertexUnnamed — fuse만 청구; 접촉-절단 인구는 자기 칸), common EMPTY ✓. red: 세그먼트 보정 off(정확 숫자 이동) · winding 팔 off(turned만 cos=-1) · 프로브 전부 python+assert. **남은 것**: 접촉-절단 · 양끝-잘림 라벨 · grouping shim · tess 부분 밴드 · STEP 호 — 각자 인구가 올 때 | ✔ 2026-08-23 |

### 남은 항목 — **비었다** (2026-08-08)

☑★★ **정정의 정정(2026-09-11, 칸 ㊷ `78d770f`): 이름과 자리도 도착했다.** 이 자리에 *"「전부」는
타입의 «내용»이고, 곡면의 «이름과 자리»는 남아 있다"* 고 적혀 있었다 — 그 남은 절반이 그날 닫혔다
(진실 enum = `Surface`, 아레나가 그것을 든다; 열린 항목 8). ⚠★★★ **`SurfaceCache.tol` 은 «안 지음»이
아니라 «동기 미확인»이다**(2026-09-13 실측): 소비자 **0**, 그리고 정점과 달리 **곡면 진실은 정확**해
(`Known` = 유리수 점 셋, 원통 = 유리수 `CylinderDef`) «참 곡면에서 떨어진 잔차»가 애초에 없다 — 발견
정점은 세 평면의 f64 해라 면에서 떠 있지만(그래서 `PointCache.tol` 이 실물, 14 소비자), 곡면은 그렇지
않다. f64 계수의 반올림은 `distance_eps` 로 **질의 때** 나온다(저장 불가). 유일한 후보는 `Through` datum
평면이 세 정점의 tol 을 유도해 무는 것인데 그나마 **하중 없는 지름길**이다. ⇒ 「곡면 실현」 칸은 필드를
짓기 전에 **정말 필요한가부터** 확인한다(대칭이 만든 유령 필드일 수 있다 — 「소비자 없는 배관 금지」).
남은 것은
그것은 타입의 내용이 아니라 **실현**의 몫이다.
⚠★★★★ 그리고 반전이 **새 항목 하나를 드러냈다**: 곡면의 push 문이 캐시를 **받는다**(`push_plane
(cache, points, motion)`) ⇒ 진실과 캐시가 어긋난 상태가 «표현 가능»하다. 형제(`push_edge`)는
**유도**한다. 최종 상태는 「문이 진실만 받는다」이고, 그 열쇠가 `realize_surface` 다(열린 항목 8).

★★★★★ **진실 타입의 «내용»이 전부 최종형에 도달했다.** `PlanePoints`·`Vertex`(⏳ 오늘 `VertexDef`)·`Surface`·
`PlaneName`·`Profile2d`·`Edge`·`FramePlacement` — 그리고 `Motion::Frame { plane }` 은 **이름 없는
평면에도 이미 정확한 정의**다. 남은 것은 타입이 아니라 **판정층의 실현**이고, 그래서 마지막 행이
이행표를 떠나 열린 항목으로 갔다(아래 표의 S5(ii)-2b 줄).

★ **«전부 ✔» 는 아니다**: 위 완료 표의 **`3b`(좌표 재생)는 여전히 ⏸** 이다 — 목적지와 `WitnessPoint` 재생의
무충돌만 확인해 둔 상태로, 재보고 보류다.

| | 항목 | 선행 |
|---|---|---|
| S5(i)-a ✔ | **datum 평면 연산** — `Operation::DatumPlane{def: DatumDef{Stated(SketchPlane) \| Offset{frame, dist}}}` + `OpOutput::DatumPlane{plane, frame}` + `Model::surface_handle_at`(좁은 읽기 접근자, S1 봉인 무손상) + `rebind` 의 `Handle<Surface>` 갈래(`LogCell::Surface`). 새 공개 생성자 없음 — 다섯 `SketchPlane` 생성자가 곧 datum 의 어휘다. 규약: 캐시 법선 `−(진술된 법선)`(근거 둘 — S9 의 781, `frame_sign`), 배치는 `Stated` 면 **무조건 `Named`**(ZX 의 `+u = +ẑ` 가 유도를 금지), `Offset` 은 push 전 3중 정규화(flip→부호 / 세계 되당김 / `dist == 0` 은 `ZeroOffset` 거절)로 «한 평면에 한 핸들»을 지킨다. ★ **앵커 실측**(`tests/plane_anchor.rs`): 기울어진 평면의 저장 `d` 는 앵커마다 최대 **22 ulp** 다르고 그 넷은 **링의 네 점**이라 이 흔들림은 datum 이전부터 있었다; `spans_exactly` 는 다섯 앵커 **전부 false** 라 판정 경로가 앵커를 읽지 않고, 최악 앵커에서도 위상 동일·부피 비트 동일·좌표 **1.1e-15**(ε 의 1/4400). 축 정렬 평면은 `d` 가 곱 하나라 앵커에 **무감각** — 그래서 «먼 원점 + world_xy» 픽스처는 아무것도 재지 못한다 | ✔ 2026-08-07 |
| S5(i)-b ✔ | **`Operation::Extrude` 의 평면 어휘 교체** — `plane: SketchPlane`(값) → `frame: SketchFrame`(핸들+배치+측정된 flip). 구현은 pad/pocket 도로 재사용(`realized_plane` → `swept_profile` → `build_prism(base_cap_surface = Some(프레임의 평면))`), `rebind` 의 `Extrude` 팔이 재고정으로(**7/7 변종이 핸들을 싣는다 — R 의 전제 소멸**), 신설 `SketchFrame::world(model, Axis)` 가 `world_zx` 예외를 한 곳에 가둔다. 호출부 112곳/4크레이트, 단일 커밋. ★ **방향은 `flip` 이 든다**: 같은 평면을 `+n`/`−n` 으로 진술하면 **한 핸들이고 두 프레임 모두 정준 `ŵ`** 를 보고하므로(실측), datum 이 «호출자가 진술한 법선» 에 대해 flip 을 재야 교체가 무행동이 된다. ★★ **C2 가 교체 전에 갈림을 잡았다**: 실현한 축을 다시 lift 하면 정규화가 필요한 축(`(0.6,0.8,0)`→원시`(3,4,0)`, `uu=25`)이 `0.6000000000000001` 로 돌아와 `exact()` 가 뒤집히고 평면이 **조용히 프레임-노드 도로로** 간다(ulp 가 아니라 다른 아레나). 수리 = `RatFrame::of_plane_frame`(유리수로 묻는다, `inv_sqrt_exact`). 내 가설 둘이 죽었다 — 원점 상쇄가 원인이 아니고(단위 축은 이동해도 살아남는다), «`|u_raw|²=1` 이니 정확» 은 **이미 단위인 축에만** 참이었다. C3 차등: 9 평면 가족 × 2 프로파일 = **18/18 아레나 동일 + 노드 수 일치**. census: `c ` 는 **기울어짐 9줄만 이동**(위상 이동 0, 최대 상대 편차 **4.4e-16** = ε 의 1/2000), `wide_planes` 불변, ★ **`seeded_hits` 는 455 불변 — 하락 예측이 틀렸다**(코퍼스에 세계 평면 extrude 가 **0개**이고 그 적중은 `add_cuboid` 의 면이었다; 카운터의 출처를 오독했다). ★ 훅이 못 본 19곳(oracle 은 전부 `#[ignore]`)을 손으로 돌려 **결함 하나를 잡았다** | ✔ 2026-08-07 |
| S5(ii)-1 ✔ | **`PlanePoints::Through([Handle<Vertex>;3])` + 생산자 `DatumDef::ThroughVertices` + 유리수 닫힘 판정** (2026-08-08, 커밋 `ba52b8b`·`3c968f9`·`5c9b3ab`). 이름은 **push 시점에 한 번** 유도한다(`Model::through_points_rat` → `plane_name_exact`) — 그래서 `frame_chain`·far cap·`base_rat`·모든 술어가 **한 줄도 안 바뀐다**. 판정 표는 같은 유도로 증인 삼각형을 **그 자리에서** 만든다(연산 하나 동안 사는 거울이므로 규칙 1 위반이 아니다). ★ **정렬은 키에만, 방향은 호출자의 정점 순서** — `dist` 가 양수 전용이라 순서가 유일한 방향 선택이고, 뒤집으면 «같은 핸들 + 반대 프레임» 이다. ★★ **`transform` 은 핸들을 그대로 두고 노드를 기록한다**(정점은 복제되지 브릿지되지 않는다); 정점이 base 를 정하고 모션이 옮기며 **둘은 더해질 뿐 곱해지지 않는다**. 대가: 그런 datum 위의 솔리드는 정확한 강체 이동에도 **항상** 노드를 얻는다. ★★★★ **가장 조용히 틀릴 뻔한 자리는 컴파일러가 못 본 곳이었다** — 변종 추가가 낸 비망라 에러는 **2개**뿐이고, `transform::points_move` 의 `let`-`else`(원통용 폴백)는 거기 없었다. 원통은 진실이 기하를 **안 들어서** «나를 것 없음» 이 참이지만 `Through` 의 진실은 기하를 **참조로 든다** ⇒ 같은 답을 하면 노드 없는 경로로 가 **캐시만 움직이고 진실은 제자리**에 남는다. `false` 를 반환하게 고쳤다 | ✔ 2026-08-08 |
| S5(ii)-2a ✔ | **증인 없는 평면의 구간·고정밀 계수** (2026-08-08, 커밋 `7acf5c5`). `plane_iv_through`·`plane_hp_through` — 동차점 셋의 **사영 join**. **차수 15 → 9**(분모 털고 아핀 외적이 아니라 3×4 의 네 소행렬식), `cramer_iv` 의 **쌍대**라 새 산술 0. ★★★★ **초안은 설계가 금지한 연산을 제안했다** — `Dvec/D` 로 나눠 아핀 좌표를 만들려 했고, 문서가 *"나누면 무리수가 되고 오차가 낀다"* 며 이미 금지한 것이었다. **타입이 이미 강제하고 있었다**: `Iv` 에는 나눗셈이 없고 `HpIv::div_exact` 는 반경 0 을 요구한다. ★★★ **배율의 부호는 값 안에서 없앤다**(실을 자리가 시그니처에 없다 — 실을 곳 없는 값은 아무도 안 쓴다). ★★ **실측이 단계를 갈랐다**: 깊이 1 은 필터가 살지만(미결 0/800, 최악 상대 반경 6.3e-10) **깊이 2 는 차수 81 이라 계수가 `f64` 밖 8/8** ⇒ 상승 전용. `NaN` 반경이 우연히 옳게 굴러가는 것에 기대지 않고 `None` 으로 명시한다. 생산자 없음(census 비트 동일), 잠금은 `three_planes_big` 과 같은 모양의 **차등 테스트** | ✔ 2026-08-08 |
| ~~S5(ii)-2b~~ | ★★★★★ **이행표를 떠났다 — 타입이 아니라 «능력» 이다**(2026-08-08). 여기 «배선» 이라 적혀 있던 것이 틀렸다: `WorkingPlaneDef::Through` 는 **생산자를 가질 수 없다**. 혼합 프레임 datum 은 이름이 없고 ⇒ `frame_chain` 이 사퇴하고 ⇒ `SketchFrame` 이 없고 ⇒ base cap 이 못 되고 ⇒ 판정 표에 **들어오지도 않는다**(`a_plane_with_no_name_cannot_host_a_sketch` 가 그 사슬을 실행한다). 여기 적혀 있던 다섯 질문(인덱스 공간·`standard_for`·interning·`tri_pt3` 전성·깊이 제한)은 전부 그 벽 뒤다. → **열린 항목 16 «이름 없는 평면의 프레임»** | — |

★★★ **S5(ii) 앞에서 «먼저 잰다» 를 했고, 그 실측이 단계를 잘랐다**(2026-08-07~08).

- **재기 전 예측이 틀렸다**: 발견 좌표의 폭은 최대 **59비트**, 127 초과 **0건**, 연산을 40회
  쌓아도 이름 폭은 **상수**다(모션은 계수를 안 건드리고 불리언은 `push_plane` 을 안 부른다).
- ★★★ **그러나 그것으로 규칙 1 을 반박할 수는 없다.** 한때 «폭이 반박됐으니 변종을 미루자» 는
  결론까지 갔고, **이 문서의 반증표가 그 추론을 이미 금지하고 있었다**(*"코퍼스 최대값이 상한"*
  ✗ — **폭은 타입에서 유도하라**). 규칙 1 은 타입 주장이다: `PlaneName` 이 `Narrow | Wide` 인
  이상 세 계수 행의 Cramer 인 발견 좌표는 **임의로 넓을 수 있고**, 그 사실을 위해 만든 것이
  `MeetPoint::Wide`(음성 대조 201비트)다. 실측이 정한 것은 **목적지가 아니라 단계**다.
- ★★ 그 대신 실측이 **두 결함**을 새로 드러냈고 그것이 1단계의 동기다: (a) 세 정점을
  **좌표로 적으면 다른 평면**이 나온다(기울어진 인구에서 **220/220**, 축정렬 음성 대조
  552/0), (b) **혼합 프레임 정점은 어떤 프레임에서도 유리수 좌표가 없다**(회전 불리언에서
  20중 12) — 이것이 «값 대신 핸들» 의 진짜 근거이고 **폭이 아니라 프레임**이다.

**답이 전역일 필요가 없다는 것**(Q5 의 per-variant 교리)은 1단계가 실행으로 확인했다 —
`Through` 는 datum 평면만 받고 `Known` 경로는 한 줄도 안 바뀐다.

★★★★ **정정**: 여기 적혀 있던 *"`Through` 는 모션 아래 «움직인 정점을 가리킴» 으로 수송된다"*
는 **틀렸다**(구현이 반증). `transform_solid` 는 정점을 **복제**하므로 원래 정점을 가리키는
datum 은 새 복사본을 추적하지 않는다. 실제 규약은 **핸들을 그대로 두고 모션 노드를 기록**하는
것이고, 정점이 base 를 정하고 모션이 옮기며 **둘은 더해질 뿐 곱해지지 않는다**. 동일성은
그래도 성립하지만 그 근거는 «정점을 따라감» 이 아니라 **interning**(같은 셋 → 같은 이름 →
같은 핸들)이다.

★ 판정층 개명은 «각 타입을 처음 만지는 단계에 얹는다» 로 적어 뒀으나 이행표가 비어 얹을 단계가
없어졌고, **자기 단계로 집행됐다**(2026-08-08, 4커밋 + census 비트 동일 ×4). 결과와 항목별
처분은 §판정 이름 규칙의 대응 목록에 있다.

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
  기록 층(이행표·열린 항목)이 자라고 그림 층은 얼어 있는 구조라, 행만 쓰면 그림이 조용히 뒤처진다 —
  실측: §캐시의 `Model` 그림은 2026-08-05 에 고정됐고 그 뒤 동일성 표 둘(8/09·8/17)이 **행으로만**
  들어왔다. 같은 사고가 곡면 그림(칸 ㊷)·동일성 그림(2026-09-12)에서 두 번 났다.
- **어느 크레이트에 두나** — 「유도된 값 + 그 산술」은 `nacre-scalar`(최하단; `nacre-cip` 이 닿아야
  하므로), 「아레나 항목의 진실」은 `nacre-topo`. 구조적 강제: **`nacre-cip` 은 `nacre-topo` 에 의존하지
  않는다** ⇒ 판정이 쓰는 타입은 전부 topo 아래에 있어야 한다(`PlaneName` 이 scalar 에 사는 이유).
  본문은 design.md 「크레이트 구조」.
- **병렬 불변식**: 병렬 구간에서 **곡면 push 금지**(재생 결정성). ⚠ **정정(2026-09-11)**: 이 규칙은
  `push_surface` 를 이름 짓고 있었는데 **그 함수는 없다**(S6b 에서 `push_surface(_with_points/
  _unrecorded)` 사망, 실측 0건). 오늘의 문은 `Model::push_plane`·`push_cylinder`(사설 `push_raw` 로
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
| 혼합 정의 datum (정점 2 + 좌표 1) | `Through` 는 핸들 3 전용, `Known` 은 값 3 전용 | `Through` 원소를 `PlanePointRef(At \| Vertex)` 로 일반화 — 평면의 점 데이터 수준이라, 기각된 `VertexDef::At`(정점 아레나)과 다른 자리다 |
| 중간 평면(midplane) | 같은 프레임 평행 면 사이는 유리수 평균이라 `Known` 으로 계산 가능 ✓. **프레임이 다른** 두 면 사이는 세계 좌표가 무리수 + 어느 정점도 안 지난다 | 필요해지면 정의 변종 `Mid([Handle<Surface>; 2])` |
| 스케치 내부의 정확한 각도(정육각형·30° 변) | 꼭짓점에 √3·tan30 — `[Rat; 2]` 로 못 적는다. 앱이 계산한 f64 의 십진수가 진실이 된다(§design 6.0 — 의도된 동작) | **정확한 각도는 좌표가 아니라 모션으로 표현된다**(회전·프레임) — 결함이 아니라 이 숫자 시스템의 정의적 성질. 여기 명시해 둔다 |

검사했고 막히지 않는 것: datum 이 참조한 솔리드의 이동(이동본 정점으로 재지시), superseded
참조(append-only), 원형 패턴(`360/n` 은 항상 유리수), 깊은 체이닝(비용만 — 실측 회전 3200회),
datum 낀 불리언(비용만 — 아래 1), 4+평면 정점(validate 몫), 오프셋의 오프셋(연산당 1단),
STEP 출력, undo/replay.

---

## 열린 항목 — 측정할 것과 알려진 절벽

0. ✔ **닫힘 (2026-08-09, `696dbd3`·`3bf46e0`·docs) — 정점 풀이의 산술 천장이 사라졌다.**
   원문(2026-08-07): `three_planes_rat` 의 여인수 전개는 덧셈·뺄셈이 `a.num·b.den ± b.num·a.den`
   을 약분 전에 만들며 넘치므로, 그 `None` 이 「점이 `Rat` 에 안 든다」와 「중간값이 넘쳤다」를
   뭉쳤다(기울어진 십진 가족에서 받은 정점의 2/3 거절, 그중 100%가 드는 점). 지은 것과 배운 것:
   - **수선은 문서의 처방 그대로**: `three_planes_big` 에서 정수 코어(`three_planes_int`)를
     분리(순수 리팩토링 커밋 — blame 격리)하고, `three_planes_rat` 에 `.or_else`(행별 분모
     걷기 → 코어 → `Narrow` 복귀). 이제 `None` = «유일한 만남 없음» 또는 «점이 정말 `Rat`
     밖» — 산술은 원인이 아니다. 비용은 decline 경로에서만.
   - ★★★★ **2026-08-08 의 «프로덕션 이득 아마 0» 정정이 반증됐다**(반증표). 그 정정은 발견
     정점만 잰 인접 명제였고, reuse 에 도달하는 인구는 **움직인 솔리드의 구성 코너**다 —
     probe 실측: wf 프레임 프리즘 코너 **8/8 이 rat 에서 넘치고 8/8 이 드는 점**.
     `solid_points` 는 첫 실패에서 솔리드 전체를 포기하므로, 십진 프레임 피연산자는 클래스
     reuse(6.3)를 통째로 잃고 있었다.
   - **실측(수선 후)**: 프레임 프리즘 `solid_points` = Some(8), 각 점의 replay 좌표가 저장
     좌표와 **비트 동일**(한 점 한 실현 도로 — `a_framed_prisms_corners_solve_for_reuse`);
     census-형태 pocket 불리언에서 새 상담 방향(도구-소유 클래스가 표적의 점을 묻는다)이
     **1/11 클래스 발화**(Empty — 도구 far cap). census 는 **비트 동일** — 발화한 계획의
     산출이 배열 산출과 비트까지 같았다(그 클래스는 어차피 결과에 면을 내지 않는다).
   - **정직한 잔여**: 벽-프레임 도구(pocket 의 도구 모양) **자신**은 여전히 decline — 그
     꼭짓점의 담체가 혼합 프레임(밑캡 = 벽의 프레임, 벽 = 자기 스케치 프레임)이라 열린 항목
     13/14 의 인구다. 이 항목의 천장과는 다른 벽.

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
3. ~~**`from_decimal` 창 밖 치수**~~ — 셌다(S3): 폴백에 기대던 프로파일은 스위트 전체에서
   **0건**(유일한 창 밖 좌표는 의도적 부정 테스트 하나였고 구성 시점 에러 잠금으로 대체).
   전 스위트 그린 + census 150줄 비트 동일이 증거.
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
8. ✔ **닫힘 (2026-09-11, 칸 ㊷ `78d770f`) — 최종 개명의 두 절반이 다 집행됐다.**
   ✔ 판정층(2026-08-08): `Mag`·`Approx`/`HpApprox`(+필드 `value`/`error`)·`WitnessPoint`·
   `WorkingPlane` — 4커밋, census 비트 동일 ×4, 항목별 처분은 §판정 이름 규칙의 대응 목록.
   ✔ **아레나 반전**(2026-09-11): `surfaces: Store<Surface>` 가 **진실**을 들고
   `surface_cache: Vec<SurfaceCache>` 가 실현을 든다. `SurfaceTruth`→`Surface` 개명 동반.
   두 접근자의 서명은 그대로이므로 `.surface(`·`surface_truth(` 호출 철자는 무변이고,
   바뀐 것은 «어느 쪽을 돌려주는가»다.

   ★★★★ **개명과 반전은 한 커밋이었다 — 갈랐던 계획이 타입 체계에 반증됐다.** 진실을
   `Surface` 로 부르는 순간 topo 자신의 `use nacre_geom::Surface` 가 **E0255** 이므로, 갈라
   가면 `Handle<Surface>` **122자리**를 `Handle<nacre_geom::Surface>` 로 적고 되돌려야 한다.
   개명은 census 를 움직일 수 없고(의미가 없는 변경), 두 `Surface` 의 혼동은 **하드 오류**
   (geom 은 tuple 변종, 진실은 struct 변종 ⇒ E0532)이므로 합쳐서 잃는 것이 없다.

   ★★★ **철자 규칙(영구)**: 두 `Surface` 가 한 파일에 오면 **맨이름은 진실**(= `Handle<Surface>`
   가 이름 짓는 것)이고 실현은 `nacre_geom::Surface` 로 적는다. 실측: 두 이름이 다 필요한 파일
   **17**, 그 안에서 경로를 붙인 geom 패턴 **99**(32는 이미 그 형태), `Handle<Surface>` **무변**.
   ⇒ 반대로 하면 한 파일에서 `Handle<Surface>` 는 진실을, 맨 `Surface::Plane(p)` 는 캐시를 뜻해
   **한 단어가 표지 없이 두 가지**가 된다.

   ☑ **`set_plane_points_for_test` 는 «다시 열지» 않고 사망했다.** 유일한 픽스처가 깊은 삼중항을
   `push_plane`(프로덕션 interning 문)으로 **먼저** 진술하고 `add_cuboid` 의 캡이 그 위로
   interning 되게 한다(그 interning 을 테스트가 단언한다) ⇒ **`Store` 에 test-gated 문도 내지
   않았다.** 「의도적으로 어긋날 수 있는」 곡면이 이제 이 저장소에 없다.

   ⚠★★★★ **감사가 잡은 것 — 반전이 «가드 하나»를 조용히 떨어뜨렸다** (2026-09-11, 같은 날 수리).
   `Store::get` 은 `cfg(debug_assertions)` 로 **cross-store 가드**를 든다(*"Handle was minted by a
   different Store"*). 반전 전 `Model::surface` 는 `self.surfaces.get(h)` 였으므로 그 가드를 공짜로
   받았는데, 캐시를 `Vec` 으로 색인하면서 **사라졌다**(실측: 남의 모델 핸들이 패닉 없이 다른 칸을
   읽었다). ⚠ 그리고 **형제 둘은 애초에 없었다** — `vertex_point`·`vertex_tol`·`edge_curve` 전부
   통과. ⇒ 「형제를 따른다」가 **가드를 잃는 방향**이었다.
   ☑ 수리는 부류 전체로: 캐시 읽기가 자기 `Store` 에 먼저 물어본다(`Model::debug_guard`, 디버그
   전용 ⇒ release 비용 0). 잠금 `a_foreign_handle_cannot_read_a_cache` 가 **다섯 문 전부**를 문다
   (가드를 `surface` 에서·`vertex_point` 에서 각각 빼 보고 그 이름으로 빨개지는 것을 확인).
   ★ 확장이 **기존 테스트를 하나도 깨지 않았다** ⇒ 그 구멍은 있었지만 아직 밟히지 않았다.

   ★★★ **실측한 계기의 분업**(위반을 심어 확인): 캐시 순열 → census 가 움직이고 패닉 ·
   **모든 유도값이 같은 진실 재진술**(평면 삼중항의 순환 회전) → census **비트 동일**이고
   `arena_sig` **초록**, 진실 다이제스트만 본다 ⇒ 그 다이제스트가 존재하는 이유가 그 한 부류다 ·
   평면 캐시를 든 원통 진실 → 변종 스윕보다 **먼저** rim circle 유도가 잡는다(캐시를 읽으므로).
   ⇒ 「변종 일치」는 검사만 되는 게 아니라 **유도에 하중이 걸려 있다**.

   ⚠ 이 항목은 *"`SurfaceCache` 실형이 **판정 통합에서** 생길 때 개명·반전을 **한 번에** 한다"*
   고 적고 있었다. 2026-09-11 이 그 조건(판정 통합 ✗ → 곡면 실현)과 그 묶음(반전 ✔ 먼저)을 둘 다
   반박해 둘로 갈랐고, 반전 쪽이 위와 같이 닫혔다. 표지는 번호가 아니라 **내용**으로 붙인다 —
   §관문 규칙이 이미 다른 뜻으로 「8b」를 쓴다. 남은 절반:

   - **「8·캐시 실형」(`SurfaceCache`) — 조건은 「판정 통합」이 아니라 「곡면 실현」이다.** 근거는 이
     문서 자신의 문장: *"`Through` 평면의 `SurfaceCache` 는 **세 정점의 실현에서 유도된다**"* ⇒
     `tol` 은 판정이 아니라 **실현**이 알게 된다. ⏸ 모양은 미정이다(§캐시의 주석 — 원통에 4계수형이
     없다). ⚠ 그리고 **소비자 유무를 그 칸에서 확인한다**: 오늘 판정은 `distance_eps`+증인 tol 로
     자기 방식이고 STEP 은 tol 을 안 쓴다 ⇒ 「소비자 없는 배관」 금지에 걸릴 수 있다.

   ⚠★★★★★ **반전이 드러낸 것 — 문이 캐시를 «받는다»(남은 절반의 진짜 목표).** 실측:
   `push_plane(cache: nacre_geom::Plane, points, motion)` · `push_cylinder(cache, def, motion)` ·
   `push_raw(truth, cache)` — **캐시를 호출자가 준다.** 캐시의 정의는 「진실에서 언제든 다시
   계산할 수 있는 것」인데, 받는 문은 **둘이 어긋난 상태를 표현 가능하게** 만든다. 그래서
   `points_coverage` 의 *"cache and truth disagree about what this surface is"* 검사가 **존재해야
   했다**. 형제는 다르다: `push_edge` 는 캐시를 `derive_edge_curve` 로 **유도**하므로 간선은 그
   어긋남이 **표현 불가**다.
   ★ 그리고 그 칸이 **모델 캐시에 처음으로 정밀도 계약을 준다** — 오늘 캐시는 «만들 때 나온 f64»
   (사다리 없음, 기울면 최대 4 ulp)인데, `realize_vertex(v, NearestF64)` 의 출력을 쓰면 «정확 반올림 +
   축별 경계 `[Mag;3]`»가 된다. 최초판 887행의 약속(*"STEP 은 정확 반올림해 낸다"*)이 그때 닫히고,
   `PointCache.tol` 의 뜻이 «잔차 하나»에서 «축별 경계 셋»으로 바뀐다(§캐시 «세 오차» 표 갱신 대상).
   ⚠★★★ **정정(2026-09-13) — 「캐시를 판정 1단의 씨앗으로」는 이득이 작다, 그리고 이유는 «캐시가 정확성을
   진다»가 아니다.** 여기 그 대가를 그렇게 적었는데 부정확했다. 캐시의 오차는 엉터리가 아니다(㊸ 뒤면
   realize 통로에서 나온 건전한 `Mag`). 진짜는 **캐시가 실현의 «출력»(나눠 반올림한 f64)이고 판정의 정확
   단계는 «입력»(유리수 정의 = 세 평면)을 필요로 한다**는 것 — `PointCache` 엔 `base:[Rat;3]` 이 없어
   `WitnessPoint::at_with_tol` 에 애초에 못 들어가고(반올림한 f64 를 유리수로 도로 들면 «다른 점»의 부호를
   정한다), 씨앗이 될 수 있는 건 1단 f64 필터뿐인데 그건 `base.to_f64()` 세 번이고 비싼 hp 는 이미
   메모돼 있다. ⇒ 성능 이득 ~0. 구조 이득은 **위 계단 2 통일**로 따로 받는다(§이름과 interning 위,
   열린 항목 23).
   ⇒ **옳은 최종 상태는 「문이 진실만 받고 캐시를 유도한다」**이고, 그러면 그 변종 검사가 *필요
   없어진다* — 그것이 곡면 실현 칸의 판정 기준이다. ⚠ 유도에는 `realize_surface` 가 필요하고
   (`Known` 의 세 `Rat` 에서 법선을 유리수로 유도하면 넘친다 ⇒ `PlaneName{Narrow|Wide}` 경유),
   그래서 반전 칸은 이것을 닫지 않았다.

9. **`Inexact` 소멸의 두 반증은 지반이 제거됐다(S6a)** — #28(축만 든 호출자)은 `from_axes`
   리프트가 닫았고(십진 진실 def + S2 이름 + S4 wide 프레임), "정확한 형태가 없는 평면
   0.07%/모델 25.8%" 실측은 **S2 이전 수치**(원인이 이름의 i128 넘침이었고 그 원인이
   죽었다)다. 스위트-내 재측정치는 `points_coverage` 관문의 0(생산 경로 전수에서 점 없는
   live 평면 face 없음); 코퍼스(OCCT 스위트) 규모의 재측정은 S6b 착수 시 한 번 더 돌려
   기록한다.
11. **노드 생략의 더 강한 게이트** — 오늘의 게이트는 `exact()`(축이 정확 유리수 직교로
   리프트되는가)다. `PlaneFrame`+`inv_sqrt_exact` 로 «실현이 정확 f64 에 떨어지는가»를 직접
   묻는 더 강한 게이트가 가능하지만, **표현식을 바꾸면 노드 인구가 움직인다** — 두 ★★ 주석의
   경고 그대로, 교체는 census 관문 동반 필수(S9 에서 기록만).
12. **정점 캐시엔 «버리고 재생» 보증이 없다(S7)** — `rebuild_edge_cache` 의 정점판을 짓지
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
   ⇒ 남은 것은 「보증이 없다」는 **사실**이지 그 근거들이 아니다. 정점판 재생(=`refine_caches`)은
   **「8·반전」 뒤 곡면 실현과 함께** 열린다 — 그때 캐시가 실현의 메모가 되고 계약이 한 문장이 된다
   (*"모든 `vertex_point(v)` 가 `realize_vertex(v, NearestF64)` 와 같아지게 한다"*).
13. **reuse 의 발견-정점 갈래는 아직 포기다(S7)** — `solid_points` 는 측정 좌표를 만나면
   `None`(→ Arrange). `nacre-cip` 의 `ImplicitPoint`(세 평면의 암시적 점)로 갈아타면 융합
   fold 의 클래스 재사용이 불리언 **결과** 피연산자에도 걸린다 — 판정 기계 교체라 행동
   변화이고, 별도 측정과 함께 여는 항목.
14. **혼합-프레임 정점은 reuse 가 답하지 못한다(S7 기록)** — 호출자가 세계 좌표로 명시한
   밑캡 위의 프레임-스케치 코너는 세 담체의 모션이 갈려(둘은 프레임, 하나는 세계) 유리수
   pullback 이 없다. 옛 `Origin` 길은 base 정점으로 답했지만 그 정점이 S7 에서 소멸했다.
   실측 인구는 `the_def_road_answers_for_the_populations_it_can_name` 의 ④ 가 핀한다.
15. ✔ **닫힘 (2026-08-09, `0bd6d29`·`378d86b`·docs) — 술어가 이름의 정수를 읽는다.**
   원문(2026-08-08): 술어가 `Wide` 계수를 못 읽어 wide 이름 평면은 정확 지름길을 전부 잃고
   톨러런스+상승으로 갔다(실측 3.2배). 지은 것과 배운 것:
   - **그릇은 `Expansion` 이 아니라 BigInt 다**(원문 처방 반증 — 반증표). scalar 에 정수 부호
     술어 셋(`int_plane_side`·`int_cmp_coord`·`int_dir_sign`, f64 쌍둥이와 같은 규약) +
     `PlaneName::coeff_ints`. wide 의 참값 잠금은 **스케일 불변성**(narrow 계수 × 200비트
     소수 = 진짜 >i128 입력, 답은 narrow-Shewchuk 판이 안다) + 음수 배율 = 방향 민감성의
     음성 대조(`plane_side` 는 j 만, `cmp` 는 무엇도, `dir_sign` 은 행마다 뒤집힌다).
   - ★★★★ **wide 이름은 «움직인 평면» 이 아니라 «세계 이름» 이다** — 계획이 두 번의 검토를
     통과하고도 품고 있던 오류를 구현 중 probe 가 반증했다(반증표). 담체가 **발견 정점**
     (사슬 없는 세계 유리수)이라 datum 의 이름은 세계를 말하고, 그 평면을 처음 만드는
     불리언의 표는 **혼합 프레임**(세계-wide 캡 × 프레임의 벽 × `FrameWide` 벽 — 뒤의 것은
     프레임 자체가 무리수라 정수 기술이 **존재하지 않는다**). 게이트가 실제로 여는 인구는
     ① **세계-게이트**(전 평면 무이동·전부 이름·하나 이상 wide), ② **공유-모션 게이트**
     (전 평면 한 사슬 — wide 면을 품은 솔리드를 통째로 돌린/거울한 2세대; 홀수 미러 패리티
     보정은 계수에서 **x 유지·y/z/d 부정** = `parity·C`, «x 만 부정» 은 점-규약이지
     계수-규약이 아니다 — 외적이 pseudovector 라 `det(C)` 가 한 번 더 붙는다).
   - ★ **③ 세계-진술 문**(2026-08-24, `e45f290` — 문이 셋이 됐다). ①②는 **평면들이 한
     프레임을 공유**할 것을 요구한다. 2세대 배열(x로 붙인 결과를 다시 y로 옮긴 몸체)은 한
     코너에 **서로 다른 사슬**의 담체를 모으므로 어느 쪽에도 안 맞았다 — 그런데 그 담체들은
     저마다 세계를 **정확히** 말할 수 있다(안 움직였거나, 사슬이 평행이동으로 접히면 유리수
     이름이 그대로 옮겨간다). `Model::world_plane_name` 이 그 진술이고, `vertex_meet` 은 그것을
     **폴백으로만** 탄다 — ①②가 `None` 을 준 자리에서만. 순서가 계약이다: 먼저 태우면 오늘
     답이 있는 코너의 점 철자가 바뀐다(폴백으로 두면 census 220행 비트 동일, 실측).
     한계: **이동한 Wide 이름**은 좁은 그릇이 없어 접을 수 없고(안 움직인 Wide 는 그대로
     세계 이름이라 문제없다), 회전·거울 사슬은 평행이동으로 안 접힌다 ⇒ 둘 다 `None`.
     답에는 leaf 가 없다 — **이름을 옮긴 것이지 점을 옮긴 게 아니므로** 결과가 이미 세계다.
     쌍둥이 게이트 `reuse::solid_points` 는 같은 규칙의 두 번째 스펠링을 들고 아직 ①②만
     안다(기권 → Arrange, 느릴 뿐 틀리지 않는다). `through_meets` 는 공유-프레임 정점과
     세계 정점이 섞인 삼중을 계속 거절한다 — 각각 정확해도 **다른 프레임**의 점이라 옳다.
     코퍼스: census `xy` 가족(x-후-y·y-후-x·3셀·회전 음성 — 산술을 고정하되 **문 없이도
     빌드된다**, 실측). ★ 문을 실제로 요구하는 행은 `xy needsworld`와 통합 잠금
     `a_two_by_two_grid_fuses`뿐이고, 그 셀의 숫자(포켓 하나 + 관통 구멍 둘)는 **못 박힌 것**이다
     — 이 인구는 **단조가 아니다**(포켓 하나 + 구멍 1~4는 빌드, 같은 포켓 + 구멍 4·7은 거절).
   - **σ 접기**: 정준 이름은 방향이 없고(interning `flipped` 실측 1,916) 저장 방향과 반대일
     수 있다 — `name_stored_ints` 가 구성 시 한 번, 증인 base 외적(계약: 이름 평면 **위**의
     정확한 세 점 — debug_assert 로 박음) 대 이름의 BigRational 부호 × `frame_sign` 으로
     접어 **`base_coeffs` 와 같은 방향**을 들게 한다 ⇒ 세 술어의 rescue 는 기존 팔에 자모
     이식, 새 부호 추론 0.
   - **실측**(`wide_datum_cost`, 단일 스레드 — ★ 두 측정을 병렬로 돌리면 전역 `climb_census`
     가 서로 오염된다): 1세대(혼합 표) narrow 110 불변·nameless 434 불변(음성 대조 ✓),
     wide **454 → 330**(−27% — 혼합 표 안의 세계-게이트 부분집합). **2세대**(전부 세계 이름
     + wide 슬라이스, 세계 도구로 재절단): narrow 대조군 139 대 wide **20** — **wide 가
     narrow 보다 7배 싸졌다**(정수 길은 wide 참여 질문에만 열리므로 — census 안전 게이트).
   - **정직한 잔여**: ① 1세대 혼합-프레임 질문(330의 대부분)은 이름 길이 **원리적으로** 못
     연다(`FrameWide` 무리수 프레임) 또는 프레임-간 이송이 필요하다(유리수 사슬 한정 —
     별도 항목감). ② narrow 의 계수-불일치 인구(2세대 대조군의 139 — f64 지름길이 깊은
     유리수 증인의 반올림에 사퇴)는 census 가 움직일 수 있는 행동 변화라 **별도 후속**.
     ③ `cmp_coord` 의 공유-모션은 범위 밖(회전축 추론은 `cancel_cmp_coord` 의 것).
   - 경로 증명: `a_starved_judge_answers_only_through_the_name` — 8비트로 굶긴 판정이
     2⁻¹⁰⁰ 분리 질문에 ±1 을 답하면 이름 길뿐이고, 이름을 빼면 같은 판정이 0 으로 돌아온다.

16. ★★★★★ **이름 없는 평면에 프레임을 주는 일 — S5(ii)-2b 가 이행표를 떠나 온 자리**
   (신설 2026-08-08). 혼합 프레임 datum 을 지원한다는 것은 **판정층 배선이 아니다.**
   벽은 둘(프레임 → 판정 표)이고, 혼합의 원인도 둘이다(정점들이 서로 다른 프레임 /
   한 정점의 담체가 갈라짐 — `pure_frames` 실측이 갈랐다).

   ### ✔ 16-1 (2026-08-09, `cdff866`..`9922b1f`) — 순수-혼합 population 이 끝-대-끝으로 열렸다

   **정점들이 서로 다른 프레임**(각 정점은 자기 프레임에서 유리수)인 쪽:
   - **cip**: `HpApprox::div`(구간 분모 — `denom_lo` 로 0에서 떼고 몫 법칙)·`inv_sqrt`(도함수
     상한, 분모가 2의 거듭제곱이라 `Mag::over` 가 정확) + **`MoveNode::FrameThrough`** — 평면의
     정의점 셋(사슬 이종 허용)을 싣고 정준 기저를 **구간으로 유도**한다(수선의 발 — 동결 규약
     그대로, 그래서 구간 나눗셈이 필요했다; «원점 = points[0]» 대안은 규약 동결이 기각).
     분기(ẑ×n vs ŷ×n)는 **고정 128비트 실현에서 한 번** 판정해 저장 — 정밀도마다 다시 정하면
     기저가 흔들려 프레임이 아니다. ★ 계획의 `negate` 필드는 검토가 죽였다: 정준 원점은
     부호·스케일 불변이라 전역 부정의 자유는 `flip` 하나뿐이다. ★★ «`MoveNode` 의 `PartialEq`
     는 소비자 없음» 은 **틀렸었다** — `shared_base` 의 전-노드 비교가 하중을 진다 →
     `WitnessPoint` 가 **정의 동등성**(base+chain, 캐시 제외)을 얻었다.
   - **topo**: 이름 키가 못 받는 진술은 **진술 키**(`surface_through_ids` — 정렬 삼중항+모션)로
     intern. 진술 동일성이지 기하 동일성이 아니다 — 후자는 규칙 6 대로 술어가 매번.
   - **ops**: 원인 3분기(공유 프레임 → 이름 길 그대로 / 각자 순수·프레임 상이 → **수용, 이름
     없음** / 담체 갈라짐 → `VerticesInMixedFrames` 유지). **검증이 push 보다 먼저**
     (거절-후-커밋 금지), 실패는 `ThroughFrameUndecided` — `CollinearVertices` 도
     `DegenerateGeometry` 도 아니다(구간이 0을 못 벗어난 것은 퇴화의 증명이 아니다).
     `frame_chain` 제3 도로 + `collect_planes` 이종-사슬 가지(`rotated=true` — 신호의 뜻은
     «정확 f64 기술 없음»). 공유 유도는 `through_witness_points` 하나 — 검증한 것과 프레임
     되는 것이 어긋날 수 없다.
   - **실측**: 프리즘이 Common 을 **통과**, 부피 = 스케치 진술 그대로(0.5), validate 클린,
     재구축 결정적(같은 진술 → 같은 아레나 슬롯). 판정 비용 **416 climbs**(wide 이름 454·narrow
     110, 고갈 0) — **Wide 와 같은 자릿수, 절벽 아님** ⇒ `SurfaceCache` 선실현 레버는 안
     당긴다. 어휘 도달률(삼중항): 프레임-위-프레임 10.7%→**39.3%**, 두 번 쌓으면 3.2%→**27.7%**;
     회전-대-세계는 0 그대로(전부 straddle — 아래). 판정층 새 기계 **0** — `plane_iv`/`plane_hp`
     는 원래 세 점을 각자 실현했고, `standard_for` 는 그 점들이 표에 있어 자동 집계했다.

   ### ✔ 16-2 (2026-08-09, `4c47747`..`5a027c9`) — 걸친 정점의 datum, join 의 첫 호출자

   **`JudgedPoint { Pure, Meet }`**: 판정 프레임의 정의점은 자기 프레임의 유리수(16-1)이거나
   **세 담체 평면의 교점**이다 — 좌표는 어디에도 없고, 담체들의 Cramer(`cramer_hp`)가 낸 동차점
   `[Dvec:D]` 로 실현되며, 나눗셈은 캐시 anchor 한 곳뿐이다(생성자가 rung 에서 그 `D` 부호까지
   증명하므로 안전). 동차점 셋의 평면이 **2a 의 사영 join(`plane_hp_through`)** — 고아로 2단계를
   기다린 그 기계가 첫 호출자를 얻었고 dead_code allow 가 떨어졌다(`plane_iv_through` 의 f64 판만
   16-3 의 필터를 기다린다).
   - **잠금 = Pure-대-Meet 차등**: 순수 정점은 자기 담체 셋의 Meet 로도 적을 수 있으므로, 같은
     세 정점의 두 표기가 한 기저를 내야 한다 — 정지·37° 회전 두 판 모두 단언(회전 판이 어려운
     쪽: 담체는 평면으로 돌고 순수 표기는 점으로 돌아, 담체 실현·Cramer·join·정규화 전부가
     맞아야만 만난다). + 걸친 정의의 tol-내 실현, 결정성, `D` 미결 거절, 한 직선 공유 담체 거절.
   - **producer**: 원인 분기 재구성 — 1패스는 프레임만 묻고(풀이는 도로가 정해진 뒤: 같은 실패가
     도로마다 다른 뜻이다 — 너무 넓은 좌표는 이름 길을 죽이지만 Meet 은 좌표에게 맞기를 요구한
     적이 없다), straddle 은 `Nameless` 로 합류. **`VerticesInMixedFrames` 의 남은 유일한
     population = 담체가 이름 없는 datum 인 정점**(깊이 — 16-3 의 질문; 인구가 생기면 라벨 분리를
     잰다). 검증은 push 전, 거절은 `ThroughFrameUndecided` 그대로.
   - **실측**: 걸친 datum 이 수용(이름 없음·진술 interning)→돌출→validate→재구축 결정성까지
     끝-대-끝. **datum 어휘의 벽이 이 층에서 완전히 열렸다** — 기하학적으로 멀쩡한 삼중항의
     분류-수용률이 전 population 100%(회전 불리언 4.9+95.1, 프레임 10.7+28.6+60.7, ×2
     3.2+24.5+72.3; blocked 0; 음성 대조 둘 다 혼합 버킷 0).

   ### ✔ 16-3 (2026-08-09, `7c49573`..) — **판정 표는 Through 변종이 필요 없었다 — 항목 닫힘**

   판정 표의 계약은 S6b 이후로 **«그 평면 위의 정확한 세 점, n_out 감김»** 이지 «면 코너» 가
   아니었고, 판정된 평면은 그런 점을 **정의상** 갖고 있다: 자기 정준 프레임의 probe
   `(0,0,0)·(1,0,0)·(0,1,0)` — 원점은 수선의 발(평면 위 정확), û·v̂ 는 평면 안, `frame_chain`
   이 후속 모션까지 붙이고, 상승은 그 정의를 임의 정밀도로 실현한다. `collect_planes` 의
   암시적-점 가지가 그 삼각형을 지으면서 **다섯 질문이 증발했다**: 인덱스 공간(표 항목이
   없으니 무의미)·`standard_for`(probe 가 표에 앉아 자동)·interning(16-1)·`tri_pt3` 전성(다시
   전체)·깊이(아래 잔여 하나). **판정층 새 기계 0, `Witness` 무변경, 술어 무변경.**
   - **실측**: 걸친-datum 프리즘의 Common **통과** — 한 몸·validate 클린·부피 = 스케치 진술
     그대로. 비용 **434 climbs**(wide 454·narrow 110, 고갈 0). 클래스 병합도 실행으로 증명:
     캡을 갈라 한 판정 surface 에 두 face → 한 클래스(같은 진술 → 같은 probe → 같은 사슬 →
     `shared_base` 의 증명된 0).
   - `ImplicitPlaneUnsupported` 는 태어난 지 한 단계 만에 도달 불가로 소멸(경계는 움직인다).
   - `plane_iv_through`(f64 필터 판)는 **소비자 없이 은퇴** — 존치 기준(은퇴 잔재 → 삭제,
     `ImprintSketch` 선례) 적용. 잠금은 이사했다: 16-2 의 Pure-대-Meet 차등이
     `plane_hp_through` 를 끝-대-끝으로 문다. 복원 경로는 dev-log.
   - ★ 정리 후보(비차단): 이름 없는 평면의 증인이 두 철자다(순수-혼합 = 해 점 / 걸침 = probe)
     — probe 가 실전에서 더 증명되면 한 철자로 접는 것을 재고, 접기 전에 비용을 잰다.

   **유일한 잔여 — 깊이**: 담체가 이름 없는 datum 인 정점(datum 위 datum 의 정점)은 producer
   가 `VerticesInMixedFrames` 로 거절한다. 인구 0; 인구가 생기면 재개하고 라벨 분리를 잰다.

17. ✔ **닫힘 (2026-08-09, `42ce674`·`c2ceb7f`·docs) — 폭은 진술의 문제였던 적이 없다.**
   원문: 담체가 프레임을 공유하고 만나는 점도 정확한 유리수인데 `Rat` 에만 안 들어가면
   producer 가 `VertexPointTooWide` 로 거절했다. 문서 처방대로 닫혔다:
   - **scalar**: `plane_name_from_meets([&MeetPoint;3])` — `plane_name_big` 의 lift-이후 몸통을
     분리(`plane_name_from_lifted`)해 폭 무관 점에서 같은 join·정준화. 잠금은 양쪽 다
     손-oracle: wide 만남 셋이 놓인 **T: 7x+11y−13z+1=0**(narrow 이름) 과 **T_w: 2⁷⁰x +
     y/5³⁰ − z + 1 = 0**(정준 정수 ~2¹⁴⁰ — wide 이름) — 유도가 자기 자신과 비교되지 않는다.
   - **topo/ops**: `through_meets`(폭 필터 없는 solve) 위에 `plane_name_through` 를 다시
     세우고(`through_points_rat` 는 전-narrow **투영**으로 유지 — 증인 base 가 `[Rat;3]` 인
     타입 사실), producer 의 거절 줄과 **`OpError::VertexPointTooWide` 변종을 삭제**
     (`ImplicitPlaneUnsupported` 선례 — 도달 불가 변종은 지운다). Wide 만남의 캐시 앵커는
     **첫 저장 정점의 세계 캐시**(정의의 자기 재생 — 8/8 잠금이 근거; wide 점의 f64 는 5⁻⁴⁰
     처럼 정확 lift 가 `i128` 을 떠나므로 재생-재현이 불가능하고, 그럴 필요도 없었다).
   - **구현이 잡은 둘째 벽(계획 1회차가 예견한 자리)**: «named ⇔ narrow 풀이» 불변식을 이
     항목이 깨므로 ① `surface_witness_triangle` 에 probe 폴백(named 인데 풀이-불가 →
     자기 프레임의 정준 probe 셋 — collect_planes 의 probe 가지와 같은 논증), ② 항목 15 의
     σ-접기(`name_stored_ints`)는 증인 base 가 이름 평면 위가 아니면 — probe 증인은
     프레임-국소 좌표라 **이름과 다른 프레임을 말한다** — debug_assert 가 아니라 **decline**
     이 옳다(`name_ints` None → 그 평면은 톨러런스 길 유지: 느리고, 틀리지 않는다).
   - **e2e 실행 증명**: wide-만남 datum 이 이름·interning(재진술 = 같은 핸들; 공유-담체
     삼중 = 놓인 평면의 기존 핸들로 intern)·프레임·스케치·**불리언**까지 통과 — named 평면의
     probe-증인 인구의 첫 실행. moved 팔(37° 공유 사슬) 포함. 판정층·collect_planes·
     frame_chain 은 예고대로 **0 변경**.
   - 코퍼스 인구는 여전히 0(census 비트 동일 ×3) — 능력이 열렸고, 픽스처는 합법적 수제
     모델(`push_plane`+`push_vertex`)로 상시 잠금이 됐다. `MeetPoint::coord_f64` 는 구현 중
     소비자가 죽어 **커밋 사이에 태어나고 은퇴했다**(미소비 기계는 남기지 않는다).

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

    실측: `nacre-tess/src/lib.rs:480` 이 `model.vertex_point(v)` 로 **f64 캐시**를 읽는다 —
    그 캐시는 자기 오차를 `PointCache.tol` 로 들고 있다(즉 **최근접 f64 가 아닐 수 있다**).
    ⇒ 화면에는 충분하지만, 내보내기(OBJ·STEP)는 **고정밀 실현 → 한 번만 반올림**을 원한다.

    ☑ **기계는 이미 있다**: `nacre-cip` 의 `HpApprox`(`add`/`sub`/`mul`/`div`/`inv_sqrt` 전부
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

    #### 설계 (2026-09-10 이어진 논의로 확정된 모양 — 아직 아무것도 안 지음)

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

18. ⏳★★★★ **필독 문서가 죽은 «타입 구조»를 현재형으로 가르친다** — **DNA 조각 ✔ (2026-09-14), 나머지 둘은 남음.**

    ⚠ **2026-09-11 의 진단표는 틀렸다**(2026-09-14 재측정이 반증). 그 표는 코드 열을 **0**으로 적었는데,
    계기가 *"코드 어디든 나오면 살아 있음"* 이라 **주석이 죽은 타입을 가려 주고** 있었다(`ImprintSketch` 가
    «코드 1곳»으로 통과했는데 그 1곳이 *"since `ImprintSketch` retired"* 라는 주석). **«선언 또는 주석 제거
    후 사용»으로 교정하면 design.md 가 이름 짓는 타입 321개 중 죽은 것이 17이 아니라 47개다.**
    ⚠ 그리고 `Constructed`/`Discovered` 는 **죽은 낱말이 아니다** — `EPS_CONSTRUCTED` 는 살아 있는 공개
    상수(사용 10곳)이고, 항목 23 이 그 낱말을 `PointCache::{Constructed, Discovered}` 변종으로 되살린다.
    **죽은 것은 «어휘»가 아니라 «논지»다.**

    ✔ **닫힌 조각 — DNA (2026-09-14).** `overview.md` 는 `CLAUDE.md` 가 *"절대 원칙, 위반 금지"* 로
    지정한 문서인데 원칙 4 본문이 죽은 enum 으로 원칙을 진술했다. **원칙(*"tolerance는 발견된 교차에만"*)은
    글자 그대로 살리고 «타입 시스템에 새긴다»는 기전 절만** 오늘(정점의 실현 캐시)로 정정하고 S7 의 반증을
    날짜와 함께 기록했다. `design.md §0` 다섯째도 **같은 문장으로**. 함께: `overview:48`(예고했던 relaxation
    3단 사다리는 짓지 않았다 — 답은 `nacre_ops::realize`) · `overview:51`(`SurfaceDef` 이름표).

    ⇒ **남은 둘**(각각 자기 칸):
    ① **`design.md §6` 의 칸 서사를 dev-log 로** — §10 이 *"design.md는 설계 규칙·불변만 유지"* 라 명령하는데
    §6 이 **1447줄(문서의 64%)·날짜 든 줄 112**로 그걸 어긴다. 이것이 문서가 2260줄이 된 구조적 원인이고,
    은퇴한 내부 이름 12개와 «낡은 주장» 43줄 중 26줄이 전부 거기 있어 **옮기면 함께 해소된다**.
    ② **`design.md §4`·§3 의 죽은 논지** — §4(273–376)는 헤더·리드·`Origin` 블록·`SurfaceDef` 블록·
    `Edge.curve`·곁표·«효과가 세 가지다»가 죽었고 `Face`·HalfEdge·Adjacency·rim 규칙은 산다. §3 은 relax
    블록(246–263) + 269. ⚠ **267 «정밀도 분업»은 글자 그대로 둔다** — 소스 포인터 13개가 인용하고 relax
    블록 밖이다. ⚠ §9 는 «현재 열»과 «그때의 근거»가 한 행에 섞였다(2071).

    ★ **규율**(이 칸이 세운 것): 스윕 결과는 «후보»지 작업 목록이 아니다 — 자리마다 **거짓 / 이미 정정된
    기록 / 날짜 붙은 기록 / 미래 표기**로 가른다. 그리고 **숫자가 아니라 «분류»가 잠금이다**(죽은 이름을
    정정하려면 그 이름을 불러야 하므로 스윕 수는 안 줄 수 있다).

    ⚠ **덤으로 드러난 것 — 두 필독 문서의 원칙 «목록»이 다르다.** overview 는 {정확기하, append-only,
    연산이력, **tolerance**, **진실/캐시 분리**}, design §0 은 {정확기하, append-only, 연산이력,
    **tessellation**, **tolerance**}. 셋만 겹친다. 그래서 `design §6 649` 의 *"(원칙 4)"* 는 오기가 아니라
    **overview 의 번호를 따른 것**이다. ⇒ 번호를 맞추려면 **어느 다섯이 DNA 인지부터** 정해야 한다 —
    별건(사용자 결정 사항).

19. ⚠★★★★ **`design.md` §3 의 `Curve::Intersection` 은 split 이전 설계의 잔재다** (2026-09-11, 기록만).

    실측: `git log -S "surfaces: [Handle<Surface>; 2]" -- docs/design.md` → **`2f84fae`
    (2026-07-07), 최초 설계 커밋, 그 뒤 한 번도 개정되지 않았다.** 진실/캐시 분리는 `4c51a5e`
    (2026-08-01), 곡면 진실은 `a7c9bd1`(2026-08-06). **그때 geom 의 `Surface` 가 진실이었고**
    (그 doc 이 스스로 *"the exact truth of a face's geometry"* 라 적고 있었다 — 이번 칸이 고쳤다),
    그래서 `Handle<Surface>` 는 「진실을 가리키는 핸들」이었다.

    ⇒ 아레나 반전 뒤 그 변종은 **영구히 불가능**하다: `Handle<T>` 는 아레나 항목을 이름 짓고, 항목은
    `Handle<MotionNode>` 를 든 topo 의 진실이므로 **topo 아래 크레이트는 그 타입을 이름 지을 수 없다**
    ⇒ geom→topo 순환. ⇒ geom 은 `nacre-store` 의존을 **영원히 갖지 않는다**.
    ☑ `Store`/`Handle` 을 최하위에 둔 **결론은 그대로 옳다** — 근거만 *"geom 도 쓴다"* 에서
    *"기하·위상을 모르는 순수 인프라이고 topo·ops·validate 가 쓴다"* 로 바뀐다(이번 칸이 `nacre-store`·
    `nacre-geom` 의 주석과 design.md 68행을 그렇게 고쳤다).

    ⇒ **그 변종의 세 조각은 오늘의 분리선 양쪽으로 갈라진다**:

    | design §3 이 적은 것 | 오늘 그것이 사는 곳 |
    |---|---|
    | *"절차적 정의(진실) = 두 곡면"* | **`Edge` 의 진실(담체)이 이미 기록한다** |
    | *"근사 스플라인 `NurbsCurve` — 진실이 아니라 캐시"* | `EdgeCache` |
    | *"`cache_err` — 캐시가 진짜 교차에서 벗어난 최대 거리"* | `EdgeCache`(§캐시의 수렴 주석) |

    ⏸ **코드 인구 0** 이므로 짓지 않는다(실측: geom 에 `Handle<` 0건, `Curve = Line | Circle`).
    M7 이 마친 교차를 실제로 만들 때, `EdgeCache` 를 `{curve, err}` 로 키우는 것이 그 자리다.

20. ⚠★★★★ **`let`-`else` 30곳이 새 곡면 종류를 «조용히» 놓친다** (2026-09-12 실측, 진단만).

    `m.surface(h)` 호출 **106**곳의 모양을 세었다:

    | 받아서 하는 일 | 수 | 구(Sphere)·원뿔이 추가되면 |
    |---|---|---|
    | `match` / `matches!` | **58** | ☑ **비망라 컴파일 오류** — 컴파일러가 전부 짚는다 |
    | `let`-`else` / `if let` 로 한 변종만 | **30** | ⚠ **조용히 `else` 로 떨어진다** |
    | enum 통째로 · enum 메서드 | 18 | 손봐야 함(캐시 변종화의 대가) |

    ★ 이미 이름이 있는 함정이다 — *"새 enum 변종은 `match` 엔 보이고 `let`-`else` 엔 안 보인다"*.
    ⚠ **캐시를 변종으로 바꿔도 이건 안 고쳐진다**(`let SurfaceCache::Plane{..} = … else` 도 똑같이
    조용하다) ⇒ **별개의 일**이다.

    ⇒ **30곳을 한 건씩 읽어 의도를 판정해야 한다**, 그리고 의도가 둘로 갈린다:

    | 그 줄의 의도 | 옳은 철자 |
    |---|---|
    | **「평면이어야만 한다」**(여기 다른 종류가 올 리 없다) | **타입 문** `plane(h) -> Option<&geom::Plane>` — 구가 와도 «평면 아님» 으로 정직하게 답한다 |
    | **「모든 곡면 종류를 다뤄야 한다」**(오늘 둘뿐이라 줄여 쓴 것) | **`match`** — 그래야 새 종류가 컴파일 오류로 온다 |

    ⏸ 오늘 인구는 **0**(곡면 종류가 둘뿐이라 아직 아무것도 안 놓친다) ⇒ 지금 짓지 않는다.
    **구·원뿔·NURBS 를 더하기 «전에» 반드시** — 그 칸의 첫 단계가 이 30곳 판정이다.

21. ⏳★★★★ **저장소 가시성을 «전부 비공개»로 — 방향은 확정, 막고 있는 것이 둘** (2026-09-12 실측).

    오늘 곡면·모션만 비공개이고 정점·간선·면·셸·솔리드는 **공개 필드**다. 같은 위험(캐시 없는
    raw push)을 **곡면은 타입으로, 나머지는 검사로** 막는 셈이고, 문서가 그것을 알고 적어 두기도
    했다 — *"A raw `edges.push` without a cache entry desyncs the two … the store stays `pub`"*.
    ⇒ 통일 방향은 **↑ 전부 비공개**다(↓ 공개로 내리면 칸 ㊷이 얻은 봉인을 잃는다).

    **밖에서 쓰는 연산 실측**(topo 밖):

    | 연산 | 수 | 문으로 대신 |
    |---|---|---|
    | `.get()` | **683** | ☑ §문의 이름의 `vertex`/`edge`/`face`/`shell`/`solid` — 기계적 |
    | `.len()` | 97 | ☑ `*_count()` (`surface_count` 선례) |
    | `.iter()` | 60 | ⚠ 아래 (b) |
    | `.push()` | **76** | ⚠ 아래 (a) |

    ☑ **정점·간선은 이미 깨끗하다** — 밖에서 raw push **0건**(둘 다 `push_vertex`/`push_edge` 를
    지난다). 처음 보인 4건은 `Tessellation` 자기 저장소였다(오탐).

    ⚠ **(a) 면·셸에는 문이 «없다».** 프로덕션 raw push **74곳**(faces 40 · shells 34 —
    `boolean.rs`·`ops.rs`·`transform.rs`): 불리언이 결과를 조립하며 아레나에 직접 민다.
    `push_face`/`push_shell` 이 아예 없어서다. ⇒ **이 항목의 본체는 「비공개로 바꾸기」가 아니라
    「면·셸의 문을 설계하기」**이고, 그 문이 지켜야 할 불변식(면↔셸의 규칙, 캐시 없음)을 먼저
    정해야 한다.

    ⚠ **(b) 전량 순회 60곳의 정당성이 미검증이다.** 곡면 저장소는 *"전량 순회 문이 의도적으로
    없다 — 아레나엔 superseded 도 있으니 소비자는 live face 를 걷는다"* 인데, 정점·간선·면은
    밖에서 전량 순회가 된다. 그 60곳이 **superseded 를 같이 읽고 있는지** 재야 한다 —
    읽고 있다면 그것은 가시성 문제가 아니라 **오늘의 결함**이다.

    ⇒ 순서: (b)를 먼저 재고(결함이면 그것부터), (a)를 설계하고, 그다음 필드를 닫는다.

    ★ **같은 칸에 얹는 것 — `VertexDef` 를 `Vertex` 로 접는다** (2026-09-13 확정). `Vertex { def }` 는
    S7 이 `point`·`origin` 을 캐시로 보내고 남은 **한 필드 껍데기**다(`f86e7e3` → `76a07b2`). `Surface` 는
    진실이 enum 그 자체인데 정점만 struct 한 겹을 더 벗겨야 한다 — `m.vertex(v).def` 가 아니라
    `m.vertex(v)` 를 바로 `match`. 실측: `VertexDef` 173 · `.def` 218 자리, 전부 기계적이고 타입이 바뀌므로
    놓친 자리는 컴파일 오류, census 무영향. 문의 이름 칸이 `vertex_point` 111곳으로 정점 자리를 어차피
    지나므로 **같은 손으로** 한다 — 따로 하면 같은 파일을 두 번 연다.

    ★ **같은 칸에 얹는 것 둘째 — `Rotation.point`·`Motion::Rotate.point` → `pivot`** (2026-09-13 확정).
    `Axis` 는 방향(X/Y/Z)뿐이라 회전엔 «어느 선»인지가 더 필요하고, 그 필드가 `point` 다 — doc 은 이미
    *"the rational **pivot** `point`"* 라 6번 부르는데 이름만 뜻을 안 말한다. 접근 30곳(테스트 제외), 순수 개명.

23. ⏳★★★★ **`Bounded` 통합 — 「값 + 경계」에 이름이 둘이고, 산술은 한쪽에만 있다** (2026-09-12 실측).

    | 오늘 | 어디 | 내용 |
    |---|---|---|
    | `pub type Bounded = (BigFloat, Mag)` | scalar (칸 ㊵가 만듦) | 튜플, 산술 없음 |
    | `pub(crate) struct HpApprox { value: BigFloat, error: Mag }` | cip `interval.rs` | **같은 내용**, 메서드 **20개** |
    | `pub(crate) struct Approx { value: f64, error: f64 }` | cip | f64 판 |

    한 개념에 이름이 둘이고, ㊵ 감사 네 라운드가 못 봤다. ⇒ **`Bounded { value: f64, error: f64 }` ·
    `HpBounded { value: BigFloat, error: Mag }`** 로 통일하고 튜플 별칭은 죽인다(문서 규칙 그대로 — 같은 것의
    다른 정밀도는 `Hp` 접두사 하나로만 다르다). 「근사(`Approx`)」가 아니라 「참값이 반경 안에 있음이
    보장됨」이 이 타입의 정체이므로 이름도 그것을 말해야 한다.

    ⚠★★★ **개명이 아니라 «이사»다.** 산술(`add`·`sub`·`mul`·`div`·`div_exact`·`inv_sqrt`·`sign` — 20개)과
    변환 헬퍼(`rat_to_hp`·`bigint_to_hp`·`ub`·`lb`)가 cip 에 살고, Rust 는 남의 타입에 inherent impl 을 못
    다니 **타입과 산술이 함께 scalar 로 내려간다.** 크레이트 규칙(«유도된 값 + 그 산술 = scalar»)과 맞고,
    ㊵가 튜플을 새로 만들어야 했던 것 자체가 그 증거다(cip 의 struct 가 `pub(crate)` 라 scalar 가 못 썼다).
    개명 자리(`Approx` 50·`HpApprox` 115)는 그 뒤의 기계적 일. census 무영향.

    ★★ **같은 수술로 판정 쪽 «cache» 어휘를 없앤다** — `WorkingPlane.cache` → `realized`,
    `WitnessPoint { coord, tol }` → `[Bounded; 3]` 꼴. §캐시가 갈라 적은 «f64 가 둘» 오해가 **이름 수준에서**
    불가능해진다. `HpPointCache` 도 이 때 다시 본다: 모델 캐시가 아니라 연산 하나짜리 실현이므로 «Cache»
    가 붙을 이름이 아니다 — `[HpBounded; 3]` 이거나 `WitnessPoint` 자신이다.

    ★★★★ **원자만이 아니라 «계단 2»도 통일된다 — «경계 지어진 점»** (2026-09-13, §캐시 수치 층 그림).
    「값+오차」가 오늘 **세 곳에 값 배열·오차 배열을 «따로»** 든다: `realize` 출력 `([f64;3],[Mag;3])` ·
    `WitnessPoint{coord:[f64;3], tol:[f64;3]}` · `PointCache{coord, tol}`. ★ `realize` 는 **내부적으로 이미
    `[Bounded;3]` 로 묶어** 든다(`Approached([Bounded;3], _)`) — «묶은 점» 타입이 한 자리에 있고 나머지가
    안 따라간 것뿐이다. 묶으면 정리이자 **건전성**이다: «참값 ∈ value±error» 가 구조로 서서 transform 이
    값만 옮기고 오차를 안 옮기는 desync 가 불가능해진다.
    ⚠ **그러나 `WitnessPoint` 를 `PointCache` 로 접지는 않는다.** `WitnessPoint = [Bounded;3] + 정의
    (`base:[Rat;3]`·chain) + hp(메모)` 로, 그 «정의» 가 정확 단계의 입력이다 — 캐시는 일부러 안 든다. 계단
    2 를 **공유**하되(둘 다 `[Bounded;3]` 를 품는다) 계단 3 은 상위집합으로 남는다. 접으면 진실/캐시 경계가
    지워진다.

24. ✔ **`Branch` → `Pierce` — 이름이 «무슨 점»인지 말한다** (2026-09-13 완료).

    최초판 Q5(2026-08-01)의 임시어 «가지 번호»가 변종 이름에 남아 있었다 — `root: QuadRoot{Lo,Hi,Double}` 이
    «어느 근»을 이미 받는데도. 이 점은 도법기하의 **관통점(piercing point)**: 두 평면의 교선이 원통을 뚫는 점.

    개명한 것 — **① 관통점 엔티티**(topo+ops+validate): `VertexDef::Pierce`, 형제 enum 넷
    (`NodeId`·`PointOn`·`OnLine`·`Located`)의 `Pierce` 변종, 헬퍼 `pierce_*`(`pierce_meet`·`pierce_name`·
    `pierce_vertex_tol`…), 산문·테스트 라벨. 그대로 둔 것 — **② 근**(scalar/cip 의 `*_branch`,
    `QuadRoot` 자체)과 **③ 그래프/제어**(`SketchError::BranchingVertex`, "code branch" 산문).
    ★ **구조가 경계를 강제했다**: scalar/cip 는 topo 에 의존하지 않아 그 `branch` 는 관통점일 수 없다(어느 근).
    그래서 개명 전 분석표가 `cmp_coord_branch` 를 ①로 넣은 것은 틀렸다 — 근 뜻이라 그대로다.
    순수 개명 · **census 비트 동일**(debug·release 양 프로파일, 398행).

25. ⏳★★★★ **원/원호의 진실은 «실현값»이 아니라 «정의»여야 한다 — 평면의 normal-vs-coefficients 와 같은 갈래**
    (2026-09-13 진단, 곡선 마일스톤 입력).

    오늘 원호는 «실현된 값»으로 저장된다: `Seg2d::Arc{center: [Rat;2], radius: Rat, ccw}` ·
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
    `radius`→`r²` 를 `Seg2d::Arc`·`CylinderDef` 에서 함께 + 중심의 정의 기반 표현 신설) — 지금 고치는 게
    아니라 그 계획의 입력이다. ★ `arc_turns`/`arc_rat` 이 갈린 것은 이것과 무관한 편의 설탕(끝점을 90°k 로
    유도 vs 받음)이고, 정의 기반이 되면 둘 다 «정의를 진술하는 한 방법»으로 자연히 정리된다.

26. ⏳★★★ **입력을 «순서 있는 고리»로 통일한다 — `Edge2d`/`Seg2d` 가 하나가 되고 `from_edges` 의
    순서 재발견이 사라진다** (2026-09-13 진단).

    오늘 스케치 조각은 **둘**이다: 입력 `Edge2d::Arc{center, radius, start, end, ccw}`(양끝을 자기가 든다) ·
    저장 `Seg2d::Arc{center, radius, ccw}`(시작=앞 꼭짓점, 끝=다음 꼭짓점). 차이는 **«양끝을 자기가 드느냐»**
    하나뿐이고, 그것이 곧 「순서 없음 ↔ 순서 있음」의 자국이다. `from_edges` 는 `Vec<Edge2d>` 를 받아 끝점
    탐색으로 **사슬을 재발견**한다 — 뭉텅이로 흩어져 들어와도 잇기 위한 층이다.

    ★ **코드캐드라 그 층이 필요 없다.** kit 의 `PenPath` 는 이미 «펜»이다(`SketchSeg::Arc` = *"an arc from
    where the pen stands"*, `close()` 가 시작점으로 되돌려 닫는다) — 순서를 **이미 안다**. 그런데 그걸
    `Vec<Edge2d>` 로 펼쳐 순서를 버리고, 커널이 `from_edges` 로 도로 찾는다. 「펜 → 조각 목록 → 순서
    재발견」의 가운데가 헛돌음이다.

    ⇒ **셋을 한다**: (a) 「순서 있는 + 호 되는」 입력 문 — `from_rings(Vec<Vec<Point2>>)` 가 절반이나 이미
    그 문인데(순서 있는 고리) **`Ring2d::polygon`, 즉 직선만** 받는다; 호 버전이 비어 있다. (b) `from_edges`
    의 순서 탐색을 지우고 그 문으로 대체(②검증 `OpenChain`·`BranchingVertex` 와 ③`classify` 겹침·구멍
    판정은 **입력 형태와 무관**하므로 그대로 산다). (c) `Seg2d` 를 `Edge2d` 로 흡수 — `start`/`end` 필드가
    소멸하고 조각 타입이 하나가 된다. ★ 살아남는 이름은 **`Edge2d`** 다: 검증 문(`arc_rat`→`radius_of`,
    `arc_turns`, `circle`)이 거기 붙어 있고 `Seg2d` 는 geom 의 검증 없는 데이터다(실측: `Edge2d::` 58곳).

    ☑ **열린 항목 25 를 쉽게 만든다**: radius→r² 를 오늘은 `Edge2d`·`Seg2d`·`CylinderDef` **셋**에서
    해야 하는데, (c) 뒤엔 **둘**이다. ⚠ 실제 코드 변경이라 자기 칸이 필요하고, kit 쪽은 «펜을 `Edge2d`
    목록으로 펼치는 단계»를 없애 펜을 커널의 고리 문에 직접 넘기는 것이 자연스럽다.
