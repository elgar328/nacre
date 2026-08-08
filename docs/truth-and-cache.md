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
   > 강하다 — 정점은 `VertexDef::ThreePlane` 으로 **좌표를 아예 안 들고**(세 핸들뿐), 모션은
   > 규칙 3 대로 **면이 든다**.
   >
   > ★★★ **그래서 «연산이 누적되면 결국 넘치지 않나» 의 답은 «아니오» 이고, 논증이 아니라
   > 실측이다**(`tests/point_width.rs`): 40회를 쌓아도 이름 폭이 **1·4·2비트로 상수**이고, 그
   > 사이 surface 는 246·484 개로 늘지만 **distinct 이름은 6(10)개 그대로** — 움직인 면은 전부
   > **같은 이름 + 모션 노드**다. 셋이 맞물린다: 불리언은 `push_plane` 을 **한 번도 부르지
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
   매번 답한다(느릴 뿐, 틀리지 않는다).
7. **`Rat` 을 넓히지 않는다(C8).** 넓어지는 것은 중간값(BigInt 이름 유도·`Expansion` 술어·
   BigFloat 상승)뿐이다.

---

## 최종 타입 — 진실

```rust
// ─── nacre-geom / nacre-topo ── 진실 (아레나, append-only) ─────────────

pub enum Surface {
    Plane {
        points: PlanePoints,
        motion: Option<Handle<MotionNode>>,   // ← SurfaceDef 를 흡수한다
    },
    // Cylinder { point: [Rat;3], dir: [Rat;3], radius: Rat, ref_dir: [Rat;3],
    //            motion: Option<Handle<MotionNode>> }        — M6, 자기 진실로
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

/// 정점 = 자기 **정의**. 좌표는 진실이 아니라 캐시다(PointCache). `Origin`·`point` 는 소멸.
/// ★ 단일형이 아니라 **두 변종**이다(S7 에서 Q3 수정 — 아래 반증표): M5 에 원통 seam 정점이
/// 실재하고 그 점은 «세 평면»으로 적을 수 없다. 변종별 불변식(Q5)이 이 구조의 근거다.
pub struct Vertex {
    pub def: VertexDef,
}

pub enum VertexDef {
    /// 세 평면의 교점 — 이름이 곧 점. (D != 0 은 좌표 재생이 생기는 자리에서 단언한다.)
    ThreePlane([Handle<Surface>; 3]),
    /// 두 곡면의 교차 «곡선» 위의 점 — M3 원통 seam(테두리 원의 θ=0). 점을 못 박는 매개
    /// 정보(원통의 `ref_dir`)는 M6 의 원통 진실과 함께 오고, 그때까지 **좌표 캐시가
    /// load-bearing** 이다(정직 기록). M6 는 변종을 더한다: `Branch{surfaces, branch}`
    /// (이차곡면 셋의 최대 8점), 원뿔 꼭짓점 `Apex(Handle)` 등.
    OnSeam([Handle<Surface>; 2]),
}

/// 곡면 집합은 「담체」를 정하고, 경계가 나머지를 정한다.
/// 정점(0차원)만 곡면 집합으로 완전히 정해진다 — 그래서 정점만 좌표를 진실에서 버린다.
pub struct Edge {
    pub surfaces: [Handle<Surface>; 2],       // 담체 (두 끝점 면집합의 교집합으로 파생 불가)
    pub vertices: [Handle<Vertex>; 2],        // 경계
}
// Store<Curve> 는 없다 — 곡선은 진실이 아니라 EdgeCache 다.

// ─── 모션 (nacre-topo — Handle 이 필요하다) ────────────────────────────

pub enum Motion {
    Rotate    { axis: Axis, point: [Rat; 3], angle: Angle },
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
  평면의 두 변종은 판정층 거울(`WorkingPlaneDef`)의 `match` 하나로 흡수된다.
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

---

## 스케치

```rust
// ─── 프로파일 — 프레임을 모르는 순수 2D. 구성 시점에 Rat 으로 확정된다. ──

pub struct Profile2d {                        // 한 재료 영역 — 현행 의미 그대로
    outer: Ring2d,
    holes: Vec<Ring2d>,
}
pub struct Ring2d { points: Vec<[Rat; 2]> }   // 지금은 직선 변만 — 호는 M6 에서 세그먼트
                                              // enum(LineTo/ArcTo)으로 확장할 자리만 남긴다
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
PlanePoints (진실 — 점 셋)
   │ 유도 (정준화: 분모 털기 → gcd → 부호 규약)
   ▼
PlaneName = Narrow([Rat;4]) | Wide([BigInt;4])    ← 저장은 이것 하나
   ├─ narrow() → Option<&[Rat;4]>    산술·프레임·지름길 — 사본이 아니라 빌려 읽는다
   └─ + 모션  → SurfaceKey           interning 표(`surface_ids`)의 키
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

/// interning 표의 **키** — 이름 + 그 이름이 진술된 모션. 같은 계수라도 Constructed(세계)와
/// Moved(모션 전 프레임)는 다른 평면이다 ⇒ 구조가 같을 때만 합친다.
pub type SurfaceKey = (PlaneName, Option<Handle<MotionNode>>);
```

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
    // 진실 — append-only
    surfaces: Store<Surface>,                 // ★ 비공개(S1 ✔) — push 는 생성자 경유만, 읽기는
                                              //   좁은 접근자(surface/surface_count/motion —
                                              //   전량 순회 없음: 아레나엔 superseded 도 있다)
    pub vertices: Store<Vertex>,
    pub edges:    Store<Edge>,
    pub motions:  Store<MotionNode>,          // interned
    // faces / shells / solids / live_solids …

    // 캐시 — 핸들 인덱스 병렬, 통째로 버리고 재생 가능. live 도달분만 lazy 채움.
    pub vertex_cache:  Vec<PointCache>,
    pub surface_cache: Vec<SurfaceCache>,
    pub edge_cache:    Vec<EdgeCache>,        // 평가 가능한 곡선 (M5 는 직선)

    // 이름 — 진실에서 유도해 하나만 저장(Narrow|Wide — §이름과 interning).
    // 곁표인 이유: Known 은 항상 있지만 무리수 모션의 Through 는 없을 수 있다.
    pub surface_name: HashMap<Handle<Surface>, PlaneName>,
    pub surface_ids:  HashMap<SurfaceKey, Handle<Surface>>,
}

// 수치 층 — 같은 것의 다른 정밀도는 Hp 접두사 하나로만 다르다. 값+tol 이 한 덩어리.
pub struct Mag { m: f64, e: i64 }                                  // f64 밖 범위의 보수적 크기
pub struct Approx        { value: f64,      error: f64 }           // 중간 스칼라
pub struct HpApprox      { value: BigFloat, error: Mag }
pub struct PointCache    { coord: Point3,        tol: Option<f64> }  // S7 실형: 측정치 유무
pub struct HpPointCache  { coord: [BigFloat; 3], tol: [Mag; 3] }
pub struct SurfaceCache  { coeffs: [f64; 4], tol: [f64; 4], inv_norm: f64 }
pub struct HpSurfaceCache{ coeffs: [BigFloat; 4], tol: [Mag; 4] }
pub struct EdgeCache     { curve: Curve }                          // 평가 가능한 담체 곡선
```

| 캐시 | 키 | 수명 | |
|---|---|---|---|
| 모델 캐시 (f64) | 핸들 인덱스 (밀집) | 모델과 같이 | live 도달분만 lazy — `trial_bound` 는 일부러 우회 |
| 고정밀 메모 | (정의, 정밀도) (희소) | **연산 하나** | 판정이 실제로 만든 점에만 |
| 각도·√ 표 | `(Angle, prec)` / `(Rat, prec)` | 프로세스(스레드) | 값이 키의 순수 함수 — `cos 37°`·`1/√(n·n)` 은 어디서나 같다. 정확 반올림이라 실현이 유일하고 tol 도 유일 |

셋은 인덱스 공간이 달라 합치지 않는다. `Through` 평면의 `SurfaceCache` 는 세 정점의 실현에서
유도된다 — 법선 = 점 차의 외적이라 점의 오차 기계를 재사용하고, 이름 실현의 tol 은 지름길에만
쓰이므로 틀려도 조용한 오답이 아니라 느려질 뿐이다.

---

## 판정 (연산 동안만 산다 — nacre-cip)

이름 규칙 셋: ① **`Working*` = 모델 타입의 판정층 쌍둥이이자 표의 뿌리**(수명 = 연산 하나 —
`Model` 에 `Working…` 이 담기면 읽는 즉시 이상해 보여야 한다). ★ 접두사는 **수명 표식이
아니라 뿌리 표식**이다 — 부품(`WitnessPoint`·`MoveNode`·캐시들)은 `Working*` 컨테이너 안에
살며 수명을 상속하므로 접두사를 반복하지 않고, 이름은 역할을 말한다(캐시 타입은 모델 쪽과
일부러 공유된다 — 반올림 사본임이 타입에 보이도록). ② **거울의 변종 이름은 진실과 같다**
(`Known`/`Through`). ③ 평면 정의를 이루는 유리수 점은 **증인 점 `WitnessPoint`** 다.
오늘 코드와의 대응(이행 중 개명): `Pt3` → `WitnessPoint`, `WorkingPoint` → `WorkingVertex`,
`PlaneGeom`/`FaceInfo` → `WorkingPlane`, `tri_pt3` → `def`, `HpIv` → `HpApprox`
(스칼라 구간 — 점 실현 묶음은 `HpPointCache`).

```rust
/// 판정용 평면 — `Surface::Plane` 의 쌍둥이. 정의를 펼쳐 들고 두 실현을 메모한다.
pub struct WorkingPlane {
    pub def: WorkingPlaneDef,                 // 정의 — 진실의 두 변종을 그대로 비춘다
    pub name: Option<PlaneName>,              // 이름 — Narrow 는 Shewchuk, Wide 는 Expansion
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

/// `PlanePoints` 의 쌍둥이 — 변종 이름까지 1:1.
pub enum WorkingPlaneDef {
    /// 증인 삼각형 — 유리수 base + 사슬. 오늘의 `tri_pt3` 그대로.
    Known([WitnessPoint; 3]),
    /// 정하는 정점 셋 — 각 정점을 **평면 표의 인덱스 셋**으로 지목한다(= `WorkingVertex`
    /// 의 저장 철자 — 표 안에 저장되므로 자기참조를 피해 인덱스로 적는다. 기계는 한 벌).
    /// ★ 유리수 닫힘 datum 의 판정은 `name` 의 **Wide 정확 경로가 먼저** 받으므로,
    /// 이 동차 상승은 무리수 모션 datum 의 **최종 심급**이다(§열린 항목 1).
    /// ★★★ **차수는 ~12 가 아니라 9 다**(확정 2026-08-08, S5(ii)-2a). 여기 적었던 ~12 는
    /// 어림이었고, 실제로는 «어떤 식을 쓰는가» 가 정한다: 분모를 털고 아핀 외적을 하면 **15**,
    /// 동차점 셋의 **사영 join**(3×4 행렬의 네 3×3 소행렬식)이면 **9**. 후자를 쓴다 —
    /// `cramer_iv`(세 평면 → 점, meet)의 **정확한 쌍대**라 `det3_iv` 네 번이고 새 산술이 없다.
    /// 기계는 `nacre-cip::kernel::frame3::{plane_iv_through, plane_hp_through}` 에 **서 있고
    /// 아직 호출자가 없다**(생산자는 2b).
    Through([[usize; 3]; 3]),
}

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
  ★★★★ **«유리수 닫힘이면 `Wide` 이름이 곧 정확 계수라 공짜» 는 오늘 참이 아니다 — 미래형이다**
  (확인 2026-08-08). 술어는 계수를 `exact_coeffs()`·`base_coeffs()` 로 묻고 **둘 다
  `narrow()` 를 읽는다** ⇒ `Wide` 이름이면 둘 다 `None` 이고 **정확 지름길이 전부 declined**,
  판정은 톨러런스+상승 경로로 간다. `PlaneName::Wide` 의 BigInt 계수를 실제로 소비하는 곳은
  `nacre-cip` 의 **`WideFrame`**(프레임 실현) 하나뿐이다 — topo 주석의 *"`Wide` carries
  identity only"* 그대로. ⇒ **술어에게 `Wide` 계수를 가르치는 일이 아직 안 됐고**, 그것이
  아래 열린 항목 15 다. 실측 비용: **상승 278 → 887**(narrow 이름 대 wide 이름, 3.2배;
  둘 다 성공하고 예산 고갈 0 — 비용이지 절벽이 아니다).
  **무리수 모션이 낀 datum 만** 실현 상승 전용이다: 정점을 자기 세 평면의 고정밀 실현에서
  동차좌표로 만들고 그 위에서 계수를 구간으로 유도한다(**차수 9** — 위 참조) — 정확성 위험이
  아니라 비용 위험이었고, ★ **이제 기계가 있어 쟀다**(S5(ii)-2a, `frame3.rs` 단위 테스트):
  - **깊이 1 은 필터가 산다.** 생성 200 사례 전부 결정, 계수 800개 중 미결 **0**, 최악 상대
    반경 **6.3e-10**. 차수 9 가 `Iv` 의 여유를 먹지 않는다 — 재기 전에는 몰랐던 것이고,
    이것이 이 단계의 진짜 관문이었다.
  - ★★★ **깊이 2 는 필터가 없다.** 담체가 또 `Through` 면 차수가 **81** 이 되고, 계수가
    `f64` 범위를 **8/8 전부** 벗어난다(고정밀 쪽은 멀쩡하다). 함수는 그때 `None` 을 돌려
    **상승으로 보낸다** — `NaN` 반경이 우연히 `sign()=None` 을 내는 것에 기대지 않는다.
    ⇒ 깊이 2 는 «느린 길» 이 아니라 **상승 전용**이다. 깊이 제한을 두어야 하는지가 2b 의
    질문이 된다(C6 이 유계를, C5 가 비순환을 이미 준다 — 남은 것은 비용뿐이다).
  - ★★ **배율의 부호는 값 안에서 없앤다.** join 은 행에 대해 다중선형이라 결과가 참 평면의
    `D0·D1·D2` 배이고, 음수면 **평면 방향이 뒤집힌다**(`frame_sign`·바깥 법선·라벨 프레임이
    전부 그 위에 있다). `Judge::plane_iv(k) -> [Iv;4]` 에 부호를 실을 자리가 없으므로 —
    실을 곳 없는 값은 아무도 안 쓰는 값이다 — 함수가 스스로 정규화한다. 생성 200 중 **172**가
    음수 `D` 를 지나므로 그 이빨은 실제로 물렸고, 정규화를 지우면 대조 테스트가 깨진다(확인).
  - `D` 가 0 을 품으면 `None` → 상승 → 안 갈라지면 이름 붙은 거절. 조용한 폴백은 없다.

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
| 3b | 좌표 재생 — 재보고 보류: 공유 프레임 풀이 + 모션 재생 = 비트 동일(8/8)로 목적지와 `Pt3` 재생의 무충돌만 확인해 둠 | ⏸ |
| S1 | `Store<Surface>`·`Store<MotionNode>` 봉인 — 좁은 접근자(`surface`/`surface_count`/`motion`) + `compile_fail,E0616` 잠금. 기록 없는 surface 는 크레이트 밖에서 표현 불가(테스트 전용 입구 `push_surface_unrecorded` 만 예외, `test-util` 게이트) | ✔ 2026-08-05 |
| S2 | **임의정밀 이름** — `PlaneName{Narrow\|Wide}` 그릇(`plane_name_big` 꼬리의 `to_i128` 분기 하나), `surface_coeffs`→`surface_name` 개명, Wide 도 intern. **동일성 팔 절단**: 점 가진 평면의 이름 실패가 공선뿐이 됨 — ★ f64 폴백 연쇄의 **프레임 팔은 S4 몫**(Wide 는 `narrow()=None` 이라 프레임·지름길을 안 연다, 반증표 그대로). 잠금: 스칼라 2(`a_plane_too_wide…`·narrow/wide 합의) + topo 1(`a_wide_plane_interns_but_opens_no_shortcut`) | ✔ 2026-08-05 |
| S4 | **프레임 팔 절단** — `Motion::Frame{plane, placement, flip}` + `FramePlacement{Canonical\|Named}`. `Canonical`(기본)은 사슬 펼침 시 유도: 좁으면 오늘의 `plane_frame_default` 그대로(비트 보존 — 유도의 이사), **넘치면 wide 도로**(`MoveNode::FrameWide` = `PlaneFrame` 의 BigInt 쌍둥이, 실현은 `HpIv` 전 구간 + f64 캐시는 128비트 실현의 좁힘). Wide 이름·`n·n` 넘침(1.6%) 인구의 프레임이 열려 §12 연쇄의 프레임 팔이 닫힘. 잠금: cip 2(on-plane·tol-bounds wide 쌍둥이) + ops 3(`a_wide_plane_hosts_a_canonical_frame`·nn-넘침 개통·**종단** `a_pad_on_a_wall_with_overflowing_squares_takes_the_exact_road`) + census `wf` 가족. ★ 계획의 "심기 + 공개 SketchFrame 통일"은 **S9 로 분리**(아래) | ✔ 2026-08-05 |
| S3 | **`Profile2d` 유리수화 + 공선 중간점 정리** — `Profile2d{outer: Ring2d, holes}`·`Ring2d{points: Vec<[Rat;2]>}`, 생성자가 리프트(창 밖 = `ProfileOutsideDecimalWindow` 구성 시점 에러) + 공선 중간점 소멸(엄격 내부만 — 중복점·스파이크는 생존해 제 이름으로 보고). `check()`·`from_rings` 분류가 진실 위 정확 술어로(`orient2d_rat` scalar + geom `_rat` 쌍둥이 — geom 이 scalar 의존 획득, §design 1 격리 규칙 준수). 잠금: 쌍둥이 프리즘 비트 동일+전 코너 3-평면 정의(Q2 ② 닫힘), 창 에러 2, 십진-이진 부호 분기(십진이 이긴다), 비인접 공선 벽 interning 재핀. census 150줄 비트 동일. ★ 실측: check 는 34ms@100점(호출당 ~1.7µs, gcd 지배 — §열린 항목 7) | ✔ 2026-08-05 |
| S6a | **점 없는 평면의 소멸** — `Inexact` 소멸(S6b 타입 교체)의 전제. `PlaneDef` 를 점 셋 단일 필드로(origin=points[0]·ref_dir=points[1]−points[0]·극성=점 순서 — 불변식이 구조, 좁은-계수 def 실패 계급 사망, `named_plane_points` 은퇴), **`from_axes` 가 축의 십진 진실을 정의로**(#28 «축만 든 호출자» 인구 개통 — 45°급 프레임의 프리즘이 `WideFrame::named_of` 로 정확 경로; 내부의 실현-기저 호출 3곳은 의도적 무-def `realized_plane` 분리 — 두-정확-기술 재발 방지), `add_cylinder` 캡 점 기록(add_cuboid 선례), 넘침-이동은 노드 기록(`motion_is_exact` 프로브에 점 수송 포함). 잠금: Named×Wide 프레임 단위 + 축-전용 기울어진 프리즘 종단 + 원통 캡 interning + 넘침-이동 + **전수 관문**(`points_coverage` — 생산 경로별 모델의 live 평면 face 전수가 점 보유). census 150줄 비트 동일 ×3회 | ✔ 2026-08-05 |
| S6b | **타입 교체** — 진실 스토어(`SurfaceTruth{Plane{points: PlanePoints::Known, motion} \| Cylinder{motion}}`)가 캐시 store 와 인덱스-평행으로 탄생, `SurfaceDef`·`surface_defs`·`surface_points`·`push_surface(_with_points/_unrecorded)`·`Violation::UndefinedSurface`·`RejectReason::{InexactSurface, CoordinateOutOfRange}` **사망**. push 는 `push_plane`(interning, flipped 는 f64 캐시 내적 그대로 — 같은 평면이라 부호 정확)·`push_cylinder(motion)`·test-util `push_plane_unregistered`/`set_plane_points_for_test`. f64 프리즘 폴백 → **이름 붙은 거절**(`PlaneWithoutExactForm`·`DistOutsideDecimalWindow` — `Swept::along` 삭제, build_prism 정확-전용). 부수 개선: 이동된 원통이 `Inexact` 강등 대신 모션 기록. ★ 구현 중 반박 1건: normal_def 의 `v = n×u` 곱이 작은-지수 전폭 법선(분모 10²¹→10⁴²)에서 넘침 — proptest 가 폴백 소멸 당일 발견, 원시 방향조차 137비트라 **기저-교차 셔플**(`w = x̂×n` + 대수 부호 `det[ẑ,x̂,n]=n₁`)로 재구성(곱 0개, 전역). 아레나 반전(캐시가 Store·진실이 Vec — `Handle<T>` 타입 매개변수가 강제)은 최종 개명 시 제자리로(§열린 항목). census 150줄 비트 동일 ×4 | ✔ 2026-08-06 |
| S9 | **공개 스케치 API 통일 + world 평면 사전 심기** — ① `Model::new()` 가 세계 축 평면 셋을 심는다(핸들 0·1·2 = XY·YZ·ZX, points 는 `axis_plane` 삼중 `[0,u,v]`, **캐시 방향은 −축** — extrude 밑캡의 감각과 일치, +축이면 실측 781 캡 flip 재도입; `#[derive(Default)]` 제거 = 무씨앗 뒷문 폐쇄, `world_plane(Axis)` 접근자, `stat seeded_hits` 반증성 다리 신설 = 실측 455). census ε-재기준 1회: 평면 digest 이동 127/150줄, **결과는 143/150 비트 동일 + 나머지 7줄도 부피·면적·centroid 전부 비트 동일**(정점 해시만 이동 — Cramer 가 사실상 스케일-불변으로 반올림, 스칼라 최대 편차 정확히 0), ERR/EMPTY·피연산자 정점 해시 문자 동일. ② 공개 `SketchFrame{plane, placement, flip}`(필드 private + 검증 생성자 — 리터럴 우회 봉쇄): `named()` 가 구성 시점 거절 `FrameOutsideDecimalWindow`·`OriginNotOnPlane`(신규 scalar `plane_residual_sign` — orient2d_rat 급 **전역**, Wide 는 BigInt 팔)·`RefDirParallelToNormal`(판정은 `WideFrame::named_of` 재사용 — 폭에 전역이라 None = 평행뿐), 이름 없는 평면 = `PlaneWithoutExactForm` 재사용. `face_sketch_frame` 신설(이음새 — face_frame 이 만들던 값을 버리지 않고 공개). ③ 내부 통일: flip 측정은 `measured_frame` 한 곳, 노드 push 는 `push_frame_node` 한 곳(extrude·face 두 도로가 한 모양, 게이트 표현식 문자 유지, census 비트 동일). ★ **`Operation` 의 평면-핸들 어휘 교체는 S5 로 유예** — replay 자기완결성: 로그 속 핸들의 합법 표적은 씨앗·기존 면·datum 뿐인데 datum op 가 S5 에야 생긴다. ★ 잠금서 확정 둘: 씨앗 intern 직접 증거(원점 상자 바닥/왼쪽/앞 + z=0 밑캡 = 씨앗 핸들, 아레나 6 유지), ZX 의 canonical 프레임은 `−x̂`(스크립트 삼중과 다름 — Named 로 말할 사례임을 잠금이 명문화) | ✔ 2026-08-06 |
| S8 | **Edge 최종형** — `Edge{surfaces: [Handle<Surface>;2], vertices: [Handle<Vertex>;2]}`: 담체 두 면(오름차순 정렬 쌍) + 경계 두 점, `curve`·`bounds: Option`·`origin` 사망. `Store<Curve>` → `edge_cache: Vec<EdgeCache>`(인덱스-평행 캐시): 유일 입구 `push_edge`(eager 파생, 퇴화 검사는 **팔별** — rim `[v,v]` 는 합법) + `rebuild_edge_cache`(«버리고 재생» 잠금이 비트 동일 증명) + `edge_curve` 접근자·`derive_edge_curve`(직선 = 끝점 through_points, rim 원 = 담체에서 — 신설 geom `line_plane`, seam = 자기-인접 `[cyl,cyl]` 잠정 표기). transform pass 2(곡선 이동) 통째 소멸. validate: 신설 `EdgeCarrierMismatch`(담체 ≠ 인접 관측, `[plane,plane]` 자기쌍 검출) + `UnboundedEdgeInLoop`·`RefKind::EdgeCurve`·`StepError::UnboundedEdge` 순삭. ★ 구현 중 발견 2건: ① **담체는 wall 로 추측하면 틀린다** — 세 평면이 한 직선을 공유하는 인구(해결된 4-평면 동시성)에서 각 면의 arrangement 는 제3의 평면을 wall 로 (옳게) 지목 — 담체는 **전 링 선-주사한 인접성**에서 읽는다(실측: debug_assert 발화가 잡음). ② «전 생산 직선 비트 동일» 주장이 이동 경로에서 반박 — pass 2 는 방향을 직접 회전, 파생은 끝점 차 재정규화라 방향 ~1 ulp(실측 2.8e-16, 직선 83/84 비트 동일, 원 최대 2.2e-16 — 직선 기하는 비관측이라 무해). ③ VertexOffCurve 의 직선 갈래는 **타입상 항진**이 됐다(끝점이 자기 직선 위) — 검사는 원(rim)으로 이빨 유지, `.max(tol_of(edge.origin))` 은 상수 `EPS_CONSTRUCTED` 로 재철자(**무-행동이 아니었다** — `Discovered{tol:0}` 정점의 하한을 edge 항이 받치고 있었음, 실측). census 전 커밋 비트 동일 | ✔ 2026-08-06 |
| S7 | **`Origin` 소멸 — 정점은 자기 정의를 들고, 좌표는 캐시가 된다** — `Vertex{def: VertexDef{ThreePlane([3]) \| OnSeam([2])}}`(Q3 수정: seam 정점이 단일형을 반박 — 반증표), `point`·`Origin`(3변종) 사망, `vertex_cache: Vec<PointCache{coord, tol: Option<f64>}>` 인덱스-평행 + 유일 입구 `push_vertex` + 접근자 `vertex_point`/`vertex_tol`. **`rebuild_vertex_cache` 는 없다**(3b ⏸ — 발견 좌표는 배열이 공들인 값 1992 중 238 이 순진 Cramer 와 다르고, seam 좌표는 load-bearing): S8 이 모서리에서 얻은 «버리고 재생» 보증은 정점엔 아직 없음을 정직 기록. 소멸한 기계: 스케치 프레임 base 정점(Q2 — 프레임 공유 세 평면의 유리수 Cramer + 사슬 재생이 저장 좌표를 **비트 동일**로 재현, 8/8 실측을 영구 잠금으로 승격)·`remap_origin`·`solid_motion`+정점용 `move_node`(면이 자기 leaf 를 든다 — 규칙 3)·한-홉 base 불변식(타입이 흡수: 중복 적용이 표현 불가)·`exact.rs::base_f64/top_f64`. reuse `solid_points` 는 def 경로로(구성=`Pt3::exact` 문자 동일, 발견=포기 문자 동일, 이동=세 이름의 checked-i128 Cramer→replay; **혼합 프레임은 정직한 decline** = 기록된 유일한 차이). 게이트 `origins_are_remappable`→`defs_are_remappable` 전 정점 확장(발화 0 + **양성 대조**), validate: `tol_of` 1식화·`VertexOffDefinition` **전 정점 확장**(+양성 대조)·신설 `VertexDefCarrierMismatch`(변종 ⇔ 담체 종류). 신설 `nacre_scalar::three_planes_rat`. census **전 커밋 비트 동일**(좌표 verbatim 이사 — 재기준 없음) | ✔ 2026-08-07 |

### 남은 항목 — **순서는 다음 계획에서** (선행 관계만 적는다)

| | 항목 | 선행 |
|---|---|---|
| S5(i)-a ✔ | **datum 평면 연산** — `Operation::DatumPlane{def: DatumDef{Stated(SketchPlane) \| Offset{frame, dist}}}` + `OpOutput::DatumPlane{plane, frame}` + `Model::surface_handle_at`(좁은 읽기 접근자, S1 봉인 무손상) + `rebind` 의 `Handle<Surface>` 갈래(`LogCell::Surface`). 새 공개 생성자 없음 — 다섯 `SketchPlane` 생성자가 곧 datum 의 어휘다. 규약: 캐시 법선 `−(진술된 법선)`(근거 둘 — S9 의 781, `frame_sign`), 배치는 `Stated` 면 **무조건 `Named`**(ZX 의 `+u = +ẑ` 가 유도를 금지), `Offset` 은 push 전 3중 정규화(flip→부호 / 세계 되당김 / `dist == 0` 은 `ZeroOffset` 거절)로 «한 평면에 한 핸들»을 지킨다. ★ **앵커 실측**(`tests/plane_anchor.rs`): 기울어진 평면의 저장 `d` 는 앵커마다 최대 **22 ulp** 다르고 그 넷은 **링의 네 점**이라 이 흔들림은 datum 이전부터 있었다; `spans_exactly` 는 다섯 앵커 **전부 false** 라 판정 경로가 앵커를 읽지 않고, 최악 앵커에서도 위상 동일·부피 비트 동일·좌표 **1.1e-15**(ε 의 1/4400). 축 정렬 평면은 `d` 가 곱 하나라 앵커에 **무감각** — 그래서 «먼 원점 + world_xy» 픽스처는 아무것도 재지 못한다 | ✔ 2026-08-07 |
| S5(i)-b ✔ | **`Operation::Extrude` 의 평면 어휘 교체** — `plane: SketchPlane`(값) → `frame: SketchFrame`(핸들+배치+측정된 flip). 구현은 pad/pocket 도로 재사용(`realized_plane` → `swept_profile` → `build_prism(base_cap_surface = Some(프레임의 평면))`), `rebind` 의 `Extrude` 팔이 재고정으로(**7/7 변종이 핸들을 싣는다 — R 의 전제 소멸**), 신설 `SketchFrame::world(model, Axis)` 가 `world_zx` 예외를 한 곳에 가둔다. 호출부 112곳/4크레이트, 단일 커밋. ★ **방향은 `flip` 이 든다**: 같은 평면을 `+n`/`−n` 으로 진술하면 **한 핸들이고 두 프레임 모두 정준 `ŵ`** 를 보고하므로(실측), datum 이 «호출자가 진술한 법선» 에 대해 flip 을 재야 교체가 무행동이 된다. ★★ **C2 가 교체 전에 갈림을 잡았다**: 실현한 축을 다시 lift 하면 정규화가 필요한 축(`(0.6,0.8,0)`→원시`(3,4,0)`, `uu=25`)이 `0.6000000000000001` 로 돌아와 `exact()` 가 뒤집히고 평면이 **조용히 프레임-노드 도로로** 간다(ulp 가 아니라 다른 아레나). 수리 = `RatFrame::of_plane_frame`(유리수로 묻는다, `inv_sqrt_exact`). 내 가설 둘이 죽었다 — 원점 상쇄가 원인이 아니고(단위 축은 이동해도 살아남는다), «`|u_raw|²=1` 이니 정확» 은 **이미 단위인 축에만** 참이었다. C3 차등: 9 평면 가족 × 2 프로파일 = **18/18 아레나 동일 + 노드 수 일치**. census: `c ` 는 **기울어짐 9줄만 이동**(위상 이동 0, 최대 상대 편차 **4.4e-16** = ε 의 1/2000), `wide_planes` 불변, ★ **`seeded_hits` 는 455 불변 — 하락 예측이 틀렸다**(코퍼스에 세계 평면 extrude 가 **0개**이고 그 적중은 `add_cuboid` 의 면이었다; 카운터의 출처를 오독했다). ★ 훅이 못 본 19곳(oracle 은 전부 `#[ignore]`)을 손으로 돌려 **결함 하나를 잡았다** | ✔ 2026-08-07 |
| S5(ii)-1 ✔ | **`PlanePoints::Through([Handle<Vertex>;3])` + 생산자 `DatumDef::ThroughVertices` + 유리수 닫힘 판정** (2026-08-08, 커밋 `ba52b8b`·`3c968f9`·`5c9b3ab`). 이름은 **push 시점에 한 번** 유도한다(`Model::through_points_rat` → `plane_name_exact`) — 그래서 `frame_chain`·far cap·`base_rat`·모든 술어가 **한 줄도 안 바뀐다**. 판정 표는 같은 유도로 증인 삼각형을 **그 자리에서** 만든다(연산 하나 동안 사는 거울이므로 규칙 1 위반이 아니다). ★ **정렬은 키에만, 방향은 호출자의 정점 순서** — `dist` 가 양수 전용이라 순서가 유일한 방향 선택이고, 뒤집으면 «같은 핸들 + 반대 프레임» 이다. ★★ **`transform` 은 핸들을 그대로 두고 노드를 기록한다**(정점은 복제되지 브릿지되지 않는다); 정점이 base 를 정하고 모션이 옮기며 **둘은 더해질 뿐 곱해지지 않는다**. 대가: 그런 datum 위의 솔리드는 정확한 강체 이동에도 **항상** 노드를 얻는다. ★★★★ **가장 조용히 틀릴 뻔한 자리는 컴파일러가 못 본 곳이었다** — 변종 추가가 낸 비망라 에러는 **2개**뿐이고, `transform::points_move` 의 `let`-`else`(원통용 폴백)는 거기 없었다. 원통은 진실이 기하를 **안 들어서** «나를 것 없음» 이 참이지만 `Through` 의 진실은 기하를 **참조로 든다** ⇒ 같은 답을 하면 노드 없는 경로로 가 **캐시만 움직이고 진실은 제자리**에 남는다. `false` 를 반환하게 고쳤다 | ✔ 2026-08-08 |
| S5(ii)-2a ✔ | **증인 없는 평면의 구간·고정밀 계수** (2026-08-08, 커밋 `7acf5c5`). `plane_iv_through`·`plane_hp_through` — 동차점 셋의 **사영 join**. **차수 15 → 9**(분모 털고 아핀 외적이 아니라 3×4 의 네 소행렬식), `cramer_iv` 의 **쌍대**라 새 산술 0. ★★★★ **초안은 설계가 금지한 연산을 제안했다** — `Dvec/D` 로 나눠 아핀 좌표를 만들려 했고, 문서가 *"나누면 무리수가 되고 오차가 낀다"* 며 이미 금지한 것이었다. **타입이 이미 강제하고 있었다**: `Iv` 에는 나눗셈이 없고 `HpIv::div_exact` 는 반경 0 을 요구한다. ★★★ **배율의 부호는 값 안에서 없앤다**(실을 자리가 시그니처에 없다 — 실을 곳 없는 값은 아무도 안 쓴다). ★★ **실측이 단계를 갈랐다**: 깊이 1 은 필터가 살지만(미결 0/800, 최악 상대 반경 6.3e-10) **깊이 2 는 차수 81 이라 계수가 `f64` 밖 8/8** ⇒ 상승 전용. `NaN` 반경이 우연히 옳게 굴러가는 것에 기대지 않고 `None` 으로 명시한다. 생산자 없음(census 비트 동일), 잠금은 `three_planes_big` 과 같은 모양의 **차등 테스트** | ✔ 2026-08-08 |
| S5(ii)-2b | **배선** — 위 기계에 생산자를 붙인다: 판정층 `WorkingPlaneDef::Through`, `Witness::tri_pt3` 의 전(total)성 제거(`plane_def` 를 지나는 **7개 술어 지점**). 오늘은 `VerticesInMixedFrames` 로 **이름 붙여 거절**한다. ★★ **선결 조건**: 이 단계가 «이름 없는 평면» 을 허용한다. ★ 내가 여기 적었던 «`frame_chain` 의 조용한 f64 강등을 먼저 없애라» 는 **틀렸다** — 그런 강등이 없다: `motion_chain` 의 `None` 소비자 넷이 전부 정직하다(`planes.rs` ×2 는 `FrameOutOfRange`, `exact.rs` 는 `PlaneWithoutExactForm`(S6b 가 f64 프리즘 폴백을 죽였다), `reuse.rs` 는 동등하게 옳은 배열 도로, `rotated_vertex.rs` 는 `Pt3Error::Downgrade`). 그 함수의 doc 주석이 낡았을 뿐이다. **진짜 선결 조건은 하나**: 이름이 없으면 **interning 도 없으므로**(§이름과 interning: *"`!=` 는 아무것도 증명하지 않는다"*) «한 평면에 한 핸들» 을 무엇이 대신 지킬지 그 단계가 답해야 한다. ★★★ **2a 가 남긴 네 질문**: ① **인덱스 공간** — `Through([[usize;3];3])` 이 가리키는 표가 무엇인가. 면 표는 **두 피연산자의 면**(`plane_index_setup` = `collect_planes(a)` ++ `collect_planes(b)`)이라 **제3 의 솔리드**에서 온 정점의 담체는 거기 없고, `WorkingPlane` 은 애초에 **클래스** 표(`dense_planes`)의 타입이다 ⇒ 표에 덧붙이거나(«synthetic faces appended later» 기계가 있다) 그 경우를 거절하거나 다른 인덱스 공간으로 적는다. ② `standard_for(&planes)` 가 표에서 정밀도 기준을 유도하므로 **표가 늘면 기준이 움직인다**. ③ 위 interning 질문. ④ `Witness::tri_pt3` 의 전성 제거 — 7개 술어 지점 + `FaceInfo`·`PlaneGeom` 두 구현체. ⑤ **깊이 제한** — 깊이 2 가 상승 전용이므로(2a 실측) `OnceCell` 메모가 «평면당 한 번» 을 실제로 보장하는지 확인하거나 제한을 둔다 | 2a ✔ — 남은 것은 배선 |

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

★ 판정층 개명(`Pt3`→`WitnessPoint`·`WorkingPoint`→`WorkingVertex`·`PlaneGeom`→`WorkingPlane`,
§판정 이름 규칙)은 별도 단계가 아니라 **각 타입을 처음 만지는 단계에 얹는다** — 기계적 개명이라
관문은 컴파일이다.

### 관문 규칙

- 기본 관문은 **위상 정확 일치 + 좌표 ε(모델 크기 상대, `2⁻⁴⁰` — 실측 여유 2000배)**, 비트
  동일은 보너스 신호. 합격/불합격이 아니라 **최대 편차 숫자를 찍는다**(누적 드리프트 감시).
- ★★ **8b — 비트 동일 관문은 그 안의 population 에 대해서만 보증한다.** 바꾸려는 것이 닿는
  population 을 대장에 먼저 넣고(17자리 `fw`·기울어진 `tp` 가족은 이미 있다), *"어느
  population 인가"* 는 추측하지 말고 계측이 이름을 대게 한다. 좌표 관문은 «답은 같은데 더 나쁜
  길로 갔다»를 원리적으로 못 보므로 **"정확 경로를 탔는가"를 직접 단언하는 테스트**를 함께 둔다.
- **병렬 불변식**: 병렬 구간에서 `push_surface` 금지(재생 결정성) — 오늘은 전 호출부가 순차,
  M6 교차 곡면에서 다시 본다(그때는 정준 키 정렬 후 일괄 push).

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

0. ★★★ **`three_planes_rat` 이 «드는 점»을 거절한다 — 타입이 아니라 풀이의 문제**
   (실측 2026-08-07). 그 함수는 유리수 여인수 전개를 `checked_*` 로 하는데, 곱셈은
   `num-rational` 이 교차약분하지만 **덧셈·뺄셈은 `a.num·b.den ± b.num·a.den` 을 만들며 약분
   전에 넘친다.** 그래서 그 `None` 은 **두 사실을 뭉친다** — 「점이 `Rat` 에 안 든다」와
   「중간값이 넘쳤다」.

   기울어진 십진 가족에서 **받은 정점의 2/3을 거절하고, 그중 100%가 `Rat` 에 드는 점**이다
   (`declined_and_wide = 0`, 최대 폭 59비트).

   ★ **고칠 수법이 이미 크레이트에 있다**: `plane_name_exact` 는
   `narrow().or_else(plane_name_big())` 이고, `plane_name_big` 은 분모를 한 번에 걷어내
   **정수 산술 + 마지막에 gcd 하나**로 풀고 **맞으면 `Narrow` 로 되돌린다**. 정점 풀이에는 그
   절반이 통째로 없었다. 쌍둥이 `three_planes_big`(+`MeetPoint`)이 이제 서 있고 차등
   proptest 로 잠겨 있으므로, **수선은 `three_planes_rat` 에 `.or_else` 한 줄**이다.
   S5(ii) 와 **독립**이다. — 아직 연결하지 않았다.
   ★★ **정정(2026-08-08)**: 여기 «`reuse.rs` 의 정확 경로를 **즉시** 넓힌다» 고 적었던 것은
   **검증되지 않았고 오늘은 아마 0 이다**. `three_planes_rat` 의 **프로덕션 호출자는
   `reuse.rs:101` 하나뿐**이고(나머지 둘은 `#[cfg(test)]`), 그 함수는 **발견 정점을 먼저
   거절한다**(`vertex_tol().is_some() ⇒ None`) — 2/3 거절을 잰 인구가 거기 도달조차 하지 않는다.
   원인을 가르는 것 자체의 옳음은 그대로다(결정적 테스트가 잠갔다). **reuse 를 실제로 넓히는
   것은 열린 항목 13**(`ImplicitPoint` 로 갈아타기)이다.

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
   한정). 느릴 뿐 틀리지 않는다.
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
8. **최종 개명은 유예됐다(S6b)** — 진실 enum 은 `SurfaceTruth`(문서의 `Surface` 이름은 아직
   geom 의 f64 캐시가 쥠), 캐시는 `Store<geom::Surface>` 그대로. `Handle<T>` 의 타입
   매개변수가 아레나 반전을 강제하고(공유 인덱스라 의미 무손실), `SurfaceCache{coeffs,tol,
   inv_norm}` 실형이 판정 통합에서 생길 때 개명·반전을 한 번에 기계적으로 한다. 판정층
   개명(`PlaneGeom`→`WorkingPlane` 등)도 같은 자리(또는 S7)로 유예.
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
   않았다: 3b(좌표 재생)가 ⏸ 이고, 발견 좌표는 배열이 공들여 만든 값(1992 중 238 이 순진
   Cramer 와 다름)이며 seam 좌표는 M6 까지 load-bearing 이다. 정점 `D≠0` 의 완전한 유리수
   단언(문서 :85 의 «공짜 단언»)도 같은 이유로 유예 — 오늘은 `push_vertex` 의 핸들 상이성
   debug_assert 까지. 좌표 재생이 생기는 자리(M6/판정 통합)에서 셋을 함께 연다.
13. **reuse 의 발견-정점 갈래는 아직 포기다(S7)** — `solid_points` 는 측정 좌표를 만나면
   `None`(→ Arrange). `nacre-cip` 의 `ImplicitPoint`(세 평면의 암시적 점)로 갈아타면 융합
   fold 의 클래스 재사용이 불리언 **결과** 피연산자에도 걸린다 — 판정 기계 교체라 행동
   변화이고, 별도 측정과 함께 여는 항목.
14. **혼합-프레임 정점은 reuse 가 답하지 못한다(S7 기록)** — 호출자가 세계 좌표로 명시한
   밑캡 위의 프레임-스케치 코너는 세 담체의 모션이 갈려(둘은 프레임, 하나는 세계) 유리수
   pullback 이 없다. 옛 `Origin` 길은 base 정점으로 답했지만 그 정점이 S7 에서 소멸했다.
   실측 인구는 `the_def_road_answers_for_the_populations_it_can_name` 의 ④ 가 핀한다.
15. ★★★ **술어가 `Wide` 계수를 못 읽는다 — §판정이 «공짜» 라고 적은 그 경로가 없다**
   (신설 2026-08-08). `exact_coeffs()`·`base_coeffs()` 가 `narrow()` 만 읽으므로 wide 이름
   평면은 정확 지름길을 전부 잃고 톨러런스+상승으로 간다. 실측 대가는 **상승 3.2배**(278→887,
   고갈 0). 고치는 방향은 `plane_name_big` 이 이미 보인 것과 같다 — 분모를 걷어낸 **정수**
   계수를 `Expansion` 에 넘기는 것이고, `PlaneName::Wide` 는 이미 그 정수를 들고 있다.
   ★ 우선순위: 절벽이 아니므로 **S5(ii)-2 를 막지 않는다**. 다만 S5(ii)-2 가 wide 이름 인구를
   키우면 이 배수가 곱해진다.

10. **감김(`oriented_ring`)은 아직 f64 다** — 배치된 3D 점의 면적벡터·법선 내적(ops).
   `check()` 가 단순성(≠0 면적)을 진실 위에서 보증하므로 지금은 건전하지만, f64 폴백 소멸
   (S6 이후)과 함께 재검할 것.
