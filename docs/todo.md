# todo

맨 위 항목이 현재 위치다. 끝난 항목은 체크하지 않고 **지운다** — 기록은 커밋이 든다.
항목은 번호가 아니라 이름을 갖는다.

## 지금

### 모션 사슬을 읽을 때 접는다

**무엇이 문제인가.** `Model::chain_translation` 은 사슬의 **첫 비이동 노드에서 사퇴**하고, 그 사퇴가
`Model::world_plane_name` 을 통해 읽는 곳 여덟 자리(호출 11회)로 퍼진다: 캐시 유도의 게이트
(`Model::derive_surface_cache`) · `Model::vertex_meet_of` 안의 `world_road` 클로저 ·
`planes::world_plane_coeffs` · `planes::pierce_corner` · `planes::disk_of` ·
`transform::pierce_line_reversed` · `combinatorics::curved_wall` · `combinatorics::pierce_name_from_def`.
구현은 아직 없다.

**거절 532 의 원인별 분해.** 곡면 캐시 유도는 5,317 에서 성립하고 532 에서 거절한다(거절 = 모션 456 ·
이동된 원통 76; 넘침 0 · 이름 없음 0 · wide 0). 그 532 를 사슬의 구성으로 쪼갠 값(거절이 나는
`Model::decline_reason` 바로 그 자리에서 분류한다 — 모델을 도는 탐침은 push 인구가 아니라 «살아남은
곡면»을 세어 명제 옆을 잰다):

| 부류 | 평면 456 | 원통 76 |
|---|---|---|
| `Frame` 노드를 든다 | **169** | 0 |
| 임의 각 회전을 든다 | **197** | **7** |
| 이동 + **사분각** 회전뿐 | **24** | **69** |
| 위 + **거울**까지 | **66** | 0 |

169+197+24+66 = 456(정확히 분할된다). 있는 기계로 열리는 인구는 **90 평면 + 69 원통 = 532 중
159(30%)** 다.

- **영구히 밖**: `Frame` 169 — 기저가 `1/√유리수` 다(`chain_fixes` 의 doc 이 못 박은 경계).
- **임의 각 197(+원통 7)**: Niven — 「무리수 모션이 낀 datum 은 이름이 없다」 항목의 인구다.

**접는 것은 아레나가 아니라 «읽기»다.** `self.motions` 를 만지는 자리는 저장소 전체에 `push` 와 `get`
둘뿐이고(수정·삭제 없음), 노드 핸들은 `SurfaceKey`·`ThroughKey`·`CylinderKey` 의 interning 열쇠라
노드를 합치면 사라지는 것이 자식이 아니라 **동일성 판정 전체**다. `chain_translation` 의 doc 이 그
규칙이다 — *"it does not license dropping the history."* 그러므로 이 일은 **노드를 하나도 만들지 않고
하나도 지우지 않는다.**

**산술은 다 있다 — 없는 것은 «노드를 차례로 적용하는 걷기»뿐이다.** `Isometry::plane_coeffs`(회전은
`try_exact_cos_sin` 경유, 피벗·이동 포함) · `point_rat` · `dir_rat` · `mirror_plane_coeffs` ·
`mirror_point_rat`. 모양은 기존 가족을 따른다(`chain_fixes` 라는 private 걷기 → 공개 얼굴 둘).
계획: **`chain_plane_coeffs` / `chain_point` / `chain_dir`** 셋이 한 걷기를 공유한다 — `Isometry`
자신의 삼총사와 같은 절단면이다. 새 스칼라 함수는 필요 없다.

- **반환형은 `Isometry` 일 수 없다.** 거울은 `det = −1` 이라 `Isometry` 에 자리가 없는데 그 인구가
  **66** 이다. `SurfaceDeriveCounts` 의 doc 이 적은 처방 *"a motion wants a chain folded to an
  `Isometry`"* 는 그 66 에 대해 틀렸다 — 이 일과 함께 고친다.
- 거울의 패리티는 **이름이 흡수한다**: `canonical_plane_coeffs` 가 첫 비영 성분을 양으로 강제하므로
  정준 이름은 방향을 말하지 않는다 — 평행이동에서 그러는 것과 같다.
- `chain_dir` 은 **거울을 만날 수 없다**: 방향을 옮기는 것은 원통뿐인데(평면은 계수로 움직인다)
  `OpError::MirrorNotPlanar` 가 미러링을 먼저 거절한다(`nacre-ops/tests/edge_carriers.rs` 가
  *"A mirrored cylinder has no population"* 으로 적어 뒀다).

**이 일은 결과를 움직일 수 있다.** `derive_surface_cache` 가 `world_plane_name` 으로 게이트하므로
**캐시 확장과 이름 확장은 쪼갤 수 없고**, `vertex_meet_of` 의 `world_road` 가 자동으로 넓어져 코너가 새
datum 도로를 얻는다(그 doc 의 *"인구를 열 수는 있어도 움직일 수는 없다"* 계약이 이 자리다).

- **음성 대조군은 트리에 있다.** topo 의 테스트 `a_corner_of_two_translation_chains_solves_in_the_world`
  가 **90°** 사슬에 대해 *"a rotated carrier has no world name to meet with"* 를 단언하고, 바로 위 주석이
  스스로를 반박한다(*"A quarter turn carries x = 0 to y = 0, so the cache is exact … what declines here
  is the rotation"*). **그 줄이 빨개지는 것이 «접기가 작동한다»의 증거다.** 다시 겨눌 때 **37° 짝을 같은
  파일에 남긴다** — 무엇이 여전히 사퇴하는지 보여 주는 대조가 없으면 운 좋은 통과와 구분되지 않는다.
- **«이름이 생긴다» ≠ «길이 열린다» — 둘째 관문.** `combinatorics::curved_count` 의 원통 관문은
  *"refuses any class that is rotated or has no narrow rational name"* 이라, 회전된 클래스는 이름을 얻고도
  계속 거절된다. `FaceInfo::rotated` 는 «세계 진술이 없다»와 «행의 `tri` 가 진실이 아니라 실현본이다»를
  **겸하므로 건드리지 않는다**(`planes::world_plane_coeffs` 의 doc: 그 플래그를 뒤집으면 공유 벽이 두 몸으로
  갈린다). 그래서 인구를 `declined_motion` 의 감소로만 세면 **과대평가**다: «이름을 얻은 곡면»과 «실제로
  답이 달라진 census 행»을 **따로** 센다.
- **폴백 순서는 계약이다.** `world_road` 는 **넓히되 승격하지 않는다** — 먼저 태우면 답이 이미 있는
  코너의 점 철자가 바뀐다.

**같은 축 회전의 각도 합성은 역량으로만 적는다.** 사슬은 앞부분만이 아니라 전체를 봐야 하고, 같은
축·같은 피벗의 회전은 각도를 더해 접을 수 있으며, 30°+60°=90° 면 유리수 표현이 되살아나 주변 이동까지
함께 접힌다 — `Angle(Rat)` 이 도 단위 유리수이고 `checked_add` 가 정확히 mod 360 축약하므로 타입
수준에서 성립한다. 그러나 코퍼스에 **같은 축·같은 피벗의 인접 회전쌍이 0** 이고 **사슬 깊이 최대 2** 다
(`chain_motion` 이 한 호출당 최대 둘(회전→이동)을 남기고, 픽스처가 한 솔리드를 두 번 회전시키지 않는다).
이것은 «픽스처가 안 한다»이지 «사용자가 못 한다»가 아니다 — 같은 축으로 두 번 돌리면 즉시 생긴다.
위 설계는 이것을 막지 않는다(노드별 적용 위에 «인접 정규화» 한 겹을 얹는 모양이라 되돌릴 것이 없다).
같은 축이라도 **피벗이 다르면** 각도만 더해지지 않는다 — 합성 피벗에 cos/sin 이 든다.

## 다음 — f64 는 실현 통로 하나로

숫자 규칙 「실현 통로는 하나다」를 코드가 아직 다 지키지 않는 자리들.

### 실현 함수가 세 갈래다

규칙은 통로 하나(`정의 → 좌표`)인데 코드에는 실현이 세 곳에 있다: `nacre_ops::realize_vertex`(정점, `realize.rs`), `nacre-judge` `frame3` 의 `WitnessPoint::realize(prec)`(판정의 증인), `nacre-ops` `exact.rs` 의 `realize(pts)`(구성). 여기에 스칼라의 `to_f64` 가족과 `_tracked`·`_memoized`·`_rounded` 변종이 붙는다. 먼저 셀 것: 각 함수의 호출처와, 같은 정의를 두 갈래가 실현할 때 비트가 같은가. 목표 모양은 스칼라 → f64 하나(`nacre-scalar`), 정의 → 좌표 하나이고, 그 밖에서 f64 좌표를 짓는 것은 가시성이나 clippy `disallowed_methods` 로 컴파일 단계에서 막는다.

### 픽스처는 제품 도로로

**무엇이 문제인가.** topo 의 `test-util` 문 `add_cuboid` 가 **804곳(51 파일)**, `add_cylinder` 가
**182곳(23 파일)** 에서 불린다. topo 는 ops 를 모르므로 둘 다 실현 없이 `PointCache::Unrealized` 로
push 하고, census 의 불리언 코퍼스 전부가 이 상자·원통으로 서 있다. 존재 이유는 «topo 층 테스트가 ops
없이 솔리드를 원했다»와 «스케치 도로보다 코퍼스가 먼저 있었다» 둘이다(`Model::new()` 의 씨앗은 세계
평면뿐, 정점은 없다).

**모양.** 테스트 지원 모듈에 **같은 서명의 `cuboid(m, min, max)`·`cylinder(m, …)`** 를 두고 안에서
`apply(Extrude)`(사각형·온전한 원 프로파일)를 부른다. 모든 픽스처가 push 깔때기
(`push_vertex_realized`)를 지나 `Bounded` 가 되고, «깔때기를 안 지난 push» 라는 사유가 제품에서 사라지며,
census 가 «제품 도로의 census»가 된다. topo 의 `add_cuboid`/`add_cylinder` 는 은퇴한다. 986곳 치환은 같은
서명이라 기계적이다.

**남기는 것 하나 — 심기용 날것 문.** validate 의 자기 테스트는 **틀린 모델**(매달린 핸들·뒤집힌 면·정의와
어긋난 캐시)을 심어야 하고 제품 도로로는 틀린 모델을 만들 수 없다. 그래서 `push_vertex`(날것 캐시)·
`push_plane_unregistered` 같은 문은 `test-util` 아래 topo 에 남는다 — 제품에 없는 것을 테스트용으로 두는
것이 정당한 유일한 경우다.

**위험은 census 가 움직이는가다.** 압출 상자의 평면 이름은 `add_cuboid` 와 같은 정준형이라 **비트 동일이
예측**이고, 원통은 `add_cylinder` 의 seam 모형과 압출 원의 seam 이 같은지 A/B 가 답해야 한다.

### 곡면의 실현과 내보내기의 수치 정밀도

**무엇이 문제인가.** 정점은 태어날 때 정의에서 실현되는데(`push_vertex_realized`), 곡면의 «방향»은
생산자가 넘긴 f64 그대로다. `nacre-step` 은 면의 곡면을 `Model::surface_cache` 에서 그대로 꺼내
`normal`·`ref_dir`·`radius` 를 내보낸다(앵커 — `AXIS2_PLACEMENT_3D` 의 location — 만 진실의 첫 점에서
나온다). 계획: `realize_surface(s, Precision)` — 코드에 없다. 정점 쪽 문
`nacre_ops::realize_vertex{,_decimal}` 과 `Precision { NearestF64, Bits(usize) }`·`Realized` 는 서 있고,
자릿수는 `Realized::to_decimal` 만 든다.

두 손잡이를 가르는 것이 이 항목의 핵심이다.

| | 무엇 | 상태 |
|---|---|---|
| **A. 세밀도** | 삼각형 개수 — `TessConfig { tol, max_angle_deg }`. `tol` 은 현 편차(sagitta), `max_angle_deg` 는 「2° 면 800px 원이 0.06px 안」 — 순전히 **보기**의 수 | 있다 |
| **B. 수치 정밀도** | 각 좌표를 몇 비트로 **실현**하나 | 정점에만 있다. 곡면의 방향·곡선 위 샘플점에는 없다 |

- **B 가 여는 것은 «방향»이다.** 평면의 단위 법선, 원통의 축과 `ref_dir` 은 전부 `normalize()`
  (`self / n2.sqrt()`)를 타므로 유리수가 아니고, STEP 이 내보내는 값이 바로 그것이다.
- **행(`raw`)은 아직 생산자의 것이다.** 저장된 `raw` 가 «실현된 진실 점들의 외적»과 같은 것은 **2,747**,
  다른 것이 **2,010** 이다(회전 사슬 456 은 답할 수 없음). 행을 진실에서 가져오는 일은 앵커(캐시 비트가
  바뀐 것 61)와 달리 **약 2,010개를 움직이는 큰 변경**이다.
- **넘치는 자리는 수선의 발이다.** 이름 유도는 `plane_name_exact` → `plane_name_big` 으로 폭을 넘는다.
  넘치는 것은 `plane_origin_projection` 이 `n·n` 을 만들어 계수를 제곱하는 길이고(코퍼스 8건), 앵커는
  그 길을 쓰지 않는다. 남은 것은 **방향의 정확 반올림**이다.
- **가장 값 하는 곳은 정점이 아니라 곡선 위 샘플점이다.** `nacre-tess` 는 `model.vertex_point(v)` 로 f64
  캐시를 읽고, 원·원통을 잘게 나눈 점은 `cos/sin` 평가라 정밀도가 그대로 드러나며 **OBJ 는 거의 전부 그
  샘플점**이다. 정점은 대부분 안 움직인다: 축정렬·불리언 인구는 캐시가 이미 최근접(48/48 · 96/96 — 두
  길이 같은 유리수를 같게 반올림하니 어긋날 수 없다), 기울어진 프레임 인구는 움직인다(0/12 최근접,
  최대 4 ulp; 36 좌표 중 29). 가르는 것은 f64 표현 가능성이 아니다 — 축정렬 48 좌표 중 dyadic 은 32 뿐인데
  48 전부 최근접이다. 폭(7비트 vs 59비트)은 «산술이 얼마나 있었나»의 대리 지표이지 기제가 아니다.
- **순서: B → 그 다음 자릿수 지정.** B 없이 「STEP 텍스트에 f64 보다 많은 자리」는 거짓말이다(f64 의
  17자리 뒤는 이진 반올림의 부산물). 거꾸로 하면 없는 정밀도를 출하한다.
- **기계는 있다.** `nacre-scalar` 의 `HpBounded`(`add`/`sub`/`mul`/`div`/`inv_sqrt` 전부 `prec` 를 받는다) ·
  `WitnessPoint::realize(prec)` · `round_to_f64`/`round_to_digits` · `sqrt_bounded`·`realize_quad`·
  `realize_seam_point`. 정확 반올림이라 실현이 유일하므로 **같은 정밀도면 답이 하나**다.
- **잠금이 공짜로 딸려 온다.** 같은 모델은 **바이트 동일**한 파일을 내고, 실현 정밀도를 **두 배로 올려도
  f64 결과가 안 바뀌어야** 한다. 그 둘이 「정밀도 기능이 아니라 정확성 고침」의 증거다.
- 내보내기는 `&Model` 서명을 지킨다 — `&mut Model` 을 받으면 이름이 약속하지 않은 일을 한다.
- 이름: 실현 문은 `realize_*` 이고 정밀도를 항상 명시한다(기본값 없음). `_at` 은 이 커널에서 위치의
  낱말(`point_at`·`normal_at`·`surface_handle_at`)이라 쓰지 않는다. 판정은 이 문을 안 부르고 자기 어휘
  (`hp_coord`·`judge_precision`·`trial_bound`)를 쓴다.

### 정점 캐시의 «버리고 재생» 보증

**무엇이 문제인가.** 간선에는 `rebuild_edge_cache` + `edge_cache_discard_and_regenerate_bit_identical` 이
있는데 정점판이 없다. 정점 `D≠0` 의 완전한 유리수 단언도 같은 이유로 없다 — `push_vertex` 의 핸들 상이성
`debug_assert` 까지다.

**지금 서 있는 계약.** 연산이 만든 모든 정점에 대해 `vertex_point(v)` 는 `realize_cache` 의 답과 비트
동일이거나, 그 길이 이름으로 거절한 것이다 — census 가 매 행에서 단언한다. 캐시는 태어날 때 실현되고,
둘째 쓰기 문은 `refine_vertex_cache` 하나다(`Ceiling` 만 `Bounded` 로 올린다; 내보내기 직전에 부르는 비싼
문이고, 뒤이은 연산은 캐시를 읽어 진실이 되는 것을 정하므로 그 뒤로는 연산하지 않는다). 남은 것은
**캐시를 통째로 버리고 정의에서 다시 세워도 비트 동일**이라는 보증이다. 실현값↔담체 캐시 거리는 최대
1.07e-14 로 validate ε 의 다섯 자릿수 아래다.

### `surface_name` 곁표 접기

**무엇이 문제인가.** `pub surface_name: HashMap<Handle<Surface>, PlaneName>` 은 `Model` 의 마지막 공개
필드다. 도착점은 `SurfaceCache` 의 평면 팔이 이름을 드는 것. 인구: 제품 17 · 테스트 53. 쓰는 자리는
`push_plane_raw` 한 곳이고, 이름은 유도된다 — *"The name is derived, so it cannot disagree with the thing
it names."*

**무엇이 위험한가.** 이 필드는 **유일하게 조용한 실패 모드**를 가진다: 이름이 있어야 할 자리에 없거나
없어야 할 자리에 생겨도 **census 가 못 본다**(그 인구를 안 든다).

**재생 잠금의 모양.** 캐시 전체 재생은 서지 않는다 — `derive_surface_cache` 가 `raw` 를
`surface_cache(h)` 에서 복사하므로 재생이 «지우려는 값»을 먼저 읽는다(막는 것은 향이 아니라 행이다).
«앵커만 버리고 재생해도 비트 동일»은 선다 — 출발점이다. 향을 점 순서에서 뽑는 길은 닫혀 있다: 비씨앗
`Known` 평면 **342개**가 자기 점 순서와 반대 방향의 캐시를 단다. 계획: 이름만 재생하는
`rebuild_surface_names()` + `surface_name_discard_and_regenerate_bit_identical`, 제약 둘 —

1. **제품 도로 전용**: `push_plane_unregistered` 는 이름을 건너뛰므로 rebuild 가 이름을 **추가**한다.
2. **오름차순 핸들 순서**: `Through` 이름이 담체 이름에 의존하는 자기참조이고, append-only 라 의존이 항상
   낮은 인덱스다.

### `SurfaceCache` 를 enum 으로

**무엇이 문제인가.** `SurfaceCache` 는 `realized: nacre_geom::Surface` 하나를 감싼 struct 다. 도착점은
`enum SurfaceCache` 와 `.plane()`/`.cylinder()` 체이닝. 근거는 «확정한 구조의 이행»이지 «재발 방지»가
아니다(막을 인구는 0 이다). 인구: 값 21 + 테스트 44, 메서드 이전 8, 통째 출구는 `transform.rs` 한 곳.
테스트 30곳의 종류 판별도 함께 지나간다.

**가장 위험한 한 줄은 `Model::flipped_against`** — 캐시 normal 로 interning 의 `flipped` 비트를 정하고
그 비트가 **스위트에서 1,916번** 발화한다. 잘못 재철자하면 면 방향이 조용히 뒤집힌다.

## 다음 — 타입

### 원·원호의 진실: 무리수 중심

**지금 참인 것.** 반지름의 진실은 r² 다: `Edge2d::Arc { center: [Rat; 2], r2: Rat, ccw }` ·
`CylinderDef { origin, dir, ref_dir, r2: BigRat }`. r 은 실현(`radius_f64`)이거나 제곱근이 유리수일 때만
`radius_exact()` 다. 「유리수 점을 r 로 짓는 자리 여섯」이 `radius_exact()` 뒤에 서 있다 — 아래 단계의
방문 목록이다.

**남은 문제 — 중심.** 기울어진 코너의 필렛은 **중심 자체가 무리수**다(모서리에서 이등분선 방향으로 r —
단위 이등분선이 무리수). kit 이 «축정렬 직각 코너가 아니면 거절»하는 진짜 이유다. 기울어진 평면이 무리수
법선이라 «점 셋(정의)»을 저장하고 실현하듯, 기울어진 원호도 «만든 방법»을 저장하고 중심을 실현해야 한다.

**설계 방향 — 원의 타입은 하나, 넓히는 것은 «수의 타입».** 「필렛-두-직선」「한쪽만 접함」「두 점 + r²」
같은 **정의 변종의 열거는 답이 아니다** — 구속이 하나 늘 때마다 커널 변종이 하나 는다.

- **다른 CAD.** Parasolid·ACIS·OCCT 는 원을 «중심 f64 + 반지름 f64» 한 타입으로 들고 정확 산술이 없다.
  유도는 커널 밖에서 한다 — 스케치 구속 해석기(SolveSpace·PlaneGCS·Onshape 자체 솔버)나 작도 명령
  (AutoCAD `CIRCLE TTR`·`FILLET`: 오프셋 원 두 개의 교점을 f64 로). «접한다»는 1e-9 안에서만 참이라
  b-rep 이 정점·간선마다 허용오차를 든다 — 이 커널이 거부한 «구성 시점의 tolerance»가 그 대가다. 정확 진영
  (CGAL 원형 커널, LEDA real, CORE Expr)은 «식을 저장하고 필요한 정밀도로 구간 평가, 분리 한계로 0 판정» —
  이 커널이 정점에서 하는 것과 같은 진영이다.
- **직선·원 작도는 제곱근으로 닫혀 있다.** 접선·필렛·아폴로니우스는 전부 «오프셋한 직선/원의 교점»이라
  직선끼리는 사칙연산, 원이 끼면 제곱근 **하나**. 유도된 호 위의 필렛은 근호가 중첩된다(차수 2·4·8…).
  그래도 «자와 컴퍼스로 작도 가능한 수» 한 부류다. 타원·스플라인은 밖이다.
- **저장 타입.** 원은 끝까지 `{ center: [N; 2], r2: N }` 이고 `N` 이 `Rat` 에서 «작도된 수의 정의»로
  넓어진다. 계획: `NumDef { Rat, Add, Sub, Mul, Div, Sqrt }` — 값이 아니라 **정의(식)**, 실현은 `HpBounded`
  구간 캐시(`Vertex` + `PointCache` 와 같은 분리). `a + b√c` 는 깊이 1 특수형이고 필렛 위의 필렛은 깊이 2 —
  두 단계의 저장 타입은 같다. «쪽»(두 교점 중 어느 것)은 정의의 일부다 — 두 중심을 잇는 방향의 좌/우,
  유리수 부호 하나.
- **경우의 수는 라이브러리에 있다.** 두 직선 접 = `meet(offset(l₁,r), offset(l₂,r))` · 직선 접 + 점 통과 =
  `meet(offset(l,r), circle(p,r²), 쪽)` · 두 원 접 = `meet(circle(c₁,(r₁±r)²), circle(c₂,(r₂±r)²), 쪽)` ·
  두 점 + r = `meet(circle(p,r²), circle(q,r²), 쪽)` · 세 점 = 수직이등분선 둘의 교점(근호 없음) · 접점 =
  수선의 발 또는 `c₁ + r₁/(r₁+r)·(x−c₁)`(근호가 새로 안 생긴다). 「필렛」「접선 호」「TTR 원」은 이 두
  연산(오프셋·교점)을 부르는 **함수**이지 커널 변종이 아니다 — «커널의 문은 고리, 펜은 kit» 과 같은 선.
- **두 단계의 실제 차이 = 부호 판정 엔진.** 커널이 수에게 묻는 것은 «이 식의 부호» 하나다. 근호 하나:
  `X + Y√c` 는 제곱해 유리수끼리 비교 — `Vertex::Pierce`/`QuadRoot` 가 쓰는 부호 탑 그대로, 새 엔진 없음.
  중첩: 제곱해 없애면 깊이마다 차수가 두 배라 폭발 ⇒ 구간 평가 + 정밀도 상승(cip 의 에스컬레이션) +
  **멈추는 규칙 = 분리 한계**(식의 차수·계수 크기에서 «0 이 아니면 |x| ≥ 이만큼», BFMSS) — cip 가
  `Undecidable` 로 돌려주는 자리가 «0 임이 증명됨»이 된다. 그 앞에 **구조적 0**: 필렛의 접점은 «수선의
  발»로 정의됐으니 직선 위에 있음은 증명이지 측정이 아니다(접하는 두 원의 판별식 0 도) — 분리 한계까지
  가는 일은 드물다.
- **경계.** 담기는 것은 «순차 작도»(앞서 만든 것에서 다음이 결정됨)로 풀리는 구속이다. 여러 구속을
  **연립**으로만 풀 수 있는 경우는 3차 이상이 나올 수 있고 그건 제곱근 밖 — 이름 붙여 거절하거나, 훨씬
  뒤에 일반 대수적 수(다항식 + 고립 구간)로. 구속 해석기는 이 위의 프런트엔드로, 구속을 «작도 순서»로
  풀어 `NumDef` 식을 내놓는다.
- **남은 계단 둘.** `NumDef` 를 세우되 깊이 1 만 허용하고 기존 부호 탑으로 판정(유리수 직선·원에 접하는
  필렛 전부) → 같은 타입에서 깊이 제한을 풀고 구간 정제 + 분리 한계. 저장 타입은 첫 계단에서 한 번
  정해지고 그 뒤 바뀌지 않는다.
- **열린 비용.** 같은 원을 다른 정의로 두 번 만들었을 때의 동일성 — 정의 해시로는 다르고 대수적 동일성은
  `sign(a−b)=0` 판정이다. 인턴은 정의 기준으로 두고 기하 동일성은 불리언의 일치 판정에 맡기는 것이 평면과
  같은 길이다.
- `arc_turns`/`arc_rat` 의 갈림은 이것과 무관한 편의 설탕(끝점을 90°k 로 유도 vs 받음)이고, 정의 기반이
  되면 둘 다 «정의를 진술하는 한 방법»으로 정리된다.

### pad/pocket 을 kit 의 설탕으로

커널의 `PadOnFace`/`PocketOnFace` 는 «면 프레임 + 압출 + Fuse/Cut» 래퍼이고 새 정확 술어가 없다 —
설탕 판별 기준(«편의 레이어는 커널 op 을 조합만 한다»)에 그대로 걸린다. 커널은 extrude + boolean 만 든다.
kit 이 두 변종을 부르는 곳은 `build.rs` 한 자리다. 제거는 `Operation` 어휘(저장된 로그의 재생)와 kit
재작성을 함께 재야 하므로 **별도 플랜**이다.

## 알려진 결함과 절벽

### 불리언이 이름 붙여 거절하는 인구

- `ObliqueCylinderCut` — 비스듬한 평면 × 원통(타원 교선)은 짓지 않았다.
- `CylinderPairContact` — 옆면이 만나는 원통 쌍(4차 교선)은 짓지 않았다.
- 접선의 선이 셋째 평면 안에 놓이는 배치는 `CoincidentNodes` 로 거절한다.
- 옆면 구멍에 자국이 닿으면 `DeclineKind::CylFaceHole` 로 거절한다.
- 스케치의 임의 각도 호와 호–호 접합(`ArcSweepNotQuarterTurn`·`ArcsMeetAtVertex`).
- 원통이 낀 입력에서는 클래스 reuse(닿을 수 없는 평면 건너뛰기)가 꺼진다 — 밴드 소속 판정의 재설계가 필요하다.

### 발행 설정이 문서의 규칙을 집행하지 않는다

`design.md` 는 `nacre-oracle` 을 「발행하지 않는 dev 전용 크레이트」라고 적지만, 워크스페이스 14개 크레이트 어느 `Cargo.toml` 에도 `publish` 키가 없다 — 전부 기본값 `true` 다. 지금 릴리스를 돌리면 오라클 하네스까지 나간다. 발행 전에 `publish = false` 를 단다.

### `loop_winding` 의 전제가 거짓인 링

사전식 최소 노드가 영역의 극점이 아닌 링(호가 최소를 넘어서는 경우)에 대한 일반 답이 없다. 그런 링이 스위트에 실재한다(`design.md` 「감김과 방향」). 같은 자리의 빚: `frame_sign` 과 옆면 `orient_sign` 의 부호를 잠그는 픽스처(bore 옆면에 구멍이 나는 것)가 없다.

### tessellate 는 첫 거절에서 모델 전체를 멈춘다

`SelfTouchingBoundary` 로 거절되는 접촉 모양(`AtEnd`, 곡선 간선 접촉, 다리 둘 이상)이 남아 있고, 한 면의 거절이 전체 메시를 막는다.

### 혼합 프레임 정점은 reuse 가 답하지 못한다

호출자가 세계 좌표로 명시한 밑캡 위의 프레임-스케치 코너는 세 담체의 모션이 갈려(둘은 프레임, 하나는
세계) 유리수 pullback 이 없다. 인구는 `the_def_road_answers_for_the_populations_it_can_name` 의 ④ 가
핀한다.

같은 부류의 실현 거절 인구: 불리언 결과 정점 4,924 중 `RealizeError::NoMeet` **346**(회전 피연산자가 섞인
결과 — 두 모션 이력·Wide 이름), 피연산자 6,988 중 `NoMeet` **192** · `NoCurvedPoint` **12**(회전된
원통의 seam: 세계에 진술 못 하는 담체) · 반사된 원통(`world_cylinder_def` 는 순수 이동만 안다). 전부 구성
폴백으로 서고 census `r` 행이 센다. 실현 문이 넓어지면 이 수가 준다 — 「모션 사슬을 읽을 때 접는다」가 그
첫 걸음이다.

### 무리수 모션이 낀 datum 은 이름이 없다

interning 이 불가하고(열쇠는 `ThroughKey` — 정점 삼중 + 모션), 동일성은 술어가 매번 판정한다.
**«느릴 뿐 틀리지 않는다»는 미측정 주장이다.** 불리언은 평면 클래스별로 셀 복합체를 만들므로, 같은 평면이
두 핸들이면 «두 클래스»가 된다 — 그것이 정말 비용뿐인지(결과 동일)는 **같은 평면을 두 진술로 넣은
불리언을 하나 만들어 census 로** 재야 안다. 인구: 임의 각 회전 사슬(「모션 사슬을 읽을 때 접는다」의 197).

### 가라앉은 좌표의 경계 0

`Realized::to_f64` 의 오차 팔이 `Mag::of(val).times(Mag::pow2(-53))` 이라 `val == 0.0` 이면 **경계도 0**
이다 ⇒ 0 이 아닌 값에 «참값 ∈ 0 ± 0» 이 실린다(`nearest_f64_big_exact` 의 `exact` 깃발은 `false` 로
정직하지만 경계는 아니다). **인구 0**: 그러려면 좌표가 2⁻¹⁰⁷⁵ 아래여야 하는데 `Rat` 이름은 i128 비율이고
가장 넓게 잰 `Wide` 이름도 168비트다. 고치려면 `to_f64` 의 오차 모형을 손대야 한다(가라앉은 값의 경계는
최소 비정규수의 반). 닿으면 `Unrepresentable` 류로 거절하는 쪽이 맞을 수도 있다 — 그때 잰다.

### 곡면 둘 이상이 만나는 점과 seam 담체

정점 어휘는 `ThreePlane`·`OnSeam`·`Pierce { planes, cylinder, root: QuadRoot }` 다 — 평면 둘 + 원통
하나(최대 2점)까지 «어느 점»을 말한다. **곡면이 둘 이상 낀 점**(원통·원통·평면은 최대 4, 이차곡면 셋은
최대 8)은 말할 변종이 없다. 가지 표기는 결정적이어야 하고, 실제 형상을 만나 변종별로 정한다.

모서리에도 같은 자리가 있다: 원통 **seam** 은 두 면의 교차가 아니라 **한 면의 매개화 이음매**라 두-면 교차
담체로 적히지 않는다. 표기는 자기-인접 `[h_cyl, h_cyl]` 이고 validate 가 «담체 동일 ⇔ 원통» 을 지킨다
(`Violation::EdgeCarrierMismatch`). 이 표기는 잠정이다 — 매개 표현은 원통의 `ref_dir` 과 함께 정한다.

## 보류 — 인구가 생기면

### `Through` 판정 비용: 깊이 2 는 필터가 없다

계기: `nacre-ops/tests/wide_datum_cost.rs`.

- **유리수 닫힘.** 비용은 조각 수가 아니라 **없는 정확 경로의 대가**다. 이름이 `Narrow` 인 datum 으로
  자르면 상승 **278**, `Wide` 면 **887**(mean 256비트, 고갈 0). 두 팔은 서로 다른 정점 삼중을 써야 하므로
  도구 위치가 다르다 — 3.2배는 계수가 아니라 자릿수로 읽는다.
- **무리수 모션.** 차수 **9**(사영 join). 깊이 1: 필터가 산다(200/200 결정, 미결 0/800, 최악 상대 반경
  6.3e-10). **깊이 2: 필터가 없다**(차수 81, 계수가 `f64` 범위 밖 8/8) ⇒ 상승 전용. 필터가 무력한 자리의
  다음 수는 `SurfaceCache` 선실현이고, 깊이 2 가 정확히 그 자리다.

### wide 프레임 실현 비용

`FrameWide` 의 축 실현은 캐시 없이 점마다 돈다(인구가 작아 수용). 스케치→돌출 벽의 정준 이름은
**~115비트에 캡**된다(십진 창이 곱을 묶는다) — `Wide` 이름의 면은 구성 경로에서 안 나오고 첫 생산자는
datum 이다. `n·n` 넘침(1.6%)은 그 경로에서 실재한다(census `wf` 가족).

### 정의만 가리키는 평면의 순회·직렬화

`Model::reachable` 은 면을 통해서만 돈다. 세계 축 평면·datum 이 가리키는 평면·이동본을 정의 경유로도
따라가야 한다. 인구는 실재한다: 씨앗 셋은 심긴 직후 어느 면도 참조하지 않는 orphan 이다(validate 는 도달
집합만 검사해 무위반). 강제하는 것은 없다 — 스냅샷 포맷이 없고, `nacre-step`·`nacre-tess` 는 surface
store 를 돌지 않으며(둘 다 면 경유), 쓰이지 않는 datum 은 어디로도 새지 않는다. **스냅샷 포맷이 생길 때**
일반화한다.

### `check()` 의 비용

`Profile2d::check` 는 진실 위에서 돈다. 볼록 링·17자리 좌표: 34ms@100점 · 3.2s@1000점 · 78s@5000점,
호출당 ~1.7µs(narrow 경로의 gcd 약분이 지배). 손 스케치(수십 점)는 ms 미만이라 수용한다. 수천 점 생성기가
실재해지면 대책 둘: (a) `Rat` 비교 기반 정확 bbox 선별(교차쌍 대부분 기각), (b) 실현 f64 + 건전 오차 한계
필터 → `Rat` 상승(CIP 필터 철학의 2D 판). **둘 다 그 인구가 생기기 전엔 짓지 않는다.**

### 노드 생략의 더 강한 게이트

게이트는 프레임의 `exact()`(축이 정확 유리수 직교로 리프트되는가)다. `PlaneFrame` + `inv_sqrt_exact` 로
«실현이 정확 f64 에 떨어지는가»를 직접 묻는 더 강한 게이트가 가능하지만, **표현식을 바꾸면 노드 인구가
움직인다** — 교체는 **census 관문 동반 필수**다.

### 병렬 효율 회복

reuse 가 추적·셀 패스를 건너뛰면서 남은 일이 169 클래스 중 27 에 몰려 병렬 효율이 4.1배 → 2.3배가 됐다.
정확성이 아니라 스케줄링 문제(균등 인덱스 분할).

### astro-float 의 in-place 연산

할당이 회전 불리언 시간의 36.5%이고 횟수는 포크 없이 못 줄인다(`add`/`mul` 이 새 `BigFloat` 를 반환).
순수 로직에 남은 가장 큰 레버(~1.4배). 포크가 아니라 upstream 기여가 맞는 형태이고, 그 전에
`det3_hp`·`cramer_hp` 의 임시값을 줄이는 값싼 버전이 있다.

### 브라우저용 순수 Rust 할당자

미측정.

### 같은 축 회전의 각도 합성

「모션 사슬을 읽을 때 접는다」에 적힌 대로, 인접쌍 인구가 생기면.

### `validate` 를 연산 직후 자동으로 돌릴 것인가

`nacre-ops` 는 `nacre-validate` 에 dev-의존만 한다 — 연산은 validate 를 부르지 않고 테스트·census 가 부른다. 디버그 빌드에서 강제하려면 의존 방향(validate → topo, ops → validate)과 비용을 먼저 잰다. validate 의 tess 검사(crack-free·출처 정합)도 없다 — 그 검사는 `nacre-tess` 의 자기 테스트에 있다.

### tess 증분 갱신

`tessellate` 는 매번 통째로 짓는다(`Tessellation` 에 무효화 집합이 없다). 면 단위 증분 갱신과 연산별 tolerance override 는 GUI 가 요구할 때.

### interning 표 셋을 하나로

`surface_ids`·`surface_through_ids`·`cylinder_ids` → 표 하나, 열쇠는 사설 깔때기에서 진실이 고른다(`Name(PlaneName, motion) | Verbatim(Surface)`). 새 생산자가 표를 잊어 조용히 중복 핸들을 만드는 실수를 구조로 막고, 구·원뿔이 와도 새 팔이 필요 없다. 전제: `Surface`·`PlanePoints` 에 `Eq, Hash`. 잴 것: 열쇠 크기(약 288 B — 문제면 `Box`). 같은 자리: validate 에 중복 곡면 검사가 없다.

### 판정층의 남은 이름과 타입

`Standard` → `ProofStandard` 개명, 좌표 없이 부호만 답하는 동차 정점(`WorkingVertex`), 고정밀 실현의 타입화(지금은 `(usize, [HpBounded; 3])` 튜플), `FaceInfo`/`WorkingPlane` 합치기(개명이 아니라 설계 작업), `plane_motion(h)` → `m.surface(h).motion()`, `edge_curve(e)` → `edge_cache(e).curve()`.

### 어휘 밖의 요구가 오면

모서리를 축으로 지목하는 회전(챔퍼), 임의 축 회전·임의 평면 미러, 곡면의 반사(`MirrorNotPlanar`), `Through` 원소를 「좌표 값 | 정점 핸들」로 일반화, midplane 정의, 정점을 가리키는 프레임 원점, 원뿔 꼭짓점 등 새 `Vertex` 변종, 구·원뿔·일반 이차곡면쌍.

### undo·체크포인트

append-only 에서 undo 는 연산별 (store 길이, 루트) 체크포인트로 O(1)이고, 로그 중간 편집은 체크포인트부터의 재생 + (연산, 입력) 메모다. 먼저 확인할 위험은 재생 후 핸들 안정성이다: 로그 속 핸들은 인덱스 어휘라 **통째 재생에서만** 성립하고 중간 수정은 하류 인덱스를 밀어낸다(`design.md` 「저장소와 동일성」).

### v2 로 미룬 것

로그 중간 편집을 위한 계보 참조(`OpRef { op, output_slot }`)와 op 로그의 소유자(`Document`), `Store` 스냅샷·직렬화 포맷, 세션 메모리 관리(compact 보다 재구축 우선), 경량 STEP 라이터와 export 시 unseam 옵션, M7 SSI 의 방법 선택.

## 정리 로드맵

거대 파일 분할은 끝났다. 남은 것은 **구현 단순화**이고 **방향은 미결이다** — 아래는 잰 것이지 계획이 아니다. 무엇을 어떻게 합칠지는 의논해서 정한다.

- **같은 일을 하는 길이 여러 갈래인 자리**(제품 호출 / 테스트 호출 수는 그때 다시 센다): `realize_def`·`realize_cache` 와 그 `_tracked`(실제 본체) · `realize_inv_sqrt` → `_rounded` → `_memoized`(각각 호출처 하나인 3단 포장) · `det3` 다섯 벌(`det3`·`_sign`·`_f64`·`_hp`, 그리고 두 크레이트에 있는 `_big`) · `point_in_mixed_ring` 의 두 크레이트 구현(`_inner` 198줄 / geom `_opt` 89줄) · `base_coeffs_rat` 세 정의 · `three_planes` 의 `_rat`·`_big` 사다리 · `boolean` 은 `boolean_with_report` 의 얇은 포장. `nacre-judge` 의 `pub fn coeff_exact`·`coeff_normal_ok` 는 테스트만 부른다(그 보증은 평면이 `exact_coeffs`·`exact_normal` 을 믿을 수 있을 때만 든다는 타입으로 옮겨 갔다). 「f64 는 실현 통로 하나로」를 컴파일 단계에서 강제하는 것이 같은 갈래다.
- **300줄이 넘는 함수 열둘**: `boolean/reconstruct.rs` 의 `reconstruct` 1,113 · `arrangement/trace_plane.rs` 의 `trace_transversal_face` 510 · `cyl_chart/census.rs` 의 `census` 496(테스트 전용) · `cyl_chart/regions.rs` 의 `walk` 485 · `planes/table.rs` 의 `collect_planes` 399 · `transform.rs` 의 `transform_solid` 382 · `boolean/coplanar.rs` 의 `merge_component` 363 · `arrangement/split_circles.rs` 의 `split_circles` 353 · `arrangement/split.rs` 의 `split_at_crossings` 351 · `boolean/grouping.rs` 의 `group_faces` 344 · `ops/datum.rs` 의 `datum_plane` 317 · `planes/cyl_gate.rs` 의 `cylinder_gate` 317. **줄 수는 신호지 규칙이 아니다**(overview 「모듈의 자리」) — 쪼갤 근거는 «안의 한 단계를 다른 호출자가 이름으로 부를 만한가»이고, 그 답을 이미 든 것은 셋이다: `split_at_crossings` 의 `timed!` 구간 넷 · `merge_component` 의 번호 매긴 절 · `reconstruct` 의 `'mat:`·`'faces:` 루프. 나머지 아홉은 긴 것뿐이다.
- **주석이 줄의 33~46%** 다(코드 / 주석: `arrangement` 5,585 / 3,014 · `combinatorics` 3,679 / 2,624 · `nacre-ops/src` 최상위 3,395 / 2,987 · `nacre-scalar` 3,130 / 2,150). `error.rs` 1,193줄의 대부분은 `RejectReason` 변종 doc 이고 그것은 사용자에게 가는 문서다 — 줄일 것은 함수 본문 안의 서사 주석이다.
- `pub(super)` 는 단계의 입구와 테스트가 이름으로 부르는 것에만 달려 있다. `frame3`(26)·`arrangement`(41) 가 가장 많고, 그 가운데 테스트만 부르는 것은 계측을 정리할 때 함께 내려간다.
- `combinatorics/names.rs` 의 `NodeId` 철자 관문(`rg 'NodeId::(ThreePlane|Pierce)'`, 「비어 있어야 한다」)은 히트 여덟을 든다. `arrangement/aliases.rs` 와 `arrangement/split.rs:39` 는 변종을 **패턴으로** 읽고, 단위 테스트 다섯이 직접 짓는다. `split.rs:65` 는 생성이되 의도된 것이다 — 그 쌍은 이미 `{wc, w}` 로 정준이라 그대로 되돌리고, `NodeId::pierce` 를 거치면 `QuadRoot::canonical` 이 쌍 순서를 다시 매겨 root 가 뒤집힐 수 있다. 관문이 그 여덟을 면제하거나, 「비어 있어야 한다」가 무엇을 금지하는지 문장이 다시 말한다.
- `tests/suite/curved_nesting.rs` 의 재연산 census·crossing census doc 은 「14 → 0」 식의 사다리 서사다. 본문의 기대표는 거의 전부 `Ok(n)` 이 됐다. 서사 속의 살아 있는 규칙(섹터는 가장 가까운 rim 노드 사이의 호 run 으로 읽는다 등)을 `design.md` 「원통」 절과 대조해 없는 것만 옮기고, doc 은 현재 표 하나로 줄인다. 같은 부류로 남은 표지: 「Stage 0」·「2b's corpus」·「rule 207」, 그리고 이름이 낡은 테스트 `a_sketch_on_a_prism_side_wall_takes_the_f64_path_today`(지금 단언하는 것은 보고된 `SketchPlane` 의 `exact()` 가 `None` 이라는 것뿐이다).
