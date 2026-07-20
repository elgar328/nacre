# winding 기반 단일 arrangement 엔진 (M5) — 설계

M5 불리언을 **케이스별 라우팅(detector + 생존표)** 에서 **문헌식 winding 단일 arrangement**(면당 세분 → sub-face를 in/out/on 분류 → op별 keep)로 통합하는 설계 문서. dev-log가 "큰 후속"(1347)으로 명명한 것의 착수 근거·아키텍처·증분 로드맵.

이 문서는 진실이 아니라 **설계 의도 + 실측 근거**다(overview.md 방법론). 각 주장은 코드·측정으로 검증하고, 어긋나면 코드가 진실.

---

## 1. 왜 (결정 + 근거)

**결정(2026-07-20, 사용자):** arrangement/ray-casting 통합 방식은 design.md가 **M7**(하이브리드 메시, design.md:426)에 두지만, M5를 그때까지 **case-routed**로 두면 **경우의 수가 폭발**한다 — 실제로 same_ground·single-shared·containment로 detector·생존표 가지가 계속 늘고 있다. 그 통합 방식을 **M5 평면으로 앞당긴다**(평면은 메시·SSI 불요라 M7 분류법을 exact 평면에 그대로 적용 가능). 케이스별 접근은 이행기 과도물이며, 통합 후 폐기된다.

**니치 갭을 계속 메꾸지 않는 이유:** detector·생존표 배관은 통합 엔진이 오면 삭제된다. 그 케이스별 코드는 새 엔진에 전달되지 않는다. 전달되는 자산 = **테스트(골든·OCCT 오라클, 영구 회귀망)** + **프리미티브(indirect predicate·point_in_solid_idx·plane-triple)** + 현재 robustness.

---

## 2. DNA 호환 (핵심 — mesh-as-truth 문제)

문헌(Cherchi/Zhou)의 "mesh arrangement"는 float 삼각형 수프를 진실로 삼아 overview.md §1("메시를 진실로 만들지 말 것")을 위반한다 — **삼각형이 진실이면.** nacre는 이 긴장을 이미 해소했다: **arrangement가 float 메시가 아니라 exact plane-triple 위에서 돈다.** 모든 arrangement 노드는 `ThreePlane` implicit point, 모든 판정은 exact `indirect_orient3d`/`cmp_coord`. 메시는 design.md:428대로 "임시 분류 도구"일 뿐. **∴ plane-exact arrangement는 DNA 합법, float-mesh-as-truth만 금지.** 회전 입력은 `tolerant` 층이 `Pt3` 재구성으로 건전 처리(design.md:25).

**append-only:** 결과는 새 Face/Shell/Solid를 push하고 `live_solids`를 옮긴다(get_mut 없음). **tolerance는 Discovered seam 정점에만**(Constructed는 exact).

---

## 3. 문헌 계보 + LGPL 위생

계보 = **Cherchi 2020(mesh arrangements)·2022(interactive/robust booleans)·Attene 2020(indirect predicates)·Requicha(boundary evaluation, in/out/on 분류)**. (Zhou 2016·Nef는 nacre 문서에 없음 — 이 계보는 Cherchi/Attene/Requicha.)

**LGPL 위생(우리 커널 라이센스 불변):** 위 참조구현 C++가 LGPL이라, 열람·번역하면 우리 코드가 오염될 수 있다. **→ 설계·구현은 논문(수학은 저작권 대상 아님)·MIT/Apache 소스(`geometry-predicates`)로만**(design.md:414 기존 규칙). 완성 후 실행-대조는 OCCT처럼 허용. **nacre 라이센스 변경 없음.**

---

## 4. 아키텍처

per-op(Fuse/Cut/Common), 두 solid A·B의 **모든 면**에 대해:

1. **면당 기여 수집(통합).** 면 F(A의)에 대한 B의 모든 기여를 한 세분 입력으로:
   - **횡단 seam 곡선** `F∩∂B` — `seam_paths_on`/`edge_crosses_face`(교차점=plane-triple `{F,Q1,Q2}`).
   - **공면 겹침** — F가 B-면과 coincident(canon 일치)면 그 겹침 경계 `coplanar_boundary_crossings`.
   현재 이 둘은 case-routed(`seam_paths_on` xor `coplanar_reconstruct`)로 **분리**돼 있다 — 통합이 1차 작업.
2. **per-face 세분.** 곡선이 F를 sub-face로 나눔. `F∩∂B`는 generically **1-manifold(단순 곡선)** → `boundary_runs`/`run_classes`/`stitch_cycles`(기존 셀 추출기) 재사용. degree-4 내부 자기교차는 드문 퇴화 → 초기 honest-reject(후속 SoS).
3. **sub-face 분류 = propagation-winding.** (F는 항상 A-경계.) B 기준:
   - **in/out(횡단):** 원본 정점을 `point_in_solid_idx`로 anchor(≥3평면라 유효) → seam 곡선 넘을 때 flip 전파(`run_classes` 일반화). **임의-점 술어 신규 불필요.**
   - **on±(공면):** F가 B-면과 coincident인 sub-face는 in/out이 아니라 "on" — 법선 부호로 on+(같은 법선)·on−(반대).
4. **keep + orient.** Requicha 규칙(예 Fuse = A_out_B ⊕ B_out_A ⊕ A_on+B)을 분류 결과에 **균일 적용**(keep/flip 테이블 lib.rs:1976 재사용, 케이스 감지 없이). `loop_winding`/`orient_seam_loop`로 CCW 방향.
5. **조립.** `LocalFace{plane_idx,loop_nodes,inner,flip}` 방출 → **`assemble_fuse_cut` 그대로**(weld-by-triple로 A-조각·B-조각이 같은 seam 정점 공유, cavity·multi-solid 분리 전부 재사용).

---

## 5. 재사용 맵 (신규 술어 0 목표)

| 필요 | 기존 (재사용) | 위치 |
|------|---------------|------|
| 교차점(edge×face=plane-triple) | `edge_crosses_face` | arrange 734 |
| seam 곡선 조립(arc/loop) | `seam_paths_on`·`seam_segments_on` | arrange 1124·112 |
| 직선 위 점 순서 | `order_along` | arrange 215 |
| 평면 위 점-in-polygon | `point_in_ring`/`point_on_ring`/`every_ray` | arrange 642·664·772 |
| 셀 추출(kept/dropped run→cycle) | `boundary_runs`/`run_classes`/`stitch_cycles` | arrange 402·491·268 |
| 3D in/out(winding-parity) anchor | `point_in_solid_idx` | arrange 861 |
| 고리 감김·방향 | `loop_winding`/`turn_at`/`orient_seam_loop` | arrange 1005·540·1051 |
| 공면 겹침 경계 | `coplanar_boundary_crossings` | lib 4102 |
| op별 keep/flip 규칙 | keep 테이블 | lib 1976 |
| **조립(정점 weld·cavity·multi-solid)** | **`assemble_fuse_cut`** | lib 3512 |
| 회전 건전성 | `tolerant` 층(`t_*`) | tolerant.rs |

**제거(통합 후):** kind별 detector들(`detect_*`)·`coplanar_result_unified`/`classify_and_emit` 생존표 분기·`overlap_fuse_cut`/`general_boolean` 이중 경로 + dispatch case-routing.
**유지:** Requicha in/out/on± keep **규칙**(작고 고정, 균일 적용). 즉 "생존표 소거"가 아니라 "**detector 소거 + 규칙 균일화**."

---

## 6. scope 축소 통찰 (2개)

- **글로벌 arrangement 그래프 불필요.** 현 transversal 엔진(`reconstruct_face_paths`)이 이미 **per-face 분해 + weld-by-triple 조립**(글로벌 상태 0, 병렬)으로 정확한 결과를 낸다. winding 엔진도 같은 per-face 패턴 — 조립이 글로벌 봉합. (문헌의 전역 arrangement 구조 불요.)
- **임의-점 in/out 분류기 불필요.** `point_in_solid_idx`는 ≥3평면 점 전용(cell 중심 분류 불가)이나, **propagation**(원본 정점 anchor + seam flip)이 이를 우회 — 새 술어 species 안 만듦.

∴ 이건 **from-scratch 일반 arrangement가 아니라, 기존 per-face 코드(횡단 `reconstruct_face_paths` + 공면 `coplanar_reconstruct`)의 통합·일반화.**

---

## 7. 스파이크 GO/NO-GO 실측 (2026-07-20)

**핵심 가설:** anchor(원본 정점 `point_in_solid_idx`) + propagation이 공면 케이스도 분류할 seed를 주는가?

**측정(same_ground `a=[0,1]³`, `b=[0.5,1.5]²×[0,1]`, a의 8정점을 b에 대해 anchor):**
- **6/8 anchor 성공**(Outside) — off-boundary 정점.
- **2/8 실패**(`NO_CLEAR_RAY`): `(1,1,0)`·`(1,1,1)` = a의 코너가 b의 공유평면(z=0/z=1) 위 + b 발자국 내부 = **"on" 경계 위**.

**해석 = GO(정제됨):** 실패한 2개는 정확히 **on 케이스**(b 경계 위, in/out이 미정의인 게 정상) — 실패가 아니라 on±이 담당할 sub-face. off-boundary 6개가 **propagation seed로 충분.** ∴ "anchor+propagate로 공면 분류 가능" 확인, 그리고 **on-boundary 정점은 anchor 대상이 아니라 on± 대상**이라는 설계 요건을 실측으로 확정.

**함의(설계):** propagation의 anchor는 **B 경계에 안 닿는 원본 정점**만 쓴다(닿는 것은 on±). 어떤 면이 그런 anchor를 하나도 못 가지면(전 정점이 B 경계 위) → 그 면 honest-reject(후속). **완전 de-risk는 M-D 병렬 검증**(전 골든 diff)이 담당.

---

## 8. 정직한 scope / 리스크

**세션-단위 벽돌이 아니라 M5 최대 구조물.** 병렬 개발(현 엔진 보존, 커토버는 M-D 검증 후). 리스크:
- (a) **propagation-winding이 생존표를 재현하는가** — M-A/M-D가 diff로 판별.
- (b) **anchor 가용성** — 어떤 면이 anchor 정점을 못 가지면 honest-reject(§7 실측이 대부분 면은 가짐을 시사).
- (c) **degeneracy** — seam flip 경로가 퇴화 통과 시 미지; propagation이 원본 정점만 anchor라 완화되나 SoS가 필요할 수 있음(design.md "공면은 규칙, 섭동은 우연 퇴화만"). 막히면 honest-reject 유지.
- **비-리스크:** degree-4 내부교차(드문 퇴화 honest-reject), 글로벌 그래프(per-face+assemble로 불요).

---

## 9. 증분 로드맵 (cold-start 금지)

- **M-A** — per-face **propagation-winding 분류** 프리미티브(원본 정점 anchor + seam flip → sub-face in/out-B), isolated + test-only, `run_classes` 재사용. 순수 횡단서 구엔진과 일치.
- **M-B**(완료) — **unified per-face 분류**: 글로벌 detector 없이 로컬 `other_classes.contains(&pi)`로 공면-vs-횡단 dispatch + 공면 on±(`coincident_overlap_faces`+`coplanar_survival`). classification-only(방출=M-C). **실측 발견:** dispatch·on±은 정확하나 same_ground **벽**의 완전 분류는 M-A 횡단엔진의 **공면-인접 프런티어**(on-boundary 정점 `no_clear_ray` + 공유-캡 위 seam `vertex_on_face_plane`) — 부분-anchor 전파 + 공면-aware seam이 필요, 별도 후속.
- **M-C**(완료) — **통합 드라이버 `winding_boolean`**: 글로벌 detector 없이 per-face dispatch로 **기존 작업기 재사용**(횡단→`reconstruct_face_paths`, 조립→`assemble_fuse_cut`; seam 빌드는 `build_seam`로 추출·공유). **첫 end-to-end** — 순수 횡단(corner·tunnel·sever)에서 `winding_boolean == boolean`(부피·manifold·면수·결정성) 검증. 공면 arm·containment는 정직거절(후속). **발견:** `reconstruct_face_paths`의 분류=M-A propagation이라 순수 횡단선 현 엔진=winding(§6대로 재작성 아닌 재사용); M-A/M-B propagation은 프런티어에서 실사용.
- **M-D** — `assemble_fuse_cut` 배선 + **병렬 검증**(전 골든 diff: same_ground·single-shared·transversal…; containment 등 현 실패 케이스도 실측).
- **M-E** — **커토버**: detector·생존표·`coplanar_result_unified`·`overlap_fuse_cut` 은퇴, design.md 최종 갱신.

각 마일스톤: isolated → 병렬 검증 → 배선, 매 단계 측정-먼저·honest-reject·전 게이트(fmt·clippy·스위트·OCCT). 막히면 현 엔진 보존.
