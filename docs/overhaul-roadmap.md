# 오버홀 통합 로드맵 (큰 틀 — 세부 셀에 앞선 북극성)

> 이 문서는 **빌드·통합 순서와 브랜치 전략, 회전 아키텍처의 고수준 매핑**을 정한다. **재설계가 아니다** —
> 설계 자체는 `design.md §4·§5·§TIP`·`overhaul.md`에 있고, 여기서는 그것을 **셀로 내려놓는 순서**를 고정한다.
> 이후 모든 세부 셀 계획은 이 로드맵에 정합해야 한다.

## 왜 이 문서

셀 단위로 반응하며 진행하다 큰 흐름 없이 우왕좌왕했다(예: "유리수 좌표를 술어가 읽는다"는 오해 → 정정).
부재했던 것은 **(A) overhaul 브랜치 → main 통합 전략**과 **(B) 회전 아키텍처 고수준 설계**다. 이게 없으면
성급한 회전 구현이 나쁜 방식으로 굳는다. 이 문서가 그 북극성이다.

## 목표 상태 (오버홀 완료 시 — 곡면 없이 M5 + 회전)

**M6(곡면) 진입 금지. M5의 직선·평면 요소 그대로 + 회전.** 완료 시 다음이 된다:

- **입력은 모두 유리수(fraction)** — 치수·각도·좌표를 유리수 타입으로 받는다. 연산 중 **자연스럽게 f64/dd로
  강등**되어 **비트 폭발 차단**. 단 **직선 이동 연속 / 동일축 회전 연속은 번들링**해 강체로 먼저 계산(1회 실현).
- **연산 로그(op-log)** 완비 — 스케치·extrude·Rotate·Translate·Boolean이 로그로 남아 결정적 replay.
- **TIP 적용** — 회전으로 생긴 무리수 좌표를 정의+tol로 들고 f64 필터→astro-float 상승으로 **조용히 안 틀리게 판정**.
- **TNP — 토대만 present, 핵심 난제는 미착수(정직)**. TNP는 **완전히 풀리는 문제 아님**(어떤 CAD도 못 풂).
  우리 방식 = typed provenance(Handle, 문자열 아님) + **위상 변경 시 감지-후-질문**(마술 auto-resolve 아님).
  - **커널에 이미 있음(공짜)**: op-log·replay·Handle 동일성 = provenance 토대. end-state "TNP 어느 정도" = 이것.
  - **제일 어려운 미착수 부분**: **파라메트릭 편집** — 로그 중간 치수를 바꿔 재-replay할 때 위상이 바뀌면 op-log
    참조("이 면에 pad")를 새 엔티티로 **재매핑 or 질문**. 이 편집-시 참조 재해소가 핵심 난제(2D toy엔 없음).
  - **이번 회전 오버홀엔 안 넣음**(회전과 직교). 어려운 편집은 **별도 미래 스레드**(experiment에서 그 난제부터
    밀어본 뒤 커널화). **v1 포함은 별도 결정.**
- **모델링 흐름**:
  - 기본 xyz 평면 중 아무 곳에 **2D 스케치**(스케치 안에 **각도 조건** 가능) → **extrude** → 솔리드.
  - 솔리드에 **회전·이동**을 줄 수 있음.
  - **두 솔리드**(각각 임의 회전·이동)를 **부울** — 우연 공면·횡단·포함(감쌈) 등 **어떤 상황에서도 안정적**.
    판정이 진짜 애매하면 **사용자에게 의도를 물음**(비인터랙티브는 기본 Reject).
  - **면 위 pad/pocket** = 공유 평면 위 2D 스케치+extrude로 **별도 솔리드 만든 뒤 부울**(특이 조건 일괄 커버).
  - 기존 솔리드 표면을 **똑같이/일부** 그리면 그 **평면·선분·점을 공유**(명시 공유 §5).
- **다중 솔리드 출력 명확화** — 부울 결과로 **여러 솔리드**가 나올 수 있고(거절 아님), **각 솔리드는 non-manifold가
  아님이 보장**된다(per-solid validate).

**★ "어떤 상황에서도 안정적"의 정직한 뜻**: 조용히 절대 안 틀린다(sound). 그러나 초월수 상등은 유한 정밀도로
결정 불가 — 진짜 애매하면 자동 해결이 아니라 **declare-0 → ask/reject**. "안정적" = "never silently wrong",
"always auto-resolve"가 아니다.

## 현 상태 (무엇이 어디에)

- **main**: M5 완성 — (5d) exactness sweep 종료, coplanar-contact 셀들. 축정렬 평면 불리언.
- **overhaul 브랜치** (main + (5d) 병합됨):
  - `nacre-scalar` — #1 `Rat`/`Angle`, #2 `orient2d_judge`+astro-float. **아직 아무도 소비 안 함(은행)**.
  - `nacre-ops` surface 공유 — #3. **Case A**(동작 무변경).
  - `docs/overhaul.md`, `experiments/exact2d`(2D 실험 검증분).

## 핵심 원리 (설계 확정 — 재논의 금지)

1. **유리수(Rat)가 사는 곳 (넓게, 정확히)**: 입력 + **유리수-순수 파생** + op-log **정의**. 유리수-순수 파생은
   유리수 `+·−·×` — **번들된 연속 직선 이동(유리수 합)**·**동일축 연속 회전(유리수 각 합)** 포함. **유리수-순수
   chain이면 최종 좌표도 유리수(exact)**. (회전 없이 직선 치수만 나오면 좌표가 유리수로 남는다.)
2. **f64/dd로 강등되는 경계 (딱 이 다섯)**: (a) **부울 등 복잡 연산의 좌표 생성**(비트 폭발 차단), (b) **뉴턴
   반복**, (c) **판정 술어**(f64 필터/dd 상승 = TIP), (d) **비트 임계 초과**(i128 한계), (e) **회전 실현**
   (cos/sin 무리수; **90°계열만 예외로 유리수 exact**). **강등 = 캐시만 f64/dd, 정의는 유리수 그대로**(재계산 가능).
3. **부울·뉴턴·술어는 유리수를 직접 안 먹는다 — f64/dd 이미지를 읽는다**: Constructed = direct(f64 좌표),
   Discovered = indirect(평면 계수), 회전 = 정의+tol로 f64 필터 → astro-float 상승(**TIP**). 즉 "유리수 좌표는
   존재하나, 판정 경로는 그 f64/dd 이미지·계수·정의를 읽는다"(유리수 자체를 술어에 넣지 않음).
4. **exactness의 출처** = "정의가 exact + 좌표/판정을 정의에서 필요 정밀도로 계산·상승"(누적 f64가 진실이 되는
   점 없음).
5. **회전 = tol>0의 유일한 원천**(M5). TIP는 tol>0가 있어야 활성(축정렬은 전부 tol 0이라 상승 미발화).

## (A) 브랜치 → main 전략

**★ 결정: 큰 일괄 병합이 아니라 작은 증분 병합(나눠서 조금씩 안전하게).** `M-rot-core`를 통째로 안 미룬다 —
**각 셀/작은 응집 조각이 안전(Case A·순수 추가·정직 거절)해지는 대로 main에 증분 병합**한다.

- **작동 방식**: 미래 회전 셀은 **작게 green**으로 짓고, 안전하면 곧 main으로 (직접 main 작업 또는 overhaul→main
  소규모 병합). 큰 diff·장기 drift 회피.
- **"안전"의 뜻**(main에 올려도 되는): (a) **Case A**(동작 무변경), 또는 (b) **순수 추가**(새 크레이트/필드, 미소비),
  또는 (c) **정직 거절로 가드된 미완 기능**(예: "회전 표현은 되나 회전 불리언은 아직 Reject" — 조용히 안 틀림).
  세 경우 모두 main의 릴리스 상태를 깨지 않는다(M5의 `Unsupported` 거절과 같은 결).
- **기존 foundation(#1–#3)+실험**: 이미 overhaul에 done. **곧 main으로 catch-up 병합**(순수 추가라 안전) 후, 이후
  회전 작업을 main에서 증분으로. (또는 overhaul 유지·증분 병합 — 어느 쪽이든 "조금씩".)
- **정기 `main→(작업브랜치)` 병합**으로 drift 방지.
- 각 단계((B))는 **여러 작은 셀**이고 각 셀이 병합 후보. 단계 경계가 곧 "이제 이만큼은 main에서 쓸 수 있다"의 지점.

## (B) 전체 시퀀스 (단계 — 각 단계는 여러 셀)

| 단계 | 내용 | 산출 | 판정? |
|---|---|---|---|
| 0 | **foundation** (완료·은행) | Rat/Angle·orient2d_judge·surface 공유 | — |
| 0.4 | **★ 다중 솔리드 (완료·독립)** | ✅ 부울이 `Vec<Solid>` 반환(축정렬 severing) + `DISCONNECTED_RESULT` 은퇴. (5c) 성분 분할 재사용. sever는 전역 validate 그대로 통과. **후속**: 경계-공유(edge-touch)·multi-outer-with-cavity·`Cut(A,A)` empty. design.md §10 기록 | — |
| 0.5 | **★ 3D 실험 (완료·GO)** | ✅ `experiments/exact3d`(격리): H-a orient3d 바운드(위반0·tightness0.121)·H-b 계수-tol·H-c indirect implicit-point(이질적 provenance, wrong0)·H-d 3축 체인 전파·H-e 축변경 번들링·aux(일반 형상 상승0%). #1 위험 해소. FINDINGS=GO. 단계 2가 `Pt3`/`orient3d_judge` 이식 | 실험 |
| 1 | **회전·이동 표현** | **1a ✅ `Transform`+유리수 이동**(뼈대: `Isometry`(nacre-scalar 최초 소비)·`transform_solid` 다중패스 복제·Discovered 정의 보존·`Plane`/`Line` translate·validate/tess/STEP·OCCT diff). **1b ✅ 축정렬 회전**(`Isometry.rotate:Option<Rotation>`·`Axis{X,Y,Z}`·`Store<Rotation>` 공유 노드·`Origin::Rotated{base,rotation}`(tol 저장 안 함·§⑦)·정점만 마킹(엣지=파생캐시)·geom 무변경(성분 회전+생성자 재구성)·90°계열 exact(Constructed·boolean 허용)·비-90° boolean `ROTATED_UNSUPPORTED` 거절(O(1) 가드)·OCCT diff). **1c ✅ 재회전 forest(체인·base=루트)**(회전된 솔리드 재회전 = 모든 회전을 노드로 체인(`Rotation.parent`)·base=비-Rotated 루트·exact 회전도 기록(forest 완전)·`solid_rotation`+균일 불변식 `debug_assert`·**always-chain**(번들 아님·가드 불요·sound)·**관찰 동작 변화 없음**=순수 stage-2 인프라·white-box `forest_probe` 테스트·OCCT diff). **1d ✅ 90°계열 정확 실현**(`Angle::cos_sin_f64`=f64 실현 진실 단일화[quadrantal→정확 0.0/±1.0]·`apply_point`/`apply_dir`이 호출·1c n0 발견 **B0 버그 해소**[boolean→exact 90° 회전이 tol-0 Discovered 정점을 무효화하던 것]·rigid offset-상쇄 잔차 0 보존·90° 정점 정확 축정렬 격자·§⑦ "90°계열 tol 0" 코드 정합; design.md §10). **defer**: **번들링**(같은-축 인접 누적각 1회 실현·인접-가드·stage 2 직전)·강등·임의 유리수 축(Rodrigues) — 단계 2 판정 전까지만 필요. **저순위 최적화(정확성 무관)**: 순수 유리수 체인(유리수 base+90°계열+유리수 피벗)을 f64가 아니라 **유리수 산술로 캐시까지 exact 실현**(원점 밖 피벗 90°도 tol 0). 정확성엔 불필요(상승/공면-여유가 이미 보장)·이득은 축정렬 경계에서 불필요한 정밀-상승 제거 → **stage 3에서 상승 빈도 실측 후 결정** | 없음 |
| 2 | **TIP 코어** (3D 실험 검증분 이식) | **2a ✅ Pt3 tol 계산 이식**(`nacre-scalar::frame3::Pt3`·⑤ 방향별 `[f64;3]` tol 확정·`|R|·old+mix` 누적·`at_with_tol` seed·임의 피벗 exact3d H-f 검증·hp ground truth). **2a-ii ✅ `orient3d_judge` 이식**(Pt3.tol 소비·필터→astro-float 상승→declare-0·`det3_bound`=§TIP② 직접분·H-a 검증 위반0·`Orient` 루트화·`JUDGE_PREC` 튜닝 노브·tol>0 경로[tol-0은 상위층이 nacre-predicates Shewchuk 라우팅]·heavy 테스트 `#[ignore]`·"층만"; design.md §10). **2b ✅ `nacre-tip` forest→Pt3 브리지**(forest **첫 소비자**·회전 정점→Pt3 조립[base=루트·chain 순회]·`Rat::try_from_f64`(f64→exact Rat)·순수-회전 가드[회전-후-이동 defer·§⑦]·`orient3d` 라우팅[tol0=predicates Shewchuk·tol>0=judge]·Discovered→2c·사유별 `TipError`; design.md §1/§10). **2c-i ✅ 간접 orient3d 이식**(`frame3` `Iv` 인터벌 동적 필터[§⑨·술어별 공식 손유도 없음]·`plane_iv`/`plane_hp`[회전 평면 계수 tol·§② 간접분·§539 따름]·`indirect_orient3d_judge`[Cramer `sign(D)·sign(M)`·나눗셈 없음·간접 전용 mag-floor declare-0·부호 관용구 2a-ii 재사용]·**H-b**[계수 tol 위반0·tightness0.410]/**H-c**[이질+near-coplanar 2코퍼스·wrong0/3341·filter/escalation 양경로 실증] `#[ignore]`·exact3d 무변경·"층만"; design.md §10). §TIP② orient3d **직접+간접 소진**. **2c-ii ✅ nacre-tip 간접 브리지**(`orient3d`에 단일-Discovered 디스패치 접기·회전 평면을 f64 계수 아닌 **면 회전 정점 3개**로 재구성[`plane_pts`·winding 무관·법선-불변]·부호 치환[V-슬롯 swap·flip]·다중 implicit/정점부족 honest defer·등가(간접=직접 두 슬롯)·unrotated exact-plane 교차검증·geom dep·"층만"·회전 seam 합성 테스트; design.md §10). **cmp-i ✅ 회전 cmp_coord 술어**(`frame`→`frame2` 리네임 chore 선행·`cramer_iv`/`cramer_hp` 추출로 orient3d와 Cramer 기계 공유[시그니처 유지·bit-identical]·`indirect_cmp_coord_judge`[두 implicit point 축좌표 `sign(Dvec_a[axis]·D_b−Dvec_b[axis]·D_a)·sign(D_a)·sign(D_b)`·나눗셈 없음·nacre-predicates 규약]·winding 불변·`Zero`=같음/declare-0·**H-g** wrong-sign 0/3177 2코퍼스 `#[ignore]`·exact3d엔 없던 새 수학 frame3 직접 검증·"층만"; design.md §10). **§TIP② 완전 소진**(orient3d 직접+간접·cmp_coord). **cmp-ii ✅ nacre-tip cmp_coord 브리지**(`cmp_coord(model,v1,v2,axis)`·두 Discovered seam의 `ThreePlane`→`plane_pts`×2→`indirect_cmp_coord_judge`·winding 무관·`NotSeam`/`PlaneUnderdetermined` defer·**독립 변환 두 큐브** heterogeneous 등가[f64 좌표 wrong-sign 0·resolved>0·swap]·unrotated exact-plane 교차검증·"층만"; design.md §10). **★ TIP 판정층 완성**(stage-3 탐색 확정: boolean 정렬이 single-implicit orient3d[2c-ii]·two-implicit cmp[cmp-ii]뿐·혼합 cmp 없음→**cmp-iii 불필요**). **후속**: stage 3(배선) | 층만 |
| 3 | **회전 불리언 (`M-rot-core`)** | 회전·이동 솔리드 seam·내외 판정이 TIP 소비 → sound(우연공면·횡단·포함). 거절 가드 은퇴. declare-0→기본 Reject. (다중 솔리드는 0.4에서 이미). **← main 병합**. **분해**(design.md §10): **3a** toleranced 술어 라우팅 층(각 평면을 면 회전 정점 세 Pt3로·`rotated` flag로 geom/frame3 분기) — **3a-i ✅** `PlaneInfo.tri_pt3`+`t_orient3d`(회전 불변 등가·ops→tip dep·"층만"; design.md §10). **3a-ii ✅** `t_cmp_coord`(좌표 오라클)+`t_plane_side`(회전 불변·orient3d_judge)·혼합-회전 버그 수정(`plane_def`가 None plane을 tri 좌표로 즉석 빌드·churn 0)·미배선; design.md §10. **3a-iii ✅** `t_plane_pair_dir_sign`(frame3 `dir_sign_judge`=cramer D·orient_sign 규약 다리·**H-i** 법선 near-coplanar wrong0/1325·퇴화 평면 skip)·**3a(라우팅 층) 완성**(4술어); design.md §10. **3b-i ✅** arrange 인덱스 기반 3종 배선(`order_along`/`side_of`→`t_orient3d`·`loop_winding`→`t_cmp_coord`·`dir_sign`/`turn_at`/`point_on_ring`/`every_ray`→`t_plane_pair_dir_sign`)·**per-predicate 파생 라우팅**(`any_rotated`=평면 `tri_pt3.is_some()`·`rotated` 스레딩 폐기)·미사용 geom import·죽은 `tri` 클로저·`tri_pt3` 필드 allow 제거·인트라-닥 링크 전체경로화·declare-0 지연의무 3d 기록(A order_along·B loop_winding lex)·unrotated 무회귀(bit-identical)·가드 유지; design.md §10. **3b-ii ✅** `t_plane_side` 배선(`edge_crosses_face`/`pierced_faces`에 정점 핸들+`model` 스레딩·arrangement 4술어 완성·straddle declare-0→VERTEX_ON_FACE_PLANE 자기보호·`plane_side` import/doc catch-22·unrotated bit-identical·direct 정점만·Discovered 끝점은 3c); design.md §10. **3c (정점-origin 계층·규모 큼·여러 하위셀)** — **3c-i ✅** frame3 `dir_orient3d_judge`(방향 orient3d `sign(det[d,x−base,y−base])`·point_in_solid ray-삼각형이 요구·`p+d` 표현불가라 방향 열로 환원·Iv-직접-열 구성상 sound·`det3_big` 추출·GT 코퍼스 wrong 0/1499·escalated 603·"층만"; design.md §10). **3c-ii ✅** toleranced `ray_triangle_cross_tol` 조립(3 `orient3d_ray`+1 orient3d+1 dir_orient3d·frame3 `orient3d_ray`=`orient3d(base,base+dir,x,y)`=dir_orient3d 인자순열·declare-0→Degenerate·정수좌표 bit-identical + 90° 등변성·"층만"; design.md §10). **3c-iii ✅** `point_in_solid_tol`(회전 내외 분류기·`ray_triangle_cross_tol`로 forward-ray winding·`face_loop_verts`→`vertex_pt3`·퇴화 in-loop `pt3_base_collinear`[base 유리수 orient2d·soundness-critical]·신선 회전 한정·회전-불변 큐브+오목 L-prism·"층만"; design.md §10). **3c-iv ✅** 라우터 `vertex_in_solid`(`solid_is_rotated`||정점 Rotated→`point_in_solid_tol`·else 현행 `point_in_solid`)로 일반 부울 2 호출부(1602 seam-free·1706 overlap_fuse_cut) 배선·dead_code allow 제거·unrotated bit-identical·dispatch 테스트·공면 감지기(3617/4330)는 공면 밀레스톤으로; design.md §10. **3c-v ✅** seam 4-plane 동시성 가드(seam 정점이 4번째 평면 위인 퇴화→FOURPLANE 거절)를 `t_orient3d`로 한 줄 배선·`planes_coplanar` 트윈-스킵 무변경(트윈=같은 Surface 공유→계수 byte-identical→exact·fuzzy는 다른-Surface뿐→공면 밀레스톤)·새 테스트 없음(가드 코드가 매 seam 정점마다 실행→기존 boolean 무회귀가 wiring 검증·회전 건전성은 t_orient3d 기존 tested)·unrotated bit-identical; design.md §10. **3c-vi ✅ (마지막 횡단 f64 자리 종료)** `is_shell_outward`(컴포넌트 outer/cavity 라벨·(5d)#5 extreme-vertex 부호) 회전-exact화: 알고리즘은 회전에도 정확(v\*=min-x 극점·`∃ 인접면 n_x<0`)이라 두 수치 단계만 exact화·**신규 술어 0**. **3c-vi-a** `component_is_outward_tol`(조립-전 `LocalFace`·`model` 무접촉): outward=`(flip?−1:1)·dir_orient3d_judge([1,0,0],tri_pt3)`(3c-i)·lex-min=노드를 `loop_triples` 동형 트리플로→`t_cmp_coord`(3b-i) argmin·`plane_def` 혼합-안전·정직 거절·4테스트·unwired. **3c-vi-b** `assemble_fuse_cut` positives를 컴포넌트별 `any_rotated` 인라인 라우팅(회전→tol·축정렬→f64 무변경)·`any_rotated` pub(crate)·라우터 추출 안 함(단일-site·dispatch가 inverted-routing 못 잡음→3c-v식 무테스트)·비회전 bit-identical·무회귀; design.md §10. **3d-i ✅ (회전 부울 첫 라이브)** `ROTATED_UNSUPPORTED` 진입 가드 은퇴(프로덕션 3줄 삭제)로 완성된 TIP 기계를 회전 입력 라이브 구동. n0 실측 **전부 DNA-safe(silent-wrong 0)**: 횡단(같은-iso) 4종 정확·회전-불변(corner Cut 2.776/Fuse 3.700·**sever 2솔리드**[3c-vi 회전 첫 라이브]·cavity), mixed(한쪽 회전)도 열림·정확, 회전 공면-접촉은 solved-or-honest-reject(Z=감지기 exact 처리·X=general_boolean이 정확히 풀거나 정직 거절, 정점-on-평면 TIP declare-0가 반올림 계수 놓친 공면 포착). 회전-불변 4종+공면 no-silent-wrong 1종 테스트, 기존 거절-단언 3곳→mixed 성공으로 전환·비회전 무회귀; design.md §10. **3d-ii ✅ (OCCT 외부 확증)** 회전 두 피연산자의 부울을 OpenCASCADE와 diff(nacre 회전 STEP export→OCCT bcut/bfuse/bcommon→Σvol/area approx)로 독립 검산. 5 테스트 전부 통과(#[ignore]·로컬 DRAWEXE): overlap Cut(Z→X full-tilt)·Fuse·Common·**sever Cut(full-tilt)→2솔리드=OCCT COMPOUND**(3c-vi 회전 outward 외부 확증)·containment→cavity=BREP_WITH_VOIDS. 회전-불변(자기 일관)을 넘어 성숙 커널과 부피/면적 일치로 **실제 정답** 확증; design.md §10. **3d-iii ✅ (적대적 스트레스)** 회전 부울 대량 굴리기(회전-불변 오라클): 48 케이스(corner/sever/containment/convex × Cut/Fuse/Common × Z43·X67·2축·**3축 Euler 임의방향**) → success 48·**SILENT_WRONG 0**·정직거절 0. DNA no-silent-wrong 경험 확증. **성능 발견**: 회전 부울 ~2s(seam)·~0.5s(seam-free)로 느림(TIP astro-float 상승)→스트레스는 `#[ignore]`(회귀 가드는 3d-i 빠른 테스트)·최적화 후보(interval tol 타이트닝); design.md §10. **성능 ✅ (§TIP⑦ hp 캐싱)** 프로파일로 병목 확정(부울당 판정 ~10k·상승 ~9%·상승당 ~2.4ms=시간 대부분·`plane_def` clone이 정의점 hp를 부울당 수십 번 재계산). `Pt3.hp: Rc<OnceCell>` lazy 공유 memoize로 정의점 astro-float 실현을 부울당 1회로 축소 → **corner 17×·sever 22×·Euler3 37×·stress 151→5.4s(28×)**; escalation 카운트 불변(순수 memoize·정확성 무관)·전 스위트/stress/OCCT green; design.md §10. **후속**: **3d-iv** main 병합·추가 성능(단계 정밀도·plane_def 참조화)·공면-접촉 계층(회전 공면·`planes_coplanar` toleranced). **핵심 지렛대**: arrange가 이미 sign-exact라 `PlaneInfo.tri` 반올림만 멈추면 통째로 sound | 활성 |
| 4 | **확장·정리 (main에서 이어감)** | 스케치 각도 조건, 참조-공유 심화(§5 엣지/정점·6 감지기 은퇴), 인터랙티브 ask-user(§6·앱), export 고정밀 재계산(§4) | — |

*(TNP 파라메트릭 편집은 이 회전 시퀀스에 **없음** — 회전과 직교인 별도 미래 스레드. 위 "명확화" 참조.)*

**★ 순서의 근거**: tol>0을 **먼저 만들고(1)** → 그걸 **sound하게 판정(2)** → **불리언에 활성(3)**. 1 없이 2는
입력이 없고, 2 없이 3은 조용히 틀린다. §TIP ①–⑨가 2–3의 내용을 이미 정의(②는 2D H4 검증; ⑤⑦은 "도입 직전
확정" = 이 단계).

## (C) 회전 아키텍처 (고수준 — 셀이 반드시 지킬 불변)

- **정의가 진리, f64는 캐시.** 회전 정점 = "부모 + 공유 회전노드(angle·유리수 축)". f64 좌표는 파생; export는
  정의에서 고정밀 재계산.
- **공유 회전 forest (§⑦).** 회전을 `Store<Rotation>` 노드로(부모 링크). **per-vertex Angle 복사 금지.** 연쇄
  회전 = 노드 사슬.
- **tol = 방향별 xyz 벡터 (§⑤).** 스칼라 금지. 회전 tol = `|R|·기존 + (|x|+|y|)·da_이산화 + 각도항(v1=0)`
  (H4 교정 — 접선형 반증, 좌표혼합 채택).
- **판정 = 필터 + 상승 (§TIP 7·8).** f64 `|det| > 오차한계`면 확정, 아니면 astro-float, floor 아래면 declare-0.
  **최악-상한(선형 합)**, 확률(RSS) 금지.
- **90°계열은 tol 0** (`try_exact_cos_sin`): 축정렬-급 exact, 상승 불요.
- **상승 층 = astro-float** (H1.5; #44 회피 = `is_positive`).
- **번들링 필수** (H4-amplification): 같은 축 유리수-각 회전은 누적각 1회 실현(증분 금지 — 지수 폭발).

## 명확화·정직한 빚 (구현 전 확정)

- **`Translate`도 필요**: "이동"은 회전과 별개 연산(또는 일반 `Transform` = 회전∘이동). 유리수 이동은 exact
  (유리수-순수), 회전만 무리수.
- **다중 솔리드**: 아래 "다중 솔리드 설계" 절 참조.
- **Transform/Copy 연산 모델**: 아래 "연산 모델 — Transform / Copy" 절 참조.
- **ask-user의 커널/앱 경계**: 커널은 declare-0에서 **Coincidence 필요 신호 + 비인터랙티브 기본 `Reject`**만.
  실제 사용자 질문 UI·site별 해소는 **앱 레이어**(GD&T·스타일처럼 커널 밖). 단계 4/앱.
- **회전 입구 둘**: (1) **솔리드 `Rotate`**(불리언 대상 명확 — 첫 순위), (2) **스케치 내 각도 조건**(각도-정의
  좌표 = `Pt2` 직접 — 단계 4). 둘 다 tol>0을 만들어 TIP를 태운다.
- **면 위 pad/pocket via 부울(목표 item)** = **이미 현 방식**(pad = extrude + Fuse; c-1..c-4). 오버홀은 공유(#3)만
  얹음 — 새로 안 지어도 됨. 회전이 붙으면 그 부울이 TIP를 타면 된다.
- **TNP = 별도 스레드**: 토대(op-log·replay·Handle provenance)는 커널에 이미 있음. 파라메트릭 편집(치수 편집→재
  replay 안정·위상변경 감지)이 남은 부분 — **회전과 직교**, experiment 재검증 후 커널화, v1 포함은 별도 결정.
  (2D H5는 toy — 커널이 더 나가 있음.)
- **non-manifold-free "보장"의 실제 뜻**: 두 매니폴드 솔리드의 부울이 **엣지/정점-만-접촉** 등 non-manifold를 낳을
  수 있다. 커널은 이를 마술로 피하는 게 아니라 **그런 구성을 정직 거절**(honest-reject)해 출력 유효성을 지킨다
  (M5 watertight 불변과 같은 결). 즉 "출력 솔리드는 항상 매니폴드 — 아니면 거절".
- **(5d)는 TIP의 tol=0 기저경로로 보존**: 축정렬 exact 술어가 회전 tol=0 케이스 그대로 = "tol=0이면 지금 술어
  그대로"(§TIP). 회전이 (5d)를 무효화하지 않고 그 위에 tol>0 층을 얹는다.

## 다중 솔리드 설계 (결정 — 독립 초기 셀, 단계 0.4)

**★ 회전과 무관·독립.** 축정렬 M5에서도 필요(현재 severing Cut을 거절 중)하고 (5c)를 재사용하니 작다 →
**회전 전에 독립 셀로 앞당겨** main 병합 가능. 설계는 아래로 확정, 코드 세부는 셀 n0가 조사.

**★ 핵심: 생각보다 작다 — (5c) 성분 분할 기계를 재사용.** 지금도 부울은 결과를 연결 성분으로 나눈다
(`assemble_fuse_cut`의 `face_components`+`is_shell_outward`). 현 로직은 "바깥향 성분 정확히 1개면 수용, 아니면
`DISCONNECTED_RESULT` 거절". 다중 솔리드 = **거절 대신 각 성분을 솔리드로 반환**. 성분 분할은 이미 완성 —
**단, cavity 배정은 outer 1개일 때만 완성**이었다: outward가 여럿이면서 cavity도 있으면 "어느 outer가 어느
cavity를 소유하는가"에 shell-scoped 판정이 필요해 미해결 → `SEVERED_WITH_CAVITY` 정직 거절·defer(아래).

- **`boolean(...) -> Result<Vec<Handle<Solid>>, BoolError>`** ✅ — 결과 솔리드 집합(결정적 순서 = 성분의 기하
  canonical 정렬 `comp_key`; 빈 결과[Cut(A,A)]는 아직 `EmptyResult`·후속).
- **`Solid { outer, cavities }` 불변** ✅ — sever는 각 성분이 cavity-free 솔리드(멀티-outer+cavity는 위 defer).
- **op-log/`OpOutput { solids: Vec }`이 Vec 기록** ✅, replay 결정적 재현.
- **다음 연산은 핸들로 지정** ✅(단일 결과=길이 1; `extrude_and_boolean`은 다중 수용 후 cap 보유 solid 반환).
- **validate**: sever 조각은 서로 다른 엣지 handle이라 **전역 validate가 그대로 통과**(이 셀은 무변경). per-solid
  validate는 **경계-공유(edge-touch) 다중 솔리드**에서만 필요 → §5 명시 공유와 함께 후속(§11-B #8).
- **`DISCONNECTED_RESULT` 은퇴** ✅ → `SEVERED_WITH_CAVITY`(도달 가능·firing 테스트)·`NO_OUTWARD_SHELL`(도달
  불가·방어 백스톱)로 대체.

**철저히, 그러나 giant-refactor 아님**: Store/Handle/op-log 기초는 이미 튼튼 — 다중 솔리드는 자연스러운 일반화
(하나→여럿). 억지 대형 리팩토링은 검증된 (5c) 기계를 흔들 위험. Vec 끝까지·per-solid validate·op-log 정합만
철저히.

## 연산 모델 — Transform / Copy (결정)

- **`Transform(A)` = supersede(제자리 변환).** A 소비, A' 살아남음. 순수 이동/회전. isometry = 유리수 이동
  ∘ Angle 회전(이동=exact, 회전=무리수, 90°계열 exact). **연속 직선 이동·동일축 회전은 번들링**(강체 1회 실현).
- **`Copy(A)` = 일반 복제(로그 노드, 새 핸들).** A 유지 + A_copy 신규. 순수 복사.
- **회전된 복사 = `Copy` 후 `Transform`** — 두 순수 연산의 합성(Transform이 복사 겸하지 않음, 직교).
- **링크 vs 스냅샷 = op-log가 이미 정함. 별도 메커니즘·슈거 없음.** DNA 원칙 3("모델=로그 replay")상 `Copy(A)`는
  로그 노드라, **상류 op 편집→재-replay하면 복사도 재파생 = 링크 동작이 로그에 내재**(공짜). 별도 link-copy 연산
  불요. **동결 스냅샷(편집해도 불변)은 미지원**(이력-없는 detached 기하 = DNA 충돌). **링크 발현은 TNP 편집과
  함께 defer** — 지금은 `Copy`가 그냥 복제.
- **기하 공유로 가볍게 하는 건 skip(payoff 낮음).** 복사의 주 용도(Copy→Transform→Boolean)는 Transform이 즉시
  diverge(변환된 새 요소 append, 공유본 버림)라 공유가 안 남는다. 대신 **로그가 가벼움**(Copy=1 op, 기하는 replay
  materialize). §5 append-only 공유는 "문맥상 같은 것"(면 위 스케치→surface 공유)용이지 통째 복사용 아님. (병목이면
  CoW 후일 고려 — M5선 불필요.)
- **Copy는 M-rot 필수 아님** — 두 솔리드 흐름은 extrude 둘. Copy는 나중 추가.

## ★ 정직한 규모 경고 — 단계 3이 가장 크다

단계 3(회전 불리언)은 "판정에 TIP를 끼우기"가 아니라 **기존 M5 불리언 기계 전체(≈6 공면 감지기 + 일반 seam +
내외 분류)를 회전 하에서 동작**시키는 일이다. 또:

- **대부분의 회전-두-솔리드 부울은 횡단(transverse) → TIP-heavy.** 서로 다른 프레임이라 §3의 "2D-프레임 exact"가
  안 통하고 TIP 상승에 의존한다.
- **§3의 2D-프레임 exact는 *같은 프레임 공유* 접촉에만** — 즉 **회전된 솔리드의 면에 스케치→pad**(공유 평면 프레임)
  가 그 exact 케이스다. 우연 공면·독립 두 솔리드는 TIP.
- 그래서 단계 3은 여러 셀로 쪼갠다: (3a) TIP를 direct/indirect 술어 상승에 배선, (3b) 감지기·seam이 회전 정의를
  읽게, (3c) declare-0 정책. (다중 솔리드·per-solid validate는 **0.4에서 이미** — 회전 불리언은 그 Vec 반환을
  그대로 씀.)

**검증(회전)**: OCCT 오라클을 **회전에도 확장** — 같은 회전·이동을 OCCT에도 적용해 부피·위상 diff. **적대적 회전
코퍼스**(근접-공면 회전, 근접-coincident, grazing) + property test(회전 불변량: 부피·불리언 대수).

## 남은 구현 리스크 (설계는 해소·구현은 셀 n0/실험이 해소)

설계 결정은 위에서 다 됐다. 아래는 **구현·검증에서 풀 리스크**와 **해소 지점**:

1. ~~**3D TIP 수학 미검증(최대)**~~ — **해소(0.5 GO)**: exact3d가 orient3d 바운드·계수-tol·indirect implicit-point·3축·축변경을 sound·타이트(일반 형상 상승0%) 실증. FINDINGS 참조.
2. **회전 표현 shape** — `Origin` 회전 변이·부모 참조(다른 솔리드 핸들)·`Store<Rotation>` 구조·tol 슬롯. `Origin`이
   f64라 `Eq/Hash` 없는 제약. → **단계 1 n0**.
3. **번들링·강등 커널 통합** — 강등 트리거 위치(값별/연산별·i128 임계)·연속 이동/동일축 회전 묶는 방식. → **단계 1 n0**.
4. **다중 솔리드 caller 리플** — `boolean` 반환 `Handle→Vec`가 **모든 호출부**(pad·pocket·extrude_and_boolean·
   테스트)를 건드림 + validate 전역/per-solid 현황. → **0.4 n0**.
5. **TNP 편집 수용 기준** — "뭐가 done인지"(편집→재replay 위상 안정 테스트) 모호. → **TNP 스레드**(defer).
6. **성능** — 3D 회전 부울의 astro-float 상승 빈도 미측정(잦으면 느림). → **단계 3 측정**.
7. **단계 순서 유연성** — 단계 1(회전 표현)은 0.5(3D 실험)와 무관 → 병렬/선행 가능. 0.4·0.5도 서로 독립.

## 열린 결정 (모두 해소)

1. **병합 방식** = **작은 증분 병합**(나눠서 조금씩) — 확정. ((A) 참조.)
2. **회전 입구** = **솔리드 `Rotate`** 먼저 — 확정. (각도-스케치는 단계 4.)
3. TIP ⑤/⑦ = 단계 2; ask-user UI = 앱; **TNP 편집 = 회전과 직교인 별도 스레드**(v1 포함 여부는 별도 결정,
   experiment 재검증 후); 다중 솔리드 = 독립 0.4·per-solid validate — 확정.

---

**즉시 다음**:
1. 이 로드맵(`docs/overhaul-roadmap.md`)을 커밋(북극성 durable화).
2. 그 다음 셀 — **둘 다 회전과 독립·병렬 가능**하니 어느 쪽부터든:
   - **0.4 다중 솔리드**(가장 쉬운 독립 승부, (5c) 재사용, 축정렬 M5에 바로 가치·main 병합), 또는
   - **0.5 3D 실험**(`experiments/exact3d` 격리, #1 TIP 수학 선검증, 커널 무변경).
3. 그 뒤: 단계 1(회전·이동 표현, `Transform`) → 단계 2(TIP 이식, 0.5 검증분) → 단계 3(회전 불리언) → … 각 셀 세부 계획.
