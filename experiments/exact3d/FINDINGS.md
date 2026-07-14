# exact3d — 3D TIP 수학 go/no-go (오버홀 단계 0.5)

**판정: GO.** 회전된 평면 기하의 3D orient3d 술어(직접·간접)가 f64 필터 + astro-float 상승으로 **sound**하게
(부호를 조용히 틀리지 않고) 동작하며, sound 바운드가 **일반 형상에서 상승을 전혀 유발하지 않을 만큼 타이트**하다
(성능도 GO). 리스크 #1(3D TIP 수학 미검증, 로드맵 최대)이 해소됐다 — 단계 2(TIP 코어)가 `Pt3`/`orient3d_judge`/
간접 필터를 nacre-scalar로 이식할 근거 확보.

이 실험은 2D(`experiments/exact2d`, → `nacre-scalar::frame`)가 못 본 **3D-특유**(overhaul.md §11-B)만 다룬다.
파운데이션(유리수·각도·exact 삼각·op-log = 2D H1/H1.5/H5)은 재검증하지 않고 `nacre-scalar`의 `Rat`/`Angle`을
재사용했다. 격리(중첩 빈 `[workspace]`)라 루트 build/test 무영향.

## 가설 표

| H | 검증 대상 | 결과 |
|---|---|---|
| **H-a** | explicit orient3d 오차 한계(`det3_bound`, 새 수학 ②) | 10000 이질적-회전 4점, **위반 0**, worst tightness **0.121** |
| **H-b** | 회전 평면 계수-tol (점 tol → 뺄셈·외적 전파) | 10000, **위반 0**, tightness **0.386** |
| **H-c** | indirect orient3d (3-평면 implicit point, 이질적 provenance + near-coplanar) | GT-stable 3264개 중 **wrong-sign 0**, 상승 1414회(definite-hard 해상 + 퇴화 declare-0) |
| **H-d** | 3축 tol 전파(회전 체인 순회 누적, §⑦) | 5000 체인(길이 2–5, x·y·z 전부), **위반 0**, tightness **0.500** |
| **H-e** | 축변경·번들링 | 90°계열 체인 **tol 0**(exact); 동일축 30스텝 incremental이 bundled 대비 tol **30.8× 증폭**(둘 다 sound → 번들 필수) |
| **aux** | 상승 빈도·속도·경로독립 | 일반 형상 상승 **0%**, near-coplanar **100%**, wrong-sign 0; 필터 **~300ns** vs 상승 **~1.2ms** |

## 핵심 결과

- **soundness(부호를 조용히 안 틀림)**: 모든 H에서 위반/wrong-sign **0**. 예측 tol(선형-합 최악 상한)이 astro-float
  512bit 진리 대비 실제 f64 오차를 한 번도 안 넘겼다. 2D처럼 후보 공식을 stress로 검증했고, 3×3 `det3_bound`와
  간접 필터 모두 첫 후보(좌표혼합·interval)가 곧바로 sound.
- **tightness/성능(진짜 3D 리스크)**: 최악-선형-합은 항상 sound 유도가 가능하므로, 관건은 "sound 바운드가 충분히
  타이트해 상승이 드문가"였다. **일반 형상에서 필터가 100% 해상(상승 0%)** — 상승은 근접-축퇴에서만 발화. 상승은
  느리지만(~4000×) 드물어 amortized 비용은 필터를 따른다. ⇒ 성능 GO.
- **간접 술어(회전 boolean이 실제 소비할 핵심)**: 3-평면 implicit point를 계수-직접 형(division 없이 `sign(D)·
  sign(M)`)으로, f64 **interval 산술 필터**(design.md ⑨ "동적 필터" — 술어별 바운드 유도 없이 구성상 sound)로
  판정하고 astro-float로 상승. **이질적 provenance**(서로 다른 회전 이력이 한 술어에 섞이는, 로드맵 최난 케이스)와
  **near-coplanar**(공유 회전이 coplanarity 보존) 두 코퍼스 모두 wrong-sign 0.
- **declare-0**: 간접 상승은 **크기 floor**(`FLOOR_K·mag·2⁻ᵖʳᵉᶜ`, 항 크기 기반)로 정밀도 이하 값을 `Zero`로 선언
  → §6 "사용자에게 질문". 진짜 축퇴(M=0)는 judge=Zero·GT=degenerate로 일치. definite-but-hard는 judge(200bit)가
  전부 해상(declined 0) — 200bit면 GT가 정하는 definite를 다 정한다.

## 실험 중 발견·교정 (정직)

- **삼각형 퇴화 버그**: near-degenerate 코퍼스 초기 구성이 삼각형 q,r,s를 공선으로 만들어 `cross=0`·`mag_m≈0`이
  돼 M 부호가 정밀도마다 진동(344 wrong). 삼각형을 비퇴화(q−s·r−s 독립)로 고쳐 해소.
- **이질적 provenance vs near-degeneracy는 한 코퍼스에서 결합 불가**: 독립 회전이 base 프레임의 near-coplanar
  구조를 파괴한다(near-degeneracy는 프레임 의존). ⇒ 두 코퍼스로 분리(이질적-generic soundness + 공유회전-near
  escalation). 회전 boolean 이식 시 유의점.
- **base→f64 반올림 tol**: `Pt3::at`가 tol 0으로 두면, 체인이 안 건드린 축이 base 반올림 오차를 tol 없이 남긴다
  (per-axis 검사에서 노출; det/계수 검사에선 숨겨짐). 실제 반올림(정확 base면 0)을 초기 tol로 실었다.
- **`sign_with_floor` 상수**: soundness 우선(선형-합)이라 넉넉히 잡아도 declare-0율만 늘 뿐 틀리지 않는다.

## 파라미터

`JUDGE_PREC=200`(상승 정밀도), GT`=512`(진리; +판정 floor로 진짜 축퇴 분리), `DA_F64=16ε`(2D 검증분 재사용),
`FLOOR_K=1e6`(declare-0 floor 상수, soundness 여유). 이식 시 튜닝 대상은 `JUDGE_PREC`·`FLOOR_K`(성능/declare-0율).

## 술어 커버리지 (완전성)

회전 boolean이 소비할 indirect 술어 매핑: **orient3d** = H-a(explicit)/H-c(indirect); **`plane_side`** =
explicit orient3d 쌍둥이 = H-a; **`indirect_cmp_coord`**(두 implicit point 축 좌표 비교) = H-b(좌표 tol)의
따름정리(뺄셈 비교라 det보다 단순). 전부 커버 또는 corollary.

## 이 실험 밖 (단계 3에서 판정 — §11-B 구조 잔여)

수치 실험은 술어 **개별** soundness/tightness만 본다. 다음은 커널 **구현**이 판정한다:
- 술어-**집합** 상호일관성(boolean 위상 coherence — 2D의 transition oracle 3D판).
- boolean이 공유-토폴로지 피연산자를 받는 3D 구현 형태·경계 공유 재구성(§11-B #6).
- per-solid vs 전역 validate(§11-B #8; 0.4가 edge-disjoint sever는 전역 통과 확인, 경계-공유는 후속).
- declare-0 ask-user 실빈도·Separate 실현율(§11-B #9·10; 실사용 측정).

## 재현

```
cd experiments/exact3d && cargo test            # 7개 전부 green(위반 0)
cd experiments/exact3d && cargo test -- --nocapture   # tightness·상승 빈도·속도 수치
```
루트에서 `cargo test`는 exact3d를 건드리지 않는다(중첩 워크스페이스 격리).
