#!/usr/bin/env python3
"""필독 문서가 «커널 코드에 없는» 식별자를 드는 자리 — 후보 스윕.

★ 이 계기는 «판정에서» 제외하지 않는다. 한 물음만 답한다 — **`crates/` 에 비주석 사용이 0인가**.
  ⚠ 후보를 «만들 때»의 필터는 둘 있다(정직하게 적는다): 불용어 STOP 과 **3자 미만 배제**.
    실측(2026-09-14): STOP 낱말은 두 문서에서 전부 crates 에 살아 있어 가리는 것이 0이고,
    길이 하한이 가리는 7개(`C`·`D0`·`D1`·`D2`·`Q1`~`Q3`)는 전부 수식·질문 기호다. 그래도 **필터이므로
    적어 둔다** — 「안 적은 제외」가 이 도구가 경계하는 바로 그것이다.
  나머지(파일 이름·의존·형제 리포·주석만 있음)는 **플래그로 붙여 사람이 가른다.**
  이유: 제외를 넣었더니 «내가 고치려던 결함 둘»을 계기가 지웠다(2026-09-14).
    · 백틱 스팬 길이 상한 → `RotNode`(스팬 112자)를 통째로 건너뜀
    · 형제 리포 검사가 느슨 → `step-io` 주석과 `nacre-kit` 의 `BranchingVertex` 가 `Branch` 를 가림
  ⇒ **분류는 사람이, 계기는 재기만.**
사용: python3 tools/deadname-sweep.py [docs/design.md docs/overview.md]  (리포 루트에서)
"""
import io,re,subprocess,glob,os,sys
W=r'[A-Za-z0-9_]'
STOP={'the','and','pub','fn','let','mut','vec','usize','f64','i8','bool','true','false','none',
      'some','self','str','rust','impl','match','for','if','else','type','use','mod','crate'}
SIB=('../step-io','../nacre-kit','../nacre-playground/wasm')  # ⚠ playground 는 wasm 크레이트만:
#   web/ 은 node_modules 를 들어 `Intersection` 같은 이름이 @types/three 에서 12건으로 잡힌다(잡음).
files={os.path.basename(p).rsplit('.',1)[0]
       for r in ('crates','tools')+SIB for p in glob.glob(r+'/**/*', recursive=True)}
toml=''.join(io.open(p,encoding='utf-8').read()
             for p in glob.glob('crates/*/Cargo.toml')+['Cargo.toml'])
def uses(name, roots, exts=('*.rs',)):
    """(비주석 사용, 주석 사용)"""
    pat=re.compile(r'(?<!'+W+r')'+re.escape(name)+r'(?!'+W+r')')
    args=['grep','-rn']+[f'--include={e}' for e in exts]+[name]+list(roots)
    out=subprocess.run(args,capture_output=True,text=True).stdout
    live=cmt=0
    for x in out.splitlines():
        p=x.split(':',2)
        if len(p)<3 or not pat.search(p[2]): continue
        if p[2].strip().startswith(('//','*','/*')): cmt+=1
        else: live+=1
    return live,cmt
def sweep(doc):
    cand={}
    for i,l in enumerate(io.open(doc,encoding='utf-8').read().split('\n')):
        for tok in re.findall(r'`([^`]+)`', l):        # 길이 상한 없음
            for n in re.findall(r'\b([A-Za-z_][A-Za-z0-9_]{2,})\b', tok):
                if n.lower() not in STOP: cand.setdefault(n,i+1)
    rows=[]
    for n,ln in sorted(cand.items(), key=lambda kv: kv[1]):
        live,cmt=uses(n,('crates',))
        if live: continue                               # ← 유일한 판정
        sl,_=uses(n,SIB,('*.rs','*.ts'))
        rows.append((ln,n,cmt,n in files,
                     bool(re.search(r'(?<![A-Za-z0-9_-])'+re.escape(n)+r'(?![A-Za-z0-9_-])',toml)),
                     sl))
    return len(cand),rows
for d in sys.argv[1:] or ['docs/design.md','docs/overview.md','docs/truth-and-cache.md']:
    tot,rows=sweep(d)
    print(f"\n## {d} — 백틱 식별자 {tot} · **crates 비주석 0건 = {len(rows)}**")
    print("| 줄 | 이름 | crates주석 | 파일명? | 의존? | 형제리포 |")
    print("|---|---|---|---|---|---|")
    for ln,n,cmt,f,dep,sl in rows:
        print(f"| {ln} | `{n}` | {cmt} | {'✔' if f else ''} | {'✔' if dep else ''} | {sl} |")
