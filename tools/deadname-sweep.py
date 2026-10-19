#!/usr/bin/env python3
"""Where a required document names an identifier the kernel's code does not have -- a candidate sweep.

It answers one question and nothing else: does a backticked name have **zero non-comment uses under
`crates/`**? What it prints is a list of candidates, not a work list; a person sorts them.

Two filters shape the candidates before that question is asked:
  * STOP, a list of Rust and English words. On the three default documents it hides nothing: every
    stop word they backtick is live under `crates/`.
  * names shorter than three characters are left out. That floor does hide names -- formula
    symbols, the tool name `rg`, the naming prefix `Hp` -- so when chasing a name of two characters
    or fewer, run with it off.
Everything else -- a file name, a dependency, a use in a sibling repository, a use only in comments
-- is flagged, not excluded, and left to the reader: exclusions hide real dead names. A cap on a
backtick span's length drops a long span with every name in it, and a sibling check that matches
loosely lets a longer name there stand in for a shorter dead one.

Usage, from the repository root (needs `rg`):
    python3 tools/deadname-sweep.py [FILE...]
With no arguments it sweeps design.md, overview.md and todo.md -- the gate's line. A source file
works too: the backticked names in its comments are asked the same question.
"""
import io,re,subprocess,glob,os,sys
W=r'[A-Za-z0-9_]'
STOP={'the','and','pub','fn','let','mut','vec','usize','f64','i8','bool','true','false','none',
      'some','self','str','rust','impl','match','for','if','else','type','use','mod','crate'}
SIB=('../step-io','../brep-to-step','../nacre-kit','../nacre-playground/wasm')  # playground: its wasm crate only --
#   web/ holds node_modules, where a name like `Intersection` matches 12 times in @types/three (noise).
# `target/` is build output, not source: the sibling repos' come to gigabytes (a sweep of the three
# documents went past half an hour walking them), and a name met only there is no use of it.
def walk(root):
    for dirpath, dirs, names in os.walk(root):
        dirs[:] = [x for x in dirs if x != 'target']
        yield from (os.path.join(dirpath, n) for n in dirs + names)
files={os.path.basename(p).rsplit('.',1)[0]
       for r in ('crates','tools')+SIB for p in walk(r)}
toml=''.join(io.open(p,encoding='utf-8').read()
             for p in glob.glob('crates/*/Cargo.toml')+['Cargo.toml'])
def uses(name, roots, exts=('*.rs',)):
    """(uses outside comments, uses in comments)"""
    pat=re.compile(r'(?<!'+W+r')'+re.escape(name)+r'(?!'+W+r')')
    # `--no-ignore`/`--hidden`: search what grep -r would, not what the ignore files admit (the
    # repo's `.ignore` re-admits `graft/`). `-F`: a name is an identifier; `pat` checks the boundary.
    args=['rg','-n','--no-heading','--no-ignore','--hidden','-F']+sum([['-g',e] for e in exts],[])+['-g','!target/','--',name]+list(roots)
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
        for tok in re.findall(r'`([^`]+)`', l):        # no cap on the span's length
            for n in re.findall(r'\b([A-Za-z_][A-Za-z0-9_]{2,})\b', tok):
                if n.lower() not in STOP: cand.setdefault(n,i+1)
    rows=[]
    for n,ln in sorted(cand.items(), key=lambda kv: kv[1]):
        live,cmt=uses(n,('crates',))
        if live: continue                               # the one judgement
        sl,_=uses(n,SIB,('*.rs','*.ts'))
        rows.append((ln,n,cmt,n in files,
                     bool(re.search(r'(?<![A-Za-z0-9_-])'+re.escape(n)+r'(?![A-Za-z0-9_-])',toml)),
                     sl))
    return len(cand),rows
for d in sys.argv[1:] or ['docs/design.md','docs/overview.md','docs/todo.md']:
    tot,rows=sweep(d)
    print(f"\n## {d} — backtick identifiers {tot} · **no non-comment use in crates: {len(rows)}**")
    print("| line | name | comment uses in crates | file name? | dependency? | sibling repos |")
    print("|---|---|---|---|---|---|")
    for ln,n,cmt,f,dep,sl in rows:
        print(f"| {ln} | `{n}` | {cmt} | {'✔' if f else ''} | {'✔' if dep else ''} | {sl} |")
