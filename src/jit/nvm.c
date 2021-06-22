#include "../core.h"
#include "../core/gstack.h"
#include "../ns.h"
#include "../utils/file.h"
#include "../utils/talloc.h"
#include "../utils/mut.h"
#include "nvm.h"

#ifndef USE_PERF
  #define USE_PERF 0 // enable writing symbols to /tmp/perf-<pid>.map
#endif
#ifndef WRITE_ASM
  #define WRITE_ASM 0 // writes on every compilation, overriding the previous; view with:
#endif                // objdump -b binary -m i386 -M x86-64,intel -D --adjust-vma=$(cat asm_off) asm_bin | tail -n+7 | sed "$(cat asm_sed)"
#ifndef CSTACK
  #define CSTACK 1
#endif
#ifdef GS_REALLOC
  #undef CSTACK
  #define CSTACK 0
#endif


// separate memory management system for executable code; isn't garbage-collected
#define  BSZ(X) (1ull<<(X))
#define BSZI(X) ((u8)(64-__builtin_clzl((X)-1ull)))
#define  MMI(X) X
#define   BN(X) mmX_##X
#define buckets mmX_buckets
#include "../opt/mm_buddyTemplate.h"
#define  MMI(X) X
#define  ALSZ  17
#define PROT PROT_READ|PROT_WRITE|PROT_EXEC
#define FLAGS MAP_NORESERVE|MAP_PRIVATE|MAP_ANON|MAP_32BIT
#include "../opt/mm_buddyTemplate.c"
static void* mmX_allocN(usz sz, u8 type) { assert(sz>=16); return mmX_allocL(BSZI(sz), type); }
#undef mmX_buckets
#undef BN
#undef BSZI
#undef BSZ


// all the instructions to be called by the generated code
#if CSTACK
  #define GA1 ,B* cStack
  #define GSP (*--cStack)
  #define GS_UPD { gStack=cStack; }
#else
  #define GA1
  #define GSP (*--gStack)
  #define GS_UPD
#endif
#define P(N) B N=GSP;
#if VM_POS
  #define POS_UPD (envCurr-1)->bcL = bc;
#else
  #define POS_UPD
#endif
#define INS NOINLINE __attribute__ ((aligned(64), hot)) // idk man
INS void i_POPS(B x) {
  dec(x);
}
INS void i_INC(Value* v) {
  ptr_inc(v);
}
INS B i_ADDU(u64 v) {
  return b(v);
}
INS B i_FN1C(B f, u32* bc GA1) { P(x) GS_UPD;POS_UPD; // TODO figure out a way to instead pass an offset in bc, so that shorter `mov`s can be used to pass it
  B r = c1(f, x);
  dec(f); return r;
}
INS B i_FN1O(B f, u32* bc GA1) { P(x) GS_UPD;POS_UPD;
  B r = isNothing(x)? x : c1(f, x);
  dec(f); return r;
}
INS B i_FN2C(B w, u32* bc GA1) { P(f)P(x) GS_UPD;POS_UPD;
  B r = c2(f, w, x);
  dec(f); return r;
}
INS B i_FN2O(B w, u32* bc GA1) { P(f)P(x) GS_UPD;POS_UPD;
  B r;
  if (isNothing(x)) { dec(w); r = x; }
  else r = isNothing(w)? c1(f, x) : c2(f, w, x);
  dec(f);
  return r;
}
INS B i_FN1Ci(B x, BB2B fi, u32* bc GA1) { GS_UPD;POS_UPD;
  return fi(b((u64)0), x);
}
INS B i_FN2Ci(B w, BBB2B fi, u32* bc GA1) { P(x) GS_UPD;POS_UPD;
  return fi(b((u64)0), w, x);
}
INS B i_FN1Oi(B x, BB2B fm, u32* bc GA1) { GS_UPD;POS_UPD;
  B r = isNothing(x)? x : fm(b((u64)0), x);
  return r;
}
INS B i_FN2Oi(B w, BB2B fm, BBB2B fd, u32* bc GA1) { P(x) GS_UPD;POS_UPD;
  if (isNothing(x)) { dec(w); return x; }
  else return isNothing(w)? fm(b((u64)0), x) : fd(b((u64)0), w, x);
}
INS B i_ARR_0() { // TODO combine with ADDI
  return inc(bi_emptyHVec);
}
INS B i_ARR_p(B el0, i64 sz GA1) { assert(sz>0);
  HArr_p r = m_harrUv(sz);
  bool allNum = isNum(el0);
  r.a[sz-1] = el0;
  for (i64 i = 1; i < sz; i++) if (!isNum(r.a[sz-i-1] = GSP)) allNum = false;
  if (allNum) {
    GS_UPD;
    return withFill(r.b, m_f64(0));
  } else return r.b;
}
INS B i_DFND_0(u32* bc, Scope* sc, Block* bl GA1) { GS_UPD;POS_UPD; return m_funBlock(bl, sc); }
INS B i_DFND_1(u32* bc, Scope* sc, Block* bl GA1) { GS_UPD;POS_UPD; return m_md1Block(bl, sc); }
INS B i_DFND_2(u32* bc, Scope* sc, Block* bl GA1) { GS_UPD;POS_UPD; return m_md2Block(bl, sc); }
INS B i_OP1D(B f, u32* bc GA1) { P(m)     GS_UPD;POS_UPD; return m1_d  (m,f  ); }
INS B i_OP2D(B f, u32* bc GA1) { P(m)P(g) GS_UPD;POS_UPD; return m2_d  (m,f,g); }
INS B i_OP2H(B m          GA1) {     P(g)                 return m2_h  (m,  g); }
INS B i_TR2D(B g          GA1) {     P(h)                 return m_atop(  g,h); }
INS B i_TR3D(B f          GA1) { P(g)P(h)                 return m_fork(f,g,h); }
INS B i_TR3O(B f          GA1) { P(g)P(h) B r;
  if (isNothing(f)) { r=m_atop(g,h); dec(f); }
  else              { r=m_fork(f,g,h); }
  return r;
}
INS B i_LOCO(u32 p, Scope* sc, u32* bc GA1) {
  B l = sc->vars[p];
  if(l.u==bi_noVar.u) { POS_UPD; GS_UPD; thrM("Reading variable before its defined"); }
  return inc(l);
}
INS B i_NOVAR(u32* bc GA1) {
  POS_UPD; GS_UPD; thrM("Reading variable before its defined");
}
INS B i_LOCU(u32 p, Scope* sc) {
  B* vars = sc->vars;
  B r = vars[p];
  vars[p] = bi_optOut;
  return r;
}
INS B i_EXTO(u32 p, Scope* sc, u32* bc GA1) {
  B l = sc->ext->vars[p];
  if(l.u==bi_noVar.u) { POS_UPD; GS_UPD; thrM("Reading variable before its defined"); }
  return inc(l);
}
INS B i_EXTU(u32 p, Scope* sc) {
  B* vars = sc->ext->vars;
  B r = vars[p];
  vars[p] = bi_optOut;
  return r;
}
INS B i_SETN(B s, Scope** pscs, u32* bc GA1) {     P(x) GS_UPD; POS_UPD; v_set(pscs, s, x, false); dec(s); return x; }
INS B i_SETU(B s, Scope** pscs, u32* bc GA1) {     P(x) GS_UPD; POS_UPD; v_set(pscs, s, x, true ); dec(s); return x; }
INS B i_SETM(B s, Scope** pscs, u32* bc GA1) { P(f)P(x) GS_UPD; POS_UPD;
  B w = v_get(pscs, s);
  B r = c2(f,w,x); dec(f);
  v_set(pscs, s, r, true); dec(s);
  return r;
}
INS B i_SETNi(B x, Scope* sc, u32 p, u32* bc GA1) { GS_UPD; POS_UPD; v_setI(sc, p, inc(x), false); return x; }
INS B i_SETUi(B x, Scope* sc, u32 p, u32* bc GA1) { GS_UPD; POS_UPD; v_setI(sc, p, inc(x), true ); return x; }
INS B i_SETMi(B f, Scope* sc, u32 p, u32* bc GA1) { P(x) GS_UPD; POS_UPD;
  B w = v_getI(sc, p);
  B r = c2(f,w,x); dec(f);
  v_setI(sc, p, inc(r), true);
  return r;
}
INS B i_FLDO(B ns, u32 p, Scope* sc GA1) { GS_UPD;
  if (!isNsp(ns)) thrM("Trying to read a field from non-namespace");
  B r = inc(ns_getU(ns, sc->body->nsDesc->nameList, p));
  dec(ns);
  return r;
}
INS B i_NSPM(B o, u32 l) {
  B a = mm_alloc(sizeof(FldAlias), t_fldAlias, ftag(OBJ_TAG));
  c(FldAlias,a)->obj = o;
  c(FldAlias,a)->p = l;
  return a;
}
INS B i_CHKV(B x, u32* bc GA1) {
  if(isNothing(x)) { POS_UPD; GS_UPD; thrM("Unexpected Nothing (·)"); }
  return x;
}
INS B i_RETD(Scope* sc GA1) { GS_UPD;
  Body* b = sc->body;
  ptr_inc(sc);
  ptr_inc(b->nsDesc);
  return m_ns(sc, b->nsDesc);
}

#undef INS
#undef P
#undef GSP
#undef GS_UPD
#undef POS_UPD
#undef GA1





#include "x86_64.h"

#if USE_PERF
#include <unistd.h>
#include "../utils/file.h"
FILE* perf_map;
u32 perfid = 0;
#endif

static void* nvm_alloc(u64 sz) {
  // void* r = mmap(NULL, sz, PROT_EXEC|PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANON|MAP_32BIT, -1, 0);
  // if (r==MAP_FAILED) thrM("JIT: Failed to allocate executable memory");
  // return r;
  TmpFile* src = mmX_allocN(fsizeof(TmpFile,a,u8,sz), t_i8arr);
  arr_shVec(tag(src,ARR_TAG), sz);
  return src->a;
}
void nvm_free(u8* ptr) {
  if (!USE_PERF) mmX_free((Value*)RFLD(ptr, TmpFile, a));
}



typedef struct SRef { B v; i32 p; } SRef;
#define SREF(V,P) ((SRef){.v=V,  .p=P})
typedef struct OptRes { u32* bc; u32* offset; B refs; } OptRes;
static OptRes opt(u32* bc0) {
  TSALLOC(SRef, stk, 8);
  TSALLOC(u8, actions, 64); // 1 per instruction; 0: nothing; 1: indicates return; 2: immediate SET; 3: immediate FN1_/FN2C; 4: FN2O; 5: replace with PUSH; 6: decrement 1 data; 10+N: ignore N data
  TSALLOC(u64, data, 64); // variable length; whatever things are needed for the specific action
  u8 rm_map[] = {10,10,10,11,12,6,6,99,99,99,11,12,13,14,15,16,17,18,19};
  #define RM(N) actions[N] = rm_map[actions[N]]
  u32* bc = bc0; usz pos = 0;
  while (true) {
    u32* sbc = bc;
    #define L64 ({ u64 r = bc[0] | ((u64)bc[1])<<32; bc+= 2; r; })
    bool ret = false;
    u8 cact = 0;
    #define S(N,I) SRef N = stk[TSSIZE(stk)-1-(I)];
    switch (*bc++) { case FN1Ci: case FN1Oi: case FN2Ci: case FN2Oi: thrM("JIT optimization: didn't already expect immediate FN__");
      case ADDU: case ADDI: cact = 0; TSADD(stk,SREF(b(L64), pos)); break;
      case LOCM: { u32 d = *bc++; u32 p = *bc++;
        TSADD(stk,SREF(tag((u64)d<<32 | (u32)p, VAR_TAG), pos));
        break;
      }
      case FN1C: case FN1O: { S(f,0)
        if (!isFun(f.v) || v(f.v)->type!=t_funBI) goto defIns;
        RM(f.p); cact = 3;
        TSADD(data, (u64) c(Fun, f.v)->c1);
        goto defIns;
      }
      case FN2C: { S(f,1)
        if (!isFun(f.v) || v(f.v)->type!=t_funBI) goto defIns;
        cact = 3; RM(f.p);
        TSADD(data, (u64) c(Fun, f.v)->c2);
        goto defIns;
      }
      case FN2O: { S(f,1)
        if (!isFun(f.v) || v(f.v)->type!=t_funBI) goto defIns;
        cact = 4; RM(f.p);
        TSADD(data, (u64) c(Fun, f.v)->c1);
        TSADD(data, (u64) c(Fun, f.v)->c2);
        goto defIns;
      }
      case OP1D: { S(f,0) S(m,1)
        if (f.p==-1 | m.p==-1) goto defIns;
        B d = m1_d(inc(m.v), inc(f.v));
        cact = 5; RM(f.p); RM(m.p);
        TSADD(data, d.u);
        TSSIZE(stk)--;
        stk[TSSIZE(stk)-1] = SREF(d, pos);
        break;
      }
      case OP2D: { S(f,0) S(m,1) S(g,2)
        if (f.p==-1 | m.p==-1 | g.p==-1) goto defIns;
        B d = m2_d(inc(m.v), inc(f.v), inc(g.v));
        cact = 5; RM(f.p); RM(m.p); RM(g.p);
        TSADD(data, d.u);
        TSSIZE(stk)-= 2;
        stk[TSSIZE(stk)-1] = SREF(d, pos);
        break;
      }
      case TR2D: { S(g,0) S(h,1)
        if (g.p==-1 | h.p==-1) goto defIns;
        B d = m_atop(inc(g.v), inc(h.v));
        cact = 5; RM(g.p); RM(h.p);
        TSADD(data, d.u);
        TSSIZE(stk)--;
        stk[TSSIZE(stk)-1] = SREF(d, pos);
        break;
      }
      case TR3D: case TR3O: { S(f,0) S(g,1) S(h,2)
        if (f.p==-1 | g.p==-1 | h.p==-1) goto defIns;
        if (isNothing(f.v)) thrM("JIT optimization: didn't expect constant ·");
        B d = m_fork(inc(f.v), inc(g.v), inc(h.v));
        cact = 5; RM(f.p); RM(g.p); RM(h.p);
        TSADD(data, d.u);
        TSSIZE(stk)-= 2;
        stk[TSSIZE(stk)-1] = SREF(d, pos);
        break;
      }
      case SETN: case SETU: case SETM: { S(s,0)
        if (!isVar(s.v)) goto defIns;
        cact = 2; RM(s.p);
        TSADD(data, s.v.u);
        TSSIZE(stk)-= SETM==*sbc? 2 : 1;
        break;
      }
      case RETN: case RETD:
        ret = true;
        cact = 1;
        goto defIns;
      default: defIns:;
        TSSIZE(stk)-= stackConsumed(sbc);
        i32 added = stackAdded(sbc);
        for (i32 i = 0; i < added; i++) TSADD(stk, SREF(bi_optOut, -1))
    }
    #undef S
    TSADD(actions, cact);
    #undef L64
    #undef RM
    if (ret) break;
    bc = nextBC(sbc);
    pos++;
  }
  TSFREE(stk);
  
  TSALLOC(u32, rbc, TSSIZE(actions));
  TSALLOC(u32, roff, TSSIZE(actions));
  B refs = inc(bi_emptyHVec);
  bc = bc0;
  u64 tpos = 0, dpos = 0;
  while (true) {
    u32* sbc = bc;
    u32* ebc = nextBC(sbc);
    #define L64 ({ u64 r = bc[0] | ((u64)bc[1])<<32; bc+= 2; r; })
    u32 ctype = actions[tpos++];
    bool ret = false;
    u32 v = *bc++;
    u64 psz = TSSIZE(rbc);
    #define A64(X) { u64 a64=(X); TSADD(rbc, (u32)a64); TSADD(rbc, a64>>32); }
    switch (ctype) { default: UD;
      case 2: assert(v==SETN|v==SETU|v==SETM);
        TSADD(rbc, v==SETN? SETNi : v==SETU? SETUi : SETMi);
        u64 d = data[dpos++];
        TSADD(rbc, (u16)(d>>32));
        TSADD(rbc, (u32)d);
        break;
      case 3: assert(v==FN1C|v==FN1O|v==FN2C);
        TSADD(rbc, v==FN1C? FN1Ci : v==FN1O? FN1Oi : FN2Ci);
        A64(data[dpos++]);
        break;
      case 4: assert(v==FN2O);
        TSADD(rbc, FN2Oi);
        A64(data[dpos++]);
        A64(data[dpos++]);
        break;
      case 5:;
        u64 on = data[dpos++]; B ob = b(on);
        TSADD(rbc, isVal(ob)? ADDI : ADDU);
        A64(on);
        if (isVal(ob)) refs = vec_add(refs, ob);
        break;
      case 6:
        dec(b(data[dpos++]));
        break;
      case 10:
      case 11:case 12:case 13:case 14:case 15:case 16:case 17:case 18:case 19:
        dpos+= ctype-10;
        break;
      case 1: ret = true; goto def2; // return
      case 0: def2:; // do nothing
        TSADDA(rbc, sbc, ebc-sbc);
    }
    u64 added = TSSIZE(rbc)-psz;
    for (i32 i = 0; i < added; i++) TSADD(roff, sbc-bc0);
    #undef A64
    if (ret) break;
    bc = ebc;
  }
  bc = bc0; pos = 0;
  TSFREE(data);
  TSFREE(actions);
  if (a(refs)->ia==0) { dec(refs); refs=m_f64(0); }
  return (OptRes){.bc = rbc, .offset = roff, .refs = refs};
}
#undef SREF
void freeOpt(OptRes o) {
  TSFREEP(o.bc);
  TSFREEP(o.offset);
}

static u32 readBytes4(u8* d) {
  return d[0] | d[1]<<8 | d[2]<<16 | d[3]<<24;
}

#define ASM_TEST 0 // make -j4 debug&&./BQN&&objdump -b binary -m i386 -M x86-64,intel --insn-width=8 -D --adjust-vma=$(cat asm_off) asm_bin | tail -n+8 | sed "$(cat asm_sed);s/\\t/ /g;s/.*: //"
#if ASM_TEST
  #undef WRITE_ASM
  #define WRITE_ASM 1
  static void write_asm(u8* p, u64 sz);
  static void asm_test() {
    ALLOC_ASM(64);
    for (int i = 0; i < 16; i++) for (int j = 0; j < 16; j++) for (int k = 0; k < 16; k++) BZHI(i,j,k);
    GET_ASM();
    write_asm(bin, ASM_SIZE);
    exit(0);
  }
#endif

#if WRITE_ASM
  static void write_asm(u8* p, u64 sz) {
    i32* rp; B r = m_i32arrv(&rp, sz);
    for (u64 i = 0; i < sz; i++) rp[i] = p[i];
    file_wBytes(m_str32(U"asm_bin"), r); dec(r);
    char off[20]; snprintf(off, 20, "%p", p);
    B o = m_str8l(off);
    file_wChars(m_str32(U"asm_off"), o); dec(o);
    B s = inc(bi_emptyCVec);
    #define F(X) AFMT("s/%p$/%p   # i_" #X "/;", i_##X, i_##X);
    F(POPS) F(INC) F(ADDU) F(FN1C) F(FN1O) F(FN2C) F(FN2O) F(FN1Ci) F(FN2Ci) F(FN1Oi) F(FN2Oi) F(ARR_0) F(ARR_p) F(DFND_0) F(DFND_1) F(DFND_2) F(OP1D) F(OP2D) F(OP2H) F(TR2D) F(TR3D) F(TR3O) F(LOCO) F(LOCU) F(EXTO) F(EXTU) F(SETN) F(SETU) F(SETM) F(FLDO) F(NSPM) F(RETD) F(SETNi) F(SETUi) F(SETMi)
    #undef F
    file_wChars(m_str32(U"asm_sed"), s); dec(s);
  }
#endif

typedef B JITFn(B* cStack, Scope** pscs, Scope* sc);
static inline i32 maxi32(i32 a, i32 b) { return a>b?a:b; }
Nvm_res m_nvm(Body* body) {
  ALLOC_ASM(64);
  TSALLOC(u32, rel, 64);
  #if ASM_TEST
    asm_test();
  #endif
  Reg r_TMP  = 3; // TODO this doesn't really need to be non-volatile
  Reg r_PSCS = 14;
  Reg r_CS   = 15;
  Reg r_SC   = 12;
  PUSH(5);
  PUSH(r_TMP);
  PUSH(r_PSCS);
  PUSH(r_CS);
  PUSH(r_SC);
  MOV(r_CS  , R_A0);
  MOV(r_PSCS, R_A1);
  MOV(r_SC  , R_A2);
  for (i32 i = 1; i < body->maxPSC+1; i++) {
    MOV8rmo(R_A2, R_A2, offsetof(Scope, psc));
    MOV8mro(r_PSCS, R_A2, i*8);
  }
  if ((u64)i_SETN != (u32)(u64)i_SETN) thrM("JIT: Refusing to run with CBQN code outside of the 32-bit address range");
  // #define CCALL(F) { IMM(r_TMP, F); CALL(r_TMP); }
  #define CCALL(F) { TSADD(rel, ASM_SIZE); CALLi(F); }
  u32* origBC = body->bc;
  OptRes optRes = opt(origBC);
  Block** blocks = body->blocks->a;
  i32 depth = 0;
  u32* bc = optRes.bc;
  while (true) {
    u32* s = bc;
    u32* n = nextBC(bc);
    u32* off = origBC + optRes.offset[s-optRes.bc];
    bool ret = false;
    #define L64 ({ u64 r = bc[0] | ((u64)bc[1])<<32; bc+= 2; r; })
    // #define LEA0(O,I,OFF) { MOV(O,I); ADDI(O,OFF); }
    #define LEA0(O,I,OFF,Q) ({ i32 o=(OFF); if(o) LEAi(O,I,o); else if (Q) MOV(O,I); o?O:I; })
    #define SPOS(R,N,Q) LEA0(R, r_CS, maxi32(0, depth+(N)-1)*sizeof(B),Q)
    #if CSTACK
      // #define INV(N,D,F) MOV(R_A##N,r_CS); ADDI(r_CS,(D)*sizeof(B)); CCALL(F)
      #define INV(N,D,F) SPOS(R_A##N, D, 1); CCALL(F)
    #else
      #define INV(N,D,F) CCALL(F) // N - stack argument number; D - expected stack delta; F - called function; TODO instrs which don't need stack (POPS and things with stack delta 1)
    #endif
    #define TOPp MOV(R_A0,R_RES)
    #define TOPs if (depth) { u8 t = SPOS(r_TMP, 0, 0); MOV8mr(t, R_RES); }
    #define LSC(R,D) { if(D) MOV8rmo(R,r_PSCS,D*8); else MOV(R,r_SC); }
    #define INCV(R) INC4mo(R, offsetof(Value,refc)); // ADD4mi(r_TMP, 1); CCALL(i_INC);
    #ifdef __BMI2__ // TODO move to runtime detection maybe
      #define INCB(R,T,U) IMM(T,0xfffffffffffffull);ADD(T,R);IMM(U,0x7fffffffffffeull);CMP(T,U);{JA(lI);MOV1l(U,0x30);BZHI(U,R,U);INCV(U);LBL1(lI);}
    #else
      #define INCB(R,T,U) IMM(T,0xfffffffffffffull);ADD(T,R);IMM(U,0x7fffffffffffeull);CMP(T,U);{JA(lI);IMM(U,0xffffffffffffull);AND(U,R);INCV(U);LBL1(lI);}
    #endif
    switch (*bc++) {
      case POPS: TOPp;
        CCALL(i_POPS);
        if (depth>1) { u8 t = SPOS(r_TMP, -1, 0); MOV8rm(R_RES, t); }
      break;
      case ADDI: TOPs; { u64 x = L64; IMM(R_RES, x); IMM(r_TMP, v(b(x))); INCV(r_TMP); break; } // (u64 v, S)
      case ADDU: TOPs; // (u64 v, S)
        #if CSTACK
          IMM(R_RES, L64);
        #else
          IMM(R_A0, L64); CCALL(i_ADDU);
        #endif
      break;
      case FN1C: TOPp; IMM(R_A1,off); INV(2,0,i_FN1C); break; // (B, u32* bc, S)
      case FN2C: TOPp; IMM(R_A1,off); INV(2,0,i_FN2C); break; // (B, u32* bc, S)
      case FN1O: TOPp; IMM(R_A1,off); INV(2,0,i_FN1O); break; // (B, u32* bc, S)
      case FN2O: TOPp; IMM(R_A1,off); INV(2,0,i_FN2O); break; // (B, u32* bc, S)
      case FN1Ci:TOPp; IMM(R_A1,L64);                 IMM(R_A2,off); INV(3,0,i_FN1Ci); break; // (B, BB2B  fm, u32* bc, S)
      case FN2Ci:TOPp; IMM(R_A1,L64);                 IMM(R_A2,off); INV(3,0,i_FN2Ci); break; // (B, BBB2B fd, u32* bc, S)
      case FN1Oi:TOPp; IMM(R_A1,L64);                 IMM(R_A2,off); INV(3,0,i_FN1Oi); break; // (B, BB2B  fm,           u32* bc, S)
      case FN2Oi:TOPp; IMM(R_A1,L64); IMM(R_A2, L64); IMM(R_A3,off); INV(4,0,i_FN2Oi); break; // (B, BB2B  fm, BBB2B fd, u32* bc, S)
      case ARRM: case ARRO:;
        u32 sz = *bc++;
        if (sz) { TOPp; IMM(R_A1, sz); INV(2,0,i_ARR_p); } // (B, i64 sz, S)
        else    { TOPs;                      CCALL(i_ARR_0); } // (S)
        break;
      case DFND: TOPs; // (u32* bc, Scope* sc, Block* bl, S)
        Block* bl = blocks[*bc++];
        u64 fn = (u64)(bl->ty==0? i_DFND_0 : bl->ty==1? i_DFND_1 : bl->ty==2? i_DFND_2 : NULL);
        if (fn==0) thrM("JIT: Bad DFND argument");
        IMM(R_A0,off); MOV(R_A1,r_SC); IMM(R_A2,bl); INV(3,1,fn);
        break;
      case OP1D: TOPp; IMM(R_A1,off); INV(2,0,i_OP1D); break; // (B, u32* bc, S)
      case OP2D: TOPp; IMM(R_A1,off); INV(2,0,i_OP2D); break; // (B, u32* bc, S)
      case OP2H: TOPp; INV(1,0,i_OP2H); break; // (B, S)
      case TR2D: TOPp; INV(1,0,i_TR2D); break; // (B, S)
      case TR3D: TOPp; INV(1,0,i_TR3D); break; // (B, S)
      case TR3O: TOPp; INV(1,0,i_TR3O); break; // (B, S)
      case LOCM: TOPs; { u64 d=*bc++; u64 p=*bc++; IMM(R_RES, tag((u64)d<<32 | (u32)p, VAR_TAG).u); } break;
      case EXTM: TOPs; { u64 d=*bc++; u64 p=*bc++; IMM(R_RES, tag((u64)d<<32 | (u32)p, EXT_TAG).u); } break;
      case LOCO: TOPs; { u64 d=*bc++; u64 p=*bc++; LSC(R_A1,d);
        MOV8rmo(R_RES,R_A1,p*8+offsetof(Scope,vars)); // read variable
        INCB(R_RES,R_A2,R_A3); // increment refcount if one's needed
        if (d) { IMM(R_A2, bi_noVar.u); CMP(R_A2,R_RES); JNE(lN); IMM(R_A0,off); INV(1,1,i_NOVAR); LBL1(lN); } // check for error
      } break;
      case EXTO: TOPs; { u64 d=*bc++; IMM(R_A0,*bc++); LSC(R_A1,d); IMM(R_A2,off); INV(3,1,i_EXTO); } break; // (u32 p, Scope* sc, u32* bc, S)
      case LOCU: TOPs; { u64 d=*bc++; IMM(R_A0,*bc++); LSC(R_A1,d);                  CCALL(i_LOCU); } break; // (u32 p, Scope* sc, S)
      case EXTU: TOPs; { u64 d=*bc++; IMM(R_A0,*bc++); LSC(R_A1,d);                  CCALL(i_EXTU); } break; // (u32 p, Scope* sc, S)
      case SETN: TOPp; MOV(R_A1,r_PSCS); IMM(R_A2,off); INV(3,0,i_SETN); break; // (B, Scope** pscs, u32* bc, S)
      case SETU: TOPp; MOV(R_A1,r_PSCS); IMM(R_A2,off); INV(3,0,i_SETU); break; // (B, Scope** pscs, u32* bc, S)
      case SETM: TOPp; MOV(R_A1,r_PSCS); IMM(R_A2,off); INV(3,0,i_SETM); break; // (B, Scope** pscs, u32* bc, S)
      case SETNi:TOPp; { u64 d=*bc++; u64 p=*bc++; LSC(R_A1,d); IMM(R_A2,p); IMM(R_A3,off); INV(4,0,i_SETNi); break; } // (B, Scope* sc, u32 p, u32* bc, S)
      case SETUi:TOPp; { u64 d=*bc++; u64 p=*bc++; LSC(R_A1,d); IMM(R_A2,p); IMM(R_A3,off); INV(4,0,i_SETUi); break; } // (B, Scope* sc, u32 p, u32* bc, S)
      case SETMi:TOPp; { u64 d=*bc++; u64 p=*bc++; LSC(R_A1,d); IMM(R_A2,p); IMM(R_A3,off); INV(4,0,i_SETMi); break; } // (B, Scope* sc, u32 p, u32* bc, S)
      case FLDO: TOPp; IMM(R_A1,*bc++); MOV(R_A2,r_SC); INV(3,0,i_FLDO); break; // (B, u32 p, Scope* sc, S)
      case NSPM: TOPp; IMM(R_A1,*bc++); CCALL(i_NSPM); break; // (B, u32 l, S)
      case CHKV: TOPp; IMM(R_A1,off); INV(2,0,i_CHKV); break; // (B, u32* bc, S)
      case RETD: MOV(R_A0,r_SC); INV(1,1,i_RETD); ret=true; break; // (Scope* sc, S); stack diff 0 is wrong, but updating it is useless
      case RETN: IMM(r_TMP, &gStack); MOV8mr(r_TMP, r_CS); ret=true; break;
      default: thrF("JIT: Unsupported bytecode %i", *s);
    }
    #undef INCB
    #undef INCV
    #undef LSC
    #undef TOPs
    #undef TOPp
    #undef INV
    #undef SPOS
    #undef L64
    if (n!=bc) thrM("JIT: Wrong parsing of bytecode");
    depth+= stackDiff(s);
    if (ret) break;
  }
  freeOpt(optRes);
  POP(r_SC);
  POP(r_CS);
  POP(r_PSCS);
  POP(r_TMP);
  POP(5);
  RET();
  #undef CCALL
  GET_ASM();
  u64 sz = ASM_SIZE;
  u8* binEx = nvm_alloc(sz);
  #if USE_PERF
    if (!perf_map) {
      B s = m_str32(U"/tmp/perf-"); AFMT("%l.map", getpid());
      perf_map = file_open(s, "open", "wa");
      print(s); printf(": map\n");
      dec(s);
    }
    u32 bcPos = body->map[0];
    // printf("JIT %d:\n", perfid);
    // vm_printPos(body->comp, bcPos, -1);
    fprintf(perf_map, "%lx %lx JIT %d: BC@%u\n", (u64)binEx, sz, perfid++, bcPos);
  #endif
  memcpy(binEx, bin, sz);
  u64 relAm = TSSIZE(rel);
  // printf("allocated at %p; i_ADDU: %p\n", binEx, i_ADDU);
  for (u64 i = 0; i < relAm; i++) {
    u8* ins = binEx+rel[i];
    u32 o = readBytes4(ins+1);
    u32 n = o-(u32)(u64)ins-5;
    memcpy(ins+1, (u8[]){BYTES4(n)}, 4);
  }
  #if WRITE_ASM
    write_asm(binEx, sz);
  #endif
  FREE_ASM();
  TSFREE(rel);
  return (Nvm_res){.p = binEx, .refs = optRes.refs};
}
B evalJIT(Body* b, Scope* sc, u8* ptr) { // doesn't consume
  u32* bc = b->bc;
  pushEnv(sc, bc);
  gsReserve(b->maxStack);
  Scope* pscs[b->maxPSC+1];
  pscs[0] = sc;
  
  B* sp = gStack;
  B r = ((JITFn*)ptr)(gStack, pscs, sc);
  if (sp!=gStack) thrM("uh oh");
  
  popEnv();
  return r;
}