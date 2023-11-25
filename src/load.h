#pragma once
extern HArr* comps_curr; // global-ish state, valid during compilation; comps_max elements
enum {
  comps_path,
  comps_args,
  comps_src,
  comps_re, // HArr of re_* values (re_mode & re_scope have unspecified values)
  comps_envPos,
  comps_max
};
enum {
  re_comp, re_compOpts, re_rt, re_glyphs, re_sysNames, re_sysVals, // compiling info
  re_mode, re_scope, // only for repl_exec
  re_max
};
#define COMPS_REF(O,N) O->a[comps_##N]
#define COMPS_CREF(N) COMPS_REF(comps_curr, N)

extern B def_sysNames, def_sysVals;
B comps_getPrimitives(void);
void comps_getSysvals(B* res);

typedef struct Block Block;
typedef struct Scope Scope;
Block* bqn_comp(B str, B path, B args); // consumes all
Block* bqn_compSc(B str, B path, B args, Scope* sc, bool repl); // consumes str,path,args
Block* bqn_compScc(B str, B path, B args, B re, Scope* sc, bool loose, bool noNS); // consumes str,path,args
B rebqn_exec(B str, B path, B args, B re); // consumes str,path,args; runs in a new environment
B repl_exec(B str, B path, B args, B re); // consumes str,path,args; uses re_mode and re_scope
void init_comp(B* new_re, B* prev_re, B prim, B sys); // doesn't consume; writes re_* compiling info into new_re