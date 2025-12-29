#pragma once
typedef ux Run;

#ifndef TOPLEVEL_GC
  #define TOPLEVEL_GC 1
#endif

#if TOPLEVEL_GC
  void run_slow();
  extern GLOBAL ux run_nesting;
  extern GLOBAL ux run_countdown;
  static void nneg_ux(ux x) { assert(x == (x<<1>>1)); }
  
  // must surround all entry points that may ever itself result in doing run_start+run_end
  //   (allocating & inc/dec are okay (if at some point finalizers are added, they'll be deferred to GC or similar, and GC will ensure its own non-zero run_nesting); getters are currently okay)
  static Run run_start() {
    nneg_ux(run_nesting);
    return ++run_nesting;
  }
  static void run_end(Run e) {
    debug_assert(e == run_nesting);
    ux now = --run_nesting;
    nneg_ux(run_countdown-1);
    if (!now) {
      assert(run_countdown > 0);
      if (RARE(!--run_countdown)) run_slow();
    }
  }
  static void run_pressure() {
    run_countdown = 1;
  }
#else
  static Run run_start() { return 0; }
  static void run_end(Run e) { }
  static void run_pressure() { }
#endif
