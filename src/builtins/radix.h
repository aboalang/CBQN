#pragma once

// Radix sorting utilities
// These are leaky macros and assume counts are c0, c1,...
// which must be adjacent in memory. Search functions will also need
// #define GRADE_UD(U,D) U
// to do the appropriate sums for an ascending radix sort.

#define RDX_SUM_1(T)                                  T s0=0;                   for(usz j=0;j<256;j++) { s0+=c0[j];                                  c0[j]=s0;                               }
#define RDX_SUM_2(T)  GRADE_UD(c1[0]=0;,)             T s0=0, s1=0;             for(usz j=0;j<256;j++) { s0+=c0[j]; s1+=c1[j];                       c0[j]=s0; c1[j]=s1;                     }
#define RDX_SUM_4(T)  GRADE_UD(c1[0]=c2[0]=c3[0]=0;,) T s0=0, s1=0, s2=0, s3=0; for(usz j=0;j<256;j++) { s0+=c0[j]; s1+=c1[j]; s2+=c2[j]; s3+=c3[j]; c0[j]=s0; c1[j]=s1; c2[j]=s2; c3[j]=s3; }

#if SINGELI_SIMD
extern void (*const si_scan_pluswrap_u8)(uint8_t* v0,uint8_t* v1,uint64_t v2,uint8_t v3);
extern void (*const si_scan_pluswrap_u32)(uint32_t* v0,uint32_t* v1,uint64_t v2,uint32_t v3);
#define RADIX_SUM_1_u8   si_scan_pluswrap_u8 (c0,c0,  256,0);
#define RADIX_SUM_1_u32  si_scan_pluswrap_u32(c0,c0,  256,0);
#define RADIX_SUM_2_u8   si_scan_pluswrap_u8 (c0,c0,2*256,0);
#define RADIX_SUM_2_u32  si_scan_pluswrap_u32(c0,c0,2*256,0);
#define RADIX_SUM_4_u8   si_scan_pluswrap_u8 (c0,c0,4*256,0);
#define RADIX_SUM_4_u32  si_scan_pluswrap_u32(c0,c0,4*256,0);
#else
#define RADIX_SUM_1_u8   RDX_SUM_1(u8)
#define RADIX_SUM_1_u32  RDX_SUM_1(u32)
#define RADIX_SUM_2_u8   RDX_SUM_2(u8)
#define RADIX_SUM_2_u32  RDX_SUM_2(u32)
#define RADIX_SUM_4_u8   RDX_SUM_4(u8)
#define RADIX_SUM_4_u32  RDX_SUM_4(u32)
#endif

#if SINGELI_SIMD && !USZ_64
#define RADIX_SUM_1_usz  si_scan_pluswrap_u32(c0,c0,  256,0);
#define RADIX_SUM_2_usz  si_scan_pluswrap_u32(c0,c0,2*256,0);
#define RADIX_SUM_4_usz  si_scan_pluswrap_u32(c0,c0,4*256,0);
#else
#define RADIX_SUM_1_usz  RDX_SUM_1(usz)
#define RADIX_SUM_2_usz  RDX_SUM_2(usz)
#define RADIX_SUM_4_usz  RDX_SUM_4(usz)
#endif

u8 radix_offsets_2_u32(usz* c0, u32* v0, usz n); // selfsearch.c
